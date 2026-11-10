# Process every listed participant whose files are present. Mirrors MATLAB
# Scripts #1–#5 (probe part), #6b, #6c, #7b and #7c (averaging). Saves:
#   results/probe_trials.csv  one row per probe: labels, next event, artifact flags
#   results/probe_erps.csv    participant × bin × signal averaged waveforms (no low-pass),
#                             with the number of accepted trials of each bin
#   results/probe_qc.csv      trials and % rejected per participant, set and detector
#
# Run:  julia --project=. scripts/01_process_probe.jl

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "..", "src", "VWMPipeline.jl"))
include(joinpath(@__DIR__, "..", "src", "ProbeERP.jl"))
using .VWMPipeline, .ProbeERP
using CSV, DataFrames, Statistics, DSP

eeg_path(sub)    = joinpath(BIDS_ROOT, sub, "eeg", "$(sub)_task-$(TASK)_eeg.vhdr")
bin_path(sub)    = joinpath(BIDS_ROOT, sub, "eeg", "$(sub)_task-$(TASK)_eeg.eeg")
events_path(sub) = joinpath(BIDS_ROOT, sub, "eeg", "$(sub)_task-$(TASK)_events.tsv")
downloaded(d) = isfile(eeg_path(d)) && isfile(bin_path(d)) &&
                filesize(bin_path(d)) > 1_000_000 && isfile(events_path(d))
subjects = filter(downloaded, SUBJECTS)
isempty(subjects) && error("No EEG files found under $BIDS_ROOT. See README: Get the data.")
# Optional: reprocess only some participants and update their rows in the existing
# results files, e.g.  ONLY=sub-04,sub-10,sub-49 julia --project=. scripts/01_process_probe.jl
const ONLY = filter(!isempty, strip.(split(get(ENV, "ONLY", ""), ',')))
if !isempty(ONLY)
    subjects = filter(in(String.(ONLY)), subjects)
    isempty(subjects) && error("None of ONLY=$(join(ONLY, ",")) found under $BIDS_ROOT")
    @info "Reprocessing only $(join(subjects, ", ")); the other participants are kept from results/"
end
missing_subs = setdiff(SUBJECTS, subjects)
isempty(missing_subs) || !isempty(ONLY) ||
    @warn "$(length(missing_subs)) of $(length(SUBJECTS)) participants not downloaded yet: $(join(missing_subs, ", "))"

BINS = read_bdf(BDF_PROBE)                       # bin => (label, codes)
params = CSV.read(ARTIFACT_PARAMETERS, DataFrame)
const DETECTORS = Dict("step" => art_step, "extreme" => art_extreme, "mwpp" => art_mwpp)
left_ch  = [l for (l, _) in CLUSTER]
right_ch = [r for (_, r) in CLUSTER]
# Bins averaged per set: N2/FN400 needs 1, 2, 5, 6 and 38–45 (Script #7c);
# the probe N2pc needs 46–61 (Script #7b). All 61 are counted for the QC table.
SET_BINS = Dict("n2_fn400" => [1, 2, 5, 6, 38, 39, 40, 41, 42, 43, 44, 45],
                "n2pc" => collect(46:61))

mkpath(RESULTS)
trials_all, erps_all, qc_all, counts_all = DataFrame(), DataFrame(), DataFrame(), DataFrame()

for sub in subjects
    id = parse(Int, sub[5:end])
    @info "Processing $sub"

    # Script #1: load, downsample to 256 Hz, bipolar EOG, DC removal, high-pass.
    data, fs = load_eeg(eeg_path(sub))
    fs_orig = fs
    if fs != TARGET_FS
        r = rationalize(TARGET_FS / fs)
        for ch in collect(keys(data))
            data[ch] = numerator(r) == 1 ? eeglab_resample(data[ch], denominator(r)) :
                                           resample(Float64.(data[ch]), r)
        end
        fs = TARGET_FS
    end
    # The ICA was computed on data referenced to the left mastoid (Script #1).
    lm = copy(data["Mastoid-left"])
    for (_, x) in data
        x .-= lm
    end
    add_bipolars!(data)                 # uncorrected HEOG/VEOG (channels 70, 71)
    remove_dc!(data)
    butter_filter!(data, fs; kind=:highpass, cutoff=HIGHPASS.cutoff, order=HIGHPASS.order)

    events = read_events(events_path(sub))
    if fs_orig != fs
        # EEGLAB's pop_resample keeps fractional latencies ((lat − 1)·256/512 + 1, e.g.
        # 101.5); how pop_epoch turns them into a sample is set by EVENT_ROUNDING
        # (scripts/config.jl; test with EVENT_ROUNDING=floor|ties_away|nearest|ceil).
        events.sample = round.(Int, (events.sample .- 1) .* (fs / fs_orig), EVENT_ROUNDING) .+ 1
    end

    # Scripts #2–#4: original ICA weights, then average-mastoid reference.
    rm_ics = get(ICS_APPLIED, id, ics_to_remove(ICS_TO_REMOVE, id))
    ica = read_ica_weights(find_ica_set(ICA_WEIGHTS_DIR, id), rm_ics)
    apply_ica_weights!(data, ica)
    rereference!(data)
    data["HEOG_corr"] = data["HEOG-left"] .- data["HEOG-right"]     # channel 68
    for ch in setdiff(keys(data), KEEP_CHANNELS)
        delete!(data, ch)
    end
    GC.gc()

    # Script #5: probe events, bins and epochs (−200 to 800 ms, baseline −200–0 ms).
    codes = event_codes(events)
    id in SWAP_RESPONSE_IDS && swap_response_codes!(codes)
    trials = label_probe_trials(codes, events.sample, PROBE_CODES)
    onsets = trials.probe_sample
    ep = Dict(ch => first(epoch(data[ch], onsets, fs; tmin=EPOCH.tmin, tmax=EPOCH.tmax,
                                baseline=BASELINE)) for ch in KEEP_CHANNELS)
    times = last(epoch(data["HEOG"], onsets[1:1], fs; tmin=EPOCH.tmin, tmax=EPOCH.tmax))
    times_ms = 1000 .* times

    # Scripts #6b/#6c: artifact detection with the participant's parameters.
    binned = trials.correct
    qc = Dict{Symbol,Any}(:subject => sub, :n_probes => nrow(trials), :n_binned => count(binned),
                          :n_correct_lenient => count(trials.correct_lenient))
    keep = Dict{String,BitVector}()
    for set in ARTIFACT_SETS
        any_art = falses(nrow(trials))
        P = params[(params.id .== id) .& (params.set .== set), :]
        nrow(P) == 0 && error("no artifact parameters for participant $id, set $set")
        for row in eachrow(P)
            det = DETECTORS[row.detector]
            flagged = falses(nrow(trials))
            for ch in split(row.channels, ';')
                flagged .|= flag_trials(ep[String(ch)], times, fs, det;
                                        period=(row.start_ms / 1000, row.end_ms / 1000),
                                        threshold=row.threshold, window_ms=row.window_ms,
                                        step_ms=row.step_ms)
            end
            trials[!, Symbol("art_", set, "_", row.name)] = flagged
            qc[Symbol("pct_", set, "_", row.name)] = round(100 * mean(flagged[binned]); digits=1)
            any_art .|= flagged
        end
        keep[set] = binned .& .!any_art
        trials[!, Symbol("keep_", set)] = keep[set]
        qc[Symbol("n_kept_", set)] = count(keep[set])
        qc[Symbol("pct_rejected_", set)] = round(100 * (1 - count(keep[set]) / count(binned)); digits=1)
    end
    push!(qc_all, NamedTuple(qc); cols=:union)

    # Accepted trials in each of the 61 bins (for the comparison with ERPLAB).
    for set in ARTIFACT_SETS, b in sort(collect(keys(BINS)))
        inbin = [c in BINS[b][2] for c in trials.code]
        push!(counts_all, (subject=sub, set=set, bin=b, n_accepted=count(inbin .& keep[set])))
    end

    # Script #7b/#7c: averages per bin (pop_averager), no low-pass.
    sig_n2   = (ep["Fz"] .+ ep["FCz"]) ./ 2          # 2ChanAverageN2(Fz,FCz)
    sig_fn4  = (ep["Fz"] .+ ep["Cz"]) ./ 2           # 2ChanAverageFN400(Fz,Cz)
    sig_left  = sum(ep[c] for c in left_ch) ./ length(left_ch)    # P3 … O1
    sig_right = sum(ep[c] for c in right_ch) ./ length(right_ch)  # P4 … O2
    signals = Dict("n2_fn400" => ["N2" => sig_n2, "FN400" => sig_fn4],
                   "n2pc" => ["left" => sig_left, "right" => sig_right])
    for set in ARTIFACT_SETS, b in SET_BINS[set]
        sel = keep[set] .& [c in BINS[b][2] for c in trials.code]
        for (name, M) in signals[set]
            avg = average_bin(M, sel)
            append!(erps_all, DataFrame(subject=sub, set=set, bin=b, signal=name, n_trials=avg.n,
                                        time_ms=times_ms, amplitude=avg.wave))
        end
    end

    trials.subject = fill(sub, nrow(trials))
    append!(trials_all, select(trials, :subject, Not(:subject)); cols=:union)
end

# With ONLY: replace these participants' rows in the existing results files.
if !isempty(ONLY)
    function merge_old(name, new)
        f = joinpath(RESULTS, name)
        isfile(f) || return new
        old = CSV.read(f, DataFrame)
        old.subject = String.(old.subject)
        return vcat(old[.!in.(old.subject, Ref(Set(subjects))), :], new; cols=:union)
    end
    trials_all = merge_old("probe_trials.csv", trials_all)
    erps_all   = merge_old("probe_erps.csv", erps_all)
    qc_all     = merge_old("probe_qc.csv", qc_all)
    counts_all = merge_old("probe_bin_counts.csv", counts_all)
    sort!(qc_all, :subject)
end
CSV.write(joinpath(RESULTS, "probe_trials.csv"), trials_all)
CSV.write(joinpath(RESULTS, "probe_erps.csv"), erps_all)
CSV.write(joinpath(RESULTS, "probe_qc.csv"), qc_all)
CSV.write(joinpath(RESULTS, "probe_bin_counts.csv"), counts_all)

# Trial counts vs ERPLAB (sum over participants, N2/FN400 artifact set, bins 1–45).
ref = joinpath(REFERENCE, "matlab_trial_counts_probe.csv")
if isfile(ref) && length(unique(counts_all.subject)) == length(ALL_IDS)   # reference = sum over all 53
    m = CSV.read(ref, DataFrame)
    j = combine(groupby(counts_all[counts_all.set .== "n2_fn400", :], :bin), :n_accepted => sum => :julia)
    c = innerjoin(m[m.bin .<= 45, :], j, on=:bin)
    c.difference = c.julia .- c.accepted_total_53
    CSV.write(joinpath(RESULTS, "trial_counts_comparison.csv"), c)
    @info "Accepted trials, bins 1–45 (N2/FN400): Julia − ERPLAB, total over bins = $(sum(c.difference)); " *
          "bin 17 (all probes): Julia $(c.julia[c.bin .== 17][1]), ERPLAB $(c.accepted_total_53[c.bin .== 17][1])"
end
@info "Done: $(length(unique(qc_all.subject))) participants; kept $(sum(qc_all.n_kept_n2_fn400)) probe trials " *
      "for N2/FN400 and $(sum(qc_all.n_kept_n2pc)) for the probe N2pc."
