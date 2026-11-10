"""
    VWMPipeline

Readable Julia pipeline for the OpenNeuro dataset ds006979
(Ramzaoui & Beck, perceptual grouping in visual working memory).

It follows the original MATLAB/EEGLAB/ERPLAB scripts where possible:

| MATLAB script | Step | Here |
| --- | --- | --- |
| 0 | BDF import, channel selection, CleanLine, practice removal | already applied to the BIDS files |
| 1 | Re-reference, bipolar HEOG/VEOG, DC removal, 0.1 Hz high-pass (order 2) | `rereference!`, `add_bipolars!`, `remove_dc!`, `butter_filter!` |
| 2–4 | ICA correction and channel interpolation | not reproduced (see README) |
| 5 | Bins time-locked to the memory array (51–58), correct trials; epochs −200 to 1000 ms | `label_trials`, `epoch` |
| 6a | Blink/saccade step functions, extreme values, moving-window peak-to-peak | `art_step`, `art_extreme`, `art_mwpp` |
| 7a | Contra/ipsi waves averaged over the posterior cluster; 20 Hz low-pass (order 8) | `cluster_contra_minus_ipsi`, `butter_filter!` |
"""
module VWMPipeline

using EDF, DSP, CSV, DataFrames, Statistics, Mmap, MultivariateStats, LinearAlgebra, Random, MAT

export load_eeg, load_brainvision, read_events, read_vmrk, event_codes,
       rereference!, add_bipolars!, remove_dc!, butter_filter!,
       label_trials, swap_response_codes!, epoch, lowpass_epochs!,
       art_step, art_extreme, art_mwpp, flag_trials,
       contra_minus_ipsi, cluster_contra_minus_ipsi, mean_amplitude,
       alpha_power, percent_change, ica_blink_correct!, break_mask,
       find_ica_set, read_ica_weights, ics_to_remove, apply_ica_weights!,
       eegfiltnew_lowpass, band_power, erplab_window, eeglab_resample, kaiser_bandpass, matlab_filtfilt, sos_matrix,
       ARRAY_CODES, PROBE_CODES

# ---------------------------------------------------------------------------
# Reading data
# ---------------------------------------------------------------------------

"""
    load_eeg(path; channels=nothing) -> (data, fs)

Read EEG from an EDF/BDF file (with EDF.jl) or a BrainVision `.vhdr` file.
Returns a `Dict` mapping channel label to signal in µV, and the sampling rate
in Hz. Pass `channels` (a list of labels) to load only those channels.
"""
function load_eeg(path::AbstractString; channels=nothing)
    endswith(lowercase(path), ".vhdr") && return load_brainvision(path; channels)
    edf = EDF.read(path)
    data = Dict{String,Vector{Float64}}()
    fs = NaN
    for sig in edf.signals
        sig isa EDF.Signal || continue                 # skip annotation channels
        label = String(strip(sig.header.label))
        channels === nothing || label in channels || continue
        data[label] = Float64.(EDF.decode(sig))
        fs = sig.header.samples_per_record / edf.header.seconds_per_record
    end
    isempty(data) && error("No matching signals found in $path")
    return data, fs
end

const BV_TYPES = Dict("IEEE_FLOAT_32" => Float32, "INT_16" => Int16, "INT_32" => Int32)

"""
    load_brainvision(vhdr; channels=nothing) -> (data, fs)

Minimal BrainVision reader: parses the `.vhdr` header, then memory-maps the
binary `.eeg` file so only the requested channels are read into memory.
Supports IEEE_FLOAT_32, INT_16 and INT_32, multiplexed or vectorized.
"""
function load_brainvision(vhdr::AbstractString; channels=nothing)
    info = Dict{String,String}()
    chans = Dict{Int,Tuple{String,Float64}}()       # index → (label, µV per unit)
    section = ""
    for raw in eachline(vhdr)
        line = strip(raw)
        (isempty(line) || startswith(line, ";")) && continue
        if startswith(line, "[")
            section = line
            continue
        end
        occursin('=', line) || continue
        key, val = strip.(split(line, '='; limit=2))
        if section == "[Channel Infos]" && startswith(key, "Ch")
            parts = split(val, ',')
            label = replace(strip(parts[1]), "\\1" => ",")
            res = length(parts) >= 3 && !isempty(strip(parts[3])) ?
                  parse(Float64, strip(parts[3])) : 1.0
            chans[parse(Int, key[3:end])] = (label, res)
        else
            info[key] = val
        end
    end
    nch = parse(Int, info["NumberOfChannels"])
    fs = 1e6 / parse(Float64, info["SamplingInterval"])  # interval is in µs
    T = BV_TYPES[get(info, "BinaryFormat", "IEEE_FLOAT_32")]

    # BIDS renames files, so fall back to the .eeg next to the .vhdr.
    datafile = joinpath(dirname(vhdr), get(info, "DataFile", ""))
    isfile(datafile) || (datafile = replace(vhdr, r"\.vhdr$" => ".eeg"))
    nsamp = filesize(datafile) ÷ (sizeof(T) * nch)
    multiplexed = uppercase(get(info, "DataOrientation", "MULTIPLEXED")) == "MULTIPLEXED"

    data = Dict{String,Vector{Float64}}()
    open(datafile) do io
        M = multiplexed ? Mmap.mmap(io, Matrix{T}, (nch, nsamp)) :
                          Mmap.mmap(io, Matrix{T}, (nsamp, nch))
        for i in 1:nch
            label, res = chans[i]
            channels === nothing || label in channels || continue
            x = multiplexed ? @view(M[i, :]) : @view(M[:, i])
            data[label] = Float64.(x) .* res
        end
    end
    isempty(data) && error("No matching channels found in $vhdr")
    return data, fs
end

"""
    read_events(path) -> DataFrame

Read a BIDS `_events.tsv` file. In this dataset `sample` is 1-based
(onset 0.625 s at 512 Hz → sample 321), so it can index Julia arrays directly.
"""
read_events(path::AbstractString) =
    CSV.read(path, DataFrame; delim='\t', missingstring="n/a")

"""
    read_vmrk(path) -> DataFrame (sample, value)

Read the markers of a BrainVision `.vmrk` file (`Mk<n>=<type>,<description>,
<position>,...`; positions are 1-based samples). The numeric part of the
description becomes `value` ("S 51" → 51); markers without a number are skipped.
"""
function read_vmrk(path::AbstractString)
    samples, values = Int[], Int[]
    for line in eachline(path)
        startswith(line, "Mk") || continue
        f = split(split(line, '=', limit=2)[end], ',')
        length(f) >= 3 || continue
        code = tryparse(Int, filter(isdigit, f[2]))
        pos = tryparse(Int, strip(f[3]))
        (code === nothing || pos === nothing) && continue
        push!(samples, pos); push!(values, code)
    end
    return DataFrame(sample=samples, value=values)
end

"Event codes as integers; non-numeric codes (e.g. \"rt\") and missing become -1."
event_codes(events::DataFrame) =
    [v isa Integer ? Int(v) : something(tryparse(Int, string(coalesce(v, ""))), -1)
     for v in events.value]

# ---------------------------------------------------------------------------
# Continuous preprocessing (MATLAB Script #1 and #4)
# ---------------------------------------------------------------------------

"""
    rereference!(data; refs, skip)

Subtract the average of the two mastoids from every channel, including the
raw EOG channels (final reference in Script #4), so that the bipolar
VEOG (VEOG-lower − Fp2) is reference-free. Bipolar channels are skipped.
BioSemi records reference-free (CMS/DRL).
"""
function rereference!(data::Dict{String,Vector{Float64}};
                      refs=("Mastoid-left", "Mastoid-right"),
                      skip=("HEOG", "VEOG"))
    ref = (data[refs[1]] .+ data[refs[2]]) ./ 2
    for (ch, x) in data
        (ch in refs || ch in skip) && continue
        x .-= ref
    end
    return data
end

"""
    add_bipolars!(data)

Bipolar EOG channels as in Script #1: HEOG = HEOG-left − HEOG-right,
VEOG = VEOG-lower − Fp2. (Uncorrected, as used for artifact detection.)
"""
function add_bipolars!(data::Dict{String,Vector{Float64}})
    data["HEOG"] = data["HEOG-left"] .- data["HEOG-right"]
    data["VEOG"] = data["VEOG-lower"] .- data["Fp2"]
    return data
end

"Remove the DC offset (mean) of every channel."
function remove_dc!(data::Dict{String,Vector{Float64}})
    for (_, x) in data
        x .-= mean(x)
    end
    return data
end

"""
    butter_filter!(data, fs; kind=:highpass, cutoff, order)

Zero-phase Butterworth filter matching ERPLAB's `pop_basicfilter` /
`pop_filterp` with `'Design','butter'`: ERPLAB's `Order` is the combined
order of the forward-backward filter, so order 2 = 12 dB/oct and order 8 =
48 dB/oct, and the cutoff is the half-amplitude point.
`kind` is `:highpass`, `:lowpass`, or `:bandpass` (with `cutoff = (lo, hi)`).
"""
function butter_filter!(data::Dict{String,Vector{Float64}}, fs::Real;
                        kind::Symbol=:highpass, cutoff, order::Integer=2)
    nyq = fs / 2
    resp = kind === :highpass ? Highpass(cutoff / nyq) :
           kind === :lowpass  ? Lowpass(cutoff / nyq) :
           kind === :bandpass ? Bandpass(cutoff[1] / nyq, cutoff[2] / nyq) :
           error("kind must be :highpass, :lowpass or :bandpass")
    f = convert(SecondOrderSections, digitalfilter(resp, Butterworth(max(1, order ÷ 2))))
    for (_, x) in data
        x .= filtfilt(f, x)
    end
    return data
end

"""
    lowpass_epochs!(E, fs; cutoff, order) -> E

Zero-phase Butterworth low-pass applied to each epoch (row of `E`) separately.
filtfilt on segments of equal length is linear, so the average of the filtered
epochs equals the filtered average: this reproduces ERPLAB's `pop_filterp` on
the averaged ERPs (Script #7a), including its behavior at the epoch edges.
Rows with NaN (epochs outside the recording) are left as they are.
"""
function lowpass_epochs!(E::AbstractMatrix, fs::Real; cutoff, order::Integer=8)
    f = convert(SecondOrderSections,
                digitalfilter(Lowpass(cutoff / (fs / 2)), Butterworth(max(1, order ÷ 2))))
    for k in axes(E, 1)
        row = @view E[k, :]
        any(isnan, row) && continue
        row .= filtfilt(f, collect(row))
    end
    return E
end

# ---------------------------------------------------------------------------
# Trial labels
# ---------------------------------------------------------------------------
# Event order in every trial: memory array (51–58) → probe 1 s later (11–42)
# → response (120 correct, 121 incorrect). The memory-array code gives the
# condition and cued side; the probe code adds change and probe type.
# (BDF_stimulus.txt time-locks to 51–58 and requires 120, ignoring 11–42.)
# Note: the dataset README and events.json label these two code sets the
# other way round; the timing in the events files shows the order above.

"Memory-array code → (condition, cued side), as in BDF_stimulus.txt bins 1–8."
const ARRAY_CODES = Dict{Int,NTuple{2,String}}(
    51 => ("2-UR", "left"),  52 => ("3-UR", "left"),
    53 => ("3-RWP", "left"), 54 => ("3-RSP", "left"),
    55 => ("2-UR", "right"), 56 => ("3-UR", "right"),
    57 => ("3-RWP", "right"), 58 => ("3-RSP", "right"),
)

"Probe code → (condition, cued side, change, probe type)."
const PROBE_CODES = Dict{Int,NTuple{4,String}}(
    11 => ("2-UR",  "left",  "change",    "unrepeated"),
    12 => ("3-UR",  "left",  "change",    "unrepeated"),
    13 => ("3-RWP", "left",  "change",    "unrepeated"),
    14 => ("3-RSP", "left",  "change",    "unrepeated"),
    17 => ("3-RWP", "left",  "change",    "repeated"),
    18 => ("3-RSP", "left",  "change",    "repeated"),
    19 => ("2-UR",  "left",  "no-change", "unrepeated"),
    20 => ("3-UR",  "left",  "no-change", "unrepeated"),
    21 => ("3-RWP", "left",  "no-change", "unrepeated"),
    22 => ("3-RSP", "left",  "no-change", "unrepeated"),
    25 => ("3-RWP", "left",  "no-change", "repeated"),
    26 => ("3-RSP", "left",  "no-change", "repeated"),
    27 => ("2-UR",  "right", "change",    "unrepeated"),
    28 => ("3-UR",  "right", "change",    "unrepeated"),
    29 => ("3-RWP", "right", "change",    "unrepeated"),
    30 => ("3-RSP", "right", "change",    "unrepeated"),
    33 => ("3-RWP", "right", "change",    "repeated"),
    34 => ("3-RSP", "right", "change",    "repeated"),
    35 => ("2-UR",  "right", "no-change", "unrepeated"),
    36 => ("3-UR",  "right", "no-change", "unrepeated"),
    37 => ("3-RWP", "right", "no-change", "unrepeated"),
    38 => ("3-RSP", "right", "no-change", "unrepeated"),
    41 => ("3-RWP", "right", "no-change", "repeated"),
    42 => ("3-RSP", "right", "no-change", "repeated"),
)
const RESPONSE_CORRECT, RESPONSE_INCORRECT = 120, 121

"""
    swap_response_codes!(codes)

Swap 120 and 121 in place, as Script #0 did for participants 21 and 22
(reversed response keys). Only needed if the BIDS files were exported before
that recoding.
"""
function swap_response_codes!(codes::AbstractVector{Int})
    for i in eachindex(codes)
        codes[i] == RESPONSE_CORRECT   ? (codes[i] = RESPONSE_INCORRECT) :
        codes[i] == RESPONSE_INCORRECT ? (codes[i] = RESPONSE_CORRECT) : nothing
    end
    return codes
end

"""
    label_trials(codes, samples) -> DataFrame

One row per trial, starting at each memory-array code (51–58): array sample
(the time-locking event for N2pc and CDA), probe sample, condition, cued
side, change, probe type, accuracy, and whether the probe code agrees with
the array code (`consistent`).
"""
function label_trials(codes::AbstractVector{Int}, samples::AbstractVector{<:Integer})
    rows = NamedTuple[]
    for i in eachindex(codes)
        haskey(ARRAY_CODES, codes[i]) || continue
        cond, cue = ARRAY_CODES[codes[i]]
        probe_sample, change, probe, correct, consistent =
            missing, missing, missing, missing, missing
        j = i + 1
        while j <= length(codes) && !haskey(ARRAY_CODES, codes[j])
            c = codes[j]
            if haskey(PROBE_CODES, c) && probe_sample === missing
                pc, pcue, change, probe = PROBE_CODES[c]
                probe_sample = samples[j]
                consistent = (pc == cond) && (pcue == cue)
            elseif c == RESPONSE_CORRECT || c == RESPONSE_INCORRECT
                correct = c == RESPONSE_CORRECT
                break
            end
            j += 1
        end
        push!(rows, (stim_sample=Int(samples[i]), probe_sample=probe_sample,
                     condition=cond, cue=cue, change=change, probe=probe,
                     correct=correct, consistent=consistent))
    end
    return DataFrame(rows)
end

# ---------------------------------------------------------------------------
# Epochs (MATLAB Script #5)
# ---------------------------------------------------------------------------

"""
    epoch(x, onsets, fs; tmin, tmax, baseline) -> (E, times)

Cut one channel into epochs (trials × time), baseline-corrected to the mean
of `baseline` (seconds). Epochs that run past the recording are NaN.
"""
function epoch(x::AbstractVector{<:Real}, onsets::AbstractVector{<:Integer}, fs::Real;
               tmin=-0.2, tmax=1.0, baseline=(-0.2, 0.0))
    # As EEGLAB's pop_epoch: [-200 1000] ms at 256 Hz gives -199.2 to 996.1 ms
    # (307 samples; the end point is excluded).
    i0, i1 = ceil(Int, tmin * fs - 1e-9), ceil(Int, tmax * fs - 1e-9) - 1
    times = collect((i0:i1) ./ fs)
    bl = baseline === nothing ? Int[] :
         findall(t -> baseline[1] <= t <= baseline[2], times)
    E = fill(NaN, length(onsets), length(times))
    for (k, s) in enumerate(onsets)
        a, b = s + i0, s + i1
        (a < 1 || b > length(x)) && continue
        E[k, :] .= @view x[a:b]
        isempty(bl) || (E[k, :] .-= mean(@view E[k, bl]))
    end
    return E, times
end

# ---------------------------------------------------------------------------
# Artifact detection (MATLAB Script #6a, ERPLAB functions)
# ---------------------------------------------------------------------------

# Sample conventions of ERPLAB 10.02, from its source:
# - test period (window2sample): each edge goes to sample round(t * fs) relative
#   to time 0 (MATLAB round, ties away from zero), clipped to the epoch. At 256 Hz
#   the N2pc saccade period −50 to 300 ms becomes −50.8 to 300.8 ms.
# - window and step (pop_artstep, pop_artmwppth): floor(ms * fs / 1000) samples
#   (200 ms → 51, 10 ms → 2, 100 ms → 25).
# - a window is flagged when the value is strictly greater than the threshold.
_mround(x) = round(Int, x, RoundNearestTiesAway)
function _period(times, p, fs)
    off = 1 - _mround(times[1] * fs)               # index of the sample at time 0
    a = clamp(_mround(p[1] * fs) + off, 1, length(times))
    b = clamp(_mround(p[2] * fs) + off, 1, length(times))
    return a:b
end
_nsamp(ms, fs) = max(1, floor(Int, ms * fs / 1000 + 1e-9))

"""
    art_step(row, times, fs; period, threshold, window_ms, step_ms) -> Bool

ERPLAB `pop_artstep`: slide a window of `window_ms` in steps of `step_ms`
across `period`; flag the trial if the difference between the mean of the
second half and the mean of the first half of any window exceeds `threshold`.
Used for blinks (VEOG) and saccades (HEOG).
"""
function art_step(row::AbstractVector, times, fs; period, threshold, window_ms, step_ms)
    idx = _period(times, period, fs)
    w, s = _nsamp(window_ms, fs), _nsamp(step_ms, fs)
    h = _mround(w / 2)                             # first half: round(51/2) = 26 samples
    for start in first(idx):s:(last(idx) - w + 1)
        a = @view row[start:start+h-1]
        b = @view row[start+h:start+w-1]
        abs(mean(b) - mean(a)) > threshold && return true
    end
    return false
end

"ERPLAB `pop_artextval`: flag if any sample in `period` is beyond ±`threshold`."
art_extreme(row::AbstractVector, times, fs; period, threshold, kwargs...) =
    any(x -> abs(x) > threshold, @view row[_period(times, period, fs)])

"""
ERPLAB `pop_artmwppth`: flag if the peak-to-peak amplitude within any moving
window (`window_ms`, `step_ms`) inside `period` exceeds `threshold`.
"""
function art_mwpp(row::AbstractVector, times, fs; period, threshold, window_ms, step_ms)
    idx = _period(times, period, fs)
    w, s = _nsamp(window_ms, fs), _nsamp(step_ms, fs)
    for start in first(idx):s:(last(idx) - w + 1)
        seg = @view row[start:start+w-1]
        maximum(seg) - minimum(seg) > threshold && return true
    end
    return false
end

"""
    flag_trials(E, times, fs, detector; params...) -> BitVector

Apply an artifact `detector` to every trial of the epoch matrix `E`.
`true` = artifact. Trials containing NaN are always flagged.
"""
flag_trials(E::AbstractMatrix, times, fs, detector::Function; params...) =
    BitVector([any(isnan, @view E[k, :]) || detector(@view(E[k, :]), times, fs; params...)
               for k in axes(E, 1)])

# ---------------------------------------------------------------------------
# Lateralized ERPs (MATLAB Script #7a)
# ---------------------------------------------------------------------------

"""
    contra_minus_ipsi(E_left, E_right, cues) -> D

Contralateral minus ipsilateral waveform per trial for one electrode pair.
"""
function contra_minus_ipsi(E_left::AbstractMatrix, E_right::AbstractMatrix,
                           cues::AbstractVector{<:AbstractString})
    D = similar(E_left)
    for k in axes(D, 1)
        D[k, :] .= cues[k] == "left" ? E_right[k, :] .- E_left[k, :] :
                                        E_left[k, :] .- E_right[k, :]
    end
    return D
end

"""
    cluster_contra_minus_ipsi(epochs, pairs, cues) -> D

Average the contra-minus-ipsi waveform over a cluster of electrode pairs
(e.g. P3/4, P5/6, P7/8, P9/10, PO7/8, PO3/4, O1/2), like the cluster of
interest in Script #7a. `epochs` maps channel label → epoch matrix.
"""
function cluster_contra_minus_ipsi(epochs::AbstractDict, pairs, cues)
    Ds = [contra_minus_ipsi(epochs[l], epochs[r], cues) for (l, r) in pairs]
    return sum(Ds) ./ length(Ds)
end

"""
    erplab_window(times, window) -> UnitRange

Sample indices of a measurement window as in ERPLAB's `pop_geterpvalues`: from
the sample closest to the window start to the sample closest to the window end
(e.g. 238.3–300 ms at 256 Hz → 238.28–300.78 ms).
"""
erplab_window(times::AbstractVector, window) =
    argmin(abs.(times .- window[1])):argmin(abs.(times .- window[2]))

"Mean amplitude of each trial within a time window (seconds), ERPLAB sample selection."
function mean_amplitude(E::AbstractMatrix, times::AbstractVector, window)
    return vec(mean(E[:, erplab_window(times, window)]; dims=2))
end

# ---------------------------------------------------------------------------
# Alpha power (exploratory; MATLAB Step1_Compute_Alpha_Power)
# ---------------------------------------------------------------------------

"""
    matlab_filter(path, fs; fpass=10.0) -> filter, kind

Read the coefficients of MATLAB's `bandpass()` filter exported to CSV
(see README): 6 columns = second-order sections [b0 b1 b2 a0 a1 a2] (IIR),
1 row or 1 column = FIR coefficients. The gain is normalized to 1 at `fpass`.
"""
function matlab_filter(path::AbstractString, fs::Real; fpass=10.0)
    C = Matrix{Float64}(CSV.read(path, DataFrame; header=false))
    if size(C, 2) == 6
        bq = [Biquad(r[1] / r[4], r[2] / r[4], r[3] / r[4], r[5] / r[4], r[6] / r[4])
              for r in eachrow(C)]
        f = SecondOrderSections(bq, 1.0)
        g = abs(freqresp(f, 2π * fpass / fs))
        return SecondOrderSections(bq, 1 / g), :iir
    else
        b = vec(C)
        return b, :fir
    end
end

"""
    alpha_power(E, fs; band=(8, 12), pad_s=0.5, design=:elliptic, ...) -> P

Time-resolved band power for every trial (trials × time), as in MATLAB's
Step1_Compute_Alpha_Power: mirror-pad each trial by `pad_s` seconds, band-pass
filter, remove the padding, then take |Hilbert|². Trials containing NaN stay NaN.
"""
function alpha_power(E::AbstractMatrix, fs::Real; band=(8.0, 12.0), pad_s=0.5, order=4,
                     hilbert_padded::Bool=false, design::Symbol=:elliptic,
                     matlab_coeffs::AbstractString="")
    nyq = fs / 2
    # :matlab = the exact filter MATLAB's bandpass() built, read from CSV; IIR is
    # applied zero-phase (filtfilt), FIR once with its delay removed, as bandpass() does.
    # :elliptic = an approximation of bandpass() defaults (prototype order 6, 0.1 dB, 60 dB).
    # :butterworth = a gentle Butterworth (prototype order÷2).
    kind = :iir
    if design === :matlab
        # MATLAB's exact bandpass() filter (sos exported with zp2sos), applied
        # with MATLAB's filtfilt: identical to bandpass(x, [8 12], fs).
        S = Matrix{Float64}(CSV.read(matlab_coeffs, DataFrame; header=false))
        size(S, 2) == 6 || error("$matlab_coeffs: expected 6 columns (sos matrix from zp2sos)")
        return band_power(E, fs, S, :matlab; pad_s=pad_s, hilbert_padded=hilbert_padded)
    else
        # Minimum-order designs for MATLAB bandpass() defaults at 256 Hz (passband
        # 8–12 Hz, stopbands below 6.8 and above 29.4 Hz, 0.1 dB ripple, 60 dB):
        # elliptic 6, Chebyshev I 8, Chebyshev II 8, Butterworth 15 (scipy *ord).
        edges, proto = design === :elliptic ? ((8.0, 12.0), Elliptic(6, 0.1, 60)) :
                       design === :cheby1   ? ((8.0, 12.0), Chebyshev1(8, 0.1)) :
                       design === :cheby2   ? ((6.8497, 13.9905), Chebyshev2(8, 60)) :
                       design === :butter15 ? ((7.7891, 12.3221), Butterworth(15)) :
                       ((band[1], band[2]), Butterworth(order ÷ 2))   # :butterworth
        f = convert(SecondOrderSections,
                    digitalfilter(Bandpass(edges[1] / nyq, edges[2] / nyq), proto))
    end
    bandfilt(xp) = kind === :iir ? filtfilt(f, xp) :
                   begin
                       d = (length(f) - 1) ÷ 2
                       y = filt(f, vcat(xp, zeros(d)))
                       y[d+1:end]
                   end
    pad = round(Int, pad_s * fs)
    P = fill(NaN, size(E))
    for k in axes(E, 1)
        x = E[k, :]
        any(isnan, x) && continue
        xp = vcat(reverse(x[1:pad]), x, reverse(x[end-pad+1:end]))
        yp = bandfilt(xp)
        # As in the MATLAB script, the Hilbert transform is applied after the
        # padding is removed, which distorts power near the epoch edges.
        # hilbert_padded=true applies it to the padded signal instead.
        P[k, :] .= hilbert_padded ? abs2.(hilbert(yp))[pad+1:end-pad] :
                                    abs2.(hilbert(yp[pad+1:end-pad]))
    end
    return P
end

"Percent change from the mean of the baseline period."
function percent_change(p::AbstractVector, times, baseline)
    b = mean(p[(times .>= baseline[1]) .& (times .<= baseline[2])])
    return 100 .* (p .- b) ./ b
end

# ---------------------------------------------------------------------------
# ICA blink correction (MATLAB Scripts #2–#4, automated)
# ---------------------------------------------------------------------------

"""
    break_mask(n, event_samples, fs; gap_s=5, pre_s=0.5, post_s=1.5) -> BitVector

Samples to keep for ICA training, as `pop_erplabDeleteTimeSegments` in Script #2:
periods of more than `gap_s` without an event code (breaks) are removed, keeping
`pre_s` before the next event and `post_s` after the previous one.
"""
function break_mask(n::Integer, ev::AbstractVector{<:Integer}, fs::Real;
                    gap_s=5.0, pre_s=0.5, post_s=1.5)
    keep = falses(n)
    ev = sort(ev)
    for (k, s) in enumerate(ev)
        a = max(1, s - round(Int, pre_s * fs))
        b = min(n, s + round(Int, post_s * fs))
        a <= b && (keep[a:b] .= true)
        if k < length(ev) && ev[k+1] - s <= gap_s * fs
            keep[max(1, s):min(n, ev[k+1])] .= true
        end
    end
    return keep
end

"""
    ica_blink_correct!(data, fs, event_samples; scalp, veog="VEOG", n_remove=1,
                       n_components=30, train_fs=100, seed=1) -> NamedTuple

Automated version of the ICA correction in Scripts #2–#4:

1. Training copy (Script #2): scalp channels band-passed 1–30 Hz (Butterworth,
   48 dB/oct, zero-phase), resampled to 100 Hz, breaks removed.
2. ICA (Script #3): FastICA with `n_components` components (the original used
   extended infomax, `runica`, on all channels).
3. Correction (Script #4): the weights are applied to the analysis data, and the
   `n_remove` components most correlated with the uncorrected bipolar VEOG are
   removed. (The original removed hand-picked blink components; `n_remove`
   is the number removed for that participant.)

Only `scalp` channels are corrected; EOG channels stay uncorrected for artifact
detection, as in the original. Returns the removed components and their
correlation with VEOG.
"""
function ica_blink_correct!(data::Dict{String,Vector{Float64}}, fs::Real,
                            event_samples::AbstractVector{<:Integer};
                            scalp::AbstractVector{<:AbstractString}, veog::AbstractString="VEOG",
                            n_remove::Integer=1, n_components::Integer=30,
                            train_fs::Real=100.0, seed::Integer=1)
    n_remove == 0 && return (removed=Int[], r=Float64[])
    nyq = fs / 2
    bp = convert(SecondOrderSections,
                 digitalfilter(Bandpass(1 / nyq, 30 / nyq), Butterworth(4)))
    n = length(data[scalp[1]])
    keep = break_mask(n, event_samples, fs)
    ratio = rationalize(train_fs / fs)
    rows = map(scalp) do ch
        y = resample(filtfilt(bp, data[ch]), ratio)
        idx = clamp.(round.(Int, (0:length(y)-1) .* (fs / train_fs)) .+ 1, 1, n)
        y[keep[idx]]
    end
    m = minimum(length, rows)
    T = reduce(vcat, (permutedims(r[1:m]) for r in rows))       # channels × samples

    M = nothing
    for k in (min(n_components, length(scalp)), 20, 15)
        try
            Random.seed!(seed)
            M = fit(ICA, T, k; maxiter=1000, tol=1e-4)
            break
        catch err
            @warn "ICA with $k components did not converge; trying fewer" exception=err
        end
    end
    M === nothing && (@warn "ICA failed: no correction applied"; return (removed=Int[], r=Float64[]))

    W, mu = M.W, M.mean                                         # sources = W' (x - mu)
    A = pinv(Matrix(W'))                                        # mixing (channels × k)
    X = reduce(vcat, (permutedims(data[ch]) for ch in scalp))
    S = W' * (X .- mu)
    v = data[veog]
    r = [cor(view(S, j, :), v) for j in axes(S, 1)]
    rm = sortperm(abs.(r); rev=true)[1:min(n_remove, length(r))]
    X .-= A[:, rm] * S[rm, :]
    for (i, ch) in enumerate(scalp)
        data[ch] .= view(X, i, :)
    end
    return (removed=rm, r=r[rm])
end

"""
    find_ica_set(root, id) -> path

Find `<id>_allERPs_optimized_weights.set` (Script #3 output): first in the
participant folder `<root>/<id>/` (the derivatives layout), otherwise anywhere under `root`.
"""
function find_ica_set(root::AbstractString, id::Integer)
    name = "$(id)_allERPs_optimized_weights.set"
    direct = joinpath(root, string(id), name)        # <root>/<id>/<id>_allERPs_...set
    isfile(direct) && return direct
    for (dir, _, files) in walkdir(root)
        name in files && return joinpath(dir, name)
    end
    error("ICA weights file $name not found under $root. See README (ICA weights).")
end

"""
    read_ica_weights(set_path, remove) -> (unmixing, winv, channels, remove)

Read the ICA solution from an EEGLAB `.set` file (MATLAB format, read with
MAT.jl): unmixing = icaweights × icasphere, its pseudo-inverse (as EEGLAB recomputes
icawinv when the weights are transferred in Script #4), and the labels of the
channels in icachansind. `remove` = the components to remove (ICs_to_Remove).
Works for both EEGLAB layouts (fields saved at the top level or in `EEG`).
"""
function read_ica_weights(set_path::AbstractString, remove::AbstractVector{<:Integer})
    get_fields(f) = begin
        ks = keys(f)
        if "icaweights" in ks
            Dict(k => read(f, k) for k in ("icaweights", "icasphere", "icawinv",
                                             "icachansind", "chanlocs"))
        else
            read(f, "EEG")
        end
    end
    E = matopen(get_fields, set_path)
    W = Matrix{Float64}(E["icaweights"]) * Matrix{Float64}(E["icasphere"])
    # Script #4 transfers only icaweights, icasphere and icachansind (pop_editset), so
    # EEGLAB recomputes icawinv = pinv(icaweights * icasphere) (eeg_checkset). Do the
    # same instead of reading the stored icawinv, which can be out of date.
    winv = pinv(W)
    winv_file = haskey(E, "icawinv") ? Matrix{Float64}(E["icawinv"]) : winv
    size(winv_file) == size(winv) && maximum(abs.(winv_file .- winv)) > 1e-6 * maximum(abs.(winv)) &&
        @warn "icawinv stored in $(basename(set_path)) differs from pinv(icaweights*icasphere); using pinv, as EEGLAB does after pop_editset"
    idx = Int.(vec(E["icachansind"]))
    isempty(idx) && (idx = collect(1:size(W, 2)))         # empty = all channels
    labels = vec(string.(E["chanlocs"]["labels"]))
    return (unmixing=W, winv=winv, channels=String.(strip.(labels[idx])),
            remove=collect(Int, remove))
end

"""
    ics_to_remove(csv, id) -> Vector{Int}

Components removed for participant `id` in the original analysis
(`scripts/ICs_to_Remove.csv`, converted from ICs_to_Remove.xlsx).
"""
function ics_to_remove(csv::AbstractString, id::Integer)
    T = CSV.read(csv, DataFrame; types=Dict(:ICs_to_remove => String))
    r = T[T.id .== id, :ICs_to_remove]
    (isempty(r) || ismissing(r[1])) && return Int[]
    return parse.(Int, split(r[1]))
end

"""
    apply_ica_weights!(data, ica) -> data

Remove ICA components as EEGLAB's `pop_subcomp` does (Script #4): activations
S = W · X on the ICA channels, then X = icawinv(:, kept) · S(kept, :).
`data` must be in the reference the ICA was computed in (left mastoid).
"""
function apply_ica_weights!(data::Dict{String,Vector{Float64}}, ica)
    isempty(ica.remove) && return data
    missing_ch = setdiff(ica.channels, collect(keys(data)))
    isempty(missing_ch) || error("ICA channels not in the data: $(join(missing_ch, ", "))")
    X = reduce(vcat, (permutedims(data[ch]) for ch in ica.channels))
    S = ica.unmixing * X
    keep = setdiff(1:size(S, 1), ica.remove)
    Xc = ica.winv[:, keep] * S[keep, :]
    for (i, ch) in enumerate(ica.channels)
        data[ch] .= view(Xc, i, :)
    end
    return data
end

"""
    eegfiltnew_lowpass(x, fs, hicutoff) -> y

Low-pass filter as EEGLAB's `pop_eegfiltnew(EEG, 'hicutoff', hicutoff)` on one
epoch: windowed-sinc FIR (Hamming), default transition band
max(min(0.25·hicutoff, fs/2 − hicutoff), 2) Hz, order 3.3/(Δf/fs) rounded up to
an even number, −6 dB cutoff at hicutoff + Δf/2, zero-phase by removing the
group delay, and the epoch padded with its first and last values (as `firfilt`).
"""
function eegfiltnew_lowpass(x::AbstractVector{<:Real}, fs::Real, hicutoff::Real)
    df = max(min(0.25 * hicutoff, fs / 2 - hicutoff), 2.0)
    m = 2 * ceil(Int, 3.3 / (df / fs) / 2)                     # filter order (even)
    fc = (hicutoff + df / 2) / fs                              # cycles per sample
    n = 0:m
    h = [k == m ÷ 2 ? 2fc : sin(2π * fc * (k - m ÷ 2)) / (π * (k - m ÷ 2)) for k in n]
    w = 0.54 .- 0.46 .* cos.(2π .* n ./ m)                     # Hamming
    b = h .* w
    b ./= sum(b)
    d = m ÷ 2
    xp = vcat(fill(x[1], d), x, fill(x[end], d))
    y = filt(b, xp)                                            # causal, delay d
    return y[2d+1:end]
end

"""
    kaiser_bandpass(fs, f1, f2; transition, attenuation=60) -> taps

Linear-phase FIR band-pass (Kaiser window, as MATLAB's minimum-order FIR
designs): cutoffs `f1`, `f2` (Hz, middle of the transition bands),
`transition` = transition width (Hz), `attenuation` in dB. Odd length.
"""
function kaiser_bandpass(fs::Real, f1::Real, f2::Real; transition::Real, attenuation::Real=60)
    dw = 2π * transition / fs
    n = ceil(Int, (attenuation - 7.95) / (2.285 * dw)) + 1
    isodd(n) || (n += 1)
    β = attenuation > 50 ? 0.1102 * (attenuation - 8.7) :
        attenuation >= 21 ? 0.5842 * (attenuation - 21)^0.4 + 0.07886 * (attenuation - 21) : 0.0
    return digitalfilter(Bandpass(f1 / (fs / 2), f2 / (fs / 2)), FIRWindow(kaiser(n, β / π)))
end

"""
    band_power(E, fs, f, mode; pad_s=0.5, hilbert_padded=false) -> P

Alpha power of every trial (rows of `E`) as in Step1_Compute_Alpha_Power.m,
for any filter `f`: mirror padding, band-pass, padding removed, |Hilbert|².
`mode`: `:matlab` (IIR sos matrix, MATLAB's filtfilt, as bandpass()),
`:zerophase` (IIR, DSP.jl filtfilt), `:causal` (IIR, one pass), or `:fir`
(FIR taps, one pass with the group delay removed).
"""
function band_power(E::AbstractMatrix, fs::Real, f, mode::Symbol; pad_s=0.5, hilbert_padded=false)
    pad = round(Int, pad_s * fs)
    P = fill(NaN, size(E))
    for k in axes(E, 1)
        x = E[k, :]
        any(isnan, x) && continue
        xp = vcat(reverse(x[1:pad]), x, reverse(x[end-pad+1:end]))
        yp = mode === :matlab ? matlab_filtfilt(f, xp) :
             mode === :zerophase ? filtfilt(f, xp) :
             mode === :causal ? filt(f, xp) :
             begin
                 d = (length(f) - 1) ÷ 2
                 filt(f, vcat(xp, zeros(d)))[d+1:end]
             end
        P[k, :] .= hilbert_padded ? abs2.(hilbert(yp))[pad+1:end-pad] :
                                    abs2.(hilbert(yp[pad+1:end-pad]))
    end
    return P
end

"""
    sos_matrix(f::SecondOrderSections) -> Matrix (sections × 6)

Second-order sections as MATLAB's sos matrix [b0 b1 b2 a0 a1 a2], with the
overall gain in the first section (as `zp2sos`).
"""
function sos_matrix(f::SecondOrderSections)
    S = reduce(vcat, ([q.b0 q.b1 q.b2 1.0 q.a1 q.a2] for q in f.biquads))
    S[1, 1:3] .*= f.g
    return S
end

# Direct form II transposed, as MATLAB's filter(b, a, x, zi), for one section.
function _df2t(b, a, x::AbstractVector, z1::Real, z2::Real)
    y = similar(x, Float64)
    @inbounds for n in eachindex(x)
        yn = b[1] * x[n] + z1
        z1 = b[2] * x[n] - a[2] * yn + z2
        z2 = b[3] * x[n] - a[3] * yn
        y[n] = yn
    end
    return y
end

"""
    eeglab_resample(x, q; fc=0.9, df=0.2) -> y

Downsampling by an integer factor `q` exactly as EEGLAB's `pop_resample` with the
Signal Processing Toolbox (Script #1: 512 → 256 Hz, q = 2; 2048 → 256 Hz, q = 8).
EEGLAB's own anti-aliasing filter: a windowed-sinc FIR (`firws`) with cutoff
`fc`·π/q, Kaiser window β = 5, order from `firwsord('kaiser', 2, df/q, 0.002)`
(66 for q = 2), normalised to unit gain at DC. The data are padded with their
first and last values (`nPad` samples, a multiple of q), filtered and decimated
as MATLAB's `resample(x, 1, q, b)` (centred kernel, so output sample k is at
input sample q·k), and the padding is removed. Output length ceil(length(x)/q).
"""
function eeglab_resample(x::AbstractVector, q::Integer; fc=0.9, df=0.2)
    q == 1 && return Float64.(x)
    fcq, dfq = fc / q, df / q
    devdb = -20 * log10(0.002)
    m = Int(ceil((1 + (devdb - 8) / (2.285 * 2π * dfq / 2)) / 2) * 2)   # firwsord, even
    n = -m ÷ 2:m ÷ 2
    b = [k == 0 ? π * fcq : sin(π * fcq * k) / k for k in n]
    β = 5.0
    w = [besseli0(β * sqrt(max(0.0, 1 - (2i / m - 1)^2))) / besseli0(β) for i in 0:m]
    b .*= w
    b ./= sum(b)
    npad = ceil(Int, (m / 2) / q) * q
    xx = vcat(fill(Float64(x[1]), npad), Float64.(x), fill(Float64(x[end]), npad))
    N = cld(length(x), q)
    y = Vector{Float64}(undef, N)
    h = m ÷ 2
    @inbounds for k in 0:N-1
        c = npad + q * k + 1                       # 1-based centre in xx
        acc = 0.0
        for j in 0:m
            acc += b[j+1] * xx[c+h-j]
        end
        y[k+1] = acc
    end
    return y
end

"Modified Bessel function of the first kind, order 0 (series; for the Kaiser window)."
function besseli0(x::Real)
    s, t, k = 1.0, 1.0, 0
    while true
        k += 1
        t *= (x / (2k))^2
        s += t
        t < 1e-17 * s && return s
    end
end

"""
    matlab_filtfilt(S, x) -> y

Zero-phase filtering as MATLAB's `filtfilt` (R2025 source, default options
'const' initial states, 'filtord' transient length, 'linear' padding) for a
second-order-section matrix `S` ([b0 b1 b2 a0 a1 a2] per row; what `bandpass()`
uses): the sections are applied one after another, each forward and backward,
each with its own odd-symmetric extension at both ends of 3 × the order of the
WHOLE filter (36 samples for 6 sections) and initial conditions from its step
response. (DSP.jl's `filtfilt` filters the whole cascade at once, which handles
the edges differently.)
"""
function matlab_filtfilt(S::AbstractMatrix, x::AbstractVector)
    efflen(c) = (m = maximum(abs, c); m == 0 ? 0 : something(findlast(!iszero, c ./ m), 0))
    order = sum(max(efflen(r[1:3]), efflen(r[4:6])) - 1 for r in eachrow(S))
    nf = max(1, 3 * order)
    nf < length(x) - 1 || error("signal too short for filtfilt (transient length $nf)")
    y = Float64.(x)
    for r in eachrow(S)
        b = r[1:3] ./ r[4]
        a = r[4:6] ./ r[4]
        zi = [1+a[2] -1.0; a[3] 1.0] \ (b[2:3] .- b[1] .* a[2:3])
        n = length(y)
        xt = vcat(2y[1] .- y[nf+1:-1:2], y, 2y[n] .- y[n-1:-1:n-nf])
        xt = _df2t(b, a, xt, zi[1] * xt[1], zi[2] * xt[1])
        xt = _df2t(b, a, reverse(xt), zi[1] * xt[end], zi[2] * xt[end])
        y = reverse(xt)[nf+1:end-nf]
    end
    return y
end

end # module
