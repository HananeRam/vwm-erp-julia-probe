# Measurements (Script #11c) and statistics (R Steps 3–7, Models 4–9). Run after 01:
#   julia --project=. scripts/03_measures_and_statistics.jl
#
# 1. Mean amplitude, peak latency and 50% onset latency of the change − no-change
#    difference waves, per participant, as pop_geterpvalues in Script #11c.
# 2. Comparison with the MATLAB measurement files (reference/measurements), per
#    participant and bin.
# 3. Repeated-measures ANOVAs of Models 4–9, run on the Julia measures and, as a
#    check of the statistics code alone, on the MATLAB measures; table with the
#    published values.
# 4. Grand averages: comparison with the MATLAB grand averages, and Figures 5 and 6
#    of the paper (same design as Figures_N2_FN400_probe.m).

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "..", "src", "VWMPipeline.jl"))
include(joinpath(@__DIR__, "..", "src", "ProbeERP.jl"))
using .VWMPipeline: lowpass_epochs!
using .ProbeERP
using Printf, Plots
include(joinpath(@__DIR__, "bin_operations.jl"))

erps, times_ms, subjects = load_probe_erps()
ids = [parse(Int, s[5:end]) for s in subjects]
@info "$(length(subjects)) participants"

# --- 1. Measurements (Script #11c) ------------------------------------------------
const WINDOWS = Dict("N2" => N2_WINDOW, "FN400" => FN400_WINDOW)
const ONSET_WINDOWS = Dict("N2" => N2_ONSET_WINDOW, "FN400" => FN400_ONSET_WINDOW)
measure_bins = vcat(first.(PROBE_BINS), first.(STIMULUS_BINS))
long = DataFrame()
ga = DataFrame()
for (s, id) in zip(subjects, ids), comp in ("N2", "FN400")
    raw = n2fn400_bins(erps, s, comp)
    lp  = n2fn400_bins(erps, s, comp; lowpass=true)
    for b in measure_bins
        push!(long, (id=id, component=comp, bin=b, measure="mean_amplitude",
                     value=meanbl(raw[b].wave, times_ms, WINDOWS[comp])))
        push!(long, (id=id, component=comp, bin=b, measure="peak_latency",
                     value=peaklatbl(lp[b].wave, times_ms, WINDOWS[comp]; npoints=PEAK_NEIGHBORHOOD)))
        push!(long, (id=id, component=comp, bin=b, measure="onset_latency",
                     value=fpeaklat(lp[b].wave, times_ms, ONSET_WINDOWS[comp];
                                    frac=ONSET_FRACTION, npoints=PEAK_NEIGHBORHOOD)))
    end
    for b in 46:56, (filt, B) in (("raw", raw), ("lpfilt", lp))
        append!(ga, DataFrame(id=id, component=comp, bin=b, filter=filt, n_trials=B[b].n,
                              time_ms=times_ms, amplitude=B[b].wave))
    end
end
mkpath(RESULTS)
CSV.write(joinpath(RESULTS, "measures.csv"), long)

# --- 2. Comparison with the MATLAB measurement files -------------------------------
const REF_FILES = Dict(
    ("N2", "mean_amplitude") => ["Mean_Amplitude_Diff_Waves_N2_Probetype.txt", "Mean_Amplitude_Diff_Waves_N2_Stimulustype.txt"],
    ("FN400", "mean_amplitude") => ["Mean_Amplitude_Diff_Waves_FN400_Probetype.txt", "Mean_Amplitude_Diff_Waves_FN400_Stimulustype.txt"],
    ("N2", "peak_latency") => ["Peak_Latency_Diff_Waves_N2_Probetype.txt", "Peak_Latency_Diff_Waves_N2_Stimulustype.txt"],
    ("FN400", "peak_latency") => ["Peak_Latency_Diff_Waves_FN400_Probetype.txt", "Peak_Latency_Diff_Waves_FN400_Stimulustype.txt"],
    ("N2", "onset_latency") => ["Onset_Latency_Diff_Waves_N2_Probetype.txt", "Onset_Latency_Diff_Waves_N2_Stimulustype.txt"],
    ("FN400", "onset_latency") => ["Onset_Latency_Diff_Waves_FN400_Probetype.txt", "Onset_Latency_Diff_Waves_FN400_Stimulustype.txt"])

"Read a pop_geterpvalues 'wide' file → long DataFrame (id, bin, value)."
function read_erplab_values(path)
    L = [split(l, '\t') for l in readlines(path) if !isempty(strip(l))]
    header = strip.(L[1])
    bins = [parse(Int, match(r"^bin(\d+)", h).captures[1]) for h in header[1:end-1]]
    out = DataFrame(id=Int[], bin=Int[], value=Float64[])
    for row in L[2:end]
        f = strip.(row)
        id = parse(Int, split(f[end], '_')[1])
        for (k, b) in enumerate(bins)
            push!(out, (id, b, f[k] == "NaN" ? NaN : parse(Float64, f[k])))
        end
    end
    return out
end

matlab = DataFrame()
for ((comp, meas), files) in REF_FILES, f in files
    p = joinpath(REFERENCE, "measurements", f)
    isfile(p) || continue
    m = read_erplab_values(p)
    m.component = fill(comp, nrow(m)); m.measure = fill(meas, nrow(m))
    append!(matlab, m[in.(m.id, Ref(Set(SUBJECT_IDS))), :])   # same participants as Julia
end
cmp_rows = DataFrame()
if nrow(matlab) > 0
    both = innerjoin(rename(long, :value => :julia), rename(matlab, :value => :matlab),
                     on=[:id, :component, :bin, :measure])
    CSV.write(joinpath(RESULTS, "measures_julia_vs_matlab.csv"), both)
    for g in groupby(both, [:component, :measure, :bin])
        ok = .!isnan.(g.julia) .& .!isnan.(g.matlab)
        d = abs.(g.julia[ok] .- g.matlab[ok])
        push!(cmp_rows, (component=g.component[1], measure=g.measure[1], bin=g.bin[1],
                         n=nrow(g), n_both=count(ok),
                         nan_julia=count(isnan, g.julia), nan_matlab=count(isnan, g.matlab),
                         nan_mismatch=count(isnan.(g.julia) .!= isnan.(g.matlab)),
                         n_identical=count(d .< 1e-6), max_abs_diff=isempty(d) ? NaN : maximum(d),
                         mean_abs_diff=isempty(d) ? NaN : mean(d),
                         r=count(ok) > 2 ? cor(g.julia[ok], g.matlab[ok]) : NaN))
    end
    sort!(cmp_rows, [:component, :measure, :bin])
    CSV.write(joinpath(RESULTS, "measures_comparison.csv"), cmp_rows)
end

# --- 3. Statistics (R Steps 3–7) -------------------------------------------------------
wide_of(df, comp, meas, bins) = begin
    d = df[(df.component .== comp) .& (df.measure .== meas) .& in.(df.bin, Ref(bins)), :]
    w = unstack(d, :id, :bin, :value)
    rename!(w, [Symbol(b) => Symbol("b$b") for b in bins if Symbol(b) in propertynames(w)])
    w
end
cols(bins) = [Symbol("b$b") for b in bins]

function models(df)
    out = DataFrame()
    add!(label, a) = push!(out, (model=label, n=a.n, F=a.F, df1=a.df1, df2=a.df2, p=a.p,
                                 eta_p2=a.eta_p2, mauchly_p=a.mauchly_p, gg_eps=a.gg_eps,
                                 gg_df1=a.gg_df1, gg_df2=a.gg_df2, gg_p=a.gg_p))
    sb = first.(STIMULUS_BINS)
    for (comp, num) in (("N2", "5"), ("FN400", "6"))
        add!("$comp amplitude, stimulus type ($(num).a)", rm_anova_gg(wide_of(df, comp, "mean_amplitude", sb), cols(sb)))
        add!("$comp peak latency, stimulus type ($(num).b)", rm_anova_gg(wide_of(df, comp, "peak_latency", sb), cols(sb)))
    end
    for (comp, num) in (("N2", "8"), ("FN400", "9"))
        add!("$comp amplitude, probe type ($(num).a)", rm_anova_gg(wide_of(df, comp, "mean_amplitude", [50, 51]), cols([50, 51])))
        add!("$comp peak latency, probe type ($(num).b)", rm_anova_gg(wide_of(df, comp, "peak_latency", [50, 51]), cols([50, 51])))
    end
    # Model 4: onset latency; participants with a missing value in any of bins 46–51
    # are excluded from both analyses (as in R Step 3).
    twobytwo = DataFrame()
    for (comp, lab) in (("N2", "4.b"), ("FN400", "4.c"))
        w = wide_of(df, comp, "onset_latency", 46:51)
        w = w[[all(!isnan, coalesce.(collect(r[cols(46:51)]), NaN)) for r in eachrow(w)], :]
        add!("$comp onset latency, probe type ($lab)", rm_anova_gg(w, cols([50, 51])))
        t = rm_anova_2x2(w, cols([46, 47, 48, 49]); names=("probe", "stimulus"))
        t.component = fill(comp, nrow(t))
        append!(twobytwo, t)
    end
    return out, twobytwo
end

anova_julia, twobytwo_julia = models(long)
CSV.write(joinpath(RESULTS, "anova.csv"), anova_julia)
CSV.write(joinpath(RESULTS, "anova_2x2_onset.csv"), twobytwo_julia)
anova_matlab, twobytwo_matlab = nrow(matlab) > 0 ? models(matlab) : (nothing, nothing)

# Descriptives (mean, SD per condition), as summarySE in R.
desc = combine(groupby(long, [:component, :measure, :bin]),
               :value => (x -> count(!isnan, x)) => :n,
               :value => (x -> mean(filter(!isnan, x))) => :mean,
               :value => (x -> std(filter(!isnan, x))) => :sd)
CSV.write(joinpath(RESULTS, "descriptives.csv"), desc)

# Comparison table: Julia, statistics on the MATLAB measures, R output.
fmtF(a) = a.mauchly_p < 0.05 ? @sprintf("F(%.2f, %.2f) = %.2f, p = %.3f, ηp² = %.3f (GG)",
                                         a.gg_df1, a.gg_df2, a.F, a.gg_p, a.eta_p2) :
                               @sprintf("F(%g, %g) = %.2f, p = %.3f, ηp² = %.3f",
                                        a.df1, a.df2, a.F, a.p, a.eta_p2)
io = IOBuffer()
println(io, "Probe-locked N2/FN400: Julia vs. MATLAB/R (revised analysis)\n")
println(io, "Julia and MATLAB lines: N = $(length(subjects))",
        isempty(EXCLUDED_IDS) ? "" : " (participant(s) $(join(EXCLUDED_IDS, ", ")) excluded, see config.jl)",
        "; MATLAB/R: N = $(length(ALL_IDS)).")
# Optional sensitivity column: ALT_RESULTS=<folder of a run with ICS=listed>.
alt_dir = get(ENV, "ALT_RESULTS", "")
alt = isempty(alt_dir) || !isfile(joinpath(alt_dir, "anova.csv")) ? nothing :
      CSV.read(joinpath(alt_dir, "anova.csv"), DataFrame)
alt2 = alt === nothing || !isfile(joinpath(alt_dir, "anova_2x2_onset.csv")) ? nothing :
       CSV.read(joinpath(alt_dir, "anova_2x2_onset.csv"), DataFrame)
alt === nothing || (alt.model = String.(alt.model))
alt2 === nothing || (alt2.component = String.(alt2.component); alt2.effect = String.(alt2.effect))
println(io, "Lines: Julia pipeline | the same statistics on the MATLAB measurement files | R output",
        alt === nothing ? "" : "\n       | Julia with the ICA components listed in ICs_to_Remove.xlsx (46, 49, 69: [1])", "\n")
for r in eachrow(anova_julia)
    pp = PAPER[r.model]
    m = anova_matlab === nothing ? nothing : anova_matlab[anova_matlab.model .== r.model, :][1, :]
    println(io, r.model)
    println(io, "  Julia : ", fmtF(r))
    m === nothing || println(io, "  MATLAB: ", fmtF(m))
    if alt !== nothing
        a = alt[alt.model .== r.model, :]
        nrow(a) > 0 && println(io, "  Julia, ICs as listed: ", fmtF(a[1, :]))
    end
    @printf(io, "  MATLAB/R (R output, N = 53): F(%g, %g) = %.2f, p = %.2f, ηp² = %.3f\n\n", pp.df1, pp.df2, pp.F, pp.p, pp.eta_p2)
end
println(io, "Exploratory 2 × 2 (probe × proximity) on onset latency")
for r in eachrow(twobytwo_julia)
    pp = PAPER_2X2[(r.component, r.effect)]
    m = twobytwo_matlab === nothing ? nothing :
        twobytwo_matlab[(twobytwo_matlab.component .== r.component) .& (twobytwo_matlab.effect .== r.effect), :][1, :]
    @printf(io, "  %-6s %-15s Julia F(1, %d) = %.2f, p = %.3f", r.component, r.effect, r.df2, r.F, r.p)
    m === nothing || @printf(io, " | MATLAB F(1, %d) = %.2f, p = %.3f", m.df2, m.F, m.p)
    if alt2 !== nothing
        a = alt2[(alt2.component .== r.component) .& (alt2.effect .== r.effect), :]
        nrow(a) > 0 && @printf(io, " | ICs as listed F(1, %g) = %.2f, p = %.3f", a.df2[1], a.F[1], a.p[1])
    end
    @printf(io, " | MATLAB/R F = %.2f, p = %.3f\n", pp.F, pp.p)
end
if nrow(cmp_rows) > 0
    println(io, "\nPer-participant measures, Julia vs. MATLAB (all bins):")
    for g in groupby(cmp_rows, [:component, :measure])
        @printf(io, "  %-6s %-15s identical %d/%d, max |diff| = %.3f, mean |diff| = %.4f, NaN mismatches = %d\n",
                g.component[1], g.measure[1], sum(g.n_identical), sum(g.n_both),
                maximum(filter(!isnan, g.max_abs_diff); init=0.0),
                mean(filter(!isnan, g.mean_abs_diff)), sum(g.nan_mismatch))
    end
end

# --- 4. Grand averages -------------------------------------------------------------
gavg = combine(groupby(ga, [:component, :bin, :filter, :time_ms]), :amplitude => mean => :amplitude)
CSV.write(joinpath(RESULTS, "grand_averages.csv"), gavg)
ref_ga = joinpath(REFERENCE, "matlab_grand_average_N2_FN400.csv")
if isfile(ref_ga) && length(subjects) == length(ALL_IDS)   # MATLAB grand averages: all 53
    m = CSV.read(ref_ga, DataFrame)
    m.filter = String.(m.filter); m.component = String.(m.component)
    j = innerjoin(rename(gavg, :amplitude => :julia), rename(select(m, Not(:bin_label)), :amplitude => :matlab),
                  on=[:component, :bin, :filter, :time_ms])
    s = combine(groupby(j, [:component, :filter]), [:julia, :matlab] => ((a, b) -> maximum(abs.(a .- b))) => :max_abs_diff_uV,
                [:julia, :matlab] => ((a, b) -> cor(a, b)) => :r)
    CSV.write(joinpath(RESULTS, "grand_average_comparison.csv"), s)
    println(io, "\nGrand averages (bins 46–56), Julia vs. MATLAB:")
    for r in eachrow(s)
        @printf(io, "  %-6s %-7s max |diff| = %.3f µV, r = %.4f\n", r.component, r.filter, r.max_abs_diff_uV, r.r)
    end
end
tc = joinpath(RESULTS, "trial_counts_comparison.csv")
if isfile(tc) && isempty(EXCLUDED_IDS)   # with exclusions, see script 04 (per participant)
    c = CSV.read(tc, DataFrame)
    @printf(io, "\nAccepted trials (N2/FN400 set, bin 17 = all probes): Julia %d, ERPLAB %d\n",
            c.julia[c.bin .== 17][1], c.accepted_total_53[c.bin .== 17][1])
end
txt = String(take!(io))
write(joinpath(RESULTS, "paper_comparison.txt"), txt)
print(txt)

# Figures 5 and 6 of the paper (Figures_N2_FN400_probe.m): change − no-change
# difference waves (20 Hz low-pass), N2 (Fz, FCz) on top and FN400 (Fz, Cz) below,
# with the measurement windows as arrows. Same axes, colours, fonts and size
# (900 × 650 px, 300 dpi) as the MATLAB figures.
lp = gavg[gavg.filter .== "lpfilt", :]
const FIG_ARROW_Y = -2.5
function paper_panel(comp, bins_labels, colors, win, chlabel, winlabel; bottom=false, legend=false)
    p = plot(xlims=(-200, 800), ylims=(-3, 1.5), yticks=-3:1:1, xticks=-200:200:800,
             framestyle=:axes, grid=false, tick_direction=:out,
             tickfontsize=12, guidefontsize=14, legendfontsize=12,
             ylabel="Change - no-change (µV)", xlabel=bottom ? "Time from probe onset (ms)" : "",
             legend=legend ? :bottomright : false, foreground_color_legend=nothing,
             background_color_legend=nothing, left_margin=6Plots.mm,
             bottom_margin=(bottom ? 4 : 1) * Plots.mm)
    vline!(p, [0]; color=:black, ls=:dash, lw=1, label="")
    hline!(p, [0]; color=:black, lw=1, label="")
    for ((b, lab), col) in zip(bins_labels, colors)
        d = lp[(lp.component .== comp) .& (lp.bin .== b), :]
        plot!(p, d.time_ms, d.amplitude; label=lab, color=col, lw=2)
    end
    # Window arrow and its label
    plot!(p, collect(win), [FIG_ARROW_Y, FIG_ARROW_Y]; color=:black, lw=1.5, label="")
    scatter!(p, [win[1]], [FIG_ARROW_Y]; marker=:ltriangle, color=:black, ms=4, msw=0, label="")
    scatter!(p, [win[2]], [FIG_ARROW_Y]; marker=:rtriangle, color=:black, ms=4, msw=0, label="")
    annotate!(p, sum(win) / 2, FIG_ARROW_Y - 0.08, text(winlabel, 12, :center, :top))
    # Channel names, top left
    annotate!(p, -180, 1.45, text(chlabel, 12, :left, :top))
    return p
end
function paper_figure(bins_labels, colors, name)
    f = plot(paper_panel("N2", bins_labels, colors, N2_WINDOW, "Fz, FCz", "N2"; legend=true),
             paper_panel("FN400", bins_labels, colors, FN400_WINDOW, "Fz, Cz", "FN400"; bottom=true),
             layout=(2, 1), size=(900, 650), dpi=300, background_color=:white)
    savefig(f, joinpath(RESULTS, name * ".png"))
    savefig(f, joinpath(RESULTS, name * ".pdf"))
end
# Figure 5: 55 = 2-UR (blue), 56 = 3-UR (orange), 53 = 3-RSP (gold), 52 = 3-RWP (purple)
stim = Dict(l => b for (b, l) in STIMULUS_BINS)
paper_figure([stim[c] => c for c in CONDITIONS], COLORS, "Figure5_N2_FN400_StimulusType")
# Figure 6: 50 = repeated (salmon), 51 = unrepeated (teal)
paper_figure([50 => "Repeated probe", 51 => "Unrepeated probe"],
             ["#F8766D", "#00BFC4"], "Figure6_N2_FN400_ProbeType")

n2pc = DataFrame()
for s in subjects
    B = n2pc_bins(erps, s)
    for (b, lab) in (1 => "Unrepeated", 2 => "Repeated", 5 => "Overall")
        append!(n2pc, DataFrame(subject=s, bin=b, label=lab, time_ms=times_ms,
                                amplitude=lowpass_erp(B[b]).wave))
    end
end
gn = combine(groupby(n2pc, [:bin, :label, :time_ms]), :amplitude => mean => :amplitude)
CSV.write(joinpath(RESULTS, "n2pc_probe_grand_average.csv"), gn)
# Probe N2pc (not in the paper: no cluster), in the same style as Figures 5 and 6.
p = plot(xlims=(-200, 800), xticks=-200:200:800, framestyle=:axes, grid=false,
         tick_direction=:out, tickfontsize=12, guidefontsize=14, legendfontsize=12,
         xlabel="Time from probe onset (ms)", ylabel="Contra - ipsi, change - no-change (µV)",
         legend=:bottomright, foreground_color_legend=nothing, background_color_legend=nothing,
         size=(900, 350), dpi=300, left_margin=6Plots.mm, bottom_margin=4Plots.mm,
         background_color=:white)
vline!(p, [0]; color=:black, ls=:dash, lw=1, label="")
hline!(p, [0]; color=:black, lw=1, label="")
for (lab, col) in (("Repeated", "#F8766D"), ("Unrepeated", "#00BFC4"), ("Overall", "#000000"))
    d = gn[gn.label .== lab, :]
    plot!(p, d.time_ms, d.amplitude; label=lab == "Overall" ? "All probes" : lab * " probe",
          color=col, lw=2)
end
savefig(p, joinpath(RESULTS, "n2pc_probe.png"))
@info "Done. See results/paper_comparison.txt"
