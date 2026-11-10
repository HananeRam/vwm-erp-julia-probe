"""
    ProbeERP

Probe-locked ERPs (N2pc, N2, FN400), written to match the MATLAB/ERPLAB 10.02
scripts of the original analysis:

- `read_bdf`, `label_probe_trials`: bin assignment as ERPLAB's binlister with
  `BDF_probe.txt` (Script #5): the probe code must be followed by 120
  (correct response); memory-array codes 51–58 are ignored
- `ErpBin`, `bin_minus`, `wavgbin`, `contra_ipsi_probe`: ERPLAB bin
  operations, including the trial counts that ERPLAB carries through them
  (Script #7b, #7c)
- `meanbl`, `peaklatbl`, `fpeaklat`: `pop_geterpvalues` measurements
  (mean amplitude, local peak latency, fractional peak onset latency;
  Script #11c), ported from ERPLAB's `geterpvalues.m` and `localpeak.m`
- `clustgnd`: cluster-based permutation test against 0, ported from the
  Mass Univariate ERP Toolbox (`clustGND` → `clust_perm1`; ERP_Cluster_Analyses
  Steps 2–4)
- `rm_anova_gg`, `rm_anova_2x2`: repeated-measures ANOVAs as `rstatix::anova_test`
  (with Mauchly's test and Greenhouse–Geisser / Huynh–Feldt corrections)
"""
module ProbeERP

using Statistics, Random, Distributions, DataFrames, LinearAlgebra

export read_bdf, label_probe_trials, ErpBin, average_bin, bin_minus, wavgbin,
       contra_ipsi_probe, closest_index, meanbl, peaklatbl, fpeaklat, localpeak,
       clustgnd, rm_anova_gg, rm_anova_2x2, PROBE_IGNORE

# ---------------------------------------------------------------------------
# Bins (Script #5, BDF_probe.txt)
# ---------------------------------------------------------------------------

"Codes ignored by binlister for the probe-locked bins (memory-array codes; Script #5)."
const PROBE_IGNORE = Set(51:58)

"""
    read_bdf(path) -> Dict{Int,Tuple{String,Set{Int}}}

Read an ERPLAB bin descriptor file of the form

    bin 1
    Change 2-UR, correct
    .{11; 15; 27; 31}{120}

Returns bin number => (label, time-locking codes). Every bin of
`BDF_probe.txt` requires the next event to be 120; this is checked.
"""
function read_bdf(path::AbstractString)
    bins = Dict{Int,Tuple{String,Set{Int}}}()
    lines = [strip(l) for l in eachline(path)]
    i = 1
    while i <= length(lines)
        m = match(r"^bin\s*(\d+)"i, lines[i])
        if m === nothing
            i += 1
            continue
        end
        b = parse(Int, m.captures[1])
        j = i + 1
        while j <= length(lines) && isempty(lines[j])
            j += 1
        end
        label = String(lines[j])
        k = j + 1
        while k <= length(lines) && !startswith(lines[k], ".")
            k += 1
        end
        braces = [String(x.captures[1]) for x in eachmatch(r"\{([^}]*)\}", lines[k])]
        length(braces) == 2 && strip(braces[2]) == "120" ||
            error("bin $b: expected .{codes}{120}, found $(lines[k])")
        codes = Set(parse.(Int, filter(!isempty, strip.(split(braces[1], ';')))))
        bins[b] = (label, codes)
        i = k + 1
    end
    return bins
end

"""
    label_probe_trials(codes, samples, probe_codes) -> DataFrame

One row per probe event (codes in `probe_codes`, the 11–42 table of
VWMPipeline). `next_code` is the next event after the probe, skipping the
memory-array codes 51–58 as binlister's 'Ignore' option in Script #5 does
(other markers are not skipped). `correct` = `next_code == 120`, the
condition of every probe bin. `correct_lenient` = a 120 occurs before the next
probe or memory array, for comparison.
"""
function label_probe_trials(codes::AbstractVector{Int}, samples::AbstractVector{<:Integer},
                            probe_codes::AbstractDict)
    rows = NamedTuple[]
    for i in eachindex(codes)
        haskey(probe_codes, codes[i]) || continue
        cond, cue, change, probe = probe_codes[codes[i]]
        j = i + 1
        while j <= length(codes) && codes[j] in PROBE_IGNORE
            j += 1
        end
        nxt = j <= length(codes) ? codes[j] : -1
        lenient = false
        j = i + 1
        while j <= length(codes) && !haskey(probe_codes, codes[j]) && !(codes[j] in PROBE_IGNORE)
            codes[j] == 120 && (lenient = true; break)
            codes[j] == 121 && break
            j += 1
        end
        push!(rows, (code=codes[i], probe_sample=Int(samples[i]), condition=cond, cue=cue,
                     change=change, probe=probe, next_code=nxt, correct=nxt == 120,
                     correct_lenient=lenient))
    end
    return DataFrame(rows)
end

# ---------------------------------------------------------------------------
# Averages and bin operations (Scripts #7b, #7c)
# ---------------------------------------------------------------------------

"An averaged waveform and its number of accepted trials (ERP.ntrials.accepted)."
struct ErpBin
    wave::Vector{Float64}
    n::Int
end

"""
    average_bin(E, sel) -> ErpBin

Average of the trials `sel` of the epoch matrix `E` (pop_averager,
'Criterion','good'). A bin without trials is a flat line with n = 0, as in ERPLAB.
"""
function average_bin(E::AbstractMatrix, sel::AbstractVector{Bool})
    n = count(sel)
    n == 0 && return ErpBin(zeros(size(E, 2)), 0)
    return ErpBin(vec(mean(E[sel, :]; dims=1)), n)
end

"`nbin = binA - binB`: ERPLAB adds the trial counts of the bins used (binoperator.m)."
bin_minus(a::ErpBin, b::ErpBin) = ErpBin(a.wave .- b.wave, a.n + b.n)

"""
    wavgbin(bins...) -> ErpBin

ERPLAB `wavgbin`: average of the bins weighted by their trial counts; if all
counts are 0, the plain sum of the waveforms (as in wavgbin.m).
"""
function wavgbin(bins::ErpBin...)
    total = sum(b.n for b in bins)
    s = sum(b.n .* b.wave for b in bins)
    return ErpBin(total > 0 ? s ./ total : sum(b.wave for b in bins), total)
end

"""
    contra_ipsi_probe(left_cue, right_cue) -> ErpBin

Probe N2pc, `N2pc_ProbeERPs_BinOps_ContraIpsi.txt` and `…ContraIpsi2.txt` followed
by the 7-channel average of `N2pc_CDA_StimulusERPs_ChanOps_ClusterInterest.txt`.
`left_cue` and `right_cue` are (L = mean of the 7 left-hemisphere channels,
R = mean of the 7 right-hemisphere channels) for the left-cue and right-cue bins:

    contra = 0.5·(left-cue bin @ right channels) + 0.5·(right-cue bin @ left channels)
    ipsi   = 0.5·(left-cue bin @ left channels)  + 0.5·(right-cue bin @ right channels)

Trial counts follow ERPLAB: contra and ipsi each n_left + n_right, their
difference twice that.
"""
function contra_ipsi_probe(left_cue::NamedTuple, right_cue::NamedTuple)
    contra = 0.5 .* left_cue.R.wave .+ 0.5 .* right_cue.L.wave
    ipsi   = 0.5 .* left_cue.L.wave .+ 0.5 .* right_cue.R.wave
    n = left_cue.L.n + right_cue.L.n
    return ErpBin(contra .- ipsi, 2n)
end

# ---------------------------------------------------------------------------
# Measurements (Script #11c, pop_geterpvalues; ERPLAB geterpvalues.m/localpeak.m)
# ---------------------------------------------------------------------------

"Index of the sample closest to `t` (ERPLAB `closest`: first index on ties)."
closest_index(times_ms::AbstractVector, t::Real) = argmin(abs.(times_ms .- t))

"Baseline value (`blvalue2`): mean between the samples closest to the baseline limits."
function baseline_value(w::AbstractVector, times_ms::AbstractVector, bl)
    a = closest_index(times_ms, max(bl[1], minimum(times_ms)))
    b = closest_index(times_ms, min(bl[2], maximum(times_ms)))
    return mean(@view w[a:b])
end

_round3(x) = isnan(x) ? NaN : round(x; digits=3)      # 'Resolution', 3

"""
    meanbl(w, times_ms, window; baseline=(-200, 0)) -> Float64

'meanbl': baseline-corrected mean amplitude between the samples closest to the
window limits.
"""
function meanbl(w::AbstractVector, times_ms::AbstractVector, window; baseline=(-200.0, 0.0))
    x = w .- baseline_value(w, times_ms, baseline)
    a, b = closest_index(times_ms, window[1]), closest_index(times_ms, window[2])
    return _round3(mean(@view x[a:b]))
end

_mround(x) = round(Int, x, RoundNearestTiesAway)

"""
    localpeak(x, t; npoints, negative, replace_abs) -> (value, index)

ERPLAB `localpeak`: a local peak is a sample more extreme than the mean of the
`npoints` samples on each side and than its two immediate neighbours; searched
in x[npoints+1 : end-npoints]. Several local peaks: the most extreme (first on
ties); all of equal value: the median position. None: the absolute peak if
`replace_abs`, otherwise NaN.
"""
function localpeak(x::AbstractVector, t::AbstractVector; npoints::Int=3,
                   negative::Bool=true, replace_abs::Bool=true)
    n = length(x)
    a, b = npoints + 1, n - npoints
    seg = @view x[a:b]
    pabs = (negative ? argmin(seg) : argmax(seg)) + npoints
    vals, pos = Float64[], Int[]
    for i in a:b
        l = mean(@view x[i-npoints:i-1])
        r = mean(@view x[i+1:i+npoints])
        if negative
            (x[i] < l && x[i] < r && x[i] < x[i-1] && x[i] < x[i+1]) && (push!(vals, x[i]); push!(pos, i))
        else
            (x[i] > l && x[i] > r && x[i] > x[i-1] && x[i] > x[i+1]) && (push!(vals, x[i]); push!(pos, i))
        end
    end
    if !isempty(vals)
        u = unique(vals)
        if length(u) == 1 && length(vals) > 1
            p = _mround(median(pos))
            return (value=u[1], index=p, local_peak=true)
        elseif length(u) > 1
            k = negative ? argmin(vals) : argmax(vals)
            return (value=vals[k], index=pos[k], local_peak=true)
        else
            return (value=vals[1], index=pos[1], local_peak=true)
        end
    end
    replace_abs && return (value=x[pabs], index=pabs, local_peak=false)
    return (value=NaN, index=0, local_peak=false)
end

# Window ± npoints samples, as geterpvalues extracts it for peak measures.
function _peak_segment(w, times_ms, window, baseline, npoints)
    x = w .- baseline_value(w, times_ms, baseline)
    a, b = closest_index(times_ms, window[1]), closest_index(times_ms, window[2])
    (a - npoints < 1 || b + npoints > length(x)) &&
        error("measurement window too close to the epoch edges for a $npoints-point local peak")
    r = (a-npoints):(b+npoints)
    return x[r], times_ms[r]
end

"""
    peaklatbl(w, times_ms, window; npoints=3, negative=true, baseline=(-200, 0))

'peaklatbl' with 'Peakpolarity','negative', 'Neighborhood',3, 'PeakReplace',
'absolute': latency (ms) of the local peak, or of the absolute peak if none.
"""
function peaklatbl(w::AbstractVector, times_ms::AbstractVector, window; npoints::Int=3,
                   negative::Bool=true, baseline=(-200.0, 0.0))
    x, t = _peak_segment(w, times_ms, window, baseline, npoints)
    p = localpeak(x, t; npoints=npoints, negative=negative, replace_abs=true)
    return _round3(t[p.index])
end

"""
    fpeaklat(w, times_ms, window; frac=0.5, npoints=3, negative=true, baseline=(-200, 0))

'fpeaklat' with 'Afraction',0.5, 'PeakOnset',1, 'Fracreplace','NaN',
'Peakreplace','absolute': going back from the (local, else absolute) peak, the
latency of the first sample that is no more extreme than `frac` × peak. NaN if
that sample lies before the window (in the npoints padding) or is not
strictly beyond the target value, as in localpeak.m.
"""
function fpeaklat(w::AbstractVector, times_ms::AbstractVector, window; frac::Real=0.5,
                  npoints::Int=3, negative::Bool=true, baseline=(-200.0, 0.0))
    x, t = _peak_segment(w, times_ms, window, baseline, npoints)
    p = localpeak(x, t; npoints=npoints, negative=negative, replace_abs=true)
    n = length(x)
    target = p.value * frac
    # localpeak.m: candidate samples strictly on the other side of the target,
    # excluding the padding before the window (negative peaks) or after it (positive).
    match = negative ? Set(i for i in 1:n if x[i] > target && i >= npoints + 1) :
                       Set(i for i in 1:n if x[i] < target && i <= n - npoints)
    i = p.index
    found = false
    while i >= 1 && i <= n
        if (negative && x[i] >= target) || (!negative && x[i] <= target)
            found = true
            break
        end
        i -= 1
    end
    (found && i in match) || return NaN
    return _round3(t[i])
end

# ---------------------------------------------------------------------------
# Cluster-based permutation test (Mass Univariate ERP Toolbox, clust_perm1)
# ---------------------------------------------------------------------------

_tscores(D) = vec(mean(D; dims=1) ./ (std(D; dims=1) ./ sqrt(size(D, 1))))

function _runs(mask::AbstractVector{Bool})
    runs = UnitRange{Int}[]
    i = 1
    while i <= length(mask)
        if mask[i]
            j = i
            while j < length(mask) && mask[j+1]
                j += 1
            end
            push!(runs, i:j)
            i = j + 1
        else
            i += 1
        end
    end
    return runs
end

# Most negative cluster mass (sum of t) among clusters with t < thresh (thresh < 0); 0 if none.
_min_neg_mass(t, thresh) = minimum((sum(t[r]) for r in _runs(t .< thresh)); init=0.0)

"""
    clustgnd(D, times_ms; window=(100, 800), nperm=1000, thresh_p=0.05, rng) -> NamedTuple

One-sample (vs 0) cluster-based permutation test on the participants × time
matrix `D`, as `clustGND(GND, bin, 'time_wind',[100 800], 'chan_hood',[],
'n_perm',1000, 'alpha',0.05, 'thresh_p',0.05)` on a single (averaged) channel:

- time points with times_ms in [window] are tested;
- t-scores above tinv(1 − thresh_p/2) or below tinv(thresh_p/2) form positive
  and negative clusters; cluster mass = sum of t;
- null distribution: the most negative cluster mass after flipping the sign of
  each participant's wave at random (0 if no negative cluster);
- two-tailed p of a positive cluster = 2 · mean(null ≤ −mass), of a negative
  cluster 2 · mean(null ≤ mass).

The cluster limits depend only on the data; the p-values vary slightly with
the random permutations (MATLAB's random stream cannot be reproduced).
"""
function clustgnd(D::AbstractMatrix, times_ms::AbstractVector; window=(100.0, 800.0),
                  nperm::Integer=1000, thresh_p::Real=0.05, alpha::Real=0.05,
                  rng::AbstractRNG=MersenneTwister(2026))
    keep = findall(t -> window[1] <= t <= window[2], times_ms)
    X = D[:, keep]
    good = [all(isfinite, X[i, :]) for i in axes(X, 1)]
    X = X[good, :]
    n = size(X, 1)
    thresh = quantile(TDist(n - 1), thresh_p / 2)          # negative
    null = zeros(nperm)
    for p in 1:nperm
        s = [rand(rng) > 0.5 ? 1.0 : -1.0 for _ in 1:n]
        null[p] = _min_neg_mass(_tscores(X .* s), thresh)
    end
    t = _tscores(X)
    pvals = ones(length(t))
    rows = NamedTuple[]
    for (sign, mask) in ((:positive, t .> -thresh), (:negative, t .< thresh))
        for r in _runs(mask)
            mass = sum(t[r])
            p = sign === :positive ? 2 * mean(null .<= -mass) : 2 * mean(null .<= mass)
            pvals[r] .= p
            push!(rows, (sign=String(sign), start_ms=times_ms[keep][first(r)],
                         end_ms=times_ms[keep][last(r)], mass=mass, p=p))
        end
    end
    clusters = isempty(rows) ? DataFrame(sign=String[], start_ms=Float64[], end_ms=Float64[],
                                         mass=Float64[], p=Float64[]) : DataFrame(rows)
    clusters.significant = clusters.p .< alpha
    return (n=n, t_threshold=-thresh, times_ms=times_ms[keep], t=t, p=pvals,
            clusters=clusters, null=null)
end

# ---------------------------------------------------------------------------
# Repeated-measures ANOVAs (R: rstatix::anova_test)
# ---------------------------------------------------------------------------

_complete(wide, cols) = begin
    X = Matrix{Float64}(coalesce.(wide[:, cols], NaN))
    X[vec(all(!isnan, X; dims=2)), :]
end

"""
    rm_anova_gg(wide, conditions) -> NamedTuple

One-way repeated-measures ANOVA on participants × conditions (complete cases),
with partial η², Mauchly's test of sphericity, and Greenhouse–Geisser and
Huynh–Feldt corrected degrees of freedom and p-values (as reported by
rstatix::anova_test; for two conditions sphericity holds and ε = 1).
"""
function rm_anova_gg(wide::DataFrame, conditions)
    X = _complete(wide, conditions)
    n, k = size(X)
    grand = mean(X)
    ss_cond = n * sum((vec(mean(X; dims=1)) .- grand) .^ 2)
    ss_subj = k * sum((vec(mean(X; dims=2)) .- grand) .^ 2)
    ss_error = sum((X .- grand) .^ 2) - ss_cond - ss_subj
    df1, df2 = k - 1, (k - 1) * (n - 1)
    F = (ss_cond / df1) / (ss_error / df2)
    p = ccdf(FDist(df1, df2), F)
    gg, hf, W, pW = 1.0, 1.0, 1.0, 1.0
    if k > 2
        # Orthonormal contrasts (Helmert), covariance of the transformed scores.
        C = zeros(k, k - 1)
        for j in 1:(k-1)
            C[1:j, j] .= 1
            C[j+1, j] = -j
            C[:, j] ./= norm(C[:, j])
        end
        S = cov(X * C)
        q = k - 1
        gg = tr(S)^2 / (q * tr(S * S))
        hf = min(1.0, (n * q * gg - 2) / (q * (n - 1 - q * gg)))
        W = det(S) / (tr(S) / q)^q
        chi2 = -((n - 1) - (2q^2 + q + 2) / (6q)) * log(W)
        dfW = q * (q + 1) / 2 - 1
        pW = ccdf(Chisq(dfW), chi2)
    end
    return (n=n, F=F, df1=df1, df2=df2, p=p, eta_p2=ss_cond / (ss_cond + ss_error),
            mauchly_W=W, mauchly_p=pW, gg_eps=gg, gg_df1=gg * df1, gg_df2=gg * df2,
            gg_p=ccdf(FDist(gg * df1, gg * df2), F), hf_eps=hf,
            hf_p=ccdf(FDist(hf * df1, hf * df2), F))
end

"""
    rm_anova_2x2(wide, cells; names=("A", "B")) -> DataFrame

Two-way repeated-measures ANOVA with two levels per factor. `cells` are the
column names in the order (a1b1, a1b2, a2b1, a2b2). Each effect has 1 df, and
its F equals the squared one-sample t of the corresponding contrast.
"""
function rm_anova_2x2(wide::DataFrame, cells; names=("A", "B"))
    X = _complete(wide, collect(cells))
    n = size(X, 1)
    contrasts = [(names[1], (X[:, 1] .+ X[:, 2]) ./ 2 .- (X[:, 3] .+ X[:, 4]) ./ 2),
                 (names[2], (X[:, 1] .+ X[:, 3]) ./ 2 .- (X[:, 2] .+ X[:, 4]) ./ 2),
                 ("$(names[1]):$(names[2])", (X[:, 1] .- X[:, 2]) .- (X[:, 3] .- X[:, 4]))]
    rows = NamedTuple[]
    for (eff, d) in contrasts
        tval = mean(d) / (std(d) / sqrt(n))
        F = tval^2
        push!(rows, (effect=eff, n=n, F=F, df1=1, df2=n - 1, p=ccdf(FDist(1, n - 1), F),
                     eta_p2=F / (F + n - 1)))
    end
    return DataFrame(rows)
end

end # module
