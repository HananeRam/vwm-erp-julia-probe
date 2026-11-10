# ERPLAB bin operations of Scripts #7b and #7c, applied to the participant averages
# saved by script 01 (results/probe_erps.csv). Included by scripts 02 and 03.

using CSV, DataFrames, Statistics

"Read results/probe_erps.csv → (Dict (subject, set, bin, signal) => ErpBin, times_ms, subjects)."
function load_probe_erps(path=joinpath(RESULTS, "probe_erps.csv"))
    df = CSV.read(path, DataFrame)
    df.subject = String.(df.subject); df.set = String.(df.set); df.signal = String.(df.signal)
    sort!(df, [:subject, :set, :bin, :signal, :time_ms])
    erps = Dict{Tuple{String,String,Int,String},ErpBin}()
    for g in groupby(df, [:subject, :set, :bin, :signal])
        erps[(g.subject[1], g.set[1], g.bin[1], g.signal[1])] = ErpBin(Vector{Float64}(g.amplitude), g.n_trials[1])
    end
    times_ms = sort(unique(df.time_ms))
    # Only the participants of the analysis (config.jl: SUBJECTS, without EXCLUDED_IDS).
    return erps, times_ms, sort(intersect(unique(df.subject), SUBJECTS))
end

"20 Hz Butterworth low-pass of an averaged waveform (pop_filterp, order 8)."
lowpass_erp(b::ErpBin; fs=TARGET_FS) =
    ErpBin(vec(lowpass_epochs!(reshape(copy(b.wave), 1, :), fs; cutoff=LOWPASS.cutoff,
                               order=LOWPASS.order)), b.n)

"""
    n2fn400_bins(erps, sub, signal; lowpass=false) -> Dict{Int,ErpBin}

Bins 46–56 of `N2_FN400_ProbeERPs_BinOps_ChangeNoChange.txt` and `…ChangeNoChange2.txt`
for `signal` "N2" or "FN400" (change minus no-change difference waves):

    46 = 38 − 42 (repeated, 3-RWP)       47 = 39 − 43 (repeated, 3-RSP)
    48 = 40 − 44 (unrepeated, 3-RWP)     49 = 41 − 45 (unrepeated, 3-RSP)
    50 = wavg(46, 47) repeated           51 = wavg(48, 49) unrepeated
    52 = wavg(46, 48) 3-RWP              53 = wavg(47, 49) 3-RSP
    54 = wavg(46, 47, 48, 49) overall    55 = 1 − 5 (2-UR)   56 = 2 − 6 (3-UR)

With `lowpass=true`, every bin is then low-passed (Script #7c filters after the
bin and channel operations).
"""
function n2fn400_bins(erps, sub, signal; lowpass::Bool=false)
    g(b) = erps[(sub, "n2_fn400", b, signal)]
    B = Dict{Int,ErpBin}()
    B[46] = bin_minus(g(38), g(42)); B[47] = bin_minus(g(39), g(43))
    B[48] = bin_minus(g(40), g(44)); B[49] = bin_minus(g(41), g(45))
    B[50] = wavgbin(B[46], B[47]);   B[51] = wavgbin(B[48], B[49])
    B[52] = wavgbin(B[46], B[48]);   B[53] = wavgbin(B[47], B[49])
    B[54] = wavgbin(B[46], B[47], B[48], B[49])
    B[55] = bin_minus(g(1), g(5));   B[56] = bin_minus(g(2), g(6))
    lowpass && (B = Dict(k => lowpass_erp(v) for (k, v) in B))
    return B
end

"""
    n2pc_bins(erps, sub) -> Dict{Int,ErpBin}

Probe N2pc (Script #7b): contra − ipsi of the 7-channel cluster
(`…ContraIpsi.txt`, `…ContraIpsi2.txt`), change − no-change (`…ContraIpsi3.txt`),
weighted averages (`…ContraIpsi4.txt`):

    1 = unrepeated (avg 3-RWP/3-RSP)   2 = repeated (avg 3-RWP/3-RSP)
    3 = 3-RWP (avg probes)             4 = 3-RSP (avg probes)
    5 = overall mean (tested in ERP_Cluster_Analyses Step 2)

Left-cue and right-cue BDF bins of each condition (Script #5):
unrepeated change 3-RWP 46/50, 3-RSP 47/51; repeated change 3-RWP 48/52,
3-RSP 49/53; the same with no-change: 54/58, 55/59, 56/60, 57/61.
"""
function n2pc_bins(erps, sub)
    side(b) = (L=erps[(sub, "n2pc", b, "left")], R=erps[(sub, "n2pc", b, "right")])
    ci(bl, br) = contra_ipsi_probe(side(bl), side(br))
    # ContraIpsi2 order: 1 unrep change RWP, 2 unrep change RSP, 3 rep change RWP,
    # 4 rep change RSP, 5–8 the same for no-change.
    C = [ci(46, 50), ci(47, 51), ci(48, 52), ci(49, 53),
         ci(54, 58), ci(55, 59), ci(56, 60), ci(57, 61)]
    D = [bin_minus(C[i], C[i+4]) for i in 1:4]      # ContraIpsi3: change − no-change
    return Dict(1 => wavgbin(D[1], D[2]), 2 => wavgbin(D[3], D[4]),
                3 => wavgbin(D[1], D[3]), 4 => wavgbin(D[2], D[4]),
                5 => wavgbin(D[1], D[2], D[3], D[4]))
end
