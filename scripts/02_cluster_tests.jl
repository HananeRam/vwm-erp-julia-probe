# Cluster-based permutation tests that defined the probe-locked time windows
# (ERP_Cluster_Analyses Steps 2–4). Run after 01:
#   julia --project=. scripts/02_cluster_tests.jl
#
# - N2:    bin 54 (change − no-change, all repeated-condition probes), Fz/FCz
# - FN400: bin 54, Fz/Cz
# - N2pc:  bin 5 of the probe N2pc (contra − ipsi, change − no-change), 7-channel cluster
# One-sample tests against 0, 100–800 ms, two-tailed, 1000 permutations.
# Saves results/cluster_tests.csv, results/cluster_tests.txt and
# results/cluster_tvalues.csv (t and p per time point, as clustGND's output file).

include(joinpath(@__DIR__, "config.jl"))
include(joinpath(@__DIR__, "..", "src", "VWMPipeline.jl"))
include(joinpath(@__DIR__, "..", "src", "ProbeERP.jl"))
using .VWMPipeline: lowpass_epochs!
using .ProbeERP
using Random, Printf, Distributions
include(joinpath(@__DIR__, "bin_operations.jl"))

erps, times_ms, subjects = load_probe_erps()
@info "$(length(subjects)) participants"

waves = Dict(
    "N2"    => reduce(vcat, permutedims(n2fn400_bins(erps, s, "N2")[54].wave) for s in subjects),
    "FN400" => reduce(vcat, permutedims(n2fn400_bins(erps, s, "FN400")[54].wave) for s in subjects),
    "N2pc"  => reduce(vcat, permutedims(n2pc_bins(erps, s)[5].wave) for s in subjects))

# Time points and p-values of the MATLAB results files (reference/cluster_tests).
function matlab_clusters(name)
    f = joinpath(REFERENCE, "cluster_tests", "results_clustGND_$(name)_probe.txt")
    isfile(f) || return nothing
    L = readlines(f)
    nums(s) = [x for x in (tryparse(Float64, w) for w in split(s)) if x !== nothing]
    t, p = nums(L[1]), nums(L[2])
    sig = p .< 0.05
    runs = UnitRange{Int}[]
    i = 1
    while i <= length(sig)
        if sig[i]
            j = i
            while j < length(sig) && sig[j+1]; j += 1; end
            push!(runs, i:j); i = j + 1
        else
            i += 1
        end
    end
    return [(t[first(r)], t[last(r)]) for r in runs]
end

rows = DataFrame()
tvals = DataFrame()
io = IOBuffer()
println(io, "Cluster-based permutation tests, probe-locked (N = $(length(subjects)); ",
        "$(N_PERMUTATIONS) permutations, seed $(CLUSTER_SEED))\n")
for name in ("N2", "FN400", "N2pc")
    res = clustgnd(waves[name], times_ms; window=CLUSTER_WINDOW, nperm=N_PERMUTATIONS,
                   rng=MersenneTwister(CLUSTER_SEED))
    sig = res.clusters[res.clusters.significant, :]
    for r in eachrow(res.clusters)
        push!(rows, (component=name, sign=r.sign, start_ms=r.start_ms, end_ms=r.end_ms,
                     mass=r.mass, p=r.p, significant=r.significant))
    end
    append!(tvals, DataFrame(component=name, time_ms=res.times_ms, t=res.t, p=res.p))
    ref = matlab_clusters(name)
    @printf(io, "%-6s Julia: %s\n", name,
            nrow(sig) == 0 ? "no significant cluster" :
            join([@sprintf("%.1f–%.1f ms (%s, p = %.3f)", r.start_ms, r.end_ms, r.sign, r.p)
                  for r in eachrow(sig)], "; "))
    ref === nothing || @printf(io, "%-6s MATLAB (N = 53): %s\n", "",
        isempty(ref) ? "no significant cluster" :
        join([@sprintf("%.1f–%.1f ms", a, b) for (a, b) in ref], "; "))
    println(io)
end
println(io, "Cluster limits depend only on the data (t threshold |t| > ",
        @sprintf("%.3f", quantile(TDist(length(subjects) - 1), 0.975)),
        "); p-values vary slightly with the random permutations.")
println(io, "Paper: N2 191.4–449.2 ms and FN400 195.3–464.8 ms (split at 300 ms into the ",
        "N2 191.4–300 ms and FN400 300–464.8 ms windows); no probe-locked N2pc cluster.")

mkpath(RESULTS)
CSV.write(joinpath(RESULTS, "cluster_tests.csv"), rows)
CSV.write(joinpath(RESULTS, "cluster_tvalues.csv"), tvals)
txt = String(take!(io))
write(joinpath(RESULTS, "cluster_tests.txt"), txt)
print(txt)
