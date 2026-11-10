# Probe-locked ERPs in visual working memory: a Julia reanalysis

A Julia pipeline that reproduces the probe-locked analyses (N2pc, N2, FN400) of my open
EEG dataset [ds006979 on OpenNeuro](https://doi.org/10.18112/openneuro.ds006979.v1.1.4)
(version 1.1.4; Ramzaoui & Beck, *Scientific Reports*, Registered Report; Stage 1
accepted, Stage 2 under review), from the BIDS data and the original ICA weights.

It is the companion of [vwm-erp-julia](https://github.com/HananeRam/vwm-erp-julia),
which reproduces the stimulus-locked analyses (N2pc, CDA, negative slow wave, alpha), and
reuses its tested preprocessing code (`src/VWMPipeline.jl`).

In the task, a single-item probe appeared 1 s after the memory array; participants
reported whether its color had changed. Probe-locked ERPs index the comparison of the
probe with memory:

- **N2** (Fz, FCz) and **FN400** (Fz, Cz): change minus no-change difference waves, by
  stimulus type (2-UR, 3-UR, 3-RWP, 3-RSP) and, in the grouped arrays, by probe type
  (repeated vs. unrepeated item);
- **probe N2pc** (14 posterior channels): contralateral minus ipsilateral, change
  minus no-change.

All analyses use correct-response trials only. The original analysis used 53 participants;
participant 3 is excluded here (N = 52, see Notes).

## What is reproduced

| MATLAB / R step | What it does | Julia |
| --- | --- | --- |
| Scripts #1–#4 | Downsampling to 256 Hz, 0.1 Hz high-pass, original ICA weights and removed components, average-mastoid reference | `src/VWMPipeline.jl` (from vwm-erp-julia), `scripts/01_process_probe.jl` |
| Script #5, `BDF_probe.txt` | Probe bins: the probe code must be followed by 120 (correct), memory-array codes 51–58 ignored; epochs −200 to 800 ms, baseline −200–0 ms | `read_bdf`, `label_probe_trials`, `epoch` |
| Script #6c | N2/FN400 artifacts: saccades (step function on HEOG, −50–800 ms, 48 µV), ±200 µV and 150 µV moving-window peak-to-peak on Fz, FCz, Cz | `scripts/artifact_parameters_probe.csv` (from the `AR_*_N2_FN400_probe.xlsx` files) |
| Script #6b | Probe N2pc artifacts: saccades during the probe (−50–350 ms, 16 µV) and in the epoch (−200–400 ms, 64 µV), ±200 µV and 150 µV peak-to-peak on the 14 posterior channels and both HEOGs | same file (`AR_*_N2pc_probe*.xlsx`) |
| Script #7c, bin operations | Change − no-change waves (bins 46–49), trial-weighted averages (50–54), 2-UR and 3-UR (55–56); 20 Hz low-pass after the bin operations | `scripts/bin_operations.jl` (`bin_minus`, `wavgbin`, with ERPLAB's trial counts) |
| Script #7b, bin operations | Probe N2pc: contra/ipsi with equal hemifield weights, contra − ipsi, change − no-change, weighted averages | `n2pc_bins`, `contra_ipsi_probe` |
| ERP_Cluster_Analyses Steps 2–4 | `clustGND`: one-sample cluster permutation tests, 100–800 ms, 1000 permutations | `clustgnd`, `scripts/02_cluster_tests.jl` |
| Script #11c | Mean amplitude (N2 191.4–300 ms, FN400 300–464.8 ms, no low-pass); negative local peak latency (3-point neighbourhood, low-passed); 50% fractional peak onset latency (N2 91.4–300 ms, FN400 200–464.8 ms) | `meanbl`, `peaklatbl`, `fpeaklat`, ported from ERPLAB's `geterpvalues.m` and `localpeak.m` |
| R Steps 3–7 | Repeated-measures ANOVAs of Models 4–9 (with Greenhouse–Geisser correction when Mauchly's test is significant) and the exploratory 2 × 2 on onset latency | `rm_anova_gg`, `rm_anova_2x2`, `scripts/03_measures_and_statistics.jl` |

## Get the data

```bash
bash scripts/download_data.sh
```

This downloads the EEG files of the 53 participants (about 37 GB) and their ICA weights
(about 5.5 GB) into `data/ds006979`. If vwm-erp-julia already downloaded some of
them, use the same folder: `DS006979=/path/to/data/ds006979 bash scripts/download_data.sh`,
and run the Julia scripts with the same `DS006979` variable.

## Run

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'     # once
julia --project=. scripts/01_process_probe.jl           # preprocessing, artifacts, averages (~1–2 min per participant)
julia --project=. scripts/02_cluster_tests.jl           # cluster tests (time windows)
julia --project=. scripts/03_measures_and_statistics.jl # measures, ANOVAs, comparison, figures
```

## Outputs (`results/`)

| File | Content |
| --- | --- |
| `probe_trials.csv` | One row per probe: condition, cue, change, probe type, next event, artifact flags, kept trials per analysis |
| `probe_qc.csv`, `probe_bin_counts.csv` | Trials binned, rejected (per detector) and kept; accepted trials per bin and participant |
| `trial_counts_comparison.csv` | Accepted trials per bin (sum over participants) vs. ERPLAB |
| `probe_erps.csv` | Participant averages per bin (no low-pass) |
| `cluster_tests.txt`, `cluster_tests.csv`, `cluster_tvalues.csv` | Cluster tests, Julia vs. MATLAB |
| `measures.csv` | Mean amplitude, peak latency and onset latency per participant and bin |
| `measures_julia_vs_matlab.csv`, `measures_comparison.csv` | The same, next to the MATLAB values, and agreement per bin |
| `anova.csv`, `anova_2x2_onset.csv`, `descriptives.csv` | Models 4–9 |
| `paper_comparison.txt` | Julia, the same statistics on the MATLAB measures, and the published values |
| `grand_averages.csv`, `grand_average_comparison.csv` | Grand averages and their agreement with the MATLAB grand averages |
| `Figure5_N2_FN400_StimulusType.png/.pdf`, `Figure6_N2_FN400_ProbeType.png/.pdf` | Figures 5 and 6 of the paper, same design as `Figures_N2_FN400_probe.m` (20 Hz low-pass, 300 dpi) |
| `n2pc_probe.png` | Probe N2pc grand average, same style (not in the paper) |

`reference/` holds the MATLAB outputs used for these comparisons: the measurement files
of Script #11c, the cluster test results (time points and p-values), the grand averages and
trial counts of Script #9c.

## Validation

Run on ds006979 version 1.1.4 with the original ICA weights, for the 52 participants of the
probe-locked analyses (participant 3 excluded, see Notes). "MATLAB" is the same statistics
applied to the MATLAB measurement files of the same 52 participants (Script #11c output).

| Step | Julia vs. MATLAB/ERPLAB |
| --- | --- |
| Bin assignment | 56,764 vs. 56,763 binned probes (1 trial) |
| Artifact rejection | accepted trials identical: 54,973 (N2/FN400) and 50,248 (probe N2pc), every participant |
| Participant averages (bins 1, 2, 5, 6, 38–45) | r = 1.0000 for 51 of 52 participants (sub-22: 0.9998); median RMS difference 0.001 µV |
| Cluster tests (N = 52) | N2 191.4–449.2 ms, FN400 195.3–464.8 ms, no probe N2pc cluster: identical to the paper (N = 53) |
| Peak latencies | 1,040 / 1,040 identical |
| Onset latencies | 950 / 955 identical (the others differ by one sample, 3.9 ms) |
| Mean amplitudes | mean difference 0.0013–0.0017 µV, max 0.065 µV (MATLAB rounds to 3 decimals) |

| Model | Julia (N = 52) | MATLAB (N = 52) | Paper (N = 53) |
| --- | --- | --- | --- |
| N2 amplitude, stimulus type (5.a) | F(3, 153) = 1.21, p = .308 | F(3, 153) = 1.20, p = .310 | F(3, 156) = 1.24, p = .30 |
| N2 peak latency, stimulus type (5.b) | F(3, 153) = 0.90, p = .442 | F(3, 153) = 0.90, p = .442 | F(3, 156) = 1.12, p = .34 |
| FN400 amplitude, stimulus type (6.a) | F(2.35, 119.82) = 0.16, p = .884 | F(2.35, 119.82) = 0.16, p = .883 | F(2.38, 123.76) = 0.18, p = .87 |
| FN400 peak latency, stimulus type (6.b) | F(3, 153) = 0.85, p = .471 | F(3, 153) = 0.85, p = .471 | F(3, 156) = 0.83, p = .48 |
| N2 amplitude, probe type (8.a) | F(1, 51) = 0.51, p = .477 | F(1, 51) = 0.51, p = .476 | F(1, 52) = 0.60, p = .44 |
| N2 peak latency, probe type (8.b) | F(1, 51) = 3.18, p = .081 | F(1, 51) = 3.18, p = .081 | F(1, 52) = 2.40, p = .13 |
| FN400 amplitude, probe type (9.a) | F(1, 51) = 1.54, p = .220 | F(1, 51) = 1.55, p = .219 | F(1, 52) = 1.59, p = .21 |
| FN400 peak latency, probe type (9.b) | F(1, 51) = 0.23, p = .636 | F(1, 51) = 0.23, p = .636 | F(1, 52) = 0.22, p = .64 |
| N2 onset latency, probe type (4.b) | F(1, 38) = 0.56, p = .457 | F(1, 38) = 0.56, p = .457 | F(1, 39) = 0.78, p = .38 |
| FN400 onset latency, probe type (4.c) | F(1, 28) = 0.52, p = .475 | F(1, 28) = 0.51, p = .480 | F(1, 29) = 0.38, p = .54 |

The exploratory 2 × 2 on onset latency is identical to MATLAB for all six effects (e.g. FN400
stimulus F(1, 28) = 3.42, p = .075). No conclusion differs from the paper.

The measurement code alone was also checked on the MATLAB waveforms (the
`…ChangeNoChange_waves_ClusterInterest` files of the dataset derivatives): mean amplitude, peak latency and onset latency
are identical to the MATLAB measurement files for all 53 participants (3,095 / 3,095
values). Applied to the MATLAB measurement files of all 53 participants, the statistics code
reproduces every published value of Models 4–9 and of the 2 × 2.

Trial counts and participant averages were compared with the ERPLAB files of the dataset
derivatives (`<id>_N2_FN400_ProbeERPs_erp_ar.erp`, `<id>_N2pc_ProbeERPs_erp_ar.erp`); the
continuous data after Scripts #1 and #4 with the intermediate MATLAB datasets, which are not
shared.

## Notes

- **Bin assignment.** As in ERPLAB, a probe is binned only if the next event (skipping
  memory-array codes) is a correct response (120). `probe_qc.csv` also reports how many
  probes a more lenient rule (any 120 before the next trial) would accept.
- **Trial counts.** ERPLAB carries trial counts through bin operations: a difference
  bin gets the sum of the counts of its bins, and `wavgbin` weights by these counts.
  The Julia bin operations do the same, so the weighted averages (bins 50–54) match.
- **Cluster tests.** The cluster limits depend only on the data. The p-values depend on the
  random permutations, which cannot reproduce MATLAB's random stream, so they vary
  slightly; the clusters of the paper are far from the significance threshold.
- **Low-pass filter.** The 20 Hz Butterworth filter (pop_filterp) is applied to the
  averaged waveforms after the bin operations, as in Scripts #7b/#7c; mean amplitudes use
  the unfiltered waveforms, latencies the filtered ones.

- **Participant 3.** In this participant's Script #0 dataset, the three unused external
  channels (EXG1, EXG7, EXG8) were not removed, so in every later, position-based step
  (Scripts #1–#6) the EOG and mastoid channels were shifted by one: the "saccade" check
  ran on VEOG-lower − HEOG-left (and rejected blinks), and the mastoid reference was
  HEOG-right and Mastoid-left. Participant 3 is therefore excluded from the probe-locked
  analyses (`EXCLUDED_IDS` in `scripts/config.jl`). It is excluded from the stimulus-locked
  analyses with either channel assignment (residual HEOG with the shifted channels; more than
  50% rejected trials and residual HEOG of about 3 µV or more with the correct ones).

- **ICA components.** For participants 46, 49 and 69, Script #4 removed other components than
  `ICs_to_Remove.xlsx` lists (46: [1 2] instead of [1]; 49: [2] instead of [1]; 69: [1 7] instead of [1]), as identified
  by reconstructing the MATLAB ICA-corrected datasets from the uncorrected data. Julia removes the
  components actually removed (`ICS_APPLIED` in `scripts/config.jl`), so that it reproduces the
  published analysis.
- **Dataset version.** Use ds006979 version 1.1.4 or later: it contains complete EEG files for
  sub-10 and sub-49, the Script #0 recording of sub-04, and the recordings used in the analysis
  for sub-46, sub-47 and sub-64.

- **Event latencies after downsampling.** EEGLAB keeps fractional event latencies after
  `pop_resample` (half of the events fall on x.5 at 256 Hz), and `pop_epoch` starts each epoch
  at the lower sample. Julia does the same (`EVENT_ROUNDING` in `scripts/config.jl`); rounding
  these halves to the nearest even sample instead shifts about a quarter of the epochs by one
  sample (3.9 ms), which changes individual peak and onset latencies.

## License

MIT. The data are CC0 (OpenNeuro ds006979).
