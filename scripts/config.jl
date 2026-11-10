# Analysis settings for the probe-locked ERPs, aligned with the MATLAB/EEGLAB/ERPLAB
# scripts (#1–#5, #6b/#6c, #7b/#7c, #11c, ERP_Cluster_Analyses Steps 2–4).
# Edit here, not in the pipeline code.

# --- Data ------------------------------------------------------------------
# Set the DS006979 environment variable, or put the dataset in data/ds006979.
const BIDS_ROOT = get(ENV, "DS006979", joinpath(@__DIR__, "..", "data", "ds006979"))
const TASK = "visualworkingmemorytask"
const RESULTS = joinpath(@__DIR__, "..", "results")
const REFERENCE = joinpath(@__DIR__, "..", "reference")      # MATLAB outputs, for validation

# --- Participants ----------------------------------------------------------
# The MATLAB probe-locked analyses used 53 participants (SUB list in Scripts #6b–#11c).
const ALL_IDS = [1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                 16, 17, 18, 19, 21, 22, 23, 24, 25, 28, 30,
                 32, 33, 34, 35, 36, 39, 40, 41, 43, 46, 47, 48,
                 49, 50, 51, 52, 53, 54, 58, 59, 61, 62, 64,
                 69, 70, 72, 73, 74, 76]
# Participant 3 is excluded from the probe-locked analyses: in its Script #0 dataset
# (3.set) the three unused external channels (EXG1, EXG7, EXG8) were not removed, so
# from position 65 on the channels were shifted by one in every later (position-based)
# step: wrong EOG channels for artifact detection, HEOG-right and Mastoid-left as the
# mastoid reference, and ICA weights computed on these channels. Participant 3 is
# excluded from the stimulus-locked analyses with either channel assignment.
# Set EXCLUDED_IDS = Int[] to run all 53.
const EXCLUDED_IDS = [3]
const SUBJECT_IDS = setdiff(ALL_IDS, EXCLUDED_IDS)
const SUBJECTS = ["sub-" * lpad(string(i), 2, '0') for i in SUBJECT_IDS]

# Script #0 swapped response codes 120/121 for participants 21 and 22. The BIDS
# files were exported after that recoding, so nothing to swap here.
const SWAP_RESPONSE_IDS = Int[]

# --- Preprocessing (Scripts #1–#4), as in vwm-erp-julia --------------------------
const TARGET_FS = 256.0                     # Hz (pop_resample in Script #1)
const HIGHPASS = (cutoff=0.1, order=2)      # Hz; ERPLAB order (12 dB/oct)
const LOWPASS  = (cutoff=20.0, order=8)     # Hz; ERPLAB order (48 dB/oct), pop_filterp in #7b/#7c

# Blink correction: the original ICA (EEGLAB weights in the dataset derivatives,
# <id>_allERPs_optimized_weights.set, searched under ICA_WEIGHTS_DIR) and the
# components removed in ICs_to_Remove.csv (Script #4).
const ICA_WEIGHTS_DIR = get(ENV, "ICA_WEIGHTS",
                            joinpath(BIDS_ROOT, "derivatives", "ReplicationPeterson2015"))
const ICS_TO_REMOVE = joinpath(@__DIR__, "ICs_to_Remove.csv")
# Components actually removed by Script #4 where they differ from ICs_to_Remove.xlsx,
# identified by reconstructing the MATLAB ..._ica_corr_reref.set files from the
# uncorrected data and the weights (residual < 1e-4): to reproduce the published analysis, Julia removes these.
const ICS_APPLIED = Dict(46 => [1, 2], 49 => [2], 69 => [1, 7])

# Event latency after downsampling 512 → 256 Hz: EEGLAB keeps (lat − 1)/2 + 1, which
# is fractional (x.5) for half of the events, and pop_epoch starts the epoch at the
# lower sample (floor). With floor, sub-01's averages equal ERPLAB's (r = 1.0000,
# RMS difference 0.001 µV, vs. 0.27 µV with round-half-to-even).
const EVENT_ROUNDING = Dict("nearest" => RoundNearest, "ties_away" => RoundNearestTiesAway,
                            "floor" => RoundDown, "ceil" => RoundUp)[get(ENV, "EVENT_ROUNDING", "floor")]

# --- Epochs and bins (Script #5) ---------------------------------------------------
# Locked to the probe (codes 11–42): −200 to 800 ms → −199.2 to 796.9 ms at 256 Hz.
const EPOCH    = (tmin=-0.2, tmax=0.8)      # seconds
const BASELINE = (-0.2, 0.0)
const BDF_PROBE = joinpath(@__DIR__, "BDF_probe.txt")

# --- Channels ----------------------------------------------------------------------
# N2pc cluster pairs (left, right), channels 21–27 and 58–64 of the MATLAB datasets.
const CLUSTER = [("P3", "P4"), ("P5", "P6"), ("P7", "P8"), ("P9", "P10"),
                 ("PO7", "PO8"), ("PO3", "PO4"), ("O1", "O2")]
# N2 = mean of Fz and FCz (channels 38, 47); FN400 = mean of Fz and Cz (38, 48)
# (N2_FN400_ProbeERPs_ChanOps_ClusterInterest.txt).
const N2_CHANNELS    = ["Fz", "FCz"]
const FN400_CHANNELS = ["Fz", "Cz"]
# Channels kept after ICA. HEOG = uncorrected bipolar HEOG (channel 70, formed
# before ICA); HEOG_corr = corrected bipolar HEOG (channel 68, after ICA).
const KEEP_CHANNELS = unique(vcat([[l, r] for (l, r) in CLUSTER]..., ["Fz", "FCz", "Cz"],
                                  ["HEOG", "HEOG_corr"]))

# --- Artifact detection (Scripts #6b and #6c) ----------------------------------------
# One row per participant, set and detector, converted from the AR_*_probe.xlsx
# spreadsheets: n2_fn400 = Script #6c (N2, FN400), n2pc = Script #6b (probe N2pc).
const ARTIFACT_PARAMETERS = joinpath(@__DIR__, "artifact_parameters_probe.csv")
const ARTIFACT_SETS = ("n2_fn400", "n2pc")

# --- Measurements (Script #11c) ------------------------------------------------------
# Mean amplitude ('meanbl') on the ERPs without the 20 Hz low-pass; peak latency
# ('peaklatbl') and 50% fractional peak onset latency ('fpeaklat') on the low-passed
# ERPs; negative peaks, 3-point neighbourhood, baseline −200–0 ms.
const N2_WINDOW          = (191.4, 300.0)   # ms
const FN400_WINDOW       = (300.0, 464.8)
const N2_ONSET_WINDOW    = (91.4, 300.0)
const FN400_ONSET_WINDOW = (200.0, 464.8)
const PEAK_NEIGHBORHOOD  = 3
const ONSET_FRACTION     = 0.5

# --- Cluster tests (ERP_Cluster_Analyses Steps 2–4) ----------------------------------
const CLUSTER_WINDOW = (100.0, 800.0)       # ms
const N_PERMUTATIONS = 1000                 # as in the MATLAB scripts
const CLUSTER_SEED   = 2026

# --- Bins used (after the bin operations of Script #7c) --------------------------------
# Stimulus type: 55 = 2-UR, 56 = 3-UR, 52 = 3-RWP, 53 = 3-RSP
# Probe type:    46–49 = repeated/unrepeated × 3-RWP/3-RSP; 50 = repeated,
#                51 = unrepeated (3-RWP and 3-RSP collapsed); 54 = overall mean.
const STIMULUS_BINS = [55 => "2-UR", 56 => "3-UR", 52 => "3-RWP", 53 => "3-RSP"]
const PROBE_BINS = [46 => "Repeated_3-RWP", 47 => "Repeated_3-RSP", 48 => "Unrepeated_3-RWP",
                    49 => "Unrepeated_3-RSP", 50 => "Repeated", 51 => "Unrepeated"]
const CONDITIONS = ["2-UR", "3-UR", "3-RSP", "3-RWP"]   # order used in the figures
const COLORS = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"]

# --- Published values (Stage 2 manuscript), for the comparison table --------------------
const PAPER = Dict(
    "N2 amplitude, stimulus type (5.a)"   => (F=1.24, df1=3.0, df2=156.0, p=0.30, eta_p2=0.023),
    "N2 peak latency, stimulus type (5.b)" => (F=1.12, df1=3.0, df2=156.0, p=0.34, eta_p2=0.021),
    "FN400 amplitude, stimulus type (6.a)" => (F=0.18, df1=2.38, df2=123.76, p=0.87, eta_p2=0.003),
    "FN400 peak latency, stimulus type (6.b)" => (F=0.83, df1=3.0, df2=156.0, p=0.48, eta_p2=0.016),
    "N2 onset latency, probe type (4.b)"   => (F=0.78, df1=1.0, df2=39.0, p=0.38, eta_p2=0.020),
    "FN400 onset latency, probe type (4.c)" => (F=0.38, df1=1.0, df2=29.0, p=0.54, eta_p2=0.013),
    "N2 amplitude, probe type (8.a)"      => (F=0.60, df1=1.0, df2=52.0, p=0.44, eta_p2=0.011),
    "N2 peak latency, probe type (8.b)"   => (F=2.40, df1=1.0, df2=52.0, p=0.13, eta_p2=0.044),
    "FN400 amplitude, probe type (9.a)"   => (F=1.59, df1=1.0, df2=52.0, p=0.21, eta_p2=0.030),
    "FN400 peak latency, probe type (9.b)" => (F=0.22, df1=1.0, df2=52.0, p=0.64, eta_p2=0.004))
# Exploratory 2 × 2 (probe × proximity) on onset latency.
const PAPER_2X2 = Dict(
    ("N2", "probe") => (F=2.22, p=0.144), ("N2", "stimulus") => (F=0.10, p=0.757),
    ("N2", "probe:stimulus") => (F=1.66, p=0.205),
    ("FN400", "probe") => (F=0.10, p=0.754), ("FN400", "stimulus") => (F=4.06, p=0.053),
    ("FN400", "probe:stimulus") => (F=0.98, p=0.330))
# Cluster tests (manuscript and results_clustGND_*_probe.txt).
const PAPER_CLUSTERS = Dict("N2" => [(191.4, 449.2)], "FN400" => [(195.3, 464.8)], "N2pc" => [])
