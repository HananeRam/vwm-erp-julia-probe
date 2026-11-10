# Analysis settings for the probe-locked ERPs, aligned with the MATLAB/EEGLAB/ERPLAB
# scripts (#1–#5, #6b/#6c, #7b/#7c, #11c, ERP_Cluster_Analyses Steps 2–4).
# Edit here, not in the pipeline code.

# --- Data ------------------------------------------------------------------
# Set the DS006979 environment variable, or put the dataset in data/ds006979.
const BIDS_ROOT = get(ENV, "DS006979", joinpath(@__DIR__, "..", "data", "ds006979"))
const TASK = "visualworkingmemorytask"
# Output folder; RESULTS_DIR=/path runs the pipeline into another folder (e.g. a sensitivity run).
const RESULTS = get(ENV, "RESULTS_DIR", joinpath(@__DIR__, "..", "results"))
const REFERENCE = joinpath(@__DIR__, "..", "reference")      # MATLAB outputs, for validation

# --- Participants ----------------------------------------------------------
# The MATLAB probe-locked analyses used 53 participants (SUB list in Scripts #6b–#11c).
const ALL_IDS = [1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                 16, 17, 18, 19, 21, 22, 23, 24, 25, 28, 30,
                 32, 33, 34, 35, 36, 39, 40, 41, 43, 46, 47, 48,
                 49, 50, 51, 52, 53, 54, 58, 59, 61, 62, 64,
                 69, 70, 72, 73, 74, 76]
# Participant 3: in its original Script #0 dataset (3.set) the three unused external
# channels (EXG1, EXG7, EXG8) were not removed, so every later (position-based) MATLAB step
# used channels shifted by one from position 65 on. 3.set was corrected and sub-03 was
# reprocessed from Script #1 (new ICA, components 1 and 3 removed); it is included (N = 53).
# Julia selects channels by label, so it needs only the new ICA weights file for sub-03.
const EXCLUDED_IDS = Int[]
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
# Optional per-participant override of ICs_to_Remove.csv (id => components). Empty: the
# components in ICs_to_Remove.xlsx are removed for every participant, as in the revised
# MATLAB analysis (Script #4 was rerun for 46, 49 and 69 with the table's components).
const ICS_APPLIED = Dict{Int,Vector{Int}}()
const ICS_MODE = get(ENV, "ICS", "applied")

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

# --- Revised MATLAB/R results (N = 53, corrected sub-03, ICs as in ICs_to_Remove.xlsx), for the comparison table --------------------
const PAPER = Dict(
    "N2 amplitude, stimulus type (5.a)"   => (F=1.12, df1=3.0, df2=156.0, p=0.34, eta_p2=0.021),
    "N2 peak latency, stimulus type (5.b)" => (F=0.76, df1=3.0, df2=156.0, p=0.52, eta_p2=0.014),
    "FN400 amplitude, stimulus type (6.a)" => (F=0.19, df1=2.38, df2=123.67, p=0.87, eta_p2=0.004),
    "FN400 peak latency, stimulus type (6.b)" => (F=1.09, df1=3.0, df2=156.0, p=0.36, eta_p2=0.020),
    "N2 onset latency, probe type (4.b)"   => (F=0.07, df1=1.0, df2=39.0, p=0.80, eta_p2=0.002),
    "FN400 onset latency, probe type (4.c)" => (F=0.94, df1=1.0, df2=28.0, p=0.34, eta_p2=0.033),
    "N2 amplitude, probe type (8.a)"      => (F=0.32, df1=1.0, df2=52.0, p=0.57, eta_p2=0.006),
    "N2 peak latency, probe type (8.b)"   => (F=3.63, df1=1.0, df2=52.0, p=0.06, eta_p2=0.065),
    "FN400 amplitude, probe type (9.a)"   => (F=1.61, df1=1.0, df2=52.0, p=0.21, eta_p2=0.030),
    "FN400 peak latency, probe type (9.b)" => (F=0.26, df1=1.0, df2=52.0, p=0.61, eta_p2=0.005))
# Exploratory 2 × 2 (probe × proximity) on onset latency.
const PAPER_2X2 = Dict(
    ("N2", "probe") => (F=3.10, p=0.086), ("N2", "stimulus") => (F=0.01, p=0.917),
    ("N2", "probe:stimulus") => (F=1.35, p=0.252),
    ("FN400", "probe") => (F=0.45, p=0.510), ("FN400", "stimulus") => (F=2.72, p=0.110),
    ("FN400", "probe:stimulus") => (F=0.21, p=0.648))
# Cluster tests (manuscript and results_clustGND_*_probe.txt).
const PAPER_CLUSTERS = Dict("N2" => [(191.4, 449.2)], "FN400" => [(195.3, 464.8)], "N2pc" => [])
