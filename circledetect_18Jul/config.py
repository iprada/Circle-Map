import os

# External Tools (Ensure these are installed and in PATH)
BWA_MEM2_PATH = "bwa-mem2"
SAMBLASTER_PATH = "samblaster"
SAMTOOLS_PATH = "samtools"

# Algorithm Parameters
MIN_MAPQ = 10           # Min MAPQ for split reads (Circle-Map default: 10)
MIN_MAPQ_DISCORDANT = 20  # Min MAPQ for discordant reads (higher = fewer false SV calls)
MIN_CLIP_LEN = 15       # Minimum soft-clip length to consider
MIN_SPLIT_SUPPORT = 3   # Min split reads per junction (filters 72% of 1-read noise)
PROB_CUTOFF = 0.95      # Posterior probability threshold
MAX_INSERT_SIZE = 1000  # Max expected insert size for proper pairs

# Markov Chain / Alignment Parameters
MATCH_SCORE = 2.0
MISMATCH_PENALTY = -3.0
GAP_OPEN = -5.0         # Markov transition cost: Match -> Gap
GAP_EXTEND = -1.0       # Markov transition cost: Gap -> Gap

# Base encoding for the vectorized realignment kernel.
# A/C/G/T -> 0/1/2/3 ; anything else (incl. N) -> 4.
BASE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
# Background nucleotide frequency is computed per-junction from the local
# reference window (see ``compute_bg_freqs`` in realigner.py), following
# Circle-Map's position-specific scoring matrix approach.
REALIGN_FLANK = 50

# --- Insert size distribution estimation (mirrors Circle-Map -s / -iq) ---
INSERT_SAMPLE_SIZE = 100_000   # Number of concordant pairs to sample
INSERT_MAPQ = 60               # MAPQ cutoff for insert-size estimation

# --- Confidence scoring (depth ratio + continuity) ---
# Depth-ratio windows (mirror Circle-Map -E / -b).
DEPTH_INSIDE_LENGTH = 100      # bp inside each breakpoint used as "inside edge"
DEPTH_EXTENSION = 200          # bp flanking each breakpoint used as "outside flank"
DEPTH_MIN_BASEQ = 20
COVERAGE_MAPQ = 20             # MAPQ filter applied to coverage via read_callback
# Continuity scoring.
CONTINUITY_MAPQ = 20
JUNCTION_WINDOW = 200          # Half-window around each breakpoint for bridging
# Discordant reads within this distance of a junction snap to that junction's
# node in the breakpoint graph (mirrors Circle-Map -K clustering_dist).
DISCORDANT_MATCH_WINDOW = 500
# Minimum discordant read pairs required to create a directed edge between two nodes.
# Direction is preserved so SCC detection finds true multi-segment circles
# (which require bidirectional connectivity). 1 allows single-pair edges.
MIN_DISCORDANT_EDGES = 1

# Output
OUTPUT_DIR = "./eccdna_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)