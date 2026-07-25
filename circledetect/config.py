import os

# External Tools (Ensure these are installed and in PATH)
BWA_MEM2_PATH = "bwa-mem2"
SAMBLASTER_PATH = "samblaster"
SAMTOOLS_PATH = "samtools"

# Algorithm Parameters
MIN_MAPQ = 20
MIN_CLIP_LEN = 15       # Minimum soft-clip length to consider
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
# Uniform background nucleotide frequency used as the null model.
BG_FREQ = 0.25
# Realignment context (one-sided); the doubled junction window is 2 * 2 * this.
REALIGN_FLANK = 50

# Output
OUTPUT_DIR = "./eccdna_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)