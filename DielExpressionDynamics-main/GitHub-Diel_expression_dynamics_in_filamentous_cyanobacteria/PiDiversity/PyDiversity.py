import subprocess
import sys
import os

# Function to create a conda environment and install necessary packages
def setup_conda_environment():
    env_name = "genomic_analysis"
    python_version = "3.8"
    
    # Check if the conda environment exists
    try:
        envs_output = subprocess.check_output(["conda", "env", "list"]).decode("utf-8")
        if env_name in envs_output:
            print(f"Conda environment '{env_name}' already exists.")
        else:
            raise Exception(f"Environment '{env_name}' not found.")
    except Exception as e:
        print(e)
        print(f"Creating conda environment '{env_name}'...")
        subprocess.run(["conda", "create", "-n", env_name, f"python={python_version}", "-y"], check=True)
        print(f"Installing necessary packages in '{env_name}'...")
        subprocess.run(["conda", "install", "-n", env_name, "-c", "bioconda", "biopython", "pysam", "pyvcf", "-y"], check=True)
    
    print(f"Conda environment '{env_name}' is ready. Please activate it using the command below and then rerun the script:")
    print(f"\n    conda activate {env_name}\n")
    sys.exit()

# Check if the script is being run inside the right environment
current_env = os.environ.get('CONDA_DEFAULT_ENV')
if current_env != "genomic_analysis":
    print(f"Current environment: {current_env}")
    setup_conda_environment()

# After environment setup, import necessary libraries
from itertools import product
import vcf  # pyvcf for VCF file handling
from Bio.Blast import NCBIXML  # Part of Biopython for BLAST XML parsing
from collections import defaultdict
from Bio import SeqIO  # Biopython for sequence input/output
from Bio.Seq import Seq  # Biopython for sequence objects
from Bio.SeqRecord import SeqRecord  # Biopython for sequence records
import csv
import pysam  # For working with BAM files

# Define input files and parameters
bam_file = "csome_s.bam"
fasta_file = "csome.fasta"
chromosome = "NC_010628.1"

# Define the window of interest within the chromosome
window_start = 6058888
window_end = 6059027

# Define output file
output_file = "Npun_locus_VP_PiDiversity.txt"

# Open the BAM file using pysam
bamfile = pysam.AlignmentFile(bam_file, "rb")

# Parse the FASTA file and process the sequence of the specified chromosome
for myseq in SeqIO.parse(fasta_file, "fasta"):
    print(myseq.id)
    print(len(myseq.seq))

    if myseq.id == chromosome:
        # Open the output file and write the header
        with open(output_file, 'w') as f:
            f.write("Position\tBase\tDiversity\tCoverage\n")

        # Loop through the references in the BAM file
        for i in bamfile.references:
            if i == chromosome:
                print(i)

                # Loop through the specified window and analyze each position
                for pileupcolumn in bamfile.pileup(i, window_start, window_end):
                    column = str()
                    variants = 0
                    readcount = 0
                    Gs = 0
                    As = 0
                    Cs = 0
                    Ts = 0

                    truepos = pileupcolumn.pos
                    truebase = myseq.seq[truepos]

                    if window_start < pileupcolumn.pos < window_end:

                        for pileupread in pileupcolumn.pileups:
                            readcount += 1

                            if pileupread.query_position is not None:
                                base = pileupread.alignment.query_sequence[pileupread.query_position]
                                if base == "A":
                                    As += 1
                                elif base == "G":
                                    Gs += 1
                                elif base == "C":
                                    Cs += 1
                                elif base == "T":
                                    Ts += 1
                                column += base

                                if truebase != base:
                                    variants += 1

                        if len(column) > 0:
                            Gfreq = Gs / len(column) if Gs > 0 else 0
                            Afreq = As / len(column) if As > 0 else 0
                            Cfreq = Cs / len(column) if Cs > 0 else 0
                            Tfreq = Ts / len(column) if Ts > 0 else 0

                        # Calculate nucleotide diversity (Pi)
                        pi = 1 - (Afreq ** 2 + Tfreq ** 2 + Gfreq ** 2 + Cfreq ** 2)
                        freqdict = {"Af": Afreq, "Tf": Tfreq, "Gf": Gfreq, "Cf": Cfreq}
                        majorBasefreq = max(Afreq, Tfreq, Gfreq, Cfreq)
                        minorsum = 1 - majorBasefreq
                        majorBasevar = max(freqdict, key=freqdict.get)

                        columnlen = len(column)
                        print("pos:", pileupcolumn.pos)
                        print("Pi: ", pi)
                        mybase = str(myseq.seq[pileupcolumn.pos])

                        print("\n")

                        # Append results to the output file
                        with open(output_file, 'a') as f:
                            f.write(f"{truepos}\t{mybase}\t{pi}\t{columnlen}\n")
