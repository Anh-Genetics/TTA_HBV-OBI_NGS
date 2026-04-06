# HBV Reference Panel for Genotyping
# =====================================
# This directory should contain the file: hbv_genotype_panel.fasta
#
# The FASTA file should include representative complete or near-complete HBV
# genome sequences for genotypes A through I (and J), with clearly labelled
# FASTA headers indicating genotype/subgenotype.
#
# Recommended header format:
#   >AY_genotype_A2_J02203  Hepatitis B virus genotype A2, complete genome
#   >AB_genotype_B1_AB073858  Hepatitis B virus genotype B1
#   etc.
#
# Recommended sources:
#   - NCBI HBV reference sequences:
#     https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/Database/nph-select.cgi?taxid=10407
#   - HBVdb (https://hbvdb.lyon.inserm.fr/HBVdb/)
#   - Published curated panels (e.g. Sunbul 2014, Kramvis 2014)
#
# Minimum recommended panel (one representative per genotype):
#   Genotype A : J02203 (A2)
#   Genotype B : AB073858 (B1)
#   Genotype C : AB014381 (C1)
#   Genotype D : X65257 (D1)
#   Genotype E : X75657
#   Genotype F : X69798 (F1)
#   Genotype G : AF160501
#   Genotype H : AY090454
#   Genotype I : FJ023669
#
# Extended panels should include multiple representatives per genotype
# for accurate subgenotype assignment.
#
# HOW TO DOWNLOAD (NCBI Entrez):
# ---------------------------------
# pip install biopython
# python -c "
# from Bio import Entrez, SeqIO
# Entrez.email = 'your@email.com'
# accessions = ['J02203','AB073858','AB014381','X65257','X75657','X69798','AF160501','AY090454','FJ023669']
# for acc in accessions:
#     handle = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
#     print(handle.read())
# " > hbv_genotype_panel.fasta
#
# NOTE: The pipeline expects this file at:
#   resources/hbv_references/hbv_genotype_panel.fasta
# Override with --hbv_ref /path/to/your.fasta

# A minimal synthetic placeholder FASTA is NOT included here to avoid
# distributing sequences of unknown licensing status.
# Please download from NCBI as instructed above.
