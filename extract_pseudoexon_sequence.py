#!/usr/bin/env python3
"""
extract_pseudoexon_sequence.py

Usage:
    python extract_pseudoexon_sequence.py chr12:32711230 --ag -50 --dg 42

Requirements:
    pip install pysam biopython
    (The fasta must be indexed with samtools faidx; pysam will use the .fai)

What it does:
    - Computes retained range: start = pos + AG, end = pos + DG
    - Fetches sequence (1-based inclusive coords)
    - Handles minus strand (when end < start) by reverse complementing
    - Translates sequence and prints summaries
"""

import argparse
import sys
from math import floor
try:
    import pysam
except Exception as e:
    sys.exit("pysam is required: pip install pysam\nError: " + str(e))

try:
    from Bio.Seq import Seq
except Exception as e:
    sys.exit("biopython is required: pip install biopython\nError: " + str(e))

# Default path
HOMO_SAPIENS_FASTA = '/stornext/Bioinf/data/lab_bahlo/ref_db/human/hg38/GATK/fasta/Homo_sapiens_assembly38.fasta'

def parse_coord(coord):
    # Accept formats like chr1:12345 or 1:12345
    if ":" not in coord:
        raise argparse.ArgumentTypeError("Coordinate must be like chr12:32711230")
    chrom, pos = coord.split(":", 1)
    try:
        pos = int(pos)
    except ValueError:
        raise argparse.ArgumentTypeError("Position must be an integer")
    return chrom, pos


def maybe_try_alt_contig(fasta, chrom, start0, end0):
    # Helper: try with given contig, else if contig missing try dropping/adding 'chr'
    try:
        seq = fasta.fetch(chrom, start0, end0)
        return seq, chrom
    except Exception:
        # try toggling 'chr' prefix
        if chrom.startswith("chr"):
            alt = chrom[3:]
        else:
            alt = "chr" + chrom
        try:
            seq = fasta.fetch(alt, start0, end0)
            return seq, alt
        except Exception as e:
            raise RuntimeError(f"Failed to fetch region {chrom}:{start0+1}-{end0} "
                               f"(also tried {alt}). Error: {e}")


def summarize_and_print(seq, chrom, start1, end1, is_minus_strand=False):
    """
    seq: uppercase nucleotide sequence (string)
    start1,end1: 1-based inclusive coordinates of seq on chrom
    is_minus_strand: whether this is from minus strand (affects display)
    """
    seq = seq.upper().replace("\n", "").replace(" ", "")
    nuc_len = len(seq)
    strand_info = " (minus strand)" if is_minus_strand else " (plus strand)"
    print(f"> Region: {chrom}:{start1}-{end1}{strand_info}  (length = {nuc_len} bp)")
    print(f"Sequence: {seq}")
    print("")

    frame = 0
    if nuc_len <= frame:
        aa_seq = ""
        full_codons = 0
    else:
        # Only translate full codons (Biopython will translate up to last full codon)
        sub = seq[frame:]
        aa_seq = str(Seq(sub).translate())
        full_codons = len(sub) // 3

    aa_len = len(aa_seq)
    contains_stop = "*" in aa_seq
    last_codon_complete = ((nuc_len - frame) % 3) == 0
    inframe_summary = "in-frame (no partial codon at end)" if last_codon_complete else "out-of-frame (partial codon at end)"
    print(f"Amino acid length  : {aa_len} aa (full codons: {full_codons})")
    print(f"In-frame?          : {inframe_summary}")
    print(f"Contains stop (*)  : {contains_stop}")
    # print AA sequence
    print(f"AA (first 60 aa)   : {aa_seq[:60]}")
    print("")


def main():
    ap = argparse.ArgumentParser(description="Extract genomic range from FASTA using AG/DG offsets and translate")
    ap.add_argument("coord", help="Genomic coord like chr12:32711230")
    ap.add_argument("--ag", type=int, required=True, help="AG offset (can be negative). start = pos + AG")
    ap.add_argument("--dg", type=int, required=True, help="DG offset (can be negative). end = pos + DG")
    args = ap.parse_args()

    chrom, pos = parse_coord(args.coord)
    start = pos + args.ag
    end = pos + args.dg

    # Check if this is minus strand (end < start)
    is_minus_strand = end < start
    
    if is_minus_strand:
        # For minus strand, swap coordinates for fetching
        fetch_start = end
        fetch_end = start
    else:
        fetch_start = start
        fetch_end = end

    if fetch_start < 1:
        sys.exit(f"Computed coordinate < 1: {fetch_start}. Check AG/DG values and input coordinate.")

    # pysam.fetch expects 0-based start, end-exclusive
    start0 = fetch_start - 1
    end0 = fetch_end  # end is inclusive in our coords, so end0 = end (end-exclusive)

    # open fasta
    try:
        fasta = pysam.FastaFile(HOMO_SAPIENS_FASTA)
    except Exception as e:
        sys.exit(f"Failed to open fasta: {e}")

    try:
        seq, used_contig = maybe_try_alt_contig(fasta, chrom, start0, end0)
    except RuntimeError as e:
        sys.exit(str(e))

    # If minus strand, reverse complement the sequence
    if is_minus_strand:
        seq = str(Seq(seq).reverse_complement())

    # print computed range explicitly
    print(f"Input: {args.coord}  AG={args.ag}  DG={args.dg}")
    if is_minus_strand:
        print(f"Detected minus strand (end < start)")
        print(f"Fetched region: {used_contig}:{fetch_start}-{fetch_end}, then reverse complemented")
        print(f"Reported as: {used_contig}:{start}-{end} (n={fetch_end-fetch_start+1} bp)\n")
    else:
        print(f"Computed retained intron range: {used_contig}:{start}-{end} (n={end-start+1} bp)\n")

    # Summarize & translate
    summarize_and_print(seq, used_contig, start, end, is_minus_strand)

    fasta.close()


if __name__ == "__main__":
    main()
