#!/usr/bin/env python3
import argparse
import random
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fragment_sequence(seq, seed, target_fragments, min_length=1):
    """
    Fragment a sequence into target_fragments pieces at random cut points.

    Args:
        seq:              The nucleotide sequence string.
        seed:             Random seed for reproducibility.
        target_fragments: Desired number of fragments.
        min_length:       Minimum contig length to retain.

    Returns:
        List of contig strings passing the min_length filter.
    """
    length = len(seq)

    # Guard: nothing to cut
    if length == 0:
        return []

    # Guard: can't cut more times than there are internal positions
    max_possible = max(1, length)
    if target_fragments > max_possible:
        print(
            f"WARNING: Sequence length ({length} bp) allows at most "
            f"{max_possible} fragment(s); requested {target_fragments}. "
            f"Adjusting.",
            file=sys.stderr,
        )
        target_fragments = max_possible

    # Single fragment — no cutting needed
    if target_fragments <= 1:
        return [seq] if length >= min_length else []

    random.seed(seed)

    # Sample (target_fragments - 1) unique cut positions from internal sites
    n_cuts = min(target_fragments - 1, length - 1)
    cut_positions = sorted(random.sample(range(1, length), n_cuts))

    contigs = []
    prev = 0
    for pos in cut_positions:
        fragment = seq[prev:pos]
        if len(fragment) >= min_length:
            contigs.append(fragment)
        prev = pos

    # Final fragment
    fragment = seq[prev:]
    if len(fragment) >= min_length:
        contigs.append(fragment)

    return contigs


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Fragment genome sequences into contigs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input",  required=True,
                        help="Input FASTA genome file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output contigs FASTA file")
    parser.add_argument("--mapping", default=None,
                        help="Path for contig→genome mapping TSV "
                             "(default: contig2bin.tsv next to output)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducible fragmentation")
    parser.add_argument("--base_fragments", type=int, default=10,
                        help="Base number of fragments per sequence")
    parser.add_argument("--fragmentation", type=float, default=1.0,
                        help="Multiplier applied to base_fragments "
                             "(>1 = more fragments, <1 = fewer)")
    parser.add_argument("--min_length", type=int, default=1,
                        help="Minimum contig length (bp) to retain")
    return parser


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    # Validate inputs
    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: Input file not found: {args.input}")

    if args.fragmentation <= 0:
        sys.exit("ERROR: --fragmentation must be > 0")

    # Resolve mapping file path
    if args.mapping:
        map_file = args.mapping
    else:
        out_dir  = os.path.dirname(os.path.abspath(args.output))
        map_file = os.path.join(out_dir, "contig2bin.tsv")

    # Compute target fragments once (seed no longer affects count)
    target_fragments = max(1, int(args.base_fragments * args.fragmentation))

    contig_index  = 1
    total_written = 0
    skipped       = 0

    with open(args.output, "w") as fasta_out, \
         open(map_file,    "w") as map_out:

        map_out.write("contig_id\toriginal_sequence\n")

        for seq_idx, record in enumerate(SeqIO.parse(args.input, "fasta")):
            seq = str(record.seq)

            if not seq:
                print(
                    f"WARNING: Skipping empty record '{record.id}'",
                    file=sys.stderr,
                )
                skipped += 1
                continue

            # Unique seed per record preserves reproducibility
            # while avoiding identical cut patterns across records
            record_seed = args.seed + seq_idx

            contigs = fragment_sequence(
                seq,
                seed=record_seed,
                target_fragments=target_fragments,
                min_length=args.min_length,
            )

            for contig_seq in contigs:
                contig_name = f"contig_{contig_index}"

                contig_record = SeqRecord(
                    Seq(contig_seq),
                    id=contig_name,
                    description=f"source={record.id} len={len(contig_seq)}",
                )
                SeqIO.write(contig_record, fasta_out, "fasta")
                map_out.write(f"{contig_name}\t{record.id}\n")

                contig_index  += 1
                total_written += 1

    print(f"Generated {total_written} contigs → {args.output}")
    if skipped:
        print(f"Skipped    {skipped} empty records", file=sys.stderr)
    print(f"Mapping file         → {map_file}")


if __name__ == "__main__":
    main()