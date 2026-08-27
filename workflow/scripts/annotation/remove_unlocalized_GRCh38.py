#!/usr/bin/env python3
"""Keep only the chromosome-level records of a GRCh38 FASTA.

The GDC GRCh38.d1.vd1 FASTA tags every header with an `rl:` (region length
class) attribute; `rl:Chromosome` marks the primary assembly, and everything
else (unlocalized scaffolds, alts, decoys, viral sequence) is dropped. Chain
files built against those extra contigs produce spurious liftover targets.

References without `rl:` tags — plain Ensembl/UCSC FASTAs, or a chromosome
subset prepared for a test run — carry no such marker. Rather than emit an
empty FASTA, the whole input is passed through and a note is written to
stderr.
"""

import sys


def _is_chromosome(header):
    """Return True/False from the header's `rl:` tag, or None when absent."""
    for tag in header.split():
        if tag.startswith("rl:"):
            return tag[3:] == "Chromosome"
    return None


def main():
    path = sys.argv[1]

    tagged = False
    with open(path) as f:
        for line in f:
            if line.startswith(">") and _is_chromosome(line.rstrip("\n")) is not None:
                tagged = True
                break

    if not tagged:
        sys.stderr.write(
            "remove_unlocalized_GRCh38.py: no 'rl:' tags in {}; "
            "passing every record through\n".format(path))

    keep = True
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                keep = True if not tagged else bool(_is_chromosome(line.rstrip("\n")))
            if keep:
                sys.stdout.write(line)


if __name__ == "__main__":
    main()
