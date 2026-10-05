#!/usr/bin/env python3

import argparse
import gzip
import ahocorasick

# Open gzip or regular files transparently
def opengz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

# Build Aho-Corasick automaton for fast multi-pattern search
def build_automaton(sites):
    A = ahocorasick.Automaton()
    for s in sites:
        A.add_word(s, s)
    A.make_automaton()
    return A

# Find earliest cutsite position in a read
def cut_position(seq, A):
    pos = None
    for end, site in A.iter(seq):
        start = end - len(site) + 1
        if pos is None or start < pos:
            pos = start
    return pos

# Expand cutsite string (e.g., "A^GCTT") into all IUPAC possibilities
IUPAC = {
    "A": ["A"],
    "T": ["T"],
    "C": ["C"],
    "G": ["G"],
    "N": ["A", "C", "G", "T"],
}

def expand_cut_sites(site_string):
    import itertools
    expanded = []
    for site in site_string.split(","):
        site = site.replace("^", "").upper()
        letters = [IUPAC.get(base, [base]) for base in site]
        for combo in itertools.product(*letters):
            expanded.append("".join(combo))
    return sorted(set(expanded))

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fastq", required=True)
    p.add_argument("--cutsite", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--rmuntrim", action="store_true")
    args = p.parse_args()

    # Expand cutsites and build automaton
    sites = expand_cut_sites(args.cutsite)
    A = build_automaton(sites)

    total = trimmed = dummy = 0

    with opengz(args.fastq) as fin, open(args.out, "wt") as fout:
        while True:
            ID = fin.readline()
            if not ID:
                break
            seq = fin.readline().rstrip()
            fin.readline()  # skip "+"
            qual = fin.readline().rstrip()

            total += 1

            cut = cut_position(seq, A)

            if cut is not None:
                trimmed += 1
                if cut == 0: # If cutsite is at position 0, replace with length-1 N read
                    dummy += 1
                    fout.write(f"{ID}N\n+\n!\n")
                else:
                    fout.write(f"{ID}{seq[:cut]}\n+\n{qual[:cut]}\n")
            elif not args.rmuntrim:
                fout.write(f"{ID}{seq}\n+\n{qual}\n")

    #print(f"Total reads: {total}")
    #print(f"Trimmed reads: {trimmed}")
    #print(f"Dummy reads (site at pos0): {dummy}")

if __name__ == "__main__":
    main()
