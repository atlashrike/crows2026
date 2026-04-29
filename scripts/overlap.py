#!/usr/bin/env python3
import argparse

def read_pos(path: str) -> set[str]:
    s = set()
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # take first whitespace-delimited token (matches your R read.table default)
            tok = line.split()[0]
            if tok:
                s.add(tok)
    return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gatk", required=True)
    ap.add_argument("--samtools", required=True)
    ap.add_argument("--angsd", required=True)
    ap.add_argument("--out", default="overlapping_variants_pos.txt")
    args = ap.parse_args()

    gatk = read_pos(args.gatk)
    sam  = read_pos(args.samtools)
    ang  = read_pos(args.angsd)

    inter = sorted(gatk & sam & ang)

    with open(args.out, "w") as out:
        out.write("LocusName\n")
        for x in inter:
            out.write(f"{x}\n")

    print(f"GATK={len(gatk)} SAMtools={len(sam)} ANGSD={len(ang)} INTERSECTION={len(inter)}")

if __name__ == "__main__":
    main()
