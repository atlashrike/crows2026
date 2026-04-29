# modified from https://github.com/mmosmond/spacetrees-ms/blob/main/athaliana/longreads/Snakefile rule polarize_hap

import sys

est_files = snakemake.input.est_sfs
pval_files = snakemake.input.est_pvals
hap_file = snakemake.input.hap
out_file = snakemake.output.polar

nucs = ['A', 'C', 'G', 'T']
nflips = 0
ntot = 0

with open(hap_file, 'r') as hap_f, open(out_file, 'w') as out:
    for chunk_index in range(4):  # 4 parts now
        sfs_path = est_files[chunk_index]
        pval_path = pval_files[chunk_index]

        with open(sfs_path, 'r') as est_f, open(pval_path, 'r') as pv_f:
            # skip 8 header lines in pvalues file
            for _ in range(8):
                next(pv_f)

            for est_line, p_line in zip(est_f, pv_f):
                hap_line = next(hap_f).rstrip('\n')

                n = [int(round(float(x))) for x in est_line.split()[0].split(',')]
                major_idx = max(range(4), key=lambda i: n[i])
                major_allele = nucs[major_idx]
                minor_idx = max((i for i in range(4) if i != major_idx),
                                key=lambda i: n[i])
                minor_allele = nucs[minor_idx]

                fields = p_line.split()
                p_major_ancestral = float(fields[2])

                ancestral = major_allele
                derived = minor_allele
                if p_major_ancestral < 0.5:
                    ancestral = minor_allele
                    derived = major_allele

                haps = hap_line.split()
                reference = haps[3]
                alternate = haps[4]

                if reference != ancestral:
                    nflips += 1
                    haps[3] = ancestral
                    haps[4] = derived
                    for i in range(5, len(haps)):
                        old_val = int(haps[i])
                        haps[i] = str(1 - old_val)

                out.write(' '.join(haps) + '\n')
                ntot += 1

print(f"Flipped {nflips} out of {ntot} alleles, fraction: {nflips/ntot if ntot else 0:.6f}")