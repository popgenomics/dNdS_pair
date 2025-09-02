# dNdS\_pair

Pairwise codon-based divergence (**dN**, **dS**, **ω = dN/dS**) between two sequences drawn from a single-gene codon alignment. Includes a companion script to run **PAML/yn00** on the same pair for cross-checks.

---

## Requirements

* Python ≥ 3.8
* [Biopython](https://biopython.org/) (`Bio.SeqIO`, `Bio.Seq`)
* POSIX shell utils for the PAML script (`bash`, `awk`, `python3`)
* Optional for comparison: **PAML** ≥ 4.10 (`yn00` available in `PATH`)

---

## Usage

### Python implementation

```bash
python3 scripts/dNdS_pair.py --fasta test/test.fas --nameA seqA --nameB seqB
```

Equivalent CLI forms:

```bash
python3 scripts/dNdS_pair.py test/test.fas seqA seqB
python3 scripts/dNdS_pair.py fasta=test/test.fas nameA=seqA nameB=seqB
```

**Input assumptions**

* `--fasta` is a **codon alignment** (single gene, multiple taxa) in FASTA.
* Sequence IDs are taken as the **first whitespace-delimited token** after `>`.
* The two sequences are truncated to a common length divisible by 3.
* Codons with `N`, `-`, or translating to a **stop** are discarded.
* If **all** shortest one-step mutation paths between two differing codons pass through a stop, that site is discarded.

**Output**

Tab-separated header + one data row. Columns:

| Column            | Meaning                                                                                                                    |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `gene`            | Input FASTA file name.                                                                                                     |
| `nameA`           | ID of the first sequence.                                                                                                  |
| `nameB`           | ID of the second sequence.                                                                                                 |
| `codons_total`    | `floor(min(len(A), len(B)) / 3)` before filtering.                                                                         |
| `codons_retained` | Codon sites analyzed after all filters.                                                                                    |
| `Ls`              | Sum of **synonymous site opportunities** across retained codons (mean `nS` of the two codons, summed).                     |
| `Ln`              | Sum of **nonsynonymous site opportunities** (`nN`) across retained codons.                                                 |
| `S_changes`       | Expected count of synonymous changes, averaged over all shortest one-nucleotide paths that avoid stops, summed over sites. |
| `N_changes`       | Expected count of nonsynonymous changes, same rule.                                                                        |
| `pS`              | Uncorrected proportion of synonymous differences: `S_changes / Ls`.                                                        |
| `pN`              | Uncorrected proportion of nonsynonymous differences: `N_changes / Ln`.                                                     |
| `dS`              | JC69-corrected synonymous divergence per site: `-3/4 * ln(1 - 4/3 * pS)` (`NA` if `pS ≥ 0.75` or undefined).               |
| `dN`              | JC69-corrected nonsynonymous divergence per site from `pN`.                                                                |
| `omega`           | `dN / dS` (`NA` if `dS` is 0 or undefined).                                                                                |

**Method (concise)**

1. Build a sense-codon table with average counts of one-step **synonymous** (`nS`) and **nonsynonymous** (`nN`) neighbors (stops excluded).
2. For each retained site, opportunities `Ls`/`Ln` = mean `nS`/`nN` across the two codons; sum over sites.
3. For codon pairs differing at `d` positions, enumerate the `d!` shortest paths (one base per step). At each step, AA unchanged ⇒ S; changed ⇒ N. Paths that hit a stop are discarded. Average S/N over valid paths; sum over sites.
4. Convert `pS`, `pN` to `dS`, `dN` with JC69; report `ω`.

**Notes**

* `Ls + Ln = 3 × codons_retained`.
* High `pS` (\~0.6–0.75) indicates saturation; JC69 becomes unstable and `dS` is set `NA`.
* Align AA first and back-translate to codons (e.g., PAL2NAL, MACSE) to avoid frameshifts/stops.

---

## PAML comparison (`run_yn00_pair.sh`)

Purpose: replicate pairwise Ka/Ks using **PAML/yn00** on the same pair to verify magnitudes of `dS`, `dN`, and `ω`.

```bash
bash scripts/run_yn00_pair.sh  test/test.fas  seqA  seqB
```

What it does:

* Extracts `seqA` and `seqB` from FASTA.
* Converts to PHYLIP (names ≤10 chars, **two spaces** before sequence), sanitizes to `ACGTN-`, trims to multiple of 3.
* Writes a minimal `yn00.ctl` and runs `yn00`.
* Produces `yn00.out` in the current directory.

Dependencies:

* **PAML** ≥ 4.10 (`yn00` in `PATH`).
* `bash`, `awk`, `python3`.

Reading the PAML output:

* Prefer **Nei–Gojobori (NG86)** or **LWL85** for comparison with the JC-style output of this script.
* `ls` reported by `yn00` should match `codons_retained`.
* In the LWL section, `L(i) sum` should equal `Ls + Ln`.

---

## Reproducibility example

```bash
# Python
python3 scripts/dNdS_pair.py --fasta test/test.fas --nameA seqA --nameB seqB \
  > results_python.tsv

# PAML sanity check
bash scripts/run_yn00_pair.sh test/test.fas seqA seqB
# Inspect yn00.out (NG86 or LWL85 sections)
```

---

## Limitations

* JC69 assumes equal base frequencies and equal rates among nucleotides; no among-site rate heterogeneity.
* Corrections are applied separately to S and N subsets at the nucleotide level.
* Sites for which **all** shortest paths pass through a stop are discarded.
* Strong saturation, misalignments, or hidden frameshifts will bias estimates.

---

## Reference

* Jukes, T. H., & Cantor, C. R. (1969). *Evolution of Protein Molecules*. In **Mammalian Protein Metabolism**, Vol. 3, 21–132. Academic Press.

---

