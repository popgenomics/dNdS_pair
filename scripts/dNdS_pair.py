#!/usr/bin/env python3
# Divergence dS/dN entre deux séquences codantes (Nei–Gojobori 1986, JC)
# Input: FASTA aligné d'un seul gène (X espèces)
# Usage: python dNdS_pair.py alignment.fasta nameA nameB

from Bio.SeqIO import parse
from Bio.Seq import Seq
import itertools, math
from os.path import basename

import argparse, textwrap, sys

def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # Normalise key=value or --key value
    norm = []
    for t in argv:
        if "=" in t and not t.startswith("-"):
            k, v = t.split("=", 1)
            norm.extend([f"--{k}", v])
        else:
            norm.append(t)

    desc = """\
    Pairwise codon-based divergence (dN, dS, ω=dN/dS) between two sequences from a single-gene codon alignment (FASTA).

    Method
      • Selects sequences with IDs exactly matching nameA and nameB.
      • Truncates both to the common length divisible by 3.
      • Precomputes, for each sense codon, the average number of synonymous (nS) and nonsynonymous (nN) one-step opportunities by enumerating all single-nucleotide neighbors; stop codons are excluded.
      • For each codon position:
          – Reject if either codon contains N or '-', or translates to a stop, or is not a valid sense codon.
          – Add opportunities Ls and Ln as the mean of nS/nN of the two codons.
          – Count expected synonymous (S_changes) and nonsynonymous (N_changes) changes by averaging over all shortest single-nucleotide mutation paths between the two codons; any path that crosses a stop is discarded.
      • Compute pS = S_changes / Ls and pN = N_changes / Ln.
      • Apply Jukes–Cantor (JC69) to obtain dS = -3/4·ln(1 - 4/3·pS) and dN analogously.
      • Report ω = dN/dS (NaN if dS is 0 or undefined).
    """
    epi = """\
    Output (tab-separated, one row):
      gene               Input FASTA file name.
      nameA              ID of the first sequence (exact match to FASTA header first token).
      nameB              ID of the second sequence (exact match to FASTA header first token).
      codons_total       Total codons before filtering = floor(min(len(A),len(B))/3).
      codons_retained    Codon sites retained after filters (no N, no '-', no stop, valid minimal paths).
      Ls                 Sum of synonymous site opportunities across retained codons (mean nS of the two codons, summed).
      Ln                 Sum of nonsynonymous site opportunities across retained codons (mean nN of the two codons, summed).
      S_changes          Expected count of synonymous changes across retained codons (averaged over minimal paths).
      N_changes          Expected count of nonsynonymous changes across retained codons (averaged over minimal paths).
      pS                 Proportion of synonymous differences = S_changes / Ls (uncorrected).
      pN                 Proportion of nonsynonymous differences = N_changes / Ln (uncorrected).
      dS                 JC-corrected synonymous divergence per site; NaN if pS ≥ 0.75 or invalid.
      dN                 JC-corrected nonsynonymous divergence per site; NaN if pN ≥ 0.75 or invalid.
      omega              dN/dS; NaN if dS = 0 or undefined.

    Notes
      • Headers are taken as the first whitespace-delimited token after '>' in FASTA.
      • Stop codons are excluded everywhere; codons whose shortest mutation paths necessarily cross a stop are skipped.
      • High pS (≈0.6–0.75) implies saturation; JC may inflate dS or become undefined.
      • All counts are on the retained codon set only.

    Examples
      python3 dNdS_pair.py ALIGN.fas Aparatermes_sp. Kalotermes_flavicollis
      python3 dNdS_pair.py --fasta ALIGN.fas --nameA Aparatermes_sp. --nameB Kalotermes_flavicollis
      python3 dNdS_pair.py fasta=ALIGN.fas nameA=Aparatermes_sp. nameB=Kalotermes_flavicollis
    """
    p = argparse.ArgumentParser(
        prog="dNdS_pair.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(desc),
        epilog=textwrap.dedent(epi)
    )

    # 
    p.add_argument("fasta", nargs="?", help="Codon alignment in FASTA (single gene, multiple taxa).")
    p.add_argument("nameA", nargs="?", help="Sequence ID for the first taxon (exact FASTA ID token).")
    p.add_argument("nameB", nargs="?", help="Sequence ID for the second taxon (exact FASTA ID token).")
    p.add_argument("--fasta", "-f", dest="fasta_opt")
    p.add_argument("--nameA", "-a", dest="nameA_opt")
    p.add_argument("--nameB", "-b", dest="nameB_opt")

    args = p.parse_args(norm)

    # 
    fasta = args.fasta_opt or args.fasta
    nameA = args.nameA_opt or args.nameA
    nameB = args.nameB_opt or args.nameB
    if not (fasta and nameA and nameB):
        p.error("Missing arguments: need FASTA, nameA, nameB.")

    class Obj: pass
    out = Obj()
    out.fasta, out.nameA, out.nameB = fasta, nameA, nameB
    return out

args = parse_args()
seqFile = args.fasta
nameA   = args.nameA
nameB   = args.nameB

# Bases et table des codons (sites synonymes/nonsynonymes par codon)
bases = ["A","T","G","C"]
codonTable = {}
for p1 in bases:
    for p2 in bases:
        for p3 in bases:
            codon = p1+p2+p3
            aa = Seq(codon).translate()
            if aa != "*":  # exclure les codons stops
                nS = 0.0; nN = 0.0
                for i, pos in enumerate([p1,p2,p3]):
                    for b in bases:
                        if b==pos: 
                            continue
                        c2 = list(codon)
                        c2[i] = b
                        aa2 = Seq("".join(c2)).translate()
                        if aa2 == aa:
                            nS += 1
                        else:
                            nN += 1
                codonTable[codon] = {
                    "aa": str(aa),
                    "nS": nS/3.0,
                    "nN": nN/3.0
                }

def is_valid_codon(c):
    return c in codonTable

def syn_nonsyn_changes_expected(c1, c2):
    """Retourne (S,N) attendus en moyennant sur les chemins minimaux
    qui ne passent pas par un codon stop. Si aucun chemin valide: None."""
    # c1 et c2 = deux codons sens, codants pour des AA (pas des Stops)
    # S = nombre attendu de changements nucleotidiques synonymes pour transformer c1 en c2 
    if c1 == c2:
        return 0.0, 0.0
    diffs = [i for i in range(3) if c1[i] != c2[i]]
    # d = distance de Hamming
    d = len(diffs)
    s_sum = 0.0
    n_sum = 0.0
    n_valid = 0
    for order in itertools.permutations(diffs, d):
        cur = list(c1)
        aa_prev = codonTable[c1]["aa"]
        s_count = 0.0
        n_count = 0.0
        valid = True
        for idx in order:
            cur[idx] = c2[idx]
            cod = "".join(cur)
            if cod not in codonTable: # stop rencontré -> chemin invalide
                valid = False
                break
            aa_new = codonTable[cod]["aa"]
            if aa_new == aa_prev:
                s_count += 1.0
            else:
                n_count += 1.0
            aa_prev = aa_new
        if valid:
            s_sum += s_count
            n_sum += n_count
            n_valid += 1
    if n_valid == 0:
        return None
    return s_sum / n_valid, n_sum / n_valid


def jc_correction(p):
    """Correction de Jukes–Cantor des proportions sur sites S ou N (pS et pN)"""
    # corrige une proportion observée p de differences nucleotidiques entre 2 sequences
    # en un nombre de substitutions d par site, corrigé pour les substitutions multiples
    # au même site nucleotidique, mais non observées.
    if p is None or math.isnan(p):
        return float('nan')
    # limites
    if p < 0:
        p = 0.0
    # is p>0.75: alors d est non definissable
    if p >= 0.75: # 1 - 4/3 p <= 0
        return float('nan')
    # correction
    return -0.75 * math.log(1.0 - (4.0/3.0)*p)

# Lire séquences
seqs = {}
#for rec in parse(fasta, "fasta"):
for rec in parse(seqFile, "fasta"):
    seqs[rec.id] = str(rec.seq).upper()

if nameA not in seqs or nameB not in seqs:
    sys.stderr.write("Séquence introuvable: vérifie nameA/nameB dans les IDs FASTA.\n")
    sys.exit(1)

sA = seqs[nameA]
sB = seqs[nameB]
L = min(len(sA), len(sB))
L -= L % 3
sA = sA[:L]
sB = sB[:L]

# Parcours par codon
Ls = 0.0
Ln = 0.0
S_changes = 0.0
N_changes = 0.0
n_codons_total = L // 3
n_codons_retained = 0

for i in range(0, L, 3):
    c1 = sA[i:i+3]
    c2 = sB[i:i+3]
    # filtre: gap/ambiguïté/stop → rejet
    if ("N" in c1) or ("N" in c2) or ("-" in c1) or ("-" in c2):
        continue
    if not (is_valid_codon(c1) and is_valid_codon(c2)):
        continue
    # opportunités: moyenne des sites S/N des deux codons
    Ls += 0.5*(codonTable[c1]["nS"] + codonTable[c2]["nS"])
    Ln += 0.5*(codonTable[c1]["nN"] + codonTable[c2]["nN"])
    # changements attendus S/N (chemins de longueur minimale)
    schg_n = syn_nonsyn_changes_expected(c1, c2)
    if schg_n is None:
        continue # écarte ce codon car tous les chemins passaient par un stop
    schg, nchg = schg_n
    S_changes += schg
    N_changes += nchg
    n_codons_retained += 1

# proportions sur sites
pS = float('nan') if Ls==0 else (S_changes / Ls)
pN = float('nan') if Ln==0 else (N_changes / Ln)

# corrections JC
dS = jc_correction(pS) if not math.isnan(pS) else float('nan')
dN = jc_correction(pN) if not math.isnan(pN) else float('nan')
omega = float('nan')
if dS is not None and not math.isnan(dS) and dS != 0:
    omega = dN / dS

# outprout
gene = basename(seqFile)

hdr = ["gene","nameA","nameB","codons_total","codons_retained","Ls","Ln","S_changes","N_changes","pS","pN","dS","dN","omega"]
vals = [
    gene, nameA, nameB,
    n_codons_total, n_codons_retained,
    f"{Ls:.6f}", f"{Ln:.6f}",
    f"{S_changes:.6f}", f"{N_changes:.6f}",
    f"{pS:.6f}" if not math.isnan(pS) else "NA",
    f"{pN:.6f}" if not math.isnan(pN) else "NA",
    f"{dS:.6f}" if not math.isnan(dS) else "NA",
    f"{dN:.6f}" if not math.isnan(dN) else "NA",
    f"{omega:.6f}" if not math.isnan(omega) else "NA"
]
print("\t".join(hdr))
print("\t".join(map(str, vals)))

