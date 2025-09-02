#!/usr/bin/env bash
# run_yn00_pair.sh ALIGN.fas nameA nameB
set -euo pipefail

ALN="$1"; A="$2"; B="$3"
WD="$(mktemp -d)"
trap 'rm -rf "$WD"' EXIT

# 1) Extraire les deux séquences par ID = 1er token du header FASTA
awk -v a="$A" -v b="$B" '
  BEGIN{keep=0}
  /^>/{
    id=$2; gsub(/^>/,"",id)            # id = 1er token après ">"
    split($0, toks, /[ \t]/); id=substr(toks[1],2)
    keep=(id==a || id==b)
  }
  keep{print}
' "$ALN" > "$WD/pair.fasta"

# Sanity: exactement 2 séquences
grep -c "^>" "$WD/pair.fasta" | grep -qx 2

# 2) FASTA -> PHYLIP (séquentiel), noms ≤10 caractères (sans espaces)
python3 - "$WD/pair.fasta" "$WD/pair.phy" << 'PY'
import sys
fasta, out = sys.argv[1], sys.argv[2]
seqs = []
with open(fasta) as fh:
    name, seq = None, []
    for line in fh:
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: seqs.append((name,''.join(seq).upper()))
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if name is not None: seqs.append((name,''.join(seq).upper()))
if len(seqs)!=2:
    sys.exit("Erreur: il faut exactement 2 séquences.")
L = min(len(seqs[0][1]), len(seqs[1][1]))
L -= L%3
seqs = [(n[:10].replace(' ','_'), s[:L]) for n,s in seqs]
if L==0:
    sys.exit("Erreur: alignement effectif de longueur 0.")
with open(out,'w') as oh:
    oh.write(f"{len(seqs)} {L}\n")
    for n,s in seqs:
        oh.write(f"{n:<10} {s}\n")
PY

# 3) yn00.ctl sans commentaires
cat > "$WD/yn00.ctl" << 'CTL'
      seqfile = pair.phy
     outfile = yn00.out
       verbose = 1
       icode = 0
   weighting = 0
   commonf3x4 = 0
       ndata = 1
CTL

# 4) Lancer yn00 en montrant les messages
( cd "$WD" && yn00 yn00.ctl )

# 5) Vérifier la sortie non vide
test -s "$WD/yn00.out"

# 6) Copier et extraire un résumé Ka/Ks/omega (méthode YN00)
cp "$WD/yn00.out" ./yn00.out

awk '
  BEGIN{meth=""}
  /Yang & Nielsen/ {meth="YN00"}
  meth=="YN00" && $1=="t" {getline; print; exit}
' "$WD/yn00.out" \
| awk 'NF{printf "YN00_summary: %s\n", $0}' || true

echo "OK: yn00.out"

