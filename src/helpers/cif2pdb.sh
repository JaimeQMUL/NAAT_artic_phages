#!/usr/bin/env bash
set -e

INPUT="$1"
ACC="$2"


OUT="uvsx/data/gmx_prep/structures/$ACC/$ACC.pdb"

echo "Converting: $ACC"

conda run -n cold maxit -input "$INPUT" -output "$OUT" -o 2

