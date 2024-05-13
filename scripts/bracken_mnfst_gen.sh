#!/bin/bash
data_dir="$1"

for file in "$data_dir"/*.br; do echo -e "$file\t$(basename -- "$file" .br)" >> bracken_manifest.tsv; done
