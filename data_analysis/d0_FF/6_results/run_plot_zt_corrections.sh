#!/usr/bin/env bash
set -euo pipefail

# Loop over jet pT ranges and rapidity bins (0..7)
BASE_DIR="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA"
pt_pairs=("5_10" "10_15" "15_20" "20_30" "30_50")
for pt in "${pt_pairs[@]}"; do
	for bin in {0..7}; do
		infile="$BASE_DIR/${pt}/tagZAllEfficiencyCorrections_bin${bin}.root"
		outprefix="${pt}_tagZAllEfficiencyCorrections_bin${bin}"
		echo "== Processing: ${infile} -> ${outprefix} =="
		if [[ -f "$infile" ]]; then
			root -l -b -q "plot_zt_corrections.C(\"${infile}\",\"${outprefix}\")"
		else
			echo "  Skipping: input file not found: ${infile}"
		fi
	done
done

echo "All requested files processed."
