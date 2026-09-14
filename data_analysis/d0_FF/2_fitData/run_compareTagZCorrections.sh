#!/usr/bin/env bash
# Run compareTagZCorrections.C with configurable arguments

BASE_DIR="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_pPb_DGaussDefault"
VAR_DIRS_CSV="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_pPb_DGaussLoosePID,/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_pPb_DGaussTightPID"          # comma-separated directories (relative or absolute)
VAR_LABELS_CSV="Default,LoosePID,TightPID"        # optional comma-separated labels for variations

echo "Running compareTagZCorrections with"
echo "  base:    $BASE_DIR"
echo "  variants: $VAR_DIRS_CSV"
echo "  labels:  $VAR_LABELS_CSV"


# Run in ROOT batch (ACLiC compiled)
root -x -l -b -q "compareTagZCorrections.C+(\"$BASE_DIR\",\"$VAR_DIRS_CSV\",\"$VAR_LABELS_CSV\")"

echo "Done. Check output folder under $BASE_DIR/TagZComparison_outputs"
