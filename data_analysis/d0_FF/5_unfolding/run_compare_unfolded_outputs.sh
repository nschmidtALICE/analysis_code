infile1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_zT_pPb_2026-02-26_weightedPP69to72/pPb_unfolded_output.root"
infile2="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_zT_pPb_2026-02-26_54MC/pPb_unfolded_output.root"

root -x -b -l -q 'compare_unfolded_outputs.C+("'"$infile1"'","'"$infile2"'","pp","pPb","cmp_pp_vs_pPb","unfolded_zT",6,true)'