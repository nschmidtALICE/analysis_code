infileMBMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS.root
infileTrigMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/55_FF_pPb_EPOS_8GeV_trigg.root
# infileMBMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_response.root
# infileTrigMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/55_FF_pPb_EPOS_8GeV_trigg_response.root


infileweights="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/output_2026-03-11/out_weights_8GeVTrigMC55.root"
root -x -l -b -q "validatepPbTrigMCweights.cpp(\"$infileMBMC\",\"$infileTrigMC\",\"$infileweights\")"