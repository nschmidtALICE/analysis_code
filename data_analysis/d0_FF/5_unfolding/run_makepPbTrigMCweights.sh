infileMBMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/73_EPOS_Pbp.txt
infileTrigMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/74_EPOS_Trigg_Pbp.txt
# infileMBMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS.root
# infileTrigMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/55_FF_pPb_EPOS_8GeV_trigg.root
# infileMBMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_response.root
# infileTrigMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/55_FF_pPb_EPOS_8GeV_trigg_response.root


outfile="out_weights_Pbp_74_73MC.root"
root -x -l -b -q "makepPbTrigMCweights.cpp(\"$infileMBMC\",\"$infileTrigMC\",\"$outfile\")"