infileppData=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/merged_response.root
infileppMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/merged_response.root
infilepPbData=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_response.root
infilepPbMC=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_response.root
outfile="out_weights_54vs69plus.root"
root -x -l -b -q "makeMCweightsUnfolding.cpp(\"$infileppData\",\"$infileppMC\",\"$infilepPbData\",\"$infilepPbMC\",\"$outfile\")"