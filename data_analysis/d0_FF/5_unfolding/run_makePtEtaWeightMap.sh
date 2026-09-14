# infilepp=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/70_72_response.root
infilepp=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/69to72_response.root
# infilepp=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/69_2018sim_D02Kpi_pthgreater15_response.root
infilepPb=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_response.root
outfile="out_weights_52vs69to72.root"

root -x -l -b -q "makePtEtaWeightMap.cpp(\"$infilepp\",\"$infilepPb\",\"$outfile\")"