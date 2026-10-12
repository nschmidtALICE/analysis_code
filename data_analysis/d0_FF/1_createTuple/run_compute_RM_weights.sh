# infilepp=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/10_pp_2016_MBMC.txt
# infilepp=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_MBMC.txt
# infilepp=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC.txt
infilepp=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_addMC29-38.txt
infilepPb=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/11_12_pPb_EPOS_Fix1a3.txt
output_basename="RM_weights_pPb"


root -x -l -b -q "compute_RM_weights.cpp+(\"${infilepp}\", \"${infilepPb}\", \"${output_basename}\")"

infilePbp=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/15_16_Pbp_EPOS_Fix1a4.txt
output_basename="RM_weights_Pbp"


root -x -l -b -q "compute_RM_weights.cpp+(\"${infilepp}\", \"${infilePbp}\", \"${output_basename}\")"

