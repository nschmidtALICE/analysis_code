#pPb
responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS_response.root"
outfile="unfolded_output_zt_lhcb_pPb.root"

measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_pPb_DGaussDefault/TagZHistograms_%s.root" #DGauss
outputFolderTag="DGaussDefault"
root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'", "'$outputFolderTag'")'
measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_pPb_DGaussTightPID/TagZHistograms_%s.root" #DGauss
outputFolderTag="DGaussTightPID"
root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'", "'$outputFolderTag'")'
measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_pPb_DGaussLoosePID/TagZHistograms_%s.root" #DGauss
outputFolderTag="DGaussLoosePID"
root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'", "'$outputFolderTag'")'


#Pbp
outfile="unfolded_output_zt_lhcb_Pbp.root"
responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/74_75_EPOS_Pbp_response.root"

measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_Pbp_DGauss/TagZHistograms_%s.root" #DGauss
outputFolderTag="DGaussDefault"
root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'", "'$outputFolderTag'")'
measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_Pbp_DGaussTightPID/TagZHistograms_%s.root" #DGauss
outputFolderTag="DGaussTightPID"
root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'", "'$outputFolderTag'")'
measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-22_Pbp_DGaussLoosePID/TagZHistograms_%s.root" #DGauss
outputFolderTag="DGaussLoosePID"
root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'", "'$outputFolderTag'")'
