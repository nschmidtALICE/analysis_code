
#pPb unfolding

# measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-21_pPb/TagZHistograms_%s.root" #CBall
measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-05-14_pPb_DGaussDefault/TagZHistograms_%s.root" #DGauss
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS_response.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/6_pp_MBMC_wMult_response_sysvar0.root"
responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/11_12_pPb_EPOS_Fix1a3_response_sysvar0.root"
outfile="unfolded_output_zt_lhcb_pPb.root"

# root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'")'

#pPb with pp mc response:
# measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-05-14_pPb_DGaussDefault/TagZHistograms_%s.root" #DGauss
# measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-07-23_pPb_DGaussDefault/TagZHistograms_%s.root" #DGauss
measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-09-14_pPb_DGaussDefault/TagZHistograms_%s.root" #DGauss
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/10_pp_2016_MBMC_response_sysvar0.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_MBMC_response_sysvar0.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/26-28_pp_2016_pth_MCs_15plus_response_multweights.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/26-28_pp_2016_pth_MCs_15plus_response_noweights.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_response_for_pPb.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_response_for_pPb_jetptweights.root"
responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_addMC29-38_response_pPbWeights.root" #contains jet-pt weights
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_response_noweights.root"
outfile="unfolded_output_zt_lhcb_pPb_withppresponse.root"

root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'")'

# measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-21_Pbp/TagZHistograms_%s.root" #CBall
# measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-05-14_Pbp_DGaussDefault/TagZHistograms_%s.root" #DGauss
# measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-07-23_Pbp_DGaussDefault/TagZHistograms_%s.root" #DGauss
measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-09-14_Pbp_DGaussDefault/TagZHistograms_%s.root" #DGauss
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/74_EPOS_Fix4_Pbp_response.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/74_75_EPOS_Pbp_response.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/15_16_Pbp_EPOS_Fix1a4_response_sysvar0.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_MBMC_response_sysvar0_for_Pbp.root" #default pp response for pPb unfolding
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_response_noweights.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_response_for_Pbp.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_response_for_Pbp_600cut.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_response_for_Pbp_jetptweights.root"
responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_addMC29-38_response_PbpWeights.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/26-28_pp_2016_pth_MCs_15plus_response_multweights.root"
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/26-28_pp_2016_pth_MCs_15plus_response_noweights.root"
outfile="unfolded_output_zt_lhcb_Pbp_withppresponse.root"

root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'")'


#pp output:
measfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-28_DGaussDefault/TagZHistograms_%s.root" #DGauss
responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/26-28_pp_2016_pth_MCs_15plus_response_sysvar0.root"
outfile="unfolded_output_zt_lhcb_pp.root"

# root -x -l -b -q unfold_d0_zt_lhcb.cpp'("'$measfile'", "'$responsefile'", "'$outfile'")'
