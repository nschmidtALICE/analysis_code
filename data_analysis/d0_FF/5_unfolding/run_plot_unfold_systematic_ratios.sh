# yield extraction systematic, DGauss vs CBall
infiledefault="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-04-22_DGauss/unfolded_output_zt_lhcb_pPb.root"
infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-04-21_CBallFix/unfolded_output_zt_lhcb_pPb.root"
#root -l -b -q 'plot_unfold_systematic_ratios.cpp+("'$infiledefault','$infilevar1'","default,CBall","D0MassFitFunction")'

infiledefault="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-04-22_DGauss/unfolded_output_zt_lhcb_Pbp.root"
infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-04-21_CBallFix/unfolded_output_zt_lhcb_Pbp.root"
#root -l -b -q 'plot_unfold_systematic_ratios.cpp+("'$infiledefault','$infilevar1'","default,CBall","D0MassFitFunction")'

# PID systematic, var 1 vs var 2
infiledefaultpPb="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-04-22_DGaussDefault/unfolded_output_zt_lhcb_pPb.root"
infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-04-22_DGaussTightPID/unfolded_output_zt_lhcb_pPb.root"
infilevar2="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-04-22_DGaussLoosePID/unfolded_output_zt_lhcb_pPb.root"
# root -l -b -q 'plot_unfold_systematic_ratios.cpp+("'$infiledefaultpPb','$infilevar1','$infilevar2'","default,TightPID,LoosePID","PIDSelection")'

infiledefaultPbp="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-04-22_DGaussDefault/unfolded_output_zt_lhcb_Pbp.root"
infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-04-22_DGaussTightPID/unfolded_output_zt_lhcb_Pbp.root"
infilevar2="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-04-22_DGaussLoosePID/unfolded_output_zt_lhcb_Pbp.root"
# root -l -b -q 'plot_unfold_systematic_ratios.cpp+("'$infiledefaultPbp','$infilevar1','$infilevar2'","default,TightPID,LoosePID","PIDSelection")'


infiledefaultPbp="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-05-14/unfolded_output_zt_lhcb_Pbp.root"
# infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_withppresponse_plots_2026-05-28_weights/unfolded_output_zt_lhcb_Pbp_withppresponse.root"
infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_withppresponse_plots_2026-06-02_jetptweights_addMC/unfolded_output_zt_lhcb_Pbp_withppresponse.root"
infilevar2="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_withppresponse_plots_2026-05-28_jetptweights/unfolded_output_zt_lhcb_Pbp_withppresponse.root"
# root -l -b -q 'plot_unfold_systematic_ratios.cpp+("'$infiledefaultPbp','$infilevar1','$infilevar2'","PbpMC,ppMCwJetPtWeightsAddMC,ppMCwJetPtWeights","PbpResponseVariationJetPtWeights")'



infiledefaultPbp="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-05-14/unfolded_output_zt_lhcb_pPb.root"
# infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_withppresponse_plots_2026-05-28_weights/unfolded_output_zt_lhcb_pPb_withppresponse.root"
infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_withppresponse_plots_2026-06-02_jetptweights_addMC/unfolded_output_zt_lhcb_pPb_withppresponse.root"
infilevar2="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_withppresponse_plots_2026-05-28_jetptweights/unfolded_output_zt_lhcb_pPb_withppresponse.root"
# root -l -b -q 'plot_unfold_systematic_ratios.cpp+("'$infiledefaultPbp','$infilevar1','$infilevar2'","pPbMC,ppMCwJetPtWeightsAddMC,ppMCwJetPtWeights","pPbResponseVariationJetPtWeights")'



infiledefaultPbp="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_withppresponse_plots_2026-06-02_jetptweights_addMC/unfolded_output_zt_lhcb_Pbp_withppresponse.root"
infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_withppresponse_plots_2026-06-04_jetptweights_addMC/unfolded_output_zt_lhcb_Pbp_withppresponse.root"
infilevar2="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_withppresponse_plots_2026-06-04/unfolded_output_zt_lhcb_Pbp_withppresponse.root"
root -l -b -q 'plot_unfold_systematic_ratios.cpp+("'$infiledefaultPbp','$infilevar1','$infilevar2'","def,firstchange,kinemeff","CodeChangesandKinematicEffPbp")'



infiledefaultpPb="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_withppresponse_plots_2026-06-02_jetptweights_addMC/unfolded_output_zt_lhcb_pPb_withppresponse.root"
infilevar1="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_withppresponse_plots_2026-06-04_jetptweights_addMC/unfolded_output_zt_lhcb_pPb_withppresponse.root"
infilevar2="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_withppresponse_plots_2026-06-04/unfolded_output_zt_lhcb_pPb_withppresponse.root"
root -l -b -q 'plot_unfold_systematic_ratios.cpp+("'$infiledefaultpPb','$infilevar1','$infilevar2'","def,firstchange,kinemeff","CodeChangesandKinematicEffpPb")'

