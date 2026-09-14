# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-04-22_DGaussDefault/unfolded_output_zt_lhcb_Pbp.root","10_15,15_20,20_30","eta0,eta1,eta2")'
# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-04-22_DGaussDefault/unfolded_output_zt_lhcb_pPb.root","10_15,15_20,20_30","eta0,eta1,eta2")'

# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-04-22_DGauss/unfolded_output_zt_lhcb_Pbp.root","10_15,15_20,20_30","eta0,eta1,eta2")' //default
# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-04-22_DGauss/unfolded_output_zt_lhcb_pPb.root","10_15,15_20,20_30","eta0,eta1,eta2")' //default


##################                              default                      ##################
# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-05-14/unfolded_output_zt_lhcb_Pbp.root","10_15,15_20,20_30","eta0,eta1,eta2")'
# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-05-14/unfolded_output_zt_lhcb_pPb.root","10_15,15_20,20_30","eta0,eta1,eta2")'
##################                              default                      ##################


# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_plots_2026-05-13/unfolded_output_zt_lhcb_Pbp.root","10_15,15_20,20_30,30_100","eta0,eta1,eta2")'
# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_plots_2026-05-13/unfolded_output_zt_lhcb_pPb.root","10_15,15_20,20_30,30_100","eta0,eta1,eta2")'


#pp output
# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pp_plots_2026-04-29/unfolded_output_zt_lhcb_pp.root","10_15,15_20,20_30","eta0,eta1,eta2")'

specialadd="" #contains jet pT weights depending on pPb or Pbp response, and also contains the latest code changes for kinematic efficiency calculation
# specialadd="_final1" #contains jet pT weights depending on pPb or Pbp response, and also contains the latest code changes for kinematic efficiency calculation
# specialadd="_jetptweights_addMC_kinEff"
# specialadd="_jetptweights_addMC"
# specialadd="_jetptweights"
# specialadd="_noweights"
root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_withppresponse_plots_2026-09-14'$specialadd'/unfolded_output_zt_lhcb_Pbp_withppresponse.root","10_15,15_20,20_30","eta0,eta1,eta2")'
root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_withppresponse_plots_2026-09-14'$specialadd'/unfolded_output_zt_lhcb_pPb_withppresponse.root","10_15,15_20,20_30","eta0,eta1,eta2")'
# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_Pbp_withppresponse_plots_2026-07-23'$specialadd'/unfolded_output_zt_lhcb_Pbp_withppresponse.root","10_15,15_20,20_30","eta0,eta1,eta2")'
# root -l -b -q 'plot_unfold_results.cpp("/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/unfolded_output_zt_lhcb_pPb_withppresponse_plots_2026-07-23'$specialadd'/unfolded_output_zt_lhcb_pPb_withppresponse.root","10_15,15_20,20_30","eta0,eta1,eta2")'
