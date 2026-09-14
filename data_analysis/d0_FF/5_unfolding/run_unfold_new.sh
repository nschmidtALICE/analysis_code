# #run unfolding for pPb data
# root -x -l -b -q unfold_new.cpp'+("pPb_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-01-26_pPb/TagZHistograms_%s.root", "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_rename_response.root")' #, 6, {"5_10", "10_15", "15_20", "20_30"}, {0,1,2,3,4,5,6,7}, false, false)'

# # run unfolding for Pbp data
# root -x -l -b -q unfold_new.cpp'+("Pbp_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-01-26_Pbp/TagZHistograms_%s.root", "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_rename_response.root")' #, 6, {"5_10", "10_15", "15_20", "20_30"}, {0,1,2,3,4,5,6,7}, false, false)'

infileresponse=/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS_response.root
# root -x -l -b -q unfold_new.cpp'+("pPb_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-03-31_pPb/TagZHistograms_%s.root", "'$infileresponse'")' #, 6, {"5_10", "10_15", "15_20", "20_30"}, {0,1,2,3,4,5,6,7}, false, false)'
root -x -l -b -q unfold_new.cpp'+("pPb_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-04-01_pPb/TagZHistograms_%s.root", "'$infileresponse'","",1.0, 6, {"10_15", "15_20", "20_30", "30_100"}, {0,1,2}, {2.5,3.0,3.5,4.0}, false, true)'

# run unfolding for Pbp data
# root -x -l -b -q unfold_new.cpp'+("Pbp_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-03-19_Pbp/TagZHistograms_%s.root", "'$infileresponse'")' #, 6, {"5_10", "10_15", "15_20", "20_30"}, {0,1,2,3,4,5,6,7}, false, false)'


# ## UNFOLDING WITH MB PLUS TRIGGERED RESPONSE
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_rename_response.root"
# responsefiletrigg="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/55_FF_pPb_EPOS_8GeV_trigg_response.root"
# verbositysetting=1
# #run unfolding for pPb data
# root -x -l -b -q unfold_new.cpp'+("pPb_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-01-26_pPb/TagZHistograms_%s.root", "'$responsefile'","'$responsefiletrigg'")' #, 6, {"5_10", "10_15", "15_20", "20_30", "30_50"}, {0,1,2,3,4,5,6,7}, false, '$verbositysetting')'

# #run unfolding for Pbp data
# root -x -l -b -q unfold_new.cpp'+("Pbp_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-01-26_Pbp/TagZHistograms_%s.root", "'$responsefile'","'$responsefiletrigg'")' #, 6, {"5_10", "10_15", "15_20", "20_30", "30_50"}, {0,1,2,3,4,5,6,7}, false, '$verbositysetting')'




## --------------------------------------------------------------------
## --------------------------------------------------------------------
## --------------------------------------------------------------------
## --------------------------------------------------------------------
##OLD OLD OLD UNFOLDING WITH PP RESPONSE
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/69to72_response.root"
# # responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/merged_response.root"
# # responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/69_2018sim_D02Kpi_pthgreater15_response.root"
# # weightfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/out_weights_52vs69to72_plots_2026-03-03/out_weights_52vs69to72.root" #2d weights based on pt and eta
# weightfile="/media/niviths/local/analysis_code/data_analysis/d0_FF/5_unfolding/out_weights_52vs69to72_plots_2026-03-09/out_weights_52vs69to72.root" #3d weights based on pt, eta and zt
# # weighthistoname="hW_ptEta_det_pPbOverPP"
# weighthistoname="hW_ptEta_d0mc_pPbOverPP"
# verbositysetting=1
# #run unfolding for pPb data
# root -x -l -b -q unfold_new.cpp'+("pPb_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-01-26_pPb/TagZHistograms_%s.root", "'$responsefile'", 6, {"5_10", "10_15", "15_20", "20_30", "30_50"}, {0,1,2,3,4,5,6,7}, false, '$verbositysetting', "'$weightfile'", "'$weighthistoname'")'

# #run unfolding for Pbp data
# root -x -l -b -q unfold_new.cpp'+("Pbp_unfolded_output.root", "/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-01-26_Pbp/TagZHistograms_%s.root", "'$responsefile'", 6, {"5_10", "10_15", "15_20", "20_30", "30_50"}, {0,1,2,3,4,5,6,7}, false, '$verbositysetting', "'$weightfile'", "'$weighthistoname'")'
# responsefile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/69to72_response.root"
