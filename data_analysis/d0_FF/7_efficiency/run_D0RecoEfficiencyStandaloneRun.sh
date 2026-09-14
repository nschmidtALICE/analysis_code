# # #int D0RecoEfficiencyStandaloneRun(
# #     // TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/51/51.root,/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/52/52.root,/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/53/53.root", //this is Pbp
# #     TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root", //this is pPb
# #                                 //  TString outputFile = "output_reco_standalone_full_Pbp.root",
# #                                  TString outputFile = "output_reco_standalone_full_pPb.root",
# #                                  double massWindow = 50.0, double minPt = 1.0,
# #                                  double minEta = 2.0, double maxEta = 4.5,
# #                                  double kaonPIDCut = 0.5, double pionPIDCut = 0.5,
# #                                  bool makePlots = true) {

# # inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/73_EPOS_Fix3_pPb.txt"
# inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS.txt"
# outputFile="output_reco_standalone_54p73_pPb.root"

# # root -x -l -b -q D0RecoEfficiencyStandaloneRun.cpp'+("'$inputFiles'","'$outputFile'",50.0,1.0,2.0,4.5,0.5,0.5,true)'

# inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/74_EPOS_Fix4_Pbp.txt"
# outputFile="output_reco_standalone_74_Pbp.root"
inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/Pbp_MC_74_75_EPOS.txt"
outputFile="output_reco_standalone_74_75_Pbp_rerun.root"
# inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/15_16_Pbp_EPOS_Fix1a4.txt"
# outputFile="output_reco_standalone_15_16_Pbp_rerun.root"

root -x -l -b -q D0RecoEfficiencyStandaloneRun.cpp'+("'$inputFiles'","'$outputFile'",35.0,1.0,2.0,4.5,0.5,0.5,true)'

# inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS.root"
# inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS.txt"
# outputFile="output_reco_standalone_54p73_pPb.root"
inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/11_12_pPb_EPOS_Fix1a3.txt"
outputFile="output_reco_standalone_11_12_pPb_rerun.root"

# root -x -l -b -q D0RecoEfficiencyStandaloneRun.cpp'+("'$inputFiles'","'$outputFile'",35.0,1.0,2.0,4.5,0.5,0.5,true)'



# inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/20251104_Pbp_EPOS_62/20251104_Pbp_EPOS_62.txt"
# outputFile="output_reco_standalone_full_Pbp_inclusive.root"

# #root -x -l -b -q D0RecoEfficiencyStandaloneRun.cpp'+("'$inputFiles'","'$outputFile'",50.0,1.0,2.0,4.5,0.5,0.5,true)'


# inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/1_pp_MC.txt"
# outputFile="output_reco_standalone_full_pp.root"

# # root -x -l -b -q D0RecoEfficiencyStandaloneRun.cpp'+("'$inputFiles'","'$outputFile'",50.0,1.0,2.0,4.5,0.5,0.5,true)'
