# #int D0RecoEfficiencyStandaloneRun(
#     // TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/51/51.root,/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/52/52.root,/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/53/53.root", //this is Pbp
#     TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root", //this is pPb
#                                 //  TString outputFile = "output_reco_standalone_full_Pbp.root",
#                                  TString outputFile = "output_reco_standalone_full_pPb.root",
#                                  double massWindow = 50.0, double minPt = 1.0,
#                                  double minEta = 2.0, double maxEta = 4.5,
#                                  double kaonPIDCut = 0.5, double pionPIDCut = 0.5,
#                                  bool makePlots = true) {

inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS.root"
outputFile="output_reco_standalone_54_pPb.root"

root -x -l -b -q D0RecoEfficiencyStandaloneRun.cpp'+("'$inputFiles'","'$outputFile'",50.0,1.0,2.0,4.5,0.5,0.5,true)'



inputFiles="/media/niviths/SSD2/lhcb_analysis_SSD/20251104_Pbp_EPOS_62/20251104_Pbp_EPOS_62.txt"
outputFile="output_reco_standalone_full_Pbp_inclusive.root"

#root -x -l -b -q D0RecoEfficiencyStandaloneRun.cpp'+("'$inputFiles'","'$outputFile'",50.0,1.0,2.0,4.5,0.5,0.5,true)'
