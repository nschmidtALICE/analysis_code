BASE_DIR="/media/niviths/local/analysis_code/data_analysis/d0_FF"

echo "Running Stage 1: Create ntuples"
cd "$BASE_DIR/1_createTuple"
    

RESPONLY_INT=0  # Set to 1 if you want response for data
DORESP_INT=1    # Set to 1 if you want response for MC
DO_JET=1
JET_MODE=1

MC_INT=0

pPbORPbp="pPb"  # Change to "Pbp" for Pbp data
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_pPb_11273plus.root"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/57_FF_pPb_DATA.root"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/9_pPb_data.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/17_pPb_data_fixedassoc.txt"
root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"
DO_JET=1
JET_MODE=2
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/20251022_pPb_Data_allD0_rerun/20251022_pPb_Data_allD0_rerun.root"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/20251016_pPb_Data_allD0/20251016_pPb_data.root"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${DO_JET}, ${JET_MODE})"

DO_JET=1
JET_MODE=1
pPbORPbp="Pbp"  # Change to "Pbp" for Pbp data
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_Pbp_1127plus.root"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/59_FF_Pbp_DATA.root"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/8_Pbp_data.txt"
root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"
DO_JET=1
JET_MODE=2
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/20251023_Pbp_Data_allD0_rerun/20251023_Pbp_Data_allD0_rerun.root"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"


RESPONLY_INT=1  # Set to 1 if you want response for data
DORESP_INT=1    # Set to 1 if you want response for MC
DO_JET=1
JET_MODE=1
MC_INT=1
pPbORPbp="pPb"  # Change to "Pbp" for Pbp data
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_EPOS_Fix_45_wMult.txt"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/11_12_pPb_EPOS_Fix1a3.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/2025_53_61_outputs/54.root"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/2025_53_61_outputs/55.root"
#root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/73_EPOS_Pbp.txt"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"


pPbORPbp="Pbp"  # Change to "Pbp" for Pbp data
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/Pbp_MC_merged.root"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/2025_53_61_outputs/53.root"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/74_EPOS_Fix4_Pbp.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/75_EPOS_Fix1_Pbp.root"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/Pbp_MC_74_75_EPOS.txt"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/15_16_Pbp_EPOS_Fix1a4.txt"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"

inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/2025_53_61_outputs/53.root"
#root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"

pPbORPbp="pp"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/68/ntuple_test-12757408.root"
#root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"

inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/69_2018sim_D02Kpi_pthgreater15.root"
#root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/70_2016sim_dijetc _15to20pth.root"
#root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/71_2016sim_dijetc_20to50pth.root"
#root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/72_2016sim_dijetc_pthgr50.root"
#root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"


#inclusive D0 ntuple creation with jet info for EPOS MC
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/20251104_Pbp_EPOS_62/20251104_Pbp_EPOS_62.txt" #"/media/niviths/SSD2/lhcb_analysis_SSD/20251104_Pbp_EPOS_62/20251104_Pbp_EPOS_62.root"
DO_JET=1
JET_MODE=2
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"



pPbORPbp="pp"
# inputFile="/media/niviths/local/tst/localRunning_MC_ntuple_test.root"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/1_pp_MC.txt"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/10_pp_2016_MBMC.txt"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"



RESPONLY_INT=0  # Set to 1 if you want response for data
DORESP_INT=1    # Set to 1 if you want response for MC
DO_JET=1
JET_MODE=1

MC_INT=0

pPbORPbp="pp"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/0_pp_Data.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/0_part_pp_Data.txt"
inputFile="/media/niviths/718CC09A0FF48654/0_pp_data_full.txt"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE})"

DO_JET=1
JET_MODE=1
MC_INT=1
RESPONLY_INT=1  # Set to 1 if you want response for data

pPbORPbp="pp"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/10_pp_2016_MBMC.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/22_pp_2017_MB_MC.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/23_pp_2018_MB_MC.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_MBMC.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/26-28_pp_2016_pth_MCs_15plus.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_addMC.txt"
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_full_2016-2018_w2016pthMC_addMC29-38.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pp_magup_2016-2018_MB_MB.txt"
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/1_pp_MC.txt"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE}, 0, false, true)"



	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_Pbp_1127plus.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/20250728_pPb_MC_output.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_Pbp_11276a7plus.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged.root"
	# "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/51/51.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/52/52.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/53/53.root"
	# "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/54/54.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/55/55.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250728_pPb_MC_output/56/56.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250609_merged/1123981.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250609_merged/1122665.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_17_MC_output_D0FF/20250514_Pbp_MC_output_D0FF.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_18_20_MC_output_D0FF/20250514_Pbp_18_20_MC_output_D0FF.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_1.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_01.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_2.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_02.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_03.root"
	# # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_04.root"
    
