BASE_DIR="/media/niviths/local/analysis_code/data_analysis/d0_FF"

echo "Running Stage 1: Create ntuples"
cd "$BASE_DIR/1_createTuple"
    
#if one argument is given, use this to set the SYSVAR variable, otherwise default to 0 (nominal)
if [ "$#" -eq 1 ]; then
	SYSVAR=$1
else
	SYSVAR=0
fi


RESPONLY_INT=0  # Set to 1 if you want response for data
DORESP_INT=1    # Set to 1 if you want response for MC
DO_JET=1
JET_MODE=1

MC_INT=0

pPbORPbp="pPb"  # Change to "Pbp" for Pbp data
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/57_FF_pPb_DATA.root"
root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE}, ${SYSVAR})"

pPbORPbp="Pbp"  # Change to "Pbp" for Pbp data
inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/59_FF_Pbp_DATA.root"
root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE}, ${SYSVAR})"

# MC_INT=1
# RESPONLY_INT=1  # Set to 1 if you want response for data
# pPbORPbp="pPb"  # Change to "Pbp" for Pbp data
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/pPb_MC_54plus73_EPOS.txt"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE}, ${SYSVAR})"

# pPbORPbp="Pbp"  # Change to "Pbp" for Pbp data
# inputFile="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/Pbp_MC_74_75_EPOS.txt"
# root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", \"${pPbORPbp}\", ${MC_INT}, ${RESPONLY_INT}, ${DORESP_INT}, ${DO_JET}, ${JET_MODE}, ${SYSVAR})"
