#!/bin/bash
# d0_analysis.sh - Script to automate D0 fragmentation function analysis steps

# Run a specific stage: ./d0_analysis.sh 1  # Run ntuple creation
# Run all stages: ./d0_analysis.sh all
# Customize parameters: ./d0_analysis.sh -p 20 -r D0 -z True 4  # Run stage 4 with custom parameters


# Define paths and default values
BASE_DIR="/media/niviths/local/analysis_code/data_analysis/d0_FF"
INPUT_FILE="/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_17_MC_output_D0FF/1109236/20250514_Pbp_MC_output_D0FF.root"
RESONANCE="D0"
ZT_MODE="True"
MC_MODE="False"
UNFOLD_VAR="1"  # Default unfolding variation

# Function to display help
# Update in the show_help function
show_help() {
    echo "D0 Fragmentation Function Analysis Pipeline"
    echo ""
    echo "Usage: $0 [options] stage"
    echo ""
    echo "Stages:"
    echo "  1 - Create ntuples (ntuplesJetAnalysis_Response_noVect.py)"
    echo "  2 - Fit data (MassFitterScript_TestReduce.py)"
    echo "  3 - Plot response matrix (plot_Response_NoVec.py)"
    echo "  4 - Plot raw yields (plot_RawSignalYields.py)"
    echo "  5 - Run unfolding (unfoldzT.C)"
    echo "  all - Run all stages in sequence"
    echo ""
    echo "Options:"
    echo "  -h, --help         Display this help message"
    echo "  -f, --file FILE    Set input ROOT file path (default: $INPUT_FILE)"
    echo "  -r, --resonance R  Set resonance (default: $RESONANCE)"
    echo "  -z, --zt MODE      Set zT mode True/False (default: $ZT_MODE)"
    echo "  -m, --mcsim        Set MC mode (default: $MC_MODE)"
    echo "  -u, --unfold VAR   Set unfolding variation 0-4 (default: 0)"
    echo ""
    echo "Example: $0 -p 10 -r D0 all"
}

# Update in the argument parsing section
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -f|--file)
            INPUT_FILE="$2"
            shift 2
            ;;
        -r|--resonance)
            RESONANCE="$2"
            shift 2
            ;;
        -z|--zt)
            ZT_MODE="$2"
            shift 2
            ;;
        -m|--mcsim)
            MC_MODE="$2"
            shift 2
            ;;
        -u|--unfold)
            UNFOLD_VAR="$2"
            shift 2
            ;;
        *)
            STAGE="$1"
            shift
            ;;
    esac
done


# Check if stage is specified
if [ -z "$STAGE" ]; then
    echo "Error: No stage specified"
    show_help
    exit 1
fi

# Function to run stage 1 - Create ntuples
run_stage1() {
    echo "Running Stage 1: Create ntuples"
    cd "$BASE_DIR/1_createTuple"
    
    # Get MC mode as integer for ROOT
    MC_INT=1
    if [ "$MC_MODE" = "False" ]; then
        MC_INT=0
    fi
    # Define input files as an array
    declare -a inputFiles=(
        "/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_pPb_11273plus.root"
        "/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_Pbp_1127plus.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_Pbp_11276a7plus.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/51/51.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/52/52.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/53/53.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250609_merged/1123981.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250609_merged/1122665.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_17_MC_output_D0FF/20250514_Pbp_MC_output_D0FF.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_18_20_MC_output_D0FF/20250514_Pbp_18_20_MC_output_D0FF.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_1.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_01.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_2.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_02.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_03.root"
        # "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF/20250514_Pbp_21_MC_output_D0FF_04.root"
    )
    
    
    # Loop through the input files and run the script for each file
    for inputFile in "${inputFiles[@]}"; do
        echo "Processing input file: $inputFile"
        # Properly escape the file path and pass the MC mode
        root -x -l -b -q "nTupleMaker.C+(\"${inputFile}\", ${MC_INT})"
        # root -x -l -b -q "nTupleCreator.C+(\"${inputFile}\", ${MC_INT})"
        # python ntuplesJetAnalysis_Response_noVect.py -f "$inputFile" -p 10 -m "$MC_INT"
   
        if [ $? -ne 0 ]; then
            echo "Error: Stage 1 failed for input file $inputFile"
            exit 1
        fi
    done
    
    echo "Stage 1 completed successfully"
}

# Function to run stage 2 - Fit data
run_stage2() {
    echo "Running Stage 2: Fit data"
    cd "$BASE_DIR/2_fitData"
    inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_Pbp_1127plus_filtered.root"
    # inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_Pbp_11276a7plus_filtered.root"
    # inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250616_newGanga_DATA/merged_pPb_11273plus_filtered.root"
    # inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250501_Pbp_data/Pbp_data_filterV1.root"
    isMCswitch="false"
    if [ "$MC_MODE" = "True" ]; then
        # inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/filtered.root"
        inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/filtered.root"
        # inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250609_merged/filtered.root"
        # inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250609_Pbp_MC_output/1122665/1122665_filtered.root"
        # inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF_filterV1_bunew.root"
        isMCswitch="true"
    fi
    root -x -l -b -q MassFitter.C'+("'$inputFileMassFit'",'$isMCswitch')'
    # python3 MassFitterScript_TestReduce.py -m "$MC_MODE"
    if [ $? -ne 0 ]; then
        echo "Error: Stage 2 failed"
        exit 1
    fi
    echo "Stage 2 completed successfully"
}

# Function to run stage 3 - Plot response matrix
run_stage3() {
    echo "Running Stage 3: Plot response matrix"
    cd "$BASE_DIR/3_makeResponseMatrix"
    # infileresp=/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_allMC.root
    # infileresp=/media/niviths/SSD2/lhcb_analysis_SSD/20250609_Pbp_MC_output/1122665/1122665_response.root
    #python plot_Response_NoVec.py
    root -x -l -b -q plotResponse.C'+()' #"'$infileresp'")'
    if [ $? -ne 0 ]; then
        echo "Error: Stage 3 failed"
        exit 1
    fi
    echo "Stage 3 completed successfully"
}

# Function to run stage 4 - Plot raw yields
run_stage4() {
    echo "Running Stage 4: Plot raw yields"
    cd "$BASE_DIR/4_plotRawYields"
    isMCswitch="false"
    if [ "$MC_MODE" = "True" ]; then
        isMCswitch="true"
    fi
    #create a list of pt ranges containing 5_10, 10_15, 15_20, 20_30, 30_50
    # PT_RANGES=("5_10")
    # PT_RANGES=("5_10" "10_15" "15_20" "20_30" "30_50")
    #loop through the pt ranges and run the script for each range
    # for PT_RANGE in "${PT_RANGES[@]}"; do
    #     echo "Running for pt range: $PT_RANGE"
    #     root -x -l -b -q plotRawYields.C'+("'$PT_RANGE'")' #, '$ZT_MODE')'
    #     if [ $? -ne 0 ]; then
    #         echo "Error: Stage 4 failed for pt range $PT_RANGE"
    #         exit 1
    #     fi
    # done
    root -x -l -b -q plotRawYields.C'+("", false, '$isMCswitch')' #, '$ZT_MODE')'
    # root -x -l -b -q plotRawYields.C'+("", true, '$isMCswitch')' #, '$ZT_MODE')'
    echo "Stage 4 completed successfully"
}

# Function to run stage 5 - Run unfolding
run_stage5() {
    echo "Running Stage 5: Unfolding"
    cd "$BASE_DIR/5_unfolding"
    
    # Define variation options: 
    # 0 = Default, 1 = SWP_prior, 2 = Flat_prior, 3 = moreHyperons, 4 = lessHyperons
    VARIATION=0
    
    if [ -n "$UNFOLD_VAR" ]; then
        VARIATION="$UNFOLD_VAR"
    fi
    
    echo "Using unfolding variation: $VARIATION"
    echo "  0 = Default"
    echo "  1 = Swap prior between prompt/non-prompt"
    echo "  2 = Flat prior"
    echo "  3 = More hyperons"
    echo "  4 = Less hyperons"
    
    # Run unfolding with the selected variation
    root -x -l -b -q unfoldzT.C'+('$VARIATION')'

    
    if [ $? -ne 0 ]; then
        echo "Error: Stage 5 (unfolding) failed"
        exit 1
    fi
    
    echo "Stage 5 completed successfully"
}

# Execute stage(s) based on input
case $STAGE in
    1)
        run_stage1
        ;;
    2)
        run_stage2
        ;;
    3)
        run_stage3
        ;;
    4)
        run_stage4
        ;;
    5)
        run_stage5
        ;;
    all)
        echo "Running complete analysis pipeline"
        run_stage1
        run_stage2
        run_stage3
        run_stage4
        run_stage5
        echo "All stages completed successfully"
        ;;
    *)
        echo "Error: Invalid stage '$STAGE'"
        show_help
        exit 1
        ;;
esac
# Exit successfully
echo "Analysis pipeline completed successfully"

exit 0
