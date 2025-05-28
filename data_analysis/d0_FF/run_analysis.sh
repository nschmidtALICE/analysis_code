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

# Function to display help
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
    echo "  all - Run all stages in sequence"
    echo ""
    echo "Options:"
    echo "  -h, --help         Display this help message"
    echo "  -f, --file FILE    Set input ROOT file path (default: $INPUT_FILE)"
    echo "  -r, --resonance R  Set resonance (default: $RESONANCE)"
    echo "  -z, --zt MODE      Set zT mode True/False (default: $ZT_MODE)"
    echo "  -m, --mcsim        Set MC mode (default: $MC_MODE)"
    echo ""
    echo "Example: $0 -p 10 -r D0 all"
}

# Parse command line arguments
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
    python ntuplesJetAnalysis_Response_noVect.py -f "$INPUT_FILE" -p 10 -m "$MC_MODE"
    if [ $? -ne 0 ]; then
        echo "Error: Stage 1 failed"
        exit 1
    fi
    echo "Stage 1 completed successfully"
}

# Function to run stage 2 - Fit data
run_stage2() {
    echo "Running Stage 2: Fit data"
    cd "$BASE_DIR/2_fitData"
    inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250501_Pbp_data/Pbp_data_filterV1.root"
    isMCswitch="false"
    if [ "$MC_MODE" = "True" ]; then
        inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_17_MC_output_D0FF/1109236/20250514_Pbp_MC_output_D0FF.root"
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
    python plot_Response_NoVec.py
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
    #create a list of pt ranges containing 5_10, 10_15, 15_20, 20_30, 30_50
    PT_RANGES=("5_10")
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
    root -x -l -b -q plotRawYields.C'+("")' #, '$ZT_MODE')'
    echo "Stage 4 completed successfully"
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
    all)
        echo "Running complete analysis pipeline"
        run_stage1
        run_stage2
        run_stage3
        run_stage4
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