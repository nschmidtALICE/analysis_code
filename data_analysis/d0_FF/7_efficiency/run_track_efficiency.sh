#!/bin/bash

# Track Reconstruction Efficiency Analysis Script
# This script compiles and runs the track efficiency analysis

echo "=== Track Reconstruction Efficiency Analysis ==="
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Configuration
INPUT_FILE="/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/51/51.root"
OUTPUT_FILE="track_efficiency_results.root"
MAKEFILE="Makefile_TrackEff"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Warning: Input file not found: $INPUT_FILE"
    echo "Using default file path in code..."
fi

# Clean previous builds
echo "Cleaning previous builds..."
make -f $MAKEFILE clean

# Build the code
echo ""
echo "Building track efficiency analysis..."
make -f $MAKEFILE all

# Check if build was successful
if [ $? -ne 0 ]; then
    echo "Error: Build failed!"
    exit 1
fi

echo ""
echo "Build successful!"

# Run the analysis
echo ""
echo "Running track efficiency analysis..."
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo ""

# Create ROOT macro to run the analysis
cat > run_track_efficiency.C << EOF
{
    // Load the compiled library
    gSystem->Load("./TrackRecoEfficiencyRun");
    
    // Run the analysis
    int result = TrackRecoEfficiencyRun(
        "$INPUT_FILE",     // input file
        "$OUTPUT_FILE",    // output file
        0.5,               // minimum pT (GeV)
        2.0,               // minimum eta
        5.0,               // maximum eta
        0.3,               // maximum ghost probability
        3.0,               // maximum track chi2
        true               // make plots
    );
    
    if (result == 0) {
        cout << "\\n=== Track efficiency analysis completed successfully! ===" << endl;
    } else {
        cout << "\\n=== Track efficiency analysis failed! ===" << endl;
    }
    
    gApplication->Terminate(result);
}
EOF

# Run with ROOT
echo "Executing analysis..."
root -l -b -q run_track_efficiency.C

# Check results
if [ $? -eq 0 ]; then
    echo ""
    echo "=== Analysis completed successfully! ==="
    echo ""
    echo "Output files generated:"
    ls -la *.root *.png *.pdf 2>/dev/null | head -10
    echo ""
    echo "Key results:"
    echo "- Track efficiency histograms saved to: $OUTPUT_FILE"
    echo "- Efficiency plots saved as PNG and PDF files"
    echo "- 2D efficiency maps (pT vs eta) for pions and kaons"
    echo "- Momentum-based efficiency maps (p vs eta)"
    echo ""
    echo "To view results:"
    echo "  root -l $OUTPUT_FILE"
    echo "  # Then browse histograms in ROOT TBrowser"
else
    echo ""
    echo "=== Analysis failed! ==="
    echo "Check the error messages above for details."
    exit 1
fi

# Clean up temporary files
rm -f run_track_efficiency.C

echo "Script completed successfully!"
