/**
 * Track Ambiguity Finder for Magnet Station Detector
 * ------------------------------------------------
 *
 * This code analyzes charged particle trajectories through a magnetic field to study
 * track reconstruction ambiguities in a multi-layer detector system. It takes simulated
 * particle data, propagates the particles through a realistic magnetic field, and compares
 * the reconstructed tracks with the known truth information.
 *
 * Primary functionality:
 * 1. Loads simulated particle data from ROOT files
 * 2. Configures magnetic field maps for accurate particle propagation
 * 3. Propagates particles through the field using Runge-Kutta integration
 * 4. Identifies detector hits corresponding to propagated particles
 * 5. Tracklets hits across detector layers into potential track candidates
 * 6. Matches hit patterns against physically valid segment combinations
 * 7. Performs linear track fitting with iterative outlier removal
 * 8. Evaluates reconstruction quality by comparing to truth information
 * 9. Produces comprehensive diagnostic plots and visualizations
 *
 * Key features:
 * - Processing of time-coincident hits across detector layers
 * - Pattern recognition based on predefined valid segment combinations
 * - Track fitting with sophisticated outlier rejection
 * - Momentum-variation studies to assess magnetic field sensitivity
 * - Comprehensive visualization of track ambiguities and quality metrics
 *
 * Input:
 * - ROOT files containing simulated particle trajectories and detector hits
 * - Magnetic field maps for particle propagation
 *
 * Output:
 * - PDF plots showing detector pattern analysis, track fitting quality,
 *   position matching accuracy, and momentum sensitivity
 * - Optional trajectory trees for detailed offline analysis
 *
 * Usage: track_ambiguity_finder([bool doTree = false])
 * - If doTree is true, creates additional ROOT file with trajectory information
 *
 * Dependencies:
 * - ROOT data analysis framework
 * - track_propagatorRK.C for particle propagation through magnetic fields
 * - FieldLoaderAndPlotter.cpp for magnetic field handling
 *
 * Author: Nicolas Schmidt (nicolas.schmidt@cern.ch / nschmidtALICE on github)
 * Date: 2025-03-05
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLegend.h>

/*this code requires the track_propagator code as well as the FieldLoaderAndPlotter code
 * the first piece of code contains the Runge-Kutta 4th order integration method for the particle propagation
 * this code can also be be run as a standalone code to simulate particle propagation through a magnetic field
 * the second piece of code contains the magnetic field map loader and plotter for the LHCb magnet
 * this code can also be run as a standalone code to visualize the magnetic field map
 */
#include "track_propagatorRK.C" // FieldLoaderAndPlotter.cpp is included in this file
#include "bar_positions_handler.h"
/**
 * Valid segment combinations for a 4-layer detector.
 * Each vector {a,b,c,d} represents a possible combination of segments
 * where a = segment index in layer 0
 *       b = segment index in layer 1
 *       c = segment index in layer 2
 *       d = segment index in layer 3
 * These combinations represent physically possible track patterns.
 */
std::vector<std::vector<int>> vec_segment_combinations = {
    {0, 0, 0, 0},
    {1, 1, 1, 0},
    {1, 1, 1, 1},
    {2, 2, 2, 1},
    {2, 2, 2, 2},
    {2, 3, 2, 2},
    {3, 3, 3, 2},
    {3, 3, 3, 3},
    {3, 4, 3, 3},
    {3, 4, 4, 3},
    {4, 4, 4, 3},
    {4, 4, 4, 4},
    {4, 5, 4, 4},
    {4, 5, 5, 4},
    {5, 5, 5, 4},
    {5, 5, 5, 5},
    {5, 6, 5, 5},
    {5, 6, 6, 5},
    {6, 6, 6, 5},
    {6, 6, 6, 6},
    {6, 7, 6, 6},
    {6, 7, 7, 6},
    {7, 7, 7, 6},
    {7, 7, 7, 7}};
    // {2, 2, 2, 2},
    // {2, 3, 2, 2},
    // {3, 3, 3, 2},
    // {3, 3, 3, 3},
    // {3, 4, 3, 3},
    // {3, 4, 4, 3},
    // {4, 4, 4, 3},
    // {4, 4, 4, 4},
    // {4, 5, 4, 4},
    // {4, 5, 5, 4},
    // {5, 5, 5, 4},
    // {5, 5, 5, 5},
    // {5, 6, 5, 5},
    // {5, 6, 6, 5},
    // {6, 6, 6, 5},
    // {6, 6, 6, 6},
    // {6, 7, 6, 6},
    // {6, 7, 7, 6},
    // {7, 7, 7, 6}};
/**
 * Module rotation angles in radians.
 * Each value corresponds to the rotation angle for a specific module.
 */
std::vector<float> vec_module_rotations = {1.2, 1.2, 1.25, 1.25, 1.25, 1.3, 1.3, 1.35, 1.35};

/**
 * Detector geometry parameters (in millimeters).
 */
// Half-width of a detector module
double module_halfwidth = 25;

// Distance between centers of adjacent bars in a layer
double layer_center_between_bars = 1.5 + 2 + 1 + 4; // = 8.5

/**
 * Position offsets for the four detector layers.
 * These values define the x-positions of each layer relative to the center.
 */
double layer_offsets[4] = {
    -module_halfwidth - layer_center_between_bars, // Layer 0: -33.5
    -module_halfwidth + layer_center_between_bars, // Layer 1: -16.5
    module_halfwidth - layer_center_between_bars,  // Layer 2: 16.5
    module_halfwidth + layer_center_between_bars   // Layer 3: 33.5
};
/**
 * Main function to find track ambiguities in the magnet station.
 *
 * This function analyzes simulated particle tracks passing through a magnetic field,
 * identifies detector hits, reconstructs track segments, and evaluates track ambiguities.
 * It produces diagnostic plots showing the reconstruction performance.
 *
 * @param doTree Flag to enable trajectory tree creation for detailed analysis
 * @return Status code (0 for success)
 */
int track_ambiguity_finder_clusterizedGhosts(bool isPbPb = true, bool doTree = false)
{
    //----------------------------------------------------------------------
    // Configuration and initialization
    //----------------------------------------------------------------------
    bool verboseoutput = false; // Set to true for detailed debugging output

    // Input simulation file
    TString inputFile = "";
    if (isPbPb)
    {
        std::cout << "Running PbPb configuration" << std::endl;
        inputFile = "/home/niviths/Downloads/magnetStationSims/20250311_PbPb_addUTinfo/20250311_PbPb_addUTinfo.root";
    }
    else
    {
        std::cout << "Running pp configuration" << std::endl;
        // inputFile = "/home/niviths/Downloads/magnetStationSims/20250311_pp_newOutput/20250311_pp_newOutput.root";
        inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250328_pp_newSegmentation/20250328_pp_newSegmentation.root";
    }

    TFile fin(inputFile);
    TTree *ntup = (TTree *)fin.Get("ntup");
    std::cout << "Input file contains " << ntup->GetEntries() << " entries" << std::endl;

    // Magnetic field configuration files
    TString q1File = "./field_inputs_LHCb/field.v5r0.c1.down.cdf";
    TString q2File = "./field_inputs_LHCb/field.v5r0.c2.down.cdf";
    TString q3File = "./field_inputs_LHCb/field.v5r0.c3.down.cdf";
    TString q4File = "./field_inputs_LHCb/field.v5r0.c4.down.cdf";

    // Load magnetic field from files
    auto field = loadLHCbMagneticField(q1File.Data(), q2File.Data(), q3File.Data(), q4File.Data());
    std::string output_file_magnet = "lhcb_field_visualization2.root";
    visualizeLHCbField(q1File.Data(), q2File.Data(), q3File.Data(), q4File.Data(), output_file_magnet);

    // Load bar positions from file
    loadBarPositions("bar_positions.txt");

    // Create particle propagator with magnet boundaries
    double x_limit = 1.5;    // 1.5 m magnet half-width
    double y_limit = 2.0;    // 2.0 m magnet half-height
    double step_size = 0.01; // Step size for propagation (m)
    ParticlePropagator propagator(field, step_size, x_limit, y_limit);

    // Create output directory for plots
    TString outputdir = "ambiguity_finder_clusterized_plots_pp/";
    if (isPbPb)
        outputdir = "ambiguity_finder_clusterized_plots_PbPb/";
    system("mkdir -p " + outputdir);

    // Set matching windows for hit pattern analysis
    double matchingWindow = 1000;      // Matching window size in mm
    int maxBarDifference = 50;         // Maximum allowed bar difference
    int maxSegmentDifferencePlus = 2;  // Maximum allowed segment difference
    int maxSegmentDifferenceMinus = 1; // Maximum allowed segment difference
    int maxTimeDifference = 5;         // Maximum allowed time difference [ns]
    double maxAngleDifference = 0.5;   // was 0.25;  // Maximum allowed angle difference [rad]
    if (isPbPb)
    {
        maxBarDifference = 50;
        maxSegmentDifferencePlus = 2;
        maxSegmentDifferenceMinus = 1;
        maxTimeDifference = 5;
        maxAngleDifference = 0.5; // was 0.3;
    }

    //----------------------------------------------------------------------
    // Set up tree reader and branches for input data
    //----------------------------------------------------------------------
    TTreeReader tree(ntup);

    // Particle properties at production
    TTreeReaderArray<float> p(tree, "p");   // Momentum magnitude (MeV/c)
    TTreeReaderArray<float> vx(tree, "vx"); // Vertex x-coordinate (mm)
    TTreeReaderArray<float> vy(tree, "vy"); // Vertex y-coordinate (mm)
    TTreeReaderArray<float> vz(tree, "vz"); // Vertex z-coordinate (mm)
    TTreeReaderArray<float> px(tree, "px"); // Momentum x-component (MeV/c)
    TTreeReaderArray<float> py(tree, "py"); // Momentum y-component (MeV/c)
    TTreeReaderArray<float> pz(tree, "pz"); // Momentum z-component (MeV/c)
    TTreeReaderArray<int> pid(tree, "pid"); // Particle ID (PDG code)
    TTreeReaderArray<int> key(tree, "key"); // Unique particle identifier

    // Magnet Station (MS) hit information
    TTreeReaderArray<float> ms_vx(tree, "ms_vx");       // MS hit x-coordinate (mm)
    TTreeReaderArray<float> ms_vy(tree, "ms_vy");       // MS hit y-coordinate (mm)
    TTreeReaderArray<float> ms_vz(tree, "ms_vz");       // MS hit z-coordinate (mm)
    TTreeReaderArray<float> ms_px(tree, "ms_px");       // MS momentum x-component (MeV/c)
    TTreeReaderArray<float> ms_py(tree, "ms_py");       // MS momentum y-component (MeV/c)
    TTreeReaderArray<float> ms_pz(tree, "ms_pz");       // MS momentum z-component (MeV/c)
    TTreeReaderArray<float> ms_time(tree, "ms_time");   // MS hit time (ns)
    TTreeReaderArray<int> ms_id(tree, "ms_id");         // MS hit particle ID
    TTreeReaderArray<int> ms_bitID(tree, "ms_segment"); // MS hit segment identifier

    // Upstream Tracker (UT) hit information
    TTreeReaderArray<float> ut_vx(tree, "ut_vx");     // UT hit x-coordinate (mm)
    TTreeReaderArray<float> ut_vy(tree, "ut_vy");     // UT hit y-coordinate (mm)
    TTreeReaderArray<float> ut_vz(tree, "ut_vz");     // UT hit z-coordinate (mm)
    TTreeReaderArray<float> ut_px(tree, "ut_px");     // UT momentum x-component (MeV/c)
    TTreeReaderArray<float> ut_py(tree, "ut_py");     // UT momentum y-component (MeV/c)
    TTreeReaderArray<float> ut_pz(tree, "ut_pz");     // UT momentum z-component (MeV/c)
    TTreeReaderArray<float> ut_time(tree, "ut_time"); // UT hit time (ns)
    TTreeReaderArray<int> nUThits(tree, "nUThits");   // Number of UT hits per track

    //----------------------------------------------------------------------
    // Track fitting and visualization containers
    //----------------------------------------------------------------------
    TGraphErrors *gSlopeFits_orig[30]; // Original track fits before outlier removal
    TF1 *funcSlopeFits_orig[30];       // Original track fit functions
    TGraphErrors *gSlopeFits_true[30]; // Track fits using only truth-matched hits
    TF1 *funcSlopeFits_true[30];       // Truth track fit functions
    TGraphErrors *gSlopeFits[30];      // Final track fits after outlier removal
    TF1 *funcSlopeFits[30];            // Final track fit functions
    int nSlopeFits = 0;                // Counter for track fits

    // Containers for clustered hits
    std::vector<int> clusterizedHits_bitID;  // Segment IDs for clustered hits
    std::vector<float> clusterizedHits_time; // Timing information for clustered hits
    std::vector<int> clusterizedHits_id;     // Particle IDs for clustered hits

    //----------------------------------------------------------------------
    // Histograms for performance analysis
    //----------------------------------------------------------------------
    // Track pattern and bar difference analysis
    TH2D *hRelativeBarDifferenceToFirstLayer = new TH2D(
        "hRelativeBarDifferenceToFirstLayer",
        "Relative bar difference to first layer (truth)",
        3, 0.5, 3.5, 21, -10.5, 10.5);

    TH1D *hRelativeBarDifferenceToFirstLayerIndiv[3];
    for (int i = 0; i < 3; i++)
    {
        hRelativeBarDifferenceToFirstLayerIndiv[i] = new TH1D(
            Form("hRelativeBarDifferenceToFirstLayerIndiv_%d", i + 1),
            Form("Relative bar difference to first layer (truth) for layer %d", i + 1),
            21, -10.5, 10.5);
    }

    // Track angle and tracklet analysis
    TH2D *hAngleDiffModules = new TH2D(
        "hAngleDiffModules",
        "Angle difference between MS tracklet and true hit direction",
        9, 0.5, 9.5, 100, -0.01, 0.4);

    TH2D *hAngleDiffFinalModules = new TH2D(
        "hAngleDiffFinalModules",
        "Angle difference between MS tracklet and extrapolated track",
        9, 0.5, 9.5, 100, -0.01, 3.2);

    TH1D *hNTracklets = new TH1D(
        "hNTracklets",
        "Number of tracklets",
        10, -0.50, 9.5);

    TH1D *hNTrackletsTrue = new TH1D(
        "hNTrackletsTrue",
        "Number of tracklets with truth",
        10, -0.50, 9.5);

    TH1D *hNTrackletsSlope = new TH1D(
        "hNTrackletsSlope",
        Form("Number of tracklets with tracklet fit (d#Phi < %1.1f rad)", maxAngleDifference),
        10, -0.50, 9.5);

    TH1D *hNTrackletsSlopeTrue = new TH1D(
        "hNTrackletsSlopeTrue",
        Form("Number of tracklets with tracklet fit (d#Phi < %1.1f rad) and truth hits", maxAngleDifference),
        10, -0.50, 9.5);

    TH1D *hNTrackletsTime = new TH1D(
        "hNTrackletsTime",
        "Number of tracklets with time cut",
        10, -0.50, 9.5);

    TH1D *hNTrackletsTimeTrue = new TH1D(
        "hNTrackletsTimeTrue",
        "Number of tracklets with time cut and truth hits",
        10, -0.50, 9.5);

    TH1D *hNTrackletsSlopeAndTimeCut = new TH1D(
        "hNTrackletsSlopeAndTimeCut",
        Form("Number of tracklets with tracklet fit (d#Phi < %1.1f rad) and time cut", maxAngleDifference),
        10, -0.50, 9.5);

    TH1D *hNTrackletsSlopeAndTimeCutTrue = new TH1D(
        "hNTrackletsSlopeAndTimeCutTrue",
        Form("Number of tracklets with tracklet fit (d#Phi < %1.1f rad) and time cut and truth hits", maxAngleDifference),
        10, -0.50, 9.5);

    // Track matching quality analysis
    double matchWindow = 1000; // Matching window size in mm
    TH2D *hMatchdYdZ = new TH2D(
        "hMatchdYdZ",
        "Matched dYdZ",
        100, -matchWindow / 5, matchWindow / 5,
        100, -matchWindow, matchWindow);

    TH2D *hMatchdXdZ = new TH2D(
        "hMatchdXdZ",
        "Matched dXdZ",
        100, -matchWindow, matchWindow,
        100, -matchWindow, matchWindow);

    TH1D *hMatchdX = new TH1D(
        "hMatchdX",
        "Matched dX",
        100, -matchWindow, matchWindow);

    TH1D *hMatchdY = new TH1D(
        "hMatchdY",
        "Matched dY",
        100, -matchWindow / 5, matchWindow / 5);

    TH1D *hMatchdZ = new TH1D(
        "hMatchdZ",
        "Matched dZ",
        100, -matchWindow, matchWindow);

    //make hMatchdX, hMatchdY, hMatchdZ for each of the 9 modules
    TH1D *hMatchdXModule[9];
    TH1D *hMatchdYModule[9];
    TH1D *hMatchdZModule[9];
    for (int i = 0; i < 9; i++)
    {
        hMatchdXModule[i] = new TH1D(
            Form("hMatchdXModule%d", i),
            Form("Matched dX for module %d", i),
            100, -matchWindow, matchWindow);
        hMatchdYModule[i] = new TH1D(
            Form("hMatchdYModule%d", i),
            Form("Matched dY for module %d", i),
            100, -matchWindow / 5, matchWindow / 5);
        hMatchdZModule[i] = new TH1D(
            Form("hMatchdZModule%d", i),
            Form("Matched dZ for module %d", i),
            100, -matchWindow, matchWindow);
    }
        //make hMatchdX, hMatchdY, hMatchdZ for each of the 4 stations
    TH1D *hMatchdXStation[4];
    TH1D *hMatchdYStation[4];
    TH1D *hMatchdZStation[4];
    for (int i = 0; i < 4; i++)
    {
        hMatchdXStation[i] = new TH1D(
            Form("hMatchdXStation%d", i),
            Form("Matched dX for station %d", i),
            100, -matchWindow, matchWindow);
        hMatchdYStation[i] = new TH1D(
            Form("hMatchdYStation%d", i),
            Form("Matched dY for station %d", i),
            100, -matchWindow / 5, matchWindow / 5);
        hMatchdZStation[i] = new TH1D(
            Form("hMatchdZStation%d", i),
            Form("Matched dZ for station %d", i),
            100, -matchWindow, matchWindow);
    }

    TH1D *hMatchdXExtrapolToTracklet = new TH1D(
        "hMatchdXExtrapolToTracklet",
        "Matched dX extrapolated to tracklet",
        100, -matchWindow, matchWindow);

    TH1D *hMatchdYExtrapolToTracklet = new TH1D(
        "hMatchdYExtrapolToTracklet",
        "Matched dY extrapolated to tracklet",
        100, -matchWindow / 5, matchWindow / 5);

    TH1D *hMatchdZExtrapolToTracklet = new TH1D(
        "hMatchdZExtrapolToTracklet",
        "Matched dZ extrapolated to tracklet",
        100, -matchWindow, matchWindow);

    TH1D *hMatchdXExtrapolToTrackletQuality = new TH1D(
        "hMatchdXExtrapolToTrackletQuality",
        "Matched dX extrapolated to tracklet with quality cut",
        100, -matchWindow, matchWindow);
    TH1D *hMatchdYExtrapolToTrackletQuality = new TH1D(
        "hMatchdYExtrapolToTrackletQuality",
        "Matched dY extrapolated to tracklet with quality cut",
        100, -matchWindow / 5, matchWindow / 5);
    TH1D *hMatchdZExtrapolToTrackletQuality = new TH1D(
        "hMatchdZExtrapolToTrackletQuality",
        "Matched dZ extrapolated to tracklet with quality cut",
        100, -matchWindow, matchWindow);

    // make the same histograms as the five above for after application of cuts
    TH2D *hMatchdYdZCut = new TH2D(
        "hMatchdYdZCut",
        "Matched dYdZ after cuts",
        100, -matchWindow / 5, matchWindow / 5,
        100, -matchWindow, matchWindow);

    TH2D *hMatchdXdZCut = new TH2D(
        "hMatchdXdZCut",
        "Matched dXdZ after cuts",
        100, -matchWindow, matchWindow,
        100, -matchWindow, matchWindow);

    TH1D *hMatchdXCut = new TH1D(
        "hMatchdXCut",
        "Matched dX after cuts",
        100, -matchWindow, matchWindow);

    TH1D *hMatchdYCut = new TH1D(
        "hMatchdYCut",
        "Matched dY after cuts",
        100, -matchWindow / 5, matchWindow / 5);

    TH1D *hMatchdZCut = new TH1D(
        "hMatchdZCut",
        "Matched dZ after cuts",
        100, -matchWindow, matchWindow);

    // Track position analysis
    TH2D *hTrueHitXZ = new TH2D(
        "hTrueHitXZ",
        "Final XZ",
        200, -3000, 3000, 200, 4700, 7500);

    TH2D *hTrueHitYZ = new TH2D(
        "hTrueHitYZ",
        "Final YZ",
        100, -3000, 3000, 100, 4500, 8500);

    TH2D *hExtrapolatedXZ = new TH2D(
        "hExtrapolatedXZ",
        "Extrapolated XZ",
        100, -3000, 3000, 100, 4500, 8500);

    TH2D *hExtrapolatedYZ = new TH2D(
        "hExtrapolatedYZ",
        "Extrapolated YZ",
        100, -3000, 3000, 100, 4500, 8500);

    TH1D *hExtrapolatedTimeVsTrue = new TH1D(
        "hExtrapolatedTimeVsTrue",
        "Extrapolated time minus true time",
        200, -5, 5);
    TH1D *hExtrapolatedTimeVsTracklet = new TH1D(
        "hExtrapolatedTimeVsTracklet",
        "Extrapolated time minus tracklet time",
        200, -5, 5);

    TH1D *hDiffBarPositionToTrueHitPositionX = new TH1D(
        "hDiffBarPositionToTrueHitPositionX",
        "Difference between bar position and true hit position X",
        200, -1000, 1000);
    TH1D *hDiffBarPositionToTrueHitPositionY = new TH1D(
        "hDiffBarPositionToTrueHitPositionY",
        "Difference between bar position and true hit position Y",
        200, -1000, 1000);
    TH1D *hDiffBarPositionToTrueHitPositionZ = new TH1D(
        "hDiffBarPositionToTrueHitPositionZ",
        "Difference between bar position and true hit position Z",
        200, -1000, 1000);

    // Momentum smearing analysis
    double maxDiffExtrapolation = 2000; // Maximum extrapolation difference in mm
    TH2D *hExtrapolationMomentumSmearingdXdZ = new TH2D(
        "hExtrapolationMomentumSmearingdXdZ",
        "Extrapolation Momentum Smearing dXdZ",
        100, -maxDiffExtrapolation, maxDiffExtrapolation,
        100, -5, maxDiffExtrapolation);

    TH2D *hExtrapolationMomentumSmearingdYdZ = new TH2D(
        "hExtrapolationMomentumSmearingdYdZ",
        "Extrapolation Momentum Smearing dYdZ",
        100, -maxDiffExtrapolation, maxDiffExtrapolation,
        100, -5, maxDiffExtrapolation);

    //----------------------------------------------------------------------
    // Set up trajectory output tree
    //----------------------------------------------------------------------
    TFile *rootFileOut = nullptr;
    TTree *treeOut = nullptr;

    // Variables to store in the tree
    int tr_particle_id;                     // Unique ID for each particle
    double tr_x, tr_y, tr_z;                // Position coordinates (m)
    double tr_px_out, tr_py_out, tr_pz_out; // Momentum components (GeV/c)
    double tr_energy;                       // Total energy (GeV)
    double tr_mass_value;                   // Particle mass (GeV/c²)
    double tr_charge_value;                 // Particle charge (e)
    int tr_step_num;                        // Step number along trajectory
    bool tr_hit_boundary;                   // Flag indicating boundary hit

    if (doTree)
    {
        rootFileOut = new TFile("propagation_multiparticle.root", "RECREATE");
        treeOut = new TTree("trajectories", "Multiple Particle Trajectories");

        // Create branches for all variables in the tree
        treeOut->Branch("particle_id", &tr_particle_id, "particle_id/I");
        treeOut->Branch("step", &tr_step_num, "step/I");
        treeOut->Branch("x", &tr_x, "x/D");
        treeOut->Branch("y", &tr_y, "y/D");
        treeOut->Branch("z", &tr_z, "z/D");
        treeOut->Branch("px", &tr_px_out, "px/D");
        treeOut->Branch("py", &tr_py_out, "py/D");
        treeOut->Branch("pz", &tr_pz_out, "pz/D");
        treeOut->Branch("energy", &tr_energy, "energy/D");
        treeOut->Branch("mass", &tr_mass_value, "mass/D");
        treeOut->Branch("charge", &tr_charge_value, "charge/D");
        treeOut->Branch("hit_boundary", &tr_hit_boundary, "hit_boundary/O");
    }
    //----------------------------------------------------------------------
    // Event loop - process each event in the input tree
    //----------------------------------------------------------------------
    int numEvt = 0; // Event counter

    // Pre-allocate memory for hit collections to avoid frequent reallocations
    clusterizedHits_bitID.reserve(10000);
    clusterizedHits_time.reserve(10000);
    clusterizedHits_id.reserve(10000);

    // Prepare hit mapping to speed up searches
    std::unordered_map<int, std::vector<size_t>> keyToHitIndices;
    std::vector<size_t> candidateHitIndices;
    candidateHitIndices.reserve(5000);
    int maxbartemp = 0;

    while (tree.Next())
    {
        // Limit the number of events processed for faster development/testing
        if (isPbPb)
        {
            if (numEvt > 10)
                break;
        }
        else
        {
            if (numEvt > 500)
                break;
        }

        std::cout << "Event " << numEvt << " with " << pid.GetSize() << " particles" << std::endl;
        numEvt++;

        // Build a mapping of particle keys to hit indices once per event
        keyToHitIndices.clear();
        for (size_t j = 0; j < ms_vz.GetSize(); ++j)
        {
            // Only include hits in valid z-range and not in support structure
            if ((ms_vz[j] >= 3000 && ms_vz[j] <= 7700) &&
                !((ms_bitID[j] >> 28 & 0x3) || (ms_bitID[j] >> 30 & 0x3)))
            {
                keyToHitIndices[ms_id[j]].push_back(j);
                //TODO remove or comment this part, it's for debugging
                int barIDtemp = ms_bitID[j] >> 22 & 0x3F;
                if(barIDtemp > maxbartemp)
                    maxbartemp = barIDtemp;
                TVector3 positionOnBar = getBarPosition(ms_bitID[j] >> 8 & 0x7, ms_bitID[j] >> 11 & 0xF, ms_bitID[j] >> 15 & 0x7, ms_bitID[j] >> 18 & 0xF, barIDtemp);
                //only fill for station 0 module 2
                if ((ms_bitID[j] >> 8 & 0x7) == 0 && (ms_bitID[j] >> 11 & 0xF) == 2)
                {
                    hDiffBarPositionToTrueHitPositionX->Fill(positionOnBar.X() - ms_vx[j]);
                    hDiffBarPositionToTrueHitPositionY->Fill(positionOnBar.Y() - ms_vy[j]);
                    hDiffBarPositionToTrueHitPositionZ->Fill(positionOnBar.Z() - ms_vz[j]);
                }
            }
        }

        //----------------------------------------------------------------------
        // Process each particle in the event
        //----------------------------------------------------------------------
        for (size_t i = 0; i < pid.GetSize(); ++i)
        {
            // Clear hit collections for this particle
            clusterizedHits_bitID.clear();
            clusterizedHits_time.clear();
            clusterizedHits_id.clear();

            // Apply selection criteria - do quick rejection tests first
            if (nUThits[i] < 3 || p[i] > 5000)
            {
                continue; // Skip particles with too few tracker hits or too high momentum
            }

            // Skip if particle has no associated hits (fast early rejection)
            if (keyToHitIndices.find(key[i]) == keyToHitIndices.end())
            {
                continue;
            }

            //----------------------------------------------------------------------
            // Set up particle properties for propagation - avoid recalculating constants
            //----------------------------------------------------------------------
            // Calculate charge based on particle type
            const double particle_charge = (pid[i] == 11 || pid[i] == 13) ? -pid[i] / fabs(pid[i]) : pid[i] / fabs(pid[i]);

            // Use pion mass as default
            const double particle_mass = 0.13957039; // Charged pion mass in GeV/c²

            // Configure propagator with particle properties
            propagator.setParticleProperties(particle_charge, particle_mass);

            //----------------------------------------------------------------------
            // Prepare initial state for propagation
            //----------------------------------------------------------------------
            // Initial position from upstream tracker (convert mm to m)
            const TVector3 initial_position(ut_vx[i] / 1000, ut_vy[i] / 1000, ut_vz[i] / 1000);

            // Initial time
            const double initial_time = ut_time[i];

            // Create momentum vector using direction and magnitude
            TVector3 direction(ut_px[i], ut_py[i], ut_pz[i]);
            direction = direction.Unit();

            // Calculate total energy from momentum and mass (convert MeV to GeV)
            const double momentum_magnitude = p[i] / 1000.0; // Convert MeV/c to GeV/c
            const double total_energy = sqrt(momentum_magnitude * momentum_magnitude +
                                             particle_mass * particle_mass);

            // Create momentum vector (convert MeV/c to GeV/c)
            const TVector3 momentum = direction * momentum_magnitude;

            // Create four-momentum vectors for all three cases at once
            TLorentzVector four_momentum(momentum.X(), momentum.Y(), momentum.Z(), total_energy);

            // Only create momentum variations if needed by propagation
            const double momentum_variation = 0.1;
            const double momentum_variation_value = momentum_variation * momentum_magnitude;
            const TVector3 momentum_variation_vector = direction * momentum_variation_value;

            //----------------------------------------------------------------------
            // Propagate particle through magnetic field - run in parallel if possible
            //----------------------------------------------------------------------
            std::vector<State> trajectory = propagator.propagate(
                initial_position, four_momentum, initial_time, 8.0, 50000);

            // debug output
            if (false)
            {
                std::cout << "Final position and time: " << trajectory.back().position.X() << " " << trajectory.back().position.Y() << " " << trajectory.back().position.Z() << " " << trajectory.back().time << std::endl;
            }

            if (trajectory.back().position.Z() > 7.7 || trajectory.back().position.Z() < 3.0)
            {
                continue;
            }

            // Skip particles that don't make it through propagation
            if (trajectory.size() < 2)
                continue;

            double final_time = trajectory.back().time;

            // Create variations for momentum uncertainty studies (±10%)
            const TVector3 momentum_plus = momentum + momentum_variation_vector;
            const TVector3 momentum_minus = momentum - momentum_variation_vector;

            TLorentzVector four_momentum_plus(
                momentum_plus.X(), momentum_plus.Y(), momentum_plus.Z(), total_energy);
            TLorentzVector four_momentum_minus(
                momentum_minus.X(), momentum_minus.Y(), momentum_minus.Z(), total_energy);

            std::vector<State> trajectory_plus = propagator.propagate(
                initial_position, four_momentum_plus, initial_time, 8.0, 50000);
            std::vector<State> trajectory_minus = propagator.propagate(
                initial_position, four_momentum_minus, initial_time, 8.0, 50000);

            // Skip if any trajectory calculation failed
            if (trajectory_plus.size() < 2 || trajectory_minus.size() < 2)
                continue;

            //----------------------------------------------------------------------
            // Record trajectory points in output tree if enabled - moved inside particle loop
            //----------------------------------------------------------------------
            if (doTree)
            {
                tr_particle_id = i;
                tr_hit_boundary = false;
                tr_mass_value = particle_mass;
                tr_charge_value = particle_charge;

                for (size_t j = 0; j < trajectory.size(); j++)
                {
                    // Extract position and momentum from trajectory point
                    tr_x = trajectory[j].position.X();
                    tr_y = trajectory[j].position.Y();
                    tr_z = trajectory[j].position.Z();
                    tr_px_out = trajectory[j].momentum.X();
                    tr_py_out = trajectory[j].momentum.Y();
                    tr_pz_out = trajectory[j].momentum.Z();

                    // Calculate total energy: E² = p²c² + m²c⁴ - cache expensive operations
                    const double mom_sqr = tr_px_out * tr_px_out + tr_py_out * tr_py_out + tr_pz_out * tr_pz_out;
                    tr_energy = sqrt(mom_sqr + particle_mass * particle_mass);
                    tr_step_num = j;

                    // Fill tree with this trajectory point
                    treeOut->Fill();
                }
            }

            //----------------------------------------------------------------------
            // Calculate final track parameters - more efficiently compute positions
            //----------------------------------------------------------------------
            // Extract final position and direction from trajectory (convert m to mm)
            const size_t last_idx = trajectory.size() - 1;
            const size_t second_last_idx = last_idx - 1;

            const TVector3 final_position(
                trajectory[last_idx].position.X() * 1000,
                trajectory[last_idx].position.Y() * 1000,
                trajectory[last_idx].position.Z() * 1000);

            // Calculate track direction from last two points
            const double x1 = trajectory[second_last_idx].position.X() * 1000;
            const double y1 = trajectory[second_last_idx].position.Y() * 1000;
            const double z1 = trajectory[second_last_idx].position.Z() * 1000;
            const double x2 = trajectory[last_idx].position.X() * 1000;
            const double y2 = trajectory[last_idx].position.Y() * 1000;
            const double z2 = trajectory[last_idx].position.Z() * 1000;
            const TVector3 final_direction(x2 - x1, y2 - y1, z2 - z1);

            // Extract final positions for momentum variations
            const TVector3 final_position_plus(
                trajectory_plus[trajectory_plus.size() - 1].position.X() * 1000,
                trajectory_plus[trajectory_plus.size() - 1].position.Y() * 1000,
                trajectory_plus[trajectory_plus.size() - 1].position.Z() * 1000);

            const TVector3 final_position_minus(
                trajectory_minus[trajectory_minus.size() - 1].position.X() * 1000,
                trajectory_minus[trajectory_minus.size() - 1].position.Y() * 1000,
                trajectory_minus[trajectory_minus.size() - 1].position.Z() * 1000);

            //----------------------------------------------------------------------
            // Match propagated track to detector hit - use precomputed hit mapping
            //----------------------------------------------------------------------
            // Initialize hit matching variables
            int segment_match = -1;          // Segment index of matched hit
            int bar_match = -1;              // Bar index of matched hit
            Point point_match(0, 0, 0);      // Position of matched hit
            float time_match = 0;            // Time of matched hit
            int bitID_match = 0;             // BitID of matched hit
            TVector3 mom_vec_match(0, 0, 0); // Momentum vector at matched hit

            // Find the earliest hit in MS with matching particle key using our mapping
            const auto &particleHits = keyToHitIndices[key[i]];
            for (size_t idx : particleHits)
            {
                // We already filtered for z-range and support structure when building the map

                // If this is first match or earlier than previous match
                if (segment_match == -1 || ms_time[idx] < time_match)
                {
                    segment_match = ms_bitID[idx] >> 18 & 0xF;
                    bar_match = ms_bitID[idx] >> 22 & 0x3F;
                    TVector3 position = getBarPosition(ms_bitID[idx] >> 8 & 0x7, ms_bitID[idx] >> 11 & 0xF, ms_bitID[idx] >> 15 & 0x7, segment_match, bar_match);
                    // convert position to a Point and store in point_match
                    point_match = Point(position.X(), position.Y(), position.Z());
                    // Create position object for this hit
                    // Point ms_pos(ms_vx[idx], ms_vy[idx], ms_vz[idx]);
                    // point_match = ms_pos;
                    time_match = ms_time[idx];
                    bitID_match = ms_bitID[idx];
                    mom_vec_match = TVector3(ms_px[idx], 0, ms_pz[idx]);
                    mom_vec_match = mom_vec_match.Unit();
                }
            }

            // Record all hits for clustering analysis - use a hash set for faster lookup
            std::unordered_set<int> bitID_set;
            for (size_t j = 0; j < ms_vz.GetSize(); ++j)
            {
                // Apply z-range cut for MS hits
                if (ms_vz[j] < 3000 || ms_vz[j] > 7700)
                    continue;

                // Skip support structure hits
                if ((ms_bitID[j] >> 28 & 0x3) || (ms_bitID[j] >> 30 & 0x3))
                    continue;

                // Use the set for faster lookup
                if (bitID_set.insert(ms_bitID[j]).second)
                {
                    // New hit - add to cluster collections
                    clusterizedHits_bitID.push_back(ms_bitID[j]);
                    clusterizedHits_time.push_back(ms_time[j]);
                    clusterizedHits_id.push_back(ms_id[j]);
                }
                else
                {
                    bool sameTDC = false;
                    // Existing hit - find and update if this hit is earlier
                    for (size_t k = 0; k < clusterizedHits_bitID.size(); ++k)
                    {
                        if (clusterizedHits_bitID[k] == ms_bitID[j])
                        {
                            // compare TDC with all other hits and require at least 2ns difference from all other hits
                            double tdc_currhit = fmod(ms_time[j] - 16, 25);
                            double tdc_existinghit = fmod(clusterizedHits_time[k] - 16, 25);
                            if (fabs(tdc_currhit - tdc_existinghit) < 2)
                            {
                                sameTDC = true;
                            }
                        }
                    }
                    if (!sameTDC)
                    {
                        // Add hit to cluster collections
                        clusterizedHits_bitID.push_back(ms_bitID[j]);
                        clusterizedHits_time.push_back(ms_time[j]);
                        clusterizedHits_id.push_back(ms_id[j]);
                    }
                }
            }

            // Skip particles with no matching MS hit
            if (segment_match == -1)
            {
                if (verboseoutput)
                {
                    std::cout << "No matching segment found" << std::endl;
                }
                continue;
            }

            // fill hExtrapolatedTimeVsTrue
            hExtrapolatedTimeVsTrue->Fill(final_time - time_match);

            // Continue with the rest of your processing...
            // [The remaining code after this point stays largely the same]
            if (verboseoutput)
            {
                // Print details of matched hit
                std::cout << "\tMatching segment found: " << segment_match
                          << " layer: " << ((bitID_match >> 15) & 0x7)
                          << " bar: " << bar_match
                          << " module: " << ((bitID_match >> 11) & 0xF)
                          << " station: " << ((bitID_match >> 8) & 0x7)
                          << " time: " << time_match << std::endl;
            }
            if (false)
            {
                std::cout << "true matched time: " << time_match << std::endl;
            }

            //----------------------------------------------------------------------
            // Fill tracking match quality histograms
            //----------------------------------------------------------------------
            //only fill the histograms below for station 0
            // if (((bitID_match >> 8) & 0x7) == 1)
            // {
            hMatchdYdZ->Fill(point_match.y - final_position.Y(),
                             point_match.z - final_position.Z());
            hMatchdXdZ->Fill(point_match.x - final_position.X(),
                             point_match.z - final_position.Z());
            hMatchdX->Fill(point_match.x - final_position.X());
            hMatchdY->Fill(point_match.y - final_position.Y());
            hMatchdZ->Fill(point_match.z - final_position.Z());
            hMatchdXModule[(bitID_match >> 11) & 0xF]->Fill(point_match.x - final_position.X());
            hMatchdYModule[(bitID_match >> 11) & 0xF]->Fill(point_match.y - final_position.Y());
            hMatchdZModule[(bitID_match >> 11) & 0xF]->Fill(point_match.z - final_position.Z());
            hMatchdXStation[(bitID_match >> 8) & 0x7]->Fill(point_match.x - final_position.X());
            hMatchdYStation[(bitID_match >> 8) & 0x7]->Fill(point_match.y - final_position.Y());
            hMatchdZStation[(bitID_match >> 8) & 0x7]->Fill(point_match.z - final_position.Z());
            // }

            // Fill position histograms
            hTrueHitXZ->Fill(point_match.x, point_match.z);
            hExtrapolatedXZ->Fill(final_position.X(), final_position.Z());
            hTrueHitYZ->Fill(point_match.y, point_match.z);
            hExtrapolatedYZ->Fill(final_position.Y(), final_position.Z());

            // Fill momentum smearing effect histograms
            double dxExtrapolSmeared = final_position_plus.X() - final_position_minus.X();
            double dzExtrapolSmeared = final_position_plus.Z() - final_position_minus.Z();
            hExtrapolationMomentumSmearingdXdZ->Fill(dxExtrapolSmeared, dzExtrapolSmeared);

            double dyExtrapolSmeared = final_position_plus.Y() - final_position_minus.Y();
            hExtrapolationMomentumSmearingdYdZ->Fill(dyExtrapolSmeared, dzExtrapolSmeared);

            //----------------------------------------------------------------------
            // Find nearby hits in same detector region
            //----------------------------------------------------------------------
            int stationMatch = (bitID_match >> 8) & 0x7; // Station ID
            int moduleMatch = (bitID_match >> 11) & 0xF; // Module ID

            // Containers for nearby hits
            std::vector<int> close_by_hits_layer;  // Layer indices
            std::vector<int> close_by_hits_bitID;  // Segment bit IDs
            std::vector<float> close_by_hits_time; // Hit times
            std::vector<bool> close_by_hits_true;  // Truth match flags

            // Loop over clustered hits to find those in same module and time window
            for (size_t j = 0; j < clusterizedHits_bitID.size(); ++j)
            {
                // Check if hit is in the same station and module
                if (((clusterizedHits_bitID[j] >> 8) & 0x7) != stationMatch ||
                    ((clusterizedHits_bitID[j] >> 11) & 0xF) != moduleMatch)
                {
                    continue;
                }

                // NOTE Check if hit is in allowed segment range and time window
                if (((clusterizedHits_bitID[j] >> 18 & 0xF) > (segment_match - maxSegmentDifferenceMinus)) &&
                    ((clusterizedHits_bitID[j] >> 18 & 0xF) < (segment_match + maxSegmentDifferencePlus)) &&
                    // NOTE use the extrapolation time instead of the true hit time here now
                    (clusterizedHits_time[j] > (final_time - maxTimeDifference)) &&
                    (clusterizedHits_time[j] < (final_time + maxTimeDifference)))
                // (clusterizedHits_time[j] > (time_match - maxTimeDifference)) &&
                // (clusterizedHits_time[j] < (time_match + maxTimeDifference)))
                {

                    // NOTE Check if hit is within 50 bars of the matched bar
                    if (abs((clusterizedHits_bitID[j] >> 22 & 0x3F) - bar_match) < maxBarDifference)
                    {
                        // Add hit to close-by collection
                        close_by_hits_layer.push_back(clusterizedHits_bitID[j] >> 15 & 0x7);
                        close_by_hits_bitID.push_back(clusterizedHits_bitID[j]);
                        close_by_hits_time.push_back(clusterizedHits_time[j]);

                        // Mark if this is a truth-matched hit
                        close_by_hits_true.push_back(clusterizedHits_id[j] == key[i]);
                    }
                }
            }

            // Skip tracks with too few nearby hits
            if (close_by_hits_layer.size() < 3)
            {
                if (verboseoutput)
                {
                    std::cout << "\t\tNot enough close by hits found: "
                              << close_by_hits_layer.size() << std::endl;
                }
                continue;
            }
            else
            {
                if (verboseoutput)
                    std::cout << "\t\tclose by hits found: " << close_by_hits_layer.size() << std::endl;
            }
            // check how many layers are hit
            std::vector<int> layers_hit(4, 0);
            for (size_t j = 0; j < close_by_hits_layer.size(); ++j)
            {
                layers_hit[close_by_hits_layer[j]]++;
            }
            int layers_hit_count = 0;
            for (size_t j = 0; j < layers_hit.size(); ++j)
            {
                if (layers_hit[j] > 0)
                    layers_hit_count++;
            }
            if (layers_hit_count < 3)
            {
                if (verboseoutput)
                    std::cout << "\t\t\tnot enough layers hit " << layers_hit_count << std::endl;
                continue;
            }
            else
            {
                if (verboseoutput)
                    std::cout << "\t\t\tlayers hit: " << layers_hit_count << std::endl;
            }
            //----------------------------------------------------------------------
            // Sort hits by layer for consistent processing
            //----------------------------------------------------------------------
            // Create indices for sorting
            std::vector<size_t> indices(close_by_hits_layer.size());
            std::iota(indices.begin(), indices.end(), 0); // Fill with 0, 1, 2, ...

            // Sort indices based on layer values (ascending order)
            std::sort(indices.begin(), indices.end(), [&close_by_hits_layer](size_t i1, size_t i2)
                      { return close_by_hits_layer[i1] < close_by_hits_layer[i2]; });

            // Create temporary copies of hit arrays
            std::vector<int> sorted_close_by_hits_layer = close_by_hits_layer;
            std::vector<int> sorted_close_by_hits_bitID = close_by_hits_bitID;
            std::vector<float> sorted_close_by_hits_time = close_by_hits_time;
            std::vector<bool> sorted_close_by_hits_true = close_by_hits_true;

            // Rearrange arrays according to sorted order
            for (size_t j = 0; j < indices.size(); ++j)
            {
                close_by_hits_layer[j] = sorted_close_by_hits_layer[indices[j]];
                close_by_hits_bitID[j] = sorted_close_by_hits_bitID[indices[j]];
                close_by_hits_time[j] = sorted_close_by_hits_time[indices[j]];
                close_by_hits_true[j] = sorted_close_by_hits_true[indices[j]];
            }

            //----------------------------------------------------------------------
            // Analyze truth matching and bar differences between layers
            //----------------------------------------------------------------------
            int layer0_true_bar = -1;
            if (verboseoutput)
            {
                // Print sorted hit information for debugging
                for (size_t j = 0; j < close_by_hits_layer.size(); ++j)
                {
                    std::cout << "\t\t\thit " << j
                              << " layer: " << close_by_hits_layer[j]
                              << " bar: " << (close_by_hits_bitID[j] >> 22 & 0x3F)
                              << " segment " << (close_by_hits_bitID[j] >> 18 & 0xF)
                              << " time: " << close_by_hits_time[j]
                              << " true: " << close_by_hits_true[j] << std::endl;
                }
            }

            // Find the true bar position in layer 0 for reference
            for (size_t j = 0; j < close_by_hits_layer.size(); ++j)
            {
                if (close_by_hits_true[j])
                {
                    if (close_by_hits_layer[j] == 0)
                    {
                        // Record the reference bar number in layer 0
                        layer0_true_bar = (close_by_hits_bitID[j] >> 22 & 0x3F);
                    }
                    else if (layer0_true_bar != -1)
                    {
                        // Calculate relative bar differences for truth-matched hits in other layers
                        int current_bar = (close_by_hits_bitID[j] >> 22 & 0x3F);
                        int bar_difference = current_bar - layer0_true_bar;

                        // Fill histograms with bar differences
                        hRelativeBarDifferenceToFirstLayer->Fill(close_by_hits_layer[j], bar_difference);

                        // Handle stations with mirror symmetry (stations 0 and 3)
                        int station = (close_by_hits_bitID[j] >> 8) & 0x7;
                        if (station == 0 || station == 3)
                        {
                            hRelativeBarDifferenceToFirstLayer->Fill(
                                close_by_hits_layer[j], layer0_true_bar - current_bar);
                        }
                        else
                        {
                            hRelativeBarDifferenceToFirstLayerIndiv[close_by_hits_layer[j] - 1]->Fill(
                                bar_difference);
                        }
                    }
                }
            }

            //----------------------------------------------------------------------
            // Group hits across layers that may form tracks
            //----------------------------------------------------------------------
            std::vector<std::vector<int>> tracklet_hits_layers;   // Layer indices for each group
            std::vector<std::vector<int>> tracklet_hits_bitIDs;   // BitIDs for each group
            std::vector<std::vector<float>> tracklet_hits_times;  // Time values for each group
            std::vector<std::vector<bool>> tracklet_hits_trues;   // Truth-matching flags for each group
            std::vector<std::vector<bool>> tracklet_hits_removed; // Flags for removing hits during filtering

            // Start grouping from hits in layer 0
            for (size_t j = 0; j < close_by_hits_layer.size(); ++j)
            {
                // Skip hits not in layer 0
                if (close_by_hits_layer[j] != 0)
                {
                    continue;
                }

                // Initialize a new hit group starting from this layer 0 hit
                std::vector<int> group_layers = {close_by_hits_layer[j]};
                std::vector<int> group_bitIDs = {close_by_hits_bitID[j]};
                std::vector<float> group_times = {close_by_hits_time[j]};
                std::vector<bool> group_trues = {close_by_hits_true[j]};
                std::vector<bool> group_removed = {false}; // Layer 0 hit is never removed
                // For each subsequent layer, find matching hits within time window
                for (int layer = 1; layer < 4; ++layer)
                {
                    for (size_t k = 0; k < close_by_hits_layer.size(); ++k)
                    {
                        // Skip if not in current target layer
                        if (close_by_hits_layer[k] != layer)
                        {
                            continue;
                        }

                        // NOTE Check time coincidence window (1ns)
                        if (abs(close_by_hits_time[j] - close_by_hits_time[k]) < maxTimeDifference / 2)
                        {
                            // Extract bar indices
                            int bar_layer0 = close_by_hits_bitID[j] >> 22 & 0x3F;
                            int bar_current = close_by_hits_bitID[k] >> 22 & 0x3F;
                            int bar_diff = abs(bar_layer0 - bar_current);

                            // NOTE Apply layer-specific constraints on bar difference
                            // NOTE These reflect physical constraints on track curvature
                            // values have been determined from true studies which set the limits
                            if (layer == 1 && bar_diff > 2)
                                continue; // Layer 1: max diff = 2
                            if (layer == 2 && bar_diff > 6)
                                continue; // Layer 2: max diff = 6
                            if (layer == 3 && bar_diff > 7)
                                continue; // Layer 3: max diff = 7

                            // Add hit to current group
                            group_layers.push_back(close_by_hits_layer[k]);
                            group_bitIDs.push_back(close_by_hits_bitID[k]);
                            group_times.push_back(close_by_hits_time[k]);
                            group_trues.push_back(close_by_hits_true[k]);
                            group_removed.push_back(true); // Initially mark for removal
                        }
                    }
                }

                // Store completed group
                tracklet_hits_layers.push_back(group_layers);
                tracklet_hits_bitIDs.push_back(group_bitIDs);
                tracklet_hits_times.push_back(group_times);
                tracklet_hits_trues.push_back(group_trues);
                tracklet_hits_removed.push_back(group_removed);
            }

            // Add this after creating all tracklets but before segment combination matching (around line 1009)

            // Create vectors to track which tracklets should be merged
            std::vector<int> tracklets_to_merge;
            std::vector<int> merge_targets;

            // Look for layer 0 hits that are adjacent to other tracklets' layer 0 hits
            for (size_t i = 0; i < tracklet_hits_layers.size(); i++)
            {
                // Skip if this tracklet has already been marked for merging
                if (std::find(tracklets_to_merge.begin(), tracklets_to_merge.end(), i) != tracklets_to_merge.end())
                    continue;

                // Extract the layer 0 bar information for this tracklet
                int current_bar = -1;
                for (size_t j = 0; j < tracklet_hits_layers[i].size(); j++)
                {
                    if (tracklet_hits_layers[i][j] == 0)
                    {
                        current_bar = (tracklet_hits_bitIDs[i][j] >> 22) & 0x3F;
                        break;
                    }
                }

                if (current_bar == -1)
                    continue; // Skip if no layer 0 hit found

                // Check if this tracklet's layer 0 hit is adjacent to another tracklet's layer 0 hit
                for (size_t k = i + 1; k < tracklet_hits_layers.size(); k++)
                {
                    // Skip if this candidate has already been marked for merging
                    if (std::find(tracklets_to_merge.begin(), tracklets_to_merge.end(), k) != tracklets_to_merge.end())
                        continue;

                    int other_bar = -1;
                    for (size_t j = 0; j < tracklet_hits_layers[k].size(); j++)
                    {
                        if (tracklet_hits_layers[k][j] == 0)
                        {
                            other_bar = (tracklet_hits_bitIDs[k][j] >> 22) & 0x3F;
                            break;
                        }
                    }

                    if (other_bar == -1)
                        continue; // Skip if no layer 0 hit found

                    // Check if bars are adjacent
                    if (abs(current_bar - other_bar) == 1)
                    {
                        if (verboseoutput)
                        {
                            std::cout << "\t\t\tFound adjacent layer 0 hits in tracklets " << i
                                      << " (bar " << current_bar << ") and " << k
                                      << " (bar " << other_bar << ")" << std::endl;
                        }

                        // Choose which tracklet should absorb the other one
                        // For example, keep the one with more hits or the one with a true hit
                        bool keep_i = tracklet_hits_layers[i].size() >= tracklet_hits_layers[k].size();

                        // // Check if either tracklet has truth-matched hits
                        // bool i_has_truth = false, k_has_truth = false;
                        // for (size_t j = 0; j < tracklet_hits_trues[i].size(); j++)
                        // {
                        //     if (tracklet_hits_trues[i][j])
                        //     {
                        //         i_has_truth = true;
                        //         break;
                        //     }
                        // }
                        // for (size_t j = 0; j < tracklet_hits_trues[k].size(); j++)
                        // {
                        //     if (tracklet_hits_trues[k][j])
                        //     {
                        //         k_has_truth = true;
                        //         break;
                        //     }
                        // }

                        // // Prefer keeping tracklet with truth-matched hits
                        // if (i_has_truth && !k_has_truth)
                        //     keep_i = true;
                        // else if (!i_has_truth && k_has_truth)
                        //     keep_i = false;

                        if (keep_i)
                        {
                            tracklets_to_merge.push_back(k);
                            merge_targets.push_back(i);
                        }
                        else
                        {
                            tracklets_to_merge.push_back(i);
                            merge_targets.push_back(k);
                        }
                        break; // Stop checking this tracklet against others
                    }
                }
            }

            // Merge tracklets
            for (size_t i = 0; i < tracklets_to_merge.size(); i++)
            {
                int source = tracklets_to_merge[i];
                int target = merge_targets[i];

                if (verboseoutput)
                {
                    std::cout << "\t\t\tMerging tracklet " << source << " into tracklet " << target << std::endl;
                }

                // Add all hits from source tracklet to target tracklet
                for (size_t j = 0; j < tracklet_hits_layers[source].size(); j++)
                {
                    int layer = tracklet_hits_layers[source][j];
                    int bitID = tracklet_hits_bitIDs[source][j];
                    float time = tracklet_hits_times[source][j];
                    bool is_true = tracklet_hits_trues[source][j];
                    bool is_removed = tracklet_hits_removed[source][j];

                    // Skip layer 0 hits since we've already decided which tracklet's layer 0 hit to keep
                    // This prevents having multiple layer 0 hits in the merged tracklet
                    if (layer == 0)
                    {
                        if (verboseoutput)
                        {
                            std::cout << "\t\t\t\tSkipping layer 0 hit from source tracklet " << source << std::endl;
                        }
                        continue;
                    }

                    // Check if this hit's layer is already in the target tracklet
                    bool duplicate = false;
                    for (size_t k = 0; k < tracklet_hits_layers[target].size(); k++)
                    {
                        if (tracklet_hits_layers[target][k] == layer &&
                            (tracklet_hits_bitIDs[target][k] >> 22 & 0x3F) == (bitID >> 22 & 0x3F))
                        {
                            duplicate = true;

                            // If we have a duplicate, keep the one with better properties (e.g., truth match)
                            if (is_true && !tracklet_hits_trues[target][k])
                            {
                                tracklet_hits_trues[target][k] = true;
                            }
                            if (!is_removed && tracklet_hits_removed[target][k])
                            {
                                tracklet_hits_removed[target][k] = false;
                            }
                            break;
                        }
                    }

                    // Add hit if not a duplicate
                    if (!duplicate)
                    {
                        tracklet_hits_layers[target].push_back(layer);
                        tracklet_hits_bitIDs[target].push_back(bitID);
                        tracklet_hits_times[target].push_back(time);
                        tracklet_hits_trues[target].push_back(is_true);
                        tracklet_hits_removed[target].push_back(is_removed);
                    }
                }

                // Mark source tracklet for removal by clearing its hits
                tracklet_hits_layers[source].clear();
                tracklet_hits_bitIDs[source].clear();
                tracklet_hits_times[source].clear();
                tracklet_hits_trues[source].clear();
                tracklet_hits_removed[source].clear();
            }

            //----------------------------------------------------------------------
            // Match hit tracklets to valid segment combinations
            //----------------------------------------------------------------------
            // Vector to store matched (group_index, combination_index) pairs
            std::vector<std::pair<int, int>> matched_combinations;

            // Compare each group against valid segment combinations
            for (size_t group_idx = 0; group_idx < tracklet_hits_layers.size(); ++group_idx)
            {
                // Skip tracklets with insufficient hits
                if (tracklet_hits_layers[group_idx].size() < 3)
                {
                    continue;
                }

                // Get segment ID from layer 0 hit in this group
                int first_segment = tracklet_hits_bitIDs[group_idx][0] >> 18 & 0xF;

                // Compare against each valid segment combination
                for (size_t comb_idx = 0; comb_idx < vec_segment_combinations.size(); ++comb_idx)
                {
                    // Skip combinations that don't match layer 0 segment
                    if (vec_segment_combinations[comb_idx][0] != first_segment)
                    {
                        continue;
                    }

                    // Extract expected segment IDs for other layers
                    int expected_segment_layer1 = vec_segment_combinations[comb_idx][1];
                    int expected_segment_layer2 = vec_segment_combinations[comb_idx][2];
                    int expected_segment_layer3 = vec_segment_combinations[comb_idx][3];

                    // Track if each layer has a hit and if it matches the combination
                    bool has_hit[3] = {false, false, false};
                    bool matches_combination[3] = {false, false, false};

                    // Check hits in the group against pattern
                    for (size_t hit_idx = 0; hit_idx < tracklet_hits_layers[group_idx].size(); ++hit_idx)
                    {
                        int layer = tracklet_hits_layers[group_idx][hit_idx];
                        int segment = tracklet_hits_bitIDs[group_idx][hit_idx] >> 18 & 0xF;

                        if (layer == 1)
                        {
                            has_hit[0] = true;
                            if (segment == expected_segment_layer1)
                            {
                                matches_combination[0] = true;
                            }
                        }
                        else if (layer == 2)
                        {
                            has_hit[1] = true;
                            if (segment == expected_segment_layer2)
                            {
                                matches_combination[1] = true;
                            }
                        }
                        else if (layer == 3)
                        {
                            has_hit[2] = true;
                            if (segment == expected_segment_layer3)
                            {
                                matches_combination[2] = true;
                            }
                        }
                    }

                    // Skip if any layer has a hit that doesn't match the expected pattern
                    if ((has_hit[0] && !matches_combination[0]) ||
                        (has_hit[1] && !matches_combination[1]) ||
                        (has_hit[2] && !matches_combination[2]))
                    {
                        continue;
                    }

                    // We found a valid match
                    if (verboseoutput)
                    {
                        std::cout << "\t\t\tgroup " << group_idx << " matches combination " << comb_idx
                                  << " with segments " << first_segment << " "
                                  << expected_segment_layer1 << " "
                                  << expected_segment_layer2 << " "
                                  << expected_segment_layer3 << std::endl;
                    }

                    matched_combinations.push_back(std::make_pair(group_idx, comb_idx));

                    // Keep hits that match the combination pattern (unmark them for removal)
                    for (size_t hit_idx = 0; hit_idx < tracklet_hits_layers[group_idx].size(); ++hit_idx)
                    {
                        int layer = tracklet_hits_layers[group_idx][hit_idx];
                        int segment = tracklet_hits_bitIDs[group_idx][hit_idx] >> 18 & 0xF;

                        if ((layer == 1 && segment == expected_segment_layer1) ||
                            (layer == 2 && segment == expected_segment_layer2) ||
                            (layer == 3 && segment == expected_segment_layer3))
                        {
                            // NOTE this boolean is important, it tracks that the hit is part of a valid combination of segments
                            tracklet_hits_removed[group_idx][hit_idx] = false;
                        }
                    }
                }
            }
            //----------------------------------------------------------------------
            // Count and analyze reconstructed track candidates
            //----------------------------------------------------------------------
            // Count tracklets with sufficient valid hits after filtering
            int validGroupCount = 0;
            int validGroupCountTrue = 0;
            for (size_t groupIdx = 0; groupIdx < tracklet_hits_layers.size(); ++groupIdx)
            {
                int validHitsInGroup = 0;
                int validHitsInGroupTrue = 0;

                if (verboseoutput)
                {
                    std::cout << "\t\t\tgroup " << groupIdx;
                    if (tracklet_hits_layers[groupIdx].size() < 3)
                    {
                        std::cout << " not enough hits in group";
                    }
                    std::cout << std::endl;
                }

                // Count hits that passed the filtering criteria
                for (size_t hitIdx = 0; hitIdx < tracklet_hits_layers[groupIdx].size(); ++hitIdx)
                {
                    if (verboseoutput)
                    {
                        std::cout << "\t\t\t\thit " << hitIdx
                                  << " layer: " << tracklet_hits_layers[groupIdx][hitIdx]
                                  << " bar: " << (tracklet_hits_bitIDs[groupIdx][hitIdx] >> 22 & 0x3F)
                                  << " segment " << (tracklet_hits_bitIDs[groupIdx][hitIdx] >> 18 & 0xF)
                                  << " time: " << tracklet_hits_times[groupIdx][hitIdx]
                                  << " true: " << tracklet_hits_trues[groupIdx][hitIdx];

                        if (tracklet_hits_removed[groupIdx][hitIdx])
                        {
                            std::cout << " removed";
                        }
                        else
                        {
                            std::cout << " kept";
                        }
                        std::cout << std::endl;
                    }

                    if (!tracklet_hits_removed[groupIdx][hitIdx])
                    {
                        validHitsInGroup++;
                        if (tracklet_hits_trues[groupIdx][hitIdx])
                        {
                            validHitsInGroupTrue++;
                        }
                    }
                }

                // Count tracklets with enough valid hits for track reconstruction
                if (validHitsInGroup >= 3)
                {
                    validGroupCount++;
                    if (validHitsInGroupTrue >= 1)
                    {
                        validGroupCountTrue++;
                    }
                }
            }

            // fill hNTracklets with the number of tracklets that follow the correct segment pattern, have at least 3 hits and are in the correct time window
            hNTracklets->Fill(validGroupCount);
            hNTrackletsTrue->Fill(validGroupCountTrue);

            //----------------------------------------------------------------------
            // Fit tracks to hit tracklets and analyze track quality
            //----------------------------------------------------------------------
            int successfulFitCount = 0;                    // Counter for tracklets with good tracklet fits
            int successfulFitCountTimeCut = 0;             // Counter for tracklets with good tracklet fits
            int successfulFitCountSlopeAndTimeCut = 0;     // Counter for tracklets with good tracklet fits
            int successfulFitCountTrue = 0;                // Counter for tracklets with good tracklet fits
            int successfulFitCountTimeCutTrue = 0;         // Counter for tracklets with good tracklet fits
            int successfulFitCountSlopeAndTimeCutTrue = 0; // Counter for tracklets with good tracklet fits

            // Process each hit group and perform track fitting
            for (size_t groupIdx = 0; groupIdx < tracklet_hits_layers.size(); ++groupIdx)
            {
                // Skip tracklets with too few hits
                if (tracklet_hits_layers[groupIdx].size() < 3)
                {
                    continue;
                }

                // Initialize track vector
                TVector3 tracklet_vec(0, 0, 0);

                // Create graph for fitting hit positions
                TGraphErrors *trackGraph = new TGraphErrors();
                bool track_contains_true_hit = false;
                double track_time = 10000000;

                // Create additional graphs for truth-matched hits if we haven't hit the limit
                if (nSlopeFits < 29)
                {
                    gSlopeFits_true[nSlopeFits] = new TGraphErrors();
                }

                // Fill the graph with valid hits
                int validHitCount = 0;
                for (size_t hitIdx = 0; hitIdx < tracklet_hits_layers[groupIdx].size(); ++hitIdx)
                {
                    // Skip hits that were filtered out
                    if (tracklet_hits_removed[groupIdx][hitIdx])
                    {
                        continue;
                    }

                    // Extract layer and bar information
                    int layer = tracklet_hits_layers[groupIdx][hitIdx];
                    int bar = tracklet_hits_bitIDs[groupIdx][hitIdx] >> 22 & 0x3F;
                    double time = tracklet_hits_times[groupIdx][hitIdx];
                    if (time < track_time)
                    {
                        track_time = time;
                    }

                    // Add point to main graph
                    // X = layer position offset, Y = bar number * scale factor
                    trackGraph->SetPoint(validHitCount, layer_offsets[layer], 5 * bar);

                    // Set appropriate error bars based on layer
                    double xError = (layer == 0) ? 1.5 : 2.5;
                    double yError = (layer == 0) ? 1.5 : 3.5;
                    trackGraph->SetPointError(validHitCount, xError, yError);

                    // For truth-matched hits, also add to the truth graph with smaller errors
                    if (tracklet_hits_trues[groupIdx][hitIdx])
                    {
                        track_contains_true_hit = true;
                        if (nSlopeFits < 29)
                        {
                            gSlopeFits_true[nSlopeFits]->SetPoint(validHitCount, layer_offsets[layer], 5 * bar);
                            gSlopeFits_true[nSlopeFits]->SetPointError(validHitCount, 0.5, 0.5);
                        }
                    }

                    validHitCount++;
                }

                // Skip fitting if not enough valid hits remain
                if (validHitCount < 3)
                {
                    delete trackGraph;
                    if (verboseoutput)
                    {
                        std::cout << "\t\t\tgroup " << groupIdx
                                  << " not enough hits in group to fit" << std::endl;
                    }
                    continue;
                }

                //----------------------------------------------------------------------
                // Perform linear track fitting
                //----------------------------------------------------------------------
                // Create fit function (y = p0 + p1*x)
                TF1 *fitFunction = new TF1("fitFunction", "pol1", 0, 7);

                // Set parameter constraints for stability
                fitFunction->SetParLimits(0, 0, 500); // Intercept limits
                fitFunction->SetParLimits(1, -1, 1);  // Slope limits

                // Perform initial fit
                trackGraph->Fit(fitFunction, "NQ"); // 'Q' = quiet mode, 'N' = no drawing

                // Save the original fit before outlier removal
                if (nSlopeFits < 29)
                {
                    gSlopeFits_orig[nSlopeFits] = new TGraphErrors(*trackGraph);
                    funcSlopeFits_orig[nSlopeFits] = new TF1(*fitFunction);
                }

                //----------------------------------------------------------------------
                // Improve fit by iteratively removing outliers
                //----------------------------------------------------------------------
                for (int layer = 1; layer < 4; ++layer)
                {
                    // Calculate expected position at this layer from fit
                    double expectedPosition = fitFunction->Eval(layer_offsets[layer]);

                    // Collect distances of points from fit line
                    std::vector<double> residuals;
                    std::vector<int> pointIndices;

                    // Find all points in this layer
                    for (int pointIdx = 0; pointIdx < trackGraph->GetN(); ++pointIdx)
                    {
                        if (abs(trackGraph->GetX()[pointIdx] - layer_offsets[layer]) < 0.1)
                        {
                            // Calculate distance from fit line
                            double x, y;
                            trackGraph->GetPoint(pointIdx, x, y);
                            double residual = abs(y - expectedPosition);

                            residuals.push_back(residual);
                            pointIndices.push_back(pointIdx);
                        }
                    }

                    // If multiple hits in this layer, remove the one furthest from fit line
                    if (residuals.size() > 1)
                    {
                        auto maxResidual = std::max_element(residuals.begin(), residuals.end());
                        int worstPointIdx = pointIndices[std::distance(residuals.begin(), maxResidual)];
                        trackGraph->RemovePoint(worstPointIdx);
                    }

                    // Refit with the updated set of points
                    fitFunction->SetParLimits(0, 0, 500);
                    fitFunction->SetParLimits(1, -1, 1);
                    trackGraph->Fit(fitFunction, "NQ");
                }

                // Save the improved fit after outlier removal
                if (nSlopeFits < 29)
                {
                    gSlopeFits[nSlopeFits] = new TGraphErrors(*trackGraph);
                    funcSlopeFits[nSlopeFits] = new TF1(*fitFunction);
                    nSlopeFits++;
                }

                // Extract fit parameters
                double intercept = fitFunction->GetParameter(0);
                double slope = fitFunction->GetParameter(1);

                if (verboseoutput)
                {
                    std::cout << "\t\t\tgroup " << groupIdx << " fit: "
                              << intercept << " + " << slope << " * t" << std::endl;
                }

                //----------------------------------------------------------------------
                // Calculate tracklet direction vector and compare with truth
                //----------------------------------------------------------------------
                // Create direction vector based on station orientation
                // Note: different stations have different orientations requiring sign flips
                tracklet_vec.SetXYZ(stationMatch == 0 || stationMatch == 2 ? -1 : 1, 0, slope);
                tracklet_vec = tracklet_vec.Unit();

                // Apply standard rotation for the station's global position
                if (stationMatch == 0 || stationMatch == 2)
                {
                    tracklet_vec.RotateY(-0.44); // -25.2 degrees
                }
                else
                {
                    tracklet_vec.RotateY(0.44); // +25.2 degrees
                }

                // Apply module-specific rotation
                if (stationMatch == 0 || stationMatch == 2)
                {
                    tracklet_vec.RotateY(vec_module_rotations[moduleMatch]);
                }
                else
                {
                    tracklet_vec.RotateY(-vec_module_rotations[moduleMatch]);
                }

                // Calculate angle between reconstructed tracklet and true momentum vector
                double angle = tracklet_vec.Angle(mom_vec_match);
                hAngleDiffModules->Fill(moduleMatch, angle);

                // calculate the angle between the tracklet and the final direction
                double angle2 = tracklet_vec.Angle(final_direction);
                hAngleDiffFinalModules->Fill(moduleMatch, angle2);

                // Calculate the time difference between the tracklet and the true time
                double time_diff = final_time - track_time;
                // std::cout << "time diff " << time_diff << " track_time: " << track_time << " final_time: " << final_time << " time_match: " << time_match << std::endl;
                hExtrapolatedTimeVsTracklet->Fill(time_diff);

                // Count tracks with good angular agreement (< 0.2 rad ≈ 11.5°)
                if (angle2 < maxAngleDifference)
                {
                    successfulFitCount++;
                    if (track_contains_true_hit)
                    {
                        successfulFitCountTrue++;
                    }
                }
                if (abs(time_diff) < maxTimeDifference / 2)
                {
                    successfulFitCountTimeCut++;
                    if (track_contains_true_hit)
                    {
                        successfulFitCountTimeCutTrue++;
                    }
                }
                if (angle2 < maxAngleDifference && abs(time_diff) < maxTimeDifference / 2)
                {
                    successfulFitCountSlopeAndTimeCut++;
                    if (track_contains_true_hit)
                    {
                        successfulFitCountSlopeAndTimeCutTrue++;
                    }
                    hMatchdYdZCut->Fill(point_match.y - final_position.Y(),
                                        point_match.z - final_position.Z());
                    hMatchdXdZCut->Fill(point_match.x - final_position.X(),
                                        point_match.z - final_position.Z());
                    hMatchdXCut->Fill(point_match.x - final_position.X());
                    hMatchdYCut->Fill(point_match.y - final_position.Y());
                    hMatchdZCut->Fill(point_match.z - final_position.Z());
                }

                // Output detailed comparison if in verbose mode
                if (verboseoutput)
                {
                    std::cout << "\t\t\tgroup " << groupIdx << " angle2: " << angle2 << std::endl;
                    std::cout << "\t\t\tgroup " << groupIdx << " momentum: "
                              << mom_vec_match.X() << " " << mom_vec_match.Y() << " " << mom_vec_match.Z() << std::endl;
                    std::cout << "\t\t\tgroup " << groupIdx << " tracklet: "
                              << tracklet_vec.X() << " " << tracklet_vec.Y() << " " << tracklet_vec.Z() << std::endl;
                }

                // Clean up
                delete trackGraph;
                delete fitFunction;

                // loop over ms hits and find corresponding hit to layer 0 hit of tracklet. get position from this hit and compare to extrapolated track position
                if (!isPbPb)
                {
                    // Find layer 0 hit of this tracklet and compare its actual position to extrapolated track position
                    for (size_t hitIdx = 0; hitIdx < tracklet_hits_layers[groupIdx].size(); ++hitIdx)
                    {
                        // Find the layer 0 hit in this tracklet
                        if (tracklet_hits_layers[groupIdx][hitIdx] == 0 && !tracklet_hits_removed[groupIdx][hitIdx])
                        {
                            int layer0_bitID = tracklet_hits_bitIDs[groupIdx][hitIdx];
                            float layer0_time = tracklet_hits_times[groupIdx][hitIdx];

                            if (verboseoutput)
                            {
                                std::cout << "\t\t\tComparing layer 0 hit position with extrapolation for tracklet " << groupIdx
                                          << " bitID: " << layer0_bitID << std::endl;
                            }

                            // Loop through MS hits to find the actual hit that corresponds to this layer 0 hit
                            for (size_t j = 0; j < ms_vz.GetSize(); ++j)
                            {
                                // Match by bitID (which is unique for each hit location)
                                if (ms_bitID[j] == layer0_bitID)
                                {
                                    // Get the actual position of the hit
                                    Point actual_hit_pos(ms_vx[j], ms_vy[j], ms_vz[j]);

                                    // Calculate position differences
                                    // double dx = actual_hit_pos.x - final_position.X();
                                    // double dy = actual_hit_pos.y - final_position.Y();
                                    // double dz = actual_hit_pos.z - final_position.Z();
                                    double dx = actual_hit_pos.x - point_match.x;
                                    double dy = actual_hit_pos.y - point_match.y;
                                    double dz = actual_hit_pos.z - point_match.z;

                                    // Fill histograms for these differences (create these histograms at the beginning of the file)
                                    hMatchdXExtrapolToTracklet->Fill(dx);
                                    hMatchdYExtrapolToTracklet->Fill(dy);
                                    hMatchdZExtrapolToTracklet->Fill(dz);

                                    // If this tracklet meets the quality criteria, record in separate histograms
                                    if (angle2 < maxAngleDifference && abs(time_diff) < maxTimeDifference / 2)
                                    {
                                        // You can create additional histograms for quality-selected tracks
                                        // Example:
                                        hMatchdXExtrapolToTrackletQuality->Fill(dx);
                                        hMatchdYExtrapolToTrackletQuality->Fill(dy);
                                        hMatchdZExtrapolToTrackletQuality->Fill(dz);
                                    }

                                    break; // No need to continue searching once we found the hit
                                }
                            }

                            break; // Only need to process one layer 0 hit per tracklet
                        }
                    }
                }
            } // end of tracklet fitting loop

            // Fill histogram with number of successful track fits
            hNTrackletsSlope->Fill(successfulFitCount);
            hNTrackletsTime->Fill(successfulFitCountTimeCut);
            hNTrackletsSlopeAndTimeCut->Fill(successfulFitCountSlopeAndTimeCut);
            // Fill histograms with number of successful track fits with at least one true hit
            hNTrackletsSlopeTrue->Fill(successfulFitCountTrue);
            hNTrackletsTimeTrue->Fill(successfulFitCountTimeCutTrue);
            hNTrackletsSlopeAndTimeCutTrue->Fill(successfulFitCountSlopeAndTimeCutTrue);
        }
    } // End of particle loop and event processing
    std::cout << "Max bar ID: " << maxbartemp << std::endl;

    //----------------------------------------------------------------------
    // Configure visualization style
    //----------------------------------------------------------------------
    gStyle->SetOptStat(0);  // Disable statistics box on plots
    gStyle->SetPadTickX(1); // Enable tick marks on X axis
    gStyle->SetPadTickY(1); // Enable tick marks on Y axis

    //----------------------------------------------------------------------
    // Create and save detector pattern analysis plots
    //----------------------------------------------------------------------
    // Plot hit pattern differences between detector layers
    TCanvas c1("c1", "Bar Position Differences", 800, 600);
    c1.SetLogz(); // Use logarithmic z-scale for better visibility of patterns
    hRelativeBarDifferenceToFirstLayer->GetXaxis()->SetTitle("Layer");
    hRelativeBarDifferenceToFirstLayer->GetYaxis()->SetTitle("Relative bar difference to first layer");
    hRelativeBarDifferenceToFirstLayer->Draw("colz");
    c1.SaveAs(Form("%shRelativeBarDifferenceToFirstLayer.pdf", outputdir.Data()));

    // Create detailed histogram showing position differences for each layer separately
    TCanvas c2("c2", "Layer-by-Layer Bar Differences", 1200, 800);
    c2.Divide(3, 1); // Create 3 panels (one for each layer after layer 0)
    for (int i = 0; i < 3; i++)
    {
        c2.cd(i + 1); // Select panel i+1
        hRelativeBarDifferenceToFirstLayerIndiv[i]->GetXaxis()->SetTitle("Relative bar difference to first layer");
        hRelativeBarDifferenceToFirstLayerIndiv[i]->GetYaxis()->SetTitle("Counts");
        hRelativeBarDifferenceToFirstLayerIndiv[i]->Draw();
    }
    c2.SaveAs(Form("%shRelativeBarDifferenceToFirstLayerIndiv.pdf", outputdir.Data()));

    //----------------------------------------------------------------------
    // Create and save track group analysis plots
    //----------------------------------------------------------------------
    // Plot distribution of track tracklets per particle
    TCanvas c3("c3", "Group Counts", 800, 600);
    hNTracklets->GetXaxis()->SetTitle("Number of tracklets");
    hNTracklets->GetYaxis()->SetTitle("Counts");
    hNTracklets->Draw();
    c3.SaveAs(Form("%shNTracklets.pdf", outputdir.Data()));

    TCanvas c3b("c3b", "Group Counts", 800, 600);
    hNTrackletsTrue->GetXaxis()->SetTitle("Number of tracklets with at least one true hit");
    hNTrackletsTrue->GetYaxis()->SetTitle("Counts");
    hNTrackletsTrue->Draw();
    c3b.SaveAs(Form("%shNTrackletsTrue.pdf", outputdir.Data()));

    // Plot distribution of successful track fits per particle
    TCanvas c3x("c3x", "Group Counts", 800, 600);
    hNTrackletsSlope->GetXaxis()->SetTitle(Form("Number of tracklets with successful tracklet fit and angle < %1.1f rad", maxAngleDifference));
    hNTrackletsSlope->GetYaxis()->SetTitle("Counts");
    hNTrackletsSlope->Draw();
    c3x.SaveAs(Form("%shNTrackletsSlope.pdf", outputdir.Data()));

    // plot with at least one true hit
    TCanvas c3c("c3c", "Group Counts", 800, 600);
    hNTrackletsSlopeTrue->GetXaxis()->SetTitle(Form("Number of tracklets with successful tracklet fit and angle < %1.1f rad and at least one true hit", maxAngleDifference));
    hNTrackletsSlopeTrue->GetYaxis()->SetTitle("Counts");
    hNTrackletsSlopeTrue->Draw();
    c3c.SaveAs(Form("%shNTrackletsSlopeTrue.pdf", outputdir.Data()));

    // Plot distribution of successful track fits per particle
    TCanvas c3y("c3y", "Group Counts", 800, 600);
    hNTrackletsTime->GetXaxis()->SetTitle("Number of tracklets with successful tracklet fit and time difference < 2.5 ns");
    hNTrackletsTime->GetYaxis()->SetTitle("Counts");
    hNTrackletsTime->Draw();
    c3y.SaveAs(Form("%shNTrackletsTime.pdf", outputdir.Data()));

    // plot with at least one true hit
    TCanvas c3d("c3d", "Group Counts", 800, 600);
    hNTrackletsTimeTrue->GetXaxis()->SetTitle("Number of tracklets with successful tracklet fit and time difference < 2.5 ns and at least one true hit");
    hNTrackletsTimeTrue->GetYaxis()->SetTitle("Counts");
    hNTrackletsTimeTrue->Draw();
    c3d.SaveAs(Form("%shNTrackletsTimeTrue.pdf", outputdir.Data()));

    // Plot distribution of successful track fits per particle
    TCanvas c3z("c3z", "Group Counts", 800, 600);
    hNTrackletsSlopeAndTimeCut->GetXaxis()->SetTitle(Form("Number of tracklets with successful tracklet fit and angle < %1.1f rad and time difference < 2.5 ns", maxAngleDifference));
    hNTrackletsSlopeAndTimeCut->GetYaxis()->SetTitle("Counts");
    hNTrackletsSlopeAndTimeCut->Draw();
    c3z.SaveAs(Form("%shNTrackletsSlopeAndTimeCut.pdf", outputdir.Data()));

    // plot with at least one true hit
    TCanvas c3e("c3e", "Group Counts", 800, 600);
    hNTrackletsSlopeAndTimeCutTrue->GetXaxis()->SetTitle(Form("Number of tracklets with successful tracklet fit and angle < %1.1f rad and time difference < 2.5 ns and at least one true hit", maxAngleDifference));
    hNTrackletsSlopeAndTimeCutTrue->GetYaxis()->SetTitle("Counts");
    hNTrackletsSlopeAndTimeCutTrue->Draw();
    c3e.SaveAs(Form("%shNTrackletsSlopeAndTimeCutTrue.pdf", outputdir.Data()));

    // plot all four above on one canvas with a legend
    TCanvas c3a("c3a", "Group Counts", 800, 600);
    hNTracklets->SetLineWidth(3);
    hNTrackletsSlope->SetLineWidth(2);
    hNTrackletsTime->SetLineWidth(2);
    hNTrackletsSlopeAndTimeCut->SetLineWidth(2);
    hNTracklets->SetLineColor(kBlack);
    hNTrackletsSlope->SetLineColor(kRed);
    hNTrackletsTime->SetLineColor(kBlue);
    hNTrackletsSlopeAndTimeCut->SetLineColor(kGreen);
    // set x axis tile
    hNTracklets->GetXaxis()->SetTitle("Number of tracklet candidates for extrapolated track");
    hNTrackletsSlopeAndTimeCut->SetLineStyle(2);
    hNTrackletsTime->SetLineStyle(2);
    hNTrackletsSlope->SetLineStyle(2);
    if (isPbPb)
    {
        hNTrackletsSlopeAndTimeCut->Draw();
        hNTrackletsSlope->Draw("same");
        hNTrackletsTime->Draw("same");
        hNTracklets->Draw("same");
    }
    else
    {
        hNTracklets->Draw();
        hNTrackletsSlope->Draw("same");
        hNTrackletsTime->Draw("same");
        hNTrackletsSlopeAndTimeCut->Draw("same");
        hNTracklets->Draw("same");
    }
    TLegend *leg = new TLegend(0.35, 0.65, 0.55, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hNTracklets, "All tracklets", "l");
    leg->AddEntry(hNTrackletsSlope, Form("#DeltaAngle < %1.1f rad", maxAngleDifference), "l");
    leg->AddEntry(hNTrackletsTime, "#Deltat < 2.5 ns", "l");
    leg->AddEntry(hNTrackletsSlopeAndTimeCut, Form("#DeltaAngle < %1.1f rad and #Deltat < 2.5 ns", maxAngleDifference), "l");
    leg->Draw();
    c3a.SaveAs(Form("%shNTrackletsAll.pdf", outputdir.Data()));
    c3a.SaveAs(Form("%shNTrackletsAll.root", outputdir.Data()));

    // Plot examples of track fits with before/after outlier removal
    TCanvas c4("c4", "Track Fitting Examples", 1800, 1600);
    c4.Divide(6, 5); // Create grid for up to 30 example fits
    for (int i = 0; i < nSlopeFits; i++)
    {
        c4.cd(i + 1); // Select panel for this fit

        // Draw original fits before outlier removal (blue)
        gSlopeFits_orig[i]->SetMarkerColor(kBlue);
        gSlopeFits_orig[i]->SetLineColor(kBlue);
        gSlopeFits_orig[i]->SetMarkerStyle(20);
        gSlopeFits_orig[i]->SetMarkerSize(2);
        gSlopeFits_orig[i]->Draw("ap"); // "a" = axis, "p" = points

        // Draw improved fits after outlier removal (red)
        gSlopeFits[i]->GetXaxis()->SetTitle("Layer offset [mm]");
        gSlopeFits[i]->GetYaxis()->SetTitle("Bar number");
        gSlopeFits[i]->SetMarkerStyle(24);
        gSlopeFits[i]->SetLineColor(kRed);
        gSlopeFits[i]->SetMarkerSize(2);
        gSlopeFits[i]->SetMarkerColor(kRed);
        gSlopeFits[i]->Draw("p,same"); // Add to existing plot

        // Draw truth-matched hits (green)
        gSlopeFits_true[i]->SetMarkerStyle(5);
        gSlopeFits_true[i]->SetMarkerColor(kGreen);
        gSlopeFits_true[i]->SetMarkerSize(2);
        gSlopeFits_true[i]->Draw("p,same"); // Add to existing plot

        // Draw fit lines for visual comparison
        funcSlopeFits[i]->SetLineColor(kRed); // Improved fit
        funcSlopeFits[i]->SetRange(-40, 40);
        funcSlopeFits[i]->SetLineStyle(2);
        funcSlopeFits[i]->SetLineWidth(1);

        funcSlopeFits_orig[i]->SetLineColor(kBlue); // Original fit
        funcSlopeFits_orig[i]->SetRange(-40, 40);
        funcSlopeFits_orig[i]->SetLineStyle(2);
        funcSlopeFits_orig[i]->SetLineWidth(1);

        // Add legend with fit parameters
        TLegend *leg = new TLegend(0.15, 0.65, 0.35, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->AddEntry((TObject *)0, Form("orig. p0: %.2f", funcSlopeFits_orig[i]->GetParameter(0)), "");
        leg->AddEntry((TObject *)0, Form("orig. p1: %.3f", funcSlopeFits_orig[i]->GetParameter(1)), "");
        leg->AddEntry((TObject *)0, Form("p0: %.2f", funcSlopeFits[i]->GetParameter(0)), "");
        leg->AddEntry((TObject *)0, Form("p1: %.3f", funcSlopeFits[i]->GetParameter(1)), "");
        leg->Draw();

        funcSlopeFits_orig[i]->Draw("same"); // Draw original fit line
        funcSlopeFits[i]->Draw("same");      // Draw improved fit line
    }
    c4.SaveAs(Form("%sslopeFits.pdf", outputdir.Data()));

    //----------------------------------------------------------------------
    // Create and save track quality analysis plots
    //----------------------------------------------------------------------
    // Plot angle difference between reconstructed tracklet and true direction
    TCanvas c5("c5", "Angular Accuracy", 800, 600);
    hAngleDiffModules->GetXaxis()->SetTitle("Module");
    hAngleDiffModules->GetYaxis()->SetTitle("Angle difference [rad]");
    hAngleDiffModules->Draw("colz");
    c5.SaveAs(Form("%shAngleDiffModules.pdf", outputdir.Data()));

    // Plot angle difference between reconstructed tracklet and final direction
    TCanvas c5x("c5x", "Angular Accuracy", 800, 600);
    hAngleDiffFinalModules->GetXaxis()->SetTitle("Module");
    hAngleDiffFinalModules->GetYaxis()->SetTitle("Angle difference [rad]");
    hAngleDiffFinalModules->Draw("colz");
    c5x.SaveAs(Form("%shAngleDiffFinalModules.pdf", outputdir.Data()));

    // Plot time difference between extrapolation and tracklet
    TCanvas c6x("c6x", "Time Difference", 800, 600);
    hExtrapolatedTimeVsTracklet->GetXaxis()->SetTitle("Time difference (extrapol - tracklets) [ns]");
    hExtrapolatedTimeVsTracklet->GetYaxis()->SetTitle("Counts");
    hExtrapolatedTimeVsTracklet->Draw();
    c6x.SaveAs(Form("%shExtrapolatedTimeVsTracklet.pdf", outputdir.Data()));

    // Plot time difference between extrapolated and true time
    TCanvas c6y("c6y", "Time Difference", 800, 600);
    hExtrapolatedTimeVsTrue->GetXaxis()->SetTitle("Time difference (extrapol - true) [ns]");
    hExtrapolatedTimeVsTrue->GetYaxis()->SetTitle("Counts");
    hExtrapolatedTimeVsTrue->Draw();
    c6y.SaveAs(Form("%shExtrapolatedTimeVsTrue.pdf", outputdir.Data()));

    //----------------------------------------------------------------------
    // Create and save track matching accuracy plots
    //----------------------------------------------------------------------
    // Plot Y-Z position matching accuracy
    TCanvas c7("c7", "Y-Z Position Matching", 800, 600);
    c7.SetLogz(); // Logarithmic scale for better visibility of distributions
    hMatchdYdZ->GetXaxis()->SetTitle("Y position difference [mm]");
    hMatchdYdZ->GetYaxis()->SetTitle("Z position difference [mm]");
    hMatchdYdZ->Draw("colz");
    c7.SaveAs(Form("%shMatchdYdZ.pdf", outputdir.Data()));

    // Plot X-Z position matching accuracy after cuts
    TCanvas c7a("c7a", "Y-Z Position Matching", 800, 600);
    c7a.SetLogz(); // Logarithmic scale for better visibility of distributions
    hMatchdYdZCut->GetXaxis()->SetTitle("Y position difference [mm]");
    hMatchdYdZCut->GetYaxis()->SetTitle("Z position difference [mm]");
    hMatchdYdZCut->Draw("colz");
    c7a.SaveAs(Form("%shMatchdYdZCut.pdf", outputdir.Data()));

    // Plot X-Z position matching accuracy
    TCanvas c8("c8", "X-Z Position Matching", 800, 600);
    c8.SetLogz();
    hMatchdXdZ->GetXaxis()->SetTitle("X position difference [mm]");
    hMatchdXdZ->GetYaxis()->SetTitle("Z position difference [mm]");
    hMatchdXdZ->Draw("colz");
    c8.SaveAs(Form("%shMatchdXdZ.pdf", outputdir.Data()));

    // Plot X position matching accuracy after cuts
    TCanvas c8a("c8a", "X-Z Position Matching", 800, 600);
    c8a.SetLogz();
    hMatchdXdZCut->GetXaxis()->SetTitle("X position difference [mm]");
    hMatchdXdZCut->GetYaxis()->SetTitle("Z position difference [mm]");
    hMatchdXdZCut->Draw("colz");
    c8a.SaveAs(Form("%shMatchdXdZCut.pdf", outputdir.Data()));

    // Plot X, Y, Z position matching distributions separately
    TCanvas c9("c9", "X Position Matching", 800, 600);
    hMatchdX->GetXaxis()->SetTitle("X position difference [mm]");
    hMatchdX->GetYaxis()->SetTitle("Counts");
    hMatchdX->Draw();
    c9.SaveAs(Form("%shMatchdX.pdf", outputdir.Data()));

    TCanvas c10("c10", "Y Position Matching", 800, 600);
    hMatchdY->GetXaxis()->SetTitle("Y position difference [mm]");
    hMatchdY->GetYaxis()->SetTitle("Counts");
    hMatchdY->Draw();
    c10.SaveAs(Form("%shMatchdY.pdf", outputdir.Data()));

    // plot all 8 hMatchdXModule, hMatchdYModules, hMatchdZModules on their own 3x3 canvas in separate pads
    TCanvas c10a1("c10a1", "X Position Matching", 800, 600);
    c10a1.Divide(3, 3);
    for (int i = 0; i < 9; i++)
    {
        c10a1.cd(i + 1);
        hMatchdXModule[i]->GetXaxis()->SetTitle("X position difference [mm]");
        hMatchdXModule[i]->GetYaxis()->SetTitle("Counts");
        hMatchdXModule[i]->Draw();
    }
    c10a1.SaveAs(Form("%shMatchdXModule.pdf", outputdir.Data()));

    TCanvas c10b("c10b", "Y Position Matching", 800, 600);
    c10b.Divide(3, 3);
    for (int i = 0; i < 9; i++)
    {
        c10b.cd(i + 1);
        hMatchdYModule[i]->GetXaxis()->SetTitle("Y position difference [mm]");
        hMatchdYModule[i]->GetYaxis()->SetTitle("Counts");
        hMatchdYModule[i]->Draw();
    }
    c10b.SaveAs(Form("%shMatchdYModule.pdf", outputdir.Data()));
    TCanvas c10c("c10c", "Z Position Matching", 800, 600);
    c10c.Divide(3, 3);
    for (int i = 0; i < 9; i++)
    {
        c10c.cd(i + 1);
        hMatchdZModule[i]->GetXaxis()->SetTitle("Z position difference [mm]");
        hMatchdZModule[i]->GetYaxis()->SetTitle("Counts");
        hMatchdZModule[i]->Draw();
        //add a gaussian fit similar to the code below
        TF1 *gaus = new TF1("gaus", "gaus", -200, 200);
        hMatchdZModule[i]->Fit(gaus, "R");
        gaus->SetLineColor(kRed);
        gaus->SetLineWidth(2);
        gaus->Draw("same");
        //add a legend with the mean value
        TLegend *leg2 = new TLegend(0.15, 0.65, 0.35, 0.85);
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->SetTextSize(0.04);
        leg2->AddEntry((TObject *)0, Form("mean: %.2f", gaus->GetParameter(1)), "");
        leg2->AddEntry((TObject *)0, Form("sigma: %.2f", gaus->GetParameter(2)), "");
        leg2->Draw();
    }
    c10c.SaveAs(Form("%shMatchdZModule.pdf", outputdir.Data()));
    
    // make a four pad plot with hMatchdXStation, hMatchdYStation, hMatchdZStation
    TCanvas c10d("c10d", "X Position Matching", 1400, 600);
    c10d.Divide(4, 1);
    for (int i = 0; i < 4; i++)
    {
        c10d.cd(i + 1);
        hMatchdXStation[i]->GetXaxis()->SetTitle("X position difference [mm]");
        hMatchdXStation[i]->GetYaxis()->SetTitle("Counts");
        hMatchdXStation[i]->Draw();
    }
    c10d.SaveAs(Form("%shMatchdXStation.pdf", outputdir.Data()));
    TCanvas c10e("c10e", "Y Position Matching", 1400, 600);
    c10e.Divide(4, 1);
    for (int i = 0; i < 4; i++)
    {
        c10e.cd(i + 1);
        hMatchdYStation[i]->GetXaxis()->SetTitle("Y position difference [mm]");
        hMatchdYStation[i]->GetYaxis()->SetTitle("Counts");
        hMatchdYStation[i]->Draw();
    }
    c10e.SaveAs(Form("%shMatchdYStation.pdf", outputdir.Data()));
    TCanvas c10f("c10f", "Z Position Matching", 1400, 600);
    c10f.Divide(4, 1);
    for (int i = 0; i < 4; i++)
    {
        c10f.cd(i + 1);
        hMatchdZStation[i]->GetXaxis()->SetTitle("Z position difference [mm]");
        hMatchdZStation[i]->GetYaxis()->SetTitle("Counts");
        hMatchdZStation[i]->Draw();
    }
    c10f.SaveAs(Form("%shMatchdZStation.pdf", outputdir.Data()));

    //plot all four stations on the same canvas
    TCanvas c10g("c10g", "Z Position Matching", 800, 600);
    hMatchdZStation[0]->SetLineColor(kRed);
    hMatchdZStation[1]->SetLineColor(kBlue);
    hMatchdZStation[2]->SetLineColor(kGreen);
    hMatchdZStation[3]->SetLineColor(kBlack);
    hMatchdZStation[0]->SetLineWidth(2);
    hMatchdZStation[1]->SetLineWidth(2);    
    hMatchdZStation[2]->SetLineWidth(2);
    hMatchdZStation[3]->SetLineWidth(2);
    hMatchdZStation[0]->GetXaxis()->SetTitle("Z position difference [mm]");
    hMatchdZStation[0]->GetYaxis()->SetTitle("Counts");
    hMatchdZStation[0]->Draw();
    hMatchdZStation[1]->Draw("same");
    hMatchdZStation[2]->Draw("same");
    hMatchdZStation[3]->Draw("same");
    TLegend *leg3 = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetTextSize(0.04);
    leg3->AddEntry(hMatchdZStation[0], "Station 0", "l");
    leg3->AddEntry(hMatchdZStation[1], "Station 1", "l");
    leg3->AddEntry(hMatchdZStation[2], "Station 2", "l");
    leg3->AddEntry(hMatchdZStation[3], "Station 3", "l");
    leg3->Draw();
    c10g.SaveAs(Form("%shMatchdZStationAll.pdf", outputdir.Data()));

    // do the same for dX and dY
    TCanvas c10h("c10h", "X Position Matching", 800, 600);
    hMatchdXStation[0]->SetLineColor(kRed);
    hMatchdXStation[1]->SetLineColor(kBlue);
    hMatchdXStation[2]->SetLineColor(kGreen);
    hMatchdXStation[3]->SetLineColor(kBlack);
    hMatchdXStation[0]->SetLineWidth(2);
    hMatchdXStation[1]->SetLineWidth(2);    
    hMatchdXStation[2]->SetLineWidth(2);
    hMatchdXStation[3]->SetLineWidth(2);
    hMatchdXStation[0]->GetXaxis()->SetTitle("X position difference [mm]");
    hMatchdXStation[0]->GetYaxis()->SetTitle("Counts");
    hMatchdXStation[0]->Draw();
    hMatchdXStation[1]->Draw("same");
    hMatchdXStation[2]->Draw("same");
    hMatchdXStation[3]->Draw("same");
    TLegend *leg4 = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->SetTextSize(0.04);
    leg4->AddEntry(hMatchdXStation[0], "Station 0", "l");
    leg4->AddEntry(hMatchdXStation[1], "Station 1", "l");
    leg4->AddEntry(hMatchdXStation[2], "Station 2", "l");
    leg4->AddEntry(hMatchdXStation[3], "Station 3", "l");
    leg4->Draw();
    c10h.SaveAs(Form("%shMatchdXStationAll.pdf", outputdir.Data()));

    TCanvas c10i("c10i", "Y Position Matching", 800, 600);
    hMatchdYStation[0]->SetLineColor(kRed);
    hMatchdYStation[1]->SetLineColor(kBlue);
    hMatchdYStation[2]->SetLineColor(kGreen);
    hMatchdYStation[3]->SetLineColor(kBlack);
    hMatchdYStation[0]->SetLineWidth(2);
    hMatchdYStation[1]->SetLineWidth(2);    
    hMatchdYStation[2]->SetLineWidth(2);
    hMatchdYStation[3]->SetLineWidth(2);
    hMatchdYStation[0]->GetXaxis()->SetTitle("Y position difference [mm]");
    hMatchdYStation[0]->GetYaxis()->SetTitle("Counts");
    hMatchdYStation[0]->Draw();
    hMatchdYStation[1]->Draw("same");
    hMatchdYStation[2]->Draw("same");
    hMatchdYStation[3]->Draw("same");
    TLegend *leg5 = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg5->SetBorderSize(0);
    leg5->SetFillStyle(0);
    leg5->SetTextSize(0.04);
    leg5->AddEntry(hMatchdYStation[0], "Station 0", "l");
    leg5->AddEntry(hMatchdYStation[1], "Station 1", "l");
    leg5->AddEntry(hMatchdYStation[2], "Station 2", "l");
    leg5->AddEntry(hMatchdYStation[3], "Station 3", "l");
    leg5->Draw();
    c10i.SaveAs(Form("%shMatchdYStationAll.pdf", outputdir.Data()));


    TCanvas c11("c11", "Z Position Matching", 800, 600);
    hMatchdZ->GetXaxis()->SetTitle("Z position difference [mm]");
    hMatchdZ->GetYaxis()->SetTitle("Counts");
    hMatchdZ->Draw();
    //add a gaussian fit around 0 to the histogram
    TF1 *gaus = new TF1("gaus", "gaus", -200, 150);
    hMatchdZ->Fit(gaus, "R");
    gaus->SetLineColor(kRed);
    gaus->SetLineWidth(2);
    gaus->Draw("same");
    //add a legend with the mean value
    TLegend *leg2 = new TLegend(0.15, 0.65, 0.35, 0.85);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry((TObject *)0, Form("mean: %.2f", gaus->GetParameter(1)), "");
    leg2->AddEntry((TObject *)0, Form("sigma: %.2f", gaus->GetParameter(2)), "");
    leg2->Draw();
    c11.SaveAs(Form("%shMatchdZ.pdf", outputdir.Data()));

    // plot hDiffBarPositionToTrueHitPositionX, hDiffBarPositionToTrueHitPositionY, hDiffBarPositionToTrueHitPositionZ
    TCanvas c12a("c12a", "X Position Matching", 800, 600);
    hDiffBarPositionToTrueHitPositionX->GetXaxis()->SetTitle("X position difference [mm]");
    hDiffBarPositionToTrueHitPositionX->GetYaxis()->SetTitle("Counts");
    hDiffBarPositionToTrueHitPositionX->Draw();
    c12a.SaveAs(Form("%shDiffBarPositionToTrueHitPositionX.pdf", outputdir.Data()));

    TCanvas c12b("c12b", "Y Position Matching", 800, 600);
    hDiffBarPositionToTrueHitPositionY->GetXaxis()->SetTitle("Y position difference [mm]");
    hDiffBarPositionToTrueHitPositionY->GetYaxis()->SetTitle("Counts");
    hDiffBarPositionToTrueHitPositionY->Draw();
    c12b.SaveAs(Form("%shDiffBarPositionToTrueHitPositionY.pdf", outputdir.Data()));

    TCanvas c12c("c12c", "Z Position Matching", 800, 600);
    hDiffBarPositionToTrueHitPositionZ->GetXaxis()->SetTitle("Z position difference [mm]");
    hDiffBarPositionToTrueHitPositionZ->GetYaxis()->SetTitle("Counts");
    hDiffBarPositionToTrueHitPositionZ->Draw();
    c12c.SaveAs(Form("%shDiffBarPositionToTrueHitPositionZ.pdf", outputdir.Data()));

    // plot hMatchdXExtrapolToTracklet
    TCanvas c11a1("c11a1", "X Position Matching", 800, 600);
    hMatchdXExtrapolToTracklet->GetXaxis()->SetTitle("X position difference [mm]");
    hMatchdXExtrapolToTracklet->GetYaxis()->SetTitle("Counts");
    hMatchdXExtrapolToTracklet->Draw();
    c11a1.SaveAs(Form("%shMatchdXExtrapolToTracklet.pdf", outputdir.Data()));

    // plot hMatchdYExtrapolToTracklet
    TCanvas c11b("c11b", "Y Position Matching", 800, 600);
    hMatchdYExtrapolToTracklet->GetXaxis()->SetTitle("Y position difference [mm]");
    hMatchdYExtrapolToTracklet->GetYaxis()->SetTitle("Counts");
    hMatchdYExtrapolToTracklet->Draw();
    c11b.SaveAs(Form("%shMatchdYExtrapolToTracklet.pdf", outputdir.Data()));

    // plot hMatchdZExtrapolToTracklet
    TCanvas c11c("c11c", "Z Position Matching", 800, 600);
    hMatchdZExtrapolToTracklet->GetXaxis()->SetTitle("Z position difference [mm]");
    hMatchdZExtrapolToTracklet->GetYaxis()->SetTitle("Counts");
    hMatchdZExtrapolToTracklet->Draw();
    c11c.SaveAs(Form("%shMatchdZExtrapolToTracklet.pdf", outputdir.Data()));

    // plot hMatchdXExtrapolToTrackletQuality
    TCanvas c11d("c11d", "X Position Matching", 800, 600);
    hMatchdXExtrapolToTrackletQuality->GetXaxis()->SetTitle("X position difference [mm]");
    hMatchdXExtrapolToTrackletQuality->GetYaxis()->SetTitle("Counts");
    hMatchdXExtrapolToTrackletQuality->Draw();
    c11d.SaveAs(Form("%shMatchdXExtrapolToTrackletQuality.pdf", outputdir.Data()));

    // plot hMatchdYExtrapolToTrackletQuality
    TCanvas c11e("c11e", "Y Position Matching", 800, 600);
    hMatchdYExtrapolToTrackletQuality->GetXaxis()->SetTitle("Y position difference [mm]");
    hMatchdYExtrapolToTrackletQuality->GetYaxis()->SetTitle("Counts");
    hMatchdYExtrapolToTrackletQuality->Draw();
    c11e.SaveAs(Form("%shMatchdYExtrapolToTrackletQuality.pdf", outputdir.Data()));

    // plot hMatchdZExtrapolToTrackletQuality
    TCanvas c11f("c11f", "Z Position Matching", 800, 600);
    hMatchdZExtrapolToTrackletQuality->GetXaxis()->SetTitle("Z position difference [mm]");
    hMatchdZExtrapolToTrackletQuality->GetYaxis()->SetTitle("Counts");
    hMatchdZExtrapolToTrackletQuality->Draw();
    c11f.SaveAs(Form("%shMatchdZExtrapolToTrackletQuality.pdf", outputdir.Data()));

    // Plot X, Y, Z position matching distributions separately after cuts
    TCanvas c9a("c9a", "X Position Matching", 800, 600);
    hMatchdXCut->GetXaxis()->SetTitle("X position difference [mm]");
    hMatchdXCut->GetYaxis()->SetTitle("Counts");
    hMatchdXCut->Draw();
    c9a.SaveAs(Form("%shMatchdXCut.pdf", outputdir.Data()));

    TCanvas c10a("c10a", "Y Position Matching", 800, 600);
    hMatchdYCut->GetXaxis()->SetTitle("Y position difference [mm]");
    hMatchdYCut->GetYaxis()->SetTitle("Counts");
    hMatchdYCut->Draw();
    c10a.SaveAs(Form("%shMatchdYCut.pdf", outputdir.Data()));

    TCanvas c11a("c11a", "Z Position Matching", 800, 600);
    hMatchdZCut->GetXaxis()->SetTitle("Z position difference [mm]");
    hMatchdZCut->GetYaxis()->SetTitle("Counts");
    hMatchdZCut->Draw();
    c11a.SaveAs(Form("%shMatchdZCut.pdf", outputdir.Data()));

    //----------------------------------------------------------------------
    // Create and save detector coverage and extrapolation plots
    //----------------------------------------------------------------------
    // Plot hit positions in X-Z projection
    TCanvas c12("c12", "True hit Positions X-Z", 800, 800);
    c12.SetLogz();
    hTrueHitXZ->GetXaxis()->SetTitle("X position [mm]");
    hTrueHitXZ->GetYaxis()->SetTitle("Z position [mm]");
    hTrueHitXZ->Draw("colz");
    c12.SaveAs(Form("%shTrueHitXZ.pdf", outputdir.Data()));

    // Plot extrapolated positions in X-Z projection
    TCanvas c13("c13", "Extrapolated Positions X-Z", 800, 600);
    c13.SetLogz();
    hExtrapolatedXZ->GetXaxis()->SetTitle("X position [mm]");
    hExtrapolatedXZ->GetYaxis()->SetTitle("Z position [mm]");
    hExtrapolatedXZ->Draw("colz");
    c13.SaveAs(Form("%shExtrapolatedXZ.pdf", outputdir.Data()));

    // Plot hit positions in Y-Z projection
    TCanvas c14("c14", "True hit Positions Y-Z", 800, 600);
    c14.SetLogz();
    hTrueHitYZ->GetXaxis()->SetTitle("Y position [mm]");
    hTrueHitYZ->GetYaxis()->SetTitle("Z position [mm]");
    hTrueHitYZ->Draw("colz");
    c14.SaveAs(Form("%shTrueHitYZ.pdf", outputdir.Data()));

    // Plot extrapolated positions in Y-Z projection
    TCanvas c15("c15", "Extrapolated Positions Y-Z", 800, 600);
    c15.SetLogz();
    hExtrapolatedYZ->GetXaxis()->SetTitle("Y position [mm]");
    hExtrapolatedYZ->GetYaxis()->SetTitle("Z position [mm]");
    hExtrapolatedYZ->Draw("colz");
    c15.SaveAs(Form("%shExtrapolatedYZ.pdf", outputdir.Data()));

    //----------------------------------------------------------------------
    // Create and save momentum sensitivity analysis plots
    //----------------------------------------------------------------------
    // Plot position sensitivity to momentum variations in X-Z
    TCanvas c16("c16", "Momentum Sensitivity X-Z", 800, 600);
    c16.SetLogz();
    hExtrapolationMomentumSmearingdXdZ->GetXaxis()->SetTitle("X difference [mm]");
    hExtrapolationMomentumSmearingdXdZ->GetYaxis()->SetTitle("Z difference [mm]");
    hExtrapolationMomentumSmearingdXdZ->Draw("colz");
    c16.SaveAs(Form("%shExtrapolationMomentumSmearingdXdZ.pdf", outputdir.Data()));

    // Plot position sensitivity to momentum variations in Y-Z
    TCanvas c17("c17", "Momentum Sensitivity Y-Z", 800, 600);
    c17.SetLogz();
    hExtrapolationMomentumSmearingdYdZ->GetXaxis()->SetTitle("Y difference [mm]");
    hExtrapolationMomentumSmearingdYdZ->GetYaxis()->SetTitle("Z difference [mm]");
    hExtrapolationMomentumSmearingdYdZ->Draw("colz");
    c17.SaveAs(Form("%shExtrapolationMomentumSmearingdYdZ.pdf", outputdir.Data()));

    //----------------------------------------------------------------------
    // Clean up resources
    //----------------------------------------------------------------------
    if (doTree)
    {
        treeOut->Write();     // Write trajectory tree to file
        rootFileOut->Close(); // Close output file
        delete rootFileOut;   // Free memory
    }

    return 0; // Return success
}