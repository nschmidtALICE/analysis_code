#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <TVector3.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TArrow.h>
#include <TStyle.h>
#include <TColor.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TLegend.h>

/**
 * LHCbMagneticField class
 * Loads the LHCb magnetic field map from four quadrant files and provides
 * trilinear interpolation to evaluate field at arbitrary positions
 */
class LHCbMagneticField {
private:
    // Structure to hold field data for a quadrant
    struct Quadrant {
        double zOffset;                     // Z offset in cm
        std::array<double, 3> Dxyz;         // Grid spacing in cm
        std::array<int, 3> Nxyz;            // Grid dimensions
        std::vector<TVector3> B;            // Field vectors
    };

    // Field map properties
    double x_min, x_max, y_min, y_max, z_min, z_max;
    double dx, dy, dz;
    int nx, ny, nz;

    // 3D storage for field data (in Tesla)
    std::vector<std::vector<std::vector<TVector3>>> field_values;

    // Conversion from CGS (Gauss) to Tesla
    static constexpr double GAUSS_TO_TESLA = 1.0e-4;

    /**
     * Parse a CDF file containing one quadrant of the field map
     * @param filename Path to CDF file
     * @param quadrant Output quadrant data structure to fill
     * @return True if parsing succeeded
     */
    bool parseQuadrantFile(const std::string& filename, Quadrant& quadrant) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open field map file: " << filename << std::endl;
            return false;
        }

        std::cout << "Reading field map from " << filename << std::endl;

        // Skip header until GEOMETRY section
        std::string line;
        bool foundGeometry = false;
        while (std::getline(file, line)) {
            if (line.find("GEOMETRY:") != std::string::npos) {
                foundGeometry = true;
                break;
            }
        }

        if (!foundGeometry) {
            std::cerr << "Error: GEOMETRY section not found in " << filename << std::endl;
            return false;
        }

        // Skip comments after GEOMETRY section
        while (std::getline(file, line)) {
            if (!line.empty() && line[0] != '#') break;
        }

        // Parse geometry parameters
        std::istringstream iss(line);
        double dx, dy, dz;
        int nx, ny, nz;
        double z_offset;
        
        if (!(iss >> dx >> dy >> dz >> nx >> ny >> nz >> z_offset)) {
            std::cerr << "Error parsing geometry line: " << line << std::endl;
            return false;
        }

        // Store parameters in quadrant
        quadrant.Dxyz = {dx, dy, dz};
        quadrant.Nxyz = {nx, ny, nz};
        quadrant.zOffset = z_offset;

        // Reserve space for field values
        const size_t dataPoints = nx * ny * nz;
        quadrant.B.reserve(dataPoints);

        // Read field values
        double x, y, z, Bx, By, Bz;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            std::istringstream dataLine(line);
            if (dataLine >> Bx >> By >> Bz) {
                // Convert Gauss to Tesla and store field vector
                quadrant.B.emplace_back(
                    Bx * GAUSS_TO_TESLA,
                    By * GAUSS_TO_TESLA, 
                    Bz * GAUSS_TO_TESLA
                );
            }
        }

        // Check if we read the expected number of data points
        if (quadrant.B.size() != dataPoints) {
            std::cerr << "Warning: Expected " << dataPoints << " data points, but read " 
                      << quadrant.B.size() << std::endl;
        }

        std::cout << "Successfully read quadrant: "
                  << "dx=" << dx << ", dy=" << dy << ", dz=" << dz
                  << ", nx=" << nx << ", ny=" << ny << ", nz=" << nz
                  << ", z_offset=" << z_offset << std::endl;

        return true;
    }

    /**
     * Combine four quadrants into a complete field map
     * Following the algorithm in MagneticFieldGridReader::fillGridFromQuadrants
     * @param quadrants Array of four quadrants (ordered: 1,2,3,4)
     * @return True if combination succeeded
     */
    bool combineQuadrants(const std::array<Quadrant, 4>& quadrants) {
        // Check that all quadrants have the same grid spacing and dimensions
        for (size_t i = 1; i < quadrants.size(); ++i) {
            if (quadrants[i].Nxyz != quadrants[0].Nxyz ||
                std::abs(quadrants[i].Dxyz[0] - quadrants[0].Dxyz[0]) > 1e-6 ||
                std::abs(quadrants[i].Dxyz[1] - quadrants[0].Dxyz[1]) > 1e-6 ||
                std::abs(quadrants[i].Dxyz[2] - quadrants[0].Dxyz[2]) > 1e-6 ||
                std::abs(quadrants[i].zOffset - quadrants[0].zOffset) > 1e-6) {
                std::cerr << "Error: Quadrant " << i << " has different parameters than quadrant 0" << std::endl;
                return false;
            }
        }

        // Get dimensions from first quadrant
        const int Nxquad = quadrants[0].Nxyz[0];
        const int Nyquad = quadrants[0].Nxyz[1];
        const int Nzquad = quadrants[0].Nxyz[2];

        // Set global grid dimensions and spacing
        nx = 2 * Nxquad - 1;  // Account for overlap
        ny = 2 * Nyquad - 1;  // Account for overlap
        nz = Nzquad;
        dx = quadrants[0].Dxyz[0];
        dy = quadrants[0].Dxyz[1];
        dz = quadrants[0].Dxyz[2];

        // Calculate grid boundaries (in cm)
        x_min = -dx * (Nxquad - 1);
        x_max = dx * (Nxquad - 1);
        y_min = -dy * (Nyquad - 1);
        y_max = dy * (Nyquad - 1);
        z_min = quadrants[0].zOffset;
        z_max = quadrants[0].zOffset + dz * (Nzquad - 1);

        // Initialize field values array
        field_values.resize(nx, std::vector<std::vector<TVector3>>(
            ny, std::vector<TVector3>(nz, TVector3(0, 0, 0))));

        // Map quadrants to proper field values based on their orientation
        for (int iz = 0; iz < Nzquad; ++iz) {
            for (int iy = 0; iy < Nyquad; ++iy) {
                for (int ix = 0; ix < Nxquad; ++ix) {
                    // Get quadrant indices
                    const int q0_index = Nxquad * (Nyquad * iz + iy) + ix;
                    
                    // Ensure we don't exceed the bounds
                    if (q0_index >= quadrants[0].B.size()) continue;

                    // 1st quadrant (positive x, positive y)
                    const TVector3& Q1 = quadrants[0].B[q0_index];
                    field_values[Nxquad + ix - 1][Nyquad + iy - 1][iz] = Q1;

                    // 2nd quadrant (negative x, positive y)
                    const TVector3& Q2 = quadrants[1].B[q0_index];
                    // Flip x component sign
                    field_values[Nxquad - ix - 1][Nyquad + iy - 1][iz] = 
                        TVector3(-Q2.X(), Q2.Y(), Q2.Z());

                    // 3rd quadrant (positive x, negative y)
                    const TVector3& Q3 = quadrants[2].B[q0_index];
                    // Flip z component sign
                    field_values[Nxquad + ix - 1][Nyquad - iy - 1][iz] = 
                        TVector3(Q3.X(), Q3.Y(), -Q3.Z());

                    // 4th quadrant (negative x, negative y)
                    const TVector3& Q4 = quadrants[3].B[q0_index];
                    // Flip x and z component signs
                    field_values[Nxquad - ix - 1][Nyquad - iy - 1][iz] = 
                        TVector3(-Q4.X(), Q4.Y(), -Q4.Z());
                }
            }
        }

        std::cout << "Field map combined successfully:" << std::endl
                  << "Grid: " << nx << "x" << ny << "x" << nz << std::endl
                  << "Bounds (cm): x=[" << x_min << "," << x_max << "], "
                  << "y=[" << y_min << "," << y_max << "], "
                  << "z=[" << z_min << "," << z_max << "]" << std::endl;

        return true;
    }

public:
    /**
     * Constructor: Load LHCb field from four quadrant files
     * @param q1_file Path to first quadrant file
     * @param q2_file Path to second quadrant file 
     * @param q3_file Path to third quadrant file
     * @param q4_file Path to fourth quadrant file
     */
    LHCbMagneticField(const std::string& q1_file, const std::string& q2_file,
                      const std::string& q3_file, const std::string& q4_file) {
        std::array<Quadrant, 4> quadrants;

        // Parse all four quadrant files
        if (!parseQuadrantFile(q1_file, quadrants[0]) ||
            !parseQuadrantFile(q2_file, quadrants[1]) ||
            !parseQuadrantFile(q3_file, quadrants[2]) ||
            !parseQuadrantFile(q4_file, quadrants[3])) {
            throw std::runtime_error("Failed to parse quadrant files");
        }

        // Combine quadrants into a single field map
        if (!combineQuadrants(quadrants)) {
            throw std::runtime_error("Failed to combine quadrants");
        }
    }

    /**
     * Get magnetic field vector at a given position using trilinear interpolation
     * @param position Position vector (in meters)
     * @return Magnetic field vector (in Tesla)
     */
    TVector3 getFieldAt(const TVector3& position) const {
        // Convert position from meters to cm for grid lookup
        double x = position.X() * 100.0;
        double y = position.Y() * 100.0;
        double z = position.Z() * 100.0;

        // Check if position is within field map bounds
        if (x < x_min || x > x_max || y < y_min || y > y_max || z < z_min || z > z_max) {
            return TVector3(0, 0, 0);
        }

        // Find grid cell containing this position
        double x_normalized = (x - x_min) / dx;
        double y_normalized = (y - y_min) / dy;
        double z_normalized = (z - z_min) / dz;

        int i0 = int(std::floor(x_normalized));
        int j0 = int(std::floor(y_normalized));
        int k0 = int(std::floor(z_normalized));

        // Clamp indices to valid range
        i0 = std::max(0, std::min(i0, nx - 2));
        j0 = std::max(0, std::min(j0, ny - 2));
        k0 = std::max(0, std::min(k0, nz - 2));

        int i1 = i0 + 1;
        int j1 = j0 + 1;
        int k1 = k0 + 1;

        // Calculate interpolation weights
        double wx = x_normalized - i0;
        double wy = y_normalized - j0;
        double wz = z_normalized - k0;

        // Perform trilinear interpolation
        TVector3 B = (1-wx)*(1-wy)*(1-wz)*field_values[i0][j0][k0] + 
                     (1-wx)*(1-wy)*wz*field_values[i0][j0][k1] +
                     (1-wx)*wy*(1-wz)*field_values[i0][j1][k0] + 
                     (1-wx)*wy*wz*field_values[i0][j1][k1] +
                     wx*(1-wy)*(1-wz)*field_values[i1][j0][k0] + 
                     wx*(1-wy)*wz*field_values[i1][j0][k1] +
                     wx*wy*(1-wz)*field_values[i1][j1][k0] + 
                     wx*wy*wz*field_values[i1][j1][k1];

        return B;
    }

    /**
     * Create visualization of the magnetic field
     * @param outputFile Name of the output ROOT file for the plots
     */
    void visualizeField(const std::string& outputFile) const {
        TFile* file = new TFile(outputFile.c_str(), "RECREATE");
        
        // Set style for better visualization
        gStyle->SetPalette(kRainBow);
        gStyle->SetOptStat(0);
        
        // Create canvases for different views
        TCanvas* cXZ = new TCanvas("cXZ", "B-field in XZ plane (y=0)", 900, 700);
        TCanvas* cYZ = new TCanvas("cYZ", "B-field in YZ plane (x=0)", 900, 700);
        TCanvas* cXY = new TCanvas("cXY", "B-field in XY plane (z=center)", 900, 700);
        TCanvas* cProfile = new TCanvas("cProfile", "B-field magnitude along Z-axis", 900, 500);
        
        // Create histograms for field magnitude
        int nBinsX = 100;
        int nBinsY = 100;
        int nBinsZ = 200;
        
        // Convert bounds from cm to m for plot labels
        double x_min_m = x_min / 100.0;
        double x_max_m = x_max / 100.0;
        double y_min_m = y_min / 100.0;
        double y_max_m = y_max / 100.0;
        double z_min_m = z_min / 100.0;
        double z_max_m = z_max / 100.0;
        
        TH2D* hXZ = new TH2D("hXZ", "B-field magnitude in XZ plane (y=0);Z [m];X [m];|B| [T]", 
                             nBinsZ, z_min_m, z_max_m, nBinsX, x_min_m, x_max_m);
        
        TH2D* hYZ = new TH2D("hYZ", "B-field magnitude in YZ plane (x=0);Z [m];Y [m];|B| [T]",
                             nBinsZ, z_min_m, z_max_m, nBinsY, y_min_m, y_max_m);
        
        double z_center_m = (z_min_m + z_max_m) / 2.0;
        TH2D* hXY = new TH2D("hXY", "B-field magnitude in XY plane (z=center);X [m];Y [m];|B| [T]",
                             nBinsX, x_min_m, x_max_m, nBinsY, y_min_m, y_max_m);
        
        TGraph* gProfile = new TGraph();
        gProfile->SetTitle("B-field magnitude along Z-axis (x=0, y=0);Z [m];|B| [T]");
        gProfile->SetLineWidth(2);
        gProfile->SetLineColor(kBlue);
        gProfile->SetMarkerStyle(20);
        gProfile->SetMarkerSize(0.5);
        
        // Fill histograms with field magnitudes
        double maxField = 0;
        std::vector<TArrow*> arrowsXZ, arrowsYZ, arrowsXY;
        
        // Sample points for vector field visualization
        int arrowsPerDim = 20; // Number of arrows in each dimension
        double arrowScale = 0.5; // Scale factor for arrow length
        
        // Fill XZ plane (y=0)
        for (int ix = 0; ix < nBinsX; ix++) {
            double x = x_min_m + ix * (x_max_m - x_min_m) / nBinsX;
            
            for (int iz = 0; iz < nBinsZ; iz++) {
                double z = z_min_m + iz * (z_max_m - z_min_m) / nBinsZ;
                
                TVector3 pos(x, 0, z);
                TVector3 field = getFieldAt(pos);
                double magnitude = field.Mag();
                
                hXZ->SetBinContent(iz+1, ix+1, magnitude);
                maxField = std::max(maxField, magnitude);
                
                // Add arrows for vector field visualization (sparse grid)
                if (ix % (nBinsX/arrowsPerDim) == 0 && iz % (nBinsZ/arrowsPerDim) == 0) {
                    TArrow* arrow = new TArrow(z, x, 
                                             z + field.Z() * arrowScale, 
                                             x + field.X() * arrowScale,
                                             0.01, ">");
                    arrow->SetLineColor(kBlack);
                    arrowsXZ.push_back(arrow);
                }
            }
        }
        
        // Fill YZ plane (x=0)
        for (int iy = 0; iy < nBinsY; iy++) {
            double y = y_min_m + iy * (y_max_m - y_min_m) / nBinsY;
            
            for (int iz = 0; iz < nBinsZ; iz++) {
                double z = z_min_m + iz * (z_max_m - z_min_m) / nBinsZ;
                
                TVector3 pos(0, y, z);
                TVector3 field = getFieldAt(pos);
                double magnitude = field.Mag();
                
                hYZ->SetBinContent(iz+1, iy+1, magnitude);
                
                // Add arrows for vector field visualization (sparse grid)
                if (iy % (nBinsY/arrowsPerDim) == 0 && iz % (nBinsZ/arrowsPerDim) == 0) {
                    TArrow* arrow = new TArrow(z, y, 
                                             z + field.Z() * arrowScale, 
                                             y + field.Y() * arrowScale,
                                             0.01, ">");
                    arrow->SetLineColor(kBlack);
                    arrowsYZ.push_back(arrow);
                }
            }
        }
        
        // Fill XY plane (z=center)
        for (int ix = 0; ix < nBinsX; ix++) {
            double x = x_min_m + ix * (x_max_m - x_min_m) / nBinsX;
            
            for (int iy = 0; iy < nBinsY; iy++) {
                double y = y_min_m + iy * (y_max_m - y_min_m) / nBinsY;
                
                TVector3 pos(x, y, z_center_m);
                TVector3 field = getFieldAt(pos);
                double magnitude = field.Mag();
                
                hXY->SetBinContent(ix+1, iy+1, magnitude);
                
                // Add arrows for vector field visualization (sparse grid)
                if (ix % (nBinsX/arrowsPerDim) == 0 && iy % (nBinsY/arrowsPerDim) == 0) {
                    TArrow* arrow = new TArrow(x, y, 
                                             x + field.X() * arrowScale, 
                                             y + field.Y() * arrowScale,
                                             0.01, ">");
                    arrow->SetLineColor(kBlack);
                    arrowsXY.push_back(arrow);
                }
            }
        }
        
        // Fill Z-axis profile
        for (int iz = 0; iz < nBinsZ; iz++) {
            double z = z_min_m + iz * (z_max_m - z_min_m) / nBinsZ;
            TVector3 field = getFieldAt(TVector3(0, 0, z));
            double magnitude = field.Mag();
            gProfile->SetPoint(iz, z, magnitude);
        }
        
        // Draw XZ plane
        cXZ->cd();
        hXZ->Draw("COLZ");
        for (auto& arrow : arrowsXZ) {
            arrow->Draw();
        }
        cXZ->Write();
        
        // Draw YZ plane
        cYZ->cd();
        hYZ->Draw("COLZ");
        for (auto& arrow : arrowsYZ) {
            arrow->Draw();
        }
        cYZ->Write();
        
        // Draw XY plane
        cXY->cd();
        hXY->Draw("COLZ");
        for (auto& arrow : arrowsXY) {
            arrow->Draw();
        }
        cXY->Write();
        
        // Draw Z profile
        cProfile->cd();
        gProfile->Draw("ALP");
        
        // Add vertical lines indicating magnet boundaries
        double magnet_start = 2.5;  // approximate magnet position start (m)
        double magnet_end = 8.0;    // approximate magnet position end (m)
        
        TLine* line1 = new TLine(magnet_start, 0, magnet_start, maxField*1.1);
        line1->SetLineColor(kRed);
        line1->SetLineStyle(2);
        line1->SetLineWidth(2);
        line1->Draw();
        
        TLine* line2 = new TLine(magnet_end, 0, magnet_end, maxField*1.1);
        line2->SetLineColor(kRed);
        line2->SetLineStyle(2);
        line2->SetLineWidth(2);
        line2->Draw();
        
        TLegend* leg = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg->AddEntry(line1, "Magnet boundaries", "l");
        leg->Draw();
        
        cProfile->Write();
        
        // Create 3D field visualization
        TCanvas* c3D = new TCanvas("c3D", "3D Field Visualization", 900, 700);
        TH3D* h3D = new TH3D("h3D", "3D Field Magnitude;X [m];Y [m];Z [m]",
                            nBinsX/2, x_min_m, x_max_m,
                            nBinsY/2, y_min_m, y_max_m,
                            nBinsZ/2, z_min_m, z_max_m);
                            
        for (int ix = 0; ix < nBinsX/2; ix++) {
            double x = x_min_m + ix * (x_max_m - x_min_m) / (nBinsX/2);
            for (int iy = 0; iy < nBinsY/2; iy++) {
                double y = y_min_m + iy * (y_max_m - y_min_m) / (nBinsY/2);
                for (int iz = 0; iz < nBinsZ/2; iz++) {
                    double z = z_min_m + iz * (z_max_m - z_min_m) / (nBinsZ/2);
                    
                    TVector3 field = getFieldAt(TVector3(x, y, z));
                    h3D->SetBinContent(ix+1, iy+1, iz+1, field.Mag());
                }
            }
        }
        
        c3D->cd();
        h3D->SetMarkerStyle(20);
        h3D->Draw("BOX2Z");
        c3D->Write();
        
        // Clean up
        file->Write();
        file->Close();
        
        std::cout << "Field visualization saved to " << outputFile << std::endl;
    }

    // Getters for field boundaries
    double getXMin() const { return x_min / 100.0; }  // in meters
    double getXMax() const { return x_max / 100.0; }
    double getYMin() const { return y_min / 100.0; }
    double getYMax() const { return y_max / 100.0; }
    double getZMin() const { return z_min / 100.0; }
    double getZMax() const { return z_max / 100.0; }
};

/**
 * Load LHCb magnetic field and return as a shared pointer
 * @param q1_file Path to first quadrant file
 * @param q2_file Path to second quadrant file
 * @param q3_file Path to third quadrant file
 * @param q4_file Path to fourth quadrant file
 * @return Shared pointer to LHCbMagneticField object
 */
std::shared_ptr<LHCbMagneticField> loadLHCbMagneticField(
    const std::string& q1_file,
    const std::string& q2_file,
    const std::string& q3_file,
    const std::string& q4_file)
{
    try {
        return std::make_shared<LHCbMagneticField>(q1_file, q2_file, q3_file, q4_file);
    }
    catch (const std::exception& e) {
        std::cerr << "Error loading LHCb magnetic field: " << e.what() << std::endl;
        return nullptr;
    }
}

/**
 * Example usage function to load field and visualize it
 */
void visualizeLHCbField(const std::string& q1_file,
                       const std::string& q2_file,
                       const std::string& q3_file,
                       const std::string& q4_file,
                       const std::string& output_file = "lhcb_field_visualization.root")
{
    auto field = loadLHCbMagneticField(q1_file, q2_file, q3_file, q4_file);
    if (field) {
        field->visualizeField(output_file);
        
        // Test field values at some sample points
        std::cout << "\nField samples:" << std::endl;
        std::cout << "Field at (0,0,3.0m): " << field->getFieldAt(TVector3(0,0,3)).Mag() << " T" << std::endl;
        std::cout << "Field at (0,0,5.0m): " << field->getFieldAt(TVector3(0,0,5)).Mag() << " T" << std::endl;
        std::cout << "Field at (0,0,7.0m): " << field->getFieldAt(TVector3(0,0,7)).Mag() << " T" << std::endl;
    }
}

/**
 * Main entry point for standalone testing
 */
int FieldLoaderAndPlotter() {
    // if (argc < 5) {
    //     std::cerr << "Usage: " << argv[0] << " q1_file q2_file q3_file q4_file [output_file]" << std::endl;
    //     return 1;
    // }
    TString q1File = "/home/niviths/Downloads/field.v5r0.c1.down.cdf";
    TString q2File = "/home/niviths/Downloads/field.v5r0.c2.down.cdf";
    TString q3File = "/home/niviths/Downloads/field.v5r0.c3.down.cdf";
    TString q4File = "/home/niviths/Downloads/field.v5r0.c4.down.cdf";
    std::string output_file = "lhcb_field_visualization.root";
    // if (argc >= 6) {
    //     output_file = argv[5];
    // }
    
    visualizeLHCbField(q1File.Data(), q2File.Data(), q3File.Data(), q4File.Data(), output_file);
    return 0;
}