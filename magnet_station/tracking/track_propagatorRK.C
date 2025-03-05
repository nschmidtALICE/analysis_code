/*
 * Charged Particle Tracking in Magnetic Fields
 *
 * This code simulates the propagation of charged particles through a magnetic field,
 * using the relativistic equations of motion and Runge-Kutta 4th order integration.
 *
 * INPUTS:
 * - Magnetic field map (ASCII format with field values)
 * - Particle properties (type, charge, momentum)
 * - Geometry boundaries (magnet aperture or detector surfaces)
 *
 * PHYSICS:
 * - Full relativistic treatment of particle motion
 * - Lorentz force equation: F = q(v × B)
 * - 4th-order Runge-Kutta integration
 *
 * OUTPUTS:
 * - ROOT file with trajectory information
 * - ASCII data file with particle tracks
 * - Visualization plots of particle trajectories
 *
 * Usage: root -l 'track_propagatorRK.C("field_map_lhcb.txt", N)'
 * where N is the number of particles to simulate
 */

// Standard C++ libraries
#include <iostream> // For input/output operations
#include <fstream>  // For file handling
#include <string>   // For string operations
#include <vector>   // For dynamic arrays
#include <cmath>    // For mathematical functions

// ROOT libraries for physics calculations and visualization
#include "TVector3.h"       // 3-vector class for positions and momenta
#include "TLorentzVector.h" // 4-vector class for relativistic calculations
#include "TFile.h"          // ROOT file handling
#include "TTree.h"          // ROOT tree for data storage
#include "TH3.h"            // 3D histograms
#include "TGraph.h"         // For plotting trajectories
#include "TCanvas.h"        // For visualization windows
#include "TMultiGraph.h"    // For multiple trajectory plots
#include "TStyle.h"         // For plot styling
#include "TBox.h"           // For drawing boundaries
#include "TLine.h"          // For drawing lines
#include "TLegend.h"        // For plot legends

/**
 * Point class - Simple container for 3D coordinates with time component
 * Used for describing geometry points and boundaries
 */
class Point
{
public:
    double x, y, z, t; // Coordinates in space (and optionally time)

    // Constructor with default value for time component
    Point(double x, double y, double z, double t = 0.0) : x(x), y(y), z(z), t(t) {}
};

/**
 * Detector surfaces definition
 * Each surface is defined by 4 corner points forming a quadrilateral in 3D space
 * Coordinates are in millimeters
 * Format: {{x1,y1,z1}, {x2,y2,z2}, {x3,y3,z3}, {x4,y4,z4}}
 * These surfaces are used to define realistic detector geometry for boundary checks
 */
std::vector<std::vector<Point>> surfaces = {
    // Station 1
    {{152.195, 99.1559, 500.266}, {165.039, 5.22261, 494.219}, {144.656, 1.04414, 516.03}, {131.811, 94.9774, 522.077}},
    // Station 2
    {{164.973, 99.1559, 529.408}, {177.817, 5.22261, 523.362}, {157.434, 1.04414, 545.073}, {144.589, 94.9774, 551.12}},
    // Station 3
    {{178.279, 99.1938, 557.082}, {191.123, 5.26059, 551.035}, {169.684, 1.00616, 571.584}, {156.84, 94.9394, 577.631}},
    // Station 4
    {{189.705, 109.082, 584.861}, {203.902, 5.26059, 578.178}, {182.462, 1.00616, 598.727}, {168.266, 104.827, 605.41}},
    // Station 5
    {{202.483, 109.082, 612.004}, {216.68, 5.26059, 605.321}, {195.24, 1.00616, 625.869}, {181.044, 104.827, 632.552}},
    // Station 6
    {{215.763, 109.114, 639.704}, {229.959, 5.29326, 633.021}, {207.518, 0.973489, 652.455}, {193.321, 104.794, 659.138}},
    // Station 7
    {{227.865, 114.058, 667.164}, {242.737, 5.29326, 660.162}, {220.296, 0.973489, 679.597}, {205.423, 109.738, 686.599}},
    // Station 8
    {{235.737, 113.079, 699.891}, {250.609, 4.31444, 692.889}, {237.98, 1.95231, 702.756}, {223.108, 110.717, 709.758}},
    // Station 9
    {{247.839, 118.023, 736.851}, {263.388, 4.31444, 729.531}, {250.758, 1.95231, 739.398}, {235.21, 115.661, 746.718}}};

/**
 * Calculate aperture half-width at a given z position
 * This function implements a z-dependent aperture that gets wider
 * as z increases (to account for beam spreading)
 * @param z Z-position in millimeters
 * @return Maximum allowed x-coordinate (aperture half-width) in millimeters
 */
double getXLimit(double z)
{
    // For positions before z=330mm, use fixed aperture of 200mm
    if (z < 330)
        return 200;
    else
        // Beyond z=330mm, aperture grows linearly
        // Equation: x = 0.4054*z + 66.22
        // Derived from points: (330mm,200mm) and (700mm,350mm)
        return 0.4054 * z + 66.22;
}

/**
 * Interpolate boundary position at a given z coordinate
 * For any z-position, finds the corresponding detector surface
 * and calculates the boundary by interpolation between corners
 * @param z Z-position in meters (will be converted to mm internally)
 * @return Point object containing boundary coordinate in meters
 */
Point interpolate(double z)
{
    // Convert input z from meters to millimeters (internal calculations use mm)
    z *= 100;

    // Default boundary value when z is outside defined surfaces range
    // Uses the getXLimit function for x and fixed y-limit of 155mm
    Point default_value(getXLimit(z) / 100, 155.0 / 100, z / 100);

    // Debug flag (set to true to enable verbose output)
    bool verbosityhere = 0;

    if (verbosityhere)
    {
        std::cout << "z = " << z << std::endl;
    }

    // Search through all defined detector surfaces
    for (const auto &surface : surfaces)
    {
        // Find min/max z-values for this surface
        double z_min = std::min({surface[0].z, surface[1].z, surface[2].z, surface[3].z});
        double z_max = std::max({surface[0].z, surface[1].z, surface[2].z, surface[3].z});

        // Check if requested z is within this surface's z-range
        if (z >= z_min && z <= z_max)
        {
            if (verbosityhere)
            {
                // Output surface index and range for debugging
                std::cout << "surface number " << &surface - &surfaces[0] << " with z_min = " << z_min << ", z_max = " << z_max << std::endl;
            }

            // Perform bilinear interpolation within the quadrilateral surface
            // First calculate interpolation parameters along two edges
            double t1 = (z - surface[0].z) / (surface[1].z - surface[0].z);
            double t2 = (z - surface[2].z) / (surface[3].z - surface[2].z);

            // Interpolate x,y along first edge (based on z-position)
            double x1 = surface[0].x + t1 * (surface[1].x - surface[0].x);
            double y1 = surface[0].y + t1 * (surface[1].y - surface[0].y);

            // Interpolate x,y along second edge (based on z-position)
            double x2 = surface[2].x + t2 * (surface[3].x - surface[2].x);
            double y2 = surface[2].y + t2 * (surface[3].y - surface[2].y);

            // Now interpolate between the two edge points
            double t = (z - surface[0].z) / (surface[2].z - surface[0].z);
            double x = x1 + t * (x2 - x1);
            double y = y1 + t * (y2 - y1);

            if (verbosityhere)
                std::cout << "x = " << x << ", y = " << y << ", z = " << z << std::endl;

            // Return interpolated position, converting back to meters
            return Point(x / 100, y / 100, z / 100);
        }
    }

    // Warning if z is outside the range of all defined surfaces
    if (verbosityhere)
        std::cout << "Warning: z = " << z / 100 << " is not covered by any surface" << std::endl;

    // Return default boundary if no matching surface is found
    return default_value;
}

/**
 * MagneticField class
 * Reads magnetic field map from file and provides interpolation
 * for field values at arbitrary positions
 */
class MagneticField
{
private:
    // Grid coordinates for the magnetic field map
    std::vector<double> x_grid, y_grid, z_grid;

    // 3D array to store field vectors at each grid point
    std::vector<std::vector<std::vector<TVector3>>> field_values;

    // Field map boundaries and grid spacing
    double x_min, x_max, y_min, y_max, z_min, z_max;
    double dx, dy, dz;

public:
    /**
     * Constructor: Load magnetic field from ASCII file
     * @param filename Path to the field map file
     * Format: First line: nx ny nz (grid dimensions)
     *         Second line: xmin xmax ymin ymax zmin zmax (field boundaries)
     *         Remaining lines: x y z Bx By Bz (field values at each point)
     */
    MagneticField(const std::string &filename)
    {
        // Open the magnetic field map file
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Could not open field map file: " << filename << std::endl;
            exit(1);
        }

        // Read header with grid information
        int nx, ny, nz; // Number of grid points in each dimension
        file >> nx >> ny >> nz;
        file >> x_min >> x_max >> y_min >> y_max >> z_min >> z_max;

        // Calculate grid spacing in each dimension
        dx = (x_max - x_min) / (nx - 1);
        dy = (y_max - y_min) / (ny - 1);
        dz = (z_max - z_min) / (nz - 1);

        // Initialize arrays for grid coordinates
        x_grid.resize(nx);
        y_grid.resize(ny);
        z_grid.resize(nz);

        // Fill grid coordinate arrays
        for (int i = 0; i < nx; i++)
            x_grid[i] = x_min + i * dx;
        for (int j = 0; j < ny; j++)
            y_grid[j] = y_min + j * dy;
        for (int k = 0; k < nz; k++)
            z_grid[k] = z_min + k * dz;

        // Initialize 3D array to hold field vectors
        // Uses nested vectors for flexibility
        field_values.resize(nx);
        for (int i = 0; i < nx; i++)
        {
            field_values[i].resize(ny);
            for (int j = 0; j < ny; j++)
            {
                field_values[i][j].resize(nz);
            }
        }

        // Read field values from file
        double x, y, z, Bx, By, Bz;
        int i, j, k;

        // Process each line in the file
        while (file >> x >> y >> z >> Bx >> By >> Bz)
        {
            // Convert coordinates to grid indices
            i = round((x - x_min) / dx);
            j = round((y - y_min) / dy);
            k = round((z - z_min) / dz);

            // Store field value if indices are valid
            if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
            {
                field_values[i][j][k] = TVector3(Bx, By, Bz);
            }
        }

        file.close();

        // Report successful loading
        std::cout << "Loaded magnetic field map with dimensions "
                  << nx << "x" << ny << "x" << nz << " from " << filename << std::endl;
    }
    /**
     * Trilinear interpolation method to get magnetic field at arbitrary position
     * @param position 3D position vector where we need the field value (in meters)
     * @return Interpolated magnetic field vector at requested position (in Tesla)
     */
    TVector3 getFieldAt(const TVector3 &position) const
    {
        double x = position.X();
        double y = position.Y();
        double z = position.Z();

        // Return zero field if position is outside the defined field region
        if (x < x_min || x > x_max || y < y_min || y > y_max || z < z_min || z > z_max)
        {
            return TVector3(0, 0, 0);
        }

        // Find grid cell containing this position
        // floor() ensures we get the lower index of the surrounding grid points
        int i0 = floor((x - x_min) / dx);
        int j0 = floor((y - y_min) / dy);
        int k0 = floor((z - z_min) / dz);

        // Safety check to ensure indices stay within array bounds
        // This prevents segmentation faults if position is near max boundary
        if (i0 >= static_cast<int>(x_grid.size() - 1))
            i0 = static_cast<int>(x_grid.size() - 2);
        if (j0 >= static_cast<int>(y_grid.size() - 1))
            j0 = static_cast<int>(y_grid.size() - 2);
        if (k0 >= static_cast<int>(z_grid.size() - 1))
            k0 = static_cast<int>(z_grid.size() - 2);

        // Calculate indices for upper corner of grid cell
        int i1 = i0 + 1;
        int j1 = j0 + 1;
        int k1 = k0 + 1;

        // Calculate interpolation weights - normalized distance within cell
        // These values determine how much each corner contributes to the interpolated result
        double wx = (x - x_grid[i0]) / dx; // 0 ≤ wx ≤ 1, weight for x direction
        double wy = (y - y_grid[j0]) / dy; // 0 ≤ wy ≤ 1, weight for y direction
        double wz = (z - z_grid[k0]) / dz; // 0 ≤ wz ≤ 1, weight for z direction

        // Perform trilinear interpolation using the 8 corners of the surrounding cell
        // For each corner, calculate its contribution based on its distance from position
        TVector3 B = (1 - wx) * (1 - wy) * (1 - wz) * field_values[i0][j0][k0] + // Weight for (i0,j0,k0)
                     (1 - wx) * (1 - wy) * wz * field_values[i0][j0][k1] +       // Weight for (i0,j0,k1)
                     (1 - wx) * wy * (1 - wz) * field_values[i0][j1][k0] +       // Weight for (i0,j1,k0)
                     (1 - wx) * wy * wz * field_values[i0][j1][k1] +             // Weight for (i0,j1,k1)
                     wx * (1 - wy) * (1 - wz) * field_values[i1][j0][k0] +       // Weight for (i1,j0,k0)
                     wx * (1 - wy) * wz * field_values[i1][j0][k1] +             // Weight for (i1,j0,k1)
                     wx * wy * (1 - wz) * field_values[i1][j1][k0] +             // Weight for (i1,j1,k0)
                     wx * wy * wz * field_values[i1][j1][k1];                    // Weight for (i1,j1,k1)

        return B; // Return interpolated magnetic field vector in Tesla
    }
};

/**
 * Calculate value of a Gaussian distribution at position x
 * Used for analytical field representation
 * @param x Position to evaluate
 * @param mu Mean of the Gaussian
 * @param sig Standard deviation of the Gaussian
 * @return Gaussian function value at x
 */
double gauss(double x, double mu, double sig)
{
    // Standard Gaussian distribution formula
    return exp(-pow((x - mu) / sig, 2.0) / 2.0) / (sqrt(2.0 * M_PI) * sig);
}

/**
 * Alternative analytical magnetic field function
 * Creates a parametrized field for testing when no field map is available
 * Field is modeled using Gaussian distributions in z direction
 * @param position 3D position vector where field is evaluated (in meters)
 * @return Magnetic field vector at requested position (in Tesla)
 */
TVector3 GetFieldValue(const TVector3 &position)
{
    TVector3 B(0.0, 0.0, 0.0); // Initialize field to zero

    // Set field to zero outside magnet aperture boundaries (±1800mm in x, ±2000mm in y)
    if (fabs(position.x() * 100) > 1800)
        B.SetX(0.);
    if (fabs(position.y() * 100) > 2000)
        B.SetY(0.);

    // Field model: Sum of two Gaussian distributions in the y-component only
    // First Gaussian: peak at z=431mm with width 82.9mm, amplitude -0.46530
    // Second Gaussian: peak at z=574mm with width 154mm, amplitude -0.864
    // Factor 1/(2π²) scales overall field strength
    B.SetY((-0.46530 * gauss(position.z() * 100, 431.0, 82.9) -
            0.864 * gauss(position.z() * 100, 574.0, 154.0)) /
           (2 * M_PI * M_PI));

    // No field in z-direction
    B.SetZ(0.);

    // Debug output (disabled by default for performance)
    std::cout << "Bx: " << B.X() << " By: " << B.Y() << " Bz: " << B.Z()
              << " at position: " << position.X() << ", " << position.Y() << ", " << position.Z() << std::endl;

    return B;
}

/**
 * Simplified constant magnetic field function
 * Creates a simple magnetic field with constant value beyond z=2m
 * Useful for basic testing of particle tracking
 * @param position 3D position vector where field is evaluated (in meters)
 * @return Magnetic field vector at requested position (in Tesla)
 */
TVector3 GetFieldValueConst(const TVector3 &position)
{
    TVector3 B(0.0, 0.0, 0.0); // Initialize field to zero

    // Simple model: constant -0.3T field in y-direction beyond z=2m (200mm)
    if (fabs(position.z() * 100) > 200)
        B.SetY(-0.3); // -0.3 Tesla in y-direction

    // No field in z-direction
    B.SetZ(0.);

    // Debug output (disabled by default for performance)
    std::cout << "Bx: " << B.X() << " By: " << B.Y() << " Bz: " << B.Z()
              << " at position: " << position.X() << ", " << position.Y() << ", " << position.Z() << std::endl;

    return B;
}

/**
 * Structure to store particle state at a single point
 * Represents a snapshot of position and momentum during trajectory
 */
struct State
{
    TVector3 position; // Position vector in meters
    TVector3 momentum; // Momentum vector in GeV/c
};

/**
 * Class for propagating charged particles through magnetic fields
 * Implements relativistic equations of motion with Runge-Kutta 4th order integration
 */
class ParticlePropagator
{
private:
    MagneticField &field; // Reference to magnetic field object
    double x_limit;       // Maximum allowed |x| position (in meters)
    double y_limit;       // Maximum allowed |y| position (in meters)
    double step_size;     // Integration step size (in meters)
    double charge;        // Particle charge (in units of elementary charge e)
    double mass;          // Particle rest mass (in GeV/c²)
    double c_light;       // Speed of light (in m/ns for time unit consistency)

    /**
     * Implementation of 4th-order Runge-Kutta method for solving particle's equation of motion
     * Updates position and momentum vectors for one time step
     * Implements fully relativistic treatment of charged particle in magnetic field
     *
     * @param position Current position vector (in meters) - updated in-place
     * @param momentum Current momentum vector (in GeV/c) - updated in-place
     * @param dt Time step for integration (in nanoseconds)
     */
    void RungeKutta4Step(TVector3 &position, TVector3 &momentum, double dt)
    {
        // q_factor represents particle charge for Lorentz force calculation
        const double q_factor = charge;

        // Calculate total relativistic energy: E² = (pc)² + (mc²)²
        double energy = sqrt(momentum.Mag2() + mass * mass); // Energy in GeV

        // ------ FIRST RK4 STEP (k1) ------
        // Calculate relativistic velocity vector: v = pc²/E
        // Divide by energy to get velocity scaled by c_light (pc/E = v/c)
        TVector3 vel1 = momentum * (c_light / energy); // Velocity in m/ns

        // Get magnetic field at current position
        TVector3 B1 = field.getFieldAt(position);

        // Calculate position change: dx = v·dt
        TVector3 k1_pos = vel1 * dt; // Position increment in meters

        // Calculate momentum change: dp = q(v×B)dt
        // Lorentz force equation: F = q(v×B), and F = dp/dt
        TVector3 k1_mom = q_factor * vel1.Cross(B1) * dt; // Momentum increment in GeV/c

        // ------ SECOND RK4 STEP (k2) ------
        // Use midpoint of first increment (position + k1/2)
        TVector3 pos2 = position + k1_pos * 0.5;
        TVector3 mom2 = momentum + k1_mom * 0.5;

        // Recalculate energy and velocity at this intermediate state
        double energy2 = sqrt(mom2.Mag2() + mass * mass);
        TVector3 vel2 = mom2 * (c_light / energy2);

        // Get field at new position
        TVector3 B2 = field.getFieldAt(pos2);

        // Calculate increments for second step
        TVector3 k2_pos = vel2 * dt;
        TVector3 k2_mom = q_factor * vel2.Cross(B2) * dt;

        // ------ THIRD RK4 STEP (k3) ------
        // Use midpoint of second increment
        TVector3 pos3 = position + k2_pos * 0.5;
        TVector3 mom3 = momentum + k2_mom * 0.5;

        // Recalculate energy and velocity
        double energy3 = sqrt(mom3.Mag2() + mass * mass);
        TVector3 vel3 = mom3 * (c_light / energy3);

        // Get field at new position
        TVector3 B3 = field.getFieldAt(pos3);

        // Calculate increments for third step
        TVector3 k3_pos = vel3 * dt;
        TVector3 k3_mom = q_factor * vel3.Cross(B3) * dt;

        // ------ FOURTH RK4 STEP (k4) ------
        // Use full step with third increment
        TVector3 pos4 = position + k3_pos;
        TVector3 mom4 = momentum + k3_mom;

        // Recalculate energy and velocity
        double energy4 = sqrt(mom4.Mag2() + mass * mass);
        TVector3 vel4 = mom4 * (c_light / energy4);

        // Get field at new position
        TVector3 B4 = field.getFieldAt(pos4);

        // Calculate increments for fourth step
        TVector3 k4_pos = vel4 * dt;
        TVector3 k4_mom = q_factor * vel4.Cross(B4) * dt;

        // Optional detailed diagnostic output (disabled by default)
        if (false) // Set to true to enable debugging
        {
            // Extensive diagnostic information about current state and RK4 increments
            // Shows position, momentum, energy, fields, etc. at each step

            std::cout << std::endl;
            std::cout << "RungeKutta4Step::Position: " << position.X() << ", " << position.Y() << ", " << position.Z() << std::endl;
            std::cout << "RungeKutta4Step::Momentum: " << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << std::endl;
            std::cout << "RungeKutta4Step::Energy: " << energy << " GeV, γ = " << energy / mass << std::endl;
            std::cout << "RungeKutta4Step::Velocity: " << vel1.Mag() << " m/ns (" << vel1.Mag() / c_light << " × c)" << std::endl;
            std::cout << "RungeKutta4Step::B1: " << B1.X() << ", " << B1.Y() << ", " << B1.Z() << " T" << std::endl;
            std::cout << "RungeKutta4Step::dt: " << dt << " ns" << std::endl;

            std::cout << "\tRungeKutta4Step::k1_pos: " << k1_pos.X() << ", " << k1_pos.Y() << ", " << k1_pos.Z()
                      << " m, scaling: " << c_light / energy << std::endl;
            std::cout << "\tRungeKutta4Step::k2_pos: " << k2_pos.X() << ", " << k2_pos.Y() << ", " << k2_pos.Z()
                      << " m, scaling: " << c_light / energy2 << std::endl;
            std::cout << "\tRungeKutta4Step::k3_pos: " << k3_pos.X() << ", " << k3_pos.Y() << ", " << k3_pos.Z()
                      << " m, scaling: " << c_light / energy3 << std::endl;
            std::cout << "\tRungeKutta4Step::k4_pos: " << k4_pos.X() << ", " << k4_pos.Y() << ", " << k4_pos.Z()
                      << " m, scaling: " << c_light / energy4 << std::endl;

            std::cout << "\tRungeKutta4Step::k1_mom: " << k1_mom.X() << ", " << k1_mom.Y() << ", " << k1_mom.Z() << " GeV/c" << std::endl;
            std::cout << "\tRungeKutta4Step::k2_mom: " << k2_mom.X() << ", " << k2_mom.Y() << ", " << k2_mom.Z() << " GeV/c" << std::endl;
            std::cout << "\tRungeKutta4Step::k3_mom: " << k3_mom.X() << ", " << k3_mom.Y() << ", " << k3_mom.Z() << " GeV/c" << std::endl;
            std::cout << "\tRungeKutta4Step::k4_mom: " << k4_mom.X() << ", " << k4_mom.Y() << ", " << k4_mom.Z() << " GeV/c" << std::endl;
        }

        // ------ FINAL RK4 UPDATE ------
        // Combine all 4 increments with standard RK4 weights (1:2:2:1)/6
        // This is the key formula for 4th-order accuracy
        position += (k1_pos + 2.0 * k2_pos + 2.0 * k3_pos + k4_pos) * (1.0 / 6.0);
        momentum += (k1_mom + 2.0 * k2_mom + 2.0 * k3_mom + k4_mom) * (1.0 / 6.0);

        // Optional energy conservation check (disabled by default)
        if (false)
        {
            // Calculate updated energy to verify conservation
            double new_energy = sqrt(momentum.Mag2() + mass * mass);
            // Show energy before and after update to check conservation

            std::cout << "\tRungeKutta4Step::New Position: " << position.X() << ", " << position.Y() << ", " << position.Z() << " m" << std::endl;
            std::cout << "\tRungeKutta4Step::New Momentum: " << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << " GeV/c" << std::endl;
            std::cout << "\tRungeKutta4Step::New Energy: " << new_energy << " GeV (Change: " << new_energy - energy << " GeV)" << std::endl;
            std::cout << "\tRungeKutta4Step::Momentum magnitude: " << momentum.Mag() << " GeV/c" << std::endl;
        }
    }

public:
    /**
     * Constructor for the ParticlePropagator class
     * @param field Reference to MagneticField object providing B field values
     * @param step_size_m Step size for numerical integration in meters (default: 0.01m)
     * @param x_lim Maximum allowed x-coordinate in meters (default: 1.0m)
     * @param y_lim Maximum allowed y-coordinate in meters (default: 1.0m)
     */
    ParticlePropagator(MagneticField &field, double step_size_m = 0.01,
                       double x_lim = 1.0, double y_lim = 1.0)
        : field(field), step_size(step_size_m), charge(1.0), mass(0.1057),
          x_limit(x_lim), y_limit(y_lim)
    {
        // Speed of light in m/ns for appropriate time units in integration
        c_light = 0.299792458;
        // Alternative value for speed of light in m/s (commented out)
        // c_light = 299792458.0;
    }

    /**
     * Set physical properties of the particle being propagated
     * @param q Charge of the particle in units of elementary charge e
     * @param m Rest mass of the particle in GeV/c²
     */
    void setParticleProperties(double q, double m)
    {
        charge = q; // in units of elementary charge
        mass = m;   // in GeV/c^2
    }

    /**
     * Set geometric limits for propagation (magnet aperture)
     * @param x_lim Maximum allowed |x| position in meters
     * @param y_lim Maximum allowed |y| position in meters
     */
    void setLimits(double x_lim, double y_lim)
    {
        x_limit = x_lim;
        y_limit = y_lim;
    }

    /**
     * Set the integration step size
     * @param step Step size in meters for numerical integration
     */
    void setStepSize(double step)
    {
        step_size = step; // in meters
    }

    /**
     * Main propagation method - tracks particle through magnetic field
     * @param initial_position Starting position vector (in meters)
     * @param initial_momentum Starting momentum vector (in GeV/c)
     * @param max_length Maximum path length to propagate (in meters)
     * @param max_steps Maximum number of integration steps to perform
     * @return Vector of State objects containing the full trajectory
     */
    std::vector<State> propagate(const TVector3 &initial_position,
                                 const TVector3 &initial_momentum,
                                 double max_length = 8.0,
                                 int max_steps = 10000)
    {
        // Create trajectory vector and add initial state
        std::vector<State> trajectory;
        trajectory.push_back({initial_position, initial_momentum});

        // Set current position and momentum to initial values
        TVector3 position = initial_position;
        TVector3 momentum = initial_momentum;

        // Debug output can be enabled by uncommenting these lines
        // std::cout << "propagate::Initial Position: " << position.X() << ", " << position.Y() << ", " << position.Z() << std::endl;
        // std::cout << "propagate::Initial Momentum: " << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << std::endl;

        // Track total path length and boundary hit status
        double total_length = 0.0;
        bool hit_boundary = false;

        // Main propagation loop
        for (int step = 0; step < max_steps && !hit_boundary; step++)
        {
            // Perform one RK4 integration step
            RungeKutta4Step(position, momentum, step_size);

            // Get boundary limits at current z-position (allows for varying aperture)
            Point interpolatedBoundary = interpolate(position.Z());

            // Check if particle hit magnet boundary
            // Original simple check: if (fabs(position.X()) >= x_limit || fabs(position.Y()) >= y_limit)
            if (fabs(position.X()) >= fabs(interpolatedBoundary.x) || fabs(position.Y()) >= fabs(interpolatedBoundary.y))
            {
                hit_boundary = true;
                // Debug output when hitting boundary (disabled by default)
                if (false)
                {
                    std::cout << "propagate::Particle hit boundary at position: ("
                              << position.X() << ", " << position.Y() << ", "
                              << position.Z() << ")" << std::endl;
                }
            }

            // Add current state to trajectory
            trajectory.push_back({position, momentum});

            // Calculate incremental path length
            if (trajectory.size() > 1)
            {
                total_length += (trajectory[trajectory.size() - 1].position -
                                 trajectory[trajectory.size() - 2].position)
                                    .Mag();
            }

            // Stop if we've reached maximum path length
            if (total_length >= max_length)
            {
                break;
            }
        }

        return trajectory;
    }

    /**
     * Overloaded propagation method accepting 4-momentum as input
     * Useful when working with relativistic particles directly from detectors/generators
     * @param initial_position Starting position vector (in meters)
     * @param four_momentum 4-momentum vector with E and p components (in GeV/c)
     * @param max_length Maximum path length to propagate (in meters)
     * @param max_steps Maximum number of integration steps to perform
     * @return Vector of State objects containing the full trajectory
     */
    std::vector<State> propagate(const TVector3 &initial_position,
                                 const TLorentzVector &four_momentum,
                                 double max_length = 10.0,
                                 int max_steps = 10000)
    {
        // Extract 3-momentum components from 4-momentum
        TVector3 momentum = four_momentum.Vect();

        // Extract mass from 4-momentum: m² = E² - p²
        mass = four_momentum.M();

        // Call the main propagate method with extracted 3-momentum
        return propagate(initial_position, momentum, max_length, max_steps);
    }
};

/**
 * Main function to demonstrate particle propagation and visualization
 * Simulates multiple particles through a magnetic field and generates output files
 * 
 * @param inputFieldASCII Path to ASCII file containing magnetic field map
 * @param nParticles Number of particles to simulate (with different momenta)
 * @param verboseInfo Whether to print detailed information during propagation
 * @return 0 on success, 1 on error
 */
int track_propagatorRK(TString inputFieldASCII = "field_map_lhcb.txt", int nParticles = 200, bool verboseInfo = false)
{
    // Validate input field file path
    if (inputFieldASCII == "")
    {
        std::cerr << "Error: No magnetic field file provided." << std::endl;
        std::cerr << "Usage: ./track_propagator-RK field_map.txt" << std::endl;
        return 1;
    }

    // Convert TString to std::string for C++ file handling
    std::string field_file = inputFieldASCII.Data();

    // Load magnetic field from file
    MagneticField field(field_file);

    // Create particle propagator with specified magnet aperture
    double x_limit = 1.5;  // 1.5 m magnet half-width
    double y_limit = 2.0;  // 2.0 m magnet half-height
    ParticlePropagator propagator(field, 0.01, x_limit, y_limit); // 0.01m step size

    // Create ROOT output file for all particle trajectories
    TFile *rootFile = new TFile("propagation_multiparticle.root", "RECREATE");

    // Create a tree to store all trajectories
    TTree *tree = new TTree("trajectories", "Multiple Particle Trajectories");

    // Variables to store in the tree
    int particle_id;       // Unique ID for each particle
    double x, y, z;        // Position coordinates (m)
    double px_out, py_out, pz_out;  // Momentum components (GeV/c)
    double energy;         // Total energy (GeV)
    double mass_value;     // Particle mass (GeV/c²)
    double charge_value;   // Particle charge (e)
    int step_num;          // Step number along trajectory
    bool hit_boundary;     // Flag indicating boundary hit

    // Create branches for all variables in the tree
    tree->Branch("particle_id", &particle_id, "particle_id/I");
    tree->Branch("step", &step_num, "step/I");
    tree->Branch("x", &x, "x/D");
    tree->Branch("y", &y, "y/D");
    tree->Branch("z", &z, "z/D");
    tree->Branch("px", &px_out, "px/D");
    tree->Branch("py", &py_out, "py/D");
    tree->Branch("pz", &pz_out, "pz/D");
    tree->Branch("energy", &energy, "energy/D");
    tree->Branch("mass", &mass_value, "mass/D");
    tree->Branch("charge", &charge_value, "charge/D");
    tree->Branch("hit_boundary", &hit_boundary, "hit_boundary/O");

    // Create ASCII output file as alternative format
    std::ofstream outfile("trajectories.dat");
    outfile << "# particle_id step x y z px py pz energy mass charge hit_boundary" << std::endl;

    // Create MultiGraph objects for combined trajectory visualization
    // These will contain all particle trajectories in different projections
    TMultiGraph *mg_xy = new TMultiGraph("mg_xy", "All Trajectories XY Projection;X (m);Y (m)");
    TMultiGraph *mg_xz = new TMultiGraph("mg_xz", "All Trajectories XZ Projection;X (m);Z (m)");
    TMultiGraph *mg_yz = new TMultiGraph("mg_yz", "All Trajectories YZ Projection;Y (m);Z (m)");

    // Define colors for different particles (cycling through these for more particles)
    int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet, kPink, kTeal, kYellow};

    // Loop through all particles to simulate
    for (int i = 0; i < nParticles; i++)
    {
        // Set standard particle properties (all pions in this case)
        double charge = 1.0;       // in units of e
        double mass = 0.13957039;  // pion mass in GeV/c²
        std::string particle_type = "pion";

        // Alternative: Vary particle types by ID (commented out)
        // switch(i % 3) {
        //     case 0: mass = 0.000511; particle_type = "electron"; break;
        //     case 1: mass = 0.13957039; particle_type = "pion"; break;
        //     case 2: mass = 0.1057; particle_type = "muon"; break;
        //     default: mass = 0.13957039; particle_type = "pion";
        // }

        // Set particle properties in the propagator
        propagator.setParticleProperties(charge, mass);
        mass_value = mass;
        charge_value = charge;

        // Create logarithmically spaced momentum values from 0.5 to 5 GeV/c
        // This ensures good coverage of both low and high momentum behaviors
        double log_min = -0.3;  // log10(0.5) ≈ -0.3
        double log_max = 0.7;   // log10(5) ≈ 0.7
        double momentum_magnitude = pow(10.0, log_min + i * ((log_max - log_min) / (nParticles - 1)));

        // Set momentum direction (mainly along z-axis with small x component)
        double angle = 0.0;  // No angle variation in this example
        // Alternative: double angle = (i * 45.0 / nParticles) * (M_PI / 180.0);
        double px = 0.1 * momentum_magnitude * sin(angle);  // Small x component (10% of total)
        double py = 0.0;                                    // No y component
        double pz = momentum_magnitude * cos(angle);        // Main momentum along z

        // Set initial position at origin (no variation in this example)
        double x_offset = 0.0;
        double y_offset = 0.0;
        // Alternative: Add position variation
        // double x_offset = (i % 3) * 0.05 - 0.05;
        // double y_offset = ((i / 3) % 3) * 0.05 - 0.05;
        
        TVector3 initial_position(x_offset, y_offset, 0.0);

        // Calculate total energy from momentum and mass: E² = p²c² + m²c⁴
        double total_energy = sqrt(px * px + py * py + pz * pz + mass * mass);
        TLorentzVector four_momentum(px, py, pz, total_energy);

        // Store initial momentum for later comparison
        double orig_px = px;
        double orig_py = py;
        double orig_pz = pz;

        // Print initial conditions if verbose mode enabled
        if (verboseInfo)
        {
            std::cout << "\n=======================================" << std::endl;
            std::cout << "Propagating particle " << i << " (" << particle_type << ")" << std::endl;
            std::cout << "Mass: " << mass << " GeV/c²" << std::endl;
            std::cout << "Initial momentum: (" << px << ", " << py << ", " << pz
                      << ") GeV/c (magnitude: " << momentum_magnitude << " GeV/c)" << std::endl;
            std::cout << "Initial position: (" << initial_position.X() << ", "
                      << initial_position.Y() << ", " << initial_position.Z() << ") m" << std::endl;
        }
        
        // Propagate particle through the magnetic field (max length 8m, up to 50k steps)
        std::vector<State> trajectory = propagator.propagate(
            initial_position, four_momentum, 8.0, 50000);

        // Create graphs for visualizing this particle's trajectory
        TGraph *graph_xy = new TGraph(trajectory.size());  // XY projection
        TGraph *graph_xz = new TGraph(trajectory.size());  // XZ projection
        TGraph *graph_yz = new TGraph(trajectory.size());  // YZ projection

        // Set visual styles for trajectory plots
        graph_xy->SetLineColor(colors[i % 10]);  // Cycle through colors
        graph_xz->SetLineColor(colors[i % 10]);
        graph_yz->SetLineColor(colors[i % 10]);
        graph_xy->SetLineWidth(2);
        graph_xz->SetLineWidth(2);
        graph_yz->SetLineWidth(2);

        // Process all points along the trajectory
        particle_id = i;
        hit_boundary = false;

        for (size_t j = 0; j < trajectory.size(); j++)
        {
            // Extract position and momentum from trajectory point
            x = trajectory[j].position.X();
            y = trajectory[j].position.Y();
            z = trajectory[j].position.Z();
            px_out = trajectory[j].momentum.X();
            py_out = trajectory[j].momentum.Y();
            pz_out = trajectory[j].momentum.Z();
            
            // Calculate total energy: E² = p²c² + m²c⁴
            energy = sqrt(px_out * px_out + py_out * py_out + pz_out * pz_out + mass * mass);
            step_num = j;

            // Check if particle hit boundary
            if (fabs(x) >= x_limit || fabs(y) >= y_limit)
            {
                hit_boundary = true;
            }

            // Add point to visualization graphs
            graph_xy->SetPoint(j, x, y);
            graph_xz->SetPoint(j, x, z);
            graph_yz->SetPoint(j, y, z);

            // Fill tree with this trajectory point
            tree->Fill();

            // Write point to ASCII file
            outfile << particle_id << " " << step_num << " "
                    << x << " " << y << " " << z << " "
                    << px_out << " " << py_out << " " << pz_out << " "
                    << energy << " " << mass << " " << charge << " "
                    << (hit_boundary ? 1 : 0) << std::endl;
        }

        // Add this particle's graphs to the multi-graphs for combined visualization
        mg_xy->Add(graph_xy, "l");  // "l" option draws lines
        mg_xz->Add(graph_xz, "l");
        mg_yz->Add(graph_yz, "l");

        // Print summary statistics if verbose mode enabled
        if (verboseInfo)
        {
            std::cout << "Final position: (" << trajectory.back().position.X() << ", "
                      << trajectory.back().position.Y() << ", " << trajectory.back().position.Z() << ") m" << std::endl;
            std::cout << "Final momentum: (" << trajectory.back().momentum.X() << ", "
                      << trajectory.back().momentum.Y() << ", " << trajectory.back().momentum.Z()
                      << ") GeV/c" << std::endl;
            std::cout << "Momentum conservation: "
                      << 100.0 * (trajectory.back().momentum.Mag() / momentum_magnitude) << "%" << std::endl;
            std::cout << "Points in trajectory: " << trajectory.size() << std::endl;
        }
        
        // Option to create individual trajectory plots for each particle
        // (Commented out to avoid creating too many files)
        // TCanvas *c_particle = new TCanvas(Form("c_particle_%d", i),
        //                                  Form("Particle %d (%s) Trajectory", i, particle_type.c_str()),
        //                                  1200, 400);
        // c_particle->Divide(3, 1);
        // c_particle->cd(1);
        // graph_xy->SetTitle(Form("Particle %d (%s) XY Projection;X (m);Y (m)", i, particle_type.c_str()));
        // graph_xy->Draw("AL");
        // c_particle->cd(2);
        // graph_xz->SetTitle(Form("Particle %d (%s) XZ Projection;X (m);Z (m)", i, particle_type.c_str()));
        // graph_xz->Draw("AL");
        // c_particle->cd(3);
        // graph_yz->SetTitle(Form("Particle %d (%s) YZ Projection;Y (m);Z (m)", i, particle_type.c_str()));
        // graph_yz->Draw("AL");
        // c_particle->SaveAs(Form("c_particle_%d.pdf", i));
        // c_particle->Write();
    }

    // Create and save XY projection plot with all trajectories
    TCanvas *c_all_xy = new TCanvas("c_all_xy", "All Trajectories XY Projection", 800, 600);
    mg_xy->Draw("A");  // "A" draws axes and all trajectories

    // Add rectangular boundary box showing magnet aperture in XY view
    TBox *boundary_xy = new TBox(-x_limit, -y_limit, x_limit, y_limit);
    boundary_xy->SetFillStyle(0);  // Transparent fill
    boundary_xy->SetLineColor(kBlack);
    boundary_xy->SetLineWidth(2);
    boundary_xy->SetLineStyle(2);  // Dashed line
    boundary_xy->Draw("same");
    
    // Save plot to PDF file
    c_all_xy->SaveAs("c_all_xy.pdf");
    c_all_xy->Write();

    // Create and save XZ projection plot with all trajectories
    TCanvas *c_all_xz = new TCanvas("c_all_xz", "All Trajectories XZ Projection", 600, 1200);
    mg_xz->Draw("A");
    c_all_xz->SaveAs("c_all_xz.pdf");
    c_all_xz->Write();

    // Create and save YZ projection plot with all trajectories
    TCanvas *c_all_yz = new TCanvas("c_all_yz", "All Trajectories YZ Projection", 800, 600);
    mg_yz->Draw("A");

    // Add boundary lines showing magnet aperture in YZ view
    TLine *boundary_yz_bottom = new TLine(-y_limit, 0, -y_limit, 8);
    TLine *boundary_yz_top = new TLine(y_limit, 0, y_limit, 8);
    boundary_yz_bottom->SetLineColor(kBlack);
    boundary_yz_top->SetLineColor(kBlack);
    boundary_yz_bottom->SetLineWidth(2);
    boundary_yz_top->SetLineWidth(2);
    boundary_yz_bottom->SetLineStyle(2);
    boundary_yz_top->SetLineStyle(2);
    boundary_yz_bottom->Draw();
    boundary_yz_top->Draw();

    c_all_yz->Write();

    // Write tree to ROOT file and close all files
    tree->Write();
    rootFile->Close();
    outfile.close();

    // Print completion message with summary
    std::cout << "\n=======================================" << std::endl;
    std::cout << "Completed propagation of " << nParticles << " particles" << std::endl;
    std::cout << "Results saved to trajectories.dat and propagation_multiparticle.root" << std::endl;

    return 0;
}