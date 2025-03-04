/*
I'll create a C++ program using ROOT classes to propagate a charged particle through a tabulated magnetic field. This requires implementing several components: magnetic field interpolation, equation of motion solver, and particle tracking.

This code implements a particle propagation system in a tabulated magnetic field using ROOT classes. The implementation includes:

### Key Components:

1. **MagneticField Class**:
   - Reads magnetic field data from an ASCII file
   - Performs trilinear interpolation to determine field at any point
   - Handles field boundaries appropriately

2. **ParticlePropagator Class**:
   - Implements 4th-order Runge-Kutta method for solving the equations of motion
   - Handles relativistic effects correctly
   - Tracks particle trajectory through the field
   - Supports both 3-momentum and 4-momentum inputs

3. **Main Program**:
   - Demonstrates usage with a sample particle (muon by default)
   - Creates visualizations of the trajectory in ROOT
   - Saves trajectory data to both ASCII and ROOT file formats

### Usage:

The program expects the magnetic field file path as a command-line argument:

```bash
./particle_propagator field_map.txt
```

The magnetic field file should have the format:
```
nx ny nz
xmin xmax ymin ymax zmin zmax
x1 y1 z1 Bx1 By1 Bz1
x2 y2 z2 Bx2 By2 Bz2
...
```

The program outputs:
1. A text file (`trajectory.dat`) with trajectory points
2. A ROOT file (`propagation.root`) containing:
   - A TTree with trajectory data
   - Three projection plots (XY, XZ, YZ) of the trajectory

You can modify the initial particle parameters (position, momentum, charge, mass) in the main function to propagate different particles.

 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TBox.h"
#include "TLine.h"
#include "TLegend.h"


class Point {
    public:
        double x, y, z, t;
        Point(double x, double y, double z, double t = 0.0) : x(x), y(y), z(z), t(t) {}
    };

std::vector<std::vector<Point>> surfaces = {
    {{152.195, 99.1559, 500.266}, {165.039, 5.22261, 494.219}, {144.656, 1.04414, 516.03}, {131.811, 94.9774, 522.077}},
    {{164.973, 99.1559, 529.408}, {177.817, 5.22261, 523.362}, {157.434, 1.04414, 545.073}, {144.589, 94.9774, 551.12}},
    {{178.279, 99.1938, 557.082}, {191.123, 5.26059, 551.035}, {169.684, 1.00616, 571.584}, {156.84, 94.9394, 577.631}},
    {{189.705, 109.082, 584.861}, {203.902, 5.26059, 578.178}, {182.462, 1.00616, 598.727}, {168.266, 104.827, 605.41}},
    {{202.483, 109.082, 612.004}, {216.68, 5.26059, 605.321}, {195.24, 1.00616, 625.869}, {181.044, 104.827, 632.552}},
    {{215.763, 109.114, 639.704}, {229.959, 5.29326, 633.021}, {207.518, 0.973489, 652.455}, {193.321, 104.794, 659.138}},
    {{227.865, 114.058, 667.164}, {242.737, 5.29326, 660.162}, {220.296, 0.973489, 679.597}, {205.423, 109.738, 686.599}},
    {{235.737, 113.079, 699.891}, {250.609, 4.31444, 692.889}, {237.98, 1.95231, 702.756}, {223.108, 110.717, 709.758}},
    {{247.839, 118.023, 736.851}, {263.388, 4.31444, 729.531}, {250.758, 1.95231, 739.398}, {235.21, 115.661, 746.718}}
};

double getXLimit(double z) {
    if(z<495) return 150;
    else return 0.043*z + 237.894 + 50; //TODO added 50 to move gap further out
}
Point interpolate(double z) {
    z*=100;
    // Default value if z is not covered by any surface
    Point default_value(getXLimit(z)/100, 155.0/100, z/100);
    bool verbosityhere = 0;
    if(verbosityhere){
        std::cout << "z = " << z << std::endl;
    }
    for (const auto& surface : surfaces) {
        double z_min = std::min({surface[0].z, surface[1].z, surface[2].z, surface[3].z});
        double z_max = std::max({surface[0].z, surface[1].z, surface[2].z, surface[3].z});


        if (z >= z_min && z <= z_max) {
            if(verbosityhere){
                std::cout << "surface number " << &surface - &surfaces[0] << " with z_min = " << z_min << ", z_max = " << z_max << std::endl;
                
            }
            // Perform bilinear interpolation within the surface
            double t1 = (z - surface[0].z) / (surface[1].z - surface[0].z);
            double t2 = (z - surface[2].z) / (surface[3].z - surface[2].z);

            double x1 = surface[0].x + t1 * (surface[1].x - surface[0].x);
            double y1 = surface[0].y + t1 * (surface[1].y - surface[0].y);

            double x2 = surface[2].x + t2 * (surface[3].x - surface[2].x);
            double y2 = surface[2].y + t2 * (surface[3].y - surface[2].y);

            double t = (z - surface[0].z) / (surface[2].z - surface[0].z);

            double x = x1 + t * (x2 - x1);
            double y = y1 + t * (y2 - y1);
            if(verbosityhere) std::cout << "x = " << x << ", y = " << y << ", z = " << z << std::endl;
            return Point(x/100, y/100, z/100);
        }
    }
    if(verbosityhere)std::cout << "Warning: z = " << z/100 << " is not covered by any surface" << std::endl;
    return default_value;
}

// Class for magnetic field interpolation from table
class MagneticField
{
private:
    std::vector<double> x_grid, y_grid, z_grid;                   // Grid coordinates
    std::vector<std::vector<std::vector<TVector3>>> field_values; // B field values at grid points
    double x_min, x_max, y_min, y_max, z_min, z_max;              // Field boundaries
    double dx, dy, dz;                                            // Grid spacing

public:
    MagneticField(const std::string &filename)
    {
        // Read magnetic field data from ASCII file
        std::ifstream file(filename);
        if (!file.is_open())
        {
            std::cerr << "Error: Could not open field map file: " << filename << std::endl;
            exit(1);
        }

        // Read header with grid information
        // Format: nx ny nz xmin xmax ymin ymax zmin zmax
        int nx, ny, nz;
        file >> nx >> ny >> nz;
        file >> x_min >> x_max >> y_min >> y_max >> z_min >> z_max;

        // Calculate grid spacing
        dx = (x_max - x_min) / (nx - 1);
        dy = (y_max - y_min) / (ny - 1);
        dz = (z_max - z_min) / (nz - 1);

        // Reserve space for grid coordinates
        x_grid.resize(nx);
        y_grid.resize(ny);
        z_grid.resize(nz);

        // Initialize grid coordinates
        for (int i = 0; i < nx; i++)
            x_grid[i] = x_min + i * dx;
        for (int j = 0; j < ny; j++)
            y_grid[j] = y_min + j * dy;
        for (int k = 0; k < nz; k++)
            z_grid[k] = z_min + k * dz;

        // Resize field values array
        field_values.resize(nx);
        for (int i = 0; i < nx; i++)
        {
            field_values[i].resize(ny);
            for (int j = 0; j < ny; j++)
            {
                field_values[i][j].resize(nz);
            }
        }

        // Read field values
        // Format: x y z Bx By Bz
        double x, y, z, Bx, By, Bz;
        int i, j, k;

        // Read every field value from the file
        while (file >> x >> y >> z >> Bx >> By >> Bz)
        {
            // Find indices
            i = round((x - x_min) / dx);
            j = round((y - y_min) / dy);
            k = round((z - z_min) / dz);

            // Check bounds
            if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
            {
                field_values[i][j][k] = TVector3(Bx, By, Bz);
            }
        }

        file.close();

        std::cout << "Loaded magnetic field map with dimensions "
                  << nx << "x" << ny << "x" << nz << " from " << filename << std::endl;
    }

    // Trilinear interpolation for field value at arbitrary position
    TVector3 getFieldAt(const TVector3 &position) const
    {
        double x = position.X();
        double y = position.Y();
        double z = position.Z();

        // Check if position is within field boundaries
        if (x < x_min || x > x_max || y < y_min || y > y_max || z < z_min || z > z_max)
        {
            // debug output
            // std::cout << "Position outside field boundaries: " << x << ", " << y << ", " << z << std::endl;
            // Return zero field outside boundaries
            return TVector3(0, 0, 0);
        }

        // Find grid indices and weights for interpolation
        int i0 = floor((x - x_min) / dx);
        int j0 = floor((y - y_min) / dy);
        int k0 = floor((z - z_min) / dz);
        
        //debug output
        // std::cout << "i0: " << i0 << " j0: " << j0 << " k0: " << k0 << std::endl;

        // Ensure we don't exceed array bounds
        if (i0 >= (x_grid.size() - 1))
            i0 = (x_grid.size() - 2);
        if (j0 >= (y_grid.size() - 1))
            j0 = (y_grid.size() - 2);
        if (k0 >= (z_grid.size() - 1))
            k0 = (z_grid.size() - 2);

        int i1 = i0 + 1;
        int j1 = j0 + 1;
        int k1 = k0 + 1;

        // Calculate interpolation weights
        double wx = (x - x_grid[i0]) / dx;
        double wy = (y - y_grid[j0]) / dy;
        double wz = (z - z_grid[k0]) / dz;

        // Trilinear interpolation
        TVector3 B = (1 - wx) * (1 - wy) * (1 - wz) * field_values[i0][j0][k0] +
                     (1 - wx) * (1 - wy) * wz * field_values[i0][j0][k1] +
                     (1 - wx) * wy * (1 - wz) * field_values[i0][j1][k0] +
                     (1 - wx) * wy * wz * field_values[i0][j1][k1] +
                     wx * (1 - wy) * (1 - wz) * field_values[i1][j0][k0] +
                     wx * (1 - wy) * wz * field_values[i1][j0][k1] +
                     wx * wy * (1 - wz) * field_values[i1][j1][k0] +
                     wx * wy * wz * field_values[i1][j1][k1];

        //debug output
        // std::cout << "Bx: " << B.X() << " By: " << B.Y() << " Bz: " << B.Z() << std::endl;

        return B;
    }
};
double gauss(double x, double mu, double sig) {
    return exp(-pow((x - mu) / sig, 2.0) / 2.0) / (sqrt(2.0 * M_PI) * sig);
}
TVector3 GetFieldValue(const TVector3 &position)
{
    TVector3 B(0.0, 0.0, 0.0);
    if (fabs(position.x()*100) > 1800) B.SetX(0.);
    if (fabs(position.y()*100) > 2000) B.SetY(0.);
    B.SetY((-0.46530 * gauss(position.z()*100, 431.0, 82.9) - 0.864 * gauss(position.z()*100, 574.0, 154.0)) / (2 * M_PI * M_PI));
    B.SetZ(0.);
    //debug output
    std::cout << "Bx: " << B.X() << " By: " << B.Y() << " Bz: " << B.Z() << " at position: " << position.X() << ", " << position.Y() << ", " << position.Z() << std::endl;
    return B;
}
TVector3 GetFieldValueConst(const TVector3 &position)
{
    TVector3 B(0.0, 0.0, 0.0);
    // if (fabs(position.x()*100) > 1800) B.SetX(0.);
    // if (fabs(position.y()*100) > 2000) B.SetY(0.);
    if (fabs(position.z()*100) > 200)B.SetY(-0.3);
    B.SetZ(0.);
    //debug output
    std::cout << "Bx: " << B.X() << " By: " << B.Y() << " Bz: " << B.Z() << " at position: " << position.X() << ", " << position.Y() << ", " << position.Z() << std::endl;
    return B;
}

// Structure to hold position and momentum
struct State {
    TVector3 position;
    TVector3 momentum;
};

// Class for propagating charged particles through magnetic field
class ParticlePropagator
{
private:
    MagneticField &field;
    double x_limit;   // Maximum absolute x coordinate
    double y_limit;   // Maximum absolute y coordinate
    double step_size; // Integration step size in meters
    double charge;    // Particle charge in units of e
    double mass;      // Particle mass in GeV/c^2
    double c_light;   // Speed of light in m/s

    // Runge-Kutta 4th order method for solving ODEs
    void RungeKutta4Step(TVector3 &position, TVector3 &momentum, double dt)
{
    // c_light is in m/ns, charge in units of e
    const double q_factor = charge;
    
    // Calculate energy and velocity with proper relativistic treatment
    double energy = sqrt(momentum.Mag2() + mass * mass);  // Total energy in GeV
    
    // Calculate k1
    TVector3 vel1 = momentum * (c_light / energy);  // v = pc/E (in m/ns)
    TVector3 B1 = field.getFieldAt(position);
    TVector3 k1_pos = vel1 * dt;  // Position change (m)
    // CORRECTED: Use velocity in Lorentz force equation F = q(v × B)
    TVector3 k1_mom = q_factor * vel1.Cross(B1) * dt;  // Momentum change (GeV/c)
    
    // Calculate k2
    TVector3 pos2 = position + k1_pos * 0.5;
    TVector3 mom2 = momentum + k1_mom * 0.5;
    double energy2 = sqrt(mom2.Mag2() + mass * mass);
    TVector3 vel2 = mom2 * (c_light / energy2);
    TVector3 B2 = field.getFieldAt(pos2);
    TVector3 k2_pos = vel2 * dt;
    // CORRECTED: Use velocity in Lorentz force equation
    TVector3 k2_mom = q_factor * vel2.Cross(B2) * dt;
    
    // Calculate k3
    TVector3 pos3 = position + k2_pos * 0.5;
    TVector3 mom3 = momentum + k2_mom * 0.5;
    double energy3 = sqrt(mom3.Mag2() + mass * mass);
    TVector3 vel3 = mom3 * (c_light / energy3);
    TVector3 B3 = field.getFieldAt(pos3);
    TVector3 k3_pos = vel3 * dt;
    // CORRECTED: Use velocity in Lorentz force equation
    TVector3 k3_mom = q_factor * vel3.Cross(B3) * dt;
    
    // Calculate k4
    TVector3 pos4 = position + k3_pos;
    TVector3 mom4 = momentum + k3_mom;
    double energy4 = sqrt(mom4.Mag2() + mass * mass);
    TVector3 vel4 = mom4 * (c_light / energy4);
    TVector3 B4 = field.getFieldAt(pos4);
    TVector3 k4_pos = vel4 * dt;
    // CORRECTED: Use velocity in Lorentz force equation
    TVector3 k4_mom = q_factor * vel4.Cross(B4) * dt;
        
        // Diagnostic output
        if (false) { // Change to false to disable debug output
            std::cout << std::endl;
            std::cout << "RungeKutta4Step::Position: " << position.X() << ", " << position.Y() << ", " << position.Z() << std::endl;
            std::cout << "RungeKutta4Step::Momentum: " << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << std::endl;
            std::cout << "RungeKutta4Step::Energy: " << energy << " GeV, γ = " << energy/mass << std::endl;
            std::cout << "RungeKutta4Step::Velocity: " << vel1.Mag() << " m/ns (" << vel1.Mag()/c_light << " × c)" << std::endl;
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
        
        // Update position and momentum using weighted average of increments
        position += (k1_pos + 2.0 * k2_pos + 2.0 * k3_pos + k4_pos) * (1.0/6.0);
        momentum += (k1_mom + 2.0 * k2_mom + 2.0 * k3_mom + k4_mom) * (1.0/6.0);
    
        
        // Calculate updated energy for diagnostics
        double new_energy = sqrt(momentum.Mag2() + mass * mass);
        
        // Debug output
        if (false) { // Change to false to disable debug output
            std::cout << "\tRungeKutta4Step::New Position: " << position.X() << ", " << position.Y() << ", " << position.Z() << " m" << std::endl;
            std::cout << "\tRungeKutta4Step::New Momentum: " << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << " GeV/c" << std::endl;
            std::cout << "\tRungeKutta4Step::New Energy: " << new_energy << " GeV (Change: " << new_energy-energy << " GeV)" << std::endl;
            std::cout << "\tRungeKutta4Step::Momentum magnitude: " << momentum.Mag() << " GeV/c" << std::endl;
        }
    }

public:
    // ParticlePropagator(MagneticField& field, double step_size_m = 0.01)
    //     : field(field), step_size(step_size_m), charge(1.0), mass(0.1057) {
    //     c_light = 299792458.0; // speed of light in m/s
    // }
    // Add this to the ParticlePropagator constructor:
    ParticlePropagator(MagneticField &field, double step_size_m = 0.01,
                       double x_lim = 1.0, double y_lim = 1.0)
        : field(field), step_size(step_size_m), charge(1.0), mass(0.1057),
          x_limit(x_lim), y_limit(y_lim)
    {
        c_light = 0.299792458;
        // c_light = 299792458.0;
    }
    void setParticleProperties(double q, double m)
    {
        charge = q; // in units of elementary charge
        mass = m;   // in GeV/c^2
    }
    void setLimits(double x_lim, double y_lim)
    {
        x_limit = x_lim;
        y_limit = y_lim;
    }
    void setStepSize(double step)
    {
        step_size = step; // in meters
    }

    // Modify the propagate method to check boundaries:
    std::vector<State> propagate(const TVector3 &initial_position,
                                    const TVector3 &initial_momentum,
                                    double max_length = 8.0,
                                    int max_steps = 10000)
                                
    {
        std::vector<State> trajectory;
        trajectory.push_back({initial_position, initial_momentum});

        TVector3 position = initial_position;
        TVector3 momentum = initial_momentum;

        //debug output
        // std::cout << "propagate::Initial Position: " << position.X() << ", " << position.Y() << ", " << position.Z() << std::endl;
        // std::cout << "propagate::Initial Momentum: " << momentum.X() << ", " << momentum.Y() << ", " << momentum.Z() << std::endl;

        double total_length = 0.0;
        bool hit_boundary = false;

        for (int step = 0; step < max_steps && !hit_boundary; step++)
        {
            RungeKutta4Step(position, momentum, step_size);
            Point interpolatedBoundary = interpolate(position.Z());
            // Check if particle hit magnet boundaries
            // if (fabs(position.X()) >= x_limit || fabs(position.Y()) >= y_limit)
            if(fabs(position.X()) >= fabs(interpolatedBoundary.x) || fabs(position.Y()) >= fabs(interpolatedBoundary.y))
            {
                hit_boundary = true;
                if (false) {
                std::cout << "propagate::Particle hit boundary at position: ("
                          << position.X() << ", " << position.Y() << ", "
                          << position.Z() << ")" << std::endl;
                }
            }

            trajectory.push_back({position, momentum});

            if (trajectory.size() > 1)
            {
                total_length += (trajectory[trajectory.size() - 1].position -
                                 trajectory[trajectory.size() - 2].position)
                                    .Mag();
            }

            if (total_length >= max_length)
            {
                break;
            }
        }

        return trajectory;
    }

    // Overloaded version that takes a TLorentzVector for momentum
    std::vector<State> propagate(const TVector3 &initial_position,
                                    const TLorentzVector &four_momentum,
                                    double max_length = 10.0,
                                    int max_steps = 10000)
    {
        // Extract 3-momentum from 4-momentum
        TVector3 momentum = four_momentum.Vect();
        // Extract mass from 4-momentum
        mass = four_momentum.M();

        return propagate(initial_position, momentum, max_length, max_steps);
    }
};

// Main function to demonstrate particle propagation
int track_propagatorRK(TString inputFieldASCII = "field_map_lhcb.txt", int nParticles = 1000, bool verboseInfo = false) 
{
    // Check command line arguments
    if (inputFieldASCII == "")
    {
        std::cerr << "Error: No magnetic field file provided." << std::endl;
        std::cerr << "Usage: ./track_propagator-RK field_map.txt" << std::endl;
        return 1;
    }

    std::string field_file = inputFieldASCII.Data();

    // Load magnetic field
    MagneticField field(field_file);

    // Create particle propagator with magnet boundaries
    double x_limit = 1.5;  // 1.5 m magnet half-width
    double y_limit = 2.0;  // 2.0 m magnet half-height
    ParticlePropagator propagator(field, 0.01, x_limit, y_limit); // 0.01m step size

    // Create ROOT output file for all particle trajectories
    TFile *rootFile = new TFile("propagation_multiparticle.root", "RECREATE");
    
    // Create a tree to store all trajectories
    TTree *tree = new TTree("trajectories", "Multiple Particle Trajectories");
    
    // Variables for tree branches
    int particle_id;
    double x, y, z, px_out, py_out, pz_out, energy, mass_value, charge_value;
    int step_num;
    bool hit_boundary;
    
    // Create branches
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

    // Create output file for text data
    std::ofstream outfile("trajectories.dat");
    outfile << "# particle_id step x y z px py pz energy mass charge hit_boundary" << std::endl;
    
    // Create TMultiGraph objects for each projection to overlay all trajectories
    TMultiGraph *mg_xy = new TMultiGraph("mg_xy", "All Trajectories XY Projection;X (m);Y (m)");
    TMultiGraph *mg_xz = new TMultiGraph("mg_xz", "All Trajectories XZ Projection;X (m);Z (m)");
    TMultiGraph *mg_yz = new TMultiGraph("mg_yz", "All Trajectories YZ Projection;Y (m);Z (m)");
    
    // Color array for different particles
    int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet, kPink, kTeal, kYellow};
    
    // Generate multiple particles with different properties
    for (int i = 0; i < nParticles; i++)
    {
        // Set particle properties - vary by particle ID
        double charge = 1.0;      // in units of e
        
        // Different particle types based on ID
        double mass;
        std::string particle_type;
        
        // Choose particle type (optional: rotate through different particle types)
        // switch(i % 3) {
        //     case 0:
        //         mass = 0.000511;  // electron
        //         particle_type = "electron";
        //         break;
        //     case 1:
        //         mass = 0.13957039; // pion
        //         particle_type = "pion";
        //         break;
        //     case 2:
        //         mass = 0.1057; // muon
        //         particle_type = "muon";
        //         break;
        //     default:
                mass = 0.13957039; // pion by default
                particle_type = "pion";
        // }
        
        propagator.setParticleProperties(charge, mass);
        mass_value = mass;
        charge_value = charge;
        
        // Vary initial momentum (e.g., log-spaced energies from 0.1 to 5 GeV)
        double log_min = -1.0;  // log10(0.1) = -1
        double log_max = 0.7;   // log10(5) ≈ 0.7
        double momentum_magnitude = pow(10.0, log_min + i * ((log_max - log_min) / (nParticles - 1)));
        
        // Small random transverse component for unique trajectories
        // double angle = (i * 45.0 / nParticles) * (M_PI / 180.0);  // Angle in radians (0 to 45 degrees)
        double angle = 0.0;
        double px = 0.1 * momentum_magnitude * sin(angle);  // Small x component
        double py = 0.0;                                    // No y component
        double pz = momentum_magnitude * cos(angle);        // Main momentum along z
        
        // Small variation in initial position (optional)
        // double x_offset = (i % 3) * 0.05 - 0.05;  // Small x variation (-0.05, 0, 0.05)
        // double y_offset = ((i / 3) % 3) * 0.05 - 0.05;  // Small y variation
        double x_offset = 0.0;
        double y_offset = 0.0;
        
        TVector3 initial_position(x_offset, y_offset, 0.0);
        
        double total_energy = sqrt(px*px + py*py + pz*pz + mass*mass);
        TLorentzVector four_momentum(px, py, pz, total_energy);
        
        // For storing original momentum
        double orig_px = px;
        double orig_py = py;
        double orig_pz = pz;
        
        if(verboseInfo){
        std::cout << "\n=======================================" << std::endl;
        std::cout << "Propagating particle " << i << " (" << particle_type << ")" << std::endl;
        std::cout << "Mass: " << mass << " GeV/c²" << std::endl;
        std::cout << "Initial momentum: (" << px << ", " << py << ", " << pz 
                  << ") GeV/c (magnitude: " << momentum_magnitude << " GeV/c)" << std::endl;
        std::cout << "Initial position: (" << initial_position.X() << ", " 
                  << initial_position.Y() << ", " << initial_position.Z() << ") m" << std::endl;
        }
        // Propagate particle
        std::vector<State> trajectory = propagator.propagate(
            initial_position, four_momentum, 8.0, 50000);
        
        // Create graphs for this particle's trajectory
        TGraph *graph_xy = new TGraph(trajectory.size());
        TGraph *graph_xz = new TGraph(trajectory.size());
        TGraph *graph_yz = new TGraph(trajectory.size());
        
        graph_xy->SetLineColor(colors[i % 10]);
        graph_xz->SetLineColor(colors[i % 10]);
        graph_yz->SetLineColor(colors[i % 10]);
        
        graph_xy->SetLineWidth(2);
        graph_xz->SetLineWidth(2);
        graph_yz->SetLineWidth(2);
        
        // Process trajectory points
        particle_id = i;
        hit_boundary = false;
        
        for (size_t j = 0; j < trajectory.size(); j++)
        {
            // Set coordinates for current point
            x = trajectory[j].position.X();
            y = trajectory[j].position.Y();
            z = trajectory[j].position.Z();
            px_out = trajectory[j].momentum.X();
            py_out = trajectory[j].momentum.Y();
            pz_out = trajectory[j].momentum.Z();
            energy = sqrt(px_out*px_out + py_out*py_out + pz_out*pz_out + mass*mass);
            step_num = j;
            
            // Check if hitting boundary
            if (fabs(x) >= x_limit || fabs(y) >= y_limit) {
                hit_boundary = true;
            }
            
            // Add point to graphs
            graph_xy->SetPoint(j, x, y);
            graph_xz->SetPoint(j, x, z);
            graph_yz->SetPoint(j, y, z);
            
            // Fill tree
            tree->Fill();
            
            // Write to text file
            outfile << particle_id << " " << step_num << " " 
                   << x << " " << y << " " << z << " " 
                   << px_out << " " << py_out << " " << pz_out << " "
                   << energy << " " << mass << " " << charge << " " 
                   << (hit_boundary ? 1 : 0) << std::endl;
        }
        
        // Add this particle's graphs to the multi-graphs
        mg_xy->Add(graph_xy, "l");
        mg_xz->Add(graph_xz, "l");
        mg_yz->Add(graph_yz, "l");
        
        if(verboseInfo){
        std::cout << "Final position: (" << trajectory.back().position.X() << ", " 
                  << trajectory.back().position.Y() << ", " << trajectory.back().position.Z() << ") m" << std::endl;
        std::cout << "Final momentum: (" << trajectory.back().momentum.X() << ", " 
                  << trajectory.back().momentum.Y() << ", " << trajectory.back().momentum.Z() 
                  << ") GeV/c" << std::endl;
        std::cout << "Momentum conservation: " 
                  << 100.0 * (trajectory.back().momentum.Mag() / momentum_magnitude) << "%" << std::endl;
        std::cout << "Points in trajectory: " << trajectory.size() << std::endl;
        }
        // // Create individual canvases for this particle (optional)
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
    
    // Draw all trajectories together
    TCanvas *c_all_xy = new TCanvas("c_all_xy", "All Trajectories XY Projection", 800, 600);
    mg_xy->Draw("A");
    
    // Add boundary box for XY view
    TBox *boundary_xy = new TBox(-x_limit, -y_limit, x_limit, y_limit);
    boundary_xy->SetFillStyle(0);
    boundary_xy->SetLineColor(kBlack);
    boundary_xy->SetLineWidth(2);
    boundary_xy->SetLineStyle(2);
    boundary_xy->Draw("same");
    c_all_xy->SaveAs("c_all_xy.pdf");
    c_all_xy->Write();
    
    TCanvas *c_all_xz = new TCanvas("c_all_xz", "All Trajectories XZ Projection", 600, 1200);
    mg_xz->Draw("A");
    
    // Add boundary lines for XZ view
    TLine *boundary_xz_left = new TLine(-x_limit, 0, -x_limit, 8);
    TLine *boundary_xz_right = new TLine(x_limit, 0, x_limit, 8);
    boundary_xz_left->SetLineColor(kBlack);
    boundary_xz_right->SetLineColor(kBlack);
    boundary_xz_left->SetLineWidth(2);
    boundary_xz_right->SetLineWidth(2);
    boundary_xz_left->SetLineStyle(2);
    boundary_xz_right->SetLineStyle(2);
    boundary_xz_left->Draw();
    boundary_xz_right->Draw();
    c_all_xz->SaveAs("c_all_xz.pdf");
    c_all_xz->Write();
    
    TCanvas *c_all_yz = new TCanvas("c_all_yz", "All Trajectories YZ Projection", 800, 600);
    mg_yz->Draw("A");
    
    // Add boundary lines for YZ view
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
    
    // // Create a TLegend with particle info
    // TCanvas *c_legend = new TCanvas("c_legend", "Particle Legend", 600, 400);
    // TLegend *legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    // legend->SetHeader("Particle Legend", "C");
    
    // for (int i = 0; i < nParticles; i++) {
    //     // Generate same particle properties as in the loop
    //     std::string particle_type;
    //     switch(i % 3) {
    //         case 0: particle_type = "electron"; break;
    //         case 1: particle_type = "pion"; break;
    //         case 2: particle_type = "muon"; break;
    //         default: particle_type = "pion";
    //     }
        
    //     double momentum_magnitude = 1.0 * pow(10.0, i * (2.0 / nParticles));
        
    //     TGraph *g = new TGraph();
    //     g->SetLineColor(colors[i % 10]);
    //     g->SetLineWidth(2);
    //     legend->AddEntry(g, Form("ID %d: %s (p = %.2f GeV/c)", 
    //                           i, particle_type.c_str(), momentum_magnitude), "l");
    // }
    
    // legend->Draw();
    // c_legend->SaveAs("c_legend.pdf");
    // c_legend->Write();
    
    // Write and close files
    tree->Write();
    rootFile->Close();
    outfile.close();

    std::cout << "\n=======================================" << std::endl;
    std::cout << "Completed propagation of " << nParticles << " particles" << std::endl;
    std::cout << "Results saved to trajectories.dat and propagation_multiparticle.root" << std::endl;

    return 0;
}
