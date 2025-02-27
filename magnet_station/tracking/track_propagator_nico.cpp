#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TH2F.h>
#include <TCanvas.h>

class Point {
public:
    double x, y, z, t;
    Point(double x, double y, double z, double t = 0.0) : x(x), y(y), z(z), t(t) {}
};

double getXLimit(double z) {
    return 0.043*z + 237.894 + 50; //TODO added 50 to move gap further out
}

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

    
Point interpolate(double z) {
    // Default value if z is not covered by any surface
    Point default_value(getXLimit(z), 155.0, z);
    bool verbosityhere = 0;
    if(z > 495 && z < 505){
        // verbosityhere = 1;
    }
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
            return Point(x, y, z);
        }
    }
    // std::cout << "Warning: z = " << z/10 << " is not covered by any surface" << std::endl;
    return default_value;
}

double gauss(double x, double mu, double sig) {
    return exp(-pow((x - mu) / sig, 2.0) / 2.0) / (sqrt(2.0 * M_PI) * sig);
}

Point GetFieldValue(const Point& point) {
    Point B(0.0, 0.0, 0.0);
    if (fabs(point.x) > 1800) B.x = 0;
    if (fabs(point.y) > 2000) B.y = 0;
    B.y = (-0.46530 * gauss(point.z, 431.0, 82.9) - 0.864 * gauss(point.z, 574.0, 154.0)) / (2 * M_PI * M_PI);
    B.z = 0;
    return B;
}


Point cross_field(Point vtx_in, Point p_in, double q, double xlim = 250, double z = 810, std::vector<Point>* track = nullptr) {
    double mom = p_in.z;
    if (mom == 0) return vtx_in;
    std::vector<double> x = {vtx_in.x, vtx_in.y, vtx_in.z, p_in.x, p_in.y, q / mom};
    double delta_z = 5.0; // mm
    std::vector<double> xold = x;
    bool hit_magnet_wall = false;
    while (x[2] < z && !hit_magnet_wall) {
        xold = x;
        double v = 0.299792458 * 1e-2 * 1000; // (GeV/c)/(T*cm)
        std::vector<double> xt = x;
        Point B = GetFieldValue(Point(xt[0], xt[1], xt[2]));
        double sr = sqrt(1.0 + xt[3] * xt[3] + xt[4] * xt[4]);
        double Ax = sr * (xt[4] * (xt[3] * B.x + B.z) - (1.0 + xt[3] * xt[3]) * B.y);
        double Ay = sr * (-xt[3] * (xt[4] * B.y + B.z) + (1.0 + xt[4] * xt[4]) * B.x);
        std::vector<double> kx1 = {delta_z * xt[3], delta_z * xt[4], delta_z, delta_z * xt[5] * v * Ax, delta_z * xt[5] * v * Ay, 0.0};
        for (size_t i = 0; i < 6; ++i) x[i] += kx1[i];
        if (track != nullptr) track->emplace_back(x[0], x[1], x[2]);
        // hit_magnet_wall = fabs(x[0]) >= xlim || fabs(x[1]) >= xlim;
        Point p_surface = interpolate(x[2]);
        // std::cout << "p_surface.x = " << p_surface.x << ", p_surface.y = " << p_surface.y << std::endl;
        hit_magnet_wall = fabs(x[0]) >=p_surface.x;// || fabs(x[1]) >= p_surface.y;
    }
    if (hit_magnet_wall) {
        z = x[2];
        std::cout << "\t\thit MS panel at z = " << z << " with x = " << x[0] << " and y = " << x[1] << std::endl;
    }

    // Linear interpolation of x, y to the requested z
    double deltaz = x[2] - xold[2];
    double deltay = x[1] - xold[1];
    double deltax = x[0] - xold[0];

    if (deltaz == 0) return Point(x[0], x[1], x[2]);

    std::vector<double> deltas(6);
    for (size_t id = 0; id < 6; ++id) deltas[id] = x[id] - xold[id];
    x[0] = deltax / deltaz * (z - x[2]) + x[0];
    x[1] = deltay / deltaz * (z - x[2]) + x[1];
    x[2] = z;
    for (size_t id = 0; id < 6; ++id) x[id] = deltas[id] / deltas[2] * (z - x[2]) + x[id];
    return Point(x[0], x[1], x[2]);
}

int track_propagator_nico() {
    const char* inputfile = "/home/niviths/Downloads/magnetStationSims/20250214_backup_all/minimumBias_MS_MagDown_3200plus.root";
    TFile fin(inputfile);
    TTree* ntup = (TTree*)fin.Get("ntup");
    std::cout << ntup->GetEntries() << std::endl;

    TString outputdir = "extraploation_plots/";
    //make output directory
    system("mkdir -p " + outputdir);

    TTreeReader tree(ntup);
    TTreeReaderArray<float> ut_vx(tree, "ut_vx");
    TTreeReaderArray<float> ut_vy(tree, "ut_vy");
    TTreeReaderArray<float> ut_vz(tree, "ut_vz");
    TTreeReaderArray<float> ut_tx(tree, "ut_tx");
    TTreeReaderArray<float> ut_ty(tree, "ut_ty");
    TTreeReaderArray<float> p(tree, "p");
    TTreeReaderArray<int> pid(tree, "pid");
    TTreeReaderArray<int>   key(tree, "key");
    TTreeReaderArray<float> ms_vx(tree, "ms_vx");
    TTreeReaderArray<float> ms_vy(tree, "ms_vy");
    TTreeReaderArray<float> ms_vz(tree, "ms_vz");
    TTreeReaderArray<int>   ms_id(tree, "ms_id");
    TTreeReaderArray<int>   ms_bitID(tree, "ms_segment");

    TH2D* hXZextrapolation = new TH2D("hXZextrapolation", "XZ extrapolation", 100, 350, 800, 100, -250, 250);
    TH2D* hYZextrapolation = new TH2D("hYZextrapolation", "YZ extrapolation", 100, 350, 800, 100, -250, 250);
    TH2D* hdZdY = new TH2D("hdZdY", "dYdZ", 100, -200, 200, 100, -200, 200);
    TH1D* hdX = new TH1D("hdX", "dX", 100, -200, 200);
    TH1D* hdY = new TH1D("hdY", "dY", 100, -200, 200);
    TH1D* hdZ = new TH1D("hdZ", "dZ", 100, -200, 200);

    int numevt = 0;
    while (tree.Next()) {
        // if (numevt > 5) break;
        std::cout << "evt " << numevt << std::endl;

        for (size_t i = 0; i < ut_vz.GetSize(); ++i) {
            if (p[i] > 5000) {
                // std::cout << "skipping track with p > 5000" << std::endl;
                continue;
            }

            double charge = pid[i] / fabs(pid[i]);
            Point part_pos(ut_vx[i] / 10, ut_vy[i] / 10, ut_vz[i] / 10); // in cm
            Point part_mom(ut_tx[i], ut_ty[i], p[i] / 1000); // in GeV

            // extrapolate to magnet
            Point part_proj = cross_field(part_pos, part_mom, charge);

            // TODO: if particle goes beyond magnet THIS NEEDS IMPROVEMENT, it cuts off particles at the magnet entrance
            // if(part_proj.z > 800 || part_proj.z < 525 || abs(part_proj.y) > 150) {
            //     // std::cout << "\t\ttrack goes beyond magnet" << std::endl;
            //     continue;
            // }
            // std::cout << "track number " << i << std::endl;
            // std::cout << "\tpart_pos " << part_pos.x << " " << part_pos.y << " " << part_pos.z << std::endl;
            // std::cout << "\tpart_proj " << part_proj.x << " " << part_proj.y << " " << part_proj.z << std::endl;

            hXZextrapolation->Fill(part_proj.z, part_proj.x);
            hYZextrapolation->Fill(part_proj.z, part_proj.y);

            //loop over ms hits to find the closest hit to the extrapolated point
            double minDist = 1e9;
            double minDist_x = 1e9;
            double minDist_y = 1e9;
            double minDist_z = 1e9;
            int minDistIndex = -1;
            for (size_t j = 0; j < ms_vz.GetSize(); ++j) {
                //check that hit is not in fiber or support
                int isfiber = (ms_bitID[j] >> 28) & 0x3; // 2 bit
                int issupp = (ms_bitID[j] >> 30) & 0x3; // 2 bit
                if (isfiber || issupp) continue;
                Point ms_pos(ms_vx[j] / 10, ms_vy[j] / 10, ms_vz[j] / 10); // in cm
                double dist = sqrt(pow(part_proj.x - ms_pos.x, 2) + pow(part_proj.y - ms_pos.y, 2) + pow(part_proj.z - ms_pos.z, 2));
                hdZdY->Fill(part_proj.z - ms_pos.z, part_proj.y - ms_pos.y);
                hdX->Fill(part_proj.x - ms_pos.x);
                hdY->Fill(part_proj.y - ms_pos.y);
                hdZ->Fill(part_proj.z - ms_pos.z);
                double dist_x = part_proj.x - ms_pos.x;
                double dist_y = part_proj.y - ms_pos.y;
                double dist_z = part_proj.z - ms_pos.z;
                if (dist < minDist) {
                    minDist = dist;
                    minDist_x = dist_x;
                    minDist_y = dist_y;
                    minDist_z = dist_z;
                    minDistIndex = j;
                }
            }

            // std::cout << "\tminDist " << minDist << " minDistIndex " << minDistIndex << std::endl;
            // std::cout << "\t\tminDist_x " << minDist_x << " minDist_y " << minDist_y << " minDist_z " << minDist_z << std::endl;
        }

        numevt++;
    }

    //plot hXZextrapolation
    TCanvas c1("c1", "c1", 800, 600);
    hXZextrapolation->Draw("colz");
    c1.SaveAs(outputdir + "hXZextrapolation.pdf");

    //plot hYZextrapolation
    TCanvas c12("c12", "c12", 800, 600);
    hYZextrapolation->Draw("colz");
    c12.SaveAs(outputdir + "hYZextrapolation.pdf");

    //plot hdZdY
    TCanvas c2("c2", "c2", 800, 600);
    hdZdY->GetXaxis()->SetTitle("dZ [cm]");
    hdZdY->GetYaxis()->SetTitle("dY [cm]");
    hdZdY->Draw("colz");
    c2.SaveAs(outputdir + "hdZdY.pdf");

    //plot hdX
    TCanvas c3("c3", "c3", 800, 600);
    hdX->GetXaxis()->SetTitle("dX [cm]");
    hdX->GetYaxis()->SetTitle("Entries");
    hdX->Draw();
    c3.SaveAs(outputdir + "hdX.pdf");

    //plot hdY
    TCanvas c4("c4", "c4", 800, 600);
    hdY->GetXaxis()->SetTitle("dY [cm]");
    hdY->GetYaxis()->SetTitle("Entries");
    hdY->Draw();
    c4.SaveAs(outputdir + "hdY.pdf");

    //plot hdZ
    TCanvas c5("c5", "c5", 800, 600);
    hdZ->GetXaxis()->SetTitle("dZ [cm]");
    hdZ->GetYaxis()->SetTitle("Entries");
    hdZ->Draw();
    c5.SaveAs(outputdir + "hdZ.pdf");

    TFile fout("output.root", "RECREATE");
    hXZextrapolation->Write();
    fout.Close();

    return 0;
}