// this is a .root macro to plot the incident angle of the particles on the Magnet Station
// it uses the actual hit position and associated momentum vector to calculate the angle
// since the detector plane is rotated for each segment, the angle is calculated with respect to the rotated plane
// contact person: Nicolas Schmidt, schmidt_n@lanl.gov

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <iostream>
#include <vector>
#include <cmath>
// default rotations
// float modPositionsx1[10] = {1675.55, 1792.38, 1909.2, 2019.78, 2136.61, 2253.43, 2363.81, 2480.64, 2603.91, 2720.73};
// float modPositionsy1[10] = {124.196, 124.196, 124.196, 123.567, 123.567, 123.567, 122.893, 122.893, 123.567, 123.567};
// float modPositionsz1[10] = {4644.63, 4920.95, 5197.27, 5469.49, 5745.81, 6022.12, 6294.66, 6570.98, 6851.08, 7127.4};
// float modPositionsx2[10] = {1508.77, 1625.6, 1742.42, 1865.5, 1982.32, 2099.15, 2222.42, 2339.24, 2449.63, 2566.45};
// float modPositionsy2[10] = {86.4716, 86.4716, 86.4716, 87.1004, 87.1004, 87.1004, 87.7749, 87.7749, 87.1004, 87.1004};
// float modPositionsz2[10] = {4891.13, 5167.44, 5443.76, 5724.18, 6000.5, 6276.81, 6556.91, 6833.23, 7105.77, 7382.09};
// float modPositionsx3[10] = {1532.1, 1648.93, 1758.87, 1882.14, 1992.08, 2108.91, 2231.98, 2341.93, 2464.79, 2580.55};
// float modPositionsy3[10] = {1062.23, 1062.23, 1111.66, 1112.34, 1161.78, 1161.78, 1162.41, 1211.84, 1212.43, 1262.4};
// float modPositionsz3[10] = {4692.04, 4968.36, 5247.58, 5527.68, 5806.91, 6083.23, 6363.64, 6642.87, 6923.59, 7207.52};

// 0.5 rotations
// float modPositionsx1[10] = {1606.39, 1723.22, 1840.04, 1956.87, 2073.7, 2190.52, 2307.35, 2424.17, 2541, 2657.82};
// float modPositionsy1[10] = {116.08, 116.08, 116.08, 116.08, 116.08, 116.08, 116.08, 116.08, 116.08, 116.08};
// float modPositionsz1[10] = {4618.94, 4895.26, 5171.58, 5447.9, 5724.21, 6000.53, 6276.85, 6553.17, 6829.49, 7105.81};
// float modPositionsx2[10] = {1577.93, 1694.76, 1811.58, 1928.41, 2045.23, 2162.06, 2278.88, 2395.71, 2512.54, 2629.36};
// float modPositionsy2[10] = {94.5871, 94.5871, 94.5871, 94.5871, 94.5871, 94.5871, 94.5871, 94.5871, 94.5871, 94.5871};
// float modPositionsz2[10] = {4916.81, 5193.13, 5469.45, 5745.77, 6022.09, 6298.4, 6574.72, 6851.04, 7127.36, 7403.68};

// float modPositionsx3[10] = {1475.63, 1592.46, 1702.4, 1819.23, 1929.17, 2046, 2162.82, 2272.77, 2389.59, 2499.54};
// float modPositionsy3[10] = {1055.41, 1055.41, 1104.85, 1104.85, 1154.29, 1154.29, 1154.29, 1203.73, 1203.73, 1253.17};
// float modPositionsz3[10] = {4674.23, 4950.54, 5229.77, 5506.09, 5785.32, 6061.64, 6337.95, 6617.18, 6893.5, 7172.73};

// 1.4 rotations
//  float modPositionsx1[10] = {1708.17, 1824.99, 1941.82, 2063.14, 2179.96, 2296.79, 2417.81, 2534.63, 2655.34, 2772.17};
//  float modPositionsy1[10] = {126.933, 126.933, 126.933, 127.205, 127.205, 127.205, 127.423, 127.423, 127.586, 127.586};
//  float modPositionsz1[10] = {4675.27, 4951.59, 5227.9, 5510.22, 5786.54, 6062.86, 6345.39, 6621.71, 6904.44, 7180.76};
//  float modPositionsx2[10] = {1476.16, 1592.98, 1709.81, 1822.14, 1938.97, 2055.79, 2168.42, 2285.25, 2398.19, 2515.02};
//  float modPositionsy2[10] = {83.7349, 83.7349, 83.7349, 83.4622, 83.4622, 83.4622, 83.2442, 83.2442, 83.0814, 83.0814};
//  float modPositionsz2[10] = {4860.49, 5136.81, 5413.12, 5683.44, 5959.76, 6236.08, 6506.19, 6782.5, 7052.41, 7328.73};
//  float modPositionsx3[10] = {1577.41, 1694.23, 1804.18, 1925.5, 2035.44, 2152.26, 2273.28, 2383.23, 2503.93, 2613.88};
//  float modPositionsy3[10] = {1066.27, 1066.27, 1115.7, 1115.98, 1165.41, 1165.41, 1165.63, 1215.07, 1215.23, 1264.67};
//  float modPositionsz3[10] = {4730.55, 5006.87, 5286.1, 5568.41, 5847.64, 6123.96, 6406.49, 6685.72, 6968.45, 7247.68};

// 1.3 rotations
//  float modPositionsx1[10] = {1698.32, 1815.15, 1931.97, 2053.86, 2170.68, 2287.51, 2409.12, 2525.95, 2647.27, 2764.09};
//  float modPositionsy1[10] = {126.226, 126.226, 126.226, 126.606, 126.606, 126.606, 126.933, 126.933, 127.205, 127.205};
//  float modPositionsz1[10] = {4663.98, 4940.3, 5216.62, 5498.46, 5774.78, 6051.09, 6333.18, 6609.5, 6891.81, 7168.13};

// float modPositionsx2[10] = {1486.01, 1602.83, 1719.66, 1831.42, 1948.24, 2065.07, 2177.11, 2293.94, 2406.27, 2523.09};
// float modPositionsy2[10] = {84.4414, 84.4414, 84.4414, 84.0616, 84.0616, 84.0616, 83.7349, 83.7349, 83.4622, 83.4622};
// float modPositionsz2[10] = {4871.77, 5148.09, 5424.41, 5695.21, 5971.53, 6247.84, 6518.4, 6794.71, 7065.04, 7341.35};

// float modPositionsx3[10] = {1567.56, 1684.39, 1794.33, 1916.22, 2026.16, 2142.99, 2264.6, 2374.54, 2495.86, 2605.8};
// float modPositionsy3[10] = {1065.56, 1065.56, 1065.56, 1115.38, 1115.38, 1115.38, 1165.14, 1165.14, 1214.85, 1214.85};
// float modPositionsz3[10] = {4719.27, 4995.59, 5274.81, 5556.65, 5835.88, 6112.2, 6394.28, 6673.51, 6955.82, 7235.05};

// 9modules 26cm spacing -1700
// float modPositionsx1[9] = {1802.76, 1908.78, 2019.97, 2125.98, 2232, 2342.92, 2448.94, 2559.57, 2665.58};
// float modPositionsy1[9] = {126.226, 126.226, 126.606, 126.606, 126.606, 126.933, 126.933, 127.205, 127.205};
// float modPositionsz1[9] = {4941.12, 5178.52, 5421.35, 5658.75, 5896.15, 6139.22, 6376.63, 6619.94, 6857.34};

// float modPositionsx2[9] = {1594.64, 1700.66, 1801.51, 1907.52, 2013.54, 2114.66, 2220.68, 2322.08, 2428.1};
// float modPositionsy2[9] = {84.4414, 84.4414, 84.0616, 84.0616, 84.0616, 83.7349, 83.7349, 83.4622, 83.4622};
// float modPositionsz2[9] = {5153.12, 5390.52, 5622.5, 5859.91, 6097.31, 6329.05, 6566.45, 6797.95, 7035.35};

// float modPositionsx3[9] = {1673.13, 1772.33, 1883.52, 1982.71, 2088.73, 2199.65, 2298.84, 2409.47, 2508.67};
// float modPositionsy3[9] = {1065.56, 1115, 1115.38, 1164.82, 1164.82, 1165.14, 1214.58, 1214.85, 1264.29};
// float modPositionsz3[9] = {4999.01, 5239.46, 5482.28, 5722.73, 5960.13, 6203.21, 6443.66, 6686.96, 6927.41};

//9modules 30cm spacing - 1800
float modPositionsx1[9] = {1802.76, 1925.09, 2052.59, 2174.92, 2297.24, 2424.47, 2546.8, 2673.74, 2796.07};
float modPositionsy1[9] = {126.226, 126.226, 126.606, 126.606, 126.606, 126.933, 126.933, 127.205, 127.205};
float modPositionsz1[9] = {4941.12, 5215.05, 5494.39, 5768.32, 6042.25, 6321.84, 6595.77, 6875.6, 7149.53};

float modPositionsx2[9] = {1594.64, 1716.97, 1834.13, 1956.46, 2078.78, 2196.21, 2318.54, 2436.26, 2558.58};
float modPositionsy2[9] = {84.4414, 84.4414, 84.0616, 84.0616, 84.0616, 83.7349, 83.7349, 83.4622, 83.4622};
float modPositionsz2[9] = {5153.12, 5427.04, 5695.55, 5969.48, 6243.41, 6511.66, 6785.59, 7053.61, 7327.54};

float modPositionsx3[9] = {1673.13, 1788.64, 1916.14, 2031.64, 2153.97, 2281.2, 2396.7, 2523.64, 2639.15};
float modPositionsy3[9] = {1065.56, 1115, 1115.38, 1164.82, 1164.82, 1165.14, 1214.58, 1214.85, 1264.29};
float modPositionsz3[9] = {4999.01, 5275.98, 5555.33, 5832.3, 6106.23, 6385.82, 6662.8, 6942.63, 7219.6};

TVector3 getVerticalVector(int module, int station)
{
    float x_orientation = 1;
    float y_orientation = 1;
    if (station == 0 || station == 2)
    {
        x_orientation = -1;
    }
    if (station == 0 || station == 1)
    {
        y_orientation = -1;
    }
    // create vectors pointing to two corners of the module
    TVector3 v1(x_orientation * modPositionsx1[module], y_orientation * modPositionsy1[module], modPositionsz1[module]); // bottom left
    TVector3 v2(x_orientation * modPositionsx2[module], y_orientation * modPositionsy2[module], modPositionsz2[module]); // bottom right
    TVector3 v3(x_orientation * modPositionsx3[module], y_orientation * modPositionsy3[module], modPositionsz3[module]); // top left

    // calculate the vector pointing from v1 to v2
    TVector3 v_hor = v2 - v1;
    TVector3 v_ver = v3 - v1;

    // calculate the normal vector of the detector plane
    TVector3 v_vertical = v_ver - v_hor;
    v_vertical.SetZ(0);

    // return the normal vector
    return v_vertical;
}

TVector3 getHorizontalVector(int module, int station)
{
    float x_orientation = 1;
    float y_orientation = 1;
    if (station == 0 || station == 2)
    {
        x_orientation = -1;
    }
    if (station == 0 || station == 1)
    {
        y_orientation = -1;
    }
    // create vectors pointing to two corners of the module
    TVector3 v1(x_orientation * modPositionsx1[module], y_orientation * modPositionsy1[module], modPositionsz1[module]); // bottom left
    TVector3 v2(x_orientation * modPositionsx2[module], y_orientation * modPositionsy2[module], modPositionsz2[module]); // bottom right
    TVector3 v3(x_orientation * modPositionsx3[module], y_orientation * modPositionsy3[module], modPositionsz3[module]); // top left

    // calculate the vector pointing from v1 to v2
    TVector3 v_hor = v2 - v1;
    TVector3 v_ver = v3 - v1;

    // calculate the normal vector of the detector plane
    TVector3 v_horizontal = v_ver - v_hor;
    v_horizontal.SetY(0);

    // v_horizontal.RotateY(-x_orientation*0.4);

    // return the normal vector
    return v_horizontal;
}

TVector3 getModulePlaneVector(int module, int station)
{
    float x_orientation = 1;
    float y_orientation = 1;
    if (station == 0 || station == 2)
    {
        x_orientation = -1;
    }
    if (station == 0 || station == 1)
    {
        y_orientation = -1;
    }
    // create vectors pointing to two corners of the module
    TVector3 v1(x_orientation * modPositionsx1[module], y_orientation * modPositionsy1[module], modPositionsz1[module]); // bottom left
    TVector3 v2(x_orientation * modPositionsx2[module], y_orientation * modPositionsy2[module], modPositionsz2[module]); // bottom right
    TVector3 v3(x_orientation * modPositionsx3[module], y_orientation * modPositionsy3[module], modPositionsz3[module]); // top left

    // calculate the vector pointing from v1 to v2
    TVector3 v_hor = v2 - v1;
    TVector3 v_ver = v3 - v1;

    // calculate the normal vector of the detector plane
    TVector3 v_plane = v_ver - v_hor;

    // return the normal vector
    return v_plane;
}

// make function that returns the y position for a given segment according to a lookup table of segment positions
// TODO adjust for new 95cm long modules
double getBarY(int segment, int layer)
{
    float modPositionsy1 = 0; // NOTE this offset is defined in parameters.xml
    // float modPositionsy1 = 100; // NOTE this offset is defined in parameters.xml
    float modPositionsy2[16] = {1050, 1050, 1050, 1050, 1050, 1050, 1150, 1150, 1150, 1150, 1150, 1250, 1250, 1250, 1250, 1250};

    // lookup table of segment positions
    double bar_positions[4][10] =
        {
            {1.5, 6, 13.5, 24, 37.5, 54, 77.5, 96},
            {1.5, 6, 12, 19.5, 30.5, 45, 65, 88.5},
            {1.5, 6, 13.5, 22.5, 34, 49.5, 71, 92},
            {3, 9, 16.5, 27, 41, 59.5, 85}};

    // std::cout << "layer: " << layer << " segment: " << segment << std::endl;
    // std::cout << "bar_positions[layer][segment]: " << bar_positions[layer][segment] << std::endl;

    // get y position for the given segment
    double y_tmp = 10 * bar_positions[layer][segment] + modPositionsy1;
    TVector3 y_vec(0, y_tmp, 0);
    // y_vec.RotateY(0.15);
    // return y position
    return y_vec.Y();
}

TVector3 getBarXYZ(int bar, int layer, int module, int segment, int station)
{

    float orientation = 1;
    if (station == 0 || station == 2)
    {
        orientation = -1;
    }
    float tiltorientation = 1;
    if (station == 0 || station == 3)
    {
        tiltorientation = -1;
    }

    // each module has 34 bars
    int nBars = 57;
    // each module has 4 layers
    int nLayers = 4;
    // each layer is separated by 10 mm
    double layerSpacing = 10;

    // create vectors pointing to two corners of the module
    TVector3 v1(modPositionsx1[module], modPositionsy1[module], modPositionsz1[module]); // bottom left
    TVector3 v2(modPositionsx2[module], modPositionsy2[module], modPositionsz2[module]); // bottom right
    TVector3 v3(modPositionsx3[module], modPositionsy3[module], modPositionsz3[module]); // top left

    // calculate the vector pointing from v1 to v2
    TVector3 v = v2 - v1;
    TVector3 v_tilt = v3 - v1;

    // calculate the length of the vector
    double length = v.Mag();
    double length_tilt = v_tilt.Mag();

    // calculate the unit vector
    TVector3 u = v.Unit();
    TVector3 u_tilt = v_tilt.Unit();

    // determine vector orthogonal to v
    TVector3 v2_norm = u.Cross(u_tilt);

    double ySeg = getBarY(segment, layer);

    // calculate the position of the bar using u and v2_norm
    //  double x = modPositionsx1[module] + u.X() * ( (bar + 0.5) * length / nBars);
    //  double z = modPositionsz1[module] + u.Z() * ( (bar + 0.5) * length / nBars);
    TVector3 pos = v1 + u * ((bar + 0.5) * length / nBars) + orientation * v2_norm * ((layer - 1.5) * layerSpacing) + ySeg * u_tilt;
    // double x = modPositionsx1[module] + zOffset + u.X() * ((bar + 0.5) * length / nBars) + orientation * v2_norm.X() * ((layer - 1.5) * layerSpacing);
    // double z = modPositionsz1[module] + zOffset + u.Z() * ((bar + 0.5) * length / nBars) + orientation * v2_norm.Z() * ((layer - 1.5) * layerSpacing);
    double x = pos.X();
    double y = pos.Y();
    double z = pos.Z();

    TVector3 barXYZ(x, y, z);

    return barXYZ;
}

// void makeangleplot_new3_moreModules(TString inputfile = "/home/niviths/Downloads/minimumBias_MS_MagDown_1006.root") //tilt
void makeangleplot_new4_10modstilted(TString inputfile = "/home/niviths/Downloads/magnetStationSims/20250227_pp_LayerNumbersFixed/20250227_pp_LayerNumbersFixed.root") // tilt
// void makeangleplot_new2(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_600plus.root") //vertical
{
    TFile fin(inputfile, "READ");
    TTreeReader tree("ntup", &fin);
    double modSpacing = 300 * TMath::Cos(0.42);
    double minZ = 3500 + 2 * 500 + 300;
    double maxZ = 7700;

    // get ms_px, ms_py, ms_pz and ms_vx, ms_vy, ms_vz, as well as p and ms_time from the tree
    TTreeReaderArray<float> ms_px(tree, "ms_px");
    TTreeReaderArray<float> ms_py(tree, "ms_py");
    TTreeReaderArray<float> ms_pz(tree, "ms_pz");
    TTreeReaderArray<float> ms_vx(tree, "ms_vx");
    TTreeReaderArray<float> ms_vy(tree, "ms_vy");
    TTreeReaderArray<float> ms_vz(tree, "ms_vz");
    TTreeReaderArray<float> ms_energy(tree, "ms_energy");
    TTreeReaderArray<int> ms_npart(tree, "ms_npart");
    TTreeReaderArray<int> ms_id(tree, "ms_id");
    TTreeReaderArray<float> ft_vx(tree, "ft_vx");
    TTreeReaderArray<float> ft_vy(tree, "ft_vy");
    TTreeReaderArray<float> ft_vz(tree, "ft_vz");
    TTreeReaderArray<int> ft_id(tree, "ft_id");
    TTreeReaderArray<int> key(tree, "key");
    TTreeReaderArray<int> parent_key(tree, "parent_key");
    TTreeReaderArray<int> nFThits(tree, "nFThits");
    TTreeReaderArray<float> p(tree, "p");
    TTreeReaderArray<int> pid(tree, "pid");
    TTreeReaderArray<float> pt(tree, "pt");
    TTreeReaderArray<float> eta(tree, "eta");
    TTreeReaderArray<float> ms_time(tree, "ms_time");
    TTreeReaderArray<int> ms_bitID(tree, "ms_segment");
    const int nModules = 9;
    // create 2D histogram of incident angle of particles on the Magnet Station
    TH2F *h2 = new TH2F("h2", "Incident angle of particles on the Magnet Station", 125, minZ, maxZ, 90, 0, 180);
    TH2F *h2_alt = new TH2F("h2_alt", "Incident angle of particles on the Magnet Station", 125, minZ, maxZ, 90, 0, 180);
    TH2F *h2_hor = new TH2F("h2_hor", "Incident angle of particles on the Magnet Station - horizontal plane", 125, minZ, maxZ, 90, 0, 180);
    TH2F *h2_hor_alt = new TH2F("h2_hor_alt", "Incident angle of particles on the Magnet Station - horizontal plane", 125, minZ, maxZ, 90, 0, 180);
    TH2F *h2_ver = new TH2F("h2_ver", "Incident angle of particles on the Magnet Station - vertical plane", 125, minZ, maxZ, 90, 0, 180);
    TH2F *h2_ver_alt = new TH2F("h2_ver_alt", "Incident angle of particles on the Magnet Station - vertical plane", 125, minZ, maxZ, 90, 0, 180);
    TH2F *h2_ver_y = new TH2F("h2_ver_y", "Incident angle of particles on the Magnet Station - vertical plane", 150, -1500, 1500, 90, 0, 180);
    TH2F *hAcc = new TH2F("hAcc", "Acceptance of particles full system", 125, 0, 20, 60, 1.8, 4.5);
    TH2F *hAcc_ft = new TH2F("hAcc_ft", "Acceptance of particles on the Forward Tracker", 125, 0, 20, 60, 1.8, 4.5);
    TH2F *hAcc_ms = new TH2F("hAcc_ms", "Acceptance of particles on the Magnet Station", 125, 0, 20, 60, 1.8, 4.5);
    TH2F *hAcc_ms_noFirstModule = new TH2F("hAcc_ms_noFirstModule", "Acceptance of particles on the Magnet Station skipping the first module", 125, 0, 20, 60, 1.8, 4.5);
    TH2F *hSeg[nModules] = {nullptr};
    TH1F *hSeg_p[nModules] = {nullptr};
    TH1F *hSeg_pt[nModules] = {nullptr};
    TH1F *hSeg_time[nModules] = {nullptr};
    for (int i = 0; i < nModules; i++)
    {
        hSeg[i] = new TH2F(Form("hSeg_%d", i), Form("Incident angle of particles on the Magnet Station module %d", i), 125, minZ, maxZ, 60, 0, 120);
        hSeg_p[i] = new TH1F(Form("hSeg_p_%d", i), Form("Momentum of particles on the Magnet Station module %d", i), 100, 0, 5000);
        hSeg_pt[i] = new TH1F(Form("hSeg_pt_%d", i), Form("Transverse momentum of particles on the Magnet Station module %d", i), 80, 0, 1600);
        hSeg_time[i] = new TH1F(Form("hSeg_time_%d", i), Form("Time of particles on the Magnet Station module %d", i), 100, 0, 300);
    }
    TH2F *hnHitsModulePerStation = new TH2F("hnHitsModulePerStation", "Number of hits per module per MC particle", 17, -0.5, 16.5, 20, -.5, 19.5);
    // TH2F *hnHitsPerParticle = new TH2F("hnHitsPerParticle", "Number of hits per MC particle", 17, -0.5, 16.5, 20, -.5, 19.5);
    TH1F *hnHitsPerParticle = new TH1F("hnHitsPerParticle", "Number of hits per MC particle", 20, -.5, 19.5);
    TH1F *hnModulesHitPerParticle = new TH1F("hnModulesHitPerParticle", "Number of modules hit per MC particle", 10, -.5, 9.5);
    TH2F *hnModulesHitPerParticleVsEta = new TH2F("hnModulesHitPerParticleVsEta", "Number of modules hit per MC particle vs Eta", 100, 1.5, 4.5, 10, -.5, 9.5);
    TH2F *hnModulesHitPerParticleVsPt = new TH2F("hnModulesHitPerParticleVsPt", "Number of modules hit per MC particle vs Pt", 100, 0, 2000, 10, -.5, 9.5);
    TH1F *hNPartPerBar = new TH1F("hNPartPerBar", "Number of particles in same bar", 10, -.5, 9.5);
    // TH2F *hnModulesHitPerParticle = new TH2F("hnModulesHitPerParticle", "Number of modules hit per MC particle vs E", 100, 0, 5000, 5, -.5, 4.5);
    TH2F *hnLayerHitsPerModule = new TH2F("hnLayerHitsPerModule", "Number of hits per layer per module", 17, -0.5, 16.5, 4, -0.5, 3.5);

    TH2F *pTvsP_lowHorAngle = new TH2F("pTvsP_lowHorAngle", "pT vs P for particles with low horizontal angle", 100, 0, 5000, 50, 0, 2500);
    TH2F *pTvsP_highHorAngle = new TH2F("pTvsP_highHorAngle", "pT vs P for particles with high horizontal angle", 100, 0, 5000, 50, 0, 2500);

    TH2F *deltaYvsDeltaZ = new TH2F("deltaYvsDeltaZ", "deltaY vs deltaZ", 100, -200, 200, 100, -20, 20);
    TH2F *deltaYvsSegment = new TH2F("deltaYvsSegment", "deltaY vs segment", 10, -0.5, 9.5, 100, -200, 200);
    float maxDelta = 80;
    TH2F *deltaXvsDeltaZ = new TH2F("deltaXvsDeltaZ", "deltaX vs deltaZ", 100, -maxDelta, maxDelta, 100, -maxDelta, maxDelta);
    TH2F *deltaXvsDeltaZPerLayer[4] = {nullptr};
    for (int i = 0; i < 4; i++)
    {
        deltaXvsDeltaZPerLayer[i] = new TH2F(Form("deltaXvsDeltaZPerLayer_%d", i), Form("deltaX vs deltaZ for layer %d", i), 100, -maxDelta, maxDelta, 100, -maxDelta, maxDelta);
    }
    TH2F *deltaXvsDeltaZPerModule[16] = {nullptr};
    for (int i = 0; i < 16; i++)
    {
        deltaXvsDeltaZPerModule[i] = new TH2F(Form("deltaXvsDeltaZPerModule_%d", i), Form("deltaX vs deltaZ for module %d", i), 100, -maxDelta, maxDelta, 100, -maxDelta, maxDelta);
    }
    // make histogram for time difference between layer 3 and layer 0 versus 1/pT
    TH2F *hTimeDiffvs1pT = new TH2F("hTimeDiffvs1pT", "Time difference between layer 3 and layer 0 vs 1/pT", 100, 0, 5.0, 100, 0, 3);
    // make histograms for x y of hits of particles which hit the ft and ms
    TH2F *hXYft = new TH2F("hXYft", "XY of hits in FT", 100, -3000, 3000, 100, -3000, 3000);
    TH2F *hXYms = new TH2F("hXYms", "XY of hits in MS", 100, -3000, 3000, 100, -3000, 3000);
    // make histograms for x y of hits of particles in the ft where the mother particle hit the ms
    TH2F *hXYftMother = new TH2F("hXYftMother", "XY of hits in FT where mother hit MS", 100, -1000, 1000, 100, -1000, 1000);

    // manual rotation of the detector plane for each module
    //  float add_rotations[nModules] = {0.0, 0.1, 0.3, 0.4, 0.5, 0.5, 0.6};
    //  float add_rotations[nModules] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1}; //used for 14
    // float add_rotations[nModules] = {0.6, 0.6, 0.6, 0.6, 0.7, 0.8, 0.9}; // used for 15,16
    // float add_rotations[nModules] = {0.7, 0.7, 0.6, 0.6, 0.7, 0.8, 0.9}; // used for 20plus
    // float add_rotations[nModules] = {0.7, 0.7, 0.7, 0.8, 0.9, 1.0, 1.0}; // used for 15,16
    // float add_rotations[nModules] = {0.8, 0.8, 0.7, 0.7, 0.7, 0.8, 0.9, 1.0, 1.0}; // used for 15,16
    // float add_rotations[nModules] = {1.1, 1.0, 0.7, 0.6, 0.7, 0.8, 0.9}; //BEST SO FAR, but unfeasible
    // float add_rotations[nModules] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    float add_rotations[nModules] = {1.2, 1.2, 1.25, 1.25, 1.25, 1.3, 1.3, 1.35, 1.35}; // 1500 10 modules with tilt
    // float add_rotations[nModules] = {0.9, 0.9, 0.9, 0.95, 0.95, 0.95, 1.0, 1.0, 1.05, 1.1}; //1500 10 modules with tilt
    // float add_rotations[nModules] = {1.0, 1.0, 1.0, 0.95, 0.95, 0.95, 0.9, 0.9, 0.95, 0.95, 1.0, 1.0, 1.05, 1.05, 1.1, 1.1}; /1200,1300
    // float add_rotations[nModules] = {0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    // float add_rotations[nModules] = {0};

    // create 7 vectors to store the average values of the incident angle for each module
    TVector3 vIncidentAverage[nModules];
    int curr_entry = 0;
    bool verbositysetting = 0;
    // loop over tree
    int maxbar = 0;
    bool onceMsg = 1;
    while (tree.Next())
    {
        // cout << "Entry: " << curr_entry << endl;
        curr_entry++;
        // if (curr_entry > 1500 || onceMsg)
        // {
        //     if (onceMsg)
        //     {
        //         onceMsg = 0;
        //         cout << "processing only 1500 events" << endl;
        //     } else {
        //         break;
        //     }
        // }
        if (curr_entry % 500 == 0)
            cout << "Entry: " << curr_entry << endl;
        if (verbositysetting)
            cout << "Entry: " << curr_entry << endl;
        if (verbositysetting)
            cout << "Number of hits: " << ms_px.GetSize() << endl;
        for (auto i = 0U; i < pt.GetSize(); ++i)
        {
            // if (parent_key[i] > -2000)
            //     continue;
            int nHitsModulePerStation[4][nModules] = {{0}};
            int nLayerHitsPerModule[4][nModules][4] = {{{0}}};
            if (p[i] > 5000)
            {
                if (verbositysetting)
                    cout << "Momentum too high" << endl;
                // continue;
            }
            else
            {
                if (verbositysetting)
                    cout << "Momentum ok" << endl;
            }
            // cout << "\tKey: " << key[i] << " particle PID: " << pid[i] << " with momentum: " << p[i] << " and eta: " << eta[i] << endl;
            bool hitInLayer[4] = {0};
            float timeInLayer[4] = {0};
            int nHitsPartSum = 0;
            int nHitsMS = 0;
            int nHitsMS_nofirstMod = 0;
            int nfthits = 0;
            // loop over size of ms_px to get all the hits per event
            for (auto j = 0U; j < ms_px.GetSize(); ++j)
            {
                if (ms_id[j] != key[i])
                    continue;
                // require hits within MS acceptance
                if (ms_vz[j] < 3000 || ms_vz[j] > maxZ)
                {
                    if (verbositysetting)
                        cout << "Hit outside of acceptance" << endl;
                    continue;
                }
                else
                {
                    if (verbositysetting)
                        cout << "Hit in acceptance" << endl;
                }
                // if((ms_vy[j]) <0)
                // // if(abs(ms_vy[j]) >100)
                // {
                //     continue;
                // }
                if (nFThits[j] > 0)
                {
                    if (verbositysetting)
                        cout << "Hit in FT" << endl;
                    continue;
                }
                else
                {
                    if (verbositysetting)
                        cout << "Hit not in FT" << endl;
                }
                // if(eta[j] < 1.5 || eta[j] > 2)
                // if(eta[j] < 2 || eta[j] > 2.5)
                // if(eta[j] < 2.5 || eta[j] > 3)
                // if(eta[j] < 3 || eta[j] > 3.5)
                // if(eta[j] < 3.5 || eta[j] > 4)
                // if(eta[j] < 4 )
                // {
                //     continue;
                // }
                // require certain momentum
                // if (p[j] > 5000)
                // {
                //     if (verbositysetting)
                //         cout << "Momentum too high" << endl;
                //     continue;
                // }
                // else
                // {
                //     if (verbositysetting)
                //         cout << "Momentum ok" << endl;
                // }

                // print layer and segment from ms_bitID
                // cout << "Segment: " << std::bitset<32>(ms_bitID[j]>>8) << endl;
                // 8 bit = 255
                // 7 bit = 127
                // 6 bit = 63
                // 5 bit = 31
                // 4 bit = 15
                // 3 bit = 7
                // <id>system:8,station:3,module:4,layer:3,segment:4,bar:6,isfiber:2,issupp:2</id>

                int station = (ms_bitID[j] >> 8) & 0x7;  // 3 bit
                int module = (ms_bitID[j] >> 11) & 0xF;  // 4 bit
                int layer = (ms_bitID[j] >> 15) & 0x7;   // 3 bit
                int segment = (ms_bitID[j] >> 18) & 0xF; // 4 bit
                int bar = (ms_bitID[j] >> 22) & 0x3F;     // 6 bits (was incorrectly using 0xFF for 8 bits)
                int isfiber = (ms_bitID[j] >> 28) & 0x3;  // 2 bits (was incorrectly shifted by 30)
                int issupp = (ms_bitID[j] >> 30) & 0x3;   // 2 bits (this was missing)
                // if(station!=0) continue;
                // cout << __LINE__ << endl;
                if (bar > maxbar)
                    maxbar = bar;
                // cout << "\t\tStation: " << station << " Module: " << module << " Layer: " << layer << " Segment: " << segment << " Bar: " << bar << " Energy: " << ms_energy[j] << endl;
                // if(layer==1 || layer==2)
                // cout << "Station: " << station << " Module: " << module << " Layer: " << layer << " Segment: " << segment << " Bar: " << bar << endl;
                // if(segment!=1) continue;
                nHitsMS++;
                if (module != 0)
                {
                    nHitsMS_nofirstMod++;
                }
                nHitsPartSum += 1;
                // hNPartPerBar->Fill(ms_npart[j]);
                hitInLayer[layer] = 1;
                if (abs(pid[i]) == 211)
                    timeInLayer[layer] = ms_time[j];
                // create a vector with the direction of the particle
                TVector3 v(ms_px[j], ms_py[j], ms_pz[j]);
                TVector3 v_hor(ms_px[j], 0, ms_pz[j]);
                TVector3 v_ver(ms_px[j], ms_py[j], 0);

                // get vector pointing to the bar from functions
                TVector3 v_bar = getBarXYZ(bar, layer, module, segment, station);

                TVector3 v_true(ms_vx[j], ms_vy[j], ms_vz[j]);

                // print v values and v_bar values
                //  cout << endl;
                //  cout << "v_true: " << v_true.X() << " " << v_true.Y() << " " << v_true.Z() << endl;
                //  cout << "v_bar: " << v_bar.X() << " " << v_bar.Y() << " " << v_bar.Z() << endl;
                // fill deltaYvsDeltaZ with the difference between v and v_bar
                deltaYvsDeltaZ->Fill(abs(v_true.Y()) - abs(v_bar.Y()), abs(v_true.Z()) - abs(v_bar.Z()));
                deltaXvsDeltaZ->Fill(abs(v_true.X()) - abs(v_bar.X()), abs(v_true.Z()) - abs(v_bar.Z()));

                deltaXvsDeltaZPerLayer[layer]->Fill(abs(v_true.X()) - abs(v_bar.X()), abs(v_true.Z()) - abs(v_bar.Z()));
                deltaXvsDeltaZPerModule[module]->Fill(abs(v_true.X()) - abs(v_bar.X()), abs(v_true.Z()) - abs(v_bar.Z()));

                deltaYvsSegment->Fill(segment, abs(v_true.Y()) - abs(v_bar.Y()));

                nHitsModulePerStation[station][module]++;
                nLayerHitsPerModule[station][module][layer]++;

                // add the vector to the average vector for the segment
                if (ms_vx[j] > 0)
                    vIncidentAverage[module] += v;

                // calculate the vector spanning the detector plane
                TVector3 v2(0, 1, 1);
                double z_rotation = 0.15; // 0.34 is about 20 degrees
                // double z_rotation = 0.34; // 0.34 is about 20 degrees
                double y_rotation = 0.4;
                if (station == 1 || station == 3)
                {
                    y_rotation = -y_rotation;
                }
                else
                {
                    y_rotation = y_rotation;
                }

                if (station == 0)
                {
                    z_rotation = z_rotation;
                }
                else if (station == 1)
                {
                    z_rotation = TMath::Pi() - z_rotation;
                }
                else if (station == 2)
                {
                    z_rotation = -z_rotation;
                }
                else if (station == 3)
                {
                    z_rotation = TMath::Pi() + z_rotation;
                }

                // add additional rotation depending on module of particle hit
                // need two cases for the different sides of the detector
                if (station == 1 || station == 3)
                    v2.RotateY(TMath::Pi() - add_rotations[module]);
                else
                    v2.RotateY(add_rotations[module]);

                // if((station==2 && module==0) || (station==3 && module==0)){
                //     cout <<  "Station: " << station << " Module: " << module << " Rotation: " << add_rotations[module] << endl;
                //     cout << v2.X() << " " << v2.Y() << " " << v2.Z() << endl;
                // }
                v2.RotateZ(z_rotation); // tilt around z axis (makes upper and lower modules)
                v2.RotateY(y_rotation); // initial tilt around y-axis (makes modules be in line with the magnet inside)

                TVector3 v2_hor(v2.X(), 0, v2.Z());
                TVector3 v2_ver(v2.X(), v2.Y(), 0);
                // calculate the angle in degrees between the two vectors
                float angle2 = v2.Angle(v) * 180 / 3.14159;
                float angle_alternative = getModulePlaneVector(module, station).Angle(v) * 180 / 3.14159;
                // if(station==0 && module==8)
                // {
                //     //print module vector and v vector
                //     TVector3 v_mod = getModulePlaneVector(module, station);
                //     //normalize v_mod
                //     v_mod = v_mod.Unit();
                //     v = v.Unit();
                //     cout << endl;
                //     cout << "M vector: " << v_mod.X() << " " << v_mod.Y() << " " << v_mod.Z() << endl;
                //     cout << "v vector: " << v.X() << " " << v.Y() << " " << v.Z() << endl;
                //     cout << "angle_alternative: " << angle_alternative << endl;
                // }
                float angle2_hor = v2_hor.Angle(v_hor) * 180 / 3.14159;
                float angle_hor_alt = getHorizontalVector(module, station).Angle(v_hor) * 180 / 3.14159;

                if (angle2_hor < 65)
                {
                    pTvsP_lowHorAngle->Fill(p[i], pt[i]);
                }
                else
                {
                    pTvsP_highHorAngle->Fill(p[i], pt[i]);
                }

                float angle_ver = v2_ver.Angle(v_ver) * 180 / 3.14159;
                float angle_ver_alt = getVerticalVector(module, station).Angle(v_ver) * 180 / 3.14159;
                // if(station==0 || station==3)
                // if(station==1 || station==2)
                // if(station==2 || station==3)
                // if (station == 0 || station == 1)
                //     angle_ver = 180 - angle_ver;

                // fill the histogram
                h2->Fill(ms_vz[j], angle2);
                h2_alt->Fill(ms_vz[j], angle_alternative);
                h2_hor->Fill(ms_vz[j], angle2_hor);
                h2_hor_alt->Fill(ms_vz[j], angle_hor_alt);
                // if(station==2||station==3)
                // if(station==0||station==1)
                h2_ver->Fill(ms_vz[j], angle_ver);
                h2_ver_alt->Fill(ms_vz[j], angle_ver_alt);
                h2_ver_y->Fill(ms_vy[j], angle_ver);
                hSeg[module]->Fill(ms_vz[j], angle2);
                hSeg_p[module]->Fill(p[j]);
                hSeg_pt[module]->Fill(pt[j]);
                hSeg_time[module]->Fill(ms_time[j]);

                // loop over ft hits and find matching hits to the current particle and fill xy histograms
                for (auto k = 0U; k < ft_id.GetSize(); ++k)
                {
                    if (ft_id[k] == key[i])
                    {
                        hXYft->Fill(ft_vx[k], ft_vy[k]);
                    }
                    if (ft_id[k] == parent_key[i])
                    {
                        hXYftMother->Fill(ft_vx[k], ft_vy[k]);
                    }
                }
            } // NOTE END LOOP OVER MS HITS
            // loop over ft hits and find matching hits to the current particle
            for (auto k = 0U; k < ft_id.GetSize(); ++k)
            {
                if (ft_id[k] == key[i])
                {
                    nfthits++;
                }
            }
            if (nHitsMS > 1 || nfthits > 1 || nHitsMS_nofirstMod > 1)
            {
                if (nHitsMS > 1)
                    hAcc_ms->Fill(1000 / (pt[i]), eta[i]);
                if (nfthits > 1)
                    hAcc_ft->Fill(1000 / (pt[i]), eta[i]);
                if (nHitsMS_nofirstMod > 1)
                    hAcc_ms_noFirstModule->Fill(1000 / (pt[i]), eta[i]);

                hAcc->Fill(1000 / (pt[i]), eta[i]);
            }
            if (timeInLayer[0] && timeInLayer[3])
            {
                hTimeDiffvs1pT->Fill(1000 / pt[i], (timeInLayer[3] - timeInLayer[0]));
            }

            if (nHitsPartSum)
                hnHitsPerParticle->Fill(nHitsPartSum);
            // fill hnHitsModulePerStation
            for (int b = 0; b < 4; b++)
            {
                int nModHits = 0;
                for (int c = 0; c < nModules; c++)
                {
                    if (nHitsModulePerStation[b][c] > 0)
                    {
                        hnHitsModulePerStation->Fill(c, nHitsModulePerStation[b][c]);
                        nModHits++;
                    }
                }
                if (nModHits)
                    hnModulesHitPerParticle->Fill(nModHits);
                if (nModHits)
                    hnModulesHitPerParticleVsEta->Fill(eta[i], nModHits);
                if (nModHits)
                    hnModulesHitPerParticleVsPt->Fill(pt[i], nModHits);
            }
            // fill hnLayerHitsPerModule
            for (int k = 0; k < 4; k++)
            {
                for (int l = 0; l < nModules; l++)
                {
                    for (int m = 0; m < 4; m++)
                    {
                        if (nLayerHitsPerModule[k][l][m] > 0)
                            hnLayerHitsPerModule->Fill(l, m);
                    }
                }
            }
            // if (!hitInLayer[0] || !hitInLayer[1] || !hitInLayer[2] || !hitInLayer[3])
            // {
            //     int nhits = 0;
            //     for (int i = 0; i < 4; i++)
            //     {
            //         if (hitInLayer[i])
            //             nhits++;
            //     }
            //     cout << "\t\tMissing hits in layer, only got " << nhits << " layers" << endl;
            // }
        } // NOTE END LOOP OVER PARTICLES
        // // loop over ms_bitID and find duplicates
        // for (auto j = 0U; j < ms_bitID.GetSize(); ++j)
        // {
        //     if(ms_bitID[j] <= 0) continue;
        //     for (auto k = j + 1; k < ms_bitID.GetSize(); ++k)
        //     {
        //         if(ms_bitID[k] <= 0) continue;
        //         if (ms_bitID[j] == ms_bitID[k])
        //         {
        //             cout << "Entry: " << curr_entry << endl;
        //             cout << "\tDuplicate: " << ms_bitID[j] << endl;

        //     int station = (ms_bitID[j] >> 8) & 0x7; //3 bit
        //     int module = (ms_bitID[j] >> 11) & 0xF; //4 bit
        //     int layer = (ms_bitID[j] >> 15) & 0x7;  //3 bit
        //     int segment = (ms_bitID[j] >> 18) & 0xF; //4 bit
        //     int bar = (ms_bitID[j] >> 22) & 0xFF; // 8 bit
        //     int isfiber = (ms_bitID[j] >> 30) & 0x3; //2 bit
        //         cout << "\t\tStation: " << station << " Module: " << module << " Layer: " << layer << " Segment: " << segment << " Bar: " << bar << endl;
        //         }
        //     }
        // }
        // loop over pt and eta to fill the acceptance histogram
        // for (auto j = 0U; j < pt.GetSize(); ++j)
        // {
        //     hAcc->Fill(1 / (pt[j] / 1000), eta[j]);
        //     // cout << "pt: " << pt[j] << " eta: " << eta[j] << endl;
        // }
    }
    cout << "Max bar: " << maxbar << endl;
    // normalize the average vectors
    for (int i = 0; i < nModules; i++)
    {
        vIncidentAverage[i].SetMag(1);
    }
    // normalice the acceptance histogram
    hAcc->Scale(1 / hAcc->GetEntries());
    hAcc_ft->Scale(1 / hAcc_ft->GetEntries());
    hAcc_ms->Scale(1 / hAcc_ms->GetEntries());
    hAcc_ms_noFirstModule->Scale(1 / hAcc_ms_noFirstModule->GetEntries());
    // print the average vectors
    for (int i = 0; i < nModules; i++)
    {
        std::cout << "module " << i + 1 << " average vector: " << vIncidentAverage[i].X() << " " << vIncidentAverage[i].Y() << " " << vIncidentAverage[i].Z() << std::endl;
    }

    // loop over x bins of h2 and fill graph with mean values along y axis
    TGraph *gMeanValues = new TGraph();
    for (int i = 1; i < h2->GetNbinsX(); i++)
    {
        float mean = 0;
        TH1D *h1 = h2->ProjectionY("h1", i, i + 1);
        mean = h1->GetMean();
        gMeanValues->SetPoint(i, h2->GetXaxis()->GetBinCenter(i), mean);
    }
    // remove 0 values from the graph for a nicer plot
    for (int i = gMeanValues->GetN(); i >= 0; i--)
    {
        if (gMeanValues->GetY()[i] < 1)
        {
            gMeanValues->RemovePoint(i);
        }
    }
    // do the same for the horizontal angle plot
    TGraph *gMeanValues_hor = new TGraph();
    for (int i = 1; i < h2_hor->GetNbinsX(); i++)
    {
        float mean = 0;
        TH1D *h1 = h2_hor->ProjectionY("h1_hor", i, i + 1);
        mean = h1->GetMean();
        gMeanValues_hor->SetPoint(i, h2_hor->GetXaxis()->GetBinCenter(i), mean);
    }
    for (int i = gMeanValues_hor->GetN(); i >= 0; i--)
    {
        if (gMeanValues_hor->GetY()[i] < 1)
        {
            gMeanValues_hor->RemovePoint(i);
        }
    }
    // do the same for alternative horizontal angle plot
    TGraph *gMeanValues_hor_alt = new TGraph();
    for (int i = 1; i < h2_hor_alt->GetNbinsX(); i++)
    {
        float mean = 0;
        TH1D *h1 = h2_hor_alt->ProjectionY("h1_hor_alt", i, i + 1);
        mean = h1->GetMean();
        gMeanValues_hor_alt->SetPoint(i, h2_hor_alt->GetXaxis()->GetBinCenter(i), mean);
    }
    for (int i = gMeanValues_hor_alt->GetN(); i >= 0; i--)
    {
        if (gMeanValues_hor_alt->GetY()[i] < 1)
        {
            gMeanValues_hor_alt->RemovePoint(i);
        }
    }

    // do the same for the vertical angle plot
    TGraph *gMeanValues_ver = new TGraph();
    for (int i = 1; i < h2_ver->GetNbinsX(); i++)
    {
        float mean = 0;
        TH1D *h1 = h2_ver->ProjectionY("h1_ver", i, i + 1);
        mean = h1->GetMean();
        gMeanValues_ver->SetPoint(i, h2_ver->GetXaxis()->GetBinCenter(i), mean);
    }
    for (int i = gMeanValues_ver->GetN(); i >= 0; i--)
    {
        if (gMeanValues_ver->GetY()[i] < 1)
        {
            gMeanValues_ver->RemovePoint(i);
        }
    }

    // plot nHitsModule
    TCanvas *c1_nHitsModule = new TCanvas("c1_nHitsModule", "c1_nHitsModule", 800, 600);
    gStyle->SetOptStat(0);
    hnHitsModulePerStation->GetXaxis()->SetTitle("Module");
    hnHitsModulePerStation->GetYaxis()->SetTitle("Number of hits");
    hnHitsModulePerStation->Draw("colz");
    c1_nHitsModule->SaveAs("hnHitsModulePerStation.pdf");

    // plot hXYft and hXYftMother
    TCanvas *c1_hXYft = new TCanvas("c1_hXYft", "c1_hXYft", 800, 600);
    gStyle->SetOptStat(0);
    hXYft->GetXaxis()->SetTitle("X [mm]");
    hXYft->GetYaxis()->SetTitle("Y [mm]");
    hXYft->Draw("colz");
    c1_hXYft->SaveAs("hXYft.pdf");

    TCanvas *c1_hXYftMother = new TCanvas("c1_hXYftMother", "c1_hXYftMother", 800, 600);
    gStyle->SetOptStat(0);
    hXYftMother->GetXaxis()->SetTitle("X [mm]");
    hXYftMother->GetYaxis()->SetTitle("Y [mm]");
    hXYftMother->Draw("colz");
    c1_hXYftMother->SaveAs("hXYftMother.pdf");

    // plot nLayerHitsPerModule
    TCanvas *c1_nLayerHitsPerModule = new TCanvas("c1_nLayerHitsPerModule", "c1_nLayerHitsPerModule", 800, 600);
    gStyle->SetOptStat(0);
    c1_nLayerHitsPerModule->SetLogz();
    hnLayerHitsPerModule->GetXaxis()->SetTitle("Module");
    hnLayerHitsPerModule->GetYaxis()->SetTitle("Layer");
    hnLayerHitsPerModule->Draw("colz");
    c1_nLayerHitsPerModule->SaveAs("hnLayerHitsPerModule.pdf");

    // plot nHitsPerParticle
    TCanvas *c1_nHitsPerParticle = new TCanvas("c1_nHitsPerParticle", "c1_nHitsPerParticle", 800, 600);
    hnHitsPerParticle->GetXaxis()->SetTitle("Number of hits");
    hnHitsPerParticle->GetYaxis()->SetTitle("Number of particles");
    hnHitsPerParticle->Draw();
    c1_nHitsPerParticle->SaveAs("hnHitsPerParticle.pdf");

    // plot nModulesHitPerParticle
    TCanvas *c1_nModulesHitPerParticle = new TCanvas("c1_nModulesHitPerParticle", "c1_nModulesHitPerParticle", 800, 600);
    hnModulesHitPerParticle->GetXaxis()->SetTitle("Number of modules hit");
    hnModulesHitPerParticle->GetYaxis()->SetTitle("Number of particles");
    hnModulesHitPerParticle->Draw();
    c1_nModulesHitPerParticle->SaveAs("hnModulesHitPerParticle.pdf");

    // plot nModulesHitPerParticleVsEta
    TCanvas *c1_nModulesHitPerParticleVsEta = new TCanvas("c1_nModulesHitPerParticleVsEta", "c1_nModulesHitPerParticleVsEta", 800, 600);
    hnModulesHitPerParticleVsEta->GetXaxis()->SetTitle("#eta");
    hnModulesHitPerParticleVsEta->GetYaxis()->SetTitle("Number of modules hit");
    hnModulesHitPerParticleVsEta->Draw("colz");
    c1_nModulesHitPerParticleVsEta->SaveAs("hnModulesHitPerParticleVsEta.pdf");

    // plot nModulesHitPerParticleVsPt
    TCanvas *c1_nModulesHitPerParticleVsPt = new TCanvas("c1_nModulesHitPerParticleVsPt", "c1_nModulesHitPerParticleVsPt", 800, 600);
    c1_nModulesHitPerParticleVsPt->SetLogz();
    hnModulesHitPerParticleVsPt->GetXaxis()->SetTitle("pT [MeV]");
    hnModulesHitPerParticleVsPt->GetYaxis()->SetTitle("Number of modules hit");
    hnModulesHitPerParticleVsPt->Draw("colz");
    c1_nModulesHitPerParticleVsPt->SaveAs("hnModulesHitPerParticleVsPt.pdf");

    // plot pTvsP_lowHorAngle
    TCanvas *c1_pTvsP_lowHorAngle = new TCanvas("c1_pTvsP_lowHorAngle", "c1_pTvsP_lowHorAngle", 800, 600);
    pTvsP_lowHorAngle->GetXaxis()->SetTitle("p [MeV]");
    pTvsP_lowHorAngle->GetYaxis()->SetTitle("pT [MeV]");
    pTvsP_lowHorAngle->Draw("colz");
    c1_pTvsP_lowHorAngle->SaveAs("pTvsP_lowHorAngle.pdf");

    // plot pTvsP_highHorAngle
    TCanvas *c1_pTvsP_highHorAngle = new TCanvas("c1_pTvsP_highHorAngle", "c1_pTvsP_highHorAngle", 800, 600);
    pTvsP_highHorAngle->GetXaxis()->SetTitle("p [MeV]");
    pTvsP_highHorAngle->GetYaxis()->SetTitle("pT [MeV]");
    pTvsP_highHorAngle->Draw("colz");
    c1_pTvsP_highHorAngle->SaveAs("pTvsP_highHorAngle.pdf");

    // plot hTimeDiffvs1pT
    TCanvas *c1_hTimeDiffvs1pT = new TCanvas("c1_hTimeDiffvs1pT", "c1_hTimeDiffvs1pT", 800, 600);
    hTimeDiffvs1pT->GetXaxis()->SetTitle("1/pT [1/MeV]");
    hTimeDiffvs1pT->GetYaxis()->SetTitle("Time difference between layer 3 and layer 0 [ns]");
    hTimeDiffvs1pT->Draw("colz");
    c1_hTimeDiffvs1pT->SaveAs("hTimeDiffvs1pT.pdf");

    // plot deltaYvsDeltaZ
    TCanvas *c1_deltaYvsDeltaZ = new TCanvas("c1_deltaYvsDeltaZ", "c1_deltaYvsDeltaZ", 800, 600);
    deltaYvsDeltaZ->GetXaxis()->SetTitle("deltaY [mm]");
    deltaYvsDeltaZ->GetYaxis()->SetTitle("deltaZ [mm]");
    deltaYvsDeltaZ->Draw("colz");
    c1_deltaYvsDeltaZ->SaveAs("deltaYvsDeltaZ.pdf");

    // plot deltaXvsDeltaZ
    TCanvas *c1_deltaXvsDeltaZ = new TCanvas("c1_deltaXvsDeltaZ", "c1_deltaXvsDeltaZ", 800, 600);
    deltaXvsDeltaZ->GetXaxis()->SetTitle("deltaX [mm]");
    deltaXvsDeltaZ->GetYaxis()->SetTitle("deltaZ [mm]");
    deltaXvsDeltaZ->Draw("colz");
    c1_deltaXvsDeltaZ->SaveAs("deltaXvsDeltaZ.pdf");

    // plot deltaYvsSegment
    TCanvas *c1_deltaYvsSegment = new TCanvas("c1_deltaYvsSegment", "c1_deltaYvsSegment", 800, 600);
    deltaYvsSegment->GetXaxis()->SetTitle("Segment");
    deltaYvsSegment->GetYaxis()->SetTitle("deltaY [mm]");
    deltaYvsSegment->Draw("colz");
    c1_deltaYvsSegment->SaveAs("deltaYvsSegment.pdf");

    // plot all 4 deltaXvsDeltaZPerLayer on one canvas
    TCanvas *c1_deltaXvsDeltaZPerLayer = new TCanvas("c1_deltaXvsDeltaZPerLayer", "c1_deltaXvsDeltaZPerLayer", 800, 600);
    c1_deltaXvsDeltaZPerLayer->Divide(2, 2);
    for (int i = 0; i < 4; i++)
    {
        c1_deltaXvsDeltaZPerLayer->cd(i + 1);
        deltaXvsDeltaZPerLayer[i]->GetXaxis()->SetTitle("deltaX [mm]");
        deltaXvsDeltaZPerLayer[i]->GetYaxis()->SetTitle("deltaZ [mm]");
        deltaXvsDeltaZPerLayer[i]->Draw("colz");
    }
    c1_deltaXvsDeltaZPerLayer->SaveAs("deltaXvsDeltaZPerLayer.pdf");

    // plot all deltaXvsDeltaZPerModule on one canvas
    TCanvas *c1_deltaXvsDeltaZPerModule = new TCanvas("c1_deltaXvsDeltaZPerModule", "c1_deltaXvsDeltaZPerModule", 800, 600);
    c1_deltaXvsDeltaZPerModule->Divide(4, 4);
    for (int i = 0; i < 16; i++)
    {
        c1_deltaXvsDeltaZPerModule->cd(i + 1);
        deltaXvsDeltaZPerModule[i]->GetXaxis()->SetTitle("deltaX [mm]");
        deltaXvsDeltaZPerModule[i]->GetYaxis()->SetTitle("deltaZ [mm]");
        deltaXvsDeltaZPerModule[i]->Draw("colz");
    }
    c1_deltaXvsDeltaZPerModule->SaveAs("deltaXvsDeltaZPerModule.pdf");

    // create a canvas and draw the histogram
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    h2->GetXaxis()->SetRangeUser(minZ, maxZ);
    h2->GetYaxis()->SetRangeUser(0, 180);
    h2->GetXaxis()->SetTitle("z [mm]");
    h2->GetYaxis()->SetTitle("angle [deg]");
    h2->Draw("colz");
    // create line at 90 degrees
    TLine *l = new TLine(minZ, 90, maxZ, 90);
    l->SetLineColor(kRed + 2);
    l->SetLineWidth(2);
    l->Draw();
    // create lines at modSpacing mm intervals for the modules
    for (int i = 1; i < nModules; i++)
    {
        TLine *l2 = new TLine(minZ + i * modSpacing, 10, minZ + i * modSpacing, 105);
        l2->SetLineColor(kGreen + 2);
        l2->SetLineWidth(2);
        l2->Draw();
    }
    // draw the mean values
    gMeanValues->SetLineColor(kBlack);
    gMeanValues->SetLineWidth(2);
    gMeanValues->Draw("same,l");
    // save the plot
    c1->SaveAs("angleplot.pdf");

    // plot h2_alt
    TCanvas *c1_alt = new TCanvas("c1_alt", "c1_alt", 800, 600);
    gStyle->SetOptStat(0);
    h2_alt->GetXaxis()->SetRangeUser(minZ, maxZ);
    h2_alt->GetYaxis()->SetRangeUser(0, 180);
    h2_alt->GetXaxis()->SetTitle("z [mm]");
    h2_alt->GetYaxis()->SetTitle("angle [deg]");
    h2_alt->Draw("colz");
    // create line at 90 degrees
    TLine *l_alt = new TLine(minZ, 90, maxZ, 90);
    l_alt->SetLineColor(kRed + 2);
    l_alt->SetLineWidth(2);
    l_alt->Draw();
    // // create lines at modSpacing mm intervals for the modules
    // for (int i = 1; i < nModules; i++)
    // {
    //     TLine *l2 = new TLine(minZ + i * modSpacing, 10, minZ + i * modSpacing, 105);
    //     l2->SetLineColor(kGreen + 2);
    //     l2->SetLineWidth(2);
    //     l2->Draw();
    // }
    // save the plot
    c1_alt->SaveAs("angleplot_alt.pdf");

    // make horizontal angle plot
    TCanvas *c1_hor = new TCanvas("c1_hor", "c1_hor", 800, 600);
    gStyle->SetOptStat(0);
    h2_hor->GetXaxis()->SetRangeUser(minZ, maxZ);
    h2_hor->GetYaxis()->SetRangeUser(0, 180);
    h2_hor->GetXaxis()->SetTitle("z [mm]");
    h2_hor->GetYaxis()->SetTitle("angle [deg]");
    h2_hor->Draw("colz");
    // create line at 90 degrees
    TLine *l_hor = new TLine(minZ, 90, maxZ, 90);
    l_hor->SetLineColor(kRed + 2);
    l_hor->SetLineWidth(2);
    l_hor->Draw();
    // create lines at modSpacing mm intervals for the modules
    for (int i = 1; i < nModules; i++)
    {
        TLine *l2 = new TLine(minZ + i * modSpacing, 10, minZ + i * modSpacing, 105);
        l2->SetLineColor(kGreen + 2);
        l2->SetLineWidth(2);
        l2->Draw();
    }
    // draw the mean values
    gMeanValues_hor->SetLineColor(kBlack);
    gMeanValues_hor->SetLineWidth(2);
    gMeanValues_hor->Draw("same,l");
    // save the plot
    c1_hor->SaveAs("angleplot_hor.pdf");

    // make vertical angle plot
    TCanvas *c1_ver = new TCanvas("c1_ver", "c1_ver", 800, 600);
    gStyle->SetOptStat(0);
    h2_ver->GetXaxis()->SetRangeUser(minZ, maxZ);
    h2_ver->GetYaxis()->SetRangeUser(0, 180);
    h2_ver->GetXaxis()->SetTitle("z [mm]");
    h2_ver->GetYaxis()->SetTitle("angle [deg]");
    h2_ver->Draw("colz");
    // create line at 90 degrees
    TLine *l_ver = new TLine(minZ, 90, maxZ, 90);
    l_ver->SetLineColor(kRed + 2);
    l_ver->SetLineWidth(2);
    l_ver->Draw();
    // create lines at modSpacing mm intervals for the modules
    for (int i = 1; i < nModules; i++)
    {
        TLine *l2 = new TLine(minZ + i * modSpacing, 10, minZ + i * modSpacing, 105);
        l2->SetLineColor(kGreen + 2);
        l2->SetLineWidth(2);
        l2->Draw();
    }
    // draw the mean values
    gMeanValues_ver->SetLineColor(kBlack);
    gMeanValues_ver->SetLineWidth(2);
    gMeanValues_ver->Draw("same,l");
    // save the plot
    c1_ver->SaveAs("angleplot_ver.pdf");

    // make vertical alternative angle plot
    TCanvas *c1_ver_alt = new TCanvas("c1_ver_alt", "c1_ver_alt", 800, 600);
    gStyle->SetOptStat(0);
    h2_ver_alt->GetXaxis()->SetRangeUser(minZ, maxZ);
    h2_ver_alt->GetYaxis()->SetRangeUser(0, 180);
    h2_ver_alt->GetXaxis()->SetTitle("z [mm]");
    h2_ver_alt->GetYaxis()->SetTitle("angle [deg]");
    h2_ver_alt->Draw("colz");
    // create line at 90 degrees
    TLine *l_ver_alt = new TLine(minZ, 90, maxZ, 90);
    l_ver_alt->SetLineColor(kRed + 2);
    l_ver_alt->SetLineWidth(2);
    l_ver_alt->Draw();
    // // create lines at modSpacing mm intervals for the modules
    // for (int i = 1; i < nModules; i++)
    // {
    //     TLine *l2 = new TLine(minZ + i * modSpacing, 10, minZ + i * modSpacing, 105);
    //     l2->SetLineColor(kGreen + 2);
    //     l2->SetLineWidth(2);
    //     l2->Draw();
    // }
    // save the plot
    c1_ver_alt->SaveAs("angleplot_ver_alt.pdf");

    // make horizontal alternative angle plot
    TCanvas *c1_hor_alt = new TCanvas("c1_hor_alt", "c1_hor_alt", 800, 600);
    gStyle->SetOptStat(0);
    h2_hor_alt->GetXaxis()->SetRangeUser(minZ, maxZ);
    h2_hor_alt->GetYaxis()->SetRangeUser(0, 180);
    h2_hor_alt->GetXaxis()->SetTitle("z [mm]");
    h2_hor_alt->GetYaxis()->SetTitle("angle [deg]");
    h2_hor_alt->Draw("colz");
    // create line at 90 degrees
    TLine *l_hor_alt = new TLine(minZ, 90, maxZ, 90);
    l_hor_alt->SetLineColor(kRed + 2);
    l_hor_alt->SetLineWidth(2);
    l_hor_alt->Draw();
    // // create lines at modSpacing mm intervals for the modules
    // for (int i = 1; i < nModules; i++)
    // {
    //     TLine *l2 = new TLine(minZ + i * modSpacing, 10, minZ + i * modSpacing, 105);
    //     l2->SetLineColor(kGreen + 2);
    //     l2->SetLineWidth(2);
    //     l2->Draw();
    // }

    // draw the mean values
    gMeanValues_hor_alt->SetLineColor(kBlack);
    gMeanValues_hor_alt->SetLineWidth(2);
    gMeanValues_hor_alt->Draw("same,l");
    // save the plot
    c1_hor_alt->SaveAs("angleplot_hor_alt.pdf");

    // make vertical angle plot
    TCanvas *c1_ver_y = new TCanvas("c1_ver_y", "c1_ver_y", 800, 600);
    gStyle->SetOptStat(0);
    h2_ver_y->GetXaxis()->SetRangeUser(-1500, 1500);
    h2_ver_y->GetYaxis()->SetRangeUser(0, 180);
    h2_ver_y->GetXaxis()->SetTitle("y [mm]");
    h2_ver_y->GetYaxis()->SetTitle("angle [deg]");
    h2_ver_y->Draw("colz");
    // create line at 90 degrees
    TLine *l_ver_y = new TLine(-1500, 90, 1500, 90);
    l_ver_y->SetLineColor(kRed + 2);
    l_ver_y->SetLineWidth(2);
    l_ver_y->Draw();
    // create lines at modSpacing mm intervals for the modules
    for (int i = 1; i < nModules; i++)
    {
        TLine *l2 = new TLine(-1500, 10, 1500, 10);
        l2->SetLineColor(kGreen + 2);
        l2->SetLineWidth(2);
        l2->Draw();
    }
    // draw the mean values
    // gMeanValues->SetLineColor(kBlack);
    // gMeanValues->SetLineWidth(2);
    // gMeanValues->Draw("same,l");
    // save the plot
    c1_ver_y->SaveAs("angleplot_ver_y.pdf");

    // make 7 panel canvas and plot each module individually
    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 1600);
    c2->Divide(4, 4);
    for (int i = 0; i < nModules; i++)
    {
        c2->cd(i + 1);
        hSeg[i]->GetXaxis()->SetTitle("z [mm]");
        hSeg[i]->GetYaxis()->SetTitle("angle [deg]");
        hSeg[i]->Draw("colz");
    }
    c2->SaveAs("angleplot_modules.pdf");

    // make 7 panel canvas and plot each modules projection individually
    TCanvas *c3 = new TCanvas("c3", "c3", 1200, 1600);
    c3->Divide(4, 4);
    TH1D *h2_projections[nModules];
    for (int i = 1; i < nModules + 1; i++)
    {
        c3->cd(i);
        h2_projections[i - 1] = h2->ProjectionY(Form("h2_projections_%d", i), h2->GetXaxis()->FindBin(minZ + (i - 1) * modSpacing), h2->GetXaxis()->FindBin(minZ + i * modSpacing));
        h2_projections[i - 1]->SetTitle(Form("module %d: Mean: %.2f", i, h2_projections[i - 1]->GetMean()));
        // normalize the projection
        h2_projections[i - 1]->Scale(1 / h2_projections[i - 1]->Integral());
        h2_projections[i - 1]->GetXaxis()->SetRangeUser(0, 120);
        h2_projections[i - 1]->Draw();
    }
    c3->SaveAs("angleplot_modules_projection.pdf");

    // plot all 7 modules projections on the same canvas with different colors
    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    TH1D *h2_projectionsDummy = (TH1D *)h2_projections[0]->Clone("h2_projectionsDummy");
    h2_projectionsDummy->GetXaxis()->SetTitle("angle [deg]");
    h2_projectionsDummy->GetYaxis()->SetTitle("counts");
    h2_projectionsDummy->SetTitle("Incident angle of particles on the Magnet Station");
    h2_projectionsDummy->SetLineColor(kBlack);
    h2_projectionsDummy->SetLineWidth(2);
    // h2_projectionsDummy->GetYaxis()->SetRangeUser(0.9, h2_projectionsDummy->GetMaximum() * 10);
    h2_projectionsDummy->Draw("l");
    TLegend *leg1 = new TLegend(0.1, 0.6, 0.4, 0.9);
    leg1->SetHeader("module");
    leg1->AddEntry(h2_projectionsDummy, Form("1: Mean %.0f", h2_projectionsDummy->GetMean()), "l");
    for (int i = 1; i < nModules; i++)
    {
        h2_projections[i]->SetLineColor(i + 1);
        h2_projections[i]->SetLineWidth(2);
        leg1->AddEntry(h2_projections[i], Form("%d: Mean %.0f", i + 1, h2_projections[i]->GetMean()), "l");
        h2_projections[i]->Draw("l,same");
    }
    leg1->Draw();
    c4->SaveAs("angleplot_modules_projection_all.pdf");

    // plot all 7 modules momentum on the same canvas with different colors
    TCanvas *c6 = new TCanvas("c6", "c6", 800, 600);
    c6->SetLogy();
    TH1F *hSeg_pDummy = (TH1F *)hSeg_p[0]->Clone("hSeg_pDummy");
    hSeg_pDummy->GetXaxis()->SetTitle("p [MeV/c]");
    hSeg_pDummy->SetTitle("Momentum of particles on the Magnet Station");
    hSeg_pDummy->Rebin(2);
    hSeg_pDummy->SetLineColor(kBlack);
    hSeg_pDummy->SetLineWidth(2);
    hSeg_pDummy->GetYaxis()->SetRangeUser(0.9, hSeg_pDummy->GetMaximum() * 10);
    hSeg_pDummy->Draw();
    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->SetHeader("module");
    leg->AddEntry(hSeg_pDummy, Form("1: Mean %.0f", hSeg_pDummy->GetMean()), "l");
    for (int i = 1; i < nModules; i++)
    {
        hSeg_p[i]->SetLineColor(i + 1);
        hSeg_p[i]->SetLineWidth(2);
        hSeg_p[i]->Rebin(2);
        leg->AddEntry(hSeg_p[i], Form("%d: Mean %.0f", i + 1, hSeg_p[i]->GetMean()), "l");
        hSeg_p[i]->Draw("same");
    }
    leg->Draw();
    c6->SaveAs("angleplot_modules_momentum_all.pdf");

    // plot all 7 modules transverse momentum on the same canvas with different colors
    TCanvas *c7 = new TCanvas("c7", "c7", 800, 600);
    c7->SetLogy();
    TH1F *hSeg_ptDummy = (TH1F *)hSeg_pt[0]->Clone("hSeg_ptDummy");
    hSeg_ptDummy->GetXaxis()->SetTitle("p_{T} [MeV/c]");
    hSeg_ptDummy->SetTitle("Transverse Momentum of particles on the Magnet Station");
    hSeg_ptDummy->Rebin(2);
    hSeg_ptDummy->SetLineColor(kBlack);
    hSeg_ptDummy->SetLineWidth(2);
    hSeg_ptDummy->GetYaxis()->SetRangeUser(0.9, hSeg_ptDummy->GetMaximum() * 10);
    hSeg_ptDummy->Draw();
    TLegend *leg2 = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg2->SetHeader("module");
    leg2->AddEntry(hSeg_ptDummy, Form("1: Mean %.0f", hSeg_ptDummy->GetMean()), "l");
    for (int i = 1; i < nModules; i++)
    {
        hSeg_pt[i]->SetLineColor(i + 1);
        hSeg_pt[i]->SetLineWidth(2);
        hSeg_pt[i]->Rebin(2);
        leg2->AddEntry(hSeg_pt[i], Form("%d: Mean %.0f", i + 1, hSeg_pt[i]->GetMean()), "l");
        hSeg_pt[i]->Draw("same");
    }
    leg2->Draw();
    c7->SaveAs("angleplot_modules_transverse_momentum_all.pdf");

    // plot all 7 modules time on the same canvas with different colors
    TCanvas *c8 = new TCanvas("c8", "c8", 800, 600);
    c8->SetLogy();
    TH1F *hSeg_timeDummy = (TH1F *)hSeg_time[0]->Clone("hSeg_timeDummy");
    hSeg_timeDummy->GetXaxis()->SetTitle("time [ns]");
    hSeg_timeDummy->GetYaxis()->SetRangeUser(0.9, hSeg_timeDummy->GetMaximum() * 10);
    hSeg_timeDummy->SetTitle("Time of particles on the Magnet Station");
    hSeg_timeDummy->Rebin(2);
    hSeg_timeDummy->SetLineColor(kBlack);
    hSeg_timeDummy->SetLineWidth(2);
    hSeg_timeDummy->Draw();
    TLegend *leg3 = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg3->SetHeader("module");
    leg3->AddEntry(hSeg_timeDummy, Form("1: Mean %.0f", hSeg_timeDummy->GetMean()), "l");
    for (int i = 1; i < nModules; i++)
    {
        hSeg_time[i]->SetLineColor(i + 1);
        hSeg_time[i]->SetLineWidth(2);
        hSeg_time[i]->Rebin(2);
        leg3->AddEntry(hSeg_time[i], Form("%d: Mean %.0f", i + 1, hSeg_time[i]->GetMean()), "l");
        hSeg_time[i]->Draw("same");
    }
    leg3->Draw();
    c8->SaveAs("angleplot_modules_time_all.pdf");

    // plot acceptance
    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
    hAcc->GetXaxis()->SetTitle("1/p_{T} [1/MeV/c]");
    hAcc->GetYaxis()->SetTitle("#eta");
    hAcc->Draw("colz");
    c5->SaveAs("angleplot_acceptance.pdf");

    // plot acceptance ft
    TCanvas *c5_ft = new TCanvas("c5_ft", "c5_ft", 800, 600);
    hAcc_ft->GetXaxis()->SetTitle("1/p_{T} [1/MeV/c]");
    hAcc_ft->GetYaxis()->SetTitle("#eta");
    hAcc_ft->Draw("colz");
    c5_ft->SaveAs("angleplot_acceptance_ft.pdf");

    // plot acceptance ms
    TCanvas *c5_ms = new TCanvas("c5_ms", "c5_ms", 800, 600);
    hAcc_ms->GetXaxis()->SetTitle("1/p_{T} [1/MeV/c]");
    hAcc_ms->GetYaxis()->SetTitle("#eta");
    hAcc_ms->Draw("colz");
    c5_ms->SaveAs("angleplot_acceptance_ms.pdf");

    // plot acceptance ms no first module
    TCanvas *c5_ms_noFirstModule = new TCanvas("c5_ms_noFirstModule", "c5_ms_noFirstModule", 800, 600);
    hAcc_ms_noFirstModule->GetXaxis()->SetTitle("1/p_{T} [1/MeV/c]");
    hAcc_ms_noFirstModule->GetYaxis()->SetTitle("#eta");
    hAcc_ms_noFirstModule->Draw("colz");
    c5_ms_noFirstModule->SaveAs("angleplot_acceptance_ms_noFirstModule.pdf");

    // calculate ratio of acceptance ms and acceptance ms no first module
    TH2F *hAccRatio = (TH2F *)hAcc_ms_noFirstModule->Clone("hAccRatio");
    hAccRatio->Divide(hAcc_ms);
    TCanvas *c5_AccRatio = new TCanvas("c5_AccRatio", "c5_AccRatio", 800, 600);
    hAccRatio->GetXaxis()->SetTitle("1/p_{T} [1/MeV/c]");
    hAccRatio->GetYaxis()->SetTitle("#eta");
    hAccRatio->Draw("colz");
    c5_AccRatio->SaveAs("angleplot_acceptance_ratio.pdf");

    // create output file and save the histogram
    TFile *fout = new TFile("angleplot.root", "RECREATE");
    h2->Write();
    fout->Close();
}