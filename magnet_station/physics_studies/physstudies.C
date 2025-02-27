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

void physstudies(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_120plus.root")
// void makeangleplot(TString inputfile = "/media/niviths/local/ms_sim/angle_studies/minimumBias_MS_MagDown_20plus.root")
{
    TFile fin(inputfile, "READ");
    TTreeReader tree("ntup", &fin);
    double maxZ = 8500;

    // get p, py, pz and ms_vx, ms_vy, ms_vz, as well as p and ms_time from the tree
    TTreeReaderArray<float> p(tree, "p");
    TTreeReaderArray<float> py(tree, "py");
    TTreeReaderArray<float> pz(tree, "pz");
    TTreeReaderArray<float> ms_vx(tree, "ms_vx");
    TTreeReaderArray<float> ms_vy(tree, "ms_vy");
    TTreeReaderArray<float> ms_vz(tree, "ms_vz");
    // TTreeReaderArray<float> p(tree, "p");
    TTreeReaderArray<float> pt(tree, "pt");
    TTreeReaderArray<float> eta(tree, "eta");
    TTreeReaderArray<int> pid(tree, "pid");
    TTreeReaderArray<float> ms_time(tree, "ms_time");
    TTreeReaderArray<int> nMShits(tree, "nMShits");
    TTreeReaderArray<int> nFThits(tree, "nFThits");
    TTreeReaderArray<int> ms_bitID(tree, "ms_segment");
    const int nModules = 9;
    // create 2D histogram of incident angle of particles on the Magnet Station
    // define a log scale for the x axis from 0.8 to 200
    double bins[500];
    for (int i = 0; i < 500; i++)
    {
        bins[i] = pow(10, 0.1 + i * 0.004);
        cout << "bin: " << bins[i] << endl;
    }
    //use the log scale for the x axis
    TH1F *h_eeAcc_FTMS = new TH1F("h_eeAcc", "Dielectron mass", 499, bins);
    TH1F *h_eeAcc_MS = new TH1F("h_eeAcc_MS", "Dielectron mass MS", 100, 0.8, 200); 
    TH1F *h_eeAcc_FT = new TH1F("h_eeAcc_FT", "Dielectron mass FT", 100, 0.8, 200); 
    // TH1F *h_eeAcc_FTMS = new TH1F("h_eeAcc_FTMS", "Dielectron mass FT+MS", 100, 0.8, 200); 
    
    
    int curr_entry = 0;
    // loop over tree
    while (tree.Next())
    {
        cout << "Entry: " << curr_entry << endl;
        curr_entry++;
        // loop over size of p to get all the hits per event
        for (auto j = 0U; j < pid.GetSize(); ++j)
        {
            // check PID to be electron
            if (abs(pid[j]) == 11)
            {
                TVector3 vP1(p[j], py[j], pz[j]);
                //loop over particles to get the second particle
                for (auto k = j + 1; k < pid.GetSize(); ++k)
                {
                    // check PID to be electron
                    if (pid[k] == -pid[j])
                    {
                        // create vector for momentum and position
                        TVector3 vP2(p[k], py[k], pz[k]);
                        // calculate the invariant mass from the two electrons
                        double mass = sqrt(2 * p[j] * p[k] * (1 - cos(vP1.Angle(vP2))));
                        cout << "mass: " << mass << endl;
                        // fill the histogram
                        h_eeAcc_FTMS->Fill(mass);
                    }
                }
            }
        }
    }

    // create canvas and plot the h_eeAcc_FTMS histogram
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    c->SetLogx();
    c->SetLogy();
    h_eeAcc_FTMS->GetXaxis()->SetTitle("m_{ee} [MeV/c^{2}]");
    h_eeAcc_FTMS->GetYaxis()->SetTitle("counts");
    h_eeAcc_FTMS->Draw();
    c->SaveAs("dielectron_mass_FTMS.pdf");
    

    // create output file and save the histogram
    TFile *fout = new TFile("dielectron.root", "RECREATE");
    h_eeAcc_FTMS->Write();
    fout->Close();
}