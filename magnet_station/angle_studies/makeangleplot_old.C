//this is a .root macro to plot the incident angle of the particles on the Magnet Station
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <iostream>
#include <vector>
#include <cmath>

void makeangleplot()
{
    TFile *fin = new TFile("minimumBias_MS_MagDown_10.root", "READ");
    TTree *tree = (TTree *)fin->Get("ntup");

    //get ms_px, ms_py, ms_pz and ms_vx, ms_vy, ms_vz, as well as p and ms_time from the tree
    std::vector<float> ms_px, ms_py, ms_pz, ms_vx, ms_vy, ms_vz, p, ms_time;
    tree->SetBranchAddress("ms_px", &ms_px);
    tree->SetBranchAddress("ms_py", &ms_py);
    tree->SetBranchAddress("ms_pz", &ms_pz);
    tree->SetBranchAddress("ms_vx", &ms_vx);
    tree->SetBranchAddress("ms_vy", &ms_vy);
    tree->SetBranchAddress("ms_vz", &ms_vz);
    tree->SetBranchAddress("p", &p);
    tree->SetBranchAddress("ms_time", &ms_time);

    //create 2D histogram of incident angle of particles on the Magnet Station
    TH2F *h2 = new TH2F("h2", "Incident angle of particles on the Magnet Station", 4000, 3500, 7500, 100, 0, 90);

    float add_rotations[7] = {0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};

    //loop over tree
    int nentries = tree->GetEntries();
    
    for(int i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        //loop over size of ms_px
        for(int j = 0; j < ms_vz.size(); j++)
        {
            if(ms_vz[j] < 3500 || ms_vz[j] > 7500)
            {
                continue;
            }
            //calculate the incident angle of the particles on the Magnet Station
            float angle = atan(sqrt(ms_px[j] * ms_px[j] + ms_py[j] * ms_py[j]) / ms_pz[j]) * 180 / 3.14159;

            //create a vector with the direction of the particle
            TVector3 v(ms_px[j], ms_py[j], ms_pz[j]);
            cout << "v: " << v.x() << " " << v.y() << " " << v.z() << endl;

            //calculate the angle between the particle and the vector spanning the detector plane
            TVector3 v2(0, 1, 1);
            v2.RotateZ(-0.34205);
            v2.RotateY(0.48);

            //add additional rotation depending on z position of particle hit (ms_vx, ms_vy, ms_vz)
            //make 7 cases for the 7 different z positions
            for(int k = 0; k < 7; k++)
            {
                if(ms_vz[j] > 3450 + k * 500 && ms_vz[j] < 3450 + (k+1) * 500)
                {
                    v2.RotateZ(add_rotations[k]);
                }
            }

            float angle2 = v.Angle(v2) * 180 / 3.14159;
            cout << "angle: " << angle << " angle2: " << angle2 << endl;

            //fill the histogram
            h2->Fill(ms_vz[j], angle2);
        }
    }

    //create a canvas and draw the histogram
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    h2->Draw("colz");
    c1->SaveAs("angleplot.png");

    //create output file and save the histogram
    TFile *fout = new TFile("angleplot.root", "RECREATE");
    h2->Write();
    fout->Close();

    fin->Close();
}