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

void makeangleplot_new(TString inputfile = "/home/niviths/Downloads/magnetStationSims/minimumBias_MS_MagDown_400plus.root")
// void makeangleplot(TString inputfile = "/media/niviths/local/ms_sim/angle_studies/minimumBias_MS_MagDown_20plus.root")
{
    TFile fin(inputfile, "READ");
    TTreeReader tree("ntup", &fin);
    double modSpacing = 460;
    double minZ = 3500+2*500;
    double maxZ = 7700;

    // get ms_px, ms_py, ms_pz and ms_vx, ms_vy, ms_vz, as well as p and ms_time from the tree
    TTreeReaderArray<float> ms_px(tree, "ms_px");
    TTreeReaderArray<float> ms_py(tree, "ms_py");
    TTreeReaderArray<float> ms_pz(tree, "ms_pz");
    TTreeReaderArray<float> ms_vx(tree, "ms_vx");
    TTreeReaderArray<float> ms_vy(tree, "ms_vy");
    TTreeReaderArray<float> ms_vz(tree, "ms_vz");
    TTreeReaderArray<float> p(tree, "p");
    TTreeReaderArray<float> pt(tree, "pt");
    TTreeReaderArray<float> eta(tree, "eta");
    TTreeReaderArray<float> ms_time(tree, "ms_time");
    TTreeReaderArray<int> ms_bitID(tree, "ms_segment");
    const int nModules = 7;
    // create 2D histogram of incident angle of particles on the Magnet Station
    TH2F *h2 = new TH2F("h2", "Incident angle of particles on the Magnet Station", 125, minZ, maxZ, 60, 0, 120);
    TH2F *hAcc = new TH2F("hAcc", "Acceptance of particles on the Magnet Station", 125, 0, 20, 60, 1.8, 4.5);
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

    // manual rotation of the detector plane for each module
    //  float add_rotations[nModules] = {0.0, 0.1, 0.3, 0.4, 0.5, 0.5, 0.6};
    //  float add_rotations[nModules] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1}; //used for 14
    // float add_rotations[nModules] = {0.6, 0.6, 0.6, 0.6, 0.7, 0.8, 0.9}; // used for 15,16
    // float add_rotations[nModules] = {0.7, 0.7, 0.6, 0.6, 0.7, 0.8, 0.9}; // used for 20plus
    // float add_rotations[nModules] = {0.7, 0.7, 0.7, 0.8, 0.9, 1.0, 1.0}; // used for 15,16
    // float add_rotations[nModules] = {0.8, 0.8, 0.7, 0.7, 0.7, 0.8, 0.9, 1.0, 1.0}; // used for 15,16
    // float add_rotations[nModules] = {1.1, 1.0, 0.7, 0.6, 0.7, 0.8, 0.9}; //BEST SO FAR, but unfeasible
    float add_rotations[nModules] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // float add_rotations[nModules] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // create 7 vectors to store the average values of the incident angle for each module
    TVector3 vIncidentAverage[nModules];
    int curr_entry = 0;
    // loop over tree
    while (tree.Next())
    {
        // cout << "Entry: " << curr_entry << endl;
        curr_entry++;
        
        // loop over size of ms_px to get all the hits per event
        for (auto j = 0U; j < ms_px.GetSize(); ++j)
        {
            // require hits within MS acceptance
            if (ms_vz[j] < 3000 || ms_vz[j] > maxZ)
            {
                continue;
            }
            // require certain momentum
            if (p[j] > 5000)
            {
                continue;
            }

            // print layer and segment from ms_bitID
            // cout << "Segment: " << std::bitset<32>(ms_bitID[j]>>8) << endl;
            //8 bit = 255
            //7 bit = 127
            //6 bit = 63
            //5 bit = 31
            //4 bit = 15
            //3 bit = 7

            int station = (ms_bitID[j] >> 8) & 0x7; //3 bit
            int module = (ms_bitID[j] >> 11) & 0xF; //4 bit
            int layer = (ms_bitID[j] >> 15) & 0x7;  //3 bit
            int segment = (ms_bitID[j] >> 18) & 0xF; //4 bit
            int bar = (ms_bitID[j] >> 22) & 0xFF; // 8 bit
            int isfiber = (ms_bitID[j] >> 30) & 0x3; //2 bit
            // cout << "Station: " << station << " Module: " << module << " Layer: " << layer << " Segment: " << segment << " Bar: " << bar << endl;

            // if(segment!=1) continue;

            // create a vector with the direction of the particle
            TVector3 v(ms_px[j], ms_py[j], ms_pz[j]);

            //add the vector to the average vector for the segment
            if(ms_vx[j] > 0) vIncidentAverage[module] += v;

            // calculate the vector spanning the detector plane
            TVector3 v2(0, 1, 1);
            if (layer == 0 || layer == 2)
            {
                v2.RotateZ(0.34205); // tilt around z axis (makes upper and lower modules)
                v2.RotateY(-0.37);   // initial tilt around y-axis (makes modules be in line with the magnet inside)
            }
            else
            {
                v2.RotateZ(-0.34205); // tilt around z axis (makes upper and lower modules)
                v2.RotateY(0.37);     // initial tilt around y-axis (makes modules be in line with the magnet inside)
            }

            // add additional rotation depending on module of particle hit
            // need two cases for the different sides of the detector
            if (layer == 0 || layer == 2)
                v2.RotateZ(-add_rotations[module]);
            else
                v2.RotateZ(add_rotations[module]);

            // calculate the angle in degrees between the two vectors
            float angle2 = v.Angle(v2) * 180 / 3.14159;

            // fill the histogram
            h2->Fill(ms_vz[j], angle2);
            hSeg[module]->Fill(ms_vz[j], angle2);
            hSeg_p[module]->Fill(p[j]);
            hSeg_pt[module]->Fill(pt[j]);
            hSeg_time[module]->Fill(ms_time[j]);
        }
        // loop over ms_bitID and find duplicates
        for (auto j = 0U; j < ms_bitID.GetSize(); ++j)
        {
            if(ms_bitID[j] <= 0) continue;
            for (auto k = j + 1; k < ms_bitID.GetSize(); ++k)
            {
                if(ms_bitID[k] <= 0) continue;
                if (ms_bitID[j] == ms_bitID[k])
                {
                    // cout << "Entry: " << curr_entry << endl;
                    // cout << "\tDuplicate: " << ms_bitID[j] << endl;

            // int station = (ms_bitID[j] >> 8) & 0x7; //3 bit
            // int module = (ms_bitID[j] >> 11) & 0xF; //4 bit
            // int layer = (ms_bitID[j] >> 15) & 0x7;  //3 bit
            // int segment = (ms_bitID[j] >> 18) & 0xF; //4 bit
            // int bar = (ms_bitID[j] >> 22) & 0xFF; // 8 bit
            // int isfiber = (ms_bitID[j] >> 30) & 0x3; //2 bit
                // cout << "\t\tStation: " << station << " Module: " << module << " Layer: " << layer << " Segment: " << segment << " Bar: " << bar << endl;
                }
            }
        }
        // loop over pt and eta to fill the acceptance histogram
        for (auto j = 0U; j < pt.GetSize(); ++j)
        {
            hAcc->Fill(1/(pt[j]/1000), eta[j]);
        }
    }
    //normalize the average vectors
    for (int i = 0; i < nModules; i++)
    {
        vIncidentAverage[i].SetMag(1);
    }
    //normalice the acceptance histogram
    hAcc->Scale(1/hAcc->GetEntries());
    //print the average vectors
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

    // create a canvas and draw the histogram
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(0);
    h2->GetXaxis()->SetRangeUser(minZ, maxZ);
    h2->GetYaxis()->SetRangeUser(0, 118);
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


    //make 7 panel canvas and plot each module individually
    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 800);
    c2->Divide(4, 2);
    for (int i = 0; i < nModules; i++)
    {
        c2->cd(i+1);
        hSeg[i]->GetXaxis()->SetTitle("z [mm]");
        hSeg[i]->GetYaxis()->SetTitle("angle [deg]");
        hSeg[i]->Draw("colz");
    }
    c2->SaveAs("angleplot_modules.pdf");

    //make 7 panel canvas and plot each modules projection individually
    TCanvas *c3 = new TCanvas("c3", "c3", 1200, 800);
    c3->Divide(4, 2);
    TH1D* h2_projections[nModules];
    for (int i = 1; i < nModules+1; i++)
    {
        c3->cd(i);
        h2_projections[i-1] = h2->ProjectionY(Form("h2_projections_%d",i), h2->GetXaxis()->FindBin(minZ + (i - 1) * modSpacing), h2->GetXaxis()->FindBin(minZ + i * modSpacing));
        h2_projections[i-1]->SetTitle(Form("module %d: Mean: %.2f", i, h2_projections[i-1]->GetMean()));
        //normalize the projection
        h2_projections[i-1]->Scale(1/h2_projections[i-1]->Integral());
        h2_projections[i-1]->GetXaxis()->SetRangeUser(0, 120);
        h2_projections[i-1]->Draw();
    }
    c3->SaveAs("angleplot_modules_projection.pdf");

    // plot all 7 modules projections on the same canvas with different colors
    TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
    TH1D* h2_projectionsDummy = (TH1D*)h2_projections[0]->Clone("h2_projectionsDummy");
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
    TH1F* hSeg_pDummy = (TH1F*)hSeg_p[0]->Clone("hSeg_pDummy");
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
    TH1F *hSeg_ptDummy = (TH1F*)hSeg_pt[0]->Clone("hSeg_ptDummy");
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
    TH1F *hSeg_timeDummy = (TH1F*)hSeg_time[0]->Clone("hSeg_timeDummy");
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


    // create output file and save the histogram
    TFile *fout = new TFile("angleplot.root", "RECREATE");
    h2->Write();
    fout->Close();
}