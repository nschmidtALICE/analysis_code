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

void makeangleplot(TString inputfile = "/home/niviths/Downloads/minimumBias_MS_MagDown_100plus.root")
// void makeangleplot(TString inputfile = "/media/niviths/local/ms_sim/angle_studies/minimumBias_MS_MagDown_20plus.root")
{
    TFile fin(inputfile, "READ");
    TTreeReader tree("ntup", &fin);

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
    TTreeReaderArray<int> ms_segment(tree, "ms_segment");

    // create 2D histogram of incident angle of particles on the Magnet Station
    TH2F *h2 = new TH2F("h2", "Incident angle of particles on the Magnet Station", 125, 3500, 7500, 60, 0, 120);
    TH2F *hAcc = new TH2F("hAcc", "Acceptance of particles on the Magnet Station", 125, 0, 20, 60, 1.8, 4.5);
    TH2F *hSeg[7] = {nullptr};
    TH1F *hSeg_p[7] = {nullptr};
    TH1F *hSeg_pt[7] = {nullptr};
    TH1F *hSeg_time[7] = {nullptr};
    for (int i = 0; i < 7; i++)
    {
        hSeg[i] = new TH2F(Form("hSeg_%d", i), Form("Incident angle of particles on the Magnet Station segment %d", i), 125, 3500, 7500, 60, 0, 120);
        hSeg_p[i] = new TH1F(Form("hSeg_p_%d", i), Form("Momentum of particles on the Magnet Station segment %d", i), 100, 0, 5000);
        hSeg_pt[i] = new TH1F(Form("hSeg_pt_%d", i), Form("Transverse momentum of particles on the Magnet Station segment %d", i), 80, 0, 1600);
        hSeg_time[i] = new TH1F(Form("hSeg_time_%d", i), Form("Time of particles on the Magnet Station segment %d", i), 100, 0, 300);
    }

    // manual rotation of the detector plane for each segment
    //  float add_rotations[7] = {0.0, 0.1, 0.3, 0.4, 0.5, 0.5, 0.6};
    //  float add_rotations[7] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1}; //used for 14
    // float add_rotations[7] = {0.6, 0.6, 0.6, 0.6, 0.7, 0.8, 0.9}; // used for 15,16
    // float add_rotations[7] = {0.7, 0.7, 0.6, 0.6, 0.7, 0.8, 0.9}; // used for 20plus
    float add_rotations[7] = {1.0, 1.0, 0.7, 0.6, 0.7, 0.8, 0.9}; // used for 15,16
    // float add_rotations[7] = {1.1, 1.0, 0.7, 0.6, 0.7, 0.8, 0.9}; //BEST SO FAR, but unfeasible
    // float add_rotations[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // create 7 vectors to store the average values of the incident angle for each segment
    TVector3 vIncidentAverage[7];

    // loop over tree
    while (tree.Next())
    {
        // loop over size of ms_px to get all the hits per event
        for (auto j = 0U; j < ms_px.GetSize(); ++j)
        {
            // require hits within MS acceptance
            if (ms_vz[j] < 3000 || ms_vz[j] > 8500)
            {
                continue;
            }
            // require certain momentum
            if (p[j] > 1000)
            {
                continue;
            }

            // print layer and segment from ms_segment
            int layer = (ms_segment[j] >> 12) & 0x0003;
            int segment = (ms_segment[j] >> 16) & 0x7;

            // if(segment!=1) continue;

            // create a vector with the direction of the particle
            TVector3 v(ms_px[j], ms_py[j], ms_pz[j]);

            //add the vector to the average vector for the segment
            if(ms_vx[j] > 0) vIncidentAverage[segment - 1] += v;

            // calculate the vector spanning the detector plane
            TVector3 v2(0, 1, 1);
            if (layer == 0 || layer == 2)
            {
                v2.RotateZ(0.34205); // tilt around z axis (makes upper and lower modules)
                v2.RotateY(-0.48);   // initial tilt around y-axis (makes modules be in line with the magnet inside)
            }
            else
            {
                v2.RotateZ(-0.34205); // tilt around z axis (makes upper and lower modules)
                v2.RotateY(0.48);     // initial tilt around y-axis (makes modules be in line with the magnet inside)
            }

            // add additional rotation depending on segment of particle hit
            // need two cases for the different sides of the detector
            if (layer == 0 || layer == 2)
                v2.RotateZ(-add_rotations[segment - 1]);
            else
                v2.RotateZ(add_rotations[segment - 1]);

            // calculate the angle in degrees between the two vectors
            float angle2 = v.Angle(v2) * 180 / 3.14159;

            // fill the histogram
            h2->Fill(ms_vz[j], angle2);
            hSeg[segment-1]->Fill(ms_vz[j], angle2);
            hSeg_p[segment-1]->Fill(p[j]);
            hSeg_pt[segment-1]->Fill(pt[j]);
            hSeg_time[segment-1]->Fill(ms_time[j]);
        }
        // loop over pt and eta to fill the acceptance histogram
        for (auto j = 0U; j < pt.GetSize(); ++j)
        {
            hAcc->Fill(1/(pt[j]/1000), eta[j]);
        }
    }
    //normalize the average vectors
    for (int i = 0; i < 7; i++)
    {
        vIncidentAverage[i].SetMag(1);
    }
    //normalice the acceptance histogram
    hAcc->Scale(1/hAcc->GetEntries());
    //print the average vectors
    for (int i = 0; i < 7; i++)
    {
        std::cout << "Segment " << i + 1 << " average vector: " << vIncidentAverage[i].X() << " " << vIncidentAverage[i].Y() << " " << vIncidentAverage[i].Z() << std::endl;
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
    h2->GetXaxis()->SetRangeUser(3550, 7049);
    h2->GetYaxis()->SetRangeUser(0, 118);
    h2->GetXaxis()->SetTitle("z [mm]");
    h2->GetYaxis()->SetTitle("angle [deg]");
    h2->Draw("colz");
    // create line at 90 degrees
    TLine *l = new TLine(3550, 90, 7049, 90);
    l->SetLineColor(kRed + 2);
    l->SetLineWidth(2);
    l->Draw();
    // create lines at 500 mm intervals for the segments
    for (int i = 1; i < 8; i++)
    {
        TLine *l2 = new TLine(3500 + i * 500, 10, 3500 + i * 500, 105);
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


    //make 7 panel canvas and plot each segment individually
    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 800);
    c2->Divide(4, 2);
    for (int i = 0; i < 7; i++)
    {
        c2->cd(i+1);
        hSeg[i]->GetXaxis()->SetTitle("z [mm]");
        hSeg[i]->GetYaxis()->SetTitle("angle [deg]");
        hSeg[i]->Draw("colz");
    }
    c2->SaveAs("angleplot_segments.pdf");

    //make 7 panel canvas and plot each segments projection individually
    TCanvas *c3 = new TCanvas("c3", "c3", 1200, 800);
    c3->Divide(4, 2);
    TH1D* h2_projections[7];
    for (int i = 1; i < 8; i++)
    {
        c3->cd(i);
        h2_projections[i-1] = h2->ProjectionY(Form("h2_projections_%d",i), h2->GetXaxis()->FindBin(3500 + (i - 1) * 500), h2->GetXaxis()->FindBin(3500 + i * 500));
        h2_projections[i-1]->SetTitle(Form("Segment %d: Mean: %.2f", i, h2_projections[i-1]->GetMean()));
        //normalize the projection
        h2_projections[i-1]->Scale(1/h2_projections[i-1]->Integral());
        h2_projections[i-1]->GetXaxis()->SetRangeUser(0, 120);
        h2_projections[i-1]->Draw();
    }
    c3->SaveAs("angleplot_segments_projection.pdf");

    // plot all 7 segments projections on the same canvas with different colors
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
    leg1->SetHeader("Segment");
    leg1->AddEntry(h2_projectionsDummy, Form("1: Mean %.0f", h2_projectionsDummy->GetMean()), "l");
    for (int i = 1; i < 7; i++)
    {
        h2_projections[i]->SetLineColor(i + 1);
        h2_projections[i]->SetLineWidth(2);
        leg1->AddEntry(h2_projections[i], Form("%d: Mean %.0f", i + 1, h2_projections[i]->GetMean()), "l");
        h2_projections[i]->Draw("l,same");
    }
    leg1->Draw();
    c4->SaveAs("angleplot_segments_projection_all.pdf");

    // plot all 7 segments momentum on the same canvas with different colors
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
    leg->SetHeader("Segment");
    leg->AddEntry(hSeg_pDummy, Form("1: Mean %.0f", hSeg_pDummy->GetMean()), "l");
    for (int i = 1; i < 7; i++)
    {
        hSeg_p[i]->SetLineColor(i + 1);
        hSeg_p[i]->SetLineWidth(2);
        hSeg_p[i]->Rebin(2);
        leg->AddEntry(hSeg_p[i], Form("%d: Mean %.0f", i + 1, hSeg_p[i]->GetMean()), "l");
        hSeg_p[i]->Draw("same");
    }
    leg->Draw();
    c6->SaveAs("angleplot_segments_momentum_all.pdf");

    // plot all 7 segments transverse momentum on the same canvas with different colors
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
    leg2->SetHeader("Segment");
    leg2->AddEntry(hSeg_ptDummy, Form("1: Mean %.0f", hSeg_ptDummy->GetMean()), "l");
    for (int i = 1; i < 7; i++)
    {
        hSeg_pt[i]->SetLineColor(i + 1);
        hSeg_pt[i]->SetLineWidth(2);
        hSeg_pt[i]->Rebin(2);
        leg2->AddEntry(hSeg_pt[i], Form("%d: Mean %.0f", i + 1, hSeg_pt[i]->GetMean()), "l");
        hSeg_pt[i]->Draw("same");
    }
    leg2->Draw();
    c7->SaveAs("angleplot_segments_transverse_momentum_all.pdf");

    // plot all 7 segments time on the same canvas with different colors
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
    leg3->SetHeader("Segment");
    leg3->AddEntry(hSeg_timeDummy, Form("1: Mean %.0f", hSeg_timeDummy->GetMean()), "l");
    for (int i = 1; i < 7; i++)
    {
        hSeg_time[i]->SetLineColor(i + 1);
        hSeg_time[i]->SetLineWidth(2);
        hSeg_time[i]->Rebin(2);
        leg3->AddEntry(hSeg_time[i], Form("%d: Mean %.0f", i + 1, hSeg_time[i]->GetMean()), "l");
        hSeg_time[i]->Draw("same");
    }
    leg3->Draw();
    c8->SaveAs("angleplot_segments_time_all.pdf");

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