#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TKey.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>

void plot_jetnConst(const char* inRoot = "jetnConst.root",
                    const char* outPrefix = "mcProd_compare") {
    TFile fin(inRoot, "READ");
    if (fin.IsZombie()) {
        std::cerr << "Failed to open input file: " << inRoot << std::endl;
        return;
    }

    std::vector<TH1*> h_jetn;
    std::vector<TH1*> h_jetPt;
    std::vector<TH1*> h_tagPt;
    std::vector<TH1*> h_tagZ;
    TIter next(fin.GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();
        if (!obj) continue;
        if (obj->InheritsFrom("TH1")) {
            TH1 *h = (TH1*)obj;
            std::string name(h->GetName());
            if (name.rfind("h_jetnConst_", 0) == 0) h_jetn.push_back(h);
            else if (name.rfind("h_jetPt_", 0) == 0) h_jetPt.push_back(h);
            else if (name.rfind("h_tagPt_", 0) == 0) h_tagPt.push_back(h);
            else if (name.rfind("h_tagZ_", 0) == 0) h_tagZ.push_back(h);
            else delete h;
        } else {
            delete obj;
        }
    }

    if (h_jetn.empty() && h_jetPt.empty() && h_tagPt.empty() && h_tagZ.empty()) {
        std::cerr << "No relevant histograms found in " << inRoot << std::endl;
        fin.Close();
        return;
    }

    // color palette
    const int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+1, kViolet, kAzure+7, kPink+1};
    const size_t ncolors = sizeof(colors) / sizeof(colors[0]);

    // create dated output directory: <outPrefix>_YYYY-MM-DD
    char datebuf[32];
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::strftime(datebuf, sizeof(datebuf), "%Y-%m-%d", &tm);
    std::string outDir = std::string(outPrefix) + "_" + datebuf;
    if (gSystem->AccessPathName(outDir.c_str()) != 0) {
        gSystem->mkdir(outDir.c_str(), true);
    }

    auto drawOverlay = [&](const std::vector<TH1*> &vec, const char *namePrefix, const char *titlePrefix, const std::string &outDir) {
        if (vec.empty()) return;
        double ymax = 0.0;
        for (size_t i = 0; i < vec.size(); ++i) {
            TH1 *h = vec[i];
            double entries = h->Integral();
            if (entries > 0) h->Scale(1.0 / entries);
            h->GetYaxis()->SetTitle("Normalized entries");
            double m = h->GetBinContent(h->GetMaximumBin());
            if (m > ymax) ymax = m;
        }
        if (ymax <= 0) ymax = 1.0;

        TCanvas *c = new TCanvas(Form("c_%s", namePrefix), titlePrefix, 900, 600);
        c->cd();
        bool useLogY = (std::string(namePrefix) == "jetPt") || (std::string(namePrefix) == "tagPt");
        if (useLogY) c->SetLogy(1);

        TH1 *h0 = vec[0];
        h0->SetLineWidth(2);
        h0->SetLineColor(colors[0 % ncolors]);
        // if log y requested, ensure a positive minimum
        if (useLogY) {
            h0->SetMinimum(1e-6);
            h0->SetMaximum(std::max(ymax * 10.0, 1e-3));
        } else {
            h0->SetMaximum(ymax * 1.2);
        }
        // if plotting jetPt, limit x-axis to 0..55 GeV for clarity
        std::string np(namePrefix);
        if (np == "jetPt") h0->GetXaxis()->SetRangeUser(0.0, 55.0);
        if (np == "tagZ") h0->GetXaxis()->SetRangeUser(0.0, 1.0);
        h0->Draw("HIST");

        TLegend *leg = new TLegend(0.45, 0.6, 0.68, 0.88);
        leg->SetBorderSize(0); leg->SetFillColor(0);
        leg->SetTextSize(0.045);
        leg->AddEntry(h0, h0->GetTitle(), "l");

        for (size_t i = 1; i < vec.size(); ++i) {
            TH1 *h = vec[i];
            h->SetLineWidth(2);
            h->SetLineColor(colors[i % ncolors]);
            if (np == "jetPt") h->GetXaxis()->SetRangeUser(0.0, 55.0);
            if (useLogY) {
                h->SetMinimum(1e-8);
            }
            h->Draw("HIST SAME");
            leg->AddEntry(h, h->GetTitle(), "l");
        }
        leg->Draw();

        c->SaveAs(Form("%s/%s_%s.png", outDir.c_str(), outPrefix, namePrefix));
        c->SaveAs(Form("%s/%s_%s.pdf", outDir.c_str(), outPrefix, namePrefix));

        TFile fout(Form("%s/%s_%s.root", outDir.c_str(), outPrefix, namePrefix), "RECREATE");
        if (!fout.IsZombie()) {
            c->Write();
            for (auto h : vec) { h->SetDirectory(&fout); h->Write(); }
            fout.Close();
        }

        // reset log state
        if (useLogY) c->SetLogy(0);
        delete leg; delete c;
    };

    gStyle->SetOptStat(0);

    drawOverlay(h_jetn, "jetnConst", "jetnConst comparison", outDir);
    drawOverlay(h_jetPt, "jetPt", "jetPt comparison", outDir);
    drawOverlay(h_tagPt, "tagPt", "tagPt comparison", outDir);
    drawOverlay(h_tagZ, "tagZ", "tagZ comparison", outDir);

    // cleanup
    // for (auto h : h_jetn) delete h;
    // for (auto h : h_jetPt) delete h;
    // for (auto h : h_tagPt) delete h;
    // for (auto h : h_tagZ) delete h;

    fin.Close();
    std::cout << "Saved comparisons to: " << outPrefix << "_*.(png,pdf,root)" << std::endl;
}
