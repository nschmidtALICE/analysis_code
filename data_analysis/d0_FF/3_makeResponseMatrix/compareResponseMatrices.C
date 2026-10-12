#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cctype>
#include <regex>
#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TLatex.h"

static std::string toLowerCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

static std::string findPtTokenNearKeywords(const std::string& histName,
                                           const std::vector<std::string>& keywords,
                                           const std::vector<std::string>& tokens,
                                           size_t window = 30) {
    const std::string hlow = toLowerCopy(histName);

    auto windowHasKeyword = [&](size_t pos) {
        const size_t start = (pos > window) ? (pos - window) : 0;
        const size_t end = std::min(hlow.size(), pos + window);
        const std::string view = hlow.substr(start, end - start);
        for (const auto& kw : keywords) {
            if (view.find(kw) != std::string::npos) return true;
        }
        return false;
    };

    std::string best;
    for (const auto& token : tokens) {
        size_t pos = hlow.find(token);
        while (pos != std::string::npos) {
            if (windowHasKeyword(pos)) {
                if (token.size() > best.size()) best = token;
            }
            pos = hlow.find(token, pos + 1);
        }
    }
    return best;
}

static std::pair<std::string,std::string> extractGenDetRanges(const std::string& histName) {
    // returns {genRange, detRange} where each is "X-Y" or empty
    const std::string hlow = toLowerCopy(histName);
    std::pair<std::string,std::string> out{"",""};
    std::smatch m;
    // patterns like gen5-10 or gen_5-10 or gen5_10 etc.
    std::regex re_gen("(gen|true|truth|mc)[-_]?(\\d+)[-_](\\d+)", std::regex_constants::icase);
    std::regex re_det("(det|reco|rec|meas|detector)[-_]?(\\d+)[-_](\\d+)", std::regex_constants::icase);
    if (std::regex_search(hlow, m, re_gen) && m.size() >= 4) {
        out.first = m[2].str() + std::string("-") + m[3].str();
    }
    if (std::regex_search(hlow, m, re_det) && m.size() >= 4) {
        out.second = m[2].str() + std::string("-") + m[3].str();
    }
    return out;
}

static void drawJetPtLabelsIfFound(const std::string& histName,
                                  double xNDC = 0.16,
                                  double yNDC = 0.88,
                                  double textSize = 0.05,
                                  int textAlign = 13) {
    auto ranges = extractGenDetRanges(histName);
    std::string genRange = ranges.first;
    std::string detRange = ranges.second;

    // fallback to token-based matching if regex didn't find explicit ranges
    if (genRange.empty() && detRange.empty()) {
        static const std::vector<std::string> tokens = {
            "5_10", "10_15", "15_20", "20_30", "30_50", "7_50",
            "5_8", "8_11", "11_15", "20_25", "25_30", "30_40", "40_60"
        };
        static const std::vector<std::string> genKeywords = {"gen", "truth", "true", "particle", "mc"};
        static const std::vector<std::string> detKeywords = {"det", "reco", "rec", "meas", "detector"};
        std::string genToken = findPtTokenNearKeywords(histName, genKeywords, tokens);
        std::string detToken = findPtTokenNearKeywords(histName, detKeywords, tokens);
        if (genToken.empty() && detToken.empty()) {
            std::cout << "\t\t[DEBUG] No gen/det pt token found in '" << histName << "'\n";
            for (const auto& t : tokens) {
                if (histName.find(t) != std::string::npos) { genToken = t; break; }
            }
        } else {
            std::cout << "\t\t[DEBUG] Found gen token '" << genToken << "' and det token '" << detToken << "' in '" << histName << "'\n";
        }
        auto formatRange = [](const std::string& token) {
            if (token.empty()) return std::string();
            std::string range = token;
            std::replace(range.begin(), range.end(), '_', '-');
            return range;
        };
        genRange = formatRange(genToken);
        detRange = formatRange(detToken);
    } else {
        std::cout << "\t\t[DEBUG] Extracted gen range '" << genRange << "' and det range '" << detRange << "' from '" << histName << "' using regex\n";
    }

    double dy = 0.065;
    double curY = yNDC;
    if (!genRange.empty()) {
        TLatex* tGen = new TLatex(xNDC, curY, Form("Gen: #it{p}_{T}^{jet}=%s (GeV/#it{c})", genRange.c_str()));
        tGen->SetNDC(); tGen->SetTextFont(42); tGen->SetTextSize(textSize); tGen->SetTextAlign(textAlign);
        tGen->Draw();
        curY -= dy;
    }
    if (!detRange.empty()) {
        TLatex* tDet = new TLatex(xNDC, curY, Form("Det: #it{p}_{T}^{jet}=%s (GeV/#it{c})", detRange.c_str()));
        tDet->SetNDC(); tDet->SetTextFont(42); tDet->SetTextSize(textSize); tDet->SetTextAlign(textAlign);
        tDet->Draw();
    }
}


// compareResponseMatrices("file1.root","file2.root","outdir")
void compareResponseMatrices(const char* file1 = "file1.root", const char* file2 = "file2.root", const char* outDir = "compare_out") {
    gROOT->SetBatch(kTRUE);

    std::string f1s(file1);
    std::string f2s(file2);
    std::string outd(outDir);

    gSystem->Exec(Form("mkdir -p %s", outd.c_str()));

    TFile *f1 = TFile::Open(f1s.c_str(), "READ");
    if (!f1 || f1->IsZombie()) { std::cerr<<"Cannot open "<<f1s<<"\n"; return; }
    TFile *f2 = TFile::Open(f2s.c_str(), "READ");
    if (!f2 || f2->IsZombie()) { std::cerr<<"Cannot open "<<f2s<<"\n"; f1->Close(); return; }

    std::string outRootName = outd + "/response_comparison.root";
    TFile outF(outRootName.c_str(), "RECREATE");

    TIter next1(f1->GetListOfKeys());
    TKey *key;
    int matched = 0;
    int processed = 0;

    while ((key = (TKey*)next1())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl) continue;
        // Only consider TH2 histograms
        if (!cl->InheritsFrom("TH2")) continue;

        std::string name = key->GetName();
        //if name is defined then print it
        if (!name.empty()) {
            std::cout<<"[DEBUG] Found TH2: "<<name<<" in file "<<f1s<<"\n";
        } else {
            std::cout<<"[DEBUG] Found unnamed TH2 in file "<<f1s<<"\n";
        }

        // Retrieve objects from both files
        TH2 *h1 = dynamic_cast<TH2*>(f1->Get(name.c_str()));
        TH2 *h2 = dynamic_cast<TH2*>(f2->Get(name.c_str()));
        if (!h1 || !h2) continue;
        //if the histogram name contains "Rapidity" then skip
        if (name.find("Rapidity") != std::string::npos) continue;
        
        matched++;
        processed++;
        std::cout<<"[DEBUG] Processing #"<<processed<<" name='"<<name<<"' h1="<<h1<<" h2="<<h2<<"\n"<<std::flush;
        // Clone to ensure independent objects
        std::unique_ptr<TH2> h1c((TH2*)h1->Clone(Form("%s_f1", name.c_str())));
        std::unique_ptr<TH2> h2c((TH2*)h2->Clone(Form("%s_f2", name.c_str())));

        // Normalize by number of entries (fall back to integral if entries==0)
        auto normalizeByEntries = [](TH2* h) {
            if (!h) return;
            double entries = h->GetEntries();
            if (entries > 0) {
                h->Scale(1.0 / entries);
                return;
            }
            // fallback to integral of bin contents
            double integral = h->Integral();
            if (integral > 0) h->Scale(1.0 / integral);
        };
        normalizeByEntries(h1c.get());
        normalizeByEntries(h2c.get());
        std::cout<<"[DEBUG] Cloned: h1c="<<h1c.get()<<" h2c="<<h2c.get()<<"\n"<<std::flush;

        // Difference: h1 - h2
        TH2 *hdiff = (TH2*)h1c->Clone(Form("%s_diff", name.c_str()));
        if (!hdiff) { std::cerr<<"[ERROR] Failed to clone hdiff for "<<name<<"\n"; continue; }
        hdiff->Add(h2c.get(), -1.0);

        // Ratio: h1 / h2 (protect division)
        TH2 *hratio = (TH2*)h1c->Clone(Form("%s_ratio", name.c_str()));
        if (!hratio) { std::cerr<<"[ERROR] Failed to clone hratio for "<<name<<"\n"; delete hdiff; continue; }
        int nx = hratio->GetNbinsX();
        int ny = hratio->GetNbinsY();
        std::cout<<"[DEBUG] hratio bins X="<<nx<<" Y="<<ny<<"\n"<<std::flush;
        for (int ix=1; ix<=nx; ++ix) {
            for (int iy=1; iy<=ny; ++iy) {
                double a = h1c->GetBinContent(ix,iy);
                double b = h2c->GetBinContent(ix,iy);
                if (b != 0) hratio->SetBinContent(ix,iy, a / b);
                else hratio->SetBinContent(ix,iy, 0);
            }
        }

        // --- Diagnostics: compare 1D true/reco shapes and per-true-bin reco stats ---
        outF.cd();

        auto normalizeShape = [](TH1* h) {
            if (!h) return;
            double integral = h->Integral();
            if (integral > 0) h->Scale(1.0 / integral);
        };

        // Projections: X = true, Y = reco
        TH1* true_f1 = h1c->ProjectionX(Form("%s_true_f1", name.c_str()));
        TH1* true_f2 = h2c->ProjectionX(Form("%s_true_f2", name.c_str()));
        TH1* reco_f1 = h1c->ProjectionY(Form("%s_reco_f1", name.c_str()));
        TH1* reco_f2 = h2c->ProjectionY(Form("%s_reco_f2", name.c_str()));

        normalizeShape(true_f1);
        normalizeShape(true_f2);
        normalizeShape(reco_f1);
        normalizeShape(reco_f2);

        if (true_f1) true_f1->Write();
        if (true_f2) true_f2->Write();
        if (reco_f1) reco_f1->Write();
        if (reco_f2) reco_f2->Write();

        // Ratio histograms for projections (file1/file2)
        TH1* true_ratio = nullptr;
        TH1* reco_ratio = nullptr;
        if (true_f1 && true_f2) {
            true_ratio = (TH1*)true_f1->Clone(Form("%s_true_ratio", name.c_str()));
            int nb = true_ratio->GetNbinsX();
            for (int ib = 1; ib <= nb; ++ib) {
                double a = true_f1->GetBinContent(ib);
                double b = true_f2->GetBinContent(ib);
                true_ratio->SetBinContent(ib, (b != 0) ? (a / b) : 0.0);
            }
            true_ratio->Write();
        }
        if (reco_f1 && reco_f2) {
            reco_ratio = (TH1*)reco_f1->Clone(Form("%s_reco_ratio", name.c_str()));
            int nb = reco_ratio->GetNbinsX();
            for (int ib = 1; ib <= nb; ++ib) {
                double a = reco_f1->GetBinContent(ib);
                double b = reco_f2->GetBinContent(ib);
                reco_ratio->SetBinContent(ib, (b != 0) ? (a / b) : 0.0);
            }
            reco_ratio->Write();
        }

        // Per-true-bin stats: mean/reco RMS and row sum (integral of P(reco|true) row)
        TH1* tmplX = h1c->ProjectionX(Form("%s_tmplX", name.c_str()));
        TH1* meanReco_f1 = (TH1*)tmplX->Clone(Form("%s_meanReco_f1", name.c_str())); meanReco_f1->Reset();
        TH1* rmsReco_f1  = (TH1*)tmplX->Clone(Form("%s_rmsReco_f1", name.c_str()));  rmsReco_f1->Reset();
        TH1* rowSum_f1   = (TH1*)tmplX->Clone(Form("%s_rowSum_f1", name.c_str()));   rowSum_f1->Reset();

        TH1* meanReco_f2 = (TH1*)tmplX->Clone(Form("%s_meanReco_f2", name.c_str())); meanReco_f2->Reset();
        TH1* rmsReco_f2  = (TH1*)tmplX->Clone(Form("%s_rmsReco_f2", name.c_str()));  rmsReco_f2->Reset();
        TH1* rowSum_f2   = (TH1*)tmplX->Clone(Form("%s_rowSum_f2", name.c_str()));   rowSum_f2->Reset();

        for (int ix = 1; ix <= nx; ++ix) {
            TH1* row1 = h1c->ProjectionY(Form("%s_row1_%d", name.c_str(), ix), ix, ix);
            if (row1) {
                meanReco_f1->SetBinContent(ix, row1->GetMean());
                rmsReco_f1->SetBinContent(ix, row1->GetRMS());
                rowSum_f1->SetBinContent(ix, row1->Integral());
                delete row1;
            }
            TH1* row2 = h2c->ProjectionY(Form("%s_row2_%d", name.c_str(), ix), ix, ix);
            if (row2) {
                meanReco_f2->SetBinContent(ix, row2->GetMean());
                rmsReco_f2->SetBinContent(ix, row2->GetRMS());
                rowSum_f2->SetBinContent(ix, row2->Integral());
                delete row2;
            }
        }

        meanReco_f1->Write(); rmsReco_f1->Write(); rowSum_f1->Write();
        meanReco_f2->Write(); rmsReco_f2->Write(); rowSum_f2->Write();

        // Save diagnostic canvases
        {
            TCanvas* cproj = new TCanvas(Form("cproj_%s", name.c_str()), Form("projections_%s", name.c_str()), 1000, 450);
            cproj->Divide(2, 1);
            cproj->cd(1);
            if (true_f1 && true_f2) {
                true_f1->SetLineColor(kBlue); true_f1->SetLineWidth(2);
                true_f2->SetLineColor(kRed);  true_f2->SetLineWidth(2);
                true_f1->SetTitle(Form("%s: True projection (shape)", name.c_str()));
                true_f1->Draw("hist");
                true_f2->Draw("hist same");
                TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(true_f1, "file1", "l");
                leg->AddEntry(true_f2, "file2", "l");
                leg->Draw();
            }
            drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.040);
            cproj->cd(2);
            if (reco_f1 && reco_f2) {
                reco_f1->SetLineColor(kBlue); reco_f1->SetLineWidth(2);
                reco_f2->SetLineColor(kRed);  reco_f2->SetLineWidth(2);
                reco_f1->SetTitle(Form("%s: Reco projection (shape)", name.c_str()));
                reco_f1->Draw("hist");
                reco_f2->Draw("hist same");
                TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(reco_f1, "file1", "l");
                leg->AddEntry(reco_f2, "file2", "l");
                leg->Draw();
            }
            drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.040);
            std::string outpng = outd + "/diag_projections_" + name + ".png";
            cproj->SaveAs(outpng.c_str());
            delete cproj;
        }
        {
            TCanvas* crat = new TCanvas(Form("crat_%s", name.c_str()), Form("projRatios_%s", name.c_str()), 1000, 450);
            crat->Divide(2, 1);
            crat->cd(1);
            if (true_ratio) {
                true_ratio->SetLineColor(kBlack);
                true_ratio->SetLineWidth(2);
                true_ratio->SetTitle(Form("%s: True proj ratio (f1/f2)", name.c_str()));
                true_ratio->Draw("hist");
            }
            drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.040);
            crat->cd(2);
            if (reco_ratio) {
                reco_ratio->SetLineColor(kBlack);
                reco_ratio->SetLineWidth(2);
                reco_ratio->SetTitle(Form("%s: Reco proj ratio (f1/f2)", name.c_str()));
                reco_ratio->Draw("hist");
            }
            drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.040);
            std::string outpng = outd + "/diag_projectionRatios_" + name + ".png";
            crat->SaveAs(outpng.c_str());
            delete crat;
        }
        {
            TCanvas* cstats = new TCanvas(Form("cstats_%s", name.c_str()), Form("stats_%s", name.c_str()), 1200, 400);
            cstats->Divide(3, 1);
            cstats->cd(1);
            meanReco_f1->SetLineColor(kBlue); meanReco_f1->SetLineWidth(2);
            meanReco_f2->SetLineColor(kRed);  meanReco_f2->SetLineWidth(2);
            meanReco_f1->SetTitle(Form("%s: Mean reco(Y) per true(X)", name.c_str()));
            meanReco_f1->Draw("hist"); meanReco_f2->Draw("hist same");
            {
                TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(meanReco_f1, "file1", "l");
                leg->AddEntry(meanReco_f2, "file2", "l");
                leg->Draw();
            }
            drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.040);
            cstats->cd(2);
            rmsReco_f1->SetLineColor(kBlue); rmsReco_f1->SetLineWidth(2);
            rmsReco_f2->SetLineColor(kRed);  rmsReco_f2->SetLineWidth(2);
            rmsReco_f1->SetTitle(Form("%s: RMS reco(Y) per true(X)", name.c_str()));
            rmsReco_f1->Draw("hist"); rmsReco_f2->Draw("hist same");
            {
                TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(rmsReco_f1, "file1", "l");
                leg->AddEntry(rmsReco_f2, "file2", "l");
                leg->Draw();
            }
            drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.040);
            cstats->cd(3);
            rowSum_f1->SetLineColor(kBlue); rowSum_f1->SetLineWidth(2);
            rowSum_f2->SetLineColor(kRed);  rowSum_f2->SetLineWidth(2);
            rowSum_f1->SetTitle(Form("%s: Row integral (sum_Y) per true(X)", name.c_str()));
            rowSum_f1->Draw("hist"); rowSum_f2->Draw("hist same");
            {
                TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->AddEntry(rowSum_f1, "file1", "l");
                leg->AddEntry(rowSum_f2, "file2", "l");
                leg->Draw();
            }
            drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.040);
            std::string outpng = outd + "/diag_stats_" + name + ".png";
            cstats->SaveAs(outpng.c_str());
            delete cstats;
        }

        delete tmplX;
        delete meanReco_f1;
        delete rmsReco_f1;
        delete rowSum_f1;
        delete meanReco_f2;
        delete rmsReco_f2;
        delete rowSum_f2;

        delete true_ratio;
        delete reco_ratio;
        delete true_f1;
        delete true_f2;
        delete reco_f1;
        delete reco_f2;

        // Write comparison histograms
        h1c->Write();
        h2c->Write();
        hdiff->Write();
        hratio->Write();

        // Create canvas and save images (allocate and delete to free resources)
        TCanvas *c = new TCanvas(Form("c_%s", name.c_str()), name.c_str(), 1200, 400);
        c->Divide(3,1);
        c->cd(1);
        gPad->SetLeftMargin(0.12);
        h1c->SetTitle(Form("%s : File1", name.c_str()));
        h1c->Draw("col");
        drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.06);
        c->cd(2);
        gPad->SetLeftMargin(0.12);
        h2c->SetTitle(Form("%s : File2", name.c_str()));
        h2c->Draw("col");
        drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.06);
        c->cd(3);
        gPad->SetLeftMargin(0.12);
        hdiff->SetTitle(Form("%s : Diff (f1 - f2)", name.c_str()));
        hdiff->Draw("col");
        drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.06);

        std::string png = outd + "/compare_" + name + ".png";
        std::cout<<"[DEBUG] Saving compare png: "<<png<<"\n"<<std::flush;
        c->SaveAs(png.c_str());
        std::cout<<"[DEBUG] Saved compare png for "<<name<<"\n"<<std::flush;

        // also save ratio canvas
        TCanvas *c2 = new TCanvas(Form("c2_%s", name.c_str()), Form("ratio_%s", name.c_str()), 600, 500);
        hratio->SetTitle(Form("%s : Ratio (f1/f2)", name.c_str()));
        hratio->Draw("col");
        drawJetPtLabelsIfFound(name, 0.16, 0.88, 0.040);
        c2->Update();
        std::string png2 = outd + "/ratio_" + name + ".png";
        std::cout<<"[DEBUG] Saving ratio png: "<<png2<<"\n"<<std::flush;
        c2->SaveAs(png2.c_str());
        std::cout<<"[DEBUG] Saved ratio png for "<<name<<"\n"<<std::flush;

        // delete canvases to free resources
        delete c;
        delete c2;
        // Cleanup histograms created by Clone but keep written ones alive in file
        delete hdiff;
        delete hratio;
    }

    std::cout<<"Matched and compared "<<matched<<" histograms. Results in "<<outRootName<<" and "<<outd<<"\n";

    outF.Close();
    f1->Close();
    f2->Close();
}
