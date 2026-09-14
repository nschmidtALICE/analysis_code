// calcTagZCorrections.C
// Load TagZHistograms_<pt>.root files and compute correction factors per tagZ bin.
// Outputs graphs into a single ROOT file named TagZCorrectionFactors_combined.root by default.

#include <TFile.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <regex>

using namespace std;

static double safeDiv(double a, double b)
{
    if (b == 0)
        return 0.0;
    return a / b;
}

static std::vector<std::string> splitCsv(const std::string &csv)
{
    std::vector<std::string> out;
    std::string s;
    for (size_t i = 0; i < csv.size(); ++i)
    {
        char c = csv[i];
        if (c == ',')
        {
            if (!s.empty())
                out.push_back(s);
            s.clear();
        }
        else if (c != ' ' && c != '\t' && c != '\r' && c != '\n')
        {
            s.push_back(c);
        }
    }
    if (!s.empty())
        out.push_back(s);
    return out;
}

void calcTagZCorrections(const std::string &inputDir = ".",
                         const std::string &filePattern = "TagZHistograms_",
                         const std::string &jetBinsCsv = "10_15,15_20,20_30,30_100",
                         const std::string &rapidityLabelsCsv = "2.5-3.0,3.0-3.5,3.5-4.0")
{
    // run ROOT in batch mode (no display)
    gROOT->SetBatch(kTRUE);

    std::filesystem::path dir(inputDir);
    if (!std::filesystem::exists(dir) || !std::filesystem::is_directory(dir))
    {
        std::cerr << "Input directory does not exist: " << inputDir << std::endl;
        return;
    }

    // Prepare output ROOT file
    std::string outName = (dir / "TagZCorrectionFactors_combined.root").string();
    TFile outF(outName.c_str(), "RECREATE");
    if (outF.IsZombie())
    {
        std::cerr << "Cannot create output file: " << outName << std::endl;
        return;
    }

    // parse provided CSVs
    std::vector<std::string> jetBins = splitCsv(jetBinsCsv);
    std::vector<std::string> rapidityLabels = splitCsv(rapidityLabelsCsv);
    int userNBins = static_cast<int>(rapidityLabels.size());
    if (jetBins.empty())
    {
        std::cerr << "No jet bins provided in jetBinsCsv=" << jetBinsCsv << std::endl;
        return;
    }

    // Containers to store histograms for all jetPt x eta bins
    std::vector<std::vector<TH1 *>> all_h_background;
    std::vector<std::vector<TH1 *>> all_h_prompt;
    std::vector<std::vector<TH1 *>> all_h_reco;
    std::vector<std::vector<TH1 *>> all_h_accept;
    std::vector<std::vector<TH1 *>> all_h_full;

    for (const auto &ptTag : jetBins)
    {
        std::string filename = filePattern + ptTag + ".root";
        std::filesystem::path fp = dir / filename;
        std::cout << "Processing: " << fp << std::endl;
        if (!std::filesystem::exists(fp))
        {
            std::cerr << "  File not found: " << fp << std::endl;
            continue;
        }
        TFile f(fp.c_str(), "READ");
        if (f.IsZombie())
        {
            std::cerr << "  Could not open: " << fp << std::endl;
            continue;
        }

        // We expect histograms with names like:
        // tagZHist_<pt>_bin<i>
        // backgroundSubtractedTagZHist_<pt>_bin<i>
        // promptSignalTagZHist_<pt>_bin<i>
        // promptSignalTagZHist_PIDWeighted_<pt>_bin<i>
        // promptSignalTagZHist_RecoWeighted_<pt>_bin<i>
        // promptSignalTagZHist_AcceptanceWeighted_<pt>_bin<i>
        // promptSignalTagZHist_FullyWeighted_<pt>_bin<i>

        // determine number of eta bins from user input (rapidities)
        int nBins = userNBins;
        if (nBins <= 0)
        {
            // try to auto-detect by probing for tagZHist_<pt>_binN
            int detect = 0;
            while (true)
            {
                std::string name = "tagZHist_" + ptTag + "_bin" + std::to_string(detect);
                TH1 *h = dynamic_cast<TH1 *>(f.Get(name.c_str()));
                if (!h)
                    break;
                ++detect;
            }
            nBins = detect;
        }
        if (nBins <= 0)
        {
            std::cerr << "  Could not determine eta bins for " << ptTag << " in file " << filename << std::endl;
            f.Close();
            continue;
        }

        // Determine tagZ binning from bin0 histogram
        TH1 *sample = dynamic_cast<TH1 *>(f.Get(("tagZHist_" + ptTag + "_bin0").c_str()));
        if (!sample)
            sample = dynamic_cast<TH1 *>(f.Get(("promptSignalTagZHist_" + ptTag + "_bin0").c_str()));
        if (!sample)
        {
            std::cerr << "  Could not find sample tagZ histogram for " << ptTag << " in file " << filename << std::endl;
            f.Close();
            continue;
        }
        int nbTagZ = sample->GetNbinsX();

        // Prepare graphs for each correction type. One graph per ptTag storing values vs tagZ (per z-bin bin content)
        // We'll create one TGraphErrors per correction type per z-bin index (i.e., for each tagZ histogram bin inside the tagZ histograms)

        // Create histograms arrays: per eta-bin create TH1 copies of 'sample' for ratios
        std::vector<TH1 *> h_background(nBins, nullptr);
        std::vector<TH1 *> h_prompt(nBins, nullptr);
        std::vector<TH1 *> h_reco(nBins, nullptr);
        std::vector<TH1 *> h_accept(nBins, nullptr);
        std::vector<TH1 *> h_full(nBins, nullptr);

        for (int i = 0; i < nBins; ++i)
        {
            // name histograms with jetpt and rapidity label if available
            std::string label = ptTag + "_eta" + std::to_string(i);
            if (i < (int)rapidityLabels.size())
                label = ptTag + "_" + rapidityLabels[i];
            TH1 *hb = (TH1 *)sample->Clone(("h_background_" + label).c_str());
            hb->Reset();
            hb->SetDirectory(nullptr);
            TH1 *hp = (TH1 *)sample->Clone(("h_prompt_" + label).c_str());
            hp->Reset();
            hp->SetDirectory(nullptr);
            TH1 *hr = (TH1 *)sample->Clone(("h_reco_" + label).c_str());
            hr->Reset();
            hr->SetDirectory(nullptr);
            TH1 *ha = (TH1 *)sample->Clone(("h_accept_" + label).c_str());
            ha->Reset();
            ha->SetDirectory(nullptr);
            TH1 *hf = (TH1 *)sample->Clone(("h_full_" + label).c_str());
            hf->Reset();
            hf->SetDirectory(nullptr);
            h_background[i] = hb;
            h_prompt[i] = hp;
            h_reco[i] = hr;
            h_accept[i] = ha;
            h_full[i] = hf;
        }

        // Loop over z bins and fill graphs from histogram bin-by-bin
        for (int i = 0; i < nBins; ++i)
        {
            std::string name_orig = "tagZHist_" + ptTag + "_bin" + std::to_string(i);
            std::string name_bkg = "backgroundSubtractedTagZHist_" + ptTag + "_bin" + std::to_string(i);
            std::string name_prompt = "promptSignalTagZHist_" + ptTag + "_bin" + std::to_string(i);
            std::string name_pid = "promptSignalTagZHist_PIDWeighted_" + ptTag + "_bin" + std::to_string(i);
            std::string name_reco = "promptSignalTagZHist_RecoWeighted_" + ptTag + "_bin" + std::to_string(i);
            std::string name_acc = "promptSignalTagZHist_AcceptanceWeighted_" + ptTag + "_bin" + std::to_string(i);
            std::string name_full = "promptSignalTagZHist_FullyWeighted_" + ptTag + "_bin" + std::to_string(i);

            TH1 *h_orig_hist = dynamic_cast<TH1 *>(f.Get(name_orig.c_str()));
            TH1 *h_bkg_hist = dynamic_cast<TH1 *>(f.Get(name_bkg.c_str()));
            TH1 *h_prompt_hist = dynamic_cast<TH1 *>(f.Get(name_prompt.c_str()));
            TH1 *h_pid_hist = dynamic_cast<TH1 *>(f.Get(name_pid.c_str()));
            TH1 *h_reco_hist = dynamic_cast<TH1 *>(f.Get(name_reco.c_str()));
            TH1 *h_acc_hist = dynamic_cast<TH1 *>(f.Get(name_acc.c_str()));
            TH1 *h_full_hist = dynamic_cast<TH1 *>(f.Get(name_full.c_str()));

            if (!h_orig_hist && !h_prompt_hist)
            {
                std::cerr << "  Missing essential histograms for bin " << i << " in " << filename << std::endl;
                continue;
            }

            TH1 *h_ref = h_orig_hist ? h_orig_hist : h_prompt_hist;
            int nb = h_ref->GetNbinsX();
            for (int b = 1; b <= nb; ++b)
            {
                double val_orig = h_orig_hist ? h_orig_hist->GetBinContent(b) : 0.0;
                double err_orig = h_orig_hist ? h_orig_hist->GetBinError(b) : 0.0;
                double val_bkg = h_bkg_hist ? h_bkg_hist->GetBinContent(b) : 0.0;
                double err_bkg = h_bkg_hist ? h_bkg_hist->GetBinError(b) : 0.0;
                double val_prompt = h_prompt_hist ? h_prompt_hist->GetBinContent(b) : 0.0;
                double err_prompt = h_prompt_hist ? h_prompt_hist->GetBinError(b) : 0.0;
                double val_reco = h_reco_hist ? h_reco_hist->GetBinContent(b) : 0.0;
                double err_reco = h_reco_hist ? h_reco_hist->GetBinError(b) : 0.0;
                double val_acc = h_acc_hist ? h_acc_hist->GetBinContent(b) : 0.0;
                double err_acc = h_acc_hist ? h_acc_hist->GetBinError(b) : 0.0;
                double val_full = h_full_hist ? h_full_hist->GetBinContent(b) : 0.0;
                double err_full = h_full_hist ? h_full_hist->GetBinError(b) : 0.0;

                // background correction = backgroundSubtracted / original
                double bg = safeDiv(val_bkg, val_orig);
                double bg_err = 0.0;
                if (val_bkg > 0 && val_orig > 0)
                    bg_err = bg * sqrt((err_bkg / val_bkg) * (err_bkg / val_bkg) + (err_orig / val_orig) * (err_orig / val_orig));

                // prompt correction = prompt / backgroundSubtracted
                double pr = safeDiv(val_prompt, val_bkg);
                double pr_err = 0.0;
                if (val_prompt > 0 && val_bkg > 0)
                    pr_err = pr * sqrt((err_prompt / val_prompt) * (err_prompt / val_prompt) + (err_bkg / val_bkg) * (err_bkg / val_bkg));

                // reconstruction efficiency estimate = original / recoWeighted (since recoWeighted = original / recoEff)
                double reco_eff = safeDiv(val_orig, val_reco);
                double reco_err = 0.0;
                if (val_orig > 0 && val_reco > 0)
                    reco_err = reco_eff * sqrt((err_orig / val_orig) * (err_orig / val_orig) + (err_reco / val_reco) * (err_reco / val_reco));

                // acceptance estimate = prompt / acceptanceWeighted
                double acc_eff = safeDiv(val_prompt, val_acc);
                double acc_err = 0.0;
                if (val_prompt > 0 && val_acc > 0)
                    acc_err = acc_eff * sqrt((err_prompt / val_prompt) * (err_prompt / val_prompt) + (err_acc / val_acc) * (err_acc / val_acc));

                // full correction (product) = original / fullyWeighted = combined_all * reco * acceptance
                double full_eff = safeDiv(val_orig, val_full);
                double full_err = 0.0;
                if (val_orig > 0 && val_full > 0)
                    full_err = full_eff * sqrt((err_orig / val_orig) * (err_orig / val_orig) + (err_full / val_full) * (err_full / val_full));

                h_background[i]->SetBinContent(b, bg);
                h_background[i]->SetBinError(b, bg_err);

                h_prompt[i]->SetBinContent(b, pr);
                h_prompt[i]->SetBinError(b, pr_err);

                h_reco[i]->SetBinContent(b, reco_eff);
                h_reco[i]->SetBinError(b, reco_err);

                h_accept[i]->SetBinContent(b, acc_eff);
                h_accept[i]->SetBinError(b, acc_err);

                h_full[i]->SetBinContent(b, full_eff);
                h_full[i]->SetBinError(b, full_err);
            }
        }

        // Write graphs into a folder within the output file
        outF.mkdir(ptTag.c_str());
        outF.cd(ptTag.c_str());
        for (int i = 0; i < nBins; ++i)
        {
            if (h_background[i])
                h_background[i]->Write();
            if (h_prompt[i])
                h_prompt[i]->Write();
            if (h_reco[i])
                h_reco[i]->Write();
            if (h_accept[i])
                h_accept[i]->Write();
            if (h_full[i])
                h_full[i]->Write();
        }

        // store pointers so we can later build a combined overview across all jetPt x eta
        all_h_background.emplace_back();
        all_h_background.back() = h_background;
        all_h_prompt.emplace_back();
        all_h_prompt.back() = h_prompt;
        all_h_reco.emplace_back();
        all_h_reco.back() = h_reco;
        all_h_accept.emplace_back();
        all_h_accept.back() = h_accept;
        all_h_full.emplace_back();
        all_h_full.back() = h_full;

        // Create an overview canvas with subpads for each jet-pt / eta bin
        int cols = std::min(4, nBins);
        int rows = (nBins + cols - 1) / cols;
        TCanvas *c = new TCanvas(("c_overview_" + ptTag).c_str(), ("Overview " + ptTag).c_str(), 1200, 800);
        c->Divide(cols, rows);

        for (int i = 0; i < nBins; ++i)
        {
            c->cd(i + 1);
            c->cd(i + 1)->SetLeftMargin(0.12);
            c->cd(i + 1)->SetBottomMargin(0.10);
            c->cd(i + 1)->SetRightMargin(0.01);
            c->cd(i + 1)->SetTopMargin(0.01);
            // determine y-range across all hists for this bin
            double yMin = 1e99, yMax = -1e99;
            auto scanRangeH = [&](TH1 *h)
            {
                if (!h)
                    return;
                int nb = h->GetNbinsX();
                for (int bi = 1; bi <= nb; ++bi)
                {
                    double val = h->GetBinContent(bi);
                    if (val == val)
                    {
                        if (val < yMin)
                            yMin = val;
                        if (val > yMax)
                            yMax = val;
                    }
                }
            };
            scanRangeH(h_background[i]);
            scanRangeH(h_prompt[i]);
            scanRangeH(h_reco[i]);
            scanRangeH(h_accept[i]);
            scanRangeH(h_full[i]);
            if (yMin > yMax)
            {
                yMin = 0;
                yMax = 1;
            }
            double ypad = (yMax - yMin) * 0.1;
            yMin -= ypad;
            yMax += ypad;

            bool firstDraw = true;
            auto drawHist = [&](TH1 *h, int color, const char *opt)
            {
                if (!h)
                    return;
                h->SetLineColor(color);
                h->SetMarkerColor(color);
                h->SetMarkerStyle(20);
                if (firstDraw)
                {
                    h->GetYaxis()->SetRangeUser(yMin, yMax);
                    h->GetXaxis()->SetTitle(sample->GetXaxis()->GetTitle());
                    h->Draw("E1");
                    firstDraw = false;
                }
                else
                {
                    h->Draw("E1same");
                }
            };

            drawHist(h_background[i], kRed, "");
            drawHist(h_prompt[i], kBlue, "");
            drawHist(h_reco[i], kGreen + 2, "");
            drawHist(h_accept[i], kMagenta, "");
            drawHist(h_full[i], kBlack, "");

            // per-pad legend
            TLegend *leg = new TLegend(0.32, 0.72, 0.92, 0.92);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.04);
            if (h_background[i])
                leg->AddEntry(h_background[i], "background / original", "lep");
            if (h_prompt[i])
                leg->AddEntry(h_prompt[i], "prompt / background", "lep");
            if (h_reco[i])
                leg->AddEntry(h_reco[i], "reco efficiency", "lep");
            if (h_accept[i])
                leg->AddEntry(h_accept[i], "acceptance", "lep");
            if (h_full[i])
                leg->AddEntry(h_full[i], "full correction", "lep");
            leg->Draw();
            // delete leg;
        }

        // write canvas and also save as PNG
        outF.cd(ptTag.c_str());
        c->Write();
        std::string pngName = (dir / ("TagZOverview_" + ptTag + ".png")).string();
        c->SaveAs(pngName.c_str());
        outF.cd();

        // delete c;

        f.Close();
        std::cout << "  Done " << filename << " -> stored graphs under " << ptTag << std::endl;
    }

    // Build combined canvas: rows = number of jetPt bins, cols = number of eta bins (use userNBins)
    int nJet = static_cast<int>(all_h_background.size());
    int nEta = userNBins;
    if (nJet > 0 && nEta > 0)
    {
        TCanvas *cAll = new TCanvas("c_all_overview", "All jetPt x eta overview", 1200, 1200);
        int cols = nEta;
        int rows = nJet;
        cAll->Divide(cols, rows);
        for (int ir = 0; ir < rows; ++ir)
        {
            for (int ic = 0; ic < cols; ++ic)
            {
                int pad = ir * cols + ic + 1;
                cAll->cd(pad);
                // cAll->cd(pad)->SetLogy(1);
                cAll->cd(pad)->SetTickx(1);
                cAll->cd(pad)->SetTicky(1);
                cAll->cd(pad)->SetRightMargin(0.01);
                cAll->cd(pad)->SetTopMargin(0.01);
                cAll->cd(pad)->SetLeftMargin(0.08);
                cAll->cd(pad)->SetBottomMargin(0.08);
                TH1 *hb = nullptr, *hp = nullptr, *hr = nullptr, *ha = nullptr, *hf = nullptr;
                if (ir < (int)all_h_background.size() && ic < (int)all_h_background[ir].size())
                    hb = all_h_background[ir][ic];
                if (ir < (int)all_h_prompt.size() && ic < (int)all_h_prompt[ir].size())
                    hp = all_h_prompt[ir][ic];
                if (ir < (int)all_h_reco.size() && ic < (int)all_h_reco[ir].size())
                    hr = all_h_reco[ir][ic];
                if (ir < (int)all_h_accept.size() && ic < (int)all_h_accept[ir].size())
                    ha = all_h_accept[ir][ic];
                if (ir < (int)all_h_full.size() && ic < (int)all_h_full[ir].size())
                    hf = all_h_full[ir][ic];

                // determine y range
                double yMin = 1e99, yMax = -1e99;
                auto scanH = [&](TH1 *h)
                { if(!h) return; for(int b=1;b<=h->GetNbinsX();++b){ double v=h->GetBinContent(b); if(v==v){ if(v<yMin) yMin=v; if(v>yMax) yMax=v; }} };
                scanH(hb);
                scanH(hp);
                scanH(hr);
                scanH(ha);
                scanH(hf);
                if (yMin > yMax)
                {
                    yMin = 0;
                    yMax = 1;
                }
                double ypad = (yMax - yMin) * 0.1;
                yMin -= ypad;
                yMax += ypad;

                bool first = true;
                auto drawH = [&](TH1 *h, int col)
                {
                    if(!h)
                        return;
                    h->SetLineColor(col);
                    h->SetMarkerColor(col);
                    if(first)
                    {
                        h->GetYaxis()->SetRangeUser(0,yMax);
                        h->SetTitle("");
                        h->Draw("E1");
                        first=false;
                    }
                    else
                    {
                        h->Draw("E1same");
                    }
                };
                drawH(hb, kRed);
                drawH(hp, kBlue);
                drawH(hr, kGreen + 2);
                drawH(ha, kMagenta);
                drawH(hf, kBlack);

                if(pad==1)
                {
                    TLegend *leg = new TLegend(0.32, 0.46, 0.92, 0.76);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextSize(0.04);
                    if (hb)
                        leg->AddEntry(hb, "background / original", "lep");
                    if (hp)
                        leg->AddEntry(hp, "prompt / background", "lep");
                    if (hr)
                        leg->AddEntry(hr, "reco efficiency", "lep");
                    if (ha)
                        leg->AddEntry(ha, "acceptance", "lep");
                    if (hf)
                        leg->AddEntry(hf, "full correction", "lep");
                    leg->Draw();
                }

                TLine* unityLine = new TLine(hb->GetXaxis()->GetXmin(), 1.0, hb->GetXaxis()->GetXmax(), 1.0);
                unityLine->SetLineStyle(2);
                unityLine->SetLineColor(kGray+2);
                unityLine->Draw("same");

                // small label: jetpt-eta
                std::string lab = (ir < (int)jetBins.size() ? jetBins[ir] : std::string(""));
                if (ic < (int)rapidityLabels.size())
                    lab += " " + rapidityLabels[ic];
                TLatex t;
                t.SetNDC();
                t.SetTextSize(0.03);
                t.DrawLatex(0.12, 0.92, lab.c_str());
            }
        }
        outF.cd();
        cAll->Write();
        std::string pngAll = (dir / "TagZOverview_all.png").string();
        cAll->SaveAs(pngAll.c_str());
        delete cAll;
    }

    outF.Close();
    std::cout << "All done. Corrections saved to " << outName << std::endl;
}

// End of file
