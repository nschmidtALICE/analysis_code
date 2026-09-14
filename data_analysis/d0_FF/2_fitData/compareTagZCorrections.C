// compareTagZCorrections.C
// Compare TagZ correction factors across systematic-variation input directories.
// Usage: compareTagZCorrections(baseDir, varDirsCsv, varLabelsCsv, jetBinsCsv, rapidityLabelsCsv)

#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TGraphErrors.h>

#include <iostream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <cmath>

static std::vector<std::string> splitCsv(const std::string &csv) {
    std::vector<std::string> out; std::string s;
    for (char c: csv) {
        if (c==',') { if(!s.empty()) out.push_back(s); s.clear(); }
        else if (c!=' ' && c!='\t' && c!='\r' && c!='\n') s.push_back(c);
    }
    if(!s.empty()) out.push_back(s);
    return out;
}

static bool file_exists_global(const std::string &p) {
    struct stat sb; return (stat(p.c_str(), &sb) == 0);
}

static double safeDiv(double a, double b){ if (b==0) return 0.; return a/b; }

// compute correction histograms for a given input directory and store per-jet,per-eta histograms
static bool loadCorrectionHists(const std::string &dir, const std::vector<std::string> &jetBins,
                                int nEta, std::vector<std::vector<TH1D*>> &out_bkg,
                                std::vector<std::vector<TH1D*>> &out_prompt,
                                std::vector<std::vector<TH1D*>> &out_reco,
                                std::vector<std::vector<TH1D*>> &out_accept,
                                std::vector<std::vector<TH1D*>> &out_full)
{
    auto join_path = [](const std::string &a, const std::string &b)->std::string {
        if (a.empty()) return b;
        if (a.back()=='/') return a + b;
        return a + '/' + b;
    };
    auto file_exists = [](const std::string &p)->bool {
        struct stat sb; return (stat(p.c_str(), &sb) == 0);
    };
    out_bkg.clear(); out_prompt.clear(); out_reco.clear(); out_accept.clear(); out_full.clear();
    for (size_t j=0;j<jetBins.size();++j) {
        const std::string &ptTag = jetBins[j];
        std::string fname = "TagZHistograms_" + ptTag + ".root";
        std::string fp = join_path(dir, fname);
        std::cout<<"[load] checking file: "<<fp<<"\n";
        if (!file_exists(fp)) { std::cerr<<"[load] file missing: "<<fp<<"\n"; return false; }
        TFile f(fp.c_str(), "READ"); if (f.IsZombie()) { std::cerr<<"[load] cannot open "<<fp<<"\n"; return false; }
        std::cout<<"[load] opened file: "<<fp<<"\n";

        // sample hist to clone binning
        TH1D *sample = dynamic_cast<TH1D*>(f.Get(("tagZHist_"+ptTag+"_bin0").c_str()));
        if (!sample) sample = dynamic_cast<TH1D*>(f.Get(("promptSignalTagZHist_"+ptTag+"_bin0").c_str()));
        if (!sample) { std::cerr<<"[load] no sample tagZ in "<<fp<<"\n"; f.Close(); return false; }
        else { std::cout<<"[load] sample hist found for ptTag="<<ptTag<<"\n"; }

        // allocate per-eta vectors
        std::vector<TH1D*> vb(nEta,nullptr), vp(nEta,nullptr), vr(nEta,nullptr), va(nEta,nullptr), vf(nEta,nullptr);

        for (int i=0;i<nEta;++i) {
            std::string label = ptTag + "_eta" + std::to_string(i);
            if (i < 1000) { /* label kept simple */ }
            TH1D *hb = (TH1D*)sample->Clone(("bkg_"+label).c_str()); hb->Reset(); hb->SetDirectory(nullptr);
            TH1D *hp = (TH1D*)sample->Clone(("pr_"+label).c_str()); hp->Reset(); hp->SetDirectory(nullptr);
            TH1D *hr = (TH1D*)sample->Clone(("reco_"+label).c_str()); hr->Reset(); hr->SetDirectory(nullptr);
            TH1D *ha = (TH1D*)sample->Clone(("acc_"+label).c_str()); ha->Reset(); ha->SetDirectory(nullptr);
            TH1D *hf = (TH1D*)sample->Clone(("full_"+label).c_str()); hf->Reset(); hf->SetDirectory(nullptr);
            vb[i]=hb; vp[i]=hp; vr[i]=hr; va[i]=ha; vf[i]=hf;

            // load histograms from file and fill ratios same as earlier macro
            std::string name_orig = "tagZHist_"+ptTag+"_bin"+std::to_string(i);
            std::string name_bkg = "backgroundSubtractedTagZHist_"+ptTag+"_bin"+std::to_string(i);
            std::string name_prompt = "promptSignalTagZHist_"+ptTag+"_bin"+std::to_string(i);
            std::string name_reco = "promptSignalTagZHist_RecoWeighted_"+ptTag+"_bin"+std::to_string(i);
            std::string name_acc = "promptSignalTagZHist_AcceptanceWeighted_"+ptTag+"_bin"+std::to_string(i);
            std::string name_full = "promptSignalTagZHist_FullyWeighted_"+ptTag+"_bin"+std::to_string(i);

            TH1D *h_orig = dynamic_cast<TH1D*>(f.Get(name_orig.c_str()));
            TH1D *h_bkg = dynamic_cast<TH1D*>(f.Get(name_bkg.c_str()));
            TH1D *h_prompt = dynamic_cast<TH1D*>(f.Get(name_prompt.c_str()));
            TH1D *h_reco = dynamic_cast<TH1D*>(f.Get(name_reco.c_str()));
            TH1D *h_acc = dynamic_cast<TH1D*>(f.Get(name_acc.c_str()));
            TH1D *h_full = dynamic_cast<TH1D*>(f.Get(name_full.c_str()));

            TH1D *h_ref = h_orig ? h_orig : h_prompt;
            if (!h_ref) { std::cout<<"[load] no reference hist for "<<ptTag<<" bin"<<i<<" (names tried: "<<name_orig<<","<<name_prompt<<")\n"; /* leave zeros */ continue; }
            // debug: print basic info about reference and related hists
            int nb = h_ref->GetNbinsX();
            double sumRef = 0; for (int bb=1; bb<=nb; ++bb) sumRef += h_ref->GetBinContent(bb);
            std::cout<<"[load] ref hist nbins="<<nb<<" integral="<<sumRef<<" for "<<ptTag<<" bin"<<i<<"\n";
            if (h_orig) { double s=0; for (int bb=1; bb<=nb; ++bb) s+=h_orig->GetBinContent(bb); std::cout<<"[load] h_orig integral="<<s<<"\n"; }
            if (h_bkg) { double s=0; for (int bb=1; bb<=nb; ++bb) s+=h_bkg->GetBinContent(bb); std::cout<<"[load] h_bkg integral="<<s<<"\n"; }
            if (h_prompt) { double s=0; for (int bb=1; bb<=nb; ++bb) s+=h_prompt->GetBinContent(bb); std::cout<<"[load] h_prompt integral="<<s<<"\n"; }
            for (int b=1;b<=nb;++b) {
                double val_orig = h_orig ? h_orig->GetBinContent(b) : 0.0;
                double val_b = h_bkg ? h_bkg->GetBinContent(b) : 0.0;
                double val_pr = h_prompt ? h_prompt->GetBinContent(b) : 0.0;
                double val_r = h_reco ? h_reco->GetBinContent(b) : 0.0;
                double val_a = h_acc ? h_acc->GetBinContent(b) : 0.0;
                double val_f = h_full ? h_full->GetBinContent(b) : 0.0;
                // compute correction definitions used earlier
                double bg = safeDiv(val_b, val_orig);
                double prc = safeDiv(val_pr, val_b);
                double reco_eff = safeDiv(val_orig, val_r);
                double acc_eff = safeDiv(val_pr, val_a);
                double full_eff = safeDiv(val_orig, val_f);
                vb[i]->SetBinContent(b, bg);
                vp[i]->SetBinContent(b, prc);
                vr[i]->SetBinContent(b, reco_eff);
                va[i]->SetBinContent(b, acc_eff);
                vf[i]->SetBinContent(b, full_eff);
            }
        }

        out_bkg.push_back(vb); out_prompt.push_back(vp); out_reco.push_back(vr); out_accept.push_back(va); out_full.push_back(vf);
        std::cout<<"[load] filled corrections for ptTag="<<ptTag<<" (nEta="<<nEta<<")\n";
        f.Close();
    }
    return true;
}

void compareTagZCorrections(const std::string &baseDir,
                            const std::string &varDirsCsv,
                            const std::string &varLabelsCsv = "",
                            const std::string &jetBinsCsv = "10_15,15_20,20_30,30_100",
                            const std::string &rapidityLabelsCsv = "2.5-3.0,3.0-3.5,3.5-4.0")
{
    gROOT->SetBatch(kTRUE);
    std::vector<std::string> jetBins = splitCsv(jetBinsCsv);
    std::vector<std::string> rapidityLabels = splitCsv(rapidityLabelsCsv);
    int nEta = (int)rapidityLabels.size();

    // parse variations
    std::vector<std::string> varDirs = splitCsv(varDirsCsv);
    std::vector<std::string> varLabels = splitCsv(varLabelsCsv);
    if (varDirs.empty()) { std::cerr<<"No variation dirs provided\n"; return; }
    // include base as first entry
    std::vector<std::string> allDirs; allDirs.push_back(baseDir);
    allDirs.insert(allDirs.end(), varDirs.begin(), varDirs.end());
    std::vector<std::string> allLabels; allLabels.push_back("default");
    if (varLabels.empty()) {
        for (size_t i=0;i<varDirs.size();++i) allLabels.push_back("var"+std::to_string(i+1));
    } else {
        allLabels.insert(allLabels.end(), varLabels.begin(), varLabels.end());
    }

    int nVar = (int)allDirs.size();

    // For each variation, load correction histograms
    std::vector<std::vector<std::vector<TH1D*>>> H_bkg(nVar), H_pr(nVar), H_reco(nVar), H_acc(nVar), H_full(nVar);
    for (int v=0; v<nVar; ++v) {
        std::vector<std::vector<TH1D*>> vb, vp, vr, va, vf;
        bool ok = loadCorrectionHists(allDirs[v], jetBins, nEta, vb, vp, vr, va, vf);
        if (!ok) { std::cerr<<"Failed loading dir: "<<allDirs[v]<<"\n"; return; }
        std::cout<<"[compare] loaded dir("<<v<<")="<<allDirs[v]<<" -> jetBins="<<vb.size()<<"\n";
        H_bkg[v] = std::move(vb); H_pr[v] = std::move(vp); H_reco[v] = std::move(vr); H_acc[v] = std::move(va); H_full[v] = std::move(vf);
    }

    // output file for comparison canvases
    auto join_path = [](const std::string &a, const std::string &b)->std::string {
        if (a.empty()) return b;
        if (a.back()=='/') return a + b;
        return a + '/' + b;
    };
    auto make_dirs = [](const std::string &p)->bool {
        std::string cmd = std::string("mkdir -p ") + p;
        int r = system(cmd.c_str());
        return (r == 0);
    };
    std::string outp = join_path(baseDir, "TagZComparison_outputs");
    if (!make_dirs(outp)) { std::cerr<<"[compare] failed to create output dir: "<<outp<<"\n"; }
    else { std::cout<<"[compare] output dir: "<<outp<<"\n"; }

    // --- Load PID TagZ correction factors (Kaon/Pion/Combined) from TagZCorrectionFactors.root in each directory
    std::vector<std::vector<std::vector<TGraphErrors*>>> gKaon_all(nVar), gPion_all(nVar), gComb_all(nVar);
    for (int v=0; v<nVar; ++v) {
        gKaon_all[v].assign(jetBins.size(), std::vector<TGraphErrors*>(nEta, nullptr));
        gPion_all[v].assign(jetBins.size(), std::vector<TGraphErrors*>(nEta, nullptr));
        gComb_all[v].assign(jetBins.size(), std::vector<TGraphErrors*>(nEta, nullptr));
        std::string pidFile = join_path(allDirs[v], "TagZCorrectionFactors.root");
        if (!file_exists_global(pidFile)) { std::cout<<"[compare] PID file not found: "<<pidFile<<"\n"; continue; }
        TFile pf(pidFile.c_str(), "READ");
        if (pf.IsZombie()) { std::cerr<<"[compare] PID file is zombie: "<<pidFile<<"\n"; pf.Close(); continue; }
        std::cout<<"[compare] opened PID file: "<<pidFile<<" (dirIndex="<<v<<")\n";
        for (size_t j=0;j<jetBins.size();++j) {
            const std::string &ptTag = jetBins[j];
            for (int i=0;i<nEta;++i) {
                std::string nameK = "tagZKaonCorrection_" + ptTag + "_bin" + std::to_string(i);
                std::string nameP = "tagZPionCorrection_" + ptTag + "_bin" + std::to_string(i);
                std::string nameC = "tagZCombinedCorrection_" + ptTag + "_bin" + std::to_string(i);
                TGraphErrors *gk = dynamic_cast<TGraphErrors*>(pf.Get(nameK.c_str()));
                TGraphErrors *gp = dynamic_cast<TGraphErrors*>(pf.Get(nameP.c_str()));
                TGraphErrors *gc = dynamic_cast<TGraphErrors*>(pf.Get(nameC.c_str()));
                if (gk) { gKaon_all[v][j][i] = (TGraphErrors*)gk->Clone(); std::cout<<"[compare][pid] loaded "<<nameK<<" points="<<gk->GetN()<<"\n"; }
                if (gp) { gPion_all[v][j][i] = (TGraphErrors*)gp->Clone(); std::cout<<"[compare][pid] loaded "<<nameP<<" points="<<gp->GetN()<<"\n"; }
                if (gc) { gComb_all[v][j][i] = (TGraphErrors*)gc->Clone(); std::cout<<"[compare][pid] loaded "<<nameC<<" points="<<gc->GetN()<<"\n"; }
            }
        }
        pf.Close();
    }

    // list of factors
    struct Factor { std::string name; decltype(H_bkg)& container; };
    // cannot use decltype easily; do individually
    std::vector<std::string> factors = {"background","prompt","reco","accept","full"};

    // colors for variations
    std::vector<int> colors = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+2};

    // For each factor create a multi-panel canvas rows=jetBins, cols=nEta with overlayed ratio curves
    for (const auto &factor : factors) {
        TCanvas *cAll = new TCanvas(("c_comp_"+factor).c_str(), ("Comparison "+factor).c_str(), 1600, 1200);
        int rows = (int)jetBins.size(); int cols = nEta; if (rows<1||cols<1) { continue; }
        cAll->Divide(cols, rows);
        for (int ir=0; ir<rows; ++ir) {
            for (int ic=0; ic<cols; ++ic) {
                int pad = ir*cols + ic + 1; cAll->cd(pad);
                // determine default and variations
                TH1D *hdef=nullptr;
                if (!H_bkg[0].empty() && ir < (int)H_bkg[0].size() && ic < (int)H_bkg[0][ir].size()) {
                    if (factor=="background") hdef = H_bkg[0][ir][ic];
                    if (factor=="prompt") hdef = H_pr[0][ir][ic];
                    if (factor=="reco") hdef = H_reco[0][ir][ic];
                    if (factor=="accept") hdef = H_acc[0][ir][ic];
                    if (factor=="full") hdef = H_full[0][ir][ic];
                }
                if (!hdef) { std::cout<<"[compare] no default hist for factor="<<factor<<" pt="<<jetBins[ir]<<" etaIndex="<<ic<<" (skipping)\n"; continue; }
                double ymin=1e99,ymax=-1e99;
                // build ratio histograms for each variation (skip index 0 which is default)
                std::vector<TH1D*> ratios; ratios.reserve(nVar-1);
                for (int v=1; v<nVar; ++v) {
                    TH1D *hvar = nullptr;
                    if (factor=="background") if (ir < (int)H_bkg[v].size()) hvar = H_bkg[v][ir][ic];
                    if (factor=="prompt") if (ir < (int)H_pr[v].size()) hvar = H_pr[v][ir][ic];
                    if (factor=="reco") if (ir < (int)H_reco[v].size()) hvar = H_reco[v][ir][ic];
                    if (factor=="accept") if (ir < (int)H_acc[v].size()) hvar = H_acc[v][ir][ic];
                    if (factor=="full") if (ir < (int)H_full[v].size()) hvar = H_full[v][ir][ic];
                    if (!hvar) { ratios.push_back(nullptr); continue; }
                    TH1D *r = (TH1D*)hvar->Clone(("ratio_"+std::to_string(v)+"_"+std::to_string(ir)+"_"+std::to_string(ic)).c_str());
                    r->SetDirectory(nullptr);
                    double int_hvar=0, int_hdef=0;
                    int nbins = hvar->GetNbinsX(); for (int bi=1; bi<=nbins; ++bi) int_hvar += hvar->GetBinContent(bi);
                    int nbdef = hdef->GetNbinsX(); for (int bi=1; bi<=nbdef; ++bi) int_hdef += hdef->GetBinContent(bi);
                    // compute ratio per-bin and assign binomial errors using N = denominator counts
                    int nbins_r = std::min(r->GetNbinsX(), hdef->GetNbinsX());
                    for (int b=1; b<=nbins_r; ++b) {
                        double vnum = hvar->GetBinContent(b);
                        double vden = hdef->GetBinContent(b);
                        double val = safeDiv(vnum, vden);
                        double err = 0.0;
                        if (vden > 0) {
                            // binomial error for efficiency-like ratio: sqrt(p*(1-p)/N)
                            double p = val;
                            err = std::sqrt(std::max(0.0, p*(1.0-p)/vden));
                        }
                        r->SetBinContent(b, val);
                        r->SetBinError(b, err);
                    }
                    // diagnostics
                    double int_r=0; int nanCount=0; int infCount=0; int nonzero=0;
                    for (int b=1;b<=r->GetNbinsX();++b) { double val=r->GetBinContent(b); if (std::isnan(val)) ++nanCount; if (std::isinf(val)) ++infCount; if (val!=0.0) ++nonzero; int_r += (val==val && std::isfinite(val)) ? val : 0.0; }
                    std::cout<<"[compare][diag] varIndex="<<v<<" factor="<<factor<<" pt="<<jetBins[ir]<<" eta="<<ic<<" nbins="<<r->GetNbinsX()
                             <<" int_hvar="<<int_hvar<<" int_hdef="<<int_hdef<<" int_r="<<int_r
                             <<" nan="<<nanCount<<" inf="<<infCount<<" nonzero="<<nonzero<<"\n";
                    // print first few bins
                    int show = std::min(5, r->GetNbinsX()); std::cout<<"[compare][diag] first_bins:"; for (int b=1;b<=show;++b) std::cout<<" "<<r->GetBinContent(b); std::cout<<"\n";
                    ratios.push_back(r);
                    // update range
                    for (int b=1;b<=r->GetNbinsX();++b) { double val=r->GetBinContent(b); if (val==val){ if (val<ymin) ymin=val; if (val>ymax) ymax=val; } }
                }
                if (ymin>ymax) { ymin=0.8; ymax=1.2; }
                // debug: how many valid ratio histograms
                int validRatios=0; for (TH1D* rr: ratios) if (rr) ++validRatios;
                std::cout<<"[compare] pad(pt="<<jetBins[ir]<<",etaIdx="<<ic<<") validRatios="<<validRatios<<" ymin="<<ymin<<" ymax="<<ymax<<"\n";
                double padYmin = ymin - 0.1*(ymax-ymin); double padYmax = ymax + 0.1*(ymax-ymin);
                // draw ratios
                bool first=true; TLegend leg(0.6,0.6,0.95,0.9); leg.SetBorderSize(0); leg.SetFillStyle(0);
                for (size_t vi=0; vi<ratios.size(); ++vi) {
                    TH1D *r = ratios[vi]; if (!r) continue;
                    int color = colors[(vi+1) % colors.size()]; r->SetLineColor(color); r->SetMarkerColor(color);
                    if (first) { r->GetYaxis()->SetRangeUser(padYmin,padYmax); r->Draw("E1"); first=false; }
                    else r->Draw("E1same");
                    std::string label = (vi+1 < (size_t)allLabels.size()) ? allLabels[vi+1] : ("var"+std::to_string(vi+1));
                    leg.AddEntry(r, label.c_str(), "lep");
                }
                // unity line
                TLine unity(hdef->GetXaxis()->GetXmin(), 1.0, hdef->GetXaxis()->GetXmax(), 1.0); unity.SetLineStyle(2); unity.SetLineColor(kGray+2); unity.Draw("same");
                if (!leg.GetListOfPrimitives()->IsEmpty()) leg.Draw();
                // label
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03); tx.DrawLatex(0.12,0.92, (jetBins[ir]+" "+ (ic < (int)rapidityLabels.size() ? rapidityLabels[ic] : std::string(""))).c_str());
                if (validRatios==0) {
                    TLatex t2; t2.SetNDC(); t2.SetTextSize(0.04); t2.SetTextColor(kRed); t2.DrawLatex(0.35,0.5,"no variation data (empty)");
                }
                // cleanup ratios
            }
        }
        std::string outname = join_path(outp, std::string("compare_")+factor+".png"); cAll->SaveAs(outname.c_str()); 
    }

    std::cout<<"Comparison plots saved to "<< outp <<"\n";

    // --- PID comparison: compute variation/default ratios for Kaon/Pion/Combined and draw multi-panel canvases
    std::vector<std::string> pidTypes = {"Kaon","Pion","Combined"};
    for (const auto &ptype : pidTypes) {
        TCanvas *cPidAll = new TCanvas((std::string("c_pid_comp_")+ptype).c_str(), (std::string("PID ratio ")+ptype).c_str(), 1600, 1200);
        int rows = (int)jetBins.size(); int cols = nEta; if (rows<1||cols<1) { continue; }
        cPidAll->Divide(cols, rows);
        for (int ir=0; ir<rows; ++ir) {
            for (int ic=0; ic<cols; ++ic) {
                int pad = ir*cols + ic + 1; cPidAll->cd(pad);
                // default graph
                TGraphErrors *gdef = nullptr;
                if (ptype=="Kaon") gdef = gKaon_all[0][ir][ic];
                if (ptype=="Pion") gdef = gPion_all[0][ir][ic];
                if (ptype=="Combined") gdef = gComb_all[0][ir][ic];
                if (!gdef) { TLatex t; t.SetNDC(); t.SetTextColor(kRed); t.DrawLatex(0.3,0.5,"no default PID graph"); continue; }
                double ymin=1e99, ymax=-1e99; TLegend leg(0.6,0.6,0.95,0.9); leg.SetBorderSize(0); leg.SetFillStyle(0);
                int colorIdx=0;
                for (int v=1; v<nVar; ++v) {
                    TGraphErrors *gvar=nullptr;
                    if (ptype=="Kaon") gvar = gKaon_all[v][ir][ic];
                    if (ptype=="Pion") gvar = gPion_all[v][ir][ic];
                    if (ptype=="Combined") gvar = gComb_all[v][ir][ic];
                    if (!gvar) continue;
                    // build ratio graph: clone var and replace y points with yvar/ydef
                    TGraphErrors *gr = (TGraphErrors*)gvar->Clone(); gr->SetName((std::string("ratio_pid_")+ptype+"_v"+std::to_string(v)+"_pt"+std::to_string(ir)+"_eta"+std::to_string(ic)).c_str());
                    int ndef = gdef->GetN(); int nvar = gvar->GetN();
                    int nmin = std::min(ndef, nvar);
                    double x_def, y_def, x_var, y_var;
                    for (int k=0; k<nmin; ++k) {
                        gdef->GetPoint(k, x_def, y_def);
                        gvar->GetPoint(k, x_var, y_var);
                        double yr = safeDiv(y_var, y_def);
                        double errr = 0.0;
                        double ey_var = gvar->GetErrorY(k);
                        double ey_def = gdef->GetErrorY(k);
                        if (y_var!=0 && y_def!=0) {
                            double rel_var = (ey_var / (y_var==0?1.0:y_var));
                            double rel_def = (ey_def / (y_def==0?1.0:y_def));
                            errr = std::abs(yr) * std::sqrt(rel_var*rel_var + rel_def*rel_def);
                        }
                        gr->SetPoint(k, x_var, yr);
                        gr->SetPointError(k, gvar->GetErrorX(k), errr);
                        if (yr<ymin) ymin=yr; if (yr>ymax) ymax=yr;
                    }
                    int color = colors[(++colorIdx) % colors.size()]; gr->SetLineColor(color); gr->SetMarkerColor(color); gr->SetMarkerStyle(20+colorIdx);
                    if (colorIdx==1) { gr->Draw("AP"); 
                        if (ptype=="Kaon") gr->GetYaxis()->SetRangeUser(0.8, 1.2);
                        else if (ptype=="Pion") gr->GetYaxis()->SetRangeUser(0.5, 1.5);
                        else gr->GetYaxis()->SetRangeUser(0.3, 1.7);}
                    else gr->Draw("Psame");
                    std::string label = (v < (int)allLabels.size()) ? allLabels[v] : (std::string("var") + std::to_string(v));
                    leg.AddEntry(gr, label.c_str(), "lep");
                }
                // unity line
                TLine* unity = new TLine(gdef->GetXaxis()->GetXmin(), 1.0, gdef->GetXaxis()->GetXmax(), 1.0); unity->SetLineStyle(2); unity->SetLineColor(kGray+2); unity->Draw("same");
                if (!leg.GetListOfPrimitives()->IsEmpty()) leg.Draw();
                TLatex* tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.03); tx->DrawLatex(0.12,0.92, (jetBins[ir]+" "+ (ic < (int)rapidityLabels.size() ? rapidityLabels[ic] : std::string(""))).c_str());
            }
        }
        std::string outpid_all = join_path(outp, std::string("pid_compare_ratio_")+ptype+".png"); cPidAll->SaveAs(outpid_all.c_str());
    }
}
// End
