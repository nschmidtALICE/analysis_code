// ROOT macro: compare unfolded distributions between two unfold_new.cpp outputs
//
// Usage (ROOT):
//   root -l -q 'compare_unfolded_outputs.C+("pp.root","pPb.root","pp","pPb","cmp", "unfolded_zT", 6, true)'
//
// This scans both ROOT files, finds matching TH1 histograms whose *base name*
// starts with histBasePrefix (default: "unfolded_zT"), optionally restricts to
// a specific iteration (suffix "_iterN"), and writes overlay+ratio plots.

#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TDirectoryFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <ctime>

namespace
{

static std::string BaseName(const std::string &p)
{
    const auto pos = p.find_last_of('/');
    return (pos == std::string::npos) ? p : p.substr(pos + 1);
}

static bool StartsWith(const std::string &s, const std::string &prefix)
{
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

static bool EndsWith(const std::string &s, const std::string &suffix)
{
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

static std::string SanitizeForFilename(std::string s)
{
    for (char &c : s)
    {
        if (c == '/' || c == '\\')
            c = '_';
        else if (!(std::isalnum((unsigned char)c) || c == '_' || c == '-' || c == '.'))
            c = '_';
    }
    // avoid very long filenames
    const size_t kMax = 200;
    if (s.size() > kMax)
        s.resize(kMax);
    return s;
}

static std::string TodayYYYYMMDD()
{
    std::time_t t = std::time(nullptr);
    std::tm *lt = std::localtime(&t);
    char buf[32] = {0};
    if (lt)
        std::strftime(buf, sizeof(buf), "%Y-%m-%d", lt);
    else
        std::snprintf(buf, sizeof(buf), "unknown-date");
    return std::string(buf);
}

struct Chi2Result
{
    double chi2 = 0.0;
    int ndf = 0;
    double redchi2 = -1.0;
};

// symmetric chi2 using (ea^2 + eb^2) in denominator
static Chi2Result ComputeChi2Sym(const TH1 *a, const TH1 *b)
{
    Chi2Result r;
    if (!a || !b)
        return r;
    const int nb = std::min(a->GetNbinsX(), b->GetNbinsX());
    for (int i = 1; i <= nb; ++i)
    {
        const double va = a->GetBinContent(i);
        const double vb = b->GetBinContent(i);
        const double ea = a->GetBinError(i);
        const double eb = b->GetBinError(i);
        const double e2 = ea * ea + eb * eb;
        if (e2 > 0)
        {
            const double d = va - vb;
            r.chi2 += d * d / e2;
            ++r.ndf;
        }
    }
    if (r.ndf > 0)
        r.redchi2 = r.chi2 / r.ndf;
    return r;
}

static void CollectTH1Paths(TDirectory *dir, const std::string &pathPrefix, std::vector<std::string> &out)
{
    if (!dir)
        return;

    TIter next(dir->GetListOfKeys());
    while (TKey *key = (TKey *)next())
    {
        const char *className = key->GetClassName();
        TClass *cl = gROOT->GetClass(className);
        if (!cl)
            continue;

        if (cl->InheritsFrom(TDirectory::Class()))
        {
            TDirectory *sub = dynamic_cast<TDirectory *>(dir->Get(key->GetName()));
            if (!sub)
                continue;
            CollectTH1Paths(sub, pathPrefix + sub->GetName() + "/", out);
        }
        else if (cl->InheritsFrom(TH1::Class()))
        {
            out.push_back(pathPrefix + key->GetName());
        }
    }
}

static void SetHistStyle(TH1 *h, int color, int marker)
{
    if (!h)
        return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(0.9);
    h->SetLineWidth(2);
}

} // namespace

void compare_unfolded_outputs(const char *fileA,
                             const char *fileB,
                             const char *labelA = "A",
                             const char *labelB = "B",
                             const char *outPrefix = "compare_unfolded",
                             const char *histBasePrefix = "unfolded_zT",
                             int iter = -1,
                             bool normalize = true)
{
    gStyle->SetOptStat(0);

    TFile *fA = TFile::Open(fileA, "READ");
    TFile *fB = TFile::Open(fileB, "READ");
    if (!fA || fA->IsZombie())
    {
        std::cerr << "ERROR: cannot open fileA: " << fileA << std::endl;
        return;
    }
    if (!fB || fB->IsZombie())
    {
        std::cerr << "ERROR: cannot open fileB: " << fileB << std::endl;
        return;
    }

    const std::string dateStr = TodayYYYYMMDD();
    const std::string outDir = std::string(outPrefix) + "_" + dateStr;
    gSystem->mkdir(outDir.c_str(), true);

    const std::string outRootPath = outDir + "/" + std::string(outPrefix) + ".root";
    TFile *fOut = TFile::Open(outRootPath.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie())
    {
        std::cerr << "ERROR: cannot create output file: " << outRootPath << std::endl;
        return;
    }

    const std::string outTxtPath = outDir + "/" + std::string(outPrefix) + "_summary.txt";
    std::ofstream summary(outTxtPath);
    summary << "# compare_unfolded_outputs summary\n";
    summary << "# fileA: " << fileA << " (" << labelA << ")\n";
    summary << "# fileB: " << fileB << " (" << labelB << ")\n";
    summary << "# histBasePrefix: " << (histBasePrefix ? histBasePrefix : "") << "\n";
    summary << "# iter: " << iter << " (use -1 for all)\n";
    summary << "# normalize: " << (normalize ? "true" : "false") << "\n\n";
    summary << "name\tintA\tintB\tintB_over_intA\tchi2_ndf\tndf\n";

    std::vector<std::string> aPaths, bPaths;
    CollectTH1Paths(fA, "", aPaths);
    CollectTH1Paths(fB, "", bPaths);

    std::set<std::string> aSet(aPaths.begin(), aPaths.end());
    std::set<std::string> bSet(bPaths.begin(), bPaths.end());

    std::vector<std::string> common;
    common.reserve(std::min(aSet.size(), bSet.size()));
    std::set_intersection(aSet.begin(), aSet.end(), bSet.begin(), bSet.end(), std::back_inserter(common));

    const std::string prefix = histBasePrefix ? std::string(histBasePrefix) : std::string();
    const std::string iterSuffix = (iter >= 1) ? (std::string("_iter") + std::to_string(iter)) : std::string();

    int nCompared = 0;
    int nMatched = 0;

    for (const std::string &path : common)
    {
        const std::string base = BaseName(path);
        if (!prefix.empty() && !StartsWith(base, prefix))
            continue;
        ++nMatched;
        if (iter >= 1 && !EndsWith(base, iterSuffix))
            continue;

        TH1 *hA = dynamic_cast<TH1 *>(fA->Get(path.c_str()));
        TH1 *hB = dynamic_cast<TH1 *>(fB->Get(path.c_str()));
        if (!hA || !hB)
            continue;

        // Clone into memory so we can scale/style safely
        TH1 *a = dynamic_cast<TH1 *>(hA->Clone((base + "_A").c_str()));
        TH1 *b = dynamic_cast<TH1 *>(hB->Clone((base + "_B").c_str()));
        if (!a || !b)
            continue;
        a->SetDirectory(nullptr);
        b->SetDirectory(nullptr);

        const double intA0 = a->Integral(1, a->GetNbinsX());
        const double intB0 = b->Integral(1, b->GetNbinsX());

        if (normalize)
        {
            if (intA0 > 0)
                a->Scale(1.0 / intA0);
            if (intB0 > 0)
                b->Scale(1.0 / intB0);
        }

        SetHistStyle(a, kBlue + 1, 20);
        SetHistStyle(b, kRed + 1, 24);

        // chi2 computed on the (possibly normalized) hists
        Chi2Result chi = ComputeChi2Sym(a, b);

        // ratio B/A
        TH1 *ratio = dynamic_cast<TH1 *>(b->Clone((base + "_ratio_BoverA").c_str()));
        ratio->SetDirectory(nullptr);
        ratio->Divide(a);
        ratio->SetLineColor(kBlack);
        ratio->SetMarkerColor(kBlack);
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(0.8);
        ratio->GetYaxis()->SetTitle("B / A");

        // Make plot
        const std::string safe = SanitizeForFilename(path);
        TCanvas *c = new TCanvas(("c_" + safe).c_str(), base.c_str(), 900, 800);
        c->cd();

        TPad *pTop = new TPad(("pTop_" + safe).c_str(), "", 0.0, 0.30, 1.0, 1.0);
        TPad *pBot = new TPad(("pBot_" + safe).c_str(), "", 0.0, 0.00, 1.0, 0.30);
        pTop->SetBottomMargin(0.02);
        pTop->SetLeftMargin(0.12);
        pTop->SetRightMargin(0.03);
        pBot->SetTopMargin(0.02);
        pBot->SetBottomMargin(0.33);
        pBot->SetLeftMargin(0.12);
        pBot->SetRightMargin(0.03);
        pTop->Draw();
        pBot->Draw();

        pTop->cd();
        a->GetYaxis()->SetTitle(normalize ? "1/N dN/dz_{T}" : "Counts");
        a->GetXaxis()->SetLabelSize(0);
        a->GetXaxis()->SetTitleSize(0);

        // y-range
        double ymax = 0.0;
        for (int i = 1; i <= a->GetNbinsX(); ++i)
            ymax = std::max(ymax, a->GetBinContent(i) + a->GetBinError(i));
        for (int i = 1; i <= b->GetNbinsX(); ++i)
            ymax = std::max(ymax, b->GetBinContent(i) + b->GetBinError(i));
        if (ymax > 0)
            a->SetMaximum(1.25 * ymax);

        a->Draw("E1");
        b->Draw("E1 SAME");

        TLegend *leg = new TLegend(0.55, 0.70, 0.94, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(a, labelA, "lep");
        leg->AddEntry(b, labelB, "lep");
        leg->Draw();

        TLatex tx;
        tx.SetNDC(true);
        tx.SetTextSize(0.032);
        tx.DrawLatex(0.14, 0.92, base.c_str());
        tx.DrawLatex(0.14, 0.87, Form("chi2/ndf = %.3f (%d)", chi.redchi2, chi.ndf));
        tx.DrawLatex(0.14, 0.82, Form("IntA=%.3g  IntB=%.3g", intA0, intB0));

        pBot->cd();
        ratio->GetXaxis()->SetTitle("z_{T}");
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetLabelSize(0.10);
        ratio->GetYaxis()->SetTitleOffset(0.55);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->SetMinimum(0.5);
        ratio->SetMaximum(1.5);
        ratio->Draw("E1");
        TLine l1(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
        l1.SetLineStyle(2);
        l1.SetLineColor(kGray + 2);
        l1.Draw();

        // Save outputs
        const std::string outPng = outDir + "/" + SanitizeForFilename(base) + ".png";
        const std::string outPdf = outDir + "/" + SanitizeForFilename(base) + ".pdf";
        c->SaveAs(outPng.c_str());
        c->SaveAs(outPdf.c_str());

        fOut->cd();
        c->Write();
        ratio->Write();

        const double intRatio = (intA0 > 0 ? intB0 / intA0 : 0.0);
        summary << base << "\t" << intA0 << "\t" << intB0 << "\t" << intRatio << "\t" << chi.redchi2 << "\t" << chi.ndf << "\n";

        delete leg;
        delete pTop;
        delete pBot;
        delete c;
        delete a;
        delete b;
        delete ratio;

        ++nCompared;
    }

    std::cout << "Output directory: " << outDir << std::endl;
    std::cout << "Wrote: " << outRootPath << std::endl;
    std::cout << "Wrote: " << outTxtPath << std::endl;
    std::cout << "Matched TH1 paths (prefix filter applied later): " << nMatched << std::endl;
    std::cout << "Compared histograms: " << nCompared << std::endl;

    summary.close();
    fOut->Close();
    fA->Close();
    fB->Close();
}
