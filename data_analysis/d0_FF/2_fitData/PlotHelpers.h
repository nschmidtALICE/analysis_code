#ifndef PLOTHELPERS_H
#define PLOTHELPERS_H

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <sstream>
#include <sys/stat.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TH1.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TF1.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

// Helper for TLegend setup
inline void setupLegend(TLegend* legend, double textSize, double margin, int font, int border, int fillStyle, int fillColor) {
    if (!legend) return;
    legend->SetTextFont(font);
    legend->SetBorderSize(border);
    legend->SetFillStyle(fillStyle);
    legend->SetFillColor(fillColor);
    legend->SetMargin(margin);
    legend->SetTextSize(textSize);
}

// Helper for canvas and pad setup
inline void setupCanvasAndPad(TCanvas* canvas, TPad* pad, double left, double bottom, double right) {
    if (pad) {
        pad->SetLeftMargin(left);
        pad->SetBottomMargin(bottom);
        pad->SetRightMargin(right);
        pad->Draw();
        pad->cd();
    } else if (canvas) {
        canvas->SetLeftMargin(left);
        canvas->SetBottomMargin(bottom);
        canvas->SetRightMargin(right);
        canvas->Draw();
        canvas->cd();
    }
}

// Helper for histogram styling
inline void styleHistogram(TH1* hist, int color, int lineWidth, int markerStyle, int markerColor, double markerSize) {
    if (!hist) return;
    hist->SetLineColor(color);
    hist->SetLineWidth(lineWidth);
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerColor(markerColor);
    hist->SetMarkerSize(markerSize);
}

// Helper for debug output
inline void debugLog(const std::string& msg) {
    std::cout << "[Plotter] " << msg << std::endl;
}

// Helper for TPad margin setup
inline void setupPadMargins(TPad* pad, double left, double bottom, double right, double top) {
    if (!pad) return;
    pad->SetLeftMargin(left);
    pad->SetBottomMargin(bottom);
    pad->SetRightMargin(right);
    pad->SetTopMargin(top);
    pad->Draw();
}

// Helper method to ensure a directory exists
inline void ensureDirectoryExists(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        std::cout << "Creating directory: " << path << std::endl;
        int result = mkdir(path.c_str(), 0755);
        if (result != 0) {
            std::cerr << "Error creating directory: " << path << std::endl;
        }
    }
}

// Helper to generate a unique identifier for plot objects
inline std::string getUniqueId() {
    std::time_t timestamp = std::time(nullptr);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(10000, 99999);
    int rand_val = dis(gen);
    std::stringstream ss;
    ss << timestamp << "_" << rand_val;
    return ss.str();
}

// Helper to add parameters to legend from RooArgSet
inline void addParameterToLegend(TLegend* legend, RooArgSet* params, const std::string& paramName,
                                const std::string& displayName, const std::string& format, bool skipIfMissing = false) {
    if (!legend || !params) return;
    RooRealVar* param = dynamic_cast<RooRealVar*>(params->find(paramName.c_str()));
    if (param || !skipIfMissing) {
        if (param) {
            legend->AddEntry((TObject*)nullptr, Form((displayName + ": " + format).c_str(),
                           param->getVal(), param->getError()), "");
        } else if (!skipIfMissing) {
            legend->AddEntry((TObject*)nullptr, (displayName + ":").c_str(), "");
        }
    }
}

#endif // PLOTHELPERS_H
