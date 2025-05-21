from json import load
import ROOT
import numpy as np

class Ratio:
    """
    Class to return the tracking ratio efficiency corrections (data/MC).
    """
    def __init__(self, path):
        """
        Load the low and high pT JSON efficiency files.
        """
        with open("{}2016_ratio_trackefflo.json".format(path)) as f: self.ratLo = load(f)
        with open("{}2016_ratio_trackeffhi.json".format(path)) as f: self.ratHi = load(f)
        
    def _idx(self, bins, pre, val, over):
        """
        Return the bin index for a given binning scheme.
        
        bins: binning map.
        pre:  prefix of the variable to find the bin.
        val:  value of the variable.
        over: return the maximum bin if above bin range.
        """
        keys = [key for key in bins if pre + "_" in key]
        for idx in range(len(keys)):
            key = "%s_%i" % (pre, idx)
            if bins[key]["min"] <= val < bins[key]["max"]: return idx
        if over and val > bins[key]["max"]: return idx
        else: return 0

    def _eff(self, xyEff, xVal, yVal, xKey, yKey, xOver, yOver):
        """
        Return the efficiency for two variables.

        xyEff: the efficiency map in terms of x and y.
        xVal:  the value for x.
        yVal:  the value for y.
        xKey:  the bin key for x.
        yKey:  the bin key for y.
        xOver: use the overflow bin for the x value.
        yOver: use the overflow bin for the y value.
        """
        info = {"value": 0.0, "error_low": 0.0}
        xBins = xyEff["binning"]
        xIdx  = self._idx(xBins, xKey + "_bin", xVal, xOver)
        xKey  = "%s_bin_%i" % (xKey, xIdx)
        yBins = xyEff["binning"][xKey]
        yIdx  = self._idx(yBins, yKey + "_bin", yVal, yOver)
        yKey  = "%s_bin_%i" % (yKey, yIdx)
        return xyEff["efficiencies"][xKey][yKey]

    def __call__(self, eta, pt, p):
        """
        Return the track correction efficiency.

        eta: pseudo-rapidity.
        pt:  transverse momentum in MeV.
        p:   momumentum in MeV.
        """
        if pt > 20000:
            return self._eff(self.ratHi, eta, pt, "eta", "pt", False, True)
        else:
            return self._eff(self.ratLo, eta, p, "eta", "p", False, True)
    #- - - - - - - - - - - - - - - - - - - - - - -
    def getCorr(self, eta, pt, p):
        
        efflist = self(eta, pt, p)
        #print("eta: {}, pT {}, p: {} | efff: {}".format(eta, pt, p,efflist["value"]))

        return efflist["value"]
    #- - - - - - - - - - - - - - - - - - - - - - -
    def getCorr_Var(self, eta, pt, p):
        
        efflist = self(eta, pt, p)
        value         = efflist["value"]
        meanStatError = 0.5*(efflist["error_low"]+efflist["error_high"])
        sysError      = efflist["systematic"]
        totErr        = (meanStatError**2+sysError**2)**0.5
        value_var     = value-totErr

        return value_var
    #- - - - - - - - - - - - - - - - - - - - - - -
    def getCorr_histo(self):
                
        xBinsJson = self.ratLo["binning"]
        #xBinsJson = self.ratHi["binning"]
        
        xBins = []
        yBins = []
        
        #-Extract the bin ranges
        for i in range(0,5):
            try:
                etaBin = xBinsJson["eta_bin_{}".format(i)]
                if len(yBins)==0: yBins.append(etaBin["min"])
                yBins.append(etaBin["max"])
                if i==0:
                    for j in range(0,5):
                        try:
                            #pTBin = etaBin["pt_bin_{}".format(j)] #high
                            pTBin = etaBin["p_bin_{}".format(j)] #-low
                            if len(xBins)==0: xBins.append(pTBin["min"])
                            xBins.append(pTBin["max"])
                        except:
                            continue
                            #print("pT_bin_{} in eta_Bin{} DOES NOT EXIST".format(j,i))
            except:
                continue
                #print("eta_bin_{} DOES NOT EXIST".format(i))

        xBinsArr = np.array(xBins)
        yBinsArr = np.array(yBins)

        #-Create histogram with set pT eta binning
        histo = ROOT.TH2F("Data_MCRatio_histo","Data_MCRatio_histo",len(xBins)-1,xBinsArr,len(yBins)-1,yBinsArr)

        #-Fill the histogram
        for i in range(0,len(xBins)):
            for j in range(0,len(yBins)):
                xVal = histo.GetXaxis().GetBinCenter(i)
                yVal = histo.GetYaxis().GetBinCenter(j)
                #eff  = self.getCorr(eta=yVal, pt=xVal, p=20000-1) #-high
                eff  = self.getCorr(eta=yVal, pt=0, p=xVal) #-low
                histo.SetBinContent(i,j,eff)
        return histo
       
if __name__ == "__main__":

    #-Create the object.
    jSonPath = "/Users/eliane/LHCb/x3872-code/filterGangaOutput/R_Data_MC/"
    ratio = Ratio(jSonPath)

    #-Testing
    #corrVal     = ratio.getCorr(eta = 2.4, pt = 5000, p = 10000)
    #corrVal_Var = ratio.getCorr_Var(eta = 2.4, pt = 5000, p = 10000)

    #-Create a histogram of the two maps
    histo = ratio.getCorr_histo()
    
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    
    myBlankHisto2 = histo.Clone("Blank")
    myBlankHisto2.Clear()
    myBlankHisto2.GetXaxis().SetNdivisions(505)
    myBlankHisto2.GetXaxis().SetTitleSize(0.05)
    myBlankHisto2.GetYaxis().SetTitleOffset(1.35)
    myBlankHisto2.GetYaxis().SetTitleSize(0.055)
    myBlankHisto2.SetYTitle("p_{T}")
    myBlankHisto2.SetXTitle("eta")
  
    
    c1 = ROOT.TCanvas("c1","c1: hist",600,600)
    c1.cd(1)
    myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
    myPad.SetLeftMargin(0.17)
    myPad.SetTopMargin(0.07)
    myPad.SetRightMargin(0.14)
    myPad.SetBottomMargin(0.15)
    myPad.Draw()
    myPad.cd()
    
    histo.GetXaxis().SetNdivisions(505)
    histo.GetXaxis().SetTitleSize(0.05)
    histo.GetYaxis().SetTitleOffset(1.35)
    histo.GetYaxis().SetTitleSize(0.055)
    #histo.SetXTitle("p_{T} [GeV]")
    histo.SetXTitle("p [GeV]")
    histo.SetYTitle("eta")
    histo.DrawCopy("colz")
    
    c1.cd()
    c1.SaveAs("./CorrHisto.pdf")

