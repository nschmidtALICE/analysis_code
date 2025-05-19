import os
import ROOT
#import parameters as pars
from array import *
import numpy
import argparse
import math
import sys
import argparse
import itertools
import math
from array import *
# put in LHCb style here...

class Plotter:
    def __init__(self,resonance,basepath="./",bin=0,binned=False, range="", name=""):
        #self.nameOfOutputFile = nameOfOutputFile
        self.resonance = resonance
        self.FitBin    = bin
        self.isBinned  = binned
        self.range     = range
        self.NameKey   = name
        if "z" in self.range:
          self.obsSelection = "zT"
        else:
          self.obsSelection = "dR"

        self.format    = "png" #pdf
        #self.format    = "pdf" #pdf
        #self.format    = "C" #pdf
        if basepath.endswith("/"):
          self.basepath  = basepath
        else:
          self.basepath  = "{}/".format(basepath)
        print("Save all histograms to: {}".format(self.basepath))
        if not os.path.exists(self.basepath):
          os.makedirs(self.basepath)
          print("make new directory: {}".format(self.basepath))
        if not os.path.exists("{}MassFits_{}/".format(self.basepath,self.obsSelection)):
          os.makedirs("{}MassFits_{}/".format(self.basepath,self.obsSelection))
        if not os.path.exists("{}TimeFits_{}/".format(self.basepath,self.obsSelection)):
          os.makedirs("{}TimeFits_{}/".format(self.basepath,self.obsSelection))

    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def extractCorParamRnd(self, idString, mass_tag_measured, myLegendList, data, data2=None, data3=None):


      #-Rnd selection correction
      canvas_Rnd     = ROOT.TCanvas("cMulti_Rnd","cMultiCorr", 800, 800)
      canvasNorm_Rnd = ROOT.TCanvas("cMultiNorm_Rnd","cMultiCorrNorm", 800, 800)
      canvasDiv_Rnd  = ROOT.TCanvas("cMultiDiv_Rnd","cMultiCorrDiv", 800, 800)

      ROOT.TGaxis.SetMaxDigits(2)
      ROOT.gStyle.SetOptStat(0)
      # Construct plot frame in 'x=mass_tag_measured'
      
      #Plot data in the frame, set histrogram properties
      #data.plotOn(frame, ROOT.RooFit.Name("datahistogram"), ROOT.RooFit.LineWidth(1), ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.MarkerStyle(20))
      h_data    = data.createHistogram("h_MassSignalRnd_data",mass_tag_measured)
      h_data.Rebin(100)
      maxList   = [h_data.GetMaximum(),0]
      if data2:
        h_data2 = data2.createHistogram("h_MassSignalRnd_data2",mass_tag_measured)
        h_data2.Rebin(100)
        maxList = [h_data.GetMaximum(),h_data2.GetMaximum()]
      if data3:
        h_data3 = data3.createHistogram("h_MassSignalRnd_data3",mass_tag_measured)
        h_data3.Rebin(100)
        maxList = [h_data.GetMaximum(),h_data2.GetMaximum(),h_data3.GetMaximum()]

   
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      canvasNorm_Rnd.cd()
      myPad2 = ROOT.TPad("myPad2Rnd", "The padRnd",0,0,1,1)
      myPad2.SetLeftMargin(0.8)
      myPad2.SetTopMargin(0.08)
      myPad2.SetRightMargin(0.1)
      myPad2.SetBottomMargin(0.15)
      myPad2.Draw()#d not workoes
      
      h_data.SetMarkerStyle(20)
      h_data.SetMarkerSize(0.5)
      h_data.SetMarkerColor(1)
      h_data.SetLineColor(ROOT.kGray+2)
      h_data.DrawNormalized("E")
      
      if data2:
        h_data2.SetMarkerStyle(20)
        h_data2.SetMarkerSize(0.5)
        h_data2.SetMarkerColor(ROOT.kRed+2)
        h_data2.SetLineColor(ROOT.kRed+2)
        h_data2.DrawNormalized("same E")
      if data3:
        h_data3.SetMarkerStyle(20)
        h_data3.SetMarkerSize(0.5)
        h_data3.SetMarkerColor(ROOT.kBlue+2)
        h_data3.SetLineColor(ROOT.kBlue+2)
        h_data3.DrawNormalized("same E")
      
      canvasNorm_Rnd.cd()
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      canvas_Rnd.cd()
      myPad = ROOT.TPad("myPadRnd", "The padRnd",0,0,1,1)
      myPad.SetLeftMargin(0.8)
      myPad.SetTopMargin(0.08)
      myPad.SetRightMargin(0.1)
      myPad.SetBottomMargin(0.15)
      myPad.Draw()#d not workoes
      
      h_data.GetYaxis().SetRangeUser(0,max(maxList)*1.3)
      h_data.SetMarkerStyle(20)
      h_data.SetMarkerSize(0.5)
      h_data.SetMarkerColor(1)
      h_data.SetLineColor(ROOT.kGray+2)
      h_data.DrawCopy("E")
      
      if data2:
        h_data2.SetMarkerStyle(20)
        h_data2.SetMarkerSize(0.5)
        h_data2.SetMarkerColor(ROOT.kRed+2)
        h_data2.SetLineColor(ROOT.kRed+2)
        h_data2.DrawCopy("same E")
      if data3:
        h_data3.SetMarkerStyle(20)
        h_data3.SetMarkerSize(0.5)
        h_data3.SetMarkerColor(ROOT.kBlue+2)
        h_data3.SetLineColor(ROOT.kBlue+2)
        h_data3.DrawCopy("same E")
        #frame.Draw("E")
      canvas_Rnd.cd()

      h_data.GetYaxis().UnZoom()
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      #-Ratio
      ROOT.TGaxis.SetMaxDigits(4)
      canvasDiv_Rnd.cd()
      
  
      if data2:
        myPad3 = ROOT.TPad("myPad3Rnd", "The pad3Rnd",0,0,1,1)
        myPad3.SetLeftMargin(0.8)
        myPad3.SetTopMargin(0.08)
        myPad3.SetRightMargin(0.1)
        myPad3.SetBottomMargin(0.15)
        myPad3.Draw()
        
        h_data2.Divide(h_data)
        h_data2.SetMarkerStyle(20)
        h_data2.SetMarkerSize(0.5)
        h_data2.SetMarkerColor(1)
        h_data2.SetLineColor(ROOT.kGray+2)
        #h_data.DrawNormalized("E")
        h_data2.DrawCopy("E")
        
        #Fit a horizontal line and save the value for the efficiency correction
        minL =h_data.GetBinCenter(0)
        maxL =h_data.GetBinCenter(h_data.GetNbinsX())
        line = ROOT.TF1("Line","[0]",minL,maxL)
        h_data2.Fit("Line","NQ","",minL,maxL)
        fitValue    = line.GetParameter(0)
        fitValueErr = line.GetParError(0)
        
        line.DrawCopy("same")
        #myLegendList[bin].AddEntry(frame," Correction Constant %2.3f#pm%2.3f" % (fitValue,fitValueErr),"")
        #myLegendList[bin].Draw()
        myLegendList[0].AddEntry(h_data2," Correction Constant %2.3f#pm%2.3f" % (fitValue,fitValueErr),"")
        myLegendList[0].Draw()
                
        canvasDiv_Rnd.cd()

      #print("will save canvas 1")
      canvas_Rnd.SaveAs("{}CorrFac{}_zT.png".format(self.basepath,idString))
      canvasDiv_Rnd.SaveAs("{}CorrFac{}_zTRatio.png".format(self.basepath,idString))
      canvasNorm_Rnd.SaveAs("{}CorrFac{}_zTNorm.png".format(self.basepath,idString))
    

      return fitValue, fitValueErr
 
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def extractCorParam(self, idString, bin, type, canvas, canvas_Norm, canvasDiv, corrValue, corrValueErr, myLegendList, mass_tag_measured, data, data2=None, data3=None):
      
      ROOT.TGaxis.SetMaxDigits(2)
      ROOT.gStyle.SetOptStat(0)
      # Construct plot frame in 'x=mass_tag_measured'
      frame = mass_tag_measured.frame(ROOT.RooFit.Bins(40))
      frame.SetTitle("Signal for {}".format(self.range))
      
      #Plot data in the frame, set histrogram properties
      #data.plotOn(frame, ROOT.RooFit.Name("datahistogram"), ROOT.RooFit.LineWidth(1), ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.MarkerStyle(20))
      h_data    = data.createHistogram("h_MassSignal{}_data".format(bin),mass_tag_measured)
      maxList   = [h_data.GetMaximum(),0]
      if data2:
        h_data2 = data2.createHistogram("h_MassSignal{}_data2".format(bin),mass_tag_measured)
        maxList = [h_data.GetMaximum(),h_data2.GetMaximum()]
      if data3:
        h_data3 = data3.createHistogram("h_MassSignal{}_data3".format(bin),mass_tag_measured)
        maxList = [h_data.GetMaximum(),h_data2.GetMaximum(),h_data3.GetMaximum()]

   
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      canvas_Norm[type].cd(bin+1)
      myPad2 = ROOT.TPad("myPad2{}".format(bin), "The pad{}".format(bin),0,0,1,1)
      myPad2.SetLeftMargin(0.8)
      myPad2.SetTopMargin(0.08)
      myPad2.SetRightMargin(0.1)
      myPad2.SetBottomMargin(0.15)
      myPad2.Draw()#d not workoes
      
      h_data.SetMarkerStyle(20)
      h_data.SetMarkerSize(0.5)
      h_data.SetMarkerColor(1)
      h_data.SetLineColor(ROOT.kGray+2)
      h_data.DrawNormalized("E")
      
      if data2:
        h_data2.SetMarkerStyle(20)
        h_data2.SetMarkerSize(0.5)
        h_data2.SetMarkerColor(ROOT.kRed+2)
        h_data2.SetLineColor(ROOT.kRed+2)
        h_data2.DrawNormalized("same E")
      if data3:
        h_data3.SetMarkerStyle(20)
        h_data3.SetMarkerSize(0.5)
        h_data3.SetMarkerColor(ROOT.kBlue+2)
        h_data3.SetLineColor(ROOT.kBlue+2)
        h_data3.DrawNormalized("same E")
      
      canvas_Norm[type].cd()
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      canvas[type].cd(bin+1)
      myPad = ROOT.TPad("myPad{}".format(bin), "The pad{}".format(bin),0,0,1,1)
      myPad.SetLeftMargin(0.8)
      myPad.SetTopMargin(0.08)
      myPad.SetRightMargin(0.1)
      myPad.SetBottomMargin(0.15)
      myPad.Draw()#d not workoes
      
      h_data.GetYaxis().SetRangeUser(0,max(maxList)*1.3)
      h_data.SetMarkerStyle(20)
      h_data.SetMarkerSize(0.5)
      h_data.SetMarkerColor(1)
      h_data.SetLineColor(ROOT.kGray+2)
      h_data.DrawCopy("E")
      
      if data2:
        h_data2.SetMarkerStyle(20)
        h_data2.SetMarkerSize(0.5)
        h_data2.SetMarkerColor(ROOT.kRed+2)
        h_data2.SetLineColor(ROOT.kRed+2)
        h_data2.DrawCopy("same E")
      if data3:
        h_data3.SetMarkerStyle(20)
        h_data3.SetMarkerSize(0.5)
        h_data3.SetMarkerColor(ROOT.kBlue+2)
        h_data3.SetLineColor(ROOT.kBlue+2)
        h_data3.DrawCopy("same E")
        #frame.Draw("E")
      canvas[type].cd()

      h_data.GetYaxis().UnZoom()
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      #- - - - - - - - - - - - - - - - - - - - - - - - -
      #-Ratio
      ROOT.TGaxis.SetMaxDigits(4)
      canvasDiv[type].cd(bin+1)
      
  
      if data2:
        myPad3 = ROOT.TPad("myPad3{}".format(bin), "The pad3{}".format(bin),0,0,1,1)
        myPad3.SetLeftMargin(0.8)
        myPad3.SetTopMargin(0.08)
        myPad3.SetRightMargin(0.1)
        myPad3.SetBottomMargin(0.15)
        myPad3.Draw()
        
        h_data2.Divide(h_data)
        h_data2.SetMarkerStyle(20)
        h_data2.SetMarkerSize(0.5)
        h_data2.SetMarkerColor(1)
        h_data2.SetLineColor(ROOT.kGray+2)
        #h_data.DrawNormalized("E")
        h_data2.DrawCopy("E")
        
        #Fit a horizontal line and save the value for the efficiency correction
        minL =h_data.GetBinCenter(0)
        maxL =h_data.GetBinCenter(h_data.GetNbinsX())
        line = ROOT.TF1("Line","[0]",minL,maxL)
        h_data2.Fit("Line","NQ","",minL,maxL)
        fitValue    = line.GetParameter(0)
        fitValueErr = line.GetParError(0)
        
        line.DrawCopy("same")
        myLegendList[bin].AddEntry(frame," Correction Constant %2.3f#pm%2.3f" % (fitValue,fitValueErr),"")
        myLegendList[bin].Draw()
        
        corrValue[type][bin]    = fitValue
        corrValueErr[type][bin] = fitValueErr
        
      if bin==(len(corrValue[type])-1):
        #print("will save canvas 1")
        canvas[type].SaveAs("{}CorrFac{}_zT.png".format(self.basepath,idString))
        #print("will save canvas 2")
        canvasDiv[type].SaveAs("{}CorrFac{}_zTRatio.png".format(self.basepath,idString))
        canvas_Norm[type].SaveAs("{}CorrFac{}_zTNorm.png".format(self.basepath,idString))
       
      canvasDiv[type].cd()

      return canvas[type], canvas_Norm[type], canvasDiv[type]
    
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def teffMap1D(self):
        pass


    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def teffMap2D(self,teff,effType,dataType,closureTF1=False,closureHisto=False):
        canvas = ROOT.TCanvas("canvas","canvas",0,0,600,400)

        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetNumberContours(100)
        ROOT.gPad.SetLogx();
        if effType == "trigger":
            teff.SetTitle("Trigger Efficiency; sqrt(pt(muon1)*pt(muon2)); eta(tag); eff")
        if effType == "probnnmu":
            if closureTF1 == False and closureHisto == False:
                teff.SetTitle("ProbNNmu Efficiency; pt(muon); eta(muon); eff")            
            elif closureTF1 == True or closureHisto == True: 
                teff.SetTitle("ProbNNmu Efficiency; pt(tag); eta(tag); eff")
        if effType == "reco":
            if closureHisto == True:
                teff.SetTitle("Muon reco * IsMuon Efficiency; pt(tag); eta(tag); eff")
        teff.Draw("colz")
        canvas.Update();

        effHist = teff.GetPaintedHistogram();
        effHist.SetMinimum(0); effHist.SetMaximum(1);
        effHist.Draw("colz")
        canvas.Update();
        canvas.SaveAs("{}/eff2DMap_{effType}Swap.pdf".format(self.basepath))
            
        return effHist
    
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def teffParam(self,effType,dataType,dtNum,dtDenom,length):
        effHistProjDenom = [0.0]*length
        effHistProjNum = [0.0]*length
        effHistProj = [0.0]*length
        effGraphProj = [0.0]*length
        func = [0.0]*length

        canvas = ROOT.TCanvas("canvas","canvas",0,0,600,400)
        canvas.Divide(2,5);
        for can in range(1,length):
            canvas.cd(can)
            ROOT.gStyle.SetOptStat(0)
            ROOT.gPad.SetLogx();

            if effType == "trigger":
                func[can] = ROOT.TF1("func{0}".format(can),"[0] - ([1]*exp(-1.0*[2]*x))",1000.0,11000.0);
                func[can].SetParLimits(0,0.3,0.9); func[can].SetParLimits(1,0.001,5.0);
                func[can].SetParLimits(2,0.000001,0.01);
            if effType == "probnnmu":
                func[can] = ROOT.TF1("func{0}".format(can),"[0] - ([1]*exp(-1.0*[2]*x))",500.0,10000.0);
                func[can].SetParLimits(0,0.8,1.2); func[can].SetParLimits(1,0.01,3.0);
                func[can].SetParLimits(2,0.0001,0.01);
                
            effHistProjNum[can] = dtNum.ProjectionX("effHistProjNum{0}".format(can),can,can)
            effHistProjDenom[can] = dtDenom.ProjectionX("effHistProjDenom{0}".format(can),can,can)
            effHistProj[can] = ROOT.TEfficiency(effHistProjNum[can],effHistProjDenom[can])

            effHistProj[can].Draw("")
            canvas.Update();
            effGraphProj[can] = effHistProj[can].GetPaintedGraph()
            effGraphProj[can].SetMinimum(0); effGraphProj[can].SetMaximum(1.1);
            if effType == "trigger": effGraphProj[can].GetXaxis().SetLimits(1000.,11000.);
            if effType == "probnnmu": effGraphProj[can].GetXaxis().SetLimits(500.,1e4);
            effGraphProj[can].GetXaxis().SetLabelSize(0.08);
            effGraphProj[can].GetYaxis().SetLabelSize(0.08);

            effGraphProj[can].Fit("func{}".format(can),"R EX0")
            effGraphProj[can].Draw("P")

            canvas.Update();
        
        canvas.SaveAs(os.path.join(self.basepath,'/effPtProj_{effType}Swap.pdf'))

        return func
        
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def individualMassFitPlot(self, sigYieldParam, extended_pdf, mass_tag_measured, data, fitTypeName):
        """
        Plot individual mass fit.
        """
        print(f"\n==== Creating mass fit plot for bin {self.FitBin} with model {fitTypeName} ====")
        
        # Add timestamp for unique identifier
        import time
        import random
        timestamp = int(time.time())
        rand_val = random.randint(10000, 99999)
        unique_id = f"{timestamp}_{rand_val}"
        
        # Check output directory
        output_dir = f"{self.basepath}MassFits_{self.obsSelection}/"
        print(f"  Output directory: {output_dir}")
        if not os.path.exists(output_dir):
            print(f"  Creating output directory: {output_dir}")
            os.makedirs(output_dir, exist_ok=True)
        
        # Check data
        print(f"  Data entries: {data.numEntries()}")
        
        # Extract all the parameters from binning Range
        binFull = mass_tag_measured.getBinning("fullRange")   # returns RooAbsBinning
        binSig  = mass_tag_measured.getBinning("signalRange") # returns RooAbsBinning
        print(f"  Mass range: {binFull.lowBound()}-{binFull.highBound()}")
        
        mass_tag_measured.setRange(binFull.lowBound(), binFull.highBound())
        frame = mass_tag_measured.frame(ROOT.RooFit.Bins(40))
        frame.SetTitle(f"Signal for {self.range}")
        
        print("  Plotting data...")
        data.plotOn(frame, ROOT.RooFit.Name("datahistogram"), 
                    ROOT.RooFit.LineWidth(1), ROOT.RooFit.MarkerSize(0.5), 
                    ROOT.RooFit.MarkerStyle(20))
        
        # Get histogram for fit assessment with unique name
        print("  Creating histogram from data...")
        h_data = data.createHistogram(f"h_MassSignal{self.FitBin}_{fitTypeName}_{unique_id}_data", mass_tag_measured)
        
        # Create a list to track objects that need cleanup
        histograms_to_cleanup = [h_data]
        
        print("  Plotting fit models...")
        try:
            # Plot total fit
            extended_pdf.plotOn(frame, ROOT.RooFit.LineWidth(2), ROOT.RooFit.Name("TotalFit"))
            
            # Plot background fit component
            extended_pdf.plotOn(frame, ROOT.RooFit.Components("bkg_pdf_ext"),
                              ROOT.RooFit.Name("BackgroundFit"),
                              ROOT.RooFit.LineStyle(ROOT.kDashed),
                              ROOT.RooFit.LineWidth(2))
            
            # Plot signal fit component
            extended_pdf.plotOn(frame, ROOT.RooFit.Components("sig_pdf_ext"),
                              ROOT.RooFit.Name("SignalFit"),
                              ROOT.RooFit.LineColor(ROOT.kRed),
                              ROOT.RooFit.LineStyle(ROOT.kDashed),
                              ROOT.RooFit.LineWidth(2))
        except Exception as e:
            print(f"  ERROR plotting fit components: {str(e)}")
            import traceback
            traceback.print_exc()
        
        # Save output file path
        output_file = f"{output_dir}Bin{self.FitBin}_{fitTypeName}{'Binned' if self.isBinned else 'Unbinned'}.{self.format}"
        print(f"  Will save plot to: {output_file}")
        
        try:
            # Draw and save plot
            print("  Creating canvas...")
            canvas = ROOT.TCanvas(f"canvas_{unique_id}", f"canvas_{unique_id}", 800*2, 400*2)
            
            # Set up pad, draw frame, etc.
            pad = ROOT.TPad(f"myPad_{unique_id}", f"The pad_{unique_id}", 0, 0, 1, 1)
            pad.SetLeftMargin(0.15)
            pad.SetBottomMargin(0.15)
            pad.Draw()
            pad.cd()
            
            frame.Draw()
            
            print(f"  Saving canvas to {output_file}")
            canvas.SaveAs(output_file)
            print("  Plot saved successfully!")
            
            # Clean up to prevent memory leaks
            canvas.Close()
            for hist in histograms_to_cleanup:
                if hist:
                    hist.SetDirectory(0)
                    
        except Exception as e:
            print(f"  ERROR saving plot: {str(e)}")
            import traceback
            traceback.print_exc()
        
        return h_data

    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def splotVals(self, resonance, extended_pdf, mass_tag_measured, data,ttree,ns,nb,fidCutString,sFile):

        splot = ROOT.RooStats.SPlot("splot","splot",data,extended_pdf,ROOT.RooArgList(ns,nb))
        data.Print("v")
        splot.Print()

        brs=[array.array("d",[0.]) for i in range(3)]
        #fvar = ROOT.TFile.Open("trialTree{}.root".format(fidCutString.translate ({ord(c): "" for c in " !@#$%^&*()[]{};:,/<>?\|`~-=_+"})), "recreate")
        fvar = ROOT.TFile.Open(sFile,"recreate")
        #removeSpecialChars = z.translate ({ord(c): "" for c in " !@#$%^&*()[]{};:,./<>?\|`~-=_+"})
        if fidCutString == "":
            Vartree = ttree.CopyTree("mass_tag_measured>{}&&mass_tag_measured<{}".format(mass_tag_measured.getMin(),mass_tag_measured.getMax()))
        else:
            Vartree = ttree.CopyTree("mass_tag_measured>{}&&mass_tag_measured<{} && {}".format(mass_tag_measured.getMin(),mass_tag_measured.getMax(),fidCutString))
        Vartree.SetName("ntuple")
        br1=Vartree.Branch("sig_sw",brs[0], "sig_sw/D")
        br2=Vartree.Branch("bkg_sw",brs[1], "bkg_sw/D")

        scale1,scale2 = [0.]*2
        #print(ttree.sumEntries())
        for i in range(int(data.sumEntries())):
        #for i in range(int(ttree.GetEntries())):
            row = data.get(i)
            brs[0][0] = row.find("sig_yield_sw").getVal()
            brs[1][0] = row.find("bkg_yield_sw").getVal()
            br1.Fill()
            br2.Fill()
            
            scale1 += brs[0][0]
            scale2 += brs[0][0]*brs[0][0]

        print(sFile)
        Vartree.Write("",1)
        fvar.Close()
        
        return scale1, scale2, (scale1/scale2)

    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def splotReweightDists(self, resonance, extended_pdf, data,ns,nb,pt_tag):

        splot = ROOT.RooStats.SPlot("splot","splot",data,extended_pdf,ROOT.RooArgList(ns,nb))
        data.Print("v")
        splot.Print()

        cdata = ROOT.TCanvas("sPlot", "sPlot demo", 400, 600);
         
        #RooDataSet *data = (RooDataSet *)ws->data("dataWithSWeights");
 
        #create weighted data set
        dataw_z = ROOT.RooDataSet(data.GetName(), data.GetTitle(), data, data.get(), "", "sig_yield_sw");
 
        #frame2 = rapidity.frame(2.0,4.5,20);
        frame2 = pt_tag.frame(2000.0,15000.,20);
        #Since the data are weighted, we use SumW2 to compute the errors.
        dataw_z.plotOn(frame2, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
 
        frame2.SetTitle("Pt tag distribution with s weights");
        frame2.Draw();
        
        cdata.SaveAs(os.path.join(self.basepath,'/sweightPttagRedPt2000.pdf'))

        tgRap = frame2.getHist()
        print(type(tgRap))

        histTFile = ROOT.TFile("sweightPttagRedPt2000.root","RECREATE")
        tgRap.Write()
        histTFile.Write()
        histTFile.Close()
        
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def organisedMassFitsPlot(self):
        pass

    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def plotMultiple(self, name, lifetime_tag, drawLogDiff, data1, data2, data3 ):
        """
        Plot background fit.
        """
        #Construct plot frame in 'x=mass_tag_measured'
        frame1 = lifetime_tag.frame(ROOT.RooFit.Bins(225))
        frame2 = lifetime_tag.frame(ROOT.RooFit.Bins(225))
 
        #print("Name of pdf: {}".format(total_pdf.GetName()))
        #Plot data in the frame, set histrogram properties
        nData=100
        if data1:
          #data1.plotOn(frame1, ROOT.RooFit.LineWidth(1), ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.MarkerColor(ROOT.kRed))
          data1.plotOn(frame1, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.LineWidth(1), ROOT.RooFit.MarkerSize(0.6), ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.MarkerColor(ROOT.kRed))#ROOT.DataError(ROOT.RooAbsData.SumW2)
          
          dummyHist1 = data1.createHistogram("d1{}{}".format(self.FitBin,name),lifetime_tag)
          nData = data1.sumEntries()
          dummyHist1C = dummyHist1.Clone("{}Clone".format(dummyHist1.GetName()))
          #dummyHist1C.LineColor(ROOT.kRed)
        if data2:
          data2.plotOn(frame1, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.LineWidth(1), ROOT.RooFit.MarkerSize(0.6), ROOT.RooFit.MarkerStyle(21),ROOT.RooFit.MarkerColor(ROOT.kBlue))
          dummyHist2 = data2.createHistogram("d2{}{}".format(self.FitBin,name),lifetime_tag)
          dummyHist2C = dummyHist2.Clone("{}Clone".format(dummyHist2.GetName()))
          #dummyHist2C.LineColor(ROOT.kBlue)
        if data3:
          data3.plotOn(frame1, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.LineWidth(1), ROOT.RooFit.MarkerSize(0.6), ROOT.RooFit.MarkerStyle(22),ROOT.RooFit.MarkerColor(ROOT.kGreen))
          dummyHist3 = data3.createHistogram("d3{}{}".format(self.FitBin,name),lifetime_tag)
          dummyHist3C = dummyHist3.Clone("{}Clone".format(dummyHist3.GetName()))
          #dummyHist3C.LineColor(ROOT.kGreen)


        myLegend1 = ROOT.TLegend(0.3,0.7,0.55,0.88)
        myLegend1.SetTextFont(42)
        myLegend1.SetBorderSize(0)
        myLegend1.SetFillStyle(0)
        myLegend1.SetFillColor(0)
        myLegend1.SetMargin(0.25)
        myLegend1.SetTextSize(0.04)
         
        myLegend2 = ROOT.TLegend(0.3,0.7,0.55,0.88)
        myLegend2.SetTextFont(42)
        myLegend2.SetBorderSize(0)
        myLegend2.SetFillStyle(0)
        myLegend2.SetFillColor(0)
        myLegend2.SetMargin(0.25)
        myLegend2.SetTextSize(0.04)
 
   
        #- - - - - - - - - - - - - -
        canvas = ROOT.TCanvas("canvasMulti{}".format(self.FitBin),"canvasMulti{}".format(self.FitBin), 800*2, 400*2)
        canvas.Divide(2)
        canvas.cd(1).SetLogy()

        myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
        myPad.SetLeftMargin(0.7)
        myPad.SetTopMargin(0.07)
        myPad.SetRightMargin(0.1)
        myPad.SetBottomMargin(0.15)
        myPad.Draw()#d not workoes
        
        max = frame1.GetMaximum()
        print("oo max in frame: {}".format(max))
        frame1.SetAxisRange(1, max*1.1, "Y");
        frame1.Draw("E")

        if data1:
          dummyHist1C.SetLineColor(ROOT.kRed)
          dummyHist1C.SetMarkerColor(ROOT.kRed)
          dummyHist1C.SetMarkerSize(0.7)
          dummyHist1C.SetMarkerStyle(20)
          myLegend1.AddEntry(dummyHist1C,"%s: %1.0f Evts"%(data1.GetName(),data1.sumEntries()),"PL")
        if data2:
          dummyHist2C.SetLineColor(ROOT.kBlue)
          dummyHist2C.SetMarkerColor(ROOT.kBlue)
          dummyHist2C.SetMarkerSize(0.7)
          dummyHist2C.SetMarkerStyle(21)
          myLegend1.AddEntry(dummyHist2C,"%s: %1.0f Evts"%(data2.GetName(),data2.sumEntries()),"PL")
        if data3:
          dummyHist3C.SetLineColor(ROOT.kGreen)
          dummyHist3C.SetMarkerColor(ROOT.kGreen)
          dummyHist3C.SetMarkerSize(0.7)
          dummyHist3C.SetMarkerStyle(22)
          myLegend1.AddEntry(dummyHist3C,"%s: %1.0f Evts"%(data3.GetName(),data3.sumEntries()),"PL")
 
        myLegend1.Draw()
  
        #- - - -
        if drawLogDiff:
          canvas.cd(2).SetLogy()
        else:
          canvas.cd(2)
        
        #dummyHist1.divide(dummyHist2)#AttributeError: 'TH1F' object has no attribute 'divide' . Bool_t Divide(const TH1* h1)
        #print("Types:")
        #print(type(dummyHist1C))
        #print(type(dummyHist2C))
        if drawLogDiff:
          yLegend="{}-{}".format(data1.GetName(),data2.GetName())
        else:
          nC1 = dummyHist1C.GetEntries()
          nC2 = dummyHist2C.GetEntries()
          if nC1>0:
            dummyHist1C.Scale(1./nC1)
          if nC2>0:
            dummyHist2C.Scale(1./nC2)
          yLegend="%s_{norm}-%s_{norm}" % (data1.GetName(),data2.GetName())

 
        self.setGraphHisto(dummyHist1C, "t [s]", yLegend)
        dummyHist1C.Add(dummyHist2C,-1)
        if drawLogDiff:
          dummyHist1C.GetYaxis().SetRangeUser(1,dummyHist1C.GetMaximum()*1.1)
          
        if data2:
          dummyHist2C.DrawCopy("E")
        if data3:
          dummyHist3C.DrawCopy("E")

        dummyHist1C.Draw("E")

        #yieldDiff=dummyHist1C.GetEntries()
        yieldDiff=dummyHist1C.Integral()
        myLegend2.AddEntry(dummyHist2C,"yield of diff: %1.0f" % (yieldDiff),"PL")
        myLegend2.Draw()

        canvas.SaveAs("{}TimeFits_{}/Mult_Bin{}_{}.pdf".format(self.basepath,self.obsSelection,name, self.FitBin))
        '''
        if dummyHist1:
          del dummyHist1
        if dummyHist2:
          del dummyHist2
        if dummyHist3:
          del dummyHist3
        '''

    def yieldGraphAbsZ(self):
        pass

    def yieldGraphAbsCosTheta(self):
        pass

    def yieldGraphNormZ(self):
        pass

    def yieldGraphNormCosTheta(self):
        pass

    def closureTest1D(self):
        pass

    def closureTest2D(self):
        pass

    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def setGraphHisto(self, Histo, Xtitel, Ytitel, border=0.1, logPlot=False):
  
      ROOT.TGaxis.SetMaxDigits(2)
      ROOT.gStyle.SetOptStat(0)
      #Set a larger y-axis range to leave space for the border
      min  = Histo.GetMinimum();
      max  = Histo.GetMaximum();
      #max  = Histo.GetBinContent(Histo.GetMaximumBin());
      range= max-min;
      #border=0.1
      if range>0:
        maxNew=min+range*(1.0+3*border);
        minNew=max-range*(1.0+border);
      else:
        maxNew=max+(-1)*range*(1.0+2*border);
        minNew=min-(-1)*range*(1.0+border);
      minNew=min;
      #minNew=0; #this is a problem for log scale
      if logPlot and minNew<1:
        minNew=1
      Histo.GetYaxis().SetRangeUser(minNew,maxNew);
      
        
      if Xtitel!="":
        Histo.GetXaxis().SetTitle(Xtitel);
      if Ytitel!="":
        Histo.GetYaxis().SetTitle(Ytitel);
      '''
      Histo.SetLineColor(1);
      Histo.SetMarkerColor(1);
      Histo.SetMarkerStyle(20);
      Histo.SetMarkerSize(0.7);
      '''
    
    #. . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def ipchi2FitPlot(self, resonance, log_ipchi2, data, total_pdf, nonprompt_pdf=None, prompt_pdf=None, 
                    background_pdf=None, prompt_yield=None, nonprompt_yield=None):
        """
        Create plot for IP chi2 fit using single Bukin functions for prompt and nonprompt components.
        """
        # Extract parameters from the fit
        RooArgSet = total_pdf.getParameters(data)
        prompt_frac = RooArgSet.getRealValue("prompt_frac", -1)
        
        # Parameters for Bukin functions
        # Prompt component
        xp_prompt = RooArgSet.getRealValue("xp_prompt", -1)
        sigma_prompt = RooArgSet.getRealValue("sigma_prompt", -1)
        xi_prompt = RooArgSet.getRealValue("xi_prompt", -1)
        rho1_prompt = RooArgSet.getRealValue("rho1_prompt", -1)
        rho2_prompt = RooArgSet.getRealValue("rho2_prompt", -1)
        
        # Non-prompt component
        xp_nonprompt = RooArgSet.getRealValue("xp_nonprompt", -1)
        sigma_nonprompt = RooArgSet.getRealValue("sigma_nonprompt", -1)
        xi_nonprompt = RooArgSet.getRealValue("xi_nonprompt", -1)
        rho1_nonprompt = RooArgSet.getRealValue("rho1_nonprompt", -1)
        rho2_nonprompt = RooArgSet.getRealValue("rho2_nonprompt", -1)
        
        # Get total signal yield
        if prompt_yield and nonprompt_yield:
            prompt_val = prompt_yield.getVal()
            nonprompt_val = nonprompt_yield.getVal()
            sig_yield = prompt_val + nonprompt_val
        else:
            sig_yield = RooArgSet.getRealValue("sig_yieldLim", -1)
            prompt_val = sig_yield * prompt_frac
            nonprompt_val = sig_yield * (1 - prompt_frac)
        
        # Get background yield if available
        if background_pdf:
            bkg_yield = RooArgSet.getRealValue("bkg_yieldLim", -1)
        
        # Create frame for plotting
        frame1 = log_ipchi2.frame(ROOT.RooFit.Bins(50))
        frame1.SetTitle(f"log(IP Chi2) Distribution for {self.range}")
        
        # Plot data
        data.plotOn(frame1, ROOT.RooFit.Name("datahistogram"), 
                  ROOT.RooFit.LineColor(ROOT.kGray+2), ROOT.RooFit.LineWidth(1), 
                  ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.MarkerStyle(20))
        
        # Plot models
        total_pdf.plotOn(frame1, ROOT.RooFit.LineColor(ROOT.kGreen+1), 
                       ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(2))
        
        # Plot prompt component
        total_pdf.plotOn(frame1, ROOT.RooFit.Components("prompt_pdf"), 
                       ROOT.RooFit.LineColor(ROOT.kBlue), 
                       ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(2))
        
        # Plot non-prompt component
        total_pdf.plotOn(frame1, ROOT.RooFit.Components("nonprompt_pdf"), 
                       ROOT.RooFit.LineColor(ROOT.kRed), 
                       ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineWidth(2))
        
        # Plot background if available
        if background_pdf:
            total_pdf.plotOn(frame1, ROOT.RooFit.Components("bkg_pdf"), 
                           ROOT.RooFit.LineColor(ROOT.kBlack), 
                           ROOT.RooFit.LineStyle(4), ROOT.RooFit.LineWidth(2))
        
        # Create legend with fit parameters
        myLegend1 = ROOT.TLegend(0.55, 0.50, 0.75, 0.88)
        myLegend1.SetTextFont(42)
        myLegend1.SetBorderSize(0)
        myLegend1.SetFillStyle(0)
        myLegend1.SetFillColor(0)
        myLegend1.SetMargin(0.25)
        myLegend1.SetTextSize(0.030)
        
        # Add Bukin function parameters to legend
        myLegend1.AddEntry(frame1, "Prompt Bukin:", "")
        myLegend1.AddEntry(frame1, f" Peak: {xp_prompt:.3f}", "")
        myLegend1.AddEntry(frame1, f" Width: {sigma_prompt:.3f}", "")
        myLegend1.AddEntry(frame1, f" Asym: {xi_prompt:.3f}", "")
        myLegend1.AddEntry(frame1, f" Tails: {rho1_prompt:.3f}, {rho2_prompt:.3f}", "")
        
        myLegend1.AddEntry(frame1, "Non-prompt Bukin:", "")
        myLegend1.AddEntry(frame1, f" Peak: {xp_nonprompt:.3f}", "")
        myLegend1.AddEntry(frame1, f" Width: {sigma_nonprompt:.3f}", "")
        myLegend1.AddEntry(frame1, f" Asym: {xi_nonprompt:.3f}", "")
        myLegend1.AddEntry(frame1, f" Tails: {rho1_nonprompt:.3f}, {rho2_nonprompt:.3f}", "")
        
        myLegend1.AddEntry(frame1, f"Prompt frac: {prompt_frac*100:.1f}%", "")
        
        # Create yield legend
        myLegend3 = ROOT.TLegend(0.25, 0.81, 0.45, 0.91)
        myLegend3.SetTextFont(42)
        myLegend3.SetBorderSize(0)
        myLegend3.SetFillStyle(0)
        myLegend3.SetFillColor(0)
        myLegend3.SetMargin(0.25)
        myLegend3.SetTextSize(0.03)
        
        myLegend3.AddEntry(frame1, f"Prompt yield: {prompt_val:.2f}", "")
        myLegend3.AddEntry(frame1, f"Non-prompt yield: {nonprompt_val:.2f}", "")
        if background_pdf:
            myLegend3.AddEntry(frame1, f"BKG yield: {bkg_yield:.2f}", "")
        
        # Create pull distribution
        dataHist = frame1.getHist("datahistogram")
        curve1 = frame1.getObject(1)
        hpull = dataHist.makePullHist(curve1, True)
        
        frame2 = log_ipchi2.frame(ROOT.RooFit.Title("Pull Distribution"), ROOT.RooFit.Bins(50))
        frame2.addPlotable(hpull, "P")
        frame2.getAttMarker().SetMarkerSize(0.5)
        frame2.getAttLine().SetLineWidth(1)
        
        # Draw fit and pull distribution
        canvasFit = ROOT.TCanvas("canvasFit", "canvasFit", 800*2, 400*2)
        canvasFit.Divide(2)
        
        # First pad with log scale for full range
        canvasFit.cd(1)
        canvasFit.GetPad(1).SetLogy()  # Set y-axis to logarithmic scale
        
        # Set appropriate y-axis range for log scale
        max_val = frame1.GetMaximum()
        frame1.SetAxisRange(0.9, max_val*5, "Y")  # Min value > 0 for log scale
        
        frame1.Draw()
        myLegend1.Draw()
        myLegend3.Draw()
        
        # Second pad for zoomed region to see asymmetry better
        canvasFit.cd(2)
        # Clone frame for zoomed view
        frame_zoom = log_ipchi2.frame(ROOT.RooFit.Title("Zoomed View"), ROOT.RooFit.Bins(50))
        data.plotOn(frame_zoom, ROOT.RooFit.Name("datahistogram_zoom"), 
                  ROOT.RooFit.LineColor(ROOT.kGray+2), ROOT.RooFit.LineWidth(1), 
                  ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.MarkerStyle(20))
        
        total_pdf.plotOn(frame_zoom, ROOT.RooFit.LineColor(ROOT.kGreen+1), 
                       ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(2))
        
        total_pdf.plotOn(frame_zoom, ROOT.RooFit.Components("prompt_pdf"), 
                       ROOT.RooFit.LineColor(ROOT.kBlue), 
                       ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(2))
        
        total_pdf.plotOn(frame_zoom, ROOT.RooFit.Components("nonprompt_pdf"), 
                       ROOT.RooFit.LineColor(ROOT.kRed), 
                       ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineWidth(2))
        
        # Focus on central region to see asymmetry
        frame_zoom.SetAxisRange(-1.5, 4, "X")
        frame_zoom.Draw()
        
        # Save plot
        output_dir = f"{self.basepath}IPChi2Fits_{self.obsSelection}/"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            
        canvasFit.SaveAs(f"{output_dir}IPChi2Fit_SingleBukin_{self.FitBin}.{self.format}")
        
        # Create pull canvas
        canvasPull = ROOT.TCanvas("canvasPull", "canvasPull", 800, 400)
        canvasPull.cd()
        ROOT.gPad.SetLeftMargin(0.15)
        frame2.GetYaxis().SetTitleOffset(1.6)
        frame2.Draw()
        canvasPull.SaveAs(f"{output_dir}IPChi2Fit_SingleBukin_Pull_{self.FitBin}.{self.format}")
        
        # Create histogram for return
        h_data = data.createHistogram(f"h_IPChi2_SingleBukin_{self.FitBin}_data", log_ipchi2)
        h_data.GetXaxis().SetRangeUser(*log_ipchi2.getRange())
        
        return h_data

