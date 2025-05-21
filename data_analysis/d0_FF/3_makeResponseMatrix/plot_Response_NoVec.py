#! /usr/bin/env python
#
# Plot Response matric from LHCb analyses
# python plotResponse.py
#lb-run DaVinci/v46r0 python plot_Response_NoVec.py

# General
import os
import sys
import argparse
import itertools
import math
from array import *
import numpy
# ROOT
import ROOT
from ROOT import gROOT
from ROOT import gStyle
from ROOT import gDirectory
from ROOT import Experimental
# from ROOT import PyROOT
from ROOT import TGraph
from math import sqrt, pow, fabs, pi, atan2, log, exp
import time
from time import process_time

import gc         #garbage collector to free memory

#from ROOT import gExperimental
ROOT.gROOT.SetBatch(True)

class RMplotter:
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__(self):
        #general variables
        self.startTime = time.process_time()
        self.currentTime = 0
        #fileFormat = "C"
        self.fileFormat = "pdf"
        #isMC=False
        self.isMC = True
        self.ColorArray = [ROOT.kViolet+5,ROOT.kAzure,ROOT.kAzure-4,ROOT.kCyan-6,ROOT.kGreen-3,ROOT.kTeal-6]
        #self.dictKey = "All"
        self.dictKey = "Psi2S"
        #self.dictKey = "X3872"
        
        Diagonal_line = ROOT.TLine();
        Diagonal_line.SetX1(0);
        Diagonal_line.SetX2(1);
        Diagonal_line.SetY1(0);
        Diagonal_line.SetY2(1);
        Diagonal_line.SetLineWidth(4);
        Diagonal_line.SetLineStyle(2);
        Diagonal_line.SetLineColor(ROOT.kGray);
        self.diagonal_line=Diagonal_line
 
        Diagonal_lineD = ROOT.TLine();
        Diagonal_lineD.SetX1(0);
        Diagonal_lineD.SetX2(1);
        Diagonal_lineD.SetY1(0);
        Diagonal_lineD.SetY2(1);
        Diagonal_lineD.SetLineWidth(4);
        Diagonal_lineD.SetLineStyle(2);
        Diagonal_lineD.SetLineColor(ROOT.kGray+2);
        self.diagonal_lineDark=Diagonal_lineD

        # if "X3872" in self.dictKey or "All" in self.dictKey :
        #self.pTBinning  = [5,10,15,20,30] #  new-future
        self.pTBinning  = [0, 5,10,15,20,30,200] #  new-future
        self.zTBinArray = [0.2,0.5,0.65,0.75,0.85,0.95,1]
        #self.zTBinArray = [0,0.5,1,1.1]
        # else:
        #   #self.pTBinning  = [0,5,10,15,20,30,40,60,200]#new-future
        #   self.pTBinning  = [0,5,10,15,20,30,40,60,200]#new-future
        #   self.zTBinArray = [0,0.3,0.4,0.5,0.6,0.75,0.84,0.9,0.94,0.97,1]
 
        input, output = self.setupInputOutput()
        
        self.fInFile    = ROOT.TFile(input)
        self.inTreeRM   = self.fInFile.Get("Response")

        self.dataSet = self.createDataSet()
        self.output  = output
        
    #------------------------------------------------------------
    #------------------------------------------------------------
    def setupInputOutput(self):
     
      inputDir="/media/niviths/local/analysis_code/data_analysis/d0_FF/3_makeResponseMatrix"

      #-MC samples
      #fFileName  = "{}/allMedley277-288.root".format(inputDir) #Bdecay
      #fFileName  = "{}/testRoot.root".format(inputDir)   # #prompt X(3872) in 2016 and 2018
      #fFileName  = "{}/allMedley289-294.root".format(inputDir)#Dir
      #fFileName  = "{}/job277_AllFilteredNV.root".format(inputDir)
      #fFileName  = "{}/AllBdecay.root".format(inputDir)
      #-Latest Data Apr-July
      fFileName  = "/media/niviths/SSD2/lhcb_analysis_SSD/20250514_Pbp_21_MC_output_D0FF_filterV1.root"
      #fFileName  = "{}/16_17_18AllBdecay.root".format(inputDir)
      #-MC with new filter. Only fiducial volume cuts applied for MC
      #fFileName  = "{}/Prompt_Aug22/all_X3872_2016_2018.root".format(inputDir)

 
      print("Found file: {}".format(fFileName))
      
      #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
      #outputDirFig = "./FiguresRM_NP_TT/"
      outputDirFig = "{}/D0_FF_ResMatr/".format(inputDir)

      if not os.path.exists(outputDirFig):
        os.makedirs(outputDirFig)
        print("make new directory: {}".format(outputDirFig))
     
      return fFileName, outputDirFig
     
    #------------------------------------------------------------
    #------------------------------------------------------------
    def createDataSet(self):
         

      #. . Dataset Parameters. .
      #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pt_jetDet    = ROOT.RooRealVar("jetPtDet","jetPtDet",0,500)
      pt_tagDet    = ROOT.RooRealVar("tagPtDet","tagPtDet",0,80)
      etaDet       = ROOT.RooRealVar("etaDet","etaDet",1.5,5.5)
      phiDet       = ROOT.RooRealVar("phiDet","phiDet",-3.5,3.5)
      tagMDet      = ROOT.RooRealVar("tagMDet","tagMDet",1.81, 1.935)
      zTDet        = ROOT.RooRealVar("zTDet","zTDet",0,1.01)
      pt_jetPart   = ROOT.RooRealVar("jetPtPart","jetPtPart",0,500)
      pt_tagPart   = ROOT.RooRealVar("tagPtPart","tagPtPart",0,80)
      etaPart      = ROOT.RooRealVar("etaPart","etaPart",1.5,5.5)
      phiPart      = ROOT.RooRealVar("phiPart","phiPart",-3.5,3.5)
      tagMPart     = ROOT.RooRealVar("tagMPart","tagMPart",1.81, 1.935)
      zTPart       = ROOT.RooRealVar("zTPart","zTPart",0,1.01)
      nConstDet    = ROOT.RooRealVar("nConstDet","nConstDet",0,50)
      nConstPart   = ROOT.RooRealVar("nConstPart","nConstPart",0,50)

      dR           = ROOT.RooRealVar("dR","dR",-1,1)

      cutVars = ROOT.RooArgSet()
      cutVars.add(pt_jetDet)
      cutVars.add(pt_tagDet)
      cutVars.add(tagMDet)
      cutVars.add(pt_jetPart)
      cutVars.add(pt_tagPart)
      cutVars.add(tagMPart)
      cutVars.add(dR)
      cutVars.add(zTDet)
      cutVars.add(zTPart)
      cutVars.add(nConstDet)
      cutVars.add(nConstPart)
  
      #cutString = "dR > -1 && MCassocType==1"
      #cutString = "dR > -1 && MCassocType==1 && nConstDet>1"
        #Cuts for RM quality
      
      MassCut             = (1.81, 1.935)
      cutString   = "tagMPart > {} && tagMPart < {} ".format(MassCut[0],MassCut[1])
      cutString += " && dR > -1 && nConstDet>1"#  && etaDet>2.5 && etaDet<4"

      data = ROOT.RooDataSet("Test1Corr", "Test1Corr", self.inTreeRM, cutVars, cutString)
      
      return data

    #------------------------------------------------------------
    #------------------------------------------------------------
    def reduceDataSet(self,jetPtMin,jetPtMax,dRCut=1):
      
      #. . Dataset Cuts. .
      #- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      # Cuts for the binning
      basicCutString = "&& jetPtDet > {} && jetPtDet < {} ".format(jetPtMin,jetPtMax)
      dRString       = "&& dR < {}".format(dRCut)

      #-necesarry cuts for each round
      if dRCut<1:
        basicCutString+=dRString
      
      #dataReduced = self.dataSet.reduce(basicCutString)
      dataReduced = self.dataSet
      
      print("------------------------")
      #check how much statistic is lost due to the reduction
      tagMDet      = ROOT.RooRealVar("tagMDet","tagMDet",1.81, 1.935)
      tagMPart     = ROOT.RooRealVar("tagMPart","tagMPart",1.81, 1.935)
      
      self.plotQAPlots(dataReduced,tagMDet,tagMPart,"Mass")

      return dataReduced
      

    #------------------------------------------------------------
    #------------------------------------------------------------
    def plotQAPlots(self, dataReduced,tagMDet,tagMPart,namekey=""):
        
        print("**plot dataset**")
        h_dataFullDet  = self.dataSet.createHistogram("h_MassSignalDetFull{}_data".format(bin),tagMDet)
        h_dataRedDet   = dataReduced.createHistogram("h_MassSignalDetRed{}_data".format(bin),tagMDet)
    
        h_dataFullPart  = self.dataSet.createHistogram("h_MassSignalPartFull{}_data".format(bin),tagMPart)
        h_dataRedPart   = dataReduced.createHistogram("h_MassSignalDetRed{}_data".format(bin),tagMPart)
    
        canvas = ROOT.TCanvas("canvas","canvas", 800*2, 400*2)
        canvas.Divide(2)
        canvas.cd(1)
        myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
        myPad.SetLeftMargin(0.7)
        myPad.SetTopMargin(0.07)
        myPad.SetRightMargin(0.1)
        myPad.SetBottomMargin(0.15)
        myPad.Draw()#d not workoes
        
        h_dataFullDet.SetLineColor(1)
        h_dataFullDet.Draw("hist")
        h_dataRedDet.SetLineColor(ROOT.kBlue)
        h_dataRedDet.DrawCopy("same hist")
     
        canvas.cd(2)
        myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
        myPad.SetLeftMargin(0.7)
        myPad.SetTopMargin(0.07)
        myPad.SetRightMargin(0.1)
        myPad.SetBottomMargin(0.15)
        myPad.Draw()#d not workoes
        
        h_dataFullPart.SetLineColor(1)
        h_dataFullPart.Draw("hist")
        h_dataRedPart.SetLineColor(ROOT.kBlue)
        h_dataRedPart.DrawCopy("same hist")

        myLegend2 = ROOT.TLegend(0.2,0.8,0.4,0.9)
        myLegend2.SetTextFont(42)
        myLegend2.SetBorderSize(0)
        myLegend2.SetFillStyle(0)
        myLegend2.SetFillColor(0)
        myLegend2.SetMargin(0.25)
        myLegend2.SetTextSize(0.04)
        myLegend2.AddEntry(h_dataFullPart,"Data Full","L")
        myLegend2.AddEntry(h_dataRedPart,"Data Reduced","L")
        myLegend2.Draw()
        canvas.SaveAs("{}{}QA{}.{}".format(self.output,self.dictKey,namekey,self.fileFormat))
 
    #------------------------------------------------------------
    #------------------------------------------------------------
    def plotRMJet(self,ptLow,ptHigh,dR=1):
      
      reducedData = self.reduceDataSet(ptLow,ptHigh)
      #reducedData = self.dataSet
    
      pt_jetDet    = ROOT.RooRealVar("jetPtDet","jetPtDet",0,100)
      pt_jetPart   = ROOT.RooRealVar("jetPtPart","jetPtPart",0,100)

      hh_data = reducedData.createHistogram("hist",pt_jetDet,ROOT.RooFit.Binning(50,0,100),ROOT.RooFit.YVar(pt_jetPart,ROOT.RooFit.Binning(50,0,100)))

      # - - - - - - - - - - - - - - - - - - - - - - -
      # Draw histogram on canvas
      # - - - - - - - - - - - - - - - - - - - - - - -
      canvas = ROOT.TCanvas("canvas","canvas", 800*2, 400*2)
      #canvas.Divide(2)
      canvas.cd(1)
      myPad = ROOT.TPad("myPad", "The pad",0,0,1,1)
      myPad.SetLeftMargin(0.7)
      myPad.SetTopMargin(0.07)
      myPad.SetRightMargin(0.1)
      myPad.SetBottomMargin(0.15)
      myPad.Draw()
      hh_data.Draw("colz")
        
      ptTag = "pT{}-{}".format(ptLow,ptHigh)
      canvas.SaveAs("{}{}Bin{}.{}".format(self.output,self.dictKey,ptTag,self.fileFormat))
  
    #------------------------------------------------------------
    #------------------------------------------------------------
    def plotRMTag(self,dR=1):
      
      reducedData = self.reduceDataSet(0,200) #Jet pT range for entire jet sample for unfoldins

      #-Variables of the RM
      pt_jetDet    = ROOT.RooRealVar("jetPtDet","jetPtDet",0,200)
      pt_jetPart   = ROOT.RooRealVar("jetPtPart","jetPtPart",0,200)
      zTDet  = ROOT.RooRealVar("zTDet","zTDet",0,1.01)
      zTPart = ROOT.RooRealVar("zTPart","zTPart",0,1.01)
             
      #-Create custom binning for the plot
      arrayDet  = numpy.array(self.zTBinArray)
      tbinsDet  = ROOT.RooBinning(len(self.zTBinArray)-1, arrayDet ,"zTBinningDet")
      arrayPart = numpy.insert(arrayDet,0,0)
      tbinsPart = ROOT.RooBinning(len(self.zTBinArray), arrayPart ,"zTBinningPart")
      
      nBinsGenLvlv    = len(self.pTBinning)-1
      nBinsDetLvlv    = len(self.pTBinning)-1
      maxCombinations = nBinsGenLvlv*nBinsDetLvlv

      hh_data     = numpy.empty((nBinsDetLvlv,nBinsGenLvlv),dtype=ROOT.TH2F)
      maxArray    = [None]*maxCombinations
      
      legGen      = [None]*nBinsGenLvlv
      legDet      = [None]*nBinsDetLvlv
      systemArray = numpy.empty((nBinsDetLvlv,nBinsGenLvlv),dtype=ROOT.TLatex)
      
      print("Max No combinations: {}".format(maxCombinations))
      print("pT Part bins: {}".format(len(self.pTBinning)-1))
      print("pT Det bins: {}".format(len(self.pTBinning)-1))
      
      # - - - - - - - - - - - - - - - - - - - - - - -
      # create all histograms for the zt responses
      # - - - - - - - - - - - - - - - - - - - - - - -
      histIndex=0
      #for indexGenerator in range(0,len(self.pTBinning)):
      for indexGenerator in range(0,len(self.pTBinning)-1):
          ptRangePart = "jetPtPart > {} && jetPtPart < {}".format(self.pTBinning[indexGenerator],self.pTBinning[indexGenerator+1])

          for indexDet in range(0,len(self.pTBinning)-1):
              ptRangeDet = "&& jetPtDet > {} && jetPtDet < {}".format(self.pTBinning[indexDet],self.pTBinning[indexDet+1])
              cuts       = ptRangePart+ptRangeDet

              #-build legend
              legGen[indexGenerator] = "p_{T}^{jet,gen}: %1.0f-%1.0f" % (self.pTBinning[indexGenerator],self.pTBinning[indexGenerator+1])
              legDet[indexDet]       = "p_{T}^{jet,det}: %1.0f-%1.0f" % (self.pTBinning[indexDet],self.pTBinning[indexDet+1])

              
              dataInBin = reducedData.reduce(ROOT.RooArgSet(zTDet,zTPart,pt_jetDet,pt_jetPart),cuts)
              #dataInBin = self.dataSet
              #-Build the histogram in the specific cur ranges for Det and Gen pt
              hh_data[indexDet,indexGenerator] = dataInBin.createHistogram("hist{}".format(histIndex),zTDet,ROOT.RooFit.Binning(tbinsDet),ROOT.RooFit.YVar(zTPart,ROOT.RooFit.Binning(tbinsPart)))
              maxArray[histIndex] = hh_data[indexDet,indexGenerator].GetMaximum()
              print(hh_data[indexDet,indexGenerator].Integral())
              histIndex+=1
              #QA
              #self.plotQAPlots(dataInBin,zTDet,zTPart,"zT")

      # - - - - - - - - - - - - - - - - - - - - - - -
      # Plot all zT responses in one plot
      # - - - - - - - - - - - - - - - - - - - - - - -
      #-Determine a global range to plot all panels in the same z-axis scale range
      #-for better comparison
      globalMax = max(maxArray)
      #globalMax = 1000
      try:
          # Check if any positive values exist in the array
          positive_values = [i for i in maxArray if i > 0]
          if positive_values:
              globalMin = min(positive_values)  # find minimum actual value
          else:
              print("Warning: No positive values found in response matrix. Using default minimum.")
              globalMin = 1e-6  # Set a reasonable default minimum
      except Exception as e:
          print(f"Error finding minimum value: {e}")
          globalMin = 1e-6  # Set a default minimum
    
      #Nx = len(self.pTBinning)-1
      #Ny = len(self.pTBinning)-1
      lMargin=0.08 #was 0.05
      rMargin=0.03 #was 0.03
      bMargin=0.07 #was 0.05
      tMargin=0.03
      
      ROOT.gStyle.SetOptStat(0);

      canvas3 = ROOT.TCanvas("C3","C3",1600,1600)
      canvas3.SetFillStyle(4000)
      pad = self.partitionCanvas(canvas3, nBinsDetLvlv, nBinsGenLvlv, lMargin, rMargin, bMargin, tMargin)
  
      print("length pads: {}".format(len(pad)))
      print("length histos: {}".format(len(hh_data)))
      #print("shape pads: {}".format(pad.shape))
      #print("shape histos: {}".format(hh_data.shape))
      panelNumber=0
      for i in range(0,nBinsDetLvlv):   #Jet pT Detector
        for j in range(0,nBinsGenLvlv): #Jet pT Generator
            canvas3.cd(0);
            
            #-Get the pads previously created.
            pname="pad_%i_%i"%(i,j)
            print(pname)
            if not pad[i][j]: continue
            pad[i][j] = ROOT.gROOT.FindObject(pname)
            pad[i][j].Draw()
            pad[i][j].SetFillStyle(4000)
            pad[i][j].SetFrameFillStyle(4000)
            pad[i][j].SetLogz()
            pad[i][j].cd()
         
            xFactor = pad[0][0].GetAbsWNDC()/pad[i][j].GetAbsWNDC();
            yFactor = pad[0][0].GetAbsHNDC()/pad[i][j].GetAbsHNDC();
  
            #-regarding nBinsGenLvlv-j-1: Invert the y axis. On the plot it increases from bottom to top
            if hh_data[i,nBinsGenLvlv-j-1]:
              print("NAme: {}".format(hh_data[i,nBinsGenLvlv-j-1].GetName()))
              hh_data[i,nBinsGenLvlv-j-1].SetTitle("")
              hFrame = self.setFrame(hh_data[i,nBinsGenLvlv-j-1], xFactor, yFactor, i, nBinsDetLvlv, j, nBinsGenLvlv)#, "p_{T}^{jet,Det}", "p_{T}^{jet,Gen}")
              hFrame.GetZaxis().SetRangeUser(globalMin,globalMax)
              hFrame.GetXaxis().SetTitle("z_{T}^{Det. lvl}")
              hFrame.GetYaxis().SetTitle("z_{T}^{Gen. lvl}")
              hFrame.DrawCopy();

              hh_data[i,nBinsGenLvlv-j-1].DrawCopy("same colz")
            
              topM  = pad[i][j].GetTopMargin()
              LeftM = pad[i][j].GetLeftMargin()
              print("top margin: {}".format(topM))

              print("x/y*0.8: {}".format(0.85-topM))
              print("y/x*1.16: {}".format(0.14-LeftM))
              #systemArray[i,Ny-j-1] = ROOT.TLatex(0.16-LeftM,0.8-topM,"#splitline{%s}{%s}" % (legGen[j],legDet[i]))
              #systemArray[i,Ny-j-1].SetNDC()
              systemArray[i,nBinsGenLvlv-j-1] = ROOT.TLatex(0.16,0.8,"#splitline{%s}{%s}" % (legGen[nBinsGenLvlv-j-1],legDet[i]))
              systemArray[i,nBinsGenLvlv-j-1].SetTextSize(0.09)
              systemArray[i,nBinsGenLvlv-j-1].Draw()
                
              if i==nBinsGenLvlv-j-1:
                self.diagonal_lineDark.Draw()
              else:
                self.diagonal_line.Draw()
            else:
              print("Error: could not find histogram at position: {}".format(panelNumber))
            panelNumber+=1
      canvas3.cd();
      ptTag = "pT{}-{}".format(5,70)
      canvas3.SaveAs("{}{}zTBin{}.{}".format(self.output,self.dictKey,ptTag,self.fileFormat))
  
  
  
    #________________________________________________________________________
    # Funtion to prepare a merged zT response canvas
    #________________________________________________________________________
    def partitionCanvas(self, C, Nx, Ny, lMargin, rMargin, bMargin, tMargin):
    
      if not C:
        return;
      #Setup Pad layout:
      #vertical space between pads
      vSpacing = 0.0;
      #vertical size of a pad
      vStep    = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
         
      #horizontal space between pads
      hSpacing = 0.0;
      #horizontal size of a pad
      hStep    = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

      pad = [[None]*Ny]*Nx
      ###print(pad)
      #colum
      for i in range(0,Nx):
          if i==0:
            hposl = 0.0;
            hposr = lMargin + hStep;
            hfactor = hposr-hposl;
            hmarl = lMargin / hfactor;
            hmarr = 0.0;
          elif i == Nx-1:
            hposl = hposr + hSpacing;
            hposr = hposl + hStep + rMargin;
            hfactor = hposr-hposl;
            hmarl = 0.0;
            hmarr = rMargin / hfactor;
          else:
            hposl = hposr + hSpacing;
            hposr = hposl + hStep;
            hmarl = 0.0;
            hmarr = 0.0;
          #row
          for j in range(0,Ny):
            ###print(j)
            if j == 0:
              vposu = 1
              vposd = 1 - vStep - tMargin;
              vfactor = vposu-vposd;
              vmard   = 0.0;
              vmaru   = tMargin / vfactor;
            elif j == Ny-1:
              vposu = vposd - vSpacing;
              vposd = vposu - vStep - bMargin;
              vfactor = vposu-vposd;
              vmard   = bMargin / vfactor;
              vmaru   = 0.0;
            else:
              vposu = vposd - vSpacing;
              vposd = vposu - vStep;
              vmard = 0.0;
              vmaru = 0.0;
           
            if vposd<0: vposd=0
            C.cd(0);
            name = "pad_%i_%i"%(i,j)
            print("print canvas name: {}".format(name))
            print("left: {}, right: {}, top: {}, bottom: {}".format(hposl,hposr,vposu,vposd))
            pad[i][j] = ROOT.TPad(name,"",hposl,vposd,hposr,vposu);
            #print("create pad: "+name)
            pad[i][j].SetLeftMargin(hmarl);
            pad[i][j].SetRightMargin(hmarr);
            pad[i][j].SetBottomMargin(vmard);
            pad[i][j].SetTopMargin(vmaru);
            pad[i][j].SetFrameBorderMode(0);
            pad[i][j].SetBorderMode(0);
            pad[i][j].SetBorderSize(0);
            pad[i][j].Draw();
      
     
      return pad

       
    #________________________________________________________________________
    # Funtion to set TH histograms to a similar style
    #________________________________________________________________________
    def setFrame(self,h , xFactor, yFactor, i, Nx, j, Ny):
        
        fname  = "frame_{}".format(h.GetName())
        hFrame = h.Clone(fname);
        
        #empty histogram
        hFrame.Reset();
        #y axis range
          
        hFrame.GetXaxis().SetRangeUser(0,1)
        hFrame.GetYaxis().SetRangeUser(0,1)

        #Format for y axis
        hFrame.GetYaxis().SetLabelFont(43);
        hFrame.GetYaxis().SetLabelSize(32);#was 16
        hFrame.GetYaxis().SetLabelOffset(0.01);#was 0.2
        hFrame.GetYaxis().SetTitleFont(43);
        hFrame.GetYaxis().SetTitleSize(35);#was 16
        hFrame.GetYaxis().SetTitleOffset(5);#was 5
        hFrame.GetYaxis().CenterTitle();
        hFrame.GetYaxis().SetNdivisions(505);
        #TICKS Y Axis
        hFrame.GetYaxis().SetTickLength(xFactor*0.04/yFactor);
        #Format for x axis
        hFrame.GetXaxis().SetLabelFont(43);
        hFrame.GetXaxis().SetLabelSize(32);#was16
        hFrame.GetXaxis().SetLabelOffset(0.01);
        hFrame.GetXaxis().SetTitleFont(43);
        hFrame.GetXaxis().SetTitleSize(38); #was 16
        hFrame.GetXaxis().SetTitleOffset(5);#was 5
        hFrame.GetXaxis().CenterTitle();
        hFrame.GetXaxis().SetNdivisions(505);
        #TICKS X Axis
        hFrame.GetXaxis().SetTickLength(yFactor*0.06/xFactor);

        #test
        #hFrame.GetZaxis().SetTickLength(xFactor*0.04/yFactor);
        #hFrame.GetZaxis().SetTickLength(0);
        hFrame.GetZaxis().SetLabelFont(43);
        hFrame.GetZaxis().SetLabelSize(32);
        hFrame.GetZaxis().SetTitle("");

        return hFrame
    #________________________________________________________________________
    # Funtion to set TH histograms to a similar style
    #________________________________________________________________________
    def setHisto(self,Histo, Xtitel, Ytitel, big=0, border=0.1):

      Histo.SetStats(0);
      Histo.SetTitle("");
       
      if big==0:
        Histo.GetYaxis().SetTitleOffset(1.4);
        Histo.GetXaxis().SetTitleOffset(1.4);
        Histo.GetXaxis().SetLabelSize(0.05);
        Histo.GetYaxis().SetLabelSize(0.05);
        Histo.GetXaxis().SetTitleSize(0.045);
        Histo.GetYaxis().SetTitleSize(0.045);
      if big==1:
        Histo.GetYaxis().SetTitleOffset(1.0);
        Histo.GetXaxis().SetTitleOffset(0.82);
        Histo.GetYaxis().SetLabelSize(0.06);
        Histo.GetXaxis().SetLabelSize(0.06);
        Histo.GetYaxis().SetTitleSize(0.07);
        Histo.GetXaxis().SetTitleSize(0.07);
      if big==2:
        Histo.GetYaxis().SetTitleOffset(1.2);
        Histo.GetXaxis().SetTitleOffset(1.2);
        Histo.GetXaxis().SetLabelSize(0.05);
        Histo.GetYaxis().SetLabelSize(0.05);
        Histo.GetXaxis().SetTitleSize(0.055);
        Histo.GetYaxis().SetTitleSize(0.055);

      Histo.GetXaxis().SetNdivisions(505);
      Histo.GetYaxis().SetNdivisions(505);
      #..make nice font
      Histo.GetXaxis().SetLabelFont(42);
      Histo.GetYaxis().SetLabelFont(42);
      Histo.GetXaxis().SetTitleFont(42);
      Histo.GetYaxis().SetTitleFont(42);
      if Xtitel!="":
        Histo.GetXaxis().SetTitle(Xtitel);
      if Ytitel!="":
        Histo.GetYaxis().SetTitle(Ytitel);

      #Check whether this is TH1 or TH2
      if Histo.InheritsFrom(ROOT.TH2.Class()):
        return
      else:
        self.setTH1Histo(Histo, Xtitel, Ytitel, big, border)

    #________________________________________________________________________
    # special part for TH1 histograms
    #________________________________________________________________________
    def setTH1Histo(self,Histo, Xtitel, Ytitel, big, border=0.1):

      #Set a larger y-axis range to leave space for the border
      min  = Histo.GetMinimum(0);
      max  = Histo.GetBinContent(Histo.GetMaximumBin());
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
      Histo.GetYaxis().SetRangeUser(minNew,maxNew);

      Histo.SetLineColor(1);
      Histo.SetMarkerColor(1);
      Histo.SetMarkerStyle(20);
      Histo.SetMarkerSize(0.7);


#________________________________________________________________________
#________________________________________________________________________
def plotResponse():
  
  plotterProcess = RMplotter()
  '''
  #plot different jet pT ranges
  plotterProcess.plotRMJet(5,10)
  plotterProcess.plotRMJet(10,15)
  plotterProcess.plotRMJet(15,20)
  plotterProcess.plotRMJet(20,30)
  plotterProcess.plotRMJet(30,40)
  plotterProcess.plotRMJet(40,70)
  plotterProcess.plotRMJet(5,70)
  '''
  
  plotterProcess.plotRMTag()
  
 

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  # Define arguments
  #parser = argparse.ArgumentParser(description="- -")

  # Parse the arguments
  #args = parser.parse_args()
 plotResponse()
