import ROOT
import pickle
import Plotter as plt
import math
import array
#import parameters as pars
import sys;
#import yaml

class Fitter:
    def __init__(self, tFile=None, resonance="", pvFile=None, rootFileName=None, isBinned=True, isMC=False, nBins=1, outfilePath=".", update=False):
        # Initialize dictionary
        self.initDictionary()
        self.TestFilename = rootFileName
        self.outfilePath  = outfilePath
        self.isBinned     = isBinned
        self.isMC         = isMC
        self.resonance    = resonance
        print("This is the MC status: {}".format(self.isMC))
        self.updateStartValues = update #this will specify if the start parameters of the fit are updated after a round of fitting

        # Add D0 specific B-decay ranges
        if "Psi" in self.resonance:
          self.BdecayMinVal =[10,20,20,20,20,15,15,15,0,0]  #for zT
          self.BdecayMaxVal =[80,80,80,80,70,70,60,50,30,30]#for zT
          
          #self.BdecayMinVal =[10,11,11,12,12,15,20,25,25] #for dR
          #self.BdecayMaxVal =[45,50,60,60,65,60,60,70,70] #for dR
        elif "X3872" in self.resonance:
          self.BdecayMinVal  =[0,0,0,0,0,0]
          self.BdecayMaxVal  =[70,65,65,60,25,15] #for zT
          #self.BdecayMinVal  =[10,11,11,20,20,30,25,30]
          #self.BdecayMaxVal  =[35,40,45,45,50,50,50,70]#for dR
        elif "D0" in self.resonance:
          self.BdecayMinVal  =[0,0,0,0,0]
          # self.BdecayMinVal  =[0,0,0,0,0,0]
          # self.BdecayMinVal = [5, 10, 10, 5, 2, 0]     # for zT
          self.BdecayMaxVal = [30, 25, 15, 10, 5]  # for zT
          # self.BdecayMaxVal = [30, 25, 20, 15, 10, 5]  # for zT
    
        #self.BdecayMinVal  =[20]
        #self.BdecayMaxVal  =[65]
  
        if len(self.BdecayMinVal)!=(nBins) or len(self.BdecayMaxVal)!=(nBins) :
          print("Error: Fix limits")
          print("number of mass bins= {}, number of limits for B-decay fraction: {}".format(nBins,len(self.BdecayMinVal)))
          # return
        
        if tFile is not None:
            print("tFiles: {}".format(tFile))
            self.tfile = ROOT.TFile(tFile)
        if pvFile is not None:
            print("pvFile: {}".format(pvFile))
            self.tfilePV = ROOT.TFile(pvFile)
        if rootFileName is not None:
            self.fInFileHisto = ROOT.TFile(rootFileName)
            self.nBins = nBins
            print("o Open data file: {}".format(self.fInFileHisto.GetName()))
            
            #If binned get histograms
            if isBinned:
              self.hTime=[]
              self.hYield=[]
              #nBins=9
              for i in range(0,nBins):
                timeHisto = self.fInFileHisto.Get("hTimeBKGsubtr_{}".format(i))
                self.hTime.append(timeHisto)
                print("o Add Histo: {}".format(timeHisto.GetName()))
                massHisto = self.fInFileHisto.Get("hMassSpectr_{}".format(i))
                #massHisto.Rebin(4)
                self.hYield.append(massHisto)
                #if i>0:
                 #self.hYield[0].Add(massHisto)
            else:
              #If unbinned get ttree
              self.inTree = self.fInFileHisto.Get("FragmNtuple")
              if self.inTree:
                print("Found ttree with name: {}".format(self.inTree.GetName()))
              else:
                print("No ttree in file")

    def initDictionary(self):
        self.massDict = {
        # Add D0 meson parameters
        "D0": {"sigma1": (0.008, 0.005, 0.012),
              "deltasigma": (1.5, 1.1, 2.5),
              "mean": (1.865, 1.860, 1.870),
              "alpha1": (2, 1, 5),
              "n": (1, 0.2, 5),
              "cb_frac": (0.5, 0.0, 0.99999),
              "pol0": (100, -5e3, 5e3),
              "pol1": (58, -2e2, 5e2),
              "pol2": (0, 0, 0),
              "mass_range": (1.81, 1.935),    # 150 MeV range
              "sig_yield": (100, 1, 5e6),
              "sig_yieldLim": (100, 1, 5e6),
              "bkg_yield": (1000, 0, 4e6),
              "bkg_yieldLim": (1000, 0, 4e6),
              "signal_region": (1.845, 1.885),  # 40 MeV range
              "SB_region": (1.820, 1.840)       # Sideband region
             }
        }
        self.ipchi2Dict = {
            # Add D0 IP chi2 parameters with Bukin function
            "SigD0": {
                "log_ipchi2_range": (-3, 5),           # Range for log(IP Chi2)
                
                # Prompt component (single Bukin parameters)
                "xp_prompt": (0.0, -0.5, 1.0),         # Peak position for prompt component
                "sigma_prompt": (0.7, 0.4, 2.0),       # Width parameter for prompt component
                "xi_prompt": (0.0, -0.5, 0.5),         # Asymmetry parameter for prompt component
                "rho1_prompt": (-0.2, -0.98, 0.98),        # Left tail parameter for prompt component
                "rho2_prompt": (0.2, -0.98, 0.98),        # Right tail parameter for prompt component
                
                # Non-prompt component (Bukin parameters)
                "xp_nonprompt": (2.5, 1.0, 4.0),       # Peak position for non-prompt component
                "sigma_nonprompt": (4.0, 0.4, 6.0),    # Width parameter for non-prompt component
                "xi_nonprompt": (0.1, 0.1, 0.1),      # More negative value creates rightward asymmetry
                "rho1_nonprompt": (-0.3, -0.8, -0.1),      # Higher positive value creates steeper left falloff                
                "rho2_nonprompt": (0.2, 0.0, 1.0),     # Right tail parameter - can remain unchanged                
                # Fraction parameters
                "prompt_frac": (0.95, 0.8, 0.999),     # Fraction of prompt component
                
                # Background parameters (unchanged)
                "bkg_param1": (0.5, 0, 1),             # Background shape parameter
                "bkg_param2": (0.5, 0, 1)              # Background shape parameter
            }
        }
        #self.cutVarsDict = {"jpsi":{"costheta":(-1.0,1.0),"pt_jet", pt_muon etc...
        #}}
     

    def updateDictionary(self, signal_pdf, data, fitFunc):

        res = self.massDict[self.resonance]

        RooArgSet = signal_pdf.getParameters(data)
        if RooArgSet:
          if fitFunc=="noSig":
            keyList = ["pol0","pol1","pol2"]
          #elif fitFunc=="DCB" or fitFunc=="DCBFixed":
          elif "DCB" in fitFunc:
            keyList = ["mean","sigma1","deltasigma","alpha1","n","cb_frac","sig_yield","bkg_yield","pol0","pol1","pol2"]
          elif fitFunc=="DGauss":
            keyList = ["mean", "sigma1","deltasigma","cb_frac","sig_yield","bkg_yield","pol0","pol1","pol2"]
          elif fitFunc=="SGauss":
            keyList = ["mean", "sigma1","sig_yield","bkg_yield","pol0","pol1","pol2"]
          else:
            return
            
          for key in keyList:
          
            newKeyVal    = RooArgSet.getRealValue(key,-1)
            #transform the tuple into a list to change its value
            tempKeyList  = list(res[key])
            #Assign a new value to the list element
            tempKeyList[0] = newKeyVal
            #Pass the list to the dictionary
            res[key] = tuple(tempKeyList)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #Update the ranges for the bdecay fraction
    def fixNAlphaValue(self, resonance , ptLimLow):
          
        res = self.massDict["{}".format(self.resonance)]
          
        if "Psi2S" in resonance:
          #These values were determined by
          #fitting the combined mass spectrum in jet pT 5-50GeV and zT: 0-1
          #The result of the MC was taken for this fixed range +-1Sigma
          #alpha and n were determined simultaniously in an unbinned fit
          startValAlpha=2.17   #2.461+-0.044  (before Apr. 22)
          minValAlpha  =startValAlpha-3*0.007
          maxValAlpha  =startValAlpha+3*0.007
          startValN=0.923  #0.554+-0.097
          minValN  =startValN-3*0.01
          maxValN  =startValN+3*0.01

        elif "X3872" in resonance:
          #These values were determined by
          #fitting the combined mass spectrum in jet pT 5-50GeV and zT: 0-1
          #The result of the MC was taken for this fixed range +-3Sigma
          #alpha and n were determined simultaniously in an unbinned fit
          startValAlpha=2.464 #1.625 (before Apr. 22)
          minValAlpha  =startValAlpha-3*0.025
          maxValAlpha  =startValAlpha+3*0.025
          startValN=0.593     #3.356 (before Apr. 22)
          minValN  =startValN-3*0.033
          maxValN  =startValN+3*0.033
       
        elif "D0" in resonance:
          # These values were determined by
          # fitting the combined mass spectrum in jet pT 5-50GeV and zT: 0-1
          # The result of the MC was taken for this fixed range +-3Sigma
          # alpha and n were determined simultaniously in an unbinned fit
          startValAlpha = 2.5
          minValAlpha   = startValAlpha - 3 * 0.05
          maxValAlpha   = startValAlpha + 3 * 0.05
          startValN     = 0.5
          minValN       = startValN - 3 * 0.05
          maxValN       = startValN + 3 * 0.05
        #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        #-transform the tuple into a list to change its value
        tempKeyList  = list(res["alpha1"])
        #-Assign a new value to the list element
        tempKeyList[0] = startValAlpha
        tempKeyList[1] = minValAlpha
        tempKeyList[2] = maxValAlpha
        #-Pass the list to the dictionary
        res["alpha1"] = tuple(tempKeyList)
        #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        #-transform the tuple into a list to change its value
        tempKeyList  = list(res["n"])
        #-Assign a new value to the list element
        tempKeyList[0] = startValN
        tempKeyList[1] = minValN
        tempKeyList[2] = maxValN
        #-Pass the list to the dictionary
        res["n"] = tuple(tempKeyList)
        #-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        if "Psi2S" in resonance:
          #-transform the tuple into a list to change its value
          tempKeyList  = list(res["cb_frac"])
          #-Assign a new value to the list element
          tempKeyList[0] = 0.4  #start val
          tempKeyList[1] = 0    #min val
          tempKeyList[2] = 0.5  #max val
          #-Pass the list to the dictionary
          res["cb_frac"] = tuple(tempKeyList)

#Update the ranges for the bdecay fraction
    def updateSBfraction(self, resonance, bin , ptLimLow):

        #print("Update b-frac limits:")
        #print("*This is the MC status: {}".format(self.isMC))

        res = self.sigLifetimeDict["Sig{}".format(self.resonance)]
        minVal = self.BdecayMinVal
        maxVal = self.BdecayMaxVal

        #-MC simulation are 100% B-decay
        if self.isMC:
            minVal  =[99,99,99,99,99,99,99,99,99,99,99]
            maxVal  =[100,100,100,100,100,100,100,100,100,100,100]
        else:
            minVal  =[0,0,0,0,0,0,0,0,0,0,0]
            maxVal  =[100,100,100,100,100,100,100,100,100,100,100]
       
        #-transform the tuple into a list to change its value
        tempKeyList  = list(res["bdecay_frac"])
        #print("**Before: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))
        #-Assign a new value to the list element
        tempKeyList[0] = 0.01*0.5*(minVal[bin]+maxVal[bin]) # 0.01 ->> change from frac into percent
        tempKeyList[1] = 0.01*minVal[bin]
        tempKeyList[2] = 0.01*maxVal[bin]
        #-Pass the list to the dictionary
        res["bdecay_frac"] = tuple(tempKeyList)
        #print("**After: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #Update the signal yield fix the yields to +- 20% of the yield from the mass fit
    def updateSigYield(self, resonance):

        print("Update SigYield limits:")
        res = self.massDict[resonance]
   
        #transform the tuple into a list to change its value
        tempKeyList  = list(res["sig_yield"])
        print("Sig yield Before: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))
        #Assign a new value to the list element
        tempKeyList[0] = tempKeyList[0]*1.0
        tempKeyList[1] = tempKeyList[0]*0.8
        if tempKeyList[1]<0:
          tempKeyList[1]=0
        tempKeyList[2] = tempKeyList[0]*1.2
        #Pass the list to the dictionary
        res["sig_yield"] = tuple(tempKeyList)
        print("Sig yield After: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))

        tempKeyList  = list(res["sig_yieldLim"])
        print("Sig yield lim Before: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))
        #Assign a new value to the list element
        tempKeyList[0] = tempKeyList[0]*1.0
        tempKeyList[1] = tempKeyList[0]*1.0
        tempKeyList[2] = tempKeyList[0]*1.0
        #Pass the list to the dictionary
        res["sig_yieldLim"] = tuple(tempKeyList)
        print("Sig lim yield After: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #Update the background yield fix the yields to +- 20% of the yield from the mass fit
    def updateBKGYield(self, resonance):

        print("Update SigYield limits:")
        res = self.massDict[resonance]
   
        #transform the tuple into a list to change its value
        tempKeyList  = list(res["bkg_yield"])
        print("BKG yield Before: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))
        #Assign a new value to the list element
        tempKeyList[0] = tempKeyList[0]*1.0
        tempKeyList[1] = tempKeyList[0]*0.8
        if tempKeyList[1]<0:
          tempKeyList[1]=0
        tempKeyList[2] = tempKeyList[0]*1.2
        #Pass the list to the dictionary
        res["bkg_yield"] = tuple(tempKeyList)
        print("BKG yield After: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))

        tempKeyList  = list(res["bkg_yieldLim"])
        print("BKG yield lim Before: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))
        #Assign a new value to the list element
        tempKeyList[0] = tempKeyList[0]*1.0
        tempKeyList[1] = tempKeyList[0]*1.0
        tempKeyList[2] = tempKeyList[0]*1.0
        #Pass the list to the dictionary
        res["bkg_yieldLim"] = tuple(tempKeyList)
        print("BKG lim yield After: {}, {}, {}".format(tempKeyList[0],tempKeyList[1],tempKeyList[2]))


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # @classmethod or @staticmethod...
    def fiducialCutString(self,**kwargs):
        """
        Define fiducial variables which want to cut on. Inputs should be tuples. Returns a string.
        """
        # check ntuple names, replace jpsi with tag eventually
        output = []
        for fidName, binVals in kwargs.items():
            if len(binVals) != 2: continue
            if binVals[0] != None:
                output += [f"{fidName} > {binVals[0]}"]
            if binVals[1] != None:
                output += [f"{fidName} < {binVals[1]}"]
        #join fiducial requirement strings together
        finalFidString = " && ".join(output)                    

        return finalFidString
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def createDataSet(self, resonance, name, fidCutString="", isMass=True, bin=0, corrVer=-1, fitVarName="ipchi2"):
        """Create a RooDataSet from the input file."""
                
        #-binned
        if self.isBinned:
            res = self.massDict[resonance]
            ipchi2_params = self.ipchi2Dict[f"Sig{resonance}"]

            print("Binned Histogram Fit")
            
            # Check if we have histograms loaded
            if not hasattr(self, 'hYield') or len(self.hYield) == 0:
                print("ERROR: No histograms found for binned mode.")
                return None
            
            if isMass:
                fitVariable = ROOT.RooRealVar("tagMass", "tagMass", *res["mass_range"])
                fitVariable.setRange("fullRange", *res["mass_range"])
                
                # Make sure bin is within range
                if bin < 0 or bin >= len(self.hYield):
                    print(f"WARNING: bin {bin} out of range (0-{len(self.hYield)-1}). Using bin 0.")
                    bin = 0
                    
                data = ROOT.RooDataHist("dataM", "dataM", ROOT.RooArgList(fitVariable), self.hYield[bin])
            else:
                # For IP chi2 in binned mode
                if not hasattr(self, 'hIPChi2'):
                    print("WARNING: No IP chi2 histograms available. Attempting to create them.")
                    self.hIPChi2 = []
                    for i in range(len(self.hYield)):
                        # Create dummy histograms
                        hist = ROOT.TH1F(f"hIPChi2_{i}", f"log(IP Chi2) {i}", 50, 
                                       ipchi2_params["log_ipchi2_range"][0], 
                                       ipchi2_params["log_ipchi2_range"][1])
                        self.hIPChi2.append(hist)
                
                log_ipchi2 = ROOT.RooRealVar("log_tag_ipchi2", "log(tag_ip_chi2)", 
                                           *ipchi2_params["log_ipchi2_range"])
                
                # Make sure bin is within range
                if bin < 0 or bin >= len(self.hIPChi2):
                    print(f"WARNING: bin {bin} out of range (0-{len(self.hIPChi2)-1}). Using bin 0.")
                    bin = 0
                    
                data = ROOT.RooDataHist("dataIPChi2", "dataIPChi2", ROOT.RooArgList(log_ipchi2), self.hIPChi2[bin])

            return data
        
        #-unbinned
        else:
            print("Unbinned Histogram Fit")
            
            # Check if inTree attribute exists
            if not hasattr(self, 'inTree') or self.inTree is None:
                print("ERROR: No TTree found in input file.")
                return None
                    
            # Get range from dictionary
            res = self.massDict[resonance]
            ipchi2_params = self.ipchi2Dict[f"Sig{resonance}"]
            mass_range = res["mass_range"]
            
            # Declare variables
            tagMass = ROOT.RooRealVar("tagMass", "tagMass", mass_range[0], mass_range[1])
            
            # Use tag_ip_chi2 for IP chi2 fits (assuming it exists in the tree)
            tag_ipchi2 = ROOT.RooRealVar("tag_ip_chi2", "tag_ip_chi2", 0, 10000)
            log_tag_ipchi2 = ROOT.RooRealVar("log_tag_ipchi2", "log_tag_ipchi2", 
                                           ipchi2_params["log_ipchi2_range"][0], 
                                           ipchi2_params["log_ipchi2_range"][1])

            # Other variables (unchanged)
            pt_jet = ROOT.RooRealVar("jetPt","jetPt",0,200)
            pt_tag = ROOT.RooRealVar("tagPt","tagPt",0,200)
            nConst = ROOT.RooRealVar("jetnConst","jetnConst",0.0,300.0)
            qValue = ROOT.RooRealVar("QValue","QValue",-2,0.5)
            dRValue = ROOT.RooRealVar("tagJetdR","tagJetdR",0.0,1.0)
            tagZ = ROOT.RooRealVar("tagZ","tagZ",0.0,1.01)

            # Create RooArgSet with all variables
            cutVars = ROOT.RooArgSet()
            cutVars.add(tagMass)
            cutVars.add(tag_ipchi2)
            cutVars.add(log_tag_ipchi2)
            cutVars.add(pt_jet)
            cutVars.add(pt_tag)
            cutVars.add(nConst)
            cutVars.add(qValue)
            cutVars.add(tagZ)
            cutVars.add(dRValue)
            
            # Add other variables (unchanged)
            distance1 = ROOT.RooRealVar("Distance1","Distance1",-10,200)
            distance2 = ROOT.RooRealVar("Distance2","Distance2",-10,200)
            distance3 = ROOT.RooRealVar("Distance3","Distance3",-10,200)
            cutVars.add(distance1)
            cutVars.add(distance2)
            cutVars.add(distance3)
        
            # Weights (unchanged)
            effWeight = ROOT.RooRealVar("EffWeight","EffWeight",0,25000)
            effWeight_0 = ROOT.RooRealVar("EffWeight_0","EffWeight_0",0,5)
            effWeight_1 = ROOT.RooRealVar("EffWeight_1","EffWeight_1",0,500)
            effWeight_2 = ROOT.RooRealVar("EffWeight_2","EffWeight_2",0,4)
            effWeight_3 = ROOT.RooRealVar("EffWeight_3","EffWeight_3",-1,2)
            effWeight_4 = ROOT.RooRealVar("EffWeight_4","EffWeight_4",0,2000)
            effWeight_Rnd = ROOT.RooRealVar("tagnRnd","tagnRnd",-1,25)
            cutVars.add(effWeight)
            cutVars.add(effWeight_0)
            cutVars.add(effWeight_1)
            cutVars.add(effWeight_2)
            cutVars.add(effWeight_3)
            cutVars.add(effWeight_4)
            cutVars.add(effWeight_Rnd)

            # Create final cut string
            finalCutString = fidCutString
            finalCutString += "&& tagPt > 2"
            finalCutString += "&& jetnConst > 1"
            
            if corrVer > -2:
                finalCutString += "&& Distance1 < 0.5 && Distance2 < 0.5 && Distance3 < 0.5"
                
            # REMOVE the log calculation from cut string - we'll handle it differently
            # finalCutString += "&& log_tag_ipchi2 = log(max(tag_ip_chi2, 0.001))"

            # First create dataset without log_tag_ipchi2 constraint
            if corrVer < 0:
                data = ROOT.RooDataSet(name, name, self.inTree, cutVars, finalCutString)
            elif corrVer == 0:
                data = ROOT.RooDataSet(name, name, self.inTree, cutVars, finalCutString, effWeight_Rnd.GetName())
            elif corrVer == 1:
                data = ROOT.RooDataSet(name, name, self.inTree, cutVars, finalCutString, effWeight.GetName())
            elif corrVer==2:
                data = ROOT.RooDataSet(name, name, self.inTree, cutVars, finalCutString, effWeight_0.GetName() )
            elif corrVer==3:
                data = ROOT.RooDataSet(name, name, self.inTree, cutVars, finalCutString, effWeight_1.GetName() )
            elif corrVer==4:
                data = ROOT.RooDataSet(name, name, self.inTree, cutVars, finalCutString, effWeight_2.GetName() )
            elif corrVer==5:
                data = ROOT.RooDataSet(name, name, self.inTree, cutVars, finalCutString, effWeight_3.GetName() )
            elif corrVer==6:
                data = ROOT.RooDataSet(name, name, self.inTree, cutVars, finalCutString, effWeight_4.GetName() )

            # # If the dataset is created successfully, create a new dataset with log(IP chi2)
            # if data and data.numEntries() > 0:
            #     # Create a new RooRealVar for the log(tag_ip_chi2) - this is a fundamental type
            #     log_tag_ipchi2_var = ROOT.RooRealVar("log_tag_ipchi2", "log(tag_ip_chi2)", 
            #                                       ipchi2_params["log_ipchi2_range"][0], 
            #                                       ipchi2_params["log_ipchi2_range"][1])
                
            #     # Create a new RooArgSet for the new dataset
            #     # Start with all variables from the original dataset
            #     newVars = ROOT.RooArgSet()
            #     iter = data.get().createIterator()
            #     var = iter.Next()
            #     while var:
            #         # Clone each variable to avoid reference issues
            #         if var.GetName() != "log_tag_ipchi2": # Skip if already exists
            #             newVars.add(var)
            #         var = iter.Next()
                
            #     # Add our new log_tag_ipchi2 variable
            #     newVars.add(log_tag_ipchi2_var)
                
            #     # Create a new dataset with the same name plus suffix
            #     weightVarName = data.weightVar().GetName() if data.isWeighted() else ""
                
            #     # Create a new empty dataset with all variables
            #     newData = ROOT.RooDataSet(f"{name}_with_log", f"{name} with log(IP chi2)",
            #                             newVars, weightVarName)
                
            #     # Loop through the original dataset and fill the new one
            #     for i in range(data.numEntries()):
            #         # Get original entry
            #         entry = data.get(i)
                    
            #         # Find the current IP chi2 value
            #         ip_chi2_val = entry.find("tag_ip_chi2").getVal() if entry.find("tag_ip_chi2") else 0.001
                    
            #         # Calculate log(IP chi2) with safety check
            #         if ip_chi2_val <= 0:
            #             log_val = ipchi2_params["log_ipchi2_range"][0]  # Use min value for invalid inputs
            #         else:
            #             log_val = math.log(ip_chi2_val)
            #             # Clip to valid range
            #             log_val = max(ipchi2_params["log_ipchi2_range"][0], 
            #                          min(ipchi2_params["log_ipchi2_range"][1], log_val))
                    
            #         # Set the log value in our variable
            #         log_tag_ipchi2_var.setVal(log_val)
                    
            #         # Copy the weight if dataset is weighted
            #         weight = data.weight() if data.isWeighted() else 1.0
                    
            #         # Add entry to new dataset
            #         newData.add(newVars, weight)
                
            #     print(f"  Created new dataset with {newData.numEntries()} entries including log(IP chi2)")
            #     return newData
            
            return data
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def massFit(self, resonance, data, fitTypeName = "DCB", bin = -1, zRange = "", splot = False, sFile = None):
        """
        Unbinned fit of mass of 5 different resonances. 
        """
        print(f"\n==== Starting massFit with {fitTypeName} model for bin {bin} ====")
        
        # Suppress RooFit messages except errors
        originalMsgLevel = ROOT.RooMsgService.instance().globalKillBelow()
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
        
        try:
            # declare mass variable, signal region and overall region
            fitRange="fullRange"
            res = self.massDict[resonance]
            fullRange   = res["mass_range"]
            signalRange = res["signal_region"]
            SBRange     = res["SB_region"]
            
            print(f"  Mass range: {fullRange}, Signal region: {signalRange}")
            print(f"  Data entries: {data.numEntries()}")
            
            # DEBUG: Print dictionary parameters at start
            print("\n=== Initial parameter values from dictionary for resonance ===")
            print(f"  Resonance: {resonance}")
            for key, value in res.items():
                if isinstance(value, tuple):
                    if len(value) >= 3:
                        print(f"  {key}: {value[0]} (range: {value[1]} to {value[2]})")
                    elif len(value) == 2:
                        print(f"  {key}: {value[0]} (range: {value[1]} to N/A)")
                    elif len(value) == 1:
                        print(f"  {key}: {value[0]} (no range specified)")
                    else:
                        print(f"  {key}: {value} (empty tuple)")
                else:
                    print(f"  {key}: {value}")
            
            mass_tag_measured = ROOT.RooRealVar("tagMass", "tagMass", *res["mass_range"])
            mass_tag_measured.setRange("fullRange",*res["mass_range"])
            mass_tag_measured.setRange("signalRange",*res["signal_region"])
            mass_tag_measured.setRange("SBleft",fullRange[0],SBRange[0])
            mass_tag_measured.setRange("SBright",SBRange[1],fullRange[1])

            # declare Signal parameters
            if "BKG" not in resonance:
                sigma1     = ROOT.RooRealVar("sigma1","sigma1",*res["sigma1"])
                deltasigma = ROOT.RooRealVar("deltasigma","deltasigma",*res["deltasigma"])
                sigma2     = ROOT.RooFormulaVar("sigma2", "sigma2", "sigma1*deltasigma", ROOT.RooArgList(sigma1,deltasigma))
                mean       = ROOT.RooRealVar("mean","mean",*res["mean"])
                alpha1     = ROOT.RooRealVar("alpha1","alpha1",*res["alpha1"])
                alpha2     = ROOT.RooFormulaVar("alpha2", "alpha2", "-1*alpha1", ROOT.RooArgList(alpha1))
                n          = ROOT.RooRealVar("n","n",*res["n"])
                cb_frac    = ROOT.RooRealVar("cb_frac","cb_frac",*res["cb_frac"])
                
                sig_yield  = ROOT.RooRealVar("sig_yield","sig_yield",*res["sig_yield"])
                
                # DEBUG: Print initial signal parameters
                print("\n=== Initial signal fit parameters ===")
                print(f"  sigma1: {sigma1.getVal()} (range: {sigma1.getMin()} to {sigma1.getMax()})")
                print(f"  deltasigma: {deltasigma.getVal()} (range: {deltasigma.getMin()} to {deltasigma.getMax()})")
                print(f"  sigma2: {sigma2.getVal()} (derived)")
                print(f"  mean: {mean.getVal()} (range: {mean.getMin()} to {mean.getMax()})")
                print(f"  alpha1: {alpha1.getVal()} (range: {alpha1.getMin()} to {alpha1.getMax()})")
                print(f"  alpha2: {alpha2.getVal()} (derived)")
                print(f"  n: {n.getVal()} (range: {n.getMin()} to {n.getMax()})")
                print(f"  cb_frac: {cb_frac.getVal()} (range: {cb_frac.getMin()} to {cb_frac.getMax()})")
                print(f"  sig_yield: {sig_yield.getVal()} (range: {sig_yield.getMin()} to {sig_yield.getMax()})")

                #double-sided crystal ball function
                if "DCB" in fitTypeName:
                    CB1_pdf    = ROOT.RooCBShape("Sig1_pdf","Sig1_pdf",mass_tag_measured,mean,sigma1,alpha1,n)
                    CB2_pdf    = ROOT.RooCBShape("Sig2_pdf","Sig2_pdf",mass_tag_measured,mean,sigma2,alpha2,n)
                    sig_pdf    = ROOT.RooAddPdf("sig_pdf", "Signal", ROOT.RooArgList(CB2_pdf,CB1_pdf), ROOT.RooArgList(cb_frac))
                    print(f"  Using Double Crystal Ball PDF (CB2*{cb_frac.getVal()} + CB1*{1-cb_frac.getVal()})")
                
                #double Gauss function
                if fitTypeName == "DGauss":
                    Gauss1_pdf    = ROOT.RooGaussian("Sig1_pdf","Sig1_pdf",mass_tag_measured,mean,sigma1)
                    Gauss2_pdf    = ROOT.RooGaussian("Sig2_pdf","Sig2_pdf",mass_tag_measured,mean,sigma2)
                    sig_pdf       = ROOT.RooAddPdf("sig_pdf", "Signal", ROOT.RooArgList(Gauss2_pdf,Gauss1_pdf), ROOT.RooArgList(cb_frac))
                    print(f"  Using Double Gaussian PDF (G2*{cb_frac.getVal()} + G1*{1-cb_frac.getVal()})")
                
                #single Gauss function
                if fitTypeName == "SGauss" or fitTypeName == "noSig":
                    Gauss1_pdf    = ROOT.RooGaussian("Gauss1_pdf","Gauss1_pdf",mass_tag_measured,mean,sigma1)
                    sig_pdf       = Gauss1_pdf
                    print(f"  Using Single Gaussian PDF (sigma: {sigma1.getVal()})")
                
                # signal pdf
                sig_pdf_ext  = ROOT.RooExtendPdf("sig_pdf_ext", "sig_pdf_ext", sig_pdf, sig_yield, fitRange)

            # Always use polynomial background regardless of fit type
            bkg_yield = ROOT.RooRealVar("bkg_yield", "bkg_yield", *res["bkg_yield"])
            
            # DEBUG: Print initial background parameters
            print("\n=== Initial background fit parameters ===")
            print(f"  bkg_yield: {bkg_yield.getVal()} (range: {bkg_yield.getMin()} to {bkg_yield.getMax()})")
            
            # Always use polynomial background
            poly0 = ROOT.RooRealVar("pol0", "pol0", *res["pol0"])
            poly1 = ROOT.RooRealVar("pol1", "pol1", *res["pol1"])
            bkg_pdf = ROOT.RooPolynomial("bkg_pdf", "bkg_pdf", mass_tag_measured, ROOT.RooArgList(poly0, poly1))
            print(f"  Using Polynomial background (p0: {poly0.getVal()}, p1: {poly1.getVal()})")
            
            # background pdf
            bkg_pdf_ext = ROOT.RooExtendPdf("bkg_pdf_ext", "bkg_pdf_ext", bkg_pdf, bkg_yield, fitRange)

            #-Build the PDF that will be fit to the data
            if fitTypeName == "noSig":
                extended_pdf = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(bkg_pdf_ext))
            else:
                extended_pdf = ROOT.RooAddPdf("model", "model", ROOT.RooArgList(sig_pdf_ext, bkg_pdf_ext))

            # Set fit print level to provide more information
            fitOptions = ROOT.RooFit.Save(True)
            fitOptions = ROOT.RooFit.PrintLevel(1)
            
            #-Fit the data
            print("\n=== Starting fit ===")
            if fitTypeName == "noSig":
                print("  Fitting in sideband regions only")
                fit_result = extended_pdf.fitTo(data, fitOptions, ROOT.RooFit.Range("SBleft,SBright"))
            else:
                print("  Fitting in full range")
                fit_result = extended_pdf.fitTo(data, fitOptions, ROOT.RooFit.Range(fitRange))
            
            if hasattr(self, 'updateStartValues') and self.updateStartValues:
                if hasattr(self, 'updateDictionary'):
                    self.updateDictionary(extended_pdf,data,fitTypeName)
                else:
                    print("Warning: updateDictionary method not available")

            # # DEBUG: Print fit status
            # print(f"\n=== Fit result status: {fit_result.status()} ===")
            # print(f"  Covariance matrix quality: {fit_result.covQual()}")
            # print(f"  EDM: {fit_result.edm()}")
            # print(f"  Number of function calls: {fit_result.numFitCtrl().ncalls}")
            
            integral_fullB = bkg_pdf_ext.createIntegral(mass_tag_measured,ROOT.RooFit.NormSet(mass_tag_measured),ROOT.RooFit.Range("fullRange"))
            integral_sigB  = bkg_pdf_ext.createIntegral(mass_tag_measured,ROOT.RooFit.NormSet(mass_tag_measured),ROOT.RooFit.Range("signalRange"))
            integral_fullS = sig_pdf_ext.createIntegral(mass_tag_measured,ROOT.RooFit.NormSet(mass_tag_measured),ROOT.RooFit.Range("fullRange"))
            integral_sigS  = sig_pdf_ext.createIntegral(mass_tag_measured,ROOT.RooFit.NormSet(mass_tag_measured),ROOT.RooFit.Range("signalRange"))

            parameters = extended_pdf.getParameters(data) #returns RooArgSet
            
            # DEBUG: Print all parameters after fit
            print("\n=== Final fit parameters ===")
            iter = parameters.createIterator()
            param = iter.Next()
            while param:
                print(f"  {param.GetName()}: {param.getVal()} ± {param.getError()}")
                param = iter.Next()
            
            #-Now some conversions:
            #-The integral integrates a shape. And thus gives no information on the yield
            #-under the integral in terms of event. This needs to be determined separatley
            #-https://root-forum.cern.ch/t/createintegral-gives-unexpected-result/32627/3
            fullRangeSYield = parameters.getRealValue("sig_yield",-1)
            fullRangeBYield = parameters.getRealValue("bkg_yield",-1)
            
            #-Scale factor SF=Nevt_full/FullInt_PDF -> Nevt_limited=FullInt_PDF*SF
            if integral_fullB:
              SfactorBKG = fullRangeBYield/integral_fullB.getVal()
            else:
              SfactorBKG = 0
            if integral_fullS:
              SfactorS   = fullRangeSYield/integral_fullS.getVal()
            else:
              SfactorS = 0
            NevtBKG_SignalRange = SfactorBKG*integral_sigB.getVal()
            NevtS_SignalRange   = SfactorS*integral_sigS.getVal()
            
            # DEBUG: Print integral calculations
            print("\n=== Yield calculations ===")
            print(f"  Full range signal yield: {fullRangeSYield}")
            print(f"  Full range background yield: {fullRangeBYield}")
            print(f"  Signal fraction in signal region: {integral_sigS.getVal()/integral_fullS.getVal():.4f}")
            print(f"  Background fraction in signal region: {integral_sigB.getVal()/integral_fullB.getVal():.4f}")
            print(f"  Signal yield in signal region: {NevtS_SignalRange:.2f}")
            print(f"  Background yield in signal region: {NevtBKG_SignalRange:.2f}")
            print(f"  S/B in signal region: {NevtS_SignalRange/NevtBKG_SignalRange if NevtBKG_SignalRange > 0 else 'Inf'}")
            print(f"  Significance (S/√(S+B)): {NevtS_SignalRange/math.sqrt(NevtS_SignalRange+NevtBKG_SignalRange) if NevtS_SignalRange+NevtBKG_SignalRange > 0 else 'NaN'}")

            #Plot fit spectra
            plot = plt.Plotter(resonance,bin=bin,range=zRange, basepath=self.outfilePath, binned= self.isBinned)
            histogram = plot.individualMassFitPlot(sig_yield, extended_pdf, mass_tag_measured, data, fitTypeName)
            
            parameters = extended_pdf.getParameters(data) #returns RooArgSet

            parameterArr    = [0]*12
            parameterErrArr = [0]*10
            
            parameterArr[0]  = parameters.getRealValue("sig_yield",-1)
            parameterArr[1]  = parameters.getRealValue("bkg_yield",-1)
            parameterArr[2]  = parameters.getRealValue("mean",-1)
            parameterArr[3]  = parameters.getRealValue("sigma1",-1)
            parameterArr[4]  = parameters.getRealValue("deltasigma",-1)
            parameterArr[5]  = parameters.getRealValue("alpha1",-1)
            parameterArr[6]  = parameters.getRealValue("n",-1)
            parameterArr[7]  = parameters.getRealValue("cb_frac",-1)
            
            # When extracting parameters, always use polynomial parameters
            parameterArr[8]  = parameters.getRealValue("pol0",-1)
            parameterArr[9]  = parameters.getRealValue("pol1",-1)
            
            # Special info, not from fitter
            parameterArr[10]  = NevtS_SignalRange
            parameterArr[11]  = NevtBKG_SignalRange
            
            if fitTypeName == "noSig":
              parameterErrArr[0]  = 0
            else:
              parameterErrArr[0] = sig_yield.getError()
            parameterErrArr[1] = bkg_yield.getError()
            parameterErrArr[2] = mean.getError()
            parameterErrArr[3] = sigma1.getError()
            parameterErrArr[4] = deltasigma.getError()
            parameterErrArr[5] = alpha1.getError()
            parameterErrArr[6] = n.getError()
            parameterErrArr[7] = cb_frac.getError()
            
            # When extracting parameter errors, always use polynomial errors
            parameterErrArr[8] = poly0.getError()
            parameterErrArr[9] = poly1.getError()
            
            # DEBUG: Print summary of fit results arrays
            print("\n=== Fit results summary (parameter arrays) ===")
            param_names = ["sig_yield", "bkg_yield", "mean", "sigma1", "deltasigma", "alpha1", "n", "cb_frac", 
                          "pol0", "pol1", "NevtS_SignalRange", "NevtBKG_SignalRange"]
            for i, name in enumerate(param_names):
                if i < len(parameterErrArr):
                    print(f"  {name}: {parameterArr[i]} ± {parameterErrArr[i]}")
                else:
                    print(f"  {name}: {parameterArr[i]}")
            
            print(f"\n==== {fitTypeName} fit completed successfully ====")
            return histogram, parameterArr, parameterErrArr

        except Exception as e:
            print(f"Error in massFit: {str(e)}")
            import traceback
            traceback.print_exc()
            return None, [0]*12, [0]*10


    def ipchi2Fit(self, resonance, data, background=None, figKey="All", bin=0, zRange=""):
        """
        Perform IP chi2 fit for a given resonance using single Bukin functions for prompt and non-prompt components.
        Each component is modeled with a single Bukin PDF for better simplicity and stability.
        """
        print(f"\n==== Starting IP chi2 fit with Bukin function for bin {bin} ====")
        
        # Suppress RooFit messages
        originalMsgLevel = ROOT.RooMsgService.instance().globalKillBelow()
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

        # Initialize a list to keep track of histograms for cleanup
        histograms_to_cleanup = []

        try:
            # Create plotting object
            plot = plt.Plotter(resonance, bin=bin, basepath=self.outfilePath, binned=self.isBinned, range=zRange, name=figKey)

            # Get dictionaries with parameters
            ipchi2_params = self.ipchi2Dict[f"Sig{resonance}"]
            mass_params = self.massDict[resonance]
            
            # Print dataset info
            print(f"  Signal dataset entries: {data.numEntries()}")
            if background:
                print(f"  Background dataset entries: {background.numEntries()}")
            
            # Create RooFit variable for log(IP Chi2)
            log_ipchi2 = ROOT.RooRealVar("log_tag_ipchi2", "log(tag_ip_chi2)", *ipchi2_params["log_ipchi2_range"])
            
            # Create variables for the model
            sig_yield = ROOT.RooRealVar("sig_yield", "sig_yield", *mass_params["sig_yield"])
            sig_yieldLim = ROOT.RooRealVar("sig_yieldLim", "sig_yieldLim", *mass_params["sig_yieldLim"])
            bkg_yieldLim = ROOT.RooRealVar("bkg_yieldLim", "bkg_yieldLim", *mass_params["bkg_yieldLim"])
            prompt_frac = ROOT.RooRealVar("prompt_frac", "prompt_frac", *ipchi2_params["prompt_frac"])
            
            # Create single Bukin for prompt component
            xp_prompt = ROOT.RooRealVar("xp_prompt", "xp_prompt", *ipchi2_params["xp_prompt"])
            sigma_prompt = ROOT.RooRealVar("sigma_prompt", "sigma_prompt", *ipchi2_params["sigma_prompt"])
            xi_prompt = ROOT.RooRealVar("xi_prompt", "xi_prompt", *ipchi2_params["xi_prompt"])
            rho1_prompt = ROOT.RooRealVar("rho1_prompt", "rho1_prompt", *ipchi2_params["rho1_prompt"])
            rho2_prompt = ROOT.RooRealVar("rho2_prompt", "rho2_prompt", *ipchi2_params["rho2_prompt"])
            
            # Create the prompt Bukin PDF
            prompt_pdf = ROOT.RooBukinPdf("prompt_pdf", "prompt_pdf", 
                                         log_ipchi2, xp_prompt, sigma_prompt, xi_prompt, 
                                         rho1_prompt, rho2_prompt)
            
            # Create non-prompt component using a single Bukin
            xp_nonprompt = ROOT.RooRealVar("xp_nonprompt", "xp_nonprompt", *ipchi2_params["xp_nonprompt"])
            sigma_nonprompt = ROOT.RooRealVar("sigma_nonprompt", "sigma_nonprompt", *ipchi2_params["sigma_nonprompt"])
            xi_nonprompt = ROOT.RooRealVar("xi_nonprompt", "xi_nonprompt", *ipchi2_params["xi_nonprompt"])
            rho1_nonprompt = ROOT.RooRealVar("rho1_nonprompt", "rho1_nonprompt", *ipchi2_params["rho1_nonprompt"])
            rho2_nonprompt = ROOT.RooRealVar("rho2_nonprompt", "rho2_nonprompt", *ipchi2_params["rho2_nonprompt"])
            
            nonprompt_pdf = ROOT.RooBukinPdf("nonprompt_pdf", "nonprompt_pdf", 
                                           log_ipchi2, xp_nonprompt, sigma_nonprompt, xi_nonprompt, 
                                           rho1_nonprompt, rho2_nonprompt)
            
            # Calculate yields for prompt and non-prompt
            prompt_yield = ROOT.RooFormulaVar("prompt_yield", "prompt_yield", "sig_yieldLim*prompt_frac", 
                                             ROOT.RooArgList(sig_yieldLim, prompt_frac))
            nonprompt_yield = ROOT.RooFormulaVar("nonprompt_yield", "nonprompt_yield", "sig_yieldLim*(1-prompt_frac)", 
                                                ROOT.RooArgList(sig_yieldLim, prompt_frac))
            
            # Signal-only model
            model = ROOT.RooAddPdf("model", "model", 
                                    ROOT.RooArgList(prompt_pdf, nonprompt_pdf),
                                    ROOT.RooArgList(prompt_yield, nonprompt_yield))
            bkg_pdf = None
                
            # Perform the fit
            print("  Performing IP chi2 fit with Bukin function...")
            cmd_args = ROOT.RooLinkedList()
            cmd_args.Add(ROOT.RooFit.Save(True))
            cmd_args.Add(ROOT.RooFit.PrintLevel(1))  # More verbose output
            
            result = model.fitTo(data, cmd_args)
            
            # Get parameters after fit
            parameters = model.getParameters(data)
            
            # Create plot
            print("  Creating IP chi2 fit plot...")
            histogram = plot.ipchi2FitPlot(resonance, log_ipchi2, data, model, 
                                         nonprompt_pdf, prompt_pdf, bkg_pdf,
                                         prompt_yield, nonprompt_yield)
            
            # Extract fit parameters for simplified Bukin model
            parameterArr = [0]*12
            parameterErrArr = [0]*12
            
            # Basic yield parameters
            parameterArr[0] = sig_yield.getVal()
            parameterArr[1] = prompt_frac.getVal()
            
            # Prompt Bukin parameters
            parameterArr[2] = xp_prompt.getVal()
            parameterArr[3] = sigma_prompt.getVal()
            parameterArr[4] = xi_prompt.getVal()
            parameterArr[5] = rho1_prompt.getVal()
            parameterArr[6] = rho2_prompt.getVal()
            
            # Non-prompt Bukin parameters
            parameterArr[7] = xp_nonprompt.getVal()
            parameterArr[8] = sigma_nonprompt.getVal()
            parameterArr[9] = xi_nonprompt.getVal()
            parameterArr[10] = rho1_nonprompt.getVal()
            parameterArr[11] = rho2_nonprompt.getVal()
            
            # Extract errors for key parameters
            parameterErrArr[0] = sig_yield.getError()
            parameterErrArr[1] = prompt_frac.getError()
            
            parameterErrArr[2] = xp_prompt.getError()
            parameterErrArr[3] = sigma_prompt.getError()
            parameterErrArr[4] = xi_prompt.getError()
            
            parameterErrArr[7] = xp_nonprompt.getError()
            parameterErrArr[8] = sigma_nonprompt.getError()
            parameterErrArr[9] = xi_nonprompt.getError()
            
            print("  IP chi2 fit with simplified Bukin function completed successfully")
            print(f"  Prompt fraction: {prompt_frac.getVal():.3f} ± {prompt_frac.getError():.3f}")
            print(f"  Prompt asymmetry parameter: xi={xi_prompt.getVal():.3f}")
            print(f"  Non-prompt asymmetry parameter: xi={xi_nonprompt.getVal():.3f}")
            
            # Clean up histograms
            for hist in histograms_to_cleanup:
                if hist:
                    hist.SetDirectory(0)
                    hist.Delete()
            
            return histogram, parameterArr, parameterErrArr
            
        except Exception as e:
            print(f"Error in ipchi2Fit with Bukin function: {str(e)}")
            import traceback
            traceback.print_exc()
            return None, [0]*12, [0]*12

