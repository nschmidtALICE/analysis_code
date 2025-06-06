
def analyze_mc_efficiency(ntuple_file):
    """
    Analyze reconstruction efficiency and purity using MC-truth information
    from the ntuple file.
    
    Parameters:
    ----------
    ntuple_file : str
        Path to the ROOT ntuple file
    """
    import ROOT
    
    # Open the ROOT file
    root_file = ROOT.TFile(ntuple_file, "READ")
    if not root_file.IsOpen():
        print(f"Could not open {ntuple_file}")
        return
    
    # Get the tree
    tree = root_file.Get("d0jets")
    if not tree:
        print("Could not find 'd0jets' tree in file")
        root_file.Close()
        return
    
    # Create histograms for efficiency studies
    h_mc_d0_pt = ROOT.TH1F("h_mc_d0_pt", "MC D0 p_{T};p_{T} [GeV/c];Entries", 100, 0, 50)
    h_reco_d0_pt = ROOT.TH1F("h_reco_d0_pt", "Reconstructed D0 p_{T};p_{T} [GeV/c];Entries", 100, 0, 50)
    
    h_mc_d0_eta = ROOT.TH1F("h_mc_d0_eta", "MC D0 #eta;#eta;Entries", 50, 1.5, 5.0)
    h_reco_d0_eta = ROOT.TH1F("h_reco_d0_eta", "Reconstructed D0 #eta;#eta;Entries", 50, 1.5, 5.0)
    
    h_mc_d0_z = ROOT.TH1F("h_mc_d0_z", "MC D0 Fragmentation Fraction;z = p_{T}^{D0}/p_{T}^{jet};Entries", 50, 0, 1)
    h_reco_d0_z = ROOT.TH1F("h_reco_d0_z", "Reconstructed D0 Fragmentation Fraction;z = p_{T}^{D0}/p_{T}^{jet};Entries", 50, 0, 1)
    
    h_d0_pt_res = ROOT.TH1F("h_d0_pt_res", "D0 p_{T} Resolution;(p_{T}^{reco} - p_{T}^{MC})/p_{T}^{MC};Entries", 100, -0.5, 0.5)
    h_d0_mass_res = ROOT.TH1F("h_d0_mass_res", "D0 Mass Resolution;(M^{reco} - M^{MC})/M^{MC};Entries", 100, -0.05, 0.05)
    
    h_origin = ROOT.TH1F("h_origin", "D0 Origin;Origin;Entries", 3, 0.5, 3.5)
    h_origin.GetXaxis().SetBinLabel(1, "Prompt")
    h_origin.GetXaxis().SetBinLabel(2, "From B")
    h_origin.GetXaxis().SetBinLabel(3, "Other")
    
    # Loop over all entries
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        
        # Process MC D0s
        for i in range(tree.mc_d0_pt.size()):
            h_mc_d0_pt.Fill(tree.mc_d0_pt[i])
            h_mc_d0_eta.Fill(tree.mc_d0_eta[i])
            h_origin.Fill(tree.mc_d0_origin[i])
            
            # Check if D0 is in a jet
            if tree.mc_d0_in_jet[i] > 0:
                h_mc_d0_z.Fill(tree.mc_d0_z[i])
            
            # Check if D0 was reconstructed
            if tree.mc_d0_matched[i] >= 0:
                reco_idx = tree.mc_d0_matched[i]
                h_reco_d0_pt.Fill(tree.mc_d0_pt[i])
                h_reco_d0_eta.Fill(tree.mc_d0_eta[i])
                
                if tree.d0_in_jet[reco_idx] > 0:
                    h_reco_d0_z.Fill(tree.d0_z[reco_idx])
        
        # Process resolution
        for i in range(tree.d0_pt_res.size()):
            h_d0_pt_res.Fill(tree.d0_pt_res[i])
        
        for i in range(tree.d0_mass_res.size()):
            h_d0_mass_res.Fill(tree.d0_mass_res[i])
    
    # Calculate and draw efficiency plots
    c1 = ROOT.TCanvas("c1", "D0 Reconstruction Efficiency", 1200, 800)
    c1.Divide(2, 2)
    
    c1.cd(1)
    h_eff_pt = ROOT.TEfficiency(h_reco_d0_pt, h_mc_d0_pt)
    h_eff_pt.SetTitle("D0 Reconstruction Efficiency vs p_{T}")
    h_eff_pt.Draw()
    
    c1.cd(2)
    h_eff_eta = ROOT.TEfficiency(h_reco_d0_eta, h_mc_d0_eta)
    h_eff_eta.SetTitle("D0 Reconstruction Efficiency vs #eta")
    h_eff_eta.Draw()
    
    c1.cd(3)
    h_eff_z = ROOT.TEfficiency(h_reco_d0_z, h_mc_d0_z)
    h_eff_z.SetTitle("D0 Reconstruction Efficiency vs z")
    h_eff_z.Draw()
    
    c1.cd(4)
    h_origin.Draw()
    
    c1.SaveAs("d0_efficiency_plots.pdf")
    
    # Create resolution plots
    c2 = ROOT.TCanvas("c2", "D0 Resolution", 900, 450)
    c2.Divide(2, 1)
    
    c2.cd(1)
    h_d0_pt_res.Draw()
    
    c2.cd(2)
    h_d0_mass_res.Draw()
    
    c2.SaveAs("d0_resolution_plots.pdf")
    
    # Print summary statistics
    print("\nMC Analysis Summary:")
    print("  Total MC D0s: {}".format(h_mc_d0_pt.GetEntries()))
    print("  Reconstructed D0s: {}".format(h_reco_d0_pt.GetEntries()))
    print("  Overall efficiency: {:.1f}%".format(
        100.0 * h_reco_d0_pt.GetEntries() / h_mc_d0_pt.GetEntries() if h_mc_d0_pt.GetEntries() > 0 else 0))
    print("  Prompt D0s: {:.1f}%".format(
        100.0 * h_origin.GetBinContent(1) / h_origin.GetEntries() if h_origin.GetEntries() > 0 else 0))
    print("  D0s from B: {:.1f}%".format(
        100.0 * h_origin.GetBinContent(2) / h_origin.GetEntries() if h_origin.GetEntries() > 0 else 0))
    print("  pT resolution: {:.2f}%".format(100.0 * h_d0_pt_res.GetStdDev()))
    print("  Mass resolution: {:.2f}%".format(100.0 * h_d0_mass_res.GetStdDev()))
    
    root_file.Close()



 if __name__ == "__main__":

    # If using MC, perform additional MC truth validation
    if Type == 'MC':
        print("\nPerforming MC truth validation analysis...")
        analyze_mc_efficiency(CustomOutputFileName)