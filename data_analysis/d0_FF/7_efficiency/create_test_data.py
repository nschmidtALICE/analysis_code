#!/usr/bin/env python3

"""
Simple test script to create a minimal ROOT file for testing the 
D0 reconstruction efficiency calculator.

This creates a dummy ROOT file with the expected tree structure
but with minimal synthetic data for testing purposes.
"""

import ROOT
import numpy as np
import random

def create_test_root_file(filename="test_mc_d0.root", n_events=100):
    """Create a test ROOT file with the expected tree structure"""
    
    print(f"Creating test MC file: {filename}")
    
    # Create ROOT file and tree
    root_file = ROOT.TFile(filename, "RECREATE")
    tree = ROOT.TTree("d0jets", "Test D0 MC data")
    
    # Define variables for the tree - use simple arrays instead of std::vector
    from array import array
    
    n_d0s = array('i', [0])
    
    # We'll use a maximum size for arrays
    MAX_D0S = 10
    MAX_DAUGHTERS = 20
    MAX_MC_D0S = 10
    
    d0_pt = ROOT.std.vector('float')()
    d0_eta = ROOT.std.vector('float')()
    d0_phi = ROOT.std.vector('float')()
    d0_mass = ROOT.std.vector('float')()
    d0_px = ROOT.std.vector('float')()
    d0_py = ROOT.std.vector('float')()
    d0_pz = ROOT.std.vector('float')()
    d0_e = ROOT.std.vector('float')()
    
    n_daughters = array('i', [0])
    dau_pid = ROOT.std.vector('int')()
    dau_pt = ROOT.std.vector('float')()
    dau_eta = ROOT.std.vector('float')()
    dau_phi = ROOT.std.vector('float')()
    dau_px = ROOT.std.vector('float')()
    dau_py = ROOT.std.vector('float')()
    dau_pz = ROOT.std.vector('float')()
    dau_d0_idx = ROOT.std.vector('int')()
    dau_pnn_k = ROOT.std.vector('float')()
    dau_pnn_pi = ROOT.std.vector('float')()
    dau_prb_ghost = ROOT.std.vector('float')()
    dau_trckChi2 = ROOT.std.vector('float')()
    
    # MC variables
    mc_d0_pt = ROOT.std.vector('float')()
    mc_d0_eta = ROOT.std.vector('float')()
    mc_d0_phi = ROOT.std.vector('float')()
    mc_d0_px = ROOT.std.vector('float')()
    mc_d0_py = ROOT.std.vector('float')()
    mc_d0_pz = ROOT.std.vector('float')()
    mc_d0_matched = ROOT.std.vector('int')()
    
    # Set up tree branches
    tree.Branch("n_d0s", n_d0s, "n_d0s/I")
    tree.Branch("d0_pt", d0_pt)
    tree.Branch("d0_eta", d0_eta)
    tree.Branch("d0_phi", d0_phi)
    tree.Branch("d0_mass", d0_mass)
    tree.Branch("d0_px", d0_px)
    tree.Branch("d0_py", d0_py)
    tree.Branch("d0_pz", d0_pz)
    tree.Branch("d0_e", d0_e)
    
    tree.Branch("n_daughters", n_daughters, "n_daughters/I")
    tree.Branch("dau_pid", dau_pid)
    tree.Branch("dau_pt", dau_pt)
    tree.Branch("dau_eta", dau_eta)
    tree.Branch("dau_phi", dau_phi)
    tree.Branch("dau_px", dau_px)
    tree.Branch("dau_py", dau_py)
    tree.Branch("dau_pz", dau_pz)
    tree.Branch("dau_d0_idx", dau_d0_idx)
    tree.Branch("dau_pnn_k", dau_pnn_k)
    tree.Branch("dau_pnn_pi", dau_pnn_pi)
    tree.Branch("dau_prb_ghost", dau_prb_ghost)
    tree.Branch("dau_trckChi2", dau_trckChi2)
    
    tree.Branch("mc_d0_pt", mc_d0_pt)
    tree.Branch("mc_d0_eta", mc_d0_eta)
    tree.Branch("mc_d0_phi", mc_d0_phi)
    tree.Branch("mc_d0_px", mc_d0_px)
    tree.Branch("mc_d0_py", mc_d0_py)
    tree.Branch("mc_d0_pz", mc_d0_pz)
    tree.Branch("mc_d0_matched", mc_d0_matched)
    
    # Generate events
    for event in range(n_events):
        if event % 10 == 0:
            print(f"Generating event {event}/{n_events}")
        
        # Clear vectors
        n_d0s[0] = 0
        d0_pt.clear()
        d0_eta.clear()
        d0_phi.clear()
        d0_mass.clear()
        d0_px.clear()
        d0_py.clear()
        d0_pz.clear()
        d0_e.clear()
        
        n_daughters[0] = 0
        dau_pid.clear()
        dau_pt.clear()
        dau_eta.clear()
        dau_phi.clear()
        dau_px.clear()
        dau_py.clear()
        dau_pz.clear()
        dau_d0_idx.clear()
        dau_pnn_k.clear()
        dau_pnn_pi.clear()
        dau_prb_ghost.clear()
        dau_trckChi2.clear()
        
        mc_d0_pt.clear()
        mc_d0_eta.clear()
        mc_d0_phi.clear()
        mc_d0_px.clear()
        mc_d0_py.clear()
        mc_d0_pz.clear()
        mc_d0_matched.clear()
        
        # Generate MC D0s (1-3 per event)
        n_mc_d0s = random.randint(1, 3)
        n_reco_d0s = 0
        
        for mc_idx in range(n_mc_d0s):
            # Generate MC D0 kinematics
            pt = random.uniform(2.5, 15.0)
            eta = random.uniform(2.0, 5.0)
            phi = random.uniform(-np.pi, np.pi)
            
            # Calculate momentum components
            px = pt * np.cos(phi)
            py = pt * np.sin(phi)
            pz = pt * np.sinh(eta)
            p = np.sqrt(px*px + py*py + pz*pz)
            
            mc_d0_pt.push_back(pt)
            mc_d0_eta.push_back(eta)
            mc_d0_phi.push_back(phi)
            mc_d0_px.push_back(px)
            mc_d0_py.push_back(py)
            mc_d0_pz.push_back(pz)
            
            # Simulate reconstruction efficiency (80% chance to be reconstructed)
            is_reconstructed = random.random() < 0.8
            
            if is_reconstructed:
                # Create reconstructed D0 with some smearing
                reco_pt = pt + random.gauss(0, 0.1)
                reco_eta = eta + random.gauss(0, 0.02)
                reco_phi = phi + random.gauss(0, 0.02)
                reco_mass = 1.86484 + random.gauss(0, 0.015)  # D0 mass with resolution
                
                reco_px = reco_pt * np.cos(reco_phi)
                reco_py = reco_pt * np.sin(reco_phi)
                reco_pz = reco_pt * np.sinh(reco_eta)
                reco_e = np.sqrt(reco_px*reco_px + reco_py*reco_py + reco_pz*reco_pz + reco_mass*reco_mass)
                
                d0_pt.push_back(reco_pt)
                d0_eta.push_back(reco_eta)
                d0_phi.push_back(reco_phi)
                d0_mass.push_back(reco_mass)
                d0_px.push_back(reco_px)
                d0_py.push_back(reco_py)
                d0_pz.push_back(reco_pz)
                d0_e.push_back(reco_e)
                
                # Link MC to reco
                mc_d0_matched.push_back(n_reco_d0s)
                
                # Create daughters for this D0
                # Kaon daughter
                k_pt = random.uniform(1.0, pt/2)
                k_eta = eta + random.gauss(0, 0.1)
                k_phi = phi + random.gauss(0, 0.1)
                k_px = k_pt * np.cos(k_phi)
                k_py = k_pt * np.sin(k_phi)
                k_pz = k_pt * np.sinh(k_eta)
                
                dau_pid.push_back(321)  # K+ PID
                dau_pt.push_back(k_pt)
                dau_eta.push_back(k_eta)
                dau_phi.push_back(k_phi)
                dau_px.push_back(k_px)
                dau_py.push_back(k_py)
                dau_pz.push_back(k_pz)
                dau_d0_idx.push_back(n_reco_d0s)
                dau_pnn_k.push_back(random.uniform(0.6, 0.95))  # Good kaon PID
                dau_pnn_pi.push_back(random.uniform(0.05, 0.3)) # Poor pion PID for kaon
                dau_prb_ghost.push_back(random.uniform(0.0, 0.2))
                dau_trckChi2.push_back(random.uniform(0.5, 2.5))
                
                # Pion daughter
                pi_pt = random.uniform(1.0, pt/2)
                pi_eta = eta + random.gauss(0, 0.1)
                pi_phi = phi + random.gauss(0, 0.1)
                pi_px = pi_pt * np.cos(pi_phi)
                pi_py = pi_pt * np.sin(pi_phi)
                pi_pz = pi_pt * np.sinh(pi_eta)
                
                dau_pid.push_back(-211)  # π- PID
                dau_pt.push_back(pi_pt)
                dau_eta.push_back(pi_eta)
                dau_phi.push_back(pi_phi)
                dau_px.push_back(pi_px)
                dau_py.push_back(pi_py)
                dau_pz.push_back(pi_pz)
                dau_d0_idx.push_back(n_reco_d0s)
                dau_pnn_k.push_back(random.uniform(0.05, 0.3))   # Poor kaon PID for pion
                dau_pnn_pi.push_back(random.uniform(0.6, 0.95))  # Good pion PID
                dau_prb_ghost.push_back(random.uniform(0.0, 0.2))
                dau_trckChi2.push_back(random.uniform(0.5, 2.5))
                
                n_reco_d0s += 1
            else:
                # No reconstruction match
                mc_d0_matched.push_back(-1)
        
        # Set counters
        n_d0s[0] = n_reco_d0s
        n_daughters[0] = len(dau_pid)
        
        # Fill the tree
        tree.Fill()
    
    # Write and close
    root_file.Write()
    root_file.Close()
    
    print(f"Test file created: {filename}")
    print(f"Generated {n_events} events with MC D0 data")
    return filename

if __name__ == "__main__":
    import sys
    
    filename = "test_mc_d0.root"
    n_events = 1000
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    if len(sys.argv) > 2:
        n_events = int(sys.argv[2])
    
    create_test_root_file(filename, n_events)
    print(f"\nYou can now test the efficiency calculator with:")
    print(f"./calculate_d0_reco_efficiency {filename} test_output.root --plot")
