#include "TrackRecoEfficiency.h"
#include <iostream>
#include <vector>
#include <mutex>
#include <atomic>
#include <cstdlib>
#include <algorithm>
// #include <omp.h>

TrackRecoEfficiency::TrackRecoEfficiency(TString inputFileName, TString outputFileName)
    : m_inputFile(nullptr), m_outputFile(nullptr), m_tree(nullptr),
      m_minPt(1.0), m_minEta(2.0), m_maxEta(5.0), 
      m_maxGhostProb(0.3), m_maxTrackChi2(3.0),
      m_requireFromD0(false), m_useParallel(true), m_nThreads(4)
{
    // Open input file
    m_inputFile = TFile::Open(inputFileName, "READ");
    if (!m_inputFile || m_inputFile->IsZombie())
    {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return;
    }

    // Create output file
    m_outputFile = new TFile(outputFileName, "RECREATE");
    if (!m_outputFile || m_outputFile->IsZombie())
    {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        return;
    }

    // Get tree from input file
    m_tree = (TTree *)m_inputFile->Get("d0jets");
    if (!m_tree)
    {
        std::cerr << "Error: Could not find tree 'd0jets' in input file" << std::endl;
        return;
    }

    // Initialize vector pointers
    mc_particle_pid = nullptr;
    mc_particle_pt = nullptr;
    mc_particle_eta = nullptr;
    mc_particle_phi = nullptr;
    mc_particle_px = nullptr;
    mc_particle_py = nullptr;
    mc_particle_pz = nullptr;
    mc_particle_e = nullptr;
    mc_particle_matched = nullptr;
    mc_particle_fromD0 = nullptr;

    track_pt = nullptr;
    track_eta = nullptr;
    track_phi = nullptr;
    track_px = nullptr;
    track_py = nullptr;
    track_pz = nullptr;
    track_pid = nullptr;
    track_prb_ghost = nullptr;
    track_chi2 = nullptr;
    track_mc_matched = nullptr;

    // Set default binning
    // Default pT bins (GeV) - fine binning for tracks
    m_ptBins = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0};

    // Default eta bins - LHCb acceptance
    m_etaBins = {2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0};

    // Default p bins (GeV)
    m_pBins = {2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 100.0};
}

TrackRecoEfficiency::~TrackRecoEfficiency()
{
    CleanUp();
}

void TrackRecoEfficiency::CleanUp()
{
    // Clean up efficiency objects
    for (auto &pair : m_efficiencyMaps)
    {
        if (pair.second)
            delete pair.second;
    }
    m_efficiencyMaps.clear();

    for (auto &pair : m_efficiencyObjects)
    {
        if (pair.second)
            delete pair.second;
    }
    m_efficiencyObjects.clear();

    for (auto &pair : m_1DHistograms)
    {
        if (pair.second)
            delete pair.second;
    }
    m_1DHistograms.clear();

    // Close files
    if (m_outputFile)
    {
        m_outputFile->Close();
        delete m_outputFile;
        m_outputFile = nullptr;
    }

    if (m_inputFile)
    {
        m_inputFile->Close();
        m_inputFile = nullptr;
    }
}

bool TrackRecoEfficiency::Initialize()
{
    if (!m_inputFile || !m_outputFile || !m_tree)
    {
        std::cerr << "Error: Files or tree not properly initialized" << std::endl;
        return false;
    }

    InitializeBranches();
    CreateHistograms();

    std::cout << "TrackRecoEfficiency initialized successfully" << std::endl;
    std::cout << "Input file: " << m_inputFile->GetName() << std::endl;
    std::cout << "Output file: " << m_outputFile->GetName() << std::endl;
    std::cout << "Tree entries: " << m_tree->GetEntries() << std::endl;

    return true;
}

void TrackRecoEfficiency::InitializeBranches()
{
    // Set branch addresses for MC truth particles
    // These branch names may need to be adjusted based on your ROOT file structure
    if (m_tree->FindBranch("mc_particle_pid"))
        m_tree->SetBranchAddress("mc_particle_pid", &mc_particle_pid);
    else if (m_tree->FindBranch("mc_dau_pid"))
        m_tree->SetBranchAddress("mc_dau_pid", &mc_particle_pid);
    
    if (m_tree->FindBranch("mc_particle_pt"))
        m_tree->SetBranchAddress("mc_particle_pt", &mc_particle_pt);
    else if (m_tree->FindBranch("mc_dau_pt"))
        m_tree->SetBranchAddress("mc_dau_pt", &mc_particle_pt);
    
    if (m_tree->FindBranch("mc_particle_eta"))
        m_tree->SetBranchAddress("mc_particle_eta", &mc_particle_eta);
    else if (m_tree->FindBranch("mc_dau_eta"))
        m_tree->SetBranchAddress("mc_dau_eta", &mc_particle_eta);
    
    if (m_tree->FindBranch("mc_particle_phi"))
        m_tree->SetBranchAddress("mc_particle_phi", &mc_particle_phi);
    else if (m_tree->FindBranch("mc_dau_phi"))
        m_tree->SetBranchAddress("mc_dau_phi", &mc_particle_phi);
    
    if (m_tree->FindBranch("mc_particle_px"))
        m_tree->SetBranchAddress("mc_particle_px", &mc_particle_px);
    else if (m_tree->FindBranch("mc_dau_px"))
        m_tree->SetBranchAddress("mc_dau_px", &mc_particle_px);
    
    if (m_tree->FindBranch("mc_particle_py"))
        m_tree->SetBranchAddress("mc_particle_py", &mc_particle_py);
    else if (m_tree->FindBranch("mc_dau_py"))
        m_tree->SetBranchAddress("mc_dau_py", &mc_particle_py);
    
    if (m_tree->FindBranch("mc_particle_pz"))
        m_tree->SetBranchAddress("mc_particle_pz", &mc_particle_pz);
    else if (m_tree->FindBranch("mc_dau_pz"))
        m_tree->SetBranchAddress("mc_dau_pz", &mc_particle_pz);
    
    if (m_tree->FindBranch("mc_particle_matched"))
        m_tree->SetBranchAddress("mc_particle_matched", &mc_particle_matched);
    else if (m_tree->FindBranch("mc_dau_matched"))
        m_tree->SetBranchAddress("mc_dau_matched", &mc_particle_matched);
    
    // Set branch addresses for reconstructed tracks
    if (m_tree->FindBranch("n_tracks"))
        m_tree->SetBranchAddress("n_tracks", &n_tracks);
    else if (m_tree->FindBranch("n_daughters"))
        m_tree->SetBranchAddress("n_daughters", &n_tracks);
    
    if (m_tree->FindBranch("track_pt"))
        m_tree->SetBranchAddress("track_pt", &track_pt);
    else if (m_tree->FindBranch("dau_pt"))
        m_tree->SetBranchAddress("dau_pt", &track_pt);
    
    if (m_tree->FindBranch("track_eta"))
        m_tree->SetBranchAddress("track_eta", &track_eta);
    else if (m_tree->FindBranch("dau_eta"))
        m_tree->SetBranchAddress("dau_eta", &track_eta);
    
    if (m_tree->FindBranch("track_phi"))
        m_tree->SetBranchAddress("track_phi", &track_phi);
    else if (m_tree->FindBranch("dau_phi"))
        m_tree->SetBranchAddress("dau_phi", &track_phi);
    
    if (m_tree->FindBranch("track_px"))
        m_tree->SetBranchAddress("track_px", &track_px);
    else if (m_tree->FindBranch("dau_px"))
        m_tree->SetBranchAddress("dau_px", &track_px);
    
    if (m_tree->FindBranch("track_py"))
        m_tree->SetBranchAddress("track_py", &track_py);
    else if (m_tree->FindBranch("dau_py"))
        m_tree->SetBranchAddress("dau_py", &track_py);
    
    if (m_tree->FindBranch("track_pz"))
        m_tree->SetBranchAddress("track_pz", &track_pz);
    else if (m_tree->FindBranch("dau_pz"))
        m_tree->SetBranchAddress("dau_pz", &track_pz);
    
    if (m_tree->FindBranch("track_pid"))
        m_tree->SetBranchAddress("track_pid", &track_pid);
    else if (m_tree->FindBranch("dau_pid"))
        m_tree->SetBranchAddress("dau_pid", &track_pid);
    
    if (m_tree->FindBranch("track_prb_ghost"))
        m_tree->SetBranchAddress("track_prb_ghost", &track_prb_ghost);
    else if (m_tree->FindBranch("dau_prb_ghost"))
        m_tree->SetBranchAddress("dau_prb_ghost", &track_prb_ghost);
    
    if (m_tree->FindBranch("track_chi2"))
        m_tree->SetBranchAddress("track_chi2", &track_chi2);
    else if (m_tree->FindBranch("dau_trckChi2"))
        m_tree->SetBranchAddress("dau_trckChi2", &track_chi2);
}

void TrackRecoEfficiency::CreateHistograms()
{
    // Create efficiency histograms for pions (pT vs eta)
    m_efficiencyMaps["pion_numerator"] = new TH2F("pion_numerator", "Pion Track Efficiency Numerator;p_{T} [GeV];#eta",
                                                  m_ptBins.size() - 1, &m_ptBins[0],
                                                  m_etaBins.size() - 1, &m_etaBins[0]);
    m_efficiencyMaps["pion_denominator"] = new TH2F("pion_denominator", "Pion Track Efficiency Denominator;p_{T} [GeV];#eta",
                                                    m_ptBins.size() - 1, &m_ptBins[0],
                                                    m_etaBins.size() - 1, &m_etaBins[0]);

    // Create efficiency histograms for kaons (pT vs eta)
    m_efficiencyMaps["kaon_numerator"] = new TH2F("kaon_numerator", "Kaon Track Efficiency Numerator;p_{T} [GeV];#eta",
                                                  m_ptBins.size() - 1, &m_ptBins[0],
                                                  m_etaBins.size() - 1, &m_etaBins[0]);
    m_efficiencyMaps["kaon_denominator"] = new TH2F("kaon_denominator", "Kaon Track Efficiency Denominator;p_{T} [GeV];#eta",
                                                    m_ptBins.size() - 1, &m_ptBins[0],
                                                    m_etaBins.size() - 1, &m_etaBins[0]);

    // Create TEfficiency objects for proper error handling
    m_efficiencyObjects["pion_efficiency"] = new TEfficiency("pion_efficiency", "Pion Track Reconstruction Efficiency;p_{T} [GeV];#eta",
                                                             m_ptBins.size() - 1, &m_ptBins[0],
                                                             m_etaBins.size() - 1, &m_etaBins[0]);

    m_efficiencyObjects["kaon_efficiency"] = new TEfficiency("kaon_efficiency", "Kaon Track Reconstruction Efficiency;p_{T} [GeV];#eta",
                                                             m_ptBins.size() - 1, &m_ptBins[0],
                                                             m_etaBins.size() - 1, &m_etaBins[0]);

    // Create momentum-based histograms if p bins are defined
    if (!m_pBins.empty())
    {
        m_efficiencyMaps["pion_numerator_p"] = new TH2F("pion_numerator_p", "Pion Track Efficiency Numerator (p);p [GeV];#eta",
                                                        m_pBins.size() - 1, &m_pBins[0],
                                                        m_etaBins.size() - 1, &m_etaBins[0]);
        m_efficiencyMaps["pion_denominator_p"] = new TH2F("pion_denominator_p", "Pion Track Efficiency Denominator (p);p [GeV];#eta",
                                                          m_pBins.size() - 1, &m_pBins[0],
                                                          m_etaBins.size() - 1, &m_etaBins[0]);

        m_efficiencyMaps["kaon_numerator_p"] = new TH2F("kaon_numerator_p", "Kaon Track Efficiency Numerator (p);p [GeV];#eta",
                                                        m_pBins.size() - 1, &m_pBins[0],
                                                        m_etaBins.size() - 1, &m_etaBins[0]);
        m_efficiencyMaps["kaon_denominator_p"] = new TH2F("kaon_denominator_p", "Kaon Track Efficiency Denominator (p);p [GeV];#eta",
                                                          m_pBins.size() - 1, &m_pBins[0],
                                                          m_etaBins.size() - 1, &m_etaBins[0]);

        m_efficiencyObjects["pion_efficiency_p"] = new TEfficiency("pion_efficiency_p", "Pion Track Reconstruction Efficiency (p);p [GeV];#eta",
                                                                   m_pBins.size() - 1, &m_pBins[0],
                                                                   m_etaBins.size() - 1, &m_etaBins[0]);

        m_efficiencyObjects["kaon_efficiency_p"] = new TEfficiency("kaon_efficiency_p", "Kaon Track Reconstruction Efficiency (p);p [GeV];#eta",
                                                                   m_pBins.size() - 1, &m_pBins[0],
                                                                   m_etaBins.size() - 1, &m_etaBins[0]);
    }

    // Create 1D efficiency histograms vs pT
    m_1DHistograms["pion_eff_vs_pt"] = new TH1F("pion_eff_vs_pt", "Pion Track Efficiency vs p_{T};p_{T} [GeV];Efficiency",
                                                m_ptBins.size() - 1, &m_ptBins[0]);
    m_1DHistograms["kaon_eff_vs_pt"] = new TH1F("kaon_eff_vs_pt", "Kaon Track Efficiency vs p_{T};p_{T} [GeV];Efficiency",
                                                m_ptBins.size() - 1, &m_ptBins[0]);

    // Create 1D efficiency histograms vs eta
    m_1DHistograms["pion_eff_vs_eta"] = new TH1F("pion_eff_vs_eta", "Pion Track Efficiency vs #eta;#eta;Efficiency",
                                                 m_etaBins.size() - 1, &m_etaBins[0]);
    m_1DHistograms["kaon_eff_vs_eta"] = new TH1F("kaon_eff_vs_eta", "Kaon Track Efficiency vs #eta;#eta;Efficiency",
                                                 m_etaBins.size() - 1, &m_etaBins[0]);

    std::cout << "Created track efficiency histograms with:" << std::endl;
    std::cout << "  pT bins: " << m_ptBins.size() - 1 << " (" << m_ptBins.front() << " - " << m_ptBins.back() << " GeV)" << std::endl;
    std::cout << "  eta bins: " << m_etaBins.size() - 1 << " (" << m_etaBins.front() << " - " << m_etaBins.back() << ")" << std::endl;
    if (!m_pBins.empty())
    {
        std::cout << "  p bins: " << m_pBins.size() - 1 << " (" << m_pBins.front() << " - " << m_pBins.back() << " GeV)" << std::endl;
    }
}

void TrackRecoEfficiency::CalculateTrackEfficiency()
{
    std::cout << "\nCalculating track reconstruction efficiency..." << std::endl;

    // Get efficiency histograms
    TH2F *h_pion_num = m_efficiencyMaps["pion_numerator"];
    TH2F *h_pion_den = m_efficiencyMaps["pion_denominator"];
    TH2F *h_kaon_num = m_efficiencyMaps["kaon_numerator"];
    TH2F *h_kaon_den = m_efficiencyMaps["kaon_denominator"];

    TEfficiency *eff_pion = m_efficiencyObjects["pion_efficiency"];
    TEfficiency *eff_kaon = m_efficiencyObjects["kaon_efficiency"];

    // Get momentum-based histograms (optional)
    TH2F *h_pion_num_p = nullptr;
    TH2F *h_pion_den_p = nullptr;
    TH2F *h_kaon_num_p = nullptr;
    TH2F *h_kaon_den_p = nullptr;
    TEfficiency *eff_pion_p = nullptr;
    TEfficiency *eff_kaon_p = nullptr;

    if (!m_pBins.empty())
    {
        h_pion_num_p = m_efficiencyMaps["pion_numerator_p"];
        h_pion_den_p = m_efficiencyMaps["pion_denominator_p"];
        h_kaon_num_p = m_efficiencyMaps["kaon_numerator_p"];
        h_kaon_den_p = m_efficiencyMaps["kaon_denominator_p"];
        eff_pion_p = m_efficiencyObjects["pion_efficiency_p"];
        eff_kaon_p = m_efficiencyObjects["kaon_efficiency_p"];
    }

    if (!h_pion_num || !h_pion_den || !h_kaon_num || !h_kaon_den)
    {
        std::cerr << "Error: Efficiency histograms not available" << std::endl;
        return;
    }

    // Loop over all events
    Long64_t nEntries = m_tree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;

    // Counters for statistics
    long totalPionsInAcceptance = 0;
    long totalKaonsInAcceptance = 0;
    long reconstructedPions = 0;
    long reconstructedKaons = 0;

    for (Long64_t entry = 0; entry < nEntries; ++entry)
    {
        m_tree->GetEntry(entry);

        if (entry % 100000 == 0)
        {
            std::cout << "Processing entry " << entry << "/" << nEntries << " (" << (100.0 * entry / nEntries) << "%)" << std::endl;
        }

        if (!mc_particle_pid || !mc_particle_pt || !mc_particle_eta)
            continue;

        // Loop over all MC truth particles
        for (size_t mc_idx = 0; mc_idx < mc_particle_pid->size(); ++mc_idx)
        {
            int pid = mc_particle_pid->at(mc_idx);
            double pt = mc_particle_pt->at(mc_idx);
            double eta = mc_particle_eta->at(mc_idx);

            // Skip if not pion or kaon
            if (abs(pid) != 211 && abs(pid) != 321)
                continue;

            // Check if particle is in acceptance
            if (!IsInAcceptance(pt, eta))
                continue;

            // Calculate momentum if needed
            double p = 0.0;
            if (!m_pBins.empty() && mc_particle_px && mc_particle_py && mc_particle_pz)
            {
                p = CalculateMomentum(mc_particle_px->at(mc_idx), 
                                    mc_particle_py->at(mc_idx), 
                                    mc_particle_pz->at(mc_idx));
            }

            // Check if this is a pion or kaon
            bool isPion = (abs(pid) == 211);
            bool isKaon = (abs(pid) == 321);

            if (isPion)
            {
                totalPionsInAcceptance++;
                h_pion_den->Fill(pt, eta);
                if (h_pion_den_p) h_pion_den_p->Fill(p, eta);
            }
            else if (isKaon)
            {
                totalKaonsInAcceptance++;
                h_kaon_den->Fill(pt, eta);
                if (h_kaon_den_p) h_kaon_den_p->Fill(p, eta);
            }

            // Check if this MC particle has a reconstructed track match
            bool hasRecoMatch = false;
            if (mc_particle_matched && mc_idx < mc_particle_matched->size())
            {
                int matchedTrackIdx = mc_particle_matched->at(mc_idx);
                if (matchedTrackIdx >= 0 && matchedTrackIdx < n_tracks)
                {
                    // Check if the reconstructed track passes quality cuts
                    if (PassesQualityCuts(matchedTrackIdx))
                    {
                        hasRecoMatch = true;
                    }
                }
            }

            // Fill numerator and efficiency objects
            if (isPion)
            {
                if (hasRecoMatch)
                {
                    reconstructedPions++;
                    h_pion_num->Fill(pt, eta);
                    if (h_pion_num_p) h_pion_num_p->Fill(p, eta);
                }
                eff_pion->Fill(hasRecoMatch, pt, eta);
                if (eff_pion_p) eff_pion_p->Fill(hasRecoMatch, p, eta);
            }
            else if (isKaon)
            {
                if (hasRecoMatch)
                {
                    reconstructedKaons++;
                    h_kaon_num->Fill(pt, eta);
                    if (h_kaon_num_p) h_kaon_num_p->Fill(p, eta);
                }
                eff_kaon->Fill(hasRecoMatch, pt, eta);
                if (eff_kaon_p) eff_kaon_p->Fill(hasRecoMatch, p, eta);
            }
        }
    }

    // Create efficiency maps
    TH2F *h_pion_eff = (TH2F *)h_pion_num->Clone("pion_efficiency_map");
    h_pion_eff->Divide(h_pion_den);
    h_pion_eff->SetTitle("Pion Track Reconstruction Efficiency;p_{T} [GeV];#eta");
    m_efficiencyMaps["pion_efficiency_map"] = h_pion_eff;

    TH2F *h_kaon_eff = (TH2F *)h_kaon_num->Clone("kaon_efficiency_map");
    h_kaon_eff->Divide(h_kaon_den);
    h_kaon_eff->SetTitle("Kaon Track Reconstruction Efficiency;p_{T} [GeV];#eta");
    m_efficiencyMaps["kaon_efficiency_map"] = h_kaon_eff;

    // Create momentum-based efficiency maps if available
    if (h_pion_num_p && h_pion_den_p)
    {
        TH2F *h_pion_eff_p = (TH2F *)h_pion_num_p->Clone("pion_efficiency_map_p");
        h_pion_eff_p->Divide(h_pion_den_p);
        h_pion_eff_p->SetTitle("Pion Track Reconstruction Efficiency;p [GeV];#eta");
        m_efficiencyMaps["pion_efficiency_map_p"] = h_pion_eff_p;
    }

    if (h_kaon_num_p && h_kaon_den_p)
    {
        TH2F *h_kaon_eff_p = (TH2F *)h_kaon_num_p->Clone("kaon_efficiency_map_p");
        h_kaon_eff_p->Divide(h_kaon_den_p);
        h_kaon_eff_p->SetTitle("Kaon Track Reconstruction Efficiency;p [GeV];#eta");
        m_efficiencyMaps["kaon_efficiency_map_p"] = h_kaon_eff_p;
    }

    std::cout << "Track reconstruction efficiency calculation completed" << std::endl;

    // Print statistics
    double pionEff = (totalPionsInAcceptance > 0) ? (double)reconstructedPions / totalPionsInAcceptance : 0.0;
    double kaonEff = (totalKaonsInAcceptance > 0) ? (double)reconstructedKaons / totalKaonsInAcceptance : 0.0;

    std::cout << "\nStatistics:" << std::endl;
    std::cout << "  Total MC pions in acceptance: " << totalPionsInAcceptance << std::endl;
    std::cout << "  Reconstructed pions: " << reconstructedPions << std::endl;
    std::cout << "  Pion track efficiency: " << pionEff * 100.0 << "% (" << reconstructedPions << "/" << totalPionsInAcceptance << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "  Total MC kaons in acceptance: " << totalKaonsInAcceptance << std::endl;
    std::cout << "  Reconstructed kaons: " << reconstructedKaons << std::endl;
    std::cout << "  Kaon track efficiency: " << kaonEff * 100.0 << "% (" << reconstructedKaons << "/" << totalKaonsInAcceptance << ")" << std::endl;
}

bool TrackRecoEfficiency::IsInAcceptance(double pt, double eta)
{
    return (pt >= m_minPt && eta >= m_minEta && eta <= m_maxEta);
}

bool TrackRecoEfficiency::PassesQualityCuts(int track_idx)
{
    if (track_idx < 0)
        return false;

    // Check ghost probability if available
    if (track_prb_ghost && track_idx < (int)track_prb_ghost->size())
    {
        if (track_prb_ghost->at(track_idx) > m_maxGhostProb)
            return false;
    }

    // Check track chi2 if available
    if (track_chi2 && track_idx < (int)track_chi2->size())
    {
        if (track_chi2->at(track_idx) > m_maxTrackChi2)
            return false;
    }

    return true;
}

double TrackRecoEfficiency::CalculateMomentum(double px, double py, double pz)
{
    return TMath::Sqrt(px * px + py * py + pz * pz);
}

std::string TrackRecoEfficiency::GetParticleName(int pid)
{
    switch (abs(pid))
    {
    case 211:
        return "pion";
    case 321:
        return "kaon";
    default:
        return "unknown";
    }
}

void TrackRecoEfficiency::SaveResults()
{
    if (!m_outputFile)
    {
        std::cerr << "Error: Output file not available" << std::endl;
        return;
    }

    std::cout << "\nSaving results to " << m_outputFile->GetName() << "..." << std::endl;

    m_outputFile->cd();

    // Save all efficiency histograms
    for (auto &pair : m_efficiencyMaps)
    {
        if (pair.second)
        {
            pair.second->Write();
            std::cout << "  Saved histogram: " << pair.first << std::endl;
        }
    }

    // Save all efficiency objects
    for (auto &pair : m_efficiencyObjects)
    {
        if (pair.second)
        {
            pair.second->Write();
            std::cout << "  Saved efficiency object: " << pair.first << std::endl;
        }
    }

    // Save all 1D histograms
    for (auto &pair : m_1DHistograms)
    {
        if (pair.second)
        {
            pair.second->Write();
            std::cout << "  Saved 1D histogram: " << pair.first << std::endl;
        }
    }

    m_outputFile->Write();
    std::cout << "Results saved successfully" << std::endl;
}

void TrackRecoEfficiency::PlotEfficiency(const std::string &particleType)
{
    std::string histName = particleType + "_efficiency_map";
    TH2F *h_eff = m_efficiencyMaps[histName];
    
    if (!h_eff)
    {
        std::cerr << "Error: Efficiency histogram '" << histName << "' not found" << std::endl;
        return;
    }

    std::cout << "Plotting " << particleType << " track efficiency..." << std::endl;

    // Set up canvas
    TCanvas *c1 = new TCanvas("c_track_eff", (particleType + " Track Efficiency").c_str(), 800, 600);
    c1->SetRightMargin(0.15);

    // Set up histogram style
    h_eff->SetStats(0);
    h_eff->GetZaxis()->SetTitle("Track Efficiency");
    h_eff->GetZaxis()->SetRangeUser(0.0, 1.0);

    // Draw efficiency map
    h_eff->Draw("COLZ");

    // Add text with overall efficiency
    TH2F *h_num = m_efficiencyMaps[particleType + "_numerator"];
    TH2F *h_den = m_efficiencyMaps[particleType + "_denominator"];
    if (h_num && h_den)
    {
        double overallEff = (double)h_num->GetEntries() / h_den->GetEntries();
        TString effText = TString::Format("Overall %s Track Efficiency: %.3f", 
                                         particleType.c_str(), overallEff);

        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.15, 0.92, effText);
    }

    // Save canvas
    TString plotName = TString::Format("%s_track_efficiency.png", particleType.c_str());
    c1->SaveAs(plotName);
    std::cout << "Saved plot: " << plotName << std::endl;

    // Also save as PDF
    plotName = TString::Format("%s_track_efficiency.pdf", particleType.c_str());
    c1->SaveAs(plotName);

    delete c1;
}

// Main function for running track efficiency calculation
int TrackRecoEfficiencyRun(TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/51/51.root", 
                          TString outputFile = "track_efficiency_output.root",
                          double minPt = 0.5, double minEta = 2.0, double maxEta = 5.0,
                          double maxGhostProb = 0.3, double maxChi2 = 3.0,
                          bool makePlots = true)
{
    std::cout << "=== Track Reconstruction Efficiency Calculator ===" << std::endl;
    std::cout << "Input file: " << inputFile << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    std::cout << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Minimum pT: " << minPt << " GeV" << std::endl;
    std::cout << "  Eta range: " << minEta << " - " << maxEta << std::endl;
    std::cout << "  Max ghost probability: " << maxGhostProb << std::endl;
    std::cout << "  Max track chi2: " << maxChi2 << std::endl;
    std::cout << "  Generate plots: " << (makePlots ? "Yes" : "No") << std::endl;
    std::cout << std::endl;

    try
    {
        // Create track efficiency calculator
        TrackRecoEfficiency calculator(inputFile, outputFile);

        // Configure parameters
        calculator.SetPtRange(minPt);
        calculator.SetEtaRange(minEta, maxEta);
        calculator.SetQualityCuts(maxGhostProb, maxChi2);

        // Set up fine binning for track studies
        std::vector<double> customPtBins;
        for (double pt = 0.5; pt < 2.0; pt += 0.1) customPtBins.push_back(pt);
        for (double pt = 2.0; pt < 5.0; pt += 0.25) customPtBins.push_back(pt);
        for (double pt = 5.0; pt < 10.0; pt += 0.5) customPtBins.push_back(pt);
        for (double pt = 10.0; pt < 20.0; pt += 1.0) customPtBins.push_back(pt);
        for (double pt = 20.0; pt <= 30.0; pt += 2.5) customPtBins.push_back(pt);

        std::vector<double> customEtaBins;
        for (double eta = 2.0; eta <= 5.0; eta += 0.1) customEtaBins.push_back(eta);

        std::vector<double> customPBins;
        for (double p = 2.0; p < 10.0; p += 1.0) customPBins.push_back(p);
        for (double p = 10.0; p < 30.0; p += 2.5) customPBins.push_back(p);
        for (double p = 30.0; p <= 100.0; p += 10.0) customPBins.push_back(p);

        calculator.SetPtBins(customPtBins);
        calculator.SetEtaBins(customEtaBins);
        calculator.SetPBins(customPBins);

        // Initialize
        if (!calculator.Initialize())
        {
            std::cerr << "Error: Failed to initialize calculator" << std::endl;
            return 1;
        }

        // Calculate track efficiency
        std::cout << "Starting track efficiency calculation..." << std::endl;
        calculator.CalculateTrackEfficiency();

        // Save results
        calculator.SaveResults();

        // Generate plots if requested
        if (makePlots)
        {
            std::cout << "\nGenerating efficiency plots..." << std::endl;
            calculator.PlotEfficiency("pion");
            calculator.PlotEfficiency("kaon");
        }

        std::cout << "\n=== Track efficiency calculation completed successfully! ===" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

// Main function
int main(int argc, char* argv[])
{
    // Default parameters
    TString inputFile = "/media/niviths/SSD2/lhcb_analysis_SSD/20250708_newMC_fixedTrueAssociation/51/51.root";
    TString outputFile = "track_efficiency_output.root";
    double minPt = 0.5;
    double minEta = 2.0;
    double maxEta = 5.0;
    double maxGhostProb = 0.3;
    double maxChi2 = 3.0;
    bool makePlots = true;

    // Parse command line arguments
    if (argc > 1) inputFile = argv[1];
    if (argc > 2) outputFile = argv[2];
    if (argc > 3) minPt = atof(argv[3]);
    if (argc > 4) minEta = atof(argv[4]);
    if (argc > 5) maxEta = atof(argv[5]);
    if (argc > 6) maxGhostProb = atof(argv[6]);
    if (argc > 7) maxChi2 = atof(argv[7]);
    if (argc > 8) makePlots = (atoi(argv[8]) != 0);

    // Run the analysis
    return TrackRecoEfficiencyRun(inputFile, outputFile, minPt, minEta, maxEta, 
                                 maxGhostProb, maxChi2, makePlots);
}
