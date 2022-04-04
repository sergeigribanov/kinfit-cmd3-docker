#include <cmath>
#include <iostream>
#include <set>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <kfcmd/hypos/Hypo2ChPions2Photons.hpp>
#include "TrPh.h"

const double TrPh::dZ_ = 20;
const double TrPh::dRho_ = 1;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;
const double TrPh::minPhEn_ = 20;
const double TrPh::maxPhEn_ = 1000;

TrPh::TrPh(TTree* tree) : kfcmd::core::TrPh(tree) {}

TrPh::~TrPh() {}

bool TrPh::cutTracks_() {
  trackIndices_.clear();
  for (int i = 0; i < nt; i++) {
    bool point = (std::fabs(tz[i]) < dZ_) && (std::fabs(trho[i]) < dRho_);
    bool dedx = (tdedx[i] > mindEdX_) && (tdedx[i] < maxdEdX_);
    bool ptot = (tptot[i] > minTPtot_) && (tptot[i] < maxTPtot_);
    if (point && dedx && ptot) trackIndices_.push_back(i);
  }
  if (trackIndices_.size() == 2) {
    int totalCharge = 0;
    for (int i = 0; i < 2; ++i) totalCharge += tcharge[trackIndices_[i]];
    return (totalCharge == 0);
  }
  return false;
}

bool TrPh::cutPhotons_() {
  photonIndices_.clear();
  for (int i = 0; i < nph; i++)
    if ((phen[i] > minPhEn_) && (phen[i] < maxPhEn_))
      photonIndices_.push_back(i);
  if (photonIndices_.size() > 1) return true;
  return false;
}

Int_t TrPh::Cut(Long64_t entry) {
  bool bvar;
  bvar = (nt > 1) && (nph > 1);
  if (!bvar) return -1;
  if (!cutTracks_()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(trackIndices_.begin(), trackIndices_.end(),
     [&charges](int i, int j) { return charges[i] < charges[j]; });
  if (!cutPhotons_()) return -1;
  return 1;
}

void TrPh::setupOutptuBranches_(TTree* tree) {
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_mgg/D");
  tree->Branch("in_mgg", &in_mgg_, "in_mgg/D");
  tree->Branch("kf_mgg", &kf_mgg_, "kf_mgg/D");
  tree->Branch("in_mpipi", &in_mpipi_, "in_mpipi/D");
  tree->Branch("kf_mpipi", &kf_mpipi_, "kf_mpipi/D");
  tree->Branch("kf_vtx_x", &kf_vtx_x_, "kf_vtx_x/D");
  tree->Branch("kf_vtx_y", &kf_vtx_y_, "kf_vtx_y/D");
  tree->Branch("kf_vtx_z", &kf_vtx_z_, "kf_vtx_z/D");
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  const std::set<std::string> sGG = {"g0", "g1"};
  const std::set<std::string> sPiPi = {"pi+", "pi-"};
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutptuBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::Hypo2ChPions2Photons hypo(2 * emeas, magneticField);
  double tchi2;
  int errCode;
  in_mgg_ = 0;
  kf_mgg_ = 0;
  in_mpipi_ = 0;
  kf_mpipi_ = 0;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    hypo.setBeamXY(xbeam, ybeam);
    kf_err_ = 1;
    if (!hypo.fillTrack("pi-", trackIndices_[0], *this)) continue;
    if (!hypo.fillTrack("pi+", trackIndices_[1], *this)) continue;
    kf_chi2_ = std::numeric_limits<double>::infinity();
    for (std::size_t iph = 0; iph + 1 < photonIndices_.size(); ++iph) {
      if (!hypo.fillPhoton("g0", photonIndices_[iph], *this)) continue;
      for (std::size_t jph = iph + 1; jph < photonIndices_.size(); ++jph) {
        if (!hypo.fillPhoton("g1", photonIndices_[jph], *this)) continue;
        hypo.optimize();
        errCode = hypo.getErrorCode();
        if (errCode != 0) continue;
        tchi2 = hypo.getChiSquare();
        if (tchi2 < kf_chi2_) {
          kf_err_ = 0;
          kf_chi2_ = tchi2;
          in_mgg_ = hypo.getInitialMomentum(sGG).M();
          kf_mgg_ = hypo.getFinalMomentum(sGG).M();
          in_mpipi_ = hypo.getInitialMomentum(sPiPi).M();
          kf_mpipi_ = hypo.getFinalMomentum(sPiPi).M();
          kf_vtx_x_ = hypo.getFinalVertexX("vtx0");
          kf_vtx_y_ = hypo.getFinalVertexY("vtx0");
          kf_vtx_z_ = hypo.getFinalVertexZ("vtx0");
        }
      }
    }
    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  delete outfl;
}
