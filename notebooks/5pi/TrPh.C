#include <cmath>
#include <iostream>
#include <kfcmd/hypos/Hypo4ChPions2Photons.hpp>
#include "TrPh.h"

const double TrPh::dZ_ = 20;
const double TrPh::dRho_ = 1;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;
const double TrPh::minPhEn_ = 20;
const double TrPh::maxPhEn_ = 1000;

const std::set<std::string> TrPh::sM3Pi_[4] = {{"pi-_0", "pi+_0", "g0", "g1"},
                                               {"pi-_0", "pi+_1", "g0", "g1"},
                                               {"pi-_1", "pi+_0", "g0", "g1"},
                                               {"pi-_1", "pi+_1", "g0", "g1"}};

const std::set<std::string> TrPh::sgg_ = {"g0", "g1"};

TrPh::TrPh(TTree *tree)
    : kfcmd::core::TrPh(tree) {}

TrPh::~TrPh() {}

bool TrPh::cutTracks_() {
  trackIndices_.clear();
  for (int i = 0; i < nt; i++) {
    bool point = (std::fabs(tz[i]) < dZ_) && (std::fabs(trho[i]) < dRho_);
    bool dedx = (tdedx[i] > mindEdX_) && (tdedx[i] < maxdEdX_);
    bool ptot = (tptot[i] > minTPtot_) && (tptot[i] < maxTPtot_);
    if (point && dedx && ptot) trackIndices_.push_back(i);
  }
  if (trackIndices_.size() == 4) {
    int totalCharge = 0;
    for (int i = 0; i < 4; ++i) totalCharge += tcharge[trackIndices_[i]];
    return (totalCharge == 0);
  }
  return false;
}

bool TrPh::cutPhotons_() {
  photonIndices_.clear();
  for (int i = 0; i < nph; i++)
    if ((phen[i] > minPhEn_) && (phen[i] < maxPhEn_))
      photonIndices_.push_back(i);
  if (photonIndices_.size() > 1)
    return true;
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

void TrPh::setupOutputBranches_(TTree* tree) {
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("in_m3pi", in_m3pi_, "in_m3pi[4]/D");
  tree->Branch("kf_m3pi", kf_m3pi_, "kf_m3pi[4]/D");
  tree->Branch("in_mgg", &in_mgg_, "in_mgg/D");
  tree->Branch("kf_mgg", &kf_mgg_, "kf_mgg/D");
  tree->Branch("vtx", vtx_, "vtx[3]/D");
  tree->Branch("kf_err_mpi0", &kf_err_mpi0_, "kf_err_mpi0/I");
  tree->Branch("kf_chi2_mpi0", &kf_chi2_mpi0_, "kf_chi2_mpi0/D");
  tree->Branch("in_m3pi_mpi0", in_m3pi_mpi0_, "in_m3pi_mpi0[4]/D");
  tree->Branch("kf_m3pi_mpi0", kf_m3pi_mpi0_, "kf_m3pi_mpi0[4]/D");
  tree->Branch("vtx_mpi0", vtx_mpi0_, "vtx_mpi0[3]/D");
}

void TrPh::fit_(kfcmd::hypos::Hypo4ChPions2Photons *hypo) {
    hypo->setBeamXY(xbeam, ybeam);
    hypo->fixVertexComponent("vtx0", xbeam, kfbase::core::VERTEX_X);
    hypo->fixVertexComponent("vtx0", ybeam, kfbase::core::VERTEX_Y);
    kf_err_ = 1;
    kf_chi2_ = std::numeric_limits<double>::infinity();
    std::vector<int> mi_perm = {0, 1};
    do {
      std::vector<int> pl_perm = {2, 3};
      do {
        if (!hypo->fillTrack("pi-_0", trackIndices_[mi_perm[0]], *this))
          continue;
        if (!hypo->fillTrack("pi-_1", trackIndices_[mi_perm[1]], *this))
          continue;
        if (!hypo->fillTrack("pi+_0", trackIndices_[pl_perm[0]], *this))
          continue;
        if (!hypo->fillTrack("pi+_1", trackIndices_[pl_perm[1]], *this))
          continue;
        for (std::size_t iph = 0; iph + 1 < photonIndices_.size(); ++iph) {
          if (!hypo->fillPhoton("g0", photonIndices_[iph], *this))
            continue;
          for (std::size_t jph = iph + 1; jph < photonIndices_.size(); ++jph) {
            if (!hypo->fillPhoton("g1", photonIndices_[jph], *this))
              continue;
            hypo->optimize();
            if (hypo->getErrorCode() != 0)
              continue;
            kf_err_ = 0;
            double tchi2 = hypo->getChiSquare();
            if (tchi2 < kf_chi2_) {
              kf_chi2_ = tchi2;
              in_mgg_ = hypo->getInitialMomentum(sgg_).M();
              kf_mgg_ = hypo->getFinalMomentum(sgg_).M();
              for (int k = 0; k < 4; ++k) {
                in_m3pi_[k] = hypo->getInitialMomentum(sM3Pi_[k]).M();
                kf_m3pi_[k] = hypo->getFinalMomentum(sM3Pi_[k]).M();
              }
              auto vtx = hypo->getFinalVertex("vtx0");
              vtx.GetXYZ(vtx_);
            }
          }
        }
      }
      while (std::next_permutation(pl_perm.begin(), pl_perm.end()));
    }
    while (std::next_permutation(mi_perm.begin(), mi_perm.end()));
}

void TrPh::fit_mpi0_(kfcmd::hypos::Hypo4ChPions2Photons *hypo) {
    hypo->setBeamXY(xbeam, ybeam);
    hypo->fixVertexComponent("vtx0", xbeam, kfbase::core::VERTEX_X);
    hypo->fixVertexComponent("vtx0", ybeam, kfbase::core::VERTEX_Y);
    kf_err_mpi0_ = 1;
    kf_chi2_mpi0_ = std::numeric_limits<double>::infinity();
    std::vector<int> mi_perm = {0, 1};
    do {
      std::vector<int> pl_perm = {2, 3};
      do {
        if (!hypo->fillTrack("pi-_0", trackIndices_[mi_perm[0]], *this))
          continue;
        if (!hypo->fillTrack("pi-_1", trackIndices_[mi_perm[1]], *this))
          continue;
        if (!hypo->fillTrack("pi+_0", trackIndices_[pl_perm[0]], *this))
          continue;
        if (!hypo->fillTrack("pi+_1", trackIndices_[pl_perm[1]], *this))
          continue;
        for (std::size_t iph = 0; iph + 1 < photonIndices_.size(); ++iph) {
          if (!hypo->fillPhoton("g0", photonIndices_[iph], *this))
            continue;
          for (std::size_t jph = iph + 1; jph < photonIndices_.size(); ++jph) {
            if (!hypo->fillPhoton("g1", photonIndices_[jph], *this))
              continue;
            hypo->optimize();
            if (hypo->getErrorCode() != 0)
              continue;
            kf_err_mpi0_ = 0;
            double tchi2 = hypo->getChiSquare();
            if (tchi2 < kf_chi2_mpi0_) {
              kf_chi2_mpi0_ = tchi2;
              for (int k = 0; k < 4; ++k) {
                in_m3pi_mpi0_[k] = hypo->getInitialMomentum(sM3Pi_[k]).M();
                kf_m3pi_mpi0_[k] = hypo->getFinalMomentum(sM3Pi_[k]).M();
              }
              auto vtx = hypo->getFinalVertex("vtx0");
              vtx.GetXYZ(vtx_mpi0_);
            }
          }
        }
      }
      while (std::next_permutation(pl_perm.begin(), pl_perm.end()));
    }
    while (std::next_permutation(mi_perm.begin(), mi_perm.end()));
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::Hypo4ChPions2Photons hypo(2.e-3 * emeas, magneticField);
  kfcmd::hypos::Hypo4ChPions2Photons hypo_mpi0(2.e-3 * emeas, magneticField);
  hypo_mpi0.enableConstraint("m-pi0-constraint");
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    fit_(&hypo);
    fit_mpi0_(&hypo_mpi0);
    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  delete outfl;
}
