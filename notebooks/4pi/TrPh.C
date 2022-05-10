#include <cmath>
#include <iostream>
#include <set>
#include <kfcmd/hypos/Hypo4ChPions.hpp>
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
  if (trackIndices_.size() == 4) {
    int totalCharge = 0;
    for (int i = 0; i < 4; ++i) totalCharge += tcharge[trackIndices_[i]];
    return (totalCharge == 0);
  }
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
  return 1;
}

void TrPh::setupOutputBranches_(TTree* tree) {
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("in_p", in_p_, "in_p[4][3]/D");
  tree->Branch("in_e", in_e_, "in_e[4]/D");
  tree->Branch("kf_p", kf_p_, "kf_p[4][3]/D");
  tree->Branch("kf_e", kf_e_, "kf_e[4]/D");
  tree->Branch("vtx", vtx_, "vtx[3]/D");
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::Hypo4ChPions hypo(2.e-3 * emeas, magneticField);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    hypo.setBeamXY(xbeam, ybeam);
    hypo.fixVertexParameter("vtx0", 0, xbeam);
    hypo.fixVertexParameter("vtx0", 1, ybeam);
    kf_err_ = 1;
    kf_chi2_ = std::numeric_limits<double>::infinity();
    std::vector<int> mi_perm = {0, 1};
    do {
      std::vector<int> pl_perm = {2, 3};
      do {
        if (!hypo.fillTrack("pi-_0", trackIndices_[mi_perm[0]], *this))
          continue;
        if (!hypo.fillTrack("pi-_1", trackIndices_[mi_perm[1]], *this))
          continue;
        if (!hypo.fillTrack("pi+_0", trackIndices_[pl_perm[0]], *this))
          continue;
        if (!hypo.fillTrack("pi+_1", trackIndices_[pl_perm[1]], *this))
          continue;
        hypo.optimize();
        if (hypo.getErrorCode() != 0)  continue;
        kf_err_ = 0;
        double tchi2 = hypo.getChiSquare();
        if (tchi2 < kf_chi2_) {
          kf_chi2_ = tchi2;
          auto in_pimi0 = hypo.getInitialMomentum("pi-_0");
          auto in_pimi1 = hypo.getInitialMomentum("pi-_1");
          auto in_pipl0 = hypo.getInitialMomentum("pi+_0");
          auto in_pipl1 = hypo.getInitialMomentum("pi+_1");
          auto kf_pimi0 = hypo.getFinalMomentum("pi-_0");
          auto kf_pimi1 = hypo.getFinalMomentum("pi-_1");
          auto kf_pipl0 = hypo.getFinalMomentum("pi+_0");
          auto kf_pipl1 = hypo.getFinalMomentum("pi+_1");
          auto vtx = hypo.getFinalVertex("vtx0");
          in_p_[0][0] = in_pimi0.Px();
          in_p_[0][1] = in_pimi0.Py();
          in_p_[0][2] = in_pimi0.Pz();
          in_p_[1][0] = in_pimi1.Px();
          in_p_[1][1] = in_pimi1.Py();
          in_p_[1][2] = in_pimi1.Pz();
          in_p_[2][0] = in_pipl0.Px();
          in_p_[2][1] = in_pipl0.Py();
          in_p_[2][2] = in_pipl0.Pz();
          in_p_[3][0] = in_pipl1.Px();
          in_p_[3][1] = in_pipl1.Py();
          in_p_[3][2] = in_pipl1.Pz();
          kf_p_[0][0] = kf_pimi0.Px();
          kf_p_[0][1] = kf_pimi0.Py();
          kf_p_[0][2] = kf_pimi0.Pz();
          kf_p_[1][0] = kf_pimi1.Px();
          kf_p_[1][1] = kf_pimi1.Py();
          kf_p_[1][2] = kf_pimi1.Pz();
          kf_p_[2][0] = kf_pipl0.Px();
          kf_p_[2][1] = kf_pipl0.Py();
          kf_p_[2][2] = kf_pipl0.Pz();
          kf_p_[3][0] = kf_pipl1.Px();
          kf_p_[3][1] = kf_pipl1.Py();
          kf_p_[3][2] = kf_pipl1.Pz();
          in_e_[0] = in_pimi0.E();
          in_e_[1] = in_pimi1.E();
          in_e_[2] = in_pipl0.E();
          in_e_[3] = in_pipl1.E();
          kf_e_[0] = kf_pimi0.E();
          kf_e_[1] = kf_pimi1.E();
          kf_e_[2] = kf_pipl0.E();
          kf_e_[3] = kf_pipl1.E();
          vtx.GetXYZ(vtx_);
        }
      } while (std::next_permutation(pl_perm.begin(), pl_perm.end()));
    } while (std::next_permutation(mi_perm.begin(), mi_perm.end()));
      out_tree->Fill();
    }
  outfl->cd();
  out_tree->Write();
  delete outfl;
}
