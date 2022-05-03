#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <set>
#include <kfcmd/hypos/Hypo2Ks_NoKsMasses.hpp>
#include "TrPh.h"

const double TrPh::dZ_ = 30;
const double TrPh::dRho_ = 30;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;
const std::set<std::string> TrPh::sKs1_ = {"pi+_1", "pi-_1"};
const std::set<std::string> TrPh::sKs2_ = {"pi+_2", "pi-_2"};

TrPh::TrPh(TTree *tree) :
  kfcmd::core::TrPh(tree) {
}

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

Int_t TrPh::Cut(Long64_t) {
  if (nt < 6) return -1;
  if (!cutTracks_()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(trackIndices_.begin(), trackIndices_.end(),
            [&charges](int i, int j) { return charges[i] < charges[j]; });
  return 1;
}

void TrPh::setupOutputBranches_(TTree *tree) {
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("in_mks1", &in_mks1_, "in_mks1/D");
  tree->Branch("in_mks2", &in_mks2_, "in_mks2/D");
  tree->Branch("kf_mks1", &kf_mks1_, "kf_mks1/D");
  tree->Branch("kf_mks2", &kf_mks2_, "kf_mks2/D");
  tree->Branch("vtx0", vtx0_, "vtx0[3]/D");
  tree->Branch("vtx1", vtx1_, "vtx1[3]/D");
  tree->Branch("vtx2", vtx2_, "vtx2[3]/D");
  tree->Branch("vtx1_dr", &vtx1_dr_, "vtx1_dr/D");
  tree->Branch("vtx2_dr", &vtx2_dr_, "vtx2_dr/D");
}

void TrPh::Loop(const std::string &outpath, double mfield) {
  if (fChain == 0)
    return;
  fChain->GetEntry(0);
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  auto out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  kfcmd::hypos::Hypo2Ks_NoKsMasses hypo(2.e-3 * emeas, mfield);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    if (jentry % 100 == 0)
      std::cout << jentry << " / " << nentries << std::endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0)
      continue;
    hypo.setBeamXY(xbeam, ybeam);
    hypo.fixVertexComponent("vtx0", xbeam, kfbase::core::VERTEX_X);
    hypo.fixVertexComponent("vtx0", ybeam, kfbase::core::VERTEX_Y);
    std::vector<int> mi_perm = {0, 1};
    kf_err_ = 1;
    kf_chi2_ = std::numeric_limits<double>::infinity();
    do {
      std::vector<int> pl_perm = {2, 3};
      do {
        if (!hypo.fillTrack("pi-_1", trackIndices_[mi_perm[0]], *this))
          continue;
        if (!hypo.fillTrack("pi-_2", trackIndices_[mi_perm[1]], *this))
          continue;
        if (!hypo.fillTrack("pi+_1", trackIndices_[pl_perm[0]], *this))
          continue;
        if (!hypo.fillTrack("pi+_2", trackIndices_[pl_perm[1]], *this))
          continue;
        Eigen::VectorXd tmpv(4);
        auto ks1P = hypo.getInitialMomentum("pi-_1") + hypo.getInitialMomentum("pi+_1");
        tmpv << ks1P.Px(), ks1P.Py(), ks1P.Pz(), 1.e-3;
        hypo.setInitialParticleParams("ks1", tmpv);
        auto ks2P = hypo.getInitialMomentum("pi-_2") +
          hypo.getInitialMomentum("pi+_2");
        tmpv << ks2P.Px(), ks2P.Py(), ks2P.Pz(), 1.e-3;
        hypo.setInitialParticleParams("ks2", tmpv);
        hypo.optimize();
        if (hypo.getErrorCode() != 0)
          continue;
        kf_err_ = 0;
        double tchi2 = hypo.getChiSquare();
        if (tchi2 < kf_chi2_) {
          kf_chi2_ = tchi2;
          in_mks1_ = hypo.getInitialMomentum(sKs1_).M();
          in_mks2_ = hypo.getInitialMomentum(sKs2_).M();
          kf_mks1_ = hypo.getFinalMomentum(sKs1_).M();
          kf_mks2_ = hypo.getFinalMomentum(sKs2_).M();
          auto vtx0 = hypo.getFinalVertex("vtx0");
          auto vtx1 = hypo.getFinalVertex("vtx1");
          auto vtx2 = hypo.getFinalVertex("vtx2");
          vtx0.GetXYZ(vtx0_);
          vtx1.GetXYZ(vtx1_);
          vtx2.GetXYZ(vtx2_);
          vtx1_dr_ = (vtx1 - vtx0).Mag();
          vtx2_dr_ = (vtx2 - vtx0).Mag();
        }
      } while (std::next_permutation(pl_perm.begin(), pl_perm.end()));
    } while (std::next_permutation(mi_perm.begin(), mi_perm.end()));
    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  outfl->Close();
  delete outfl;
}
