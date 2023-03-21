#include <cmath>
#include <iostream>
#include "TrPh.h"

const double TrPh::dZ_ = 20;
const double TrPh::dRho_ = 1;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;

const std::set<std::string> TrPh::sM3Pi_ = {"pi-_0", "pi+_0", "pi0"};

const std::set<std::string> TrPh::spich_ = {"pi-_0", "pi-_1", "pi+_0", "pi+_1"};

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
  tree->Branch("in_m3pi", &in_m3pi_, "in_m3pi/D");
  tree->Branch("kf_m3pi", &kf_m3pi_, "kf_m3pi/D");
  tree->Branch("kf_vtx_m3pi", &kf_vtx_m3pi_, "kf_vtx_m3pi/D");
  tree->Branch("vtx", vtx_, "vtx[3]/D");
}

void TrPh::fit_(kfcmd::hypos::Hypo4ChPionsLostPi0_EtaMass *hypo,
		kfcmd::hypos::Hypo4PiVertex* hypo_vertex,
		kfcmd::hypos::Hypo4ChPionsLostPi0_EtaMass *hypo_tmp) {
  hypo->setBeamXY(xbeam, ybeam);
  // hypo->fixVertexParameter("vtx0", 0, xbeam);
  // hypo->fixVertexParameter("vtx0", 1, ybeam);

  hypo_vertex->setBeamXY(xbeam, ybeam);
  hypo_vertex->fixVertexParameter("vtx0", 0, xbeam);
  hypo_vertex->fixVertexParameter("vtx0", 1, ybeam);

  hypo_tmp->setBeamXY(xbeam, ybeam);
  hypo_tmp->fixVertexParameter("vtx0", 0, xbeam);
  hypo_tmp->fixVertexParameter("vtx0", 1, ybeam);

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

      hypo_vertex->fillTrack("pi-_0", trackIndices_[mi_perm[0]], *this);
      hypo_vertex->fillTrack("pi-_1", trackIndices_[mi_perm[1]], *this);
      hypo_vertex->fillTrack("pi+_0", trackIndices_[pl_perm[0]], *this);
      hypo_vertex->fillTrack("pi+_1", trackIndices_[pl_perm[1]], *this);
      hypo_vertex->optimize();

      hypo_tmp->fixVertexParameter("vtx0", 2, hypo_vertex->getFinalVertex("vtx0").Z());
      hypo_tmp->fillTrack("pi-_0", trackIndices_[mi_perm[0]], *this);
      hypo_tmp->fillTrack("pi-_1", trackIndices_[mi_perm[1]], *this);
      hypo_tmp->fillTrack("pi+_0", trackIndices_[pl_perm[0]], *this);
      hypo_tmp->fillTrack("pi+_1", trackIndices_[pl_perm[1]], *this);
      
      hypo->updateInitialParams();
      auto p_pich = hypo->getInitialMomentum(spich_);
      Eigen::VectorXd tmp_pi0_pars(3);
      tmp_pi0_pars << -p_pich.Px(), -p_pich.Py(), -p_pich.Pz();
      hypo->setInitialParticleParams("pi0", tmp_pi0_pars); 
      hypo->optimize();

      hypo_tmp->updateInitialParams();
      auto p_pich_tmp = hypo_tmp->getInitialMomentum(spich_);
      Eigen::VectorXd tmp_pi0_pars_tmp(3);
      tmp_pi0_pars_tmp << -p_pich_tmp.Px(), -p_pich_tmp.Py(), -p_pich_tmp.Pz();
      hypo_tmp->setInitialParticleParams("pi0", tmp_pi0_pars_tmp); 
      hypo_tmp->optimize();
      
      if (kf_err_ != 0 && kf_err_ != 2) {
        kf_err_ = hypo->getErrorCode();
      }
      double tmp_chi2 = hypo->getChiSquare();
      if (tmp_chi2 < kf_chi2_ && hypo->getErrorCode() == 0) {
        kf_err_ = hypo->getErrorCode();
        kf_chi2_ = tmp_chi2;
        in_m3pi_ = hypo->getInitialMomentum(sM3Pi_).M();
        kf_m3pi_ = hypo->getFinalMomentum(sM3Pi_).M();
	kf_vtx_m3pi_ = hypo_tmp->getInitialMomentum(sM3Pi_).M();
	auto vtx = hypo->getFinalVertex("vtx0");
        vtx.GetXYZ(vtx_);
      }
    } while (std::next_permutation(pl_perm.begin(), pl_perm.end()));
  } while (std::next_permutation(mi_perm.begin(), mi_perm.end()));
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::Hypo4ChPionsLostPi0_EtaMass hypo(2.e-3 * emeas, magneticField, 100);
  kfcmd::hypos::Hypo4PiVertex hypo_vertex(2.e-3 * emeas, magneticField, 100);
  kfcmd::hypos::Hypo4ChPionsLostPi0_EtaMass hypo_tmp(2.e-3 * emeas, magneticField, 100);
  Eigen::MatrixXd icov = Eigen::MatrixXd::Zero(3, 3);
  icov(0, 0) = 1. / (6.e-2 * 6.e-2);
  icov(1, 1) = 1. / (6.e-2 * 6.e-2);
  hypo.setVertexInverseCovarianceMatrix("vtx0", icov);

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    fit_(&hypo, &hypo_vertex, &hypo_tmp);
    out_tree->Fill();
  }
  outfl->cd();
 out_tree->Write();
  delete outfl;
}
