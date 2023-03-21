#include <TVector3.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <set>
#include "TrPh.h"

const std::set<std::string> TrPh::sgg_ = {"g0", "g1"};
const std::set<std::string> TrPh::spipi_ = {"pi+", "pi-"};
const double TrPh::dZ_ = 30;
const double TrPh::dRho_ = 30;
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
  if (trackIndices_.size() != 3) return false;
  int totalCharge = 0;
  for (int i = 0; i < 3; ++i) totalCharge += tcharge[trackIndices_[i]];
  if (totalCharge == 1) {
    hypo_ = 0;
    return true;
  }
  if (totalCharge == -1) {
    hypo_ = 1;
    return true;
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

TVector3 TrPh::getVTX0Sim() const {
  TVector3 result;
  for (int i = 0; i < nsim; ++i) {
    if (simorig[i] == 0) {
      result = TVector3(simvtx[i], simvty[i], simvtz[i]);
      break;
    }
  }
  return result;
}

TVector3 TrPh::getVTX1Sim(int hypo) const {
  TVector3 result;
  for (int i = 0; i < nsim; ++i) {
     if (hypo == 0 && simorig[i] == 211) {
       result = TVector3(simvtx[i], simvty[i], simvtz[i]);
       break;
     }
     if (hypo == 1 && simorig[i] == -211) {
       result = TVector3(simvtx[i], simvty[i], simvtz[i]);
       break;
     }
  }
  return result;
}

void TrPh::setupOutputBranches_(TTree *tree) {
  tree->Branch("evnum", &evnum, "evnum/I");
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("hypo", &hypo_, "hypo/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("in_mgg", &in_mgg_, "in_mgg/D");
  tree->Branch("kf_mgg", &kf_mgg_, "kf_mgg/D");
  tree->Branch("in_mpipi", &in_mpipi_, "in_mpipi/D");
  tree->Branch("kf_mpipi", &kf_mpipi_, "kf_mpipi/D");
  tree->Branch("kf_vtx0", kf_vtx0_, "kf_vtx0[3]/D");
  tree->Branch("kf_vtx1", kf_vtx1_, "kf_vtx1[3]/D");
  tree->Branch("kf_ct_out_pimi", &kf_ct_out_pimi_, "kf_ct_out_pimi/D");
  tree->Branch("kf_ct_in_pimi", &kf_ct_in_pimi_, "kf_ct_in_pimi/D");
  tree->Branch("kf_ct_out_pipl", &kf_ct_out_pipl_, "kf_ct_out_pipl/D");
  tree->Branch("kf_ct_in_pipl", &kf_ct_in_pipl_, "kf_ct_in_pipl/D");
  tree->Branch("kf_ct_out_mu", &kf_ct_out_mu_, "kf_ct_out_mu/D");
  tree->Branch("kf_vtx_dr", &kf_vtx_dr_, "kf_vtx_dr/D");
  tree->Branch("sim_vtx_dr", &sim_vtx_dr_, "sim_vtx_dr/D");
  tree->Branch("sim_vtx0", sim_vtx0_, "sim_vtx0[3]/D");
  tree->Branch("sim_vtx1", sim_vtx1_, "sim_vtx1[3]/D");
}

void TrPh::fit_(kfcmd::hypos::Hypo2ChPions2Photons_PiPlusDecay *hypo) {
  hypo->setBeamXY(xbeam, ybeam);
  hypo->fixVertexParameter("vtx0", 0, xbeam);
  hypo->fixVertexParameter("vtx0", 1, ybeam);
  kf_err_ = 1;
  kf_chi2_ = std::numeric_limits<double>::infinity();
  std::vector<int> perm = {1, 2};
  if (!hypo->fillTrack("pi-", trackIndices_[0], *this)) {
    return;
  }
  do {
    if (!hypo->fillTrack("pi+", trackIndices_[perm[0]], *this)) {
      continue;
    }
    if (!hypo->fillTrack("mu+", trackIndices_[perm[1]], *this)) {
      continue;
    }
    for (std::size_t iph = 0; iph + 1 < photonIndices_.size(); ++iph) {
      if (!hypo->fillPhoton("g0", photonIndices_[iph], *this)) continue;
      for (std::size_t jph = iph + 1; jph < photonIndices_.size(); ++jph) {
        if (!hypo->fillPhoton("g1", photonIndices_[jph], *this)) continue;
        hypo->updateInitialParams();
        auto nu_p = hypo->getInitialMomentum("mu+").Vect() -
          hypo->getInitialMomentum("pi+").Vect();
        Eigen::VectorXd tmp_pars(3);
        tmp_pars << nu_p.Mag(), nu_p.Theta(), nu_p.Phi();
        hypo->setInitialParticleParams("nu", tmp_pars);
        Eigen::VectorXd tmp_pars_pipl = hypo->getParticleInitialParams("pi+");
        // Eigen::VectorXd tmp_pars_mupl = hypo->getParticleInitialParams("mu+");
        tmp_pars_pipl(6) = 20.; // std::sqrt(std::pow(tmp_pars_mupl(2), 2) + std::pow(tmp_pars_mupl(4), 2));
        hypo->setInitialParticleParams("pi+", tmp_pars_pipl);
        hypo->optimize();
        if (hypo->getErrorCode() != 0) continue;
        double tchi2 = hypo->getChiSquare();
        if (tchi2 >= kf_chi2_) continue;
        kf_err_ = 0;
        kf_chi2_ = tchi2;
        in_mgg_ = hypo->getInitialMomentum(sgg_).M();
        kf_mgg_ = hypo->getFinalMomentum(sgg_).M();
        in_mpipi_ = hypo->getInitialMomentum(spipi_).M();
        kf_mpipi_ = hypo->getFinalMomentum(spipi_).M();
        auto vtx0 = hypo->getFinalVertex("vtx0");
        auto vtx1 = hypo->getFinalVertex("vtx1");
        vtx0.GetXYZ(kf_vtx0_);
        vtx1.GetXYZ(kf_vtx1_);
        Eigen::VectorXd tpars = hypo->getParticleFinalParams("pi+");
        kf_ct_out_pipl_ = tpars(5);
        kf_ct_in_pipl_ = tpars(6);
        kf_ct_out_mu_ = hypo->getParticleFinalParams("mu+")(5);
        kf_vtx_dr_ = (vtx1 - vtx0).Mag();
        auto sim_vtx0 = getVTX0Sim();
        auto sim_vtx1 = getVTX1Sim(0);
        sim_vtx0.GetXYZ(sim_vtx0_);
        sim_vtx1.GetXYZ(sim_vtx1_);
        sim_vtx_dr_ = (sim_vtx1 - sim_vtx0).Mag();
      }
    }
  } while (std::next_permutation(perm.begin(), perm.end()));
}

void TrPh::fit_(kfcmd::hypos::Hypo2ChPions2Photons_PiMinusDecay *hypo) {
  hypo->setBeamXY(xbeam, ybeam);
  hypo->fixVertexParameter("vtx0", 0, xbeam);
  hypo->fixVertexParameter("vtx0", 1, ybeam);
  kf_err_ = 1;
  kf_chi2_ = std::numeric_limits<double>::infinity();
  std::vector<int> perm = {0, 1};
  if (!hypo->fillTrack("pi+", trackIndices_[2], *this)) {
    return;
  }
  do {
    if (!hypo->fillTrack("pi-", trackIndices_[perm[0]], *this)) {
      continue;
    }
    if (!hypo->fillTrack("mu-", trackIndices_[perm[1]], *this)) {
      continue;
    }
    for (std::size_t iph = 0; iph + 1 < photonIndices_.size(); ++iph) {
      if (!hypo->fillPhoton("g0", photonIndices_[iph], *this)) continue;
      for (std::size_t jph = iph + 1; jph < photonIndices_.size(); ++jph) {
        if (!hypo->fillPhoton("g1", photonIndices_[jph], *this)) continue;
        hypo->updateInitialParams();
        auto nu_p = hypo->getInitialMomentum("mu-").Vect() -
          hypo->getInitialMomentum("pi-").Vect();
        Eigen::VectorXd tmp_pars(3);
        tmp_pars << nu_p.Mag(), nu_p.Theta(), nu_p.Phi();
        hypo->setInitialParticleParams("nu", tmp_pars);
        Eigen::VectorXd tmp_pars_pimi = hypo->getParticleInitialParams("pi-");
        // Eigen::VectorXd tmp_pars_mumi = hypo->getParticleInitialParams("mu-");
        tmp_pars_pimi(6) = 20.; // std::sqrt(std::pow(tmp_pars_mumi(2), 2) +
                                //      std::pow(tmp_pars_mumi(4), 2));
        hypo->setInitialParticleParams("pi-", tmp_pars_pimi);
        hypo->optimize();
        if (hypo->getErrorCode() != 0) continue;
        double tchi2 = hypo->getChiSquare();
        if (tchi2 >= kf_chi2_) continue;
        kf_err_ = 0;
        kf_chi2_ = tchi2;
        in_mgg_ = hypo->getInitialMomentum(sgg_).M();
        kf_mgg_ = hypo->getFinalMomentum(sgg_).M();
        in_mpipi_ = hypo->getInitialMomentum(spipi_).M();
        kf_mpipi_ = hypo->getFinalMomentum(spipi_).M();
        auto vtx0 = hypo->getFinalVertex("vtx0");
        auto vtx1 = hypo->getFinalVertex("vtx1");
        vtx0.GetXYZ(kf_vtx0_);
        vtx1.GetXYZ(kf_vtx1_);
        Eigen::VectorXd tpars = hypo->getParticleFinalParams("pi-");
        kf_ct_out_pimi_ = tpars(5);
        kf_ct_in_pimi_ = tpars(6);
        kf_ct_out_mu_ = hypo->getParticleFinalParams("mu-")(5);
        kf_vtx_dr_ = (vtx1 - vtx0).Mag();
        auto sim_vtx0 = getVTX0Sim();
        auto sim_vtx1 = getVTX1Sim(1);
        sim_vtx0.GetXYZ(sim_vtx0_);
        sim_vtx1.GetXYZ(sim_vtx1_);
        sim_vtx_dr_ = (sim_vtx1 - sim_vtx0).Mag();
      }
    }
  } while (std::next_permutation(perm.begin(), perm.end()));
}

void TrPh::Loop(const std::string &outpath, double magneticField) {
  if (fChain == 0) return;
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::Hypo2ChPions2Photons_PiPlusDecay hypo_plus(2.e-3 * emeas, magneticField, 50);
  kfcmd::hypos::Hypo2ChPions2Photons_PiMinusDecay hypo_minus(2.e-3 * emeas, magneticField, 50);
  double tchi2;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    switch (hypo_) {
    case 0:
      fit_(&hypo_plus);
      break;
    case 1:
      fit_(&hypo_minus);
      break;
    default:
      break;
    }
    out_tree->Fill();
  }
  outfl->cd();
  out_tree->Write();
  delete outfl;
}
