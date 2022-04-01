#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <set>

#include <TStopwatch.h>

#include <kfcmd/hypos/Hypo3ChPionsKPlus.hpp>
#include <kfcmd/hypos/Hypo3ChPionsKMinus.hpp>

#include "TrPh.h"

const double TrPh::dZ_ = 30;
const double TrPh::dRho_ = 30;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;

TrPh::TrPh(TTree *tree) :
  kfcmd::core::TrPh(tree) {
}

TrPh::~TrPh() {
}

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
  if (nt < 4) return -1;
  if (!cutTracks_()) return -1;
  std::vector<Int_t> charges(nt);
  std::copy(tcharge, tcharge + nt, charges.begin());
  std::sort(trackIndices_.begin(), trackIndices_.end(),
            [&charges](int i, int j) { return charges[i] < charges[j]; });
  return 1;
}

void TrPh::setupOutptuBranches_(TTree* tree) {
  tree->Branch("kf_mks", &kf_mks_, "kf_mks/D");
  tree->Branch("in_mks", &in_mks_, "in_mks/D");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");
  tree->Branch("kf_vtx0_x", &kf_vtx0_x_, "kf_vtx0_x/D");
  tree->Branch("kf_vtx0_y", &kf_vtx0_y_, "kf_vtx0_y/D");
  tree->Branch("kf_vtx0_z", &kf_vtx0_z_, "kf_vtx0_z/D");
  tree->Branch("kf_vtx1_x", &kf_vtx1_x_, "kf_vtx1_x/D");
  tree->Branch("kf_vtx1_y", &kf_vtx1_y_, "kf_vtx1_y/D");
  tree->Branch("kf_vtx1_z", &kf_vtx1_z_, "kf_vtx1_z/D");
  tree->Branch("kf_vtx_dr", &kf_vtx_dr_, "kf_vtx_dr/D");
  tree->Branch("kf_vtx_drho", &kf_vtx_drho_, "kf_vtx_drho/D");
  tree->Branch("kf_dedx_vtx0_K", &kf_dedx_vtx0_K_, "kf_dedx_vtx0_K/D");
  tree->Branch("kf_dedx_vtx0_pi", &kf_dedx_vtx0_pi_, "kf_dedx_vtx0_pi/D");
  tree->Branch("kf_dedx_vtx1_pi", kf_dedx_vtx1_pi_, "kf_dedx_vtx1_pi[2]/D");
  tree->Branch("kf_p_vtx0_K", &kf_p_vtx0_K_, "kf_p_vtx0_K/D");
  tree->Branch("kf_p_vtx0_pi", &kf_p_vtx0_pi_, "kf_p_vtx0_pi/D");
  tree->Branch("kf_p_vtx1_pi", kf_p_vtx1_pi_, "kf_p_vtx1_pi[2]/D");
}

int TrPh::getStatus_() const {
  return status_;
}

void TrPh::setStatus_(int status) {
  status_ = status;
}

void TrPh::printStatus_(int entry_number, int nentries) {
  int status =  (100.0 * entry_number) / nentries;
  if (getStatus_() == status) return;
  setStatus_(status);
  timePerEntry_.Stop();
  std::cout << " [STATUS : " <<
    std::setfill('0') << std::setw(2) <<
    getStatus_()  << "%]" <<
    std::setw(0) <<
    "\tCPU TIME: " <<
    std::setprecision(3) <<
    std::fixed <<
    timePerEntry_.CpuTime() <<
    "\tREAL TIME: " <<
    timePerEntry_.RealTime() << std::endl;
  timePerEntry_.Start();
}

void TrPh::printSummary_(int* npassed, int ncutted, int nentries) {
  time_.Stop();
  std::cout << "_________________________________________" << std::endl;
  std::cout << "TOTAL NUMBER OF EVENTS BEFORE CUTS:\t" << nentries << std::endl;
  std::cout << "NUBER OF EVENTS TO THE INPUT OF KINFIT:\t" << ncutted << std::endl;
  std::cout << "NUMBER OF EVENTS PASSED e+e- --> KsK+pi- HYPO:\t" << npassed[0] << std::endl;
  std::cout << "NUMBER OF EVENTS PASSED e+e- --> KsK-pi+ HYPO:\t" << npassed[1] << std::endl;
  std::cout << "NUMBER OF EVENTS PASSED e+e- --> KsK+-pi-+ HYPOS:\t" << npassed[0] + npassed[1] << std::endl;
  std::cout << "TOTAL CPU TIME :\t" << time_.CpuTime() << std::endl;
  std::cout << "TOTAL REAL TIME :\t" << time_.RealTime() << std::endl;
  std::cout << "_________________________________________" << std::endl;
}

void TrPh::Loop(const std::string& outpath, double magneticField) {
  if (fChain == 0) return;
  time_.Start();
  std::set<std::string> sKs = {"pi+_0", "pi-_0"};
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutptuBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::Hypo3ChPionsKPlus hypo_plus(2 * emeas, magneticField);
  hypo_plus.setBeamXY(xbeam, ybeam);
  hypo_plus.fixVertexComponent("vtx0", xbeam, kfbase::core::VERTEX_X);
  hypo_plus.fixVertexComponent("vtx0", ybeam, kfbase::core::VERTEX_Y);
  kfcmd::hypos::Hypo3ChPionsKMinus hypo_minus(2 * emeas, magneticField);
  hypo_minus.setBeamXY(xbeam, ybeam);
  hypo_minus.fixVertexComponent("vtx0", xbeam, kfbase::core::VERTEX_X);
  hypo_minus.fixVertexComponent("vtx0", ybeam, kfbase::core::VERTEX_Y);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  double kf_chi2_plus;
  double tchi2_plus;
  double kf_chi2_minus;
  double tchi2_minus;
  bool flag_plus;
  bool flag_minus;
  int errCode;
  double v_kf_mks_plus = 0;
  double v_in_mks_plus = 0;
  double v_kf_mks_minus = 0;
  double v_in_mks_minus = 0;

  double v_dedx_vtx0_k_plus = 0;
  double v_p_vtx0_k_plus = 0;
  double v_dedx_vtx0_k_minus = 0;
  double v_p_vtx0_k_minus = 0;

  double v_dedx_vtx1_pi_plus[2] = {0., 0.};
  double v_p_vtx1_pi_plus[2] = {0., 0.};
  double v_dedx_vtx1_pi_minus[2] = {0., 0.};
  double v_p_vtx1_pi_minus[2] = {0., 0.};

  double v_dedx_vtx0_pi_plus = 0;
  double v_p_vtx0_pi_plus = 0;
  double v_dedx_vtx0_pi_minus = 0;
  double v_p_vtx0_pi_minus = 0;

  TVector3 v_vtx0_plus;
  TVector3 v_vtx1_plus;
  TVector3 v_vtx0_minus;
  TVector3 v_vtx1_minus;
  int ncutted = 0;
  int npassed[2] = {0, 0};
  timePerEntry_.Start();
  for (Long64_t jentry = 0; jentry < nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (Cut(ientry) < 0) continue;
    ncutted++;
    hypo_plus.setBeamXY(xbeam, ybeam);
    hypo_plus.fixVertexComponent("vtx0", xbeam, kfbase::core::VERTEX_X);
    hypo_plus.fixVertexComponent("vtx0", ybeam, kfbase::core::VERTEX_Y);
    hypo_minus.setBeamXY(xbeam, ybeam);
    hypo_minus.fixVertexComponent("vtx0", xbeam, kfbase::core::VERTEX_X);
    hypo_minus.fixVertexComponent("vtx0", ybeam, kfbase::core::VERTEX_Y);
    flag_plus = false;
    flag_minus = false;
    kf_chi2_plus = std::numeric_limits<double>::infinity();
    kf_chi2_minus = std::numeric_limits<double>::infinity();
    for (int im = 0; im < 2; ++im) {
      if (!hypo_plus.fillTrack("pi-_0", trackIndices_[im], *this)) continue;
      if (!hypo_plus.fillTrack("pi-_1", trackIndices_[1 - im], *this)) continue;
      for (int ip = 2; ip < 4; ++ip) {
        if (!hypo_plus.fillTrack("pi+_0", trackIndices_[ip], *this)) continue;
        if (!hypo_plus.fillTrack("k+", trackIndices_[5 - ip], *this)) continue;
        hypo_plus.optimize();
        errCode = hypo_plus.getErrorCode();
        if (errCode != 0) continue;
        tchi2_plus = hypo_plus.getChiSquare();
        if (tchi2_plus < kf_chi2_plus) {
          flag_plus = true;
          kf_chi2_plus = tchi2_plus;
          v_in_mks_plus = hypo_plus.getInitialMomentum(sKs).M();
          v_kf_mks_plus = hypo_plus.getFinalMomentum(sKs).M();
          v_vtx0_plus = hypo_plus.getFinalVertex("vtx0");
          v_vtx1_plus = hypo_plus.getFinalVertex("vtx1");

          v_dedx_vtx0_k_plus = tdedx[trackIndices_[5 - ip]];
          v_p_vtx0_k_plus = hypo_plus.getFinalMomentum("k+").P();
          v_dedx_vtx1_pi_plus[0] = tdedx[trackIndices_[im]];
          v_dedx_vtx1_pi_plus[1] = tdedx[trackIndices_[ip]];
          v_p_vtx1_pi_plus[0] = hypo_plus.getFinalMomentum("pi-_0").P();
          v_p_vtx1_pi_plus[1] = hypo_plus.getFinalMomentum("pi+_0").P();
          v_dedx_vtx0_pi_plus = tdedx[trackIndices_[1 - im]];
          v_p_vtx0_pi_plus = hypo_plus.getFinalMomentum("pi-_1").P();
        }
      }
    }

    for (int im = 0; im < 2; ++im) {
      if (!hypo_minus.fillTrack("pi-_0", trackIndices_[im], *this)) continue;
      if (!hypo_minus.fillTrack("k-", trackIndices_[1 - im], *this)) continue;
      for (int ip = 2; ip < 4; ++ip) {
        if (!hypo_minus.fillTrack("pi+_0", trackIndices_[ip], *this)) continue;
        if (!hypo_minus.fillTrack("pi+_1", trackIndices_[5 - ip], *this)) continue;
        hypo_minus.optimize();
        errCode = hypo_minus.getErrorCode();
        if (errCode != 0) continue;
        tchi2_minus = hypo_minus.getChiSquare();
        if (tchi2_minus < kf_chi2_minus) {
          flag_minus = true;
          kf_chi2_minus = tchi2_minus;
          v_in_mks_minus = hypo_minus.getInitialMomentum(sKs).M();
          v_kf_mks_minus = hypo_minus.getFinalMomentum(sKs).M();
          v_vtx0_minus = hypo_minus.getFinalVertex("vtx0");
          v_vtx1_minus = hypo_minus.getFinalVertex("vtx1");

          v_dedx_vtx0_k_minus = tdedx[trackIndices_[1 - im]];
          v_p_vtx0_k_minus = hypo_minus.getFinalMomentum("k-").P();
          v_dedx_vtx1_pi_minus[0] = tdedx[trackIndices_[im]];
          v_dedx_vtx1_pi_minus[1] = tdedx[trackIndices_[ip]];
          v_p_vtx1_pi_minus[0] = hypo_minus.getFinalMomentum("pi-_0").P();
          v_p_vtx1_pi_minus[1] = hypo_minus.getFinalMomentum("pi+_0").P();
          v_dedx_vtx0_pi_minus = tdedx[trackIndices_[5 - ip]];
          v_p_vtx0_pi_minus = hypo_minus.getFinalMomentum("pi+_1").P();
        }
      }
    }

    if (flag_plus && flag_minus) {
      if (kf_chi2_plus > kf_chi2_minus) {
        flag_plus = false;
      } else {
        flag_minus = false;
      }
    }

    if (flag_plus) {
      npassed[0]++;
      kf_mks_ = v_kf_mks_plus;
      in_mks_ = v_in_mks_plus;
      kf_chi2_ = kf_chi2_plus;
      kf_vtx0_x_ = v_vtx0_plus.X();
      kf_vtx0_y_ = v_vtx0_plus.Y();
      kf_vtx0_z_ = v_vtx0_plus.Z();
      kf_vtx1_x_ = v_vtx1_plus.X();
      kf_vtx1_y_ = v_vtx1_plus.Y();
      kf_vtx1_z_ = v_vtx1_plus.Z();
      kf_vtx_dr_ = (v_vtx1_minus - v_vtx0_minus).Mag();
      kf_vtx_drho_ = (v_vtx1_minus - v_vtx0_minus).Perp();
      kf_dedx_vtx0_K_ = v_dedx_vtx0_k_plus;
      kf_p_vtx0_K_ = v_p_vtx0_k_plus;
      kf_dedx_vtx0_pi_ = v_dedx_vtx0_pi_plus;
      kf_p_vtx0_pi_ = v_p_vtx0_pi_plus;
      std::copy(v_dedx_vtx1_pi_plus, v_dedx_vtx1_pi_plus + 2, kf_dedx_vtx1_pi_);
      std::copy(v_p_vtx1_pi_plus, v_p_vtx1_pi_plus + 2, kf_p_vtx1_pi_);
    }

    if (flag_minus) {
      npassed[1]++;
      kf_mks_ = v_kf_mks_minus;
      in_mks_ = v_in_mks_minus;
      kf_chi2_ = kf_chi2_minus;
      kf_vtx0_x_ = v_vtx0_minus.X();
      kf_vtx0_y_ = v_vtx0_minus.Y();
      kf_vtx0_z_ = v_vtx0_minus.Z();
      kf_vtx1_x_ = v_vtx1_minus.X();
      kf_vtx1_y_ = v_vtx1_minus.Y();
      kf_vtx1_z_ = v_vtx1_minus.Z();
      kf_vtx_dr_ = (v_vtx1_minus - v_vtx0_minus).Mag();
      kf_vtx_drho_ = (v_vtx1_minus - v_vtx0_minus).Perp();
      kf_dedx_vtx0_K_ = v_dedx_vtx0_k_minus;
      kf_p_vtx0_K_ = v_p_vtx0_k_minus;
      kf_dedx_vtx0_pi_ = v_dedx_vtx0_pi_minus;
      kf_p_vtx0_pi_ = v_p_vtx0_pi_minus;
      std::copy(v_dedx_vtx1_pi_minus, v_dedx_vtx1_pi_minus + 2, kf_dedx_vtx1_pi_);
      std::copy(v_p_vtx1_pi_minus, v_p_vtx1_pi_minus + 2, kf_p_vtx1_pi_);
    }

    out_tree->Fill();
    printStatus_(jentry, nentries);
  }
  printSummary_(npassed, ncutted, nentries);
  outfl->cd();
  out_tree->Write();
  outfl->Close();
}
