#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <set>

#include <TStopwatch.h>

#include "TrPh.h"

const double TrPh::dZ_ = 30;
const double TrPh::dRho_ = 30;
const double TrPh::mindEdX_ = 0;
const double TrPh::maxdEdX_ = 15000;
const double TrPh::minTPtot_ = 5;
const double TrPh::maxTPtot_ = 1000;
const std::set<std::string> TrPh::decaypds_ = {"pi-_1", "pi+_1"};

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

double TrPh::getAngle_(const kfhypos::HypoKsKPlusPiMinus_NoKsMass &hypo) const {
  auto dv = hypo.getFinalVertex("vtx1") - hypo.getFinalVertex("vtx0");
  auto p = hypo.getFinalMomentum(decaypds_).Vect();
  return dv.Angle(p);
}

double TrPh::getAngle_(const kfhypos::HypoKsKMinusPiPlus_NoKsMass &hypo) const {
  auto dv = hypo.getFinalVertex("vtx1") - hypo.getFinalVertex("vtx0");
  auto p = hypo.getFinalMomentum(decaypds_).Vect();
  return dv.Angle(p);
}

void TrPh::setupOutputBranches_(TTree* tree) {
  tree->Branch("evnum", &evnum, "evnum/I");
  tree->Branch("kf_err", &kf_err_, "kf_err/I");
  tree->Branch("kf_chi2", &kf_chi2_, "kf_chi2/D");

  tree->Branch("kf_ks_px", &kf_ks_px_, "kf_ks_px/D");
  tree->Branch("kf_ks_py", &kf_ks_py_, "kf_ks_py/D");
  tree->Branch("kf_ks_pz", &kf_ks_pz_, "kf_ks_pz/D");
  tree->Branch("kf_ks_pe", &kf_ks_pe_, "kf_ks_pe/D");

  tree->Branch("kf_ks_vc_x", &kf_ks_vc_x_, "kf_ks_vc_x/D");
  tree->Branch("kf_ks_vc_y", &kf_ks_vc_y_, "kf_ks_vc_y/D");
  tree->Branch("kf_ks_vc_z", &kf_ks_vc_z_, "kf_ks_vc_z/D");

  tree->Branch("in_tks", &in_tks_, "in_tks/D");
  tree->Branch("kf_tks", &kf_tks_, "kf_tks/D");
  tree->Branch("kf_mks", &kf_mks_, "kf_mks/D");
  tree->Branch("in_mks", &in_mks_, "in_mks/D");
  tree->Branch("sim_vtx_dr", &sim_vtx_dr_, "sim_vtx_dr/D");
  tree->Branch("sim_vtx_drho", &sim_vtx_drho_, "sim_vtx_drho/D");
  tree->Branch("kf_ks_decay_prod_angle",
               &kf_ks_decay_prod_angle_,
               "kf_ks_decay_prod_angle/D");
  tree->Branch("kf_hypo", &kf_hypo_, "kf_hypo/I");
  tree->Branch("sim_hypo", &sim_hypo_, "sim_hypo/I");
  tree->Branch("sim_ee_vtx", sim_ee_vtx_, "sim_ee_vtx[3]/D");
  tree->Branch("sim_ks_vtx", sim_ks_vtx_, "sim_ks_vtx[3]/D");
  tree->Branch("kf_vtx0", kf_vtx0_, "kf_vtx0[3]/D");
  tree->Branch("kf_vtx1", kf_vtx1_, "kf_vtx1[3]/D");
  tree->Branch("kf_vtx_dr", &kf_vtx_dr_, "kf_vtx_dr/D");
  tree->Branch("kf_vtx_drho", &kf_vtx_drho_, "kf_vtx_drho/D");
  tree->Branch("kf_dedx_vtx0_K", &kf_dedx_vtx0_K_, "kf_dedx_vtx0_K/D");
  tree->Branch("kf_dedx_vtx0_pi", &kf_dedx_vtx0_pi_, "kf_dedx_vtx0_pi/D");
  tree->Branch("kf_dedx_vtx1_pi", kf_dedx_vtx1_pi_, "kf_dedx_vtx1_pi[2]/D");
  tree->Branch("kf_p_vtx0_K", &kf_p_vtx0_K_, "kf_p_vtx0_K/D");
  tree->Branch("kf_p_vtx0_pi", &kf_p_vtx0_pi_, "kf_p_vtx0_pi/D");
  tree->Branch("kf_p_vtx1_pi", kf_p_vtx1_pi_, "kf_p_vtx1_pi[2]/D");

  tree->Branch("kf_ct_out_pipl1", &kf_ct_out_pipl1_, "kf_ct_out_pipl1/D");
  tree->Branch("kf_ct_out_pimi1", &kf_ct_out_pimi1_, "kf_ct_out_pimi1/D");
  tree->Branch("kf_ct_out_pich0", &kf_ct_out_pich0_, "kf_ct_out_pich0/D");
  tree->Branch("kf_ct_out_kch0", &kf_ct_out_kch0_, "kf_ct_out_kch0/D");
}

void TrPh::fillSimInfo_() {
  sim_hypo_ = -1;
  for (int i = 0; i < nsim; ++i) {
    switch (simorig[i]) {
    case 0:
      sim_ee_vtx_[0] = simvtx[i];
      sim_ee_vtx_[1] = simvty[i];
      sim_ee_vtx_[2] = simvtz[i];
      switch (simtype[i]) {
      case 321:
        sim_hypo_ = 0;
        break;
      case -321:
        sim_hypo_ = 1;
        break;
      default:
        break;
      }
      break;
    case 310:
      sim_ks_vtx_[0] = simvtx[i];
      sim_ks_vtx_[1] = simvty[i];
      sim_ks_vtx_[2] = simvtz[i];
      break;
    default:
      break;
    }
  }
  auto dvtx = TVector3(sim_ks_vtx_) - TVector3(sim_ee_vtx_);
  sim_vtx_dr_ = dvtx.Mag();
  sim_vtx_drho_ = dvtx.Perp();
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
  auto outfl = TFile::Open(outpath.c_str(), "recreate");
  TTree* out_tree = new TTree("kf_data", "");
  setupOutputBranches_(out_tree);
  fChain->GetEntry(0);
  kfcmd::hypos::HypoKsKPlusPiMinus_NoKsMass hypo_plus(2.e-3 * emeas, magneticField);
  kfcmd::hypos::HypoKsKMinusPiPlus_NoKsMass hypo_minus(2.e-3 * emeas, magneticField);
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  double kf_ks_prod_angle_plus = 0;
  double kf_ks_prod_angle_minus = 0;

  double in_tks_plus = 0;
  double kf_tks_plus = 0;
  double in_tks_minus = 0;
  double kf_tks_minus = 0;

  double kf_ks_px_plus = 0;
  double kf_ks_py_plus = 0;
  double kf_ks_pz_plus = 0;
  double kf_ks_pe_plus = 0;
  double kf_ks_vc_x_plus = 0;
  double kf_ks_vc_y_plus = 0;
  double kf_ks_vc_z_plus = 0;

  double kf_ks_px_minus = 0;
  double kf_ks_py_minus = 0;
  double kf_ks_pz_minus = 0;
  double kf_ks_pe_minus = 0;
  double kf_ks_vc_x_minus = 0;
  double kf_ks_vc_y_minus = 0;
  double kf_ks_vc_z_minus = 0;

  int kf_err_plus;
  int kf_err_minus;

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
    kf_err_plus = 1;
    kf_err_minus = 1;
    kf_err_ = 1;
    kf_hypo_ = -1;
    fillSimInfo_();
    ncutted++;
    hypo_plus.setBeamXY(xbeam, ybeam);
    hypo_plus.fixVertexParameter("vtx0", 0, xbeam);
    hypo_plus.fixVertexParameter("vtx0", 1, ybeam);
    hypo_minus.setBeamXY(xbeam, ybeam);
    hypo_minus.fixVertexParameter("vtx0", 0, xbeam);
    hypo_minus.fixVertexParameter("vtx0", 1, ybeam);
    flag_plus = false;
    flag_minus = false;
    kf_chi2_plus = std::numeric_limits<double>::infinity();
    kf_chi2_minus = std::numeric_limits<double>::infinity();
    for (int im = 0; im < 2; ++im) {
      if (!hypo_plus.fillTrack("pi-_1", trackIndices_[im], *this)) continue;
      if (!hypo_plus.fillTrack("pi-_0", trackIndices_[1 - im], *this)) continue;
      for (int ip = 2; ip < 4; ++ip) {
        if (!hypo_plus.fillTrack("pi+_1", trackIndices_[ip], *this)) continue;
        if (!hypo_plus.fillTrack("k+", trackIndices_[5 - ip], *this)) continue;
        hypo_plus.updateInitialParams();
        auto ks_im = hypo_plus.getInitialMomentum("pi-_1") +
          hypo_plus.getInitialMomentum("pi+_1");
        Eigen::VectorXd tmpv(4);
        tmpv << ks_im.Px(), ks_im.Py(), ks_im.Pz(), 1.e-3;
        hypo_plus.setInitialParticleParams("ks", tmpv);
        hypo_plus.optimize();
        errCode = hypo_plus.getErrorCode();
        if (errCode != 0) {
             if (errCode != 1) kf_err_ = errCode; 
            continue;
        }
        kf_err_plus = 0;
        tchi2_plus = hypo_plus.getChiSquare();
        if (tchi2_plus < kf_chi2_plus) {
          flag_plus = true;
          kf_chi2_plus = tchi2_plus;

          in_tks_plus = hypo_plus.getParticleInitialParams("ks")(3);
          kf_tks_plus = hypo_plus.getParticleFinalParams("ks")(3);
          auto ksP = hypo_plus.getFinalMomentum("ks");
          kf_ks_px_plus = ksP.Px();
          kf_ks_py_plus = ksP.Py();
          kf_ks_pz_plus = ksP.Pz();
          kf_ks_pe_plus = ksP.E();
          v_vtx0_plus = hypo_plus.getFinalVertex("vtx0");
          v_vtx1_plus = hypo_plus.getFinalVertex("vtx1");
          auto kf_ks_vc = v_vtx0_plus + ksP.Vect() * kf_tks_plus - v_vtx1_plus;
          kf_ks_vc_x_plus = kf_ks_vc.X();
          kf_ks_vc_y_plus = kf_ks_vc.Y();
          kf_ks_vc_z_plus = kf_ks_vc.Z();
          kf_ks_prod_angle_plus = getAngle_(hypo_plus);
          v_in_mks_plus = hypo_plus.getInitialMomentum(decaypds_).M();
          v_kf_mks_plus = hypo_plus.getFinalMomentum(decaypds_).M();
          v_dedx_vtx0_k_plus = tdedx[trackIndices_[5 - ip]];
          v_p_vtx0_k_plus = hypo_plus.getFinalMomentum("k+").P();
          v_dedx_vtx1_pi_plus[0] = tdedx[trackIndices_[im]];
          v_dedx_vtx1_pi_plus[1] = tdedx[trackIndices_[ip]];
          v_p_vtx1_pi_plus[0] = hypo_plus.getFinalMomentum("pi-_1").P();
          v_p_vtx1_pi_plus[1] = hypo_plus.getFinalMomentum("pi+_1").P();
          v_dedx_vtx0_pi_plus = tdedx[trackIndices_[1 - im]];
          v_p_vtx0_pi_plus = hypo_plus.getFinalMomentum("pi-_0").P();
          kf_ct_out_pipl1_ = hypo_plus.getParticleFinalParams("pi+_1")(5);
          kf_ct_out_pimi1_ = hypo_plus.getParticleFinalParams("pi-_1")(5);
          kf_ct_out_pich0_ = hypo_plus.getParticleFinalParams("pi-_0")(5);
          kf_ct_out_kch0_ = hypo_plus.getParticleFinalParams("k+")(5);
        }
      }
    }

    for (int im = 0; im < 2; ++im) {
      if (!hypo_minus.fillTrack("pi-_1", trackIndices_[im], *this)) continue;
      if (!hypo_minus.fillTrack("k-", trackIndices_[1 - im], *this)) continue;
      for (int ip = 2; ip < 4; ++ip) {
        if (!hypo_minus.fillTrack("pi+_1", trackIndices_[ip], *this)) continue;
        if (!hypo_minus.fillTrack("pi+_0", trackIndices_[5 - ip], *this)) continue;
        hypo_minus.updateInitialParams();
        auto ks_im = hypo_minus.getInitialMomentum("pi-_1") +
                     hypo_minus.getInitialMomentum("pi+_1");
        Eigen::VectorXd tmpv(4);
        tmpv << ks_im.Px(), ks_im.Py(), ks_im.Pz(), 1.e-3;
        hypo_minus.setInitialParticleParams("ks", tmpv);
        hypo_minus.optimize();
        errCode = hypo_minus.getErrorCode();
        if (errCode != 0) {
            if (errCode != 1) kf_err_ = errCode;
            continue;
        }
        kf_err_minus = 0;
        tchi2_minus = hypo_minus.getChiSquare();
        if (tchi2_minus < kf_chi2_minus) {
          flag_minus = true;
          kf_chi2_minus = tchi2_minus;

          auto ksP = hypo_minus.getFinalMomentum("ks");
          kf_ks_px_minus = ksP.Px();
          kf_ks_py_minus = ksP.Py();
          kf_ks_pz_minus = ksP.Pz();
          kf_ks_pe_minus = ksP.E();
          in_tks_minus = hypo_minus.getParticleInitialParams("ks")(3);
          kf_tks_minus = hypo_minus.getParticleFinalParams("ks")(3);
          v_vtx0_minus = hypo_minus.getFinalVertex("vtx0");
          v_vtx1_minus = hypo_minus.getFinalVertex("vtx1");
          auto kf_ks_vc = v_vtx0_minus + ksP.Vect() * kf_tks_minus - v_vtx1_minus;
          kf_ks_vc_x_minus = kf_ks_vc.X();
          kf_ks_vc_y_minus = kf_ks_vc.Y();
          kf_ks_vc_z_minus = kf_ks_vc.Z();
          kf_ks_prod_angle_minus = getAngle_(hypo_minus);
          v_in_mks_minus = hypo_minus.getInitialMomentum(decaypds_).M();
          v_kf_mks_minus = hypo_minus.getFinalMomentum(decaypds_).M();
          v_dedx_vtx0_k_minus = tdedx[trackIndices_[1 - im]];
          v_p_vtx0_k_minus = hypo_minus.getFinalMomentum("k-").P();
          v_dedx_vtx1_pi_minus[0] = tdedx[trackIndices_[im]];
          v_dedx_vtx1_pi_minus[1] = tdedx[trackIndices_[ip]];
          v_p_vtx1_pi_minus[0] = hypo_minus.getFinalMomentum("pi-_1").P();
          v_p_vtx1_pi_minus[1] = hypo_minus.getFinalMomentum("pi+_1").P();
          v_dedx_vtx0_pi_minus = tdedx[trackIndices_[5 - ip]];
          v_p_vtx0_pi_minus = hypo_minus.getFinalMomentum("pi+_0").P();
          kf_ct_out_pipl1_ = hypo_minus.getParticleFinalParams("pi+_1")(5);
          kf_ct_out_pimi1_ = hypo_minus.getParticleFinalParams("pi-_1")(5);
          kf_ct_out_pich0_ = hypo_minus.getParticleFinalParams("pi+_0")(5);
          kf_ct_out_kch0_ = hypo_minus.getParticleFinalParams("k-")(5);
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
      kf_err_ = kf_err_plus;
      kf_hypo_ = 0;
      npassed[0]++;
      kf_mks_ = v_kf_mks_plus;
      in_mks_ = v_in_mks_plus;
      kf_chi2_ = kf_chi2_plus;
      in_tks_ = in_tks_plus;
      kf_tks_ = kf_tks_plus;

      kf_ks_px_ = kf_ks_px_plus;
      kf_ks_py_ = kf_ks_py_plus;
      kf_ks_pz_ = kf_ks_pz_plus;
      kf_ks_pe_ = kf_ks_pe_plus;
      kf_ks_vc_x_ = kf_ks_vc_x_plus;
      kf_ks_vc_y_ = kf_ks_vc_y_plus;
      kf_ks_vc_z_ = kf_ks_vc_z_plus;

      kf_ks_decay_prod_angle_ = kf_ks_prod_angle_plus;
      v_vtx0_plus.GetXYZ(kf_vtx0_);
      v_vtx1_plus.GetXYZ(kf_vtx1_);
      kf_vtx_dr_ = (v_vtx1_plus - v_vtx0_plus).Mag();
      kf_vtx_drho_ = (v_vtx1_plus - v_vtx0_plus).Perp();
      kf_dedx_vtx0_K_ = v_dedx_vtx0_k_plus;
      kf_p_vtx0_K_ = v_p_vtx0_k_plus;
      kf_dedx_vtx0_pi_ = v_dedx_vtx0_pi_plus;
      kf_p_vtx0_pi_ = v_p_vtx0_pi_plus;
      std::copy(v_dedx_vtx1_pi_plus, v_dedx_vtx1_pi_plus + 2, kf_dedx_vtx1_pi_);
      std::copy(v_p_vtx1_pi_plus, v_p_vtx1_pi_plus + 2, kf_p_vtx1_pi_);
    }

    if (flag_minus) {
      kf_err_ = kf_err_minus;
      kf_hypo_ = 1;
      npassed[1]++;
      kf_mks_ = v_kf_mks_minus;
      in_mks_ = v_in_mks_minus;
      kf_chi2_ = kf_chi2_minus;
      in_tks_ = in_tks_minus;
      kf_tks_ = kf_tks_minus;

      kf_ks_px_ = kf_ks_px_minus;
      kf_ks_py_ = kf_ks_py_minus;
      kf_ks_pz_ = kf_ks_pz_minus;
      kf_ks_pe_ = kf_ks_pe_minus;
      kf_ks_vc_x_ = kf_ks_vc_x_minus;
      kf_ks_vc_y_ = kf_ks_vc_y_minus;
      kf_ks_vc_z_ = kf_ks_vc_z_minus;

      kf_ks_decay_prod_angle_ = kf_ks_prod_angle_minus;
      v_vtx0_minus.GetXYZ(kf_vtx0_);
      v_vtx1_minus.GetXYZ(kf_vtx1_);
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
  delete outfl;
}
