//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 24 21:50:46 2019 by ROOT version 6.18/04
// from TTree tr_ph/Tree with the non-collinear events
// found on file: tr_ph_run021171.root
//////////////////////////////////////////////////////////

#ifndef _KFCmdEtaPiPi_EtaTo2Gamma_PiDecay_TrPh_H_
#define _KFCmdEtaPiPi_EtaTo2Gamma_PiDecay_TrPh_H_

#include <vector>
#include <kfcmd/core/TrPh.hpp>
#include <kfcmd/hypos/Hypo2ChPions2Photons_PiPlusDecay.hpp>
#include <kfcmd/hypos/Hypo2ChPions2Photons_PiMinusDecay.hpp>
#include <TVector3.h>

class TrPh : public kfcmd::core::TrPh {
 public:
  TrPh(TTree* tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t entry) override final;
  virtual void Loop(const std::string& outpath,
                    double magneticField = 1.3) override final;

 private:
  bool cutTracks_();
  bool cutPhotons_();
  TVector3 getVTX0Sim() const;
  TVector3 getVTX1Sim(int) const;
  void setupOutputBranches_(TTree *);
  void fit_(kfcmd::hypos::Hypo2ChPions2Photons_PiPlusDecay *);
  void fit_(kfcmd::hypos::Hypo2ChPions2Photons_PiMinusDecay *);
  static const std::set<std::string> sgg_;
  static const std::set<std::string> spipi_;
  static const double dZ_;
  static const double dRho_;
  static const double mindEdX_;
  static const double maxdEdX_;
  static const double minTPtot_;
  static const double maxTPtot_;
  static const double minPhEn_;
  static const double maxPhEn_;
  int kf_err_;
  double kf_chi2_;
  double in_mgg_;
  double kf_mgg_;
  double in_mpipi_;
  double kf_mpipi_;
  double kf_vtx0_[3];
  double kf_vtx1_[3];
  double kf_ct_out_pipl_;
  double kf_ct_in_pipl_;
  double kf_ct_out_pimi_;
  double kf_ct_in_pimi_;
  double kf_ct_out_mu_;
  double kf_vtx_dr_;
  double sim_vtx_dr_;
  double sim_vtx0_[3];
  double sim_vtx1_[3];
  // 0 - pi+ decay, 1 - pi- decay
  int hypo_;
  std::vector<std::size_t> trackIndices_;
  std::vector<std::size_t> photonIndices_;
};
#endif
