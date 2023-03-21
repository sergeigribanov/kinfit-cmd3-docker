#ifndef _KFCmd_KsChKaonChPion_TrPh_
#define _KFCmd_KsChKaonChPion_TrPh_

#include <kfcmd/core/TrPh.hpp>
#include <kfcmd/hypos/HypoKsKMinusPiPlus_NoKsMass.hpp>
#include <kfcmd/hypos/HypoKsKPlusPiMinus_NoKsMass.hpp>

namespace kfhypos = kfcmd::hypos;

class TrPh : public kfcmd::core::TrPh {
public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t entry) override final;
  virtual void Loop(const std::string&, double magneticField = 1.3) override final;
private:
  bool cutTracks_();
  double getAngle_(const kfhypos::HypoKsKPlusPiMinus_NoKsMass &) const;
  double getAngle_(const kfhypos::HypoKsKMinusPiPlus_NoKsMass &) const;
  void setupOutputBranches_(TTree*);
  void fillSimInfo_();
  void printStatus_(int, int);
  void printSummary_(int*, int, int);
  int getStatus_() const;
  void setStatus_(int);
  static const double dZ_;
  static const double dRho_;
  static const double mindEdX_;
  static const double maxdEdX_;
  static const double minTPtot_;
  static const double maxTPtot_;
  static const std::set<std::string> decaypds_;
  int status_;
  int numiters_;
  int kf_err_;
  double kf_chi2_;
  double kf_mks_;
  double in_mks_;
  double kf_vtx0_[3];
  double kf_vtx1_[3];
  double kf_vtx_dr_;
  double kf_vtx_drho_;
  double kf_dedx_vtx0_K_;
  double kf_dedx_vtx0_pi_;
  double kf_dedx_vtx1_pi_[2];
  double kf_p_vtx0_K_;
  double kf_p_vtx0_pi_;
  double kf_p_vtx1_pi_[2];
  double kf_ks_decay_prod_angle_;
  double in_tks_;
  double kf_tks_;

  double kf_ks_px_;
  double kf_ks_py_;
  double kf_ks_pz_;
  double kf_ks_pe_;

  double kf_ks_vc_x_;
  double kf_ks_vc_y_;
  double kf_ks_vc_z_;

  int kf_hypo_; // -1 - undefined, 0 - KsK+Pi-, 1 - KsK-Pi+
  int sim_hypo_; // -1 - undefined, 0 - KsK+Pi-, 1 - KsK-Pi+
  double sim_ee_vtx_[3];
  double sim_ks_vtx_[3];
  double sim_vtx_dr_;
  double sim_vtx_drho_;

  double kf_ct_out_pipl1_;
  double kf_ct_out_pimi1_;
  double kf_ct_out_pich0_;
  double kf_ct_out_kch0_;

  TStopwatch time_;
  TStopwatch timePerEntry_;
  std::vector<std::size_t> trackIndices_;
};

#endif
