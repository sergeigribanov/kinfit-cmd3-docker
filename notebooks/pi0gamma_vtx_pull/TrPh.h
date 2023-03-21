#ifndef _KFCmd_3Gamma_TrPh_
#define _KFCmd_3Gamma_TrPh_

#include <set>
#include <string>
#include <vector>
#include <kfcmd/core/TrPh.hpp>
#include "Hypo3PhotonsCustom.hpp"

class TrPh : public kfcmd::core::TrPh {
public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t energy) override final;
  virtual void Loop(const std::string&, double mfield = 1.3) override final;
private:
  bool cutPhotons_();
  void setupOutputBranches_(TTree *);
  void fit_(Hypo3PhotonsCustom*);
  void fillSimInfo_();
  std::vector<std::size_t> photonIndices_;
  static const double min_energy_;
  static const std::set<std::string> s_phpair0_;
  static const std::set<std::string> s_phpair1_;
  static const std::set<std::string> s_phpair2_;
  static const std::set<std::string> s_all_photons_;
  double kf_err_;
  double kf_chi2_;
  double in_mgg_[3];
  double kf_mgg_[3];
  double kf_vtx_[3];
  double in_total_p_[4];
  double kf_total_p_[4];
  double sim_ee_vtx_z_;
  double sigma_vtx0_z_;
};

#endif
