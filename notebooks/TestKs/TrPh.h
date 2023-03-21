#ifndef _KFCmd_KsChKaonChPion_TrPh_
#define _KFCmd_KsChKaonChPion_TrPh_

#include <kfcmd/core/TrPh.hpp>


class TrPh : public kfcmd::core::TrPh {
public:
  TrPh(TTree *tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t entry) override final;
  virtual void Loop(const std::string&, double magneticField = 1.3) override final;
private:
  void setupOutputBranches_(TTree*);
  void fillSimInfo_();
  int trph_is_ks_;
  int sim_is_ks_;
  int kf_err_;
  int niter_;
  double kf_chi2_;
  double kf_mks_;
};

#endif
