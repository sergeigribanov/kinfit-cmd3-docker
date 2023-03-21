#ifndef _KFCmdGaussSim_TrPh_EtaTo2Gamma_TrPh_H_
#define _KFCmdGaussSim_TrPh_EtaTo2Gamma_TrPh_H_
#include <vector>
#include <TRandom.h>
#include <TH1F.h>
#include <TH2F.h>
#include <kfcmd/core/TrPh.hpp>
#include "Hypo3PhotonsCustom.hpp"

namespace kfcore = kfcmd::core;

typedef struct {
  Eigen::VectorXd parGamma0;
  Eigen::VectorXd parGamma1;
  Eigen::VectorXd parGamma2;
  Eigen::VectorXd parVtx0;
  Eigen::MatrixXd invCovGamma0;
  Eigen::MatrixXd invCovGamma1;
  Eigen::MatrixXd invCovGamma2;
  Eigen::MatrixXd invCovVtx0;
} InputParams;

typedef struct Histograms {
  TH1F *chi2Hist;
  TH1F *qHist;
  TH1F *sigma_z;
  TH1F *vtx_dz;
  TH1F *vtx_dz_pull;
  TH2F *chi2_vs_q;
  Histograms();
  virtual ~Histograms();
  void write(const std::string&);
} Histograms;

class TrPh : public kfcmd::core::TrPh {
 public:
  TrPh(TTree* tree = 0);
  virtual ~TrPh();
  virtual Int_t Cut(Long64_t entry) override final;
  virtual void Loop(const std::string& outpath,
                    double magneticField = 1.3) override final;
  void setEntry(int);
  void setNEvents(int);
private:
  bool cutTracks_();
  bool cutPhotons_();
  bool fitOnce_(Hypo3PhotonsCustom *);
  // void fillSimInfo_();
  void fillParams_(const Hypo3PhotonsCustom &);
  void refit_(Hypo3PhotonsCustom*, Histograms*);
  int entry_;
  int nevents_;
  double kf_vtx_[3];
  double sim_ee_vtx_[3];
  double sigma_z_vtx_;
  std::vector<std::size_t> photonIndices_;
  TRandom rnd_;
  InputParams params_;
  static const double min_energy_;
  static const std::set<std::string> s_phpair0_;
  static const std::set<std::string> s_phpair1_;
  static const std::set<std::string> s_phpair2_;
  static const std::set<std::string> s_all_photons_;
};

#endif
