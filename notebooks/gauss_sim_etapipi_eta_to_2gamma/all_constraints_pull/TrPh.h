#ifndef _KFCmdGaussSim_TrPh_EtaTo2Gamma_TrPh_H_
#define _KFCmdGaussSim_TrPh_EtaTo2Gamma_TrPh_H_
#include <vector>
#include <TRandom.h>
#include <TH1F.h>
#include <TH2F.h>
#include <kfcmd/core/TrPh.hpp>
#include <kfcmd/hypos/Hypo2ChPions2Photons.hpp>

namespace kfcore = kfcmd::core;
namespace kfhypos = kfcmd::hypos;

typedef struct {
  Eigen::VectorXd parPiMi;
  Eigen::VectorXd parPiPl;
  Eigen::VectorXd parGamma0;
  Eigen::VectorXd parGamma1;
  Eigen::MatrixXd invCovPiMi;
  Eigen::MatrixXd invCovPiPl;
  Eigen::MatrixXd invCovGamma0;
  Eigen::MatrixXd invCovGamma1;
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
  bool fitOnce_(kfhypos::Hypo2ChPions2Photons *);
  // void fillSimInfo_();
  void fillParams_(const kfhypos::Hypo2ChPions2Photons &);
  void refit_(kfhypos::Hypo2ChPions2Photons*, Histograms*);
  static const double dZ_;
  static const double dRho_;
  static const double mindEdX_;
  static const double maxdEdX_;
  static const double minTPtot_;
  static const double maxTPtot_;
  static const double minPhEn_;
  static const double maxPhEn_;
  int entry_;
  int nevents_;
  double kf_vtx_[3];
  double sim_ee_vtx_[3];
  double sigma_z_vtx_;
  std::vector<std::size_t> trackIndices_;
  std::vector<std::size_t> photonIndices_;
  TRandom rnd_;
  InputParams params_;
};

#endif
