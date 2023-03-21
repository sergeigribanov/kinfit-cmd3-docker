#ifndef _KFCmdGaussSim_KsChKaonChPion_TrPh_
#define _KFCmdGaussSim_KsChKaonChPion_TrPh_
#include <vector>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom.h>
#include <kfcmd/core/TrPh.hpp>
#include <kfcmd/hypos/HypoKsKMinusPiPlus_NoKsMass.hpp>
#include <kfcmd/hypos/HypoKsKPlusPiMinus_NoKsMass.hpp>


namespace kfhypos = kfcmd::hypos;

typedef struct {
  int chMinus0_track_index;
  int chPlus0_track_index;
  int piMinus1_track_index;
  int piPlus1_track_index;
  Eigen::VectorXd parChMinus0;
  Eigen::VectorXd parChPlus0;
  Eigen::VectorXd parPiMinus1;
  Eigen::VectorXd parPiPlus1;
  Eigen::MatrixXd invCovChMinus0;
  Eigen::MatrixXd invCovChPlus0;
  Eigen::MatrixXd invCovPiMinus1;
  Eigen::MatrixXd invCovPiPlus1;
} InputParams;

typedef struct Histograms {
  TH1F *chi2Hist;
  TH1F *qHist;
  TH1F *sigma_z_vtx0;
  TH1F *sigma_x_vtx1;
  TH1F *sigma_y_vtx1;
  TH1F *sigma_z_vtx1;
  TH1F *vtx0_dz;
  TH1F *vtx1_dx;
  TH1F *vtx1_dy;
  TH1F *vtx1_dz;
  TH1F *vtx0_dz_pull;
  TH1F *vtx1_dx_pull;
  TH1F *vtx1_dy_pull;
  TH1F *vtx1_dz_pull;
  TH2F *chi2_vs_q;
  Histograms();
  virtual ~Histograms();
  void write(const std::string &);
} Histograms;

  class TrPh : public kfcmd::core::TrPh {
  public:
    TrPh(TTree *tree = 0);
    virtual ~TrPh();
    virtual Int_t Cut(Long64_t entry) override final;
    virtual void Loop(const std::string &,
                      double magneticField = 1.3) override final;
    void setEntry(int);
    void setNEvents(int);

  private:
    bool cutTracks_();
    bool detectSimHypo_();
    bool fitOnceKPlus_(kfcmd::hypos::HypoKsKPlusPiMinus_NoKsMass*);
    bool fitOnceKMinus_(kfcmd::hypos::HypoKsKMinusPiPlus_NoKsMass*);
    void fillParamsKPlus_(const kfcmd::hypos::HypoKsKPlusPiMinus_NoKsMass &);
    void fillParamsKMinus_(const kfcmd::hypos::HypoKsKMinusPiPlus_NoKsMass&);
    void refitKPlus_(kfcmd::hypos::HypoKsKPlusPiMinus_NoKsMass*, Histograms*);
    void refitKMinus_(kfcmd::hypos::HypoKsKMinusPiPlus_NoKsMass*, Histograms*);
    static const double dZ_;
    static const double dRho_;
    static const double mindEdX_;
    static const double maxdEdX_;
    static const double minTPtot_;
    static const double maxTPtot_;
    int entry_;
    int nevents_;
    int sim_hypo_;
    double kf_vtx0_[3];
    double kf_vtx1_[3];
    double sim_ee_vtx_[3];
    double sim_ks_vtx_[3];
    double sigma_z_vtx0_;
    double sigma_x_vtx1_;
    double sigma_y_vtx1_;
    double sigma_z_vtx1_;
    std::vector<std::size_t> trackIndices_;
    TRandom rnd_;
    InputParams params_;
  };

#endif
