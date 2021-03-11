#include "Plotting_Header.h"
// #include "/home/joshua/PCG_Software/Plotting/Plotting_Class.h"
#include "TGraphErrors.h"

// Fit to describe the NCellEffi
Double_t func1(Double_t *x, Double_t *par)
{
Double_t a = par[0];
Double_t b = par[1];
Double_t c = par[2];

Double_t res1 = TMath::Exp(a - TMath::Power((x[0]-b)/c,2));
Double_t res = 1-par[3]*TMath::Exp(- res1);

return res;
}

class Effi{

public:
  Effi();
  Effi(TString input, TString Period = "13TeV", TString suffix ="png");
  ~Effi();
  void FillHistos();
  float calcErr(float Emc, float Errmc, float Edata, float Errdata);
  void GetEffiHists(TH2F *hdata2d, TH2F *hMC2d, TH1D *&hEffiData,  TH1D *&hEffiMC, TH1D *&hRatio, TH1D *&hCorr, bool doSmooth = false);
  void GetTrueHists(TH2F* h2, TH1D *&h, TString name = "h");
  void FillCorrHistos();
  void LoadTB();
  void GetHistReweighted();
  void SetOtherHistos();
  void DoSBSubtraction();

  TH2F* SubtractSidebandBack(TH2F* hMInv, TH2F* hMInvBack, TH2F* hNCellSB, TH2F* hNCell);


  void SetPlotting();

  // Effi Plots
  void PlotEffi_AllClusAndTB();    // Effi Plot with all clusters from P2
  void PlotEffi_AllGammaAndTB();   // Effi Plot with clusters in pi0 mass range on right side of peak
  void PlotEffi_RWGammaAndTB();    // Effi Plot with clusters in pi0 mass range on right side of peak + data reweighted
  void PlotEffi_SBGammaAndTB();    // Effi Plot with clusters in pi0 mass range on right side of peak + data reweighted
  void PlotEffi_RWGammaHighAndTB_MCCheck();    // Effi Plot with clusters in pi0 mass range on right side of peak + data reweighted MC check
  void PlotEffi_RWGammaHighAndTB();    // Effi Plot with clusters in pi0 mass range on right side of peak + data reweighted
  void PlotEffi_TrueAndReweighted();    // Effi Plot with clusters in pi0 mass range on right side of peak + data reweighted + MC reweighted
  void PlotEffi_HighAllMethods();    // Effi Plot with clusters in pi0 mass range only higher gamma with RW and SB subtraction
  void PlotEffi_HighAllMethodsMCSBSub();    // Effi Plot with clusters in pi0 mass range only higher gamma with RW and SB subtraction
  void PlotEffi_HighAndWideRWSBSub();    // Effi Plot with clusters in pi0 mass range (wide and only high ) with RW and SB sub
  void PlotEffi_HighAndWideRWSBSub_Light();    // Effi Plot with clusters in pi0 mass range (wide and only high ) with RW and SB sub
  void PlotEffi_HighAndWideRWSBSub_Smoothed();    // Effi Plot with clusters in pi0 mass range (wide smoothed vs unsmoothed ) with RW and SB sub
  void PlotEffi_WideOnlySBSub();    // Effi Plot with clusters in pi0 mass range (wide smoothed vs unsmoothed ) with RW and SB sub
  void PlotEffi_WideSBSubGammaAndConv();    // Effi Plot with clusters in pi0 mass range wide, photons and conversions with reweighting
  void PlotEffi_HighAndWide_Light();    // Effi Plot with clusters in pi0 mass range (wide and only high ) without RW and SB sub
  void PlotEffi_WideSmoothed();    // Effi Plot with clusters in pi0 mass range (wide and smoothed ) without RW and SB sub
  void PlotEffi_Electrons();    // Effi Plot with clusters in pi0 mass range (wide and smoothed ) without RW and SB sub

  // Ratio Plots
  void PlotRatio_GammasAndRW();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_GammasAndRWLeft();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_GammasAndRWWide();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_GammasAndConvRWWide();   // Plot ratio with TB + gammas and conversions in wide range
  void PlotRatio_GammasAndRWHigh();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_GammasAndRWHighMCSBSub();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_SBGammasAndRW();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  // void PlotRatio_GammasAndRWWide();   // Plot ratio with TB + all clusters + gammas (wide) + gammas (wide) reweighted and SB sub

  // Corr Plots
  void PlotCorr_GammasAndRW();   // Plot correction factor with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotCorr_GammasAndRW_Left();   // Plot correction factor with TB + all clusters + gammas (left) + gammas (left) reweighted
  void PlotCorr_GammasAndRW_Wide();   // Plot correction factor with TB + all clusters + gammas (left) + gammas (left) reweighted
  void PlotCorr_GammasAndConv_Wide();   // Plot correction factor with TB gammas and conversions in wide range
  void PlotCorr_GammasAndRW_High();   // Plot correction factor with TB + all clusters + gammas (left) + gammas (left) reweighted
  void PlotCorr_HighAllMethods();   // Plot correction factor with TB + all clusters + gammas (left) + gammas (left) reweighted
  void PlotCorr_Electrons();   // Plot correction factor with TB + electrons

  // Fitting of CorrPlots
  void FitCorr_GammasWide();      // Fitting of correction with

  // Purity
  void PlotPurity();
  void PlotPurity2();
  // ExampleBin
  void PlotExampleBin();
  void PlotExampleBinCells();
  void ScaleTo(TH1D *h, TH1D *& h2, float down = 0.2, float up = 0.3);

  void PlotNCellRatios();
  void PlotFraction1cell();
  void PlotNCellDistr();

  TFile* fcontrolOut = nullptr;

private:

  TFile* fdata = nullptr;

  TString fPeriod = "";

  TString fsuffix = "png";

  float PlotenergyHigh = 8.;

  TString sEnergy = "pp #sqrt{s} = 13 TeV";

  // Test Beam
  TGraphErrors *grTB_data = nullptr;
  TGraphErrors *grTB_MC = nullptr;
  TGraphErrors *grTB_Ratio = nullptr;
  TGraphErrors *grTB_Corr = nullptr;

  // All clusters
  TH2F* hNCellVsETMNL_data = nullptr;
  TH2F* hNCellVsETMNL_MC = nullptr;

  // Different clusterizer
  TH2F* hNCellVsETMNLS500A100_data = nullptr;
  TH2F* hNCellVsETMNLS500A100_MC = nullptr;
  TH2F* hNCellVsETMNLS500A105_data = nullptr;
  TH2F* hNCellVsETMNLS500A105_MC = nullptr;
  TH2F* hNCellVsETMNLS500A110_data = nullptr;
  TH2F* hNCellVsETMNLS500A110_MC = nullptr;

  // TH2F* hNCellVsENL_data = nullptr;
  // Gammas
  TH2F* hNCellVsEGammasNL_data = nullptr;
  TH2F* hNCellVsEGammasNL_MC = nullptr;
  TH2F* hNCellVsETrueGammasNL_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElec_MC = nullptr;
  // Gammas Left
  TH2F* hNCellVsEGammasNLLeft_data = nullptr;
  TH2F* hNCellVsEGammasNLLeft_MC = nullptr;
  TH2F* hNCellVsETrueGammasNLLeft_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElecLeft_MC = nullptr;
  // Gammas Wide
  TH2F* hNCellVsEGammasNLWide_data = nullptr;
  TH2F* hNCellVsEGammasNLWide_MC = nullptr;
  TH2F* hNCellVsETrueGammasNLWide_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElecWide_MC = nullptr;
  // Gammas High
  TH2F* hNCellVsEGammasNLHigh_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_MC = nullptr;
  TH2F* hNCellVsETrueGammasNLHigh_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElecHigh_MC = nullptr;
  // Gammas High sideband subtracted
  TH2F* hNCellVsEGammasNLHigh_SBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_MCSBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_SBSub_MC = nullptr;
  // Gammas Wide sideband subtracted
  TH2F* hNCellVsEGammasNLWide_SBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLWide_MCSBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLWide_SBSub_MC = nullptr;
  // Gammas SideBand
  TH2F* hNCellVsEGammasNLSB_data = nullptr;
  TH2F* hNCellVsEGammasNLSB_MC = nullptr;
  TH2F* hNCellVsETrueGammasNLSB_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElecSB_MC = nullptr;
  // Gammas SideBand High
  TH2F* hNCellVsEGammasNLSBHigh_data = nullptr;
  TH2F* hNCellVsEGammasNLSBHigh_MC = nullptr;
  TH2F* hNCellVsETrueGammasNLSBHigh_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElecSBHigh_MC = nullptr;
  // Gammas Low
  TH2F* hNCellVsEGammasNLLow_data = nullptr;
  TH2F* hNCellVsEGammasNLLow_MC = nullptr;
  TH2F* hNCellVsETrueGammasNLLow_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElecLow_MC = nullptr;

  TH2F* hNCellVsEGammasNL_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLLeft_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLWide_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLWide_RW_SBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_SBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_MCSBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLSB_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLSBHigh_RW_data = nullptr;
  // conversions only
  TH2F* hNCellVsEConversionsNLWide_RW_data = nullptr;

  TH2F* hNCellVsEGammasNL_RW_MC = nullptr;
  TH2F* hNCellVsEGammasNLLeft_RW_MC = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_MC = nullptr;
  TH2F* hNCellVsEGammasNLWide_RW_MC = nullptr;
  TH2F* hNCellVsEGammasNLWide_RW_SBSub_MC = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_SBSub_MC = nullptr;
  TH2F* hNCellVsEGammasNLSBHigh_RW_MC = nullptr;

  // electron clusters
  TH2F* hNCellVsEelecNL_data = nullptr;
  TH2F* hNCellVsEelecNL_MC = nullptr;

  // TH2F* hNCellVsTrueEelecNL_data = nullptr;

  // all clusters
  TH1D* hNCell_AllClus_Effi_data = nullptr;
  TH1D* hNCell_AllClus_Effi_MC = nullptr;
  TH1D* hNCell_AllClus_Effi_Ratio = nullptr;
  TH1D* hNCell_AllClus_Effi_Corr = nullptr;

  // gammas
  TH1D* hNCell_Gammas_Effi_data = nullptr;
  TH1D* hNCell_Gammas_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_Effi_Corr = nullptr;

  TH1D* hNCell_TrueGammas_Effi_MC = nullptr;

  // true electrons
  TH1D* hNCell_TrueElec_Effi_MC = nullptr;

  // gammas left
  TH1D* hNCell_GammasLeft_Effi_data = nullptr;
  TH1D* hNCell_GammasLeft_Effi_MC = nullptr;
  TH1D* hNCell_GammasLeft_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasLeft_Effi_Corr = nullptr;

  // gammas wide
  TH1D* hNCell_GammasWide_Effi_data = nullptr;
  TH1D* hNCell_GammasWide_Effi_MC = nullptr;
  TH1D* hNCell_GammasWide_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasWide_Effi_Corr = nullptr;

  // gammas wide
  TH1D* hNCell_GammasWide_Smoothed_Effi_data = nullptr;
  TH1D* hNCell_GammasWide_Smoothed_Effi_MC = nullptr;
  TH1D* hNCell_GammasWide_Smoothed_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasWide_Smoothed_Effi_Corr = nullptr;

  // gammas SideBand
  TH1D* hNCell_GammasSB_Effi_data = nullptr;
  TH1D* hNCell_GammasSB_Effi_MC = nullptr;
  TH1D* hNCell_GammasSB_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasSB_Effi_Corr = nullptr;

  // gammas SideBand high
  TH1D* hNCell_GammasSBHigh_Effi_data = nullptr;
  TH1D* hNCell_GammasSBHigh_Effi_MC = nullptr;
  TH1D* hNCell_GammasSBHigh_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasSBHigh_Effi_Corr = nullptr;

  // gammas highestClus
  TH1D* hNCell_GammasHigh_Effi_data = nullptr;
  TH1D* hNCell_GammasHigh_Effi_MC = nullptr;
  TH1D* hNCell_GammasHigh_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasHigh_Effi_Corr = nullptr;

  // gammas lowestClus
  TH1D* hNCell_GammasLow_Effi_data = nullptr;
  TH1D* hNCell_GammasLow_Effi_MC = nullptr;
  TH1D* hNCell_GammasLow_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasLow_Effi_Corr = nullptr;


  // gammas reweighted
  TH1D* hNCell_Gammas_RW_Effi_data = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_Corr = nullptr;

  // gammas reweighted Pure MC
  TH1D* hNCell_Gammas_RW_Effi_MC_consistency = nullptr;
  // TH1D* hNCell_Gammas_RW_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_Ratio_consistency = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_Corr_consistency = nullptr;

  // gammas reweighted left
  TH1D* hNCell_GammasLeft_RW_Effi_data = nullptr;
  TH1D* hNCell_GammasLeft_RW_Effi_MC = nullptr;
  TH1D* hNCell_GammasLeft_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasLeft_RW_Effi_Corr = nullptr;

  // gammas reweighted wide
  TH1D* hNCell_GammasWide_RW_Effi_data = nullptr;
  TH1D* hNCell_GammasWide_RW_Effi_MC = nullptr;
  TH1D* hNCell_GammasWide_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasWide_RW_Effi_Corr = nullptr;

  // conversions reweighted wide
  TH1D* hNCell_ConversionsWide_RW_Effi_data = nullptr;
  TH1D* hNCell_ConversionsWide_RW_Effi_MC = nullptr;
  TH1D* hNCell_ConversionsWide_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_ConversionsWide_RW_Effi_Corr = nullptr;

  // gammas reweighted Sideband
  TH1D* hNCell_GammasSB_RW_Effi_data = nullptr;
  TH1D* hNCell_GammasSB_RW_Effi_MC = nullptr;
  TH1D* hNCell_GammasSB_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasSB_RW_Effi_Corr = nullptr;

  // gammas reweighted Sideband highest cluster
  TH1D* hNCell_GammasSBHigh_RW_Effi_data = nullptr;
  TH1D* hNCell_GammasSBHigh_RW_Effi_MC = nullptr;
  TH1D* hNCell_GammasSBHigh_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasSBHigh_RW_Effi_Corr = nullptr;

  // gammas reweighted highest clus
  TH1D* hNCell_GammasHigh_RW_Effi_data = nullptr;
  TH1D* hNCell_GammasHigh_RW_Effi_MC = nullptr;
  TH1D* hNCell_GammasHigh_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasHigh_RW_Effi_Corr = nullptr;

  // gammas reweighted highest clus Sideband subtracted
  TH1D* hNCell_GammasHigh_RW_SBSub_Effi_data = nullptr;
  TH1D* hNCell_GammasHigh_RW_SBSub_Effi_MC = nullptr;
  TH1D* hNCell_GammasHigh_RW_SBSub_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasHigh_RW_SBSub_Effi_Corr = nullptr;

  // gammas reweighted highest clus Sideband subtracted
  TH1D* hNCell_GammasHigh_RW_MCSBSub_Effi_data = nullptr;
  TH1D* hNCell_GammasHigh_RW_MCSBSub_Effi_MC = nullptr;
  TH1D* hNCell_GammasHigh_RW_MCSBSub_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasHigh_RW_MCSBSub_Effi_Corr = nullptr;

  // gammas reweighted both clus Sideband subtracted
  TH1D* hNCell_GammasWide_RW_SBSub_Effi_data = nullptr;
  TH1D* hNCell_GammasWide_RW_SBSub_Effi_MC = nullptr;
  TH1D* hNCell_GammasWide_RW_SBSub_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasWide_RW_SBSub_Effi_Corr = nullptr;

  // gammas both clus Sideband subtracted
  TH1D* hNCell_GammasWide_SBSub_Effi_data = nullptr;
  TH1D* hNCell_GammasWide_SBSub_Effi_MC = nullptr;
  TH1D* hNCell_GammasWide_SBSub_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasWide_SBSub_Effi_Corr = nullptr;

  // gammas reweighted both clus Sideband subtracted SMOOTHED
  TH1D* hNCell_GammasWide_RW_SBSub_Smoothed_Effi_data = nullptr;
  TH1D* hNCell_GammasWide_RW_SBSub_Smoothed_Effi_MC = nullptr;
  TH1D* hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Corr = nullptr;

  // electron histograms
  TH1D* hNCell_Electrons_Effi_data = nullptr;
  TH1D* hNCell_Electrons_Effi_MC = nullptr;
  TH1D* hNCell_Electrons_Effi_Ratio = nullptr;
  TH1D* hNCell_Electrons_Effi_Corr = nullptr;



  TH1D* hGammaPurity = nullptr;
  TH1D* hGammaPurityLeft = nullptr;
  TH1D* hGammaPurityWide = nullptr;
  TH1D* hGammaPurityHigh = nullptr;
  TH1D* hGammaPuritySB = nullptr;
  TH1D* hElecPurity = nullptr;
  TH1D* hElecPurityLeft = nullptr;
  TH1D* hElecPurityWide = nullptr;
  TH1D* hElecPurityHigh = nullptr;
  TH1D* hElecPuritySB = nullptr;
  TH1D* hUnPurity = nullptr;
  TH1D* hUnPurityLeft = nullptr;
  TH1D* hUnPurityWide = nullptr;
  TH1D* hUnPurityHigh = nullptr;
  TH1D* hUnPuritySB = nullptr;

  TH1D* hTrueGammas = nullptr;
  TH1D* hTrueElec = nullptr;
  TH1D* hTrueGammasLeft = nullptr;
  TH1D* hTrueElecLeft = nullptr;
  TH1D* hTrueGammasWide = nullptr;
  TH1D* hTrueElecWide = nullptr;
  TH1D* hTrueGammasHigh = nullptr;
  TH1D* hTrueElecHigh = nullptr;
  TH1D* hTrueGammasSB = nullptr;
  TH1D* hTrueElecSB = nullptr;


  // InvMass
  TH2F* hInvMassVsPt_MC = nullptr;
  TH2F* hInvMassVsPtGG = nullptr;
  TH2F* hInvMassVsPtGC = nullptr;
  TH2F* hInvMassVsPtCC = nullptr;
  TH2F* hInvMassVsPt_data = nullptr;
  TH2F* hInvMassVsPtBack_data = nullptr;
  TH2F* hInvMassVsPtBack_MC = nullptr;

  TH1D* hInvMassVsPt_Slice = nullptr;
  TH1D* hInvMassVsPtGG_Slice = nullptr;
  TH1D* hInvMassVsPtGC_Slice = nullptr;
  TH1D* hInvMassVsPtCC_Slice = nullptr;
  TH1D* hInvMassBackVsPt_Slice = nullptr;
  TH1D* hInvMassVsPt_data_Slice = nullptr;

  TH2F* hInvMassVsPtHigh_MC = nullptr;
  TH2F* hInvMassVsPtHighGamma_MC = nullptr;
  TH2F* hInvMassVsPtHighElec_MC = nullptr;
  TH2F* hInvMassVsPtHigh1cell_MC = nullptr;
  TH2F* hInvMassVsPtHigh2cell_MC = nullptr;
  TH2F* hInvMassVsPtHigh3cell_MC = nullptr;
  TH2F* hInvMassVsPtHigh_data = nullptr;
  TH2F* hInvMassVsPtHigh1cell_data = nullptr;
  TH2F* hInvMassVsPtHigh2cell_data = nullptr;
  TH2F* hInvMassVsPtHigh3cell_data = nullptr;

  // for those histos, low is loaded and high added to result in wide
  TH2F* hInvMassVsPtWide_MC = nullptr;
  TH2F* hInvMassVsPtWideGamma_MC = nullptr;
  TH2F* hInvMassVsPtWideElec_MC = nullptr;
  TH2F* hInvMassVsPtWide_data = nullptr;


  TH2F* hInvMassVsHighGammaPtBack_data = nullptr;
  TH2F* hInvMassVsHighGammaPtBack_MC = nullptr;

  TH2F* hInvMassVsWideGammaPtBack_data = nullptr;
  TH2F* hInvMassVsWideGammaPtBack_MC = nullptr;

  TH1D* hInvMassVsPtHigh_Slice_MC = nullptr;
  TH1D* hInvMassVsPtHighGamma_Slice_MC = nullptr;
  TH1D* hInvMassVsPtHighElec_Slice_MC = nullptr;
  TH1D* hInvMassVsPtHigh1cell_Slice_MC = nullptr;
  TH1D* hInvMassVsPtHigh2cell_Slice_MC = nullptr;
  TH1D* hInvMassVsPtHigh3cell_Slice_MC = nullptr;
  TH1D* hInvMassVsPtHigh_Slice_data = nullptr;
  TH1D* hInvMassVsPtHigh1cell_Slice_data = nullptr;
  TH1D* hInvMassVsPtHigh2cell_Slice_data = nullptr;
  TH1D* hInvMassVsPtHigh3cell_Slice_data = nullptr;


  // PlottingHists
  TH2F* hDummyEffi = nullptr;
  TH2F* hDummyRatio = nullptr;
  TH2F* hDummyCorr = nullptr;
  TH2F* hDummyPurity = nullptr;
  TH2F* hDummyMInv = nullptr;
  TH2F* hDummyNCellRatio = nullptr;
  TH2F* hDummyNCellFraction = nullptr;

  TCanvas* hDummyCan = nullptr;

  Float_t textSizeSinglePad = 0.044;
  Float_t textSizeLabelsRel = 0.044;

  const int nBinsE = 20;
  Double_t arrEbins[21] = {0.0, 0.7, 0.8, 0.9, 1.0,   1.2, 1.4, 1.6, 1.8, 2.0,
                               2.4, 2.8, 3.2, 3.6, 4.0,   4.5, 5.0, 6.0, 7.0, 8.0,
                               10.};


};



//----------------------------
Effi::Effi(){

}
//----------------------------
Effi::~Effi(){

}
//----------------------------
Effi::Effi(TString input, TString Period, TString suffix){
  fdata = TFile::Open(input);
  fPeriod = Period;
  fsuffix = suffix;
  fcontrolOut = new TFile("fcontrolOut.root", "Recreate");
}

void Effi::LoadTB(){

  TFile TB("TBEffi.root");
  TGraph* tmpTBData = (TGraph*) TB.Get("gData1X_Ecall100MeV");
  TGraph* tmpTBMC = (TGraph*) TB.Get("gMCe1_G3_Eell100MeV");

  double *xData = tmpTBData->GetX();
  double *yData = tmpTBData->GetY();
  double *xMC = tmpTBMC->GetX();
  double *yMC = tmpTBMC->GetY();
  Double_t ey[8] = {0.04, 0.03, 0.03, 0.02, 0.015, 0.007, 0.001, 0.001};
  Double_t ex[8] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
  Double_t eMC[8] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001};


  grTB_data = new TGraphErrors(8, xData, yData, ex, ey);
  grTB_MC = new TGraphErrors(8, xMC, yMC, eMC, eMC);

  for(int i = 0; i < grTB_data->GetN(); ++i){
    grTB_data->SetPointError(i, 0.001, ey[i]);
  }



  grTB_Ratio = (TGraphErrors*) grTB_MC->Clone("grTBRatio");

  // SetPointError

  for(int i = 0; i < grTB_Ratio->GetN(); ++i){
    grTB_Ratio->SetPointY(i, grTB_data->Eval(grTB_MC->GetPointX(i)) / grTB_MC->GetPointY(i));
    grTB_Ratio->SetPointError(i, 0.001, ey[i]/grTB_MC->GetPointY(i));
  }

  grTB_Corr = (TGraphErrors*) grTB_MC->Clone("grTBEffi");
  for(int i = 0; i < grTB_data->GetN(); ++i){
    float CF = -1;
    float CFerr = 0.001;
    if(grTB_MC->GetPointY(i) < 1){
      CF = 1 - (( 1 - grTB_data->Eval(grTB_MC->GetPointX(i)))/(1 - grTB_MC->GetPointY(i)));
      CFerr = -ey[i]/(1-grTB_MC->GetPointY(i));
    } else{
      CF = -1;
      CFerr = 0.001;
    }
    grTB_Corr->SetPointY(i, CF);
    grTB_Corr->SetPointError(i, 0.0001, CFerr);
  }
}

void Effi::FillHistos(){
  hNCellVsETMNL_data = (TH2F*) fdata->Get("hNCellVsETMNL_data");
  hNCellVsETMNL_data ->Sumw2();
  hNCellVsETMNL_data ->SetDirectory(0);
  hNCellVsETMNL_MC = (TH2F*) fdata->Get("hNCellVsETMNL_MC");
  hNCellVsETMNL_MC->Sumw2();
  hNCellVsETMNL_MC->SetDirectory(0);

  hNCellVsETMNLS500A100_data = (TH2F*) fdata->Get("hNCellVsETMNLS500A100_data");
  hNCellVsETMNLS500A100_data ->Sumw2();
  hNCellVsETMNLS500A100_data ->SetDirectory(0);
  hNCellVsETMNLS500A100_MC = (TH2F*) fdata->Get("hNCellVsETMNLS500A100_MC");
  hNCellVsETMNLS500A100_MC ->Sumw2();
  hNCellVsETMNLS500A100_MC ->SetDirectory(0);

  hNCellVsETMNLS500A105_data = (TH2F*) fdata->Get("hNCellVsETMNLS500A105_data");
  hNCellVsETMNLS500A105_data ->Sumw2();
  hNCellVsETMNLS500A105_data ->SetDirectory(0);
  hNCellVsETMNLS500A105_MC = (TH2F*) fdata->Get("hNCellVsETMNLS500A105_MC");
  hNCellVsETMNLS500A105_MC->Sumw2();
  hNCellVsETMNLS500A105_MC->SetDirectory(0);

  hNCellVsETMNLS500A110_data = (TH2F*) fdata->Get("hNCellVsETMNLS500A110_data");
  hNCellVsETMNLS500A110_data ->Sumw2();
  hNCellVsETMNLS500A110_data ->SetDirectory(0);
  hNCellVsETMNLS500A110_MC = (TH2F*) fdata->Get("hNCellVsETMNLS500A110_MC");
  hNCellVsETMNLS500A110_MC ->Sumw2();
  hNCellVsETMNLS500A110_MC ->SetDirectory(0);

  hNCellVsEGammasNL_data = (TH2F*) fdata->Get("hNCellVsEGammasNL_data");
  hNCellVsEGammasNL_data ->Sumw2();
  hNCellVsEGammasNL_data ->SetDirectory(0);
  hNCellVsEGammasNL_MC = (TH2F*) fdata->Get("hNCellVsEGammasNL_MC");
  hNCellVsEGammasNL_MC ->Sumw2();
  hNCellVsEGammasNL_MC ->SetDirectory(0);
  hNCellVsETrueGammasNL_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNL_MC");
  hNCellVsETrueGammasNL_MC ->Sumw2();
  hNCellVsETrueGammasNL_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElec_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElec_MC");
  hNCellVsEGammasNLTrueElec_MC ->Sumw2();
  hNCellVsEGammasNLTrueElec_MC ->SetDirectory(0);

  hNCellVsEGammasNLLeft_data = (TH2F*) fdata->Get("hNCellVsEGammasNLLeft_data");
  hNCellVsEGammasNLLeft_data ->Sumw2();
  hNCellVsEGammasNLLeft_data ->SetDirectory(0);
  hNCellVsEGammasNLLeft_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLLeft_MC");
  hNCellVsEGammasNLLeft_MC ->Sumw2();
  hNCellVsEGammasNLLeft_MC ->SetDirectory(0);
  hNCellVsETrueGammasNLLeft_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLLeft_MC");
  hNCellVsETrueGammasNLLeft_MC ->Sumw2();
  hNCellVsETrueGammasNLLeft_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElecLeft_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecLeft_MC");
  hNCellVsEGammasNLTrueElecLeft_MC ->Sumw2();
  hNCellVsEGammasNLTrueElecLeft_MC ->SetDirectory(0);

  hNCellVsEGammasNLWide_data = (TH2F*) fdata->Get("hNCellVsEGammasNLWide_data");
  hNCellVsEGammasNLWide_data ->Sumw2();
  hNCellVsEGammasNLWide_data ->SetDirectory(0);
  hNCellVsEGammasNLWide_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLWide_MC");
  hNCellVsEGammasNLWide_MC->Sumw2();
  hNCellVsEGammasNLWide_MC->SetDirectory(0);
  hNCellVsETrueGammasNLWide_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLWide_MC");
  hNCellVsETrueGammasNLWide_MC ->Sumw2();
  hNCellVsETrueGammasNLWide_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElecWide_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecWide_MC");
  hNCellVsEGammasNLTrueElecWide_MC ->Sumw2();
  hNCellVsEGammasNLTrueElecWide_MC ->SetDirectory(0);

  hNCellVsEGammasNLSB_data = (TH2F*) fdata->Get("hNCellVsEGammasNLSideBand_data");
  hNCellVsEGammasNLSB_data ->Sumw2();
  hNCellVsEGammasNLSB_data ->SetDirectory(0);
  hNCellVsEGammasNLSB_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLSideBand_MC");
  hNCellVsEGammasNLSB_MC ->Sumw2();
  hNCellVsEGammasNLSB_MC ->SetDirectory(0);
  hNCellVsETrueGammasNLSB_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLSideBand_MC");
  hNCellVsETrueGammasNLSB_MC ->Sumw2();
  hNCellVsETrueGammasNLSB_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElecSB_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecSideBand_MC");
  hNCellVsEGammasNLTrueElecSB_MC->Sumw2();
  hNCellVsEGammasNLTrueElecSB_MC->SetDirectory(0);

  hNCellVsEGammasNLSBHigh_data = (TH2F*) fdata->Get("hNCellVsEGammasNLSideBandOnlyHighClus_data");
  hNCellVsEGammasNLSBHigh_data ->Sumw2();
  hNCellVsEGammasNLSBHigh_data ->SetDirectory(0);
  hNCellVsEGammasNLSBHigh_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLSideBandOnlyHighClus_MC");
  hNCellVsEGammasNLSBHigh_MC ->Sumw2();
  hNCellVsEGammasNLSBHigh_MC ->SetDirectory(0);
  hNCellVsETrueGammasNLSBHigh_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLSideBandOnlyHighClus_MC");
  hNCellVsETrueGammasNLSBHigh_MC ->Sumw2();
  hNCellVsETrueGammasNLSBHigh_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElecSBHigh_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecSideBandOnlyHighClus_MC");
  hNCellVsEGammasNLTrueElecSBHigh_MC ->Sumw2();
  hNCellVsEGammasNLTrueElecSBHigh_MC ->SetDirectory(0);

  hNCellVsEGammasNLHigh_data = (TH2F*) fdata->Get("hNCellVsEGammasNLOnlyHighClus_data");
  hNCellVsEGammasNLHigh_data ->Sumw2();
  hNCellVsEGammasNLHigh_data ->SetDirectory(0);
  hNCellVsEGammasNLHigh_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLOnlyHighClus_MC");
  hNCellVsEGammasNLHigh_MC ->Sumw2();
  hNCellVsEGammasNLHigh_MC ->SetDirectory(0);
  hNCellVsETrueGammasNLHigh_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLOnlyHighClus_MC");
  hNCellVsETrueGammasNLHigh_MC ->Sumw2();
  hNCellVsETrueGammasNLHigh_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElecHigh_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecOnlyHighClus_MC");
  hNCellVsEGammasNLTrueElecHigh_MC->Sumw2();
  hNCellVsEGammasNLTrueElecHigh_MC->SetDirectory(0);

  hNCellVsEGammasNLLow_data = (TH2F*) fdata->Get("hNCellVsEGammasNLOnlyLowClus_data");
  hNCellVsEGammasNLLow_data ->Sumw2();
  hNCellVsEGammasNLLow_data ->SetDirectory(0);
  hNCellVsEGammasNLLow_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLOnlyLowClus_MC");
  hNCellVsEGammasNLLow_MC ->Sumw2();
  hNCellVsEGammasNLLow_MC ->SetDirectory(0);
  hNCellVsETrueGammasNLLow_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLOnlyLowClus_MC");
  hNCellVsETrueGammasNLLow_MC ->Sumw2();
  hNCellVsETrueGammasNLLow_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElecLow_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecOnlyLowClus_MC");
  hNCellVsEGammasNLTrueElecLow_MC ->Sumw2();
  hNCellVsEGammasNLTrueElecLow_MC ->SetDirectory(0);



  hNCellVsEelecNL_data = (TH2F*) fdata->Get("hNCellVsEelecNL_data");
  hNCellVsEelecNL_data ->Sumw2();
  hNCellVsEelecNL_data ->SetDirectory(0);
  hNCellVsEelecNL_MC = (TH2F*) fdata->Get("hNCellVsTrueEelecNL_MC");
  hNCellVsEelecNL_MC ->Sumw2();
  hNCellVsEelecNL_MC ->SetDirectory(0);


  // int reb = 4;
  // hNCellVsETMNL_data->RebinX(reb);
  // hNCellVsETMNL_MC->RebinX(reb);
  //
  // hNCellVsEGammasNL_data->RebinX(reb);
  // hNCellVsEGammasNL_MC->RebinX(reb);
  // hNCellVsETrueGammasNL_MC->RebinX(reb);
  // hNCellVsEGammasNLTrueElec_MC->RebinX(reb);
  //
  // hNCellVsEGammasNLLeft_data->RebinX(reb);
  // hNCellVsEGammasNLLeft_MC->RebinX(reb);
  // hNCellVsETrueGammasNLLeft_MC->RebinX(reb);
  // hNCellVsEGammasNLTrueElecLeft_MC->RebinX(reb);

}


// -------------------------------------------------
// Steering Efficiency histos calculation
// -------------------------------------------------
void Effi::FillCorrHistos(){
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // all clusters
  GetEffiHists(hNCellVsETMNL_data, hNCellVsETMNL_MC, hNCell_AllClus_Effi_data, hNCell_AllClus_Effi_MC, hNCell_AllClus_Effi_Ratio, hNCell_AllClus_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas
  GetEffiHists(hNCellVsEGammasNL_data, hNCellVsEGammasNL_MC, hNCell_Gammas_Effi_data, hNCell_Gammas_Effi_MC, hNCell_Gammas_Effi_Ratio, hNCell_Gammas_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas reweighted
  GetEffiHists(hNCellVsEGammasNL_RW_data, hNCellVsETrueGammasNL_MC, hNCell_Gammas_RW_Effi_data, hNCell_Gammas_RW_Effi_MC, hNCell_Gammas_RW_Effi_Ratio, hNCell_Gammas_RW_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // gammas reweighted MC
  GetEffiHists(hNCellVsEGammasNL_RW_MC, hNCellVsETrueGammasNL_MC, hNCell_Gammas_RW_Effi_MC_consistency, hNCell_Gammas_RW_Effi_MC_consistency, hNCell_Gammas_RW_Effi_Ratio_consistency, hNCell_Gammas_RW_Effi_Corr_consistency);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // gammas left
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  GetEffiHists(hNCellVsEGammasNLLeft_data, hNCellVsEGammasNLLeft_MC, hNCell_GammasLeft_Effi_data, hNCell_GammasLeft_Effi_MC, hNCell_GammasLeft_Effi_Ratio, hNCell_GammasLeft_Effi_Corr);
  // gammas left reweighted
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  GetEffiHists(hNCellVsEGammasNLLeft_RW_data, hNCellVsETrueGammasNLLeft_MC, hNCell_GammasLeft_RW_Effi_data, hNCell_GammasLeft_RW_Effi_MC, hNCell_GammasLeft_RW_Effi_Ratio, hNCell_GammasLeft_RW_Effi_Corr);

  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas Wide
  // should be used for corrections on gammas + electrons

  GetEffiHists(hNCellVsEGammasNLWide_data, hNCellVsEGammasNLWide_MC, hNCell_GammasWide_Smoothed_Effi_data, hNCell_GammasWide_Smoothed_Effi_MC, hNCell_GammasWide_Smoothed_Effi_Ratio, hNCell_GammasWide_Smoothed_Effi_Corr, true);
  GetEffiHists(hNCellVsEGammasNLWide_data, hNCellVsEGammasNLWide_MC, hNCell_GammasWide_Effi_data, hNCell_GammasWide_Effi_MC, hNCell_GammasWide_Effi_Ratio, hNCell_GammasWide_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas wide reweighted
  GetEffiHists(hNCellVsEGammasNLWide_RW_data, /*hNCellVsETrueGammasNLWide_MC*/hNCellVsEGammasNLWide_RW_MC, hNCell_GammasWide_RW_Effi_data, hNCell_GammasWide_RW_Effi_MC, hNCell_GammasWide_RW_Effi_Ratio, hNCell_GammasWide_RW_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // conversions wide reweighted
  GetEffiHists(hNCellVsEConversionsNLWide_RW_data, hNCellVsEGammasNLTrueElecWide_MC, hNCell_ConversionsWide_RW_Effi_data, hNCell_ConversionsWide_RW_Effi_MC, hNCell_ConversionsWide_RW_Effi_Ratio, hNCell_ConversionsWide_RW_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;


  // gammas Side band
  GetEffiHists(hNCellVsEGammasNLSB_data, hNCellVsEGammasNLSB_MC, hNCell_GammasSB_Effi_data, hNCell_GammasSB_Effi_MC, hNCell_GammasSB_Effi_Ratio, hNCell_GammasSB_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas Side band reweighted
  GetEffiHists(hNCellVsEGammasNLSB_RW_data, hNCellVsETrueGammasNLSB_MC, hNCell_GammasSB_RW_Effi_data, hNCell_GammasSB_RW_Effi_MC, hNCell_GammasSB_RW_Effi_Ratio, hNCell_GammasSB_RW_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // gammas Side band highest cluster
  GetEffiHists(hNCellVsEGammasNLSBHigh_data, hNCellVsEGammasNLSBHigh_MC, hNCell_GammasSBHigh_Effi_data, hNCell_GammasSBHigh_Effi_MC, hNCell_GammasSBHigh_Effi_Ratio, hNCell_GammasSBHigh_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas Side band reweighted highest cluster
  GetEffiHists(hNCellVsEGammasNLSBHigh_RW_data, hNCellVsETrueGammasNLSBHigh_MC, hNCell_GammasSBHigh_RW_Effi_data, hNCell_GammasSBHigh_RW_Effi_MC, hNCell_GammasSBHigh_RW_Effi_Ratio, hNCell_GammasSBHigh_RW_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // gammas High
  GetEffiHists(hNCellVsEGammasNLHigh_data, hNCellVsEGammasNLHigh_MC, hNCell_GammasHigh_Effi_data, hNCell_GammasHigh_Effi_MC, hNCell_GammasHigh_Effi_Ratio, hNCell_GammasHigh_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas High reweighted
  GetEffiHists(hNCellVsEGammasNLHigh_RW_data, hNCellVsETrueGammasNLHigh_MC, hNCell_GammasHigh_RW_Effi_data, hNCell_GammasHigh_RW_Effi_MC, hNCell_GammasHigh_RW_Effi_Ratio, hNCell_GammasHigh_RW_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas High reweighted sideband subtracted
  GetEffiHists(hNCellVsEGammasNLHigh_RW_SBSub_data, hNCellVsETrueGammasNLHigh_MC, hNCell_GammasHigh_RW_SBSub_Effi_data, hNCell_GammasHigh_RW_SBSub_Effi_MC, hNCell_GammasHigh_RW_SBSub_Effi_Ratio, hNCell_GammasHigh_RW_SBSub_Effi_Corr);
  // GetEffiHists(hNCellVsEGammasNLHigh_RW_SBSub_data, hNCellVsEGammasNLHigh_RW_SBSub_MC, hNCell_GammasHigh_RW_SBSub_Effi_data, hNCell_GammasHigh_RW_SBSub_Effi_MC, hNCell_GammasHigh_RW_SBSub_Effi_Ratio, hNCell_GammasHigh_RW_SBSub_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas High reweighted  MC sideband subtracted
  GetEffiHists(hNCellVsEGammasNLHigh_RW_MCSBSub_data, hNCellVsEGammasNLHigh_RW_SBSub_MC, hNCell_GammasHigh_RW_MCSBSub_Effi_data, hNCell_GammasHigh_RW_MCSBSub_Effi_MC, hNCell_GammasHigh_RW_MCSBSub_Effi_Ratio, hNCell_GammasHigh_RW_MCSBSub_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // gammas low and high wide reweighted sideband subtracted
  GetEffiHists(hNCellVsEGammasNLWide_RW_SBSub_data, /*hNCellVsETrueGammasNLWide_MC*/hNCellVsEGammasNLWide_RW_SBSub_MC, hNCell_GammasWide_RW_SBSub_Effi_data, hNCell_GammasWide_RW_SBSub_Effi_MC, hNCell_GammasWide_RW_SBSub_Effi_Ratio, hNCell_GammasWide_RW_SBSub_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas wide only sideband subtracted
  GetEffiHists(hNCellVsEGammasNLWide_SBSub_data, hNCellVsEGammasNLWide_SBSub_MC, hNCell_GammasWide_SBSub_Effi_data, hNCell_GammasWide_SBSub_Effi_MC, hNCell_GammasWide_SBSub_Effi_Ratio, hNCell_GammasWide_SBSub_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
  // gammas low and high wide reweighted sideband subtracted SMOOTHED
  GetEffiHists(hNCellVsEGammasNLWide_RW_SBSub_data, hNCellVsETrueGammasNLWide_MC, hNCell_GammasWide_RW_SBSub_Smoothed_Effi_data, hNCell_GammasWide_RW_SBSub_Smoothed_Effi_MC, hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Ratio, hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Corr, true);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // gammas Low
  GetEffiHists(hNCellVsEGammasNLLow_data, hNCellVsEGammasNLLow_MC, hNCell_GammasLow_Effi_data, hNCell_GammasLow_Effi_MC, hNCell_GammasLow_Effi_Ratio, hNCell_GammasLow_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // electron hists
  GetEffiHists(hNCellVsEelecNL_data, hNCellVsEelecNL_MC, hNCell_Electrons_Effi_data, hNCell_Electrons_Effi_MC, hNCell_Electrons_Effi_Ratio, hNCell_Electrons_Effi_Corr);
  cout<<"FillCorrHistos: "<<__LINE__<<endl;


  // true gamma in pi0 range?
  GetTrueHists(hNCellVsETrueGammasNL_MC, hNCell_TrueGammas_Effi_MC, "trueGammas");
  cout<<"FillCorrHistos: "<<__LINE__<<endl;

  // true elec in pi0 range?
  GetTrueHists(hNCellVsEGammasNLTrueElec_MC, hNCell_TrueElec_Effi_MC, "trueElecs");
  cout<<"FillCorrHistos: "<<__LINE__<<endl;
}


// -------------------------------------------------
// Steering for sideband subtraction
// -------------------------------------------------
void Effi::DoSBSubtraction(){
  // data high gamma
  cout<<"hInvMassVsPtHigh_data->Integrl(): "<<hNCellVsEGammasNLSBHigh_data->Integral()<<endl;
  fcontrolOut->cd();
  hNCellVsEGammasNLSBHigh_data->Write();
  hNCellVsEGammasNLHigh_data->Write();
  hNCellVsEGammasNLHigh_SBSub_data = SubtractSidebandBack( hInvMassVsPtHigh_data, hInvMassVsHighGammaPtBack_data , hNCellVsEGammasNLSBHigh_data , hNCellVsEGammasNLHigh_data);
  // hNCellVsEGammasNLHigh_SBSub_data = (TH2F*) hNCellVsEGammasNLHigh_data->Clone("hNCellVsEGammasNLHigh_SBSub_data");//SubtractSidebandBack( hInvMassVsPtHigh_data, hInvMassVsHighGammaPtBack_data , hNCellVsEGammasNLSBHigh_data , hNCellVsEGammasNLHigh_data);
  // MC high gamma
  cout<<"hInvMassVsPtHigh_MC->Integrl(): "<<hNCellVsEGammasNLSBHigh_MC->Integral()<<endl;
  hInvMassVsPtHigh_MC->Write();


  hNCellVsEGammasNLHigh_SBSub_MC = SubtractSidebandBack( hInvMassVsPtHigh_MC, hInvMassVsHighGammaPtBack_MC , hNCellVsEGammasNLSBHigh_MC , hNCellVsEGammasNLHigh_MC);
  // hNCellVsEGammasNLHigh_SBSub_MC = (TH2F*) hNCellVsEGammasNLHigh_MC->Clone("hNCellVsEGammasNLHigh_SBSub_MC");
  // data high gamma
  hNCellVsEGammasNLHigh_MCSBSub_data = SubtractSidebandBack( hInvMassVsPtHigh_MC, hInvMassVsHighGammaPtBack_MC , hNCellVsEGammasNLSBHigh_MC , hNCellVsEGammasNLHigh_MC);
  // data wide gamma
  hNCellVsEGammasNLWide_SBSub_data = SubtractSidebandBack( hInvMassVsPtWide_data, hInvMassVsWideGammaPtBack_data , hNCellVsEGammasNLSB_data , hNCellVsEGammasNLWide_data);
  // MC wide gamma
  hNCellVsEGammasNLWide_SBSub_MC = SubtractSidebandBack( hInvMassVsPtWide_MC, hInvMassVsWideGammaPtBack_MC , hNCellVsEGammasNLSB_MC , hNCellVsEGammasNLWide_MC);

}

// -------------------------------------------------
// SUBTRACT SIDEBAND FROM SIGNAL REGION
// Make Sure the sideband and the peak use the same method or gamma selection!
// -------------------------------------------------
TH2F* Effi::SubtractSidebandBack(TH2F* hMInv, TH2F* hMInvBack, TH2F* hNCellSB, TH2F* hNCell){


  // STrategy:
  // Fill a histo with a pT dependent scaling factors
  // this only works if the inv mass histo is only filled with the pT of the higher photon!
  TFile fMassPos("fMassPos.root");
  TF1 *funcMassPos = nullptr;
  if(fMassPos.IsOpen()){
    TString name = "MassPos8TeV_func";
    if(fPeriod.Contains("13TeVLowB")) name = "MassPos13TeVLowB_func";
    else if(fPeriod.Contains("13TeVNomB")) name = "MassPos13TeVNomB_func";
    else if(fPeriod.Contains("8TeV")) name = "MassPos8TeV_func";
    funcMassPos = (TF1*) fMassPos.Get(name);
  }

  TH2F* hScale = (TH2F*) hNCell->Clone("hScale");
  hScale->Reset();
  for(int i = 1; i <= hNCell->GetNbinsX(); ++i){
    float xValLow = hNCell->GetXaxis()->GetBinLowEdge(i) + 0.001;
    float xValUp = hNCell->GetXaxis()->GetBinUpEdge(i) - 0.001;

    // get inv. mass distribution for signal and background
    // scale background to same evt.
    TH1D* hMInv_Proj = (TH1D*) hMInv->ProjectionX("hMInv_Proj", hMInv->GetYaxis()->FindBin(xValLow), hMInv->GetYaxis()->FindBin(xValUp));
    TH1D* hMInvBack_Proj = (TH1D*) hMInvBack->ProjectionX("hMInvBack_Proj", hMInvBack->GetYaxis()->FindBin(xValLow), hMInvBack->GetYaxis()->FindBin(xValUp));
    cout<<"hMInv_Proj->GetBinContent(20): "<<hMInv_Proj->GetBinContent(20)<<endl;
    float scaleFac = hMInv_Proj->Integral(hMInv_Proj->FindBin(0.17), hMInv_Proj->FindBin(0.3)) / hMInvBack_Proj->Integral(hMInvBack_Proj->FindBin(0.17), hMInvBack_Proj->FindBin(0.3));
    hMInvBack_Proj->Scale(scaleFac);

    // ranges for background and signal make that better??
    float rangeSB[2] = {0.135 + 0.05, 0.135 + 0.2};
    float rangeSignal[2] = {0.08, 0.15};
    cout<<"\n\n Before Mass Pos \n\n";
    if(funcMassPos){
      cout<<"\n\n In Mass Pos \n\n";
      rangeSB[0] = funcMassPos->Eval(0.5*(xValLow + xValUp )) + 0.05;
      rangeSB[1] = funcMassPos->Eval(0.5*(xValLow + xValUp )) + 0.2;
      rangeSignal[0] = funcMassPos->Eval(0.5*(xValLow + xValUp )) - 0.0;
      rangeSignal[1] = funcMassPos->Eval(0.5*(xValLow + xValUp )) + 0.02;
    }

    // calculate the integral of the background in signal range to estimate the amount of background
    // Do the same for the same evt. in the sideband region
    // calculate the ratio between the background and the sideband to later scale the sideband contribution
    float backInt = hMInvBack_Proj->Integral(hMInvBack_Proj->FindBin(rangeSignal[0]), hMInvBack_Proj->FindBin(rangeSignal[1]));
    float SBInt = hMInv_Proj->Integral(hMInv_Proj->FindBin(rangeSB[0]), hMInv_Proj->FindBin(rangeSB[1]));
    cout<<"hMInv_Proj: "<<hMInvBack_Proj->GetBinContent(20)<<endl;
    cout<<"backInt: "<<backInt<<"      SBInt: "<<SBInt<<"   backInt/SBInt: "<<backInt/SBInt<<endl;
    for(int iy = 1; iy < hScale->GetNbinsY(); ++iy){
      if(isnan(backInt/SBInt)){
        hScale->SetBinContent(i, iy, 1);
      }else {
        hScale->SetBinContent(i, iy, backInt/SBInt);
      }
    }
  }



  // Scale background to same evt.

  cout<<__LINE__<<endl;
  // scale NCell vs E distribution in the sideband region to the Signal region
  // this is now the scaled NCell vs. E sideband contribution that has to be subtracted from the NCell vs E distribution in the signal region
  hNCellSB->Multiply(hScale);
  cout<<"hNCellSB->Integral(): "<<hNCellSB->Integral()<<endl;
  fcontrolOut->cd();
  hScale->Write();
  hNCellSB->Write(Form("%s_control",hNCellSB->GetName()));

  // subtract SB NCell distribution from Signal NCell distribution
  // This will also remove some photons!
  TH2F* hSBSubtracted = (TH2F*) hNCell->Clone("hSBSubtracted");
  hSBSubtracted->SetDirectory(0);
  hSBSubtracted->Add(hNCellSB, -1);
  cout<<__LINE__<<endl;
  return hSBSubtracted;

}

void Effi::SetOtherHistos(){

  // right side of peak
  hGammaPurity = (TH1D*) hNCellVsETrueGammasNL_MC->ProjectionX("hGammaPurity");
  hTrueGammas = (TH1D*) hNCellVsETrueGammasNL_MC->ProjectionX("hTrueGammas");
  hTrueElec = (TH1D*) hNCellVsEGammasNLTrueElec_MC->ProjectionX("hTrueElec");
  hElecPurity = (TH1D*) hNCellVsEGammasNLTrueElec_MC->ProjectionX("hElecPurity");
  hGammaPurity->Divide(hTrueGammas, (TH1D*) hNCellVsEGammasNL_MC->ProjectionX("hGammasForPurity"), 1, 1, "B");
  hElecPurity->Divide(hTrueElec, (TH1D*) hNCellVsEGammasNL_MC->ProjectionX("hElecForPurity"), 1, 1, "B");
  hUnPurity = (TH1D*) hGammaPurity->Clone("hUnPurity");
  hUnPurity->Add(hElecPurity);
  for(int i = 1; i <= hUnPurity->GetNbinsX(); ++i){ hUnPurity->SetBinContent(i, 1- hUnPurity->GetBinContent(i));}

  // left side of peak
  hGammaPurityLeft = (TH1D*) hNCellVsETrueGammasNLLeft_MC->ProjectionX("hGammaPurityLeft");
  hTrueGammasLeft = (TH1D*) hNCellVsETrueGammasNLLeft_MC->ProjectionX("hTrueGammasLeft");
  hTrueElecLeft = (TH1D*) hNCellVsEGammasNLTrueElecLeft_MC->ProjectionX("hTrueElecLeft");
  hElecPurityLeft = (TH1D*) hNCellVsEGammasNLTrueElecLeft_MC->ProjectionX("hElecPurityLeft");
  hGammaPurityLeft->Divide(hTrueGammasLeft, (TH1D*) hNCellVsEGammasNLLeft_MC->ProjectionX("hGammasForPurityLeft"), 1, 1, "B");
  hElecPurityLeft->Divide(hTrueElecLeft, (TH1D*) hNCellVsEGammasNLLeft_MC->ProjectionX("hElecForPurityLeft"), 1, 1, "B");

  hUnPurityLeft = (TH1D*) hGammaPurityLeft->Clone("hUnPurityLeft");
  hUnPurityLeft->Add(hElecPurityLeft);
  for(int i = 1; i <= hUnPurityLeft->GetNbinsX(); ++i){ hUnPurityLeft->SetBinContent(i, 1- hUnPurityLeft->GetBinContent(i));}

  // both sides of peak
  hGammaPurityWide = (TH1D*) hNCellVsETrueGammasNLWide_MC->ProjectionX("hGammaPurityWide");
  hTrueGammasWide = (TH1D*) hNCellVsETrueGammasNLWide_MC->ProjectionX("hTrueGammasWide");
  hTrueElecWide = (TH1D*) hNCellVsEGammasNLTrueElecWide_MC->ProjectionX("hTrueElecWide");
  hElecPurityWide = (TH1D*) hNCellVsEGammasNLTrueElecWide_MC->ProjectionX("hElecPurityWide");
  hGammaPurityWide->Divide(hTrueGammasWide, (TH1D*) hNCellVsEGammasNLWide_MC->ProjectionX("hGammasForPurityWide"), 1, 1, "B");
  hElecPurityWide->Divide(hTrueElecWide, (TH1D*) hNCellVsEGammasNLWide_MC->ProjectionX("hElecForPurityWide"), 1, 1, "B");

  hUnPurityWide = (TH1D*) hGammaPurityWide->Clone("hUnPurityWide");
  hUnPurityWide->Add(hElecPurityWide);
  for(int i = 1; i <= hUnPurityWide->GetNbinsX(); ++i){ hUnPurityWide->SetBinContent(i, 1- hUnPurityWide->GetBinContent(i));}


  // both sides of peak, high gamma
  hGammaPurityHigh = (TH1D*) hNCellVsETrueGammasNLHigh_MC->ProjectionX("hGammaPurityHigh");
  hTrueGammasHigh = (TH1D*) hNCellVsETrueGammasNLHigh_MC->ProjectionX("hTrueGammasHigh");
  hTrueElecHigh = (TH1D*) hNCellVsEGammasNLTrueElecHigh_MC->ProjectionX("hTrueElecHigh");
  hElecPurityHigh = (TH1D*) hNCellVsEGammasNLTrueElecHigh_MC->ProjectionX("hElecPurityHigh");
  hGammaPurityHigh->Divide(hTrueGammasHigh, (TH1D*) hNCellVsEGammasNLHigh_MC->ProjectionX("hGammasForPurityHigh"), 1, 1, "B");
  hElecPurityHigh->Divide(hTrueElecHigh, (TH1D*) hNCellVsEGammasNLHigh_MC->ProjectionX("hElecForPurityHigh"), 1, 1, "B");


  // get integral of signal and background (THIS IS NOT PT DEPENDENT WHICH IS NOT CORRECT)
  // TH1D* hMInv_Proj = (TH1D*) hMInv->ProjectionX("hMInv_Proj");
  // TH1D* hMInvBack_Proj = (TH1D*) hMInvBack->ProjectionX("hMInvBack_Proj");
  // float scaleFacMInv = hMInv_Proj->Integral(hMInv_Proj->FindBin(0.17), hMInv_Proj->FindBin(0.3)) / hMInvBack_Proj->Integral(hMInvBack_Proj->FindBin(0.17), hMInvBack_Proj->FindBin(0.3));
  // hMInvBack_Proj->Scale(scaleFacMInv);
  // cout<<"scaleFacMInv: "<<scaleFacMInv<<endl;
  // float intSignalBack = hMInv_Proj->Integral(hMInv_Proj->FindBin(0.05), hMInv_Proj->FindBin(0.15));
  // float intBack = hMInvBack_Proj->Integral(hMInvBack_Proj->FindBin(0.05), hMInvBack_Proj->FindBin(0.15));

  // fraction of electrons in SB


  hUnPurityHigh = (TH1D*) hGammaPurityHigh->Clone("hUnPurityHigh");
  hUnPurityHigh->Add(hElecPurityHigh);
  for(int i = 1; i <= hUnPurityHigh->GetNbinsX(); ++i){ hUnPurityHigh->SetBinContent(i, 1- hUnPurityHigh->GetBinContent(i));}



  // Sideband, SB gamma
  hGammaPuritySB = (TH1D*) hNCellVsETrueGammasNLSB_MC->ProjectionX("hGammaPuritySB");
  hTrueGammasSB = (TH1D*) hNCellVsETrueGammasNLSB_MC->ProjectionX("hTrueGammasSB");
  hTrueElecSB = (TH1D*) hNCellVsEGammasNLTrueElecSB_MC->ProjectionX("hTrueElecSB");
  hElecPuritySB = (TH1D*) hNCellVsEGammasNLTrueElecSB_MC->ProjectionX("hElecPuritySB");
  hGammaPuritySB->Divide(hTrueGammasSB, (TH1D*) hNCellVsEGammasNLSB_MC->ProjectionX("hGammasForPuritySB"), 1, 1, "B");
  hElecPuritySB->Divide(hTrueElecSB, (TH1D*) hNCellVsEGammasNLSB_MC->ProjectionX("hElecForPuritySB"), 1, 1, "B");

  hUnPuritySB = (TH1D*) hGammaPuritySB->Clone("hUnPuritySB");
  hUnPuritySB->Add(hElecPuritySB);
  for(int i = 1; i <= hUnPuritySB->GetNbinsX(); ++i){ hUnPuritySB->SetBinContent(i, 1- hUnPuritySB->GetBinContent(i));}


  hInvMassVsPt_MC = (TH2F*) fdata->Get("hInvMassVsPt_MC");
  hInvMassVsPt_MC->Sumw2();
  hInvMassVsPtBack_MC = (TH2F*) fdata->Get("hInvMassVsPtBack_MC");
  hInvMassVsPtBack_MC->Sumw2();
  hInvMassVsPt_data = (TH2F*) fdata->Get("hInvMassVsPt_data");
  hInvMassVsPt_data->Sumw2();
  hInvMassVsPtGG = (TH2F*) fdata->Get("hInvMassVsPtGG_MC");
  hInvMassVsPtGG->Sumw2();
  hInvMassVsPtGC = (TH2F*) fdata->Get("hInvMassVsPtGC_MC");
  hInvMassVsPtGC->Sumw2();
  hInvMassVsPtCC = (TH2F*) fdata->Get("hInvMassVsPtCC_MC");
  hInvMassVsPtCC->Sumw2();
  float minpt = 2.;
  float maxpt = 3.;
  hInvMassVsPt_Slice = (TH1D*) hInvMassVsPt_MC->ProjectionX("hInvMassVsPt_Slice", hInvMassVsPt_MC->GetYaxis()->FindBin(minpt), hInvMassVsPt_MC->GetYaxis()->FindBin(maxpt));
  hInvMassVsPt_data_Slice = (TH1D*) hInvMassVsPt_data->ProjectionX("hInvMassVsPt_data_Slice", hInvMassVsPt_MC->GetYaxis()->FindBin(minpt), hInvMassVsPt_MC->GetYaxis()->FindBin(maxpt));
  hInvMassVsPtGG_Slice = (TH1D*) hInvMassVsPtGG->ProjectionX("hInvMassVsPtGG_Slice", hInvMassVsPtGG->GetYaxis()->FindBin(minpt), hInvMassVsPtGG->GetYaxis()->FindBin(maxpt));
  hInvMassVsPtGC_Slice = (TH1D*) hInvMassVsPtGC->ProjectionX("hInvMassVsPtGC_Slice", hInvMassVsPtGC->GetYaxis()->FindBin(minpt), hInvMassVsPtGC->GetYaxis()->FindBin(maxpt));
  hInvMassVsPtCC_Slice = (TH1D*) hInvMassVsPtCC->ProjectionX("hInvMassVsPtCC_Slice", hInvMassVsPtCC->GetYaxis()->FindBin(minpt), hInvMassVsPtCC->GetYaxis()->FindBin(maxpt));
  hInvMassBackVsPt_Slice = (TH1D*) hInvMassVsPtBack_MC->ProjectionX("hInvMassBackVsPt_Slice", hInvMassVsPtBack_MC->GetYaxis()->FindBin(minpt), hInvMassVsPtBack_MC->GetYaxis()->FindBin(maxpt));

  // load background hists
  hInvMassVsPtBack_data = (TH2F*) fdata->Get("hInvMassVsPtBack_data");
  hInvMassVsPtBack_data->Sumw2();
  hInvMassVsPtBack_data = (TH2F*) fdata->Get("hInvMassVsPtBack_data");
  hInvMassVsPtBack_data->Sumw2();
  hInvMassVsPtBack_MC = (TH2F*) fdata->Get("hInvMassVsPtBack_MC");
  hInvMassVsPtBack_MC->Sumw2();


  hInvMassVsPtHigh_MC = (TH2F*) fdata->Get("hInvMassVsPtHigh_MC");
  hInvMassVsPtHigh_MC->Sumw2();
  hInvMassVsHighGammaPtBack_MC = (TH2F*) fdata->Get("hInvMassVsHighGammaPtBack_MC");
  hInvMassVsHighGammaPtBack_MC->Sumw2();
  hInvMassVsPtHighGamma_MC = (TH2F*) fdata->Get("hInvMassVsPtHighGamma_MC");
  hInvMassVsPtHighGamma_MC->Sumw2();
  hInvMassVsPtHighElec_MC = (TH2F*) fdata->Get("hInvMassVsPtHighElec_MC");
  hInvMassVsPtHighElec_MC ->Sumw2();
  hInvMassVsPtHigh1cell_MC = (TH2F*) fdata->Get("hInvMassVsPtHigh1cell_MC");
  hInvMassVsPtHigh1cell_MC->Sumw2();
  hInvMassVsPtHigh2cell_MC = (TH2F*) fdata->Get("hInvMassVsPtHigh2cell_MC");
  hInvMassVsPtHigh2cell_MC ->Sumw2();
  hInvMassVsPtHigh3cell_MC = (TH2F*) fdata->Get("hInvMassVsPtHigh3cell_MC");
  hInvMassVsPtHigh3cell_MC ->Sumw2();
  hInvMassVsPtHigh_data = (TH2F*) fdata->Get("hInvMassVsPtHigh_data");
  hInvMassVsPtHigh_data->Sumw2();
  hInvMassVsHighGammaPtBack_data = (TH2F*) fdata->Get("hInvMassVsHighGammaPtBack_data");
  hInvMassVsHighGammaPtBack_data->Sumw2();
  hInvMassVsPtHigh1cell_data = (TH2F*) fdata->Get("hInvMassVsPtHigh1cell_data");
  hInvMassVsPtHigh1cell_data ->Sumw2();
  hInvMassVsPtHigh2cell_data = (TH2F*) fdata->Get("hInvMassVsPtHigh2cell_data");
  hInvMassVsPtHigh2cell_data->Sumw2();
  hInvMassVsPtHigh3cell_data = (TH2F*) fdata->Get("hInvMassVsPtHigh3cell_data");
  hInvMassVsPtHigh3cell_data->Sumw2();

  // histos for low and high = wide
  hInvMassVsPtWide_MC = (TH2F*) fdata->Get("hInvMassVsPtLow_MC");
  hInvMassVsPtWide_MC->Sumw2();
  hInvMassVsWideGammaPtBack_MC = (TH2F*) fdata->Get("hInvMassVsLowGammaPtBack_MC");
  hInvMassVsWideGammaPtBack_MC ->Sumw2();
  hInvMassVsPtWideGamma_MC = (TH2F*) fdata->Get("hInvMassVsPtLowGamma_MC");
  hInvMassVsPtWideGamma_MC ->Sumw2();
  hInvMassVsPtWideElec_MC = (TH2F*) fdata->Get("hInvMassVsPtLowElec_MC");
  hInvMassVsPtWideElec_MC ->Sumw2();
  hInvMassVsPtWide_data = (TH2F*) fdata->Get("hInvMassVsPtLow_data");
  hInvMassVsPtWide_data ->Sumw2();
  hInvMassVsWideGammaPtBack_data = (TH2F*) fdata->Get("hInvMassVsLowGammaPtBack_data");
  hInvMassVsWideGammaPtBack_data->Sumw2();

  // add high stuff
  hInvMassVsPtWide_MC->Add(hInvMassVsPtHigh_MC);
  hInvMassVsWideGammaPtBack_MC->Add(hInvMassVsHighGammaPtBack_MC);
  hInvMassVsPtWideGamma_MC->Add(hInvMassVsPtHighGamma_MC);
  hInvMassVsPtWideElec_MC->Add(hInvMassVsPtHighElec_MC);
  hInvMassVsPtWide_data->Add(hInvMassVsPtHigh_data);
  hInvMassVsWideGammaPtBack_data->Add(hInvMassVsHighGammaPtBack_data);




  hInvMassVsPtHigh_Slice_MC = (TH1D*)  hInvMassVsPtHigh_MC->ProjectionX("hInvMassVsPtHigh_Slice_MC", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHighGamma_Slice_MC = (TH1D*)  hInvMassVsPtHighGamma_MC->ProjectionX("hInvMassVsPtHighGamma_Slice_MC", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHighElec_Slice_MC = (TH1D*)  hInvMassVsPtHighElec_MC->ProjectionX("hInvMassVsPtHighElec_Slice_MC", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHigh1cell_Slice_MC = (TH1D*)  hInvMassVsPtHigh1cell_MC->ProjectionX("hInvMassVsPtHigh1cell_Slice_MC", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHigh2cell_Slice_MC = (TH1D*)  hInvMassVsPtHigh2cell_MC->ProjectionX("hInvMassVsPtHigh2cell_Slice_MC", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHigh3cell_Slice_MC = (TH1D*)  hInvMassVsPtHigh3cell_MC->ProjectionX("hInvMassVsPtHigh3cell_Slice_MC", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHigh_Slice_data = (TH1D*)  hInvMassVsPtHigh_data->ProjectionX("hInvMassVsPtHigh_Slice_data", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHigh1cell_Slice_data = (TH1D*)  hInvMassVsPtHigh1cell_data->ProjectionX("hInvMassVsPtHigh1cell_Slice_data", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHigh2cell_Slice_data = (TH1D*)  hInvMassVsPtHigh2cell_data->ProjectionX("hInvMassVsPtHigh2cell_Slice_data", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));
  hInvMassVsPtHigh3cell_Slice_data = (TH1D*)  hInvMassVsPtHigh3cell_data->ProjectionX("hInvMassVsPtHigh3cell_Slice_data", hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.001), hInvMassVsPtHigh_MC->GetYaxis()->FindBin(2.999));




}


// ------------------------------------------
// Get NCell efficiency hists from 2d NCellVs E distributions
// ------------------------------------------
void Effi::GetEffiHists(TH2F *hdata2d, TH2F *hMC2d, TH1D *&hEffiData,  TH1D *&hEffiMC, TH1D *&hRatio, TH1D *&hCorr, bool doSmooth){

  TH1D* hdataNCell1 = (TH1D*) hdata2d->ProjectionX("hdataNCell1", hdata2d->GetYaxis()->FindBin(1.5), hdata2d->GetYaxis()->FindBin(20));
  TH1D* hdataNCellOnly1 = (TH1D*) hdata2d->ProjectionX("hdataNCellOnly1", hdata2d->GetYaxis()->FindBin(1.5), hdata2d->GetYaxis()->FindBin(1.5));
  TH1D* hdataNCell2 = (TH1D*) hdata2d->ProjectionX("hdataNCell2", hdata2d->GetYaxis()->FindBin(2.5), hdata2d->GetYaxis()->FindBin(20));
  TH1D* hMCNCell1 = (TH1D*) hMC2d->ProjectionX("hMCNCell1", hMC2d->GetYaxis()->FindBin(1.5), hMC2d->GetYaxis()->FindBin(20));
  TH1D* hMCNCellOnly1 = (TH1D*) hMC2d->ProjectionX("hMCNCellOnly1", hMC2d->GetYaxis()->FindBin(1.5), hMC2d->GetYaxis()->FindBin(1.5));
  TH1D* hMCNCell2 = (TH1D*) hMC2d->ProjectionX("hMCNCell2", hMC2d->GetYaxis()->FindBin(2.5), hMC2d->GetYaxis()->FindBin(20));


  hEffiData = (TH1D*) hdataNCell2->Clone("hEffidata");
  hEffiData->Divide(hdataNCell2, hdataNCell1, 1, 1, "B");

  hEffiMC = (TH1D*) hMCNCell2->Clone("hEffiMC");
  hEffiMC->Divide(hMCNCell2, hMCNCell1, 1, 1, "B");

  // smoothing on effi histos
  if(doSmooth){
    // logistisches wachstum
    // TF1 * func = new TF1("func", "[2]/(1+exp(-[0]*(x - [1]))) - [2] + 1", 0.7, 3);
    TF1 * func = new TF1("func", func1, 1., 8, 4);
    func->SetParameter(0, 2.8);
    func->SetParameter(1, 9.5);
    func->SetParameter(2, 9.9);
    func->SetParameter(3, 970);
    // TF1 * func = new TF1("func", "-(x-[0])^([1])", 0.7, 5);
    // func->SetParameters(1, 1);
    // func->SetParLimits(1, -100, -1);
    // func->SetParameters(0.1, 1);
    // func->SetParLimits(1, -100, -1);
    // fit data and take values from fit
    // TFile outfile("fNCellEffi_InPi0Range.root", "Recreate");
    // outfile.cd();
    // hEffiData->Write("data");
    // hEffiMC->Write("MC");

    hEffiData->Fit(func, "MR0");
    // func->Write("func_data");
    for(int i = 5; i <= hEffiData->GetNbinsX(); ++i){
      hEffiData->SetBinContent(i, func->Eval(hEffiData->GetBinCenter(i)));
    }

    hEffiMC->Fit(func, "MR0");
    // func->Write("func_MC");
    for(int i = 5; i <= hEffiMC->GetNbinsX(); ++i){
      hEffiMC->SetBinContent(i, func->Eval(hEffiMC->GetBinCenter(i)));
    }
    // outfile.Close();
  }


  hRatio = (TH1D*) hEffiData->Clone("hEffiRatio");
  hRatio->Divide(hEffiMC);


  // Calculate NCell Effi
  hCorr = (TH1D*) hRatio->Clone("NCellEffi");

  for(int i = 1; i <= hCorr->GetNbinsX(); ++i){
    float CF = -1;
    float err = 0;
    if(hEffiMC->GetBinContent(i) < 1){
      CF = 1 - (( 1 - hEffiData->GetBinContent(i))/(1 - hEffiMC->GetBinContent(i)));
      err = calcErr(hEffiMC->GetBinContent(i), hEffiMC->GetBinError(i), hEffiData->GetBinContent(i), hEffiData->GetBinError(i));
    }
    cout<<CF<<endl;
    hCorr->SetBinContent(i, CF);
    hCorr->SetBinError(i, err);
  }

}

// ------------------------------------------
// Error calculation using gaussian error propagation for NCell effi
// ------------------------------------------
float Effi::calcErr(float Emc, float Errmc, float Edata, float Errdata){
  float err1 = Errmc / (1-Edata) ;
  float err2 =  ((1 - Emc) * Errdata ) / ((1-Edata)*(1-Edata) );
  return sqrt(err1*err1 + err2*err2);
}

// ------------------------------------------
// GetTrueHists
// ------------------------------------------
void Effi::GetTrueHists(TH2F* h2, TH1D *&h, TString name ){
  TH1D* hNCell1 = (TH1D*) h2->ProjectionX("hNCell1", h2->GetYaxis()->FindBin(1.5), h2->GetYaxis()->FindBin(20));
  TH1D* hNCell2 = (TH1D*) h2->ProjectionX("hNCell2", h2->GetYaxis()->FindBin(2.5), h2->GetYaxis()->FindBin(20));
  h = (TH1D*) hNCell2->Clone(name);
  h->Divide(hNCell2, hNCell1, 1, 1, "B");
}



// ------------------------------------------
// Purity correction/Reweighting
// ------------------------------------------
void Effi::GetHistReweighted(){
  hNCellVsEGammasNL_RW_data = (TH2F*) hNCellVsEGammasNL_data->Clone("hNCellVsEGammasNL_RW_data");
  // get electron fraction
  TH2F* hElecFrac = (TH2F*) hNCellVsEGammasNLTrueElec_MC->Clone("hElecFrac");
  hElecFrac->Divide(hNCellVsEGammasNL_MC);
  hElecFrac->Multiply(hNCellVsEGammasNL_data);
  hNCellVsEGammasNL_RW_data->Add(hElecFrac, -1);


  hNCellVsEGammasNLLeft_RW_data = (TH2F*) hNCellVsEGammasNLLeft_data->Clone("hNCellVsEGammasNLLeft_RW_data");
  // get electron fraction
  TH2F* hElecFracLeft = (TH2F*) hNCellVsEGammasNLTrueElecLeft_MC->Clone("hElecFracLeft");
  hElecFracLeft->Divide(hNCellVsEGammasNLLeft_MC);
  hElecFracLeft->Multiply(hNCellVsEGammasNLLeft_data);
  hNCellVsEGammasNLLeft_RW_data->Add(hElecFracLeft, -1);


  hNCellVsEGammasNLWide_RW_data = (TH2F*) hNCellVsEGammasNLWide_data->Clone("hNCellVsEGammasNLWide_RW_data");
  // get electron fraction
  TH2F* hElecFracWide = (TH2F*) hNCellVsEGammasNLTrueElecWide_MC->Clone("hElecFracWide");
  hElecFracWide->Divide(hNCellVsEGammasNLWide_MC);
  hElecFracWide->Multiply(hNCellVsEGammasNLWide_data);
  hNCellVsEGammasNLWide_RW_data->Add(hElecFracWide, -1);


  hNCellVsEGammasNLWide_RW_MC = (TH2F*) hNCellVsEGammasNLWide_MC->Clone("hNCellVsEGammasNLWide_RW_MC");
  hNCellVsEGammasNLWide_RW_MC->Add(hNCellVsEGammasNLTrueElecWide_MC, -1);



  hNCellVsEGammasNLWide_RW_SBSub_data = (TH2F*) hNCellVsEGammasNLWide_SBSub_data->Clone("hNCellVsEGammasNLWide_RW_SBSub_data");
  // get electron fraction
  // was there from before
  TH2F* hElecFracWide2 = (TH2F*) hNCellVsEGammasNLTrueElecWide_MC->Clone("hElecFracWide2");
  // hElecFracWide2->Divide(hNCellVsEGammasNLWide_SBSub_MC);
  // hElecFracWide2->Multiply(hNCellVsEGammasNLWide_SBSub_data);
  hElecFracWide2->Scale(hNCellVsEGammasNLWide_SBSub_data->Integral()/hNCellVsEGammasNLWide_SBSub_MC->Integral());
  hNCellVsEGammasNLWide_RW_SBSub_data->Add(hElecFracWide2, -1);
  //
  // TFile *fout = new TFile("fStrangeHists.root", "recreate");
  // fout->cd();
  // hNCellVsEGammasNLWide_data->Write("hNCellVsEGammasNLWide_SBSub_data");
  // hNCellVsEGammasNLTrueElecWide_MC->Write("hNCellVsEGammasNLTrueElecWide_MC");
  // hNCellVsEGammasNLWide_MC->Write("hNCellVsEGammasNLWide_MC");
  // hNCellVsEGammasNLWide_RW_SBSub_data->Write("hNCellVsEGammasNLWide_RW_SBSub_data");
  // fout->Close();


  hNCellVsEGammasNLWide_RW_SBSub_MC = (TH2F*) hNCellVsEGammasNLWide_SBSub_MC->Clone("hNCellVsEGammasNLWide_RW_SBSub_MC");
  // get electron fraction
  TH2F* hElecFracWideMC = (TH2F*) hNCellVsEGammasNLTrueElecWide_MC->Clone("hElecFracWideMC");
  // hElecFracWideMC->Divide(hNCellVsEGammasNLWide_SBSub_MC);
  // hElecFracWideMC->Multiply(hNCellVsEGammasNLWide_SBSub_MC);
  hNCellVsEGammasNLWide_RW_SBSub_MC->Add(hElecFracWideMC, -1);


// reweight to get the conversion contribution
  hNCellVsEConversionsNLWide_RW_data = (TH2F*) hNCellVsEGammasNLWide_data->Clone("hNCellVsEConversionsNLWide_RW_data");
  // get photon fraction
  TH2F* hGammaFracWide = (TH2F*) hNCellVsETrueGammasNLWide_MC->Clone("hGammaFracWide");
  hGammaFracWide->Divide(hNCellVsEGammasNLWide_MC);
  hGammaFracWide->Multiply(hNCellVsEGammasNLWide_data);
  hNCellVsEConversionsNLWide_RW_data->Add(hGammaFracWide, -1);


  hNCellVsEGammasNLSB_RW_data = (TH2F*) hNCellVsEGammasNLSB_data->Clone("hNCellVsEGammasNLSB_RW_data");
  // get electron fraction
  TH2F* hElecFracSB = (TH2F*) hNCellVsEGammasNLTrueElecSB_MC->Clone("hElecFracSB");
  hElecFracSB->Divide(hNCellVsEGammasNLSB_MC);
  hElecFracSB->Multiply(hNCellVsEGammasNLSB_data);
  hNCellVsEGammasNLSB_RW_data->Add(hElecFracSB, -1);



  hNCellVsEGammasNLSBHigh_RW_data = (TH2F*) hNCellVsEGammasNLSBHigh_data->Clone("hNCellVsEGammasNLSBHigh_RW_data");
  // get electron fraction
  TH2F* hElecFracSBHigh = (TH2F*) hNCellVsEGammasNLTrueElecSBHigh_MC->Clone("hElecFracSBHigh");
  hElecFracSBHigh->Divide(hNCellVsEGammasNLSBHigh_MC);
  hElecFracSBHigh->Multiply(hNCellVsEGammasNLSBHigh_data);
  hNCellVsEGammasNLSBHigh_RW_data->Add(hElecFracSBHigh, -1);

  // ------------------------------------------------
  // Special case for Sideband subtracted case (hNCellVsEGammasNLHigh_data) has sideband contribution subtracted!!!
  // special purity correction needed! with corrected fractions

  hNCellVsEGammasNLHigh_RW_data = (TH2F*) hNCellVsEGammasNLHigh_data->Clone("hNCellVsEGammasNLHigh_RW_data");
  // get electron fraction
  TH2F* hElecFracHigh = (TH2F*) hNCellVsEGammasNLTrueElecHigh_MC->Clone("hElecFracHigh");
  hElecFracHigh->Divide(hNCellVsEGammasNLHigh_MC);
  hElecFracHigh->Multiply(hNCellVsEGammasNLHigh_data);
  hNCellVsEGammasNLHigh_RW_data->Add(hElecFracHigh, -1);



  hNCellVsEGammasNLHigh_RW_SBSub_data = (TH2F*) hNCellVsEGammasNLHigh_SBSub_data->Clone("hNCellVsEGammasNLHigh_RW_SBSub_data");
  // get electron fraction
  TH2F* hElecFracHigh_SBSub = (TH2F*) hNCellVsEGammasNLTrueElecHigh_MC->Clone("hElecFracHigh_SBSub");
  hElecFracHigh_SBSub->Divide(hNCellVsEGammasNLHigh_MC);
  hElecFracHigh_SBSub->Multiply(hNCellVsEGammasNLHigh_data);
  hNCellVsEGammasNLHigh_RW_SBSub_data->Add(hElecFracHigh_SBSub, -1);


  hNCellVsEGammasNLHigh_RW_MCSBSub_data = (TH2F*) hNCellVsEGammasNLHigh_MCSBSub_data->Clone("hNCellVsEGammasNLHigh_RW_MCSBSub_data");
  // get electron fraction
  TH2F* hElecFracHigh3 = (TH2F*) hNCellVsEGammasNLTrueElecHigh_MC->Clone("hElecFracHigh3");
  hElecFracHigh3->Divide(hNCellVsEGammasNLHigh_MC);
  hElecFracHigh3->Multiply(hNCellVsEGammasNLHigh_data);
  hNCellVsEGammasNLHigh_RW_MCSBSub_data->Add(hElecFracHigh3, -1);


  hNCellVsEGammasNLHigh_RW_SBSub_MC = (TH2F*) hNCellVsEGammasNLHigh_SBSub_MC->Clone("hNCellVsEGammasNLHigh_RW_SBSub_MC");
  // get electron fraction
  TH2F* hElecFracHighMC = (TH2F*) hNCellVsEGammasNLTrueElecHigh_MC->Clone("hElecFracHighMC");
  hElecFracHighMC->Divide(hNCellVsEGammasNLHigh_MC);
  hElecFracHighMC->Multiply(hNCellVsEGammasNLHigh_MC);
  hNCellVsEGammasNLHigh_RW_SBSub_MC->Add(hElecFracHighMC, -1);
  //



  // hNCellVsEGammasNLWide_RW_MCSBSub_data = (TH2F*) hNCellVsEGammasNLWide_MCSBSub_data->Clone("hNCellVsEGammasNLWide_RW_MCSBSub_data");
  // // get electron fraction
  // TH2F* hElecFracHigh = (TH2F*) hNCellVsEGammasNLTrueElecHigh_MC->Clone("hElecFracHigh");
  // hElecFracWide->Divide(hNCellVsEGammasNLWide_MC);
  // hElecFracWide->Multiply(hNCellVsEGammasNLWide_data);
  // hNCellVsEGammasNLWide_RW_MCSBSub_data->Add(hElecFracWide, -1);




  hNCellVsEGammasNL_RW_MC = (TH2F*) hNCellVsEGammasNL_MC->Clone("hNCellVsEGammasNL_RW_MC");
  // get electron fraction
  TH2F* hElecFracMC = (TH2F*) hNCellVsEGammasNLTrueElec_MC->Clone("hElecFracMC");
  hElecFracMC->Divide(hNCellVsEGammasNL_MC);
  hElecFracMC->Multiply(hNCellVsEGammasNL_MC);
  hNCellVsEGammasNL_RW_MC->Add(hElecFracMC, -1);


  hNCellVsEGammasNLLeft_RW_MC = (TH2F*) hNCellVsEGammasNL_MC->Clone("hNCellVsEGammasNLLeft_RW_MC");
  // get electron fraction
  TH2F* hElecFracLeftMC = (TH2F*) hNCellVsEGammasNLTrueElecLeft_MC->Clone("hElecFracLeftMC");
  hElecFracLeftMC->Divide(hNCellVsEGammasNLLeft_MC);
  hElecFracLeftMC->Multiply(hNCellVsEGammasNLLeft_MC);
  hNCellVsEGammasNLLeft_RW_MC->Add(hElecFracLeftMC, -1);


}

void Effi::ScaleTo(TH1D *h, TH1D *& h2, float down, float up){
  float inth2 = h2->Integral(h2->FindBin(down+ 0.0001), h2->FindBin(up - 0.0001));
  float inth = h->Integral(h->FindBin(down+ 0.0001), h->FindBin(up - 0.0001));
  h->Scale(inth2/inth);

}


// ------------------------------------------
// Plotting
// ------------------------------------------
void Effi::SetPlotting(){


  if(fPeriod.Contains("13TeV")) sEnergy = "pp #sqrt{s} = 13 TeV";
  if(fPeriod.Contains("8TeV")) sEnergy = "pp #sqrt{s} = 8 TeV";

  StyleSettingsPaper();


  hDummyEffi    = new TH2F("hDummyEffi","hDummyEffi",1000,0.5, PlotenergyHigh,1000,0., 1.4);
  SetStyleHistoTH2ForGraphs(hDummyEffi, "#it{E} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyRatio    = new TH2F("hDummyRatio","hDummyRatio",1000,0.5, PlotenergyHigh,1000,0.9, 1.4);
  SetStyleHistoTH2ForGraphs(hDummyRatio, "#it{E} (GeV)","#varepsilon_{MC}/#varepsilon_{data}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyCorr    = new TH2F("hDummyCorr","hDummyCorr",1000,0.5, PlotenergyHigh,1000,0., 1.4);
  SetStyleHistoTH2ForGraphs(hDummyCorr, "#it{E} (GeV)","#rho", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyPurity    = new TH2F("hDummyPurity","hDummyPurity",1000,0.5, PlotenergyHigh,1000,0., 1.1);
  SetStyleHistoTH2ForGraphs(hDummyPurity, "#it{E} (GeV)","P", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyMInv    = new TH2F("hDummyMInv","hDummyMInv",1000,0.01, 0.4,1000,0., 10000);
  SetStyleHistoTH2ForGraphs(hDummyMInv, "#it{M}_{inv} (GeV/c^{2})","counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyNCellRatio    = new TH2F("hDummyNCellRatio","hDummyNCellRatio",1000,0.7, PlotenergyHigh,1000,0.75, 1.45);
  SetStyleHistoTH2ForGraphs(hDummyNCellRatio, "#it{E}_{clus} (GeV/c^{2})","N_{clus}^{data}/N_{clus}^{MC}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyNCellFraction    = new TH2F("hDummyNCellFraction","hDummyNCellFraction",1000,0.4, PlotenergyHigh,1000,0., 1.);
  SetStyleHistoTH2ForGraphs(hDummyNCellFraction, "#it{E}_{clus} (GeV/c^{2})","N_{clus}^{Ncell = 1}/N_{clus}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);
  // gPad->SetLogx();
  hDummyCan = new TCanvas("Can", "", 1200, 1000);
  DrawPaperCanvasSettings(hDummyCan, 0.1, 0.003, 0.003, 0.1);


  DrawSetMarkerTGraphErr(grTB_data, 21, 1.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_MC, 25, 2.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_Ratio, 21, 1.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_Corr, 21, 1.5, kBlue + 2, kBlue + 2, 2);

  DrawSetMarker(hNCell_AllClus_Effi_data, 20, 2, kRed + 2, kRed + 2);
  DrawSetMarker(hNCell_AllClus_Effi_MC, 8, 3, kOrange + 2, kOrange + 2); //hNCell_AllClus_Effi_MC->SetFillColor(kOrange + 2);
  DrawSetMarker(hNCell_AllClus_Effi_Ratio, 20, 2, kRed + 2, kRed + 2);
  DrawSetMarker(hNCell_AllClus_Effi_Corr, 20, 2, kRed + 2, kRed + 2);

  int colGamma = kMagenta + 2;
  DrawSetMarker(hNCell_Gammas_Effi_data, 33, 2.4, colGamma, colGamma);
  DrawSetMarker(hNCell_Gammas_Effi_MC, 8, 3.4, colGamma, colGamma); //hNCell_Gammas_Effi_MC->SetFillColor(kMagenta - 6);
  DrawSetMarker(hNCell_Gammas_Effi_Ratio, 33, 2.4, colGamma, colGamma);
  DrawSetMarker(hNCell_Gammas_Effi_Corr, 33, 2.4, colGamma, colGamma);

  // true gamma
  DrawSetMarker(hNCell_TrueGammas_Effi_MC, 27, 2.4, kSpring - 8, kSpring - 8);
  // true elec
  DrawSetMarker(hNCell_TrueElec_Effi_MC, 2, 2.4, kSpring + 4, kSpring + 4);

  int colGammaLeft = kCyan + 2;
  DrawSetMarker(hNCell_GammasLeft_Effi_data, 29, 2, colGammaLeft, colGammaLeft);
  DrawSetMarker(hNCell_GammasLeft_Effi_MC, 7, 3, colGammaLeft, colGammaLeft); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasLeft_Effi_Ratio, 29, 2, colGammaLeft, colGammaLeft);
  DrawSetMarker(hNCell_GammasLeft_Effi_Corr, 29, 2, colGammaLeft, colGammaLeft);

  int colGammaWide = kGreen + 2;
  DrawSetMarker(hNCell_GammasWide_Effi_data, 34, 2, colGammaWide, colGammaWide);
  DrawSetMarker(hNCell_GammasWide_Effi_MC, 8, 3, colGammaWide, colGammaWide); //hNCell_GammasWide_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasWide_Effi_Ratio, 34, 2, colGammaWide, colGammaWide);
  DrawSetMarker(hNCell_GammasWide_Effi_Corr, 34, 2, colGammaWide, colGammaWide);

  int colGammaWideSmoothed = kMagenta - 7;
  DrawSetMarker(hNCell_GammasWide_Smoothed_Effi_data, 34, 2, colGammaWideSmoothed, colGammaWideSmoothed);
  DrawSetMarker(hNCell_GammasWide_Smoothed_Effi_MC, 8, 3, colGammaWideSmoothed, colGammaWideSmoothed); //hNCell_GammasWide_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasWide_Smoothed_Effi_Ratio, 34, 2, colGammaWideSmoothed, colGammaWideSmoothed);
  DrawSetMarker(hNCell_GammasWide_Smoothed_Effi_Corr, 34, 2, colGammaWideSmoothed, colGammaWideSmoothed);

  int colGammaHigh = kYellow  - 2;
  DrawSetMarker(hNCell_GammasHigh_Effi_data, 21, 2, colGammaHigh, colGammaHigh);
  DrawSetMarker(hNCell_GammasHigh_Effi_MC, 9, 3, colGammaHigh, colGammaHigh); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasHigh_Effi_Ratio, 21, 2, colGammaHigh, colGammaHigh);
  DrawSetMarker(hNCell_GammasHigh_Effi_Corr, 21, 2, colGammaHigh, colGammaHigh);
  //
  // int colGammaLow = kPink;
  // DrawSetMarker(hNCell_GammasLow_Effi_data, 31, 2, colGammaLow, colGammaLow);
  // DrawSetMarker(hNCell_GammasLow_Effi_MC, 2, 3, colGammaLow, colGammaLow); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  // DrawSetMarker(hNCell_GammasLow_Effi_Ratio, 31, 2, colGammaLow, colGammaLow);
  // DrawSetMarker(hNCell_GammasLow_Effi_Corr, 31, 2, colGammaLow, colGammaLow);

  int colGammaLow = kPink;
  DrawSetMarker(hNCell_GammasLow_Effi_data, 31, 2, colGammaLow, colGammaLow);
  DrawSetMarker(hNCell_GammasLow_Effi_MC, 2, 3, colGammaLow, colGammaLow); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasLow_Effi_Ratio, 31, 2, colGammaLow, colGammaLow);
  DrawSetMarker(hNCell_GammasLow_Effi_Corr, 31, 2, colGammaLow, colGammaLow);

  int colGammaRWHigh = kPink  - 2;
  DrawSetMarker(hNCell_GammasHigh_RW_Effi_data, 21, 2, colGammaRWHigh, colGammaRWHigh);
  DrawSetMarker(hNCell_GammasHigh_RW_Effi_MC, 9, 3, colGammaRWHigh + 4, colGammaRWHigh + 4); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasHigh_RW_Effi_Ratio, 21, 2, colGammaRWHigh, colGammaRWHigh);
  DrawSetMarker(hNCell_GammasHigh_RW_Effi_Corr, 21, 2, colGammaRWHigh, colGammaRWHigh);


  int colGammaRWHighSBSub = kGreen + 2;
  DrawSetMarker(hNCell_GammasHigh_RW_SBSub_Effi_data, 25, 2, colGammaRWHighSBSub, colGammaRWHighSBSub);
  DrawSetMarker(hNCell_GammasHigh_RW_SBSub_Effi_MC, 9, 3, colGammaRWHighSBSub, colGammaRWHighSBSub); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasHigh_RW_SBSub_Effi_Ratio, 25, 2, colGammaRWHighSBSub, colGammaRWHighSBSub);
  DrawSetMarker(hNCell_GammasHigh_RW_SBSub_Effi_Corr, 25, 2, colGammaRWHighSBSub, colGammaRWHighSBSub);


  int colGammaRWHighMCSBSub = kOrange + 2;
  DrawSetMarker(hNCell_GammasHigh_RW_MCSBSub_Effi_data, 21, 2, colGammaRWHighMCSBSub, colGammaRWHighMCSBSub);
  DrawSetMarker(hNCell_GammasHigh_RW_MCSBSub_Effi_MC, 9, 3, colGammaRWHighMCSBSub, colGammaRWHighMCSBSub); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasHigh_RW_MCSBSub_Effi_Ratio, 21, 2, colGammaRWHighMCSBSub, colGammaRWHighMCSBSub);
  DrawSetMarker(hNCell_GammasHigh_RW_MCSBSub_Effi_Corr, 21, 2, colGammaRWHighMCSBSub, colGammaRWHighMCSBSub);


  int colGammaRWWideSBSub = kGray + 2;
  DrawSetMarker(hNCell_GammasWide_RW_SBSub_Effi_data, 25, 2, colGammaRWWideSBSub, colGammaRWWideSBSub);
  DrawSetMarker(hNCell_GammasWide_RW_SBSub_Effi_MC, 8, 3, colGammaRWWideSBSub, colGammaRWWideSBSub); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasWide_RW_SBSub_Effi_Ratio, 25, 2, colGammaRWWideSBSub, colGammaRWWideSBSub);
  DrawSetMarker(hNCell_GammasWide_RW_SBSub_Effi_Corr, 25, 2, colGammaRWWideSBSub, colGammaRWWideSBSub);


  int colGammaWideSBSub = kBlue - 7;
  DrawSetMarker(hNCell_GammasWide_SBSub_Effi_data, 25, 2, colGammaWideSBSub, colGammaWideSBSub);
  DrawSetMarker(hNCell_GammasWide_SBSub_Effi_MC, 8, 3, colGammaWideSBSub, colGammaWideSBSub); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasWide_SBSub_Effi_Ratio, 25, 2, colGammaWideSBSub, colGammaWideSBSub);
  DrawSetMarker(hNCell_GammasWide_SBSub_Effi_Corr, 25, 2, colGammaWideSBSub, colGammaWideSBSub);



  int colGammaRWWideSBSub_Smoothed = kGreen + 3;
  DrawSetMarker(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_data, 28, 2, colGammaRWWideSBSub_Smoothed, colGammaRWWideSBSub_Smoothed);
  DrawSetMarker(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_MC, 8, 3, colGammaRWWideSBSub_Smoothed, colGammaRWWideSBSub_Smoothed); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Ratio, 28, 2, colGammaRWWideSBSub_Smoothed, colGammaRWWideSBSub_Smoothed);
  DrawSetMarker(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Corr, 28, 2, colGammaRWWideSBSub_Smoothed, colGammaRWWideSBSub_Smoothed);


  int colGammaRW = kBlack;
  DrawSetMarker(hNCell_Gammas_RW_Effi_data, 33, 2.4, colGammaRW, colGammaRW);
  DrawSetMarker(hNCell_Gammas_RW_Effi_MC, 7, 3.4, colGammaRW, colGammaRW); //hNCell_Gammas_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_Gammas_RW_Effi_Ratio, 33, 2.4, colGammaRW, colGammaRW);
  DrawSetMarker(hNCell_Gammas_RW_Effi_Corr, 33, 2.4, colGammaRW, colGammaRW);

  // int colGammaRW = kBlack;
  DrawSetMarker(hNCell_GammasSB_Effi_data, 33, 2.4, colGammaRW, colGammaRW);
  DrawSetMarker(hNCell_GammasSB_Effi_MC, 7, 3.4, colGammaRW, colGammaRW); //hNCell_Gammas_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_GammasSB_Effi_Ratio, 33, 2.4, colGammaRW, colGammaRW);
  DrawSetMarker(hNCell_GammasSB_Effi_Corr, 33, 2.4, colGammaRW, colGammaRW);

  int colGammaRWLeft = kPink + 2;
  DrawSetMarker(hNCell_GammasLeft_RW_Effi_data, 33, 2.4, colGammaRWLeft, colGammaRWLeft);
  DrawSetMarker(hNCell_GammasLeft_RW_Effi_MC, 7, 3.4, colGammaRWLeft, colGammaRWLeft); //hNCell_GammasLeft_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_GammasLeft_RW_Effi_Ratio, 33, 2.4, colGammaRWLeft, colGammaRWLeft);
  DrawSetMarker(hNCell_GammasLeft_RW_Effi_Corr, 33, 2.4, colGammaRWLeft, colGammaRWLeft);

  int colGammaRWWide = kPink + 2;
  DrawSetMarker(hNCell_GammasWide_RW_Effi_data, 33, 2.4, colGammaRWWide, colGammaRWWide);
  DrawSetMarker(hNCell_GammasWide_RW_Effi_MC, 7, 3.4, colGammaRWWide, colGammaRWWide); //hNCell_GammasWide_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_GammasWide_RW_Effi_Ratio, 33, 2.4, colGammaRWWide, colGammaRWWide);
  DrawSetMarker(hNCell_GammasWide_RW_Effi_Corr, 33, 2.4, colGammaRWWide, colGammaRWWide);

  int colConvRWWide = kBlue + 7;
  DrawSetMarker(hNCell_ConversionsWide_RW_Effi_data, 28, 2.4, colConvRWWide, colConvRWWide);
  DrawSetMarker(hNCell_ConversionsWide_RW_Effi_MC, 5, 3.4, colConvRWWide, colConvRWWide); //hNCell_GammasWide_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_ConversionsWide_RW_Effi_Ratio, 28, 2.4, colConvRWWide, colConvRWWide);
  DrawSetMarker(hNCell_ConversionsWide_RW_Effi_Corr, 28, 2.4, colConvRWWide, colConvRWWide);

  int colElectrons = kBlue + 7;
  DrawSetMarker(hNCell_Electrons_Effi_data, 28, 2.4, colElectrons, colElectrons);
  DrawSetMarker(hNCell_Electrons_Effi_MC, 5, 3.4, colConvRWWide, colElectrons); //hNCell_GammasWide_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_Electrons_Effi_Ratio, 28, 2.4, colElectrons, colElectrons);
  DrawSetMarker(hNCell_Electrons_Effi_Corr, 28, 2.4, colElectrons, colElectrons);


  DrawSetMarker(hGammaPurity, 1, 3.4, kRed + 2, kRed + 2);
  DrawSetMarker(hElecPurity, 1, 3.4, kBlack, kBlack);
  DrawSetMarker(hUnPurity, 1, 3.4, kGreen + 2, kGreen + 2);

  DrawSetMarker(hGammaPurityLeft, 8, 3.4, kRed + 2, kRed + 2);
  DrawSetMarker(hElecPurityLeft, 8, 3.4, kBlack, kBlack);
  DrawSetMarker(hUnPurityLeft, 8, 3.4, kGreen + 2, kGreen + 2);

  DrawSetMarker(hGammaPurityWide, 2, 3.4, kRed + 2, kRed + 2);
  DrawSetMarker(hElecPurityWide, 2, 3.4, kBlack, kBlack);
  DrawSetMarker(hUnPurityWide, 2, 3.4, kGreen + 2, kGreen + 2);


  DrawSetMarker(hGammaPurityHigh, 8, 3.4, kRed + 2, kRed + 2);
  DrawSetMarker(hElecPurityHigh, 8, 3.4, kBlack, kBlack);
  DrawSetMarker(hUnPurityHigh, 8, 3.4, kGreen + 2, kGreen + 2);

  DrawSetMarker(hGammaPuritySB, 1, 3.4, kRed + 2, kRed + 2);
  DrawSetMarker(hElecPuritySB, 1, 3.4, kBlack, kBlack);
  DrawSetMarker(hUnPuritySB, 1, 3.4, kGreen + 2, kGreen + 2);

  // Inv Mass stuff
  DrawSetMarker(hInvMassVsPt_Slice, 20, 1.4, kBlack, kBlack);
  DrawSetMarker(hInvMassBackVsPt_Slice, 24, 1.4, kBlack, kBlack);
  DrawSetMarker(hInvMassVsPt_data_Slice, 20, 2., kBlue, kBlue);
  DrawSetMarker(hInvMassVsPtGG_Slice, 2, 3.4, kRed + 2, kRed + 2);
  DrawSetMarker(hInvMassVsPtGC_Slice, 7, 3.4, kGreen + 2, kGreen + 2);
  DrawSetMarker(hInvMassVsPtCC_Slice, 8, 3.4, kOrange + 2, kOrange + 2);


  DrawSetMarker(hInvMassVsPtHigh_Slice_data, 20, 2.4, kBlack, kBlack);
  DrawSetMarker(hInvMassVsPtHigh1cell_Slice_data, 8, 3.4, kOrange + 2, kOrange + 2);
  DrawSetMarker(hInvMassVsPtHigh2cell_Slice_data, 9, 3.4, kGreen + 2, kGreen + 2);
  DrawSetMarker(hInvMassVsPtHigh3cell_Slice_data, 8, 3.4, kBlue + 2, kBlue + 2);
  DrawSetMarker(hInvMassVsPtHigh_Slice_MC, 24, 2.4, kRed, kRed);
  DrawSetMarker(hInvMassVsPtHigh1cell_Slice_MC, 8, 3.4, kOrange + 2, kOrange + 2);
  DrawSetMarker(hInvMassVsPtHigh2cell_Slice_MC, 9, 3.4, kGreen + 2, kGreen + 2);

  // MC reweighted for consistency check
  DrawSetMarker(hNCell_Gammas_RW_Effi_MC_consistency, 33, 2.4, kGreen , kGreen);
}

void Effi::PlotEffi_AllClusAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_AllClus_Effi_data->Draw("same,p");
  hNCell_AllClus_Effi_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.65, 0.16, 0.95, 0.34, 40);
  leg->AddEntry(hNCell_AllClus_Effi_data, "P2, data, all clus.", "p");
  leg->AddEntry(hNCell_AllClus_Effi_MC, "P2, MC, all clus.", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd("0.12 < M_{#gamma#gamma} < 0.14 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_AllClusAndTB.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}


void Effi::PlotEffi_AllGammaAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_Gammas_Effi_data->Draw("same,p");
  hNCell_Gammas_Effi_MC->Draw("same,histc");
  hNCell_TrueGammas_Effi_MC->Draw("same,histc");
  hNCell_TrueElec_Effi_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.65, 0.16, 0.95, 0.34, 40);
  leg->AddEntry(hNCell_Gammas_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_Gammas_Effi_MC, "P2, MC, #gamma", "l");
  leg->AddEntry(hNCell_TrueGammas_Effi_MC, "P2, MC, true #gamma", "l");
  leg->AddEntry(hNCell_TrueElec_Effi_MC, "P2, MC, true e^{#pm}", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotEffi_RWGammaAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  // hNCell_Gammas_RW_Effi_MC_consistency->Draw("same,p");
  hNCell_Gammas_RW_Effi_data->Draw("same,p");
  hNCell_Gammas_Effi_data->Draw("same,p");
  // hNCell_Gammas_Effi_MC->Draw("same,E3");
  hNCell_TrueGammas_Effi_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.65, 0.16, 0.95, 0.34, 40);
  leg->AddEntry(hNCell_Gammas_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_data, "P2, data, #gamma RW", "p");
  // leg->AddEntry(hNCell_Gammas_RW_Effi_MC, "P2, MC, #gamma", "l");
  leg->AddEntry(hNCell_TrueGammas_Effi_MC, "P2, MC, true #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_Reweighted.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotEffi_SBGammaAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");

  // hNCell_Gammas_RW_Effi_MC_consistency->Draw("same,p");
  hNCell_GammasSB_Effi_MC->Draw("same,histc");
  hNCell_GammasSB_Effi_data->Draw("same,p");
  // hNCell_Gammas_Effi_MC->Draw("same,E3");

  TLegend *leg = GetAndSetLegend2(0.65, 0.16, 0.95, 0.34, 40);
  leg->AddEntry(hNCell_GammasSB_Effi_data, "P2, data, SB #gamma", "p");
  leg->AddEntry(hNCell_GammasSB_Effi_MC, "P2, MC, SB #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" only higher #gamma taken!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_Sideband.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotEffi_RWGammaHighAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  // hNCell_Gammas_RW_Effi_MC_consistency->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_MC->Draw("same,histc");
  hNCell_GammasHigh_Effi_data->Draw("same,p");
  hNCell_GammasHigh_Effi_MC->Draw("same,histc");
  hNCell_TrueGammas_Effi_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.65, 0.16, 0.95, 0.4, 40);
  leg->AddEntry(hNCell_GammasHigh_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_MC, "P2, MC, #gamma RW", "l");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_data, "P2, data, #gamma RW", "p");
  leg->AddEntry(hNCell_GammasHigh_Effi_MC, "P2, MC, #gamma", "l");
  leg->AddEntry(hNCell_TrueGammas_Effi_MC, "P2, MC, true #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" only higher #gamma taken!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaHighAndTB_Reweighted.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotEffi_RWGammaHighAndTB_MCCheck(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  // hNCell_Gammas_RW_Effi_MC_consistency->Draw("same,p");
  // hNCell_GammasHigh_RW_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_MC->Draw("same,p");
  // hNCell_GammasHigh_Effi_data->Draw("same,p");
  hNCell_GammasHigh_Effi_MC->Draw("same,p");
  hNCell_TrueGammas_Effi_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.65, 0.16, 0.95, 0.4, 30);
  // leg->AddEntry(hNCell_GammasHigh_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_MC, "P2, MC, #gamma RW", "p");
  // leg->AddEntry(hNCell_GammasHigh_RW_Effi_data, "P2, data, #gamma RW", "p");
  leg->AddEntry(hNCell_GammasHigh_Effi_MC, "P2, MC, #gamma", "p");
  leg->AddEntry(hNCell_TrueGammas_Effi_MC, "P2, MC, true #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" only higher #gamma taken!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaHighAndTB_Reweighted_MCConsistency.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotEffi_TrueAndReweighted(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_Gammas_RW_Effi_MC_consistency->Draw("same,p");
  hNCell_Gammas_RW_Effi_data->Draw("same,p");
  hNCell_Gammas_Effi_data->Draw("same,p");
  // hNCell_Gammas_Effi_MC->Draw("same,E3");
  hNCell_TrueGammas_Effi_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.65, 0.16, 0.95, 0.34, 40);
  leg->AddEntry(hNCell_Gammas_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_data, "P2, data, #gamma RW", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_MC_consistency, "P2, MC, #gamma RW", "p");
  leg->AddEntry(hNCell_TrueGammas_Effi_MC, "P2, MC, true #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}}< M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_Reweighted_consistency.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotEffi_HighAllMethods(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_GammasHigh_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_SBSub_Effi_data->Draw("same,p");
  hNCell_GammasHigh_Effi_MC->Draw("same,histc");
  hNCell_GammasHigh_RW_Effi_MC->Draw("same,histc");
  hNCell_GammasHigh_RW_SBSub_Effi_MC->Draw("same,histc");

  // hNCell_Gammas_Effi_MC->Draw("same,E3");

  TLegend *leg = GetAndSetLegend2(0.55, 0.16, 0.95, 0.44, 40);
  leg->AddEntry(hNCell_GammasHigh_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_data, "P2, data, #gamma RW", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_data, "P2, data, #gamma RW + SB sub", "p");
  leg->AddEntry(hNCell_GammasHigh_Effi_MC, "P2, MC, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_MC, "P2, MC, #gamma RW", "l");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_MC, "P2, MC, #gamma RW SB sub", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}}< M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" only high #gamma selected!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_SBSub.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotEffi_HighAllMethodsMCSBSub(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_GammasHigh_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_MCSBSub_Effi_data->Draw("same,p");
  // hNCell_GammasWide_RW_SBSub_Effi_data->Draw("same,p");
  hNCell_GammasHigh_Effi_MC->Draw("same,histc");
  hNCell_GammasHigh_RW_SBSub_Effi_MC->Draw("same,histc");

  // hNCell_Gammas_Effi_MC->Draw("same,E3");

  TLegend *leg = GetAndSetLegend2(0.55, 0.16, 0.95, 0.44, 40);
  leg->AddEntry(hNCell_GammasHigh_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_data, "P2, data, #gamma RW", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_MCSBSub_Effi_data, "P2, data, #gamma RW + MC SB sub", "p");
  leg->AddEntry(hNCell_GammasHigh_Effi_MC, "P2, MC, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_MC, "P2, MC, #gamma RW", "l");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_MC, "P2, MC, #gamma RW SB sub", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" only high #gamma selected!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_MCSidebandSub.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotEffi_HighAndWideRWSBSub(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_GammasHigh_Effi_data->Draw("same,p");
  hNCell_GammasWide_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_data->Draw("same,p");
  hNCell_GammasWide_RW_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_SBSub_Effi_data->Draw("same,p");
  hNCell_GammasWide_RW_SBSub_Effi_data->Draw("same,p");

  hNCell_GammasHigh_RW_SBSub_Effi_MC->Draw("same,histc");
  hNCell_GammasWide_RW_SBSub_Effi_MC->Draw("same,histc");

  // hNCell_Gammas_Effi_MC->Draw("same,E3");

  TLegend *leg = GetAndSetLegend2(0.55, 0.16, 0.95, 0.44, 40);
  leg->AddEntry(hNCell_GammasHigh_Effi_data, "P2, data, high #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_data, "P2, data, high #gamma RW", "p");
  leg->AddEntry(hNCell_GammasWide_RW_Effi_data, "P2, data, #gamma RW", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_data, "P2, data, high #gamma RW SB sub", "p");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_data, "P2, data, #gamma RW SB sub", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_MC, "P2, MC, high #gamma RW SB sub", "l");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_MC, "P2, MC, #gamma RW SB sub", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_HighAndWideComp.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_WideSBSubGammaAndConv(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");

  hNCell_GammasWide_Effi_data->Draw("same,p");
  hNCell_GammasWide_RW_Effi_data->Draw("same,p");
  hNCell_ConversionsWide_RW_Effi_data->Draw("same,p");

  hNCell_GammasWide_Effi_MC->Draw("same,histc");
  hNCell_GammasWide_RW_Effi_MC->Draw("same,histc");
  hNCell_ConversionsWide_RW_Effi_MC->Draw("same,histc");
  // hNCell_TrueGammas_Effi_MC->Draw("same,p");
  // hNCell_TrueElec_Effi_MC->Draw("same,p");

  cout<<"hNCell_GammasWide_RW_Effi_data->GetBinContent(7) = "<<hNCell_GammasWide_RW_Effi_data->GetBinContent(7)<<endl;
  cout<<"hNCell_ConversionsWide_RW_Effi_data->GetBinContent(7) = "<<hNCell_ConversionsWide_RW_Effi_data->GetBinContent(7)<<endl;
  cout<<"hNCell_GammasWide_RW_Effi_MC->GetBinContent(7) = "<<hNCell_GammasWide_RW_Effi_MC->GetBinContent(7)<<endl;
  cout<<"hNCell_ConversionsWide_RW_Effi_MC->GetBinContent(7) = "<<hNCell_ConversionsWide_RW_Effi_MC->GetBinContent(7)<<endl;

  // hNCell_Gammas_Effi_MC->Draw("same,E3");

  TLegend *leg = GetAndSetLegend2(0.55, 0.16, 0.95, 0.44, 40);
  leg->AddEntry(hNCell_GammasWide_Effi_data, "P2, data, #gamma + #gamma_{conv}, wide, SB Sub", "p");
  leg->AddEntry(hNCell_GammasWide_RW_Effi_data, "P2, data, #gamma, wide, SB Sub", "p");
  leg->AddEntry(hNCell_GammasWide_RW_Effi_MC, "P2, MC, #gamma, wide, SB Sub", "l");
  leg->AddEntry(hNCell_GammasWide_Effi_MC, "P2, MC, #gamma + #gamma_{conv}, wide, SB Sub", "l");
  leg->AddEntry(hNCell_ConversionsWide_RW_Effi_data, "P2, data, #gamma_{conv}, wide, SB Sub", "p");
  leg->AddEntry(hNCell_ConversionsWide_RW_Effi_MC, "P2, MC, #gamma_{conv}, wide, SB Sub", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndConv_WideComp.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_HighAndWideRWSBSub_Smoothed(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");

  hNCell_GammasWide_RW_SBSub_Effi_data->Draw("same,p");
  hNCell_GammasWide_RW_SBSub_Smoothed_Effi_data->Draw("same,p");

  hNCell_GammasWide_RW_SBSub_Effi_MC->Draw("same,histc");
  hNCell_GammasWide_RW_SBSub_Smoothed_Effi_MC->Draw("same,histc");

  // hNCell_Gammas_Effi_MC->Draw("same,E3");

  TLegend *leg = GetAndSetLegend2(0.55, 0.16, 0.95, 0.44, 40);

  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_data, "P2, data, #gamma RW SB sub", "p");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_MC, "P2, MC, #gamma RW SB sub", "l");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_data, "P2, data, #gamma RW SB sub, smoothed", "p");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_MC, "P2, MC, #gamma RW SB sub, smoothed", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_SmoothedComp.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_WideOnlySBSub(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");

  hNCell_GammasWide_SBSub_Effi_data->Draw("same,p");
  hNCell_GammasWide_RW_SBSub_Effi_data->Draw("same,p");

  hNCell_GammasWide_SBSub_Effi_MC->Draw("same,histc");
  hNCell_GammasWide_RW_SBSub_Effi_MC->Draw("same,histc");

  // hNCell_Gammas_Effi_MC->Draw("same,E3");

  TLegend *leg = GetAndSetLegend2(0.55, 0.16, 0.95, 0.44, 40);

  leg->AddEntry(hNCell_GammasWide_SBSub_Effi_data, "P2, data, #gamma SB sub", "p");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_data, "P2, data, #gamma RW SB sub", "p");
  leg->AddEntry(hNCell_GammasWide_SBSub_Effi_MC, "P2, MC, #gamma SB sub", "l");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_MC, "P2, MC, #gamma RW SB sub", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_OnlySBSub.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_HighAndWideRWSBSub_Light(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  // hNCell_GammasHigh_Effi_data->Draw("same,p");
  // hNCell_GammasWide_Effi_data->Draw("same,p");
  // hNCell_GammasHigh_RW_Effi_data->Draw("same,p");
  // hNCell_GammasWide_RW_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_SBSub_Effi_data->Draw("same,p");
  hNCell_GammasWide_RW_SBSub_Effi_data->Draw("same,p");

  hNCell_GammasHigh_RW_SBSub_Effi_MC->Draw("same,histc");
  hNCell_GammasWide_RW_SBSub_Effi_MC->Draw("same,histc");

  // hNCell_Gammas_Effi_MC->Draw("same,E3");

  TLegend *leg = GetAndSetLegend2(0.45, 0.16, 0.95, 0.36, 40);
  // leg->AddEntry(hNCell_GammasHigh_Effi_data, "P2, data, high #gamma", "p");
  // leg->AddEntry(hNCell_GammasWide_Effi_data, "P2, data, #gamma", "p");
  // leg->AddEntry(hNCell_GammasHigh_RW_Effi_data, "P2, data, high #gamma RW", "p");
  // leg->AddEntry(hNCell_GammasWide_RW_Effi_data, "P2, data, #gamma RW", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_data, "P2, data, high #gamma RW SB sub", "p");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_data, "P2, data, #gamma RW SB sub", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_MC, "P2, MC, high #gamma RW SB sub", "l");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_MC, "P2, MC, #gamma RW SB sub", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_HighAndWideComp_Light.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotEffi_HighAndWide_Light(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_GammasHigh_Effi_data->Draw("same,p");
  hNCell_GammasLow_Effi_data->Draw("same,p");
  hNCell_GammasWide_Effi_data->Draw("same,p");
  hNCell_Gammas_Effi_data->Draw("same,p");
  hNCell_GammasLeft_Effi_data->Draw("same,p");
  hNCell_GammasHigh_Effi_MC->Draw("same,histc");
  hNCell_GammasLow_Effi_MC->Draw("same,histc");
  hNCell_GammasWide_Effi_MC->Draw("same,histc");
  hNCell_Gammas_Effi_MC->Draw("same,histc");
  hNCell_GammasLeft_Effi_MC->Draw("same,histc");


  TLegend *leg = GetAndSetLegend2(0.45, 0.16, 0.95, 0.46, 40);
  leg->AddEntry(hNCell_GammasWide_Effi_data, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_Effi_MC, "P2, MC, #gamma", "l");
  leg->AddEntry(hNCell_GammasHigh_Effi_data, "P2, data, high #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_Effi_MC, "P2, MC, high #gamma", "l");
  leg->AddEntry(hNCell_GammasLow_Effi_data, "P2, data, low #gamma", "p");
  leg->AddEntry(hNCell_GammasLow_Effi_MC, "P2, MC, low #gamma", "l");
  leg->AddEntry(hNCell_Gammas_Effi_data, "P2, data, right #gamma", "p");
  leg->AddEntry(hNCell_Gammas_Effi_MC, "P2, MC, right #gamma", "l");
  leg->AddEntry(hNCell_GammasLeft_Effi_data, "P2, data, left #gamma", "p");
  leg->AddEntry(hNCell_GammasLeft_Effi_MC, "P2, MC, left #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_Comp_Light.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_WideSmoothed(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_GammasWide_Effi_data->Draw("same,p");
  hNCell_GammasWide_Effi_MC->Draw("same,histc");
  hNCell_GammasWide_Smoothed_Effi_data->Draw("same,p");
  hNCell_GammasWide_Smoothed_Effi_MC->Draw("same,histc");


  TLegend *leg = GetAndSetLegend2(0.45, 0.16, 0.95, 0.46, 40);
  leg->AddEntry(hNCell_GammasWide_Effi_data, "P2, data", "p");
  leg->AddEntry(hNCell_GammasWide_Effi_MC, "P2, MC", "l");
  leg->AddEntry(hNCell_GammasWide_Smoothed_Effi_data, "P2, data, smoothed", "p");
  leg->AddEntry(hNCell_GammasWide_Smoothed_Effi_MC, "P2, MC, smoothed", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_GammaAndTB_Smoothed.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_Electrons(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_Electrons_Effi_data->Draw("same,p");
  hNCell_Electrons_Effi_MC->Draw("same,histc");


  TLegend *leg = GetAndSetLegend2(0.45, 0.16, 0.95, 0.46, 40);
  leg->AddEntry(hNCell_Electrons_Effi_data, "P2, data, electrons", "p");
  leg->AddEntry(hNCell_Electrons_Effi_MC, "P2, MC, electrons", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("electron clusters selected via. track matching",0.15,0.92,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Effi_Electrons.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotRatio_GammasAndRW(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_RW_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.65, 0.66, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_Gammas_Effi_Ratio, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_Ratio, "P2, data, RW #gamma", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}}< M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Ratio_GammaAndTB_Reweighted.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotRatio_SBGammasAndRW(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_Effi_Ratio->Draw("same,p");
  hNCell_GammasSB_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.65, 0.66, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_Gammas_Effi_Ratio, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasSB_Effi_Ratio, "P2, data, SB #gamma", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}}< M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Ratio_SBGammaAndTB.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotRatio_GammasAndRWLeft(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_GammasLeft_Effi_Ratio->Draw("same,p");
  hNCell_GammasLeft_RW_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.65, 0.66, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasLeft_Effi_Ratio, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasLeft_RW_Effi_Ratio, "P2, data, RW #gamma", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Ratio_GammaAndTB_Reweighted_Left.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

//
// void Effi::PlotRatio_GammasAndRWWide(){
//   hDummyCan->cd();
//   hDummyRatio->Draw();
//   DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
//   grTB_Ratio->Draw("same,p");
//
//   hNCell_AllClus_Effi_Ratio->Draw("same,p");
//   hNCell_GammasWide_Effi_Ratio->Draw("same,p");
//   hNCell_GammasWide_RW_Effi_Ratio->Draw("same,p");
//
//   TLegend *leg = GetAndSetLegend2(0.65, 0.66, 0.95, 0.84, 40);
//   leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
//   leg->AddEntry(hNCell_GammasWide_Effi_Ratio, "P2, data, #gamma", "p");
//   leg->AddEntry(hNCell_GammasWide_RW_Effi_Ratio, "P2, data, RW #gamma", "p");
//   leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
//   leg->Draw("same");
//
//   drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
//   drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
//   drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
//   // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
//   hDummyCan->SaveAs(Form("%s/%s/Ratio_GammaAndTB_Reweighted_Wide.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
// }

void Effi::PlotRatio_GammasAndRWHigh(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_GammasHigh_Effi_Ratio->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_Ratio->Draw("same,p");
  hNCell_GammasHigh_RW_SBSub_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.65, 0.66, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasHigh_Effi_Ratio, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_Ratio, "P2, data, RW #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_Ratio, "P2, data, RW, SB sub #gamma", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" only higher #gamma",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Ratio_GammaAndTB_Reweighted_High.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotRatio_GammasAndRWWide(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_GammasWide_Effi_Ratio->Draw("same,p");
  hNCell_GammasWide_Smoothed_Effi_Ratio->Draw("same,p");
  hNCell_GammasWide_RW_Effi_Ratio->Draw("same,p");
  hNCell_GammasWide_RW_SBSub_Effi_Ratio->Draw("same,p");
  hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.5, 0.6, 0.95, 0.82, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasWide_Effi_Ratio, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_Smoothed_Effi_Ratio, "P2, data, #gamma, smoothed", "p");
  leg->AddEntry(hNCell_GammasWide_RW_Effi_Ratio, "P2, data, RW #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_Ratio, "P2, data, RW, SB sub #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Ratio, "P2, data, RW, SB sub #gamma, smoothed", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} < M_{#gamma#gamma} - 0.05 < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" only higher #gamma",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Ratio_GammaAndTB_Reweighted_SBSub_Wide.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotRatio_GammasAndConvRWWide(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_GammasWide_RW_Effi_Ratio->Draw("same,p");
  hNCell_ConversionsWide_RW_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.65, 0.66, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_GammasWide_RW_Effi_Ratio, "P2, data, RW #gamma", "p");
  leg->AddEntry(hNCell_ConversionsWide_RW_Effi_Ratio, "P2, data, RW #gamma_{conv}", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} < M_{#gamma#gamma} - 0.05 < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" only higher #gamma",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Ratio_GammaAndConv_Reweighted_Wide.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotRatio_GammasAndRWHighMCSBSub(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_GammasHigh_Effi_Ratio->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_Ratio->Draw("same,p");
  hNCell_GammasHigh_RW_MCSBSub_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.55, 0.55, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasHigh_Effi_Ratio, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_Ratio, "P2, data, RW #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_MCSBSub_Effi_Ratio, "P2, data, RW, MC SB sub #gamma", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" only higher #gamma",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Ratio_GammaAndTB_Reweighted_High_MCSBSub.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotCorr_GammasAndRW(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_AllClus_Effi_Corr->Draw("same,p");
  hNCell_Gammas_Effi_Corr->Draw("same,p");
  hNCell_Gammas_RW_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Corr, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_Gammas_Effi_Corr, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_Corr, "P2, data, RW #gamma", "p");
  leg->AddEntry(grTB_Corr, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Corr_GammaAndTB_Reweighted.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotCorr_GammasAndRW_Left(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_AllClus_Effi_Corr->Draw("same,p");
  hNCell_GammasLeft_Effi_Corr->Draw("same,p");
  hNCell_GammasLeft_RW_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Corr, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasLeft_Effi_Corr, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasLeft_RW_Effi_Corr, "P2, data, RW #gamma", "p");
  leg->AddEntry(grTB_Corr, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Corr_GammaLeftAndTB_Reweighted.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotCorr_GammasAndRW_Wide(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_AllClus_Effi_Corr->Draw("same,p");
  hNCell_GammasWide_Effi_Corr->Draw("same,p");
  hNCell_GammasWide_Smoothed_Effi_Corr->Draw("same,p");
  hNCell_GammasWide_RW_Effi_Corr->Draw("same,p");
  hNCell_GammasWide_RW_SBSub_Effi_Corr->Draw("same,p");
  // hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Corr, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasWide_Effi_Corr, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_Smoothed_Effi_Corr, "P2, data, #gamma, smoothed", "p");
  leg->AddEntry(hNCell_GammasWide_RW_Effi_Corr, "P2, data, RW #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_RW_SBSub_Effi_Corr, "P2, data, RW SB sub. #gamma", "p");
  // leg->AddEntry(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Corr, "P2, data, RW SB sub. #gamma, smoothed", "p");
  leg->AddEntry(grTB_Corr, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Corr_GammaWideAndTB_Reweighted_SBSub.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}



void Effi::PlotCorr_GammasAndConv_Wide(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_GammasWide_RW_Effi_Corr->Draw("same,p");
  hNCell_ConversionsWide_RW_Effi_Corr->Draw("same,p");
  hNCell_AllClus_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
  leg->AddEntry(hNCell_GammasWide_RW_Effi_Corr, "P2, data, RW #gamma", "p");
  leg->AddEntry(hNCell_ConversionsWide_RW_Effi_Corr, "P2, data, RW #gamma_{conv}", "p");
  // leg->AddEntry(hNCell_GammasWide_RW_SBSub_Smoothed_Effi_Corr, "P2, data, RW SB sub. #gamma, smoothed", "p");
  leg->AddEntry(grTB_Corr, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Corr_GammaAndConv_Wide_Reweighted.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotCorr_GammasAndRW_High(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_AllClus_Effi_Corr->Draw("same,p");
  hNCell_GammasHigh_Effi_Corr->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Corr, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasHigh_Effi_Corr, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_Corr, "P2, data, RW #gamma", "p");
  leg->AddEntry(grTB_Corr, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" Only higher #gamma selected",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Corr_GammaHighAndTB_Reweighted.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotCorr_HighAllMethods(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_GammasHigh_Effi_Corr->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_Corr->Draw("same,p");
  hNCell_GammasHigh_RW_SBSub_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
  leg->AddEntry(hNCell_GammasHigh_Effi_Corr, "P2, ", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_Effi_Corr, "P2, RW", "p");
  leg->AddEntry(hNCell_GammasHigh_RW_SBSub_Effi_Corr, "P2, RW SB sub", "p");
  leg->AddEntry(grTB_Corr, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" Only higher #gamma selected",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Corr_GammaHighAndTB_AllMethods.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}

void Effi::PlotCorr_Electrons(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_Electrons_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
  leg->AddEntry(hNCell_Electrons_Effi_Corr, "P2, e^{#pm}", "p");
  leg->AddEntry(grTB_Corr, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/%s/Corr_Electrons.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotPurity(){

  hDummyCan->cd();
  hDummyPurity->Draw();



  hGammaPurity->Draw("same,hist,c");
  hGammaPurityLeft->Draw("same,hist,c");
  hGammaPurityWide->Draw("same,hist,c");
  hElecPurity->Draw("same,hist,c");
  hElecPurityLeft->Draw("same,hist,c");
  hElecPurityWide->Draw("same,hist,c");
  hUnPurityLeft->Draw("same,hist,c");
  hUnPurity->Draw("same,hist,c");
  hUnPurityWide->Draw("same,hist,c");

  TLegend *leg = GetAndSetLegend2(0.15, 0.75, 0.35, 0.95, 40);
  leg->AddEntry(hGammaPurity, "#gamma, right side", "l");
  leg->AddEntry(hElecPurity, "e^{#pm}, right side", "l");
  leg->AddEntry(hGammaPurityLeft, "#gamma, left side", "l");
  leg->AddEntry(hElecPurityLeft, "e^{#pm}, left side", "l");
  leg->AddEntry(hGammaPurityWide, "#gamma, both sides", "l");
  leg->AddEntry(hElecPurityWide, "e^{#pm}, both sides", "l");
  leg->Draw("same");
  TLegend *leg2 = GetAndSetLegend2(0.45, 0.86, 0.65, 0.95, 40);
  leg2->AddEntry(hUnPurity, "contamination", "l");
  leg2->AddEntry(hUnPurityLeft, "contamination, left side", "l");
  leg2->AddEntry(hUnPurityWide, "contamination, both sides", "l");
  leg2->Draw("same");

  hDummyCan->SaveAs(Form("%s/%s/Purity.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}


void Effi::PlotPurity2(){

  hDummyCan->cd();
  hDummyPurity->Draw();



  hGammaPuritySB->Draw("same,hist,c");
  hGammaPurityHigh->Draw("same,hist,c");
  hGammaPurityWide->Draw("same,hist,c");

  hElecPuritySB->Draw("same,hist,c");
  hElecPurityHigh->Draw("same,hist,c");
  hElecPurityWide->Draw("same,hist,c");

  hUnPurityHigh->Draw("same,hist,c");
  hUnPuritySB->Draw("same,hist,c");
  hUnPurityWide->Draw("same,hist,c");

  TLegend *leg = GetAndSetLegend2(0.15, 0.75, 0.35, 0.95, 40);
  leg->AddEntry(hGammaPurityHigh, "#gamma, high wide", "l");
  leg->AddEntry(hElecPurityHigh, "e^{#pm}, high wide", "l");
  leg->AddEntry(hGammaPuritySB, "#gamma, SB", "l");
  leg->AddEntry(hElecPuritySB, "e^{#pm}, SB", "l");
  leg->AddEntry(hGammaPurityWide, "#gamma, both sides", "l");
  leg->AddEntry(hElecPurityWide, "e^{#pm}, both sides", "l");
  leg->Draw("same");
  TLegend *leg2 = GetAndSetLegend2(0.45, 0.86, 0.65, 0.95, 40);
  leg2->AddEntry(hUnPurityHigh, "contamination, high", "l");
  leg2->AddEntry(hUnPuritySB, "contamination, SB", "l");
  leg2->AddEntry(hUnPurityWide, "contamination, both sides", "l");
  leg2->Draw("same");

  hDummyCan->SaveAs(Form("%s/%s/Purity2.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}



void Effi::PlotExampleBin(){

  hDummyCan->cd();
  hDummyMInv->GetYaxis()->SetRangeUser(0., hInvMassVsPt_Slice->GetMaximum()*1.3);
  hDummyMInv->Draw();

  float intBack = hInvMassBackVsPt_Slice->Integral(hInvMassBackVsPt_Slice->FindBin(0.2), hInvMassBackVsPt_Slice->FindBin(0.3));
  float intSignal = hInvMassVsPt_Slice->Integral(hInvMassVsPt_Slice->FindBin(0.2), hInvMassVsPt_Slice->FindBin(0.3));
  hInvMassBackVsPt_Slice->Scale(intSignal/intBack);

  hInvMassVsPt_Slice->Draw("same,p");
  hInvMassBackVsPt_Slice->Draw("same,p");
  hInvMassVsPtGG_Slice->Draw("same,histc");
  hInvMassVsPtGC_Slice->Draw("same,histc");
  hInvMassVsPtCC_Slice->Draw("same,histc");

  TH1F* hMCTrueSum = (TH1F*) hInvMassVsPtGG_Slice->Clone("hMCTrueSum");
  hMCTrueSum->Add(hInvMassVsPtGC_Slice);
  hMCTrueSum->Add(hInvMassVsPtCC_Slice);
  DrawSetMarker(hMCTrueSum, 2, 3.4, kBlue + 2, kBlue + 2);

  hMCTrueSum->Draw("same,histc");

  TFile fmasspos("fMassPos.root");
  TString name = "MassPos8TeV_func";
  cout<<"fPeriod: "<<fPeriod<<endl;
  if(fPeriod.Contains("13TeVLowB")) name = "MassPos13TeVLowB_func";
  else if(fPeriod.Contains("13TeVNomB")) name = "MassPos13TeVNomB_func";
  else if(fPeriod.Contains("8TeV")) name = "MassPos8TeV_func";

  cout<<"name: "<<name<<endl;
  TF1* func = (TF1*) fmasspos.Get(name);


  TF1 *fgaus = new TF1("fgaus", "gaus(0)", 0.1, 0.15);
  fgaus->SetParameters(100, 0.13, 0.05);
  hInvMassVsPt_Slice->Fit(fgaus, "MQR0");

  float mean = func->Eval(2.4);

  cout<<"peak pos from global fit: "<<mean<<"   from gaus: "<<fgaus->GetParameter(1)<<endl;

  TLegend *leg = GetAndSetLegend2(0.55, 0.74, 0.95, 0.95, 40);
  leg->AddEntry(hInvMassVsPt_Slice, "MC, all cluster pairs", "p");
  leg->AddEntry(hInvMassBackVsPt_Slice, "MC, scaled bck.", "p");
  leg->AddEntry(hInvMassVsPtGG_Slice, "true #gamma #gamma", "l");
  leg->AddEntry(hInvMassVsPtGC_Slice, "true #gamma e^{#pm}", "l");
  leg->AddEntry(hInvMassVsPtCC_Slice, "true e^{#pm} e^{#pm}", "l");
  leg->AddEntry(hMCTrueSum, "true Signal", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("2 < p_{T} < 3 GeV/c",0.15,0.87,textSizeLabelsRel,kFALSE);

  DrawLines(mean, mean, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kRed + 1, 2);
  DrawLines(mean - 0.05, mean - 0.05, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);
  DrawLines(mean + 0.02, mean + 0.02, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);


  hDummyCan->SaveAs(Form("%s/%s/ExampleBin.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));



/////////////////////////////////////////////

  // hInvMassVsPt_data_Slice->Scale( hInvMassVsPt_Slice->Integral(0.05, 0.2) / hInvMassVsPt_data_Slice->Integral(0.05, 0.2));
  hDummyMInv->Draw();

  hInvMassVsPt_Slice->Draw("same");
  hInvMassVsPt_data_Slice->Draw("same,p");

  TLegend *leg2 = GetAndSetLegend2(0.7, 0.8, 0.95, 0.95, 40);
  leg2->AddEntry(hInvMassVsPt_Slice, "MC, all cluster pairs", "l");
  leg2->AddEntry(hInvMassVsPt_data_Slice, "data, all cluster pairs", "l");
  leg2->Draw("same");

  drawLatexAdd(sEnergy,0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("2 < p_{T} < 3 GeV/c",0.15,0.87,textSizeLabelsRel,kFALSE);

  DrawLines(mean, mean, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kRed + 1, 2);
  DrawLines(mean - 0.05, mean - 0.05, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);
  DrawLines(mean + 0.02, mean + 0.02, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);

  hDummyCan->SaveAs(Form("%s/%s/ExampleBinWithData.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}


void Effi::PlotExampleBinCells(){

  hDummyCan->cd();
  hDummyMInv->GetYaxis()->SetRangeUser(0., hInvMassVsPtHigh_Slice_data->GetMaximum()*1.3);
  hDummyMInv->DrawCopy();

  ScaleTo(hInvMassVsPtHigh_Slice_data, hInvMassVsPtHigh_Slice_MC);
  hInvMassVsPtHigh_Slice_MC->Draw("same,p");
  hInvMassVsPtHigh_Slice_data->Draw("same");


  TFile fmasspos("fMassPos.root");
  TString name = "MassPos8TeV_func";
  cout<<"fPeriod: "<<fPeriod<<endl;
  if(fPeriod.Contains("13TeVLowB")) name = "MassPos13TeVLowB_func";
  else if(fPeriod.Contains("13TeVNomB")) name = "MassPos13TeVNomB_func";
  else if(fPeriod.Contains("8TeV")) name = "MassPos8TeV_func";
  TF1* func = (TF1*) fmasspos.Get(name);

  float mean = func->Eval(2.4);


  TLegend *leg = GetAndSetLegend2(0.7, 0.8, 0.95, 0.95, 40);
  leg->AddEntry(hInvMassVsPtHigh_Slice_MC, "MC, only high clus", "l");
  leg->AddEntry(hInvMassVsPtHigh_Slice_data, "data, only high clus", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("2 < p_{T, #gamma} < 3 GeV/c",0.15,0.87,textSizeLabelsRel,kFALSE);

  DrawLines(mean, mean, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kRed + 1, 2);
  DrawLines(mean - 0.05, mean - 0.05, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);
  DrawLines(mean + 0.02, mean + 0.02, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);


  hDummyCan->SaveAs(Form("%s/%s/ExampleBinHigh.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));



  //-----------------------------------------

  hDummyMInv->GetYaxis()->SetRangeUser(0., hInvMassVsPtHigh_Slice_data->GetMaximum()*1.3);
  hDummyMInv->DrawCopy();

  // ScaleTo(hInvMassVsPtHigh_Slice_data, hInvMassVsPtHigh_Slice_MC);
  hInvMassVsPtHigh_Slice_data->Draw("same");
  hInvMassVsPtHigh1cell_Slice_data->Draw("same");
  hInvMassVsPtHigh2cell_Slice_data->Draw("same");




  TLegend *leg2 = GetAndSetLegend2(0.7, 0.8, 0.95, 0.95, 40);
  leg2->AddEntry(hInvMassVsPtHigh_Slice_data, "data, only high clus", "l");
  leg2->AddEntry(hInvMassVsPtHigh1cell_Slice_data, "data, 1 cells", "l");
  leg2->AddEntry(hInvMassVsPtHigh2cell_Slice_data, "data, 2 cells", "l");

  leg2->Draw("same");

  drawLatexAdd(sEnergy,0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("2 < p_{T, #gamma} < 3 GeV/c",0.15,0.87,textSizeLabelsRel,kFALSE);

  DrawLines(mean, mean, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kRed + 1, 2);
  DrawLines(mean - 0.05, mean - 0.05, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);
  DrawLines(mean + 0.02, mean + 0.02, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);


  hDummyCan->SaveAs(Form("%s/%s/ExampleBinHighCells.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));


//-------------------------------------------
  hDummyMInv->GetYaxis()->SetRangeUser(0., hInvMassVsPtHigh_Slice_MC->GetMaximum()*1.3);
  hDummyMInv->DrawCopy();

  // ScaleTo(hInvMassVsPtHigh_Slice_data, hInvMassVsPtHigh_Slice_MC);
  hInvMassVsPtHigh_Slice_MC->Draw("same");
  hInvMassVsPtHigh1cell_Slice_MC->Draw("same");
  hInvMassVsPtHigh2cell_Slice_MC->Draw("same");




  TLegend *leg3 = GetAndSetLegend2(0.7, 0.8, 0.95, 0.95, 40);
  leg3->AddEntry(hInvMassVsPtHigh_Slice_MC, "MC, only high clus", "l");
  leg3->AddEntry(hInvMassVsPtHigh1cell_Slice_MC, "MC, 1 cells", "l");
  leg3->AddEntry(hInvMassVsPtHigh2cell_Slice_MC, "MC, 2 cells", "l");

  leg3->Draw("same");

  drawLatexAdd(sEnergy,0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("2 < p_{T, #gamma} < 3 GeV/c",0.15,0.87,textSizeLabelsRel,kFALSE);

  DrawLines(mean, mean, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kRed + 1, 2);
  DrawLines(mean - 0.05, mean - 0.05, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);
  DrawLines(mean + 0.02, mean + 0.02, 0, hInvMassVsPt_Slice->GetMaximum()*0.3 , 2, kGray+2, 2);


  hDummyCan->SaveAs(Form("%s/%s/ExampleBinHighCellsMC.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));


}

void Effi::PlotNCellRatios(){

  // Normal threshold
  TH1D* hS500A100Data_all = (TH1D*) hNCellVsETMNLS500A100_data->ProjectionX("hS500A100Data_all", 2, 2);
  TH1D* hS500A100MC_all = (TH1D*) hNCellVsETMNLS500A100_MC->ProjectionX("hS500A100MC_all", 2, 2);

  TH1D* hS500A100Data_2cell = (TH1D*) hNCellVsETMNLS500A100_data->ProjectionX("hS500A100Data_2cell", 3, 10);
  TH1D* hS500A100MC_2cell = (TH1D*) hNCellVsETMNLS500A100_MC->ProjectionX("hS500A100MC_2cell", 3, 10);

  TH1D* hS500A100Ratio_all = (TH1D*) hS500A100Data_all->Clone("hS500A100Ratio_all");
  hS500A100Ratio_all->Divide(hS500A100MC_all);
  DrawSetMarker(hS500A100Ratio_all, 1, 3., kBlue + 2, kBlue + 2);

  TH1D* hS500A100Ratio_2cell = (TH1D*) hS500A100Data_2cell->Clone("hS500A100Ratio_2cell");
  hS500A100Ratio_2cell->Divide(hS500A100MC_2cell);
  DrawSetMarker(hS500A100Ratio_2cell, 8, 3., kRed + 2, kRed + 2);

  // S500A105
  TH1D* hS500A105Data_all = (TH1D*) hNCellVsETMNLS500A105_data->ProjectionX("hS500A105Data_all", 2, 2);
  TH1D* hS500A105MC_all = (TH1D*) hNCellVsETMNLS500A105_MC->ProjectionX("hS500A105MC_all", 2, 2);

  TH1D* hS500A105Data_2cell = (TH1D*) hNCellVsETMNLS500A105_data->ProjectionX("hS500A105Data_2cell", 3, 10);
  TH1D* hS500A105MC_2cell = (TH1D*) hNCellVsETMNLS500A105_MC->ProjectionX("hS500A105MC_2cell", 3, 10);

  TH1D* hS500A105Ratio_all = (TH1D*) hS500A105Data_all->Clone("hS500A105Ratio_all");
  hS500A105Ratio_all->Divide(hS500A105MC_all);
  DrawSetMarker(hS500A105Ratio_all, 1, 3., kBlack, kBlack);

  TH1D* hS500A105Ratio_2cell = (TH1D*) hS500A105Data_2cell->Clone("hS500A105Ratio_2cell");
  hS500A105Ratio_2cell->Divide(hS500A105MC_2cell);
  DrawSetMarker(hS500A105Ratio_2cell, 8, 3., kBlack, kBlack);

  hDummyCan->cd();
  hDummyNCellRatio->Draw("");
  hS500A100Ratio_all->Draw("same");
  hS500A100Ratio_2cell->Draw("same");
  // hS500A105Ratio_all->Draw("same");
  hS500A105Ratio_2cell->Draw("same");

  TLegend *leg = GetAndSetLegend2(0.5, 0.15, 0.9, 0.35, 40);
  leg->AddEntry(hS500A100Ratio_all, "N_{cell} = 1, S500A100", "l");
  leg->AddEntry(hS500A100Ratio_2cell, "N_{cell} #geq 2, S500A100", "l");
  // leg->AddEntry(hS500A105Ratio_all, "N_{cell} = 1, S500A105", "l");
  leg->AddEntry(hS500A105Ratio_2cell, "N_{cell} #geq 2, S500A105", "l");
  leg->Draw("same");

  hDummyCan->SaveAs(Form("%s/%s/RatioNCells.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}




void Effi::PlotFraction1cell(){


  // Normal threshold
  TH1D* hS500A100Data_all = (TH1D*) hNCellVsETMNLS500A100_data->ProjectionX("hS500A100Data_all");
  TH1D* hS500A100Data_2cell = (TH1D*) hNCellVsETMNLS500A100_data->ProjectionX("hS500A100Data_2cell", 3, 10);
  hS500A100Data_2cell->Divide(hS500A100Data_all);

  TH1D* hS500A100MC_all = (TH1D*) hNCellVsETMNLS500A100_MC->ProjectionX("hS500A100MC_all");
  TH1D* hS500A100MC_2cell = (TH1D*) hNCellVsETMNLS500A100_MC->ProjectionX("hS500A100MC_2cell", 3, 10);
  hS500A100MC_2cell->Divide(hS500A100MC_all);

  // 500 105
  TH1D* hS500A105Data_all = (TH1D*) hNCellVsETMNLS500A105_data->ProjectionX("hS500A105Data_all");
  TH1D* hS500A105Data_2cell = (TH1D*) hNCellVsETMNLS500A105_data->ProjectionX("hS500A105Data_2cell", 3, 10);
  hS500A105Data_2cell->Divide(hS500A105Data_all);

  TH1D* hS500A105MC_all = (TH1D*) hNCellVsETMNLS500A105_MC->ProjectionX("hS500A105MC_all");
  TH1D* hS500A105MC_2cell = (TH1D*) hNCellVsETMNLS500A105_MC->ProjectionX("hS500A105MC_2cell", 3, 10);
  hS500A105MC_2cell->Divide(hS500A105MC_all);

  // 500 110
  TH1D* hS500A110Data_all = (TH1D*) hNCellVsETMNLS500A110_data->ProjectionX("hS500A110Data_all");
  TH1D* hS500A110Data_2cell = (TH1D*) hNCellVsETMNLS500A110_data->ProjectionX("hS500A110Data_2cell", 3, 10);
  hS500A110Data_2cell->Divide(hS500A110Data_all);

  TH1D* hS500A110MC_all = (TH1D*) hNCellVsETMNLS500A110_MC->ProjectionX("hS500A110MC_all");
  TH1D* hS500A110MC_2cell = (TH1D*) hNCellVsETMNLS500A110_MC->ProjectionX("hS500A110MC_2cell", 3, 10);
  hS500A110MC_2cell->Divide(hS500A110MC_all);


  DrawSetMarker(hS500A100Data_2cell, 20, 2., kRed+2, kRed+2);
  DrawSetMarker(hS500A100MC_2cell, 8, 2., kRed+2, kRed+2);

  DrawSetMarker(hS500A105Data_2cell, 20, 2., kBlack , kBlack);
  DrawSetMarker(hS500A105MC_2cell, 8, 2., kBlack , kBlack);

  DrawSetMarker(hS500A110Data_2cell, 20, 2., kGreen + 2, kGreen + 2);
  DrawSetMarker(hS500A110MC_2cell, 8, 2., kGreen + 2 , kGreen + 2);


  hDummyCan->cd();
  hDummyNCellFraction->Draw("");
  hS500A100Data_2cell->Draw("same,p");
  hS500A100MC_2cell->Draw("same,hist");
  hS500A105Data_2cell->Draw("same,p");
  hS500A105MC_2cell->Draw("same,hist");
  hS500A110Data_2cell->Draw("same,p");
  hS500A110MC_2cell->Draw("same,hist");

  TLegend *leg = GetAndSetLegend2(0.5, 0.15, 0.9, 0.35, 40);
  leg->AddEntry(hS500A100Data_2cell, "data, S500A100", "p");
  leg->AddEntry(hS500A100MC_2cell, "MC, S500A100", "l");
  leg->AddEntry(hS500A105Data_2cell, "data, S500A105", "p");
  leg->AddEntry(hS500A105MC_2cell, "MC, S500A105", "l");
  leg->AddEntry(hS500A110Data_2cell, "data, S500A110", "p");
  leg->AddEntry(hS500A110MC_2cell, "MC, S500A110", "l");
  leg->Draw("same");

  hDummyCan->SaveAs(Form("%s/%s/NCell1Fraction.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}



void Effi::PlotNCellDistr(){


  // Normal threshold
  TH1D* hS500A100Data_all = (TH1D*) hNCellVsETMNLS500A100_data->ProjectionX("hS500A100Data_all");
  TH1D* hS500A100Data_2cell = (TH1D*) hNCellVsETMNLS500A100_MC->ProjectionX("hS500A100MC_all", 3, 10);
  TH1D* hS500A100MC_all = (TH1D*) hNCellVsETMNLS500A100_MC->ProjectionX("hS500A100MC_all");
  TH1D* hS500A100MC_2cell = (TH1D*) hNCellVsETMNLS500A100_MC->ProjectionX("hS500A100MC_2cell", 3, 10);


  DrawSetMarker(hS500A100Data_all, 20, 2., kRed+2, kRed+2);
  DrawSetMarker(hS500A100Data_2cell, 8, 2., kRed+2, kRed+2);

  DrawSetMarker(hS500A100MC_all, 20, 2., kBlack , kBlack);
  DrawSetMarker(hS500A100MC_2cell, 8, 2., kBlack , kBlack);


  hDummyCan->cd();
  // hDummyCan->SetLogy();
  hS500A100Data_all->Draw("p");
  hS500A100Data_2cell->Draw("same,hist");
  hS500A100MC_all->Draw("same,p");
  hS500A100MC_2cell->Draw("same,hist");

  TLegend *leg = GetAndSetLegend2(0.5, 0.15, 0.9, 0.35, 40);
  leg->AddEntry(hS500A100Data_all, "data, all cells", "p");
  leg->AddEntry(hS500A100Data_2cell, "data N_{#geq 2}", "l");
  leg->AddEntry(hS500A100MC_all, "MC, all cells", "p");
  leg->AddEntry(hS500A100MC_2cell, "MC N_{#geq 2}", "l");
  leg->Draw("same");

  hDummyCan->SaveAs(Form("%s/%s/NCellDistribution.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));

}


void Effi::FitCorr_GammasWide(){
    if(!hNCell_GammasWide_Effi_Corr){
      cout<<"hNCell_GammasWide_Effi_Corr does not exit... returning\n";
    }
    TF1 *funcGaus = new TF1("funcGaus", "gaus(0)", 0, 5);
    TF1 *funcLinear = new TF1("funcLinear", "pol1(0)", 0, 5);
    hNCell_GammasWide_Effi_Corr->Fit(funcGaus, "MR0");
    hNCell_GammasWide_Effi_Corr->Fit(funcLinear, "MR0");
    cout<<"###############################################"<<endl;
    cout<<"#### Printing parameters for gaussian fit #####"<<endl;
    cout<<"### Applied to gammas + conv. in wide range ###"<<endl;
    for(int i = 0; i < 3; ++i){
      cout<<"param["<<i<<"] = "<<funcGaus->GetParameter(i)<<endl;
    }
    cout<<"###############################################"<<endl;
    cout<<"##### Printing parameters for linear fit ######"<<endl;
    cout<<"### Applied to gammas + conv. in wide range ###"<<endl;
    for(int i = 0; i < 2; ++i){
      cout<<"param["<<i<<"] = "<<funcLinear->GetParameter(i)<<endl;
    }
    cout<<"###############################################"<<endl;

    hDummyCan->cd();
    hDummyCorr->Draw();
    DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
    hNCell_GammasWide_Effi_Corr->Draw("same,p");
    funcGaus->Draw("same");
    funcLinear->SetLineColor(kOrange + 2);
    funcLinear->Draw("same");

    TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
    leg->AddEntry(hNCell_GammasWide_Effi_Corr, "P2 data", "p");
    leg->AddEntry(funcGaus, "gaussian fit", "l");
    leg->AddEntry(funcLinear, "linear fit", "l");
    leg->Draw("same");

    drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);

    hDummyCan->SaveAs(Form("%s/%s/Fit_GammasWide_Gaus.%s", fPeriod.Data(), fsuffix.Data(), fsuffix.Data()));




}




void PlottingNCellEffi(){

  {

    gSystem->Exec("mkdir -p 13TeVNomB/pdf");
    gSystem->Exec("mkdir -p 13TeVNomB/png");
    // Effi plot("Pi0Tagging_13TeV_nom_02_28_TBNL.root", "13TeVNomB", "pdf");
    Effi plot("Pi0Tagging_13TeV_nom_03_10_noIso_TBNL.root", "13TeVNomB", "pdf");
    plot.FillHistos();
    plot.SetOtherHistos();
    plot.DoSBSubtraction();
    plot.GetHistReweighted();
    plot.FillCorrHistos();
    cout<<__LINE__<<endl;
    plot.LoadTB();
    cout<<__LINE__<<endl;
    plot.SetPlotting();
    plot.PlotEffi_AllClusAndTB();
    plot.PlotEffi_AllGammaAndTB();
    plot.PlotEffi_RWGammaAndTB();
    plot.PlotEffi_SBGammaAndTB();
    plot.PlotEffi_RWGammaHighAndTB();
    plot.PlotEffi_RWGammaHighAndTB_MCCheck();
    plot.PlotEffi_TrueAndReweighted();
    plot.PlotRatio_GammasAndRW();
    plot.PlotRatio_GammasAndRWLeft();
    plot.PlotRatio_GammasAndRWWide();
    plot.PlotRatio_GammasAndConvRWWide();
    plot.PlotRatio_GammasAndRWHigh();
    plot.PlotRatio_GammasAndRWHighMCSBSub();
    // plot.PlotRatio_GammasAndRWWide();
    plot.PlotRatio_SBGammasAndRW();
    plot.PlotEffi_HighAllMethods();
    plot.PlotEffi_HighAllMethodsMCSBSub();
    plot.PlotEffi_HighAndWideRWSBSub();
    plot.PlotEffi_HighAndWideRWSBSub_Light();
    plot.PlotEffi_HighAndWideRWSBSub_Smoothed();
    plot.PlotEffi_WideOnlySBSub();
    plot.PlotEffi_WideSBSubGammaAndConv();
    plot.PlotEffi_HighAndWide_Light();
    plot.PlotEffi_WideSmoothed();
    plot.PlotEffi_Electrons();
    plot.PlotCorr_GammasAndRW();
    plot.PlotCorr_GammasAndRW_Left();
    plot.PlotCorr_GammasAndRW_Wide();
    plot.PlotCorr_GammasAndConv_Wide();
    plot.PlotCorr_GammasAndRW_High();
    plot.PlotCorr_HighAllMethods();
    plot.PlotCorr_Electrons();
    plot.PlotPurity();
    plot.PlotPurity2();
    plot.PlotExampleBin();
    plot.PlotExampleBinCells();
    plot.PlotNCellRatios();
    plot.PlotFraction1cell();
    plot.PlotNCellDistr();
    plot.FitCorr_GammasWide();
  }

  // gSystem->Exec("mkdir -p 13TeVLowB/pdf");
  // gSystem->Exec("mkdir -p 13TeVLowB/png");
  // Effi plot("Pi0Tagging_13TeV_nom_02_25_TBNL.root", "13TeVLowB", "pdf");
  // plot.FillHistos();
  // plot.SetOtherHistos();
  // plot.DoSBSubtraction();
  // plot.GetHistReweighted();
  // plot.FillCorrHistos();
  // cout<<__LINE__<<endl;
  // plot.LoadTB();
  // cout<<__LINE__<<endl;
  // plot.SetPlotting();
  // plot.PlotEffi_AllClusAndTB();
  // plot.PlotEffi_AllGammaAndTB();
  // plot.PlotEffi_RWGammaAndTB();
  // plot.PlotEffi_SBGammaAndTB();
  // plot.PlotEffi_RWGammaHighAndTB();
  // plot.PlotEffi_RWGammaHighAndTB_MCCheck();
  // plot.PlotEffi_TrueAndReweighted();
  // plot.PlotRatio_GammasAndRW();
  // plot.PlotRatio_GammasAndRWLeft();
  // plot.PlotRatio_GammasAndRWWide();
  // plot.PlotRatio_GammasAndRWHigh();
  // plot.PlotRatio_GammasAndRWHighMCSBSub();
  // // plot.PlotRatio_GammasAndRWWide();
  // plot.PlotRatio_SBGammasAndRW();
  // plot.PlotEffi_HighAllMethods();
  // plot.PlotEffi_HighAllMethodsMCSBSub();
  // plot.PlotEffi_HighAndWideRWSBSub();
  // plot.PlotEffi_HighAndWideRWSBSub_Light();
  // plot.PlotCorr_GammasAndRW();
  // plot.PlotCorr_GammasAndRW_Left();
  // plot.PlotCorr_GammasAndRW_Wide();
  // plot.PlotCorr_GammasAndRW_High();
  // plot.PlotCorr_HighAllMethods();
  // plot.PlotPurity();
  // plot.PlotPurity2();
  // plot.PlotExampleBin();
  // plot.PlotExampleBinCells();
  // plot.PlotNCellRatios();
  // plot.PlotFraction1cell();
  // plot.PlotNCellDistr();

  // {
  //   gSystem->Exec("mkdir 8TeV");
  //   Effi plot("Pi0Tagging_8TeV_02_08.root", "8TeV");
  //   plot.FillHistos();
  //   plot.SetOtherHistos();
  //   plot.DoSBSubtraction();
  //   plot.GetHistReweighted();
  //   plot.FillCorrHistos();
  //   plot.LoadTB();
  //   plot.SetPlotting();
  //   plot.PlotEffi_AllClusAndTB();
  //   plot.PlotEffi_AllGammaAndTB();
  //   plot.PlotEffi_RWGammaAndTB();
  //   plot.PlotEffi_SBGammaAndTB();
  //   plot.PlotEffi_RWGammaHighAndTB();
  //   plot.PlotEffi_RWGammaHighAndTB_MCCheck();
  //   plot.PlotEffi_TrueAndReweighted();
  //   plot.PlotRatio_GammasAndRW();
  //   plot.PlotRatio_GammasAndRWLeft();
  //   plot.PlotRatio_GammasAndRWWide();
  //   plot.PlotRatio_GammasAndRWHigh();
  //   plot.PlotRatio_GammasAndRWHighMCSBSub();
  //   plot.PlotRatio_SBGammasAndRW();
  //   plot.PlotEffi_HighAllMethods();
  //   plot.PlotEffi_HighAllMethodsMCSBSub();
  //   plot.PlotCorr_GammasAndRW();
  //   plot.PlotCorr_GammasAndRW_Left();
  //   plot.PlotCorr_GammasAndRW_Wide();
  //   plot.PlotCorr_GammasAndRW_High();
  //   plot.PlotCorr_HighAllMethods();
  //   plot.PlotPurity();
  //   plot.PlotPurity2();
  //   plot.PlotExampleBin();
  //   plot.PlotExampleBinCells();
  //   plot.PlotNCellRatios();
  //   plot.PlotFraction1cell();
  //   plot.PlotNCellDistr();
  // }
}
