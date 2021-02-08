#include "Plotting_Header.h"
// #include "/home/joshua/PCG_Software/Plotting/Plotting_Class.h"
#include "TGraph.h"


class Effi{

public:
  Effi();
  Effi(TString input, TString Period = "13TeV");
  ~Effi();
  void FillHistos();
  float calcErr(float Emc, float Errmc, float Edata, float Errdata);
  void GetEffiHists(TH2F *hdata2d, TH2F *hMC2d, TH1D *&hEffiData,  TH1D *&hEffiMC, TH1D *&hRatio, TH1D *&hCorr);
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

  // Ratio Plots
  void PlotRatio_GammasAndRW();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_GammasAndRWLeft();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_GammasAndRWWide();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_GammasAndRWHigh();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_GammasAndRWHighMCSBSub();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotRatio_SBGammasAndRW();   // Plot ratio with TB + all clusters + gammas (right) + gammas (right) reweighted

  // Corr Plots
  void PlotCorr_GammasAndRW();   // Plot correction factor with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotCorr_GammasAndRW_Left();   // Plot correction factor with TB + all clusters + gammas (left) + gammas (left) reweighted
  void PlotCorr_GammasAndRW_Wide();   // Plot correction factor with TB + all clusters + gammas (left) + gammas (left) reweighted
  void PlotCorr_GammasAndRW_High();   // Plot correction factor with TB + all clusters + gammas (left) + gammas (left) reweighted
  void PlotCorr_HighAllMethods();   // Plot correction factor with TB + all clusters + gammas (left) + gammas (left) reweighted

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

private:

  TFile* fdata = nullptr;

  TString fPeriod = "";

  TString sEnergy = "pp #sqrt{s} = 13 TeV";

  // Test Beam
  TGraph *grTB_data = nullptr;
  TGraph *grTB_MC = nullptr;
  TGraph *grTB_Ratio = nullptr;
  TGraph *grTB_Corr = nullptr;

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

  TH2F* hNCellVsEGammasNL_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLLeft_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLWide_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_SBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_MCSBSub_data = nullptr;
  TH2F* hNCellVsEGammasNLSB_RW_data = nullptr;
  TH2F* hNCellVsEGammasNLSBHigh_RW_data = nullptr;

  TH2F* hNCellVsEGammasNL_RW_MC = nullptr;
  TH2F* hNCellVsEGammasNLLeft_RW_MC = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_MC = nullptr;
  TH2F* hNCellVsEGammasNLHigh_RW_SBSub_MC = nullptr;
  TH2F* hNCellVsEGammasNLSBHigh_RW_MC = nullptr;

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

  // gammas reweighted left
  TH1D* hNCell_GammasWide_RW_Effi_data = nullptr;
  TH1D* hNCell_GammasWide_RW_Effi_MC = nullptr;
  TH1D* hNCell_GammasWide_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasWide_RW_Effi_Corr = nullptr;

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

  TH2F* hInvMassVsHighGammaPtBack_data = nullptr;
  TH2F* hInvMassVsHighGammaPtBack_MC = nullptr;

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
Effi::Effi(TString input, TString Period){
  fdata = TFile::Open(input);
  fPeriod = Period;
}

void Effi::LoadTB(){

  TFile TB("NCellEfficiencies.root");
  grTB_data = (TGraph*) TB.Get("Effi_Data_TB15");

  // for 8TeV scale it back
  // if(fPeriod.Contains("8TeV")) {
  for(int i = 0; i < grTB_data->GetN(); ++i){
    cout<<"grTB_data->GetPointX(i)"<<grTB_data->GetPointX(i)<<endl;
    grTB_data->SetPointX(i, grTB_data->GetPointX(i)*(1./1.015));
  }
  // }

  grTB_MC = (TGraph*) TB.Get("Effi_MC_TB");

  grTB_Ratio = (TGraph*) grTB_MC->Clone("grTBRatio");
  for(int i = 0; i < grTB_Ratio->GetN(); ++i){
    grTB_Ratio->SetPointY(i, grTB_data->Eval(grTB_MC->GetPointX(i)) / grTB_MC->GetPointY(i));
  }

  grTB_Corr = (TGraph*) grTB_MC->Clone("grTBEffi");
  for(int i = 0; i < grTB_data->GetN(); ++i){
    float CF = -1;
    if(grTB_MC->GetPointY(i) < 1){
      CF = 1 - (( 1 - grTB_data->Eval(grTB_MC->GetPointX(i)))/(1 - grTB_MC->GetPointY(i)));
    } else{
      CF = -1;
    }
    grTB_Corr->SetPointY(i, CF);
  }
}

void Effi::FillHistos(){
  hNCellVsETMNL_data = (TH2F*) fdata->Get("hNCellVsETMNL_data");
  hNCellVsETMNL_MC = (TH2F*) fdata->Get("hNCellVsETMNL_MC");

  hNCellVsETMNLS500A100_data = (TH2F*) fdata->Get("hNCellVsETMNLS500A100_data");
  hNCellVsETMNLS500A100_MC = (TH2F*) fdata->Get("hNCellVsETMNLS500A100_MC");

  hNCellVsETMNLS500A105_data = (TH2F*) fdata->Get("hNCellVsETMNLS500A105_data");
  hNCellVsETMNLS500A105_MC = (TH2F*) fdata->Get("hNCellVsETMNLS500A105_MC");

  hNCellVsETMNLS500A110_data = (TH2F*) fdata->Get("hNCellVsETMNLS500A110_data");
  hNCellVsETMNLS500A110_MC = (TH2F*) fdata->Get("hNCellVsETMNLS500A110_MC");

  hNCellVsEGammasNL_data = (TH2F*) fdata->Get("hNCellVsEGammasNL_data");
  hNCellVsEGammasNL_MC = (TH2F*) fdata->Get("hNCellVsEGammasNL_MC");
  hNCellVsETrueGammasNL_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNL_MC");
  hNCellVsEGammasNLTrueElec_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElec_MC");

  hNCellVsEGammasNLLeft_data = (TH2F*) fdata->Get("hNCellVsEGammasNLLeft_data");
  hNCellVsEGammasNLLeft_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLLeft_MC");
  hNCellVsETrueGammasNLLeft_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLLeft_MC");
  hNCellVsEGammasNLTrueElecLeft_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecLeft_MC");

  hNCellVsEGammasNLWide_data = (TH2F*) fdata->Get("hNCellVsEGammasNLWide_data");
  hNCellVsEGammasNLWide_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLWide_MC");
  hNCellVsETrueGammasNLWide_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLWide_MC");
  hNCellVsEGammasNLTrueElecWide_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecWide_MC");

  hNCellVsEGammasNLSB_data = (TH2F*) fdata->Get("hNCellVsEGammasNLSideBand_data");
  hNCellVsEGammasNLSB_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLSideBand_MC");
  hNCellVsETrueGammasNLSB_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLSideBand_MC");
  hNCellVsEGammasNLTrueElecSB_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecSideBand_MC");

  hNCellVsEGammasNLSBHigh_data = (TH2F*) fdata->Get("hNCellVsEGammasNLSideBandOnlyHighClus_data");
  hNCellVsEGammasNLSBHigh_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLSideBandOnlyHighClus_MC");
  hNCellVsETrueGammasNLSBHigh_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLSideBandOnlyHighClus_MC");
  hNCellVsEGammasNLTrueElecSBHigh_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecSideBandOnlyHighClus_MC");

  hNCellVsEGammasNLHigh_data = (TH2F*) fdata->Get("hNCellVsEGammasNLOnlyHighClus_data");
  hNCellVsEGammasNLHigh_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLOnlyHighClus_MC");
  hNCellVsETrueGammasNLHigh_MC = (TH2F*) fdata->Get("hNCellVsETrueGammasNLOnlyHighClus_MC");
  hNCellVsEGammasNLTrueElecHigh_MC = (TH2F*) fdata->Get("hNCellVsEGammasNLTrueElecOnlyHighClus_MC");


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
  // all clusters
  GetEffiHists(hNCellVsETMNL_data, hNCellVsETMNL_MC, hNCell_AllClus_Effi_data, hNCell_AllClus_Effi_MC, hNCell_AllClus_Effi_Ratio, hNCell_AllClus_Effi_Corr);
  // gammas
  GetEffiHists(hNCellVsEGammasNL_data, hNCellVsEGammasNL_MC, hNCell_Gammas_Effi_data, hNCell_Gammas_Effi_MC, hNCell_Gammas_Effi_Ratio, hNCell_Gammas_Effi_Corr);
  // gammas reweighted
  GetEffiHists(hNCellVsEGammasNL_RW_data, hNCellVsETrueGammasNL_MC, hNCell_Gammas_RW_Effi_data, hNCell_Gammas_RW_Effi_MC, hNCell_Gammas_RW_Effi_Ratio, hNCell_Gammas_RW_Effi_Corr);

  // gammas reweighted MC
  GetEffiHists(hNCellVsEGammasNL_RW_MC, hNCellVsETrueGammasNL_MC, hNCell_Gammas_RW_Effi_MC_consistency, hNCell_Gammas_RW_Effi_MC_consistency, hNCell_Gammas_RW_Effi_Ratio_consistency, hNCell_Gammas_RW_Effi_Corr_consistency);

  // gammas left
  GetEffiHists(hNCellVsEGammasNLLeft_data, hNCellVsEGammasNLLeft_MC, hNCell_GammasLeft_Effi_data, hNCell_GammasLeft_Effi_MC, hNCell_GammasLeft_Effi_Ratio, hNCell_GammasLeft_Effi_Corr);
  // gammas left reweighted
  GetEffiHists(hNCellVsEGammasNLLeft_RW_data, hNCellVsETrueGammasNLLeft_MC, hNCell_GammasLeft_RW_Effi_data, hNCell_GammasLeft_RW_Effi_MC, hNCell_GammasLeft_RW_Effi_Ratio, hNCell_GammasLeft_RW_Effi_Corr);

  // gammas Wide
  GetEffiHists(hNCellVsEGammasNLWide_data, hNCellVsEGammasNLWide_MC, hNCell_GammasWide_Effi_data, hNCell_GammasWide_Effi_MC, hNCell_GammasWide_Effi_Ratio, hNCell_GammasWide_Effi_Corr);
  // gammas wide reweighted
  GetEffiHists(hNCellVsEGammasNLWide_RW_data, hNCellVsETrueGammasNLWide_MC, hNCell_GammasWide_RW_Effi_data, hNCell_GammasWide_RW_Effi_MC, hNCell_GammasWide_RW_Effi_Ratio, hNCell_GammasWide_RW_Effi_Corr);

  // gammas Side band
  GetEffiHists(hNCellVsEGammasNLSB_data, hNCellVsEGammasNLSB_MC, hNCell_GammasSB_Effi_data, hNCell_GammasSB_Effi_MC, hNCell_GammasSB_Effi_Ratio, hNCell_GammasSB_Effi_Corr);
  // gammas Side band reweighted
  GetEffiHists(hNCellVsEGammasNLSB_RW_data, hNCellVsETrueGammasNLSB_MC, hNCell_GammasSB_RW_Effi_data, hNCell_GammasSB_RW_Effi_MC, hNCell_GammasSB_RW_Effi_Ratio, hNCell_GammasSB_RW_Effi_Corr);

  // gammas Side band highest cluster
  GetEffiHists(hNCellVsEGammasNLSBHigh_data, hNCellVsEGammasNLSBHigh_MC, hNCell_GammasSBHigh_Effi_data, hNCell_GammasSBHigh_Effi_MC, hNCell_GammasSBHigh_Effi_Ratio, hNCell_GammasSBHigh_Effi_Corr);
  // gammas Side band reweighted highest cluster
  GetEffiHists(hNCellVsEGammasNLSBHigh_RW_data, hNCellVsETrueGammasNLSBHigh_MC, hNCell_GammasSBHigh_RW_Effi_data, hNCell_GammasSBHigh_RW_Effi_MC, hNCell_GammasSBHigh_RW_Effi_Ratio, hNCell_GammasSBHigh_RW_Effi_Corr);

  // gammas High
  GetEffiHists(hNCellVsEGammasNLHigh_data, hNCellVsEGammasNLHigh_MC, hNCell_GammasHigh_Effi_data, hNCell_GammasHigh_Effi_MC, hNCell_GammasHigh_Effi_Ratio, hNCell_GammasHigh_Effi_Corr);
  // gammas High reweighted
  GetEffiHists(hNCellVsEGammasNLHigh_RW_data, hNCellVsETrueGammasNLHigh_MC, hNCell_GammasHigh_RW_Effi_data, hNCell_GammasHigh_RW_Effi_MC, hNCell_GammasHigh_RW_Effi_Ratio, hNCell_GammasHigh_RW_Effi_Corr);
  // gammas High reweighted sideband subtracted
  GetEffiHists(hNCellVsEGammasNLHigh_RW_SBSub_data, hNCellVsEGammasNLHigh_RW_SBSub_MC, hNCell_GammasHigh_RW_SBSub_Effi_data, hNCell_GammasHigh_RW_SBSub_Effi_MC, hNCell_GammasHigh_RW_SBSub_Effi_Ratio, hNCell_GammasHigh_RW_SBSub_Effi_Corr);

  // gammas High reweighted  MC sideband subtracted
  GetEffiHists(hNCellVsEGammasNLHigh_RW_MCSBSub_data, hNCellVsEGammasNLHigh_RW_SBSub_MC, hNCell_GammasHigh_RW_MCSBSub_Effi_data, hNCell_GammasHigh_RW_MCSBSub_Effi_MC, hNCell_GammasHigh_RW_MCSBSub_Effi_Ratio, hNCell_GammasHigh_RW_MCSBSub_Effi_Corr);

  // true gamma in pi0 range?
  GetTrueHists(hNCellVsETrueGammasNL_MC, hNCell_TrueGammas_Effi_MC, "trueGammas");

  // true elec in pi0 range?
  GetTrueHists(hNCellVsEGammasNLTrueElec_MC, hNCell_TrueElec_Effi_MC, "trueElecs");
}


// -------------------------------------------------
// Steering for sideband subtraction
// -------------------------------------------------
void Effi::DoSBSubtraction(){
  // data high gamma
  hNCellVsEGammasNLHigh_SBSub_data = SubtractSidebandBack( hInvMassVsPtHigh_data, hInvMassVsHighGammaPtBack_data , hNCellVsEGammasNLSBHigh_data , hNCellVsEGammasNLHigh_data);
  // MC high gamma
  hNCellVsEGammasNLHigh_SBSub_MC = SubtractSidebandBack( hInvMassVsPtHigh_MC, hInvMassVsHighGammaPtBack_MC , hNCellVsEGammasNLSBHigh_MC , hNCellVsEGammasNLHigh_MC);
  // data high gamma
  hNCellVsEGammasNLHigh_MCSBSub_data = SubtractSidebandBack( hInvMassVsPtHigh_data, hInvMassVsHighGammaPtBack_data , hNCellVsEGammasNLSBHigh_MC , hNCellVsEGammasNLHigh_data);

}

// -------------------------------------------------
// SUBTRACT SIDEBAND FROM SIGNAL REGION
// Make Sure the sideband and the peak use the same method or gamma selection!
// -------------------------------------------------
TH2F* Effi::SubtractSidebandBack(TH2F* hMInv, TH2F* hMInvBack, TH2F* hNCellSB, TH2F* hNCell){


  // STrategy:
  // Fill a histo with a pT dependent scaling factors
  // this only works if the inv mass isto is only filled with the pT of the higher photon!
  TFile fMassPos("fMassPos.root");
  TGraph *grMassPos = nullptr;
  if(fMassPos.IsZombie()){
    TString name = (fPeriod.Contains("13") ? "MassPos13TeV_func" : "MassPos8TeV_func");
    grMassPos = (TGraph*) fMassPos.Get(name);
  }
  // MassPos8TeV_func
  TH2F* hScale = (TH2F*) hNCell->Clone("hScale");
  hScale->Reset();
  for(int i = 1; i <= hNCell->GetNbinsX(); ++i){
    float xValLow = hNCell->GetXaxis()->GetBinLowEdge(i) + 0.001;
    float xValUp = hNCell->GetXaxis()->GetBinUpEdge(i) - 0.001;

    TH1D* hMInv_Proj = (TH1D*) hMInv->ProjectionX("hMInv_Proj", hMInv->GetYaxis()->FindBin(xValLow), hMInv->GetYaxis()->FindBin(xValUp));
    TH1D* hMInvBack_Proj = (TH1D*) hMInvBack->ProjectionX("hMInvBack_Proj", hMInvBack->GetYaxis()->FindBin(xValLow), hMInvBack->GetYaxis()->FindBin(xValUp));
    cout<<"hMInv_Proj->GetBinContent(20): "<<hMInv_Proj->GetBinContent(20)<<endl;
    float scaleFac = hMInv_Proj->Integral(hMInv_Proj->FindBin(0.17), hMInv_Proj->FindBin(0.3)) / hMInvBack_Proj->Integral(hMInvBack_Proj->FindBin(0.17), hMInvBack_Proj->FindBin(0.3));
    hMInvBack_Proj->Scale(scaleFac);

    // ranges for background and signal make that better??
    float rangeSB[2] = {0.135 + 0.05, 0.135 + 0.2};
    float rangeSignal[2] = {0.08, 0.15};
    if(grMassPos){
      rangeSB[0] = grMassPos->Eval(0.5*(xValLow + xValUp )) + 0.05;
      rangeSB[1] = grMassPos->Eval(0.5*(xValLow + xValUp )) + 0.2;
      rangeSignal[0] = grMassPos->Eval(0.5*(xValLow + xValUp )) - 0.02;
      rangeSignal[1] = grMassPos->Eval(0.5*(xValLow + xValUp )) + 0.05;
    }

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
  // scale NCell vs E distribution to the Signal region
  hNCellSB->Multiply(hScale);
  // hNCellSB->Scale(backInt/SBInt);
  cout<<__LINE__<<endl;
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
  hInvMassVsPt_data = (TH2F*) fdata->Get("hInvMassVsPt_data");
  hInvMassVsPtGG = (TH2F*) fdata->Get("hInvMassVsPtGG_MC");
  hInvMassVsPtGC = (TH2F*) fdata->Get("hInvMassVsPtGC_MC");
  hInvMassVsPtCC = (TH2F*) fdata->Get("hInvMassVsPtCC_MC");
  float minpt = 2.;
  float maxpt = 3.;
  hInvMassVsPt_Slice = (TH1D*) hInvMassVsPt_MC->ProjectionX("hInvMassVsPt_Slice", hInvMassVsPt_MC->GetYaxis()->FindBin(minpt), hInvMassVsPt_MC->GetYaxis()->FindBin(maxpt));
  hInvMassVsPt_data_Slice = (TH1D*) hInvMassVsPt_data->ProjectionX("hInvMassVsPt_data_Slice", hInvMassVsPt_MC->GetYaxis()->FindBin(minpt), hInvMassVsPt_MC->GetYaxis()->FindBin(maxpt));
  hInvMassVsPtGG_Slice = (TH1D*) hInvMassVsPtGG->ProjectionX("hInvMassVsPtGG_Slice", hInvMassVsPtGG->GetYaxis()->FindBin(minpt), hInvMassVsPtGG->GetYaxis()->FindBin(maxpt));
  hInvMassVsPtGC_Slice = (TH1D*) hInvMassVsPtGC->ProjectionX("hInvMassVsPtGC_Slice", hInvMassVsPtGC->GetYaxis()->FindBin(minpt), hInvMassVsPtGC->GetYaxis()->FindBin(maxpt));
  hInvMassVsPtCC_Slice = (TH1D*) hInvMassVsPtCC->ProjectionX("hInvMassVsPtCC_Slice", hInvMassVsPtCC->GetYaxis()->FindBin(minpt), hInvMassVsPtCC->GetYaxis()->FindBin(maxpt));

  // load background hists
  hInvMassVsPtBack_data = (TH2F*) fdata->Get("hInvMassVsPtBack_data");
  hInvMassVsPtBack_MC = (TH2F*) fdata->Get("hInvMassVsPtBack_MC");


  hInvMassVsPtHigh_MC = (TH2F*) fdata->Get("hInvMassVsPtHigh_MC");
  hInvMassVsHighGammaPtBack_MC = (TH2F*) fdata->Get("hInvMassVsHighGammaPtBack_MC");
  hInvMassVsPtHighGamma_MC = (TH2F*) fdata->Get("hInvMassVsPtHighGamma_MC");
  hInvMassVsPtHighElec_MC = (TH2F*) fdata->Get("hInvMassVsPtHighElec_MC");
  hInvMassVsPtHigh1cell_MC = (TH2F*) fdata->Get("hInvMassVsPtHigh1cell_MC");
  hInvMassVsPtHigh2cell_MC = (TH2F*) fdata->Get("hInvMassVsPtHigh2cell_MC");
  hInvMassVsPtHigh3cell_MC = (TH2F*) fdata->Get("hInvMassVsPtHigh3cell_MC");
  hInvMassVsPtHigh_data = (TH2F*) fdata->Get("hInvMassVsPtHigh_data");
  hInvMassVsHighGammaPtBack_data = (TH2F*) fdata->Get("hInvMassVsHighGammaPtBack_data");
  hInvMassVsPtHigh1cell_data = (TH2F*) fdata->Get("hInvMassVsPtHigh1cell_data");
  hInvMassVsPtHigh2cell_data = (TH2F*) fdata->Get("hInvMassVsPtHigh2cell_data");
  hInvMassVsPtHigh3cell_data = (TH2F*) fdata->Get("hInvMassVsPtHigh3cell_data");

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
void Effi::GetEffiHists(TH2F *hdata2d, TH2F *hMC2d, TH1D *&hEffiData,  TH1D *&hEffiMC, TH1D *&hRatio, TH1D *&hCorr){

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
  // TH2F* hElecFracHigh = (TH2F*) hNCellVsEGammasNLTrueElecHigh_MC->Clone("hElecFracHigh");
  hElecFracHigh->Divide(hNCellVsEGammasNLHigh_MC);
  hElecFracHigh->Multiply(hNCellVsEGammasNLHigh_data);
  hNCellVsEGammasNLHigh_RW_SBSub_data->Add(hElecFracHigh, -1);


  hNCellVsEGammasNLHigh_RW_MCSBSub_data = (TH2F*) hNCellVsEGammasNLHigh_MCSBSub_data->Clone("hNCellVsEGammasNLHigh_RW_MCSBSub_data");
  // get electron fraction
  // TH2F* hElecFracHigh = (TH2F*) hNCellVsEGammasNLTrueElecHigh_MC->Clone("hElecFracHigh");
  hElecFracHigh->Divide(hNCellVsEGammasNLHigh_MC);
  hElecFracHigh->Multiply(hNCellVsEGammasNLHigh_data);
  hNCellVsEGammasNLHigh_RW_MCSBSub_data->Add(hElecFracHigh, -1);


  hNCellVsEGammasNLHigh_RW_SBSub_MC = (TH2F*) hNCellVsEGammasNLHigh_SBSub_MC->Clone("hNCellVsEGammasNLHigh_RW_SBSub_MC");
  // get electron fraction
  // TH2F* hElecFracHigh = (TH2F*) hNCellVsEGammasNLTrueElecHigh_MC->Clone("hElecFracHigh");
  // hElecFracHigh->Divide(hNCellVsEGammasNLHigh_MC);
  // hElecFracHigh->Multiply(hNCellVsEGammasNLHigh_data);
  hNCellVsEGammasNLHigh_RW_SBSub_MC->Add(hElecFracHigh, -1);
  //
  // TH1D* hElecPurityRecalc = hElecPurityHigh->Clone("hElecPurityRecalc");
  // hElecPurityRecalc->Scale();
  // hElecPuritySB



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


  hDummyEffi    = new TH2F("hDummyEffi","hDummyEffi",1000,0.5, 10,1000,0., 1.4);
  SetStyleHistoTH2ForGraphs(hDummyEffi, "#it{E} (GeV)","#varepsilon", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyRatio    = new TH2F("hDummyRatio","hDummyRatio",1000,0.5, 10,1000,0.9, 1.4);
  SetStyleHistoTH2ForGraphs(hDummyRatio, "#it{E} (GeV)","#varepsilon_{MC}/#varepsilon_{data}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyCorr    = new TH2F("hDummyCorr","hDummyCorr",1000,0.5, 10,1000,0., 1.4);
  SetStyleHistoTH2ForGraphs(hDummyCorr, "#it{E} (GeV)","#rho", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyPurity    = new TH2F("hDummyPurity","hDummyPurity",1000,0.5, 10,1000,0., 1.1);
  SetStyleHistoTH2ForGraphs(hDummyPurity, "#it{E} (GeV)","P", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyMInv    = new TH2F("hDummyMInv","hDummyMInv",1000,0.01, 0.4,1000,0., 10000);
  SetStyleHistoTH2ForGraphs(hDummyMInv, "#it{M}_{inv} (GeV/c^{2})","counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyNCellRatio    = new TH2F("hDummyNCellRatio","hDummyNCellRatio",1000,0.7, 10,1000,0.75, 1.45);
  SetStyleHistoTH2ForGraphs(hDummyNCellRatio, "#it{E}_{clus} (GeV/c^{2})","N_{clus}^{data}/N_{clus}^{MC}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyNCellFraction    = new TH2F("hDummyNCellFraction","hDummyNCellFraction",1000,0.4, 10,1000,0., 1.);
  SetStyleHistoTH2ForGraphs(hDummyNCellFraction, "#it{E}_{clus} (GeV/c^{2})","N_{clus}^{Ncell = 1}/N_{clus}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);
  // gPad->SetLogx();
  hDummyCan = new TCanvas("Can", "", 1200, 1000);
  DrawPaperCanvasSettings(hDummyCan, 0.1, 0.003, 0.003, 0.1);


  DrawSetMarkerTGraph(grTB_data, 21, 1.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraph(grTB_MC, 25, 2.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraph(grTB_Ratio, 21, 1.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraph(grTB_Corr, 21, 1.5, kBlue + 2, kBlue + 2, 2);

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

  int colGammaLeft = kGreen + 2;
  DrawSetMarker(hNCell_GammasLeft_Effi_data, 34, 2, colGammaLeft, colGammaLeft);
  DrawSetMarker(hNCell_GammasLeft_Effi_MC, 8, 3, colGammaLeft, colGammaLeft); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasLeft_Effi_Ratio, 34, 2, colGammaLeft, colGammaLeft);
  DrawSetMarker(hNCell_GammasLeft_Effi_Corr, 34, 2, colGammaLeft, colGammaLeft);

  int colGammaWide = kGreen + 2;
  DrawSetMarker(hNCell_GammasWide_Effi_data, 34, 2, colGammaWide, colGammaWide);
  DrawSetMarker(hNCell_GammasWide_Effi_MC, 8, 3, colGammaWide, colGammaWide); //hNCell_GammasWide_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasWide_Effi_Ratio, 34, 2, colGammaWide, colGammaWide);
  DrawSetMarker(hNCell_GammasWide_Effi_Corr, 34, 2, colGammaWide, colGammaWide);

  int colGammaHigh = kYellow  - 2;
  DrawSetMarker(hNCell_GammasHigh_Effi_data, 21, 2, colGammaHigh, colGammaHigh);
  DrawSetMarker(hNCell_GammasHigh_Effi_MC, 9, 3, colGammaHigh, colGammaHigh); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasHigh_Effi_Ratio, 21, 2, colGammaHigh, colGammaHigh);
  DrawSetMarker(hNCell_GammasHigh_Effi_Corr, 21, 2, colGammaHigh, colGammaHigh);

  int colGammaRWHigh = kPink  - 2;
  DrawSetMarker(hNCell_GammasHigh_RW_Effi_data, 21, 2, colGammaRWHigh, colGammaRWHigh);
  DrawSetMarker(hNCell_GammasHigh_RW_Effi_MC, 9, 3, colGammaRWHigh + 4, colGammaRWHigh + 4); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasHigh_RW_Effi_Ratio, 21, 2, colGammaRWHigh, colGammaRWHigh);
  DrawSetMarker(hNCell_GammasHigh_RW_Effi_Corr, 21, 2, colGammaRWHigh, colGammaRWHigh);


  int colGammaRWHighSBSub = kGreen + 2;
  DrawSetMarker(hNCell_GammasHigh_RW_SBSub_Effi_data, 21, 2, colGammaRWHighSBSub, colGammaRWHighSBSub);
  DrawSetMarker(hNCell_GammasHigh_RW_SBSub_Effi_MC, 9, 3, colGammaRWHighSBSub, colGammaRWHighSBSub); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasHigh_RW_SBSub_Effi_Ratio, 21, 2, colGammaRWHighSBSub, colGammaRWHighSBSub);
  DrawSetMarker(hNCell_GammasHigh_RW_SBSub_Effi_Corr, 21, 2, colGammaRWHighSBSub, colGammaRWHighSBSub);


  int colGammaRWHighMCSBSub = kOrange + 2;
  DrawSetMarker(hNCell_GammasHigh_RW_MCSBSub_Effi_data, 21, 2, colGammaRWHighMCSBSub, colGammaRWHighMCSBSub);
  DrawSetMarker(hNCell_GammasHigh_RW_MCSBSub_Effi_MC, 9, 3, colGammaRWHighMCSBSub, colGammaRWHighMCSBSub); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_GammasHigh_RW_MCSBSub_Effi_Ratio, 21, 2, colGammaRWHighMCSBSub, colGammaRWHighMCSBSub);
  DrawSetMarker(hNCell_GammasHigh_RW_MCSBSub_Effi_Corr, 21, 2, colGammaRWHighMCSBSub, colGammaRWHighMCSBSub);


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
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Effi_AllClusAndTB.pdf", fPeriod.Data()));

}


void Effi::PlotEffi_AllGammaAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Effi_GammaAndTB.pdf", fPeriod.Data()));
}


void Effi::PlotEffi_RWGammaAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Effi_GammaAndTB_Reweighted.pdf", fPeriod.Data()));
}


void Effi::PlotEffi_SBGammaAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Effi_GammaAndTB_Sideband.pdf", fPeriod.Data()));
}


void Effi::PlotEffi_RWGammaHighAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Effi_GammaHighAndTB_Reweighted.pdf", fPeriod.Data()));
}

void Effi::PlotEffi_RWGammaHighAndTB_MCCheck(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Effi_GammaHighAndTB_Reweighted_MCConsistency.pdf", fPeriod.Data()));
}


void Effi::PlotEffi_TrueAndReweighted(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Effi_GammaAndTB_Reweighted_consistency.pdf", fPeriod.Data()));
}


void Effi::PlotEffi_HighAllMethods(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Effi_GammaAndTB_SBSub.pdf", fPeriod.Data()));
}

void Effi::PlotEffi_HighAllMethodsMCSBSub(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_GammasHigh_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_Effi_data->Draw("same,p");
  hNCell_GammasHigh_RW_MCSBSub_Effi_data->Draw("same,p");
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
  drawLatexAdd(" M_{#pi^{0}}< M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" only high #gamma selected!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/Effi_GammaAndTB_MCSidebandSub.pdf", fPeriod.Data()));
}

void Effi::PlotRatio_GammasAndRW(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Ratio_GammaAndTB_Reweighted.pdf", fPeriod.Data()));
}

void Effi::PlotRatio_SBGammasAndRW(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Ratio_SBGammaAndTB.pdf", fPeriod.Data()));
}


void Effi::PlotRatio_GammasAndRWLeft(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Ratio_GammaAndTB_Reweighted_Left.pdf", fPeriod.Data()));
}


void Effi::PlotRatio_GammasAndRWWide(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_GammasWide_Effi_Ratio->Draw("same,p");
  hNCell_GammasWide_RW_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.65, 0.66, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasWide_Effi_Ratio, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_RW_Effi_Ratio, "P2, data, RW #gamma", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05 < M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/Ratio_GammaAndTB_Reweighted_Wide.pdf", fPeriod.Data()));
}

void Effi::PlotRatio_GammasAndRWHigh(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Ratio_GammaAndTB_Reweighted_High.pdf", fPeriod.Data()));
}

void Effi::PlotRatio_GammasAndRWHighMCSBSub(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Ratio_GammaAndTB_Reweighted_High_MCSBSub.pdf", fPeriod.Data()));
}


void Effi::PlotCorr_GammasAndRW(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Corr_GammaAndTB_Reweighted.pdf", fPeriod.Data()));
}

void Effi::PlotCorr_GammasAndRW_Left(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Corr_GammaLeftAndTB_Reweighted.pdf", fPeriod.Data()));
}


void Effi::PlotCorr_GammasAndRW_Wide(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_AllClus_Effi_Corr->Draw("same,p");
  hNCell_GammasWide_Effi_Corr->Draw("same,p");
  hNCell_GammasWide_RW_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.5, 0.35, 0.7, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Corr, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_GammasWide_Effi_Corr, "P2, data, #gamma", "p");
  leg->AddEntry(hNCell_GammasWide_RW_Effi_Corr, "P2, data, RW #gamma", "p");
  leg->AddEntry(grTB_Corr, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(" M_{#pi^{0}} - 0.05< M_{#gamma#gamma} < M_{#pi^{0}} + 0.02 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s/Corr_GammaWideAndTB_Reweighted.pdf", fPeriod.Data()));
}


void Effi::PlotCorr_GammasAndRW_High(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Corr_GammaHighAndTB_Reweighted.pdf", fPeriod.Data()));
}


void Effi::PlotCorr_HighAllMethods(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, 5., 1,1, 2, kGray+2, 2);
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
  hDummyCan->SaveAs(Form("%s/Corr_GammaHighAndTB_AllMethods.pdf", fPeriod.Data()));
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

  hDummyCan->SaveAs(Form("%s/Purity.pdf", fPeriod.Data()));

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

  hDummyCan->SaveAs(Form("%s/Purity2.pdf", fPeriod.Data()));

}



void Effi::PlotExampleBin(){

  hDummyCan->cd();
  hDummyMInv->GetYaxis()->SetRangeUser(0., hInvMassVsPt_Slice->GetMaximum()*1.3);
  hDummyMInv->Draw();

  hInvMassVsPt_Slice->Draw("same,p");
  hInvMassVsPtGG_Slice->Draw("same");
  hInvMassVsPtGC_Slice->Draw("same");
  hInvMassVsPtCC_Slice->Draw("same");

  TH1F* hMCTrueSum = (TH1F*) hInvMassVsPtGG_Slice->Clone("hMCTrueSum");
  hMCTrueSum->Add(hInvMassVsPtGC_Slice);
  hMCTrueSum->Add(hInvMassVsPtCC_Slice);
  DrawSetMarker(hMCTrueSum, 2, 3.4, kBlue + 2, kBlue + 2);

  hMCTrueSum->Draw("same");

  TFile fmasspos("fMassPos.root");
  TString name = "MassPos8TeV_func";
  if(fPeriod.Contains("13TeV")) name = "MassPos13TeV_func";
  TF1* func = (TF1*) fmasspos.Get(name);


  TF1 *fgaus = new TF1("fgaus", "gaus(0)", 0.1, 0.15);
  fgaus->SetParameters(100, 0.13, 0.05);
  hInvMassVsPt_Slice->Fit(fgaus, "MQR0");

  float mean = func->Eval(2.4);

  cout<<"peak pos from global fit: "<<mean<<"   from gaus: "<<fgaus->GetParameter(1)<<endl;

  TLegend *leg = GetAndSetLegend2(0.7, 0.8, 0.95, 0.95, 40);
  leg->AddEntry(hInvMassVsPt_Slice, "MC, all cluster pairs", "l");
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


  hDummyCan->SaveAs(Form("%s/ExampleBin.pdf", fPeriod.Data()));



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

  hDummyCan->SaveAs(Form("%s/ExampleBinWithData.pdf", fPeriod.Data()));

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
  if(fPeriod.Contains("13TeV")) name = "MassPos13TeV_func";
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


  hDummyCan->SaveAs(Form("%s/ExampleBinHigh.pdf", fPeriod.Data()));



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


  hDummyCan->SaveAs(Form("%s/ExampleBinHighCells.pdf", fPeriod.Data()));


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


  hDummyCan->SaveAs(Form("%s/ExampleBinHighCellsMC.pdf", fPeriod.Data()));


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

  hDummyCan->SaveAs(Form("%s/RatioNCells.pdf", fPeriod.Data()));

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

  hDummyCan->SaveAs(Form("%s/NCell1Fraction.pdf", fPeriod.Data()));

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
  hDummyCan->SetLogy();
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

  hDummyCan->SaveAs(Form("%s/NCellDistribution.pdf", fPeriod.Data()));

}




void PlottingNCellEffi(){

  {
    gSystem->Exec("mkdir 13TeV");
    Effi plot("Pi0Tagging_13TeV_02_08.root", "13TeV");
    plot.FillHistos();
    plot.SetOtherHistos();
    plot.DoSBSubtraction();
    plot.GetHistReweighted();
    plot.FillCorrHistos();
    plot.LoadTB();
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
    plot.PlotRatio_GammasAndRWHigh();
    plot.PlotRatio_GammasAndRWHighMCSBSub();
    plot.PlotRatio_SBGammasAndRW();
    plot.PlotEffi_HighAllMethods();
    plot.PlotEffi_HighAllMethodsMCSBSub();
    plot.PlotCorr_GammasAndRW();
    plot.PlotCorr_GammasAndRW_Left();
    plot.PlotCorr_GammasAndRW_Wide();
    plot.PlotCorr_GammasAndRW_High();
    plot.PlotCorr_HighAllMethods();
    plot.PlotPurity();
    plot.PlotPurity2();
    plot.PlotExampleBin();
    plot.PlotExampleBinCells();
    plot.PlotNCellRatios();
    plot.PlotFraction1cell();
    plot.PlotNCellDistr();
  }
  // {
  //   gSystem->Exec("mkdir 8TeV");
  //   Effi plot("Pi0Tagging_8TeV_02_04.root", "8TeV");
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
