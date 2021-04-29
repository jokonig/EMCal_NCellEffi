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
  void SetAddition(TString tmp)       {fMethod = tmp;};
  void SetSpecialName(TString tmp)       {fSpecialName = tmp; };
  float calcErr(float Emc, float Errmc, float Edata, float Errdata);
  void GetEffiHists(TH2F *hdata2d, TH2F *hMC2d, TH1D *&hEffiData,  TH1D *&hEffiMC, TH1D *&hRatio, TH1D *&hCorr, bool doSmooth = false);
  void GetTrueHists(TH2F* h2, TH1D *&h, TString name = "h");
  void FillCorrHistos();
  void LoadTB();
  void GetHistReweighted();
  void SetOtherHistos();
  void DoSBSubtraction();
  void SetPurityHistosAfterSBSub();

  TH2F* SubtractSidebandBack(TH2F* hMInv, TH2F* hMInvBack, TH2F* hNCellSB, TH2F* hNCell);

  TH2F* SubtractHadronsFromMC(TH2F* hMInv, TH2F* hMInvBack, TH2F* hNCellSB, TH2F* hNCell);


  void SetPlotting();

  // Effi Plots
  void PlotEffi_TB();    // Effi Plot with all clusters from P2
  void PlotEffi_AllClusAndTB();    // Effi Plot with all clusters from P2
  void PlotEffi_AllClusAndTBAndTrueGamma();    // Effi Plot with all clusters from P2 and true gammas
  void PlotEffi_Wide();            // Effi Plot with clusters in pi0 mass range
  void PlotEffi_WideWithTrue();            // Effi Plot with clusters in pi0 mass range with true gammas
  void PlotEffi_WideWithRW();            // Effi Plot with clusters in pi0 mass range with true gammas
  void PlotEffi_WideWithRWSBSub();            // Effi Plot with clusters in pi0 mass range with true gammas
  void PlotEffi_TrueVsRecE();            // Effi Plot with clusters in pi0 mass range with true gammas


  // Ratio Plots
  void PlotRatio_Wide();   // Plot ratio with TB + wide clusters with and without RW and SB sub
  void PlotRatio_TBAndAll();   // Plot ratio with TB + all clus
  void PlotRatio_ConvMod();   // Plot ratio with TB + wide clusters with and without RW and SB sub

  // Corr Plots
  void PlotCorr_Wide();   // Plot correction factor with TB + all clusters + gammas (right) + gammas (right) reweighted
  void PlotCorr_TBAndAll();   // Plot correction factor with TB + all clusters

  // Fitting of CorrPlots
  void FitCorr_GammasWide();      // Fitting of correction with

  void PlotMCClosure();      // Fitting of correction with
  void PlotDataClosure();      // Fitting of correction with

  void PlotTrueVsRecE();

  // Purity
  void PlotPurity();
  // ExampleBin
  void PlotExampleBin();

  void PlotNCellVsE();

  void PlotEffiSources();

  void ScaleTo(TH1D *h, TH1D *& h2, float down = 0.2, float up = 0.3);

  void WriteToFile();

  TFile* fcontrolOut = nullptr;

private:

  TFile* fdata = nullptr;

  TFile *fMassPos = nullptr;

  TString fPeriod = "";

  TString fMethod = "";

  TString fSpecialName = "";

  TString fsuffix = "png";

  float PlotenergyHigh = 8.;

  TString sEnergy = "pp #sqrt{s} = 13 TeV";

  TString StrSelectedGamma(){
    if(fMethod.Contains("Low")) return "lower clus. selected";
    if(fMethod.Contains("High")) return "higher clus. selected";
    else return "both clus. selected";
  }
  TString StrSelectedRange(){
    if(fMethod.Contains("Wide")) return "0.09 < M_{inv} < 0.17 GeV/#it{c}^{2}";
    if(fMethod.Contains("Left")) return "M_{#pi^{0}} - 0.05 < M_{inv} < M_{#pi^{0}} GeV/#it{c}^{2}";
    if(fMethod.Contains("Right")) return "M_{#pi^{0}} < M_{inv} < M_{#pi^{0}} + 0.02 GeV/#it{c}^{2}";
    else return "0.09 < M_{inv} < 0.17 GeV/#it{c}^{2}";
  }
  std::array<double, 2> getRangeSignal(float pT = 0){
    if(!fMassPos) fMassPos = TFile::Open("fMassPos.root");
    TF1 *func = (TF1*) fMassPos->Get(Form("MassPos%s_func", fPeriod.Data()));
    if(fMethod.Contains("Right")){
      return {func->Eval(pT), func->Eval(pT) + 0.02};
    }
    else if(fMethod.Contains("Left")){
      return {func->Eval(pT) - 0.05, func->Eval(pT)};
    }
    return {0.09, 0.17};
  }

  // Test Beam
  TGraphErrors *grTB_data = nullptr;
  TGraphErrors *grTB_MC = nullptr;
  TGraphErrors *grTB_Ratio = nullptr;
  TGraphErrors *grTB_Corr = nullptr;

  // All clusters
  TH2F* hNCellVsETMNL_data = nullptr;
  TH2F* hNCellVsETMNL_MC = nullptr;

  TH2F* hNCellVsEGammasNLTrue_TrueE_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrue_RecE_MC = nullptr;

  // Gammas Wide
  TH2F* hNCellVsEGammasNL_data = nullptr;
  TH2F* hNCellVsEGammasNL_MC = nullptr;
  TH2F* hNCellVsETrueGammasNL_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElec_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueHadrons_MC = nullptr;

  // Gammas Wide sideband subtracted
  TH2F* hNCellVsEGammasNL_SBSub_data = nullptr;
  TH2F* hNCellVsEGammasNL_SBSub_MC = nullptr;
  TH2F* hNCellVsETrueGammasNL_SBSub_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElec_SBSub_MC = nullptr;

  // MC sideband subtracted
  TH2F* hNCellVsEGammasNL_MCSBSub_data = nullptr;
  TH2F* hNCellVsEGammasNL_MCSBSub_MC = nullptr;
  // MC sideband subtracted and reweighted
  TH2F* hNCellVsEGammasNL_MCSBSub_RW_data = nullptr;
  TH2F* hNCellVsEGammasNL_MCSBSub_RW_MC = nullptr;

  // Gammas Wide reweighted
  TH2F* hNCellVsEGammasNL_RW_data = nullptr;
  TH2F* hNCellVsEGammasNL_RW_MC = nullptr;

  // Gammas Wide reweighted but 10% different conversion probability
  TH2F* hNCellVsEGammasNL_RW_ConvMod10_data = nullptr;

  // Gammas Wide RW + sideband subtracted
  TH2F* hNCellVsEGammasNL_RW_SBSub_data = nullptr;
  TH2F* hNCellVsEGammasNL_RW_MCSBSub_data = nullptr;
  TH2F* hNCellVsEGammasNL_RW_SBSub_MC = nullptr;

  // Gammas SideBand
  TH2F* hNCellVsEGammasNLSB_data = nullptr;
  TH2F* hNCellVsEGammasNLSB_MC = nullptr;
  TH2F* hNCellVsETrueGammasNLSB_MC = nullptr;
  TH2F* hNCellVsEGammasNLTrueElecSB_MC = nullptr;

  // conversions only
  TH2F* hNCellVsEConversionsNL_RW_data = nullptr;
  TH2F* hNCellVsEConversionsNL_RW_MC = nullptr;


  // TH2F* hNCellVsTrueEelecNL_data = nullptr;


  // all clusters
  TH1D* hNCell_AllClus_Effi_data = nullptr;
  TH1D* hNCell_AllClus_Effi_MC = nullptr;
  TH1D* hNCell_AllClus_Effi_Ratio = nullptr;
  TH1D* hNCell_AllClus_Effi_Corr = nullptr;


  TH1D* hNCell_TrueGammas_TrueE_Effi_data = nullptr;
  TH1D* hNCell_TrueGammas_TrueE_Effi_MC = nullptr;
  TH1D* hNCell_TrueGammas_TrueE_Effi_Ratio = nullptr;
  TH1D* hNCell_TrueGammas_TrueE_Effi_Corr = nullptr;

  TH1D* hNCell_TrueGammas_RecE_Effi_data = nullptr;
  TH1D* hNCell_TrueGammas_RecE_Effi_MC = nullptr;
  TH1D* hNCell_TrueGammas_RecE_Effi_Ratio = nullptr;
  TH1D* hNCell_TrueGammas_RecE_Effi_Corr = nullptr;


  // gammas wide
  TH1D* hNCell_Gammas_Effi_data = nullptr;
  TH1D* hNCell_Gammas_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_Effi_Corr = nullptr;

  // gammas wide
  TH1D* hNCell_Gammas_Smoothed_Effi_data = nullptr;
  TH1D* hNCell_Gammas_Smoothed_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_Smoothed_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_Smoothed_Effi_Corr = nullptr;

  // True gammas wide
  TH1D* hNCell_TrueGammas_data = nullptr;
  TH1D* hNCell_TrueGammas_MC = nullptr;
  TH1D* hNCell_TrueGammas_Ratio = nullptr;
  TH1D* hNCell_TrueGammas_Corr = nullptr;

  // True gammas Sideband
  TH1D* hNCell_TrueGammasSB_data = nullptr;
  TH1D* hNCell_TrueGammasSB_MC = nullptr;
  TH1D* hNCell_TrueGammasSB_Ratio = nullptr;
  TH1D* hNCell_TrueGammasSB_Corr = nullptr;

  // gammas SideBand
  TH1D* hNCell_GammasSB_Effi_data = nullptr;
  TH1D* hNCell_GammasSB_Effi_MC = nullptr;
  TH1D* hNCell_GammasSB_Effi_Ratio = nullptr;
  TH1D* hNCell_GammasSB_Effi_Corr = nullptr;


  // gammas reweighted Pure MC
  TH1D* hNCell_Gammas_RW_Effi_MC_consistency = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_MC_consistency2 = nullptr;
  // TH1D* hNCell_Gammas_RW_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_Ratio_consistency = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_Corr_consistency = nullptr;

  // gammas reweighted wide
  TH1D* hNCell_Gammas_RW_Effi_data = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_RW_Effi_Corr = nullptr;

  // gammas SB sub wide
  TH1D* hNCell_Gammas_SBSub_Effi_data = nullptr;
  TH1D* hNCell_Gammas_SBSub_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_SBSub_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_SBSub_Effi_Corr = nullptr;

  // gammas reweighted SB sub wide
  TH1D* hNCell_Gammas_RW_SBSub_Effi_data = nullptr;
  TH1D* hNCell_Gammas_RW_SBSub_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_RW_SBSub_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_RW_SBSub_Effi_Corr = nullptr;

  // gammas reweighted MC hadrons sub wide
  TH1D* hNCell_Gammas_RW_MCHadSub_Effi_data = nullptr;
  TH1D* hNCell_Gammas_RW_MCHadSub_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_RW_MCHadSub_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_RW_MCHadSub_Effi_Corr = nullptr;

  // conversions reweighted wide
  TH1D* hNCell_Conversions_RW_Effi_data = nullptr;
  TH1D* hNCell_Conversions_RW_Effi_MC = nullptr;
  TH1D* hNCell_Conversions_RW_Effi_Ratio = nullptr;
  TH1D* hNCell_Conversions_RW_Effi_Corr = nullptr;

  // gammas reweighted wide 10% diff conv. probability
  TH1D* hNCell_Gammas_RW_ConvMod10_Effi_data = nullptr;
  TH1D* hNCell_Gammas_RW_ConvMod10_Effi_MC = nullptr;
  TH1D* hNCell_Gammas_RW_ConvMod10_Effi_Ratio = nullptr;
  TH1D* hNCell_Gammas_RW_ConvMod10_Effi_Corr = nullptr;



  TH2F* hNCellVsEAllClusTrueGamma = nullptr;
  TH2F* hNCellVsEAllClusTrueElec = nullptr;
  TH2F* hNCellVsEAllClusTrueHadrons = nullptr;


  // InvMass
  TH2F* hInvMassVsPt_MC = nullptr;
  TH2F* hInvMassVsPtGG = nullptr;
  TH2F* hInvMassVsPtGC = nullptr;
  TH2F* hInvMassVsPtCC = nullptr;
  TH2F* hInvMassVsPt_data = nullptr;


  TH1D* hInvMassVsPt_Slice = nullptr;
  TH1D* hInvMassVsPtGG_Slice = nullptr;
  TH1D* hInvMassVsPtGC_Slice = nullptr;
  TH1D* hInvMassVsPtCC_Slice = nullptr;
  TH1D* hInvMassBackVsPt_Slice = nullptr;
  TH1D* hInvMassVsPt_data_Slice = nullptr;


  // for those histos, low is loaded and high added to result in wide
  TH2F* hInvMassVsClusterPt_MC = nullptr;
  TH2F* hInvMassVsClusterPt_data = nullptr;

  TH2F* hInvMassVsClusterPtGamma_MC = nullptr;
  TH2F* hInvMassVsPtGamma_MC = nullptr;
  TH2F* hInvMassVsPtElec_MC = nullptr;

  TH2F* hInvMassVsClusterPtBack_data = nullptr;
  TH2F* hInvMassVsClusterPtBack_MC = nullptr;


  TH1D* hGammaPurity = nullptr;
  TH1D* hGammaPuritySB = nullptr;
  TH1D* hGammaPuritySBSub = nullptr;
  TH1D* hElecPurity = nullptr;
  TH1D* hElecPuritySB = nullptr;
  TH1D* hElecPuritySBSub = nullptr;
  TH1D* hUnPurity = nullptr;
  TH1D* hUnPuritySB = nullptr;
  TH1D* hUnPuritySBSub = nullptr;

  TH1D* hTrueGammas = nullptr;
  TH1D* hTrueElec = nullptr;
  TH1D* hTrueGammasSB = nullptr;
  TH1D* hTrueElecSB = nullptr;


  // PlottingHists
  TH2F* hDummyEffi = nullptr;
  TH2F* hDummyRatio = nullptr;
  TH2F* hDummyCorr = nullptr;
  TH2F* hDummyPurity = nullptr;
  TH2F* hDummyMInv = nullptr;
  TH2F* hDummyMInvRatio = nullptr;
  TH2F* hDummyMassPos = nullptr;
  TH2F* hDummyMCClosure = nullptr;
  TH2F* hDummyNCellRatio = nullptr;
  TH2F* hDummyNCellFraction = nullptr;
  TH2F* hDummyERecVsETrue = nullptr;
  TH2F* hDummyERecVsETrue2 = nullptr;
  TH2F* hDummyFraction = nullptr;

  TH2F* hRecDivTrueE = nullptr;
  TH2F* hRecDivTrueEOneCell = nullptr;
  TH2F* hRecDivTrueETwoCell = nullptr;
  TH2F* hRecDivTrueEThreeCell = nullptr;

  TCanvas* hDummyCan = nullptr;
  TCanvas* hDummyCan2d = nullptr;

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

  TFile TB("Fig35right.root");
  grTB_data = (TGraphErrors*) TB.Get("gData");
  TGraph* tmpTBMC = (TGraph*) TB.Get("gG3");

  // double *xData = tmpTBData->GetX();
  // double *yData = tmpTBData->GetY();
  double *xMC = tmpTBMC->GetX();
  double *yMC = tmpTBMC->GetY();
  // Double_t ey[8] = {0.04, 0.03, 0.03, 0.02, 0.015, 0.007, 0.001, 0.001};
  Double_t ex[8] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
  Double_t eMC[8] = {0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001};


  // grTB_data = new TGraphErrors(8, xData, yData, ex, ey);
  grTB_MC = new TGraphErrors(8, xMC, yMC, eMC, eMC);

  for(int i = 0; i < grTB_MC->GetN(); ++i){
    grTB_MC->SetPointError(i, 0.001, eMC[i]);
  }



  grTB_Ratio = (TGraphErrors*) grTB_MC->Clone("grTBRatio");

  // SetPointError

  for(int i = 0; i < grTB_Ratio->GetN(); ++i){
    grTB_Ratio->SetPointY(i, grTB_data->Eval(grTB_MC->GetPointX(i)) / grTB_MC->GetPointY(i));
    grTB_Ratio->SetPointError(i, 0.001, grTB_data->GetErrorY(i)/grTB_MC->GetPointY(i));
  }

  grTB_Corr = (TGraphErrors*) grTB_MC->Clone("grTBEffi");
  for(int i = 0; i < grTB_data->GetN(); ++i){
    float CF = -1;
    float CFerr = 0.001;
    if(grTB_MC->GetPointY(i) < 1){
      CF = 1 - (( 1 - grTB_data->Eval(grTB_MC->GetPointX(i)))/(1 - grTB_MC->GetPointY(i)));
      CFerr = grTB_data->GetErrorY(i)/(1-grTB_MC->GetPointY(i));
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


  hNCellVsEGammasNL_data = (TH2F*) fdata->Get(Form("hNCellVsEGammasNL%s_data", fMethod.Data()));
  hNCellVsEGammasNL_data ->Sumw2();
  hNCellVsEGammasNL_data ->SetDirectory(0);
  hNCellVsEGammasNL_MC = (TH2F*) fdata->Get(Form("hNCellVsEGammasNL%s_MC", fMethod.Data()));
  hNCellVsEGammasNL_MC->Sumw2();
  hNCellVsEGammasNL_MC->SetDirectory(0);
  hNCellVsETrueGammasNL_MC = (TH2F*) fdata->Get(Form("hNCellVsETrueGammasNL%s_MC", fMethod.Data()));
  hNCellVsETrueGammasNL_MC ->Sumw2();
  hNCellVsETrueGammasNL_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElec_MC = (TH2F*) fdata->Get(Form("hNCellVsEGammasNLTrueElec%s_MC", fMethod.Data()));
  hNCellVsEGammasNLTrueElec_MC ->Sumw2();
  hNCellVsEGammasNLTrueElec_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueHadrons_MC = (TH2F*) hNCellVsEGammasNL_MC->Clone("hNCellVsEGammasNLTrueHadrons_MC");
  hNCellVsEGammasNLTrueHadrons_MC->Add(hNCellVsETrueGammasNL_MC, -1);
  hNCellVsEGammasNLTrueHadrons_MC->Add(hNCellVsEGammasNLTrueElec_MC, -1);
  hNCellVsEGammasNLTrueHadrons_MC->Sumw2();

  TString isHighLow = "";
  if(fMethod.Contains("High")) isHighLow = "High";
  else if(fMethod.Contains("Low")) isHighLow = "Low";
  hNCellVsEGammasNLSB_data = (TH2F*) fdata->Get(Form("hNCellVsEGammasNLSideBand%s_data", isHighLow.Data()));
  hNCellVsEGammasNLSB_data ->Sumw2();
  hNCellVsEGammasNLSB_data ->SetDirectory(0);
  hNCellVsEGammasNLSB_MC = (TH2F*) fdata->Get(Form("hNCellVsEGammasNLSideBand%s_MC", isHighLow.Data()));
  hNCellVsEGammasNLSB_MC ->Sumw2();
  hNCellVsEGammasNLSB_MC ->SetDirectory(0);
  hNCellVsETrueGammasNLSB_MC = (TH2F*) fdata->Get(Form("hNCellVsETrueGammasNLSideBand%s_MC", isHighLow.Data()));
  hNCellVsETrueGammasNLSB_MC ->Sumw2();
  hNCellVsETrueGammasNLSB_MC ->SetDirectory(0);
  hNCellVsEGammasNLTrueElecSB_MC = (TH2F*) fdata->Get(Form("hNCellVsEGammasNLTrueElecSideBand%s_MC", isHighLow.Data()));
  hNCellVsEGammasNLTrueElecSB_MC->Sumw2();
  hNCellVsEGammasNLTrueElecSB_MC->SetDirectory(0);

  hRecDivTrueE = (TH2F*) fdata->Get("hRecDivTrueE_MC");
  hRecDivTrueEOneCell = (TH2F*) fdata->Get("hRecDivTrueEOneCell_MC");
  hRecDivTrueETwoCell = (TH2F*) fdata->Get("hRecDivTrueETwoCell_MC");
  hRecDivTrueEThreeCell = (TH2F*) fdata->Get("hRecDivTrueEThreeCell_MC");



  // true gammas, electrons, hadrons
  hNCellVsEAllClusTrueGamma = (TH2F*) fdata->Get("hNCellVsETrueGamma_MC");
  hNCellVsEAllClusTrueElec = (TH2F*) fdata->Get("hNCellVsETrueElec_MC");
  hNCellVsEAllClusTrueHadrons = (TH2F*) fdata->Get("hNCellVsETrueHadr_MC");

}


// -------------------------------------------------
// Steering Efficiency histos calculation
// -------------------------------------------------
void Effi::FillCorrHistos(){
  // all clusters
  GetEffiHists(hNCellVsETMNL_data, hNCellVsETMNL_MC, hNCell_AllClus_Effi_data, hNCell_AllClus_Effi_MC, hNCell_AllClus_Effi_Ratio, hNCell_AllClus_Effi_Corr);

  // wide clusters
  GetEffiHists(hNCellVsEGammasNL_data, hNCellVsEGammasNL_MC, hNCell_Gammas_Effi_data, hNCell_Gammas_Effi_MC, hNCell_Gammas_Effi_Ratio, hNCell_Gammas_Effi_Corr);

  // wide reweighted
  GetEffiHists(hNCellVsEGammasNL_RW_data, hNCellVsEGammasNL_RW_MC, hNCell_Gammas_RW_Effi_data, hNCell_Gammas_RW_Effi_MC, hNCell_Gammas_RW_Effi_Ratio, hNCell_Gammas_RW_Effi_Corr);

  // wide SB subtracted
  GetEffiHists(hNCellVsEGammasNL_SBSub_data, hNCellVsEGammasNL_SBSub_MC, hNCell_Gammas_SBSub_Effi_data, hNCell_Gammas_SBSub_Effi_MC, hNCell_Gammas_SBSub_Effi_Ratio, hNCell_Gammas_SBSub_Effi_Corr);

  // wide SB subtracted + RW
  cout<<"here it comes \n";
  GetEffiHists(hNCellVsEGammasNL_RW_SBSub_data, hNCellVsEGammasNL_RW_SBSub_MC, hNCell_Gammas_RW_SBSub_Effi_data, hNCell_Gammas_RW_SBSub_Effi_MC, hNCell_Gammas_RW_SBSub_Effi_Ratio, hNCell_Gammas_RW_SBSub_Effi_Corr);

  // hadrons from MC subtracted and reweighted
  GetEffiHists(hNCellVsEGammasNL_MCSBSub_RW_data, hNCellVsEGammasNL_MCSBSub_RW_MC, hNCell_Gammas_RW_MCHadSub_Effi_data, hNCell_Gammas_RW_MCHadSub_Effi_MC, hNCell_Gammas_RW_MCHadSub_Effi_Ratio, hNCell_Gammas_RW_MCHadSub_Effi_Corr);

  // true photons
  GetEffiHists(hNCellVsETrueGammasNL_MC, hNCellVsETrueGammasNL_MC, hNCell_TrueGammas_data, hNCell_TrueGammas_MC, hNCell_TrueGammas_Ratio, hNCell_TrueGammas_Corr);

  // true photons Sideband
  GetEffiHists(hNCellVsETrueGammasNLSB_MC, hNCellVsETrueGammasNLSB_MC, hNCell_TrueGammasSB_data, hNCell_TrueGammasSB_MC, hNCell_TrueGammasSB_Ratio, hNCell_TrueGammasSB_Corr);


  // wide reweighted conv. mod 10%
  GetEffiHists(hNCellVsEGammasNL_RW_ConvMod10_data, hNCellVsEGammasNL_RW_MC, hNCell_Gammas_RW_ConvMod10_Effi_data, hNCell_Gammas_RW_ConvMod10_Effi_MC, hNCell_Gammas_RW_ConvMod10_Effi_Ratio, hNCell_Gammas_RW_ConvMod10_Effi_Corr);

  // true photons vs true energy
  GetEffiHists(hNCellVsEGammasNLTrue_TrueE_MC, hNCellVsEGammasNLTrue_TrueE_MC, hNCell_TrueGammas_TrueE_Effi_data, hNCell_TrueGammas_TrueE_Effi_MC, hNCell_TrueGammas_TrueE_Effi_Ratio, hNCell_TrueGammas_TrueE_Effi_Corr);
  // true photons vs rec energy
  GetEffiHists(hNCellVsEGammasNLTrue_RecE_MC, hNCellVsEGammasNLTrue_RecE_MC, hNCell_TrueGammas_RecE_Effi_data, hNCell_TrueGammas_RecE_Effi_MC, hNCell_TrueGammas_RecE_Effi_Ratio, hNCell_TrueGammas_RecE_Effi_Corr);
}


// -------------------------------------------------
// Steering for sideband subtraction
// -------------------------------------------------
void Effi::DoSBSubtraction(){

  // data wide gamma
  hNCellVsEGammasNL_SBSub_data = SubtractSidebandBack( hInvMassVsClusterPt_data, hInvMassVsClusterPtBack_data , hNCellVsEGammasNLSB_data , hNCellVsEGammasNL_data);
  // MC wide gamma
  hNCellVsEGammasNL_SBSub_MC = SubtractSidebandBack( hInvMassVsClusterPt_MC, hInvMassVsClusterPtBack_MC , hNCellVsEGammasNLSB_MC , hNCellVsEGammasNL_MC);
  // true electrons
  hNCellVsEGammasNLTrueElec_SBSub_MC = SubtractSidebandBack( hInvMassVsClusterPt_MC, hInvMassVsClusterPtBack_MC , hNCellVsEGammasNLTrueElecSB_MC , hNCellVsEGammasNLTrueElec_MC);
  // true photons
  hNCellVsETrueGammasNL_SBSub_MC = SubtractSidebandBack( hInvMassVsClusterPt_MC, hInvMassVsClusterPtBack_MC , hNCellVsETrueGammasNLSB_MC , hNCellVsETrueGammasNL_MC);

  // MC true hadron subtraction
  hNCellVsEGammasNL_MCSBSub_data = SubtractHadronsFromMC(hInvMassVsClusterPtBack_data, hInvMassVsClusterPtBack_MC, hNCellVsEGammasNLTrueHadrons_MC, hNCellVsEGammasNL_data);
  // MC true hadron subtraction on MC
  cout<<"MC only"<<endl;
  hNCellVsEGammasNL_MCSBSub_MC = SubtractHadronsFromMC(hInvMassVsClusterPtBack_MC, hInvMassVsClusterPtBack_MC, hNCellVsEGammasNLTrueHadrons_MC, hNCellVsEGammasNL_MC);

}

// -------------------------------------------------
// SUBTRACT SIDEBAND FROM SIGNAL REGION
// Make Sure the sideband and the peak use the same method or gamma selection!
// -------------------------------------------------
TH2F* Effi::SubtractSidebandBack(TH2F* hMInv_tmp, TH2F* hMInvBack_tmp, TH2F* hNCellSB_tmp, TH2F* hNCell_tmp){

  TH2F* hMInv = (TH2F*) hMInv_tmp->Clone("hMInv");
  TH2F* hMInvBack = (TH2F*) hMInvBack_tmp->Clone("hMInvBack");
  TH2F* hNCellSB = (TH2F*) hNCellSB_tmp->Clone("hNCellSB");
  TH2F* hNCell = (TH2F*) hNCell_tmp->Clone("hNCell");

  // STrategy:
  // Fill a histo with a pT dependent scaling factors
  // this only works if the inv mass histo is only filled with the pT of the higher photon!


  TH2F* hScale = (TH2F*) hNCell->Clone("hScale");
  hScale->Reset();
  // loop over energy bins
  for(int i = 1; i <= hNCell->GetNbinsX(); ++i){
    float xValLow = hNCell->GetXaxis()->GetBinLowEdge(i) + 0.001;
    float xValUp = hNCell->GetXaxis()->GetBinUpEdge(i) - 0.001;

    // get inv. mass distribution for signal and background
    // scale background to same evt.
    TH1D* hMInv_Proj = (TH1D*) hMInv->ProjectionX(Form("hMInv_Proj%i", i), hMInv->GetYaxis()->FindBin(xValLow), hMInv->GetYaxis()->FindBin(xValUp));
    TH1D* hMInvBack_Proj = (TH1D*) hMInvBack->ProjectionX(Form("hMInvBack_Proj%i", i), hMInvBack->GetYaxis()->FindBin(xValLow), hMInvBack->GetYaxis()->FindBin(xValUp));
    float scaleFac = hMInv_Proj->Integral(hMInv_Proj->FindBin(0.2), hMInv_Proj->FindBin(0.3)) / hMInvBack_Proj->Integral(hMInvBack_Proj->FindBin(0.2), hMInvBack_Proj->FindBin(0.3));
    hMInvBack_Proj->Scale(scaleFac);

    // ranges for background and signal
    std::array<double, 2> rangeSB = {0.2, 0.3};
    std::array<double, 2> rangeSignal = getRangeSignal((xValLow + xValUp)*0.5);


    // calculate the integral of the background in signal range to estimate the amount of background
    // Do the same for the same evt. in the sideband region
    // calculate the ratio between the background and the sideband to later scale the sideband contribution
    float backInt = hMInvBack_Proj->Integral(hMInvBack_Proj->FindBin(rangeSignal[0]), hMInvBack_Proj->FindBin(rangeSignal[1]));
    float SBInt = hMInv_Proj->Integral(hMInv_Proj->FindBin(rangeSB[0]), hMInv_Proj->FindBin(rangeSB[1]));
    // cout<<"hMInv_Proj: "<<hMInvBack_Proj->GetBinContent(20)<<endl;
    // cout<<"backInt: "<<backInt<<"      SBInt: "<<SBInt<<"   backInt/SBInt: "<<backInt/SBInt<<endl;
    // loop over ncell bins
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
  hScale->Write(Form("%s_hScale",hNCellSB->GetName()));
  hNCellSB->Write(Form("%s_control",hNCellSB->GetName()));

  // subtract SB NCell distribution from Signal NCell distribution
  // This will also remove some photons!
  TH2F* hSBSubtracted = (TH2F*) hNCell->Clone("hSBSubtracted");
  hSBSubtracted->SetDirectory(0);
  hSBSubtracted->Add(hNCellSB, -1);
  cout<<__LINE__<<endl;
  return hSBSubtracted;

}


// -------------------------------------------------
// SUBTRACT HADRON CONTRIBUTION WITH MC INFO FROM SIGNAL REGION
// Make Sure the sideband and the peak use the same method or gamma selection!
// -------------------------------------------------
TH2F* Effi::SubtractHadronsFromMC(TH2F* hMInvBack_tmp, TH2F* hMInvBackMC_tmp, TH2F* hNCell_TrueHadr, TH2F* hNCell_tmp){


  TH2F* hNCell_TrueHadr_tmp = (TH2F*) hNCell_TrueHadr->Clone("hNCell_TrueHadr_tmp");
  TH2F* hScale = (TH2F*) hNCell_TrueHadr->Clone("hScale");
  hScale->Reset();
  // loop over energy bins
  for(int i = 1; i <= hNCell_TrueHadr->GetNbinsX(); ++i){
    float xValLow = hNCell_TrueHadr->GetXaxis()->GetBinLowEdge(i) + 0.001;
    float xValUp = hNCell_TrueHadr->GetXaxis()->GetBinUpEdge(i) - 0.001;

    std::array<double, 2> rangeSignal = getRangeSignal((xValLow + xValUp)*0.5);
    // get inv mass background for data and MC
    // get a scale factor for the MC to fit the data
    TH1D* hMInvBack_Proj = (TH1D*) hMInvBack_tmp->ProjectionX(Form("hMInvBack_Proj%i", i), hMInvBack_tmp->GetYaxis()->FindBin(xValLow), hMInvBack_tmp->GetYaxis()->FindBin(xValUp));
    TH1D* hMInvBackMC_Proj = (TH1D*) hMInvBackMC_tmp->ProjectionX(Form("hMInvBackMC_Proj%i", i), hMInvBackMC_tmp->GetYaxis()->FindBin(xValLow), hMInvBackMC_tmp->GetYaxis()->FindBin(xValUp));
    float intBackMC = hMInvBackMC_Proj->Integral(hMInvBackMC_Proj->FindBin(rangeSignal[0]), hMInvBackMC_Proj->FindBin(rangeSignal[1]));
    float intBack = hMInvBack_Proj->Integral(hMInvBack_Proj->FindBin(rangeSignal[0]), hMInvBack_Proj->FindBin(rangeSignal[1]));
    cout<<"intBackMC: "<<intBackMC<<"  intBack: "<<intBack<<endl;
    float scaleFac =  (intBackMC != 0) ? intBack / intBackMC : 1;
    cout<<"scaleFac: "<<scaleFac<<endl;
    if(isnan(scaleFac)) scaleFac = 1;
    for(int b = 1; b < hScale->GetNbinsY(); ++b){
      hScale->SetBinContent(i, b, scaleFac);
    }
  }


   // scale mc true ????
   hNCell_TrueHadr_tmp->Multiply(hScale);

   TH2F* hSBSubtracted = (TH2F*) hNCell_tmp->Clone("hSBSubtracted");
   hSBSubtracted->SetDirectory(0);
   hSBSubtracted->Add(hNCell_TrueHadr_tmp, -1);
   return hSBSubtracted;
}

void Effi::SetOtherHistos(){

  // both sides of peak
  hGammaPurity = (TH1D*) hNCellVsETrueGammasNL_MC->ProjectionX("hGammaPurity");
  hTrueGammas = (TH1D*) hNCellVsETrueGammasNL_MC->ProjectionX("hTrueGammas");
  hTrueElec = (TH1D*) hNCellVsEGammasNLTrueElec_MC->ProjectionX("hTrueElec");
  hElecPurity = (TH1D*) hNCellVsEGammasNLTrueElec_MC->ProjectionX("hElecPurity");
  hGammaPurity->Divide(hTrueGammas, (TH1D*) hNCellVsEGammasNL_MC->ProjectionX("hGammasForPurity"), 1, 1, "B");
  hElecPurity->Divide(hTrueElec, (TH1D*) hNCellVsEGammasNL_MC->ProjectionX("hElecForPurity"), 1, 1, "B");

  hUnPurity = (TH1D*) hGammaPurity->Clone("hUnPurity");
  hUnPurity->Add(hElecPurity);
  for(int i = 1; i <= hUnPurity->GetNbinsX(); ++i){ hUnPurity->SetBinContent(i, 1- hUnPurity->GetBinContent(i));}


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



  TString isHighLow = "Both";
  if(fMethod.Contains("Low")) isHighLow = "Low";
  if(fMethod.Contains("High")) isHighLow = "High";

  hInvMassVsClusterPt_data = (TH2F*) fdata->Get(Form("hInvMassVsPt%s_data", isHighLow.Data()));
  hInvMassVsClusterPt_data ->Sumw2();
  hInvMassVsClusterPtBack_data = (TH2F*) fdata->Get(Form("hInvMassVsGammaPtBack%s_data", isHighLow.Data()));
  hInvMassVsClusterPtBack_data->Sumw2();

  hInvMassVsClusterPt_MC = (TH2F*) fdata->Get(Form("hInvMassVsPt%s_MC", isHighLow.Data()));
  hInvMassVsClusterPt_MC ->Sumw2();
  hInvMassVsClusterPtBack_MC = (TH2F*) fdata->Get(Form("hInvMassVsGammaPtBack%s_MC", isHighLow.Data()));
  hInvMassVsClusterPtBack_MC->Sumw2();

  // get MC information photons, electrons

  hInvMassVsPtGamma_MC = (TH2F*) fdata->Get(Form("hInvMassVsPtGamma%s_MC", isHighLow.Data()));
  hInvMassVsPtElec_MC = (TH2F*) fdata->Get(Form("hInvMassVsPtElec%s_MC", isHighLow.Data()));

  // all true gammas for true energy
  hNCellVsEGammasNLTrue_TrueE_MC = (TH2F*) fdata->Get(Form("hNCellVsEGammasNLTrue_TrueE_MC"));
  // all true gammas for rec energy
  hNCellVsEGammasNLTrue_RecE_MC = (TH2F*) fdata->Get(Form("hNCellVsEGammasNLTrue_RecE_MC"));


}

//---------------------------------------------
// Get purity histograms after SB subtraction
//---------------------------------------------
void Effi::SetPurityHistosAfterSBSub(){
  if(!hNCellVsETrueGammasNL_SBSub_MC){
    cout<<"Warning!! hNCellVsETrueGammasNL_SBSub_MC not there..."<<endl;
    cout<<"Exit SetPurityHistosAfterSBSub, set these histos first!"<<endl;
    return;
  }

  hGammaPuritySBSub = (TH1D*) hNCellVsETrueGammasNL_SBSub_MC->ProjectionX("hGammaPuritySBSub");
  hElecPuritySBSub = (TH1D*) hNCellVsEGammasNLTrueElec_SBSub_MC->ProjectionX("hElecPuritySBSub");
  TH1D* hTrueGammasSBSub = (TH1D*) hNCellVsETrueGammasNL_SBSub_MC->ProjectionX("hTrueGammasSBSub");
  TH1D* hTrueElecSBSub = (TH1D*) hNCellVsEGammasNLTrueElec_SBSub_MC->ProjectionX("hTrueElecSBSub");
  hGammaPuritySBSub->Divide(hTrueGammasSBSub, (TH1D*) hNCellVsEGammasNL_SBSub_MC->ProjectionX("hGammasForPuritySBSub"), 1, 1, "B");
  hElecPuritySBSub->Divide(hTrueElecSBSub, (TH1D*) hNCellVsEGammasNL_SBSub_MC->ProjectionX("hGammasForPuritySBSub"), 1, 1, "B");

  hUnPuritySBSub = (TH1D*) hGammaPuritySBSub->Clone("hUnPuritySBSub");
  hUnPuritySBSub->Add(hElecPuritySBSub);
  for(int i = 1; i <= hUnPuritySBSub->GetNbinsX(); ++i){ hUnPuritySBSub->SetBinContent(i, 1- hUnPuritySBSub->GetBinContent(i));}


}


// ------------------------------------------
// Get NCell efficiency hists from 2d NCellVs E distributions
// ------------------------------------------
void Effi::GetEffiHists(TH2F *hdata2d, TH2F *hMC2d, TH1D *&hEffiData,  TH1D *&hEffiMC, TH1D *&hRatio, TH1D *&hCorr, bool doSmooth){

  TH1D* hdataNCell1 = (TH1D*) hdata2d->ProjectionX("hdataNCell1", hdata2d->GetYaxis()->FindBin(1.5), hdata2d->GetYaxis()->FindBin(19));

  TH1D* hdataNCellOnly1 = (TH1D*) hdata2d->ProjectionX("hdataNCellOnly1", hdata2d->GetYaxis()->FindBin(1.5), hdata2d->GetYaxis()->FindBin(1.5));
  TH1D* hdataNCell2 = (TH1D*) hdata2d->ProjectionX("hdataNCell2", hdata2d->GetYaxis()->FindBin(2.5), hdata2d->GetYaxis()->FindBin(19));
  TH1D* hMCNCell1 = (TH1D*) hMC2d->ProjectionX("hMCNCell1", hMC2d->GetYaxis()->FindBin(1.5), hMC2d->GetYaxis()->FindBin(19));
  TH1D* hMCNCellOnly1 = (TH1D*) hMC2d->ProjectionX("hMCNCellOnly1", hMC2d->GetYaxis()->FindBin(1.5), hMC2d->GetYaxis()->FindBin(1.5));
  TH1D* hMCNCell2 = (TH1D*) hMC2d->ProjectionX("hMCNCell2", hMC2d->GetYaxis()->FindBin(2.5), hMC2d->GetYaxis()->FindBin(19));


  hEffiData = (TH1D*) hdataNCell2->Clone("hEffidata");
  hEffiData->Divide(hdataNCell2, hdataNCell1, 1, 1, "B");
  for(int i = 1; i <= hEffiData->GetNbinsX(); ++i){
    if(isnan(hEffiData->GetBinContent(i))) hEffiData->SetBinContent(i, 1);
  }


  hEffiMC = (TH1D*) hMCNCell2->Clone("hEffiMC");
  hEffiMC->Divide(hMCNCell2, hMCNCell1, 1, 1, "B");
  for(int i = 1; i <= hEffiMC->GetNbinsX(); ++i){
    if(isnan(hEffiMC->GetBinContent(i))) hEffiMC->SetBinContent(i, 1);
  }

  // smoothing on effi histos
  if(doSmooth){
    // logistisches wachstum
    // TF1 * func = new TF1("func", "[2]/(1+exp(-[0]*(x - [1]))) - [2] + 1", 0.7, 3);
    TF1 * func = new TF1("func", func1, 1., 8, 4);
    func->SetParameter(0, 2.8);
    func->SetParameter(1, 9.5);
    func->SetParameter(2, 9.9);
    func->SetParameter(3, 970);

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
  TH2F* hElecFracWide = (TH2F*) hNCellVsEGammasNLTrueElec_MC->Clone("hElecFracWide");
  hElecFracWide->Divide(hNCellVsEGammasNL_MC);
  hElecFracWide->Multiply(hNCellVsEGammasNL_data);
  hNCellVsEGammasNL_RW_data->Add(hElecFracWide, -1);


  hNCellVsEGammasNL_RW_MC = (TH2F*) hNCellVsEGammasNL_MC->Clone("hNCellVsEGammasNL_RW_MC");
  hNCellVsEGammasNL_RW_MC->Add(hNCellVsEGammasNLTrueElec_MC, -1);


  hNCellVsEGammasNL_RW_ConvMod10_data = (TH2F*) hNCellVsEGammasNL_data->Clone("hNCellVsEGammasNL_RW_ConvMod10_data");
  // get electron fraction
  TH2F* hElecFracWideConvMod10 = (TH2F*) hNCellVsEGammasNLTrueElec_MC->Clone("hElecFracWideConvMod10");
  hElecFracWideConvMod10->Divide(hNCellVsEGammasNL_MC);
  hElecFracWideConvMod10->Multiply(hNCellVsEGammasNL_data);
  hNCellVsEGammasNL_RW_ConvMod10_data->Add(hElecFracWide, -0.9);


  hNCellVsEGammasNL_RW_SBSub_data = (TH2F*) hNCellVsEGammasNL_SBSub_data->Clone("hNCellVsEGammasNL_RW_SBSub_data");
  // get electron fraction
  // was there from before
  TH2F* hElecFracWide2 = (TH2F*) hNCellVsEGammasNLTrueElec_SBSub_MC->Clone("hElecFracWide2");
  hElecFracWide2->Divide(hNCellVsEGammasNL_SBSub_MC);
  hElecFracWide2->Multiply(hNCellVsEGammasNL_SBSub_data);
  // hElecFracWide2->Scale(hNCellVsEGammasNL_SBSub_data->Integral()/hNCellVsEGammasNL_SBSub_MC->Integral());
  hNCellVsEGammasNL_RW_SBSub_data->Add(hElecFracWide2, -1);


  hNCellVsEGammasNL_RW_SBSub_MC = (TH2F*) hNCellVsEGammasNL_SBSub_MC->Clone("hNCellVsEGammasNL_RW_SBSub_MC");
  // get electron fraction
  TH2F* hElecFracWideMC = (TH2F*) hNCellVsEGammasNLTrueElec_SBSub_MC->Clone("hElecFracWideMC");
  // hElecFracWideMC->Divide(hNCellVsEGammasNL_SBSub_MC);
  // hElecFracWideMC->Multiply(hNCellVsEGammasNL_SBSub_MC);
  hNCellVsEGammasNL_RW_SBSub_MC->Add(hElecFracWideMC, -1);


// reweight to get the conversion contribution
  hNCellVsEConversionsNL_RW_data = (TH2F*) hNCellVsEGammasNL_data->Clone("hNCellVsEConversionsNL_RW_data");
  // get photon fraction
  TH2F* hGammaFracWide = (TH2F*) hNCellVsETrueGammasNL_MC->Clone("hGammaFracWide");
  hGammaFracWide->Divide(hNCellVsEGammasNL_MC);
  hGammaFracWide->Multiply(hNCellVsEGammasNL_data);
  hNCellVsEConversionsNL_RW_data->Add(hGammaFracWide, -1);


  // Reweight the MC hadron subtracted ones
  hNCellVsEGammasNL_MCSBSub_RW_data = (TH2F*) hNCellVsEGammasNL_MCSBSub_data->Clone("hNCellVsEGammasNL_MCSBSub_RW_data");
  TH2F* hElecFracWideMC3 = (TH2F*) hNCellVsEGammasNLTrueElec_MC->Clone("hElecFracWideMC3");
  hElecFracWideMC3->Divide(hNCellVsEGammasNL_MCSBSub_MC);
  hElecFracWideMC3->Multiply(hNCellVsEGammasNL_MCSBSub_data);
  hNCellVsEGammasNL_MCSBSub_RW_data->Add(hElecFracWideMC3, -1);

  hNCellVsEGammasNL_MCSBSub_RW_MC = (TH2F*) hNCellVsEGammasNL_MCSBSub_MC->Clone("hNCellVsEGammasNL_MCSBSub_RW_MC");
  TH2F* hElecFracWideMC4 = (TH2F*) hNCellVsEGammasNLTrueElec_MC->Clone("hElecFracWideMC4");
  hNCellVsEGammasNL_MCSBSub_RW_MC->Add(hElecFracWideMC4, -1);


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

  gSystem->Exec(Form("mkdir -p %s_%s/%s/%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data()));


  if(fPeriod.Contains("13TeV")) sEnergy = "pp #sqrt{s} = 13 TeV";
  if(fPeriod.Contains("8TeV")) sEnergy = "pp #sqrt{s} = 8 TeV";

  StyleSettingsPaper();


  hDummyEffi    = new TH2F("hDummyEffi","hDummyEffi",1000,0.5, PlotenergyHigh,1000,0.1, 1.4);
  SetStyleHistoTH2ForGraphs(hDummyEffi, "#it{E} (GeV)","#nu", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyRatio    = new TH2F("hDummyRatio","hDummyRatio",1000,0.5, PlotenergyHigh,1000,0.9, 1.4);
  SetStyleHistoTH2ForGraphs(hDummyRatio, "#it{E} (GeV)","#nu_{MC}/#nu_{data}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyCorr    = new TH2F("hDummyCorr","hDummyCorr",1000,0.5, PlotenergyHigh,1000,0., 1.4);
  SetStyleHistoTH2ForGraphs(hDummyCorr, "#it{E} (GeV)","#rho", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyPurity    = new TH2F("hDummyPurity","hDummyPurity",1000,0.5, PlotenergyHigh,1000,-0.05, 1.2);
  SetStyleHistoTH2ForGraphs(hDummyPurity, "#it{E} (GeV)","P", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyMInv    = new TH2F("hDummyMInv","hDummyMInv",1000,0.01, 0.4,1000,0., 50000);
  SetStyleHistoTH2ForGraphs(hDummyMInv, "#it{M}_{inv} (GeV/c^{2})","counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyMInvRatio    = new TH2F("hDummyMInvRatio","hDummyMInvRatio",1000,0.01, 0.4,1000,0., 20);
  SetStyleHistoTH2ForGraphs(hDummyMInvRatio, "#it{M}_{inv} (GeV/c^{2})","ratio to rot. back.", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyMassPos    = new TH2F("hDummyMassPos","hDummyMassPos",1000,0., 20,1000,0.8, 1.2);
  SetStyleHistoTH2ForGraphs(hDummyMassPos, "E_{clus} (GeV)" ,"#it{M}_{inv}/#it{M}_{#pi^{0}}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyERecVsETrue    = new TH2F("hDummyERecVsETrue","hDummyERecVsETrue",1000,0.8, 1.2,1000,0., 1.);
  SetStyleHistoTH2ForGraphs(hDummyERecVsETrue, "#it{E}_{rec.}/#it{E}_{true}","normalized counts", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyERecVsETrue2    = new TH2F("hDummyERecVsETrue2","hDummyERecVsETrue2",1000,0., 10,1000,0.9, 1.15);
  SetStyleHistoTH2ForGraphs(hDummyERecVsETrue2,"#it{E}_{rec} (GeV)", "#it{E}_{rec.}/#it{E}_{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyMCClosure    = new TH2F("hDummyMCClosure","hDummyMCClosure",1000,0.5, PlotenergyHigh,1000, 0.8, 1.4);
  SetStyleHistoTH2ForGraphs(hDummyMCClosure, "#it{E} (GeV)","MC_{rec.}/MC_{true}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyNCellRatio    = new TH2F("hDummyNCellRatio","hDummyNCellRatio",1000,0.7, PlotenergyHigh,1000,0.75, 1.45);
  SetStyleHistoTH2ForGraphs(hDummyNCellRatio, "#it{E}_{clus} (GeV)","N_{clus}^{data}/N_{clus}^{MC}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);

  hDummyNCellFraction    = new TH2F("hDummyNCellFraction","hDummyNCellFraction",1000,0.4, PlotenergyHigh,1000,0., 1.);
  SetStyleHistoTH2ForGraphs(hDummyNCellFraction, "#it{E}_{clus} (GeV)","N_{clus}^{Ncell = 1}/N_{clus}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1., 510, 510);

  hDummyFraction    = new TH2F("hDummyFraction","hDummyFraction",1000,0.4, PlotenergyHigh,1000,0., 1.);
  SetStyleHistoTH2ForGraphs(hDummyFraction, "#it{E}_{clus} (GeV)","N_{clus}^{specie}/N_{clus}^{all}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1., 510, 510);
  // gPad->SetLogx();
  hDummyCan = new TCanvas("Can", "", 1200, 1000);
  DrawPaperCanvasSettings(hDummyCan, 0.1, 0.003, 0.003, 0.1);
  hDummyCan2d = new TCanvas("hDummyCan2d", "", 1200, 1000);
  DrawPaperCanvasSettings(hDummyCan2d, 0.11, 0.11, 0.003, 0.11);


  DrawSetMarkerTGraphErr(grTB_data, 21, 1.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_MC, 25, 2.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_Ratio, 21, 1.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_Corr, 21, 1.5, kBlue + 2, kBlue + 2, 2);


  int colAllClus = kRed + 2;
  DrawSetMarker(hNCell_AllClus_Effi_data, 20, 2, colAllClus, colAllClus);
  DrawSetMarker(hNCell_AllClus_Effi_MC, 2, 3, colAllClus, colAllClus); //hNCell_Gammas_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_AllClus_Effi_Ratio, 20, 2, colAllClus, colAllClus);
  DrawSetMarker(hNCell_AllClus_Effi_Corr, 20, 2, colAllClus, colAllClus);

  //
  // // true gamma
  DrawSetMarker(hNCell_TrueGammas_MC, 7, 2.8, kBlack, kBlack);
  DrawSetMarker(hNCell_TrueGammasSB_MC, 7, 2.8, kRed + 7,  kRed + 7);


  DrawSetMarker(hNCell_TrueGammas_TrueE_Effi_MC, 7, 2.8, kOrange + 7, kOrange + 7);
  DrawSetMarker(hNCell_TrueGammas_RecE_Effi_MC, 7, 2.8, kRed + 2, kRed + 2);
  // // true elec
  // DrawSetMarker(hNCell_TrueElec_Effi_MC, 2, 2.4, kSpring + 4, kSpring + 4);

  //
  int colGammaWide = kGreen + 2;
  DrawSetMarker(hNCell_Gammas_Effi_data, 34, 2, colGammaWide, colGammaWide);
  DrawSetMarker(hNCell_Gammas_Effi_MC, 8, 3, colGammaWide, colGammaWide); //hNCell_Gammas_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_Gammas_Effi_Ratio, 34, 2, colGammaWide, colGammaWide);
  DrawSetMarker(hNCell_Gammas_Effi_Corr, 34, 2, colGammaWide, colGammaWide);

  int colGammaRWWideSBSub = kOrange + 2;
  DrawSetMarker(hNCell_Gammas_RW_SBSub_Effi_data, 28, 2, colGammaRWWideSBSub, colGammaRWWideSBSub);
  DrawSetMarker(hNCell_Gammas_RW_SBSub_Effi_MC, 8, 3, colGammaRWWideSBSub, colGammaRWWideSBSub); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_Gammas_RW_SBSub_Effi_Ratio, 28, 2, colGammaRWWideSBSub, colGammaRWWideSBSub);
  DrawSetMarker(hNCell_Gammas_RW_SBSub_Effi_Corr, 28, 2, colGammaRWWideSBSub, colGammaRWWideSBSub);

  int colGammaRWWideHadronsSub = kGray + 2;
  DrawSetMarker(hNCell_Gammas_RW_MCHadSub_Effi_data, 29, 2.5, colGammaRWWideHadronsSub, colGammaRWWideHadronsSub);
  DrawSetMarker(hNCell_Gammas_RW_MCHadSub_Effi_MC, 9, 3, colGammaRWWideHadronsSub, colGammaRWWideHadronsSub); //hNCell_GammasLeft_Effi_MC->SetFillColor(kGreen - 6);
  DrawSetMarker(hNCell_Gammas_RW_MCHadSub_Effi_Ratio, 29, 2.5, colGammaRWWideHadronsSub, colGammaRWWideHadronsSub);
  DrawSetMarker(hNCell_Gammas_RW_MCHadSub_Effi_Corr, 29, 2.5, colGammaRWWideHadronsSub, colGammaRWWideHadronsSub);


  int colGammaRWWide = kPink + 2;
  DrawSetMarker(hNCell_Gammas_RW_Effi_data, 33, 2.4, colGammaRWWide, colGammaRWWide);
  DrawSetMarker(hNCell_Gammas_RW_Effi_MC, 7, 3.4, colGammaRWWide, colGammaRWWide); //hNCell_Gammas_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_Gammas_RW_Effi_Ratio, 33, 2.4, colGammaRWWide, colGammaRWWide);
  DrawSetMarker(hNCell_Gammas_RW_Effi_Corr, 33, 2.4, colGammaRWWide, colGammaRWWide);


  int colGammaSBWide = kCyan + 2;
  DrawSetMarker(hNCell_Gammas_SBSub_Effi_data, 27, 2.4, colGammaSBWide, colGammaSBWide);
  DrawSetMarker(hNCell_Gammas_SBSub_Effi_MC, 7, 3.4, colGammaSBWide, colGammaSBWide); //hNCell_Gammas_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_Gammas_SBSub_Effi_Ratio, 27, 2.4, colGammaSBWide, colGammaSBWide);
  DrawSetMarker(hNCell_Gammas_SBSub_Effi_Corr, 27, 2.4, colGammaSBWide, colGammaSBWide);

  int colGammaRWWideConvMod = kYellow + 2;
  DrawSetMarker(hNCell_Gammas_RW_ConvMod10_Effi_data, 33, 2.4, colGammaRWWideConvMod, colGammaRWWideConvMod);
  DrawSetMarker(hNCell_Gammas_RW_ConvMod10_Effi_MC, 7, 3.4, colGammaRWWideConvMod, colGammaRWWideConvMod); //hNCell_Gammas_RW_Effi_MC->SetFillColor(colGammaRW);
  DrawSetMarker(hNCell_Gammas_RW_ConvMod10_Effi_Ratio, 33, 2.4, colGammaRWWideConvMod, colGammaRWWideConvMod);
  DrawSetMarker(hNCell_Gammas_RW_ConvMod10_Effi_Corr, 33, 2.4, colGammaRWWideConvMod, colGammaRWWideConvMod);
}




void Effi::PlotEffi_TB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");

  TLegend *leg = GetAndSetLegend2(0.5, 0.16, 0.95, 0.16 + 0.05*4, 40);
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd("EMCal test beam",0.7,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd("0.12 < M_{#gamma#gamma} < 0.14 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_TB.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}



void Effi::PlotEffi_AllClusAndTB(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_AllClus_Effi_data->Draw("same,p");
  hNCell_AllClus_Effi_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.5, 0.16, 0.95, 0.16 + 0.05*4, 40);
  leg->AddEntry(hNCell_AllClus_Effi_data, "P2, data, all clus.", "p");
  leg->AddEntry(hNCell_AllClus_Effi_MC, "P2, MC, all clus.", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd("0.12 < M_{#gamma#gamma} < 0.14 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_AllClusAndTB.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_AllClusAndTBAndTrueGamma(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_AllClus_Effi_data->Draw("same,p");
  hNCell_AllClus_Effi_MC->Draw("same,histc");
  hNCell_TrueGammas_MC->Draw("same,histc");
  // hNCell_TrueElec_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.5, 0.16, 0.95, 0.16 + 0.05*5, 40);
  leg->AddEntry(hNCell_AllClus_Effi_data, "P2, data, all clus.", "p");
  leg->AddEntry(hNCell_AllClus_Effi_MC, "P2, MC, all clus.", "l");
  leg->AddEntry(hNCell_TrueGammas_MC, "P2, MC, true #gamma", "l");
  // leg->AddEntry(hNCell_TrueElec_MC, "P2, MC, true e^{#pm}", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd("0.12 < M_{#gamma#gamma} < 0.14 (GeV/#it{c}^{2})",0.15,0.87,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_AllClusAndTBAndTrue.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}


void Effi::PlotEffi_Wide(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_Gammas_Effi_data->Draw("same,p");
  hNCell_Gammas_Effi_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.5, 0.16, 0.95, 0.16 + 0.05*4, 40);
  leg->AddEntry(hNCell_Gammas_Effi_data, "P2, data", "p");
  leg->AddEntry(hNCell_Gammas_Effi_MC, "P2, MC", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_Wide.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotMCClosure(){
  hDummyCan->cd();
  hDummyMCClosure->SetTitle("");
  hDummyMCClosure->GetYaxis()->SetTitle("MC_{rec.}/MC_{val.}");
  hDummyMCClosure->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);

  TH1F* hMCClosure1 = (TH1F*) hNCell_Gammas_RW_SBSub_Effi_MC->Clone("hMCClosure1");
  hMCClosure1->Divide(hNCell_TrueGammas_MC);
  DrawSetMarker(hMCClosure1, 28, 2, kOrange + 2, kOrange + 2);
  hMCClosure1->Draw("same,p");
  TH1F* hMCClosure2 = (TH1F*) hNCell_Gammas_RW_Effi_MC->Clone("hMCClosure2");
  hMCClosure2->Divide(hNCell_TrueGammas_MC);
  DrawSetMarker(hMCClosure2, 33, 2, kPink + 2, kPink + 2);
  hMCClosure2->Draw("same,p");
  TH1F* hMCClosure3 = (TH1F*) hNCell_Gammas_Effi_MC->Clone("hMCClosure3");
  hMCClosure3->Divide(hNCell_TrueGammas_MC);
  DrawSetMarker(hMCClosure3, 33, 2, kGreen + 2, kGreen + 2);
  hMCClosure3->Draw("same,p");
  TH1F* hMCClosure4 = (TH1F*) hNCell_Gammas_SBSub_Effi_MC->Clone("hMCClosure4");
  hMCClosure4->Divide(hNCell_TrueGammas_MC);
  DrawSetMarker(hMCClosure4, 25, 2, kGray + 2, kGray + 2);
  hMCClosure4->Draw("same,p");
  TH1F* hMCClosure5 = (TH1F*) hNCell_Gammas_RW_MCHadSub_Effi_MC->Clone("hMCClosure5");
  hMCClosure5->Divide(hNCell_TrueGammas_MC);
  DrawSetMarker(hMCClosure5, 34, 2, kGray + 2, kGray + 2);
  hMCClosure5->Draw("same,p");

  TGraphErrors* grTB_MC_ratio = (TGraphErrors*) grTB_MC->Clone("grTB_MC_ratio");
  for(int i = 0; i < hNCell_TrueGammas_MC->GetNbinsX(); ++i){
    cout<<hNCell_TrueGammas_MC->GetBinCenter(i)<<" -> "<<hNCell_TrueGammas_MC->GetBinContent(i)<<endl;
  }
  for(int i = 0; i < grTB_MC_ratio->GetN(); ++i){
    cout<<grTB_MC_ratio->GetPointX(i)<<" = "<<hNCell_TrueGammas_MC->Interpolate(grTB_MC_ratio->GetPointX(i))<<" --> "<<grTB_MC_ratio->GetPointY(i)<<endl;
    grTB_MC_ratio->SetPointY(i, grTB_MC->GetPointY(i) / hNCell_TrueGammas_MC->Interpolate(grTB_MC_ratio->GetPointX(i)));
    int bin = hNCell_TrueGammas_MC->FindBin(grTB_MC_ratio->GetPointX(i));
    float binerr = hNCell_TrueGammas_MC->GetBinError(bin);
    float binCont = hNCell_TrueGammas_MC->Interpolate(grTB_MC_ratio->GetPointX(i));
    grTB_MC_ratio->SetPointError(i, 0.001, binerr*grTB_MC->GetPointY(i) / (binCont*binCont)); //0.02 is systematic
    cout<<"err: "<<binerr*grTB_MC->GetPointY(i) / (binCont*binCont)<<endl;
  }
  grTB_MC_ratio->Draw("same,p");


  TLegend *leg = GetAndSetLegend2(0.4, 0.13, 0.95, 0.35, 40);
  leg->AddEntry(hMCClosure3, "P2, rec. = in M_{#pi^{0}}", "p");
  leg->AddEntry(hMCClosure2, "P2, rec. = in M_{#pi^{0}} + RW", "p");
  leg->AddEntry(hMCClosure4, "P2, rec. = in M_{#pi^{0}} + SB sub", "p");
  leg->AddEntry(hMCClosure1, "P2, rec. = in M_{#pi^{0}} + SB sub + RW", "p");
  leg->AddEntry(hMCClosure5, "P2, rec. = in M_{#pi^{0}} + h^{#pm} sub + RW", "p");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("MC closure",0.7,0.87,textSizeLabelsRel,kTRUE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.82,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.77,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/MC_Closure.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}


void Effi::PlotDataClosure(){
  hDummyCan->cd();
  hDummyMCClosure->SetTitle("");
  hDummyMCClosure->GetYaxis()->SetTitle("X / Fit to TB");
  hDummyMCClosure->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);

  TF1 *func = new TF1("func", func1, 0.5, 10, 4);
  func->SetParameter(0, 2.8);
  func->SetParameter(1, 9.5);
  func->SetParameter(2, 9.9);
  func->SetParameter(3, 970);

  grTB_data->Fit(func, "MR0");

  // build Ratio To Fit

  TGraphErrors* grTB_data_ratio = (TGraphErrors*) grTB_data->Clone("grTB_data_ratio");
  for(int i = 0; i < grTB_data_ratio->GetN(); ++i) grTB_data_ratio->SetPointY(i, grTB_data_ratio->GetPointY(i)/ func->Eval(grTB_data_ratio->GetPointX(i)));
  grTB_data_ratio->Draw("same,P");

  TH1F* hNCell_Gammas_Effi_data_ratio = (TH1F*) hNCell_Gammas_Effi_data->Clone("hNCell_Gammas_Effi_data_ratio");
  hNCell_Gammas_Effi_data_ratio->Divide(func, 1);
  hNCell_Gammas_Effi_data_ratio->Draw("same");
  TH1F* hNCell_Gammas_SBSub_Effi_data_ratio = (TH1F*) hNCell_Gammas_SBSub_Effi_data->Clone("hNCell_Gammas_SBSub_Effi_data_ratio");
  hNCell_Gammas_SBSub_Effi_data_ratio->Divide(func,1);
  hNCell_Gammas_SBSub_Effi_data_ratio->Draw("same");
  TH1F* hNCell_Gammas_RW_Effi_data_ratio = (TH1F*) hNCell_Gammas_RW_Effi_data->Clone("hNCell_Gammas_RW_Effi_data_ratio");
  hNCell_Gammas_RW_Effi_data_ratio->Divide(func, 1);
  hNCell_Gammas_RW_Effi_data_ratio->Draw("same");
  TH1F* hNCell_Gammas_RW_SBSub_Effi_data_ratio = (TH1F*) hNCell_Gammas_RW_SBSub_Effi_data->Clone("hNCell_Gammas_RW_SBSub_Effi_data_ratio");
  hNCell_Gammas_RW_SBSub_Effi_data_ratio->Divide(func, 1);
  hNCell_Gammas_RW_SBSub_Effi_data_ratio->Draw("same");
  TH1F* hNCell_Gammas_RW_MCHadSub_Effi_data_ratio = (TH1F*) hNCell_Gammas_RW_MCHadSub_Effi_data->Clone("hNCell_Gammas_RW_MCHadSub_Effi_data_ratio");
  hNCell_Gammas_RW_MCHadSub_Effi_data_ratio->Divide(func, 1);
  hNCell_Gammas_RW_MCHadSub_Effi_data_ratio->Draw("same");

  TLegend *leg = GetAndSetLegend2(0.4, 0.6, 0.95, 0.8, 40);
  leg->AddEntry(hNCell_Gammas_Effi_data_ratio, "P2, rec. = in M_{#pi^{0}}", "p");
  leg->AddEntry(hNCell_Gammas_SBSub_Effi_data_ratio, "P2, rec. = in M_{#pi^{0}} + RW", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_data_ratio, "P2, rec. = in M_{#pi^{0}} + SB sub", "p");
  leg->AddEntry(hNCell_Gammas_RW_SBSub_Effi_data_ratio, "P2, rec. = in M_{#pi^{0}} + SB sub + RW", "p");
  leg->AddEntry(hNCell_Gammas_RW_MCHadSub_Effi_data_ratio, "P2, rec. = in M_{#pi^{0}} + h^{#pm} sub + RW", "p");
  leg->AddEntry(grTB_data_ratio, "TB", "p");
  leg->Draw("same");

  // func->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("data closure",0.7,0.87,textSizeLabelsRel,kTRUE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Data_Closure.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_WideWithTrue(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_Gammas_Effi_data->Draw("same,p");
  hNCell_Gammas_Effi_MC->Draw("same,histc");
  hNCell_TrueGammas_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.5, 0.16, 0.95, 0.16 + 0.05*5, 40);
  leg->AddEntry(hNCell_Gammas_Effi_data, "P2, data", "p");
  leg->AddEntry(hNCell_Gammas_Effi_MC, "P2, MC", "l");
  leg->AddEntry(hNCell_TrueGammas_MC, "P2, MC, true #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_WithTrue.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_WideWithRW(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_Gammas_Effi_data->Draw("same,p");
  hNCell_Gammas_Effi_MC->Draw("same,histc");
  hNCell_Gammas_RW_Effi_data->Draw("same,p");
  hNCell_Gammas_RW_Effi_MC->Draw("same,histc");
  hNCell_TrueGammas_MC->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.5, 0.16, 0.95, 0.16 + 0.05*7, 40);
  leg->AddEntry(hNCell_Gammas_Effi_data, "P2, data", "p");
  leg->AddEntry(hNCell_Gammas_Effi_MC, "P2, MC", "l");
  leg->AddEntry(hNCell_Gammas_RW_Effi_data, "P2, data, RW", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_MC, "P2, MC, RW", "l");
  leg->AddEntry(hNCell_TrueGammas_MC, "P2, MC, true #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_RW.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotEffi_WideWithRWSBSub(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_Gammas_Effi_data->Draw("same,p");
  hNCell_Gammas_Effi_MC->Draw("same,histc");
  // hNCell_Gammas_RW_Effi_data->Draw("same,p");
  // hNCell_Gammas_RW_Effi_MC->Draw("same,histc");
  hNCell_Gammas_RW_SBSub_Effi_data->Draw("same,p");
  hNCell_Gammas_RW_SBSub_Effi_MC->Draw("same,histc");
  hNCell_TrueGammas_MC->Draw("same,histc");

  TGraph gr(hNCell_Gammas_RW_SBSub_Effi_data);
  gr.Print();

  TLegend *leg = GetAndSetLegend2(0.5, 0.16, 0.95, 0.16 + 0.05*7, 40);
  leg->AddEntry(hNCell_Gammas_Effi_data, "P2, data", "p");
  leg->AddEntry(hNCell_Gammas_Effi_MC, "P2, MC", "l");
  leg->AddEntry(hNCell_Gammas_RW_SBSub_Effi_data, "P2, data, RW+SBSub", "p");
  leg->AddEntry(hNCell_Gammas_RW_SBSub_Effi_MC, "P2, MC, RW+SBSub", "l");
  leg->AddEntry(hNCell_TrueGammas_MC, "P2, MC, true #gamma", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_RW_SB.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
  hDummyCan->SetLogy();
  hDummyCan->SetLogx();
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_RW_SB_Log.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
  hDummyCan->SetLogy(0);
  hDummyCan->SetLogx(0);

}


void Effi::PlotEffi_TrueVsRecE(){
  hDummyCan->cd();
  hDummyEffi->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_data->Draw("same,p");
  grTB_MC->Draw("same,l");
  hNCell_TrueGammas_TrueE_Effi_MC->Draw("same,histc");
  hNCell_TrueGammas_RecE_Effi_MC->Draw("same,histc");


  TLegend *leg = GetAndSetLegend2(0.5, 0.16, 0.95, 0.16 + 0.05*7, 40);
  leg->AddEntry(hNCell_TrueGammas_TrueE_Effi_MC, "P2, true #gamma, true E (no tagging)", "l");
  leg->AddEntry(hNCell_TrueGammas_RecE_Effi_MC, "P2, true #gamma (no tagging)", "l");
  leg->AddEntry(grTB_data, "TB, data, e^{#pm} (B=0T)", "p");
  leg->AddEntry(grTB_MC, "TB, MC, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_TrueE.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
  hDummyCan->SetLogy();
  hDummyCan->SetLogx();
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_TrueE_Log.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
  hDummyCan->SetLogy(0);
  hDummyCan->SetLogx(0);

}




void Effi::PlotRatio_Wide(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_RW_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_SBSub_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_RW_SBSub_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_RW_MCHadSub_Effi_Ratio->Draw("same,p");
  // hNCell_Gammas_RW_ConvMod10_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.55, 0.6, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_Gammas_Effi_Ratio, "P2, data", "p");
  leg->AddEntry(hNCell_Gammas_SBSub_Effi_Ratio, "P2, data, SB", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_Ratio, "P2, data, RW", "p");
  leg->AddEntry(hNCell_Gammas_RW_SBSub_Effi_Ratio, "P2, data, RW+SB", "p");
  leg->AddEntry(hNCell_Gammas_RW_MCHadSub_Effi_Ratio, "P2, data, RW+h^{#pm}sub", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.82,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Ratio.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotRatio_TBAndAll(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  // hNCell_Gammas_RW_ConvMod10_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.55, 0.6, 0.95, 0.84, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Ratio_TBAndAll.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotRatio_ConvMod(){
  hDummyCan->cd();
  hDummyRatio->Draw();
  DrawLines(0.5, PlotenergyHigh, 1,1, 2, kGray+2, 2);
  grTB_Ratio->Draw("same,p");

  hNCell_AllClus_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_RW_Effi_Ratio->Draw("same,p");
  // hNCell_Gammas_SBSub_Effi_Ratio->Draw("same,p");
  // hNCell_Gammas_RW_SBSub_Effi_Ratio->Draw("same,p");
  // hNCell_Gammas_RW_MCHadSub_Effi_Ratio->Draw("same,p");
  hNCell_Gammas_RW_ConvMod10_Effi_Ratio->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.65, 0.7, 0.95, 0.94, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Ratio, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_Gammas_Effi_Ratio, "P2, data", "p");
  // leg->AddEntry(hNCell_Gammas_SBSub_Effi_Ratio, "P2, data, SB", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_Ratio, "P2, data, RW", "p");
  leg->AddEntry(hNCell_Gammas_RW_ConvMod10_Effi_Ratio, "P2, data, RW+10% conv", "p");
  // leg->AddEntry(hNCell_Gammas_RW_MCHadSub_Effi_Ratio, "P2, data, RW+h^{#pm}sub", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.82,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.77,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Ratio_ConvMod.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotCorr_Wide(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1, 1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_AllClus_Effi_Corr->Draw("same,p");
  hNCell_Gammas_Effi_Corr->Draw("same,p");
  hNCell_Gammas_RW_Effi_Corr->Draw("same,p");
  hNCell_Gammas_SBSub_Effi_Corr->Draw("same,p");
  hNCell_Gammas_RW_SBSub_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.12, 0.76, 0.4, 0.96, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Corr, "P2, data, all clus", "p");
  leg->AddEntry(hNCell_Gammas_Effi_Corr, "P2, data", "p");
  leg->AddEntry(hNCell_Gammas_SBSub_Effi_Corr, "P2, data, SB", "p");
  leg->AddEntry(hNCell_Gammas_RW_Effi_Corr, "P2, data, RW", "p");
  leg->AddEntry(hNCell_Gammas_RW_SBSub_Effi_Corr, "P2, data, RW+SB", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.55,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.55,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.55,0.82,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.55,0.77,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Corr.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
  // cout<<Form("%s_%s/%s/Corr.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data())<<endl;
}

void Effi::PlotCorr_TBAndAll(){
  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1, 1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_AllClus_Effi_Corr->Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.15, 0.63, 0.4, 0.82, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Corr, "P2, data, all clus", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Corr_TBAndAll.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::FitCorr_GammasWide(){

  cout<<"Fitting SBsub and RW correction wth pol1:\n";
  TF1 *funcpol1 = new TF1("funcpol1", "[0] + [1]*x", 0.5, 10);
  TF1 *funcpol1_TB = new TF1("funcpol1_TB", "[0] + [1]*x", 0.5, 10);
  TF1 *funcGaus_All = new TF1("funcGaus_All", "gaus(0)", 0.5, 3);

  hNCell_Gammas_RW_SBSub_Effi_Corr->Fit(funcpol1, "MR0");
  funcpol1->SetLineColor(kOrange + 2);
  funcpol1->SetLineWidth(3);
  grTB_Corr->Fit(funcpol1_TB, "MR0");
  funcpol1_TB->SetLineColor(kBlue - 2);
  funcpol1_TB->SetLineWidth(3);
  hNCell_AllClus_Effi_Corr->Fit(funcGaus_All, "MR0");
  funcGaus_All->SetLineColor(kRed + 2);
  funcGaus_All->SetLineWidth(3);

  hDummyCan->cd();
  hDummyCorr->Draw();
  DrawLines(0.5, PlotenergyHigh, 1, 1, 2, kGray+2, 2);
  grTB_Corr->Draw("same,p");

  hNCell_AllClus_Effi_Corr->Draw("same,p");
  // hNCell_Gammas_Effi_Corr->Draw("same,p");
  // hNCell_Gammas_RW_Effi_Corr->Draw("same,p");
  // hNCell_Gammas_SBSub_Effi_Corr->Draw("same,p");
  hNCell_Gammas_RW_SBSub_Effi_Corr->Draw("same,p");

  funcpol1->Draw("same");
  funcpol1_TB->Draw("same");
  funcGaus_All->Draw("same");

  TLegend *leg = GetAndSetLegend2(0.15, 0.63, 0.4, 0.82, 40);
  leg->AddEntry(hNCell_AllClus_Effi_Corr, "P2, data, all clus", "p");
  // leg->AddEntry(hNCell_Gammas_Effi_Corr, "P2, data", "p");
  // leg->AddEntry(hNCell_Gammas_SBSub_Effi_Corr, "P2, data, SB", "p");
  // leg->AddEntry(hNCell_Gammas_RW_Effi_Corr, "P2, data, RW", "p");
  leg->AddEntry(hNCell_Gammas_RW_SBSub_Effi_Corr, "P2, data, RW+SB", "p");
  leg->AddEntry(grTB_Ratio, "TB, e^{#pm} (B=0T)", "l");
  leg->Draw("same");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("#gamma identified via. #pi^{0} tagging",0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.87,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedRange(),0.15,0.82,textSizeLabelsRel,kFALSE);
  // drawLatexAdd(" #gamma reweighted!",0.15,0.82,textSizeLabelsRel,kFALSE);
  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Corr_Fit.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));



}


void Effi::PlotExampleBin(){

  gSystem->Exec(Form("mkdir %s_%s/%s/%s/ExampleBins", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data()));
  hDummyCan->cd();

  std::vector<float> MInvBins = {0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.};

  TH1F* hMassPosData = new TH1F("hMassPosData", "", MInvBins.size() -1 , &MInvBins[0]);
  TH1F* hMassPosMC = new TH1F("hMassPosMC", "", MInvBins.size() -1 , &MInvBins[0]);

  for(int i = 0; i < MInvBins.size() - 1; ++i){

    float binLow = MInvBins.at(i) + 0.001;
    float binHigh = MInvBins.at(i+1) - 0.001;


    TH1F* hSE = (TH1F*) hInvMassVsClusterPt_MC->ProjectionX("sliceSameEvt", hInvMassVsClusterPt_MC->GetYaxis()->FindBin(binLow), hInvMassVsClusterPt_MC->GetYaxis()->FindBin(binHigh));
    TH1F* hSE_data = (TH1F*) hInvMassVsClusterPt_data->ProjectionX("sliceSameEvtData", hInvMassVsClusterPt_data->GetYaxis()->FindBin(binLow), hInvMassVsClusterPt_data->GetYaxis()->FindBin(binHigh));
    TH1F* hTrueGamma = (TH1F*) hInvMassVsPtGamma_MC->ProjectionX("sliceGamma", hInvMassVsPtGamma_MC->GetYaxis()->FindBin(binLow), hInvMassVsPtGamma_MC->GetYaxis()->FindBin(binHigh));
    TH1F* hTrueElec = (TH1F*) hInvMassVsPtElec_MC->ProjectionX("sliceElec", hInvMassVsPtElec_MC->GetYaxis()->FindBin(binLow), hInvMassVsPtElec_MC->GetYaxis()->FindBin(binHigh));
    TH1F* hBack = (TH1F*) hInvMassVsClusterPtBack_MC->ProjectionX("sliceBack", hInvMassVsClusterPtBack_MC->GetYaxis()->FindBin(binLow), hInvMassVsClusterPtBack_MC->GetYaxis()->FindBin(binHigh));
    TH1F* hContamination = (TH1F*) hSE->Clone("hContamination");
    hContamination->Add(hTrueGamma, -1);
    hContamination->Add(hTrueElec, -1);

    hDummyMInv->GetYaxis()->SetRangeUser(0., hSE->GetMaximum()*1.4);
    hDummyMInv->Draw();

    hBack->Scale(hSE->Integral(hSE->GetXaxis()->FindBin(0.2), hSE->GetXaxis()->FindBin(0.3))  /  hBack->Integral(hBack->GetXaxis()->FindBin(0.2), hBack->GetXaxis()->FindBin(0.3)));
    hSE_data->Scale(hSE->Integral(hSE->GetXaxis()->FindBin(0.2), hSE->GetXaxis()->FindBin(0.3))  /  hSE_data->Integral(hSE_data->GetXaxis()->FindBin(0.2), hSE_data->GetXaxis()->FindBin(0.3)));
    DrawSetMarker(hSE, 20, 2, kBlack, kBlack);
    DrawSetMarker(hTrueGamma, 2, 3, kRed + 2, kRed + 2);
    DrawSetMarker(hTrueElec, 8, 3, kGreen + 2, kGreen + 2);
    DrawSetMarker(hContamination, 9, 3, kBlue + 2, kBlue + 2);
    DrawSetMarker(hBack, 24, 2, kGray + 2, kGray + 2);

    DrawSetMarker(hSE_data, 24, 2, kOrange + 2, kOrange + 2);

    hSE->Draw("same, p");
    hTrueGamma->Draw("same,histc");
    hTrueElec->Draw("same,histc");
    hContamination->Draw("same,histc");
    hBack->Draw("same,p");


    drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(StrSelectedGamma(),0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(Form("%.01f < E (clus) #leq %.01f GeV", binLow, binHigh),0.15,0.92 - 0.05,textSizeLabelsRel,kFALSE);

    DrawLines(0.09, 0.09, 0., hSE->GetMaximum() * 0.3, 2, 2, kGray + 2 );
    DrawLines(0.17, 0.17, 0., hSE->GetMaximum() * 0.3, 2, 2, kGray + 2 );

    TLegend *leg = GetAndSetLegend2(0.67, 0.82-5*0.05, 0.95, 0.82, 40);
    leg->AddEntry(hSE, "same evt.", "p");
    leg->AddEntry(hBack, "rot. back.", "p");
    leg->AddEntry(hTrueGamma, "true #gamma clus", "l");
    leg->AddEntry(hTrueElec, "true e^{#pm} clus", "l");
    leg->AddEntry(hContamination, "hadron clus.", "l");
    leg->Draw("same");

    hDummyCan->SaveAs(Form("%s_%s/%s/%s/ExampleBins/ExampleBin_%i.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), i,  fsuffix.Data()));

    hSE_data->Draw("same,p");
    leg->AddEntry(hSE_data, "same evt. data", "p");
    leg->Draw("same");
    hDummyCan->SaveAs(Form("%s_%s/%s/%s/ExampleBins/ExampleBin_wData_%i.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), i,  fsuffix.Data()));

    // Plot ratio of background to different components
    hDummyMInvRatio->GetYaxis()->SetRangeUser(0., 2.);
    hDummyMInvRatio->Draw("grid");

    TH1F* hTrueGammaRatio = (TH1F*) hTrueGamma->Clone(Form("hTrueGammaRatio_%i", i));
    TH1F* hTrueElecRatio = (TH1F*) hTrueElec->Clone(Form("hTrueElecRatio_%i", i));
    TH1F* hContaminationRatio = (TH1F*) hContamination->Clone(Form("hContaminationRatio_%i", i));

    hTrueGammaRatio->Divide(hBack);
    hTrueElecRatio->Divide(hBack);
    hContaminationRatio->Divide(hBack);

    hTrueGammaRatio->Draw("same,histc");
    hTrueElecRatio->Draw("same,histc");
    hContaminationRatio->Draw("same,histc");

    leg->Draw("same");

    drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(StrSelectedGamma(),0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(Form("%.01f < p_{T}(clus) #leq %.01f GeV/#it{c}", binLow, binHigh),0.15,0.92 - 0.05,textSizeLabelsRel,kFALSE);


    hDummyCan->SaveAs(Form("%s_%s/%s/%s/ExampleBins/ExampleBin_Ratio_%i.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), i,  fsuffix.Data()));

    TF1* fgaus = new TF1("fgaus", "gaus(0)", 0.09, 0.15);
    hSE->Fit(fgaus, "MQR0");
    float mean = fgaus->GetParameter(1);
    float sig = fgaus->GetParameter(2);
    fgaus->SetRange(mean - sig, mean + sig);
    hSE->Fit(fgaus, "MQR0");
    hMassPosMC->SetBinContent(i+1, fgaus->GetParameter(1)/0.135);
    hMassPosMC->SetBinError(i+1, fgaus->GetParError(1)/0.135);

    fgaus->SetRange(0.09, 0.15);
    hSE_data->Fit(fgaus, "MQR0");
    mean = fgaus->GetParameter(1);
    sig = fgaus->GetParameter(2);
    fgaus->SetRange(mean - sig, mean + sig);
    hSE_data->Fit(fgaus, "MQR0");
    hMassPosData->SetBinContent(i+1, fgaus->GetParameter(1)/0.135);
    hMassPosData->SetBinError(i+1, fgaus->GetParError(1)/0.135);

  }

  hDummyCan->cd();
  DrawSetMarker(hMassPosData, 20, 2, kBlack, kBlack);
  DrawSetMarker(hMassPosMC, 24, 2, kRed + 1, kRed + 1);
  hDummyMassPos->Draw();
  hMassPosData->Draw("same,pe");
  hMassPosMC->Draw("same,pe");

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);

  TLegend *leg = GetAndSetLegend2(0.15, 0.93-2*0.05, 0.3, 0.93, 40);
  leg->AddEntry(hMassPosData, "data", "p");
  leg->AddEntry(hMassPosMC, "MC", "p");
  leg->Draw();

  hDummyCan->SaveAs(Form("%s_%s/%s/%s/MassPos.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(),  fsuffix.Data()));

}


void Effi::PlotTrueVsRecE(){
  gSystem->Exec(Form("mkdir %s_%s/%s/%s/TrueVsRecE", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data()));

  std::vector<float> EBins = {0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.};
  TH1D* hMeanMassPos = new TH1D("hMeanMassPos", "", EBins.size() -1, &EBins[0]);
  TH1D* hMeanMassPosOneCell = new TH1D("hMeanMassPosOneCell", "", EBins.size() -1, &EBins[0]);
  TH1D* hMeanMassPosTwoCell = new TH1D("hMeanMassPosTwoCell", "", EBins.size() -1, &EBins[0]);
  TH1D* hMeanMassPosThreeCell = new TH1D("hMeanMassPosThreeCell", "", EBins.size() -1, &EBins[0]);


  for(int i = 0; i < EBins.size() - 1; ++i){
    hDummyCan->cd();
    float binLow = EBins.at(i) + 0.001;
    float binHigh = EBins.at(i+1) - 0.001;

    TH1F* hRecDivTrueE_plot = (TH1F*) hRecDivTrueE->ProjectionY(Form("hRecDivTrueE%i",i), hRecDivTrueE->GetXaxis()->FindBin(binLow), hRecDivTrueE->GetXaxis()->FindBin(binHigh));
    TH1F* hRecDivTrueEOneCell_plot = (TH1F*) hRecDivTrueEOneCell->ProjectionY(Form("hRecDivTrueEOneCell%i",i), hRecDivTrueE->GetXaxis()->FindBin(binLow), hRecDivTrueE->GetXaxis()->FindBin(binHigh));
    TH1F* hRecDivTrueETwoCell_plot = (TH1F*) hRecDivTrueETwoCell->ProjectionY(Form("hRecDivTrueETwoCell%i",i), hRecDivTrueE->GetXaxis()->FindBin(binLow), hRecDivTrueE->GetXaxis()->FindBin(binHigh));
    TH1F* hRecDivTrueEThreeCell_plot = (TH1F*) hRecDivTrueEThreeCell->ProjectionY(Form("hRecDivTrueEThreeCell%i",i), hRecDivTrueE->GetXaxis()->FindBin(binLow), hRecDivTrueE->GetXaxis()->FindBin(binHigh));

    DrawSetMarker(hRecDivTrueE_plot, 2, 3, kRed + 2, kRed + 2);
    DrawSetMarker(hRecDivTrueEOneCell_plot, 8, 3, kGreen + 2, kGreen + 2);
    DrawSetMarker(hRecDivTrueETwoCell_plot, 9, 3, kBlue + 2, kBlue + 2);
    DrawSetMarker(hRecDivTrueEThreeCell_plot, 24, 2, kGray + 2, kGray + 2);


    hRecDivTrueE_plot->Scale(1./hRecDivTrueE_plot->Integral());
    hRecDivTrueEOneCell_plot->Scale(1./hRecDivTrueEOneCell_plot->Integral());
    hRecDivTrueETwoCell_plot->Scale(1./hRecDivTrueETwoCell_plot->Integral());
    hRecDivTrueEThreeCell_plot->Scale(1./hRecDivTrueEThreeCell_plot->Integral());
    hDummyERecVsETrue->GetYaxis()->SetRangeUser(0., hRecDivTrueE_plot->GetMaximum()*1.5);
    hDummyERecVsETrue->Draw();
    hRecDivTrueE_plot->Draw("same,histc");
    hRecDivTrueEOneCell_plot->Draw("same,histc");
    hRecDivTrueETwoCell_plot->Draw("same,histc");
    hRecDivTrueEThreeCell_plot->Draw("same,histc");

    TF1 *func = new TF1("func", "gaus(0)", 0.8, 1.2);
    hRecDivTrueE_plot->Fit(func, "MQR0");
    DrawLines(func->GetParameter(1), func->GetParameter(1), 0, hRecDivTrueE_plot->GetMaximum()*0.2, 3, kRed + 2, 2);
    hMeanMassPos->SetBinContent(i+1, func->GetParameter(1));
    hMeanMassPos->SetBinError(i+1, func->GetParError(1));

    hRecDivTrueEOneCell_plot->Fit(func, "MQR0");
    DrawLines(func->GetParameter(1), func->GetParameter(1), 0, hRecDivTrueE_plot->GetMaximum()*0.2, 3, kGreen + 2, 2);
    hMeanMassPosOneCell->SetBinContent(i+1, func->GetParameter(1));
    hMeanMassPosOneCell->SetBinError(i+1, func->GetParError(1));

    hRecDivTrueETwoCell_plot->Fit(func, "MQR0");
    DrawLines(func->GetParameter(1), func->GetParameter(1), 0, hRecDivTrueE_plot->GetMaximum()*0.2, 3, kBlue + 2, 2);
    hMeanMassPosTwoCell->SetBinContent(i+1, func->GetParameter(1));
    hMeanMassPosTwoCell->SetBinError(i+1, func->GetParError(1));

    hRecDivTrueEThreeCell_plot->Fit(func, "MQR0");
    DrawLines(func->GetParameter(1), func->GetParameter(1), 0, hRecDivTrueE_plot->GetMaximum()*0.2, 3, kGray + 2, 2);
    hMeanMassPosThreeCell->SetBinContent(i+1, func->GetParameter(1));
    hMeanMassPosThreeCell->SetBinError(i+1, func->GetParError(1));

    TLegend *leg = GetAndSetLegend2(0.15, 0.82-4*0.05, 0.3, 0.82, 40);
    leg->AddEntry(hRecDivTrueE_plot, "all clus", "l");
    leg->AddEntry(hRecDivTrueEOneCell_plot, "1 cell clus", "l");
    leg->AddEntry(hRecDivTrueETwoCell_plot, "2 cell clus", "l");
    leg->AddEntry(hRecDivTrueEThreeCell_plot, "3 cell clus", "l");
    leg->Draw();

    drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd("true #gamma clus",0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd(Form("%.01f < E_(clus) #leq %.01f GeV", binLow, binHigh),0.15,0.92 - 0.05,textSizeLabelsRel,kFALSE);


    hDummyCan->SaveAs(Form("%s_%s/%s/%s/TrueVsRecE/TrueVsRecE_%i.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), i, fsuffix.Data()));


  }

  DrawSetMarker(hMeanMassPos, 20, 3, kRed + 2, kRed + 2);
  DrawSetMarker(hMeanMassPosOneCell, 24, 3, kGreen + 2, kGreen + 2);
  DrawSetMarker(hMeanMassPosTwoCell, 27, 3, kBlue + 2, kBlue + 2);
  DrawSetMarker(hMeanMassPosThreeCell, 34, 2, kGray + 2, kGray + 2);
  hDummyCan->cd();

  hDummyERecVsETrue2->Draw();
  hMeanMassPos->Draw("same,pe");
  hMeanMassPosOneCell->Draw("same,pe");
  hMeanMassPosTwoCell->Draw("same,pe");
  hMeanMassPosThreeCell->Draw("same,pe");

  TLegend *leg = GetAndSetLegend2(0.15, 0.9-4*0.05, 0.3, 0.9, 40);
  leg->AddEntry(hMeanMassPos, "all clus", "l");
  leg->AddEntry(hMeanMassPosOneCell, "1 cell clus", "l");
  leg->AddEntry(hMeanMassPosTwoCell, "2 cell clus", "l");
  leg->AddEntry(hMeanMassPosThreeCell, "3 cell clus", "l");
  leg->Draw();

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("true #gamma clus",0.15,0.92,textSizeLabelsRel,kFALSE);

  hDummyCan->SaveAs(Form("%s_%s/%s/%s/TrueVsRecE.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(),fsuffix.Data()));


}


void Effi::PlotPurity(){
  hDummyCan->cd();

  hDummyPurity->Draw();


  DrawSetMarker(hGammaPurity, 1, 3, kRed + 3, kRed + 3);
  DrawSetMarker(hElecPurity, 1, 3, kGreen + 3, kGreen + 3);
  DrawSetMarker(hUnPurity, 1, 3, kBlue + 3, kBlue + 3);

  DrawSetMarker(hGammaPuritySB, 7, 3, kRed+ 3, kRed+ 3);
  DrawSetMarker(hElecPuritySB, 7, 3, kGreen+ 3, kGreen+ 3);
  DrawSetMarker(hUnPuritySB, 7, 3, kBlue+ 3, kBlue + 3);

  DrawSetMarker(hGammaPuritySBSub, 2, 3, kRed+ 3, kRed+ 3);
  DrawSetMarker(hElecPuritySBSub, 2, 3, kGreen+ 3, kGreen+ 3);
  DrawSetMarker(hUnPuritySBSub, 2, 3, kBlue+ 3, kBlue + 3);

  hGammaPurity->Draw("same,histc");
  hElecPurity->Draw("same,histc");
  hUnPurity->Draw("same,histc");
  hGammaPuritySB->Draw("same,histc");
  hElecPuritySB->Draw("same,histc");
  hUnPuritySB->Draw("same,histc");
  hGammaPuritySBSub->Draw("same,histc");
  hElecPuritySBSub->Draw("same,histc");
  hUnPuritySBSub->Draw("same,histc");

  TLegend *leg = GetAndSetLegend2(0.2, 0.72, 0.45, 0.83, 40);
  leg->AddEntry(hGammaPurity, "#gamma Signal", "l");
  leg->AddEntry(hElecPurity, "e^{#pm} Signal", "l");
  leg->AddEntry(hUnPurity, "h^{#pm} Signal", "l");

  TLegend *leg2 = GetAndSetLegend2(0.45, 0.72, 0.7, 0.83, 40);
  leg2->AddEntry(hGammaPuritySB, "#gamma SB", "l");
  leg2->AddEntry(hElecPuritySB, "e^{#pm} SB", "l");
  leg2->AddEntry(hUnPuritySB, "h^{#pm} SB", "l");

  TLegend *leg3 = GetAndSetLegend2(0.7, 0.72, 0.95, 0.83, 40);
  leg3->AddEntry(hGammaPuritySBSub, "#gamma SB sub.", "l");
  leg3->AddEntry(hElecPuritySBSub, "e^{#pm} SB sub.", "l");
  leg3->AddEntry(hUnPuritySBSub, "h^{#pm} SB sub.", "l");

  leg->Draw("same");
  leg2->Draw("same");
  leg3->Draw("same");

  DrawLines(0.5, PlotenergyHigh, 0, 0, 2, kGray + 2, 2);

  drawLatexAdd(sEnergy,0.7,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("Signal: " + StrSelectedRange(),0.15,0.92 - 0.05,textSizeLabelsRel,kFALSE);

  hDummyCan->SaveAs(Form("%s_%s/%s/%s/Purity.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

}

void Effi::PlotNCellVsE(){

  hDummyCan2d->cd();

  SetStyleHistoTH2ForGraphs(hNCellVsEGammasNL_data, "#it{E}_{clus} (GeV)","N_{Ncells}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);
  hNCellVsEGammasNL_data->Draw("colz");
  hNCellVsEGammasNL_data->GetXaxis()->SetMoreLogLabels(1);
  hDummyCan2d->SetLogx();
  // hDummyCan2d->SetLogy();
  hDummyCan2d->SetLogz();
  drawLatexAdd(sEnergy,0.6,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd(StrSelectedGamma(),0.15,0.92,textSizeLabelsRel,kFALSE);
  drawLatexAdd("Signal: " + StrSelectedRange(),0.15,0.92 - 0.05,textSizeLabelsRel,kFALSE);
  hDummyCan2d->SaveAs(Form("%s_%s/%s/%s/NCellVsE.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::PlotEffiSources(){

    if(!hNCellVsEAllClusTrueGamma){cout<<"Error loading hNCellVsEAllClusTrueGamma..."; return;}

    TH1D*  hEffiTrueGammas= nullptr;
    TH1D*  hEffiTrueElec= nullptr;
    TH1D*  hEffiTrueHadrons= nullptr;
    GetTrueHists(hNCellVsEAllClusTrueGamma, hEffiTrueGammas);
    GetTrueHists(hNCellVsEAllClusTrueElec, hEffiTrueElec);
    GetTrueHists(hNCellVsEAllClusTrueHadrons, hEffiTrueHadrons);

    DrawSetMarker(hEffiTrueGammas, 20, 1.9, kRed + 1, kRed + 1);
    DrawSetMarker(hEffiTrueElec, 34, 1.9, kBlue + 1, kBlue + 1);
    DrawSetMarker(hEffiTrueHadrons, 34, 1.9, kGray + 2, kGray + 2);

    hDummyCan->cd();
    hDummyEffi->Draw();

    hEffiTrueGammas->Draw("same,pe");
    hEffiTrueElec->Draw("same,pe");
    hEffiTrueHadrons->Draw("same,pe");
    grTB_MC->Draw("same,l");



    TLegend *leg = GetAndSetLegend2(0.5, 0.82, 0.7, 0.96, 40);
    leg->AddEntry(hEffiTrueGammas, "true #gamma", "p");
    leg->AddEntry(hEffiTrueElec, "true e^{#pm}", "p");
    leg->AddEntry(hEffiTrueHadrons, "true hadrons", "p");
    leg->AddEntry(grTB_MC , "TB e^{#pm} (B=0T), no material", "l");
    leg->Draw("same");

    drawLatexAdd(sEnergy,0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd("rec. E_{clus} used for P2 MC",0.15,0.87,textSizeLabelsRel,kFALSE);

    DrawLines(0.5, PlotenergyHigh, 1, 1, 2, kGray + 2, 2);

    hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_TrueContributions.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));

    TH1D*  hEffiTrueGammasAbs= (TH1D*) hNCellVsEAllClusTrueGamma->ProjectionX("hEffiTrueGammasAbs", 1, 10);
    TH1D*  hEffiTrueElecAbs= (TH1D*) hNCellVsEAllClusTrueElec->ProjectionX("hEffiTrueElecAbs", 1, 10);
    TH1D*  hEffiTrueHadronsAbs= (TH1D*) hNCellVsEAllClusTrueHadrons->ProjectionX("hEffiTrueHadronsAbs", 1, 10);

    TH1D*  hEffiAllClus = (TH1D*) hEffiTrueHadronsAbs->Clone("hEffiAllClus");
    hEffiAllClus->Add(hEffiTrueElecAbs);
    hEffiAllClus->Add(hEffiTrueGammasAbs);

    hEffiTrueGammasAbs->Divide(hEffiAllClus);
    hEffiTrueElecAbs->Divide(hEffiAllClus);
    hEffiTrueHadronsAbs->Divide(hEffiAllClus);

    DrawSetMarker(hEffiTrueGammasAbs, 20, 3.5, kRed + 1, kRed + 1);
    DrawSetMarker(hEffiTrueElecAbs, 34, 3.5, kBlue + 1, kBlue + 1);
    DrawSetMarker(hEffiTrueHadronsAbs, 34, 3.5, kGray + 2, kGray + 2);


    hDummyCan->cd();
    hDummyFraction->Draw();
    hEffiTrueGammasAbs->Draw("same,hist");
    hEffiTrueElecAbs->Draw("same,hist");
    hEffiTrueHadronsAbs->Draw("same,hist");

    TLegend *leg2 = GetAndSetLegend2(0.7, 0.82, 0.95, 0.93, 40);
    leg2->AddEntry(hEffiTrueGammas, "true #gamma", "l");
    leg2->AddEntry(hEffiTrueElec, "true e^{#pm}", "l");
    leg2->AddEntry(hEffiTrueHadrons, "true hadrons", "l");
    leg2->Draw("same");

    drawLatexAdd(sEnergy,0.15,0.92,textSizeLabelsRel,kFALSE);
    drawLatexAdd("rec. E_{clus} used for P2 MC",0.15,0.87,textSizeLabelsRel,kFALSE);

    hDummyCan->SaveAs(Form("%s_%s/%s/%s/Effi_TrueContributions_Fraction.%s", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data(), fsuffix.Data()));
}


void Effi::WriteToFile(){

  TFile fout(Form("%s_%s/%s/%s/histos.root", fPeriod.Data(), fMethod.Data(), fSpecialName.Data(), fsuffix.Data()), "Recreate");
  fout.cd();

  hNCell_Gammas_RW_SBSub_Effi_data->Write("hNCell_Gammas_RW_SBSub_Effi_data", TObject::kOverwrite);
  hNCell_Gammas_RW_SBSub_Effi_MC->Write("hNCell_Gammas_RW_SBSub_Effi_MC", TObject::kOverwrite);
  hNCell_Gammas_RW_SBSub_Effi_Ratio->Write("hNCell_Gammas_RW_SBSub_Effi_Ratio", TObject::kOverwrite);
  hNCell_Gammas_RW_SBSub_Effi_Corr->Write("hNCell_Gammas_RW_SBSub_Effi_Corr", TObject::kOverwrite);

  hNCell_AllClus_Effi_data->Write("hNCell_AllClus_Effi_data", TObject::kOverwrite);
  hNCell_AllClus_Effi_MC->Write("hNCell_AllClus_Effi_MC", TObject::kOverwrite);
  hNCell_AllClus_Effi_Ratio->Write("hNCell_AllClus_Effi_Ratio", TObject::kOverwrite);
  hNCell_AllClus_Effi_Corr->Write("hNCell_AllClus_Effi_Corr", TObject::kOverwrite);

  fout.Close();
}



void PlottingNCellEffi_Methods(TString period = "13TeVNomB"){
  TString fSpecialName = "";
  {
    std::vector<TString> Method = {"Wide"};
    // std::vector<TString> Method = {"Wide", "Low", "High", "Left"};
    for(auto &i : Method){

      fSpecialName = "Pi0Tagging_13TeV_nom_04_26_WithTRD_WithBorderCells_1cellFT";

      // Effi plot("Pi0Tagging_13TeV_nom_02_28_TBNL.root", "13TeVNomB", "pdf");
      // Effi plot("Pi0Tagging_13TeV_nom_03_11_Iso02_TBNL_noScale.root", "13TeVNomB", "pdf"); // high stat with std. settings
      Effi plot("rootFiles/Pi0Tagging_13TeV_nom_04_26_WithTRD_WithBorderCells_1cellFT.root", "13TeVNomB", "pdf"); // high stat wo fine tuning MC scale (1.00) data scale (1.015)
      // Effi plot("rootFiles/Pi0Tagging_13TeV_nom_03_17_Iso02_TBNL_noScale_MCScaled_1cellFT.root", "13TeVNomB", "pdf"); // high stat wo fine tuning MC scale (1.00) data scale (1.015)
      // Effi plot("Pi0Tagging_13TeV_nom_03_11_Iso02_TBNL_NoFT.root", "13TeVNomB", "pdf"); // no Fine tuning
      plot.SetSpecialName(fSpecialName);
      plot.SetAddition(i);
      plot.FillHistos();
      plot.SetOtherHistos();
      plot.DoSBSubtraction();
      plot.SetPurityHistosAfterSBSub();
      plot.GetHistReweighted();
      plot.FillCorrHistos();
      cout<<__LINE__<<endl;
      plot.LoadTB();
      cout<<__LINE__<<endl;
      plot.SetPlotting();
      plot.PlotEffi_TB();
      plot.PlotEffi_AllClusAndTB();
      plot.PlotEffi_AllClusAndTBAndTrueGamma();
      plot.PlotEffi_Wide();
      plot.PlotEffi_WideWithTrue();
      plot.PlotEffi_WideWithRW();
      plot.PlotEffi_WideWithRWSBSub();
      plot.PlotEffi_TrueVsRecE();

      plot.PlotRatio_Wide();
      plot.PlotRatio_TBAndAll();
      plot.PlotRatio_ConvMod();

      plot.PlotCorr_TBAndAll();
      plot.PlotCorr_Wide();

      plot.FitCorr_GammasWide();

      plot.PlotMCClosure();
      plot.PlotDataClosure();

      plot.PlotExampleBin();
      plot.PlotPurity();
      plot.PlotNCellVsE();
      plot.PlotTrueVsRecE();

      plot.PlotEffiSources();

      plot.WriteToFile();
    }
  }


//   {
//   // else if(period.Contains("8TeV")){
//   std::vector<TString> Method = {"Wide", "Low", "High", "Left"};
//   for(auto &i : Method){
//     // gSystem->Exec(Form("mkdir -p 8TeV_%s/pdf", i.Data()));
//     // gSystem->Exec(Form("mkdir -p 8TeV_%s/png", i.Data()));
//     // Effi plot("Pi0Tagging_13TeV_nom_02_28_TBNL.root", "13TeVNomB", "pdf");
//     Effi plot("rootFiles/Pi0Tagging_8TeV_03_12_NoTRD_NoShaper_TBNL.root", "8TeV", "pdf"); // high stat with std. settings
//     // Effi plot("Pi0Tagging_13TeV_nom_03_11_Iso02_TBNL_NoFT.root", "13TeVNomB", "pdf"); // no Fine tuning
//     plot.SetAddition(i);
//     plot.FillHistos();
//     plot.SetOtherHistos();
//     plot.DoSBSubtraction();
//     plot.GetHistReweighted();
//     plot.FillCorrHistos();
//     cout<<__LINE__<<endl;
//     plot.LoadTB();
//     cout<<__LINE__<<endl;
//     plot.SetPlotting();
//     plot.PlotEffi_TB();
//     plot.PlotEffi_AllClusAndTB();
//     plot.PlotEffi_Wide();
//     plot.PlotEffi_WideWithTrue();
//     plot.PlotEffi_WideWithRW();
//     plot.PlotEffi_WideWithRWSBSub();
//
//     plot.PlotRatio_Wide();
//
//     plot.PlotCorr_TBAndAll();
//     plot.PlotCorr_Wide();
//
//     plot.PlotMCClosure();
//     plot.PlotDataClosure();
//
//     plot.PlotExampleBin();
//     plot.PlotPurity();
//   }
// }

}
