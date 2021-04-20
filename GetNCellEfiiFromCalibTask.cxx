#include "header.h"
#include "Plotting_Header.h"
#include "/home/joshua/PCG_Software/Plotting/Plotting_Class.h"

float calcErr(float Emc, float Errmc, float Edata, float Errdata);
void GetEffiHists(TH2F *hdata2d, TH2F *hMC2d, TH1D *&hEffiData,  TH1D *&hEffiMC, TH1D *&hRatio, TH1D *&hCorr);
void LoadTB();

TGraphErrors* grTB_data;
TGraphErrors* grTB_MC;
TGraphErrors* grTB_Ratio;
TGraphErrors* grTB_Corr;

void GetNCellEfiiFromCalibTask(TString method = "PCMEMC"){

TFile *fData = TFile::Open(Form("/media/joshua/external_drive/Data/2021_04_04_pp13TeV_NCell/%s_data.root", method.Contains("PCM") ? "ConvCaloCalibration_1_7" : "ConvCaloCalibration_2_46"));
TFile *fMC = TFile::Open(Form("/media/joshua/external_drive/Data/2021_04_04_pp13TeV_NCell/%s_MC.root", method.Contains("PCM") ? "ConvCaloCalibration_1_7" : "ConvCaloCalibration_2_46"));

TH2F* hDataTBScale = (TH2F*) GetHistFromFile(fData, "ESD_EVsNcell_Pi0Tagged", method.Contains("PCM") ? "00010113_00200009327000008250400000_411793706f030000000_0163103100000010" : "2_00010113_411793706f030000000_01631031000000d0");
TH2F* hDataTBScaleFT = (TH2F*) GetHistFromFile(fData, "ESD_EVsNcell_Pi0Tagged", method.Contains("PCM") ? "00010113_00200009327000008250400000_411790106f030000000_0163103100000010" : "2_00010113_411790106f030000000_01631031000000d0");
TH2F* hDataTBNoScale = (TH2F*) GetHistFromFile(fData, "ESD_EVsNcell_Pi0Tagged", method.Contains("PCM") ? "00010113_00200009327000008250400000_411790406f030000000_0163103100000010" : "2_00010113_411790406f030000000_01631031000000d0");
TH2F* hMCTB = (TH2F*) GetHistFromFile(fMC, "ESD_EVsNcell_Pi0Tagged", method.Contains("PCM") ? "00010113_00200009327000008250400000_411793706f030000000_0163103100000010" : "2_00010113_411793706f030000000_01631031000000d0");
TH2F* hMCTBFT = (TH2F*) GetHistFromFile(fMC, "ESD_EVsNcell_Pi0Tagged", method.Contains("PCM") ? "00010113_00200009327000008250400000_411790106f030000000_0163103100000010" : "2_00010113_411790106f030000000_01631031000000d0");


TH1D* hData_Effi_TBScale = nullptr;
TH1D* hMC_Effi_TB = nullptr;
TH1D* h_Ratio_TBScale = nullptr;
TH1D* h_Corr_TBScale = nullptr;

TH1D* hData_Effi_TBScaleFT = nullptr;
TH1D* hMC_Effi_TBFT = nullptr;
TH1D* h_Ratio_TBScaleFT = nullptr;
TH1D* h_Corr_TBScaleFT = nullptr;

TH1D* hData_Effi_TBNoScale = nullptr;
TH1D* h_Ratio_TBNoScale = nullptr;
TH1D* h_Corr_TBNoScale = nullptr;


LoadTB();

GetEffiHists(hDataTBScale, hMCTB, hData_Effi_TBScale, hMC_Effi_TB, h_Ratio_TBScale, h_Corr_TBScale);
GetEffiHists(hDataTBScaleFT, hMCTBFT, hData_Effi_TBScaleFT, hMC_Effi_TBFT, h_Ratio_TBScaleFT, h_Corr_TBScaleFT);
GetEffiHists(hDataTBNoScale, hMCTB, hData_Effi_TBNoScale, hMC_Effi_TB, h_Ratio_TBNoScale, h_Corr_TBNoScale);


DrawSetMarker(hData_Effi_TBScale, 20, 2, kCyan + 2, kCyan + 2);
DrawSetMarker(hMC_Effi_TB, 20, 2, kCyan + 2, kCyan + 2);
DrawSetMarker(h_Ratio_TBScale, 20, 2, kCyan + 2, kCyan + 2);
DrawSetMarker(h_Corr_TBScale, 20, 2, kCyan + 2, kCyan + 2);

DrawSetMarker(hData_Effi_TBScaleFT, 28, 2, kGreen+3, kGreen+3);
DrawSetMarker(hMC_Effi_TBFT, 28, 2, kGreen+3, kGreen+3);
DrawSetMarker(h_Ratio_TBScaleFT, 28, 2, kGreen+3, kGreen+3);
DrawSetMarker(h_Corr_TBScaleFT, 28, 2, kGreen+3, kGreen+3);


DrawSetMarker(hData_Effi_TBNoScale, 27, 3, kRed + 2, kRed + 2);
DrawSetMarker(h_Ratio_TBNoScale, 27, 3, kRed + 2, kRed + 2);
DrawSetMarker(h_Corr_TBNoScale, 27, 3, kRed + 2, kRed + 2);



Float_t textSizeSinglePad = 0.044;
Float_t textSizeLabelsRel = 0.044;

TH2F* hDummyRatio    = new TH2F("hDummyRatio","hDummyRatio",1000,0.5, 8,1000,0.9, 1.4);
SetStyleHistoTH2ForGraphs(hDummyRatio, "#it{E} (GeV)","#nu_{MC}/#nu_{data}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);
hDummyRatio->SetStats(0);

TH2F* hDummyEffi    = new TH2F("hDummyEffi","hDummyEffi",1000,0.5, 8,1000,0.1, 1.4);
SetStyleHistoTH2ForGraphs(hDummyEffi, "#it{E} (GeV)","#nu", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);
hDummyEffi->SetStats(0);

TH2F* hDummyCorr    = new TH2F("hDummyCorr","hDummyCorr",1000,0.5, 8,1000,0., 1.4);
SetStyleHistoTH2ForGraphs(hDummyCorr, "#it{E} (GeV)","#rho", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 1.1,1.2, 510, 510);
hDummyCorr->SetStats(0);


TCanvas *hDummyCan = new TCanvas("Can", "", 1200, 1000);
DrawPaperCanvasSettings(hDummyCan, 0.1, 0.003, 0.003, 0.1);

hDummyCan->cd();

hDummyEffi->Draw();

grTB_data->Draw("same,pe");
grTB_MC->Draw("same,l");
hData_Effi_TBScale->Draw("same,pe");
hData_Effi_TBScaleFT->Draw("same,pe");
hData_Effi_TBNoScale->Draw("same,pe");
hMC_Effi_TB->Draw("same,histc");
hMC_Effi_TBFT->Draw("same,histc");

gSystem->Exec("mkdir PlotsCalibTrain");
hDummyCan->SaveAs(Form("PlotsCalibTrain/Effi_%s.pdf", method.Data()));


hDummyRatio->Draw();

grTB_Ratio->Draw("same,pe");
h_Ratio_TBScale->Draw("same,pe");
h_Ratio_TBScaleFT->Draw("same,pe");
h_Ratio_TBNoScale->Draw("same,pe");
hDummyCan->SaveAs(Form("PlotsCalibTrain/Ratio_%s.pdf", method.Data()));

hDummyCorr->Draw();

grTB_Corr->Draw("same,pe");
h_Corr_TBScale->Draw("same,pe");
h_Corr_TBScaleFT->Draw("same,pe");
h_Corr_TBNoScale->Draw("same,pe");
hDummyCan->SaveAs(Form("PlotsCalibTrain/Corr_%s.pdf", method.Data()));

}

// ------------------------------------------
// Get NCell efficiency hists from 2d NCellVs E distributions
// ------------------------------------------
void GetEffiHists(TH2F *hdata2d, TH2F *hMC2d, TH1D *&hEffiData,  TH1D *&hEffiMC, TH1D *&hRatio, TH1D *&hCorr){

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
float calcErr(float Emc, float Errmc, float Edata, float Errdata){
  float err1 = Errmc / (1-Edata) ;
  float err2 =  ((1 - Emc) * Errdata ) / ((1-Edata)*(1-Edata) );
  return sqrt(err1*err1 + err2*err2);
}


void LoadTB(){

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
