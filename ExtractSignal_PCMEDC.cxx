#include "/home/joshua/PCG_Software/Plotting/Plotting_Class.h"
#include "header.h"
#include "Plotting_Header.h"


TGraphErrors* grTB_data;
TGraphErrors* grTB_MC;
TGraphErrors* grTB_Ratio;
TGraphErrors* grTB_Corr;

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

float GetYield(TH2F* hInvMass, TH2F* hInvMassBack, float pTLow, float pTUp, double &err, TString name = ""){
  TH1F* hSignal = (TH1F*) hInvMass->ProjectionX("hSignal", hInvMass->GetYaxis()->FindBin(pTLow + 0.001), hInvMass->GetYaxis()->FindBin(pTUp - 0.001));
  TH1F* hBack = (TH1F*) hInvMassBack->ProjectionX("hBack", hInvMassBack->GetYaxis()->FindBin(pTLow + 0.001), hInvMassBack->GetYaxis()->FindBin(pTUp - 0.001));


  /// scale back to Signal
  float scSig = hSignal->Integral(hSignal->FindBin(0.2001), hSignal->FindBin(0.3001));
  cout<<"scSig "<<scSig<<endl;
  float scBck = hBack->Integral(hSignal->FindBin(0.2001), hSignal->FindBin(0.3001));
  hBack->Scale(scSig/scBck);

  hSignal->Add(hBack);

  TFile fout("extractSignal.root", "Update");
  hSignal->Write(Form("hSignal_%s_%f", name.Data(), pTLow));
  fout.Close();

  TF1 fgaus("fgaus", "gaus(0)", 0.09, 0.18);
  fgaus.SetParameters( 100, 0.13, 0.05);
  hSignal->Fit(&fgaus, "MQR");
  cout<<" fgaus.GetParameter(1): "<<fgaus.GetParameter(1)<<endl;
  fgaus.SetRange(fgaus.GetParameter(1) - 2*fgaus.GetParameter(2), fgaus.GetParameter(1) + 2*fgaus.GetParameter(2));
  hSignal->Fit(&fgaus, "MQR");
  return hSignal-> IntegralAndError(hSignal->FindBin(fgaus.GetParameter(1) - 3*fgaus.GetParameter(2)), hSignal->FindBin( fgaus.GetParameter(1) + 3*fgaus.GetParameter(2)), err);

}


void ExtractSignal_PCMEDC(){

  LoadTB();

  DrawSetMarkerTGraphErr(grTB_data, 21, 1.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_MC, 25, 2.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_Ratio, 21, 1.5, kBlue + 2, kBlue + 2, 2);
  DrawSetMarkerTGraphErr(grTB_Corr, 21, 1.5, kBlue + 2, kBlue + 2, 2);

  TFile fout("extractSignal_PCMEDC.root", "Recreate");
  fout.Close();


  //  stuff from train
  TFile fNCell2_data("/media/joshua/external_drive/Data/2020_11_12_pp13TeV_NCell/GammaConvCalo_2131_data_2018.root");
  TFile fNCell2_MC("/media/joshua/external_drive/Data/2020_11_12_pp13TeV_NCell/GammaConvCalo_2131_MC_2018.root");
  TFile NoNCell_data("/media/joshua/external_drive/Data/2020_11_12_pp13TeV_NCell/GammaConvCalo_2131_data_2018.root");
  TFile NoNCell_MC("/media/joshua/external_drive/Data/2020_11_12_pp13TeV_NCell/GammaConvCalo_2131_MC_2018.root");


  TFile fNellFile("/media/joshua/external_drive/Data/2021_03_18_pp13TeV/MC_complete/GammaCalo_2062_MC.root");


  TH2F* hInvMassVsPt_NCell2_data = (TH2F*) GetHistFromFile(&fNCell2_data, "ESD_Mother_InvMass_E_Calib", "00010113_0dm00009f9730000dge0404000_411790109fe32220000_0163103100000010");
  TH2F* hInvMassVsPt_NCell2_MC = (TH2F*) GetHistFromFile(&fNCell2_MC, "ESD_Mother_InvMass_E_Calib", "00010113_0dm00009f9730000dge0404000_411790109fe32220000_0163103100000010");
  TH2F* hInvMassVsPt_Back_NCell2_data = (TH2F*) GetHistFromFile(&fNCell2_data, "ESD_Background_InvMass_E_Calib", "00010113_0dm00009f9730000dge0404000_411790109fe32220000_0163103100000010");
  TH2F* hInvMassVsPt_Back_NCell2_MC = (TH2F*) GetHistFromFile(&fNCell2_MC, "ESD_Background_InvMass_E_Calib", "00010113_0dm00009f9730000dge0404000_411790109fe32220000_0163103100000010");


  TH2F* hInvMassVsPt_NCell1_data = (TH2F*) GetHistFromFile(&fNCell2_data, "ESD_Mother_InvMass_E_Calib", "00010113_0dm00009f9730000dge0404000_411790109fe30220000_0163103100000010");
  TH2F* hInvMassVsPt_NCell1_MC = (TH2F*) GetHistFromFile(&fNCell2_MC, "ESD_Mother_InvMass_E_Calib", "00010113_0dm00009f9730000dge0404000_411790109fe30220000_0163103100000010");
  TH2F* hInvMassVsPt_Back_NCell1_data = (TH2F*) GetHistFromFile(&fNCell2_data, "ESD_Background_InvMass_E_Calib", "00010113_0dm00009f9730000dge0404000_411790109fe30220000_0163103100000010");
  TH2F* hInvMassVsPt_Back_NCell1_MC = (TH2F*) GetHistFromFile(&fNCell2_MC, "ESD_Background_InvMass_E_Calib", "00010113_0dm00009f9730000dge0404000_411790109fe30220000_0163103100000010");


  TH2F* Cluster_NCell_VS_E = (TH2F*) GetHistFromFile(&fNellFile, "ClusterEnergyVsNCells_beforeQA 411792109fe30220000", "00010113_411792109fe30220000_0r631031000000d0");
  // return;



  const int numBins = 17;
  const double pTArr[numBins + 1] = {0.7, 0.8, 0.9, 1., 1.2,    1.4, 1.6, 1.8, 2.0, 2.5,
                                     3.0, 4.0, 6.0, 8.0, 10.,   12., 16, 20};
  TH1F* hRaw_NCell2_data =  new TH1F("hRaw_NCell2_data", "hRaw_NCell2_data", numBins, pTArr);
  TH1F* hRaw_NCell2_MC =  new TH1F("hRaw_NCell2_MC", "hRaw_NCell2_MC", numBins, pTArr);
  TH1F* hRaw_NCell1_data =  new TH1F("hRaw_NCell1_data", "hRaw_NCell1_data", numBins, pTArr);
  TH1F* hRaw_NCell1_MC =  new TH1F("hRaw_NCell1_MC", "hRaw_NCell1_MC", numBins, pTArr);

  double err;
  for(int i = 0; i < numBins - 1; ++i){

    hRaw_NCell2_data->SetBinContent(i+1, GetYield(hInvMassVsPt_NCell2_data, hInvMassVsPt_Back_NCell2_data, pTArr[i], pTArr[i+1], err, ""));
    hRaw_NCell2_data->SetBinError(i+1, err);

    hRaw_NCell2_MC->SetBinContent(i+1, GetYield(hInvMassVsPt_NCell2_MC, hInvMassVsPt_Back_NCell2_MC, pTArr[i], pTArr[i+1], err, ""));
    hRaw_NCell2_MC->SetBinError(i+1, err);

    hRaw_NCell1_data->SetBinContent(i+1, GetYield(hInvMassVsPt_NCell1_data, hInvMassVsPt_Back_NCell1_data, pTArr[i], pTArr[i+1], err, ""));
    hRaw_NCell1_data->SetBinError(i+1, err);

    hRaw_NCell1_MC->SetBinContent(i+1, GetYield(hInvMassVsPt_NCell1_MC, hInvMassVsPt_Back_NCell1_MC, pTArr[i], pTArr[i+1], err, ""));
    hRaw_NCell1_MC->SetBinError(i+1, err);
  }


  TH1F* hCorrected_NCell2 = (TH1F*) hRaw_NCell2_data->Clone("hCorrected_NCell2");
  hCorrected_NCell2->Divide(hRaw_NCell2_MC);

  TH1F* hCorrected_NCell1 = (TH1F*) hRaw_NCell1_data->Clone("hCorrected_NCell1");
  hCorrected_NCell1->Divide(hRaw_NCell1_MC);


  // build the double ratio
  TH1F* hCorrectionFact = (TH1F*) hCorrected_NCell2->Clone("hCorrectionFact");
  hCorrectionFact->Divide(hCorrected_NCell1);

  // Plotting
  StyleSettingsPaper();
  TCanvas *Can = new TCanvas("Can", "", 1200, 1000);
  DrawPaperCanvasSettings(Can, 0.1, 0.003, 0.003, 0.12);
  TH2F* hDummy    = new TH2F("hDummy","hDummy",1000,0.4, 20,1000,-0.2, 1.4);
  SetStyleHistoTH2ForGraphs(hDummy, "#it{E}_{#gamma, calo} (GeV)","ratio to NCell #geq 2", 0.85*0.044,0.044, 0.85*0.044,0.044, 1.1,1.2, 510, 510);
  Can->cd();
  hDummy->GetYaxis()->SetRangeUser(-0.2, 0.4);
  hDummy->Draw();


  DrawSetMarker(hCorrectionFact, 25, 2, kBlack, kBlack);

  DrawLines(0.5, 20, 0,0, 2, kGray+2, 2);
  for(int i = 1; i <= hCorrectionFact->GetNbinsX(); ++i) hCorrectionFact->SetBinContent(i, hCorrectionFact->GetBinContent(i) - 1);
  hCorrectionFact->Draw("same,pe");


  gSystem->Exec("mkdir Plots_PCMEDC_corr");
  Can->SaveAs("Plots_PCMEDC_corr/Correction.pdf");
  Can->SaveAs("Plots_PCMEDC_corr/Correction.png");


  TH1F* hNCells2_Vs_E = (TH1F*) Cluster_NCell_VS_E->ProjectionX("NCell2", 3, 20);
  TH1F* hNCells1_Vs_E = (TH1F*) Cluster_NCell_VS_E->ProjectionX("NCell1", 2, 2);

  TH1F* hNCells2_Vs_E_Reb = (TH1F*) hNCells2_Vs_E->Rebin(numBins, "hNCells2_Vs_E" ,pTArr);
  TH1F* hNCells1_Vs_E_Reb = (TH1F*) hNCells1_Vs_E->Rebin(numBins, "hNCells2_Vs_E" ,pTArr);
  // hNCells1_Vs_E->RebinAxis(19,pTArr);

  TH1F* hRho = (TH1F*) hNCells2_Vs_E_Reb->Clone("hRho");
  hRho->Multiply(hCorrectionFact);
  hRho->Divide(hNCells1_Vs_E_Reb);


  TFile fExternal("/home/joshua/PCG_Software/EMCal_NCellEffi/13TeVNomB_Wide/Pi0Tagging_13TeV_nom_03_25_NoTRD_1cellFT/pdf/histos.root");
  TH1F* hAllClus = (TH1F*) fExternal.Get("hNCell_AllClus_Effi_Corr");
  TH1F* hTaggedClus = (TH1F*) fExternal.Get("hNCell_Gammas_RW_SBSub_Effi_Corr");

  SetStyleHistoTH2ForGraphs(hDummy, "#it{E}_{#gamma, calo} (GeV)","#rho", 0.85*0.044,0.044, 0.85*0.044,0.044, 1.1,1.2, 510, 510);
  Can->cd();
  hDummy->GetYaxis()->SetRangeUser(-0.05, 1.4);
  hDummy->GetXaxis()->SetRangeUser(0., 8.);

  hDummy->Draw();

  DrawSetMarker(hRho, 33, 2.9, kGreen+2, kGreen+2);
  DrawSetMarker(hAllClus, 20, 2, kRed+2, kRed+2);
  DrawSetMarker(hTaggedClus, 27, 2.9, kOrange+2, kOrange+2);
  hRho->Draw("same");
  grTB_Corr->Draw("same,pe");
  hAllClus->Draw("same,pe");
  hTaggedClus->Draw("same,pe");

  TLegend *leg = GetAndSetLegend2(0.15, 0.92, 0.5, 0.92 - 0.05*3, 40);
  leg->AddEntry(hRho, "from PCM-EDC #pi^{0} spectrum", "p");
  leg->AddEntry(grTB_Corr, "TB", "p");
  leg->AddEntry(hAllClus, "All clus", "p");
  leg->AddEntry(hTaggedClus, "Tagged clus.", "p");
  leg->Draw("same");

  gSystem->Exec("mkdir Plots_PCMEDC_corr");
  Can->SaveAs("Plots_PCMEDC_corr/Rho.pdf");




}
