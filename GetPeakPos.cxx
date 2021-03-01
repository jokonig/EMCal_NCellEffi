#include "Plotting_Header.h"


void GetPeakPos(TString period = "8TeV"){

  TString name = "Pi0Tagging_8TeV_02_24_ReverseTBNL.root";
  if(period.Contains("13TeVLowB")) name = "Pi0Tagging_13TeV_low_02_24_TBNL.root";
  if(period.Contains("13TeVNomB")) name = "Pi0Tagging_13TeV_nom_02_25_TBNL.root";
  TFile *fdata = TFile::Open(name, "Update");

  TH2F* hdata = (TH2F*) fdata->Get("hInvMassVsPt_data");
  TH2F* hMC = (TH2F*) fdata->Get("hInvMassVsPt_MC");

  double xdata[100];
  double ydata[100];
  double xdataErr[100];
  double ydataErr[100];

  std::vector<float> vec = {1.4, 1.5, 1.6, 1.7, 1.8,   1.9, 2.0, 2.2, 2.4, 2.6,
                            2.8, 3.0, 3.2, 3.4, 3.6,   3.8, 4.0, 4.4, 4.8, 5.2,
                            5.6, 6.0, 8.0, 10.0, 14.0, 20.0};

  if(period.Contains("13TeVLowB"))gSystem->Exec("mkdir PeakPos13TeVLowB");
  else if(period.Contains("13TeVNomB"))gSystem->Exec("mkdir PeakPos13TeVNomB");
  else gSystem->Exec("mkdir PeakPos8TeV");
  for(int i = 0; i < vec.size() - 1; ++i){
    TH1D* hprojData = (TH1D*) hdata->ProjectionX(Form("hdata_%i", i), hdata->GetYaxis()->FindBin(vec.at(i) + 0.001), hdata->GetYaxis()->FindBin(vec.at(i+1) - 0.001));

    TF1* fgaus = new TF1("fgaus", "gaus(0)", 0.1, 0.15);
    fgaus->SetParLimits(1, 0.11, 0.15);
    hprojData->Fit(fgaus, "M0R");
    float mean = fgaus->GetParameter(1);
    float sig = fgaus->GetParameter(2);
    fgaus->SetRange(mean - 2*sig, mean + 2*sig);
    hprojData->Fit(fgaus, "M0R");
    mean = fgaus->GetParameter(1);
    sig = fgaus->GetParameter(2);
    fgaus->SetRange(mean - sig, mean + sig);
    hprojData->Fit(fgaus, "M0R");
    xdata[i] =  (vec.at(i)+vec.at(i+1)) / 2.;
    xdataErr[i] =  0.001;
    ydata[i] =  fgaus->GetParameter(1);
    ydataErr[i] =  fgaus->GetParError(1);
    TCanvas Can("Can", "", 500, 500);
    Can.cd();
    hprojData->GetXaxis()->SetRangeUser(0., 0.25);
    hprojData->Draw();
    fgaus->Draw("same");
    drawLatexAdd(Form("%.01f < E < %.01f GeV", vec.at(i), vec.at(i+1)), 0.1, 0.95, 0.035);
    if(period.Contains("13TeVLowB"))Can.SaveAs(Form("PeakPos13TeVLowB/DataBin%i.png", i));
    else if(period.Contains("13TeVNomB"))Can.SaveAs(Form("PeakPos13TeVNomB/DataBin%i.png", i));
    else Can.SaveAs(Form("PeakPos8TeV/DataBin%i.png", i));

  }

  TGraphErrors grData(vec.size() - 1, xdata, ydata, xdataErr, ydataErr);


  double xMC[100];
  double yMC[100];
  double xMCErr[100];
  double yMCErr[100];
  for(int i = 0; i < vec.size() - 1; ++i){
    TH1D* hprojMC = (TH1D*) hMC->ProjectionX(Form("hMC_%i", i), hdata->GetYaxis()->FindBin(vec.at(i) + 0.001), hdata->GetYaxis()->FindBin(vec.at(i+1) - 0.001));

    TF1* fgaus = new TF1("fgaus", "gaus(0)", 0.1, 0.15);
    fgaus->SetParLimits(1, 0.11, 0.15);
    hprojMC->Fit(fgaus, "M0R");
    float mean = fgaus->GetParameter(1);
    float sig = fgaus->GetParameter(2);
    fgaus->SetRange(mean - 2*sig, mean + 2*sig);
    hprojMC->Fit(fgaus, "MR0");
    mean = fgaus->GetParameter(1);
    sig = fgaus->GetParameter(2);
    fgaus->SetRange(mean - sig, mean + sig);
    hprojMC->Fit(fgaus, "MR0");
    xMC[i] =  (vec.at(i)+vec.at(i+1)) / 2.;
    xMCErr[i] =  0.001;
    yMC[i] =  fgaus->GetParameter(1);
    yMCErr[i] =  fgaus->GetParError(1);
    TCanvas Can("Can", "", 500, 500);
    Can.cd();
    hprojMC->GetXaxis()->SetRangeUser(0., 0.25);
    hprojMC->Draw();

    fgaus->Draw("same");
    drawLatexAdd(Form("%.01f < E < %.01f GeV", vec.at(i), vec.at(i+1)), 0.1, 0.95, 0.035);
    if(period.Contains("13TeVLowB"))Can.SaveAs(Form("PeakPos13TeVLowB/MCBin%i.png", i));
    else if(period.Contains("13TeVNomB"))Can.SaveAs(Form("PeakPos13TeVNomB/MCBin%i.png", i));
    else Can.SaveAs(Form("PeakPos8TeV/MCBin%i.png", i));
  }

  TGraphErrors grMC(vec.size() - 1, xMC, yMC, xMCErr, yMCErr);

  TGraph *grRatio = (TGraph*) grData.Clone("ratio");
  for(int i = 0; i < grMC.GetN(); ++i){
    grRatio->SetPointY(i, grData.GetPointY(i)/grMC.GetPointY(i));
  }

  grRatio->Print();


  TF1 funcData("funcData", "[0] + [1]*TMath::Power(x,[2])", 1.6,6);
  funcData.SetParameters(1, 1, 1);
  TF1 funcMC("funcMC", "[0] + [1]*TMath::Power(x,[2])", 1.6,6);
  funcMC.SetParameters(1, 1, 1);
  grData.Fit(&funcData, "MQR0");
  grMC.Fit(&funcMC, "MQR0");
  //
  // TFile fMassPos("fMassPos.root", "Update");
  // fMassPos.cd();
  // if(period.Contains("13TeVLowB")){
  //   grData.Write("MassPos13TeVLowB_data", TObject::kOverwrite);
  //   grMC.Write("MassPos13TeVLowB_MC", TObject::kOverwrite);
  //   funcData.Write("MassPos13TeVLowB_func", TObject::kOverwrite);
  // } else if(period.Contains("13TeVNomB")){
  //   grData.Write("MassPos13TeVNomB_data", TObject::kOverwrite);
  //   grMC.Write("MassPos13TeVNomB_MC", TObject::kOverwrite);
  //   funcData.Write("MassPos13TeVNomB_func", TObject::kOverwrite);
  // } else if (period.Contains("8TeV")) {
  //   grData.Write("MassPos8TeV_data", TObject::kOverwrite);
  //   grMC.Write("MassPos8TeV_MC", TObject::kOverwrite);
  //   funcData.Write("MassPos8TeV_func", TObject::kOverwrite);
  // }
  // fMassPos.Close();



  TCanvas Can("Can", "", 1000, 1000);
  DrawPaperCanvasSettings(&Can, 0.14, 0.003, 0.003, 0.1);
  Can.cd();
  DrawSetMarkerTGraph(&grData, 20, 2, kBlack, kBlack);
  grData.SetTitle("");
  grData.GetYaxis()->SetTitle("M_{#gamma#gamma} (GeV/#it{c}^{2})");
  grData.GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
  grData.GetYaxis()->SetRangeUser(0.12, 0.15);
  DrawSetMarkerTGraph(&grMC, 24, 2, kRed, kRed);
  grData.GetXaxis()->SetRangeUser(0, 8);
  grData.Draw("AP");
  funcData.SetLineColor(1);
  funcData.Draw("same");
  funcMC.Draw("same");
  grMC.SetLineColor(kRed);
  grMC.Draw("same,p");

  TLegend *leg = GetAndSetLegend2(0.17, 0.78, 0.4, 0.9, 40);
  leg->AddEntry(&grData,"data", "p");
  leg->AddEntry(&grMC,"MC", "p");
  leg->Draw("same");

  if(period.Contains("13TeVLowB")) {
    drawLatexAdd("pp #sqrt{s} = 13 TeV (B=0.2T)", 0.17,0.92,0.044,kFALSE);
    Can.SaveAs(Form("PeakPos13TeVLowB/PeakPos%s.png", period.Data()));
  }
  else if(period.Contains("13TeVNomB")) {
    drawLatexAdd("pp #sqrt{s} = 13 TeV  (B=0.5T)", 0.17,0.92,0.044,kFALSE);
    Can.SaveAs(Form("PeakPos13TeVNomB/PeakPos%s.png", period.Data()));
  }
  else {
    drawLatexAdd("pp #sqrt{s} = 8 TeV", 0.17,0.92,0.044,kFALSE);
    Can.SaveAs(Form("PeakPos8TeV/PeakPosReverseNL%s.png", period.Data()));
  }

}
