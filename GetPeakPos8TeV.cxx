


void GetPeakPos(){

  TFile *fdata = TFile::Open("Pi0Tagging_8TeV_01_28.root", "Update");

  TH2F* hdata = (TH2F*) fdata->Get("hInvMassVsPt_data");
  TH2F* hMC = (TH2F*) fdata->Get("hInvMassVsPt_MC");

  double xdata[100];
  double ydata[100];

  std::vector<float> vec = {1.4, 1.5, 1.6, 1.7, 1.8,   1.9, 2.0, 2.2, 2.4, 2.6,
                            2.8, 3.0, 3.2, 3.4, 3.6,   3.8, 4.0, 4.4, 4.8, 5.2,
                            5.6, 6.0, 6.5, 7.0, 7.5,   8.0, 8.5, 9.0, 9.5, 10.0,
                            12.0, 14.0, 20.0};

  for(int i = 0; i < vec.size() - 1; ++i){
    TH1D* hprojData = (TH1D*) hdata->ProjectionX(Form("hdata_%i", i), hdata->GetYaxis()->FindBin(vec.at(i) + 0.001), hdata->GetYaxis()->FindBin(vec.at(i+1) - 0.001));

    TF1* fgaus = new TF1("fgaus", "gaus(0)", 0.1, 0.15);
    fgaus->SetParLimits(1, 0.11, 0.15);
    hprojData->Fit(fgaus, "M0R");
    float mean = fgaus->GetParameter(1);
    float sig = fgaus->GetParameter(2);
    fgaus->SetRange(mean - 2*sig, mean + 2*sig);
    hprojData->Fit(fgaus, "M0R");
    xdata[i] =  (vec.at(i)+vec.at(i+1)) / 2.;
    ydata[i] =  fgaus->GetParameter(1);
    TCanvas Can("Can", "", 500, 500);
    Can.cd();
    hprojData->Draw();
    fgaus->Draw("same");
    Can.SaveAs(Form("DataBin%i.png", i));

  }

  TGraph grData(vec.size() - 1, xdata, ydata);


  double xMC[100];
  double yMC[100];
  for(int i = 0; i < vec.size() - 1; ++i){
    TH1D* hprojMC = (TH1D*) hMC->ProjectionX(Form("hMC_%i", i), hdata->GetYaxis()->FindBin(vec.at(i) + 0.001), hdata->GetYaxis()->FindBin(vec.at(i+1) - 0.001));

    TF1* fgaus = new TF1("fgaus", "gaus(0)", 0.1, 0.15);
    fgaus->SetParLimits(1, 0.11, 0.15);
    hprojMC->Fit(fgaus, "M0R");
    float mean = fgaus->GetParameter(1);
    float sig = fgaus->GetParameter(2);
    fgaus->SetRange(mean - 2*sig, mean + 2*sig);
    hprojMC->Fit(fgaus, "MR0");
    xMC[i] =  (vec.at(i)+vec.at(i+1)) / 2.;
    yMC[i] =  fgaus->GetParameter(1);
    TCanvas Can("Can", "", 500, 500);
    Can.cd();
    hprojMC->Draw();
    fgaus->Draw("same");
    Can.SaveAs(Form("MCBin%i.png", i));
  }

  TGraph grMC(vec.size() - 1, xMC, yMC);

  TGraph *grRatio = (TGraph*) grData.Clone("ratio");
  for(int i = 0; i < grMC.GetN(); ++i){
    grRatio->SetPointY(i, grData.GetPointY(i)/grMC.GetPointY(i));
  }

  grRatio->Print();


  TF1 funcData("funcData", "[0] + [1]*TMath::Power(x,[2])", 1.4,7);
  TF1 funcMC("funcMC", "[0] + [1]*TMath::Power(x,[2])", 1.4,7);
  grData.Fit(&funcData, "MQR0");
  grMC.Fit(&funcMC, "MQR0");

  TFile fMassPos("fMassPos.root", "Update");
  fMassPos.cd();
  grData.Write("MassPos8TeV_data", TObject::kOverwrite);
  grMC.Write("MassPos8TeV_MC", TObject::kOverwrite);
  funcData.Write("MassPos8TeV_func", TObject::kOverwrite);
  fMassPos.Close();



  TCanvas Can("Can", "", 1000, 1000);
  Can.cd();
  grData.Draw();
  funcData.SetLineColor(1);
  funcData.Draw("same");
  funcMC.Draw("same");
  grMC.SetLineColor(kRed);
  grMC.Draw("same");
  Can.SaveAs("PeakPos.png");

}
