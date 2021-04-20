#include "header.h"
#include "Plotting_Header.h"

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


void ExtractSignal(){
  TFile fout("extractSignal.root", "Recreate");
  fout.Close();


  //  stuff from train
  TFile fNCell2("/home/joshua/PCG_Software/AnalysisSoftware_Clones/NCellEffi_2021_01_11/00010113_411792106fe32220000_0r631031000000d0/13TeV/Pi0_data_GammaConvV1Correction_00010113_411792106fe32220000_0r631031000000d0.root");
  TFile NoNCell("/home/joshua/PCG_Software/AnalysisSoftware_Clones/NCellEffi_2021_01_11/00010113_411792106fe30220000_0r631031000000d0/13TeV/Pi0_data_GammaConvV1Correction_00010113_411792106fe30220000_0r631031000000d0.root");
  TH1D* hNCell2 = (TH1D*) fNCell2.Get("CorrectedYieldNormEff");
  hNCell2->Rebin(2);
  TH1D* hNoNCell = (TH1D*) NoNCell.Get("CorrectedYieldNormEff");
  hNoNCell->Rebin(2);
  TH1D* hNoNCell_Ratio = (TH1D*) hNCell2->Clone("hNoNCell_Ratio");
  hNoNCell_Ratio->Divide(hNCell2, hNoNCell, 1, 1, "B");


  TFile *f = TFile::Open("InvMassCorrected_13TeV_nom_03_25_TBNL.root");

  std::vector<TH2F*> vecSignal;
  std::vector<TH2F*> vecBack;

  int nMethods = 10;
  // NCell >= 2 cells

  // 0 -> normal NCell >=2 cut
  // 1 -> No NCell cut (useless here)
  // 2 -> all cluster
  // 3 -> TB
  // 4 -> Gamma and elec applied on all clusters
  // 5 -> Gamma and elec but seperated ()TB for gamma and elec for elec
  // 6 ->Gamma and elec but seperated ()GammaAndElec for gamma and elec for elec
  for(int i = 0; i < nMethods; ++i){
    vecSignal.push_back( (TH2F*) f->Get(Form("hInvMassVsPt_%i_MC", i)));
    vecSignal.at(vecSignal.size() - 1)->Sumw2();
    vecBack.push_back( (TH2F*) f->Get(Form("hInvMassVsPtBack_%i_MC", i)));
    vecBack.at(vecBack.size() - 1)->Sumw2();
  }


  const int numBins = 14;
  const double pTArr[numBins + 1] = {1.4, 1.6, 1.8, 2.0, 2.5,     3.0, 3.5, 4.0, 5.0, 6.0,
                                     8.0, 10., 12., 16, 20};
  std::vector<TH1F*> vecHistYields;
  vecHistYields.clear();
  vecHistYields.push_back( new TH1F("hRawYieldNCell2", "hRawYieldNCell2", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNoNCell", "hRawYieldNoNCell", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNCEAllClus", "hRawYieldNCEAllClus", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNCETB", "hRawYieldNCETB", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNCEGammaAndElecGaus", "hRawYieldNCEGammaAndElecGaus", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNCETBAndElec", "hRawYieldNCETBAndElec", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNCEGammaAndElec", "hRawYieldNCEGammaAndElec", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNCEGammaAndElec2", "hRawYieldNCEGammaAndElecAll2", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNCEGammaAndElec3", "hRawYieldNCEGammaAndElecAll3", numBins, pTArr));
  vecHistYields.push_back( new TH1F("hRawYieldNCEGammaLin", "hRawYieldNCEGammaLin", numBins, pTArr));


  for(int i = 0; i < numBins - 1; ++i){
    for(int j = 0; j < nMethods; ++j){
      double err = 0;
      float yield = GetYield(vecSignal.at(j), vecBack.at(j), pTArr[i], pTArr[i+1], err);
      vecHistYields.at(j)->SetBinContent(i+1, yield);
      vecHistYields.at(j)->SetBinError(i+1, err);
    }

  }

  // make ratios to std ncell = 2
  for(int i = 1; i < vecHistYields.size(); ++i){
    vecHistYields.at(i)->Divide(vecHistYields.at(i), vecHistYields.at(0));
    // vecHistYields.at(i)->Divide(vecHistYields.at(i), vecHistYields.at(0), 1, 1, "B");
  }


  // Plotting
  StyleSettingsPaper();
  TCanvas *Can = new TCanvas("Can", "", 1200, 1000);
  DrawPaperCanvasSettings(Can, 0.1, 0.003, 0.003, 0.12);
  TH2F* hDummy    = new TH2F("hDummy","hDummy",1000,0.4, 20,1000,0.9, 1.4);
  SetStyleHistoTH2ForGraphs(hDummy, "#it{p}_{T} (GeV/#it{c})","ratio to NCell #geq 2", 0.85*0.044,0.044, 0.85*0.044,0.044, 1.1,1.2, 510, 510);
  Can->cd();
  hDummy->Draw();

  DrawSetMarker(hNoNCell_Ratio, 25, 2, kBlack, kBlack);
  // DrawSetMarker(vecHistYields.at(1), 24, 3, kBlue + 2, kBlue + 2);
  DrawSetMarker(vecHistYields.at(2), 20, 2, kRed + 2, kRed + 2);
  DrawSetMarker(vecHistYields.at(3), 24, 2, kBlue + 2, kBlue + 2);
  DrawSetMarker(vecHistYields.at(9), 28, 3, kOrange + 2, kOrange + 2);
  // DrawSetMarker(vecHistYields.at(5), 34, 2, kBlack, kBlack);
  // DrawSetMarker(vecHistYields.at(9), 34, 2, kCyan + 2, kCyan + 2);

  DrawLines(0.5, 20, 1,1, 2, kGray+2, 2);
  // vecHistYields.at(1)->Draw("same,pe");
  hNoNCell_Ratio->Draw("same,pe");
  vecHistYields.at(2)->Draw("same,pe");
  vecHistYields.at(3)->Draw("same,pe");
  vecHistYields.at(9)->Draw("same,pe");
  // vecHistYields.at(5)->Draw("same,pe");
  // vecHistYields.at(6)->Draw("same,pe");

  TLegend *leg = GetAndSetLegend2(0.3, 0.92, 0.95, 0.92 - 0.05*3, 40);
  leg->AddEntry(hNoNCell_Ratio, "no N_{cell} cut (from train)", "p");
  leg->AddEntry(vecHistYields.at(2), "N_{cell} func: All clus", "p");
  leg->AddEntry(vecHistYields.at(3), "N_{cell} func: TB", "p");
  leg->AddEntry(vecHistYields.at(9), "N_{cell} func: linear fit to tagged #gamma", "p");
  // leg->AddEntry(vecHistYields.at(5), "N_{cell} func: TB for #gamma and elec for e^{#pm}", "p");
  // leg->AddEntry(vecHistYields.at(6), "N_{cell} func: GammaAndElec3", "p");
  leg->Draw("same");

  gSystem->Exec("mkdir Plots_Spectra");
  Can->SaveAs("Plots_Spectra/Comparison.pdf");
  Can->SaveAs("Plots_Spectra/Comparison.png");


}
