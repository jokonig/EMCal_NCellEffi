#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TColor.h"
#include "TFile.h"
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <string>

//__________________________________________________________________________________________________________

void SetPlotStyle() {
	const Int_t nRGBs = 5;
	const Int_t nCont = 255;

	Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[nRGBs] = { 0.31, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[nRGBs]  = { 0.51, 1., 0.12, 0.00, 0.00};

	TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);
	gStyle->SetNumberContours(nCont);
}



void StyleSettingsPaper( TString format = ""){
    //gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);   //show day and time
    gStyle->SetOptStat(0);  //show statistic
    gStyle->SetPalette(1,0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetTitleFillColor(0);
    gStyle->SetTextSize(0.5);
    gStyle->SetLabelSize(0.03,"xyz");
    gStyle->SetLabelOffset(0.002,"xyz");
    gStyle->SetTitleFontSize(0.04);
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(0.7,"x");
    gStyle->SetCanvasColor(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(1);

    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.09);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadLeftMargin(0.13);


    TGaxis::SetMaxDigits(3);
    gErrorIgnoreLevel=kError;

    if (format.CompareTo("eps") == 0 ||format.CompareTo("pdf") == 0  ) gStyle->SetLineScalePS(1);
    SetPlotStyle();
}
//
//  TCanvas *canvasExample  = new TCanvas("canvasExample","",0,0,1300,850);
//  Example:     DrawPaperCanvasSettings( canvasExample, 0.085, 0.01, 0.01, 0.105);
//__________________________________________________________________________________________________________
void DrawPaperCanvasSettings(
    TCanvas* c1,
    Double_t leftMargin,
    Double_t rightMargin,
    Double_t topMargin,
    Double_t bottomMargin
){
    c1->SetTickx();
    c1->SetTicky();
    c1->SetGridx(0);
    c1->SetGridy(0);
    c1->SetLogy(0);
    c1->SetLeftMargin(leftMargin);
    c1->SetRightMargin(rightMargin);
    c1->SetTopMargin(topMargin);
    c1->SetBottomMargin(bottomMargin);
    c1->SetFillColor(0);
}
//
//     Example
//     Double_t minY                                 = 0.91;
//     Double_t maxY                                 = 1.039;
//     Double_t minX                                 = 0.27;
//     Double_t maxX                                 = 2.99e2;
//     Double_t textSizeSinglePad                    = 0.05;
//     Double_t textSizeLabelsPixel                  = 35;
//     Double_t textSizeLabelsRel                    = 35./canvasheight;
//
// it is often best to use a dummy histogram for the style settings
//     TH2F * histExampleDummy    = new TH2F("histExampleDummy","histExampleDummy",1000,minX, maxX,1000,minY, maxY);
// set the style of the dummy (dummy, x title, y title, x label size, x title size, y label size, y title size, x title offset, y title offset, ndivisions x, ndivisions y)
//     SetStyleHistoTH2ForGraphs(histExampleDummy, "#it{E}_{rec} (GeV)","#it{E}_{rec}/#it{E}_{in}", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81, 510, 510);
//
//   histExampleDummy->GetXaxis()->SetNoExponent(); // when SetLogx() is used, one can think about adding more labels
//   histExampleDummy->GetXaxis()->SetMoreLogLabels(kTRUE);
//   histExampleDummy->DrawCopy();
//
//__________________________________________________________________________________________________________
void SetStyleHistoTH2ForGraphs(
    TH2* histo,
    TString XTitle,
    TString YTitle,
    Size_t xLableSize,
    Size_t xTitleSize,
    Size_t yLableSize,
    Size_t yTitleSize,
    Float_t xTitleOffset    = 1,
    Float_t yTitleOffset    = 1,
    Int_t xNDivisions       = 510,
    Int_t yNDivisions       = 510,
    Font_t textFontLabel    = 42,
    Font_t textFontTitle    = 62
){
    histo->SetXTitle(XTitle);
    histo->SetYTitle(YTitle);
    histo->SetTitle("");

    histo->GetXaxis()->SetLabelFont(textFontLabel);
    histo->GetYaxis()->SetLabelFont(textFontLabel);
    histo->GetXaxis()->SetTitleFont(textFontTitle);
    histo->GetYaxis()->SetTitleFont(textFontTitle);

    histo->GetXaxis()->SetLabelSize(xLableSize);
    histo->GetXaxis()->SetTitleSize(xTitleSize);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(yLableSize);
    histo->GetYaxis()->SetTitleSize(yTitleSize);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);
}

//
//  draw a line to guide the eye (good for some comparisons)
//  min x, max x, min y, max y, line width, line color, line style)
//  DrawLines(minX, maxX, 1., 1., 1, kGray+2, 7);
//
//__________________________________________________________________________________________________________
void DrawLines(
    Float_t startX, Float_t endX,
    Float_t startY, Float_t endY,
    Float_t linew, Float_t lineColor = 4, Style_t lineStyle = 1
){
    TLine * l1 = new TLine (startX,startY,endX,endY);
    l1->SetLineColor(lineColor);
    l1->SetLineWidth(linew);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
}
// function to set the plotting style (takes graph, markerstyle, markersize, and two times the color (for marker and lines))
// Example
//
//  Color_t  colorData          = kBlack;
//  Style_t  markerStyleData    = 20; // do not use asymmetric markers (like triangles [numbers 23, 26, 32 are forbidden!])
//  Size_t   markerSize         = 1.5; //1.5 or 2 make reasonable sizes on most paper plots
//   DrawSetMarker(  histoData,   30,               markerSize,     kOrange-8,     kOrange-8);
//
//__________________________________________________________________________________________________________
void DrawSetMarker(
    TH1* histo1,
    Style_t markerStyle,
    Size_t markerSize,
    Color_t markerColor,
    Color_t lineColor
) {
    histo1->SetMarkerStyle(markerStyle);
    histo1->SetMarkerSize(markerSize);
    histo1->SetLineWidth(markerSize);
    histo1->SetMarkerColor(markerColor);
    histo1->SetLineColor(lineColor);
    if(markerStyle < 10)histo1->SetLineStyle(markerStyle);
    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
}
//
// function to set the plotting style (takes graph, markerstyle, markersize, and two times the color (for marker and lines))
//  Color_t  colorData          = kBlack;
//  Style_t  markerStyleData    = 20; // do not use asymmetric markers (like triangles [numbers 23, 26, 32 are forbidden!])
//  Size_t   markerSize         = 1.5; //1.5 or 2 make reasonable sizes on most paper plots
//
//  DrawSetMarkerTGraphErr(graphData,        markerStyleData,  markerSize,     colorData,     colorData);
//
//__________________________________________________________________________________________________________
void DrawSetMarkerTGraphErr(
    TGraphErrors* graph,
    Style_t markerStyle,
    Size_t markerSize,
    Color_t markerColor,
    Color_t lineColor,
    Width_t lineWidth       = 1,
    Bool_t boxes            = kFALSE,
    Color_t fillColor       = 0,
    Bool_t isHollow         = kFALSE
) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            if (!isHollow){
                graph->SetFillStyle(1001);
            } else {
                graph->SetFillStyle(0);
            }
        } else {
            graph->SetFillStyle(0);
        }
    }
}
//__________________________________________________________________________________________________________
void DrawSetMarkerTGraph(
    TGraph* graph,
    Style_t markerStyle,
    Size_t markerSize,
    Color_t markerColor,
    Color_t lineColor,
    Width_t lineWidth       = 1,
    Bool_t boxes            = kFALSE,
    Color_t fillColor       = 0,
    Bool_t isHollow         = kFALSE
) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            if (!isHollow){
                graph->SetFillStyle(1001);
            } else {
                graph->SetFillStyle(0);
            }
        } else {
            graph->SetFillStyle(0);
        }
    }
}


//__________________________________________________________________________________________________________
void DrawSetMarkerTGraphAsym(
    TGraphAsymmErrors* graph,
    Style_t markerStyle,
    Size_t markerSize,
    Color_t markerColor,
    Color_t lineColor,
    Width_t lineWidth   =1,
    Bool_t boxes        = kFALSE,
    Color_t fillColor   = 0,
    Bool_t isHollow     = kFALSE
) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(markerSize);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
    if (boxes){
        graph->SetFillColor(fillColor);
        if (fillColor!=0){
            if (!isHollow){
                graph->SetFillStyle(1001);
            } else {
                graph->SetFillStyle(0);
            }
        } else {
            graph->SetFillStyle(0);
        }
    }
}




//__________________________________________________________________________________________________________
void SetStyleTLatex(
    TLatex* text,
    Size_t textSize,
    Width_t lineWidth,
    Color_t textColor = 1,
    Font_t textFont = 42,
    Bool_t kNDC = kTRUE,
    Short_t align = 11
){
    if (kNDC) {text->SetNDC();}
    text->SetTextFont(textFont);
    text->SetTextColor(textColor);
    text->SetTextSize(textSize);
    text->SetLineWidth(lineWidth);
    text->SetTextAlign(align);
}


//
//     Example
//     Double_t textSizeLabelsRel    = 35./canvasheight;
//     drawLatexAdd("ALICE Performance",0.15,0.92,textSizeLabelsRel,kFALSE);
//     drawLatexAdd("Xe-Xe, #sqrt{#it{s}_{NN}} = 90 TeV",0.15,0.92-textSizeLabelsRel,textSizeLabelsRel,kFALSE);
//     drawLatexAdd("e^{#pm} rec. with EMCal",0.15,0.92-2*textSizeLabelsRel,textSizeLabelsRel,kFALSE);
//
//__________________________________________________________________________________________________________

void drawLatexAdd(TString latextext, Double_t textcolumn, Double_t textrow, Double_t textSizePixel,Bool_t setFont = kFALSE, Bool_t setFont2 = kFALSE, Bool_t alignRight = kFALSE, Color_t textcolor = kBlack){
    TLatex *latexDummy                  = new TLatex(textcolumn ,textrow,latextext);
    SetStyleTLatex( latexDummy, textSizePixel,4);
    if(setFont)
        latexDummy->SetTextFont(62);
    if(setFont2)
        latexDummy->SetTextFont(43);
    if(alignRight)
        latexDummy->SetTextAlign(31);
    latexDummy->SetTextColor(textcolor);
    latexDummy->Draw();
}



//__________________________________________________________________________________________________________
void DrawGammaSetMarkerTF1(
    TF1* fit1,
    Style_t lineStyle,
    Size_t lineWidth,
    Color_t lineColor
) {
    fit1->SetLineColor(lineColor);
    fit1->SetLineStyle(lineStyle);
    fit1->SetLineWidth(lineWidth);
}

//__________________________________________________________________________________________________________

TLegend *GetAndSetLegend2(
    Double_t positionX,
    Double_t positionY,
    Double_t positionXRight,
    Double_t positionYUp,
    Size_t textSize,
    Int_t columns               = 1,
    TString header              = "",
    Font_t textFont             = 43,
    Double_t margin             = 0
){
    TLegend *legend = new TLegend(positionX,positionY,positionXRight,positionYUp);
    legend->SetNColumns(columns);
    legend->SetLineColor(0);
    legend->SetLineWidth(0);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetLineStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(textFont);
    legend->SetTextSize(textSize);
    if (margin != 0) legend->SetMargin(margin);
    if (header.CompareTo("")!= 0) legend->SetHeader(header);
    return legend;
}




struct LegLeft{

	LegLeft(double xmin, double ymin, double xmax, double ymax, Size_t size){
		leg = GetAndSetLegend2(xmin, ymin, xmax, ymax, size);
		leg->SetNColumns(2);
		leg->SetTextAlign(32);
	}

	TLegend* Leg()			{return leg;};

	void AddEntry(TObject* h, TString s = "", TString p = ""){
		leg->AddEntry((TObject*)0,s,"");
		leg->AddEntry(h,"#color[0]{i}",p);
	}

	void Draw(TString s = ""){
		leg->Draw(s);
	}

	TLegend *leg = nullptr;

};
