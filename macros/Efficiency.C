#include <string>
#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"

using namespace std;

TCanvas *canv;
TString outname = "prova.root";
bool first;
/////////////////////////

// USE: root -b -L effRej.C

/////////////////////////

TGraphErrors* EffRej(TString file1, TString file2, TString la1, TString la2, TString nplot1, TString nplot2, TString xname, int reb =1, TString fix="", float CutValue = -999)
{
  TFile f1( file1, "READ");
  //TFile f2( file2, "READ");

  TString n1 = nplot1;
  TString n2 = nplot2;
  TH1D* h1= (TH1D*) f1.Get( n1 )->Clone("h1");
  TH1D* h2= (TH1D*) f1.Get( n2 )->Clone("h2");
  h1->Rebin(reb);
  h2->Rebin(reb);
//  h1->SetAxisRange(0,200);
//  h2->SetAxisRange(0,200);

    
  TLegend *l = new TLegend(0.45,0.60, 0.85,0.85);
  //  TLegend *l = new TLegend(0.15,0.60, 0.45,0.85);
  l->AddEntry(h1, la1, "lp");
  l->AddEntry(h2, la2, "lp");
  l->SetShadowColor(0);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->SetTextFont(22);
  l->SetTextSize(0.04);

  canv = new TCanvas();
  canv->cd();
  gStyle->SetOptStat(0);

  h1->SetMarkerStyle(23);
  h2->SetMarkerStyle(22);
  h2->SetLineColor(2);
  h2->SetMarkerColor(2);
  h1->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->SetTitle(xname);
  
  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  
  float max = (h1->GetMaximum() >  h2->GetMaximum()) ? h1->GetMaximum() : h2->GetMaximum();
  h1->GetYaxis()->SetRangeUser(0, max*1.3);

  if ( xname.Contains("ntrk") )h1->GetXaxis()->SetRangeUser(0,20);
  if ( xname.Contains("Mass") )h1->GetXaxis()->SetRangeUser(0,7);
  h1->SetTitle("");
  h2->SetTitle("");
  //  h1->DrawNormalized();
  //  h2->DrawNormalized("same");
  h1->Draw();
  h2->Draw("same");
  l->Draw("same");
  nplot1.ReplaceAll(" ","_");
  canv->Print("eff/" + nplot1 + fix + ".pdf"); 
  
  // Now we can work on the efficiency vs rejection plot
 
  //  TH1D * heff = h1->Clone();
  TGraphErrors *effplot = new TGraphErrors;
  int maxbin = h1->GetNbinsX();
  int startBin = 0;
  if ( xname.Contains("MET") ) startBin = h1->FindBin(25.);
  if ( xname.Contains("WmT") ) startBin = h1->FindBin(40.);
  

  int cutBin = -1;
  if (CutValue != 999) cutBin = h1->FindBin(CutValue);

  for (int i = 1; i< maxbin-1; i++) {
    //    if ( xname.Contains("MET") && i<startBin ) continue;
    //    if ( xname.Contains("WmT") && i<startBin ) continue;
    if (i < startBin) continue;
    float sigeff = h1->Integral(startBin, i) / h1->Integral(startBin,maxbin);
    float bkgeff = h2->Integral(startBin, i) / h2->Integral(startBin,maxbin);
    if (sigeff > 0.001 && sigeff < 0.999 && bkgeff > 0.001 && bkgeff < 0.999) {
      effplot->SetPoint(i, bkgeff, sigeff);
      if ( i == cutBin && cutBin != -1) effplot->SetPointError(i, 0, 1);
    }
  }
  effplot->GetXaxis()->SetTitle("Fake rate");
  effplot->GetYaxis()->SetTitle("Efficiency");
  effplot->GetXaxis()->SetRangeUser(0,1);
  effplot->GetYaxis()->SetRangeUser(0,1);
  //  effplot->SetAxisRange(0,1);
  effplot->SetMarkerSize(5);
  effplot->SetMarkerStyle(7);
  canv->SetGridx();
  canv->SetGridy();
  effplot->Draw("AP");
  canv->Print("eff/" + nplot1 + fix + "Eff.pdf");
  //  delete canv;
  //  delete l   ;
  return effplot;
}

void EfficiencyIsoCut()
{
  TString filename1 = "../looper/output/out_ZeePlusTT536_40to60_BasicLoop.root";
  TString filename2 = "../looper/output/out_ZeePlusTT536_40to60_AODarea.root";	       
  TString filename3 = "../looper/output/out_ZeePlusTT536_40to60_DeltaBetaSimple.root";
  TString filename4 = "../looper/output/out_ZeePlusTT536_40to60_DeltaBetaWeights.root";
  TString filename5 = "../looper/output/out_ZeePlusTT536_40to60_DeltaBetaWeights01.root";
  TString label1    = "True";
  TString label2    = "Fake";
  TString EE = "";
  //void EffRej(TString file1, TString file2, TString la1, TString la2, TString nplot1, TString nplot2, TString xname, int reb =1, TString fix="")
  //  EffRej(filename1, filename2 , label1, label2, "Muon 1JetIncl 1 B-tag MET C4" , "MET",1);
  TGraphErrors* h_1 = EffRej(filename1, filename1 , label1, label2, "h_true"+EE+"_el_iso_cor" ,"h_fake"+EE+"_el_iso_cor" , "RelIso",1, "uncorr"); 
  TGraphErrors* h_2 = EffRej(filename2, filename2 , label1, label2, "h_true"+EE+"_el_iso_cor" ,"h_fake"+EE+"_el_iso_cor" , "RelIso",1, "area");
  TGraphErrors* h_3 = EffRej(filename3, filename3 , label1, label2, "h_true"+EE+"_el_iso_cor" ,"h_fake"+EE+"_el_iso_cor" , "RelIso",1, "DB");
  TGraphErrors* h_4 = EffRej(filename4, filename4 , label1, label2, "h_true"+EE+"_el_iso_cor" ,"h_fake"+EE+"_el_iso_cor" , "RelIso",1, "DBW");

  canv = new TCanvas();
  canv->cd();
  gStyle->SetOptStat(0);
  //  h_1->GetYaxis()->SetRangeUser(0.4,1);
  h_1->GetXaxis()->SetLimits(0.0,0.6);


  h_1->SetMarkerColor(kRed);
  h_2->SetMarkerColor(kBlue);
  h_3->SetMarkerColor(kBlack);
  h_4->SetMarkerColor(kGreen);

  h_1->SetLineColor(kRed);
  h_2->SetLineColor(kBlue);
  h_3->SetLineColor(kBlack);
  h_4->SetLineColor(kGreen);

  h_1->SetMarkerStyle(20);
  h_2->SetMarkerStyle(21);
  h_3->SetMarkerStyle(22);
  h_4->SetMarkerStyle(23);

  h_1->SetMarkerSize(1);
  h_2->SetMarkerSize(1);
  h_3->SetMarkerSize(1);
  h_4->SetMarkerSize(1);
  

  TLegend* leg = new TLegend(0.32,0.20,0.70,0.35);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->AddEntry(h_1,"PFCand loop, uncorr","PL");
  leg->AddEntry(h_2,"PFCand loop, area corr","PL");
  leg->AddEntry(h_3,"PFCand loop, simple DeltaBeta","PL");
  leg->AddEntry(h_4,"PFCand loop, DeltaBetaWeights","PL");

  h_1->Draw("AP");
  h_2->Draw("Psame");
  h_3->Draw("Psame");
  h_4->Draw("Psame");
  leg->Draw();
  canv->Print("eff/CompareEff"+EE+"_4060.pdf");



  return;
}

void EfficiencySieieCut()
{
  TString filename1 = "../looper/output/out_ZeePlusTT700.root";
  TString filename2 = "../looper/output/out_ZeePlusTT700patch.root";	       
//  TString filename3 = "../looper/output/out_ZeePlusTT536_40to60_DeltaBetaSimple.root";
//  TString filename4 = "../looper/output/out_ZeePlusTT536_40to60_DeltaBetaWeights.root";
//  TString filename5 = "../looper/output/out_ZeePlusTT536_40to60_DeltaBetaWeights01.root";
  TString label1    = "True";
  TString label2    = "Fake";
  TString EE = "";
  TString cut = "M_";
  //TString variable = "sieie";
  //TString variable = "dEtaIn";
  //TString variable = "dPhiIn";
  //TString variable = "hOverE";
  //TString variable = "d0corr";
  TString variable = "z0corr";
  //TString variable = "ooemoop";
  float cutValue = -999;
  if (EE=="EE") cutValue = -999;
  //void EffRej(TString file1, TString file2, TString la1, TString la2, TString nplot1, TString nplot2, TString xname, int reb =1, TString fix="")
  //  EffRej(filename1, filename2 , label1, label2, "Muon 1JetIncl 1 B-tag MET C4" , "MET",1);
  TGraphErrors* h_1 = EffRej(filename1, filename1 , label1, label2, "h_"+cut+"true"+EE+"_el_"+variable ,"h_"+cut+"fake"+EE+"_el_"+variable , variable,1, "rel7"     , cutValue); 
  TGraphErrors* h_2 = EffRej(filename2, filename2 , label1, label2, "h_"+cut+"true"+EE+"_el_"+variable ,"h_"+cut+"fake"+EE+"_el_"+variable , variable,1, "rel7patch", cutValue);
//  TGraphErrors* h_3 = EffRej(filename3, filename3 , label1, label2, "h_true"+EE+"_el_iso_cor" ,"h_fake"+EE+"_el_iso_cor" , "RelIso",1, "DB");
//  TGraphErrors* h_4 = EffRej(filename4, filename4 , label1, label2, "h_true"+EE+"_el_iso_cor" ,"h_fake"+EE+"_el_iso_cor" , "RelIso",1, "DBW");

  canv = new TCanvas();
  canv->cd();
  gStyle->SetOptStat(0);
  //  h_1->GetYaxis()->SetRangeUser(0.4,1);



  h_1->SetMarkerColor(kRed);
  h_2->SetMarkerColor(kBlue);
//  h_3->SetMarkerColor(kBlack);
//  h_4->SetMarkerColor(kGreen);

  h_1->SetLineColor(kRed);
  h_2->SetLineColor(kBlue);
//  h_3->SetLineColor(kBlack);
//  h_4->SetLineColor(kGreen);

  h_1->SetMarkerStyle(20);
  h_2->SetMarkerStyle(21);
//  h_3->SetMarkerStyle(22);
//  h_4->SetMarkerStyle(23);

  h_1->SetMarkerSize(1);
  h_2->SetMarkerSize(1);
//  h_3->SetMarkerSize(1);
//  h_4->SetMarkerSize(1);
  

  TLegend* leg = new TLegend(0.62,0.20,0.90,0.35);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->AddEntry(h_1,"Rel7","PL");
  leg->AddEntry(h_2,"Rel7 5x5","PL");
//  leg->AddEntry(h_3,"PFCand loop, simple DeltaBeta","PL");
//  leg->AddEntry(h_4,"PFCand loop, DeltaBetaWeights","PL");

  h_1->Draw("AP");
  h_1->GetXaxis()->SetLimits(0,1);
  h_1->GetXaxis()->SetRangeUser(0,1);
  h_1->Draw("AP");
  h_2->Draw("Psame");

//  h_3->Draw("Psame");
//  h_4->Draw("Psame");
  leg->Draw();
  canv->Print("eff/Compare_"+variable+"Eff_"+cut+EE+".pdf");



  return;
}

void Efficiency()
{
  //  EfficiencyIsoCut();
  EfficiencySieieCut();
}
