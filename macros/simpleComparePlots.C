
//[cms-tas01] /tas/gzevi/ElectronID/ElectronLooperCMS2/macros > root -l 'simpleComparePlots.C("h_nvtx_eff", "Isolation Efficiency vs NPV")' 

void simpleComparePlots(TString histoName, TString title){

  gStyle->SetOptStat(0);

  //  TFile* f_00eta09 = TFile::Open(Form("trktest_hp/histo_algo%02d_00eta09.root",algo));


  TFile* f_1 = TFile::Open("../looper/output/out_ZeePlusTT536_10to20_BasicLoop.root");
  TFile* f_2 = TFile::Open("../looper/output/out_ZeePlusTT536_10to20_AODarea.root");
  TFile* f_3 = TFile::Open("../looper/output/out_ZeePlusTT536_10to20_DeltaBetaSimple.root");
  TFile* f_4 = TFile::Open("../looper/output/out_ZeePlusTT536_10to20_DeltaBetaWeights.root");
  TFile* f_5 = TFile::Open("../looper/output/out_ZeePlusTT536_10to20_DeltaBetaWeights01.root");


  TH1F* h_1 = (TH1F*) f_1->Get(histoName);
  TH1F* h_2 = (TH1F*) f_2->Get(histoName);
  TH1F* h_3 = (TH1F*) f_3->Get(histoName);
  TH1F* h_4 = (TH1F*) f_4->Get(histoName);
  TH1F* h_5 = (TH1F*) f_5->Get(histoName);

  h_1->SetMarkerColor(kRed);
  h_2->SetMarkerColor(kBlue);
  h_3->SetMarkerColor(kBlack);
  h_4->SetMarkerColor(kGreen);
  h_5->SetMarkerColor(kViolet);

  h_1->SetLineColor(kRed);
  h_2->SetLineColor(kBlue);
  h_3->SetLineColor(kBlack);
  h_4->SetLineColor(kGreen);
  h_5->SetLineColor(kViolet);

  h_1->SetMarkerStyle(20);
  h_2->SetMarkerStyle(21);
  h_3->SetMarkerStyle(22);
  h_4->SetMarkerStyle(23);
  h_5->SetMarkerStyle(24);

  TCanvas c;


//  h_00eta09->Scale(1./h_00eta09->Integral());
//  h_09eta16->Scale(1./h_09eta16->Integral());
//  h_16eta25->Scale(1./h_16eta25->Integral());

  h_1->SetTitle(title);
  h_1->GetXaxis()->SetTitle(histoName);
  h_1->GetYaxis()->SetTitle("");
  const float minA[5] = {h_1->GetMinimum(), h_2->GetMinimum(), h_3->GetMinimum(), h_4->GetMinimum(), h_5->GetMinimum()};
  const float maxA[5] = {h_1->GetMaximum(), h_2->GetMaximum(), h_3->GetMaximum(), h_4->GetMaximum(), h_5->GetMaximum()};
  float min = TMath::MinElement(5, minA);
  float max = TMath::MaxElement(5, maxA);
  h_1->GetYaxis()->SetRangeUser(min*0.8, max*1.2);

  if ( !histoName.Contains("_eff") && !histoName.Contains("_fr") )  {

    //    c.SetLogy();
    h_1->GetYaxis()->SetTickLength(0.01);
    int rebin = 1;
    if (histoName.Contains("_fake")) rebin = 2;
    h_1->Rebin(rebin);
    h_2->Rebin(rebin);
    h_3->Rebin(rebin);
    h_4->Rebin(rebin);
    h_5->Rebin(rebin);
    if (min < 1) h_1->GetYaxis()->SetRangeUser(1, rebin*max);
    if ( histoName.Contains("_true")) h_1->GetXaxis()->SetRangeUser(0, 1);
    h_1->Draw("hP");
    h_2->Draw("samehP");
    h_3->Draw("samehP");
    h_4->Draw("samehP");
    //h_5->Draw("samehP");
  }
  else {
    h_1->Draw("P");
    h_2->Draw("sameP");
    h_3->Draw("sameP");
    h_4->Draw("sameP");
    h_5->Draw("sameP");
  }

  TLegend* leg = new TLegend(0.22,0.75,0.60,0.89);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  //  leg->SetNColumns(3);
 if ( !histoName.Contains("_eff") && !histoName.Contains("_fr") )  {
   //  leg->AddEntry(h_1,"from AOD, uncorr","PL");
  leg->AddEntry(h_1,"PFCand loop, uncorr","PL");
  leg->AddEntry(h_2,"PFCand loop, area corr","PL");
  leg->AddEntry(h_3,"PFCand loop, simple DeltaBeta","PL");
  leg->AddEntry(h_4,"PFCand loop, DeltaBetaWeights","PL");
 }
 else {
   //  leg->AddEntry(h_1,"from AOD, uncorr, <0.15","PL");
  leg->AddEntry(h_1,"PFCand loop, uncorr, <0.15","PL");
  leg->AddEntry(h_2,"PFCand loop, area corr, <0.15","PL");
  leg->AddEntry(h_3,"PFCand loop, simple DeltaBeta, <0.15","PL");
  leg->AddEntry(h_4,"PFCand loop, DeltaBetaWeights, <0.15","PL");
  leg->AddEntry(h_5,"PFCand loop, DeltaBetaWeights, <0.10","PL");
 }
  leg->Draw();

  c.SaveAs("compare/"+histoName+"_1020.pdf");

}

/*
root -b -q comparePlots.C\(\"nlayer_strip_3d\",10\)
for plot in {algo,isHighPurity,pT,eta,phi,dxy,dz,nlayer_strip_3d,nlayer_tot_3d,nlayer_strip,nlayer_pixel}; do for algo in {9,10}; do root -b -q comparePlots.C\(\"${plot}\",${algo}\); done; done
*/
