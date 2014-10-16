
void makeEffMergedFiles(const TString& infile) {
  TFile *file = TFile::Open(infile, "UPDATE");
  TH1D * eff = file->Get("h_ptb_eff")->Clone("eff");
  eff->Divide(h_ptb_num, h_ptb_den, 1.0, 1.0, "B");
  eff->Write();
  TH1D * effW2 = file->Get("h_ptb_numW2")->Clone("effW2");
  effW2->Divide(h_ptb_numW2, h_ptb_denW2, 1.0, 1.0, "B");
  effW2->Write();
  TH1D * eff05 = file->Get("h_ptb05_eff")->Clone("eff05");
  eff05->Divide(h_ptb05_num, h_ptb05_den, 1.0, 1.0, "B");
  eff05->Write();
  TH1D * effiso = file->Get("h_ptbiso_eff")->Clone("effiso");
  effiso->Divide(h_ptbiso_num, h_ptbiso_den, 1.0, 1.0, "B");
  effiso->Write();
  file->Close();
}
