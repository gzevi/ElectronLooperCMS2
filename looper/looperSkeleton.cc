#include "dilepStudyLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include <bitset>
#include "TGraph.h"

#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
//#include "../Tools/pfjetMVAtools.h"

#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/muonSelections.h"
#include "../CORE/jetSelections.h"
#include "../CORE/metSelections.h"
#include "../CORE/jetSmearingTools.h"
#include "../CORE/jetcorr/JetCorrectionUncertainty.h"
#include "../CORE/MITConversionUtilities.h"

//#include "CMS2.h" // This can be a different CMS2.h from the one in ../CORE and ../Tools, which have to be based on the full CMS2.h 
// --> actually not! it always seems to load the one in ../CORE (this is done by jetSelectons and jetSmearingTools) --> so need to comment out jetSelections and jetSmearingTools
// --> also keep in mind that ../Tools has to keep its own CMS2.h, otherwise it won't compile, and same for ../CORE
// --> To check for conflict with other CMS2.h, just comment out the #ifndef and #endif lines in the local CMS2.h (this will show where the conflicting definitions are)
// --> Actually can't get to run in this configuration. CMS2.h has to be included by someone above, maybe to be compiled? not clear. otherwise can't find any variables at all...
// --> So we're stuck with adding the variables by hand to ../CORE/CMS2.h... ??? YES. 
// --> Actually the local CMS2.h and .cc are the ones to commit! At the beginning of a new setup they will be copied to ../CORE
 
bool verbose              = false;
bool doTenPercent         = false;
bool doLowerPtThresh      = false;
bool doEM                 = false;
bool doME                 = false;
bool requireTrigMatch     = true;
bool doSS                 = false;
bool doOS                 = false;
bool m_miniAOD            = true;
      bool is25ns = false;

bool  m_newISO = false;
bool  m_newID = false;

using namespace std;
using namespace tas;

//--------------------------------------------------------------------

int getMotherIndex(int motherid){
  for(unsigned int i = 0; i < genps_id().size() ; i++){
    if( motherid == genps_id().at(i) ) return i;
  }
  
  return -1;
}

//--------------------------------------------------------------------

float dilepStudyLooper::dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 ) { 
  
  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();
  
  return sqrt(dphi*dphi + deta*deta);
  
}
//--------------------------------------------------------------------

float dilepStudyLooper::dRbetweenVectors2(LorentzVector vec1, LorentzVector vec2 ) { 
  
  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();
  
  return (dphi*dphi + deta*deta);
  
}
//--------------------------------------------------------------------
float dilepStudyLooper::Mt( LorentzVector p4, float met, float met_phi )
{
    return sqrt( 2*met*( p4.pt() - ( p4.Px()*cos(met_phi) + p4.Py()*sin(met_phi) ) ) );
}
//--------------------------------------------------------------------
dilepStudyLooper::dilepStudyLooper()
{
  
  std::cout << " construct " << std::endl;
  g_createTree   = false;
  initialized = false;
}

//--------------------------------------------------------------------

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
  already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

// void dilepStudyLooper::InitBaby(){
// }

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void dilepStudyLooper::closeOutput()
{
  outFile->cd();
  //  outTree->Write();
  outFile->Write();
  outFile->Close();
  delete outFile;
}

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

int dilepStudyLooper::ScanChain(TChain* chain, const TString& prefix, int sign, int isocortype)

{
  
  //  cout << "ciao " << isData << endl;
  bool isData = false;
  bool isEE = false;
  bool isEM = false;
  bool isMM = false;
  if( prefix.Contains("data") || prefix.Contains("2012") 
     || prefix.Contains("dimu") || prefix.Contains("diel")
     || prefix.Contains("mueg") ){
    cout << "DATA!!!" << endl;
    isData       = true;
    doTenPercent = false;
    if (prefix.Contains("DoubleEle")) {
      isEE = true;
      std::cout << "DoubleElectron data" << std::endl;
    }
    else if (prefix.Contains("DoubleMu")) {
      isMM = true;
      std::cout << "DoubleMu data" << std::endl;
    }
    else if (prefix.Contains("MuEG")) {
      isEM = true;
      std::cout << "MuEG data" << std::endl;
    }
  }
  
  cout << "IS DATA: " << isData << endl;

  
  if( doTenPercent ) cout << "Processing 10% of MC" << endl;
  

  // -----------------------------------------------------------
  
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  
  makeOutput(prefix);
  //  if(g_createTree) makeTree(prefix, doFakeApp, frmode);
  
  //  BookHistos(prefix);
  BookHistos("h");
  
  cout << " done with initialization "  << endl;
  
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  unsigned int nEventsPreReco = 0;
  unsigned int nEventsPass = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;
  
  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;
  
  // test chain
  if (!chain)
  {
    throw std::invalid_argument("at::ScanChain: chain is NULL!");
  }
  if (chain->GetListOfFiles()->GetEntries()<1)
  {
    throw std::invalid_argument("at::ScanChain: chain has no files!");
  }
  if (not chain->GetFile())
  {
    throw std::invalid_argument("at::ScanChain: chain has no files or file path is invalid!");
  }

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());
    
    cout << currentFile->GetTitle() << endl;
    
    if (!f || f->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain is invalid or corrupt: %s", currentFile->GetTitle()));
    }
    
    // get the trees in each file
    // TTree *tree = (TTree*)f->Get("Events");
    TTree *tree = dynamic_cast<TTree*>(f->Get("Events"));
    if (!tree || tree->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain has an invalid TTree or is corrupt: %s", currentFile->GetTitle()));
    }
    
    //Matevz
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);
    
    cms2.Init(tree);
    
    unsigned int nEntries = tree->GetEntries();
    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;
      /////////      cout << nEventsTotal << endl;
      
      if( doTenPercent ){
        if( !(nEventsTotal%10==0) ) continue;
      }
      
      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
        
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
          
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }
      
      //Matevz
      tree->LoadTree(z);
      
      cms2.GetEntry(z);

      //onlyAOD if( evt_ww_rho_vor() != evt_ww_rho_vor() ){
      //onlyAOD   cout << "Skipping event with rho = nan!!!" << endl;
      //onlyAOD   continue;
      //onlyAOD }
      
      //      InitBaby();
      
      isdata_ = isData ? 1 : 0;
      
      if( verbose ){
        cout << "-------------------------------------------------------"   << endl;
        cout << "Event " << z                                               << endl;
        cout << "File  " << currentFile->GetTitle()                         << endl;
        //cout << evt_dataset().at(0) << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
        cout << "-------------------------------------------------------"   << endl;
      }
      
      //      TString datasetname(evt_dataset().at(0));

      TString filename = currentFile->GetTitle();

      
      if (!isData) {
	for(unsigned int idx = 0; idx < genps_id().size(); idx++) {
	  
	  int pid = abs(genps_id().at(idx));          
	  if(genps_status().at(idx) != 3 && genps_status().at(idx) != 22 && genps_status().at(idx) != 23) continue; // this shouldn't happen, since only status 3 are in CMS2          
	  if(pid != 11) continue;
	  if ( genps_p4().at(idx).pt() < 10.) continue;
	  if ( fabs(genps_p4().at(idx).eta()) > 2.5) continue;

	}
      }
        
    ++nEventsPass;
    
    //      outTree->Fill();
    
  } // entries
  
  delete f;
} // currentFile


cout << endl;
cout << "Sample: " << prefix << endl;
cout << endl;
cout << "Processed events: " << nEventsTotal << endl;

// outFile->cd(); // Make sure histograms get written out

// hSetEff["h_nvtx_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_nvtx_num_true"] , hSetEff["h_nvtx_den_true"] , "h_nvtx_eff", "Eff ;# vtxs");
// hSetEff["h_nvtx_fr"]  = (TH1F*) MakeEfficiencyPlot(hSetEff["h_nvtx_num_fake"] , hSetEff["h_nvtx_den_fake"] , "h_nvtx_fr", "Eff ;# vtxs");
// hSetEff["h_eta_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_eta_num_true"] , hSetEff["h_eta_den_true"] , "h_eta_eff", "Eff ;eta");
// hSetEff["h_eta_fr"]  = (TH1F*) MakeEfficiencyPlot(hSetEff["h_eta_num_fake"] , hSetEff["h_eta_den_fake"] , "h_eta_fr", "Eff ;eta");
// hSetEff["h_pt_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_pt_num_true"] , hSetEff["h_pt_den_true"] , "h_pt_eff", "Eff ;pt");
// hSetEff["h_pt_fr"]  = (TH1F*) MakeEfficiencyPlot(hSetEff["h_pt_num_fake"] , hSetEff["h_pt_den_fake"] , "h_pt_fr", "Eff ;pt");

closeOutput();
//  if(g_createTree) closeTree();

already_seen.clear();

if (nEventsChain != nEventsTotal) 
std::cout << "ERROR: number of events from files (" << nEventsChain 
<< ") is not equal to total number of processed events (" << nEventsTotal << ")" << std::endl;

return 0;

}


//--------------------------------------------------------------------

TH1F* dilepStudyLooper::MakeEfficiencyPlot(TH1F* num_hist, TH1F* den_hist, const std::string& name, const std::string& title) // from PAC/RootTools/MiscTools
{
  // check that hists are valid
  if (!num_hist || !den_hist)
    {
      throw runtime_error("rt::MakeEfficiencyPlot: one of the Histograms are NULL");
    }
  
  // verify that all histograms have same binning
  if (den_hist->GetNbinsX() != num_hist->GetNbinsX()) 
    {
      throw runtime_error("rt::MakeEfficiencyPlot: Histograms must have same number of bins");
    }
  
  // get the new histogram
  TH1F* temp = (TH1F*) num_hist->Clone(name.c_str());
  temp->SetTitle(title.empty() ? name.c_str() : title.c_str());
  temp->Reset();
  if (!temp->GetSumw2N())
    {
      temp->Sumw2();
    }
  
  // Do the calculation
  temp->Divide(num_hist, den_hist, 1.0, 1.0, "B");
  
  // Done
  return temp;
}

void dilepStudyLooper::BookHistos(const TString& prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;
  
  // TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  // rootdir->cd();
  if (outFile) outFile->cd();
  
  const int max_nvtx = 40;
  const int nbins_el = 8;
  const int nbins_mu = 8;
  
  h_nvtx = new TH1F(Form("%s_nvtx",prefix.Data()),";N(vtx)",max_nvtx,0,max_nvtx);

  //  hSetEff["h_pt_num_fake"] = new TH1F(Form("%s_pt_num_fake",prefix.Data()),";pt", 9, 10, 100);

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()

 
void dilepStudyLooper::makeOutput(const TString& prefix){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  
  TString tpsuffix = "";
  if( doTenPercent ) tpsuffix = "_tenPercent";
  
  outFile   = new TFile(Form("output/%s%s.root",prefix.Data(),tpsuffix.Data()), "RECREATE");
  //  outFile   = new TFile(Form("output/%s/%s_smallTree%s%s.root",g_version,prefix.Data(),frsuffix,tpsuffix), "RECREATE");iso03_pf2012_ch
  //  outFile   = new TFile("baby.root","RECREATE");
  outFile->cd();
  //  outTree = new TTree("t","Tree");
  
  //Set branch addresses
  //variables must be declared in dilepStudyLooper.h
}

void dilepStudyLooper::labelAxis(TH2F* h, int axis, int lep) {
  
  TAxis* h_axis;
  if (axis == 0) h_axis = h->GetXaxis();
  else if (axis == 1) h_axis = h->GetYaxis();
  else {
    std::cout << "labelAxis: didn't recognize axis: " << axis << ", returning.." << std::endl;
    return;
  }
  
  // electrons
  if (lep == 11) {
    h_axis->SetBinLabel(1,"All");
    h_axis->SetBinLabel(2,"Cand");
    h_axis->SetBinLabel(3,"WWDenom");
    h_axis->SetBinLabel(4,"Loose");
    h_axis->SetBinLabel(5,"Med");
    h_axis->SetBinLabel(6,"Iso");
    h_axis->SetBinLabel(7,"L+Iso");
    h_axis->SetBinLabel(8,"M+Iso");
  }
  
  // muons
  else if (lep == 13) {
    h_axis->SetBinLabel(1,"All");
    h_axis->SetBinLabel(2,"Cand");
    h_axis->SetBinLabel(3,"Tight");
    h_axis->SetBinLabel(4,"Iso < 1.0");
    h_axis->SetBinLabel(5,"Iso < 0.7");
    h_axis->SetBinLabel(6,"Iso < 0.4");
    h_axis->SetBinLabel(7,"Iso < 0.15");
    h_axis->SetBinLabel(8,"T+Iso");
  }
  
  else {
    std::cout << "labelAxis: didn't recognize lep: " << lep << ", returning.." << std::endl;
    return;
  }
  
  return;
}


//----------------------------

float dilepStudyLooper::getdphi( float phi1 , float phi2 ){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}


