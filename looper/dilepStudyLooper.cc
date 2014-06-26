#include "dilepStudyLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include <bitset>

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
// --> So we're stuck with adding the variables by hand to ../CORE/CMS2.h... ???
 
bool verbose              = false;
bool doTenPercent         = false;
bool doLowerPtThresh      = false;
bool doEM                 = false;
bool doME                 = false;
bool requireTrigMatch     = true;
bool doSS                 = false;
bool doOS                 = false;
bool m_miniAOD            = true;

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
  

  
  const TString trigname_mm = "HLT_Mu17_Mu8_v";
  const TString trigname_mmtk = "HLT_Mu17_TkMu8_v";
  const TString trigname_me = "HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
  const TString trigname_em = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
  const TString trigname_ee = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";

  const TString trigname_e8  = "Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
  const TString trigname_e17 = "Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
  const TString trigname_e8jet  = "Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v";
  const TString trigname_e17jet = "Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v";
  
  //------------------------------------------------------------------------------------------------------
  // set json\, vertex reweighting function and msugra cross section files
  //------------------------------------------------------------------------------------------------------
  
  if( !initialized ){
    
    //set json
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );
    
    //    if( prefix.Contains("ttall_massivebin") ) 
    //    set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer12MC_PUS10_19fb_Zselection.root",true);
    
    initialized = true;
  }
  
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
  int nSkip_els_conv_dist = 0;
  int nTrueEle = 0;
  int nPromptEle = 0;
  int nNonPromptEle = 0;
  int nFakeEle = 0;
  int nPromptEleL = 0;
  int nNonPromptEleL = 0;
  int nFakeEleL = 0;
  int nPromptEleM = 0;
  int nNonPromptEleM = 0;
  int nFakeEleM = 0;
  int nPromptEleT = 0;
  int nNonPromptEleT = 0;
  int nFakeEleT = 0;
  int nFOEle = 0; // This is the basis for the Fake Rate. One could use all Reco Electrons, but ideally should use something that's not biased by trigger

  int nPromptEleBin[6] = {};
  int nPromptEleLBin[6] = {};
  int nPromptEleMBin[6] = {};
  int nPromptEleTBin[6] = {};
  int nFakeEleBin[6] = {};
  int nFakeEleLBin[6] = {};
  int nFakeEleMBin[6] = {};
  int nFakeEleTBin[6] = {};

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
        cout << evt_dataset().at(0) << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
        cout << "-------------------------------------------------------"   << endl;
      }
      
      TString datasetname(evt_dataset().at(0));
      
      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------
      
      //      if( !cleaning_goodVertexApril2011() )                          continue;
      if( isData && !goodrun(evt_run(), evt_lumiBlock()) ) continue;

      //---------------------
      // skip duplicates
      //---------------------
      
      if( isData ) {
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){
          continue;
        }
      }
      TString filename = currentFile->GetTitle();
      if (filename.Contains("CMS3") || filename.Contains("miniAOD")) m_miniAOD = true;
      else m_miniAOD = false;

      //-------------------------------------
      // skip events with bad els_conv_dist
      //-------------------------------------
// onlyAOD      
// onlyAOD      bool skipEvent = false;
// onlyAOD      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
// onlyAOD        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
// onlyAOD          skipEvent = true;
// onlyAOD        }
// onlyAOD        if( els_sigmaIEtaIEta().at(iEl) != els_sigmaIEtaIEta().at(iEl) ){
// onlyAOD          skipEvent = true;
// onlyAOD        }
// onlyAOD        if( els_sigmaIEtaIEtaSC().at(iEl) != els_sigmaIEtaIEtaSC().at(iEl) ){
// onlyAOD          skipEvent = true;
// onlyAOD        }
// onlyAOD      }
// onlyAOD
// onlyAOD      if( skipEvent ){
// onlyAOD        nSkip_els_conv_dist++;
// onlyAOD        continue;
// onlyAOD      }

      //---------------------------------------------
      // count vertices
      //---------------------------------------------
      int nvtx = 0;
      for (size_t v = 0; v < cms2.vtxs_position().size(); ++v){
        if(isGoodVertex(v)) ++nvtx;
      }

      //---------------------------------------------
      // gen selection for DY MC
      //---------------------------------------------
      
      //if (isDYMC) {
      //  int nel = 0;
      //  int nmu = 0;
      //  // require 2 gen muons, to study Z->mm
      //  for (unsigned int igen = 0; igen < genps_p4().size(); ++igen) {
      //    if (genps_status().at(igen) != 3) continue;
      //    int id = genps_id().at(igen);
      //    if ((abs(id) != 11) && (abs(id) != 13)) continue;
      //    if (abs(genps_id_mother().at(igen)) != 23) continue;
      //    if (abs(id) == 11) ++nel;
      //    else if (abs(id) == 13) ++nmu;
      //  }
      //  if (nmu < 2) continue;
      //}

      //---------------------------------------------
      // trigger selection
      //---------------------------------------------
      bool pass_trig_mm        = passUnprescaledHLTTriggerPattern(trigname_mm.Data())                                     ? 1 : 0;
      bool pass_trig_mmtk      = passUnprescaledHLTTriggerPattern(trigname_mmtk.Data())                                   ? 1 : 0;
      bool pass_trig_me        = passUnprescaledHLTTriggerPattern(trigname_me.Data()) ? 1 : 0;
      bool pass_trig_em        = passUnprescaledHLTTriggerPattern(trigname_em.Data()) ? 1 : 0;
      bool pass_trig_ee        = passUnprescaledHLTTriggerPattern(trigname_ee.Data()) ? 1 : 0;
      bool pass_trig_e8        = passHLTTriggerPattern(trigname_e8.Data()) ? 1 : 0;
      bool pass_trig_e17       = passHLTTriggerPattern(trigname_e17.Data()) ? 1 : 0;
      bool pass_trig_e8jet     = passHLTTriggerPattern(trigname_e8jet.Data()) ? 1 : 0;
      bool pass_trig_e17jet    = passHLTTriggerPattern(trigname_e17jet.Data()) ? 1 : 0;
      pass_trig_ee = ( pass_trig_e8 || pass_trig_e8jet || pass_trig_e17 || pass_trig_e17jet); // patch to work on single electron trigger
      
      // require the event to pass one of these triggers (only for data)
      if (isData) {
        if (isEE && !pass_trig_ee) continue;
        if (isEM) {
          bool pass = false;
          if (doEM && pass_trig_em) pass = true;
          if (doME && pass_trig_me) pass = true;
          if (!pass) continue;
        }
        if (isMM && !pass_trig_mm && !pass_trig_mmtk) continue;
        if (!pass_trig_mm && !pass_trig_mmtk && !pass_trig_me && !pass_trig_em && !pass_trig_ee) continue;
      }

      // INFO: Print list of passed triggers
      //cout<<"passed triggers: ";
      //for( unsigned int itrig = 0 ; itrig < cms2.hlt_trigNames().size() ; ++itrig ){
      //	TString name(hlt_trigNames().at(itrig));
      //	if( passHLTTrigger( name )) cout<<name<<" ";
      //}
      //cout<<endl;


      //---------------------------------------------
      // reco electron selection
      //---------------------------------------------
      
      float pt_lead_el = -1.;
      int idx_lead_el = -1;
      float pt_subl_el = -1.;
      int idx_subl_el = -1;
      
      float pt_lead_denom_el = -1.;
      int idx_lead_denom_el = -1;
      float pt_subl_denom_el = -1.;
      int idx_subl_denom_el = -1;
      
      float pt_lead_mediso_el = -1.;
      int idx_lead_mediso_el = -1;
      float pt_subl_mediso_el = -1.;
      int idx_subl_mediso_el = -1;

      dilepStudyLooper::electron eleDencuts; // should correspond to a dilepton trigger (not WP80, which is too tight)
      eleDencuts.sieie = 0.011;
      eleDencuts.dEtaIn = 0.01;
      eleDencuts.dPhiIn = 0.15;
      eleDencuts.hOverE = 0.10;
      eleDencuts.d0corr = 999;
      eleDencuts.z0corr = 999;
      eleDencuts.ooemoop = 0.05;
      eleDencuts.iso_cor = 999;
      eleDencuts.pfiso_ch = 999;
      eleDencuts.pfiso_em = 999;
      eleDencuts.pfiso_nh = 999;
      eleDencuts.detiso_ch = 0.2;
      eleDencuts.detiso_em = 0.3; /// should be 0.2, but shifted by PU
      eleDencuts.detiso_nh = 0.2;
      eleDencuts.valid_pixelHits = 0;
      eleDencuts.lost_pixelhits = 999;
      eleDencuts.vtxFitConversion = false; // false means don't check
      eleDencuts.chi2n = 999;

      dilepStudyLooper::electron eleDencutsEE = eleDencuts;
      eleDencutsEE.sieie = 0.031;
      eleDencutsEE.dEtaIn = 0.01;
      eleDencutsEE.dPhiIn = 0.10;
      eleDencutsEE.hOverE = 0.075;


      dilepStudyLooper::electron eleWP80cuts;
      eleWP80cuts.sieie = 0.01;
      eleWP80cuts.dEtaIn = 0.007;
      eleWP80cuts.dPhiIn = 0.06;
      eleWP80cuts.hOverE = 0.10;
      eleWP80cuts.d0corr = 999;
      eleWP80cuts.z0corr = 999;
      eleWP80cuts.ooemoop = 0.05;
      eleWP80cuts.iso_cor = 999;
      eleWP80cuts.pfiso_ch = 999;
      eleWP80cuts.pfiso_em = 999;
      eleWP80cuts.pfiso_nh = 999;
      eleWP80cuts.detiso_ch = 0.05;
      eleWP80cuts.detiso_em = 0.15;
      eleWP80cuts.detiso_nh = 0.1;
      eleWP80cuts.valid_pixelHits = 0;
      eleWP80cuts.lost_pixelhits = 999;
      eleWP80cuts.vtxFitConversion = false; // false means don't check
      eleWP80cuts.chi2n = 999;

      dilepStudyLooper::electron eleWP80cutsEE = eleWP80cuts;
      eleWP80cutsEE.sieie = 0.03;
      eleWP80cutsEE.dPhiIn = 0.03;
      eleWP80cutsEE.hOverE = 0.05;
      eleWP80cutsEE.valid_pixelHits = 3;
      eleWP80cutsEE.pfiso_em = 0.1;
      eleWP80cutsEE.chi2n = 1.5;
      
      dilepStudyLooper::electron eleLcuts;
      eleLcuts.sieie = 0.01;
      eleLcuts.dEtaIn = 0.007;
      eleLcuts.dPhiIn = 0.15;
      eleLcuts.hOverE = 0.12;
      eleLcuts.d0corr = 0.02;
      eleLcuts.z0corr = 0.2;
      eleLcuts.ooemoop = 0.05;
      eleLcuts.iso_cor = 0.15;
      eleLcuts.pfiso_ch = 999;
      eleLcuts.pfiso_em = 999;
      eleLcuts.pfiso_nh = 999;
      eleLcuts.detiso_ch = 999;
      eleLcuts.detiso_em = 999;
      eleLcuts.detiso_nh = 999;
      eleLcuts.valid_pixelHits = 0;
      eleLcuts.lost_pixelhits = 1;
      eleLcuts.vtxFitConversion = true; // false means don't check
      eleLcuts.chi2n = 999;

      dilepStudyLooper::electron eleLcutsEE = eleLcuts;
      eleLcutsEE.sieie = 0.03;
      eleLcutsEE.dEtaIn = 0.009;
      eleLcutsEE.dPhiIn = 0.10;
      eleLcutsEE.hOverE = 0.10;
      eleLcutsEE.iso_cor = 0.15; // should be 0.10 for <20 GeV

      dilepStudyLooper::electron eleMcuts;
      eleMcuts.sieie = 0.01;
      eleMcuts.dEtaIn = 0.004;
      eleMcuts.dPhiIn = 0.06;
      eleMcuts.hOverE = 0.12;
      eleMcuts.d0corr = 0.02;
      eleMcuts.z0corr = 0.1;
      eleMcuts.ooemoop = 0.05;
      eleMcuts.iso_cor = 0.15;
      eleMcuts.pfiso_ch = 999;
      eleMcuts.pfiso_em = 999;
      eleMcuts.pfiso_nh = 999;
      eleMcuts.detiso_ch = 999;
      eleMcuts.detiso_em = 999;
      eleMcuts.detiso_nh = 999;
      eleMcuts.valid_pixelHits = 0;
      eleMcuts.lost_pixelhits = 1;
      eleMcuts.vtxFitConversion = true; // false means don't check
      eleMcuts.chi2n = 999;

      dilepStudyLooper::electron eleMcutsEE = eleMcuts;
      eleMcutsEE.sieie = 0.03;
      eleMcutsEE.dEtaIn = 0.007;
      eleMcutsEE.dPhiIn = 0.03;
      eleMcutsEE.hOverE = 0.10;
      eleMcutsEE.iso_cor = 0.15; // should be 0.10 for <20 GeV

      dilepStudyLooper::electron eleTcuts;
      eleTcuts.sieie = 0.01;
      eleTcuts.dEtaIn = 0.004;
      eleTcuts.dPhiIn = 0.03;
      eleTcuts.hOverE = 0.12;
      eleTcuts.d0corr = 0.02;
      eleTcuts.z0corr = 0.1;
      eleTcuts.ooemoop = 0.05;
      eleTcuts.iso_cor = 0.1;
      eleTcuts.pfiso_ch = 999;
      eleTcuts.pfiso_em = 999;
      eleTcuts.pfiso_nh = 999;
      eleTcuts.detiso_ch = 999;
      eleTcuts.detiso_em = 999;
      eleTcuts.detiso_nh = 999;
      eleTcuts.valid_pixelHits = 0;
      eleTcuts.lost_pixelhits = 0;
      eleTcuts.vtxFitConversion = true; // false means don't check
      eleTcuts.chi2n = 999;

      dilepStudyLooper::electron eleTcutsEE = eleTcuts;
      eleTcutsEE.sieie = 0.03;
      eleTcutsEE.dEtaIn = 0.005;
      eleTcutsEE.dPhiIn = 0.02;
      eleTcutsEE.hOverE = 0.10;
      eleTcutsEE.iso_cor = 0.1; // should be 0.10 for <20 GeV
      
      if (!isData) {
	for(unsigned int idx = 0; idx < genps_id().size(); idx++) {
	  
	  int pid = abs(genps_id().at(idx));          
	  if(genps_status().at(idx) != 3 && genps_status().at(idx) != 22 && genps_status().at(idx) != 23) continue; // this shouldn't happen, since only status 3 are in CMS2          
	  if(pid != 11) continue;
	  if ( genps_p4().at(idx).pt() < 10.) continue;
	  if ( fabs(genps_p4().at(idx).eta()) > 2.5) continue;
	  nTrueEle++;
	}
      }
      
      for (unsigned int iel = 0; iel < els_p4().size(); ++iel) { 
        LorentzVector& el_p4 = els_p4().at(iel);
        float pt = el_p4.pt();
        float eta = els_etaSC().at(iel);
        if (!isData && (pt < 10 || pt > 50)) continue;
        if (isData && (pt < 10 || pt > 35)) continue;
	// Quick check for single electron:
	int neleDen = 1;
	for (unsigned int jel = 0; jel < els_p4().size(); ++jel) {
	  if (!isData) continue;
	  if (iel==jel) continue;
	  if (els_p4().at(iel).pt() < 10 ) continue;
	  dilepStudyLooper::electron eleStruct2; 
	  fillElectronStructure(jel, eleStruct2,  false, false, 0, false, false, false);
	  if (passElectronCuts(eleStruct2, eleDencuts, eleDencutsEE)) neleDen++;
	  //cout<<"neleDen = "<<neleDen<<endl;
	}
	if (neleDen > 1) continue;

	//if (eta < 0 ) continue;
	//	if (isData && evt_pfmet() > 20) continue;	
	//	if (isData && Mt( el_p4, evt_pfmet(), evt_pfmetphi() ) > 25) continue; // need to calculate mt	
	//if (pt < 20) continue;
	//	if (pt < 10 || pt > 20) continue;
	if (fabs(eta) > 2.5) continue;
        h_el_pt->Fill(pt,1);
        h_el_eta->Fill(eta,1);

	int effBin = 0;
	if (fabs(eta) <= 1.479) {
	  if (pt < 15.) effBin = 0;
	  else if (pt < 20.) effBin = 2;
	  else effBin = 4;	  
	}
	else {
	  if (pt < 15.) effBin = 1;
	  else if (pt < 20.) effBin = 3;
	  else effBin = 5;
	}

	nFOEle++;
	

	bool ecalDriven = (cms2.els_type()[iel] & (1<<ISECALDRIVEN));
	bool trackerDriven = (cms2.els_type()[iel] & (1<<ISTRACKERDRIVEN));
	int seedCode = 0;
	if (ecalDriven) seedCode = 1;
	if (trackerDriven) seedCode = 2;
	if (ecalDriven && trackerDriven) seedCode = 3;
	

	bool useMap = false;
	bool doPFCandLoop = true;
	int useDeltaBetaWeights = 1; // 0: no, 1: 1/DR^2, 2: log(pt^2/DR^2)
	bool DeltaBetaSimple = false; 
	float isocut = 0.15;
	bool areaCorrection = false;
	//if ( prefix.Contains("700") || prefix.Contains("CMS3") ) useMap = true;

	// don't use map
	if (isocortype==0) { doPFCandLoop = false; useDeltaBetaWeights = 0; DeltaBetaSimple = false; areaCorrection = false;} // uncorrected
	if (isocortype==1) { doPFCandLoop = false; useDeltaBetaWeights = 0; DeltaBetaSimple = false; areaCorrection = true; }// area
	if (isocortype==2) { doPFCandLoop = true; useDeltaBetaWeights = 0; DeltaBetaSimple = true; areaCorrection = false;  }// DeltaBeta
	if (isocortype==3) { doPFCandLoop = true; useDeltaBetaWeights = 1; DeltaBetaSimple = false; areaCorrection = false;  }// DeltaBetaWeights DR
	if (isocortype==4) { doPFCandLoop = true; useDeltaBetaWeights = 1; DeltaBetaSimple = false; areaCorrection = false; isocut = 0.10; }// DeltaBetaWeights DR. cut at 0.1
	if (isocortype==5) { doPFCandLoop = true; useDeltaBetaWeights = 0; DeltaBetaSimple = false; areaCorrection = false; }// uncorrected, do loop
	if (isocortype==6) { doPFCandLoop = true; useDeltaBetaWeights = 2; DeltaBetaSimple = false; areaCorrection = false; }// DeltaBetaWeight log pT DR
	// use map
	if (isocortype==10) { useMap = true; doPFCandLoop = false; useDeltaBetaWeights = 0; DeltaBetaSimple = false; areaCorrection = false;} // uncorrected
	if (isocortype==11) { useMap = true; doPFCandLoop = false; useDeltaBetaWeights = 0; DeltaBetaSimple = false; areaCorrection = true; }// area
	if (isocortype==12) { useMap = true; doPFCandLoop = true; useDeltaBetaWeights = 0; DeltaBetaSimple = true; areaCorrection = false;  }// DeltaBeta
	if (isocortype==13) { useMap = true; doPFCandLoop = true; useDeltaBetaWeights = 1; DeltaBetaSimple = false; areaCorrection = false;  }// DeltaBetaWeights DR
	if (isocortype==14) { useMap = true; doPFCandLoop = true; useDeltaBetaWeights = 1; DeltaBetaSimple = false; areaCorrection = false; isocut = 0.10; }// DeltaBetaWeights DR. cut at 0.1
	if (isocortype==15) { useMap = true; doPFCandLoop = true; useDeltaBetaWeights = 0; DeltaBetaSimple = false; areaCorrection = false; }// uncorrected, do loop
	if (isocortype==16) { useMap = true; doPFCandLoop = true; useDeltaBetaWeights = 2; DeltaBetaSimple = false; areaCorrection = false; }// DeltaBetaWeight log pT DR


        // Fill the electron structure
        dilepStudyLooper::electron eleStruct; 
	fillElectronStructure(iel, eleStruct,  useMap, doPFCandLoop, useDeltaBetaWeights, DeltaBetaSimple, areaCorrection, true);
        
        // Fill "hall"              
        fillElectronQuantities(hSet1, eleStruct);
        
        // Truth matching: t=true (status 3), f=fake (including non-prompt)
        // Can't separate non-prompt because we only have status 3 particles in CMS2 ntuple
        bool truthmatched = false;
	bool fakematched = false;
	bool nonpromptmatched = false;
	if (!isData) {
	  for(unsigned int idx = 0; idx < genps_id().size(); idx++) {
	    
	    int pid = abs(genps_id().at(idx));          
	    if(genps_status().at(idx) != 3 && genps_status().at(idx) != 22 && genps_status().at(idx) != 23) continue; // this shouldn't happen, since only status 3 are in CMS2          
	    if(pid != 11) continue;
	    float dr = dRbetweenVectors(genps_p4().at(idx), el_p4);
	    h_el_truthmatch_DR->Fill(dr, 1);
	    if ( dr > 0.02 ) continue;
	    else truthmatched = true;
	  }
	  if (truthmatched) {nPromptEle++; nPromptEleBin[effBin]++;}
	  else if (abs(els_mc_id().at(iel) == 11)) { nNonPromptEle++; nonpromptmatched = true; }
	  else { nFakeEle++; nFakeEleBin[effBin]++; fakematched = true; }
	}
	if (isData) fakematched = true;
	//	cout<<" els_mc_id "<<els_mc_id().at(iel)<<", els_mc3_id "<<els_mc3_id().at(iel)<<", els_mc_motherid "<<els_mc_motherid().at(iel)<<", els_mc_motherid "<<els_mc3_motherid().at(iel)<<endl;
        if (fabs(eta) <= 1.479)   truthmatched ? fillElectronQuantities(hSet2, eleStruct) : ( fakematched ? fillElectronQuantities(hSet2f, eleStruct) : fillElectronQuantities(hSet2np, eleStruct) );
        else   truthmatched ? fillElectronQuantities(hSet2E, eleStruct) : ( fakematched ? fillElectronQuantities(hSet2Ef, eleStruct) : fillElectronQuantities(hSet2Enp, eleStruct) );
        

	//	if (fabs(eta) <= 1.479 && truthmatched ) cout<< "sieie "<<eleStruct.sieie<<endl;

        // Now start with the selections
        // First plot the N-1 distributions for the WP80 selections. 
        // Plot all variables for electrons passing all cuts 2->N

	// After denominator cuts
	if (passElectronCuts(eleStruct, eleDencuts, eleDencutsEE)) {
	  if (fabs(eta) <= 1.479) truthmatched ? fillElectronQuantities(hSet7, eleStruct) : fillElectronQuantities(hSet7f, eleStruct);
	  else truthmatched ? fillElectronQuantities(hSet7E, eleStruct) : fillElectronQuantities(hSet7Ef, eleStruct);
	}
	
	// N-1
        if (fabs(eta) <= 1.479) truthmatched ? fillElectronQuantitiesN1(hSet3, eleStruct, eleWP80cuts) : fillElectronQuantitiesN1(hSet3f, eleStruct, eleWP80cuts); 
        else truthmatched ? fillElectronQuantitiesN1(hSet3E, eleStruct, eleWP80cutsEE) : fillElectronQuantitiesN1(hSet3Ef, eleStruct, eleWP80cutsEE); 
        

	// Run Without Isolation Cut
      //eleLcuts.iso_cor = 999;
      //eleLcutsEE.iso_cor = 999;
      //eleMcuts.iso_cor = 999;
      //eleMcutsEE.iso_cor = 999;
      //eleTcuts.iso_cor = 999;
      //eleTcutsEE.iso_cor = 999;
      eleDencuts.iso_cor = 999;
      eleDencutsEE.iso_cor = 999;

      // Fill N-1 plots
      // LOOSE
      if (fabs(eta) <= 1.479) truthmatched ? fillElectronQuantitiesN1(hSet4, eleStruct, eleLcuts) :  ( fakematched ? fillElectronQuantitiesN1(hSet4f, eleStruct, eleLcuts) : fillElectronQuantitiesN1(hSet4np, eleStruct, eleLcuts) );
      else truthmatched ? fillElectronQuantitiesN1(hSet4E, eleStruct, eleLcutsEE) :  ( fakematched ? fillElectronQuantitiesN1(hSet4Ef, eleStruct, eleLcutsEE) : fillElectronQuantitiesN1(hSet4Enp, eleStruct, eleLcutsEE) );
      // MEDIUM
      if (fabs(eta) <= 1.479) truthmatched ? fillElectronQuantitiesN1(hSet5, eleStruct, eleMcuts) :  ( fakematched ? fillElectronQuantitiesN1(hSet5f, eleStruct, eleMcuts) : fillElectronQuantitiesN1(hSet5np, eleStruct, eleMcuts) );
      else truthmatched ? fillElectronQuantitiesN1(hSet5E, eleStruct, eleMcutsEE) :  ( fakematched ? fillElectronQuantitiesN1(hSet5Ef, eleStruct, eleMcutsEE) : fillElectronQuantitiesN1(hSet5Enp, eleStruct, eleMcutsEE) );
      // TIGHT
      if (fabs(eta) <= 1.479) truthmatched ? fillElectronQuantitiesN1(hSet6, eleStruct, eleTcuts) :  ( fakematched ? fillElectronQuantitiesN1(hSet6f, eleStruct, eleTcuts) : fillElectronQuantitiesN1(hSet6np, eleStruct, eleTcuts) );
      else truthmatched ? fillElectronQuantitiesN1(hSet6E, eleStruct, eleTcutsEE) :  ( fakematched ? fillElectronQuantitiesN1(hSet6Ef, eleStruct, eleTcutsEE) : fillElectronQuantitiesN1(hSet6Enp, eleStruct, eleTcutsEE) );
        
      // Instead of MEDIUM, fill only if ncluster == 1
//      if ( seedCode == 3 ) {
//	if (fabs(eta) <= 1.479)   truthmatched ? fillElectronQuantities(hSet5, eleStruct) : ( fakematched ? fillElectronQuantities(hSet5f, eleStruct) : fillElectronQuantities(hSet5np, eleStruct) );
//	else   truthmatched ? fillElectronQuantities(hSet5E, eleStruct) : ( fakematched ? fillElectronQuantities(hSet5Ef, eleStruct) : fillElectronQuantities(hSet5Enp, eleStruct) );
//      }



      // Fill Plots After Cuts (not N-1), and counters
	if (passElectronCuts(eleStruct, eleLcuts, eleLcutsEE)) {
	  if (truthmatched) {nPromptEleL++; nPromptEleLBin[effBin]++;}
	  else if (nonpromptmatched)  nNonPromptEleL++;
	  else {nFakeEleL++; nFakeEleLBin[effBin]++;}
//	  if (fabs(eta) <= 1.479) truthmatched ? fillElectronQuantities(hSet4, eleStruct) : ( fakematched ? fillElectronQuantities(hSet4f, eleStruct) : fillElectronQuantities(hSet4np, eleStruct) );
//	  else truthmatched ? fillElectronQuantities(hSet4E, eleStruct) : ( fakematched ? fillElectronQuantities(hSet4Ef, eleStruct) : fillElectronQuantities(hSet4np, eleStruct) );
	}

	if (passElectronCuts(eleStruct, eleMcuts, eleMcutsEE)) {
	  if (truthmatched) {nPromptEleM++; nPromptEleMBin[effBin]++;}
	  else if (nonpromptmatched)  nNonPromptEleM++;
	  else {nFakeEleM++; nFakeEleMBin[effBin]++;}
//	  if (fabs(eta) <= 1.479) truthmatched ? fillElectronQuantities(hSet5, eleStruct) : ( fakematched ? fillElectronQuantities(hSet5f, eleStruct) : fillElectronQuantities(hSet5np, eleStruct) );
//	  else truthmatched ? fillElectronQuantities(hSet5E, eleStruct) : ( fakematched ? fillElectronQuantities(hSet5Ef, eleStruct) : fillElectronQuantities(hSet5Enp, eleStruct) );
	}

	if (passElectronCuts(eleStruct, eleTcuts, eleTcutsEE)) {
	  if (truthmatched) {nPromptEleT++; nPromptEleTBin[effBin]++;}
	  else if (nonpromptmatched)  nNonPromptEleT++;
	  else {nFakeEleT++; nFakeEleTBin[effBin]++;}
//	  if (fabs(eta) <= 1.479) truthmatched ? fillElectronQuantities(hSet6, eleStruct) : ( fakematched ? fillElectronQuantities(hSet6f, eleStruct) : fillElectronQuantities(hSet6np, eleStruct) );
//	  else truthmatched ? fillElectronQuantities(hSet6E, eleStruct) : ( fakematched ? fillElectronQuantities(hSet6Ef, eleStruct) : fillElectronQuantities(hSet6Enp, eleStruct) );
	}

	//Debug: plot after each cut, one cut at the time
//	float isoCut = eleLcuts.iso_cor;
//	electron e = eleStruct;
//	if ( e.pt < 20 && fabs(e.eta) > 1.479 && eleLcuts.iso_cor == 0.10) isoCut = 0.07;
//	if ( e.pt < 20 && fabs(e.eta) > 1.479 && eleLcuts.iso_cor == 0.15) isoCut = 0.10; 
//	bool passVtx = eleLcuts.vtxFitConversion ? !e.vtxFitConversion : true;
//	if (fabs(eta) <= 1.479 && !truthmatched ) {
//	  if ( e.sieie < eleLcuts.sieie ) 	      fillElectronQuantities(	hSetCut1, eleStruct);
//	  if ( fabs(e.dEtaIn) < eleLcuts.dEtaIn )     fillElectronQuantities(	hSetCut2, eleStruct);
//	  if ( fabs(e.dPhiIn) < eleLcuts.dPhiIn )     fillElectronQuantities(	hSetCut3, eleStruct);
//	  if ( e.hOverE < eleLcuts.hOverE ) 	      fillElectronQuantities(	hSetCut4, eleStruct);
//	  if ( fabs(e.d0corr) < eleLcuts.d0corr )     fillElectronQuantities(	hSetCut5, eleStruct);
//	  if ( fabs(e.z0corr) < eleLcuts.z0corr )     fillElectronQuantities(	hSetCut6, eleStruct);
//	  if ( fabs(e.ooemoop) < eleLcuts.ooemoop )   fillElectronQuantities(	hSetCut7, eleStruct);
//	  if ( e.iso_cor < isoCut  ) 	              fillElectronQuantities(	hSetCut8, eleStruct);
//	  if (passVtx)                                fillElectronQuantities(	hSetCut9, eleStruct);
//	}
  
  
	//////////////////////
	// EFFICIENCY PLOTS //
	// Fill "pass" and "all" nvtx plots for "truthmatched" and "Fakes"
	// Dave Evans: "Denominator: Reco Electron (Signal) FO (Fakes)" --> so we should use different denominators? Yes because the meanings are different. Bkg Eff is always with respect to some arbitrary quantity that makes sense
	bool passDenEFF = true && truthmatched;
	//	bool passDenFR = passElectronCuts(eleStruct, eleDencuts, eleDencutsEE) && pt < 30 && fakematched; // just use RECO to start with. Later move to a more "electron-like" object
	bool passDenFR = true && fakematched; // just use RECO to start with. Later move to a more "electron-like" object
	bool passNum = passElectronCuts(eleStruct, eleMcuts, eleMcutsEE);
	//	bool passNum = eleStruct.iso_cor < isocut;
	TString TorF = truthmatched ? "true" : (fakematched ? "fake" : "");
	if ( passDenEFF || passDenFR ) {
	  hSetEff[("h_nvtx_den_"+TorF)]->Fill(nvtx);
	  hSetEff[("h_eta_den_"+TorF)]->Fill(fabs(eta));
	  hSetEff[("h_pt_den_"+TorF)]->Fill(pt);
	  if (passNum) { 
	    hSetEff[("h_nvtx_num_"+TorF)]->Fill(nvtx);
	    hSetEff[("h_eta_num_"+TorF)]->Fill(fabs(eta));
	    hSetEff[("h_pt_num_"+TorF)]->Fill(pt);
	  }
	}


//        // loose ID, no iso
//        bool pass_loose = passElectronSelection_ZMet2012_v3_NoIso( iel,true,true,false);
//        // med ID, no iso
//        bool pass_med = passElectronSelection_Stop2012_v3_NoIso( iel,true,true,false);
//        bool pass_iso = bool(eleStruct.iso_cor < 0.15);
        
        // require match to trigger object
        //        if (requireTrigMatch) {
        //          bool matched = false;
        //          if (isEE) {
        //            matched = objectPassTrigger(els_p4().at(iel),trigname_ee,82);
        //          } else if (isEM) {
        //            if (doEM) matched |= objectPassTrigger(els_p4().at(iel),trigname_em,82);
        //            if (doME) matched |= objectPassTrigger(els_p4().at(iel),trigname_me,82);
        //          }
        //          if (!matched) continue;
        //        }
        
        
      }
      

      //---------------------------------------------
      // putting together event selections
      //---------------------------------------------
      
      h_nvtx->Fill(nvtx);
//      
//      if (isEE && pass_trig_ee) {
//        for (unsigned int ilead=0; ilead < el_lead_flags.size(); ++ilead) {
//          for (unsigned int isubl=0; isubl < el_subl_flags.size(); ++isubl) {
//            h_ee_events->Fill(el_lead_flags.at(ilead),el_subl_flags.at(isubl));
//          }
//        }
//      }
      
      // plot hlt obj kinematics
//      int trigidx = findTriggerIndex(triggerName(trigname_em));
//      std::vector<LorentzVector> trigp4 = hlt_trigObjs_p4()[trigidx];
//      std::vector<int> trigid = hlt_trigObjs_id()[trigidx];
//      
//      LorentzVector el_hlt;
//      LorentzVector mu_hlt;
//      for (unsigned int ihlt = 0; ihlt < trigp4.size(); ++ihlt){
//        if (trigid.at(ihlt) == 82) { // HLT code for electron
//          el_hlt = trigp4.at(ihlt);
//          h_em_el_hlt_pt->Fill(trigp4.at(ihlt).pt());
//          h_em_el_hlt_eta->Fill(trigp4.at(ihlt).eta());
//          if (pt_lead_el < 0.) {
//            h_em_el_hlt_noreco_pt->Fill(trigp4.at(ihlt).pt());
//            h_em_el_hlt_noreco_eta->Fill(trigp4.at(ihlt).eta());
//          }
//        }
//      } // loop on hlt objects
      

        
    ++nEventsPass;
    
    //      outTree->Fill();
    
  } // entries
  
  delete f;
} // currentFile

if( nSkip_els_conv_dist > 0 )
cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

cout << endl;
cout << "Sample: " << prefix << endl;
cout << endl;
cout << "Processed events: " << nEventsTotal << endl;
cout << "Events before reco cuts: " << nEventsPreReco << endl;
cout << "Passed events: " << nEventsPass << endl;
cout << "Gen Electrons: " << nTrueEle << endl;
cout << "FO Electrons: " << nFOEle << endl;
 cout << "All:  \tTruthMatch: " << nPromptEle << "\t NonPrompt: " << nNonPromptEle<< "\t Fake: " << nFakeEle <<"\t Eff: " << 1.*nPromptEle/nTrueEle <<"\t FakeEff: " << 1.*(nFakeEle+nNonPromptEle)/(nFakeEle+nNonPromptEle) << endl;
cout << "Loose: \tTruthMatch: " << nPromptEleL << "\t NonPrompt: " << nNonPromptEleL<< "\t Fake: " << nFakeEleL <<"\t Eff: " << 1.*nPromptEleL/nTrueEle <<"\t FakeEff: " << 1.*(nFakeEleL+nNonPromptEleL)/(nFakeEle+nNonPromptEle) << endl;
cout << "Medium:\tTruthMatch: " << nPromptEleM << "\t NonPrompt: " << nNonPromptEleM<< "\t Fake: " << nFakeEleM <<"\t Eff: " << 1.*nPromptEleM/nTrueEle <<"\t FakeEff: " << 1.*(nFakeEleM+nNonPromptEleM)/(nFakeEle+nNonPromptEle) << endl;
cout << "Tight: \tTruthMatch: " << nPromptEleT << "\t NonPrompt: " << nNonPromptEleT<< "\t Fake: " << nFakeEleT <<"\t Eff: " << 1.*nPromptEleT/nTrueEle <<"\t FakeEff: " << 1.*(nFakeEleT+nNonPromptEleT)/(nFakeEle+nNonPromptEle) << endl;
 cout << "EffVReco: "<< 1.*nPromptEle/nPromptEle <<", "<< 1.*nPromptEleL/nPromptEle <<", "<< 1.*nPromptEleM/nPromptEle <<", "<< 1.*nPromptEleT/nPromptEle <<endl;
 cout << "FR_Fake: "<<1.*(nFakeEle)/(nFakeEle) <<", "<< 1.*(nFakeEleL)/(nFakeEle) <<", "<< 1.*(nFakeEleM)/(nFakeEle) <<", "<< 1.*(nFakeEleT)/(nFakeEle) <<endl;
 cout << "FR_NonPrompt: "<<1.*(nNonPromptEle)/(nNonPromptEle) <<", "<< 1.*(nNonPromptEleL)/(nNonPromptEle) <<", "<< 1.*(nNonPromptEleM)/(nNonPromptEle) <<", "<< 1.*(nNonPromptEleT)/(nNonPromptEle) <<endl;
 cout << "FR_Tot: "<<1.*(nFakeEle+nNonPromptEle)/(nFakeEle+nNonPromptEle) <<", "<< 1.*(nFakeEleL+nNonPromptEleL)/(nFakeEle+nNonPromptEle) <<", "<< 1.*(nFakeEleM+nNonPromptEleM)/(nFakeEle+nNonPromptEle) <<", "<< 1.*(nFakeEleT+nNonPromptEleT)/(nFakeEle+nNonPromptEle) <<endl;

 cout << " float EffVReco1015[3] = {"<<   1.*nPromptEleLBin[0]/nPromptEleBin[0] <<", "<< 1.*nPromptEleMBin[0]/nPromptEleBin[0] <<", "<< 1.*nPromptEleTBin[0]/nPromptEleBin[0] <<"}"<<endl;
 cout << " float FR_Fake1015[3] = {"<<     1.*(nFakeEleLBin[0])/(nFakeEleBin[0]) <<", "<< 1.*(nFakeEleMBin[0])/(nFakeEleBin[0]) <<", "<< 1.*(nFakeEleTBin[0])/(nFakeEleBin[0]) <<"}"<<endl;
 cout << " float EffVReco1015EE[3] = {"<< 1.*nPromptEleLBin[1]/nPromptEleBin[1] <<", "<< 1.*nPromptEleMBin[1]/nPromptEleBin[1] <<", "<< 1.*nPromptEleTBin[1]/nPromptEleBin[1] <<"}"<<endl;
 cout << " float FR_Fake1015EE[3] = {"<<   1.*(nFakeEleLBin[1])/(nFakeEleBin[1]) <<", "<< 1.*(nFakeEleMBin[1])/(nFakeEleBin[1]) <<", "<< 1.*(nFakeEleTBin[1])/(nFakeEleBin[1]) <<"}"<<endl;
 cout << " float EffVReco1520[3] = {"<<   1.*nPromptEleLBin[2]/nPromptEleBin[2] <<", "<< 1.*nPromptEleMBin[2]/nPromptEleBin[2] <<", "<< 1.*nPromptEleTBin[2]/nPromptEleBin[2] <<"}"<<endl;
 cout << " float FR_Fake1520[3] = {"<<     1.*(nFakeEleLBin[2])/(nFakeEleBin[2]) <<", "<< 1.*(nFakeEleMBin[2])/(nFakeEleBin[2]) <<", "<< 1.*(nFakeEleTBin[2])/(nFakeEleBin[2]) <<"}"<<endl;
 cout << " float EffVReco1520EE[3] = {"<< 1.*nPromptEleLBin[3]/nPromptEleBin[3] <<", "<< 1.*nPromptEleMBin[3]/nPromptEleBin[3] <<", "<< 1.*nPromptEleTBin[3]/nPromptEleBin[3] <<"}"<<endl;
 cout << " float FR_Fake1520EE[3] = {"<<   1.*(nFakeEleLBin[3])/(nFakeEleBin[3]) <<", "<< 1.*(nFakeEleMBin[3])/(nFakeEleBin[3]) <<", "<< 1.*(nFakeEleTBin[3])/(nFakeEleBin[3]) <<"}"<<endl;
 cout << " float EffVReco20[3] = {"<<     1.*nPromptEleLBin[4]/nPromptEleBin[4] <<", "<< 1.*nPromptEleMBin[4]/nPromptEleBin[4] <<", "<< 1.*nPromptEleTBin[4]/nPromptEleBin[4] <<"}"<<endl;
 cout << " float FR_Fake20[3] = {"<<       1.*(nFakeEleLBin[4])/(nFakeEleBin[4]) <<", "<< 1.*(nFakeEleMBin[4])/(nFakeEleBin[4]) <<", "<< 1.*(nFakeEleTBin[4])/(nFakeEleBin[4]) <<"}"<<endl;
 cout << " float EffVReco20EE[3] = {"<<   1.*nPromptEleLBin[5]/nPromptEleBin[5] <<", "<< 1.*nPromptEleMBin[5]/nPromptEleBin[5] <<", "<< 1.*nPromptEleTBin[5]/nPromptEleBin[5] <<"}"<<endl;
 cout << " float FR_Fake20EE[3] = {"<<     1.*(nFakeEleLBin[5])/(nFakeEleBin[5]) <<", "<< 1.*(nFakeEleMBin[5])/(nFakeEleBin[5]) <<", "<< 1.*(nFakeEleTBin[5])/(nFakeEleBin[5]) <<"}"<<endl;
cout << endl;



 outFile->cd(); // Make sure histograms get written out
 hSetEff["h_nvtx_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_nvtx_num_true"] , hSetEff["h_nvtx_den_true"] , "h_nvtx_eff", "Eff ;# vtxs");
 hSetEff["h_nvtx_fr"]  = (TH1F*) MakeEfficiencyPlot(hSetEff["h_nvtx_num_fake"] , hSetEff["h_nvtx_den_fake"] , "h_nvtx_fr", "Eff ;# vtxs");
 hSetEff["h_eta_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_eta_num_true"] , hSetEff["h_eta_den_true"] , "h_eta_eff", "Eff ;eta");
 hSetEff["h_eta_fr"]  = (TH1F*) MakeEfficiencyPlot(hSetEff["h_eta_num_fake"] , hSetEff["h_eta_den_fake"] , "h_eta_fr", "Eff ;eta");
 hSetEff["h_pt_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_pt_num_true"] , hSetEff["h_pt_den_true"] , "h_pt_eff", "Eff ;pt");
 hSetEff["h_pt_fr"]  = (TH1F*) MakeEfficiencyPlot(hSetEff["h_pt_num_fake"] , hSetEff["h_pt_den_fake"] , "h_pt_fr", "Eff ;pt");

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


  
  h_el_pt = new TH1F(Form("%s_el_pt",prefix.Data()),";electron p_{T} [GeV]",100,0.,100.);
  h_el_eta = new TH1F(Form("%s_el_eta",prefix.Data()),";electron eta",100,-3.,3.);
  h_el_truthmatch_DR = new TH1F(Form("%s_el_truthmatch_DR",prefix.Data()),";electron truth DR",500,0.,5.);
  h_el_lead_pt = new TH1F(Form("%s_el_lead_pt",prefix.Data()),";leading electron p_{T} [GeV]",100,0.,100.);
  h_el_subl_pt = new TH1F(Form("%s_el_subl_pt",prefix.Data()),";subleading electron p_{T} [GeV]",100,0.,100.);
  h_el_lead_eta = new TH1F(Form("%s_el_lead_eta",prefix.Data()),";leading electron #eta",100,-3.,3.);
  h_el_subl_eta = new TH1F(Form("%s_el_subl_eta",prefix.Data()),";subleading electron #eta",100,-3.,3.);
  h_ee_mll = new TH1F(Form("%s_ee_mll",prefix.Data()),";m_{ee} [GeV]",150,0.,150.);
  h_ee_denom_mll = new TH1F(Form("%s_ee_denom_mll",prefix.Data()),";m_{ee} [GeV]",150,0.,150.);
  h_ee_mediso_mll = new TH1F(Form("%s_ee_mediso_mll",prefix.Data()),";m_{ee} [GeV]",150,0.,150.);
      
  h_em_el_hlt_pt = new TH1F(Form("%s_em_el_hlt_pt",prefix.Data()),";HLT electron p_{T} [GeV]",100,0.,100.);
  h_em_mu_hlt_pt = new TH1F(Form("%s_em_mu_hlt_pt",prefix.Data()),";HLT muon p_{T} [GeV]",100,0.,100.);
  h_em_el_hlt_eta = new TH1F(Form("%s_em_el_hlt_eta",prefix.Data()),";HLT electron #eta",100,-3.,3.);
  h_em_mu_hlt_eta = new TH1F(Form("%s_em_mu_hlt_eta",prefix.Data()),";HLT muon #eta",100,-3.,3.);
  h_em_hlt_mll = new TH1F(Form("%s_em_hlt_mll",prefix.Data()),";HLT m_{e#mu} [GeV]",150,0.,150.);
  h_em_hlt_dr = new TH1F(Form("%s_em_hlt_dr",prefix.Data()),";HLT #DeltaR(e,#mu)",600,0.,6.);
    
  h_em_el_hlt_noreco_pt = new TH1F(Form("%s_em_el_hlt_noreco_pt",prefix.Data()),";HLT electron p_{T} [GeV]",100,0.,100.);
  h_em_mu_hlt_noreco_pt = new TH1F(Form("%s_em_mu_hlt_noreco_pt",prefix.Data()),";HLT muon p_{T} [GeV]",100,0.,100.);
  h_em_el_hlt_noreco_eta = new TH1F(Form("%s_em_el_hlt_noreco_eta",prefix.Data()),";HLT electron #eta",100,-3.,3.);
  h_em_mu_hlt_noreco_eta = new TH1F(Form("%s_em_mu_hlt_noreco_eta",prefix.Data()),";HLT muon #eta",100,-3.,3.);
  h_em_hlt_noreco_mll = new TH1F(Form("%s_em_hlt_noreco_mll",prefix.Data()),";HLT m_{e#mu} [GeV]",150,0.,150.);
  h_em_hlt_noreco_dr = new TH1F(Form("%s_em_hlt_noreco_dr",prefix.Data()),";HLT #DeltaR(e,#mu)",600,0.,6.);
    
  h_ee_events = new TH2F(Form("%s_ee_events",prefix.Data()),";leading e;subleading e",nbins_el,-0.5,nbins_el-0.5,nbins_el,-0.5,nbins_el-0.5);
   
  labelAxis(h_ee_events,0,11); labelAxis(h_ee_events,1,11); 
   
  bookElectronHistos(hSet1, "h_all");
  bookElectronHistos(hSet2, "h_true");bookElectronHistos(hSet2E, "h_trueEE");
  bookElectronHistos(hSet2f, "h_fake");  bookElectronHistos(hSet2Ef, "h_fakeEE");   bookElectronHistos(hSet2np, "h_np");  bookElectronHistos(hSet2Enp, "h_npEE"); 
  bookElectronHistos(hSet3, "h_WP80_true");  bookElectronHistos(hSet3E, "h_WP80_trueEE");
  bookElectronHistos(hSet3f, "h_WP80_fake");  bookElectronHistos(hSet3Ef, "h_WP80_fakeEE");   bookElectronHistos(hSet3np, "h_WP80_np");  bookElectronHistos(hSet3Enp, "h_WP80_npEE");
  bookElectronHistos(hSet4, "h_L_true");  bookElectronHistos(hSet4E, "h_L_trueEE");
  bookElectronHistos(hSet4f, "h_L_fake");  bookElectronHistos(hSet4Ef, "h_L_fakeEE");  bookElectronHistos(hSet4np, "h_L_np");  bookElectronHistos(hSet4Enp, "h_L_npEE");
  bookElectronHistos(hSet5, "h_M_true");  bookElectronHistos(hSet5E, "h_M_trueEE");
  bookElectronHistos(hSet5f, "h_M_fake");  bookElectronHistos(hSet5Ef, "h_M_fakeEE");  bookElectronHistos(hSet5np, "h_M_np");  bookElectronHistos(hSet5Enp, "h_M_npEE");
  bookElectronHistos(hSet6, "h_T_true");  bookElectronHistos(hSet6E, "h_T_trueEE");
  bookElectronHistos(hSet6f, "h_T_fake");  bookElectronHistos(hSet6Ef, "h_T_fakeEE");  bookElectronHistos(hSet6np, "h_T_np");  bookElectronHistos(hSet6Enp, "h_T_npEE");
  bookElectronHistos(hSet7, "h_Den_true");  bookElectronHistos(hSet7E, "h_Den_trueEE");
  bookElectronHistos(hSet7f, "h_Den_fake");  bookElectronHistos(hSet7Ef, "h_Den_fakeEE");  bookElectronHistos(hSet7np, "h_Den_np");  bookElectronHistos(hSet7Enp, "h_Den_npEE");

  //bookElectronHistos(hSetCut1, "h_cut1"); 
  //bookElectronHistos(hSetCut2, "h_cut2"); 
  //bookElectronHistos(hSetCut3, "h_cut3"); 
  //bookElectronHistos(hSetCut4, "h_cut4"); 
  //bookElectronHistos(hSetCut5, "h_cut5"); 
  //bookElectronHistos(hSetCut6, "h_cut6"); 
  //bookElectronHistos(hSetCut7, "h_cut7"); 
  //bookElectronHistos(hSetCut8, "h_cut8"); 
  //bookElectronHistos(hSetCut9, "h_cut9"); 

  hSetEff["h_nvtx_den_true"] = new TH1F(Form("%s_nvtx_den_true",prefix.Data()),";N(vtx)",max_nvtx/2,0,max_nvtx);
  hSetEff["h_nvtx_den_fake"] = new TH1F(Form("%s_nvtx_den_fake",prefix.Data()),";N(vtx)",max_nvtx/2,0,max_nvtx);
  hSetEff["h_nvtx_num_true"] = new TH1F(Form("%s_nvtx_num_true",prefix.Data()),";N(vtx)",max_nvtx/2,0,max_nvtx);
  hSetEff["h_nvtx_num_fake"] = new TH1F(Form("%s_nvtx_num_fake",prefix.Data()),";N(vtx)",max_nvtx/2,0,max_nvtx);
  hSetEff["h_eta_den_true"] = new TH1F(Form("%s_eta_den_true",prefix.Data()),";Eta",10, 0, 2.5);
  hSetEff["h_eta_den_fake"] = new TH1F(Form("%s_eta_den_fake",prefix.Data()),";Eta",10, 0, 2.5);
  hSetEff["h_eta_num_true"] = new TH1F(Form("%s_eta_num_true",prefix.Data()),";Eta",10, 0, 2.5);
  hSetEff["h_eta_num_fake"] = new TH1F(Form("%s_eta_num_fake",prefix.Data()),";Eta",10, 0, 2.5);
  hSetEff["h_pt_den_true"] = new TH1F(Form("%s_pt_den_true",prefix.Data()),";pt", 9, 10, 100);
  hSetEff["h_pt_den_fake"] = new TH1F(Form("%s_pt_den_fake",prefix.Data()),";pt", 9, 10, 100);
  hSetEff["h_pt_num_true"] = new TH1F(Form("%s_pt_num_true",prefix.Data()),";pt", 9, 10, 100);
  hSetEff["h_pt_num_fake"] = new TH1F(Form("%s_pt_num_fake",prefix.Data()),";pt", 9, 10, 100);

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

//--------------------------------------------------------------------

float dilepStudyLooper::electronPFiso(const unsigned int index, const bool cor) {
  
  float pt     = cms2.els_p4().at(index).pt();
  float etaAbs = fabs(cms2.els_etaSC().at(index));
  
  // get effective area
  float AEff = 0.;
  if (etaAbs <= 1.0) AEff = 0.10;
  else if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.12;
  else if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.085;
  else if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.11;
  else if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.12;
  else if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.12;
  else if (etaAbs > 2.4) AEff = 0.13;
  
//  float pfiso_ch = cms2.els_iso03_pf2012ext_ch().at(index); // CMSSW sequence run by CMS2
//  float pfiso_em = cms2.els_iso03_pf2012ext_em().at(index); // CMSSW sequence run by CMS2
//  float pfiso_nh = cms2.els_iso03_pf2012ext_nh().at(index); // CMSSW sequence run by CMS2
  float pfiso_ch = cms2.els_iso03_pf2012_ch().at(index); // Calculated by CMS2 ElectronMaker
  float pfiso_em = cms2.els_iso03_pf2012_em().at(index); // Calculated by CMS2 ElectronMaker
  float pfiso_nh = cms2.els_iso03_pf2012_nh().at(index); // Calculated by CMS2 ElectronMaker
//  float pfiso_ch = cms2.els_pfChargedHadronIso().at(index); // From AOD (wrong cone in 53X)
//  float pfiso_em = cms2.els_pfPhotonIso().at(index); // From AOD (wrong cone in 53X)
//  float pfiso_nh = cms2.els_pfNeutralHadronIso().at(index); // From AOD (wrong cone in 53X)
  
  // rho
  float rhoPrime = m_miniAOD ? std::max(cms2.evt_fixgrid_rho_ctr(), float(0.0)) : std::max(cms2.evt_kt6pf_foregiso_rho(), float(0.0));
  float pfiso_n = pfiso_em + pfiso_nh;
  if (cor) pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, float(0.0));
  float pfiso = (pfiso_ch + pfiso_n) / pt;
  
  return pfiso;
}

//--------------------------------------------------------------------

float dilepStudyLooper::muonPFiso(const unsigned int imu, const bool cor) {
  float chiso = cms2.mus_isoR03_pf_ChargedHadronPt().at(imu);
  float nhiso = cms2.mus_isoR03_pf_NeutralHadronEt().at(imu);
  float emiso = cms2.mus_isoR03_pf_PhotonEt().at(imu);
  float deltaBeta = cms2.mus_isoR03_pf_PUPt().at(imu);
  float pt = cms2.mus_p4().at(imu).pt();
  
  float absiso = chiso + nhiso + emiso;
  if (cor) absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return (absiso / pt);
  
}

//--------------------------------------------------------------------

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

//--------------------------------------------------------------------
// WW definitions - taken from HWW2012CORE
//--------------------------------------------------------------------

bool dilepStudyLooper::ElectronFOIdV4(unsigned int i) {
  
	float pt = cms2.els_p4().at(i).pt();
	float etaSC = cms2.els_etaSC().at(i);
  
	if (fabs(etaSC)<1.479) {
		if (cms2.els_sigmaIEtaIEta().at(i)>0.01		||
        fabs(cms2.els_dEtaIn().at(i))>0.007 	||
        fabs(cms2.els_dPhiIn().at(i))>0.15 		||
        cms2.els_hOverE().at(i)>0.12 			||
        cms2.els_tkIso().at(i)/pt>0.2 			||
        (cms2.els_ecalIso().at(i) - 1.0)/pt>0.2 ||
        cms2.els_hcalIso().at(i)/pt>0.2 ) return false;
	} else {
		if (cms2.els_sigmaIEtaIEta().at(i)>0.03		|| 
        fabs(cms2.els_dEtaIn().at(i))>0.009 	||
        fabs(cms2.els_dPhiIn().at(i))>0.10 		|| 
        cms2.els_hOverE().at(i)>0.10 			||
        cms2.els_tkIso().at(i)/pt>0.2 			||
        cms2.els_ecalIso().at(i)/pt>0.2 		||
        cms2.els_hcalIso().at(i)/pt>0.2 ) return false;
	}
  
  // MIT conversion
	// onlyAOD	if ( isFromConversionMIT(i) ) return false;
	// conversion rejection - hit based
	if ( cms2.els_exp_innerlayers().at(i) > 0 ) return false;
	
	return true;
} 

bool dilepStudyLooper::ElectronFOV4(unsigned int i){
  return ww_elBase(i) && ElectronFOIdV4(i) && ww_eld0PV(i) && ww_eldZPV(i);
}

bool dilepStudyLooper::ww_elBase(unsigned int index){
  if (cms2.els_p4().at(index).pt() < ptthresh_low) return false;
  if (fabs(cms2.els_p4().at(index).eta()) > 2.5) return false;
  return true;
}

double dilepStudyLooper::dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
//  cout<<"(vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt()"<<endl;
//  cout<<"("<<vtx.z()<<"-"<<pv.z()<<") - (("<<vtx.x()<<"-"<<pv.x()<<")*"<<p4.x()<<"+("<<vtx.y()<<"-"<<pv.y()<<")*"<<p4.y()<<")/"<<p4.pt()<<" * "<<p4.z()<<"/"<<p4.pt()<<endl;
//  cout<<(vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt()<<endl;
  return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}

bool dilepStudyLooper::ww_eld0PV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return false;
  double dxyPV = cms2.els_d0()[index]-
  cms2.vtxs_position()[vtxIndex].x()*sin(cms2.els_trk_p4()[index].phi())+
  cms2.vtxs_position()[vtxIndex].y()*cos(cms2.els_trk_p4()[index].phi());
  return fabs(dxyPV) < 0.02;
}

bool dilepStudyLooper::ww_eldZPV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return false;
  // double dzPV = cms2.els_z0corr()[index]-cms2.vtxs_position()[iMax].z();
  double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
  return fabs(dzpv)<0.1;
}

int dilepStudyLooper::primaryVertex() {
  return 0;
}


//--------------------------------------------------------------------
// object-trigger matching: taken from SingleLepton2012/looper/stopUtils.cc
//--------------------------------------------------------------------

int dilepStudyLooper::findTriggerIndex(const TString& trigName)
{
  vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
  if(found_it != end_it) return found_it - begin_it;
  return -1;
}

//--------------------------------------------------------------------

TString dilepStudyLooper::triggerName(const TString& triggerPattern){
  
  //-------------------------------------------------------
  // get exact trigger name corresponding to given pattern
  //-------------------------------------------------------
  
  bool    foundTrigger  = false;
  TString exact_hltname = "";
  
  for( unsigned int itrig = 0 ; itrig < hlt_trigNames().size() ; ++itrig ){
    if( TString( hlt_trigNames().at(itrig) ).Contains( triggerPattern ) ){
      foundTrigger  = true;
      exact_hltname = hlt_trigNames().at(itrig);
      break;
    }
  }
  
  if( !foundTrigger) return "TRIGGER_NOT_FOUND";
  
  return exact_hltname;
  
}


bool dilepStudyLooper::objectPassTrigger(const LorentzVector &obj, const TString& trigname, int type, float drmax ){
  
  if (type != 82 && type != 83) {
    cout << __FILE__ << " " << __LINE__ << " Error! invalid HLT object type: " << type << endl;
    return false;
  }
  
  TString exact_trigname = triggerName( trigname );
  
  if( exact_trigname.Contains("TRIGGER_NOT_FOUND") ){
    cout << __FILE__ << " " << __LINE__ << " Error! couldn't find trigger name " << trigname << endl;
    return false;
  }
  
  int trigidx = findTriggerIndex(exact_trigname);
  std::vector<LorentzVector> trigp4 = hlt_trigObjs_p4()[trigidx];
  std::vector<int> trigid = hlt_trigObjs_id()[trigidx];
  
  if( trigp4.size() == 0 ) return false;
  
  for (unsigned int i = 0; i < trigp4.size(); ++i){
    if (trigid.at(i) != type) continue;
    float dr = dRbetweenVectors(trigp4[i], obj);
    if( dr < drmax ) return true;
  }
  
  return false;
}

//----------------------------

float dilepStudyLooper::getdphi( float phi1 , float phi2 ){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

//-------------------------------------------------------------------

// based on muonIsoValuePF2012 from muonSelections.cc

isovals dilepStudyLooper::muonChIsoValuePF2012 (const unsigned int imu, const float R, const int ivtx)
{
  
  // isolation sums
  isovals vals;
  vals.chiso00 = 0.;
  vals.chiso04 = 0.;
  vals.chiso07 = 0.;
  vals.chiso10 = 0.;
  
  // loop on pfcandidates
  for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ++ipf) {
    
    // skip electrons and muons
    const int particleId = abs(cms2.pfcands_particleId()[ipf]);
    if (particleId == 11)    continue;
    if (particleId == 13)    continue;
    
    // deltaR between electron and cadidate
    const float dR = dRbetweenVectors(cms2.pfcands_p4()[ipf], cms2.mus_p4()[imu]);
    if (dR > R)              continue;
    
    // charged hadrons closest vertex
    // should be the primary vertex
    if (cms2.pfcands_charge().at(ipf) != 0) {
      //        if (particleId == 211 || particleId == 321 || particleId == 2212 || particleId == 999211) {
      // onlyAOD if (cms2.pfcands_vtxidx().at(ipf) != ivtx) continue;
      if (dR < 0.0001)
        continue;
      
      float pt = cms2.pfcands_p4()[ipf].pt();
      vals.chiso00 += pt;
      if (pt > 0.4) vals.chiso04 += pt;
      if (pt > 0.7) vals.chiso07 += pt;
      if (pt > 1.0) vals.chiso10 += pt;
    }
    
  } // loop on pf cands
  
  return vals;
}


void dilepStudyLooper::bookElectronHistos(std::map<std::string, TH1F*> & hSet, TString prefix ){

  TH1F * h_el_pt = new TH1F(Form("%s_el_pt",prefix.Data()),";electron pt",40,0.,100);
  TH1F * h_el_sieie = new TH1F(Form("%s_el_sieie",prefix.Data()),";electron sieie",100,0.,0.05);
  TH1F * h_el_dEtaIn = new TH1F(Form("%s_el_dEtaIn",prefix.Data()),";electron dEtaIn",100,-0.02,0.02);
  TH1F * h_el_dPhiIn = new TH1F(Form("%s_el_dPhiIn",prefix.Data()),";electron dPhiIn",100,-0.2,0.2);
  TH1F * h_el_hOverE = new TH1F(Form("%s_el_hOverE",prefix.Data()),";electron hOverE",100,0.,0.5);
  TH1F * h_el_d0corr = new TH1F(Form("%s_el_d0corr",prefix.Data()),";electron d0corr",100,-0.2,0.2);
  TH1F * h_el_z0corr = new TH1F(Form("%s_el_z0corr",prefix.Data()),";electron z0corr",100,-1.,1.);
  TH1F * h_el_ooemoop = new TH1F(Form("%s_el_ooemoop",prefix.Data()),";electron ooemoop",100,-0.5,0.5);
  TH1F * h_el_iso_cor = new TH1F(Form("%s_el_iso_cor",prefix.Data()),";electron iso_cor",100,0.,2.);
  TH1F * h_el_iso_corsize = new TH1F(Form("%s_el_iso_corzise",prefix.Data()),";electron iso_corsize",100,-1.,1.);
  TH1F * h_el_pfiso_chPU = new TH1F(Form("%s_el_pfiso_chPU",prefix.Data()),";electron pfiso_chPU",100,0.,2.);
  TH1F * h_el_pfiso_ch = new TH1F(Form("%s_el_pfiso_ch",prefix.Data()),";electron pfiso_ch",100,0.,2.);
  TH1F * h_el_pfiso_em = new TH1F(Form("%s_el_pfiso_em",prefix.Data()),";electron pfiso_em",100,0.,2.);
  TH1F * h_el_pfiso_nh = new TH1F(Form("%s_el_pfiso_nh",prefix.Data()),";electron pfiso_nh",100,0.,2.);
  TH1F * h_el_detiso_ch = new TH1F(Form("%s_el_detiso_ch",prefix.Data()),";electron detiso_ch",100,0.,2.);
  TH1F * h_el_detiso_em = new TH1F(Form("%s_el_detiso_em",prefix.Data()),";electron detiso_em",100,0.,2.);
  TH1F * h_el_detiso_nh = new TH1F(Form("%s_el_detiso_nh",prefix.Data()),";electron detiso_nh",100,0.,2.);
  TH1F * h_el_valid_pixelHits = new TH1F(Form("%s_el_valid_pixelHits",prefix.Data()),";electron valid_pixelHits",6,0,6);
  TH1F * h_el_lost_pixelhits = new TH1F(Form("%s_el_lost_pixelhits",prefix.Data()),";electron lost_pixelhits",4,0,4);
  TH1F * h_el_vtxFitConversion = new TH1F(Form("%s_el_vtxFitConversion",prefix.Data()),";electron vtxFitConversion",3,0,3);
  TH1F * h_el_chi2n = new TH1F(Form("%s_el_chi2n",prefix.Data()),";electron chi2n",100,0,10);
  TH1F * h_el_fbrem     = new TH1F(Form("%s_el_fbrem"    , prefix.Data()), ";electron fbrem    ", 100, -1.5, 1.5);
  TH1F * h_el_dEtaOut   = new TH1F(Form("%s_el_dEtaOut"  , prefix.Data()), ";electron dEtaOut  ", 100, -0.25, 0.25);
  TH1F * h_el_etaWidth  = new TH1F(Form("%s_el_etaWidth" , prefix.Data()), ";electron etaWidth ", 100, 0, 0.2);
  TH1F * h_el_phiWidth  = new TH1F(Form("%s_el_phiWidth" , prefix.Data()), ";electron phiWidth ", 100, 0, 0.4);
  TH1F * h_el_e1x5e5x5  = new TH1F(Form("%s_el_e1x5e5x5" , prefix.Data()), ";electron e1x5e5x5 ", 100, 0, 1.2);
  TH1F * h_el_r9        = new TH1F(Form("%s_el_r9"       , prefix.Data()), ";electron r9       ", 100, 0, 1.2);
  TH1F * h_el_sieieSC   = new TH1F(Form("%s_el_sieieSC"  , prefix.Data()), ";electron sieieSC  ", 100, 0, 0.2);
  TH1F * h_el_eoverpIn  = new TH1F(Form("%s_el_eoverpIn" , prefix.Data()), ";electron eoverpIn ", 100, 0, 20);
  TH1F * h_el_eoverpOut = new TH1F(Form("%s_el_eoverpOut", prefix.Data()), ";electron eoverpOut", 100, 0, 20);
  TH1F * h_el_psOverRaw = new TH1F(Form("%s_el_psOverRaw", prefix.Data()), ";electron psOverRaw", 100, 0, 0.5);
  TH1F * h_el_seed      = new TH1F(Form("%s_el_seed"     , prefix.Data()), ";electron seed     ", 5, 0, 5);
  TH1F * h_el_ncluster  = new TH1F(Form("%s_el_ncluster" , prefix.Data()), ";electron ncluster ", 15, 0, 15);
  hSet.insert(std::pair<std::string, TH1F*> ("pt",h_el_pt));
  hSet.insert(std::pair<std::string, TH1F*> ("sieie",h_el_sieie));
  hSet.insert(std::pair<std::string, TH1F*> ("dEtaIn",h_el_dEtaIn));   
  hSet.insert(std::pair<std::string, TH1F*> ("dPhiIn",h_el_dPhiIn));   
  hSet.insert(std::pair<std::string, TH1F*> ("hOverE",h_el_hOverE));
  hSet.insert(std::pair<std::string, TH1F*> ("d0corr",h_el_d0corr));
  hSet.insert(std::pair<std::string, TH1F*> ("z0corr",h_el_z0corr));
  hSet.insert(std::pair<std::string, TH1F*> ("ooemoop",h_el_ooemoop));
  hSet.insert(std::pair<std::string, TH1F*> ("iso_cor",h_el_iso_cor));
  hSet.insert(std::pair<std::string, TH1F*> ("iso_corsize",h_el_iso_corsize));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_chPU",h_el_pfiso_chPU));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_ch",h_el_pfiso_ch));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_em",h_el_pfiso_em));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_nh",h_el_pfiso_nh));
  hSet.insert(std::pair<std::string, TH1F*> ("detiso_ch",h_el_detiso_ch));
  hSet.insert(std::pair<std::string, TH1F*> ("detiso_em",h_el_detiso_em)); 
  hSet.insert(std::pair<std::string, TH1F*> ("detiso_nh",h_el_detiso_nh));
  hSet.insert(std::pair<std::string, TH1F*> ("valid_pixelHits",h_el_valid_pixelHits));
  hSet.insert(std::pair<std::string, TH1F*> ("lost_pixelhits",h_el_lost_pixelhits));
  hSet.insert(std::pair<std::string, TH1F*> ("vtxFitConversion",h_el_vtxFitConversion));
  hSet.insert(std::pair<std::string, TH1F*> ("chi2n",h_el_chi2n));
  hSet.insert(std::pair<std::string, TH1F*> ("fbrem"     , h_el_fbrem      ));
  hSet.insert(std::pair<std::string, TH1F*> ("dEtaOut"   , h_el_dEtaOut    ));
  hSet.insert(std::pair<std::string, TH1F*> ("etaWidth"  , h_el_etaWidth   ));
  hSet.insert(std::pair<std::string, TH1F*> ("phiWidth"  , h_el_phiWidth   ));
  hSet.insert(std::pair<std::string, TH1F*> ("e1x5e5x5"  , h_el_e1x5e5x5   ));
  hSet.insert(std::pair<std::string, TH1F*> ("r9"        , h_el_r9         ));
  hSet.insert(std::pair<std::string, TH1F*> ("sieieSC"   , h_el_sieieSC    ));
  hSet.insert(std::pair<std::string, TH1F*> ("eoverpIn"  , h_el_eoverpIn   ));
  hSet.insert(std::pair<std::string, TH1F*> ("eoverpOut" , h_el_eoverpOut  ));
  hSet.insert(std::pair<std::string, TH1F*> ("psOverRaw" , h_el_psOverRaw  ));
  hSet.insert(std::pair<std::string, TH1F*> ("seed"      , h_el_seed       ));
  hSet.insert(std::pair<std::string, TH1F*> ("ncluster"  , h_el_ncluster   ));
  return;
}


void dilepStudyLooper::fillElectronQuantities(std::map<std::string, TH1F*> & hSet, electron e ){

  fillUnderOverFlow( hSet["pt"]           , e.pt, 1);
  fillUnderOverFlow( hSet["sieie"]           , e.sieie, 1);
  fillUnderOverFlow( hSet["dEtaIn"]          , e.dEtaIn, 1);
  fillUnderOverFlow( hSet["dPhiIn"]          , e.dPhiIn, 1);
  fillUnderOverFlow( hSet["hOverE"]          , e.hOverE, 1);
  fillUnderOverFlow( hSet["d0corr"]          , e.d0corr, 1);
  fillUnderOverFlow( hSet["z0corr"]          , e.z0corr, 1);
  fillUnderOverFlow( hSet["ooemoop"]         , e.ooemoop, 1);
  fillUnderOverFlow( hSet["iso_cor"]         , e.iso_cor, 1);
  fillUnderOverFlow( hSet["iso_corsize"]     , e.iso_cor - e.iso_uncor, 1);
  fillUnderOverFlow( hSet["pfiso_chPU"]       , e.pfiso_chPU, 1);
  fillUnderOverFlow( hSet["pfiso_ch"]         , e.pfiso_ch, 1);
  fillUnderOverFlow( hSet["pfiso_em"]         , e.pfiso_em, 1);
  fillUnderOverFlow( hSet["pfiso_nh"]         , e.pfiso_nh, 1);
  fillUnderOverFlow( hSet["detiso_ch"]         , e.detiso_ch, 1);
  fillUnderOverFlow( hSet["detiso_em"]         , e.detiso_em, 1);
  fillUnderOverFlow( hSet["detiso_nh"]         , e.detiso_nh, 1);
  fillUnderOverFlow( hSet["valid_pixelHits"] , e.valid_pixelHits, 1);
  fillUnderOverFlow( hSet["lost_pixelhits"]  , e.lost_pixelhits, 1);
  fillUnderOverFlow( hSet["vtxFitConversion"], e.vtxFitConversion, 1);  
  fillUnderOverFlow( hSet["chi2n"], e.chi2n, 1);  
  fillUnderOverFlow( hSet["fbrem"]        , e.fbrem     , 1);
  fillUnderOverFlow( hSet["dEtaOut"]      , e.dEtaOut   , 1);
  fillUnderOverFlow( hSet["etaWidth"]     , e.etaWidth  , 1);
  fillUnderOverFlow( hSet["phiWidth"]     , e.phiWidth  , 1);
  fillUnderOverFlow( hSet["e1x5e5x5"]     , e.e1x5e5x5  , 1);
  fillUnderOverFlow( hSet["r9"]	          , e.r9	      , 1);
  fillUnderOverFlow( hSet["sieieSC"]      , e.sieieSC   , 1);
  fillUnderOverFlow( hSet["eoverpIn"]     , e.eoverpIn  , 1);
  fillUnderOverFlow( hSet["eoverpOut"]    , e.eoverpOut , 1);
  fillUnderOverFlow( hSet["psOverRaw"]    , e.psOverRaw , 1);
  fillUnderOverFlow( hSet["seed"]         , e.seed      , 1);
  fillUnderOverFlow( hSet["ncluster"]     , e.ncluster  , 1);	  



//  hSet["sieie"]->Fill(e.sieie, 1);
//  hSet["dEtaIn"]->Fill(e.dEtaIn, 1);
//  hSet["dPhiIn"]->Fill(e.dPhiIn, 1);
//  hSet["hOverE"]->Fill(e.hOverE, 1);
//  hSet["d0corr"]->Fill(e.d0corr, 1);
//  hSet["z0corr"]->Fill(e.z0corr, 1);
//  hSet["ooemoop"]->Fill(e.ooemoop, 1);
//  hSet["iso_cor"]->Fill(e.iso_cor, 1);
//  hSet["valid_pixelHits"]->Fill(e.valid_pixelHits, 1);
//  hSet["lost_pixelhits"]->Fill(e.lost_pixelhits, 1);
//  hSet["vtxFitConversion"]->Fill(e.vtxFitConversion, 1);  
  
  return;
}

void dilepStudyLooper::fillElectronQuantitiesN1(std::map<std::string, TH1F*> & hSet, electron e, electron cut){
  ULong64_t p = 0;
  fillElectronCutsResult ( e, cut, p );
    //cout<<(bitset<64>) p <<endl;
  ULong64_t all = 0;
  for (uint i = 0; i <= 17; i++) all |= 1ll<<i;
				   //(1ll<<0 | 1ll<<1 | 1ll<<2 | 1ll<<3 | 1ll<<4 | 1ll<<5 | 1ll<<6 | 1ll<<7 | 1ll<<8 | 1ll<<9 | 1ll<<10);
  int i = 0;
  
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["sieie"]->Fill(e.sieie, 1);                         i++;  
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["dEtaIn"]->Fill(e.dEtaIn, 1);			   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["dPhiIn"]->Fill(e.dPhiIn, 1);			   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["hOverE"]->Fill(e.hOverE, 1);			   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["d0corr"]->Fill(e.d0corr, 1);			   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["z0corr"]->Fill(e.z0corr, 1);			   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["ooemoop"]->Fill(e.ooemoop, 1);			   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["iso_cor"]->Fill(e.iso_cor, 1);			   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["pfiso_ch"]->Fill(e.pfiso_ch, 1);		       i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["pfiso_em"]->Fill(e.pfiso_em, 1);		       i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["pfiso_nh"]->Fill(e.pfiso_nh, 1);		      i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["detiso_ch"]->Fill(e.detiso_ch, 1);		       i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["detiso_em"]->Fill(e.detiso_em, 1);		       i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["detiso_nh"]->Fill(e.detiso_nh, 1);		      i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["valid_pixelHits"]->Fill(e.valid_pixelHits, 1);	   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["lost_pixelhits"]->Fill(e.lost_pixelhits, 1);	   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["vtxFitConversion"]->Fill(e.vtxFitConversion, 1);   i++;
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) hSet["chi2n"]->Fill(e.chi2n, 1);                         i++;
  if ( (p & all) == all ) {
   fillUnderOverFlow( hSet["pt"]           , e.pt, 1);
    fillUnderOverFlow( hSet["fbrem"]        , e.fbrem     , 1);
    fillUnderOverFlow( hSet["dEtaOut"]      , e.dEtaOut   , 1);
    fillUnderOverFlow( hSet["etaWidth"]     , e.etaWidth  , 1);
    fillUnderOverFlow( hSet["phiWidth"]     , e.phiWidth  , 1);
    fillUnderOverFlow( hSet["e1x5e5x5"]     , e.e1x5e5x5  , 1);
    fillUnderOverFlow( hSet["r9"]	          , e.r9	      , 1);
    fillUnderOverFlow( hSet["sieieSC"]      , e.sieieSC   , 1);
    fillUnderOverFlow( hSet["eoverpIn"]     , e.eoverpIn  , 1);
    fillUnderOverFlow( hSet["eoverpOut"]    , e.eoverpOut , 1);
    fillUnderOverFlow( hSet["psOverRaw"]    , e.psOverRaw , 1);
    fillUnderOverFlow( hSet["seed"]         , e.seed      , 1);
    fillUnderOverFlow( hSet["ncluster"]     , e.ncluster  , 1);	  
    fillUnderOverFlow( hSet["pfiso_chPU"]     , e.pfiso_chPU , 1);
  }


  return;
}

void dilepStudyLooper::fillElectronCutsResult(electron e, electron cut, ULong64_t & pass) {
  int i=0;
  if ( e.sieie < cut.sieie ) pass   |= 1ll<<i; i++;
  if ( fabs(e.dEtaIn) < cut.dEtaIn ) pass |= 1ll<<i; i++;
  if ( fabs(e.dPhiIn) < cut.dPhiIn ) pass |= 1ll<<i; i++;
  if ( e.hOverE < cut.hOverE ) pass |= 1ll<<i; i++;
  if ( fabs(e.d0corr) < cut.d0corr ) pass |= 1ll<<i; i++;
  if ( fabs(e.z0corr) < cut.z0corr ) pass |= 1ll<<i; i++;
  if ( fabs(e.ooemoop) < cut.ooemoop ) pass |=   1ll<<i; i++;
  float isoCut = cut.iso_cor;
  if ( e.pt < 20 && fabs(e.eta) > 1.479 && cut.iso_cor == 0.10) isoCut = 0.07;
  if ( e.pt < 20 && fabs(e.eta) > 1.479 && cut.iso_cor == 0.15) isoCut = 0.10; 
  if ( e.iso_cor < isoCut  ) pass |=   1ll<<i; i++;
  if ( e.pfiso_ch < cut.pfiso_ch ) pass |= 1ll<<i; i++;
  if ( e.pfiso_em < cut.pfiso_em ) pass |= 1ll<<i; i++;
  if ( e.pfiso_nh < cut.pfiso_nh ) pass |= 1ll<<i; i++;
  if ( e.detiso_ch < cut.detiso_ch ) pass |= 1ll<<i; i++;
  if ( e.detiso_em < cut.detiso_em ) pass |= 1ll<<i; i++;
  if ( e.detiso_nh < cut.detiso_nh ) pass |= 1ll<<i; i++;
  if ( e.valid_pixelHits >= cut.valid_pixelHits ) pass |= 1ll<<i; i++;
  if ( e.lost_pixelhits <= cut.lost_pixelhits ) pass |= 1ll<<i; i++;
  bool passVtx = cut.vtxFitConversion ? !e.vtxFitConversion : true;
  if (passVtx) pass |= 1ll<<i;  i++;
  if ( e.chi2n <= cut.chi2n ) pass |= 1ll<<i;
  return;
}



bool dilepStudyLooper::passElectronCuts(electron eleStruct, electron cut, electron cutEE) {
  ULong64_t p = 0;
  ULong64_t pEE = 0;
  ULong64_t all1 = 0;
  bool pass = false;
  for (uint i = 0; i <= 17; i++) all1 |= 1ll<<i; //this is  (ULong64_t) 262143, for reference
  fillElectronCutsResult ( eleStruct, cut, p );
  fillElectronCutsResult ( eleStruct, cutEE, pEE );
  if ( (fabs(eleStruct.eta) <= 1.479 && (p & all1) == all1) || (fabs(eleStruct.eta) > 1.479 && (pEE & all1) == all1) ) pass = true;
  return pass;
}
void dilepStudyLooper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

void  dilepStudyLooper::electronPFiso2(float &pfiso_ch, float &pfiso_em, float &pfiso_nh, float &pfiso_chPU, const float R, const unsigned int iel, const int ivtx, bool useMap, int useDeltaBetaWeights) {
  // Mostly taken from electronIsoValuePF2012 in electronSelections.cc

  pfiso_ch = 0;
  pfiso_em = 0;
  pfiso_nh = 0;
  pfiso_chPU = 0;

  for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ++ipf) {
    
    // skip electrons and muons
    const int particleId = abs(cms2.pfcands_particleId()[ipf]);
    if (particleId == 11)    continue;
    if (particleId == 13)    continue;
//    if ( cms2.pfcands_mva_nothing_gamma()[ipf] > 0.99 && cand.superCluster().isNonnull() 
//	  && cms2.pfcands_superClusterRef().isNonnull() 
//	  && cand.superCluster() == it->superClusterRef())   continue; // Remove photons with same supercluster (can't do with CMS2 ntuples)
//    if (particleId == 22 && cms2.pfcands_mva_nothing_gamma()[ipf] > 0) continue; // Shortcut to remove photons (can't do with Slim CMS2 ntuples)
    
    // deltaR between electron and cadidate
    const float dR = dRbetweenVectors(cms2.pfcands_p4()[ipf], cms2.els_p4().at(iel));
    if (dR > R)              continue;

    // Endcap Vetoes
    if (!(cms2.els_fiduciality()[iel] & (1<<ISEB)) && !useMap) { // only apply these vetoes if we don't use the PFCand Map
      if (particleId == 211 && dR <= 0.015)   continue;
      if (particleId == 22  && dR <= 0.08)    continue;
    } 
    
    // Charged hadron vertex matching
    bool isPU = false;
    if (particleId == 211) {
      if (m_miniAOD) {
	if ( cms2.pfcands_fromPV()[ipf] < 2) isPU = true;
      }
      else {
	int pfVertexIndex = cms2.pfcands_vtxidx().at(ipf); 
	if (pfVertexIndex != ivtx) isPU = true;
      }
//      cout<<"charghed pfcand with fromPV "<< (int) cms2.pfcands_fromPV()[ipf] ;
//      cout<<". <1 ? "<< (cms2.pfcands_fromPV()[ipf] < 1);
//      cout<<". <2 ? "<< (cms2.pfcands_fromPV()[ipf] < 2);
//      cout<<". <3 ? "<< (cms2.pfcands_fromPV()[ipf] < 3)<<endl;
      
    }

    // Use Map (only REL 7)
    bool skipcand = false;
    if (useMap) {
      std::vector<int> pfCandMatchedToEle = els_PFCand_idx().at(iel);
      for ( unsigned int i = 0; i < pfCandMatchedToEle.size(); i++ ) {
	if (pfCandMatchedToEle[i] == (int) ipf) skipcand = true;
      }
    }
    if (skipcand) continue;

    LorentzVector & p4 = cms2.pfcands_p4()[ipf];
    float weight = 1.0;
    if (useDeltaBetaWeights>0 && (particleId == 22 || particleId == 130) ) {      
      float sumPU = useDeltaBetaWeights == 2 ? 1. : 0;
      float sumNPU = useDeltaBetaWeights == 2 ? 1. : 0;
      float tempSum = 0;
      //      cout<<"Found neutral PFCand with pT "<<p4.pt()<<endl;
      for (unsigned int jpf = 0; jpf < cms2.pfcands_p4().size(); ++jpf) { // loop over PFCands
	if (abs(cms2.pfcands_particleId()[jpf]) != 211)  continue; // only keep charged
	LorentzVector & p4_2 = cms2.pfcands_p4()[jpf];
	if ( fabs(p4.Eta() - p4_2.Eta()) > 0.5 ) continue; // only keep close-by (avoid DR calculation with all the event!)
	if ( std::min(::fabs(p4.Phi() - p4_2.Phi()), 2 * M_PI - fabs(p4.Phi() - p4_2.Phi())) > 0.5 ) continue; 
	float dR2 = dRbetweenVectors2(p4, p4_2);
	if (dR2 > 0.25) continue; // use a circle instead of square
	//	cout<<"Nearby PFCand at dR "<<dR2<<" is ";
	bool fromPV = true;
	if (m_miniAOD) {
	  if (cms2.pfcands_fromPV().at(jpf) < 2 ) fromPV = false;
	} 
	else {
	  if (cms2.pfcands_vtxidx().at(jpf) == ivtx ) fromPV = false;
	}

	if (useDeltaBetaWeights==1) { // 1/DR^2 metric
	  if (fromPV ) {sumNPU += 1./dR2; /*cout<<"NPU"<<endl;*/} 
	  else {sumPU += 1./dR2; /*cout<<"PU"<<endl;*/}
	}
	if (useDeltaBetaWeights==2) { // log(pt/DR) metric
	  tempSum= 1.*p4_2.pt()*p4_2.pt()/dR2; //cout<<"tempSum "<<tempSum<<endl;
	  if (fromPV  && tempSum > 0) {sumNPU *= tempSum; /*cout<<"NPU"<<endl;*/}
	  else if (tempSum > 0) {sumPU *= tempSum; /*cout<<"PU"<<endl;*/}
	}	
      }
      //cout<<"sumPU "<<sumPU<<", sumNPU "<<sumNPU<<endl;
      if (useDeltaBetaWeights==2) {sumPU=0.5*log(sumPU); sumNPU=0.5*log(sumNPU);  /*cout<<"after log: sumPU "<<sumPU<<", sumNPU "<<sumNPU<<endl;*/}
      weight = 1.0 * sumNPU / (sumNPU+sumPU);
      //cout<<"Weight is "<<weight<<endl;
    }
    
    
    if (particleId == 211 && !isPU)      pfiso_ch += p4.pt();
    if (particleId == 211 &&  isPU)      pfiso_chPU += p4.pt();
    if (particleId == 22)       pfiso_em += p4.pt() * weight;
    if (particleId == 130)      pfiso_nh += p4.pt() * weight;


  }// end loop on the candidates

  return;
}

void dilepStudyLooper::fillElectronStructure( const unsigned int iel, electron & eleStruct, bool useMap, bool doPFCandLoop, int useDeltaBetaWeights, bool DeltaBetaSimple, bool areaCorrection, bool fillPFiso) {
	float ch = 0;
	float em = 0;
	float nh = 0;
	float chPU = 0;
	float pt = els_p4().at(iel).pt();
	eleStruct.pt = els_p4().at(iel).pt() ;
        eleStruct.eta = els_etaSC().at(iel);
        eleStruct.sieie = m_miniAOD ? els_sigmaIEtaIEta_full5x5().at(iel) : els_sigmaIEtaIEta().at(iel);
        eleStruct.dEtaIn = els_dEtaIn().at(iel);
        eleStruct.dPhiIn = els_dPhiIn().at(iel);
        eleStruct.hOverE = els_hOverE().at(iel);
        eleStruct.d0corr = els_d0corr().at(iel); // for miniAOD move to els_dxyPV
        eleStruct.z0corr = dzPV(els_vertex_p4()[iel], els_trk_p4()[iel], vtxs_position()[0]);//els_z0corr().at(iel); // for miniAOD move to els_dzPV
        eleStruct.ooemoop = (1.0/els_ecalEnergy().at(iel)) - (els_eOverPIn().at(iel)/els_ecalEnergy().at(iel)) ;
	eleStruct.iso_uncor = (fillPFiso && !doPFCandLoop) ? electronPFiso(iel, false) : 0;
	if (fillPFiso) {
	  if ( doPFCandLoop ) electronPFiso2(ch, em, nh, chPU, 0.3, iel, firstGoodVertex(), useMap, useDeltaBetaWeights); 	
	  if (!doPFCandLoop) eleStruct.iso_cor = electronPFiso(iel, areaCorrection);// false = NO PILEUP AREA CORRECTION
	  else {
	    if (DeltaBetaSimple) eleStruct.iso_cor = ( ch + std::max(0.0, nh + em - 0.5 * chPU) ) / pt;
	    else eleStruct.iso_cor = (ch + em + nh) / pt; 
	    //	cout<<"eleStruct.iso_cor "<<eleStruct.iso_cor<<" and pt "<<pt<<endl;
	  }	
	}

//	eleStruct.pfiso_ch = prefix.Contains("700") ? ch/pt : cms2.els_iso03_pf2012ext_ch().at(iel) / pt; //els_pfChargedHadronIso().at(iel) / pt;
//	eleStruct.pfiso_em = prefix.Contains("700") ? em/pt : cms2.els_iso03_pf2012ext_em().at(iel) / pt; //els_pfNeutralHadronIso().at(iel) / pt;
//	eleStruct.pfiso_nh = prefix.Contains("700") ? nh/pt : cms2.els_iso03_pf2012ext_nh().at(iel) / pt; //els_pfPhotonIso().at(iel) / pt;
	eleStruct.pfiso_chPU = chPU/pt;  //ComponentsOfPFIso
	eleStruct.pfiso_ch = ch/pt;  //ComponentsOfPFIso
	eleStruct.pfiso_em = em/pt;  //ComponentsOfPFIso
	eleStruct.pfiso_nh = nh/pt;  //ComponentsOfPFIso
	eleStruct.detiso_ch = els_tkIso().at(iel)/pt;    // Detector Iso
	eleStruct.detiso_em = els_ecalIso().at(iel)/pt;  // Detector Iso
	eleStruct.detiso_nh = els_hcalIso().at(iel)/pt;  // Detector Iso
        eleStruct.valid_pixelHits = els_valid_pixelhits().at(iel);
        eleStruct.lost_pixelhits = els_lost_pixelhits().at(iel);
        eleStruct.vtxFitConversion = m_miniAOD ? !els_conv_vtx_flag().at(iel) : isMITConversion(iel, 0,   1e-6,   2.0,   true,  false);
	eleStruct.chi2n = els_chi2().at(iel) / els_ndof().at(iel);
	eleStruct.fbrem = els_fbrem().at(iel);
	eleStruct.dEtaOut = els_dEtaOut().at(iel);
	eleStruct.etaWidth = els_etaSCwidth().at(iel);
	eleStruct.phiWidth = els_phiSCwidth().at(iel);
	eleStruct.e1x5e5x5 =  (els_e5x5().at(iel)) !=0. ? 1.-(els_e1x5().at(iel) / els_e5x5().at(iel)) : -1. ;
	eleStruct.r9 = els_r9().at(iel); // e3x3(seed) / ele.superCluster()->rawEnergy();
	eleStruct.sieieSC = m_miniAOD ? 0 : els_sigmaIEtaIEtaSC().at(iel);
	eleStruct.eoverpIn = els_eOverPIn().at(iel);
	eleStruct.eoverpOut = els_eOverPOut().at(iel);

}


