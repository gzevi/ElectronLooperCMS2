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
  int nPromptEleVBin[6] = {};
  int nPromptEleLBin[6] = {};
  int nPromptEleMBin[6] = {};
  int nPromptEleTBin[6] = {};
  int nFakeEleBin[6] = {};
  int nFakeEleVBin[6] = {};
  int nFakeEleLBin[6] = {};
  int nFakeEleMBin[6] = {};
  int nFakeEleTBin[6] = {};
  int nNPEleBin[6] = {};
  int nNPEleVBin[6] = {};
  int nNPEleLBin[6] = {};
  int nNPEleMBin[6] = {};
  int nNPEleTBin[6] = {};
  int nPromptEleVBinNoIso[6] = {};
  int nPromptEleLBinNoIso[6] = {};
  int nPromptEleMBinNoIso[6] = {};
  int nPromptEleTBinNoIso[6] = {};
  int nFakeEleVBinNoIso[6] = {};
  int nFakeEleLBinNoIso[6] = {};
  int nFakeEleMBinNoIso[6] = {};
  int nFakeEleTBinNoIso[6] = {};
  int nNPEleVBinNoIso[6] = {};
  int nNPEleLBinNoIso[6] = {};
  int nNPEleMBinNoIso[6] = {};
  int nNPEleTBinNoIso[6] = {};


  int CountVstdIDstdISO[6][3] = {};//stdIDstdISO // 6 regions (pt/eta), 3 types (truth, nonprompt, fake)
  int CountLstdIDstdISO[6][3] = {};
  int CountMstdIDstdISO[6][3] = {};
  int CountTstdIDstdISO[6][3] = {};
  int CountVnewIDstdISO[6][3] = {};//newIDstdISO
  int CountLnewIDstdISO[6][3] = {};
  int CountMnewIDstdISO[6][3] = {};
  int CountTnewIDstdISO[6][3] = {};
  int CountVnewIDnewISO[6][3] = {};//newIDnewISO
  int CountLnewIDnewISO[6][3] = {};
  int CountMnewIDnewISO[6][3] = {};
  int CountTnewIDnewISO[6][3] = {};
  int CountVstdIDnoISO[6][3] = {};//stdIDnoISO
  int CountLstdIDnoISO[6][3] = {};
  int CountMstdIDnoISO[6][3] = {};
  int CountTstdIDnoISO[6][3] = {};
  int CountVnewIDnoISO[6][3] = {};//newIDnoISO  
  int CountLnewIDnoISO[6][3] = {};
  int CountMnewIDnoISO[6][3] = {};
  int CountTnewIDnoISO[6][3] = {};


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
      if (filename.Contains("CMS3") || filename.Contains("miniAOD") || prefix.Contains("CMS3")) m_miniAOD = true;
      else m_miniAOD = false;

      if (prefix.Contains("20bx25")) is25ns = true;

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
  
      dilepStudyLooper::electron eleVcuts; //VETO
      eleVcuts.sieie = 0.01;
      eleVcuts.dEtaIn =0.007;
      eleVcuts.dPhiIn = 0.8; //wow
      eleVcuts.hOverE = 0.15;
      eleVcuts.d0corr = 0.04;
      eleVcuts.z0corr = 0.2;
      eleVcuts.ooemoop = 999; //wow
      eleVcuts.iso_cor = 0.15;
      eleVcuts.pfiso_ch = 999;
      eleVcuts.pfiso_em = 999;
      eleVcuts.pfiso_nh = 999;
      eleVcuts.detiso_ch = 999;
      eleVcuts.detiso_em = 999;
      eleVcuts.detiso_nh = 999;
      eleVcuts.valid_pixelHits = 0;
      eleVcuts.lost_pixelhits = 3;
      eleVcuts.vtxFitConversion = false; // false means don't check
      eleVcuts.chi2n = 999;

      dilepStudyLooper::electron eleVcutsEE = eleVcuts; //VETO
      eleVcutsEE.sieie = 0.03;
      eleVcutsEE.dEtaIn = 0.01;
      eleVcutsEE.dPhiIn = 0.7; //wow
      eleVcutsEE.hOverE = 999; //wow
      eleVcutsEE.iso_cor = 0.15;//0.15; // should be 0.10 for <20 GeV    
      
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
      eleTcuts.dEtaIn =0.004;
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
        if (!isData && (pt < 20 || pt > 50)) continue;
        if (isData && (pt < 10 || pt > 35)) continue;
	//if (pt < 20) continue;
	//	if (pt < 10 || pt > 20) continue;
	if (fabs(eta) > 2.5) continue;
	//dilepStudyLooper::electron eleStructTemp; 
	//fillElectronStructure(iel, eleStructTemp,  false, false, 0, false, false, false);
	//if ( !passElectronCuts(eleStructTemp, eleDencuts, eleDencutsEE)) continue;
	// Quick check for single electron:
	int neleDen = 1;
	for (unsigned int jel = 0; jel < els_p4().size(); ++jel) {
	  if (!isData) continue;
	  if (iel==jel) continue;
	  if (els_p4().at(iel).pt() < 10 ) continue;
	  dilepStudyLooper::electron eleStruct2; 
	  fillElectronStructure(jel, eleStruct2,  false, false, 0, false, false, false, false);
	  if (passElectronCuts(eleStruct2, eleDencuts, eleDencutsEE)) neleDen++;
	  //cout<<"neleDen = "<<neleDen<<endl;
	}
	if (neleDen > 1) continue;

	//if (eta < 0 ) continue;
	//	if (isData && evt_pfmet() > 20) continue;	
	//	if (isData && Mt( el_p4, evt_pfmet(), evt_pfmetphi() ) > 25) continue; // need to calculate mt	

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


        // Truth matching: t=true (status 3), f=fake (including non-prompt)
        // Can't separate non-prompt because we only have status 3 particles in CMS2 ntuple
        bool truthmatched = false;
	bool fakematched = false;
	bool nonpromptmatched = false;
	float status1E = 0;
	float status3E = 0;
	float truthE = 0;
	if (!isData) {
	  float bestmatch = 999;
	  int bestmatchMother = 999;
	  //cout<<"Look for truth for electron with pt "<<pt<<endl;
	  for(unsigned int idx = 0; idx < genps_id().size(); idx++) {	    
	    int pid = abs(genps_id().at(idx));    
	    if(pid != 11) continue;      
	    float dr = dRbetweenVectors(genps_p4().at(idx), el_p4);
	    //if ( dr < 0.3 ) cout<<"Found electron in record with mother "<<genps_id_mother().at(idx)<<" and dr "<<dr<<endl;
	    if(genps_status().at(idx) != 3 && genps_status().at(idx) != 22 && genps_status().at(idx) != 23  && abs(genps_id_mother().at(idx)) != 24) continue; // trying everything
	    
	    h_el_truthmatch_DR->Fill(dr, 1);
	    if ( dr > 0.02 ) {
	      if (dr < bestmatch ) {bestmatch = dr; bestmatchMother=genps_id_mother().at(idx);}
	      continue;
	    }
	    else {
	      truthmatched = true;
//	if (abs(els_mc_id().at(iel)) != 11) {
//	  cout<<"pt " << el_p4.pt()<<". TRUTHMATCHED with a PDG "<<genps_id().at(idx)<<" status "<<genps_status().at(idx)<<" with mother "<<genps_id_mother().at(idx)<<" with pt "<<genps_p4().at(idx).pt();
//	  cout<<". The lepton has mc_id "<< els_mc_id().at(iel)<<" mc_mother_id "<<els_mc_motherid().at(iel)<<" mcdr "<< els_mcdr().at(iel)<< " mc pt "<< els_mc_p4().at(iel).pt()<<" mc3_id "<<els_mc3_id().at(iel)<<" mc3_mother_id "<<els_mc3_motherid().at(iel) <<endl; 
//	}
	      if (genps_status().at(idx) == 1) status1E = genps_p4().at(idx).E();
	      else status3E = genps_p4().at(idx).E();
	      //cout<<"Matched raw "<< els_eSCRaw().at(iel)<<" with status1 "<<status1pt<<" and status3 "<<status3pt<<endl;
	    }
	  }
	  if (truthmatched) {nPromptEle++; nPromptEleBin[effBin]++; truthE = (status3E > 1.) ? status3E : status1E;}	  
	  else if (abs(els_mc_id().at(iel)) == 11) { nNonPromptEle++; nNPEleBin[effBin]++; nonpromptmatched = true; 
	    //cout<<"NONPROMPT: mc_id "<< els_mc_id().at(iel)<<" mc_mother_id "<<els_mc_motherid().at(iel)<<" mcdr "<< els_mcdr().at(iel)<< " mc3_id "<<els_mc3_id().at(iel)<<" mc3_mother_id "<<els_mc3_motherid().at(iel)<<". bestmatch is "<<bestmatch<<" with mother "<<bestmatchMother <<endl;     
    } 
	  else if (bestmatch > 0.2) { nFakeEle++; nFakeEleBin[effBin]++; fakematched = true; }
	  // Classify NONPROMPT as FAKE
	  if (nonpromptmatched == true) fakematched = true;
	}
	if (isData) fakematched = true;



        // Fill the electron structure
        dilepStudyLooper::electron eleStruct; 
	fillElectronStructure(iel, eleStruct,  useMap, doPFCandLoop, useDeltaBetaWeights, DeltaBetaSimple, areaCorrection, true, truthmatched, truthE);
	//if ( fabs(eleStruct.z0corr) < 1 ) continue;        // STUDY DZ TAILS
        // Fill "hall"              
        fillElectronQuantities(hSet1, eleStruct);        


	fillTrueFakeHistos( eta, truthmatched, fakematched, nonpromptmatched, hSet2, hSet2f, hSet2np, hSet2E, hSet2Ef, hSet2Enp, eleStruct);
	if (pt > 10 && pt < 20) fillTrueFakeHistos( eta, truthmatched, fakematched, nonpromptmatched, hSet8, hSet8f, hSet8np, hSet8E, hSet8Ef, hSet8Enp, eleStruct);
	if (pt > 20) 	fillTrueFakeHistos( eta, truthmatched, fakematched, nonpromptmatched, hSet9, hSet9f, hSet9np, hSet9E, hSet9Ef, hSet9Enp, eleStruct);

 
        // Now start with the selections
        // First plot the N-1 distributions for the WP80 selections. 
        // Plot all variables for electrons passing all cuts 2->N

	// After denominator cuts
	if (passElectronCuts(eleStruct, eleDencuts, eleDencutsEE)) {
	  fillTrueFakeHistos( eta, truthmatched, fakematched, nonpromptmatched, hSet7, hSet7f, hSet7np, hSet7E, hSet7Ef, hSet7Enp, eleStruct);
	}
	
	// N-1
	fillTrueFakeHistosN1( eta, truthmatched, fakematched, nonpromptmatched, hSet3, hSet3f, hSet3np, hSet3E, hSet3Ef, hSet3Enp, eleStruct, eleWP80cuts, eleWP80cutsEE);

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

      fillTrueFakeHistosN1( eta, truthmatched, fakematched, nonpromptmatched, hSet4, hSet4f, hSet4np, hSet4E, hSet4Ef, hSet4Enp,   eleStruct, eleLcuts, eleLcutsEE);
      fillTrueFakeHistosN1( eta, truthmatched, fakematched, nonpromptmatched, hSet5, hSet5f, hSet5np, hSet5E, hSet5Ef, hSet5Enp,   eleStruct, eleLcuts, eleLcutsEE);
      fillTrueFakeHistosN1( eta, truthmatched, fakematched, nonpromptmatched, hSet6, hSet6f, hSet6np, hSet6E, hSet6Ef, hSet6Enp,   eleStruct, eleLcuts, eleLcutsEE);

      // Fill Plots After Cuts (not N-1), and counters


      // Plots after ID and ISO cuts //stdIDstdISO
      // Full Standard
      if (passElectronCuts(eleStruct, eleVcuts, eleVcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountVstdIDstdISO);
      if (passElectronCuts(eleStruct, eleLcuts, eleLcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountLstdIDstdISO);
      if (passElectronCuts(eleStruct, eleMcuts, eleMcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountMstdIDstdISO);
      if (passElectronCuts(eleStruct, eleTcuts, eleTcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountTstdIDstdISO);

      // New ID and std ISO //newIDstdISO
      eleVcuts.dEtaIn =   0.012;                 // These are the preliminary CSA14 ID cuts   //eleVcuts.dEtaIn =   0.007; // These are the standard ID cuts     
      eleVcutsEE.dEtaIn = (is25ns ? 0.015 : 0.022); 					      //eleVcutsEE.dEtaIn = 0.01;					       
      eleLcuts.dEtaIn =   0.012; 				                  	      //eleLcuts.dEtaIn =   0.007;				       
      eleLcutsEE.dEtaIn = (is25ns ? 0.014 : 0.021);					      //eleLcutsEE.dEtaIn = 0.009;				       
      eleMcuts.dEtaIn =   0.009;					                      //eleMcuts.dEtaIn =   0.004;				       
      eleMcutsEE.dEtaIn = (is25ns ? 0.012 : 0.019);					      //eleMcutsEE.dEtaIn = 0.007;				       
      eleTcuts.dEtaIn =   0.009;		                 			      //eleTcuts.dEtaIn =   0.004;				       
      eleTcutsEE.dEtaIn = (is25ns ? 0.010 : 0.017);					      //eleTcutsEE.dEtaIn = 0.005;       
      eleLcutsEE.hOverE = 0.12;
      eleMcutsEE.hOverE = 0.12;
      eleTcutsEE.hOverE = 0.12;
      eleVcuts.sieie = 0.013;
      eleLcuts.sieie = 0.013;
      eleVcutsEE.sieie = 0.033;
      eleLcutsEE.sieie = 0.033;
      eleMcutsEE.sieie = 0.031;
      eleTcutsEE.sieie = 0.031;
                                
      if (passElectronCuts(eleStruct, eleVcuts, eleVcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountVnewIDstdISO);
      if (passElectronCuts(eleStruct, eleLcuts, eleLcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountLnewIDstdISO);
      if (passElectronCuts(eleStruct, eleMcuts, eleMcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountMnewIDstdISO);
      if (passElectronCuts(eleStruct, eleTcuts, eleTcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountTnewIDstdISO);

      // New ID and new ISO //newIDnewISO
      eleVcuts.iso_cor =   (is25ns ? 0.23 : 0.20); // These are the preliminary CSA14 ISO cuts   //eleVcuts.iso_cor =   0.15; // These are the standard ISO cuts   
      eleVcutsEE.iso_cor = (is25ns ? 0.25 : 0.20);						 //eleVcutsEE.iso_cor = 0.15; // should be 0.10 for <20 GeV        
      eleLcuts.iso_cor =   (is25ns ? 0.23 : 0.20);						 //eleLcuts.iso_cor =   0.15;				      
      eleLcutsEE.iso_cor = (is25ns ? 0.25 : 0.20);						 //eleLcutsEE.iso_cor = 0.15; // should be 0.10 for <20 GeV	      
      eleMcuts.iso_cor =   (is25ns ? 0.23 : 0.20);						 //eleMcuts.iso_cor =   0.15;				      
      eleMcutsEE.iso_cor = (is25ns ? 0.21 : 0.18);						 //eleMcutsEE.iso_cor = 0.15; // should be 0.10 for <20 GeV	      
      eleTcuts.iso_cor =   (is25ns ? 0.18 : 0.15);						 //eleTcuts.iso_cor =   0.1;					      
      eleTcutsEE.iso_cor = (is25ns ? 0.16 : 0.13);						 //eleTcutsEE.iso_cor = 0.1; // should be 0.10 for <20 GeV         
      if (passElectronCuts(eleStruct, eleVcuts, eleVcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountVnewIDnewISO);
      if (passElectronCuts(eleStruct, eleLcuts, eleLcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountLnewIDnewISO);
      if (passElectronCuts(eleStruct, eleMcuts, eleMcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountMnewIDnewISO);
      if (passElectronCuts(eleStruct, eleTcuts, eleTcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountTnewIDnewISO);

      // Turn ISO off: New ID only //newIDnoISO  
      eleVcuts.iso_cor   = 999;
      eleVcutsEE.iso_cor = 999;
      eleLcuts.iso_cor   = 999;
      eleLcutsEE.iso_cor = 999;
      eleMcuts.iso_cor   = 999;
      eleMcutsEE.iso_cor = 999;
      eleTcuts.iso_cor   = 999;
      eleTcutsEE.iso_cor = 999;
      if (passElectronCuts(eleStruct, eleVcuts, eleVcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountVnewIDnoISO);
      if (passElectronCuts(eleStruct, eleLcuts, eleLcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountLnewIDnoISO);
      if (passElectronCuts(eleStruct, eleMcuts, eleMcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountMnewIDnoISO);
      if (passElectronCuts(eleStruct, eleTcuts, eleTcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountTnewIDnoISO);               



      // Turn ISO off: Standard ID only //stdIDnoISO  
      eleVcuts.dEtaIn =   0.007; // These are the standard ID cuts   
      eleVcutsEE.dEtaIn = 0.01;				       
      eleLcuts.dEtaIn =   0.007;				       
      eleLcutsEE.dEtaIn = 0.009;				       
      eleMcuts.dEtaIn =   0.004;				       
      eleMcutsEE.dEtaIn = 0.007;				       
      eleTcuts.dEtaIn =   0.004;				       
      eleTcutsEE.dEtaIn = 0.005;   
      eleLcutsEE.hOverE = 0.10;
      eleMcutsEE.hOverE = 0.10;
      eleTcutsEE.hOverE = 0.10;
      eleVcuts.sieie = 0.01;
      eleLcuts.sieie = 0.01;
      eleVcutsEE.sieie = 0.03;
      eleLcutsEE.sieie = 0.03;
      eleMcutsEE.sieie = 0.03;
      eleTcutsEE.sieie = 0.03;               
      if (passElectronCuts(eleStruct, eleVcuts, eleVcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountVstdIDnoISO);
      if (passElectronCuts(eleStruct, eleLcuts, eleLcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountLstdIDnoISO);
      if (passElectronCuts(eleStruct, eleMcuts, eleMcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountMstdIDnoISO);
      if (passElectronCuts(eleStruct, eleTcuts, eleTcutsEE)) fillCounters(effBin, truthmatched, fakematched, nonpromptmatched, CountTstdIDnoISO);

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
cout << "Reco truthmatch: " << nPromptEle << endl;
cout << "Reco non-prompt: " << nNonPromptEle << endl;
cout << "Reco fake: " << nFakeEle << endl;

 outFile->cd(); // Make sure histograms get written out

 // std::map<std::string, TGraph*> hGraph;

//stdIDstdISO
//newIDstdISO
//newIDnewISO
//stdIDnoISO	
//newIDnoISO  

  for (unsigned int b = 0; b < 6; b++) {
    int tot = nPromptEleBin[b];
    int Vet = CountVstdIDstdISO[b][0];
    int Los = CountLstdIDstdISO[b][0];
    int Med = CountMstdIDstdISO[b][0];
    int Tig = CountTstdIDstdISO[b][0];
    //cout<<"eff "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double eff[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    tot = nFakeEleBin[b] + nNPEleBin[b];
    Vet = CountVstdIDstdISO[b][1] + CountVstdIDstdISO[b][2];
    Los = CountLstdIDstdISO[b][1] + CountLstdIDstdISO[b][2];
    Med = CountMstdIDstdISO[b][1] + CountMstdIDstdISO[b][2];
    Tig = CountTstdIDstdISO[b][1] + CountTstdIDstdISO[b][2];
    //cout<<"fake "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double fake[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    TGraph * g_0  = new TGraph(4, fake, eff);
    TString gname = "stdIDstdISO_" + TString::Itoa(b, 10) + "_ROC";
    g_0->SetNameTitle(gname, gname);
    g_0->Write();
  }
  for (unsigned int b = 0; b < 6; b++) {
    int tot = nPromptEleBin[b];
    int Vet = CountVnewIDstdISO[b][0];
    int Los = CountLnewIDstdISO[b][0];
    int Med = CountMnewIDstdISO[b][0];
    int Tig = CountTnewIDstdISO[b][0];
    //cout<<"eff "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double eff[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    tot = nFakeEleBin[b] + nNPEleBin[b];
    Vet = CountVnewIDstdISO[b][1] + CountVnewIDstdISO[b][2];
    Los = CountLnewIDstdISO[b][1] + CountLnewIDstdISO[b][2];
    Med = CountMnewIDstdISO[b][1] + CountMnewIDstdISO[b][2];
    Tig = CountTnewIDstdISO[b][1] + CountTnewIDstdISO[b][2];
    //cout<<"fake "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double fake[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    TGraph * g_0  = new TGraph(4, fake, eff);
    TString gname = "newIDstdISO_" + TString::Itoa(b, 10) + "_ROC";
    g_0->SetNameTitle(gname, gname);
    g_0->Write();
  }  for (unsigned int b = 0; b < 6; b++) {
    int tot = nPromptEleBin[b];
    int Vet = CountVnewIDnewISO[b][0];
    int Los = CountLnewIDnewISO[b][0];
    int Med = CountMnewIDnewISO[b][0];
    int Tig = CountTnewIDnewISO[b][0];
    //cout<<"eff "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double eff[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    tot = nFakeEleBin[b] + nNPEleBin[b];
    Vet = CountVnewIDnewISO[b][1] + CountVnewIDnewISO[b][2];
    Los = CountLnewIDnewISO[b][1] + CountLnewIDnewISO[b][2];
    Med = CountMnewIDnewISO[b][1] + CountMnewIDnewISO[b][2];
    Tig = CountTnewIDnewISO[b][1] + CountTnewIDnewISO[b][2];
    //cout<<"fake "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double fake[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    TGraph * g_0  = new TGraph(4, fake, eff);
    TString gname = "newIDnewISO_" + TString::Itoa(b, 10) + "_ROC";
    g_0->SetNameTitle(gname, gname);
    g_0->Write();
  }  for (unsigned int b = 0; b < 6; b++) {
    int tot = nPromptEleBin[b];
    int Vet = CountVstdIDnoISO[b][0];
    int Los = CountLstdIDnoISO[b][0];
    int Med = CountMstdIDnoISO[b][0];
    int Tig = CountTstdIDnoISO[b][0];
    //cout<<"eff "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double eff[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    tot = nFakeEleBin[b] + nNPEleBin[b];
    Vet = CountVstdIDnoISO[b][1] + CountVstdIDnoISO[b][2];
    Los = CountLstdIDnoISO[b][1] + CountLstdIDnoISO[b][2];
    Med = CountMstdIDnoISO[b][1] + CountMstdIDnoISO[b][2];
    Tig = CountTstdIDnoISO[b][1] + CountTstdIDnoISO[b][2];
    //cout<<"fake "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double fake[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    TGraph * g_0  = new TGraph(4, fake, eff);
    TString gname = "stdIDnoISO_" + TString::Itoa(b, 10) + "_ROC";
    g_0->SetNameTitle(gname, gname);
    g_0->Write();
  }  for (unsigned int b = 0; b < 6; b++) {
    int tot = nPromptEleBin[b];
    int Vet = CountVnewIDnoISO[b][0];
    int Los = CountLnewIDnoISO[b][0];
    int Med = CountMnewIDnoISO[b][0];
    int Tig = CountTnewIDnoISO[b][0];
    //cout<<"eff "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double eff[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    tot = nFakeEleBin[b] + nNPEleBin[b];
    Vet = CountVnewIDnoISO[b][1] + CountVnewIDnoISO[b][2];
    Los = CountLnewIDnoISO[b][1] + CountLnewIDnoISO[b][2];
    Med = CountMnewIDnoISO[b][1] + CountMnewIDnoISO[b][2];
    Tig = CountTnewIDnoISO[b][1] + CountTnewIDnoISO[b][2];
    //cout<<"fake "<<b<<" "<<tot<<" "<<Vet<<" "<<Los<<" "<<Med<<" "<<Tig<<endl;
    double fake[4] = { 1.0* Vet/tot,  1.0* Los/tot,  1.0* Med/tot,  1.0* Tig/tot};
    TGraph * g_0  = new TGraph(4, fake, eff);
    TString gname = "newIDnoISO_" + TString::Itoa(b, 10) + "_ROC";
    g_0->SetNameTitle(gname, gname);
    g_0->Write();
  }


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

  bookElectronHistos(hSet8, "h_true1020");bookElectronHistos(hSet8E, "h_true1020EE");
  bookElectronHistos(hSet8f, "h_fake1020");  bookElectronHistos(hSet8Ef, "h_fake1020EE");   bookElectronHistos(hSet8np, "h_np1020");  bookElectronHistos(hSet8Enp, "h_np1020EE"); 
  bookElectronHistos(hSet9, "h_true20");bookElectronHistos(hSet9E, "h_true20EE");
  bookElectronHistos(hSet9f, "h_fake20");  bookElectronHistos(hSet9Ef, "h_fake20EE");   bookElectronHistos(hSet9np, "h_np20");  bookElectronHistos(hSet9Enp, "h_np20EE"); 


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
  TH1F * h_el_dEtaIn = new TH1F(Form("%s_el_dEtaIn",prefix.Data()),";electron dEtaIn",80,-0.04,0.04);
  TH1F * h_el_dPhiIn = new TH1F(Form("%s_el_dPhiIn",prefix.Data()),";electron dPhiIn",100,-0.2,0.2);
  TH1F * h_el_hOverE = new TH1F(Form("%s_el_hOverE",prefix.Data()),";electron hOverE",100,0.,0.5);
  TH1F * h_el_hOverEBC = new TH1F(Form("%s_el_hOverEBC",prefix.Data()),";electron hOverEBC",100,0.,0.5);
  TH1F * h_el_d0corr = new TH1F(Form("%s_el_d0corr",prefix.Data()),";electron d0corr",100,-0.2,0.2);
  TH1F * h_el_z0corr = new TH1F(Form("%s_el_z0corr",prefix.Data()),";electron z0corr",100,-1.,1.);
  TH1F * h_el_ooemoop = new TH1F(Form("%s_el_ooemoop",prefix.Data()),";electron ooemoop",100,-0.5,0.5);
  TH1F * h_el_iso_cor = new TH1F(Form("%s_el_iso_cor",prefix.Data()),";electron iso_cor",55,0.,1.1);
  TH1F * h_el_iso_corsize = new TH1F(Form("%s_el_iso_corzise",prefix.Data()),";electron iso_corsize",55,0.,1.1);
  TH1F * h_el_pfiso_chPU = new TH1F(Form("%s_el_pfiso_chPU",prefix.Data()),";electron pfiso_chPU",55,0.,1.1);
  TH1F * h_el_pfiso_ch = new TH1F(Form("%s_el_pfiso_ch",prefix.Data()),";electron pfiso_ch",55,0.,1.1);
  TH1F * h_el_pfiso_em = new TH1F(Form("%s_el_pfiso_em",prefix.Data()),";electron pfiso_em",55,0.,1.1);
  TH1F * h_el_pfiso_nh = new TH1F(Form("%s_el_pfiso_nh",prefix.Data()),";electron pfiso_nh",55,0.,1.1);
  TH1F * h_el_detiso_ch = new TH1F(Form("%s_el_detiso_ch",prefix.Data()),";electron detiso_ch",55,0.,1.1);
  TH1F * h_el_detiso_em = new TH1F(Form("%s_el_detiso_em",prefix.Data()),";electron detiso_em",55,0.,1.1);
  TH1F * h_el_detiso_nh = new TH1F(Form("%s_el_detiso_nh",prefix.Data()),";electron detiso_nh",55,0.,1.1);
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
  TH1F * h_el_eRawOverTrue  = new TH1F(Form("%s_el_eRawOverTrue" , prefix.Data()), ";electron eRawOverTrue ", 100, 0, 2);
  TH1F * h_el_eRawMinusTrue  = new TH1F(Form("%s_el_eRawMinusTrue" , prefix.Data()), ";electron eRawMinusTrue ", 100, -50, 50);
  TH1F * h_el_eRawMtrueOtrue  = new TH1F(Form("%s_el_eRawMtrueOtrue" , prefix.Data()), ";electron eRawMtrueOtrue ", 100, -1, 1);
  TH1F * h_el_pfiso_ch_pt = new TH1F(Form("%s_el_pfiso_ch_pt", prefix.Data()), "; electron ch_pt", 100, 0, 10);
  TH1F * h_el_pfiso_em_pt = new TH1F(Form("%s_el_pfiso_em_pt", prefix.Data()), "; electron em_pt", 100, 0, 10);
  TH1F * h_el_pfiso_nh_pt = new TH1F(Form("%s_el_pfiso_nh_pt", prefix.Data()), "; electron nh_pt", 100, 0, 10);
  TH1F * h_el_pfiso_chPU_pt = new TH1F(Form("%s_el_pfiso_ch_ptPU", prefix.Data()), "; electron ch_ptPU", 100, 0, 10);
  TH1F * h_el_pfiso_el_pt = new TH1F(Form("%s_el_pfiso_el_pt", prefix.Data()), "; electron el_pt", 100, 0, 10);
  TH1F * h_el_pfiso_mu_pt = new TH1F(Form("%s_el_pfiso_mu_pt", prefix.Data()), "; electron mu_pt", 100, 0, 10);
  TH1F * h_el_pfiso_ch_N  = new TH1F(Form("%s_el_pfiso_ch_N" , prefix.Data()), "; electron ch_N", 20, 0, 20);
  TH1F * h_el_pfiso_em_N  = new TH1F(Form("%s_el_pfiso_em_N" , prefix.Data()), "; electron em_N", 20, 0, 20);
  TH1F * h_el_pfiso_nh_N  = new TH1F(Form("%s_el_pfiso_nh_N" , prefix.Data()), "; electron nh_N", 20, 0, 20);
  TH1F * h_el_pfiso_chPU_N  = new TH1F(Form("%s_el_pfiso_ch_NPU" , prefix.Data()), "; electron ch_NPU", 20, 0, 20);
  TH1F * h_el_pfiso_el_N  = new TH1F(Form("%s_el_pfiso_el_N" , prefix.Data()), "; electron el_N", 20, 0, 20);
  TH1F * h_el_pfiso_mu_N  = new TH1F(Form("%s_el_pfiso_mu_N" , prefix.Data()), "; electron mu_N", 20, 0, 20);
  hSet.insert(std::pair<std::string, TH1F*> ("pt",h_el_pt));
  hSet.insert(std::pair<std::string, TH1F*> ("sieie",h_el_sieie));
  hSet.insert(std::pair<std::string, TH1F*> ("dEtaIn",h_el_dEtaIn));   
  hSet.insert(std::pair<std::string, TH1F*> ("dPhiIn",h_el_dPhiIn));   
  hSet.insert(std::pair<std::string, TH1F*> ("hOverE",h_el_hOverE));
  hSet.insert(std::pair<std::string, TH1F*> ("hOverEBC",h_el_hOverEBC));
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
  hSet.insert(std::pair<std::string, TH1F*> ("eRawOverTrue"  , h_el_eRawOverTrue   ));
  hSet.insert(std::pair<std::string, TH1F*> ("eRawMinusTrue"  , h_el_eRawMinusTrue   ));
  hSet.insert(std::pair<std::string, TH1F*> ("eRawMtrueOtrue"  , h_el_eRawMtrueOtrue   ));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_ch_pt",h_el_pfiso_ch_pt));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_em_pt",h_el_pfiso_em_pt));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_nh_pt",h_el_pfiso_nh_pt));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_chPU_pt",h_el_pfiso_chPU_pt));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_el_pt",h_el_pfiso_el_pt));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_mu_pt",h_el_pfiso_mu_pt));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_ch_N",h_el_pfiso_ch_N));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_em_N",h_el_pfiso_em_N));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_nh_N",h_el_pfiso_nh_N));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_chPU_N",h_el_pfiso_chPU_N));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_el_N",h_el_pfiso_el_N));
  hSet.insert(std::pair<std::string, TH1F*> ("pfiso_mu_N",h_el_pfiso_mu_N));

  return;
}


void dilepStudyLooper::fillElectronQuantities(std::map<std::string, TH1F*> & hSet, electron e ){

  fillUnderOverFlow( hSet["pt"]           , e.pt, 1);
  fillUnderOverFlow( hSet["sieie"]           , e.sieie, 1);
  fillUnderOverFlow( hSet["dEtaIn"]          , e.dEtaIn, 1);
  fillUnderOverFlow( hSet["dPhiIn"]          , e.dPhiIn, 1);
  fillUnderOverFlow( hSet["hOverE"]          , e.hOverE, 1);
  fillUnderOverFlow( hSet["hOverEBC"]        , e.hOverEBC, 1);
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
  fillUnderOverFlow( hSet["eRawOverTrue"] , e.eRawOverTrue  , 1);	  
  fillUnderOverFlow( hSet["eRawMinusTrue"] , e.eRawMinusTrue  , 1);	  
  fillUnderOverFlow( hSet["eRawMtrueOtrue"] , e.eRawMtrueOtrue  , 1);	  



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
  if ( (p & (all & ~(1ll<<i))) == (all & ~(1ll<<i)) ) {hSet["hOverE"]->Fill(e.hOverE, 1);  fillUnderOverFlow( hSet["hOverEBC"]     , e.hOverEBC , 1);}			   i++;
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
    fillUnderOverFlow( hSet["eRawOverTrue"] , e.eRawOverTrue  , 1);	  
    fillUnderOverFlow( hSet["eRawMinusTrue"] , e.eRawMinusTrue  , 1);	  
    fillUnderOverFlow( hSet["eRawMtrueOtrue"] , e.eRawMtrueOtrue  , 1);	  
    fillUnderOverFlow( hSet["pfiso_chPU"]   , e.pfiso_chPU , 1);
  
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
  bool passVtx = cut.vtxFitConversion ? e.vtxFitConversion : true;
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

void  dilepStudyLooper::electronPFiso2(float &pfiso_ch, float &pfiso_em, float &pfiso_nh, float &pfiso_chPU, const float R, const unsigned int iel, const int ivtx, bool useMap, int useDeltaBetaWeights, bool fillPlots, hMAP & hSet) {
  // Mostly taken from electronIsoValuePF2012 in electronSelections.cc

  pfiso_ch = 0;
  pfiso_em = 0;
  pfiso_nh = 0;
  pfiso_chPU = 0;
  int n_ch = 0, n_em = 0, n_chPU = 0, n_nh = 0, n_el = 0, n_mu = 0;
  for (unsigned int ipf = 0; ipf < cms2.pfcands_p4().size(); ++ipf) {
    
    // skip electrons and muons
    const int particleId = abs(cms2.pfcands_particleId()[ipf]);

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

    LorentzVector & p4 = cms2.pfcands_p4()[ipf];
    // Fill some plots
    if (fillPlots) {
      if (particleId == 211 && !isPU)  { fillUnderOverFlow( hSet["pfiso_ch_pt"]   , p4.pt(), 1); n_ch++; }
      if (particleId == 211 &&  isPU)  { fillUnderOverFlow( hSet["pfiso_chPU_pt"] , p4.pt(), 1); n_chPU++; }
      if (particleId == 22)            { fillUnderOverFlow( hSet["pfiso_em_pt"]   , p4.pt(), 1); n_em++; }
      if (particleId == 130)           { fillUnderOverFlow( hSet["pfiso_nh_pt"]   , p4.pt(), 1); n_nh++; }
      if (particleId == 11)            { fillUnderOverFlow( hSet["pfiso_el_pt"]   , p4.pt(), 1); n_el++; }
      if (particleId == 13)            { fillUnderOverFlow( hSet["pfiso_mu_pt"]   , p4.pt(), 1); n_mu++; }
    }
    if (particleId == 11)    continue;
    if (particleId == 13)    continue;
    // Use Map (only REL 7)
    bool skipcand = false;
    if (useMap) {
      std::vector<int> pfCandMatchedToEle = els_PFCand_idx().at(iel);
      for ( unsigned int i = 0; i < pfCandMatchedToEle.size(); i++ ) {
	if (pfCandMatchedToEle[i] == (int) ipf) skipcand = true;
      }
    }
    if (skipcand) continue;

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
  if (fillPlots) {
    fillUnderOverFlow( hSet["pfiso_ch_N"]   ,  n_ch  , 1);
    fillUnderOverFlow( hSet["pfiso_chPU_N"] ,  n_chPU, 1);
    fillUnderOverFlow( hSet["pfiso_em_N"]   ,  n_em  , 1);
    fillUnderOverFlow( hSet["pfiso_nh_N"]   ,  n_nh  , 1);
    fillUnderOverFlow( hSet["pfiso_el_N"]   ,  n_el  , 1);
    fillUnderOverFlow( hSet["pfiso_mu_N"]   ,  n_mu  , 1);
  }


  return;
}

void dilepStudyLooper::fillElectronStructure( const unsigned int iel, electron & eleStruct, bool useMap, bool doPFCandLoop, int useDeltaBetaWeights, bool DeltaBetaSimple, bool areaCorrection, bool fillPFiso, bool truth, float truthE) {
	float ch = 0;
	float em = 0;
	float nh = 0;
	float chPU = 0;
	float pt = els_p4().at(iel).pt();
	float eta = els_etaSC().at(iel);
	eleStruct.pt = els_p4().at(iel).pt() ;
        eleStruct.eta = els_etaSC().at(iel);
        eleStruct.sieie = m_miniAOD ? els_sigmaIEtaIEta_full5x5().at(iel) : els_sigmaIEtaIEta().at(iel);
        eleStruct.dEtaIn = els_etaSC().at(iel) > 0 ? els_dEtaIn().at(iel) : -1.*els_dEtaIn().at(iel);
        eleStruct.dPhiIn = els_dPhiIn().at(iel);
        eleStruct.hOverE = els_hOverE().at(iel);
        eleStruct.hOverEBC = els_hOverEBC().at(iel);
        eleStruct.d0corr = els_dxyPV().at(iel); //els_d0corr().at(iel); // for miniAOD move to els_dxyPV
        eleStruct.z0corr = els_dzPV().at(iel); //dzPV(els_vertex_p4()[iel], els_trk_p4()[iel], vtxs_position()[0]);//els_z0corr().at(iel); // for miniAOD move to els_dzPV
        eleStruct.ooemoop = (1.0/els_ecalEnergy().at(iel)) - (els_eOverPIn().at(iel)/els_ecalEnergy().at(iel)) ;
	eleStruct.iso_uncor = (fillPFiso && !doPFCandLoop) ? electronPFiso(iel, false) : 0;
	if (fillPFiso) {
	  if (!doPFCandLoop) eleStruct.iso_cor = electronPFiso(iel, areaCorrection);// false = NO PILEUP AREA CORRECTION
	  else {
	    
	    if (m_miniAOD && useDeltaBetaWeights==0 && useMap==false && is25ns ) { // if everything is standard, just take the standard values
	      ch = els_pfChargedHadronIso().at(iel);
	      nh = els_pfNeutralHadronIso().at(iel);
	      em = els_pfPhotonIso().at(iel);
	      chPU = els_pfPUIso().at(iel);
	    } 
	    else  {
	      bool  fillPlots = true;
	      if ( truth ) {
		if (fabs(eta) <= 1.479) electronPFiso2(ch, em, nh, chPU, 0.3, iel, firstGoodVertex(), useMap, useDeltaBetaWeights, fillPlots, hSet2);
		else electronPFiso2(ch, em, nh, chPU, 0.3, iel, firstGoodVertex(), useMap, useDeltaBetaWeights, fillPlots, hSet2E);
	      }
	      else {
		if (fabs(eta) <= 1.479) electronPFiso2(ch, em, nh, chPU, 0.3, iel, firstGoodVertex(), useMap, useDeltaBetaWeights, fillPlots, hSet2f);
		else electronPFiso2(ch, em, nh, chPU, 0.3, iel, firstGoodVertex(), useMap, useDeltaBetaWeights, fillPlots, hSet2Ef);
	      }

	    //cout<<"ch "<<ch<<" "<<els_pfChargedHadronIso().at(iel)<<endl;
	    //cout<<"nh "<<nh<<" "<<els_pfNeutralHadronIso().at(iel)<<endl;
	    //cout<<"em "<<em<<" "<<els_pfPhotonIso().at(iel)<<endl;
	    //cout<<"PU "<<chPU<<" "<<els_pfPUIso().at(iel)<<endl;
	    //cout<<endl;
	    }
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
	eleStruct.eRawOverTrue =  (truthE > 0) ? els_eSCRaw().at(iel) / truthE : -999;
	eleStruct.eRawMinusTrue =  (truthE > 0) ? els_eSCRaw().at(iel) - truthE : -999;
	eleStruct.eRawMtrueOtrue =  (truthE > 0) ? (els_eSCRaw().at(iel) - truthE) / truthE : -999;

}

void dilepStudyLooper::fillTrueFakeHistos( float eta, bool truth, bool fake, bool nonprompt, hMAP & hSet, hMAP & hSetf, hMAP & hSetnp, hMAP & hSetE, hMAP & hSetEf, hMAP & hSetEnp, electron & eleStruct){

  if (fabs(eta) <= 1.479)   {
    if (truth) fillElectronQuantities(hSet, eleStruct);
    else if (fake) fillElectronQuantities(hSetf, eleStruct);
    else if (nonprompt) fillElectronQuantities(hSetnp, eleStruct);
  }
  else {
    if (truth) fillElectronQuantities(hSetE, eleStruct);
    else if (fake) fillElectronQuantities(hSetEf, eleStruct);
    else if (nonprompt) fillElectronQuantities(hSetEnp, eleStruct);		 
  }
}

void dilepStudyLooper::fillTrueFakeHistosN1( float eta, bool truth, bool fake, bool nonprompt, hMAP & hSet, hMAP & hSetf, hMAP & hSetnp, hMAP & hSetE, hMAP & hSetEf, hMAP & hSetEnp, electron & eleStruct,  electron cut, electron cutEE){

  if (fabs(eta) <= 1.479)   {
    if (truth) fillElectronQuantitiesN1(hSet, eleStruct, cut);
    else if (fake) fillElectronQuantitiesN1(hSetf, eleStruct, cut);
    else if (nonprompt) fillElectronQuantitiesN1(hSetnp, eleStruct, cut);
  }
  else {
    if (truth) fillElectronQuantitiesN1(hSetE, eleStruct, cutEE);
    else if (fake) fillElectronQuantitiesN1(hSetEf, eleStruct, cutEE);
    else if (nonprompt) fillElectronQuantitiesN1(hSetEnp, eleStruct, cutEE);		 
  }

}
void dilepStudyLooper::fillCounters( int effBin,  bool truth, bool fake, bool nonprompt, int CountV[][3]) { 

  if (truth) CountV[effBin][0]++;
  else if (nonprompt)  CountV[effBin][1]++;
  else if (fake) CountV[effBin][2]++;
  
}
