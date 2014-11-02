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
//const int nptbins = 9;
//const float ptbins[] = {0.00, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 110.0, 150};
const int nptbins = 11;
const float ptbins[] = {0.00, 20.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 110.0, 150};

const int nhtbins = 11;
const float htbins[] = {0.00, 100.0, 175, 250., 325., 400.0, 475.0, 550.0, 625.0, 700.0, 800.0, 900.0};

const int nptbins2 = 7;
const float ptbins2[] = {30.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 150};

const int nhtbins2 = 5;
const float htbins2[] = {0.00, 150.0, 250., 350., 500.0, 900.0};

const int nQ2bins = 6;
const float Q2bins[] = {0.00, 150.0, 250., 350., 500.0, 700.0, 1500};

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
      bool is_qcd = false;
      if (filename.Contains("MuEnriched")) is_qcd = true;
      bool is_bjet = false;
      if (filename.Contains("BJet")) is_bjet = true;
      
      if (!isData) {
	// First, find the muon from b
	LorentzVector b1p4;
	LorentzVector b2p4;
	LorentzVector muFromBp4;
	std::vector<LorentzVector> qgp4;
	bool foundb1 = false;
	bool foundb2 = false;
	bool foundMuFromB = false;
	int idxb1 = -1;
	int idxb2 = -1;
	int idxMuFromB = -1;
	float status1rho = 0;
	float status1rho25 = 0;
	for(unsigned int idx = 0; idx < genps_id().size(); idx++) {	 
	  int pid = abs(genps_id().at(idx));
	  int status = genps_status().at(idx);
	  // total pt at status 1
	  if (status==1 && fabs(genps_p4().at(idx).eta()) < 1.2 && !(pid==12 || pid==14|| pid==16)) status1rho += genps_p4().at(idx).pt();
	  if (status==1 && fabs(genps_p4().at(idx).eta()) < 2.5 && !(pid==12 || pid==14|| pid==16)) status1rho25 += genps_p4().at(idx).pt();

	  if (genps_p4().at(idx).pt() < 10) continue;
	  if (pid==5 && status==3) {
	    if (!foundb1) {b1p4 = genps_p4().at(idx); foundb1 = true; idxb1 = idx;}
	    else if (!foundb2) {b2p4 = genps_p4().at(idx); foundb2 = true; idxb2 = idx;}	    
	  }
	  if (pid==13 && status==1 && !foundMuFromB) {
	    if ( idIsBeauty( genps_id_mother().at(idx)) && genps_p4().at(idx).pt() > 20 && fabs(genps_p4().at(idx).eta()) < 0.9 ) { // ETA RESTRICTED
	      muFromBp4 = genps_p4().at(idx);
	      foundMuFromB = true;
	      idxMuFromB = idx;
	    }
	  }
	  if ( status==3 ) { // get the other status3 quarks/gluons
	    if ((pid<5 || pid==21) && fabs(genps_p4().at(idx).eta()) < 2 ) {
	      qgp4.push_back(genps_p4().at(idx));
	    }
	  }
	} // loop over Gen
	
	if (!foundMuFromB) continue;
	// which b is the right one?
	LorentzVector bP4;
	int bIdx = -1;
	float dRb1 = dRbetweenVectors2(b1p4, muFromBp4 ) ;
	float dRb2 = dRbetweenVectors2(b2p4, muFromBp4 ) ;
	//cout<<"found muon from b, has DR2 "<<dRb1<<" and "<<dRb2<<endl;
	if (dRb1 <= dRb2 && dRb1 < 0.04) {bP4 = b1p4; bIdx  = idxb1;}
	if (dRb1 >= dRb2 && dRb2 < 0.04) {bP4 = b2p4; bIdx  = idxb2;}
	//cout<<"assigned to "<<bIdx<<endl;
	if (bIdx < 0) continue;
	if (bP4.pt() < 30 || bP4.pt() > 200) continue;

	// Check DR between mu and others
	bool bIsIso = true;
	int otherBidx = (bIdx == idxb1) ? idxb2 : idxb1;
	if (otherBidx>=0) if (dRbetweenVectors2(muFromBp4, genps_p4().at(otherBidx)) < 0.64 ) bIsIso = false;
	for (unsigned int i = 0; i < qgp4.size(); i++) {
	  if (dRbetweenVectors2(muFromBp4, qgp4.at(i)) < 0.64 ) bIsIso = false;
	}
	// Use genJets for DR with muon
	// Can't do this for ttbar_pythia6_8tev sample
	bool bIsIsoGenJet = true;
	bool bIsIsoGenJet40 = true;
	if (filename.Contains("ttbar_pythia6_8tev")) bIsIsoGenJet = false;
	else {
	  for(unsigned int idx = 0; idx < genjets_p4().size(); idx++) {
	    float pt = genjets_p4().at(idx).pt();
	    if (pt < 20) continue;
	    float dR = dRbetweenVectors2(genjets_p4().at(idx), muFromBp4 );
	    if (dR < 0.64) { bIsIsoGenJet = false; }
	    if (dR < 0.64 && pt > 40) { bIsIsoGenJet40 = false; continue;}

	  }
	}

	
	// Next calculate muon isolation
	// Remove all status1 within 0.4 of status3 b/q/g. Remove status3 leptons. Only the UE is left
	float iso = 0;
	float iso05 = 0;
	float rhoFromHard = 0;
	for(unsigned int idx = 0; idx < genps_id().size(); idx++) {	  
	  int pid = abs(genps_id().at(idx));      
	  if (pid==12 || pid==14|| pid==16) continue;
	  float ipt =  genps_p4().at(idx).pt();
	  LorentzVector ip4 =  genps_p4().at(idx);
	  int status = genps_status().at(idx);
	  // Energy from hard scattering
	  if (fabs(genps_p4().at(idx).eta()) < 1.2) {
	    if (status==3 && (pid==11 || pid==13 || pid==15)) rhoFromHard += ipt;
	    if (status == 1) {
        if (idxb1!=-1) if (dRbetweenVectors2(b1p4, ip4) < 0.16 ) rhoFromHard += ipt;
        if (idxb2!=-1) if (dRbetweenVectors2(b2p4, ip4) < 0.16 ) rhoFromHard += ipt;
        for (unsigned int iqg = 0; iqg < qgp4.size(); iqg++) {
          if (dRbetweenVectors2(qgp4.at(iqg), ip4) < 0.16 ) rhoFromHard += ipt;
        }
      }
	  }	 
	  // calculate muon isolation
	  if (pid!=13 && status==1) {
      float dR2 = dRbetweenVectors2(muFromBp4,  ip4 ) ;
      if (dR2 > 0.09) continue;
      iso += ipt;
      if (genps_p4().at(idx).pt()>0.5)  iso05 += ipt;
    }
	}

	float bPt = bP4.pt();
	double evt_weight = 1;
	if (is_qcd) {
          double p0 = -0.31816;
          double p1 = 0.0411647;
          double p2 = -0.00164925;
          double p3 = 2.20514e-05;
          double wForMC3PT = p0 + p1 * bPt + p2 * bPt * bPt + p3 * bPt * bPt * bPt;
          evt_weight *= wForMC3PT;
        }

	float rhoUE = status1rho - rhoFromHard;

	double evt_weight2 = evt_weight;
	if (is_qcd) {
	  if (rhoUE < 10) evt_weight2 *= 1.14;
	  else if (rhoUE < 300) {
	    double p0 = 0.735926;
	    double p1 = 0.00315761;
	    double wForUE = p0 + p1 * rhoUE;
	    evt_weight2 *= wForUE;
	  }
	  else evt_weight2 *= 1;
        }
	
	// Here we want b-jets to look like QCD. We need to cut some tails. Keep 80-450.
	double evt_weightHT25_B = 1;
	if (status1rho25>80 && status1rho25<450 && is_bjet) {
	  double p0 = 213.621;
	  double p1 = -3.00715;
	  double p2 = 0.0169485;
	  double p3 = -4.3368e-05;
	  double p4 = 4.21988e-08;
	  evt_weightHT25_B *= p0 + p1 * status1rho25 + p2 * status1rho25 * status1rho25 + p3 * status1rho25 * status1rho25 * status1rho25 + p4 * status1rho25*status1rho25*status1rho25*status1rho25;
	}
	else if (status1rho25<80 || status1rho25>450) evt_weightHT25_B = 0.; // Got rid of below 80 and above 450 for everyone

	// Then we want qcd (and bjets) to look like ttbar
	double evt_weightHT25 = evt_weightHT25_B;
	if (status1rho25>150 && status1rho25<450 && (is_bjet || is_qcd)) {
	  double constant = 2.01193e+02;
	  double MPV = 7.60896e+02;
	  double Sigma = 1.71578e+02;
	  evt_weightHT25 *= TMath::Landau(status1rho25, MPV, Sigma)*constant;
	}
	else if (status1rho25<150 || status1rho25>450) evt_weightHT25 = 0.; 

	// Now we do the same for the normal status1rho
	// Here we want b-jets to look like QCD. We need to cut some tails. Keep 80-450.
	double evt_weightHT_B = 1;
	if (status1rho>50 && status1rho<350 && is_bjet) {
	  double p0 = 114.453;
	  double p1 = -1.92313;
	  double p2 = 0.0134909;
	  double p3 = -4.42116e-05;
	  double p4 = 5.70836e-08;
	  evt_weightHT_B *= p0 + p1 * status1rho + p2 * status1rho * status1rho + p3 * status1rho * status1rho * status1rho + p4 * status1rho*status1rho*status1rho*status1rho;
	}
	else if (status1rho<50 || status1rho>350) evt_weightHT_B = 0.; 

	// Then we want qcd (and bjets) to look like ttbar
	double evt_weightHT = evt_weightHT_B;
	if (status1rho>100 && status1rho<350 && (is_bjet || is_qcd)) {
	  double constant = 5.76853e+02;
	  double MPV = 8.02710e+02;
	  double Sigma = 2.05853e+02;
	  evt_weightHT *= TMath::Landau(status1rho, MPV, Sigma)*constant;
	}
	else if (status1rho<100 || status1rho>350) evt_weightHT = 0.; 
	
	hSetEff["h_status1rho"]->Fill(status1rho, evt_weight);
	hSetEff["h_status1rhoW1"]->Fill(status1rho, 1);
	hSetEff["h_status1rho25W1"]->Fill(status1rho25, 1);
	hSetEff["h_status1rho25WHT"]->Fill(status1rho25, evt_weightHT25);
	hSetEff["h_status1rhoWHT"]->Fill(status1rho, evt_weightHT);
	hSetEff["h_status1rho25WHTB"]->Fill(status1rho25, evt_weightHT25_B);
	hSetEff["h_status1rhoWHTB"]->Fill(status1rho, evt_weightHT_B);
	hSetEff["h_rhoFromHard"]->Fill(rhoFromHard, evt_weight);
	hSetEff["h_rhoUE"]->Fill(rhoUE, evt_weight);
	hSetEff["h_rhoUEW2"]->Fill(rhoUE, evt_weight2);


	float muPt = muFromBp4.pt();
	float reliso = 1.0*iso/muPt;
	float reliso05 = 1.0*iso05/muPt;
	hSetEff["h_bquark_pt"]->Fill(bPt, evt_weight);
	hSetEff["h_bquark_pt_w1"]->Fill(bPt, 1);
	hSetEff["h_bquark_eta"]->Fill(bP4.eta(), evt_weight);
	hSetEff["h_bquark_eta_w1"]->Fill(bP4.eta(), 1);
	hSetEff["h_bquark_p"]->Fill(bP4.P(), evt_weight);
	hSetEff["h_bquark_p_w1"]->Fill(bP4.P(), 1);

	//cout<<"Found muon with pt "<<muPt<<" and reliso "<<reliso<<" and bPt "<<bPt<<endl;
	if (bPt > 40 && bPt < 50 ) {
	  hSetEff["h_iso_4050"]->Fill(reliso, evt_weight);
	  //cout<<"filled! "<<reliso<<endl;
	}

	//float scale = (pdfinfo_x1() + pdfinfo_x2())*4000;

	if (bPt > 60 && bPt < 80 )  hSetEff["h_iso_6080"]->Fill(reliso, evt_weight);

	if (reliso < 1.0) hSetEff["h_ptb_den"]->Fill(bPt, evt_weight);
	if (reliso < 0.3) hSetEff["h_ptb_num"]->Fill(bPt, evt_weight);
	
	// HT categories to show that QCD FR depends on HT
	if (status1rho < 150) { 
	  if (reliso < 1.0) hSetEff["h_ptb_den_ht1"]->Fill(bPt, evt_weight);
	  if (reliso < 0.3) hSetEff["h_ptb_num_ht1"]->Fill(bPt, evt_weight);
	}
	else if (status1rho > 150 && status1rho < 220) {
	  if (reliso < 1.0) hSetEff["h_ptb_den_ht2"]->Fill(bPt, evt_weight);
	  if (reliso < 0.3) hSetEff["h_ptb_num_ht2"]->Fill(bPt, evt_weight);
	}
	else if (status1rho > 220) {
	  if (reliso < 1.0) hSetEff["h_ptb_den_ht3"]->Fill(bPt, evt_weight);
	  if (reliso < 0.3) hSetEff["h_ptb_num_ht3"]->Fill(bPt, evt_weight);
	}
	if (reliso < 1.0) hSetEff2D["h_ptb_htS1_den"]->Fill(bPt, status1rho, 1);
	if (reliso < 0.3) hSetEff2D["h_ptb_htS1_num"]->Fill(bPt, status1rho, 1);
	if (reliso < 1.0) hSetEff2D["h_ptb_htS3_den"]->Fill(bPt, rhoFromHard, 1);
	if (reliso < 0.3) hSetEff2D["h_ptb_htS3_num"]->Fill(bPt, rhoFromHard, 1);
//	if (reliso < 1.0) hSetEff2D["h_ptb_x1x2_den"]->Fill(bPt, scale, 1);
//	if (reliso < 0.3) hSetEff2D["h_ptb_x1x2_num"]->Fill(bPt, scale, 1);
	if (reliso < 1.0) hSetEff["h_htS1_den"]->Fill(status1rho, evt_weight);
	if (reliso < 0.3) hSetEff["h_htS1_num"]->Fill(status1rho, evt_weight);
	if (reliso < 1.0) hSetEff["h_htS3_den"]->Fill(rhoFromHard, evt_weight);
	if (reliso < 0.3) hSetEff["h_htS3_num"]->Fill(rhoFromHard, evt_weight);
	if (reliso < 1.0) hSetEff["h_ptb_denW1"]->Fill(bPt, 1);
	if (reliso < 0.3) hSetEff["h_ptb_numW1"]->Fill(bPt, 1);
	if (reliso < 1.0) hSetEff["h_ptb_denWHT"]->Fill(bPt, evt_weightHT);
	if (reliso < 0.3) hSetEff["h_ptb_numWHT"]->Fill(bPt, evt_weightHT);
	if (reliso < 1.0) hSetEff["h_ptb_denWHT25"]->Fill(bPt, evt_weightHT25);
	if (reliso < 0.3) hSetEff["h_ptb_numWHT25"]->Fill(bPt, evt_weightHT25);
	if (reliso < 1.0) hSetEff["h_ptb_denWHTB"]->Fill(bPt, evt_weightHT_B);
	if (reliso < 0.3) hSetEff["h_ptb_numWHTB"]->Fill(bPt, evt_weightHT_B);
	if (reliso < 1.0) hSetEff["h_ptb_denWHTB25"]->Fill(bPt, evt_weightHT25_B);
	if (reliso < 0.3) hSetEff["h_ptb_numWHTB25"]->Fill(bPt, evt_weightHT25_B);
	if (bIsIso) { 
	  if (reliso < 1.0) hSetEff["h_ptbiso_den"]->Fill(bPt, evt_weight);
	  if (reliso < 0.3) hSetEff["h_ptbiso_num"]->Fill(bPt, evt_weight);
	}
	if (bIsIsoGenJet) { 
	  if (reliso < 1.0) hSetEff["h_ptbisojet_den"]->Fill(bPt, evt_weight);
	  if (reliso < 0.3) hSetEff["h_ptbisojet_num"]->Fill(bPt, evt_weight);
	}
	if (bIsIsoGenJet40) { 
	  if (reliso < 1.0) hSetEff["h_ptbisojet40_den"]->Fill(bPt, evt_weight);
	  if (reliso < 0.3) hSetEff["h_ptbisojet40_num"]->Fill(bPt, evt_weight);
	}
	if (reliso05 < 1.0) hSetEff["h_ptb05_den"]->Fill(bPt, evt_weight);
	if (reliso05 < 0.3) hSetEff["h_ptb05_num"]->Fill(bPt, evt_weight);

      } // isData
        
    ++nEventsPass;
    
    //      outTree->Fill();
    
  } // entries
  
  delete f;
} // currentFile


cout << endl;
cout << "Sample: " << prefix << endl;
cout << endl;
cout << "Processed events: " << nEventsTotal << endl;

 outFile->cd(); // Make sure histograms get written out

 hSetEff["h_ptb_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_num"] , hSetEff["h_ptb_den"] , "h_ptb_eff", "Eff ;pTb");
 hSetEff["h_ptb_eff_ht1"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_num_ht1"] , hSetEff["h_ptb_den_ht1"] , "h_ptb_eff", "Eff ;pTb");
 hSetEff["h_ptb_eff_ht2"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_num_ht2"] , hSetEff["h_ptb_den_ht2"] , "h_ptb_eff", "Eff ;pTb");
 hSetEff["h_ptb_eff_ht3"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_num_ht3"] , hSetEff["h_ptb_den_ht3"] , "h_ptb_eff", "Eff ;pTb");
 hSetEff2D["h_ptb_htS1_eff"] = (TH2F*) MakeEfficiencyPlot2D(hSetEff2D["h_ptb_htS1_num"] , hSetEff2D["h_ptb_htS1_den"] , "h_ptb_htS1_eff", "Eff ;pTb ;HT S1");
 hSetEff2D["h_ptb_htS3_eff"] = (TH2F*) MakeEfficiencyPlot2D(hSetEff2D["h_ptb_htS3_num"] , hSetEff2D["h_ptb_htS3_den"] , "h_ptb_htS3_eff", "Eff ;pTb ;HT S3");
 // hSetEff2D["h_ptb_x1x2_eff"] = (TH2F*) MakeEfficiencyPlot2D(hSetEff2D["h_ptb_x1x2_num"] , hSetEff2D["h_ptb_x1x2_den"] , "h_ptb_x1x2_eff", "Eff ;pTb ;x1+x2");
 hSetEff["h_htS1_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_htS1_num"] , hSetEff["h_htS1_den"] , "h_htS1_eff", "Eff ;SumPtStatus1");
 hSetEff["h_htS3_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_htS3_num"] , hSetEff["h_htS3_den"] , "h_htS3_eff", "Eff ;SumPtStatus3");
 hSetEff["h_ptb_effW1"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_numW1"] , hSetEff["h_ptb_denW1"] , "h_ptb_effW1", "Eff ;pTb");
 hSetEff["h_ptb_effWHT"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_numWHT"] , hSetEff["h_ptb_denWHT"] , "h_ptb_effWHT", "Eff ;pTb");
 hSetEff["h_ptb_effWHT25"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_numWHT25"] , hSetEff["h_ptb_denWHT25"] , "h_ptb_effWHT25", "Eff ;pTb");
 hSetEff["h_ptb_effWHTB"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_numWHTB"] , hSetEff["h_ptb_denWHTB"] , "h_ptb_effWHTB", "Eff ;pTb");
 hSetEff["h_ptb_effWHTB25"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb_numWHTB25"] , hSetEff["h_ptb_denWHTB25"] , "h_ptb_effWHTB25", "Eff ;pTb");
 hSetEff["h_ptbiso_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptbiso_num"] , hSetEff["h_ptbiso_den"] , "h_ptbiso_eff", "Eff ;pTb");
 hSetEff["h_ptbisojet_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptbisojet_num"] , hSetEff["h_ptbisojet_den"] , "h_ptbisojet_eff", "Eff ;pTb");
 hSetEff["h_ptbisojet40_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptbisojet40_num"] , hSetEff["h_ptbisojet40_den"] , "h_ptbisojet40_eff", "Eff ;pTb");
 hSetEff["h_ptb05_eff"] = (TH1F*) MakeEfficiencyPlot(hSetEff["h_ptb05_num"] , hSetEff["h_ptb05_den"] , "h_ptb05_eff", "Eff ;pTb");
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

TH2*  dilepStudyLooper::MakeEfficiencyPlot2D(TH2* num_hist, TH2* den_hist, const std::string& name, const std::string& title)
{
  // check that hists are valid
  if (!num_hist || !den_hist)
    {
      throw runtime_error("rt::MakeEfficiencyPlot2D: one of the Histograms are NULL");
    }
  
  // verify that all histograms have same binning
  if ((den_hist->GetNbinsX() != num_hist->GetNbinsX()) && (den_hist->GetNbinsY() != num_hist->GetNbinsY())) 
    {
      throw runtime_error("rt::MakeEfficiencyPlot2D: Histograms must have same number of bins");
    }
  
  // get the new histogram
  TH2F* temp = dynamic_cast<TH2F*>(num_hist->Clone(name.c_str()));
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
  TH1::SetDefaultSumw2();

  h_nvtx = new TH1F(Form("%s_nvtx",prefix.Data()),";N(vtx)",max_nvtx,0,max_nvtx);

  hSetEff["h_iso_4050"] = new TH1F(Form("%s_iso_4050",prefix.Data()),";GenIso", 20, 0, 1);
  hSetEff["h_iso_6080"] = new TH1F(Form("%s_iso_6080",prefix.Data()),";GenIso", 20, 0, 1);

  //  hSetEff["h_FRvsBpt"] = new TH1F(Form("%s_FRvsBpt",prefix.Data()),";GenIso", 20, 0, 1);


  hSetEff["h_ptb_num"] = new TH1F(Form("%s_ptb_num",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_den"] = new TH1F(Form("%s_ptb_den",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_num_ht1"] = new TH1F(Form("%s_ptb_num_ht1",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_den_ht1"] = new TH1F(Form("%s_ptb_den_ht1",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_num_ht2"] = new TH1F(Form("%s_ptb_num_ht2",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_den_ht2"] = new TH1F(Form("%s_ptb_den_ht2",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_num_ht3"] = new TH1F(Form("%s_ptb_num_ht3",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_den_ht3"] = new TH1F(Form("%s_ptb_den_ht3",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff2D["h_ptb_htS1_num"] = new TH2F(Form("%s_ptb_htS1_num",prefix.Data()),";pt ;HT S1", nptbins2, ptbins2, nhtbins2, htbins2);
  hSetEff2D["h_ptb_htS1_den"] = new TH2F(Form("%s_ptb_htS1_den",prefix.Data()),";pt ;HT S1", nptbins2, ptbins2, nhtbins2, htbins2);
  hSetEff2D["h_ptb_htS3_num"] = new TH2F(Form("%s_ptb_htS3_num",prefix.Data()),";pt ;HT S3", nptbins2, ptbins2, nhtbins2, htbins2);
  hSetEff2D["h_ptb_htS3_den"] = new TH2F(Form("%s_ptb_htS3_den",prefix.Data()),";pt ;HT S3", nptbins2, ptbins2, nhtbins2, htbins2);
//  hSetEff2D["h_ptb_x1x2_num"] = new TH2F(Form("%s_ptb_x1x2_num",prefix.Data()),";pt ;x1+x2", nptbins2, ptbins2, nhtbins2, htbins2);
//  hSetEff2D["h_ptb_x1x2_den"] = new TH2F(Form("%s_ptb_x1x2_den",prefix.Data()),";pt ;x1+x2", nptbins2, ptbins2, nhtbins2, htbins2);
  hSetEff["h_htS1_num"] = new TH1F(Form("%s_htS1_num",prefix.Data()),";htS1", nhtbins, htbins);
  hSetEff["h_htS1_den"] = new TH1F(Form("%s_htS1_den",prefix.Data()),";htS1", nhtbins, htbins);
  hSetEff["h_htS3_num"] = new TH1F(Form("%s_htS3_num",prefix.Data()),";htS3", nhtbins, htbins);
  hSetEff["h_htS3_den"] = new TH1F(Form("%s_htS3_den",prefix.Data()),";htS3", nhtbins, htbins);
  hSetEff["h_ptb_numW1"] = new TH1F(Form("%s_ptb_numW1",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_denW1"] = new TH1F(Form("%s_ptb_denW1",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_numWHT"] = new TH1F(Form("%s_ptb_numWHT",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_denWHT"] = new TH1F(Form("%s_ptb_denWHT",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_numWHT25"] = new TH1F(Form("%s_ptb_numWHT25",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_denWHT25"] = new TH1F(Form("%s_ptb_denWHT25",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_numWHTB"] = new TH1F(Form("%s_ptb_numWHTB",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_denWHTB"] = new TH1F(Form("%s_ptb_denWHTB",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_numWHTB25"] = new TH1F(Form("%s_ptb_numWHTB25",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb_denWHTB25"] = new TH1F(Form("%s_ptb_denWHTB25",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptbiso_num"] = new TH1F(Form("%s_ptbiso_num",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptbiso_den"] = new TH1F(Form("%s_ptbiso_den",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptbisojet_num"] = new TH1F(Form("%s_ptbisojet_num",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptbisojet_den"] = new TH1F(Form("%s_ptbisojet_den",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptbisojet40_num"] = new TH1F(Form("%s_ptbisojet40_num",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptbisojet40_den"] = new TH1F(Form("%s_ptbisojet40_den",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb05_num"] = new TH1F(Form("%s_ptb05_num",prefix.Data()),";pt", nptbins, ptbins);
  hSetEff["h_ptb05_den"] = new TH1F(Form("%s_ptb05_den",prefix.Data()),";pt", nptbins, ptbins);


  hSetEff["h_bquark_pt"] = new TH1F(Form("%s_bquark_pt",prefix.Data()),";pt", 60, 0, 300);
  hSetEff["h_bquark_pt_w1"] = new TH1F(Form("%s_bquark_ptW1",prefix.Data()),";pt", 60, 0, 300);
  hSetEff["h_bquark_eta"] = new TH1F(Form("%s_bquark_eta",prefix.Data()),";eta", 30, -1.5, 1.5);
  hSetEff["h_bquark_eta_w1"] = new TH1F(Form("%s_bquark_etaW1",prefix.Data()),";eta", 30, -1.5, 1.5);
  hSetEff["h_bquark_p"] = new TH1F(Form("%s_bquark_p",prefix.Data()),";p", 60, 0, 300);
  hSetEff["h_bquark_p_w1"] = new TH1F(Form("%s_bquark_pW1",prefix.Data()),";p", 60, 0, 300);

  hSetEff["h_status1rho"] = new TH1F(Form("%s_status1rho",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_status1rhoW1"] = new TH1F(Form("%s_status1rhoW1",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_status1rho25W1"] = new TH1F(Form("%s_status1rho25W1",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_status1rho25WHT"] = new TH1F(Form("%s_status1rho25WHT",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_status1rhoWHT"] = new TH1F(Form("%s_status1rhoWHT",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_status1rho25WHTB"] = new TH1F(Form("%s_status1rho25WHTB",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_status1rhoWHTB"] = new TH1F(Form("%s_status1rhoWHTB",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_rhoFromHard"] = new TH1F(Form("%s_rhoFromHard",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_rhoUE"] = new TH1F(Form("%s_rhoUE",prefix.Data()),";SumPt", 100, 0, 1000);
  hSetEff["h_rhoUEW2"] = new TH1F(Form("%s_rhoUEW2",prefix.Data()),";SumPt", 100, 0, 1000);


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


 bool dilepStudyLooper::idIsBeauty(int id) {
  id = abs(id);
  if (
      id == 5       ||
      id == 511     ||
      id == 521     ||
      id == 10511   ||
      id == 10521   ||
      id == 513     ||
      id == 523     ||
      id == 10513   ||
      id == 10523   ||
      id == 20513   ||
      id == 20523   ||
      id == 515     ||
      id == 525     ||
      id == 531     ||
      id == 10531   ||
      id == 533     ||
      id == 10533   ||
      id == 20533   ||
      id == 535     ||
      id == 541     ||
      id == 10541   ||
      id == 543     ||
      id == 10543   ||
      id == 20543   ||
      id == 545     ||
      id == 551     ||
      id == 10551   ||
      id == 100551  ||
      id == 110551  ||
      id == 200551  ||
      id == 210551  ||
      id == 553     ||
      id == 10553   ||
      id == 20553   ||
      id == 30553   ||
      id == 100553  ||
      id == 110553  ||
      id == 120553  ||
      id == 130553  ||
      id == 200553  ||
      id == 210553  ||
      id == 220553  ||
      id == 300553  ||
      id == 9000553 ||
      id == 9010553 ||
      id == 555     ||
      id == 10555   ||
      id == 20555   ||
      id == 100555  ||
      id == 110555  ||
      id == 120555  ||
      id == 200555  ||
      id == 557     ||
      id == 100557  ||
      id == 5122    || 
      id == 5112    ||
      id == 5212    ||
      id == 5222    ||
      id == 5114    ||
      id == 5214    ||
      id == 5224    ||
      id == 5132    ||
      id == 5232    ||
      id == 5312    ||
      id == 5322    ||
      id == 5314    ||
      id == 5324    ||
      id == 5332    ||
      id == 5334    ||
      id == 5142    ||
      id == 5242    ||
      id == 5412    ||
      id == 5422    ||
      id == 5414    ||
      id == 5424    ||
      id == 5342    ||
      id == 5432    ||
      id == 5434    ||
      id == 5442    ||
      id == 5444    ||
      id == 5512    ||
      id == 5522    ||
      id == 5514    ||
      id == 5524    ||
      id == 5532    ||
      id == 5534    ||
      id == 5542    ||
      id == 5544    ||
      id == 5554 
      ) {
    return true;
  }
  else return false;
}

