#ifndef dilepStudyLooper_h
#define dilepStudyLooper_h

#include <vector>
#include <list> 
#include <string>
#include <map>
#include <set>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

#include "TChain.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
typedef std::map<std::string, TH1F*> hMAP;

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

struct isovals {
  float chiso00;
  float chiso04;
  float chiso07;
  float chiso10;
};


class dilepStudyLooper
{
public: 
  dilepStudyLooper();
  ~dilepStudyLooper() {}

  int  ScanChain(TChain *chain, const TString& prefix = "", int sign = 0, int isocortype = 0 );
  void BookHistos (const TString& prefix);
  void InitBaby();
  float electronPFiso(const unsigned int index, const bool cor = false);
  void  electronPFiso2(float &pfiso_ch, float &pfiso_em, float &pfiso_nh,  float &pfiso_chPU, const float R, const unsigned int iel, const int ivtx, bool useMap, int useDeltaBetaWeights, bool fillPlots, hMAP & hSet);

  float muonPFiso(const unsigned int imu, const bool cor = false);
  float dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 );
  float dRbetweenVectors2(LorentzVector vec1, LorentzVector vec2 );
  float Mt( LorentzVector p4, float met, float met_phi );
  void labelAxis(TH2F* h, int axis, int lep);
  
  //
  // WW Electron Id
  //
  bool   ww_elBase(unsigned int i);
  bool   ww_eld0PV(unsigned int i);
  bool   ww_eldZPV(unsigned int i);
  bool ElectronFOV4(unsigned int i);
  bool ElectronFOIdV4(unsigned int i);
  double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv);
  int primaryVertex();

  // Set globals
  void set_createTree   (bool  b)    { g_createTree   = b; }
  void set_version      (const char* v)    { g_version      = v; }
  void set_json         (const char* v)    { g_json         = v; }        

  // Baby ntuple methods
  void makeOutput (const TString& prefix);
  void closeOutput ();

  // trigger-object matching
  int findTriggerIndex(const TString& trigName);
  TString triggerName(const TString& triggerPattern);
  bool objectPassTrigger(const LorentzVector &obj, const TString& trigname, int type, float drmax = 0.1 );

  // utils
  float getdphi( float phi1 , float phi2 );
  isovals muonChIsoValuePF2012 (const unsigned int imu, const float R = 0.3, const int ivtx = 0);

  struct electron {
    float pt;
    float eta;
    float sieie;
    float dEtaIn;
    float dPhiIn;
    float hOverE;
    float hOverEBC;
    float d0corr;
    float z0corr;
    float ooemoop;
    float iso_cor;
    float iso_uncor;
    float pfiso_chPU;
    float pfiso_ch;
    float pfiso_em;
    float pfiso_nh;
    float detiso_ch;
    float detiso_em;
    float detiso_nh;
    float valid_pixelHits;
    float lost_pixelhits;
    bool vtxFitConversion;   
    float chi2n;
    float fbrem    ;
    float dEtaOut  ;
    float etaWidth ;
    float phiWidth ;
    float e1x5e5x5 ;
    float r9       ;
    float sieieSC  ;
    float eoverpIn ;
    float eoverpOut;
    float psOverRaw;
    int   seed;
    int   ncluster;
    float eRawOverTrue;
    float eRawMinusTrue;
    float eRawMtrueOtrue;
  };

  struct effRejCounter {
    int t; //truth
    int np; //non-prompt
    int f; //fake
  };

  void fillElectronStructure( const unsigned int iel, electron & eleStruct, bool useMap, bool doPFCandLoop, int useDeltaBetaWeights, bool DeltaBetaSimple, bool areaCorrection, bool fillPFiso, bool truth, float truthpt = 0);
  void fillElectronQuantities(std::map<std::string, TH1F*> & hSet, electron e );
  void bookElectronHistos(std::map<std::string, TH1F*> & hSet, TString prefix );
  void fillElectronQuantitiesN1(std::map<std::string, TH1F*> & hSet, electron e, electron cut);
  void fillElectronCutsResult(electron e, electron cut, ULong64_t & pass);
  void fillUnderOverFlow(TH1F *h1, float value, float weight);
  bool passElectronCuts(electron eleStruct, electron cut, electron cutEE);
  TH1F* MakeEfficiencyPlot(TH1F* num_hist, TH1F* den_hist, const std::string& name, const std::string& title);

  void fillTrueFakeHistos( float eta, bool truth, bool fake, bool nonprompt, hMAP & hSet, hMAP & hSetf, hMAP & hSetnp, hMAP & hSetE, hMAP & hSetEf, hMAP & hSetEnp, electron & eleStruct);
  void fillTrueFakeHistosN1( float eta, bool truth, bool fake, bool nonprompt, hMAP & hSet, hMAP & hSetf, hMAP & hSetnp, hMAP & hSetE, hMAP & hSetEf, hMAP & hSetEnp, electron & eleStruct,  electron cut, electron cutEE);
  //  void fillCounters( int effBin,  bool truth, bool fake, bool nonprompt, effRejCounter * CountV);
  void fillCounters( int effBin,  bool truth, bool fake, bool nonprompt, int CountV[][3]); 



private:

  // Globals
  bool  g_createTree;
  const char* g_version;
  const char* g_json;      
  bool initialized;
  bool isdata_;
  TFile* outFile;

  float ptthresh_high;
  float ptthresh_low;

  // histograms
  TH1F* h_nvtx;

  TH1F* h_el_pt;
  TH1F* h_el_eta;
  TH1F* h_el_truthmatch_DR;
  TH1F* h_el_lead_pt;
  TH1F* h_el_subl_pt;
  TH1F* h_el_lead_eta;
  TH1F* h_el_subl_eta;
  TH1F* h_ee_mll;
  TH1F* h_ee_denom_mll;
  TH1F* h_ee_mediso_mll;

  TH1F* h_em_el_hlt_pt;
  TH1F* h_em_mu_hlt_pt;
  TH1F* h_em_el_hlt_eta;
  TH1F* h_em_mu_hlt_eta;
  TH1F* h_em_hlt_dr;
  TH1F* h_em_hlt_mll;

  TH1F* h_em_el_hlt_noreco_pt;
  TH1F* h_em_mu_hlt_noreco_pt;
  TH1F* h_em_el_hlt_noreco_eta;
  TH1F* h_em_mu_hlt_noreco_eta;
  TH1F* h_em_hlt_noreco_dr;
  TH1F* h_em_hlt_noreco_mll;

  TH2F* h_ee_events;
  
//  TH1F* h_el_sieie;
//  TH1F* h_el_dEtaIn;
//  TH1F* h_el_dPhiIn;

    std::map<std::string, TH1F*> hSet1;
    std::map<std::string, TH1F*> hSet2;
    std::map<std::string, TH1F*> hSet3;
    std::map<std::string, TH1F*> hSet4;
    std::map<std::string, TH1F*> hSet5;
    std::map<std::string, TH1F*> hSet6;
  std::map<std::string, TH1F*> hSet2f;
  std::map<std::string, TH1F*> hSet3f;
  std::map<std::string, TH1F*> hSet4f;
  std::map<std::string, TH1F*> hSet5f;
  std::map<std::string, TH1F*> hSet6f;
  std::map<std::string, TH1F*> hSet2np;
  std::map<std::string, TH1F*> hSet3np;
  std::map<std::string, TH1F*> hSet4np;
  std::map<std::string, TH1F*> hSet5np;
  std::map<std::string, TH1F*> hSet6np;
  std::map<std::string, TH1F*> hSet2E;
  std::map<std::string, TH1F*> hSet3E;
  std::map<std::string, TH1F*> hSet4E;
  std::map<std::string, TH1F*> hSet5E;
  std::map<std::string, TH1F*> hSet6E;
  std::map<std::string, TH1F*> hSet2Ef;
  std::map<std::string, TH1F*> hSet3Ef;
  std::map<std::string, TH1F*> hSet4Ef;
  std::map<std::string, TH1F*> hSet5Ef;
  std::map<std::string, TH1F*> hSet6Ef;
  std::map<std::string, TH1F*> hSet2Enp;
  std::map<std::string, TH1F*> hSet3Enp;
  std::map<std::string, TH1F*> hSet4Enp;
  std::map<std::string, TH1F*> hSet5Enp;
  std::map<std::string, TH1F*> hSet6Enp;

    std::map<std::string, TH1F*> hSet7;
    std::map<std::string, TH1F*> hSet7f;
    std::map<std::string, TH1F*> hSet7np;
    std::map<std::string, TH1F*> hSet7E;
    std::map<std::string, TH1F*> hSet7Ef;
    std::map<std::string, TH1F*> hSet7Enp;
    std::map<std::string, TH1F*> hSet8;
    std::map<std::string, TH1F*> hSet8f;
    std::map<std::string, TH1F*> hSet8np;
    std::map<std::string, TH1F*> hSet8E;
    std::map<std::string, TH1F*> hSet8Ef;
    std::map<std::string, TH1F*> hSet8Enp;
    std::map<std::string, TH1F*> hSet9;
    std::map<std::string, TH1F*> hSet9f;
    std::map<std::string, TH1F*> hSet9np;
    std::map<std::string, TH1F*> hSet9E;
    std::map<std::string, TH1F*> hSet9Ef;
    std::map<std::string, TH1F*> hSet9Enp;


std::map<std::string, TH1F*> hSetCut1;  
std::map<std::string, TH1F*> hSetCut2;
std::map<std::string, TH1F*> hSetCut3;
std::map<std::string, TH1F*> hSetCut4;
std::map<std::string, TH1F*> hSetCut5;
std::map<std::string, TH1F*> hSetCut6;
std::map<std::string, TH1F*> hSetCut7;
std::map<std::string, TH1F*> hSetCut8;
std::map<std::string, TH1F*> hSetCut9;

std::map<TString, TH1F*> hSetEff;

 enum EgammaElectronType {
   ISECALENERGYCORRECTED,  // if false, the electron "ecalEnergy" is just the supercluster energy 
   ISMOMENTUMCORRECTED,    // has E-p combination been applied
   ISECALDRIVEN,
   ISTRACKERDRIVEN
 };

enum EgammaFiduciality {
  ISEB,
  ISEBEEGAP,
  ISEE,
  ISEEGAP,
  ISEBETAGAP,
  ISEBPHIGAP,
  ISEEDEEGAP,
  ISEERINGGAP,
  ISGAP
};


};

#endif
