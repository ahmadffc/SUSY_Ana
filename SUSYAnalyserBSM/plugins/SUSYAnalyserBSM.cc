// -*- C++ -*-
//
// Package:    SUSY_Ana/SUSYAnalyserBSM
// Class:      SUSYAnalyserBSM
//
/**\class SUSYAnalyserBSM SUSYAnalyserBSM.cc SUSY_Ana/SUSYAnalyserBSM/plugins/SUSYAnalyserBSM.cc
 Description: [one line class summary]
  Implementation:
  [Notes on implementation]
*/

// Original Author:  Muhammad Ahmad
//         Created:  Tue, 16 Feb 2016 05:06:00 GMT
//
//

        


#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"


#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// ECAL spike cleaning - Swiss cross
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/Common/interface/ValueMap.h"





#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
//
// class declaration
//

using namespace edm;
using namespace std;

class SUSYAnalyserBSM : public edm::EDAnalyzer {
   public:
      explicit SUSYAnalyserBSM(const edm::ParameterSet&);
      ~SUSYAnalyserBSM();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
     virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob();

      // ----------member data ---------------------------
      edm::Service<TFileService> fs_;  

      void createHistogram(const std::string& folderName);
      bool check(std::string process, std::string pCheck);
      void init(const edm::TriggerResults &, const edm::TriggerNames & HLTNames);

      std::map<std::string, unsigned int> prescales;
      std::map<std::string, unsigned int> prescale_counter; 
      std::map<std::string, unsigned int> trigger_indices;

      edm::InputTag muoLabel_;
      edm::InputTag jetLabel_;
      edm::InputTag tauLabel_;
      edm::InputTag metLabel_;
      edm::InputTag rechitBLabel_;
      edm::InputTag rechitELabel_;
      edm::InputTag pvSrc_;
      edm::InputTag triggerLabel_;
      edm::InputTag filterLabel_;
      edm::InputTag rhoLabel_;
      

              
      bool isMCBH;
      bool DEBUG_;


      edm::EDGetTokenT<EcalRecHitCollection> ebRecHitsToken_;
      edm::EDGetTokenT<EcalRecHitCollection> eeRecHitsToken_;

      reco::TrackBase::TrackQuality _trackQuality;
 
      std::vector< std::map<std::string,TH1*> > histos_; 
      std::vector<std::string> cutNames_;
     
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_; 
      edm::EDGetToken eleLabelToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;
      
      //Electron Decisions
      edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
      //Photon Decisions
      edm::EDGetToken phoLabelToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > phoLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > phoTightIdMapToken_; 
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      //TTree
      TTree* tree;
      float JetE[25];
      float JetPx[25];
      float JetPy[25];
      float JetPz[25];
      float JetPt[25];
      float JetEt[25];
      float JetEta[25];
      float JetPhi[25];      
      
      float EleE[25];
      float ElePx[25];
      float ElePy[25];
      float ElePz[25];
      float ElePt[25];
      float EleEt[25];
      float EleEta[25];
      float ElePhi[25]; 

      float PhE[25];
      float PhPx[25];
      float PhPy[25];
      float PhPz[25];
      float PhPt[25];
      float PhEt[25];
      float PhEta[25];
      float PhPhi[25];
      
      float MuE[25];
      float MuPx[25];
      float MuPy[25];
      float MuPz[25];
      float MuPt[25];
      float MuEt[25];
      float MuEta[25];
      float MuPhi[25];    
 


     float miniIso[25];
     float Rvalue;
//      float Iso_hist;
      float ST;
      float HT;
      float mBH;
      float Met;
      float MetPx;
      float MetPy;
      float MetPhi;      
      float Sphericity;
      float JetArr[4];      
      float EleArr[4];
      float MuArr[4];      
      float PhArr[4];
      int NPV;
      int NTracks;
      int NJets;
      int NElectrons;
      int NPhotons;      
      int NMuons; 
      bool NoScrap;     
      
      int Multiplicity; 
      bool isLeptonPhoton;
      bool isEleChannel;
      bool isMuChannel;
      bool isPhChannel;

      //
     
      int runno;
      int evtno;
      int lumiblock;
      int isRealData;
      float muon_d0;
      int No_of_Muons=0;  
      int No_of_Muons_2=0;                                                                              
  float LeadingArr[4];    

  //HLT info (Jets and MET)
  bool firedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1;
  bool firedHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1;
  bool firedHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1;
  bool firedHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1;
  bool firedHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1;
  bool firedHLT_DoubleMu8_Mass8_PFHT300_v1;
  bool firedHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v1;
  bool firedHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v1;

  bool firedHLT_TripleMu_12_10_5_v1;
  bool firedHLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1; 
  bool firedHLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1;
  bool firedHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1;

  bool firedHLT_Ele8_CaloIdM_TrackIdM_PFJet30_v1;
  bool firedHLT_Ele12_CaloIdM_TrackIdM_PFJet30_v1;
  bool firedHLT_Ele18_CaloIdM_TrackIdM_PFJet30_v1;
  bool firedHLT_Ele23_CaloIdM_TrackIdM_PFJet30_v1;
  bool firedHLT_Ele33_CaloIdM_TrackIdM_PFJet30_v1;
  bool firedHLT_Mu8_v1;
  bool firedHLT_Mu17_v1;
  bool firedHLT_Mu24_v1;
  bool firedHLT_Mu34_v1;

  //TODO 
  double Reliso_el;
  double Reliso_mu;
  int pass_eleID_medium;      
  
  Float_t rho_;

  };


  // constants, enums and typedefs

  // static data member definitions

  // constructors and destructor
  SUSYAnalyserBSM::SUSYAnalyserBSM(const edm::ParameterSet& iConfig):
  muoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag")),
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag")),
  tauLabel_(iConfig.getUntrackedParameter<edm::InputTag>("tauTag")),
  metLabel_(iConfig.getUntrackedParameter<edm::InputTag>("metTag")),
  pvSrc_(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertex")),
  triggerLabel_(iConfig.getUntrackedParameter<edm::InputTag>("triggerTag")),  
  filterLabel_(iConfig.getUntrackedParameter<edm::InputTag>("filterTag")),  
  rhoLabel_(iConfig.getUntrackedParameter<edm::InputTag>("rho_lable")),
  isMCBH(iConfig.getUntrackedParameter<bool>("MCLabel",false)),
  DEBUG_(iConfig.getUntrackedParameter<bool>("DEBUG",false)),
  
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  phoLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoLooseIdMap"))),
  phoMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoMediumIdMap"))),
  phoTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoTightIdMap"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")))
  {
   //now do what ever initialization is needed
   beamSpotToken_            = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
   eleLabelToken_            = mayConsume<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronTag"));
   vtxMiniAODToken_          = mayConsume<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("verticesMiniAOD"));
   conversionsMiniAODToken_  = mayConsume< reco::ConversionCollection >(iConfig.getParameter<edm::InputTag>("conversionsMiniAOD"));
   phoLabelToken_            = mayConsume<edm::View<reco::Photon> >(iConfig.getParameter<edm::InputTag>("photonTag"));
 }


SUSYAnalyserBSM::~SUSYAnalyserBSM()

{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SUSYAnalyserBSM::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   // reset all variables:
   std::fill(std::begin( JetE   ),  std::end( JetE   ),   0.    );
   std::fill(std::begin( JetPx  ),  std::end( JetPx  ),   0.    );
   std::fill(std::begin( JetPy  ),  std::end( JetPy  ),   0.    );
   std::fill(std::begin( JetPz  ),  std::end( JetPz  ),   0.    );
   std::fill(std::begin( JetPt  ),  std::end( JetPt  ),   0.    );
   std::fill(std::begin( JetEt  ),  std::end( JetEt  ),   0.    );
   std::fill(std::begin( JetEta ),  std::end( JetEta ),  -9999. );
   std::fill(std::begin( JetPhi ),  std::end( JetPhi ),  -9999. );
   std::fill(std::begin( EleE   ),  std::end( EleE   ),   0.    );
   std::fill(std::begin( ElePx  ),  std::end( ElePx  ),   0.    );
   std::fill(std::begin( ElePy  ),  std::end( ElePy  ),   0.    );
   std::fill(std::begin( ElePz  ),  std::end( ElePz  ),   0.    );
   std::fill(std::begin( ElePt  ),  std::end( ElePt  ),   0.    );
   std::fill(std::begin( EleEt  ),  std::end( EleEt  ),   0.    );
   std::fill(std::begin( EleEta ),  std::end( EleEta ),  -9999. );
   std::fill(std::begin( ElePhi ),  std::end( ElePhi ),  -9999. );
   std::fill(std::begin( MuE    ),  std::end( MuE    ),   0.    );
   std::fill(std::begin( MuPx   ),  std::end( MuPx   ),   0.    );
   std::fill(std::begin( MuPy   ),  std::end( MuPy   ),   0.    );
   std::fill(std::begin( MuPz   ),  std::end( MuPz   ),   0.    );
   std::fill(std::begin( MuPt   ),  std::end( MuPt   ),   0.    );
   std::fill(std::begin( MuEt   ),  std::end( MuEt   ),   0.    );
   std::fill(std::begin( MuEta  ),  std::end( MuEta  ),  -9999. );
   std::fill(std::begin( MuPhi  ),  std::end( MuPhi  ),  -9999. );
   std::fill(std::begin( PhE    ),  std::end( PhE    ),   0.    );
   std::fill(std::begin( PhPx   ),  std::end( PhPx   ),   0.    );
   std::fill(std::begin( PhPy   ),  std::end( PhPy   ),   0.    );
   std::fill(std::begin( PhPz   ),  std::end( PhPz   ),   0.    );
   std::fill(std::begin( PhPt   ),  std::end( PhPt   ),   0.    );
   std::fill(std::begin( PhEt   ),  std::end( PhEt   ),   0.    );
   std::fill(std::begin( PhEta  ),  std::end( PhEta  ),  -9999. );
   std::fill(std::begin( PhPhi  ),  std::end( PhPhi  ),  -9999. );
   std::fill(std::begin( miniIso  ),  std::end( miniIso  ),  0. );
   using namespace edm;

   runno = iEvent.id().run();
   evtno  = iEvent.id().event();
   lumiblock = iEvent.luminosityBlock();
   isRealData = iEvent.isRealData();
   
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   // first: get all objects from the event.
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  // All PF Candidate for alternate isolation
  
   // Get the beam spot
   edm::Handle<reco::BeamSpot> theBeamSpot;
   iEvent.getByToken(beamSpotToken_,theBeamSpot);  

   edm::Handle<edm::View<reco::GsfElectron> > electrons;
   iEvent.getByToken(eleLabelToken_, electrons);
  
  
   //Get the conversions collection
   edm::Handle<reco::ConversionCollection> conversions;
   iEvent.getByToken(conversionsMiniAODToken_, conversions);	

   // Get the electron ID data from the event stream.
   //   // Note: this implies that the VID ID modules have been run upstream.
   // If you need more info, check with the EGM group.
   edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
   iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
   
   edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
   iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);

   edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
   iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);

   edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
   iEvent.getByToken(eleTightIdMapToken_ ,tight_id_decisions);

  edm::Handle<edm::ValueMap<bool> > loose_id_decisions_ph;
  iEvent.getByToken(phoLooseIdMapToken_ ,loose_id_decisions_ph);
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions_ph;
  iEvent.getByToken(phoMediumIdMapToken_,medium_id_decisions_ph);
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions_ph;
  iEvent.getByToken(phoTightIdMapToken_ ,tight_id_decisions_ph);

   edm::Handle<double> rhoHandle;
   iEvent.getByLabel(rhoLabel_,rhoHandle);
   
   if(rhoHandle.isValid()) {
   rho_ = *(rhoHandle.product());
   }

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  iEvent.getByLabel(muoLabel_,muonHandle);
  const edm::View<pat::Muon> & muons = *muonHandle;   

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  iEvent.getByLabel(jetLabel_,jetHandle);
  const edm::View<pat::Jet> & jets = *jetHandle;
   
   edm::Handle<edm::View<pat::MET> > metHandle;
   iEvent.getByLabel(metLabel_,metHandle);
   const edm::View<pat::MET> & mets = *metHandle;
   
   edm::Handle<edm::View<reco::Photon> > photons;
   iEvent.getByToken(phoLabelToken_,photons);


   edm::Handle<edm::View<pat::Tau> > tauHandle;
   iEvent.getByLabel(tauLabel_,tauHandle);
   
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerPrescales_, triggerPrescales);
   // Swiss cross - ECAL cleaning
   edm::View<pat::Photon>::const_iterator photon;
   Handle<EcalRecHitCollection> Brechit;//barrel
   Handle<EcalRecHitCollection> Erechit;//endcap
   iEvent.getByLabel(rechitBLabel_,Brechit);
   iEvent.getByLabel(rechitELabel_,Erechit);

   edm::ESHandle<CaloTopology> pTopology;
   edm::ESHandle<CaloTopology> theCaloTopo_;
   iSetup.get<CaloTopologyRecord>().get(theCaloTopo_);
   
   // Loop over all objects and evaluate BH properties: ST, mBH, multiplicity ...
   ST=0.;
   HT=0.;
//   Iso_hist =0.; 
   //Lorentz vectors for jets, egamma, muons, leptons, and all objects
   math::XYZTLorentzVectorF pJet;
   math::XYZTLorentzVectorF pEle;      
   math::XYZTLorentzVectorF pMu;
   math::XYZTLorentzVectorF pPh;
   math::XYZTLorentzVectorF pLep;
   math::XYZTLorentzVectorF pObj;           

   //Beam spot
   math::XYZPoint bs;
       
   // MET
   math::XYZTLorentzVectorF pBH;
   pBH += mets[0].p4();
   Met = mets[0].pt();
   MetPhi = mets[0].phi();
  
   ST += Met;
   
   MetPx = mets[0].px();
   MetPy = mets[0].py();
      
   //Sphericity
   float sumPx2=mets[0].px()*mets[0].px();
   float sumPy2=mets[0].py()*mets[0].py();
   float sumPxPy=mets[0].px()*mets[0].py();   
   
   //Leading objects
   std::vector<double> leadingJets (4,-1);
   std::vector<double> leadingElectrons (4,-1);
   std::vector<double> leadingPhotons (4,-1);
   std::vector<double> leadingMuons (4,-1);
   std::vector<double> leadingObjects (4,-1);
   std::vector<double> leadingLeptons (4,-1);
      
   float leadingElePt = 0;
   float leadingMuPt  = 0;
   float leadingPhPt  = 0;

   int jetcnt         = -1;
   int elecnt         = -1;
   int mucnt          = -1;
   int phcnt          = -1;
   int ngoodmuons     = 0;
   int ngoodphotons   = 0;
   int ngoodelectrons = 0;
   int ngoodjets      = 0;   
   int nPVcount       = 0;
   
   //------ Primary Vertices------- 
   edm::Handle< reco::VertexCollection > PVCollection; 
   iEvent.getByLabel(pvSrc_, PVCollection);
   const reco::Vertex & vertex_ = PVCollection->front();
   if (PVCollection->empty()) return; // skip the event if no PV found

   if (iEvent.getByLabel(pvSrc_, PVCollection )) {
     for (reco::VertexCollection::const_iterator pv = PVCollection->begin(); pv != PVCollection->end(); ++pv) {
       //--- vertex selection
       //std::cout<<" PV "<<pv->x()<<" "<<pv->y()<<" "pv->z()<<std::endl;
       //std::cout<<"PV parameters "<<pv->isFake()<<" "<<fabs(pv->z())<<" "<<pv->position().Rho()<<std::endl;
       //std::cout<<"PV position "<<pv->position()<<endl;
       if (!pv->isFake() && pv->ndof() > 4 && fabs(pv->z()) <= 24. && pv->position().Rho() <= 2.) ++nPVcount;       
     }
   }
   
   // No scraping
   bool noscrap = true;
   
   // HLT results   
   firedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1 = false;
   firedHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1 = false;
   firedHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1 = false;
   firedHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1 = false;
   firedHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1 = false;
   firedHLT_DoubleMu8_Mass8_PFHT300_v1 = false;
   firedHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v1 = false;
   firedHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v1  = false;

   firedHLT_TripleMu_12_10_5_v1  = false;
   firedHLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1   = false; 
   firedHLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1  = false;
   firedHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1  = false;

   firedHLT_Ele8_CaloIdM_TrackIdM_PFJet30_v1  = false;
   firedHLT_Ele12_CaloIdM_TrackIdM_PFJet30_v1  = false;
   firedHLT_Ele18_CaloIdM_TrackIdM_PFJet30_v1  = false;
   firedHLT_Ele23_CaloIdM_TrackIdM_PFJet30_v1  = false;
   firedHLT_Ele33_CaloIdM_TrackIdM_PFJet30_v1  = false;
   firedHLT_Mu8_v1  = false;
   firedHLT_Mu17_v1  = false;
   firedHLT_Mu24_v1  = false;
   firedHLT_Mu34_v1  = false;

  
   TriggerResults tr;
   Handle<TriggerResults> h_trigRes;
   iEvent.getByLabel(triggerLabel_, h_trigRes);
   tr = *h_trigRes;

   // MET filter results   
  
   std::vector<string> triggerList;
   Service<service::TriggerNamesService> tns;
   bool foundNames = tns->getTrigPaths(*h_trigRes, triggerList);
   if (!foundNames) std::cout << "Could not get trigger names!\n";
   if (tr.size()!=triggerList.size()) std::cout << "ERROR: length of names and paths not the same: " 
                                                << triggerList.size() << "," << tr.size() << endl;
   // dump trigger list at first event
   for (unsigned int i=0; i< tr.size(); i++) {
   //  std::cout << "["<<i<<"] = " << triggerList[i]<<setw(40)<<
   //  ": Prescale " << triggerPrescales->getPrescaleForIndex(i) << ": " << (tr[i].accept() ? "Event Passed" : "Event Failed") << endl;
     if ( !tr[i].accept() == 1 ) continue;
     if( triggerList[i] == "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1")  { firedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1  = true; }
     if( triggerList[i] == "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1")  { firedHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1  = true; }
     if( triggerList[i] == "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1")  { firedHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1  = true; }
     if( triggerList[i] == "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1")  { firedHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1  = true; }
     if( triggerList[i] == "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1")  { firedHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1  = true; }
     if( triggerList[i] == "HLT_DoubleMu8_Mass8_PFHT300_v1")  { firedHLT_DoubleMu8_Mass8_PFHT300_v1  = true; }
     if( triggerList[i] == "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v1")  { firedHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v1  = true; }
     if( triggerList[i] == "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v1")  { firedHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v1  = true; }

     if( triggerList[i] == "HLT_TripleMu_12_10_5_v1")  { firedHLT_TripleMu_12_10_5_v1  = true; }
     if( triggerList[i] == "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1")  { firedHLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1  = true; } 
     if( triggerList[i] == "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1")  { firedHLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1  = true; }
     if( triggerList[i] == "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1")  { firedHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1  = true; }

     if( triggerList[i] == "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v1")  { firedHLT_Ele8_CaloIdM_TrackIdM_PFJet30_v1  = true; }
     if( triggerList[i] == "HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v1")  { firedHLT_Ele12_CaloIdM_TrackIdM_PFJet30_v1  = true; }
     if( triggerList[i] == "HLT_Ele18_CaloIdM_TrackIdM_PFJet30_v1")  { firedHLT_Ele18_CaloIdM_TrackIdM_PFJet30_v1  = true; }
     if( triggerList[i] == "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v1")  { firedHLT_Ele23_CaloIdM_TrackIdM_PFJet30_v1  = true; }
     if( triggerList[i] == "HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v1")  { firedHLT_Ele33_CaloIdM_TrackIdM_PFJet30_v1  = true; }
     if( triggerList[i] == "HLT_Mu8_v1")  { firedHLT_Mu8_v1  = true; }
     if( triggerList[i] == "HLT_Mu17_v1")  { firedHLT_Mu17_v1  = true; }
     if( triggerList[i] == "HLT_Mu24_v1")  { firedHLT_Mu24_v1  = true; }
     if( triggerList[i] == "HLT_Mu34_v1")  { firedHLT_Mu34_v1  = true; }


   }
  
//cout<<"ibefore jet loop 1st Jet pt ="<<JetPt[0]<<" ,2nd Jet pt = "<<JetPt[1]<<endl;

  for(edm::View<pat::Jet>::const_iterator jet = jets.begin(); jet!=jets.end(); ++jet){     
//jetcnt++;
    //cout << " pt of jets === " << jet->pt() << endl; 
      // Loose Jet ID (equivalent to previous medium)
    if((
    jet->neutralHadronEnergyFraction()    <  0.99 &&  // 0.90 for tight
    jet->neutralEmEnergyFraction()        <  0.99 &&  // 0.90 for tight
    jet->numberOfDaughters()              >  1
    )                                             && 
    //((
    //abs(jet->eta())                       <= 2.4  && 
    //jet->chargedHadronEnergyFraction()    >  0    && 
    //jet->chargedMultiplicity()            >  0    && 
    //jet->chargedEmEnergyFraction()        <  0.99
    //)                                             || 
    //abs(jet->eta())                       >  2.4
    //)                                             && 
    abs(jet->eta())                       < 2.4  &&
    jet->pt()                             >  10
    ) {
     jetcnt++;  
/*    cout << "   evtno = " << evtno <<
        "|  jetcnt = " << jetcnt <<
        "|  pt = " << jet->pt() <<
        "|  eta = " << jet->eta() <<
        "|  NHEF = " << jet->neutralHadronEnergyFraction() <<
        "|  NEEF =  "<< jet->neutralEmEnergyFraction() <<
        "|  NOD  === " << jet->numberOfDaughters() <<
        "|  "<<  endl;
*/
    pBH += jet->p4();
    ST  += jet->et();
    HT  += jet->pt();
    sumPx2  += jet->px()*jet->px();
    sumPy2  += jet->py()*jet->py();
    sumPxPy += jet->px()*jet->py();      
       
    leadingJets.push_back(jet->pt());
    leadingObjects.push_back(jet->pt());

    pJet += jet->p4();
    pObj += jet->p4(); 
         
    JetE[jetcnt]   = jet->energy();
    JetPx[jetcnt]  = jet->px();
    JetPy[jetcnt]  = jet->py();
    JetPz[jetcnt]  = jet->pz();     
    JetPt[jetcnt]  = jet->pt();
    JetEt[jetcnt]  = jet->et();
    JetEta[jetcnt] = jet->eta();      
    JetPhi[jetcnt] = jet->phi();
    ++ngoodjets;
   }//JetID
  }   
//cout<<"1st Jet pt ="<<JetPt[0]<<" ,2nd Jet pt = "<<JetPt[1]<<endl;
   
   for (size_t i = 0; i < electrons->size(); ++i){
   const auto e = electrons->ptrAt(i);
    //if(!e->gsfTrack()) continue;
    ++elecnt;
  
    // Electron Medium ID
    if(
    e->pt()           		  >  20.    &&
    fabs(e->eta())    		  <  2.5    &&
    (*medium_id_decisions)[e]      == 1      
    ){  
   
    pBH += e->p4();
    ST  += e->et();
    sumPx2 += e->px()*e->px();
    sumPy2 += e->py()*e->py();
    sumPxPy += e->px()*e->py();
    if (e->pt() > leadingElePt) leadingElePt = e->pt();

    leadingElectrons.push_back(e->pt());
    leadingLeptons.push_back(e->pt());
    leadingObjects.push_back(e->pt());                       
    pEle += e->p4();
    pLep += e->p4();      
    pObj += e->p4(); 
    EleE[elecnt]   = e->energy();
    ElePx[elecnt]  = e->px();
    ElePy[elecnt]  = e->py();
    ElePz[elecnt]  = e->pz();     
    ElePt[elecnt]  = e->pt();
    EleEt[elecnt]  = e->et();
    EleEta[elecnt] = e->eta();      
    ElePhi[elecnt] = e->phi();
    //(*loose_id_decisions)[e];
    //(*tight_id_decisions)[e];
    ++ngoodelectrons;
     }//eleID
    }  
   

        
  for (size_t i = 0; i < photons->size(); ++i){
  const auto ph = photons->ptrAt(i);
  ++phcnt;

 if(
  ph->pt()                        >   20       &&
  abs(ph->eta())                  <   2.4      &&
  (*medium_id_decisions_ph)[ph]   ==  1
  ){
  
  //If not a spike, increment # photons
  ++ngoodphotons;

  pBH += ph->p4();
  ST += ph->et();
  sumPx2 += ph->px()*ph->px();
  sumPy2 += ph->py()*ph->py();
  sumPxPy += ph->px()*ph->py();
  if (ph->pt() > leadingPhPt) leadingPhPt = ph->pt();

  leadingPhotons.push_back(ph->pt());
  leadingObjects.push_back(ph->pt());      

  pPh += ph->p4();    
  pObj += ph->p4();      

  PhE[ngoodphotons-1] = ph->energy();
  PhPx[ngoodphotons-1] = ph->px();
  PhPy[ngoodphotons-1] = ph->py();
  PhPz[ngoodphotons-1] = ph->pz();
  PhPt[ngoodphotons-1] = ph->pt();
  PhEt[ngoodphotons-1] = ph->et();
  PhEta[ngoodphotons-1] = ph->eta();
  //(*loose_id_decisions_ph)[ph];
  //(*tight_id_decisions_ph)[ph];      
  PhPhi[ngoodphotons-1] = ph->phi();
  }
  
}


edm::Handle< double > rrho_;
  iEvent.getByLabel("fixedGridRhoFastjetCentralNeutral", rrho_); // Central rho recommended for SUSY
  double rrho = *rrho_;

  rrho =  rrho+1;

//  cout << "value of rho_ === "<< rrho << endl;



  for(edm::View<pat::Muon>::const_iterator mu = muons.begin(); mu!=muons.end(); ++mu){
    ++mucnt;
  
     //double MuEta = mu->eta();
if(!mu->innerTrack().isNull()){


 float EffectiveArea = 0.0;

        if (fabs(mu->eta()) >= 0.0 && fabs(mu->eta()) < 0.8 ) EffectiveArea = 0.0735;
        if (fabs(mu->eta()) >= 0.8 && fabs(mu->eta()) < 1.3 ) EffectiveArea = 0.0619;
        if (fabs(mu->eta()) >= 1.3 && fabs(mu->eta()) < 2.0 ) EffectiveArea = 0.0465;
        if (fabs(mu->eta()) >= 2.0 && fabs(mu->eta()) < 2.2 )   EffectiveArea = 0.0433;
        if (fabs(mu->eta()) >= 2.2 && fabs(mu->eta()) <= 2.5 )   EffectiveArea = 0.0577;



     Rvalue = (10.0 / min(max(mu->pt(),50.0),200.0)) * (10.0 / min(max(mu->pt(),50.0),200.0)) ;

     miniIso[mucnt-1]= (mu->pfIsolationR04().sumChargedHadronPt + max(0., mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt - rrho * Rvalue * EffectiveArea ))/mu->pt();

  float product = EffectiveArea * rrho;

if (evtno == 13222)
{

cout << "Aaaa = "<<mu->pfIsolationR04().sumChargedHadronPt<< ", BBB = " << mu->pfIsolationR04().sumNeutralHadronEt + mu->pfIsolationR04().sumPhotonEt << ", Cccc = " << rrho * Rvalue * EffectiveArea << endl;
cout << "   evtno = " << evtno <<
	"|  Mun["<<mucnt<<"]" <<
	"|  miniIso== " << miniIso[mucnt-1] <<
        "|  eta = " << mu->eta() << 
	"|  EffectiveArea = " << EffectiveArea <<
	"|  rho =  "<< rrho<<
	"|  product  === " << product <<
	"|  "<<  endl; 
}

   if( mu->pt() > 5  &&  
     fabs(mu->eta()) < 2.4 &&
     fabs(mu->innerTrack()->dxy(vertex_.position())) < 0.05 &&
     fabs(mu->innerTrack()->dz(vertex_.position())) < 0.1 ) 
   {
    ++ngoodmuons;

//  float iso = gamma.chargedHadronIso;
//  iso = max(iso - rho*ea, (float)0.);

  // EffectiveArea=EffectiveArea+1.0;

//   for (size_t i=0; i<cutNames_.size(); ++i){
//     std::map<std::string,TH1*>& histo = histos_[i];
//     histo["Iso_hist"]->Fill(miniIso);
//   }

  //cout << "value of rho_ ==222222222222= "<< rrho << endl;

//  cout << "eta = " << mu->eta() << ", EffectiveArea = " << EffectiveArea << ",rho =  "<< rrho<< endl;

//  cout << "  EffectiveArea * rrho === " << product << endl;

//     cout <<" value of Isolation === "<< miniIso <<endl; 

//    cout<< " vtx_dxy  +++== " << fabs(mu->innerTrack()->dxy(vertex_.position())) << " vtx_dz  +++== " << fabs(mu->innerTrack()->dz(vertex_.position())) <<endl;




    pBH += mu->p4();
    ST += mu->et();
    sumPx2 += mu->px()*mu->px();
    sumPy2 += mu->py()*mu->py();
    sumPxPy += mu->px()*mu->py();
    if (mu->pt() > leadingMuPt) leadingMuPt = mu->pt();

    leadingMuons.push_back(mu->pt());
    leadingLeptons.push_back(mu->pt());
    leadingObjects.push_back(mu->pt());

    pMu += mu->p4();
    pLep += mu->p4();      
    pObj += mu->p4(); 

    MuE[ngoodmuons-1]   = mu->energy();
    MuPx[ngoodmuons-1]  = mu->px();
    MuPy[ngoodmuons-1]  = mu->py();
    MuPz[ngoodmuons-1]  = mu->pz(); 
    MuPt[ngoodmuons-1]  = mu->pt();
    MuEt[ngoodmuons-1]  = mu->et();
    MuEta[ngoodmuons-1] = mu->eta();      
    MuPhi[ngoodmuons-1] = mu->phi(); 
    }//muonID



  }


}  

//    if (ngoodmuons >0) No_of_Muons = No_of_Muons +1;
//    cout << evtno<< " --  Hellooo inside   "  << ngoodmuons <<  "    Number of Muons in the Event  are ==== " << No_of_Muons <<endl;
   
   //Sorting
   std::sort(leadingJets.begin(), leadingJets.end());
   std::sort(leadingElectrons.begin(), leadingElectrons.end());
   std::sort(leadingPhotons.begin(), leadingPhotons.end());
   std::sort(leadingMuons.begin(), leadingMuons.end());
   std::sort(leadingLeptons.begin(), leadingLeptons.end());
   std::sort(leadingObjects.begin(), leadingObjects.end());

   std::reverse(leadingJets.begin(), leadingJets.end());
   std::reverse(leadingElectrons.begin(), leadingElectrons.end());
   std::reverse(leadingPhotons.begin(), leadingPhotons.end());
   std::reverse(leadingMuons.begin(), leadingMuons.end());
   std::reverse(leadingLeptons.begin(), leadingLeptons.end());
   std::reverse(leadingObjects.begin(), leadingObjects.end());

   for (int i=0; i<4; ++i){
     JetArr[i] = leadingJets[i];
     EleArr[i] = leadingElectrons[i];
     PhArr[i] = leadingPhotons[i];
     MuArr[i] = leadingMuons[i];
     LeadingArr[i] = leadingObjects[i];
   } 

   float trace = sumPx2 + sumPy2;
   float det = sumPx2*sumPy2 - sumPxPy*sumPxPy;
   float lambda2 = (trace - sqrt(trace*trace - 4*det))/2.0;
   Sphericity = 2*lambda2/trace;
   
   mBH = 0.;
   
   mBH = pBH.M();
   
   //h_norm -> Fill(ST);
   for (size_t i=0; i<cutNames_.size(); ++i){
     std::map<std::string,TH1*>& histo = histos_[i];       
     histo["ST"]->Fill(ST);
     histo["HT"]->Fill(HT);
   }   
   
   NPV = nPVcount;
   NTracks = 0;
   NJets = ngoodjets;
   NElectrons = ngoodelectrons;
   NMuons = ngoodmuons; 	//not using muons.size() here as |dxy(bs)| < 0.2 might remove cosmics
   NPhotons = ngoodphotons; 	//not using photons.size() here as Swiss cross might remove spikes
   Multiplicity = ngoodjets + ngoodelectrons + ngoodphotons + ngoodmuons;
   
   // Attempts to reduce the size of output file and reject noise and unused events (Mult < 2)
   //if (isMCBH == false && Multiplicity < 2) return;
          
   NoScrap = noscrap; 
   isLeptonPhoton = (ngoodelectrons + ngoodmuons + ngoodphotons > 0);
   
   // Classify event by leading lepton/photon
   isEleChannel = (isLeptonPhoton && (leadingElePt > leadingMuPt) && (leadingElePt > leadingPhPt));
   isMuChannel = (isLeptonPhoton && (leadingMuPt > leadingElePt) && (leadingMuPt > leadingPhPt));
   isPhChannel = (isLeptonPhoton && (leadingPhPt > leadingElePt) && (leadingPhPt > leadingMuPt));
   
   //Filling objects "averaged" angles
   
   
   tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void SUSYAnalyserBSM::beginJob()
{
  cutNames_.push_back("No_Cut");
  
  tree = fs_->make<TTree>("t","t");
  //h_norm = new TH1F("h_norm","",500,0,50000);

  tree->Branch("NJets",&NJets,"NJets/I");
  tree->Branch("JetE",&JetE,"JetE[25]/F");
  tree->Branch("JetPx",&JetPx,"JetPx[25]/F");
  tree->Branch("JetPy",&JetPy,"JetPy[25]/F");
  tree->Branch("JetPz",&JetPz,"JetPz[25]/F");
  tree->Branch("JetPt",&JetPt,"JetPt[25]/F");
  tree->Branch("JetEt",&JetEt,"JetEt[25]/F");
  tree->Branch("JetEta",&JetEta,"JetEta[25]/F");
  tree->Branch("JetPhi",&JetPhi,"JetPhi[25]/F");

  tree->Branch("EleE" ,&EleE, "EleE[25]/F");
  tree->Branch("ElePx",&ElePx,"ElePx[25]/F");
  tree->Branch("ElePy",&ElePy,"ElePy[25]/F");
  tree->Branch("ElePz",&ElePz,"ElePz[25]/F");
  tree->Branch("ElePt",&ElePt,"ElePt[25]/F");
  tree->Branch("EleEt",&EleEt,"EleEt[25]/F");
  tree->Branch("EleEta",&EleEta,"EleEta[25]/F");
  tree->Branch("ElePhi",&ElePhi,"ElePhi[25]/F");
 
  tree->Branch("PhE",&PhE,"PhE[25]/F");
  tree->Branch("PhPx",&PhPx,"PhPx[25]/F");
  tree->Branch("PhPy",&PhPy,"PhPy[25]/F");
  tree->Branch("PhPz",&PhPz,"PhPz[25]/F");
  tree->Branch("PhPt",&PhPt,"PhPt[25]/F");
  tree->Branch("PhEt",&PhEt,"PhEt[25]/F");
  tree->Branch("PhEta",&PhEta,"PhEta[25]/F");
  tree->Branch("PhPhi",&PhPhi,"PhPhi[25]/F");

  tree->Branch("MuE",&MuE,"MuE[25]");
  tree->Branch("MuPx",&MuPx,"MuPx[25]");
  tree->Branch("MuPy",&MuPy,"MuPy[25]");
  tree->Branch("MuPz",&MuPz,"MuPz[25]");
  tree->Branch("MuPt",&MuPt,"MuPt[25]");
  tree->Branch("MuEt",&MuEt,"MuEt[25]");
  tree->Branch("MuEta",&MuEta,"MuEta[25]");
  tree->Branch("MuPhi",&MuPhi,"MuPhi[25]");
 
  tree->Branch("miniIso",&miniIso,"miniIso[25]"); 
 
  tree->Branch("ST",&ST,"ST/F");
  tree->Branch("HT",&HT,"HT/F");
  tree->Branch("mBH",&mBH,"mBH/F");
  tree->Branch("Met",&Met,"Met/F");
  tree->Branch("MetPx",&MetPx,"MetPx/F");
  tree->Branch("MetPy",&MetPy,"MetPy/F");
  tree->Branch("MetPhi",&MetPhi,"MetPhi/F");   
  tree->Branch("Sphericity", &Sphericity, "Sphericity/F");
  tree->Branch("Jet",JetArr,"JetArr[4]");   
  tree->Branch("Ele",EleArr,"EleArr[4]");
  tree->Branch("Mu",MuArr,"MuArr[4]");
  tree->Branch("Ph",PhArr,"PhArr[4]");   
  tree->Branch("NPV",&NPV,"NPV/I");   
  tree->Branch("NTracks",&NTracks,"NTracks/I");   
  tree->Branch("NElectrons",&NElectrons,"NElectrons/I");
  tree->Branch("NPhotons",&NPhotons,"NPhotons/I");
  tree->Branch("NMuons",&NMuons,"NMuons/I");
  tree->Branch("NoScrap", &NoScrap, "NoScrap/B");     
  
  tree->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
  tree->Branch("isLeptonPhoton", &isLeptonPhoton, "isLeptonPhoton/B");   
  tree->Branch("isEleChannel", &isEleChannel, "isEleChannel/B");
  tree->Branch("isMuChannel", &isMuChannel, "isMuChannel/B"); 
  tree->Branch("isPhChannel", &isPhChannel, "isPhChannel/B");
  
 
  tree->Branch("runno",&runno,"runno/I");
  tree->Branch("evtno",&evtno,"evtno/I"); 
  tree->Branch("lumiblock",&lumiblock,"lumiblock/I");               
  tree->Branch("isRealData",&isRealData,"isRealData/I");
  tree->Branch("muon_d0",&muon_d0,"muon_d0/F");                                                      
  
  tree->Branch("firedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1",&firedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1,"firedHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1/O");
  tree->Branch("firedHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1",&firedHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1,"firedHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1/O");
  tree->Branch("firedHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1",&firedHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1,"firedHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1/O");
  tree->Branch("firedHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1",&firedHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1,"firedHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1/O");
  tree->Branch("firedHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1",&firedHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1,"firedHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v1/O");
  tree->Branch("firedHLT_DoubleMu8_Mass8_PFHT300_v1",&firedHLT_DoubleMu8_Mass8_PFHT300_v1,"firedHLT_DoubleMu8_Mass8_PFHT300_v1/O");
  tree->Branch("firedHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v1",&firedHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v1,"firedHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v1/O");
  tree->Branch("firedHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v1",&firedHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v1,"firedHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v1/O");

  tree->Branch("firedHLT_TripleMu_12_10_5_v1",&firedHLT_TripleMu_12_10_5_v1,"firedHLT_TripleMu_12_10_5_v1/O");
  tree->Branch("firedHLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1",&firedHLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1,"firedHLT_DiMu9_Ele9_CaloIdL_TrackIdL_v1/O"); 
  tree->Branch("firedHLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1",&firedHLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1,"firedHLT_Mu8_DiEle12_CaloIdL_TrackIdL_v1/O");
  tree->Branch("firedHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1",&firedHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1,"firedHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v1/O");

  tree->Branch("firedHLT_Ele8_CaloIdM_TrackIdM_PFJet30_v1",&firedHLT_Ele8_CaloIdM_TrackIdM_PFJet30_v1,"firedHLT_Ele8_CaloIdM_TrackIdM_PFJet30_v1/O");
  tree->Branch("firedHLT_Ele12_CaloIdM_TrackIdM_PFJet30_v1",&firedHLT_Ele12_CaloIdM_TrackIdM_PFJet30_v1,"firedHLT_Ele12_CaloIdM_TrackIdM_PFJet30_v1/O");
  tree->Branch("firedHLT_Ele18_CaloIdM_TrackIdM_PFJet30_v1",&firedHLT_Ele18_CaloIdM_TrackIdM_PFJet30_v1,"firedHLT_Ele18_CaloIdM_TrackIdM_PFJet30_v1/O");
  tree->Branch("firedHLT_Ele23_CaloIdM_TrackIdM_PFJet30_v1",&firedHLT_Ele23_CaloIdM_TrackIdM_PFJet30_v1,"firedHLT_Ele23_CaloIdM_TrackIdM_PFJet30_v1/O");
  tree->Branch("firedHLT_Ele33_CaloIdM_TrackIdM_PFJet30_v1",&firedHLT_Ele33_CaloIdM_TrackIdM_PFJet30_v1,"firedHLT_Ele33_CaloIdM_TrackIdM_PFJet30_v1/O");
  tree->Branch("firedHLT_Mu8_v1",&firedHLT_Mu8_v1,"firedHLT_Mu8_v1/O");
  tree->Branch("firedHLT_Mu17_v1",&firedHLT_Mu17_v1,"firedHLT_Mu17_v1/O");
  tree->Branch("firedHLT_Mu24_v1",&firedHLT_Mu24_v1,"firedHLT_Mu24_v1/O");
  tree->Branch("firedHLT_Mu34_v1",&firedHLT_Mu34_v1,"firedHLT_Mu34_v1/O");

  
    
  for (size_t i=0; i<cutNames_.size(); ++i)
    createHistogram(cutNames_[i]);
  
}

void SUSYAnalyserBSM::createHistogram(const std::string& folderName){
  TFileDirectory subDir = fs_->mkdir(folderName);
  std::map<std::string, TH1*> container;
  
  container["ST"] = subDir.make<TH1F>("ST", "ST", 500, 0, 100000);
  container["HT"] = subDir.make<TH1F>("HT", "HT", 500, 0, 100000); 
  histos_.push_back(container);          
}

bool SUSYAnalyserBSM::check(std::string process, std::string pCheck) {
  bool value= false;
  
  if (process == pCheck) value = true;
  
  return value;
}


// ------------ method called once each job just after ending the event loop  ------------
void SUSYAnalyserBSM::endJob()
{
  //h_norm->Write();
}

void
SUSYAnalyserBSM::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//The following says we do not know what parameters are allowed so do no validation
// Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
   }


//define this as a plug-in
DEFINE_FWK_MODULE(SUSYAnalyserBSM);
