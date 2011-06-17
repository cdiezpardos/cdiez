/* \class WMuNuAnalyzer
 *  WMuNu Analyzer modified to obtain the "MC truth" efficiency for a single muon
 *  Specially devoted to MET & Bckg studies
 *  \author M.I. Josa
 */

#include "master.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
  

#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#include <map>
#include <vector>

#define MZ 91.1876
#define MMuon 0.105658

class TFile;
class TH1F;
class TH2F;
double MassMuon=0.105658369;

class WMuNu : public edm::EDAnalyzer {
public:
      typedef math::XYZPoint Point;
  WMuNu(const edm::ParameterSet& pset);
  virtual ~WMuNu();
  virtual void beginJob();
  bool GoodEWKMuon(reco::MuonRef Mu);


 const std::vector<reco::MuonRef> selectMuons(edm::Handle<reco::MuonCollection> muonCollection);
 bool compMomentum(reco::MuonRef m1, reco::MuonRef m2);
 double computeMass(reco::MuonRef m1, reco::MuonRef m2);
 double computeMass(reco::MuonRef m1, reco::MuonRef m2,double dpt1, double dpt2);
/*
          const reco::TrackToTrackMap * tevMap_default;
          const reco::TrackToTrackMap * tevMap_1stHit;
          const reco::TrackToTrackMap * tevMap_picky;
          const reco::TrackToTrackMap * tevMap_dyt;
*/

  bool FoundZ(edm::Handle<reco::MuonCollection>& muonCollection, edm::Handle<reco::TrackCollection>& trackCollection);
  virtual void endJob();
  virtual void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
private:
  edm::InputTag muonTag_;
  edm::InputTag metTag_;
  edm::InputTag jetTag_;

  edm::InputTag triggerSummaryLabel_;
  edm::InputTag hltFilterTag_;

  bool useOnlyGlobalMuons_;
  bool useMuonQualityCuts_;
  double ptCut_;
  double etaCut_;
  bool isTotIso_;
  double isoCut03_;
  double acopCut_;

// Histograms
  std::map<std::string,TH1D*> h1_;
  std::map<std::string,TH2D*> h2_;

  std::ostringstream ss;

  // Root output file


  unsigned int numberOfEvents;
  unsigned int numberOfMuons;
  unsigned int Nneutrinos;
};

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;

/// Constructor
WMuNu::WMuNu(const ParameterSet& pset) :
      muonTag_(pset.getUntrackedParameter<edm::InputTag> ("MuonTag", edm::InputTag("muons"))),
      useOnlyGlobalMuons_(pset.getUntrackedParameter<bool>("UseOnlyGlobalMuons", true)),
      useMuonQualityCuts_(pset.getUntrackedParameter<bool>("UseMuonQualityCuts")),
      ptCut_(pset.getUntrackedParameter<double>("PtCut", 20.)),
      etaCut_(pset.getUntrackedParameter<double>("EtaCut", 2.1)),
      isTotIso_(pset.getUntrackedParameter<bool>("IsTotIso", true)),
      isoCut03_(pset.getUntrackedParameter<double>("IsoCut03", 0.15)),
      acopCut_(pset.getUntrackedParameter<double>("AcopCut", 999999.)),
      hltFilterTag_(pset.getParameter<edm::InputTag>("hltL3FilterLabel")),
      triggerSummaryLabel_(pset.getParameter<edm::InputTag>("triggerSummaryLabel"))


{
  LogDebug("WMuNuAnalyzer")<<" WMuNuAnalyzer constructor called";
}

/// Destructor
WMuNu::~WMuNu(){
}

void WMuNu::beginJob(){
  // Create output files

  edm::Service<TFileService> fs;


  numberOfEvents = 0;
  numberOfMuons = 0;

  h2_["hEtaMu1Mu2_gen"]=fs->make<TH2D>("hEtaMu1Mu2_gen","#eta Mu1 vs #eta Mu2",100,-3.5,3.5,100,-3.5,3.5);
  h2_["hEtaPtMu_gen"]=fs->make<TH2D>("hEtaPtMu_gen","#eta vs #pt gen; #eta; p_{t} GeV",70,-2.1,2.1, 200,0.,200.);
  h2_["hEtaPtMu_sel"]=fs->make<TH2D>("hEtaPtMu_sel","#eta vs #pt sel; #eta; p_{t} GeV",70,-2.1,2.1,200,0.,200.);

  h1_["hNMu"]=fs->make<TH1D>("NMu","Nb. muons in the event",10,0.,10.);
  h1_["hNMu_trigger"]=fs->make<TH1D>("NMu_trigger","Nb. muons after trigger",10,0.,10.);
  h1_["hNMu_sel"]=fs->make<TH1D>("NMu_sel","Nb. muons in the event",10,0.,10.);

  h1_["hPtMu"]=fs->make<TH1D>("ptMu","Pt mu",1000,0.,1000.);
  h1_["hPtMu_trigger"]=fs->make<TH1D>("ptMu_trigger","Pt mu",1000,0.,1000.);
  h1_["hPtMu_gen"]=fs->make<TH1D>("ptMu_gen","Pt mu",1000,0.,1000.);
  h1_["hPtMu_sel"]=fs->make<TH1D>("ptMu_sel","Pt mu",1000,0.,1000.);
  h1_["hPtMu_PtEtacut"]=fs->make<TH1D>("ptMu_PtEtacut","Pt mu",1000,0.,1000.);
  h1_["hPtMu_aisl"]=fs->make<TH1D>("ptMu_aisl","Pt mu",1000,0.,1000.);
  h1_["hPtMu_qual"]=fs->make<TH1D>("ptMu_qual","Pt mu",1000,0.,1000.);
  h1_["hPtMu_reco"]=fs->make<TH1D>("ptMu_reco","Pt mu",1000,0.,1000.);
          double xbin[7]={-2.1,-1.2,-0.9,0.,0.9,1.2,2.1};
  h1_["hEtaMu"]=fs->make<TH1D>("etaMu","Eta mu",6,xbin);
  h1_["hEtaMu_gen"]=fs->make<TH1D>("etaMu_gen","Eta mu",6,xbin);
  h1_["hEtaMu_sel"]=fs->make<TH1D>("etaMu_sel","Eta mu",6,xbin);
  h1_["hEtaMu_trigger"]=fs->make<TH1D>("etaMu_trigger","Eta mu",6,xbin);
  h1_["hEtaMu_PtEtacut"]=fs->make<TH1D>("etaMu_PtEtacut","Eta mu",6,xbin);
  h1_["hEtaMu_aisl"]=fs->make<TH1D>("etaMu_aisl","Eta mu",6,xbin);
  h1_["hEtaMu_qual"]=fs->make<TH1D>("etaMu_qual","Eta mu",6,xbin);
  h1_["hEtaMu_reco"]=fs->make<TH1D>("etaMu_reco","Eta mu",6,xbin);


  h1_["hPhiMu"]=fs->make<TH1D>("phiMu","Pt mu",100,-TMath::Pi(), TMath::Pi());
  h1_["hPhiMu_gen"]=fs->make<TH1D>("phiMu_gen","Pt mu",100,-TMath::Pi(), TMath::Pi());
  h1_["hPhiMu_sel"]=fs->make<TH1D>("phiMu_sel","Pt mu",100,-TMath::Pi(), TMath::Pi());
  h1_["hPhiMu_trigger"]=fs->make<TH1D>("phiMu_trigger","Pt mu",100,-TMath::Pi(), TMath::Pi());
  h1_["hPhiMu_PtEtacut"]=fs->make<TH1D>("phiMu_PtEtacut","Pt mu",100,-TMath::Pi(), TMath::Pi());
  h1_["hPhiMu_aisl"]=fs->make<TH1D>("phiMu_aisl","Pt mu",100,-TMath::Pi(), TMath::Pi());
  h1_["hPhiMu_qual"]=fs->make<TH1D>("phiMu_qual","Pt mu",100,-TMath::Pi(), TMath::Pi());
  h1_["hPhiMu_reco"]=fs->make<TH1D>("phiMu_reco","Pt mu",100,-TMath::Pi(), TMath::Pi());


}

void WMuNu::endJob(){
  LogVerbatim("") << "WMuNuAnalyzer>>> FINAL PRINTOUTS -> BEGIN";
  LogVerbatim("") << "WMuNuAnalyzer>>> Number of analyzed events= " << numberOfEvents;
  LogVerbatim("") << "WMuNuAnalyzer>>> Number of analyzed muons= " << numberOfMuons;
  LogVerbatim("") << "WMuNuAnalyzer>>> FINAL PRINTOUTS -> END";
}

bool WMuNu::GoodEWKMuon(reco::MuonRef Mu){
     bool goodMuon=true;

     int validmuonhits=Mu->globalTrack()->hitPattern().numberOfValidMuonHits();
     int    nMatches = Mu->numberOfMatches();
//     double d0 = Mu->globalTrack()->dxy(beamSpotHandle->position());
     double chi2=Mu->globalTrack()->normalizedChi2();
     double numberOfValidHits=Mu->globalTrack()->numberOfValidHits(); // el criterio de calidad los pide en el silicon track, no en el global


     goodMuon = (chi2<10.)&&(numberOfValidHits>=11)&&(validmuonhits>0)&&(nMatches>1);

         
    if (!Mu->isGlobalMuon()||!Mu->isTrackerMuon()){goodMuon=false;}    
    return goodMuon;
}

const std::vector<reco::MuonRef> WMuNu::selectMuons(Handle<reco::MuonCollection> muonCollection){
  std::vector<reco::MuonRef> muones;
 

 for (unsigned int i=0; i<muonCollection->size(); i++) {
    reco::MuonRef mu(muonCollection,i);
    //CORTE CALIDAD
    if ( (useOnlyGlobalMuons_ && !mu->isGlobalMuon()) ) {cout << "no es global" << endl; continue;}
    if ( (useMuonQualityCuts_ && !GoodEWKMuon(mu)) ) {cout<<"no calidad"<<endl; continue;}


    //CORTE CINEMATICO
    if (mu->pt()<ptCut_ || fabs(mu->eta())>etaCut_) {cout<<"corte cinematico"<<endl;continue;}
    //CORTE AISLAMIENTO
    double ptsum = mu->isolationR03().sumPt;
    double cal=mu->isolationR03().emEt+mu->isolationR03().hadEt;
    if (isTotIso_==true){
      if ( (cal+ptsum)/mu->pt()> isoCut03_ ) continue; }   // Depends on what kind of iso you are applying, this changes
    else {     if ( ptsum/mu->pt()>isoCut03_ ) {cout<<"isolation"<<endl;continue;  }}
    
    muones.push_back(mu);
  }
  return muones;
}



double WMuNu::computeMass(reco::MuonRef m1, reco::MuonRef m2) {
  const math::XYZTLorentzVector Muon4_1 (m1->px(), m1->py() , m1->pz(), m1->p());
  const math::XYZTLorentzVector Muon4_2 (m2->px(), m2->py() , m2->pz(), m2->p());
  const math::XYZTLorentzVector myZCand=Muon4_1+Muon4_2;
  return myZCand.mass();
}

void WMuNu::analyze(const Event & event, const EventSetup& eventSetup){

   double pt_sel[15];
   double eta_sel[15];
   double phi_sel[15];
   bool event_sel = true;


   // Get the Collections from file
   Handle<reco::MuonCollection> muonCollection;
   if (event.getByLabel(muonTag_, muonCollection)) {
      LogTrace("")<<"Reconstructed Muon tracks: " << muonCollection->size() << endl;
   } else {
      LogTrace("") << ">>> Muon collection does not exist !!!";
      return;
   }
/*
  Handle<reco::TrackCollection> trackCollection;
   if (!event.getByLabel("generalTracks", trackCollection)) {
      LogTrace("") << ">>> Track collection does not exist !!!";
      return;
   }

*/
   Handle<TriggerResults> triggerResults;
    if (event.getByLabel(InputTag("TriggerResults::HLT"), triggerResults)) {
    }   else{
                    LogError("") << ">>> TRIGGER collection does not exist !!!";
                    return;
    }
    const edm::TriggerNames trigNames = event.triggerNames(*triggerResults);
    edm::Handle<TriggerEvent> triggerObj;
     event.getByLabel(triggerSummaryLabel_,triggerObj);
     if(!triggerObj.isValid()) {
          LogTrace("") << "Summary HLT objects not found, "
             "skipping event";
             return;
     }


   Handle<reco::GenParticleCollection> genParticles;
   if(!event.getByLabel( "prunedGenParticles", genParticles )){
   LogTrace("")<<"genParticles dont exist!!!"<<endl;
   return;
   }
/*   reco::BeamSpot vertexBeamSpot;
//   reco::BeamSpot offlineBeamSpot;
   edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
   if(event.getByType(recoBeamSpotHandle)){
   vertexBeamSpot = *recoBeamSpotHandle;

   }
   else { LogTrace("") << ">>> BEAMSPOT collection does not exist !!!";
      return; 
   }
*/
/*//------------------
   bool fired = false;
   int itrig1 = trigNames.triggerIndex("HLT_Mu15");
   if (triggerResults->accept(itrig1)){  fired = true;}
   else{ LogTrace("")<<" Event did not fire the HLT_Mu15 trigger!!"; return;}
///------------------ */
      // Beam spot
      Handle<reco::BeamSpot> beamSpotHandle;
      if (!event.getByLabel(InputTag("offlineBeamSpot"), beamSpotHandle)) {
            LogTrace("") << ">>> No beam spot found !!!";
            return;
      }
/*

          edm::Handle<reco::TrackToTrackMap> tevMapH_default;
          edm::Handle<reco::TrackToTrackMap> tevMapH_1stHit;
          edm::Handle<reco::TrackToTrackMap> tevMapH_picky;
          edm::Handle<reco::TrackToTrackMap> tevMapH_dyt;


          event.getByLabel("tevMuons", "default", tevMapH_default);
          tevMap_default = tevMapH_default.product();

          event.getByLabel("tevMuons", "firstHit", tevMapH_1stHit);
          tevMap_1stHit = tevMapH_1stHit.product();

          event.getByLabel("tevMuons", "picky", tevMapH_picky);
          tevMap_picky = tevMapH_picky.product();

          event.getByLabel("tevMuons", "dyt", tevMapH_dyt);
          tevMap_dyt = tevMapH_dyt.product();

          TrackToTrackMap::const_iterator iTeV_default;
          TrackToTrackMap::const_iterator iTeV_1stHit;
          TrackToTrackMap::const_iterator iTeV_picky;
          TrackToTrackMap::const_iterator iTeV_dyt;
*/


   numberOfEvents++;
   numberOfMuons+=muonCollection->size();


    float etagen= -99., ptgen =-99., phigen =-99.;
    bool genfound = false;


      //loop sobre muones generados
      double eta_GEN[15]={0}, pt_GEN[15]={0}, phi_GEN[15]={0},  px_GEN[15]={0}, py_GEN[15]={0}, pz_GEN[15]={0}, p_GEN[15]={0};
      unsigned int ngen = 0;
        for(unsigned int i = 0; i < genParticles->size(); ++i) {
            const GenParticle& p = (*genParticles)[i];
            int id = p.pdgId();
            int st = p.status();
   
            if (st==3) continue;
            if (abs(id)!=13) continue;

            int nmothers = p.numberOfMothers();
            //printf("Boson Id %d numberOfMothers %d\n", id, nmothers);
            if (nmothers!=1) continue;
            const Candidate * mom = p.mother();
            int momId=mom->pdgId();

            while(abs(momId)==13 && mom->numberOfMothers()==1){
                  mom = mom->mother();
                  momId=mom->pdgId();
            }
            //if(abs(momId) !=23 ) continue;

            ptgen = p.pt();
            etagen =p.eta();

            if (ptgen<20.  || fabs(etagen)>etaCut_) continue;
            ngen++; //if(ngen>2) continue;
            genfound = true;
            eta_GEN[ngen-1]=etagen;
            phi_GEN[ngen-1]=p.phi();
            //cout <<" ptgen "<<ptgen<<" "<< pt_GEN[ngen-1]<<endl;
            //break;
            h1_["hPtMu_gen"]->Fill(ptgen);
            if (ptgen>=20.) {
            h1_["hEtaMu_gen"]->Fill(etagen);
            h1_["hPhiMu_gen"]->Fill(p.phi()); }
            h2_["hEtaPtMu_gen"]->Fill(etagen,ptgen);
     
//Empieza el loop sobre variables reconstruidas--------------------

  //variables
//  vector<MuonRef> muones;
  double Zm_aux=0;
//  MuonRef muon1, muon2;
  
  //Z mass = 91.1876GeV
  double numSel = 0;  
  for (unsigned int i=0; i<muonCollection->size(); i++) {
    bool muon_sel = true;
    reco::MuonRef mu1(muonCollection,i);

    Geom::Phi<double> deltaphi(p.phi()-mu1->phi());
    double  delta= sqrt(pow(deltaphi.value(),2)+pow(etagen-mu1->eta(),2));
    if (delta > 0.5) continue;
    numSel++;

    h1_["hEtaMu"]->Fill(mu1->eta());
    h1_["hPtMu"]->Fill(mu1->pt());
    h1_["hPhiMu"]->Fill(mu1->phi());

//      double pt = (*cocktail).pt();
      double pt = mu1->pt();
      double eta = mu1->eta();


      if (pt<ptCut_ || fabs(eta)>etaCut_) {muon_sel = false; continue;}
      LogTrace("") << "\t... eta= " << eta;
      LogTrace("") << "\t... pt= " << pt << " GeV";

      if (muon_sel == false) continue;
      //muones que pasan cortes calidad+ son globales + cortes pt, eta
      h1_["hPtMu_PtEtacut"]->Fill(ptgen);
      h1_["hEtaMu_PtEtacut"]->Fill(etagen);
      h1_["hPhiMu_PtEtacut"]->Fill(p.phi());


     if( !mu1->isGlobalMuon() || !mu1->isTrackerMuon() ) continue;
     double d0 = mu1->globalTrack()->dxy(beamSpotHandle->position());
     double chi2=mu1->globalTrack()->normalizedChi2();
     double numberOfValidHits=mu1->globalTrack()->hitPattern().numberOfValidTrackerHits(); // el criterio de calidad los pide en el silicon track, no en el global

     int pixelHits = mu1->globalTrack()->hitPattern().numberOfValidPixelHits();

     int validmuonhits=mu1->globalTrack()->hitPattern().numberOfValidMuonHits();
     int    nMatches = mu1->numberOfMatches();
      if((numberOfValidHits <11)||validmuonhits<1 || nMatches <= 1 || pixelHits <1 ) continue;




      //histogramas de muones Globales+corte calidad
      h1_["hPtMu_reco"]->Fill(ptgen);
      h1_["hPhiMu_reco"]->Fill(p.phi());
      h1_["hEtaMu_reco"]->Fill(etagen);
      LogTrace("") << "\t... pt= " << pt << " GeV";

      if(fabs(d0)>=0.2 || chi2>=10.) continue;   
      h1_["hPtMu_qual"]->Fill(ptgen);
      h1_["hPhiMu_qual"]->Fill(p.phi());
      h1_["hEtaMu_qual"]->Fill(etagen);



//CORTE AISLAMIENTO
    double ptsum = mu1->isolationR03().sumPt;
    double cal=mu1->isolationR03().emEt+mu1->isolationR03().hadEt;
    if (isTotIso_==true){
      if ( (cal+ptsum)/mu1->pt()> isoCut03_ ) continue; }   // Depends on what kind of iso you are applying, this changes
    else {     if ( ptsum/mu1->pt()>isoCut03_ ) {continue;  }}

      //muones todo lo anterior y aislados
      h1_["hPtMu_aisl"]->Fill( ptgen);
      h1_["hPhiMu_aisl"]->Fill( p.phi());
      h1_["hEtaMu_aisl"]->Fill( etagen);

      bool haytrigger = false;
      double deltatrigmax = 1.;
      
      const trigger::TriggerObjectCollection & toc(triggerObj->getObjects());
      for ( size_t ia = 0; ia < triggerObj->sizeFilters(); ++ ia) {
//      cout << triggerObj->filterTag(ia) <<endl;
      if( triggerObj->filterTag(ia)  == hltFilterTag_) {
    //  cout << triggerObj->filterTag(ia) <<endl;
      const trigger::Keys & k = triggerObj->filterKeys(ia);
       for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
       double pttrig= toc[*ki].pt();
       double etatrig=toc[*ki].eta();
       double phitrig=toc[*ki].phi();
      Geom::Phi<double> deltaphi(phitrig-mu1->phi());
       if(sqrt(pow(etatrig-mu1->eta(),2)+pow(deltaphi.value(),2)) < 0.2 ) {haytrigger=true; break;}
       }
      }
      }
      if (haytrigger==false) continue;

//        h1_["hPtMu_sel"]->Fill(mu1->pt());

         h1_["hPtMu_sel"]->Fill(ptgen);
         if (p.pt()>=ptCut_) {
         h1_["hEtaMu_sel"]->Fill(etagen);
         h1_["hPhiMu_sel"]->Fill(p.phi());}
         h2_["hEtaPtMu_sel"]->Fill(etagen,ptgen);
         h1_["hNMu_sel"]->Fill(numSel);

  }//for i
  }//end generator loop
  
}


DEFINE_FWK_MODULE(WMuNu);

