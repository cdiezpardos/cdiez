/* class ZMuMuTrackingAnalyzer *  ZMuMuTracking Analyzer modified t
 include *many* more plots than the standard 
 *  Specially devoted to MET & Bckg studies
 *  \author C. Diez
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

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMinuit.h"
#include "TFitterMinuit.h"
#include "Minuit2/FCNBase.h"
//#include "utils.h"
#include "TH1D.h"
#include "TTree.h"

#include <map>
#include <vector>

#define MZ 91.1876
#define MMuon 0.105658
//--------------------------------------------------------------------
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include <fstream>
#include <iostream>

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"

#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"


class TFile;
class TH1F;
class TH2F;
class TH2D;
//double MassMuon=0.105658369;

class ZMuMuTracking : public edm::EDAnalyzer {
public:
      typedef math::XYZPoint Point;
  ZMuMuTracking(const edm::ParameterSet& pset);
  virtual ~ZMuMuTracking();
	  virtual void beginJob();
	  bool GoodEWKMuon(reco::MuonRef Mu);
	     // PathInfoCollection hltPaths_;

	 const std::vector<reco::MuonRef> selectMuons2(edm::Handle<reco::MuonCollection> muonCollection);

	  bool FoundZ(edm::Handle<reco::MuonCollection>& muonCollection, edm::Handle<reco::TrackCollection>& trackCollection);
	  virtual void endJob();
	  virtual void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
	private:

  TTree * Trackingtree;
  double pt, eta, Isolation;
  int PassCut, charge;


	  edm::InputTag muonTag_;
	  edm::InputTag jetTag_;

	//------------------------
	  edm::InputTag triggerSummaryLabel_;
	  edm::InputTag hltFilterTag_;
          edm::InputTag hltFilterTagMu11_;
          edm::InputTag hltFilterTagMu15_;

	  bool useOnlyGlobalMuons_;
	  bool useMuonQualityCuts_;
	  double ptCut_;
	  double etaCut_;
	  bool isTotIso_;
	  double isoCut03_;
	  double massTMin_;
	  double massTMax_;
	  double acopCut_;
/*
          const reco::TrackToTrackMap * tevMap_default;
          const reco::TrackToTrackMap * tevMap_1stHit;
          const reco::TrackToTrackMap * tevMap_picky;
          const reco::TrackToTrackMap * tevMap_dyt;
*/

	// Histograms
	  std::map<std::string,TH1D*> h1_;
	  std::map<std::string,TH2D*> h2_;



	  unsigned int numberOfEvents;
	  unsigned int numberOfMuons;

         TFile* ff;
         TH1D* h; 
         TH1D* hMassajuste;
         TH1D* hZ_MC;
//         TFitterMinuit *minimizer;
	};

	using namespace std;
	using namespace edm;
	using namespace reco;
	using namespace trigger;
	/// Constructor
	ZMuMuTracking::ZMuMuTracking(const ParameterSet& pset) :
	      muonTag_(pset.getUntrackedParameter<edm::InputTag> ("MuonTag", edm::InputTag("muons"))),
	      useOnlyGlobalMuons_(pset.getUntrackedParameter<bool>("UseOnlyGlobalMuons", true)),
	      useMuonQualityCuts_(pset.getUntrackedParameter<bool>("UseMuonQualityCuts")),
	      ptCut_(pset.getUntrackedParameter<double>("PtCut", 20.)),
	      etaCut_(pset.getUntrackedParameter<double>("EtaCut", 2.1)),
	      isTotIso_(pset.getUntrackedParameter<bool>("IsTotIso", true)),
	      isoCut03_(pset.getUntrackedParameter<double>("IsoCut03", 0.15)),
	      acopCut_(pset.getUntrackedParameter<double>("AcopCut", 999999.)),
	      triggerSummaryLabel_(pset.getParameter<edm::InputTag>("triggerSummaryLabel")),
              hltFilterTag_(pset.getParameter<edm::InputTag>("hltL3FilterLabel")),
              hltFilterTagMu11_(pset.getParameter<edm::InputTag>("hltL3FilterLabelMu11")),
              hltFilterTagMu15_(pset.getParameter<edm::InputTag>("hltL3FilterLabelMu15"))
         //     hltFilterTag_(pset.getUntrackedParameter<std::vector <std::string> >("hltL3FilterLabel")),


	{
	  LogDebug("ZMuMuTrackingAnalyzer")<<" ZMuMuTrackingAnalyzer constructor called";
	}

	/// Destructor
	ZMuMuTracking::~ZMuMuTracking(){
	}

	void ZMuMuTracking::beginJob(){
	  // Create output files

	  edm::Service<TFileService> fs;
          Trackingtree = fs->make<TTree>("Trackingtree","check pt-Zmumu relation");
          Trackingtree->Branch("Isolation",&Isolation,"Isolation/D");
          Trackingtree->Branch("ptmuon",&pt,"pt/D");
          Trackingtree->Branch("eta",&eta,"eta/D");
          Trackingtree->Branch("charge",&charge,"charge/I"); 
          Trackingtree->Branch("PassCut",&PassCut,"PassCut/I");
 
         ff = TFile::Open("Example.root", "RECREATE");
         ff->cd();
         h = new TH1D("h","",100,40.,140.);
         hMassajuste = new TH1D("hMassajuste","",100,40.,140.);
	 hZ_MC = new TH1D("hZ_MC","",100,40.,140.);


	  numberOfEvents = 0;
	  numberOfMuons = 0;
          double xbin[7]={-2.1,-1.2,-0.9,0.,0.9,1.2,2.1};
	  h1_["hPtMu_DENtr"]=fs->make<TH1D>("ptMu_DENtr","Pt mu",100,0.,100.);
	  h1_["hPhiMu_DENtr"]=fs->make<TH1D>("phiMu_DENtr","Phi mu",24,-TMath::Pi(), TMath::Pi());
          h1_["hEtaMu_DENtr"]=fs->make<TH1D>("etaMu_DENtr","Eta mu",6,xbin);
//	  h1_["hEtaMu_DENtr"]=fs->make<TH1D>("etaMu_DENtr","Eta mu",100,-2.5,2.5);
	  h1_["hPtMu_NUMtr"]=fs->make<TH1D>("ptMu_NUMtr","Pt mu",100,0.,100.);
	  h1_["hPhiMu_NUMtr"]=fs->make<TH1D>("phiMu_NUMtr","Phi mu",24,-TMath::Pi(), TMath::Pi());
         h1_["hEtaMu_NUMtr"]=fs->make<TH1D>("etaMu_NUMtr","Eta mu",6,xbin);

	}

	void ZMuMuTracking::endJob(){
	  LogVerbatim("") << "ZMuMuTrackingAnalyzer>>> FINAL PRINTOUTS -> BEGIN";
	  LogVerbatim("") << "ZMuMuTrackingAnalyzer>>> Number of analyzed events= " << numberOfEvents;
	  LogVerbatim("") << "ZMuMuTrackingAnalyzer>>> Number of analyzed muons= " << numberOfMuons;
	  LogVerbatim("") << "ZMuMuTrackingAnalyzer>>> FINAL PRINTOUTS -> END";

	}

	bool ZMuMuTracking::GoodEWKMuon(reco::MuonRef Mu){
	     bool goodMuon=true;
	     
	     int validmuonhits=Mu->globalTrack()->hitPattern().numberOfValidMuonHits();
             int pixelHits = Mu->globalTrack()->hitPattern().numberOfValidPixelHits(); 
	     double chi2=Mu->globalTrack()->normalizedChi2();
	     double numberOfValidHits=Mu->globalTrack()->hitPattern().numberOfValidTrackerHits();	 
            int    nMatches = Mu->numberOfMatches();


	    goodMuon = (chi2<10.)&&(numberOfValidHits>=11)&&(validmuonhits>0)&&(nMatches>1)&&(pixelHits>0);
	    if (!Mu->isGlobalMuon()||!Mu->isTrackerMuon()) goodMuon=false;    
	    return goodMuon;
	}

	//subrutina para seleccionar muones-pt, eta ,(iso) como preseleccion para la eficiencia de tracker
	const std::vector<reco::MuonRef> ZMuMuTracking::selectMuons2(Handle<reco::MuonCollection> muonCollection){
	  std::vector<reco::MuonRef> muones;
	 
	  for (unsigned int i=0; i<muonCollection->size(); i++) {
	    reco::MuonRef mu1(muonCollection,i);

	//      if(!mu1->isStandAloneMuon()) continue;
	    //CORTE CINEMATICO
	// pt
	      double pt = mu1->pt();
	      LogTrace("") << "\t... pt= " << pt << " GeV";

	      if (pt<ptCut_) continue;
	// eta
	      double eta = mu1->eta();
	      LogTrace("") << "\t... eta= " << eta;
	      if (fabs(eta)>etaCut_) continue;
/*
	//CORTE AISLAMIENTO
	    double ptsum = mu1->isolationR03().sumPt;
	    double cal=mu1->isolationR03().emEt+mu1->isolationR03().hadEt;
	    if (isTotIso_==true){
	      if ( (cal+ptsum)/mu1->pt()> isoCut03_ ) continue; }   // Depends on what kind of iso you are applying, this changes
	      else {     if ( ptsum/mu1->pt()> isoCut03_ ) {continue;  }
	    }
*/	    
	    muones.push_back(mu1);    
	  }//for i

	  return muones;
	}


	void ZMuMuTracking::analyze(const Event & event, const EventSetup& eventSetup){

	   // Get the Collections from file
	   Handle<reco::MuonCollection> muonCollection;
	   if (event.getByLabel(muonTag_, muonCollection)) {
	      LogTrace("")<<"Reconstructed Muon tracks: " << muonCollection->size() << endl;
	   } else {
	      LogTrace("") << ">>> Muon collection does not exist !!!";
	      return;
	   }

	  Handle<reco::TrackCollection> trackCollection;
	   if (!event.getByLabel("generalTracks", trackCollection)) {
	      LogTrace("") << ">>> Track collection does not exist !!!";
	      return;
	   }


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

	   numberOfEvents++;
	   numberOfMuons+=muonCollection->size();

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

	LogTrace("")<<">>> Number of muons in the Event: "<<muonCollection->size() <<" Starting Loop over Muons!!";

	reco::MuonRef SelectedMuon; 

	//variables
	std::vector<reco::TrackRef> trazas;
	double Zm_aux=0;
	TrackRef track1;

	std::vector<reco::MuonRef> muones;
	reco::MuonRef muon1, muon2;

//-----------------------------------
//Tracking efficiency---------------------muons with pt>25,eta<2, iso?
//NOTA: Tal como esta corresponde a lo que os ensenie, si tengo un standAloneMuon "probe", le pido que sea tambien TrackerMuon, PERO NO que este asociado a un traza--- modificarlo, anyadiendo la referencia a la traza... 

    muones = selectMuons2(muonCollection);
    for (std::vector<reco::MuonRef>::const_iterator m1=muones.begin(); m1!=muones.end(); ++m1){

      muon1= (*m1);
    //CORTE CALIDAD
    if ( (useOnlyGlobalMuons_ && !muon1->isGlobalMuon()) ) continue;
    if ( (useMuonQualityCuts_ && !GoodEWKMuon(muon1)) ) continue; 
    double d0 = muon1->globalTrack()->dxy(beamSpotHandle->position()); if(fabs(d0)>=0.02) continue;

            double ptsum = muon1->isolationR03().sumPt;
            double cal=muon1->isolationR03().emEt+muon1->isolationR03().hadEt;
            if (isTotIso_==true){
              if ( (cal+ptsum)/muon1->pt()> isoCut03_ ) continue; }   // Depends on what kind of iso you are applying, this changes
              else {     if ( ptsum/muon1->pt()> isoCut03_ ) {continue;  }
            }


      bool haytrigger = false;
      const trigger::TriggerObjectCollection & toc(triggerObj->getObjects());
      for ( size_t ia = 0; ia < triggerObj->sizeFilters(); ++ ia) {
      if( triggerObj->filterTag(ia)  == hltFilterTag_ ){
//|| triggerObj->filterTag(ia)  == hltFilterTagMu11_ || triggerObj->filterTag(ia)  == hltFilterTagMu15_ ) {
      const trigger::Keys & k = triggerObj->filterKeys(ia);
       for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
       double etatrig=toc[*ki].eta();
       double phitrig=toc[*ki].phi();
       if(sqrt(pow(etatrig-muon1->eta(),2)+pow(phitrig-muon1->phi(),2)) < 0.2 ) {haytrigger=true; break;}
       }
      }
      }
      if (haytrigger != true) continue;//cout<<"tag  "<<tag<<endl;*/
         pt=eta=Isolation=-999; PassCut = charge = 0;
         for (std::vector<reco::MuonRef>::const_iterator m2=muones.begin(); m2!=muones.end(); ++m2){
          muon2= (*m2);
          if(muon2==muon1 || muon2->charge()*muon1->charge() !=-1.) continue;



          const math::XYZTLorentzVector Muon4_1 (muon1->px(), muon1->py() , muon1->pz(), muon1->p());
	  const math::XYZTLorentzVector Muon4_2 (muon2->px(), muon2->py() , muon2->pz(), muon2->p());
	  const math::XYZTLorentzVector myZCand=Muon4_1+Muon4_2;
	  double Zm_aux =  myZCand.mass();
          if ( Zm_aux<60. ||  Zm_aux>120.) continue;


        if (!muon2->isStandAloneMuon()) continue;
           double cal=muon2->isolationR03().emEt+muon2->isolationR03().hadEt;
            if (isTotIso_==true){
              Isolation = (cal+ptsum)/muon2->pt(); }
//              if ( (cal+ptsum)/muon2->pt()> isoCut03_ ) continue; }   // Depends on what kind of iso you are applying, this changes
              else {
              Isolation = ptsum/muon2->pt();
//if ( ptsum/muon2->pt()> isoCut03_ ) {continue;  }
            }
        if (Isolation < 0.15) {
        h1_["hPtMu_DENtr"]->Fill(muon2->pt());
        h1_["hPhiMu_DENtr"]->Fill(muon2->phi());
        h1_["hEtaMu_DENtr"]->Fill(muon2->eta());}
        pt = muon2->pt();
        eta = muon2->eta();
        charge = muon2->charge();

//referencia a la traza:
         for (unsigned int i=0; i != trackCollection->size(); i++) {
            TrackRef t1(trackCollection,i);

        reco::TrackRef trackerTrackRef1 = t1; //ref

        reco::TrackRef trackTR = muon2->innerTrack(); //ref
        if (trackTR == trackerTrackRef1){ 
//*******************
           unsigned int tracitas = 0;
           for (unsigned int ii=0; ii != trackCollection->size(); ii++) {
             TrackRef trackp(trackCollection,ii);
             if(trackp == trackTR) continue;
               if (trackp->pt() < 0.5) continue;
               Geom::Phi<double> deltaphi(trackTR->phi()-trackp->phi());

               if(sqrt(pow(deltaphi.value(),2)+pow(trackTR->eta()-trackp->eta(),2)) <0.3) tracitas++;
               }

            if (tracitas>5) continue;

//******************

           PassCut++;
           if (Isolation < 0.15) {
           h1_["hPtMu_NUMtr"]->Fill(muon2->pt());
           h1_["hPhiMu_NUMtr"]->Fill(muon2->phi());
           h1_["hEtaMu_NUMtr"]->Fill(muon2->eta());}
        continue;
        }
       if (PassCut >0 ) continue; 
//        Trackingtree->Fill();
        }
        Trackingtree->Fill();
        }
    }
//-----------------------------------

}



DEFINE_FWK_MODULE(ZMuMuTracking);

