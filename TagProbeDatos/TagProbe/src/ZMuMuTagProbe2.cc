/* class ZMuMuTagProbeAnalyzer *  ZMuMuTagProbe Analyzer modified t
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
//#include "TFitterMinuit.h"
//#include "Minuit2/FCNBase.h"
////#include "utils.h"
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

//#include "DQMServices/Core/interface/MonitorElement.h"

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
double MassMuon=0.105658369;

class ZMuMuTagProbe : public edm::EDAnalyzer {
public:
      typedef math::XYZPoint Point;
  ZMuMuTagProbe(const edm::ParameterSet& pset);
  virtual ~ZMuMuTagProbe();
	  virtual void beginJob();
	  bool GoodEWKMuon(reco::MuonRef Mu);
	     // PathInfoCollection hltPaths_;

	 const std::vector<reco::MuonRef> selectMuons(edm::Handle<reco::MuonCollection> muonCollection);
	 const std::vector<reco::MuonRef> selectMuons2(edm::Handle<reco::MuonCollection> muonCollection);
	 const std::vector<reco::TrackRef> selectTracks(edm::Handle<reco::TrackCollection> trackCollection);
	 bool compMomentum(reco::MuonRef m1, reco::MuonRef m2);
	 double computeMass(reco::TrackRef m1, reco::TrackRef m2);
         double REWEIGHTZ(double pt, double ZM);

	  bool FoundZ(edm::Handle<reco::MuonCollection>& muonCollection, edm::Handle<reco::TrackCollection>& trackCollection);
	  virtual void endJob();
	  virtual void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
	private:

  TTree * ProbeMuonVar;
  double massBoson, pt, eta, weightM, dxy, dz, dzError, dzDiff, NormChi2, Isolation;
  int Allqualcuts, TrigMatch, PassAllSelec, charge, isMuon, isGlobalMuon, isGlbTrkMuon, Numtrazas;
//  int MuonHits, PixelHits, TrackerHits, NumberMatches; 


	  edm::InputTag muonTag_;
	  edm::InputTag jetTag_;

	//------------------------
	  edm::InputTag triggerSummaryLabel_;
	  edm::InputTag hltFilterTag_;
          edm::InputTag hltFilterTagMu11_;
          edm::InputTag hltFilterTagMu15_;
       //   std::vector<std::string> hltFilterTag_;

	  bool useOnlyGlobalMuons_;
	  bool useMuonQualityCuts_;
	  double ptCut_;
	  double etaCut_;
	  bool isTotIso_;
	  double isoCut03_;
	  double massTMin_;
	  double massTMax_;
	  double acopCut_;
	  int WCharge_;
          double abc;
          double slope;
          double parZ;
          double abc2;
          double slope2;
          double parZ2;
          bool REWEIGHT;
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
         //TFitterMinuit *minimizer2;
	};

	using namespace std;
	using namespace edm;
	using namespace reco;
	using namespace trigger;
        void minuitFunction(int& nDim, double* gin, double& result, double par[], int flg);
        double myFunction(double x, double y, double z);
	/// Constructor
	ZMuMuTagProbe::ZMuMuTagProbe(const ParameterSet& pset) :
	      muonTag_(pset.getUntrackedParameter<edm::InputTag> ("MuonTag", edm::InputTag("muons"))),
	      useOnlyGlobalMuons_(pset.getUntrackedParameter<bool>("UseOnlyGlobalMuons", true)),
	      useMuonQualityCuts_(pset.getUntrackedParameter<bool>("UseMuonQualityCuts")),
	      ptCut_(pset.getUntrackedParameter<double>("PtCut", 10.)),
	      etaCut_(pset.getUntrackedParameter<double>("EtaCut", 2.1)),
	      isTotIso_(pset.getUntrackedParameter<bool>("IsTotIso", true)),
	      isoCut03_(pset.getUntrackedParameter<double>("IsoCut03", 0.15)),
	      acopCut_(pset.getUntrackedParameter<double>("AcopCut", 999999.)),
	      triggerSummaryLabel_(pset.getParameter<edm::InputTag>("triggerSummaryLabel")),
              hltFilterTag_(pset.getParameter<edm::InputTag>("hltL3FilterLabel")),
              hltFilterTagMu11_(pset.getParameter<edm::InputTag>("hltL3FilterLabelMu11")),
              hltFilterTagMu15_(pset.getParameter<edm::InputTag>("hltL3FilterLabelMu15")),
         //     hltFilterTag_(pset.getUntrackedParameter<std::vector <std::string> >("hltL3FilterLabel")),

              slope(pset.getUntrackedParameter<double>("SLOPE", -0.301)),
              parZ(pset.getUntrackedParameter<double>("ZPAR", 0.991)),
              abc(pset.getUntrackedParameter<double>("ABC",4.52)),
              slope2(pset.getUntrackedParameter<double>("SLOPE2", -0.301)),
              parZ2(pset.getUntrackedParameter<double>("ZPAR2", 0.991)),
              abc2(pset.getUntrackedParameter<double>("ABC2",4.52)),
              REWEIGHT(pset.getUntrackedParameter<bool>("doReweight",true))

	{
	  LogDebug("ZMuMuTagProbeAnalyzer")<<" ZMuMuTagProbeAnalyzer constructor called";
	}

	/// Destructor
	ZMuMuTagProbe::~ZMuMuTagProbe(){
	}

	void ZMuMuTagProbe::beginJob(){
	  // Create output files

	  edm::Service<TFileService> fs;
  //bool Allqualcuts, GlobalAndTracker;
  int MuonHits, PixelHits, TrackerHits, NumberMatches;

  ProbeMuonVar  = fs->make<TTree>("ProbeMuonVar","probe muon variables"); 
  ProbeMuonVar->Branch("dz",&dz,"dz/D");
  ProbeMuonVar->Branch("dzError",&dzError,"dzError/D");
  ProbeMuonVar->Branch("dzDiff",&dzDiff,"dzDiff/D");

  ProbeMuonVar->Branch("ptmuon",&pt,"pt/D");
  ProbeMuonVar->Branch("etamuon",&eta,"eta/D");
  ProbeMuonVar->Branch("charge",&charge,"charge/I");
  ProbeMuonVar->Branch("massBoson",&massBoson,"massBoson/D");
  ProbeMuonVar->Branch("weightM",&weightM,"weightM/D");
  ProbeMuonVar->Branch("Allqualcuts",&Allqualcuts,"Allqualcuts/I");
  ProbeMuonVar->Branch("Isolation",&Isolation,"Isolation/D");
  ProbeMuonVar->Branch("dxy",&dxy,"dxy/D");
  ProbeMuonVar->Branch("NormChi2",&NormChi2,"NormChi2/D");
  ProbeMuonVar->Branch("TrigMatch",&TrigMatch,"TrigMatch/I");
  ProbeMuonVar->Branch("PassAllSelec",&PassAllSelec,"PassAllSelec/I");
  ProbeMuonVar->Branch("isMuon",&isMuon,"isMuon/I");
  ProbeMuonVar->Branch("isGlobalMuon",&isGlobalMuon,"isGlobalMuon/I");
  ProbeMuonVar->Branch("isGlbTrkMuon",&isGlbTrkMuon,"isGlbTrkMuon/I");
  ProbeMuonVar->Branch("Numtrazas",&Numtrazas,"Numtrazas/I");

         ff = TFile::Open("Example.root", "RECREATE");
         ff->cd();
         h = new TH1D("h","",100,40.,140.);
         hMassajuste = new TH1D("hMassajuste","",100,40.,140.);
	 hZ_MC = new TH1D("hZ_MC","",100,40.,140.);


	  numberOfEvents = 0;
	  numberOfMuons = 0;
	  
	  h1_["hPtMu_DEN"]=fs->make<TH1D>("ptMu_DEN","Pt mu",300,0.,300.);
	  h1_["hPhiMu_DEN"]=fs->make<TH1D>("phiMu_DEN","Phi mu",100,-TMath::Pi(), TMath::Pi());
          double xbin[7]={-2.1,-1.2,-0.9,0.,0.9,1.2,2.1};
//	  h1_["hEtaMu_DEN"]=fs->make<TH1D>("etaMu_DEN","Eta mu",100,-2.5,2.5);
          h1_["hEtaMu_DEN"]=fs->make<TH1D>("etaMu_DEN","Eta mu",6,xbin);

	  h1_["hPtMu_NUM"]=fs->make<TH1D>("ptMu_NUM","Pt mu",300,0.,300.);
	  h1_["hPhiMu_NUM"]=fs->make<TH1D>("phiMu_NUM","Phi mu",100,-TMath::Pi(), TMath::Pi());
//	  h1_["hEtaMu_NUM"]=fs->make<TH1D>("etaMu_NUM","Eta mu",100,-2.5,2.5);
          h1_["hEtaMu_NUM"]=fs->make<TH1D>("etaMu_NUM","Eta mu",6,xbin);

	  h2_["hEtaPtMu_NUM"]=fs->make<TH2D>("EtaPtMu_NUM","Pt vs Eta mu Probe sel",100,-2.5,2.5,100,0.,100.);
          h2_["hEtaPtMu_DEN"]=fs->make<TH2D>("EtaPtMu_DEN","Pt vs Eta mu Probe gen",100,-2.5,2.5,100,0.,100.);


	  h1_["hNMu"]=fs->make<TH1D>("NMu","Nb. muons in the event",10,0.,10.);
	  h1_["hNMu_trigger"]=fs->make<TH1D>("NMu_trigger","Nb. muons after trigger",10,0.,10.);

	  h1_["hPtMu"]=fs->make<TH1D>("ptMu","Pt mu",300,0.,300.);
	  h1_["hPtMu_trigger"]=fs->make<TH1D>("ptMu_trigger","Pt mu",300,0.,300.);
	  h1_["hEtaMu"]=fs->make<TH1D>("etaMu","Eta mu",50,-2.5,2.5);


	  h1_["hMmumu_DEN"]=fs->make<TH1D>("Inv_mass_mumu_DEN","inv. mass mumu DEN",100, 40.,140.);
          h1_["hMmumu_DENless35"]=fs->make<TH1D>("Inv_mass_mumu_DENless35","inv. mass mumu DEN pt<35",100, 40.,140.);
	  h1_["hMmumu"]=fs->make<TH1D>("Inv_mass_mumu","inv. mass mumu",100, 40.,140.);
	  h2_["hAcopl"]=fs->make<TH2D>("Acopl","phi1 vs phi2",50,-TMath::Pi(), TMath::Pi(),50,-TMath::Pi(), TMath::Pi());
          h1_["acop"]=fs->make<TH1D>("acop","#Delta #phi",50,0, TMath::Pi());


 
           h1_["NumTr_DEN"]=fs->make<TH1D>("NumTr","NumTr",50,0.,50.);
           h1_["NumTr_muon"]=fs->make<TH1D>("NumTr_muon","NumTr",50,0.,50.);
           h1_["NumTr_NUM"]=fs->make<TH1D>("NumTr_NUM","NumTr",50,0.,50.);
 
           h1_["hPtMu_DENm"]=fs->make<TH1D>("ptMu_DENm","Pt mu",300,0.,300.);
           h1_["hPhiMu_DENm"]=fs->make<TH1D>("phiMu_DENm","Phi mu",100,-TMath::Pi(), TMath::Pi());
           h1_["hEtaMu_DENm"]=fs->make<TH1D>("etaMu_DENm","Eta mu",6,xbin);
 
            h1_["hEtaMu"]=fs->make<TH1D>("etaMu","Eta mu",6,xbin);
            h1_["hEtaMu_numreco"]=fs->make<TH1D>("etaMu_numreco","etaMu_numreco",6,xbin);
           h1_["probemuon_trackerhits"]=fs->make<TH1D>("probemuonhits","#hits",50,0,50);
           h1_["probe_trackerhits"]=fs->make<TH1D>("probehits","#hits",50,0,50); 
           h1_["DeltaR"]=fs->make<TH1D>("DeltaR","DeltaR",100,0,4);
 
           h1_["hMmumu_DENm"]=fs->make<TH1D>("Inv_mass_mumu_DENm","inv. mass mumu DEN",100, 40.,140.);
           h1_["hMmumu_DENless35m"]=fs->make<TH1D>("Inv_mass_mumu_DENless35m","inv. mass mumu DEN pt<35",100, 40.,140.);
 
           h1_["hMmumu_DENreco"]=fs->make<TH1D>("Inv_mass_mumu_DENreco","inv. mass mumu DEN reco",100, 40.,140.);
           h1_["hMmumu_DENless35reco"]=fs->make<TH1D>("Inv_mass_mumu_DENless35reco","inv. mass mumu DEN pt<35 reco",100, 40.,140.);
 
           h1_["hMmumu_DENiso"]=fs->make<TH1D>("Inv_mass_mumu_DENiso","inv. mass mumu DEN",100, 40.,140.);
           h1_["hMmumu_DENless35iso"]=fs->make<TH1D>("Inv_mass_mumu_DENless35iso","inv. mass mumu DEN pt<35",100, 40.,140.);
 
           h1_["hptMu_deniso"]=fs->make<TH1D>("ptdenisoeff","Pt mu den iso eff.",300,0.,300.);
           h1_["hptMu_dentrigger"]=fs->make<TH1D>("ptdentriggereff","Pt mu den trigger eff.",300,0.,300.);

           h1_["hptMu_denisoCKT"]=fs->make<TH1D>("ptdenisoeffcock","Pt mu den iso eff.",300,0.,300.);
           h1_["hptMu_dentriggerCKT"]=fs->make<TH1D>("ptdentriggereffcock","Pt mu den trigger eff.",300,0.,300.);

           h1_["hptMu_NUMCKT"]=fs->make<TH1D>("ptMu_NUM_cock","Pt mu COCKTAIL ",300,0.,300.);

           h1_["hPtMu_DENQual"]=fs->make<TH1D>("ptMu_DENQual","Pt mu",300,0.,300.);
           h1_["hPhiMu_DENQual"]=fs->make<TH1D>("phiMu_DENQual","Phi mu",100,-TMath::Pi(), TMath::Pi());
           h1_["hEtaMu_DENQual"]=fs->make<TH1D>("etaMu_DENQual","Eta mu",6,xbin);


	}

	void ZMuMuTagProbe::endJob(){
	  LogVerbatim("") << "ZMuMuTagProbeAnalyzer>>> FINAL PRINTOUTS -> BEGIN";
	  LogVerbatim("") << "ZMuMuTagProbeAnalyzer>>> Number of analyzed events= " << numberOfEvents;
	  LogVerbatim("") << "ZMuMuTagProbeAnalyzer>>> Number of analyzed muons= " << numberOfMuons;
	  LogVerbatim("") << "ZMuMuTagProbeAnalyzer>>> FINAL PRINTOUTS -> END";
	}

	bool ZMuMuTagProbe::GoodEWKMuon(reco::MuonRef Mu){
	     bool goodMuon=true;
	     
	     int validmuonhits=Mu->globalTrack()->hitPattern().numberOfValidMuonHits();
             int pixelHits = Mu->globalTrack()->hitPattern().numberOfValidPixelHits(); 
	     double numberOfValidHits=Mu->globalTrack()->hitPattern().numberOfValidTrackerHits();	 
            // d0, chi2, nhits quality cuts
            int    nMatches = Mu->numberOfMatches();//Mu->numberOfMatches();

            goodMuon = (numberOfValidHits>=11)&&(validmuonhits>0)&&(nMatches>1)&&(pixelHits>0);  
//	    goodMuon = (chi2<10.)&&(numberOfValidHits>=11)&&(validmuonhits>0)&&(nMatches>1)&&(pixelHits>0);
	    if (!Mu->isGlobalMuon()||!Mu->isTrackerMuon()) goodMuon=false;    
	    return goodMuon;
	}
	//subrutina para seleccionar buenos muones = pt, eta, iso cut, quality+global   
	const std::vector<reco::MuonRef> ZMuMuTagProbe::selectMuons(Handle<reco::MuonCollection> muonCollection){
	  std::vector<reco::MuonRef> muones;

	  for (unsigned int i=0; i<muonCollection->size(); i++) {
	    reco::MuonRef mu1(muonCollection,i);

	    //CORTE CALIDAD
	    //if(!mu1->isStandAloneMuon()) continue;
	    if ( (useOnlyGlobalMuons_ && !mu1->isGlobalMuon()) ) continue;// {cout << "no es global" << endl; continue;}
	    if ( (useMuonQualityCuts_ && !GoodEWKMuon(mu1)) ) continue; //{cout<<"no calidad"<<endl; continue;}

	    //CORTE CINEMATICO
	// pt
	      double pt = mu1->innerTrack()->pt();
	      LogTrace("") << "\t... pt= " << pt << " GeV";

	      if (pt<ptCut_) continue;
	// eta
	      double eta = mu1->innerTrack()->eta();
	      LogTrace("") << "\t... eta= " << eta;
	      if (fabs(eta)>etaCut_) continue;

	    muones.push_back(mu1);
	  }//for i

	  return muones;
	}

	//subrutina para seleccionar muones-pt, eta ,(iso) como preseleccion para la eficiencia de tracker
	const vector<reco::MuonRef> ZMuMuTagProbe::selectMuons2(Handle<reco::MuonCollection> muonCollection){
	  vector<reco::MuonRef> muones;
	 
	  for (unsigned int i=0; i<muonCollection->size(); i++) {
	    reco::MuonRef mu1(muonCollection,i);

	      if(!mu1->isStandAloneMuon()) continue;
	    //CORTE CINEMATICO
	// pt
	      double pt = mu1->pt();
	      LogTrace("") << "\t... pt= " << pt << " GeV";

	      if (pt<ptCut_) continue;
	// eta
	      double eta = mu1->eta();
	      LogTrace("") << "\t... eta= " << eta;
	      if (fabs(eta)>etaCut_) continue;

	//CORTE AISLAMIENTO
	    double ptsum = mu1->isolationR03().sumPt;
	    double cal=mu1->isolationR03().emEt+mu1->isolationR03().hadEt;
	    if (isTotIso_==true){
	      if ( (cal+ptsum)/mu1->pt()> isoCut03_ ) continue; }   // Depends on what kind of iso you are applying, this changes
	      else {     if ( ptsum/mu1->pt()> isoCut03_ ) {continue;  }
	    }
	    
	    muones.push_back(mu1);    
	  }//for i

	  return muones;
	}

	//Subrutina para seleccionar trazas con pt>25 GeV, abs(eta)<2., tales que 70< M_{mumu} <110
	const std::vector<reco::TrackRef> ZMuMuTagProbe::selectTracks(Handle<reco::TrackCollection> trackCollection){
	  std::vector<reco::TrackRef> trazas;
	  double Zm_aux = 0;

	  for (unsigned int i=0; i != trackCollection->size(); i++) {
	    TrackRef mu1(trackCollection,i);
	    //CORTE CINEMATICO
	// pt
	      double pt = mu1->pt();
	      LogTrace("") << "\t... pt= " << pt << " GeV";

	      if (pt<ptCut_) continue;
	// eta
	      double eta = mu1->eta();
	      LogTrace("") << "\t... eta= " << eta;
	      if (fabs(eta)>etaCut_) continue;
	       int numSel=0;
	       for (unsigned int j=0; j != trackCollection->size(); j++) {
	       TrackRef mu2(trackCollection,j);

	       if ( mu2->pt()<ptCut_) continue;
	       if (fabs(mu2->eta())>etaCut_) continue;
	       if( mu1->charge() * mu2->charge() != -1 || mu1==mu2) continue;
		const math::XYZTLorentzVector Muon4_1 (mu1->px(), mu1->py() , mu1->pz(), mu1->p());
		const math::XYZTLorentzVector Muon4_2 (mu2->px(), mu2->py() , mu2->pz(), mu2->p());
		const math::XYZTLorentzVector myZCand=Muon4_1+Muon4_2;
		Zm_aux=myZCand.mass();
		if ( Zm_aux<60. ||  Zm_aux>120.) continue;
		++numSel;
		if (numSel>1) {cout<<"repetida traza"<<endl; }       
		trazas.push_back(mu1);
		break;

	       }//for j
	  }//for i

	  return trazas;
	}

	double ZMuMuTagProbe::computeMass(TrackRef m1, TrackRef m2) {
	  const math::XYZTLorentzVector Muon4_1 (m1->px(), m1->py() , m1->pz(), m1->p());
	  const math::XYZTLorentzVector Muon4_2 (m2->px(), m2->py() , m2->pz(), m2->p());
	  const math::XYZTLorentzVector myZCand=Muon4_1+Muon4_2;
	  return myZCand.mass();
	}
	bool compMomentum(reco::MuonRef m1, reco::MuonRef m2){
	  return ( m1->pt() > m2->pt() );
	}



double MASSZBIN[61]={60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120};

// sin pile - up muon probe charge positive

double MASSZless35[60]={738., 692., 822., 816., 803., 820., 901., 870., 918., 1032., 1042., 1049., 1223., 1220., 1275., 1364., 1502., 1522., 1783., 1809., 2066., 2238., 2483., 3017., 3603., 4535., 6000., 8499., 12908., 18949., 23517., 22637., 16554., 10022., 5807., 3240., 2130., 1397., 1043., 640., 547., 484., 343., 288., 275., 220., 183., 158., 126., 104., 118., 91., 93., 85., 95., 73., 55., 50., 53., 49.};

double MASSZ[60]={67., 89., 131., 153., 185., 242., 273., 319., 427., 480., 569., 691., 786., 946., 1095., 1303., 1442., 1749., 2188., 2487., 3111., 3797., 4943., 6156., 8326., 11223., 16616., 25139., 40518., 61343., 77888., 79394., 63737., 42582., 25897., 15708., 10087., 7109., 5041., 3876., 3061., 2472., 2101., 1737., 1423., 1320., 1200., 1046., 943., 845., 753., 698., 604., 614., 527., 504., 431., 408., 400., 359.};


double ZMuMuTagProbe::REWEIGHTZ(double pt, double Zmass){

double weight = 1.0;
int BinZmass=-5;
const int nint = 60;
const int siz = 1;

  for (int i = 0; i < nint; i++ ){
     if(Zmass >= 60. + siz*i && Zmass < 60 + siz*(i+1) ) {BinZmass=i; break;}
  }
  if (pt>=35.) weight = parZ*MASSZ[BinZmass]*(36.1/1238.1)/(parZ*MASSZ[BinZmass]*36.1/1238.1+(abc+slope*(siz*BinZmass+60+siz*0.5)));
  else weight = parZ2*MASSZless35[BinZmass]*(36.1/1238.1)/(parZ2*MASSZless35[BinZmass]*36.1/1238.1+(abc2+slope2*(siz*BinZmass+60+siz*0.5)));

return weight;

}



	void ZMuMuTagProbe::analyze(const Event & event, const EventSetup& eventSetup){

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

	MuonRef SelectedMuon; 

	//variables
	std::vector<reco::TrackRef> trazas;
	double Zm_aux=0;
	TrackRef track1;

	std::vector<reco::MuonRef> muones;
	MuonRef muon1, muon2;
	trazas = selectTracks(trackCollection);
	muones = selectMuons(muonCollection);
	//1er loop sobre las trazas-->tag?
	for(std::vector<reco::TrackRef>::const_iterator t1=trazas.begin(); t1!=trazas.end(); ++t1){
	track1= (*t1);
   

	int tag=0; //if(tag=0) 
	//loop sobre muones--> MATCH CON LA REFERENCIA!!! 
	reco::TrackRef trackerTrackRef1 = track1; //ref
        if(trazas.size() >2 )  break;
	bool haytrigger = false;
        double dz1 = -1; double dzError1 = -1 ;

	for (std::vector<reco::MuonRef>::const_iterator m1=muones.begin(); m1!=muones.end(); ++m1){
	muon1= (*m1);


        double d0 = muon1->globalTrack()->dxy(beamSpotHandle->position()); if(fabs(d0)>=0.2) continue;
        double chi2=muon1->globalTrack()->normalizedChi2(); if(chi2>=10) continue;
	reco::TrackRef trackTR = muon1->innerTrack(); //ref
	if (trackTR != trackerTrackRef1) continue; //ref
        dz1 = muon1->globalTrack()->dz(beamSpotHandle->position());
        dzError1 = muon1->globalTrack()->dzError();

       cout << muon1->pt() << endl;
             double ptsum = muon1->isolationR03().sumPt;
             double cal=muon1->isolationR03().emEt+muon1->isolationR03().hadEt;
             if (isTotIso_==true){
               if ( (cal+ptsum)/muon1->pt()> isoCut03_ ) continue; }   // Depends on what kind of iso you are applying, this changes
               else {     if ( ptsum/muon1->pt()> isoCut03_ ) {continue;  }
             }


	const trigger::TriggerObjectCollection & toc(triggerObj->getObjects());
	for ( size_t ia = 0; ia < triggerObj->sizeFilters(); ++ ia) {
	  if( triggerObj->filterTag(ia)  == hltFilterTag_ ) {
          const trigger::Keys & k = triggerObj->filterKeys(ia);
	    for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
	    double etatrig=toc[*ki].eta();
	    double phitrig=toc[*ki].phi();
	    Geom::Phi<double> deltaphi(phitrig-muon1->phi());
            if(sqrt(pow(etatrig-muon1->eta(),2)+pow(deltaphi.value(),2)) < 0.2 ) {haytrigger=true; ++tag;}
	    }
	  }
        }
	}//for m1--------------

	if(tag>1) cout<<"HAY DOS TAGS"<<endl;
	if(tag < 1) continue;
	int probe=0;

	//2. loop sobre las trazas, si tenemos tag, cogemos las variables del  probe (de la traza!)
        TrackRef track2, ttrack2, trackp;
        Zm_aux=0.;
         int tracitas = 0;
          double deltamass=100.;     
          for (std::vector<reco::TrackRef>::const_iterator tt2=trazas.begin(); tt2!=trazas.end(); ++tt2){
          ttrack2= (*tt2);
          if(track1->charge() * ttrack2->charge() != -1) continue;
          double Zm_aux2 = computeMass(track1, ttrack2);
          if ( Zm_aux2<60.0||  Zm_aux2>120.) continue;
           if (abs(MZ-Zm_aux2) < deltamass){
          track2=ttrack2;
          deltamass = abs(MZ-Zm_aux2);
          Zm_aux = Zm_aux2;
          
         }

           for (unsigned int i=0; i != trackCollection->size(); i++) {
             TrackRef trackp(trackCollection,i);
             if(trackp == track2) continue;
              if (trackp->pt() < 0.5) continue;
//               if (trackp->pt() < 1.0) continue;
               Geom::Phi<double> deltaphi(track2->phi()-trackp->phi());

               h1_["DeltaR"]->Fill(sqrt(pow(deltaphi.value(),2)+pow(track2->eta()-trackp->eta(),2)));
               if(sqrt(pow(deltaphi.value(),2)+pow(track2->eta()-trackp->eta(),2)) <0.3) tracitas++;
               }
               h1_["probe_trackerhits"]->Fill(track2->hitPattern().numberOfValidTrackerHits());

          }

        if(!track2)  continue;
        bool PassCut = true;
        PassAllSelec = 0; TrigMatch= 0; Allqualcuts = 0; charge = 0;
        pt = eta = massBoson = weightM = NormChi2 = dxy = Isolation =  dz = dzError = dzDiff = -999.;

        Numtrazas = tracitas;
                 h1_["NumTr_DEN"]->Fill(tracitas); 
          //       if (tracitas>5) continue;

                 double weightetapt =  1.;//Wetapt;         
                 double dz2 =  track2->dz(beamSpotHandle->position());
                 double dzError2 = track2->dzError();
                 dz = dz2;
                 dzError = dzError2;
                 dzDiff = abs(dz2-dz1)/sqrt(pow(dzError2,2)+pow(dzError1,2));


                if (REWEIGHT==true) weightetapt = REWEIGHTZ(track2->pt(),Zm_aux); //cout << " weightetapt " <<  weightetapt <<endl;
                h1_["hPtMu_DEN"]->Fill(track2->pt(),weightetapt);

		h1_["hPhiMu_DEN"]->Fill(track2->phi(),weightetapt);//relleno con traza no con muon!!!
		h1_["hEtaMu_DEN"]->Fill(track2->eta(),weightetapt);
                if (track2->pt()<35) h1_["hMmumu_DENless35"]->Fill(Zm_aux,weightetapt); 
                else h1_["hMmumu_DEN"]->Fill(Zm_aux,weightetapt);
                
                h2_["hEtaPtMu_DEN"]->Fill(track2->eta(),track2->pt());
		h->Fill(Zm_aux);
                pt = track2->pt();
                eta = track2->eta(); 
                charge = track2->charge();
                massBoson = Zm_aux;
                weightM = weightetapt;

           weightetapt = 1.0; isMuon = 0; isGlobalMuon = 0; isGlbTrkMuon = 0; 
           for (unsigned int i=0; i<muonCollection->size(); i++) {
             MuonRef mu1(muonCollection,i);
            

             double ptsum2= mu1->isolationR03().sumPt;
             double cal2= mu1->isolationR03().emEt+mu1->isolationR03().hadEt;
             Isolation = (cal2+ptsum2)/mu1->pt();



             reco::TrackRef trackTR2 = mu1->innerTrack();
             reco::TrackRef trackerTrackRef2 = track2;
             if(trackTR2 != trackerTrackRef2) continue;
               isMuon = 1; 
             if(!mu1->isGlobalMuon()) continue ; isGlobalMuon = 1;
             if(!mu1->isGlobalMuon() || !mu1->isTrackerMuon()) continue; 
             isGlbTrkMuon = 1;   
                 h1_["NumTr_muon"]->Fill(tracitas);
                 h1_["hPtMu_DENm"]->Fill(track2->pt(),weightetapt);

               h1_["hPhiMu_DENm"]->Fill(track2->phi(),weightetapt);//relleno con traza no con muon!!!
               h1_["hEtaMu_DENm"]->Fill(track2->eta(),weightetapt);
                 if (track2->pt()<35) h1_["hMmumu_DENless35m"]->Fill(Zm_aux,weightetapt);
                 else h1_["hMmumu_DENm"]->Fill(Zm_aux,weightetapt);
                 h1_["probemuon_trackerhits"]->Fill(track2->hitPattern().numberOfValidTrackerHits());

           }

            

 //               continue;
		 //2. loop sobre los muones-> Ver si el "probe" pasa la seleccion, de nuevo falta referencia traza-muon!!

		 reco::TrackRef trackerTrackRef2 = track2;
		 for (vector<MuonRef>::const_iterator m2=muones.begin(); m2!=muones.end(); ++m2){
		 muon2= (*m2);
                 reco::TrackRef trackTR2 = muon2->innerTrack(); //ref
                 if (trackTR2 != trackerTrackRef2) continue; //PassCut =false; //refa
                 if(PassCut==true){
                 Allqualcuts = 1; 
                 h1_["hPtMu_DENQual"]->Fill(track2->pt(),weightetapt);
                 h1_["hPhiMu_DENQual"]->Fill(track2->phi(),weightetapt);//relleno con traza no con muon!!!
                 h1_["hEtaMu_DENQual"]->Fill(track2->eta(),weightetapt);
                 }
                 


                 double chi22=muon2->globalTrack()->normalizedChi2(); if(chi22>=10) PassCut =false;
                 double d02 = muon2->globalTrack()->dxy(beamSpotHandle->position()); if(fabs(d02)>=0.2) PassCut =false;
                 dxy = d02;
                 NormChi2 = chi22;
/*
           reco::TrackRef cocktail =  muon::tevOptimized(muon2->combinedMuon(), muon2->track(),
           *tevMap_default, *tevMap_1stHit, *tevMap_picky);
 */         
                 if(PassCut==true){
                 h1_["hEtaMu_numreco"]->Fill(track2->eta());
                 if (track2->pt()<35) h1_["hMmumu_DENless35reco"]->Fill(Zm_aux,weightetapt);
                 else h1_["hMmumu_DENreco"]->Fill(Zm_aux,weightetapt);

                 h1_["hptMu_deniso"]->Fill(track2->pt());
   //              h1_["hptMu_denisoCKT"]->Fill((*cocktail).pt());
                 }

             double ptsum2= muon2->isolationR03().sumPt;
             double cal2=muon2->isolationR03().emEt+muon2->isolationR03().hadEt;
             if (isTotIso_==true){
               Isolation = (cal2+ptsum2)/muon2->pt();  
               if ( (cal2+ptsum2)/muon2->pt()> isoCut03_ )  PassCut =false; }   // Depends on what kind of iso you are applying, this changes
             else {
               Isolation = ptsum2/muon2->pt();
               if ( ptsum2/muon2->pt()> isoCut03_ ) PassCut =false;  }
                          
         if(PassCut==true){ 
         h1_["hNMu"]->Fill(track2->eta());
         h1_["hEtaMu"]->Fill(track2->eta());
         h1_["hptMu_dentrigger"]->Fill(track2->pt());
//         h1_["hptMu_dentriggerCKT"]->Fill((*cocktail).pt());
         if (track2->pt()<35) h1_["hMmumu_DENless35iso"]->Fill(Zm_aux,weightetapt);
         else h1_["hMmumu_DENiso"]->Fill(Zm_aux,weightetapt);
         } 

		 bool haytrigger2 = false;
		 int probe=0;
		 const trigger::TriggerObjectCollection & toc2(triggerObj->getObjects());
		 for ( size_t ia = 0; ia < triggerObj->sizeFilters(); ++ ia) {
		  if( triggerObj->filterTag(ia)  == hltFilterTag_ ) {
// || triggerObj->filterTag(ia)  == hltFilterTagMu11_ || triggerObj->filterTag(ia)  == hltFilterTagMu15_ ) {
		  const trigger::Keys & k = triggerObj->filterKeys(ia);
		   for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
		   double etatrig=toc2[*ki].eta();
		   double phitrig=toc2[*ki].phi();
		   Geom::Phi<double> deltaphi(phitrig-muon2->phi()); 
		   if (sqrt(pow(etatrig-muon2->eta(),2)+pow(deltaphi.value(),2)) < 0.20) {haytrigger2=true; break;}
		    }
		   }
		  }//for size_t ia...
                  TrigMatch = haytrigger2; 
		  if (haytrigger2 == false) PassCut =false;
                  
		  ++probe; 
 
               weightetapt = 1.;	
                 if(PassCut == true){
                 
                 PassAllSelec = 1;                

//                 h1_["hptMu_NUMCKT"]->Fill((*cocktail).pt());
                 h2_["hEtaPtMu_NUM"]->Fill(muon2->innerTrack()->eta(),muon2->innerTrack()->pt());
		 h1_["hPtMu_NUM"]->Fill(muon2->innerTrack()->pt(),weightetapt);
		 h1_["hPhiMu_NUM"]->Fill(muon2->innerTrack()->phi(),weightetapt);
		 h1_["hEtaMu_NUM"]->Fill(muon2->innerTrack()->eta(),weightetapt);

		 h1_["hMmumu"]->Fill(Zm_aux);

                 h1_["NumTr_NUM"]->Fill(tracitas);
                 }
		 }//for mu2
   ProbeMuonVar->Fill();
		if(probe>1)  cout <<"HAY DOS PROBES"<<endl;

  }  //for track i---


}



DEFINE_FWK_MODULE(ZMuMuTagProbe);

