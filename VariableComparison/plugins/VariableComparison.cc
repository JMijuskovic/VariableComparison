#include "../interface/VariableComparison.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

/* Collaborating Class Header */
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"


#include "TLorentzVector.h"
#include "TFile.h"
#include <vector>
#include <cmath>
#include "TROOT.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TBranch.h"

/* C++ Headers */
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
using namespace edm;
using namespace pat;
using namespace reco;

VariableComparison::VariableComparison(const::ParameterSet& pSet)
{

  //initialise parameters:
  parameters= pSet;
  
  //declare consumes:
  theMuonCollectionLabel_ = consumes<edm::View<pat::Muon> >  (parameters.getParameter<edm::InputTag>("MuonCollection"));
  theVertexLabel_          = consumes<edm::View<reco::Vertex> >         (parameters.getParameter<edm::InputTag>("VertexLabel"));
  theBeamSpotLabel_        = consumes<reco::BeamSpot>       (parameters.getParameter<edm::InputTag>("BeamSpotLabel"));
  trigger_          = consumes<edm::View<pat::TriggerObjectStandAlone> >  (parameters.getParameter<edm::InputTag>("triggerObjects"));
  trgresultsToken_ = consumes<edm::TriggerResults> (parameters.getParameter<edm::InputTag>("triggerResults"));

  /*sample parameters*/
  noTrigger = parameters.getParameter<bool>("noTrigger");

  /*tag and probe selection*/
  pair_minInvMass = parameters.getParameter<double>("pair_minInvMass");
  pair_maxInvMass = parameters.getParameter<double>("pair_maxInvMass");     
  tag_minPt = parameters.getParameter<double>("tag_minPt");      
  tag_muonID = parameters.getParameter<string>("tag_muonID");
  tag_isoCut = parameters.getParameter<double>("tag_isoCut");
  tag_hltDrCut = parameters.getParameter<double>("tag_hltDrCut");
  tag_hltFilter = parameters.getParameter<string>("tag_hltFilter");
  probe_minPt = parameters.getParameter<double>("probe_minPt");
  probe_trackerIsoCut = parameters.getParameter<double>("probe_trackerIsoCut");      
  probe_muonIDs = parameters.getParameter<string>("probe_muonIDs");
  hlt_path = parameters.getParameter<string>("hlt_path"); 
  etaMin= parameters.getParameter<vector<double>> ("etaMin");
  etaMax=parameters.getParameter<vector<double>> ("etaMax");
}


VariableComparison::~VariableComparison(){
}

void VariableComparison::beginJob(){
  
  edm::Service<TFileService> fileService;


  TFileDirectory time = fileService->mkdir( "Time" );

  UnbSTAmuonTime = time.make<TH1F> ("UnbSTAmuonTime", "UnbSTAmuonTime;time (ns);# entries",400.,-200.,200.);
  UnbSTAmuonTimeBarrel = time.make<TH1F> ("UnbSTAmuonTimeBarrel", "UnbSTAmuonTimeBarrel;time(ns); # entries",400.,-200.,200.);
  UnbSTAmuonTimeEndcap = time.make<TH1F> ("UnbSTAmuonTimeEndcap", "UnbSTAmuonTimeEndcap;time(ns); # entries",400.,-200.,200.);
  STAmuonTime= time.make<TH1F> ("STAmuonTime", "StaMuonTime;time (ns);# entries", 400, -200., 200.);
  STAmuonTimeBarrel = time.make<TH1F>("STAmuonTimeBarrel","StaMuonTimeBarrel;time (ns);# entries", 400, -200., 200.);
  STAmuonTimeEndcap = time.make<TH1F>("STAmuonTimeEndcap", "StaMuonTimeEndcap;time (ns);# entries", 400, -200., 200.);


  TFileDirectory control = fileService->mkdir("Control");

  invMass = control.make<TH1F> ("invMass", "Invariant Mass; mass (GeV); # entries",100.,0.,200.);
  invMassInRange = control.make<TH1F> ("invMassRange", "Invariant Mass in range; mass(GeV); # entries",100.,0.,200.);
  dilepPt = control.make<TH1F> ("dilepPt", "Di lepton pt; p_{T}(GeV); # entries",100.,0.,200.);
  nVertices = control.make<TH1F> ("nVertices", "nVertices; # Vertices; # entries",180.,0.,180.);

  TFileDirectory kinematical = fileService->mkdir("Kinematical variables");
  
  for (unsigned int iEtaRegion=0; iEtaRegion<6; iEtaRegion++){
    char histname [3];
    double *etaCutMin = etaMin.data();
    double *etaCutMax = etaMax.data();

    char pt[]="ProbePt";
    sprintf(histname, "%s minEta %f maxEta %f",pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbePt.push_back(kinematical.make<TH1F> (histname, "ProbePt; p_{T}(GeV); # entries",75.,0.,150.));
    char pt_pt[]="ProbePt_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f",pt_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbePt_vsPt.push_back(kinematical.make<TProfile>(histname,"ProbePt_vsPt; p_{T} (GeV); p_{T} (GeV)",50.,0.,150.,0.,150.));
    char pt_eta[]="ProbePt_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f",pt_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbePt_vsEta.push_back(kinematical.make<TProfile>(histname,"ProbePt_vsEta; #eta; p_{T} (GeV)",48.,-2.4,2.4,0.,150.));
    char pt_phi[]="ProbePt_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f",pt_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbePt_vsPhi.push_back(kinematical.make<TProfile>(histname,"ProbePt_vsPhi; #phi; p_{T} (GeV)",48.,-TMath::Pi(),TMath::Pi(),0.,150.));

    char eta[]="ProbeEta";
    sprintf(histname, "%s minEta %f maxEta %f",eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbeEta.push_back(kinematical.make<TH1F> (histname, "ProbeEta; #eta; # entries",48.,-2.4,2.4));
    char eta_pt[]="ProbeEta_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f",eta_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbeEta_pt.push_back(kinematical.make<TProfile> (histname, "ProbeEta_vsPt; p_{T};#eta",50.,0.,150.,-2.4,2.4));
    char eta_eta[]="ProbeEta_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f",eta_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbeEta_eta.push_back(kinematical.make<TProfile> (histname, "ProbeEta_vsEta; #eta; #eta ",48.,-2.4,2.4,-2.4,2.4));
    char eta_phi[]="ProbeEta_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f",eta_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbeEta_phi.push_back(kinematical.make<TProfile> (histname, "ProbeEta_vsPhi; #phi; #eta",48.,-TMath::Pi(),TMath::Pi(),-2.4,2.4));

    char phi[]="ProbePhi";
    sprintf(histname, "%s minEta %f maxEta %f",phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbePhi.push_back(kinematical.make<TH1F> (histname, "ProbePhi; #phi; # entries",48.,-TMath::Pi(),TMath::Pi()));
    char phi_pt[]="ProbePhi_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f",phi_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbePhi_pt.push_back(kinematical.make<TProfile> (histname, "ProbePhi_vsPt; p_{T} (GeV); #phi",50.,0.,150.,-TMath::Pi(),TMath::Pi()));
    char phi_eta[]="ProbePhi_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f",phi_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbePhi_eta.push_back(kinematical.make<TProfile> (histname, "ProbePhi_vsEta; #eta; #phi",48.,-2.4,2.4,-TMath::Pi(),TMath::Pi()));
    char phi_phi[]="ProbePhi_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f",phi_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ProbePhi_phi.push_back(kinematical.make<TProfile> (histname, "ProbePhi_vsPhi; #phi; #phi ",48.,-TMath::Pi(),TMath::Pi(),-TMath::Pi(),TMath::Pi()));
  }

  TFileDirectory ID = fileService->mkdir ("ID");
  for (unsigned int iEtaRegion=0; iEtaRegion<6; iEtaRegion++){
    char histname [3];
    double *etaCutMin = etaMin.data();
    double *etaCutMax = etaMax.data();

    char nHitsGLB[]="NHitsGLB";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsGLB,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsGLB.push_back(ID.make<TH1F> (histname, "NHitsGLB; # hits; # entries", 80, 0., 80.));
    char nHitsGLB_pt[]="NHitsGLB_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsGLB_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsGLB_pt.push_back(ID.make<TProfile> (histname, "NHitsGLB; p_{T}; # hits",50.,0.,150., 0., 80.));
    char nHitsGLB_eta[]="NHitsGLB_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsGLB_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsGLB_eta.push_back(ID.make<TProfile> (histname, "NHitsGLB; #eta; # hits",48.,-2.4,2.4, 0., 80.));
    char nHitsGLB_phi[]="NHitsGLB_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsGLB_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsGLB_phi.push_back(ID.make<TProfile> (histname, "NHitsGLB; #phi; # hits",48.,-TMath::Pi(),TMath::Pi(), 0., 80.));

    char nHitsTRK[]="NHitsTRK";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsTRK,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsTRK.push_back(ID.make<TH1F> (histname, "NHitsTRK; # hits; # entries", 40, 0., 40.));
    char nHitsTRK_pt[]="NHitsTRK_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsTRK_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsTRK_pt.push_back(ID.make<TProfile> (histname, "NHitsTRK_vsPt; p_{T} (GeV); # hits", 50.,0.,150., 0., 40.));
    char nHitsTRK_eta[]="NHitsTRK_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsTRK_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsTRK_eta.push_back(ID.make<TProfile> (histname, "NHitsTRK_eta; #eta; # hits", 48.,-2.4,2.4, 0., 40.));
    char nHitsTRK_phi[]="NHitsTRK_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsTRK_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsTRK_phi.push_back(ID.make<TProfile> (histname, "NHitsTRK_pt; #phi ; # hits", 48.,-TMath::Pi(),TMath::Pi(), 0., 40.));

    char nHitsSTA[]="NHitsSTA";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsSTA,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsSTA.push_back(ID.make<TH1F> (histname, "NHitsSTA; # hits; # entries", 60, 0., 60.));
    char nHitsSTA_pt[]="NHitsSTA_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsSTA_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsSTA_pt.push_back(ID.make<TProfile> (histname, "NHitsSTA_vsPt; p_{T} (GeV); # hits", 50.,0.,150., 0., 60.));
    char nHitsSTA_eta[]="NHitsSTA_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsSTA_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsSTA_eta.push_back(ID.make<TProfile> (histname, "NHitsSTA_vsEta; #eta; # hits", 48.,-2.4,2.4, 0., 60.));
    char nHitsSTA_phi[]="NHitsSTA_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f",nHitsSTA_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NHitsSTA_phi.push_back(ID.make<TProfile> (histname, "NHitsSTA_vsPhi; #phi ; # hits", 48.,-TMath::Pi(),TMath::Pi(), 0., 60.));

    char chi2GLB[]=" Chi2GLB";
    sprintf(histname, "%s minEta %f maxEta %f", chi2GLB,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Chi2GLB.push_back(ID.make<TH1F>(histname, "Chi2GLB; chi2/ndof; # entries", 50, 0., 100.));
    char chi2GLB_pt[]=" Chi2GLB_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", chi2GLB_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Chi2GLB_pt.push_back(ID.make<TProfile>(histname, "Chi2GLB_vsPt; p_{T} (GeV); chi2/ndof", 50.,0.,150., 0., 100.));
    char chi2GLB_eta[]=" Chi2GLB_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", chi2GLB_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Chi2GLB_eta.push_back(ID.make<TProfile>(histname, "Chi2GLB_vsEta; #eta; chi2/ndof", 48.,-2.4,2.4, 0., 100.));
    char chi2GLB_phi[]=" Chi2GLB_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", chi2GLB_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Chi2GLB_phi.push_back(ID.make<TProfile>(histname, "Chi2GLB_vsPhi; #phi; chi2/ndof",  48.,-TMath::Pi(),TMath::Pi(), 0., 100.));

    char chi2TRK[]=" Chi2TRK";
    sprintf(histname, "%s minEta %f maxEta %f", chi2TRK,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Chi2TRK.push_back(ID.make<TH1F>(histname, "Chi2TRK; chi2/ndof; # entries", 100, 0., 50.));
    char chi2TRK_pt[]=" Chi2TRK_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", chi2TRK_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Chi2TRK_pt.push_back(ID.make<TProfile>(histname, "Chi2TRK_vsPt; p_{T} (GeV); chi2/ndof", 50.,0.,150., 0., 100.));
    char chi2TRK_eta[]=" Chi2TRK_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", chi2TRK_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Chi2TRK_eta.push_back(ID.make<TProfile>(histname, "Chi2TRK_vsEta; #eta; chi2/ndof", 48.,-2.4,2.4, 0., 100.));
    char chi2TRK_phi[]=" Chi2TRK_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", chi2TRK_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Chi2TRK_phi.push_back(ID.make<TProfile>(histname, "Chi2TRK_vsPhi; #phi; chi2/ndof",  48.,-TMath::Pi(),TMath::Pi(), 0., 100.));

    char nMatchedStation[]="NMatchedStation";
    sprintf(histname, "%s minEta %f maxEta %f", nMatchedStation,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NMatchedStation.push_back(ID.make<TH1F>(histname,"NMatchedStation; # stations; # entries", 10, 0., 10.));
    char nMatchedStation_pt[]="NMatchedStation_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", nMatchedStation_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NMatchedStation_pt.push_back(ID.make<TProfile>(histname,"NMatchedStation_vsPt; p_{T} (GeV); # stations", 50.,0.,150., 0., 10.));
    char nMatchedStation_eta[]="NMatchedStation_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", nMatchedStation_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NMatchedStation_eta.push_back(ID.make<TProfile>(histname,"NMatchedStation_vsEta; #eta ; # stations", 48.,-2.4,2.4, 0., 10.));
    char nMatchedStation_phi[]="NMatchedStation_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", nMatchedStation_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NMatchedStation_phi.push_back(ID.make<TProfile>(histname,"NMatchedStation_vsPhi; #phi; # stations", 48.,-TMath::Pi(),TMath::Pi(), 0., 10.));

    char nMuonValidHitsGLB[]="NMuonValidHitsGLB";
    sprintf(histname, "%s minEta %f maxEta %f", nMuonValidHitsGLB,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NMuonValidHitsGLB.push_back(ID.make<TH1F>(histname,"NMuonValidHitsGLB; # hits; # entries", 60, 0., 60.));
    char nMuonValidHitsGLB_pt[]="NMuonValidHitsGLB_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", nMuonValidHitsGLB_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NMuonValidHitsGLB_pt.push_back(ID.make<TProfile>(histname,"NMuonValidHitsGLB_vsPt; p_{T} (GeV); # hits", 50.,0.,150., 0., 60.));
    char nMuonValidHitsGLB_eta[]="NMuonValidHitsGLB_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", nMuonValidHitsGLB_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NMuonValidHitsGLB_eta.push_back(ID.make<TProfile>(histname,"NMuonValidHitsGLB_vsEta; #eta ; # hits", 48.,-2.4,2.4, 0., 60.));
    char nMuonValidHitsGLB_phi[]="NMuonValidHitsGLB_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", nMuonValidHitsGLB_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NMuonValidHitsGLB_phi.push_back(ID.make<TProfile>(histname,"NMuonValidHitsGLB_vsPhi; #phi; # hits", 48.,-TMath::Pi(),TMath::Pi(), 0., 60.));

    char pixelHitsTRK[]="PixelHitsTRK";
    sprintf(histname, "%s minEta %f maxEta %f", pixelHitsTRK,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PixelHitsTRK.push_back(ID.make<TH1F>("PixelHitsTRK","PixelHitsTRK; # hits; # entries", 10, 0., 10.));
    char pixelHitsTRK_pt[]="PixelHitsTRK_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", pixelHitsTRK_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PixelHitsTRK_pt.push_back(ID.make<TProfile>(histname,"PixelHitsTRK_vsPt; p_{T} (GeV); # hits", 50.,0.,150., 0., 10.));
    char pixelHitsTRK_eta[]="PixelHitsTRK_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", pixelHitsTRK_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PixelHitsTRK_eta.push_back(ID.make<TProfile>(histname,"PixelHitsTRK_vsEta; #eta ; # hits", 48.,-2.4,2.4, 0., 10.));
    char pixelHitsTRK_phi[]="PixelHitsTRK_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", pixelHitsTRK_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PixelHitsTRK_phi.push_back(ID.make<TProfile>(histname,"PixelHitsTRK_vsPhi; #phi; # hits", 48.,-TMath::Pi(),TMath::Pi(), 0., 10.));

    char pixelLayersTRK[]="PixelLayersTRK";
    sprintf(histname, "%s minEta %f maxEta %f", pixelLayersTRK,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PixelLayersTRK.push_back(ID.make<TH1F>("PixelLayersTRK","PixelLayersTRK; # layers; # entries", 10, 0., 10.));
    char pixelLayersTRK_pt[]="PixelLayersTRK_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", pixelLayersTRK_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PixelLayersTRK_pt.push_back(ID.make<TProfile>(histname,"PixelLayersTRK_vsPt; p_{T} (GeV); # layers", 50.,0.,150., 0., 10.));
    char pixelLayersTRK_eta[]="PixelLayersTRK_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", pixelLayersTRK_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PixelLayersTRK_eta.push_back(ID.make<TProfile>(histname,"PixelLayersTRK_vsEta; #eta ; # layers", 48.,-2.4,2.4, 0., 10.));
    char pixelLayersTRK_phi[]="PixelLayersTRK_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", pixelLayersTRK_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PixelLayersTRK_phi.push_back(ID.make<TProfile>(histname,"PixelLayersTRK_vsPhi; #phi; # layers", 48.,-TMath::Pi(),TMath::Pi(), 0., 10.));

    char trackerLayersTRK[]="TrackerLayersTRK";
    sprintf(histname, "%s minEta %f maxEta %f", trackerLayersTRK,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrackerLayersTRK.push_back(ID.make<TH1F>("TrackerLayersTRK","TrackerLayersTRK; # layers; # entries", 30, 0., 30.));
    char trackerLayersTRK_pt[]="TrackerLayersTRK_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", trackerLayersTRK_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrackerLayersTRK_pt.push_back(ID.make<TProfile>(histname,"TrackerLayersTRK_vsPt; p_{T} (GeV); # layers", 50.,0.,150., 0., 30.));
    char trackerLayersTRK_eta[]="TrackerLayersTRK_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", trackerLayersTRK_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrackerLayersTRK_eta.push_back(ID.make<TProfile>(histname,"TrackerLayersTRK_vsEta; #eta ; # layers", 48.,-2.4,2.4, 0., 30.));
    char trackerLayersTRK_phi[]="TrackerLayersTRK_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", trackerLayersTRK_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrackerLayersTRK_phi.push_back(ID.make<TProfile>(histname,"TrackerLayersTRK_vsPhi; #phi; # layers", 48.,-TMath::Pi(),TMath::Pi(), 0., 30.));
 
    char hitFractionTRK[]="HitFractionTRK";
    sprintf(histname, "%s minEta %f maxEta %f", hitFractionTRK,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    HitFractionTRK.push_back(ID.make<TH1F>(histname, "HitFractionTRK; fraction; # entries", 20, 0., 1.));
    char hitFractionTRK_pt[]="HitFractionTRK_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", hitFractionTRK_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    HitFractionTRK_pt.push_back(ID.make<TProfile>(histname,"HitFractionTRK_vsPt; p_{T} (GeV); fraction", 50.,0.,150., 0., 1.));
    char hitFractionTRK_eta[]="HitFractionTRK_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", hitFractionTRK_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    HitFractionTRK_eta.push_back(ID.make<TProfile>(histname,"HitFractionTRK_vsEta; #eta ; fraction ", 48.,-2.4,2.4, 0., 1.));
    char hitFractionTRK_phi[]="HitFractionTRK_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", hitFractionTRK_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    HitFractionTRK_phi.push_back(ID.make<TProfile>(histname,"HitFractionTRK_vsPhi; #phi; fraction", 48.,-TMath::Pi(),TMath::Pi(), 0., 1.));

    char trkStaChi2[]="TrkStaChi2";
    sprintf(histname, "%s minEta %f maxEta %f", trkStaChi2,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrkStaChi2.push_back(ID.make<TH1F>("TrkStaChi2", "TrkStaChi2; chi2; # entries", 50, 0., 100.));
    char trkStaChi2_pt[]="TrkStaChi2_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", trkStaChi2_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrkStaChi2_pt.push_back(ID.make<TProfile>(histname,"TrkStaChi2_vsPt; p_{T} (GeV); chi2", 50.,0.,150., 0., 100.));
    char trkStaChi2_eta[]="TrkStaChi2_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", trkStaChi2_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrkStaChi2_eta.push_back(ID.make<TProfile>(histname,"TrkStaChi2_vsEta; #eta ; chi2 ", 48.,-2.4,2.4, 0., 100.));
    char trkStaChi2_phi[]="TrkStaChi2_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", trkStaChi2_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrkStaChi2_phi.push_back(ID.make<TProfile>(histname,"TrkStaChi2_vsPhi; #phi; chi2", 48.,-TMath::Pi(),TMath::Pi(), 0., 100.));

    char trkKink[]="TrkKink";
    sprintf(histname, "%s minEta %f maxEta %f", trkKink,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrkKink.push_back(ID.make<TH1F>("TrkKink", "TrkKink; prob.; # entries", 100, 0., 250.));
    char trkKink_pt[]="TrkKink_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", trkKink_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrkKink_pt.push_back(ID.make<TProfile>(histname,"TrkKink_vsPt; p_{T} (GeV); prob.", 50.,0.,150., 0., 250.));
    char trkKink_eta[]="TrkKink_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", trkKink_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrkKink_eta.push_back(ID.make<TProfile>(histname,"TrkKink_vsEta; #eta ; prob ", 48.,-2.4,2.4, 0., 250.));
    char trkKink_phi[]="TrkKink_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", trkKink_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrkKink_phi.push_back(ID.make<TProfile>(histname,"TrkKink_vsPhi; #phi; prob", 48.,-TMath::Pi(),TMath::Pi(), 0., 250.));

    char segmentComp[]="SegmentComp";
    sprintf(histname, "%s minEta %f maxEta %f", segmentComp,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    SegmentComp.push_back(ID.make<TH1F>("SegmentComp" ,"SegmentComp; prob.; segmentComp", 100, 0., 1.));
    char segmentComp_pt[]="SegmentComp_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", segmentComp_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    SegmentComp_pt.push_back(ID.make<TProfile>(histname,"SegmentComp_vsPt; p_{T} (GeV); prob.", 50.,0.,150., 0., 1.));
    char segmentComp_eta[]="SegmentComp_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", segmentComp_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    SegmentComp_eta.push_back(ID.make<TProfile>(histname,"SegmentComp_vsEta; #eta ; prob ", 48.,-2.4,2.4, 0., 1.));
    char segmentComp_phi[]="SegmentComp_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", segmentComp_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    SegmentComp_phi.push_back(ID.make<TProfile>(histname,"SegmentComp_vsPhi; #phi; prob", 48.,-TMath::Pi(),TMath::Pi(), 0., 1.));

    char  dxy []=" Dxy ";
    sprintf(histname, "%s minEta %f maxEta %f",  dxy ,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Dxy.push_back(ID.make<TH1F>("Dxy", "Dxy; d_{xy} (cm); # entries", 50,-0.5,0.5));
    char dxy_pt[]="Dxy_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", dxy_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Dxy_pt.push_back(ID.make<TProfile>(histname,"Dxy_vsPt; p_{T} (GeV); d_{xy} (cm)", 50.,0.,150., -0.5, 0.5));
    char dxy_eta[]="Dxy_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", dxy_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Dxy_eta.push_back(ID.make<TProfile>(histname,"Dxy_vsEta; #eta ; d_{xy} (cm)", 48.,-2.4,2.4,-0.5,0.5));
    char dxy_phi[]="Dxy_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", dxy_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Dxy_phi.push_back(ID.make<TProfile>(histname,"Dxy_vsPhi; #phi; d_{xy} (cm)", 48.,-TMath::Pi(),TMath::Pi(),-0.5,0.5));

    char  dz []=" Dz ";
    sprintf(histname, "%s minEta %f maxEta %f",  dz ,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Dz.push_back(ID.make<TH1F>("Dz" , "Dz; d_{z} (cm); # entries", 60,-1.5,1.5));
    char dz_pt[]="Dz_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", dz_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Dz_pt.push_back(ID.make<TProfile>(histname,"Dz_vsPt; p_{T} (GeV); d_{z} (cm)", 50.,0.,150., -1.5, 1.5));
    char dz_eta[]="Dz_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", dz_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Dz_eta.push_back(ID.make<TProfile>(histname,"Dz_vsEta; #eta ; d_{z} (cm)", 48.,-2.4,2.4,-1.5,1.5));
    char dz_phi[]="Dz_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", dz_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    Dz_phi.push_back(ID.make<TProfile>(histname,"Dxy_vsPhi; #phi; d_{z} (cm)", 48.,-TMath::Pi(),TMath::Pi(),-1.5,1.5));

    char qOverPtTrkSTa[]="qOverPtTrkSta";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkSTa ,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkSta.push_back(ID.make<TH1F>("qOverPtTrkSta", "qOverPtTrkSta; q/p_{T}^{sta} - q/p_{T}^{trk}; # entries", 50,-0.05,0.05)); 
    char qOverPtTrkSTa_pt[]="qOverPtTrkSt_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkSTa_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkSta_pt.push_back(ID.make<TProfile>(histname,"qOverPtTrkSta_vsPt; p_{T} (GeV); q/p_{T}^{sta} - q/p_{T}^{trk}", 50.,0.,150., -0.05, 0.05));
    char qOverPtTrkSTa_eta[]="qOverPtTrkSta_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkSTa_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkSta_eta.push_back(ID.make<TProfile>(histname,"qOverPtTrkSta_vsEta; #eta ; q/p_{T}^{sta} - q/p_{T}^{trk}", 48.,-2.4,2.4,-0.05,0.05));
    char qOverPtTrkSTa_phi[]="qOverPtTrkSta_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkSTa_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkSta_phi.push_back(ID.make<TProfile>(histname,"qOverPtTrkSta_vsPhi; #phi; q/p_{T}^{sta} - q/p_{T}^{trk}", 48.,-TMath::Pi(),TMath::Pi(),-0.05,0.05));

    char qOverPtTrkSTa200[]="qOverPtTrkSta200";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkSTa200 ,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkSta200.push_back(ID.make<TH1F>("qOverPtTrkSta200", "qOverPtTrkSta200; q/p_{T}^{sta} - q/p_{T}^{trk}; # entries", 50,-0.2,0.2));
    char qOverPtTrkSTa200_pt[]="qOverPtTrkSta200_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkSTa200_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkSta200_pt.push_back(ID.make<TProfile>(histname,"qOverPtTrkSta200_vsPt; p_{T} (GeV); q/p_{T}^{sta} - q/p_{T}^{trk}", 50.,0.,150., -0.2, 0.2));
    char qOverPtTrkSTa200_eta[]="qOverPtTrkSta200_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkSTa200_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkSta200_eta.push_back(ID.make<TProfile>(histname,"qOverPtTrkSta200_vsEta; #eta ; q/p_{T}^{sta} - q/p_{T}^{trk}", 48.,-2.4,2.4,-0.2,0.2));
    char qOverPtTrkSTa200_phi[]="qOverPtTrkSta200_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkSTa200_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkSta200_phi.push_back(ID.make<TProfile>(histname,"qOverPtTrkSta200_vsPhi; #phi; q/p_{T}^{sta} - q/p_{T}^{trk}", 48.,-TMath::Pi(),TMath::Pi(),-0.2,0.2));

    char qOverPtTrkGlB[]="qOverPtTrkGlb";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkGlB ,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkGlb.push_back(ID.make<TH1F>("qOverPtTrkGlb", "qOverPtTrkGlb; q/p_{T}^{glb} - q/p_{T}^{trk}; # entries", 50,-0.03,0.03));
    char qOverPtTrkGlB_pt[]="qOverPtTrkGlb_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkGlB_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkGlb_pt.push_back(ID.make<TProfile>(histname,"qOverPtTrkGlb_vsPt; p_{T} (GeV); q/p_{T}^{glb} - q/p_{T}^{trk}", 50.,0.,150., -0.03, 0.03));
    char qOverPtTrkGlB_eta[]="qOverPtTrkGlb_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkGlB_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkGlb_eta.push_back(ID.make<TProfile>(histname,"qOverPtTrkGlb_vsEta; #eta ; q/p_{T}^{glb} - q/p_{T}^{trk}", 48.,-2.4,2.4,-0.03,0.03));
    char qOverPtTrkGlB_phi[]="qOverPtTrkGlb_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkGlB_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkGlb_phi.push_back(ID.make<TProfile>(histname,"qOverPtTrkGlb_vsPhi; #phi; q/p_{T}^{glb} - q/p_{T}^{trk}", 48.,-TMath::Pi(),TMath::Pi(),-0.03,0.03));

    char qOverPtTrkGlB200[]="qOverPtTrkGlb200";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkGlB200 ,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkGlb200.push_back(ID.make<TH1F>("qOverPtTrkGlb200", "qOverPtTrkGlb200; q/p_{T}^{glb} - q/p_{T}^{trk}; # entries", 50,-0.2,0.2)); 
    char qOverPtTrkGlB200_pt[]="qOverPtTrkGlb200_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkGlB200_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkGlb200_pt.push_back(ID.make<TProfile>(histname,"qOverPtTrkGlb200_vsPt; p_{T} (GeV); q/p_{T}^{glb} - q/p_{T}^{trk}", 50.,0.,150., -0.2, 0.2));
    char qOverPtTrkGlB200_eta[]="qOverPtTrkGlb200_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkGlB200_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkGlb200_eta.push_back(ID.make<TProfile>(histname,"qOverPtTrkGlb200_vsEta; #eta ; q/p_{T}^{glb} - q/p_{T}^{trk}", 48.,-2.4,2.4,-0.2,0.2));
    char qOverPtTrkGlB200_phi[]="qOverPtTrkGlb200_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtTrkGlB200_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtTrkGlb200_phi.push_back(ID.make<TProfile>(histname,"qOverPtTrkGlb200_vsPhi; #phi; q/p_{T}^{glb} - q/p_{T}^{trk}", 48.,-TMath::Pi(),TMath::Pi(),-0.2,0.2));  

    char qOverPtGlB[]="qOverPtGlobal";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtGlB ,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtGlobal.push_back(ID.make<TH1F>("qOverPtGlb","qOverPtGlb; q/p_{T}^{glb}; # entries",50,-0.2,0.2));
    char qOverPtGLB_pt[]="qOverPtGlobal_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtGLB_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtGlobal_pt.push_back(ID.make<TProfile>(histname,"qOverPtGlobal_vsPt; p_{T} (GeV); q/p_{T}^{glb}", 50.,0.,150., -0.2, 0.2));
    char qOverPtGLB_eta[]="qOverPtGlobal_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtGLB_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtGlobal_eta.push_back(ID.make<TProfile>(histname,"qOverPtGlobal_vsEta; #eta ; q/p_{T}^{glb}", 48.,-2.4,2.4,-0.2,0.2));
    char qOverPtGLB_phi[]="qOverPtGlobal_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f", qOverPtGLB_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    qOverPtGlobal_phi.push_back(ID.make<TProfile>(histname,"qOverPtGlobal_vsPhi; #phi; q/p_{T}^{glb}", 48.,-TMath::Pi(),TMath::Pi(),-0.2,0.2));

  }


  TFileDirectory efficiency=fileService->mkdir("Efficiency");
  TFileDirectory medium=efficiency.mkdir("Medium");
  Medium_Numerator_eta = medium.make<TH1F>("Medium_Numerator_eta", "Medium_Numerator_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_Numerator_pt = medium.make<TH1F>("Medium_Numerator_pt", "Medium_Numerator_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_Numerator_phi= medium.make<TH1F>("Medium_Numerator_phi", "Medium_Numerator_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_Step0_eta = medium.make<TH1F>("Medium_Step0_eta", "Medium_Step0_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_Step0_pt = medium.make<TH1F>("Medium_Step0_pt", "Medium_Step0_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_Step0_phi = medium.make<TH1F>("Medium_Step0_phi", "Medium_Step0_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_Step1_eta = medium.make<TH1F> ("Medium_Step1_eta", "Medium_Step1_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_Step1_pt = medium.make<TH1F>("Medium_Step1_pt", "Medium_Step1_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_Step1_phi= medium.make<TH1F>("Medium_Step1_phi", "Medium_Step1_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_Step2_eta = medium.make<TH1F>("Medium_Step2_eta", "Medium_Step2_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_Step2_pt = medium.make<TH1F>("Medium_Step2_pt", "Medium_Step2_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_Step2_phi= medium.make<TH1F>("Medium_Step2_phi", "Medium_Step2_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_Step1not2_eta = medium.make<TH1F>("Medium_Step1not2_eta", "Medium_Step1not2_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_Step1not2_pt = medium.make<TH1F>("Medium_Step1not2_pt", "Medium_Step1not2_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_Step1not2_phi = medium.make<TH1F>("Medium_Step1not2_phi", "Medium_Step1not2_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_Step2not1_eta = medium.make<TH1F>("Medium_Step2not1_eta", "Medium_Step2not1_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_Step2not1_pt = medium.make<TH1F>("Medium_Step2not1_pt", "Medium_Step2not1_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_Step2not1_phi= medium.make<TH1F>("Medium_Step2not1_phi", "Medium_Step2not1_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_isGlobal_eta = medium.make<TH1F>("Medium_isGlobal_eta", "Medium_isGlobal_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_isGlobal_pt = medium.make<TH1F>("Medium_isGlobal_pt", "Medium_isGlobal_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_isGlobal_phi= medium.make<TH1F>("Medium_isGlobal_phi", "Medium_isGlobal_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_glbNormChi2_eta = medium.make<TH1F>("Medium_glbNormChi2_eta", "Medium_glbNormChi2_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_glbNormChi2_pt = medium.make<TH1F>("Medium_glbNormChi2_pt", "Medium_glbNormChi2_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_glbNormChi2_phi = medium.make<TH1F>("Medium_glbNormChi2_phi", "Medium_glbNormChi2_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_trkStaChi2_eta = medium.make<TH1F>("Medium_trkStaChi2_eta", "Medium_trkStaChi2_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_trkStaChi2_pt = medium.make<TH1F>("Medium_trkStaChi2_pt", "Medium_trkStaChi2_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_trkStaChi2_phi = medium.make<TH1F>("Medium_trkStaChi2_phi", "Medium_trkStaChi2_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_trkKink_eta = medium.make<TH1F>("Medium_trkKink_eta", "Medium_trkKink_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_trkKink_pt = medium.make<TH1F>("Medium_trkKink_pt", "Medium_trkKink_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_trkKink_phi = medium.make<TH1F>("Medium_trkKink_phi", "Medium_trkKink; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Medium_muSegmCompL_eta = medium.make<TH1F>("Medium_muSegmCompL_eta", "Medium_muSegmCompL_eta; #eta; # entries", 48, -2.4, 2.4);
  Medium_muSegmCompL_pt = medium.make<TH1F>("Medium_muSegmCompL_pt", "Medium_muSegmComL_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Medium_muSegmCompL_phi = medium.make<TH1F>("Medium_muSegmCompL_phi", "MEdium_muSegmCompL_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  TFileDirectory tight=efficiency.mkdir("Tight");
  Tight_Numerator_eta =tight.make<TH1F>("Tight_Numerator_eta", "Tight_Numerator_eta; #eta; # entries", 48, -2.4, 2.4);
  Tight_Numerator_pt = tight.make<TH1F>("Tight_Numerator_pt", "Tight_Numerator_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_Numerator_phi = tight.make<TH1F>("Tight_Numerator_phi", "Tight_Numerator_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_isGlobal_eta = tight.make<TH1F>("Tight_isGlobal_eta", "Tight_isGlobal_eta; #eta; # entries", 48, -2.4, 2.4);
  Tight_isGlobal_pt = tight.make<TH1F>("Tight_isGlobal_pt", "Tight_isGlobal_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_isGlobal_phi = tight.make<TH1F>("Tight_isGlobal_phi","Tight_isGlobal_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_isPF_eta = tight.make<TH1F>("Tight_isPF_eta", "Tight_isPF_eta; #eta;# entries", 48, -2.4, 2.4);
  Tight_isPF_pt = tight.make<TH1F>("Tight_isPF_pt", "Tight_isPF_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_isPF_phi = tight.make<TH1F>("Tight_isPF_phi", "Tight_isPF_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_glbNormChi2_eta = tight.make<TH1F>("Tight_glbNormChi2_eta", "Tight_glbNormChi2_eta; #eta; # entries", 48, -2.4, 2.4);
  Tight_glbNormChi2_pt = tight.make<TH1F>("Tight_glbNormChi2_pt", "Tight_glbNormChi2_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_glbNormChi2_phi = tight.make<TH1F>("Tight_glbNormChi2_phi", "Tight_glbNormChi2_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_glbMuonValidHits_eta = tight.make<TH1F>("Tight_glbMuonValidHits_eta", "Tight_glbMuonValidHits; #eta; # entries", 48, -2.4, 2.4);
  Tight_glbMuonValidHits_pt = tight.make<TH1F>("Tight_glbMuonValidHits_pt", "Tight_glbMuonValidHits_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_glbMuonValidHits_phi = tight.make<TH1F>("Tight_glbMuonValidHits_phi", "Tight_glbMuonValidHits_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_trkMuonMatchedStations_eta =tight.make<TH1F>("Tight_trkMuonMatchedStations_eta","Tight_trkMuonMatchedStations_eta; #eta; # entries", 48, -2.4, 2.4);
  Tight_trkMuonMatchedStations_pt = tight.make<TH1F>("Tight_trkMuonMatchedStations_pt", "Tight_trkMuonMatchedStations_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_trkMuonMatchedStations_phi = tight.make<TH1F>("Tight_trkMuonMatchedStations_phi", "Tight_trkMuonMatchedStations_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_dxy_eta = tight.make<TH1F>("Tight_dxy_eta", "Tight_dxy_eta; #eta; # entries", 48, -2.4, 2.4);
  Tight_dxy_pt = tight.make<TH1F>("Tight_dxy_pt", "Tight_dxy_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_dxy_phi = tight.make<TH1F>("Tight_dxy_phi", "Tight_dxy_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_dz_eta = tight.make<TH1F>("Tight_dz_eta", "Tight_dz_eta; #eta; # entries", 48, -2.4, 2.4);
  Tight_dz_pt = tight.make<TH1F>("Tight_dz_pt", "Tight_dz_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_dz_phi = tight.make<TH1F>("Tight_dz_phi", "Tight_dz_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_trkPixelValidHits_eta= tight.make<TH1F>("Tight_trkPixelValidHits_eta", "Tight_trkPixelValidHits_eta; #eta; # entries", 48, -2.4, 2.4);
  Tight_trkPixelValidHits_pt = tight.make<TH1F>("Tight_trkPixelValidHits_pt", "Tight_trkPixelValidHits_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_trkPixelValidHits_phi= tight.make<TH1F>("Tight_trkPixelValidHits_phi", "Tight_trkPixelValidHits_phi; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());
  Tight_trkTrackerLayerWithMeas_eta = tight.make<TH1F>("Tight_trkTrackerLayerWithMeas_eta", "Tight_trkTrackerLayerWithMeas_eta; #eta; # entries", 48, -2.4, 2.4);
  Tight_trkTrackerLayerWithMeas_pt = tight.make<TH1F>("Tight_trkTrackerLayerWithMeas_pt", "Tight_trkTrackerLAyerWithMeas_pt; p_{T} (GeV); # entries", 75, 0., 150.);
  Tight_trkTrackerLayerWithMeas_phi = tight.make<TH1F>("Tight_trkTrackerLayerWithMeas_phi", "Tight_trkTrackerLayerWithMeas; #phi; # entries", 48, -TMath::Pi(), TMath::Pi());


  TFileDirectory isolation=fileService->mkdir("Isolation");
  for (unsigned int iEtaRegion=0; iEtaRegion<6; iEtaRegion++){
    char histname [3];
    double *etaCutMin = etaMin.data();
    double *etaCutMax = etaMax.data();
    
    char photon[]="PhotonIso";
    sprintf(histname, "%s minEta %f maxEta %f",photon,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PhotonIso.push_back(isolation.make<TH1F> (histname, "PhotonIso; photon iso 0.4; # entries",50.,0.,10.));
    char photon_pt[]="PhotonIso_pt";
    sprintf(histname, "%s minEta %f maxEta %f",photon_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PhotonIso_pt.push_back(isolation.make<TProfile> (histname, "PhotonIso_vsPt; p_{T} (GeV);photon iso 0.4",50.,0.,150.,0.,10.));
    char photon_eta[]="PhotonIso_eta";
    sprintf(histname, "%s minEta %f maxEta %f",photon_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PhotonIso_eta.push_back(isolation.make<TProfile> (histname, "PhotonIso_vsEta; #eta; photon iso 0.4",48.,-2.4,2.4,0.,10.));
    char photon_phi[]="PhotonIso_phi";
    sprintf(histname, "%s minEta %f maxEta %f",photon_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    PhotonIso_phi.push_back(isolation.make<TProfile> (histname, "PhotonIso_vsPhi; p_{T} (GeV);photon iso 0.4",48.,-TMath::Pi(),TMath::Pi(),0.,10.));

    char hadPU[]="ChHadIsoPU";
    sprintf(histname, "%s minEta %f maxEta %f",hadPU,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ChHadIsoPU.push_back(isolation.make<TH1F> (histname, "ChHadIsoPU; PU Charged had. iso 0.4; # entries",50.,0.,10.));
    char hadPU_pt[]="ChHadIsoPU_pt";
    sprintf(histname, "%s minEta %f maxEta %f",hadPU_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ChHadIsoPU_pt.push_back(isolation.make<TProfile> (histname, "ChHadIsoPU_vsPt;  p_{T} (GeV);PU Charged had. iso 0.4",50.,0.,150.,0.,10.));
    char hadPU_eta[]="ChHadIsoPU_eta";
    sprintf(histname, "%s minEta %f maxEta %f",hadPU_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ChHadIsoPU_eta.push_back(isolation.make<TProfile> (histname, "ChHadIsoPU_vsEta; #eta; PU Charged had. iso 0.4",48.,-2.4,2.4,0.,10.));
    char hadPU_phi[]="ChHadIsoPU_phi";
    sprintf(histname, "%s minEta %f maxEta %f",hadPU_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ChHadIsoPU_phi.push_back(isolation.make<TProfile> (histname, "ChHadIsoPU_vsPhi; #phi; PU Charged had. iso 0.4",48.,-TMath::Pi(),TMath::Pi(),0.,10.));

    char neut[]="NeutralIso";
    sprintf(histname, "%s minEta %f maxEta %f",neut,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NeutralIso.push_back(isolation.make<TH1F> (histname, "NeutralIso; neutral had. iso 0.4; # entries",50,0,10));
    char neut_pt[]="NeutralIso_pt";
    sprintf(histname, "%s minEta %f maxEta %f",neut_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NeutralIso_pt.push_back(isolation.make<TProfile> (histname, "NeutralIso_vsPt; p{T} (GeV)",50.,0.,150.,0,10));
    char neut_eta[]="NeutralIso_eta";
    sprintf(histname, "%s minEta %f maxEta %f",neut_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NeutralIso_eta.push_back(isolation.make<TProfile> (histname, "NeutralIso_vsEta; #eta; neutral had. iso 0.4",48.,-2.4,2.4,0,10));
    char neut_phi[]="NeutralIso_phi";
    sprintf(histname, "%s minEta %f maxEta %f",neut_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    NeutralIso_phi.push_back(isolation.make<TProfile> (histname, "NeutralIso_vsPhi; #phi; neutral had. iso 0.4",48.,-TMath::Pi(),TMath::Pi(),0,10));

    char had[]="ChHadIso";
    sprintf(histname, "%s minEta %f maxEta %f",had,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ChHadIso.push_back(isolation.make<TH1F> (histname, "Charged Hadron Isolation;charged had. iso 0.4; # entries",50.,0.,10));
    char had_pt[]="ChHadIso_pt";
    sprintf(histname, "%s minEta %f maxEta %f",had_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ChHadIso_pt.push_back(isolation.make<TProfile> (histname, "ChargedHadronIsolation_vsPt; p{T} (GeV); charged had. iso 0.4",50.,0.,150.,0.,10));
    char had_eta[]="ChHadIso_eta";
    sprintf(histname, "%s minEta %f maxEta %f",had_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ChHadIso_eta.push_back(isolation.make<TProfile> (histname, "ChargedHadronIsolation_vsEta; #eta ; charged had. iso 0.4",48., -2.4, 2.4,0.,10));
    char had_phi[]="ChHadIso_phi";
    sprintf(histname, "%s minEta %f maxEta %f",had_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    ChHadIso_phi.push_back(isolation.make<TProfile> (histname, "ChargedHadronIsolation_vsPhi; #phi ; charged had. iso 0.4",48.,-TMath::Pi(),TMath::Pi(),0.,10));

    char dBeta[]="DBetaRelIso";
    sprintf(histname, "%s minEta %f maxEta %f",dBeta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    DBetaRelIso.push_back(isolation.make<TH1F>(histname , "DBetaRelIso;PFIso 0.4 (dBeta);# entries", 50,0.,2.));
    char dBeta_pt[]="DBetaRelIso_pt";
    sprintf(histname, "%s minEta %f maxEta %f",dBeta_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    DBetaRelIso_pt.push_back(isolation.make<TProfile>(histname , "DBetaRelIso_vsPt; p{T} (GeV); PFIso 0.4 (dBeta)", 50.,0.,150.,0.,2.));
    char dBeta_eta[]="DBetaRelIso_eta";
    sprintf(histname, "%s minEta %f maxEta %f",dBeta_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    DBetaRelIso_eta.push_back(isolation.make<TProfile>(histname , "DBetaRelIso_vsEta; #eta ; PFIso 0.4 (dBeta)", 48., -24, 2.4,0.,2.));
    char dBeta_phi[]="DBetaRelIso_phi";
    sprintf(histname, "%s minEta %f maxEta %f",dBeta_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    DBetaRelIso_phi.push_back(isolation.make<TProfile>(histname , "DBetaRelIso_vsPhi; #phi ; PFIso 0.4 (dBeta)", 48.,-TMath::Pi(),TMath::Pi(),0.,2.));

    char track[]="TrackerIso";
    sprintf(histname, "%s minEta %f maxEta %f",track,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrackerIso.push_back(isolation.make<TH1F>(histname  , "TrackerIso;tracker. iso 0.3;# entries", 50,0.,10.));
    char track_pt[]="TrackerIso_pt";
    sprintf(histname, "%s minEta %f maxEta %f",track_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrackerIso_pt.push_back(isolation.make<TProfile>(histname  , "TrackerIso_vsPt; p{T} (GeV);tracker. iso 0.3", 50.,0.,150.,0.,10.));
    char track_eta[]="TrackerIso_pt";
    sprintf(histname, "%s minEta %f maxEta %f",track_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrackerIso_eta.push_back(isolation.make<TProfile>(histname  , "TrackerIso_vsEta; #eta ; tracker. iso 0.3", 48., -2.4, 2.4,0.,10.));
    char track_phi[]="TrackerIso_phi";
    sprintf(histname, "%s minEta %f maxEta %f",track_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    TrackerIso_phi.push_back(isolation.make<TProfile>(histname  , "TrackerIso_vsPhi; #phi; tracker. iso 0.3",48.,-TMath::Pi(),TMath::Pi(),0.,10.));

    char emcal[]="EMCalIso";
    sprintf(histname, "%s minEta %f maxEta %f",emcal,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    EMCalIso.push_back(isolation.make<TH1F>(histname , "EMCalIso;EMCal. iso 0.3;# entries", 50,0.,10.));
    char emcal_pt[]="EMCalIso_pt";
    sprintf(histname, "%s minEta %f maxEta %f",emcal_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    EMCalIso_pt.push_back(isolation.make<TProfile>(histname , "EMCalIso_vsPt; p{T} (GeV); EMCal iso 0.3", 50.,0.,150.,0.,10.));
    char emcal_eta[]="EMCalIso_eta";
    sprintf(histname, "%s minEta %f maxEta %f",emcal_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    EMCalIso_eta.push_back(isolation.make<TProfile>(histname , "EMCalIso_vsEta; #eta; EMCal. iso 0.3",48.,-2.4,+2.4,0.,10.));
    char emcal_phi[]="EMCalIso_phi";
    sprintf(histname, "%s minEta %f maxEta %f",emcal_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    EMCalIso_phi.push_back(isolation.make<TProfile>(histname , "EMCalIso_vsPhi; #phi; EMCal. iso 0.3",48.,-TMath::Pi(),TMath::Pi(),0.,10.));

    char hcal[]="HCalIso";
    sprintf(histname, "%s minEta %f maxEta %f",hcal,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    HCalIso.push_back(isolation.make<TH1F>(histname, "HCalIso;HCal. iso 0.3;# entries", 50,0.,10.));
    char hcal_pt[]="HCalIso_pt";
    sprintf(histname, "%s minEta %f maxEta %f",hcal_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    HCalIso_pt.push_back(isolation.make<TProfile>(histname, "HCalIso_vsPt; p{T} (GeV); HCal. iso 0.3", 50.,0.,150.,0.,10.));
    char hcal_eta[]="HCalIso_eta";
    sprintf(histname, "%s minEta %f maxEta %f",hcal_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    HCalIso_eta.push_back(isolation.make<TProfile>(histname, "HCalIso_vsEta; #eta; HCal. iso 0.3",48.,-2.4,2.4,0.,10.));
    char hcal_phi[]="HCalIso_phi";
    sprintf(histname, "%s minEta %f maxEta %f",hcal_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    HCalIso_phi.push_back(isolation.make<TProfile>(histname, "HCalIso_vsPhi; #phi ; HCal. iso 0.3",48.,-TMath::Pi(),TMath::Pi(),0.,10.));

  }

  TFileDirectory momentum = fileService->mkdir("Momentum variables");
  for (unsigned int iEtaRegion=0; iEtaRegion<6; iEtaRegion++){
    char histname [3];
    double *etaCutMin = etaMin.data();
    double *etaCutMax = etaMax.data();

    char Mass[] ="GoodMuMass";
    sprintf(histname, "%s minEta %f maxEta %f",Mass,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMass_.push_back(momentum.make<TH1F>(histname,"GoodMuMass; invariant mass(Gev); # entries",30.,85.,115.));

    char MassPlus[] ="GoodMuMassPlus";
    sprintf(histname, "%s minEta %f maxEta %f",MassPlus,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMassPlus.push_back(momentum.make<TH1F>(histname,"GoodMuMassPlus; invariant mass (Gev); #entries",20.,86.5,96.5));
    char MassPlus_pt[] ="GoodMuMassPlus_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f",MassPlus_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMassPlus_pt.push_back(momentum.make<TProfile>(histname,"GoodMuMassPlus_vsPt; p_{T} (GeV); invariant mass (Gev)",50.,0.,150.,86.5,96.5));
    char MassPlus_eta[] ="GoodMuMassPlus_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f",MassPlus_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMassPlus_eta.push_back(momentum.make<TProfile>(histname,"GoodMuMassPlus_vsEta; #eta ; invariant mass (Gev)",48.,-2.4,2.4,86.5,96.5));
    char MassPlus_phi[] ="GoodMuMassPlus_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f",MassPlus_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMassPlus_phi.push_back(momentum.make<TProfile>(histname,"GoodMuMassPlus_vsPhi; #phi ; invariant mass (Gev)",48.,-TMath::Pi(),TMath::Pi(),86.5,96.5));

    char MassMinus[] ="GoodMuMassMinus";
    sprintf(histname, "%s minEta %f maxEta %f",MassMinus,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMassMinus.push_back(momentum.make<TH1F>(histname,"GoodMuMassMinus; invariant mass (Gev); #entries",20.,86.5,96.5));
    char MassMinus_pt[] ="GoodMuMassMinus_vsPt";
    sprintf(histname, "%s minEta %f maxEta %f",MassMinus_pt,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMassMinus_pt.push_back(momentum.make<TProfile>(histname,"GoodMuMassMinus_vsPt; p_{T} (GeV);invariant mass (GeV)",50.,0.,150.,86.5,96.5));
    char MassMinus_eta[] ="GoodMuMassMinus_vsEta";
    sprintf(histname, "%s minEta %f maxEta %f",MassMinus_eta,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMassMinus_eta.push_back(momentum.make<TProfile>(histname,"GoodMuMassMinus_vsEta; #eta ;invariant mass (GeV)",48.,-2.4,2.4,86.5,96.5));
    char MassMinus_phi[] ="GoodMuMassMinus_vsPhi";
    sprintf(histname, "%s minEta %f maxEta %f",MassMinus_phi,etaCutMin[iEtaRegion],etaCutMax[iEtaRegion]);
    GoodMuMassMinus_phi.push_back(momentum.make<TProfile>(histname,"GoodMuMassMinus_vsPhi; #phi;invariant mass (GeV)",48.,-TMath::Pi(),TMath::Pi(),86.5,96.5));

  }

} 

void VariableComparison::endJob(){
}

void VariableComparison::analyze(const edm::Event & iEvent,const edm::EventSetup& iSetup)
{

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(theMuonCollectionLabel_, muons);
  
  edm::Handle<edm::View<pat::TriggerObjectStandAlone>>  triggerObjects;
  iEvent.getByToken(trigger_, triggerObjects);

  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);

  edm::Handle<edm::View<reco::Vertex> > vertex;
  iEvent.getByToken(theVertexLabel_, vertex);


  /*look for a primary vertex*/
  reco::Vertex::Point posVtx;
  reco::Vertex::Error errVtx;
  unsigned int theIndexOfThePrimaryVertex = 999.;


  if (vertex.isValid()){
    for (unsigned int ind=0; ind<vertex->size(); ++ind) {
      if ( (*vertex)[ind].isValid() && !((*vertex)[ind].isFake()) ) {
        theIndexOfThePrimaryVertex = ind;
	break;
      }
    }
  }

  int nVtx = 0;
  edm::View<reco::Vertex>::const_iterator vertexIt  = vertex->begin();
  edm::View<reco::Vertex>::const_iterator vertexEnd = vertex->end();

  for (; vertexIt != vertexEnd; ++vertexIt) {
      const reco::Vertex& vertex = *vertexIt;
      if (!vertex.isValid()) continue;
      ++nVtx;
     
  }

  if (theIndexOfThePrimaryVertex<100) {
    posVtx = ((*vertex)[theIndexOfThePrimaryVertex]).position();
    errVtx = ((*vertex)[theIndexOfThePrimaryVertex]).error();
  } 
  

  /*primary vertex not found - use BeamSpot instead*/
  else {
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    iEvent.getByToken(theBeamSpotLabel_,recoBeamSpotHandle);
    reco::BeamSpot bs = *recoBeamSpotHandle;
   
    posVtx = bs.position();
    errVtx(0,0) = bs.BeamWidthX();
    errVtx(1,1) = bs.BeamWidthY();
    errVtx(2,2) = bs.sigmaZ();
  }

  const reco::Muon::ArbitrationType arbitrationType = reco::Muon::SegmentAndTrackArbitration; 
  const reco::Vertex vtx(posVtx,errVtx);

  if(!muons.isValid()) return;

  TLorentzVector Mu_tag, Mu_probe;
 
  float charge = 99.;
  float InvMass = -99.;

  double *EtaCutMin = etaMin.data();
  double *EtaCutMax = etaMax.data();


  for (edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon) {

    if (muon->isTimeValid()) {
      if (muon->isStandAloneMuon() && ((fabs(muon->eta()) < 1.2  && (muon->outerTrack()->numberOfValidHits()) > 20) ||
				(fabs(muon->eta()) >= 1.2 && (muon->outerTrack()->numberOfValidHits() > 11))))
        STAmuonTime->Fill(muon->time().timeAtIpInOut);

      if (muon->isStandAloneMuon() && fabs(muon->eta()) < 0.9  && (muon->outerTrack()->numberOfValidHits()) > 20)
        STAmuonTimeBarrel->Fill(muon->time().timeAtIpInOut);
      
      if (muon->isStandAloneMuon() && fabs(muon->eta()) < 1.2  && (muon->outerTrack()->numberOfValidHits()) > 11)  
        STAmuonTimeEndcap->Fill(muon->time().timeAtIpInOut);
    }
  }

  /*tag muon selection*/
  for (edm::View<pat::Muon>::const_iterator mu_tag=muons->begin(); mu_tag!=muons->end(); ++mu_tag){

    bool pathHasFired = false;
    if( !trigResults.failedToGet() ) {
      int N_Triggers = trigResults->size();
      const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);

      for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
        if (trigResults.product()->accept(i_Trig)) {
          TString TrigPath =trigName.triggerName(i_Trig);
          if (TrigPath == hlt_path) {
            pathHasFired = true;
          }
        }
      }
    }

    bool hltDr_match = false;
    const edm::TriggerNames & names = iEvent.triggerNames(*trigResults);
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      obj.unpackPathNames(names);
      obj.unpackFilterLabels(iEvent,*trigResults);
      double Deta = mu_tag->eta() - obj.eta();
      double Dphi = mu_tag->phi() - obj.phi();  
      if  (Dphi >   TMath::Pi()) Dphi -= 2*TMath::Pi();
      else if (Dphi <= -TMath::Pi()) Dphi += 2*TMath::Pi();
      double deltaR =sqrt(Deta*Deta + Dphi*Dphi);
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
        string myfillabl=obj.filterLabels()[h];
        if (deltaR < tag_hltDrCut && myfillabl == tag_hltFilter ) {
        hltDr_match= true;
        }
      }
    }

    if (!hltDr_match) continue;
    if (!noTrigger && !pathHasFired) continue;
    if ( not(mu_tag->pt() > tag_minPt) )  continue;
    if ( fabs(mu_tag->eta()) > 2.4 ) continue;
    if (not(((mu_tag->pfIsolationR04().sumChargedHadronPt+ std::max(0.,mu_tag->pfIsolationR04().sumPhotonEt+mu_tag->pfIsolationR04().sumNeutralHadronEt - 0.5*(mu_tag->pfIsolationR04().sumPUPt))) / (mu_tag->pt()))< tag_isoCut)) continue;
    if (tag_muonID == "LOOSE"  && !mu_tag->isLooseMuon()) continue;
    if (tag_muonID == "MEDIUM" && !mu_tag->isMediumMuon()) continue;
    if (tag_muonID == "TIGHT" && !mu_tag->isTightMuon(vtx)) continue;
    if (tag_muonID == "HIGHPT" && !mu_tag->isHighPtMuon(vtx)) continue;
    if (tag_muonID == "SOFT" && !mu_tag->isSoftMuon(vtx)) continue;


    /*probe muon selection*/
    for (edm::View<pat::Muon>::const_iterator mu_probe=muons->begin(); mu_probe!=muons->end(); ++mu_probe){
      if (mu_tag==mu_probe) continue;
      if ( not (mu_probe->pt() > probe_minPt ) ) continue;
      if (not (((mu_probe->isolationR03().sumPt)/mu_probe->pt()) < probe_trackerIsoCut)) continue;
      if ( fabs(mu_probe->eta()) > 2.4 ) continue;
      reco::TrackRef recoBestTrack_tag = mu_tag->muonBestTrack();
      reco::TrackRef recoBestTrack_probe = mu_probe->muonBestTrack();
      Mu_tag.SetPxPyPzE(recoBestTrack_tag->px(), recoBestTrack_tag->py(),recoBestTrack_tag->pz(), recoBestTrack_tag->p());
      Mu_probe.SetPxPyPzE(recoBestTrack_probe->px(), recoBestTrack_probe->py(),recoBestTrack_probe->pz(), recoBestTrack_probe->p());
   
      charge= recoBestTrack_tag->charge()*recoBestTrack_probe->charge();
      
      if(charge>0) continue;

      if(mu_probe->isGlobalMuon() || (mu_probe->isTrackerMuon() && mu_probe->numberOfMatches(arbitrationType)>0)){
        if((fabs(mu_probe->eta()) < 1.2  && (mu_probe->isGlobalMuon() && (mu_probe->outerTrack()->numberOfValidHits()) > 20)) || (fabs(mu_probe->eta()) >= 1.2 && (mu_probe->isGlobalMuon() &&(mu_probe->outerTrack()->numberOfValidHits())> 11)))
            UnbSTAmuonTime->Fill(mu_probe->time().timeAtIpInOut);

        if(fabs(mu_probe->eta()) < 0.9  && (mu_probe->isGlobalMuon()&&(mu_probe->outerTrack()->numberOfValidHits()) > 20))
            UnbSTAmuonTimeBarrel->Fill(mu_probe->time().timeAtIpInOut);
		
        if(fabs(mu_probe->eta()) > 1.2  && (mu_probe->isGlobalMuon() && (mu_probe->outerTrack()->numberOfValidHits()) > 11)) 
          UnbSTAmuonTimeEndcap->Fill(mu_probe->time().timeAtIpInOut);

        InvMass=(Mu_tag+Mu_probe).M();
        invMass->Fill(InvMass);
        if((InvMass > pair_minInvMass) &&  (InvMass < pair_maxInvMass) ) {
          invMassInRange->Fill(InvMass);
          dilepPt->Fill((Mu_tag+Mu_probe).Pt());
          nVertices->Fill(nVtx);


          if(mu_probe->isMediumMuon()) {
            Medium_Numerator_eta->Fill(mu_probe->eta());
            Medium_Numerator_pt->Fill(mu_probe->pt());
	    Medium_Numerator_phi->Fill(mu_probe->phi());
          }

	  if((mu_probe->innerTrack()->validFraction()) > 0.8 && mu_probe->isLooseMuon()) {
	    Medium_Step0_eta->Fill(mu_probe->eta());
	    Medium_Step0_pt->Fill(mu_probe->pt());
	    Medium_Step0_phi->Fill(mu_probe->phi());
	  
	    if(mu_probe->isGlobalMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 3. && (mu_probe->combinedQuality().chi2LocalPosition) < 12. && (mu_probe->combinedQuality().trkKink) < 20. && (mu_probe->segmentCompatibility( arbitrationType)) > 0.303) {
	      Medium_Step1_eta->Fill(mu_probe->eta());
	      Medium_Step1_pt->Fill(mu_probe->pt());
	      Medium_Step1_phi->Fill(mu_probe->phi());
		    
              if((mu_probe->segmentCompatibility( arbitrationType)) <= 0.451) { 
	        Medium_Step1not2_eta->Fill(mu_probe->eta());
	        Medium_Step1not2_pt->Fill(mu_probe->pt());
	        Medium_Step1not2_phi->Fill(mu_probe->phi());
	      }
            }

            else if((mu_probe->segmentCompatibility( arbitrationType)) > 0.451){
              Medium_Step2not1_eta->Fill(mu_probe->eta());
	      Medium_Step2not1_pt->Fill(mu_probe->pt());
	      Medium_Step2not1_phi->Fill(mu_probe->phi());
            }

			  
	    if(mu_probe->isGlobalMuon() &&  (mu_probe->combinedQuality().chi2LocalPosition)< 12. && (mu_probe->combinedQuality().trkKink) < 20. && (mu_probe->segmentCompatibility( arbitrationType)) > 0.303 ){
	      Medium_glbNormChi2_eta->Fill(mu_probe->eta());
	      Medium_glbNormChi2_pt->Fill(mu_probe->pt());
	      Medium_glbNormChi2_phi->Fill(mu_probe->phi());
	    }
			  
	    if(mu_probe->isGlobalMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 3. && (mu_probe->combinedQuality().trkKink) < 20. && (mu_probe->segmentCompatibility( arbitrationType))  > 0.303 ) {
	      Medium_trkStaChi2_eta->Fill(mu_probe->eta());
              Medium_trkStaChi2_pt->Fill(mu_probe->pt());
	      Medium_trkStaChi2_phi->Fill(mu_probe->phi());
	    }
			  
	    if(mu_probe->isGlobalMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 3. && (mu_probe->combinedQuality().chi2LocalPosition)< 12. && (mu_probe->segmentCompatibility( arbitrationType)) > 0.303 ){
	      Medium_trkKink_eta->Fill(mu_probe->eta());
	      Medium_trkKink_pt->Fill(mu_probe->pt());
	      Medium_trkKink_phi->Fill(mu_probe->phi());
	    }
		  
    	    if(mu_probe->isGlobalMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 3. && (mu_probe->combinedQuality().chi2LocalPosition) < 12. && (mu_probe->combinedQuality().trkKink) < 20.) {
	      Medium_muSegmCompL_eta->Fill(mu_probe->eta());
	      Medium_muSegmCompL_pt->Fill(mu_probe->pt());
	      Medium_muSegmCompL_phi->Fill(mu_probe->phi());
            }

            if((mu_probe->segmentCompatibility( arbitrationType)) > 0.451){
	      Medium_Step2_eta->Fill(mu_probe->eta());
	      Medium_Step2_pt->Fill(mu_probe->pt());
	      Medium_Step2_phi->Fill(mu_probe->phi());
	    }
          }

          if(mu_probe->isTightMuon(vtx)){
	    Tight_Numerator_eta->Fill(mu_probe->eta());
            Tight_Numerator_pt->Fill(mu_probe->pt());
	    Tight_Numerator_phi->Fill(mu_probe->phi());
          }
		      
		      		      
	  if(mu_probe->isGlobalMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 10. && (mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits()) > 0 && (mu_probe->numberOfMatchedStations()) > 1 && fabs(mu_probe->muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu_probe->muonBestTrack()->dz(vtx.position())) < 0.5 && (mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits())> 0 && (mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement()) > 5){
	    Tight_isPF_eta->Fill(mu_probe->eta());
            Tight_isPF_pt->Fill(mu_probe->pt());
	    Tight_isPF_phi->Fill(mu_probe->phi());
	  }
		     
  	  if(mu_probe->isGlobalMuon() && mu_probe->isPFMuon() && (mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits()) > 0 && (mu_probe->numberOfMatchedStations()) > 1 && fabs(mu_probe->muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu_probe->muonBestTrack()->dz(vtx.position())) < 0.5 && (mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits())> 0 && (mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement()) > 5){
	    Tight_glbNormChi2_eta->Fill(mu_probe->eta());
	    Tight_glbNormChi2_pt->Fill(mu_probe->pt());
	    Tight_glbNormChi2_phi->Fill(mu_probe->phi());
          }

          if(mu_probe->isGlobalMuon() && mu_probe->isPFMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 10. &&  (mu_probe->numberOfMatchedStations()) > 1 && fabs(mu_probe->muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu_probe->muonBestTrack()->dz(vtx.position())) < 0.5 && (mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits())> 0 && (mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement()) > 5){
	    Tight_glbMuonValidHits_eta->Fill(mu_probe->eta());
	    Tight_glbMuonValidHits_pt->Fill(mu_probe->pt());
            Tight_glbMuonValidHits_phi->Fill(mu_probe->phi());
          }
		      			      
	  if(mu_probe->isGlobalMuon() && mu_probe->isPFMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 10. && (mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits()) > 0  && fabs(mu_probe->muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu_probe->muonBestTrack()->dz(vtx.position())) < 0.5 && (mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits())> 0 && (mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement()) > 5){
            Tight_trkMuonMatchedStations_eta->Fill(mu_probe->eta());
	    Tight_trkMuonMatchedStations_pt->Fill(mu_probe->pt());
	    Tight_trkMuonMatchedStations_phi->Fill(mu_probe->phi());
	  }
	      			      
	  if(mu_probe->isGlobalMuon() && mu_probe->isPFMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 10. && (mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits()) > 0 && (mu_probe->numberOfMatchedStations()) > 1 && fabs(mu_probe->muonBestTrack()->dz(vtx.position())) < 0.5 && (mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits())> 0 && (mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement()) > 5){
	    Tight_dxy_eta->Fill(mu_probe->eta());
            Tight_dxy_pt->Fill(mu_probe->pt());
	    Tight_dxy_phi->Fill(mu_probe->phi());
          }
		      
	  if(mu_probe->isGlobalMuon() && mu_probe->isPFMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 10. && (mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits()) > 0 && (mu_probe->numberOfMatchedStations()) > 1 && fabs(mu_probe->muonBestTrack()->dxy(vtx.position())) < 0.2  && (mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits())> 0 && (mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement()) > 5){
	    Tight_dz_eta->Fill(mu_probe->eta());
	    Tight_dz_pt->Fill(mu_probe->pt());
	    Tight_dz_phi->Fill(mu_probe->phi());
          }
		      			      
	  if(mu_probe->isGlobalMuon() && mu_probe->isPFMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 10. && (mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits()) > 0 && (mu_probe->numberOfMatchedStations()) > 1 && fabs(mu_probe->muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu_probe->muonBestTrack()->dz(vtx.position())) < 0.5 && (mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement()) > 5){
	    Tight_trkPixelValidHits_eta->Fill(mu_probe->eta());
	    Tight_trkPixelValidHits_pt->Fill(mu_probe->pt());
	    Tight_trkPixelValidHits_phi->Fill(mu_probe->phi());
          }
		      			      
	  if(mu_probe->isGlobalMuon() && mu_probe->isPFMuon() && (mu_probe->globalTrack()->normalizedChi2()) < 10. && (mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits()) > 0 && (mu_probe->numberOfMatchedStations()) > 1 && fabs(mu_probe->muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(mu_probe->muonBestTrack()->dz(vtx.position())) < 0.5 && (mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits())> 0 ){
	    Tight_trkTrackerLayerWithMeas_eta->Fill(mu_probe->eta());
	    Tight_trkTrackerLayerWithMeas_pt->Fill(mu_probe->pt());
            Tight_trkTrackerLayerWithMeas_phi->Fill(mu_probe->phi());
          }

          if (probe_muonIDs == "LOOSE"  && !mu_probe->isLooseMuon()) continue;
          if (probe_muonIDs == "MEDIUM" && !mu_probe->isMediumMuon()) continue;
          if (probe_muonIDs == "TIGHT" &&  !mu_probe->isTightMuon(vtx)) continue;
          if (probe_muonIDs == "HIGHPT" && !mu_probe->isHighPtMuon(vtx)) continue;
          if (probe_muonIDs == "SOFT" && !mu_probe->isSoftMuon(vtx)) continue;

          if(((mu_probe->pfIsolationR04().sumChargedHadronPt+ std::max(0.,mu_probe->pfIsolationR04().sumPhotonEt+mu_probe->pfIsolationR04().sumNeutralHadronEt - 0.5*(mu_probe->pfIsolationR04().sumPUPt)))/(mu_probe->pt())) < 0.25 && InvMass > 60 && InvMass <115){
            for (unsigned int iEtaRegion=0; iEtaRegion<6; iEtaRegion++){
              if(recoBestTrack_tag->eta()>EtaCutMin[iEtaRegion] && recoBestTrack_tag->eta()<EtaCutMax[iEtaRegion] && recoBestTrack_probe->eta()>EtaCutMin[iEtaRegion] && recoBestTrack_probe->eta()<EtaCutMax[iEtaRegion]){

                GoodMuMass_[iEtaRegion]-> Fill(InvMass);
                  if (InvMass > 86.5 && InvMass < 96.5) {      
                  if ((mu_probe->charge())>0){
                    GoodMuMassPlus[iEtaRegion]->Fill(InvMass);
                    GoodMuMassPlus_pt[iEtaRegion]->Fill((Mu_tag+Mu_probe).Pt(),InvMass);
                    GoodMuMassPlus_eta[iEtaRegion]->Fill((Mu_tag+Mu_probe).Eta(),InvMass);
                    GoodMuMassPlus_phi[iEtaRegion]->Fill((Mu_tag+Mu_probe).Phi(),InvMass);
                  }
                  else{
                    GoodMuMassMinus[iEtaRegion]->Fill(InvMass); 
                    GoodMuMassMinus_pt[iEtaRegion]->Fill((Mu_tag+Mu_probe).Pt(),InvMass);
                    GoodMuMassMinus_eta[iEtaRegion]->Fill((Mu_tag+Mu_probe).Eta(),InvMass);
                    GoodMuMassMinus_phi[iEtaRegion]->Fill((Mu_tag+Mu_probe).Phi(),InvMass);
                  }
                }
              }
            }
          }
        
       
          for (unsigned int iEtaRegion=0; iEtaRegion<6; iEtaRegion++){
            if(recoBestTrack_tag->eta()>EtaCutMin[iEtaRegion] && recoBestTrack_tag->eta()<EtaCutMax[iEtaRegion] && recoBestTrack_probe->eta()>EtaCutMin[iEtaRegion] && recoBestTrack_probe->eta()<EtaCutMax[iEtaRegion]){
 
              if (mu_probe->isGlobalMuon()){
                NHitsGLB[iEtaRegion]->Fill(mu_probe->globalTrack()->numberOfValidHits());
                NHitsGLB_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->globalTrack()->numberOfValidHits());
                NHitsGLB_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->globalTrack()->numberOfValidHits());
                NHitsGLB_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->globalTrack()->numberOfValidHits());
                Chi2GLB[iEtaRegion]->Fill(mu_probe->globalTrack()->normalizedChi2());
                Chi2GLB_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->globalTrack()->normalizedChi2());
                Chi2GLB_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->globalTrack()->normalizedChi2());
                Chi2GLB_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->globalTrack()->normalizedChi2());
                NMuonValidHitsGLB[iEtaRegion]->Fill(mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits());
                NMuonValidHitsGLB_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits());
                NMuonValidHitsGLB_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits());
                NMuonValidHitsGLB_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->globalTrack()->hitPattern().numberOfValidMuonHits());
                TrkStaChi2[iEtaRegion]->Fill(mu_probe->combinedQuality().chi2LocalPosition);       
                TrkStaChi2_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->combinedQuality().chi2LocalPosition);
                TrkStaChi2_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->combinedQuality().chi2LocalPosition);
                TrkStaChi2_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->combinedQuality().chi2LocalPosition);
                TrkKink[iEtaRegion]->Fill(mu_probe->combinedQuality().trkKink);
                TrkKink_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->combinedQuality().trkKink);
                TrkKink_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->combinedQuality().trkKink);
                TrkKink_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->combinedQuality().trkKink);
              }

              if ((mu_probe->isTrackerMuon()  && mu_probe->numberOfMatches(arbitrationType)>0) && mu_probe->isGlobalMuon()){
                SegmentComp[iEtaRegion]->Fill(mu_probe->segmentCompatibility( arbitrationType));      
                SegmentComp_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->segmentCompatibility( arbitrationType));
                SegmentComp_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->segmentCompatibility( arbitrationType));
                SegmentComp_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->segmentCompatibility( arbitrationType));
                NHitsTRK[iEtaRegion]->Fill(mu_probe->innerTrack()->numberOfValidHits() );
                NHitsTRK_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->innerTrack()->numberOfValidHits() );
                NHitsTRK_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->innerTrack()->numberOfValidHits() );
                NHitsTRK_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->innerTrack()->numberOfValidHits() );
                Chi2TRK[iEtaRegion]->Fill(mu_probe->innerTrack()->normalizedChi2());
                Chi2TRK_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->innerTrack()->normalizedChi2());
                Chi2TRK_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->innerTrack()->normalizedChi2());
                Chi2TRK_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->innerTrack()->normalizedChi2());
                PixelHitsTRK[iEtaRegion]->Fill(mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits());
                PixelHitsTRK_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits());
                PixelHitsTRK_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits());
                PixelHitsTRK_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->innerTrack()->hitPattern().numberOfValidPixelHits());
                PixelLayersTRK[iEtaRegion]->Fill(mu_probe->innerTrack()->hitPattern().pixelLayersWithMeasurement());
                PixelLayersTRK_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->innerTrack()->hitPattern().pixelLayersWithMeasurement());
                PixelLayersTRK_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->innerTrack()->hitPattern().pixelLayersWithMeasurement());
                PixelLayersTRK_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->innerTrack()->hitPattern().pixelLayersWithMeasurement());
                TrackerLayersTRK[iEtaRegion]->Fill(mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement()); 
                TrackerLayersTRK_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement());
                TrackerLayersTRK_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement());
                TrackerLayersTRK_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->innerTrack()->hitPattern().trackerLayersWithMeasurement());
                HitFractionTRK[iEtaRegion]->Fill(mu_probe->innerTrack()->validFraction());   
                HitFractionTRK_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->innerTrack()->validFraction());   
                HitFractionTRK_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->innerTrack()->validFraction());   
                HitFractionTRK_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->innerTrack()->validFraction());   
                Dxy[iEtaRegion]->Fill(mu_probe->innerTrack()->dxy(vtx.position()));
                Dxy_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->innerTrack()->dxy(vtx.position()));
                Dxy_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->innerTrack()->dxy(vtx.position()));
                Dxy_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->innerTrack()->dxy(vtx.position()));
                Dz[iEtaRegion]->Fill(mu_probe->innerTrack()->dz(vtx.position()));
                Dz_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->innerTrack()->dz(vtx.position()));
                Dz_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->innerTrack()->dz(vtx.position()));
                Dz_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->innerTrack()->dz(vtx.position()));		
              }

              if (mu_probe->isStandAloneMuon() || mu_probe->isGlobalMuon()) {
                NHitsSTA[iEtaRegion]->Fill(mu_probe->outerTrack()->numberOfValidHits());
                NHitsSTA_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->outerTrack()->numberOfValidHits());
                NHitsSTA_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->outerTrack()->numberOfValidHits());
                NHitsSTA_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->outerTrack()->numberOfValidHits());
              }

              if (mu_probe->isTrackerMuon() && mu_probe->numberOfMatches(arbitrationType)>0 ) {
                NMatchedStation[iEtaRegion]->Fill(mu_probe->numberOfMatchedStations());
                NMatchedStation_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->numberOfMatchedStations());
                NMatchedStation_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->numberOfMatchedStations());
                NMatchedStation_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->numberOfMatchedStations());
              }

              if (mu_probe->isGlobalMuon()){
                float qOverPtTrk = (mu_probe->innerTrack()->charge()) / (mu_probe->innerTrack()->pt());
                float qOverPtSta = (mu_probe->outerTrack()->charge())/ (mu_probe->outerTrack()->pt());
                float qOverPtGlb = (mu_probe->globalTrack()->charge()) /(mu_probe->globalTrack()->pt());
 
                qOverPtGlobal[iEtaRegion]->Fill(qOverPtGlb);
                qOverPtGlobal_pt[iEtaRegion]->Fill(Mu_probe.Pt(),qOverPtGlb);
                qOverPtGlobal_eta[iEtaRegion]->Fill(Mu_probe.Eta(),qOverPtGlb);
                qOverPtGlobal_phi[iEtaRegion]->Fill(Mu_probe.Phi(),qOverPtGlb);
                qOverPtTrkSta[iEtaRegion]->Fill(qOverPtSta - qOverPtTrk); 
                qOverPtTrkSta_pt[iEtaRegion]->Fill(Mu_probe.Pt(),qOverPtSta - qOverPtTrk);
                qOverPtTrkSta_eta[iEtaRegion]->Fill(Mu_probe.Eta(),qOverPtSta - qOverPtTrk);
                qOverPtTrkSta_phi[iEtaRegion]->Fill(Mu_probe.Phi(),qOverPtSta - qOverPtTrk);
                qOverPtTrkGlb[iEtaRegion]->Fill(qOverPtGlb - qOverPtTrk); 
                qOverPtTrkGlb_pt[iEtaRegion]->Fill(Mu_probe.Pt(),qOverPtGlb - qOverPtTrk);
                qOverPtTrkGlb_eta[iEtaRegion]->Fill(Mu_probe.Eta(),qOverPtGlb - qOverPtTrk);
                qOverPtTrkGlb_phi[iEtaRegion]->Fill(Mu_probe.Phi(),qOverPtGlb - qOverPtTrk);

                if ((mu_probe->innerTrack()->pt()) > 200){
                  qOverPtTrkSta200[iEtaRegion]->Fill(qOverPtSta - qOverPtTrk);
                  qOverPtTrkSta200_pt[iEtaRegion]->Fill(Mu_probe.Pt(),qOverPtSta - qOverPtTrk);
                  qOverPtTrkSta200_eta[iEtaRegion]->Fill(Mu_probe.Eta(),qOverPtSta - qOverPtTrk);
                  qOverPtTrkSta200_phi[iEtaRegion]->Fill(Mu_probe.Phi(),qOverPtSta - qOverPtTrk);
                  qOverPtTrkGlb200[iEtaRegion]->Fill(qOverPtGlb - qOverPtTrk); 
                  qOverPtTrkGlb200_pt[iEtaRegion]->Fill(Mu_probe.Pt(),qOverPtGlb - qOverPtTrk);
                  qOverPtTrkGlb200_eta[iEtaRegion]->Fill(Mu_probe.Eta(),qOverPtGlb - qOverPtTrk);
                  qOverPtTrkGlb200_phi[iEtaRegion]->Fill(Mu_probe.Phi(),qOverPtGlb - qOverPtTrk);
                }
              }
     
              if (probe_muonIDs == "LOOSE"  && !mu_probe->isLooseMuon()) continue;
              if (probe_muonIDs == "MEDIUM" && !mu_probe->isMediumMuon()) continue;
              if (probe_muonIDs == "TIGHT" &&  !mu_probe->isTightMuon(vtx)) continue;
              if (probe_muonIDs == "HIGHPT" && !mu_probe->isHighPtMuon(vtx)) continue;
              if (probe_muonIDs == "SOFT" && !mu_probe->isSoftMuon(vtx)) continue;

              ProbePt[iEtaRegion]->Fill(Mu_probe.Pt());
              ProbePt_vsPt[iEtaRegion]->Fill(Mu_probe.Pt(),Mu_probe.Pt());
              ProbePt_vsEta[iEtaRegion]->Fill(Mu_probe.Eta(),Mu_probe.Pt());
              ProbePt_vsPhi[iEtaRegion]->Fill(Mu_probe.Phi(),Mu_probe.Pt());
              ProbeEta[iEtaRegion]->Fill(Mu_probe.Eta());
              ProbeEta_pt[iEtaRegion]->Fill(Mu_probe.Pt(),Mu_probe.Eta());
              ProbeEta_eta[iEtaRegion]->Fill(Mu_probe.Eta(),Mu_probe.Eta());
              ProbeEta_phi[iEtaRegion]->Fill(Mu_probe.Phi(),Mu_probe.Eta());
              ProbePhi[iEtaRegion]->Fill(Mu_probe.Phi());        
              ProbePhi_pt[iEtaRegion]->Fill(Mu_probe.Pt(),Mu_probe.Phi());
              ProbePhi_eta[iEtaRegion]->Fill(Mu_probe.Eta(),Mu_probe.Phi());
              ProbePhi_phi[iEtaRegion]->Fill(Mu_probe.Phi(),Mu_probe.Phi());

              PhotonIso[iEtaRegion]->Fill(mu_probe->pfIsolationR04().sumPhotonEt);
              PhotonIso_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->pfIsolationR04().sumPhotonEt);
              PhotonIso_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->pfIsolationR04().sumPhotonEt);
              PhotonIso_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->pfIsolationR04().sumPhotonEt);
              ChHadIso[iEtaRegion]->Fill(mu_probe->pfIsolationR04().sumChargedHadronPt);
              ChHadIso_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->pfIsolationR04().sumChargedHadronPt);
              ChHadIso_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->pfIsolationR04().sumChargedHadronPt);
              ChHadIso_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->pfIsolationR04().sumChargedHadronPt);
              ChHadIsoPU[iEtaRegion]->Fill(mu_probe->pfIsolationR04().sumPUPt);
              ChHadIsoPU_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->pfIsolationR04().sumPUPt);
              ChHadIsoPU_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->pfIsolationR04().sumPUPt);
              ChHadIsoPU_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->pfIsolationR04().sumPUPt);
              NeutralIso[iEtaRegion]->Fill(mu_probe->pfIsolationR04().sumNeutralHadronEt);
              NeutralIso_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->pfIsolationR04().sumNeutralHadronEt);
              NeutralIso_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->pfIsolationR04().sumNeutralHadronEt);
              NeutralIso_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->pfIsolationR04().sumNeutralHadronEt);
              DBetaRelIso[iEtaRegion]->Fill((mu_probe->pfIsolationR04().sumChargedHadronPt+ std::max(0.,mu_probe->pfIsolationR04().sumPhotonEt+mu_probe->pfIsolationR04().sumNeutralHadronEt - 0.5*(mu_probe->pfIsolationR04().sumPUPt)))/(mu_probe->pt()));
              DBetaRelIso_pt[iEtaRegion]->Fill(Mu_probe.Pt(),(mu_probe->pfIsolationR04().sumChargedHadronPt+ std::max(0.,mu_probe->pfIsolationR04().sumPhotonEt+mu_probe->pfIsolationR04().sumNeutralHadronEt - 0.5*(mu_probe->pfIsolationR04().sumPUPt)))/(mu_probe->pt()));
              DBetaRelIso_eta[iEtaRegion]->Fill(Mu_probe.Eta(),(mu_probe->pfIsolationR04().sumChargedHadronPt+ std::max(0.,mu_probe->pfIsolationR04().sumPhotonEt+mu_probe->pfIsolationR04().sumNeutralHadronEt - 0.5*(mu_probe->pfIsolationR04().sumPUPt)))/(mu_probe->pt()));
              DBetaRelIso_phi[iEtaRegion]->Fill(Mu_probe.Eta(),(mu_probe->pfIsolationR04().sumChargedHadronPt+ std::max(0.,mu_probe->pfIsolationR04().sumPhotonEt+mu_probe->pfIsolationR04().sumNeutralHadronEt - 0.5*(mu_probe->pfIsolationR04().sumPUPt)))/(mu_probe->pt()));
              TrackerIso[iEtaRegion]->Fill(mu_probe->isolationR03().sumPt);
              TrackerIso_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->isolationR03().sumPt);
              TrackerIso_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->isolationR03().sumPt);
              TrackerIso_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->isolationR03().sumPt);
              EMCalIso[iEtaRegion]->Fill(mu_probe->isolationR03().emEt);
              EMCalIso_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->isolationR03().emEt);
              EMCalIso_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->isolationR03().emEt);
              EMCalIso_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->isolationR03().emEt);
              HCalIso[iEtaRegion]->Fill(mu_probe->isolationR03().hadEt);
              HCalIso_pt[iEtaRegion]->Fill(Mu_probe.Pt(),mu_probe->isolationR03().hadEt);
              HCalIso_eta[iEtaRegion]->Fill(Mu_probe.Eta(),mu_probe->isolationR03().hadEt);
              HCalIso_phi[iEtaRegion]->Fill(Mu_probe.Phi(),mu_probe->isolationR03().hadEt);             
            }                     
          }
        }
      }
    }
  }
}
DEFINE_FWK_MODULE(VariableComparison);  
