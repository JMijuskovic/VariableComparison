#ifndef VariableComparison_h
#define VariableComparison_h

#include <memory>
#include <fstream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class TH1F;
class TH2F;

class VariableComparison : public edm::EDAnalyzer {
 public:
  /* Constructor */ 
  VariableComparison(const edm::ParameterSet& pSet);
  
  /* Destructor */ 
  ~VariableComparison() override ;
  
  /* Operations */ 
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  void beginJob() override;
  void endJob() override;

 private:
  MuonServiceProxy* theService;
  edm::ParameterSet parameters;

  double pair_minInvMass;
  double pair_maxInvMass;
  std::string hlt_path;
  std::string tag_hltFilter;
  double tag_hltDrCut;
  double tag_minPt;
  double tag_isoCut;
  std::string tag_muonID;
  double probe_minPt;
  double probe_trackerIsoCut;
  std::string probe_muonIDs;
  bool noTrigger;

  std::string EtaName[3];
  double EtaCutMin;
  double EtaCutMax;
  std::vector<double> etaMin;
  std::vector<double> etaMax;

  TH1F *UnbSTAmuonTime;
  TH1F *UnbSTAmuonTimeBarrel;
  TH1F *UnbSTAmuonTimeEndcap;
  TH1F *STAmuonTime;
  TH1F *STAmuonTimeBarrel;
  TH1F *STAmuonTimeEndcap;
  TH1F *invMass;
  TH1F *invMassInRange;
  TH1F *dilepPt;
  TH1F *nVertices;
  std::vector<TH1F *> PhotonIso;
  std::vector<TH1F *> ChHadIso;
  std::vector<TH1F *> ChHadIsoPU;
  std::vector<TH1F *> NeutralIso;
  std::vector<TH1F *> DBetaRelIso;
  std::vector<TH1F *> TrackerIso;  
  std::vector<TH1F *> EMCalIso;
  std::vector<TH1F *> HCalIso;
  std::vector<TH1F *> ProbePt;
  std::vector<TH1F *> ProbeEta;
  std::vector<TH1F *> ProbePhi;
  std::vector<TH1F *> GoodMuMass_;
  std::vector<TH1F *> GoodMuMassPlus;
  std::vector<TH1F *> GoodMuMassMinus;
  std::vector<TH1F *> NHitsGLB;
  std::vector<TH1F *> NHitsTRK;
  std::vector<TH1F *> NHitsSTA;
  std::vector<TH1F *> Chi2GLB;
  std::vector<TH1F *> Chi2TRK;
  std::vector<TH1F *> NMatchedStation;
  std::vector<TH1F *> NMuonValidHitsGLB;
  std::vector<TH1F *> PixelHitsTRK;
  std::vector<TH1F *> PixelLayersTRK;
  std::vector<TH1F *> TrackerLayersTRK;
  std::vector<TH1F *> HitFractionTRK;
  std::vector<TH1F *> TrkStaChi2;
  std::vector<TH1F *> TrkKink;
  std::vector<TH1F *> SegmentComp;
  std::vector<TH1F *> Dxy;
  std::vector<TH1F *> Dz;
  std::vector<TH1F *> qOverPtTrkSta;
  std::vector<TH1F *> qOverPtTrkSta200;
  std::vector<TH1F *> qOverPtTrkGlb;
  std::vector<TH1F *> qOverPtTrkGlb200;
  std::vector<TH1F *> qOverPtGlobal;
  TH1F *Medium_Numerator_eta;
  TH1F *Medium_Numerator_pt;
  TH1F *Medium_Numerator_phi;
  TH1F *Medium_Step0_eta;
  TH1F *Medium_Step0_pt;
  TH1F *Medium_Step0_phi;
  TH1F *Medium_Step1_eta;
  TH1F *Medium_Step1_pt;
  TH1F *Medium_Step1_phi;
  TH1F *Medium_Step2_eta;
  TH1F *Medium_Step2_pt;
  TH1F *Medium_Step2_phi;
  TH1F *Medium_Step1not2_eta;
  TH1F *Medium_Step1not2_pt;
  TH1F *Medium_Step1not2_phi;
  TH1F *Medium_Step2not1_eta;
  TH1F *Medium_Step2not1_pt;
  TH1F *Medium_Step2not1_phi;
  TH1F *Medium_isGlobal_eta;
  TH1F *Medium_isGlobal_pt;
  TH1F *Medium_isGlobal_phi;
  TH1F *Medium_glbNormChi2_eta;
  TH1F *Medium_glbNormChi2_pt;
  TH1F *Medium_glbNormChi2_phi;
  TH1F *Medium_trkStaChi2_eta;
  TH1F *Medium_trkStaChi2_pt;
  TH1F *Medium_trkStaChi2_phi;
  TH1F *Medium_trkKink_eta;
  TH1F *Medium_trkKink_pt;
  TH1F *Medium_trkKink_phi;
  TH1F *Medium_muSegmCompL_eta;
  TH1F *Medium_muSegmCompL_pt;
  TH1F *Medium_muSegmCompL_phi;
  TH1F *Tight_Numerator_eta;
  TH1F *Tight_Numerator_pt;
  TH1F *Tight_Numerator_phi;
  TH1F *Tight_isGlobal_eta;
  TH1F *Tight_isGlobal_pt;
  TH1F *Tight_isGlobal_phi;
  TH1F *Tight_isPF_eta;
  TH1F *Tight_isPF_pt;
  TH1F *Tight_isPF_phi;
  TH1F *Tight_glbNormChi2_eta;
  TH1F *Tight_glbNormChi2_pt;
  TH1F *Tight_glbNormChi2_phi;
  TH1F *Tight_glbMuonValidHits_eta;
  TH1F *Tight_glbMuonValidHits_pt;
  TH1F *Tight_glbMuonValidHits_phi;
  TH1F *Tight_trkMuonMatchedStations_eta;
  TH1F *Tight_trkMuonMatchedStations_pt;
  TH1F *Tight_trkMuonMatchedStations_phi;
  TH1F *Tight_dxy_eta;
  TH1F *Tight_dxy_pt;
  TH1F *Tight_dxy_phi;
  TH1F *Tight_dz_eta;
  TH1F *Tight_dz_pt;
  TH1F *Tight_dz_phi;
  TH1F *Tight_trkPixelValidHits_eta;
  TH1F *Tight_trkPixelValidHits_pt;
  TH1F *Tight_trkPixelValidHits_phi;
  TH1F *Tight_trkTrackerLayerWithMeas_eta;
  TH1F *Tight_trkTrackerLayerWithMeas_pt;
  TH1F *Tight_trkTrackerLayerWithMeas_phi;

  std::vector<TProfile *> PhotonIso_pt;
  std::vector<TProfile *> ChHadIso_pt;
  std::vector<TProfile *> ChHadIsoPU_pt;
  std::vector<TProfile *> NeutralIso_pt;
  std::vector<TProfile *> DBetaRelIso_pt;
  std::vector<TProfile *> TrackerIso_pt;
  std::vector<TProfile *> EMCalIso_pt;
  std::vector<TProfile *> HCalIso_pt;
  std::vector<TProfile *> PhotonIso_eta;
  std::vector<TProfile *> ChHadIso_eta;
  std::vector<TProfile *> ChHadIsoPU_eta;
  std::vector<TProfile *> NeutralIso_eta;
  std::vector<TProfile *> DBetaRelIso_eta;
  std::vector<TProfile *> TrackerIso_eta;
  std::vector<TProfile *> EMCalIso_eta;
  std::vector<TProfile *> HCalIso_eta;
  std::vector<TProfile *> PhotonIso_phi;
  std::vector<TProfile *> ChHadIso_phi;
  std::vector<TProfile *> ChHadIsoPU_phi;
  std::vector<TProfile *> NeutralIso_phi;
  std::vector<TProfile *> DBetaRelIso_phi;
  std::vector<TProfile *> TrackerIso_phi;
  std::vector<TProfile *> EMCalIso_phi;
  std::vector<TProfile *> HCalIso_phi;
  std::vector<TProfile *> ProbePt_vsPt;
  std::vector<TProfile *> ProbePt_vsEta;
  std::vector<TProfile *> ProbePt_vsPhi;
  std::vector<TProfile *> ProbeEta_pt;
  std::vector<TProfile *> ProbeEta_eta;
  std::vector<TProfile *> ProbeEta_phi;
  std::vector<TProfile *> ProbePhi_pt;
  std::vector<TProfile *> ProbePhi_eta;
  std::vector<TProfile *> ProbePhi_phi;
  std::vector<TProfile *> GoodMuMassMinus_pt;
  std::vector<TProfile *> GoodMuMassMinus_eta;
  std::vector<TProfile *> GoodMuMassMinus_phi;
  std::vector<TProfile *> GoodMuMassPlus_pt;
  std::vector<TProfile *> GoodMuMassPlus_eta;
  std::vector<TProfile *> GoodMuMassPlus_phi;
  std::vector<TProfile *> NHitsGLB_pt;
  std::vector<TProfile *> NHitsGLB_eta;
  std::vector<TProfile *> NHitsGLB_phi;
  std::vector<TProfile *> NHitsTRK_pt;
  std::vector<TProfile *> NHitsTRK_eta;
  std::vector<TProfile *> NHitsTRK_phi;
  std::vector<TProfile *> NHitsSTA_pt;
  std::vector<TProfile *> NHitsSTA_eta;
  std::vector<TProfile *> NHitsSTA_phi;
  std::vector<TProfile *> Chi2GLB_pt;
  std::vector<TProfile *> Chi2GLB_eta;
  std::vector<TProfile *> Chi2GLB_phi;
  std::vector<TProfile *> Chi2TRK_pt;
  std::vector<TProfile *> Chi2TRK_eta;
  std::vector<TProfile *> Chi2TRK_phi;
  std::vector<TProfile *> NMatchedStation_pt;
  std::vector<TProfile *> NMatchedStation_eta;
  std::vector<TProfile *> NMatchedStation_phi;
  std::vector<TProfile *> NMuonValidHitsGLB_pt;
  std::vector<TProfile *> NMuonValidHitsGLB_eta;
  std::vector<TProfile *> NMuonValidHitsGLB_phi;
  std::vector<TProfile *> PixelHitsTRK_pt;
  std::vector<TProfile *> PixelHitsTRK_eta;
  std::vector<TProfile *> PixelHitsTRK_phi;
  std::vector<TProfile *> PixelLayersTRK_pt;
  std::vector<TProfile *> PixelLayersTRK_eta;
  std::vector<TProfile *> PixelLayersTRK_phi;
  std::vector<TProfile *> TrackerLayersTRK_pt;
  std::vector<TProfile *> TrackerLayersTRK_eta;
  std::vector<TProfile *> TrackerLayersTRK_phi;
  std::vector<TProfile *> HitFractionTRK_pt;
  std::vector<TProfile *> HitFractionTRK_eta;
  std::vector<TProfile *> HitFractionTRK_phi;
  std::vector<TProfile *> TrkStaChi2_pt;
  std::vector<TProfile *> TrkStaChi2_eta;
  std::vector<TProfile *> TrkStaChi2_phi;
  std::vector<TProfile *> TrkKink_pt;
  std::vector<TProfile *> TrkKink_eta;
  std::vector<TProfile *> TrkKink_phi;
  std::vector<TProfile *> SegmentComp_pt;
  std::vector<TProfile *> SegmentComp_eta;
  std::vector<TProfile *> SegmentComp_phi;
  std::vector<TProfile *> Dxy_pt;
  std::vector<TProfile *> Dxy_eta;
  std::vector<TProfile *> Dxy_phi;
  std::vector<TProfile *> Dz_pt;
  std::vector<TProfile *> Dz_eta;
  std::vector<TProfile *> Dz_phi;
  std::vector<TProfile *> qOverPtTrkSta_pt;
  std::vector<TProfile *> qOverPtTrkSta_eta;
  std::vector<TProfile *> qOverPtTrkSta_phi;
  std::vector<TProfile *> qOverPtTrkSta200_pt;
  std::vector<TProfile *> qOverPtTrkSta200_eta;
  std::vector<TProfile *> qOverPtTrkSta200_phi;
  std::vector<TProfile *> qOverPtTrkGlb_pt;
  std::vector<TProfile *> qOverPtTrkGlb_eta;
  std::vector<TProfile *> qOverPtTrkGlb_phi;
  std::vector<TProfile *> qOverPtTrkGlb200_pt;
  std::vector<TProfile *> qOverPtTrkGlb200_eta;
  std::vector<TProfile *> qOverPtTrkGlb200_phi;
  std::vector<TProfile *> qOverPtGlobal_pt;
  std::vector<TProfile *> qOverPtGlobal_eta;
  std::vector<TProfile *> qOverPtGlobal_phi;


  edm::EDGetTokenT<edm::View<pat::Muon> >   theMuonCollectionLabel_;
  edm::EDGetTokenT<edm::View<reco::Vertex> >            theVertexLabel_;
  edm::EDGetTokenT<reco::BeamSpot>          theBeamSpotLabel_;
  edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>>   trigger_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;


};

#endif
