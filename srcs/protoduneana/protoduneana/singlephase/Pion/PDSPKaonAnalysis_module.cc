////////////////////////////////////////////////////////////////////////
// Class:       PDSPKaonAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        PDSPKaonAnalysis_module.cc
//
// Generated at Thu Dec  3 12:05:55 2020 by Mousam Rai,,,Mousam.Rai@warwick.ac.uk, using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "art_root_io/TFileService.h"
#include <TTree.h>
#include <vector>
#include <string>
#include "nusimdata/SimulationBase/MCParticle.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"
#include "protoduneana/Utilities/ProtoDUNEShowerUtils.h"
#include "protoduneana/Utilities/ProtoDUNEBeamlineUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"  //--------DOM'S FUNCTION-----------
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h" //for SCE correction
#include "larsim/MCCheater/ParticleInventoryService.h" //for pi_serv
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include <algorithm>
#include <tgmath.h>

namespace analysis {
  class PDSPKaonAnalysis;
}


class analysis::PDSPKaonAnalysis : public art::EDAnalyzer {
public:
  explicit PDSPKaonAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPKaonAnalysis(PDSPKaonAnalysis const&) = delete;
  PDSPKaonAnalysis(PDSPKaonAnalysis&&) = delete;
  PDSPKaonAnalysis& operator=(PDSPKaonAnalysis const&) = delete;
  PDSPKaonAnalysis& operator=(PDSPKaonAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
//  bool sortcol( const std::vector<int>& v1, const std::vector<int>& v2 );   //<=============THIS ONE, DOM
//void calculate_tier0_quatities
//void calculate_tier1_quatities
std::vector<const simb::MCParticle*> apply_true_particle_quality_cuts(std::vector<const simb::MCParticle*> inputVector);
int count_unwanted_particles(std::vector<const simb::MCParticle*> inputVector);

  TTree *fTree;

//  unsigned int fNPFParticles;
//  unsigned int fNPrimaries;
//  int fNPrimaryDaughters;
  std::string fPFParticleLabel;
  std::string fTruthLabel;
  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fGeneratorTag;
  std::string fBeamModuleLabel;
  std::string fHitTag;

  unsigned int fEventID;

  int fTrueBeamParticlePDGCode;

  float fcompleteness;
  float fpurity;

  float fRecoBeamParticleStartX;
  float fRecoBeamParticleStartY;
  float fRecoBeamParticleStartZ;

  float fBeamInst_startVertex_X_SCE_corrected;
  float fBeamInst_startVertex_Y_SCE_corrected;
  float fBeamInst_startVertex_Z_SCE_corrected;

  float fBeamInst_startVertex_dr_SCE_corrected; // need to code
  
  int fDoesTrueBeamParticleExist; // need to code
  int fIsBeamParticleReconstructed;

  float fPrimaryRecoPfpLength;
  float fPrimaryBeamParticleLength;
  float fBeam_length_by_traj_points;
  float fTrueBeamLengthVersion3;

  int fTrueBeamParticleNhits;
  int fRecoBeamParticleNhits;

  float fTrueBeamParticleEnergy;//
  int fRecoBeamParticlePDGCode;//

  int fbestMatchedMCParticleFromPFParticlePdgCode;//

  float fTrueBeamParticleStartX; //
  float fTrueBeamParticleStartY; //
  float fTrueBeamParticleStartZ; //

//  float fTrueBeamParticleStartX_SCE_corrected; // need to code
//  float fTrueBeamParticleStartY_SCE_corrected; // need to code
//  float fTrueBeamParticleStartZ_SCE_corrected; // need to code




  float fTrueBeamParticleEndX;
  float fTrueBeamParticleEndY;
  float fTrueBeamParticleEndZ;

  float fTrueBeamParticleStartPx;
  float fTrueBeamParticleStartPy;
  float fTrueBeamParticleStartPz;

  float fRecoBeamParticleStartPx;
  float fRecoBeamParticleStartPy;
  float fRecoBeamParticleStartPz;

  float fTrueBeamStartDirX;
  float fTrueBeamStartDirY;
  float fTrueBeamStartDirZ;

  float fRecoBeamStartDirX;
  float fRecoBeamStartDirY;
  float fRecoBeamStartDirZ;

  float fRecoBeamParticleInteractionX;
  float fRecoBeamParticleInteractionY;
  float fRecoBeamParticleInteractionZ;

  float fTrueBeamParticleInteractionX;
  float fTrueBeamParticleInteractionY;
  float fTrueBeamParticleInteractionZ;

  float beam_inst_X;
  float beam_inst_Y;
  float beam_inst_Z;

  int fNumberOfReconstructedBeamParticle;
  bool fMCHasBI;
  bool beam_inst_valid;

  int fReco_beam_pfp_topology;

  int fTier0MCParticleHitsSize;
  int fTier0recoBeamParticleHitsSize;
  int fSharedTier0RecoTrueHitsSize;
  int fTier0RecoMCParticleMatch;
  float fStartVertexDr;
  float fInteractionVertexDr;

  protoana::ProtoDUNETruthUtils fProtoDUNETruthUtils;
  protoana::ProtoDUNEPFParticleUtils fProtoDUNEPFParticleUtils;
  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;

//---------------------------Tier 1-----------------------------------------

  float fTier1TrueBeamParticleEnergy;

  int fTier1RecoBeamParticlePDGCode;
  int fTier1IsBeamParticleReconstructed;
  int fTier1TrueBeamParticlePDGCode;
  int fTier1bestMatchedMCParticleFromPFParticlePdgCode;
  float fTier1purity;
  float fTier1completeness;

  float fTier1TrueBeamParticleStartX;
  float fTier1TrueBeamParticleStartY;
  float fTier1TrueBeamParticleStartZ;

//  float fTrueBeamParticleStartX_SCE_corrected;
//  float fTrueBeamParticleStartY_SCE_corrected;
//  float fTrueBeamParticleStartZ_SCE_corrected;

  float fTier1BeamInst_startVertex_X_SCE_corrected;
  float fTier1BeamInst_startVertex_Y_SCE_corrected;
  float fTier1BeamInst_startVertex_Z_SCE_corrected;


  float fTier1TrueBeamParticleEndX;
  float fTier1TrueBeamParticleEndY;
  float fTier1TrueBeamParticleEndZ;

  float fTier1TrueBeamParticleStartPx;
  float fTier1TrueBeamParticleStartPy;
  float fTier1TrueBeamParticleStartPz;

  int  fTier1TrueBeamParticleNhits;

  float fTier1RecoBeamParticleStartX;
  float fTier1RecoBeamParticleStartY;
  float fTier1RecoBeamParticleStartZ;

  float fTier1RecoBeamParticleStartPx;
  float fTier1RecoBeamParticleStartPy;
  float fTier1RecoBeamParticleStartPz;

  float fTier1TrueBeamStartDirX;
  float fTier1TrueBeamStartDirY;
  float fTier1TrueBeamStartDirZ;

  float fTier1RecoBeamStartDirX;
  float fTier1RecoBeamStartDirY;
  float fTier1RecoBeamStartDirZ;

  float fTier1RecoBeamParticleInteractionX;
  float fTier1RecoBeamParticleInteractionY;
  float fTier1RecoBeamParticleInteractionZ;

  float fTier1TrueBeamParticleInteractionX;
  float fTier1TrueBeamParticleInteractionY;
  float fTier1TrueBeamParticleInteractionZ;

  int fTier1RecoBeamParticleNhits;

  float tier1beam_inst_X;
  float tier1beam_inst_Y;
  float tier1beam_inst_Z;

/*  int fNumberOfReconstructedBeamParticle;
  bool fMCHasBI;
  bool beam_inst_valid;*/

  float fTier1PrimaryRecoPfpLength;
  float fTier1PrimaryBeamParticleLength;
  float fTier1Beam_length_by_traj_points;
  int fTier1Reco_beam_pfp_topology;

  int fTier1MCParticleHitsSize;
  int fTier1recoBeamParticleHitsSize;
  int fSharedTier1RecoTrueHitsSize;
  int fTier1RecoMCParticleMatch;

  int fTier1MCParticleNumber;
  int fTier1RecoParticleNumber;
  int fTier1MCRecoMatchedNumber;
  float fTier1TrueBeamParticleLength;
  float fTier1StartVertexDr;

//-------------------------------------------------------------------
//  protoana::ProtoDUNEBeamlineUtils fBeamlineUtils;
};


analysis::PDSPKaonAnalysis::PDSPKaonAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  //,

//    fBeamlineUtils(p.get< fhicl::ParameterSet >("BeamlineUtils"))

  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
    fPFParticleLabel    = p.get<std::string>("PFParticleLabel");
    fTruthLabel         = p.get<std::string>("TruthLabel");
    fTrackLabel         = p.get<std::string>("TrackLabel");
    fShowerLabel        = p.get<std::string>("ShowerLabel");
    fGeneratorTag       = p.get<std::string>("GeneratorTag");
    fHitTag             = p.get<std::string>("HitTag");
    fBeamModuleLabel    = p.get<std::string>("BeamModuleLabel");
}

void analysis::PDSPKaonAnalysis::analyze(art::Event const& e)
{
    std::cout << "====================== Starting Kaon Analyser ======================" << std::endl;

  // Implementation of required member function here.

    fEventID = e.id().event();
//    fTruePrimaryEnergy = -999.;
    fTrueBeamParticleEnergy                     = -9999.;
    fTrueBeamParticlePDGCode                    = -9999;
    fpurity                                     = -1.f;
    fcompleteness                               = -1.f;

    fTrueBeamParticleStartX                     = -9999;
    fTrueBeamParticleStartY                     = -9999;
    fTrueBeamParticleStartZ                     = -9999;

//    fTrueBeamParticleStartX_SCE_corrected       = -9999;
//    fTrueBeamParticleStartY_SCE_corrected       = -9999;
//    fTrueBeamParticleStartZ_SCE_corrected       = -9999;

    fBeamInst_startVertex_X_SCE_corrected       = -9999;
    fBeamInst_startVertex_Y_SCE_corrected       = -9999;
    fBeamInst_startVertex_Z_SCE_corrected       = -9999;

    fTrueBeamParticleEndX                       = -9999;
    fTrueBeamParticleEndY                       = -9999;
    fTrueBeamParticleEndZ                       = -9999;

    beam_inst_X                                 = -9999;
    beam_inst_Y                                 = -9999;
    beam_inst_Z                                 = -9999;
    beam_inst_valid                             =  true;

    fTrueBeamParticleStartPx                    = -9999;
    fTrueBeamParticleStartPy                    = -9999;
    fTrueBeamParticleStartPz                    = -9999;

    fTrueBeamParticleNhits                      = -9999;
    fRecoBeamParticleNhits                      = -9999;

    fRecoBeamParticleStartX                     = -9999;
    fRecoBeamParticleStartY                     = -9999;
    fRecoBeamParticleStartZ                     = -9999;

    fRecoBeamParticleStartPx                    = -9999;
    fRecoBeamParticleStartPy                    = -9999;
    fRecoBeamParticleStartPz                    = -9999;

    fTrueBeamStartDirX                          = -9999;
    fTrueBeamStartDirY                          = -9999;
    fTrueBeamStartDirZ                          = -9999;

    fRecoBeamStartDirX                          = -9999;
    fRecoBeamStartDirY                          = -9999;
    fRecoBeamStartDirZ                          = -9999;
            
    fRecoBeamParticleInteractionX               = -9999;
    fRecoBeamParticleInteractionY               = -9999;
    fRecoBeamParticleInteractionZ               = -9999;

    fTrueBeamParticleInteractionX               = -9999;
    fTrueBeamParticleInteractionY               = -9999;
    fTrueBeamParticleInteractionZ               = -9999;

    fNumberOfReconstructedBeamParticle          = 1;

    double tpcActiveXLow                        = -358.5282287595;
    double tpcActiveXHigh                       = 358.5282287595;
    double tpcActiveYLow                        = 0;
    double tpcActiveYHigh                       = 603.8612670895;
    double tpcActiveZLow                        = 0;
    double tpcActiveZHigh                       = 695.2862548825;

    fPrimaryRecoPfpLength                       = -9999;
    fPrimaryBeamParticleLength                  = -9999;
    fBeam_length_by_traj_points                 = -9999;

    fReco_beam_pfp_topology                     = -9999;

    fTier0MCParticleHitsSize                    = -9999;
    fTier0recoBeamParticleHitsSize              = -9999;
    fSharedTier0RecoTrueHitsSize                = -9999;
    fTier0RecoMCParticleMatch                   = 0;
    fDoesTrueBeamParticleExist                  = -9999;
    fStartVertexDr                              = -9999;
    fInteractionVertexDr                        = -9999;
//------------------------------Tier 1--------------------------------------
    fTier1TrueBeamParticleEnergy                     = -9999.;
    fTier1TrueBeamParticlePDGCode                    = -9999;
    fTier1purity                                     = -1.f;
    fTier1completeness                               = -1.f;

    fTier1TrueBeamParticleStartX                     = -9999;
    fTier1TrueBeamParticleStartY                     = -9999;
    fTier1TrueBeamParticleStartZ                     = -9999;

//    fTrueBeamParticleStartX_SCE_corrected       = -9999;
//    fTrueBeamParticleStartY_SCE_corrected       = -9999;
//    fTrueBeamParticleStartZ_SCE_corrected       = -9999;

    fTier1BeamInst_startVertex_X_SCE_corrected       = -9999;
    fTier1BeamInst_startVertex_Y_SCE_corrected       = -9999;
    fTier1BeamInst_startVertex_Z_SCE_corrected       = -9999;

    fTier1TrueBeamParticleEndX                       = -9999;
    fTier1TrueBeamParticleEndY                       = -9999;
    fTier1TrueBeamParticleEndZ                       = -9999;

    tier1beam_inst_X                                 = -9999;
    tier1beam_inst_Y                                 = -9999;
    tier1beam_inst_Z                                 = -9999;
//    beam_inst_valid                             =  true;

    fTier1TrueBeamParticleStartPx                    = -9999;
    fTier1TrueBeamParticleStartPy                    = -9999;
    fTier1TrueBeamParticleStartPz                    = -9999;

    fTier1TrueBeamParticleNhits                      = -9999;
    fTier1RecoBeamParticleNhits                      = -9999;

    fTier1RecoBeamParticleStartX                     = -9999;
    fTier1RecoBeamParticleStartY                     = -9999;
    fTier1RecoBeamParticleStartZ                     = -9999;

    fTier1RecoBeamParticleStartPx                    = -9999;
    fTier1RecoBeamParticleStartPy                    = -9999;
    fTier1RecoBeamParticleStartPz                    = -9999;

    fTier1TrueBeamStartDirX                          = -9999;
    fTier1TrueBeamStartDirY                          = -9999;
    fTier1TrueBeamStartDirZ                          = -9999;

    fTier1RecoBeamStartDirX                          = -9999;
    fTier1RecoBeamStartDirY                          = -9999;
    fTier1RecoBeamStartDirZ                          = -9999;
            
    fTier1RecoBeamParticleInteractionX               = -9999;
    fTier1RecoBeamParticleInteractionY               = -9999;
    fTier1RecoBeamParticleInteractionZ               = -9999;

    fTier1TrueBeamParticleInteractionX               = -9999;
    fTier1TrueBeamParticleInteractionY               = -9999;
    fTier1TrueBeamParticleInteractionZ               = -9999;

//    fTier1NumberOfReconstructedBeamParticle          = 1;

    fTier1PrimaryRecoPfpLength                       = -9999;
    fTier1PrimaryBeamParticleLength                  = -9999;
    fTier1Beam_length_by_traj_points                 = -9999;

    fTier1Reco_beam_pfp_topology                     = -9999;

    fTier1MCParticleHitsSize                    = -9999;
    fTier1recoBeamParticleHitsSize              = -9999;
    fSharedTier1RecoTrueHitsSize                = -9999;
    fTier1RecoMCParticleMatch                   = 0;

    fTier1MCParticleNumber                            = -9999;
    fTier1RecoParticleNumber                          = -9999; 
    fTier1MCRecoMatchedNumber                   = -9999;

    fTier1TrueBeamParticleLength                = -9999;
    fTier1StartVertexDr                         = -9999;


//------------------------------------------------------------------

//==================Accessing Truth Info Block=============================

/*========To do list=======================================================

    rename beam_inst_vertex into start_vertex_XYZ
=========================================================================*/

//For tier0, count associated hits with MCParticle

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    const simb::MCParticle* true_beam_particle = 0x0;

    std::vector< const recob::Hit* > tier0MCParticleHits;

    std::vector < const simb::MCParticle*> tier1TrueParticles;
    std::vector <const recob::PFParticle* > tier1RecoParticles;
    std::vector < const simb::MCParticle*> tier1TrueParticlesQ1;

    const sim::ParticleList & plist = pi_serv->ParticleList();

    if(!e.isRealData())
    {
        auto mcTruths = e.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
        
        true_beam_particle = fProtoDUNETruthUtils.GetGeantGoodParticle((*mcTruths)[0],e);
        fDoesTrueBeamParticleExist = 1;
        tier0MCParticleHits = fProtoDUNETruthUtils.GetMCParticleHits(clockData, *true_beam_particle, e, "hitpdune" );
    }
    
//    std::cout << "tier0MCParticleHits size: " << tier0MCParticleHits.size() << std::endl;
    fTier0MCParticleHitsSize = tier0MCParticleHits.size();

//For tier0, count associated hits with reco beam Particle

    art::ValidHandle<std::vector<recob::PFParticle>> recoParticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    std::vector<const recob::PFParticle*> beamParticles = fProtoDUNEPFParticleUtils.GetPFParticlesFromBeamSlice(e,fPFParticleLabel);
//    std::cout << "beamParticles.size(): " << beamParticles.size() << std::endl;

    std::vector< const recob::Hit* > tier0recoBeamParticleHits;

    if(beamParticles.size() == 0)
    {
        std::cerr << "We found no beam particles for this event... moving on" << std::endl;
        return;
    }

    if(beamParticles.size() == 1)
    {
        fIsBeamParticleReconstructed = 1;

        for(const recob::PFParticle* particle : beamParticles)
        {
            tier0recoBeamParticleHits = fProtoDUNEPFParticleUtils.GetPFParticleHits(*particle, e, fPFParticleLabel);
//            std::cout << "tier0recoBeamParticleHits size: " << tier0recoBeamParticleHits.size() << std::endl;
        }
    }

    else
        std::cout << "More than one reconstructed beam particle!" << std::endl;

    fTier0recoBeamParticleHitsSize = tier0recoBeamParticleHits.size();



//SharedHits block    


    std::vector< const recob::Hit* > sharedTier0RecoTrueHits = fProtoDUNETruthUtils.GetSharedHits(clockData, *true_beam_particle, *beamParticles[0], e, fPFParticleLabel);
//    std::cout << "sharedTier0RecoTrueHits size: " << sharedTier0RecoTrueHits.size() << std::endl;

    fSharedTier0RecoTrueHitsSize = sharedTier0RecoTrueHits.size();

//Reconstruction efficiency criteria

    if( ((tier0MCParticleHits.size()) != 0) && ((tier0recoBeamParticleHits.size()) != 0) && (((static_cast<float>(sharedTier0RecoTrueHits.size()))/(static_cast<float>(tier0MCParticleHits.size()))) >= 0.5) && (((static_cast<float>(sharedTier0RecoTrueHits.size()))/(static_cast<float>(tier0recoBeamParticleHits.size()))) >= 0.5) )
        fTier0RecoMCParticleMatch = 1;

//calculating quatities
    
    if(fTier0RecoMCParticleMatch == 1)
    {
        //Truth quatities
        fTrueBeamParticleEnergy     = true_beam_particle->E();
        fTrueBeamParticlePDGCode    = true_beam_particle->PdgCode();
        fTrueBeamParticleStartX     = true_beam_particle->Position(0).X();
        fTrueBeamParticleStartY     = true_beam_particle->Position(0).Y();
        fTrueBeamParticleStartZ     = true_beam_particle->Position(0).Z();
        fTrueBeamParticleEndX       = true_beam_particle->EndX();
        fTrueBeamParticleEndY       = true_beam_particle->EndY();
        fTrueBeamParticleEndZ       = true_beam_particle->EndZ();
        fTrueBeamParticleStartPx    = true_beam_particle->Px();
        fTrueBeamParticleStartPy    = true_beam_particle->Py();
        fTrueBeamParticleStartPz    = true_beam_particle->Pz();
        fTrueBeamStartDirX          = (fTrueBeamParticleStartPx)/(true_beam_particle->P());
        fTrueBeamStartDirY          = (fTrueBeamParticleStartPy)/(true_beam_particle->P());
        fTrueBeamStartDirZ          = (fTrueBeamParticleStartPz)/(true_beam_particle->P());
        fTrueBeamLengthVersion3     = true_beam_particle->Trajectory().TotalLength();
        std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;
        auto beamHandle = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("generator");

        if (beamHandle.isValid())
            art::fill_ptr_vector(beamVec, beamHandle);

        else
            std::cout << "invalid beam handle" << std::endl;

        const beam::ProtoDUNEBeamEvent beamEvent = *(beamVec.at(0));

        auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
        int nTracks = beamEvent.GetBeamTracks().size();

        if( nTracks > 0 )
        {
            beam_inst_X = beamEvent.GetBeamTracks()[0].Trajectory().End().X();
            beam_inst_Y = beamEvent.GetBeamTracks()[0].Trajectory().End().Y();
            beam_inst_Z = beamEvent.GetBeamTracks()[0].Trajectory().End().Z();

            fBeamInst_startVertex_X_SCE_corrected = beam_inst_X-SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).X();
            fBeamInst_startVertex_Y_SCE_corrected = beam_inst_Y+SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).Y();
            fBeamInst_startVertex_Z_SCE_corrected = beam_inst_Z+SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).Z();
        }

        float true_beam_particle_SCE_end_X = fTrueBeamParticleEndX-SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleEndX,fTrueBeamParticleEndY,fTrueBeamParticleEndZ)).X();
        float true_beam_particle_SCE_end_Y = fTrueBeamParticleEndY+SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleEndX,fTrueBeamParticleEndY,fTrueBeamParticleEndZ)).Y();
        float true_beam_particle_SCE_end_Z = fTrueBeamParticleEndZ+SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleEndX,fTrueBeamParticleEndY,fTrueBeamParticleEndZ)).Z();

        fTrueBeamParticleInteractionX = true_beam_particle_SCE_end_X;       
        fTrueBeamParticleInteractionY = true_beam_particle_SCE_end_Y;       
        fTrueBeamParticleInteractionZ = true_beam_particle_SCE_end_Z;       

        fBeam_length_by_traj_points = sqrt(((true_beam_particle_SCE_end_X-fBeamInst_startVertex_X_SCE_corrected)*(true_beam_particle_SCE_end_X-fBeamInst_startVertex_X_SCE_corrected)) + ((true_beam_particle_SCE_end_Y-fBeamInst_startVertex_Y_SCE_corrected)*(true_beam_particle_SCE_end_Y-fBeamInst_startVertex_Y_SCE_corrected)) + ((true_beam_particle_SCE_end_Z-fBeamInst_startVertex_Z_SCE_corrected)*(true_beam_particle_SCE_end_Z-fBeamInst_startVertex_Z_SCE_corrected)));

        fTrueBeamParticleNhits = fProtoDUNETruthUtils.GetMCParticleHits( clockData, *true_beam_particle, e, fHitTag ).size();

        fPrimaryBeamParticleLength = fProtoDUNETruthUtils.GetMCParticleLengthInTPCActiveVolume(*true_beam_particle, tpcActiveXLow, tpcActiveXHigh, tpcActiveYLow, tpcActiveYHigh, tpcActiveZLow, tpcActiveZHigh);

        //Reco Quantities
        
        fRecoBeamParticleNhits = fProtoDUNEPFParticleUtils.GetPFParticleHits( *beamParticles[0], e, fPFParticleLabel).size();

        const recob::Track* thisTrack = fProtoDUNEPFParticleUtils.GetPFParticleTrack(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);
        const recob::Shower* thisShower = fProtoDUNEPFParticleUtils.GetPFParticleShower(*beamParticles[0],e,fPFParticleLabel,fShowerLabel);

        if(thisTrack != 0x0)
        {
            fReco_beam_pfp_topology = 1;
            fPrimaryRecoPfpLength = thisTrack->Length();

            const TVector3 reco_primary_start_vertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);

            fRecoBeamParticleStartX = reco_primary_start_vertex.X();
            fRecoBeamParticleStartY = reco_primary_start_vertex.Y();
            fRecoBeamParticleStartZ = reco_primary_start_vertex.Z();

            const TVector3 reco_primary_interaction_vertex = fProtoDUNEPFParticleUtils.GetPFParticleSecondaryVertex(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);

            fRecoBeamParticleInteractionX = reco_primary_interaction_vertex.X();
            fRecoBeamParticleInteractionY = reco_primary_interaction_vertex.Y();
            fRecoBeamParticleInteractionZ = reco_primary_interaction_vertex.Z();

            fRecoBeamStartDirX = (thisTrack->StartDirection()).X();
            fRecoBeamStartDirY = (thisTrack->StartDirection()).Y();
            fRecoBeamStartDirZ = (thisTrack->StartDirection()).Z();

            float interactionVertexDx = fRecoBeamParticleInteractionX - fTrueBeamParticleInteractionX;
            float interactionVertexDy = fRecoBeamParticleInteractionY - fTrueBeamParticleInteractionY;
            float interactionVertexDz = fRecoBeamParticleInteractionZ - fTrueBeamParticleInteractionZ;

            fInteractionVertexDr = sqrt(((interactionVertexDx)*(interactionVertexDx))+((interactionVertexDy)*(interactionVertexDy))+((interactionVertexDz)*(interactionVertexDz)));
        }

        if(thisShower != 0x0) 
        {
            fReco_beam_pfp_topology = 0;
            fPrimaryRecoPfpLength = thisShower->Length();

            const TVector3 reco_primary_start_vertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);

            fRecoBeamParticleStartX = reco_primary_start_vertex.X();
            fRecoBeamParticleStartY = reco_primary_start_vertex.Y();
            fRecoBeamParticleStartZ = reco_primary_start_vertex.Z();

            fRecoBeamStartDirX = (thisShower->Direction()).X();
            fRecoBeamStartDirY = (thisShower->Direction()).Y();
            fRecoBeamStartDirZ = (thisShower->Direction()).Z();

            fInteractionVertexDr = 0;
        }

        const simb::MCParticle* bestMatchedMcParticleFromPFParticle;
        bestMatchedMcParticleFromPFParticle = (fProtoDUNETruthUtils.GetMCParticleFromPFParticle(clockData, *beamParticles[0], e, fPFParticleLabel));
        fbestMatchedMCParticleFromPFParticlePdgCode = bestMatchedMcParticleFromPFParticle->PdgCode();
    
        fpurity = fProtoDUNETruthUtils.GetPurity(clockData, *beamParticles[0], e, fPFParticleLabel); 
//        std::cout << "purity: " << fpurity << std::endl; 

        fcompleteness = fProtoDUNETruthUtils.GetCompleteness(clockData, *beamParticles[0], e, "pandora", "hitpdune");
//        std::cout << "completeness: " << fcompleteness << std::endl;

        float startVertexDx = fRecoBeamParticleStartX - fBeamInst_startVertex_X_SCE_corrected;
        float startVertexDy = fRecoBeamParticleStartY - fBeamInst_startVertex_Y_SCE_corrected;
        float startVertexDz = fRecoBeamParticleStartZ - fBeamInst_startVertex_Z_SCE_corrected;

        fStartVertexDr = sqrt(((startVertexDx)*(startVertexDx))+((startVertexDy)*(startVertexDy))+((startVertexDz)*(startVertexDz)));

    }
    fTree->Fill();

    if (fTier0RecoMCParticleMatch == 1)
    {
        //For Tier 1 true particles

        int number_true_tier1_daughter = true_beam_particle->NumberDaughters();

        for (int i = 0; i < number_true_tier1_daughter; i++)
        {   
            auto part = plist[ true_beam_particle->Daughter(i) ]; 
//            std::cout << "particle pdgcode: " << part->PdgCode() << std::endl;
//            const std::vector < const recob::Hit*> mcHitsVector = fProtoDUNETruthUtils.GetMCParticleHits(clockData, *part, e, fHitTag);
//            int nMCHits = mcHitsVector.size();
//            std::cout << "nMCHits: " << nMCHits << std::endl;
      
            tier1TrueParticles.push_back(part); 
        }

/*        int test = tier1TrueParticles.size();
        std::cout << "tier1TrueParticles size: " << test << std::endl;
//------------------------------------------------------------------------------------------------------------------
        for (int i = 0; i < test; i++)
        {
            
            if ( ((tier1TrueParticles[i]->PdgCode()) != 111) && ((tier1TrueParticles[i]->PdgCode()) != 2112) && ((tier1TrueParticles[i]->PdgCode()) < 1000000000))
                tier1TrueParticlesQ1.push_back(tier1TrueParticles[i]);
        
            else
            {
                std::cout << "------------------granddaughters------------------------" << std::endl;
                int number_true_granddaughter = tier1TrueParticles[i]->NumberDaughters();
                std::cout << "number_true_granddaughter: " << number_true_granddaughter << std::endl;
                for (int j = 0; j < number_true_granddaughter; j++)
                {
                    auto granddaughter = plist[ tier1TrueParticles[i]->Daughter(j)];
                    tier1TrueParticlesQ1.push_back(granddaughter);
                }
                std::cout << "--------------------------------------------------------" << std::endl;
//                tier1TrueParticlesQ1.erase(i);
            }
        }

        std::cout << "tier1TrueParticlesQ1 size: " << tier1TrueParticlesQ1.size() << std::endl;
        

        test2 = apply_true_particle_quality_cuts(tier1TrueParticles);*/
        tier1TrueParticlesQ1 = apply_true_particle_quality_cuts(tier1TrueParticles);
//        std::cout << "tier1TrueParticlesQ1: " << tier1TrueParticlesQ1.size() << std::endl;

        int numUnwantedParticles = count_unwanted_particles(tier1TrueParticles);

        while (numUnwantedParticles != 0)
        {
            tier1TrueParticlesQ1 = apply_true_particle_quality_cuts(tier1TrueParticles);
            numUnwantedParticles = count_unwanted_particles(tier1TrueParticlesQ1);
            
            if (numUnwantedParticles != 0)
            {
                tier1TrueParticles = tier1TrueParticlesQ1;
                tier1TrueParticlesQ1.empty();
            }
        }

//        std::cout << "tier1TrueParticlesQ1.size(): " << tier1TrueParticlesQ1.size() << std::endl;
//        int new_size = tier1TrueParticlesQ1.size();
/*        for (int i = 0; i < new_size; i++)
        {
            int pdgcode = tier1TrueParticlesQ1[i]->PdgCode();
            const std::vector < const recob::Hit*> mcHitsVector = fProtoDUNETruthUtils.GetMCParticleHits(clockData, *tier1TrueParticlesQ1[i], e, fHitTag);
            int nMCHits = mcHitsVector.size();
            std::cout << "pdgcode: " << pdgcode << std::endl;
            std::cout << "nMCHits: " << nMCHits << std::endl;
        }*/
//------------------------------------------------------------------------------------------------------------------


        //For tier 1 reco particles

        for(const recob::PFParticle* particle : beamParticles)
        {
            for(const int daughterID : particle->Daughters())
            {
                const recob::PFParticle *daughterParticle = &(recoParticles->at(daughterID));
                tier1RecoParticles.push_back(daughterParticle);
            }
        }

        int array_row_size = tier1TrueParticlesQ1.size();
//        std::cout << "array_row_size: " << array_row_size << std::endl;
        int array_col_size = tier1RecoParticles.size();
        fTier1RecoParticleNumber = array_col_size;
         
//        std::cout << "array_col_size: " << array_col_size << std::endl;

        std::vector< std::vector<int> > twoDTrueRecoSharedHitsMatrix;
        
        for (int i = 0; i < array_row_size; i++)
        {
//            std::vector<const recob::Hit*> trueHits = fProtoDUNETruthUtils.GetMCParticleHits(clockData, *tier1TrueParticles[i], e, fHitTag)

            for (int j = 0; j < array_col_size; j++)
            {
                std::vector<int> v1;
//                std::vector<const recob::Hit*> recoHits = fProtoDUNEPFParticleUtils.GetPFParticleHits(*tier1RecoParticles[j], e, fPFParticleLabel);
                std::vector< const recob::Hit* > hit_matched_vector = fProtoDUNETruthUtils.GetSharedHits(clockData, *tier1TrueParticlesQ1[i], *tier1RecoParticles[j], e, fPFParticleLabel);

                int hit_matched = hit_matched_vector.size();
                if (hit_matched != 0)
                {
                    v1 = {i, j, hit_matched};
                    twoDTrueRecoSharedHitsMatrix.push_back(v1);
                }            
            }
            
        }

        int row = twoDTrueRecoSharedHitsMatrix.size();
        fTier1MCParticleNumber = row;
//        int col = twoDTrueRecoSharedHitsMatrix[0].size();
        std::vector< std::vector<int> > twoDTrueRecoSharedHitsSortedMatchedMatrix;

//       for (int i = 0; i < row; i++)
//            std::cout << "shared hits: " << twoDTrueRecoSharedHitsMatrix[i][0] << " , " << twoDTrueRecoSharedHitsMatrix[i][1] << " , " << twoDTrueRecoSharedHitsMatrix[i][2] << std::endl;

        std::sort(twoDTrueRecoSharedHitsMatrix.begin(), twoDTrueRecoSharedHitsMatrix.end(), [] ( const std::vector<int>& v1, const std::vector<int>& v2 )->bool
        {
            return v1[2] > v2[2];
        }); //<========twoDTrueRecoSharedHitsMatrix is the vector<vector<int>>

//        for (int i = 0; i < row; i++)
//            std::cout << "SORTED shared hits: " << twoDTrueRecoSharedHitsMatrix[i][0] << " , " << twoDTrueRecoSharedHitsMatrix[i][1] << " , " << twoDTrueRecoSharedHitsMatrix[i][2] << std::endl;

//        for (int i = 0; i < row; i++)

/*        for (int i = 0; i < row; i++)
        {
            if ( (twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0] != twoDTrueRecoSharedHitsMatrix[i][0]) && (twoDTrueRecoSharedHitsSortedMatchedMatrix[i] != twoDTrueRecoSharedHitsMatrix[i]))
                twoDTrueRecoSharedHitsSortedMatchedMatrix.push_back(twoDTrueRecoSharedHitsMatrix[i])
        }
            
        for (int i = 0; i < row; i++)
            std::cout << "SORTED matched shared hits: " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0] << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][1] << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][2] << std::endl;*/

        std::vector<int> taken_true_particle_index;
        std::vector<int> taken_reco_particle_index;

        for (int i = 0; i < row; i++)
        {
            
            int true_count = std::count (taken_true_particle_index.begin(), taken_true_particle_index.end(), twoDTrueRecoSharedHitsMatrix[i][0]);
            int reco_count = std::count (taken_reco_particle_index.begin(), taken_reco_particle_index.end(), twoDTrueRecoSharedHitsMatrix[i][1]);

            if ( (true_count == 0) && (reco_count == 0) )
            {
                taken_true_particle_index.push_back(twoDTrueRecoSharedHitsMatrix[i][0]);
                taken_reco_particle_index.push_back(twoDTrueRecoSharedHitsMatrix[i][1]);
                twoDTrueRecoSharedHitsSortedMatchedMatrix.push_back(twoDTrueRecoSharedHitsMatrix[i]);
            }
//            std::cout << "taken_true_particle_index size: " << taken_true_particle_index.size() << std::endl;
//            std::cout << "taken_reco_particle_index size: " << taken_reco_particle_index.size() << std::endl;

//            std::cout << "true_count: " << true_count << std::endl;
//            std::cout << "reco_count: " << reco_count << std::endl;

        }

//        for (int i = 0; i < row; i++)
//            std::cout << "SORTED matched shared hits: " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][1] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][2] << std::endl;

//                fpurity = fProtoDUNETruthUtils.GetPurity(clockData, *beamParticles[0], e, fPFParticleLabel); 
//        fcompleteness = fProtoDUNETruthUtils.GetCompleteness(clockData, *beamParticles[0], e, "pandora", "hitpdune");
        int matched_row = twoDTrueRecoSharedHitsSortedMatchedMatrix.size();
        fTier1MCRecoMatchedNumber = matched_row;
//        std::cout << "matched_row size: " << matched_row << std::endl;
        for (int i = 0; i < matched_row; i++)
        {
            auto trueParticle = tier1TrueParticlesQ1[twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0]];
            auto recoParticle = tier1RecoParticles[twoDTrueRecoSharedHitsSortedMatchedMatrix[i][1]];

            fTier1TrueBeamParticleNhits = fProtoDUNETruthUtils.GetMCParticleHits( clockData, *trueParticle, e, fHitTag ).size();
            fTier1RecoBeamParticleNhits = fProtoDUNEPFParticleUtils.GetPFParticleHits( *recoParticle, e, fPFParticleLabel).size();
            auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
//--------------------------------true tier 1 quantities--------------------------------------------------------------------------
            fTier1TrueBeamParticleLength     = trueParticle->Trajectory().TotalLength();
            fTier1TrueBeamParticleEndX       = trueParticle->EndX();
            fTier1TrueBeamParticleEndY       = trueParticle->EndY();
            fTier1TrueBeamParticleEndZ       = trueParticle->EndZ();
            fTier1TrueBeamParticlePDGCode    = trueParticle->PdgCode();
            fTier1TrueBeamParticleStartX     = trueParticle->Position(0).X();
            fTier1TrueBeamParticleStartY     = trueParticle->Position(0).Y();
            fTier1TrueBeamParticleStartZ     = trueParticle->Position(0).Z();

//            tier1beam_inst_X = trueParticle->Trajectory().X(trueParticle->Trajectory().size());
//            tier1beam_inst_Y = trueParticle->Trajectory().Y(trueParticle->Trajectory().size());
//            tier1beam_inst_Z = trueParticle->Trajectory().Z(trueParticle->Trajectory().size());

            tier1beam_inst_X = fTier1TrueBeamParticleStartX;
            tier1beam_inst_Y = fTier1TrueBeamParticleStartY;
            tier1beam_inst_Z = fTier1TrueBeamParticleStartZ;

            fTier1BeamInst_startVertex_X_SCE_corrected = tier1beam_inst_X-SCE->GetPosOffsets(geo::Point_t(tier1beam_inst_X,tier1beam_inst_Y,tier1beam_inst_Z)).X();
            fTier1BeamInst_startVertex_Y_SCE_corrected = tier1beam_inst_Y+SCE->GetPosOffsets(geo::Point_t(tier1beam_inst_X,tier1beam_inst_Y,tier1beam_inst_Z)).Y();
            fTier1BeamInst_startVertex_Z_SCE_corrected = tier1beam_inst_Z+SCE->GetPosOffsets(geo::Point_t(tier1beam_inst_X,tier1beam_inst_Y,tier1beam_inst_Z)).Z();
        

            float tier1true_beam_particle_SCE_end_X = fTier1TrueBeamParticleEndX-SCE->GetPosOffsets(geo::Point_t(fTier1TrueBeamParticleEndX,fTier1TrueBeamParticleEndY,fTier1TrueBeamParticleEndZ)).X();
            float tier1true_beam_particle_SCE_end_Y = fTier1TrueBeamParticleEndY+SCE->GetPosOffsets(geo::Point_t(fTier1TrueBeamParticleEndX,fTier1TrueBeamParticleEndY,fTier1TrueBeamParticleEndZ)).Y();
            float tier1true_beam_particle_SCE_end_Z = fTier1TrueBeamParticleEndZ+SCE->GetPosOffsets(geo::Point_t(fTier1TrueBeamParticleEndX,fTier1TrueBeamParticleEndY,fTier1TrueBeamParticleEndZ)).Z();

            fTier1TrueBeamParticleInteractionX = tier1true_beam_particle_SCE_end_X;       
            fTier1TrueBeamParticleInteractionY = tier1true_beam_particle_SCE_end_Y;       
            fTier1TrueBeamParticleInteractionZ = tier1true_beam_particle_SCE_end_Z;       

            fTier1Beam_length_by_traj_points = sqrt(((tier1true_beam_particle_SCE_end_X-fTier1BeamInst_startVertex_X_SCE_corrected)*(tier1true_beam_particle_SCE_end_X-fTier1BeamInst_startVertex_X_SCE_corrected)) + ((tier1true_beam_particle_SCE_end_Y-fTier1BeamInst_startVertex_Y_SCE_corrected)*(tier1true_beam_particle_SCE_end_Y-fTier1BeamInst_startVertex_Y_SCE_corrected)) + ((tier1true_beam_particle_SCE_end_Z-fTier1BeamInst_startVertex_Z_SCE_corrected)*(tier1true_beam_particle_SCE_end_Z-fTier1BeamInst_startVertex_Z_SCE_corrected)));


//--------------------------------reco tier 1 quantities--------------------------------------------------------------------------
            fTier1purity = fProtoDUNETruthUtils.GetPurity(clockData, *recoParticle, e, fPFParticleLabel);
            fTier1completeness = fProtoDUNETruthUtils.GetCompleteness(clockData, *recoParticle, e, "pandora", "hitpdune");
            
            const recob::Track* thisTrack = fProtoDUNEPFParticleUtils.GetPFParticleTrack(*recoParticle,e,fPFParticleLabel,fTrackLabel); 
            const recob::Shower* thisShower = fProtoDUNEPFParticleUtils.GetPFParticleShower(*recoParticle,e,fPFParticleLabel,fShowerLabel);

            if(thisTrack != 0x0)
            {
                fTier1Reco_beam_pfp_topology = 1;
                fTier1PrimaryRecoPfpLength = thisTrack->Length();

                const TVector3 reco_primary_start_vertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(*recoParticle,e,fPFParticleLabel,fTrackLabel);

                fTier1RecoBeamParticleStartX = reco_primary_start_vertex.X();
                fTier1RecoBeamParticleStartY = reco_primary_start_vertex.Y();
                fTier1RecoBeamParticleStartZ = reco_primary_start_vertex.Z();

                const TVector3 reco_primary_interaction_vertex = fProtoDUNEPFParticleUtils.GetPFParticleSecondaryVertex(*recoParticle,e,fPFParticleLabel,fTrackLabel);

                fTier1RecoBeamParticleInteractionX = reco_primary_interaction_vertex.X();
                fTier1RecoBeamParticleInteractionY = reco_primary_interaction_vertex.Y();
                fTier1RecoBeamParticleInteractionZ = reco_primary_interaction_vertex.Z();

                fTier1RecoBeamStartDirX = (thisTrack->StartDirection()).X();
                fTier1RecoBeamStartDirY = (thisTrack->StartDirection()).Y();
                fTier1RecoBeamStartDirZ = (thisTrack->StartDirection()).Z();
            }

            if(thisShower != 0x0) 
            {
                fTier1Reco_beam_pfp_topology = 0;
                fTier1PrimaryRecoPfpLength = thisShower->Length();

                const TVector3 reco_primary_start_vertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);

                fTier1RecoBeamParticleStartX = reco_primary_start_vertex.X();
                fTier1RecoBeamParticleStartY = reco_primary_start_vertex.Y();
                fTier1RecoBeamParticleStartZ = reco_primary_start_vertex.Z();

                fTier1RecoBeamStartDirX = (thisShower->Direction()).X();
                fTier1RecoBeamStartDirY = (thisShower->Direction()).Y();
                fTier1RecoBeamStartDirZ = (thisShower->Direction()).Z();
            }

        float tier1StartVertexDx = fTier1RecoBeamParticleStartX - fTier1BeamInst_startVertex_X_SCE_corrected;
        float tier1StartVertexDy = fTier1RecoBeamParticleStartY - fTier1BeamInst_startVertex_Y_SCE_corrected;
        float tier1StartVertexDz = fTier1RecoBeamParticleStartZ - fTier1BeamInst_startVertex_Z_SCE_corrected;

        fTier1StartVertexDr = sqrt(((tier1StartVertexDx)*(tier1StartVertexDx))+((tier1StartVertexDy)*(tier1StartVertexDy))+((tier1StartVertexDz)*(tier1StartVertexDz)));
            fTree->Fill();
        }
        
    }
}

void analysis::PDSPKaonAnalysis::beginJob()
{
  // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("tree","Analyser Output Tree");

  // Truth Tree Branches    
    fTree->Branch("eventID",&fEventID,"eventID/i");

    fTree->Branch("trueBeamParticleEnergy",&fTrueBeamParticleEnergy,"trueBeamParticleEnergy/F");
    fTree->Branch("trueBeamParticlePDGCode",&fTrueBeamParticlePDGCode,"trueBeamParticlePDGCode/I");

    fTree->Branch("trueBeamParticleStartX",&fTrueBeamParticleStartX,"trueBeamParticleStartX/F");
    fTree->Branch("trueBeamParticleStartY",&fTrueBeamParticleStartY,"trueBeamParticleStartY/F");
    fTree->Branch("trueBeamParticleStartZ",&fTrueBeamParticleStartZ,"trueBeamParticleStartZ/F");

    fTree->Branch("trueBeamParticleEndX",&fTrueBeamParticleEndX,"trueBeamParticleEndX/F");
    fTree->Branch("trueBeamParticleEndY",&fTrueBeamParticleEndY,"trueBeamParticleEndY/F");
    fTree->Branch("trueBeamParticleEndZ",&fTrueBeamParticleEndZ,"trueBeamParticleEndZ/F");


    fTree->Branch("beam_inst_X",&beam_inst_X,"beam_inst_X/F");
    fTree->Branch("beam_inst_Y",&beam_inst_Y,"beam_inst_Y/F");
    fTree->Branch("beam_inst_Z",&beam_inst_Z,"beam_inst_Z/F");

    fTree->Branch("trueBeamParticleStartPx",&fTrueBeamParticleStartPx,"trueBeamParticleStartPx/F");
    fTree->Branch("trueBeamParticleStartPy",&fTrueBeamParticleStartPy,"trueBeamParticleStartPy/F");
    fTree->Branch("trueBeamParticleStartPz",&fTrueBeamParticleStartPz,"trueBeamParticleStartPz/F");

  // Reco Tree Branches   
    fTree->Branch("bestMatchedMCParticleFromPFParticlePdgCode",&fbestMatchedMCParticleFromPFParticlePdgCode,"bestMatchedMCParticleFromPFParticlePdgCode/I");
    fTree->Branch("isBeamParticleReconstructed",&fIsBeamParticleReconstructed,"isBeamParticleReconstructed/i");


    fTree->Branch("recoBeamParticleStartPx",&fRecoBeamParticleStartPx,"recoBeamParticleStartPx/F");
    fTree->Branch("recoBeamParticleStartPy",&fRecoBeamParticleStartPy,"recoBeamParticleStartPy/F");
    fTree->Branch("recoBeamParticleStartPz",&fRecoBeamParticleStartPz,"recoBeamParticleStartPz/F");


    fTree->Branch("numberOfReconstructedBeamParticle",&fNumberOfReconstructedBeamParticle,"numberOfReconstructedBeamParticle/I");
    fTree->Branch("primaryBeamParticleLength",&fPrimaryBeamParticleLength,"primaryBeamParticleLength/F");
    fTree->Branch("trueBeamParticleNhits", &fTrueBeamParticleNhits, "trueBeamParticleNhits/I");

  // Interested quantities

    fTree->Branch("reco_beam_pfp_topology",&fReco_beam_pfp_topology,"reco_beam_pfp_topology/I");//

    fTree->Branch("purity",&fpurity,"purity/F");//
    fTree->Branch("completeness",&fcompleteness,"completeness/F");//

    fTree->Branch("recoBeamParticleStartX",&fRecoBeamParticleStartX,"recoBeamParticleStartX/F");//
    fTree->Branch("recoBeamParticleStartY",&fRecoBeamParticleStartY,"recoBeamParticleStartY/F");//
    fTree->Branch("recoBeamParticleStartZ",&fRecoBeamParticleStartZ,"recoBeamParticleStartZ/F");//

    fTree->Branch("beamInst_startVertex_X_SCE_corrected",&fBeamInst_startVertex_X_SCE_corrected,"beamInst_startVertex_X_SCE_corrected/F");//
    fTree->Branch("beamInst_startVertex_Y_SCE_corrected",&fBeamInst_startVertex_Y_SCE_corrected,"beamInst_startVertex_Y_SCE_corrected/F");//
    fTree->Branch("beamInst_startVertex_Z_SCE_corrected",&fBeamInst_startVertex_Z_SCE_corrected,"beamInst_startVertex_Z_SCE_corrected/F");//

    fTree->Branch("recoBeamParticleInteractionX",&fRecoBeamParticleInteractionX,"recoBeamParticleInteractionX/F");//
    fTree->Branch("recoBeamParticleInteractionY",&fRecoBeamParticleInteractionY,"recoBeamParticleInteractionY/F");//
    fTree->Branch("recoBeamParticleInteractionZ",&fRecoBeamParticleInteractionZ,"recoBeamParticleInteractionZ/F");//

    fTree->Branch("trueBeamParticleInteractionX",&fTrueBeamParticleInteractionX,"TrueBeamParticleInteractionX/F");//
    fTree->Branch("trueBeamParticleInteractionY",&fTrueBeamParticleInteractionY,"TrueBeamParticleInteractionY/F");//
    fTree->Branch("trueBeamParticleInteractionZ",&fTrueBeamParticleInteractionZ,"TrueBeamParticleInteractionZ/F");//

    fTree->Branch("recoBeamParticleNhits", &fRecoBeamParticleNhits, "recoBeamParticleNhits/I");

    fTree->Branch("trueBeamStartDirX",&fTrueBeamStartDirX,"trueBeamStartDirX/F");//
    fTree->Branch("trueBeamStartDirY",&fTrueBeamStartDirY,"trueBeamStartDirY/F");//
    fTree->Branch("trueBeamStartDirZ",&fTrueBeamStartDirZ,"trueBeamStartDirZ/F");//

    fTree->Branch("recoBeamStartDirX",&fRecoBeamStartDirX,"recoBeamStartDirX/F");//
    fTree->Branch("recoBeamStartDirY",&fRecoBeamStartDirY,"recoBeamStartDirY/F");//
    fTree->Branch("recoBeamStartDirZ",&fRecoBeamStartDirZ,"recoBeamStartDirZ/F");//

    fTree->Branch("primaryRecoPfpLength",&fPrimaryRecoPfpLength,"primaryRecoPfpLength/F");//
    fTree->Branch("beam_length_by_traj_points",&fBeam_length_by_traj_points,"beam_length_by_traj_points/F");//

    fTree->Branch("fTier0MCParticleHitsSize",&fTier0MCParticleHitsSize,"fTier0MCParticleHitsSize/I");//
    fTree->Branch("fTier0recoBeamParticleHitsSize",&fTier0recoBeamParticleHitsSize,"fTier0recoBeamParticleHitsSize/I");//
    fTree->Branch("fSharedTier0RecoTrueHitsSize",&fSharedTier0RecoTrueHitsSize,"fSharedTier0RecoTrueHitsSize/I");//
    fTree->Branch("fTier0RecoMCParticleMatch",&fTier0RecoMCParticleMatch,"fTier0RecoMCParticleMatch/I");//

    fTree->Branch("tier1purity",&fTier1purity,"tier1purity/F");//
    fTree->Branch("tier1completeness",&fTier1completeness,"tier1completeness/F");//
    fTree->Branch("tier1beamInst_startVertex_X_SCE_corrected",&fTier1BeamInst_startVertex_X_SCE_corrected,"tier1beamInst_startVertex_X_SCE_corrected/F");//
    fTree->Branch("tier1beamInst_startVertex_Y_SCE_corrected",&fTier1BeamInst_startVertex_Y_SCE_corrected,"tier1beamInst_startVertex_Y_SCE_corrected/F");//
    fTree->Branch("tier1beamInst_startVertex_Z_SCE_corrected",&fTier1BeamInst_startVertex_Z_SCE_corrected,"tier1beamInst_startVertex_Z_SCE_corrected/F");//
    fTree->Branch("tier1beam_length_by_traj_points",&fTier1Beam_length_by_traj_points,"tier1beam_length_by_traj_points/F");//
    fTree->Branch("tier1recoBeamParticleStartX",&fTier1RecoBeamParticleStartX,"tier1recoBeamParticleStartX/F");//
    fTree->Branch("tier1recoBeamParticleStartY",&fTier1RecoBeamParticleStartY,"tier1recoBeamParticleStartY/F");//
    fTree->Branch("tier1recoBeamParticleStartZ",&fTier1RecoBeamParticleStartZ,"tier1recoBeamParticleStartZ/F");//
    fTree->Branch("tier1recoBeamParticleInteractionX",&fTier1RecoBeamParticleInteractionX,"tier1recoBeamParticleInteractionX/F");//
    fTree->Branch("tier1recoBeamParticleInteractionY",&fTier1RecoBeamParticleInteractionY,"tier1recoBeamParticleInteractionY/F");//
    fTree->Branch("tier1recoBeamParticleInteractionZ",&fTier1RecoBeamParticleInteractionZ,"tier1recoBeamParticleInteractionZ/F");//
    fTree->Branch("tier1recoBeamStartDirX",&fTier1RecoBeamStartDirX,"tier1recoBeamStartDirX/F");//
    fTree->Branch("tier1recoBeamStartDirY",&fTier1RecoBeamStartDirY,"tier1recoBeamStartDirY/F");//
    fTree->Branch("tier1recoBeamStartDirZ",&fTier1RecoBeamStartDirZ,"tier1recoBeamStartDirZ/F");//
    fTree->Branch("tier1primaryRecoPfpLength",&fTier1PrimaryRecoPfpLength,"tier1primaryRecoPfpLength/F");//
    fTree->Branch("fTier1MCParticleNumber",&fTier1MCParticleNumber,"fTier1MCParticleNumber/I");//
    fTree->Branch("fTier1RecoParticleNumber",&fTier1RecoParticleNumber,"fTier1RecoParticleNumber/I");//
    fTree->Branch("fTier1MCRecoMatchedNumber",&fTier1MCRecoMatchedNumber,"fTier1MCRecoMatchedNumber/I");//
    fTree->Branch("fTier1TrueBeamParticleLength",&fTier1TrueBeamParticleLength,"fTier1TrueBeamParticleLength/F");//
    fTree->Branch("fTier1TrueBeamParticleInteractionX",&fTier1TrueBeamParticleInteractionX,"fTier1TrueBeamParticleInteractionX/F");//
    fTree->Branch("fTier1TrueBeamParticleInteractionY",&fTier1TrueBeamParticleInteractionY,"fTier1TrueBeamParticleInteractionY/F");//fix
    fTree->Branch("fTier1TrueBeamParticleInteractionZ",&fTier1TrueBeamParticleInteractionZ,"fTier1TrueBeamParticleInteractionZ/F");//fix
    fTree->Branch("fTier1RecoBeamParticleNhits",&fTier1RecoBeamParticleNhits,"fTier1RecoBeamParticleNhits/I");//
    fTree->Branch("fTier1TrueBeamParticleNhits",&fTier1TrueBeamParticleNhits,"fTier1TrueBeamParticleNhits/I");//
    fTree->Branch("fTier1TrueBeamParticlePDGCode",&fTier1TrueBeamParticlePDGCode,"fTier1TrueBeamParticlePDGCode/I");//
    fTree->Branch("fStartVertexDr",&fStartVertexDr,"fStartVertexDr/F");
    fTree->Branch("fInteractionVertexDr",&fInteractionVertexDr,"fInteractionVertexDr/F");
    fTree->Branch("fTrueBeamLengthVersion3",&fTrueBeamLengthVersion3,"fTrueBeamLengthVersion3/F");
    fTree->Branch("fTier1StartVertexDr",&fTier1StartVertexDr,"fTier1StartVertexDr/F");
}

void analysis::PDSPKaonAnalysis::endJob()
{
  // Implementation of optional member function here.
}

/*bool analysis::PDSPKaonAnalysis::sortcol( const std::vector<int>& v1, const std::vector<int>& v2 )    //<=============THIS ONE, DOM
{
    return v1[2] > v2[2];
}*/
std::vector<const simb::MCParticle*> analysis::PDSPKaonAnalysis::apply_true_particle_quality_cuts(std::vector<const simb::MCParticle*> inputVector)
{
//    std::cout << "testing apply_true_particle_quality_cuts function" << std::endl;
    int size = inputVector.size();
    std::vector<const simb::MCParticle*> outputVector;
    const sim::ParticleList & plist = pi_serv->ParticleList();
    for (int i = 0; i < size; i++)
    {
            
        if ( ((inputVector[i]->PdgCode()) != 111) && ((inputVector[i]->PdgCode()) != 2112) && ((inputVector[i]->PdgCode()) < 1000000000))
            outputVector.push_back(inputVector[i]);
        
        else
        {
//            std::cout << "------------------granddaughters------------------------" << std::endl;
            int number_true_granddaughter = inputVector[i]->NumberDaughters();
//            std::cout << "number_true_granddaughter: " << number_true_granddaughter << std::endl;
            for (int j = 0; j < number_true_granddaughter; j++)
            {
                auto granddaughter = plist[ inputVector[i]->Daughter(j)];
                outputVector.push_back(granddaughter);
            }
//            std::cout << "--------------------------------------------------------" << std::endl;
            }
        }        

    return outputVector;

}

int analysis::PDSPKaonAnalysis::count_unwanted_particles(std::vector<const simb::MCParticle*> inputVector)
{
//    std::cout << "testing count_unwanted_particles function" << std::endl;
    int size = inputVector.size();
    int count = 0;
    for (int i = 0; i < size; i++)
    {
        if ( ((inputVector[i]->PdgCode()) == 111) || ((inputVector[i]->PdgCode()) == 2112) || ((inputVector[i]->PdgCode()) > 1000000000))
            ++count;
    }

//    std::cout << "count: " << count << std::endl;

    return count;
}


DEFINE_ART_MODULE(analysis::PDSPKaonAnalysis)
