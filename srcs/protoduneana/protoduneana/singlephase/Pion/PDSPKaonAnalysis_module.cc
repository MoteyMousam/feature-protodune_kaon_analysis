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

#include "TMath.h"
#include "TVirtualFitter.h"
#include "TPolyLine3D.h"
#include "Math/Vector3D.h"
#include "TGraph2D.h"
#include "TROOT.h"

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
//double distance2(double x, double y, double z, double * p);
//void line(double t, double * p, double & x, double & y, double & z);
//void SumDistance2(int &, double *, double & sum, double * par, int);
//TVector3 FitLine(const std::vector<TVector3> & input);

  TTree *fTree;
  TTree *fTestTree;

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
  std::string fCalorimetryTagSCE;

  unsigned int fEventID;
  int fTier0MCPDGCode;
  float fTier0Completeness;
  float fTier0Purity;
  float fTier0MCStartVertexX;
  float fTier0MCStartVertexY;
  float fTier0MCStartVertexZ;
  float fTier0RecoStartVertexX_SCE_corrected;
  float fTier0RecoStartVertexY_SCE_corrected;
  float fTier0RecoStartVertexZ_SCE_corrected;
  float fTier0StartVertexDr;
  float fTier0MCInteractionVertexX;
  float fTier0MCInteractionVertexY;
  float fTier0MCInteractionVertexZ;
  float fTier0RecoInteractionVertexX_SCE_corrected;
  float fTier0RecoInteractionVertexY_SCE_corrected;
  float fTier0RecoInteractionVertexZ_SCE_corrected;
  float fTier0InteractionVertexDr;
  float fTier0MCStartDirectionX;
  float fTier0MCStartDirectionY;
  float fTier0MCStartDirectionZ;
  float fTier0RecoStartDirectionX_SCE_corrected;
  float fTier0RecoStartDirectionY_SCE_corrected;
  float fTier0RecoStartDirectionZ_SCE_corrected;
  float fTier0MCLengthByTrajPoints;
  float fTier0RecoLength;
  int fTier0MCParticleHitsSize;
  int fTier0RecoMCParticleMatch;




  float fRecoBeamParticleStartX;
  float fRecoBeamParticleStartY;
  float fRecoBeamParticleStartZ;

  float fRecoBeamParticleInteractionX;
  float fRecoBeamParticleInteractionY;
  float fRecoBeamParticleInteractionZ;

  float fBeamInst_startVertex_X_SCE_corrected;
  float fBeamInst_startVertex_Y_SCE_corrected;
  float fBeamInst_startVertex_Z_SCE_corrected;

  float fBeamInst_startVertex_dr_SCE_corrected; // need to code
  
  int fDoesTrueBeamParticleExist; // need to code
  int fDoesRecoBeamParticleExist;


  float fPrimaryBeamParticleLength;

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





  int fNumberOfReconstructedBeamParticle;
  bool fMCHasBI;
  bool beam_inst_valid;

  int fReco_beam_pfp_topology;


//  float fTier0recoBeamParticleHitsSize;
//  float fSharedTier0RecoTrueHitsSize;

  float fStartVertexDr;



  protoana::ProtoDUNETruthUtils fProtoDUNETruthUtils;
  protoana::ProtoDUNEPFParticleUtils fProtoDUNEPFParticleUtils;
  protoana::ProtoDUNETrackUtils fProtoDUNETrackUtils;
  protoana::ProtoDUNEShowerUtils fProtoDUNEShowerUtils;

  art::ServiceHandle< cheat::ParticleInventoryService > pi_serv;

//---------------------------Tier 1-----------------------------------------
  int fTier1MCPDGCode;
  float fTier1Completeness;
  float fTier1Purity;
//  float fTier1MCStartVertexX;
//  float fTier1MCStartVertexY;
//  float fTier1MCStartVertexZ;
  float fTier1RecoStartVertexX_SCE_corrected;
  float fTier1RecoStartVertexY_SCE_corrected;
  float fTier1RecoStartVertexZ_SCE_corrected;
  float fTier1StartVertexDr;
//  float fTier1MCInteractionVertexX;
//  float fTier1MCInteractionVertexY;
//  float fTier1MCInteractionVertexZ;
////  float fTier1RecoInteractionVertexX_SCE_corrected;
//  float fTier1RecoInteractionVertexY_SCE_corrected;
//  float fTier1RecoInteractionVertexZ_SCE_corrected;
//  float fTier1InteractionVertexDr;
  float fTier1MCStartDirectionX;
  float fTier1MCStartDirectionY;
  float fTier1MCStartDirectionZ;
  float fTier1RecoStartDirectionX_SCE_corrected;
  float fTier1RecoStartDirectionY_SCE_corrected;
  float fTier1RecoStartDirectionZ_SCE_corrected;
  float fTier1MCLengthByTrajPoints;
  float fTier1RecoLength;
  int fTier1MCParticleHitsSize;
  int fTier1RecoMCParticleMatch;


  float fTier1TrueBeamParticleEnergy;

  int fTier1RecoBeamParticlePDGCode;
  int fTier1IsBeamParticleReconstructed;
  int fTier1TrueBeamParticlePDGCode;
  int fTier1bestMatchedMCParticleFromPFParticlePdgCode;

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

//  float fTier1PrimaryRecoPfpLength;
  float fTier1PrimaryBeamParticleLength;
  float fTier1Beam_length_by_traj_points;
  int fTier1Reco_beam_pfp_topology;

  int fTier1recoBeamParticleHitsSize;
  int fSharedTier1RecoTrueHitsSize;

  int fTier1MCParticleNumber;
  int fTier1RecoParticleNumber;
  int fTier1MCRecoMatchedNumber;
  float fTier1TrueBeamParticleLength;

  int fMCRecoMatchedHits;

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
    fCalorimetryTagSCE  = p.get<std::string>("CalorimetryTagSCE");
}

void analysis::PDSPKaonAnalysis::analyze(art::Event const& e)
{
    std::cout << "====================== Starting Kaon Analyser ======================" << std::endl;

  // Implementation of required member function here.

    fEventID = e.id().event();
//    fTruePrimaryEnergy = -999.;
    fTrueBeamParticleEnergy                     = -9999.;
    fTier0MCPDGCode                             = -9999;
    fTier0Purity                                = -1.f;
    fTier0Completeness                          = -1.f;

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

    fTier0MCStartVertexX                        = -9999;
    fTier0MCStartVertexY                        = -9999;
    fTier0MCStartVertexZ                        = -9999;
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

    fTier0MCInteractionVertexX               = -9999;
    fTier0MCInteractionVertexY               = -9999;
    fTier0MCInteractionVertexZ               = -9999;

    fNumberOfReconstructedBeamParticle          = 1;

    double tpcActiveXLow                        = -358.5282287595;
    double tpcActiveXHigh                       = 358.5282287595;
    double tpcActiveYLow                        = 0;
    double tpcActiveYHigh                       = 603.8612670895;
    double tpcActiveZLow                        = 0;
    double tpcActiveZHigh                       = 695.2862548825;

    fTier0RecoLength                            = -9999;
    fPrimaryBeamParticleLength                  = -9999;
    fTier0MCLengthByTrajPoints                  = -9999;

    fReco_beam_pfp_topology                     = -9999;

    fTier0MCParticleHitsSize                    = -9999;
//    fTier0recoBeamParticleHitsSize              = -9999;
//    fSharedTier0RecoTrueHitsSize                = -9999;
    fTier0RecoMCParticleMatch                   = 0;
    fDoesTrueBeamParticleExist                  = -9999;
    fStartVertexDr                              = -9999;
    fTier0RecoStartVertexX_SCE_corrected        = -9999;
    fTier0RecoStartVertexY_SCE_corrected        = -9999;
    fTier0RecoStartVertexZ_SCE_corrected        = -9999;

    fTier0RecoInteractionVertexX_SCE_corrected  = -9999;
    fTier0RecoInteractionVertexY_SCE_corrected  = -9999;
    fTier0RecoInteractionVertexZ_SCE_corrected  = -9999;
    fTier0StartVertexDr                         = -9999;
    fTier0InteractionVertexDr                   = -9999;
  
    fTier0RecoStartDirectionX_SCE_corrected     = -9999; 
    fTier0RecoStartDirectionY_SCE_corrected     = -9999; 
    fTier0RecoStartDirectionZ_SCE_corrected     = -9999; 

    fTier0MCStartDirectionX                     = -9999;
    fTier0MCStartDirectionY                     = -9999;
    fTier0MCStartDirectionZ                     = -9999;
//------------------------------Tier 1--------------------------------------
    fTier1TrueBeamParticleEnergy                     = -9999.;
    fTier1TrueBeamParticlePDGCode                    = -9999;
    fTier1Purity                                     = -1.f;
    fTier1Completeness                               = -1.f;

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

    fTier1RecoLength                                 = -9999;
    fTier1PrimaryBeamParticleLength                  = -9999;
    fTier1Beam_length_by_traj_points                 = -9999;

    fTier1Reco_beam_pfp_topology                     = -9999;

    fTier1MCParticleHitsSize                        = -9999;
    fTier1recoBeamParticleHitsSize                  = -9999;
    fSharedTier1RecoTrueHitsSize                    = -9999;
    fTier1RecoMCParticleMatch                       = 0;

    fTier1MCParticleNumber                          = -9999;
    fTier1RecoParticleNumber                        = -9999; 
    fTier1MCRecoMatchedNumber                       = -9999;
    fTier1TrueBeamParticleLength                    = -9999;
    fTier1StartVertexDr                             = -9999;
    fMCRecoMatchedHits                              = -9999;

    fTier1RecoBeamParticleStartX                    = -9999;
    fTier1RecoBeamParticleStartY                    = -9999;
    fTier1RecoBeamParticleStartZ                    = -9999;

    fTier1MCStartDirectionX                         = -9999;
    fTier1MCStartDirectionY                         = -9999;
    fTier1MCStartDirectionZ                         = -9999;
    fTier1RecoStartDirectionX_SCE_corrected         = -9999;
    fTier1RecoStartDirectionY_SCE_corrected         = -9999;
    fTier1RecoStartDirectionZ_SCE_corrected         = -9999;
    fTier1MCLengthByTrajPoints                      = -9999;
    fTier1MCPDGCode                                 = -9999;
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
    std::vector < const recob::PFParticle*> tier1RecoParticles;
    std::vector < const simb::MCParticle*> tier1TrueParticlesQ1;

    const sim::ParticleList & plist = pi_serv->ParticleList();

    if(!e.isRealData())
    {
        auto mcTruths = e.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
        
        true_beam_particle = fProtoDUNETruthUtils.GetGeantGoodParticle((*mcTruths)[0],e);
        fDoesTrueBeamParticleExist = 1;
        tier0MCParticleHits = fProtoDUNETruthUtils.GetMCParticleHits(clockData, *true_beam_particle, e, "hitpdune" );
    }
    
    else
        fDoesTrueBeamParticleExist = 0;        
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
        fDoesRecoBeamParticleExist = 1;

        for(const recob::PFParticle* particle : beamParticles)
        {
            tier0recoBeamParticleHits = fProtoDUNEPFParticleUtils.GetPFParticleHits(*particle, e, fPFParticleLabel);
//            std::cout << "tier0recoBeamParticleHits size: " << tier0recoBeamParticleHits.size() << std::endl;
        }
    }

    else
        std::cout << "More than one reconstructed beam particle!" << std::endl;

    int fTier0recoBeamParticleHitsSize = tier0recoBeamParticleHits.size();



//SharedHits block    


    std::vector< const recob::Hit* > sharedTier0RecoTrueHits = fProtoDUNETruthUtils.GetSharedHits(clockData, *true_beam_particle, *beamParticles[0], e, fPFParticleLabel);
//    std::cout << "sharedTier0RecoTrueHits size: " << sharedTier0RecoTrueHits.size() << std::endl;

    float fSharedTier0RecoTrueHitsSize = sharedTier0RecoTrueHits.size();
    
    float testPurity = fProtoDUNETruthUtils.GetPurity(clockData, *beamParticles[0], e, fPFParticleLabel); 
    float recoPurity = fSharedTier0RecoTrueHitsSize/fTier0recoBeamParticleHitsSize;
    float testCompleteness = fProtoDUNETruthUtils.GetCompleteness(clockData, *beamParticles[0], e, "pandora", "hitpdune");
    float recoCompleteness = fSharedTier0RecoTrueHitsSize/fTier0MCParticleHitsSize;
//Reconstruction efficiency criteria
//    std::cout << "recoPurity: " << recoPurity << ", recoCompleteness: " << recoCompleteness << std::endl;

    if( ((tier0MCParticleHits.size()) != 0) && (recoPurity >= 0.1) && (recoCompleteness >= 0.5) )
        fTier0RecoMCParticleMatch = 1;

//calculating quatities
//    std::cout << "1" << std::endl;      
//    std::cout << "fTier0RecoMCParticleMatch: " << fTier0RecoMCParticleMatch << std::endl;
    if(fTier0RecoMCParticleMatch == 1)
    {
        //Truth quatities
        fTrueBeamParticleEnergy     = true_beam_particle->E();
        fTier0MCPDGCode             = true_beam_particle->PdgCode();
        fTrueBeamParticleStartX     = true_beam_particle->Position(0).X();
        fTrueBeamParticleStartY     = true_beam_particle->Position(0).Y();
        fTrueBeamParticleStartZ     = true_beam_particle->Position(0).Z();
        fTrueBeamParticleEndX       = true_beam_particle->EndX();
        fTrueBeamParticleEndY       = true_beam_particle->EndY();
        fTrueBeamParticleEndZ       = true_beam_particle->EndZ();
        fTrueBeamParticleStartPx    = true_beam_particle->Px();
        fTrueBeamParticleStartPy    = true_beam_particle->Py();
        fTrueBeamParticleStartPz    = true_beam_particle->Pz();
        fTier0MCStartDirectionX     = (fTrueBeamParticleStartPx)/(true_beam_particle->P());
        fTier0MCStartDirectionY     = (fTrueBeamParticleStartPy)/(true_beam_particle->P());
        fTier0MCStartDirectionZ     = (fTrueBeamParticleStartPz)/(true_beam_particle->P());
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
//        std::cout << "2" << std::endl;
        if( nTracks > 0 )
        {
            fTier0MCStartVertexX = beamEvent.GetBeamTracks()[0].Trajectory().End().X();
            fTier0MCStartVertexY = beamEvent.GetBeamTracks()[0].Trajectory().End().Y();
            fTier0MCStartVertexZ = beamEvent.GetBeamTracks()[0].Trajectory().End().Z();

//            fBeamInst_startVertex_X_SCE_corrected = beam_inst_X-SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).X();
//            fBeamInst_startVertex_Y_SCE_corrected = beam_inst_Y+SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).Y();
//            fBeamInst_startVertex_Z_SCE_corrected = beam_inst_Z+SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).Z();
        }

/*        float true_beam_particle_SCE_end_X = fTrueBeamParticleEndX-SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleEndX,fTrueBeamParticleEndY,fTrueBeamParticleEndZ)).X();
        float true_beam_particle_SCE_end_Y = fTrueBeamParticleEndY+SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleEndX,fTrueBeamParticleEndY,fTrueBeamParticleEndZ)).Y();
        float true_beam_particle_SCE_end_Z = fTrueBeamParticleEndZ+SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleEndX,fTrueBeamParticleEndY,fTrueBeamParticleEndZ)).Z();*/

        fTier0MCInteractionVertexX = fTrueBeamParticleEndX;       
        fTier0MCInteractionVertexY = fTrueBeamParticleEndY;       
        fTier0MCInteractionVertexZ = fTrueBeamParticleEndZ;       

        fTier0MCLengthByTrajPoints = sqrt(((fTier0MCInteractionVertexX-fTier0MCStartVertexX)*(fTier0MCInteractionVertexX-fTier0MCStartVertexX)) + ((fTier0MCInteractionVertexY-fTier0MCStartVertexY)*(fTier0MCInteractionVertexY-fTier0MCStartVertexY)) + ((fTier0MCInteractionVertexZ-fTier0MCStartVertexZ)*(fTier0MCInteractionVertexZ-fTier0MCStartVertexZ)));

        fTrueBeamParticleNhits = fProtoDUNETruthUtils.GetMCParticleHits( clockData, *true_beam_particle, e, fHitTag ).size();

        fPrimaryBeamParticleLength = fProtoDUNETruthUtils.GetMCParticleLengthInTPCActiveVolume(*true_beam_particle, tpcActiveXLow, tpcActiveXHigh, tpcActiveYLow, tpcActiveYHigh, tpcActiveZLow, tpcActiveZHigh);

//        std::vector < const recob::Hit* > trueBeamHits = fProtoDUNETruthUtils.GetMCParticleHits( clockData, *true_beam_particle, e, fHitTag );
//        art::FindManyP<recob::SpacePoint> spFromHits(trueBeamHits, e, fHitTag);

//        std::cout << " X: " << (trueBeamHits[0]->get())->position().X() << std::endl;

        //Reco Quantities
        
        fRecoBeamParticleNhits = fProtoDUNEPFParticleUtils.GetPFParticleHits( *beamParticles[0], e, fPFParticleLabel).size();

        const recob::Track* thisTrack = fProtoDUNEPFParticleUtils.GetPFParticleTrack(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);
        const recob::Shower* thisShower = fProtoDUNEPFParticleUtils.GetPFParticleShower(*beamParticles[0],e,fPFParticleLabel,fShowerLabel);

        if(thisTrack != 0x0)
        {
            fReco_beam_pfp_topology = 1;
            fTier0RecoLength = thisTrack->Length();
//            std::cout << "3" << std::endl;
            auto calo = fProtoDUNETrackUtils.GetRecoTrackCalorimetry(*thisTrack, e, fTrackLabel, fCalorimetryTagSCE);

            const TVector3 reco_primary_start_vertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);

            fRecoBeamParticleStartX = reco_primary_start_vertex.X();
            fRecoBeamParticleStartY = reco_primary_start_vertex.Y();
            fRecoBeamParticleStartZ = reco_primary_start_vertex.Z();

            fTier0RecoStartVertexX_SCE_corrected = fRecoBeamParticleStartX+SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleStartX,fRecoBeamParticleStartY,fRecoBeamParticleStartZ)).X();
            fTier0RecoStartVertexY_SCE_corrected = fRecoBeamParticleStartY-SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleStartX,fRecoBeamParticleStartY,fRecoBeamParticleStartZ)).Y();
            fTier0RecoStartVertexZ_SCE_corrected = fRecoBeamParticleStartZ-SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleStartX,fRecoBeamParticleStartY,fRecoBeamParticleStartZ)).Z();

            float startVertexDx = fTier0RecoStartVertexX_SCE_corrected - fTier0MCStartVertexX;
            float startVertexDy = fTier0RecoStartVertexY_SCE_corrected - fTier0MCStartVertexY;
            float startVertexDz = fTier0RecoStartVertexZ_SCE_corrected - fTier0MCStartVertexZ;
            
            fTier0StartVertexDr = sqrt(((startVertexDx)*(startVertexDx))+((startVertexDy)*(startVertexDy))+((startVertexDz)*(startVertexDz)));
            

            const TVector3 reco_primary_interaction_vertex = fProtoDUNEPFParticleUtils.GetPFParticleSecondaryVertex(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);

            fRecoBeamParticleInteractionX = reco_primary_interaction_vertex.X();
            fRecoBeamParticleInteractionY = reco_primary_interaction_vertex.Y();
            fRecoBeamParticleInteractionZ = reco_primary_interaction_vertex.Z();

            fTier0RecoInteractionVertexX_SCE_corrected = fRecoBeamParticleInteractionX+SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleInteractionX,fRecoBeamParticleInteractionY,fRecoBeamParticleInteractionZ)).X();
            fTier0RecoInteractionVertexY_SCE_corrected = fRecoBeamParticleInteractionY-SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleInteractionX,fRecoBeamParticleInteractionY,fRecoBeamParticleInteractionZ)).Y();
            fTier0RecoInteractionVertexZ_SCE_corrected = fRecoBeamParticleInteractionZ-SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleInteractionX,fRecoBeamParticleInteractionY,fRecoBeamParticleInteractionZ)).Z();

            fRecoBeamStartDirX = (thisTrack->StartDirection()).X();
            fRecoBeamStartDirY = (thisTrack->StartDirection()).Y();
            fRecoBeamStartDirZ = (thisTrack->StartDirection()).Z();

            float interactionVertexDx = fTier0RecoInteractionVertexX_SCE_corrected - fTier0MCInteractionVertexX;
            float interactionVertexDy = fTier0RecoInteractionVertexY_SCE_corrected - fTier0MCInteractionVertexY;
            float interactionVertexDz = fTier0RecoInteractionVertexZ_SCE_corrected - fTier0MCInteractionVertexZ;


            fTier0InteractionVertexDr = sqrt(((interactionVertexDx)*(interactionVertexDx))+((interactionVertexDy)*(interactionVertexDy))+((interactionVertexDz)*(interactionVertexDz)));
//            std::cout << "4" << std::endl;

//-------------------------------------------------------------------------------------------------------------------

            bool found_calo = false;
            size_t index = 0;
            for ( index = 0; index < calo.size(); ++index)
            {
                if (calo[index].PlaneID().Plane == 2)
                {
                    found_calo = true;
                    break; 
                }
            }

            if (found_calo)
            {
                auto theXYZPoints = calo[index].XYZ();
                std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b){return (a.Z() < b.Z());});

                if (theXYZPoints.size())
                {
                    TVector3 dir((theXYZPoints.back().X() - theXYZPoints[0].X()), (theXYZPoints.back().Y() - theXYZPoints[0].Y()), (theXYZPoints.back().Z() - theXYZPoints[0].Z()));
                    fTier0RecoStartDirectionX_SCE_corrected = dir.Unit().X();
                    fTier0RecoStartDirectionY_SCE_corrected = dir.Unit().Y();
                    fTier0RecoStartDirectionZ_SCE_corrected = dir.Unit().Z();
                }

                else
                {
                    fTier0RecoStartDirectionX_SCE_corrected = -9999;
                    fTier0RecoStartDirectionY_SCE_corrected = -9999;
                    fTier0RecoStartDirectionZ_SCE_corrected = -9999;
                }

                if (theXYZPoints.size() > 1)
                {
                    TVector3 start_p1(theXYZPoints[0].X(), theXYZPoints[0].Y(), theXYZPoints[0].Z());
                    TVector3 start_p2(theXYZPoints[1].X(), theXYZPoints[1].Y(), theXYZPoints[1].Z());
                    TVector3 start_diff = start_p2 - start_p1;

                    fTier0RecoStartDirectionX_SCE_corrected = start_diff.Unit().X();
                    fTier0RecoStartDirectionY_SCE_corrected = start_diff.Unit().Y();
                    fTier0RecoStartDirectionZ_SCE_corrected = start_diff.Unit().Z();
                }

                else
                {
                    fTier0RecoStartDirectionX_SCE_corrected = -9999;
                    fTier0RecoStartDirectionY_SCE_corrected = -9999;
                    fTier0RecoStartDirectionZ_SCE_corrected = -9999;
                }

/*                if (theXYZPoints.size() > 2)
                {
                    std::vector<TVector3> input;
                    for (size_t iP = 0; iP < 3; ++iP) 
                    {
                        input.push_back(TVector3(theXYZPoints[iP].X(), theXYZPoints[iP].Y(), theXYZPoints[iP].Z()));
                    }

                    TVector3 startDiff = FitLine(input);

                    fTier0RecoStartDirectionX_SCE_corrected = startDiff.Unit().X();
                    fTier0RecoStartDirectionY_SCE_corrected = startDiff.Unit().Y();
                    fTier0RecoStartDirectionZ_SCE_corrected = startDiff.Unit().Z();
                }
                
                else
                {
                    fTier0RecoStartDirectionX_SCE_corrected = -9999;
                    fTier0RecoStartDirectionY_SCE_corrected = -9999;
                    fTier0RecoStartDirectionZ_SCE_corrected = -9999;
                }*/
            }
//---------------------------------------------------------------------------------------------------------------------
        }

        if(thisShower != 0x0) 
        {
            fReco_beam_pfp_topology = 0;
            fTier0RecoLength = thisShower->Length();

            auto calo = fProtoDUNEShowerUtils.GetRecoShowerCalorimetry(*thisShower, e, fShowerLabel, fCalorimetryTagSCE);

            const TVector3 reco_primary_start_vertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);
//            showerEnd = showerVertex + showerDir.unit()*ShowerLength;
//            SCE showerEnd and showerVertex



            fRecoBeamParticleStartX = reco_primary_start_vertex.X();
            fRecoBeamParticleStartY = reco_primary_start_vertex.Y();
            fRecoBeamParticleStartZ = reco_primary_start_vertex.Z();

            fTier0RecoStartVertexX_SCE_corrected = fRecoBeamParticleStartX+SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleStartX,fRecoBeamParticleStartY,fRecoBeamParticleStartZ)).X();
            fTier0RecoStartVertexY_SCE_corrected = fRecoBeamParticleStartY-SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleStartX,fRecoBeamParticleStartY,fRecoBeamParticleStartZ)).Y();
            fTier0RecoStartVertexZ_SCE_corrected = fRecoBeamParticleStartZ-SCE->GetPosOffsets(geo::Point_t(fRecoBeamParticleStartX,fRecoBeamParticleStartY,fRecoBeamParticleStartZ)).Z();

            float startVertexDx = fTier0RecoStartVertexX_SCE_corrected - fTier0MCStartVertexX;
            float startVertexDy = fTier0RecoStartVertexY_SCE_corrected - fTier0MCStartVertexY;
            float startVertexDz = fTier0RecoStartVertexZ_SCE_corrected - fTier0MCStartVertexZ;
            
            fTier0StartVertexDr = sqrt(((startVertexDx)*(startVertexDx))+((startVertexDy)*(startVertexDy))+((startVertexDz)*(startVertexDz)));

            fRecoBeamStartDirX = (thisShower->Direction()).X();
            fRecoBeamStartDirY = (thisShower->Direction()).Y();
            fRecoBeamStartDirZ = (thisShower->Direction()).Z();

            fTier0InteractionVertexDr = 0;

//-------------------------------------------------------------------------------------------------------------------

            bool found_calo = false;
            size_t index = 0;
            for ( index = 0; index < calo.size(); ++index)
            {
                if (calo[index].PlaneID().Plane == 2)
                {
                    found_calo = true;
                    break; 
                }
            }

            if (found_calo)
            {
                auto theXYZPoints = calo[index].XYZ();
                std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b){return (a.Z() < b.Z());});

                if (theXYZPoints.size())
                {
                    TVector3 dir((theXYZPoints.back().X() - theXYZPoints[0].X()), (theXYZPoints.back().Y() - theXYZPoints[0].Y()), (theXYZPoints.back().Z() - theXYZPoints[0].Z()));
                    fTier0RecoStartDirectionX_SCE_corrected = dir.Unit().X();
                    fTier0RecoStartDirectionY_SCE_corrected = dir.Unit().Y();
                    fTier0RecoStartDirectionZ_SCE_corrected = dir.Unit().Z();
                }

                else
                {
                    fTier0RecoStartDirectionX_SCE_corrected = -9999;
                    fTier0RecoStartDirectionY_SCE_corrected = -9999;
                    fTier0RecoStartDirectionZ_SCE_corrected = -9999;
                }

                if (theXYZPoints.size() > 1)
                {
                    TVector3 start_p1(theXYZPoints[0].X(), theXYZPoints[0].Y(), theXYZPoints[0].Z());
                    TVector3 start_p2(theXYZPoints[1].X(), theXYZPoints[1].Y(), theXYZPoints[1].Z());
                    TVector3 start_diff = start_p2 - start_p1;

                    fTier0RecoStartDirectionX_SCE_corrected = start_diff.Unit().X();
                    fTier0RecoStartDirectionY_SCE_corrected = start_diff.Unit().Y();
                    fTier0RecoStartDirectionZ_SCE_corrected = start_diff.Unit().Z();
                }

                else
                {
                    fTier0RecoStartDirectionX_SCE_corrected = -9999;
                    fTier0RecoStartDirectionY_SCE_corrected = -9999;
                    fTier0RecoStartDirectionZ_SCE_corrected = -9999;
                }

/*                if (theXYZPoints.size() > 2)
                {
                    std::vector<TVector3> input;
                    for (size_t iP = 0; iP < 3; ++iP) 
                    {
                        input.push_back(TVector3(theXYZPoints[iP].X(), theXYZPoints[iP].Y(), theXYZPoints[iP].Z()));
                    }

                    TVector3 startDiff = FitLine(input);

                    fTier0RecoStartDirectionX_SCE_corrected = startDiff.Unit().X();
                    fTier0RecoStartDirectionY_SCE_corrected = startDiff.Unit().Y();
                    fTier0RecoStartDirectionZ_SCE_corrected = startDiff.Unit().Z();
                }
                
                else
                {
                    fTier0RecoStartDirectionX_SCE_corrected = -9999;
                    fTier0RecoStartDirectionY_SCE_corrected = -9999;
                    fTier0RecoStartDirectionZ_SCE_corrected = -9999;
                }*/
            }
//---------------------------------------------------------------------------------------------------------------------
//                    std::cout << "5" << std::endl;
        }

        const simb::MCParticle* bestMatchedMcParticleFromPFParticle;
        bestMatchedMcParticleFromPFParticle = (fProtoDUNETruthUtils.GetMCParticleFromPFParticle(clockData, *beamParticles[0], e, fPFParticleLabel));
        fbestMatchedMCParticleFromPFParticlePdgCode = bestMatchedMcParticleFromPFParticle->PdgCode();
    
        fTier0Purity = testPurity;
//        std::cout << "purity: " << fpurity << std::endl; 

        fTier0Completeness = testCompleteness;
//        std::cout << "completeness: " << fcompleteness << std::endl;

/*        float startVertexDx = fRecoBeamParticleStartX - fBeamInst_startVertex_X_SCE_corrected;
        float startVertexDy = fRecoBeamParticleStartY - fBeamInst_startVertex_Y_SCE_corrected;
        float startVertexDz = fRecoBeamParticleStartZ - fBeamInst_startVertex_Z_SCE_corrected;

        fStartVertexDr = sqrt(((startVertexDx)*(startVertexDx))+((startVertexDy)*(startVertexDy))+((startVertexDz)*(startVertexDz)));*/

/*        std::cout << "---------------------------------------TIER 0 PARTICLE-------------------------------------- " << std::endl; 
        std::cout << "purity: " << fpurity << std::endl; 
        std::cout << "completeness: " << fcompleteness << std::endl;
        std::cout << "fSharedTier0RecoTrueHitsSize: " << fSharedTier0RecoTrueHitsSize << std::endl; 
        std::cout << "fTrueBeamParticleNhits: " << fTrueBeamParticleNhits << std::endl; 
        std::cout << "fRecoBeamParticleNhits: " << fRecoBeamParticleNhits << std::endl; 
        std::cout << "-------------------------------------------------------------------------------------------- " << std::endl; */
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
        fTier1MCParticleNumber = array_row_size;
//        std::cout << "array_row_size: " << array_row_size << std::endl;
        int array_col_size = tier1RecoParticles.size();
        fTier1RecoParticleNumber = array_col_size;
         
//        std::cout << "array_col_size: " << array_col_size << std::endl;

        std::vector< std::vector<int> > twoDTrueRecoSharedHitsMatrix;
        std::vector<int> tier1TrueParticleMinHit;
        for (int i = 0; i < array_row_size; i++)
        {
            for (int j = 0; j < array_col_size; j++)
            {
                std::vector<int> v1;
                std::vector<const recob::Hit*> trueHits = fProtoDUNETruthUtils.GetMCParticleHits(clockData, *tier1TrueParticlesQ1[i], e, fHitTag);
//                std::vector<const recob::Hit*> recoHits = fProtoDUNEPFParticleUtils.GetPFParticleHits(*tier1RecoParticles[j], e, fPFParticleLabel);
                std::vector< const recob::Hit* > hit_matched_vector = fProtoDUNETruthUtils.GetSharedHits(clockData, *tier1TrueParticlesQ1[i], *tier1RecoParticles[j], e, fPFParticleLabel);
//                float purity = fProtoDUNETruthUtils.GetPurity(clockData, *tier1RecoParticles[j], e, fPFParticleLabel);
//                float completeness = fProtoDUNETruthUtils.GetCompleteness(clockData, *tier1RecoParticles[j], e, "pandora", "hitpdune");
                int hit_matched = hit_matched_vector.size();
              
//                if ( (hit_matched >= 10) && (purity >= 0.1) && (completeness >= 0.5) )
                if (trueHits.size() > 14)
                {
                    v1 = {i, j, hit_matched, 0}; 
                    twoDTrueRecoSharedHitsMatrix.push_back(v1);
                    tier1TrueParticleMinHit.push_back(i);
                }            
            }
            
        }

        int row = twoDTrueRecoSharedHitsMatrix.size();

//        std::cout << "tier1TrueParticleMinHit size: " << tier1TrueParticleMinHit.size() << std::endl;
//        int col = twoDTrueRecoSharedHitsMatrix[0].size();
        std::vector< std::vector<int> > twoDTrueRecoSharedHitsSortedMatchedMatrix;

//       for (int i = 0; i < row; i++)
//            std::cout << "shared hits: " << twoDTrueRecoSharedHitsMatrix[i][0] << " , " << twoDTrueRecoSharedHitsMatrix[i][1] << " , " << twoDTrueRecoSharedHitsMatrix[i][2] << std::endl;

        std::sort(twoDTrueRecoSharedHitsMatrix.begin(), twoDTrueRecoSharedHitsMatrix.end(), [] ( const std::vector<int>& v1, const std::vector<int>& v2 )->bool
        {
            return v1[2] > v2[2];
        }); //<========twoDTrueRecoSharedHitsMatrix is the vector<vector<int>>

/*        for (int i = 0; i < row; i++)
            std::cout << "SORTED shared hits: " << twoDTrueRecoSharedHitsMatrix[i][0] << " , " << twoDTrueRecoSharedHitsMatrix[i][1] << " , " << twoDTrueRecoSharedHitsMatrix[i][2] << " , " << twoDTrueRecoSharedHitsMatrix[i][3] << std::endl;*/

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
//            std::cout << "i: " << i << ", true_count: " << true_count << ", reco_count: " << reco_count << std::endl;
//            float testPurity = fProtoDUNETruthUtils.GetPurity(clockData, *tier1RecoParticles[twoDTrueRecoSharedHitsMatrix[i][1]], e, fPFParticleLabel);
//            float testCompleteness = fProtoDUNETruthUtils.GetCompleteness(clockData, *tier1RecoParticles[twoDTrueRecoSharedHitsMatrix[i][1]], e, "pandora", "hitpdune");

            std::vector<const recob::Hit*> trueTempHits = fProtoDUNETruthUtils.GetMCParticleHits(clockData, *tier1TrueParticlesQ1[twoDTrueRecoSharedHitsMatrix[i][0]], e, fHitTag);
            std::vector<const recob::Hit*> recoTempHits = fProtoDUNEPFParticleUtils.GetPFParticleHits(*tier1RecoParticles[twoDTrueRecoSharedHitsMatrix[i][1]], e, fPFParticleLabel);

            float trueTempHitsSize = trueTempHits.size();
            float recoTempHitsSize = recoTempHits.size();

            float testPurity = (twoDTrueRecoSharedHitsMatrix[i][2])/(recoTempHitsSize);
            float testCompleteness = (twoDTrueRecoSharedHitsMatrix[i][2])/(trueTempHitsSize);

//            std::cout << "testPurity: " << testPurity << std::endl;
//            std::cout << "testCompleteness: " << testCompleteness << std::endl;

            if ( (true_count == 0) && (reco_count == 0) && (testPurity >= 0.1) && (testCompleteness >= 0.5) )
            {
                taken_true_particle_index.push_back(twoDTrueRecoSharedHitsMatrix[i][0]);
                taken_reco_particle_index.push_back(twoDTrueRecoSharedHitsMatrix[i][1]);
//                std::cout << "first pair inserted: " << twoDTrueRecoSharedHitsMatrix[i][0] << ", " << twoDTrueRecoSharedHitsMatrix[i][1] << ", " << twoDTrueRecoSharedHitsMatrix[i][2] << ", " << "1" << std::endl; 
                std::vector<int> tempVec = {twoDTrueRecoSharedHitsMatrix[i][0], twoDTrueRecoSharedHitsMatrix[i][1], twoDTrueRecoSharedHitsMatrix[i][2], 1};
                twoDTrueRecoSharedHitsSortedMatchedMatrix.push_back(tempVec);
            }

            else
            {
//                std::cout << "second pair inserted: " << twoDTrueRecoSharedHitsMatrix[i][0] << ", " << twoDTrueRecoSharedHitsMatrix[i][1] << ", " << twoDTrueRecoSharedHitsMatrix[i][2] << ", " << twoDTrueRecoSharedHitsMatrix[i][3] << std::endl; 
                std::vector<int> tempVec2 = {twoDTrueRecoSharedHitsMatrix[i][0], twoDTrueRecoSharedHitsMatrix[i][1], twoDTrueRecoSharedHitsMatrix[i][2], twoDTrueRecoSharedHitsMatrix[i][3]};
                twoDTrueRecoSharedHitsSortedMatchedMatrix.push_back(tempVec2);
            }
//            std::cout << "taken_true_particle_index size: " << taken_true_particle_index.size() << std::endl;
//            std::cout << "taken_reco_particle_index size: " << taken_reco_particle_index.size() << std::endl;

//            std::cout << "true_count: " << true_count << std::endl;
//            std::cout << "reco_count: " << reco_count << std::endl;

        }
        
/*z        for (int i = 0; i < row; i++)
            std::cout << "shared hits: " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][1] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][2] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][3] << std::endl;

        std::cout << "-------done---------" << std::endl;*/

        std::sort(twoDTrueRecoSharedHitsSortedMatchedMatrix.begin(), twoDTrueRecoSharedHitsSortedMatchedMatrix.end(), [] ( const std::vector<int>& v1, const std::vector<int>& v2 )->bool
        {
            return v1[3] > v2[3];
        });        

/*        for (int i = 0; i < row; i++)
            std::cout << "SORTED shared hits: " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][1] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][2] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][3] << std::endl;
        
        std::cout << "-------done---------" << std::endl;*/

        std::vector<int> unique_true_particles;
        std::vector< std::vector<int> > finalMatchedMatrix;

        for (int i = 0; i < row; i++)
        {
            int count = std::count (unique_true_particles.begin(), unique_true_particles.end(), twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0]);

            if (count == 0)
            {
                unique_true_particles.push_back(twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0]);
                std::vector<int> tempVec3 = {twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0], twoDTrueRecoSharedHitsSortedMatchedMatrix[i][1], twoDTrueRecoSharedHitsSortedMatchedMatrix[i][2], twoDTrueRecoSharedHitsSortedMatchedMatrix[i][3]};
                finalMatchedMatrix.push_back(tempVec3);
            }
        }

/*        for (int i = 0; i < row; i++)
            std::cout << "SORTED unique shared hits: " << finalMatchedMatrix[i][0] << " , " << finalMatchedMatrix[i][1] << " , " << finalMatchedMatrix[i][2] << " , " << finalMatchedMatrix[i][3] << std::endl;
        
        std::cout << "-------done---------" << std::endl;*/
        




//        for (int i = 0; i < row; i++)
//            std::cout << "SORTED matched shared hits: " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][0] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][1] << " , " << twoDTrueRecoSharedHitsSortedMatchedMatrix[i][2] << std::endl;

//                fpurity = fProtoDUNETruthUtils.GetPurity(clockData, *beamParticles[0], e, fPFParticleLabel); 
//        fcompleteness = fProtoDUNETruthUtils.GetCompleteness(clockData, *beamParticles[0], e, "pandora", "hitpdune");
//array_col_size

//        std::cout << "1" << std::endl;
        int matched_row = finalMatchedMatrix.size();
        fTier1MCRecoMatchedNumber = matched_row;
    //        std::cout << "matched_row size: " << matched_row << std::endl;
        for (int i = 0; i < matched_row; i++) //todo check matched_row
        {
            if (finalMatchedMatrix[i][3] == 1)
            {
                fTier1RecoMCParticleMatch = 1;
                auto trueParticle = tier1TrueParticlesQ1[finalMatchedMatrix[i][0]];
                auto recoParticle = tier1RecoParticles[finalMatchedMatrix[i][1]];
                fMCRecoMatchedHits = finalMatchedMatrix[i][2];
                fTier1MCParticleHitsSize = fProtoDUNETruthUtils.GetMCParticleHits( clockData, *trueParticle, e, fHitTag ).size();
                fTier1RecoBeamParticleNhits = fProtoDUNEPFParticleUtils.GetPFParticleHits( *recoParticle, e, fPFParticleLabel).size();
                auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    //--------------------------------true tier 1 quantities--------------------------------------------------------------------------
                fTier1TrueBeamParticleLength     = trueParticle->Trajectory().TotalLength();
                fTier1TrueBeamParticleEndX       = trueParticle->EndX();
                fTier1TrueBeamParticleEndY       = trueParticle->EndY();
                fTier1TrueBeamParticleEndZ       = trueParticle->EndZ();
                fTier1MCPDGCode    = trueParticle->PdgCode();
                fTier1TrueBeamParticleStartX     = trueParticle->Position(0).X();
                fTier1TrueBeamParticleStartY     = trueParticle->Position(0).Y();
                fTier1TrueBeamParticleStartZ     = trueParticle->Position(0).Z();
                fTier1TrueBeamParticleStartPx    = trueParticle->Px();
                fTier1TrueBeamParticleStartPy    = trueParticle->Py();
                fTier1TrueBeamParticleStartPz    = trueParticle->Pz();
                fTier1MCStartDirectionX          = (fTier1TrueBeamParticleStartPx)/(trueParticle->P());
                fTier1MCStartDirectionY          = (fTier1TrueBeamParticleStartPy)/(trueParticle->P());
                fTier1MCStartDirectionZ          = (fTier1TrueBeamParticleStartPz)/(trueParticle->P());

    //            tier1beam_inst_X = trueParticle->Trajectory().X(trueParticle->Trajectory().size());
    //            tier1beam_inst_Y = trueParticle->Trajectory().Y(trueParticle->Trajectory().size());
    //            tier1beam_inst_Z = trueParticle->Trajectory().Z(trueParticle->Trajectory().size());

//                tier1beam_inst_X = fTier1TrueBeamParticleStartX;
//                tier1beam_inst_Y = fTier1TrueBeamParticleStartY;
//                tier1beam_inst_Z = fTier1TrueBeamParticleStartZ;

//                fTier1BeamInst_startVertex_X_SCE_corrected = tier1beam_inst_X-SCE->GetPosOffsets(geo::Point_t(tier1beam_inst_X,tier1beam_inst_Y,tier1beam_inst_Z)).X();
//                fTier1BeamInst_startVertex_Y_SCE_corrected = tier1beam_inst_Y+SCE->GetPosOffsets(geo::Point_t(tier1beam_inst_X,tier1beam_inst_Y,tier1beam_inst_Z)).Y();
//                fTier1BeamInst_startVertex_Z_SCE_corrected = tier1beam_inst_Z+SCE->GetPosOffsets(geo::Point_t(tier1beam_inst_X,tier1beam_inst_Y,tier1beam_inst_Z)).Z();
            

//                float tier1true_beam_particle_SCE_end_X = fTier1TrueBeamParticleEndX-SCE->GetPosOffsets(geo::Point_t(fTier1TrueBeamParticleEndX,fTier1TrueBeamParticleEndY,fTier1TrueBeamParticleEndZ)).X();
//                float tier1true_beam_particle_SCE_end_Y = fTier1TrueBeamParticleEndY+SCE->GetPosOffsets(geo::Point_t(fTier1TrueBeamParticleEndX,fTier1TrueBeamParticleEndY,fTier1TrueBeamParticleEndZ)).Y();
//                float tier1true_beam_particle_SCE_end_Z = fTier1TrueBeamParticleEndZ+SCE->GetPosOffsets(geo::Point_t(fTier1TrueBeamParticleEndX,fTier1TrueBeamParticleEndY,fTier1TrueBeamParticleEndZ)).Z();

                fTier1TrueBeamParticleInteractionX = fTier1TrueBeamParticleEndX;       
                fTier1TrueBeamParticleInteractionY = fTier1TrueBeamParticleEndY;       
                fTier1TrueBeamParticleInteractionZ = fTier1TrueBeamParticleEndZ;       

                fTier1MCLengthByTrajPoints = sqrt(((fTier1TrueBeamParticleEndX-fTier1TrueBeamParticleStartX)*(fTier1TrueBeamParticleEndX-fTier1TrueBeamParticleStartX)) + ((fTier1TrueBeamParticleEndY-fTier1TrueBeamParticleStartY)*(fTier1TrueBeamParticleEndY-fTier1TrueBeamParticleStartY)) + ((fTier1TrueBeamParticleEndZ-fTier1TrueBeamParticleStartZ)*(fTier1TrueBeamParticleEndZ-fTier1TrueBeamParticleStartZ)));




//        fTier0MCLengthByTrajPoints = sqrt(((fTier0MCInteractionVertexX-fTier0MCStartVertexX)*(fTier0MCInteractionVertexX-fTier0MCStartVertexX)) + ((fTier0MCInteractionVertexY-fTier0MCStartVertexY)*(fTier0MCInteractionVertexY-fTier0MCStartVertexY)) + ((fTier0MCInteractionVertexZ-fTier0MCStartVertexZ)*(fTier0MCInteractionVertexZ-fTier0MCStartVertexZ)));
    //--------------------------------reco tier 1 quantities--------------------------------------------------------------------------
                fTier1Purity = fProtoDUNETruthUtils.GetPurity(clockData, *recoParticle, e, fPFParticleLabel);
                fTier1Completeness = fProtoDUNETruthUtils.GetCompleteness(clockData, *recoParticle, e, "pandora", "hitpdune");
                
                const recob::Track* thisTrack = fProtoDUNEPFParticleUtils.GetPFParticleTrack(*recoParticle,e,fPFParticleLabel,fTrackLabel); 
                const recob::Shower* thisShower = fProtoDUNEPFParticleUtils.GetPFParticleShower(*recoParticle,e,fPFParticleLabel,fShowerLabel);

                if(thisTrack != 0x0)
                {
                    fTier1Reco_beam_pfp_topology = 1;
                    fTier1RecoLength = thisTrack->Length();
                    auto calo = fProtoDUNETrackUtils.GetRecoTrackCalorimetry(*thisTrack, e, fTrackLabel, fCalorimetryTagSCE);
                    const TVector3 reco_primary_start_vertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(*recoParticle,e,fPFParticleLabel,fTrackLabel);

                    fTier1RecoBeamParticleStartX = reco_primary_start_vertex.X();
                    fTier1RecoBeamParticleStartY = reco_primary_start_vertex.Y();
                    fTier1RecoBeamParticleStartZ = reco_primary_start_vertex.Z();

                    fTier1RecoStartVertexX_SCE_corrected = fTier1RecoBeamParticleStartX+SCE->GetPosOffsets(geo::Point_t(fTier1RecoBeamParticleStartX,fTier1RecoBeamParticleStartY,fTier1RecoBeamParticleStartZ)).X();
                    fTier1RecoStartVertexY_SCE_corrected = fTier1RecoBeamParticleStartY-SCE->GetPosOffsets(geo::Point_t(fTier1RecoBeamParticleStartX,fTier1RecoBeamParticleStartY,fTier1RecoBeamParticleStartZ)).Y();
                    fTier1RecoStartVertexZ_SCE_corrected = fTier1RecoBeamParticleStartZ-SCE->GetPosOffsets(geo::Point_t(fTier1RecoBeamParticleStartX,fTier1RecoBeamParticleStartY,fTier1RecoBeamParticleStartZ)).Z();

                    const TVector3 reco_primary_interaction_vertex = fProtoDUNEPFParticleUtils.GetPFParticleSecondaryVertex(*recoParticle,e,fPFParticleLabel,fTrackLabel);

                    fTier1RecoBeamParticleInteractionX = reco_primary_interaction_vertex.X();
                    fTier1RecoBeamParticleInteractionY = reco_primary_interaction_vertex.Y();
                    fTier1RecoBeamParticleInteractionZ = reco_primary_interaction_vertex.Z();

//                    fTier1RecoBeamStartDirX = (thisTrack->StartDirection()).X();
//                    fTier1RecoBeamStartDirY = (thisTrack->StartDirection()).Y();
//                    fTier1RecoBeamStartDirZ = (thisTrack->StartDirection()).Z();
//-------------------------------------------------------------------------------------------------------------------

            bool found_calo = false;
            size_t index = 0;
            for ( index = 0; index < calo.size(); ++index)
            {
                if (calo[index].PlaneID().Plane == 2)
                {
                    found_calo = true;
                    break; 
                }
            }

            if (found_calo)
            {
                auto theXYZPoints = calo[index].XYZ();
                std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b){return (a.Z() < b.Z());});

                if (theXYZPoints.size())
                {
                    TVector3 dir((theXYZPoints.back().X() - theXYZPoints[0].X()), (theXYZPoints.back().Y() - theXYZPoints[0].Y()), (theXYZPoints.back().Z() - theXYZPoints[0].Z()));
                    fTier1RecoStartDirectionX_SCE_corrected = dir.Unit().X();
                    fTier1RecoStartDirectionY_SCE_corrected = dir.Unit().Y();
                    fTier1RecoStartDirectionZ_SCE_corrected = dir.Unit().Z();
                }

                else
                {
                    fTier1RecoStartDirectionX_SCE_corrected = -9999;
                    fTier1RecoStartDirectionY_SCE_corrected = -9999;
                    fTier1RecoStartDirectionZ_SCE_corrected = -9999;
                }

                if (theXYZPoints.size() > 1)
                {
                    TVector3 start_p1(theXYZPoints[0].X(), theXYZPoints[0].Y(), theXYZPoints[0].Z());
                    TVector3 start_p2(theXYZPoints[1].X(), theXYZPoints[1].Y(), theXYZPoints[1].Z());
                    TVector3 start_diff = start_p2 - start_p1;

                    fTier1RecoStartDirectionX_SCE_corrected = start_diff.Unit().X();
                    fTier1RecoStartDirectionY_SCE_corrected = start_diff.Unit().Y();
                    fTier1RecoStartDirectionZ_SCE_corrected = start_diff.Unit().Z();
                }

                else
                {
                    fTier1RecoStartDirectionX_SCE_corrected = -9999;
                    fTier1RecoStartDirectionY_SCE_corrected = -9999;
                    fTier1RecoStartDirectionZ_SCE_corrected = -9999;
                }

/*                if (theXYZPoints.size() > 2)
                {
                    std::vector<TVector3> input;
                    for (size_t iP = 0; iP < 3; ++iP) 
                    {
                        input.push_back(TVector3(theXYZPoints[iP].X(), theXYZPoints[iP].Y(), theXYZPoints[iP].Z()));
                    }

                    TVector3 startDiff = FitLine(input);

                    fTier1RecoStartDirectionX_SCE_corrected = startDiff.Unit().X();
                    fTier0RecoStartDirectionY_SCE_corrected = startDiff.Unit().Y();
                    fTier0RecoStartDirectionZ_SCE_corrected = startDiff.Unit().Z();
                }
                
                else
                {
                    fTier1RecoStartDirectionX_SCE_corrected = -9999;
                    fTier0RecoStartDirectionY_SCE_corrected = -9999;
                    fTier0RecoStartDirectionZ_SCE_corrected = -9999;
                }*/
            }
//---------------------------------------------------------------------------------------------------------------------

                }

                if(thisShower != 0x0) 
                {
                    fTier1Reco_beam_pfp_topology = 0;
                    fTier1RecoLength = thisShower->Length();
                    auto calo = fProtoDUNEShowerUtils.GetRecoShowerCalorimetry(*thisShower, e, fShowerLabel, fCalorimetryTagSCE);
                    const TVector3 reco_primary_start_vertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(*beamParticles[0],e,fPFParticleLabel,fTrackLabel);

                    fTier1RecoBeamParticleStartX = reco_primary_start_vertex.X();
                    fTier1RecoBeamParticleStartY = reco_primary_start_vertex.Y();
                    fTier1RecoBeamParticleStartZ = reco_primary_start_vertex.Z();

                    fTier1RecoStartVertexX_SCE_corrected = fTier1RecoBeamParticleStartX+SCE->GetPosOffsets(geo::Point_t(fTier1RecoBeamParticleStartX,fTier1RecoBeamParticleStartY,fTier1RecoBeamParticleStartZ)).X();
                    fTier1RecoStartVertexY_SCE_corrected = fTier1RecoBeamParticleStartY-SCE->GetPosOffsets(geo::Point_t(fTier1RecoBeamParticleStartX,fTier1RecoBeamParticleStartY,fTier1RecoBeamParticleStartZ)).Y();
                    fTier1RecoStartVertexZ_SCE_corrected = fTier1RecoBeamParticleStartZ-SCE->GetPosOffsets(geo::Point_t(fTier1RecoBeamParticleStartX,fTier1RecoBeamParticleStartY,fTier1RecoBeamParticleStartZ)).Z();

//                    fTier1RecoBeamStartDirX = (thisShower->Direction()).X();
//                    fTier1RecoBeamStartDirY = (thisShower->Direction()).Y();
//                    fTier1RecoBeamStartDirZ = (thisShower->Direction()).Z();

//-------------------------------------------------------------------------------------------------------------------

            bool found_calo = false;
            size_t index = 0;
            for ( index = 0; index < calo.size(); ++index)
            {
                if (calo[index].PlaneID().Plane == 2)
                {
                    found_calo = true;
                    break; 
                }
            }

            if (found_calo)
            {
                auto theXYZPoints = calo[index].XYZ();
                std::sort(theXYZPoints.begin(), theXYZPoints.end(), [](auto a, auto b){return (a.Z() < b.Z());});

                if (theXYZPoints.size())
                {
                    TVector3 dir((theXYZPoints.back().X() - theXYZPoints[0].X()), (theXYZPoints.back().Y() - theXYZPoints[0].Y()), (theXYZPoints.back().Z() - theXYZPoints[0].Z()));
                    fTier1RecoStartDirectionX_SCE_corrected = dir.Unit().X();
                    fTier1RecoStartDirectionY_SCE_corrected = dir.Unit().Y();
                    fTier1RecoStartDirectionZ_SCE_corrected = dir.Unit().Z();
                }

                else
                {
                    fTier1RecoStartDirectionX_SCE_corrected = -9999;
                    fTier1RecoStartDirectionY_SCE_corrected = -9999;
                    fTier1RecoStartDirectionZ_SCE_corrected = -9999;
                }

                if (theXYZPoints.size() > 1)
                {
                    TVector3 start_p1(theXYZPoints[0].X(), theXYZPoints[0].Y(), theXYZPoints[0].Z());
                    TVector3 start_p2(theXYZPoints[1].X(), theXYZPoints[1].Y(), theXYZPoints[1].Z());
                    TVector3 start_diff = start_p2 - start_p1;

                    fTier1RecoStartDirectionX_SCE_corrected = start_diff.Unit().X();
                    fTier1RecoStartDirectionY_SCE_corrected = start_diff.Unit().Y();
                    fTier1RecoStartDirectionZ_SCE_corrected = start_diff.Unit().Z();
                }

                else
                {
                    fTier1RecoStartDirectionX_SCE_corrected = -9999;
                    fTier1RecoStartDirectionY_SCE_corrected = -9999;
                    fTier1RecoStartDirectionZ_SCE_corrected = -9999;
                }

/*                if (theXYZPoints.size() > 2)
                {
                    std::vector<TVector3> input;
                    for (size_t iP = 0; iP < 3; ++iP) 
                    {
                        input.push_back(TVector3(theXYZPoints[iP].X(), theXYZPoints[iP].Y(), theXYZPoints[iP].Z()));
                    }

                    TVector3 startDiff = FitLine(input);

                    fTier1RecoStartDirectionX_SCE_corrected = startDiff.Unit().X();
                    fTier0RecoStartDirectionY_SCE_corrected = startDiff.Unit().Y();
                    fTier0RecoStartDirectionZ_SCE_corrected = startDiff.Unit().Z();
                }
                
                else
                {
                    fTier1RecoStartDirectionX_SCE_corrected = -9999;
                    fTier0RecoStartDirectionY_SCE_corrected = -9999;
                    fTier0RecoStartDirectionZ_SCE_corrected = -9999;
                }*/
            }
//---------------------------------------------------------------------------------------------------------------------

                }

            float tier1StartVertexDx = fTier1RecoStartVertexX_SCE_corrected - fTier0MCInteractionVertexX;
            float tier1StartVertexDy = fTier1RecoStartVertexY_SCE_corrected - fTier0MCInteractionVertexY;
            float tier1StartVertexDz = fTier1RecoStartVertexZ_SCE_corrected - fTier0MCInteractionVertexZ;

            fTier1StartVertexDr = sqrt(((tier1StartVertexDx)*(tier1StartVertexDx))+((tier1StartVertexDy)*(tier1StartVertexDy))+((tier1StartVertexDz)*(tier1StartVertexDz))); 
            fTestTree->Fill();
            }
        
            else if (finalMatchedMatrix[i][3] == 0)
            {
                fTier1RecoMCParticleMatch = 0;
                auto trueParticle = tier1TrueParticlesQ1[finalMatchedMatrix[i][0]];
                fTier1MCParticleHitsSize = fProtoDUNETruthUtils.GetMCParticleHits( clockData, *trueParticle, e, fHitTag ).size();
                fTier1MCPDGCode = trueParticle->PdgCode();
                fTestTree->Fill();
            }
        }
    }
}

void analysis::PDSPKaonAnalysis::beginJob()
{
  // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("tree","Analyser Output Tree");
    fTestTree = tfs->make<TTree>("testtree","Analyser Output Tree");

  // Tier0 Tree Branches    
    fTree->Branch("fEventID",&fEventID,"fEventID/I");
    fTree->Branch("fTier0MCPDGCode",&fTier0MCPDGCode,"fTier0MCPDGCode/I");
    fTree->Branch("fTier0Completeness",&fTier0Completeness,"fTier0Completeness/F");//
    fTree->Branch("fTier0Purity",&fTier0Purity,"fTier0Purity/F");//
    fTree->Branch("fTier0MCStartVertexX",&fTier0MCStartVertexX,"fTier0MCStartVertexX/F");
    fTree->Branch("fTier0MCStartVertexY",&fTier0MCStartVertexY,"fTier0MCStartVertexY/F");
    fTree->Branch("fTier0MCStartVertexZ",&fTier0MCStartVertexZ,"fTier0MCStartVertexZ/F");
    fTree->Branch("fTier0RecoStartVertexX_SCE_corrected",&fTier0RecoStartVertexX_SCE_corrected,"fTier0RecoStartVertexX_SCE_corrected/F");//
    fTree->Branch("fTier0RecoStartVertexY_SCE_corrected",&fTier0RecoStartVertexY_SCE_corrected,"fTier0RecoStartVertexY_SCE_corrected/F");//
    fTree->Branch("fTier0RecoStartVertexZ_SCE_corrected",&fTier0RecoStartVertexZ_SCE_corrected,"fTier0RecoStartVertexZ_SCE_corrected/F");//
    fTree->Branch("fTier0StartVertexDr",&fTier0StartVertexDr,"fTier0StartVertexDr/F");//
    fTree->Branch("fTier0MCInteractionVertexX",&fTier0MCInteractionVertexX,"fTier0MCInteractionVertexX/F");//
    fTree->Branch("fTier0MCInteractionVertexY",&fTier0MCInteractionVertexY,"fTier0MCInteractionVertexY/F");//
    fTree->Branch("fTier0MCInteractionVertexZ",&fTier0MCInteractionVertexZ,"fTier0MCInteractionVertexZ/F");//
    fTree->Branch("fTier0RecoInteractionVertexX_SCE_corrected",&fTier0RecoInteractionVertexX_SCE_corrected,"fTier0RecoInteractionVertexX_SCE_corrected/F");//
    fTree->Branch("fTier0RecoInteractionVertexY_SCE_corrected",&fTier0RecoInteractionVertexY_SCE_corrected,"fTier0RecoInteractionVertexY_SCE_corrected/F");//
    fTree->Branch("fTier0RecoInteractionVertexZ_SCE_corrected",&fTier0RecoInteractionVertexZ_SCE_corrected,"fTier0RecoInteractionVertexZ_SCE_corrected/F");//
    fTree->Branch("fTier0InteractionVertexDr",&fTier0InteractionVertexDr,"fTier0InteractionVertexDr/F");//
    fTree->Branch("fTier0MCStartDirectionX",&fTier0MCStartDirectionX,"fTier0MCStartDirectionX/F");//
    fTree->Branch("fTier0MCStartDirectionY",&fTier0MCStartDirectionY,"fTier0MCStartDirectionY/F");//
    fTree->Branch("fTier0MCStartDirectionZ",&fTier0MCStartDirectionZ,"fTier0MCStartDirectionZ/F");//
    fTree->Branch("fTier0RecoStartDirectionX_SCE_corrected",&fTier0RecoStartDirectionX_SCE_corrected,"fTier0RecoStartDirectionX_SCE_corrected/F");//
    fTree->Branch("fTier0RecoStartDirectionY_SCE_corrected",&fTier0RecoStartDirectionY_SCE_corrected,"fTier0RecoStartDirectionY_SCE_corrected/F");//
    fTree->Branch("fTier0RecoStartDirectionZ_SCE_corrected",&fTier0RecoStartDirectionZ_SCE_corrected,"fTier0RecoStartDirectionZ_SCE_corrected/F");//
    fTree->Branch("fTier0MCLengthByTrajPoints",&fTier0MCLengthByTrajPoints,"fTier0MCLengthByTrajPoints/F");//
    fTree->Branch("fTier0RecoLength",&fTier0RecoLength,"fTier0RecoLength/F");//
    fTree->Branch("fTier0MCParticleHitsSize",&fTier0MCParticleHitsSize,"fTier0MCParticleHitsSize/I");//
    fTree->Branch("fTier0RecoMCParticleMatch",&fTier0RecoMCParticleMatch,"fTier0RecoMCParticleMatch/I");//

//    fTree->Branch("doesTrueBeamParticleExist",&fDoesTrueBeamParticleExist,"doesTrueBeamParticleExist/I");
//    fTree->Branch("trueBeamParticleStartX",&fTrueBeamParticleStartX,"trueBeamParticleStartX/F");
//    fTree->Branch("trueBeamParticleStartY",&fTrueBeamParticleStartY,"trueBeamParticleStartY/F");
//    fTree->Branch("trueBeamParticleStartZ",&fTrueBeamParticleStartZ,"trueBeamParticleStartZ/F");

//    fTree->Branch("trueBeamParticleEndX",&fTrueBeamParticleEndX,"trueBeamParticleEndX/F");
//    fTree->Branch("trueBeamParticleEndY",&fTrueBeamParticleEndY,"trueBeamParticleEndY/F");
//    fTree->Branch("trueBeamParticleEndZ",&fTrueBeamParticleEndZ,"trueBeamParticleEndZ/F");




//    fTree->Branch("trueBeamParticleStartPx",&fTrueBeamParticleStartPx,"trueBeamParticleStartPx/F");
//    fTree->Branch("trueBeamParticleStartPy",&fTrueBeamParticleStartPy,"trueBeamParticleStartPy/F");
//    fTree->Branch("trueBeamParticleStartPz",&fTrueBeamParticleStartPz,"trueBeamParticleStartPz/F");

  // Reco Tree Branches   
//    fTree->Branch("bestMatchedMCParticleFromPFParticlePdgCode",&fbestMatchedMCParticleFromPFParticlePdgCode,"bestMatchedMCParticleFromPFParticlePdgCode/I");
//    fTree->Branch("doesRecoBeamParticleExist",&fDoesRecoBeamParticleExist,"doesRecoBeamParticleExist/I");


//    fTree->Branch("recoBeamParticleStartPx",&fRecoBeamParticleStartPx,"recoBeamParticleStartPx/F");
//    fTree->Branch("recoBeamParticleStartPy",&fRecoBeamParticleStartPy,"recoBeamParticleStartPy/F");
//    fTree->Branch("recoBeamParticleStartPz",&fRecoBeamParticleStartPz,"recoBeamParticleStartPz/F");
//

//    fTree->Branch("numberOfReconstructedBeamParticle",&fNumberOfReconstructedBeamParticle,"numberOfReconstructedBeamParticle/I");
//    fTree->Branch("primaryBeamParticleLength",&fPrimaryBeamParticleLength,"primaryBeamParticleLength/F");
//    fTree->Branch("trueBeamParticleNhits", &fTrueBeamParticleNhits, "trueBeamParticleNhits/I");

  // Interested quantities

//    fTree->Branch("reco_beam_pfp_topology",&fReco_beam_pfp_topology,"reco_beam_pfp_topology/I");//




 //   fTree->Branch("recoBeamParticleStartX",&fRecoBeamParticleStartX,"recoBeamParticleStartX/F");//
//    fTree->Branch("recoBeamParticleStartY",&fRecoBeamParticleStartY,"recoBeamParticleStartY/F");//
//fTree->Branch("recoBeamParticleStartZ",&fRecoBeamParticleStartZ,"recoBeamParticleStartZ/F");//

//    fTree->Branch("beamInst_startVertex_X_SCE_corrected",&fBeamInst_startVertex_X_SCE_corrected,"beamInst_startVertex_X_SCE_corrected/F");//
//    fTree->Branch("beamInst_startVertex_Y_SCE_corrected",&fBeamInst_startVertex_Y_SCE_corrected,"beamInst_startVertex_Y_SCE_corrected/F");//
//    fTree->Branch("beamInst_startVertex_Z_SCE_corrected",&fBeamInst_startVertex_Z_SCE_corrected,"beamInst_startVertex_Z_SCE_corrected/F");//
//
//    fTree->Branch("recoBeamParticleInteractionX",&fRecoBeamParticleInteractionX,"recoBeamParticleInteractionX/F");//
//    fTree->Branch("recoBeamParticleInteractionY",&fRecoBeamParticleInteractionY,"recoBeamParticleInteractionY/F");//
//    fTree->Branch("recoBeamParticleInteractionZ",&fRecoBeamParticleInteractionZ,"recoBeamParticleInteractionZ/F");//


//
//    fTree->Branch("recoBeamParticleNhits", &fRecoBeamParticleNhits, "recoBeamParticleNhits/I");

//    fTree->Branch("trueBeamStartDirX",&fTrueBeamStartDirX,"trueBeamStartDirX/F");//
//    fTree->Branch("trueBeamStartDirY",&fTrueBeamStartDirY,"trueBeamStartDirY/F");//
//    fTree->Branch("trueBeamStartDirZ",&fTrueBeamStartDirZ,"trueBeamStartDirZ/F");//

//    fTree->Branch("recoBeamStartDirX",&fRecoBeamStartDirX,"recoBeamStartDirX/F");//
//    fTree->Branch("recoBeamStartDirY",&fRecoBeamStartDirY,"recoBeamStartDirY/F");//
//    fTree->Branch("recoBeamStartDirZ",&fRecoBeamStartDirZ,"recoBeamStartDirZ/F");//





//    fTree->Branch("fTier0recoBeamParticleHitsSize",&fTier0recoBeamParticleHitsSize,"fTier0recoBeamParticleHitsSize/F");//
//    fTree->Branch("fSharedTier0RecoTrueHitsSize",&fSharedTier0RecoTrueHitsSize,"fSharedTier0RecoTrueHitsSize/F");//

    fTestTree->Branch("fTier1MCPDGCode",&fTier1MCPDGCode,"fTier1MCPDGCode/I");//
    fTestTree->Branch("tier1Completeness",&fTier1Completeness,"tier1Completeness/F");//
    fTestTree->Branch("tier1Purity",&fTier1Purity,"tier1Purity/F");//
    fTestTree->Branch("fTier1RecoStartVertexX_SCE_corrected",&fTier1RecoStartVertexX_SCE_corrected,"fTier1RecoStartVertexX_SCE_corrected/F");
    fTestTree->Branch("fTier1RecoStartVertexY_SCE_corrected",&fTier1RecoStartVertexY_SCE_corrected,"fTier1RecoStartVertexY_SCE_corrected/F");
    fTestTree->Branch("fTier1RecoStartVertexZ_SCE_corrected",&fTier1RecoStartVertexZ_SCE_corrected,"fTier1RecoStartVertexZ_SCE_corrected/F");
    fTestTree->Branch("fTier1StartVertexDr",&fTier1StartVertexDr,"fTier1StartVertexDr/F");
    fTestTree->Branch("fTier1MCStartDirectionX",&fTier1MCStartDirectionX,"fTier1MCStartDirectionX/F");
    fTestTree->Branch("fTier1MCStartDirectionY",&fTier1MCStartDirectionY,"fTier1MCStartDirectionY/F");
    fTestTree->Branch("fTier1MCStartDirectionZ",&fTier1MCStartDirectionZ,"fTier1MCStartDirectionZ/F");
    fTestTree->Branch("fTier1RecoStartDirectionX_SCE_corrected",&fTier1RecoStartDirectionX_SCE_corrected,"fTier1RecoStartDirectionX_SCE_corrected/F");
    fTestTree->Branch("fTier1RecoStartDirectionY_SCE_corrected",&fTier1RecoStartDirectionY_SCE_corrected,"fTier1RecoStartDirectionY_SCE_corrected/F");
    fTestTree->Branch("fTier1RecoStartDirectionZ_SCE_corrected",&fTier1RecoStartDirectionZ_SCE_corrected,"fTier1RecoStartDirectionZ_SCE_corrected/F");
    fTestTree->Branch("fTier1MCLengthByTrajPoints",&fTier1MCLengthByTrajPoints,"fTier1MCLengthByTrajPoints/F");
    fTestTree->Branch("fTier1RecoLength",&fTier1RecoLength,"fTier1RecoLength/F");//
    fTestTree->Branch("fTier1MCParticleHitsSize",&fTier1MCParticleHitsSize,"fTier1MCParticleHitsSize/I");//
    fTestTree->Branch("fTier1RecoMCParticleMatch",&fTier1RecoMCParticleMatch,"fTier1RecoMCParticleMatch/I");//

//    fTestTree->Branch("tier1beamInst_startVertex_X_SCE_corrected",&fTier1BeamInst_startVertex_X_SCE_corrected,"tier1beamInst_startVertex_X_SCE_corrected/F");//
//    fTestTree->Branch("tier1beamInst_startVertex_Y_SCE_corrected",&fTier1BeamInst_startVertex_Y_SCE_corrected,"tier1beamInst_startVertex_Y_SCE_corrected/F");//
//    fTestTree->Branch("tier1beamInst_startVertex_Z_SCE_corrected",&fTier1BeamInst_startVertex_Z_SCE_corrected,"tier1beamInst_startVertex_Z_SCE_corrected/F");//
//    fTestTree->Branch("tier1beam_length_by_traj_points",&fTier1Beam_length_by_traj_points,"tier1beam_length_by_traj_points/F");//
//    fTestTree->Branch("tier1recoBeamParticleStartX",&fTier1RecoBeamParticleStartX,"tier1recoBeamParticleStartX/F");//
//    fTestTree->Branch("tier1recoBeamParticleStartY",&fTier1RecoBeamParticleStartY,"tier1recoBeamParticleStartY/F");//
//    fTestTree->Branch("tier1recoBeamParticleStartZ",&fTier1RecoBeamParticleStartZ,"tier1recoBeamParticleStartZ/F");//
//    fTestTree->Branch("tier1recoBeamParticleInteractionX",&fTier1RecoBeamParticleInteractionX,"tier1recoBeamParticleInteractionX/F");//
//    fTestTree->Branch("tier1recoBeamParticleInteractionY",&fTier1RecoBeamParticleInteractionY,"tier1recoBeamParticleInteractionY/F");//
//    fTestTree->Branch("tier1recoBeamParticleInteractionZ",&fTier1RecoBeamParticleInteractionZ,"tier1recoBeamParticleInteractionZ/F");//
//    fTestTree->Branch("tier1trueBeamStartDirX",&fTier1TrueBeamParticleStartPx,"tier1trueBeamStartDirX/F");//
//    fTestTree->Branch("tier1trueBeamStartDirY",&fTier1TrueBeamParticleStartPy,"tier1trueBeamStartDirY/F");//
//    fTestTree->Branch("tier1trueBeamStartDirZ",&fTier1TrueBeamParticleStartPz,"tier1trueBeamStartDirZ/F");//
//    fTestTree->Branch("tier1recoBeamStartDirX",&fTier1RecoBeamStartDirX,"tier1recoBeamStartDirX/F");//
//    fTestTree->Branch("tier1recoBeamStartDirY",&fTier1RecoBeamStartDirY,"tier1recoBeamStartDirY/F");//
//    fTestTree->Branch("tier1recoBeamStartDirZ",&fTier1RecoBeamStartDirZ,"tier1recoBeamStartDirZ/F");//

//    fTestTree->Branch("fTier1MCParticleNumber",&fTier1MCParticleNumber,"fTier1MCParticleNumber/I");//
//    fTestTree->Branch("fTier1RecoParticleNumber",&fTier1RecoParticleNumber,"fTier1RecoParticleNumber/I");//
//    fTestTree->Branch("fTier1MCRecoMatchedNumber",&fTier1MCRecoMatchedNumber,"fTier1MCRecoMatchedNumber/I");//
//    fTestTree->Branch("fTier1TrueBeamParticleLength",&fTier1TrueBeamParticleLength,"fTier1TrueBeamParticleLength/F");//
 //   fTestTree->Branch("fTier1TrueBeamParticleInteractionX",&fTier1TrueBeamParticleInteractionX,"fTier1TrueBeamParticleInteractionX/F");//
////    fTestTree->Branch("fTier1TrueBeamParticleInteractionY",&fTier1TrueBeamParticleInteractionY,"fTier1TrueBeamParticleInteractionY/F");//fix
//    fTestTree->Branch("fTier1TrueBeamParticleInteractionZ",&fTier1TrueBeamParticleInteractionZ,"fTier1TrueBeamParticleInteractionZ/F");//fix
//    fTestTree->Branch("fTier1RecoBeamParticleNhits",&fTier1RecoBeamParticleNhits,"fTier1RecoBeamParticleNhits/I");//
//    fTestTree->Branch("fTier1TrueBeamParticleNhits",&fTier1TrueBeamParticleNhits,"fTier1TrueBeamParticleNhits/I");//
//    fTestTree->Branch("fTier1TrueBeamParticlePDGCode",&fTier1TrueBeamParticlePDGCode,"fTier1TrueBeamParticlePDGCode/I");//
//    fTestTree->Branch("fStartVertexDr",&fStartVertexDr,"fStartVertexDr/F");
//
//    fTestTree->Branch("fTrueBeamLengthVersion3",&fTrueBeamLengthVersion3,"fTrueBeamLengthVersion3/F");
//
//    fTestTree->Branch("fMCRecoMatchedHits",&fMCRecoMatchedHits,"fMCRecoMatchedHits/I");

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
            
//        if ( ((inputVector[i]->PdgCode()) != 111) && ((inputVector[i]->PdgCode()) != 2112) && ((inputVector[i]->PdgCode()) < 1000000000))
        if ( ((inputVector[i]->PdgCode()) != 111) && ((inputVector[i]->PdgCode()) < 1000000000))
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

        if ((inputVector[i]->PdgCode()) == 2112)
            outputVector.pop_back();
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
//        if ( ((inputVector[i]->PdgCode()) == 111) || ((inputVector[i]->PdgCode()) > 1000000000))

    }

//    std::cout << "count: " << count << std::endl;

    return count;
}

/*void analysis::PDSPKaonAnalysis::line(double t, double *p,
                                 double &x, double &y, double &z) {
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

double analysis::PDSPKaonAnalysis::distance2(double x,double y,double z, double *p) {
   // distance line point is D= | (xp-x0) cross  ux |
   // where ux is direction of line and x0 is a point in the line (like t = 0)
   ROOT::Math::XYZVector xp(x,y,z);
   ROOT::Math::XYZVector x0(p[0], p[2], 0.);
   ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1.);
   ROOT::Math::XYZVector u = (x1-x0).Unit();
   double d2 = ((xp-x0).Cross(u)) .Mag2();
   return d2;
}


void analysis::PDSPKaonAnalysis::SumDistance2(int &, double *, double & sum, double * par, int) {
   // the TGraph must be a global variable
   TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );
   assert(gr != 0);
   double * x = gr->GetX();
   double * y = gr->GetY();
   double * z = gr->GetZ();
   int npoints = gr->GetN();
   sum = 0;
   for (int i  = 0; i < npoints; ++i) {
      double d = distance2(x[i],y[i],z[i],par);
      sum += d;
   }
}

TVector3 analysis::PDSPKaonAnalysis::FitLine(const std::vector<TVector3> & input) 
{
    TGraph2D * gr = new TGraph2D();

    for (size_t i = 0; i < input.size(); ++i) 
    {
        gr->SetPoint(i, input[i].X(), input[i].Y(), input[i].Z());
    }

    TVirtualFitter * min = TVirtualFitter::Fitter(0,4);
    min->SetObjectFit(gr);
    min->SetFCN(SumDistance2);

    double arglist[10];
    
    arglist[0] = -1;
    min->ExecuteCommand("SET PRINT",arglist,1);
    

    double pStart[4] = {1,1,1,1};
    min->SetParameter(0,"x0",pStart[0],0.01,0,0);
    min->SetParameter(1,"Ax",pStart[1],0.01,0,0);
    min->SetParameter(2,"y0",pStart[2],0.01,0,0);
    min->SetParameter(3,"Ay",pStart[3],0.01,0,0);

    arglist[0] = 1000; // number of function calls 
    arglist[1] = 0.001; // tolerance 
    min->ExecuteCommand("MIGRAD", arglist, 2);

      // get fit parameters
    double parFit[4];

    for (int i = 0; i < 4; ++i) 
    {
        parFit[i] = min->GetParameter(i);
    }

    double startX1, startY1, startZ1;
    double startX2, startY2, startZ2;
    line(0, parFit, startX1, startY1, startZ1);
    line(1, parFit, startX2, startY2, startZ2);
      
    TVector3 diff(startX2 - startX1, startY2 - startY1, startZ2 - startZ1);

    delete gr;
    delete min;
    return diff;
}*/

DEFINE_ART_MODULE(analysis::PDSPKaonAnalysis)
