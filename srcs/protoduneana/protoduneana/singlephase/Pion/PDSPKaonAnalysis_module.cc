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

#include "larevt/SpaceChargeServices/SpaceChargeService.h" //for SCE correction


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
  TTree *fTree;
  unsigned int fEventID;
//  unsigned int fNPFParticles;
//  unsigned int fNPrimaries;
//  int fNPrimaryDaughters;
  std::string fPFParticleLabel;
  std::string fTruthLabel;
  std::string fTrackLabel;
  std::string fGeneratorTag;
  std::string fBeamModuleLabel;
  std::string fHitTag;

  float fPrimaryRecoPfpLength;
//  float fTrueKaonEnergy;
//  float fTruePrimaryEnergy;
  float fTrueBeamParticleEnergy;

  int fRecoBeamParticlePDGCode;
  int fIsBeamParticleReconstructed;
  int fTrueBeamParticlePDGCode;
  int fbestMatchedMCParticleFromPFParticlePdgCode;
  float fpurity;
  float fcompleteness;

  float fTrueBeamParticleStartX;
  float fTrueBeamParticleStartY;
  float fTrueBeamParticleStartZ;

//  float fTrueBeamParticleStartX_SCE_corrected;
//  float fTrueBeamParticleStartY_SCE_corrected;
//  float fTrueBeamParticleStartZ_SCE_corrected;

  float fBeamInstX_SCE_corrected;
  float fBeamInstY_SCE_corrected;
  float fBeamInstZ_SCE_corrected;


  float fTrueBeamParticleEndX;
  float fTrueBeamParticleEndY;
  float fTrueBeamParticleEndZ;

  float fTrueBeamParticleStartPx;
  float fTrueBeamParticleStartPy;
  float fTrueBeamParticleStartPz;

  int  fTrueBeamParticleNhits;

  float fRecoBeamParticleStartX;
  float fRecoBeamParticleStartY;
  float fRecoBeamParticleStartZ;

  float fRecoBeamParticleStartPx;
  float fRecoBeamParticleStartPy;
  float fRecoBeamParticleStartPz;

  int fRecoBeamParticleNhits;
  float beam_inst_X;
  float beam_inst_Y;
  float beam_inst_Z;

  int fNumberOfReconstructedBeamParticle;
  bool fMCHasBI;
  bool beam_inst_valid;

  protoana::ProtoDUNETruthUtils fProtoDUNETruthUtils;
  protoana::ProtoDUNEPFParticleUtils fProtoDUNEPFParticleUtils;
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
    fTrueBeamParticleEnergy = -9999.;
    fTrueBeamParticlePDGCode = 0;
    fpurity = -1.f;
    fcompleteness = -1.f;

    fTrueBeamParticleStartX                     = -9999;
    fTrueBeamParticleStartY                     = -9999;
    fTrueBeamParticleStartZ                     = -9999;

//    fTrueBeamParticleStartX_SCE_corrected       = -9999;
//    fTrueBeamParticleStartY_SCE_corrected       = -9999;
//    fTrueBeamParticleStartZ_SCE_corrected       = -9999;

    fBeamInstX_SCE_corrected                    = -9999;
    fBeamInstY_SCE_corrected                    = -9999;
    fBeamInstZ_SCE_corrected                    = -9999;

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

    fNumberOfReconstructedBeamParticle = -1;

//==================Accessing Truth Info Block=============================

/*========To do list=======================================================

    rename beam_inst_vertex into start_vertex_XYZ
=========================================================================*/
/*    if (!e.isRealData())
    {
        art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);

//        std::cout << "mcParticles: " << mcParticles->size() << std::endl;

        if(mcParticles.isValid())
        {
            for(unsigned int t = 0; t < mcParticles->size(); ++t)
            {
                const simb::MCParticle trueParticle = mcParticles->at(t);
                if(trueParticle.Process() == "primary" && fProtoDUNEPFParticleUtils.IsBeamParticle(particle, e, fPFParticleLabel) == 1)
                {
                    fTruePrimaryEnergy = trueParticle.E();
                }
            }
        }
    }*/

//=========================================================================
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    const simb::MCParticle* true_beam_particle = 0x0;

    std::vector<art::Ptr<beam::ProtoDUNEBeamEvent>> beamVec;

    if(!e.isRealData())
    {
        auto mcTruths = e.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorTag);
        
        true_beam_particle = fProtoDUNETruthUtils.GetGeantGoodParticle((*mcTruths)[0],e);

        fTrueBeamParticleEnergy = true_beam_particle->E();
        fTrueBeamParticlePDGCode = true_beam_particle->PdgCode();

        fTrueBeamParticleStartX     = true_beam_particle->Position(0).X();
        fTrueBeamParticleStartY     = true_beam_particle->Position(0).Y();
        fTrueBeamParticleStartZ     = true_beam_particle->Position(0).Z();

        fTrueBeamParticleEndX     = true_beam_particle->EndX();
        fTrueBeamParticleEndY     = true_beam_particle->EndY();
        fTrueBeamParticleEndZ     = true_beam_particle->EndZ();

//        std::cout << "fTrueBeamParticleStartZ: " << fTrueBeamParticleStartZ << std::endl;
//        std::cout << "fTrueBeamParticleEndZ: " << fTrueBeamParticleEndZ << std::endl;

        fTrueBeamParticleStartPx    = true_beam_particle->Px(); // incorrect possibly
        fTrueBeamParticleStartPy    = true_beam_particle->Py(); // incorrect possibly
        fTrueBeamParticleStartPz    = true_beam_particle->Pz(); // incorrect possibly
//        std::cout << "1" << std::endl;
/*        auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
        fTrueBeamParticleStartX_SCE_corrected = fTrueBeamParticleStartX-SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleStartX,fTrueBeamParticleStartY,fTrueBeamParticleStartZ)).X();
        fTrueBeamParticleStartY_SCE_corrected = fTrueBeamParticleStartY+SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleStartX,fTrueBeamParticleStartY,fTrueBeamParticleStartZ)).Y();
        fTrueBeamParticleStartZ_SCE_corrected = fTrueBeamParticleStartZ+SCE->GetPosOffsets(geo::Point_t(fTrueBeamParticleStartX,fTrueBeamParticleStartY,fTrueBeamParticleStartZ)).Z();*/
//        std::cout << "2" << std::endl;

        auto beamHandle = e.getValidHandle<std::vector<beam::ProtoDUNEBeamEvent>>("generator");

        if (beamHandle.isValid())
        {
            art::fill_ptr_vector(beamVec, beamHandle);
        }

        else
            std::cout << "invalid beam handle" << std::endl;

//        std::cout << "3" << std::endl;
        const beam::ProtoDUNEBeamEvent beamEvent = *(beamVec.at(0));
//        const beam::ProtoDUNEBeamEvent beamEvent = fBeamlineUtils.GetBeamEvent(e);
//        std::cout << "3.5" << std::endl;
//        std::cout << "fBeamlineUtils.IsGoodBeamlineTrigger(e): " << fBeamlineUtils.IsGoodBeamlineTrigger(e) << std::endl;
/*        if (e.isRealData() && !fBeamlineUtils.IsGoodBeamlineTrigger(e)) 
        {
            beam_inst_valid = false;
            return;
        }*/
//        std::cout << "4" << std::endl;
        int nTracks = beamEvent.GetBeamTracks().size();
//        std::cout << "nTracks: " << nTracks << std::endl;
        if( nTracks > 0 )
        {
            beam_inst_X = beamEvent.GetBeamTracks()[0].Trajectory().End().X();
            beam_inst_Y = beamEvent.GetBeamTracks()[0].Trajectory().End().Y();
            beam_inst_Z = beamEvent.GetBeamTracks()[0].Trajectory().End().Z();

//            std::cout << "true mc length: " << beamEvent.GetBeamTracks()[0].Length() << std::endl;

//            std::cout << "beam_inst_X: " << beam_inst_X << std::endl;
//            std::cout << "beam_inst_Y: " << beam_inst_Y << std::endl;
//            std::cout << "beam_inst_Z: " << beam_inst_Z << std::endl;

            auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
            fBeamInstX_SCE_corrected = beam_inst_X-SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).X();
            fBeamInstY_SCE_corrected = beam_inst_Y+SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).Y();
            fBeamInstZ_SCE_corrected = beam_inst_Z+SCE->GetPosOffsets(geo::Point_t(beam_inst_X,beam_inst_Y,beam_inst_Z)).Z();

//            std::cout << "fBeamInstX_SCE_corrected: " << fBeamInstX_SCE_corrected << std::endl;
//            std::cout << "fBeamInstY_SCE_corrected: " << fBeamInstY_SCE_corrected << std::endl;
//            std::cout << "fBeamInstZ_SCE_corrected: " << fBeamInstZ_SCE_corrected << std::endl;
        }
//        std::cout << "5" << std::endl;

//        std::cout << "fTrueBeamParticleEnergy: " << fTrueBeamParticleEnergy << std::endl;
//        std::cout << "fTrueBeamParticlePDGCode: " << fTrueBeamParticlePDGCode << std::endl;

//        std::cout << "fTrueBeamParticleStartX: " << fTrueBeamParticleStartX << std::endl;
//        std::cout << "fTrueBeamParticleStartY: " << fTrueBeamParticleStartY << std::endl;
//        std::cout << "fTrueBeamParticleStartZ: " << fTrueBeamParticleStartZ << std::endl;

//        std::cout << "fTrueBeamParticleStartPx: " << fTrueBeamParticleStartPx << std::endl;
//        std::cout << "fTrueBeamParticleStartPy: " << fTrueBeamParticleStartPy << std::endl;
//        std::cout << "fTrueBeamParticleStartPz: " << fTrueBeamParticleStartPz << std::endl;

        fTrueBeamParticleNhits = fProtoDUNETruthUtils.GetMCParticleHits( clockData, *true_beam_particle, e, fHitTag ).size();
    }

//==================Accessing Reco Info Block==============================

/*========To do list=======================================================

    search for kaons


=========================================================================*/

//    fNPrimaries = 0;

    art::ValidHandle<std::vector<recob::PFParticle>> recoParticles = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
//    std::cout << "recoParticles: " << recoParticles->size() << std::endl;

    std::vector<recob::PFParticle> recoBeamParticles;
    std::vector<recob::PFParticle> recoCosmicMuons;
    std::vector<recob::PFParticle> recoBeamParticles_2;

/*    std::cout << "recoBeamParticle before: " << recoBeamParticles.size() << std::endl;
    std::cout << "recoBeamParticles_2 before: " << recoBeamParticles_2.size() << std::endl;
    std::cout << "recoCosmicMuons before: " << recoCosmicMuons.size() << std::endl;*/

    for(unsigned int k = 0; k < recoParticles->size(); ++k)
    {
        const recob::PFParticle particle = recoParticles->at(k);

        if(fProtoDUNEPFParticleUtils.IsBeamParticle(particle, e, fPFParticleLabel) == 1)
        {
            //std::cout << "Is particle primary? " << particle.IsPrimary() << std::endl;
//            fRecoBeamParticlePDGCode = particle.Self();
//            std::cout << "recoBeamParticles_2 PDGCode written into tree: " << fRecoBeamParticlePDGCode << std::endl;

            recoBeamParticles_2.push_back(particle);
        }

        if(abs(particle.PdgCode()) != 13)         
            recoBeamParticles.push_back(particle);

        else
            recoCosmicMuons.push_back(particle);
    }

/*    std::cout << "recoBeamParticles after: " << recoBeamParticles.size() << std::endl;*/
//    std::cout << "recoBeamParticles_2 after: " << recoBeamParticles_2.size() << std::endl;
    fNumberOfReconstructedBeamParticle = recoBeamParticles_2.size();

    fIsBeamParticleReconstructed = 0;    

    if (recoBeamParticles_2.size() == 1)
    {
        fIsBeamParticleReconstructed = 1;
//        fRecoBeamParticlePDGCode = recoBeamParticles_2[0]->Self(); 
        fRecoBeamParticlePDGCode = recoBeamParticles_2[0].PdgCode();  

        const simb::MCParticle* bestMatchedMcParticleFromPFParticle;

        bestMatchedMcParticleFromPFParticle = (fProtoDUNETruthUtils.GetMCParticleFromPFParticle(clockData, recoBeamParticles_2[0], e, fPFParticleLabel));
        fbestMatchedMCParticleFromPFParticlePdgCode = bestMatchedMcParticleFromPFParticle->PdgCode();

        fpurity = fProtoDUNETruthUtils.GetPurity(clockData, recoBeamParticles_2[0], e, fPFParticleLabel); 
        std::cout << "purity: " << fpurity << std::endl; 

        fcompleteness = fProtoDUNETruthUtils.GetCompleteness(clockData, recoBeamParticles_2[0], e, "pandora", "hitpdune");
        std::cout << "completeness: " << fcompleteness << std::endl; 

//        std::cout << "bestMatchedMCParticleFromPFParticlePdgCode: " << fbestMatchedMCParticleFromPFParticlePdgCode << std::endl;

//        std::cout << "recoBeamParticles_2 PDGCode: " << recoBeamParticles_2[0] << std::endl;
//        std::cout << "recoBeamParticles_2 PDGCode written into tree: " << fRecoBeamParticlePDGCode << std::endl;
        
        auto recoBeamParticleVertex = fProtoDUNEPFParticleUtils.GetPFParticleVertex(recoBeamParticles_2[0], e, fPFParticleLabel, fTrackLabel);
    
        fRecoBeamParticleStartX = recoBeamParticleVertex.X();
        fRecoBeamParticleStartY = recoBeamParticleVertex.Y();
        fRecoBeamParticleStartZ = recoBeamParticleVertex.Z();


//        std::cout << "fRecoBeamParticleStartX: " << fRecoBeamParticleStartX << std::endl;
//        std::cout << "fRecoBeamParticleStartY: " << fRecoBeamParticleStartY << std::endl;
//        std::cout << "fRecoBeamParticleStartZ: " << fRecoBeamParticleStartZ << std::endl;

//        fRecoBeamParticleStartPx = recoBeamParticles_2[0]->Px();
//        fRecoBeamParticleStartPy = recoBeamParticles_2[0]->Py();
//        fRecoBeamParticleStartPz = recoBeamParticles_2[0]->Pz();

        const std::vector< art::Ptr< recob::Hit > > beamPFP_hits = fProtoDUNEPFParticleUtils.GetPFParticleHits_Ptrs( recoBeamParticles_2[0], e, fPFParticleLabel );

        fRecoBeamParticleNhits = beamPFP_hits.size();
        std::cout << "fRecoBeamParticleNhits: " << fRecoBeamParticleNhits << std::endl;
        
    }

    else
    {
        std::cout << "The number of reconstructed beam particle is not 1." << std::endl;
        std::cout << "Number of reconstructed beam particle is " << recoBeamParticles_2.size() << std::endl;
    }

//    std::cout << "fIsBeamParticleReconstructed: " << fIsBeamParticleReconstructed << std::endl;

//    std::cout << "recoCosmicMuons after: " << recoCosmicMuons.size() << std::endl;*/


    if(recoParticles.isValid())
    {
        const art::FindManyP<recob::Track> findTracks(recoParticles, e, fTrackLabel);
//        std::cout << "findTracks: " << findTracks.size() << std::endl;

        fPrimaryRecoPfpLength = -1.0;//not really kaon length right now. More like primary particle length

        for(unsigned int p = 0; p < recoParticles->size(); ++p)
        {
            const recob::PFParticle particle = recoParticles->at(p);

            if(particle.IsPrimary() && fProtoDUNEPFParticleUtils.IsBeamParticle(particle, e, fPFParticleLabel) == 1)
            {
                const std::vector<art::Ptr<recob::Track>>pfpTracks = findTracks.at(p);

                if(pfpTracks.size() == 1 && particle.PdgCode() == 211)
                {
                    art::Ptr<recob::Track> thisTrack = pfpTracks.at(0);
                    fPrimaryRecoPfpLength = thisTrack->Length();
//                    std::cout << "fPrimaryRecoPfpLength: " << fPrimaryRecoPfpLength <<  std::endl;
                }
            }
        }
    }

    fTree->Fill();

}

void analysis::PDSPKaonAnalysis::beginJob()
{
  // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("tree","Analyser Output Tree");

  // Truth Tree Branches    
    fTree->Branch("eventID",&fEventID,"eventID/i");
//    fTree->Branch("nPFParticles",&fNPFParticles,"nPFParticles/i");
//    fTree->Branch("nPrimaries",&fNPrimaries,"nPrimaries/i");
//    fTree->Branch("nPrimaryDaughters",&fNPrimaryDaughters,"nPrimaryDaughters/i");
//    fTree->Branch("trueKaonEnergy",&fTrueKaonEnergy,"trueKaonEnergy/F");
//    fTree->Branch("truePrimaryEnergy",&fTruePrimaryEnergy,"truePrimaryEnergy/F");
    fTree->Branch("trueBeamParticleEnergy",&fTrueBeamParticleEnergy,"trueBeamParticleEnergy/F");
    fTree->Branch("trueBeamParticlePDGCode",&fTrueBeamParticlePDGCode,"trueBeamParticlePDGCode/I");

//    fTree->Branch("kaonLength",&fKaonLength,"kaonLength/F");

    fTree->Branch("trueBeamParticleStartX",&fTrueBeamParticleStartX,"trueBeamParticleStartX/F");
    fTree->Branch("trueBeamParticleStartY",&fTrueBeamParticleStartY,"trueBeamParticleStartY/F");
    fTree->Branch("trueBeamParticleStartZ",&fTrueBeamParticleStartZ,"trueBeamParticleStartZ/F");

    fTree->Branch("trueBeamParticleEndX",&fTrueBeamParticleEndX,"trueBeamParticleEndX/F");
    fTree->Branch("trueBeamParticleEndY",&fTrueBeamParticleEndY,"trueBeamParticleEndY/F");
    fTree->Branch("trueBeamParticleEndZ",&fTrueBeamParticleEndZ,"trueBeamParticleEndZ/F");

//    fTree->Branch("trueBeamParticleStartX_SCE_corrected",&fTrueBeamParticleStartX_SCE_corrected,"trueBeamParticleStartX_SCE_corrected/F");
//    fTree->Branch("trueBeamParticleStartY_SCE_corrected",&fTrueBeamParticleStartY_SCE_corrected,"trueBeamParticleStartY_SCE_corrected/F");
//    fTree->Branch("trueBeamParticleStartZ_SCE_corrected",&fTrueBeamParticleStartZ_SCE_corrected,"trueBeamParticleStartZ_SCE_corrected/F");

    fTree->Branch("beam_inst_X",&beam_inst_X,"beam_inst_X/F");
    fTree->Branch("beam_inst_Y",&beam_inst_Y,"beam_inst_Y/F");
    fTree->Branch("beam_inst_Z",&beam_inst_Z,"beam_inst_Z/F");

    fTree->Branch("beamInstX_SCE_corrected",&fBeamInstX_SCE_corrected,"beamInstX_SCE_corrected/F");
    fTree->Branch("beamInstY_SCE_corrected",&fBeamInstY_SCE_corrected,"beamInstY_SCE_corrected/F");
    fTree->Branch("beamInstZ_SCE_corrected",&fBeamInstZ_SCE_corrected,"beamInstZ_SCE_corrected/F");

    fTree->Branch("trueBeamParticleStartPx",&fTrueBeamParticleStartPx,"trueBeamParticleStartPx/F");
    fTree->Branch("trueBeamParticleStartPy",&fTrueBeamParticleStartPy,"trueBeamParticleStartPy/F");
    fTree->Branch("trueBeamParticleStartPz",&fTrueBeamParticleStartPz,"trueBeamParticleStartPz/F");

    fTree->Branch("trueBeamParticleNhits", &fTrueBeamParticleNhits, "trueBeamParticleNhits/I");
    fTree->Branch("recoBeamParticleNhits", &fRecoBeamParticleNhits, "recoBeamParticleNhits/I");
//    fTree->Branch("beam_inst_X",&beam_inst_X,"beam_inst_X/F");
//    fTree->Branch("beam_inst_Y",&beam_inst_Y,"beam_inst_Y/F");
//    fTree->Branch("beam_inst_Z",&beam_inst_Z,"beam_inst_Z/F");
//    fTree->Branch("beam_inst_valid", &beam_inst_valid,"beam_inst_X/I");

  // Reco Tree Branches   
    fTree->Branch("bestMatchedMCParticleFromPFParticlePdgCode",&fbestMatchedMCParticleFromPFParticlePdgCode,"bestMatchedMCParticleFromPFParticlePdgCode/I");
    fTree->Branch("purity",&fpurity,"purity/F");
    fTree->Branch("completeness",&fcompleteness,"completeness/F");
    fTree->Branch("isBeamParticleReconstructed",&fIsBeamParticleReconstructed,"isBeamParticleReconstructed/i");
    fTree->Branch("primaryRecoPfpLength",&fPrimaryRecoPfpLength,"primaryRecoPfpLength/F"); 

    fTree->Branch("recoBeamParticleStartX",&fRecoBeamParticleStartX,"recoBeamParticleStartX/F");
    fTree->Branch("recoBeamParticleStartY",&fRecoBeamParticleStartY,"recoBeamParticleStartY/F");
    fTree->Branch("recoBeamParticleStartZ",&fRecoBeamParticleStartZ,"recoBeamParticleStartZ/F");

    fTree->Branch("recoBeamParticleStartPx",&fRecoBeamParticleStartPx,"recoBeamParticleStartPx/F");
    fTree->Branch("recoBeamParticleStartPy",&fRecoBeamParticleStartPy,"recoBeamParticleStartPy/F");
    fTree->Branch("recoBeamParticleStartPz",&fRecoBeamParticleStartPz,"recoBeamParticleStartPz/F");


    fTree->Branch("numberOfReconstructedBeamParticle",&fNumberOfReconstructedBeamParticle,"numberOfReconstructedBeamParticle/I");
}

void analysis::PDSPKaonAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(analysis::PDSPKaonAnalysis)
