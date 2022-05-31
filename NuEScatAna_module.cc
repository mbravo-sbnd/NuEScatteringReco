/***************************************************************************
 * Analyzer module created for Nu-e elastic scattering analysis. 
 * It makes use of objects from generation, geant4 propagation and reco1 
 * simulation steps.
 * Output: root file called analysisOutput.root, with a TDirectory called 
 * ana. ana contains the TTree analysisOutputTree, with the following 
 * branches:
 * - EventID (int): ID of the event
   - EventSubRun (int): Sub-Run of the event
   - EventRun (int): Run of the event
   - POT (double): number of protons on target to generate the sample
   - TrueIntMode (int): MC interaction mode of the event
   - TrueIntType (int): MC interaction type of the event
   - TrueIntVerX (double): X position of the MC primary interaction vertex 
      (in cm, SBND coordinate system)
   - TrueIntVerY (double): Y position of the MC primary interaction vertex 
      (in cm, SBND coordinate system)
   - TrueIntVerZ (double): Z position of the MC primary interaction vertex 
      (in cm, SBND coordinate system)
    - trueParts (std::vector<int>): MC PDG number of all simulated particles
    - truePartsX (std::vector<double>):MC X position of the origin of all 
      simulated particles (in cm, SBND coordinate system)
    - truePartsY (std::vector<double>):MC Y position of the origin of all 
      simulated particles (in cm, SBND coordinate system)
    - truePartsZ (std::vector<double>):MC Z position of the origin of all 
      simulated particles (in cm, SBND coordinate system)
    - truePartsPx (std::vector<double> ): MC X component of the original 
      momentum of all simulated particles (GeV/c) 
    - truePartsPy (std::vector<double>): MC Y component of the original 
      momentum of all simulated particles (GeV/c) 
    - truePartsPz (std::vector<double>): MC Z component of the original 
      momentum of all simulated particles (GeV/c) 
    - truePartsE (std::vector<double>): MC energy of the original 
      momentum of all simulated particles (GeV) 
    - truePartsStatus (std::vector<int>): MC status of all the simulated 
      particles (0=original, 1=from an interaction)
    - G4EDepos (std::vector<double>): simulated energy depositions (MeV)
    - G4EDeposPDG (std::vector<int>): PDG codes of the particles that made the 
      energy depositions
    - G4EDeposT (std::vector<double>): times of the energy depositions (ns)
    - G4EDeposMidX (std::vector<double>): X middle points of energy depositions
      (in cm, SBND coordinate system)
    - G4EDeposMidY (std::vector<double>): Y middle points of energy depositions
      (in cm, SBND coordinate system)
    - G4EDeposMidZ (std::vector<double>): Z middle points of energy depositions
      (in cm, SBND coordinate system)
    - G4EDeposStartX (std::vector<double>): X start points of energy depositions
      (in cm, SBND coordinate system)
    - G4EDeposStartY ( std::vector<double>): Y start points of energy depositions
      (in cm, SBND coordinate system)
    - G4EDeposStartZ (std::vector<double>): Z start points of energy depositions
      (in cm, SBND coordinate system)
    - G4EDeposEndX (std::vector<double>): X final points of energy depositions
      (in cm, SBND coordinate system)
    - G4EDeposEndY (std::vector<double>): Y final points of energy depositions
      (in cm, SBND coordinate system)
    - G4EDeposEndZ (std::vector<double>): Z final points of energy depositions
      (in cm, SBND coordinate system)
    - nTotalHits (int): total number of hits in the event
    - spX (std::vector<double>): X positions of the reconstructed space-points
      (in cm, SBND coordinate system)
    - spY (std::vector<double>): Y positions of the reconstructed space-points
      (in cm, SBND coordinate system)
    - spZ (std::vector<double>): Z positions of the reconstructed space-points
      (in cm, SBND coordinate system)
    - spHitsPTime (std::vector<double>): peak time of the hit associated to 
      each space-point (in ns)
    - spHitsInteg (std::vector<double>): integral of the hit associated to
      each space-point (in ADC x tick, 1 tick = 2 ns)
    - spHitsPlane (std::vector<int>): plane of the hit associated to each 
      space-point (0, 1 = induction planes, 2 = collection plane)
    - spHitsWire (std::vector<int>): wire ID of the hit associated to each 
      space point
    - spHitsTPC (std::vector<int>): TPC ID of the hit associated to each 
      space point
 * 
 * *************************************************************************/

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

#include "larcore/Geometry/Geometry.h"

#include "lardata/RecoBaseProxy/ProxyBase/withCollectionProxy.h"
#include "lardata/RecoBaseProxy/ProxyBase/withAssociated.h"
#include "lardata/RecoBaseProxy/ProxyBase/withParallelData.h"
#include "lardata/RecoBaseProxy/ProxyBase/withZeroOrOne.h"
#include "lardata/RecoBaseProxy/ProxyBase/getCollection.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"

#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/RecoBase/MCSFitResult.h"


#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larreco/SpacePointSolver/TripletFinder.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventory.h"
#include "larsim/MCCheater/ParticleInventoryService.h"


#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "larevt/SpaceCharge/SpaceCharge.h"

#include "sbndcode/RecoUtils/RecoUtils.h"


#include "TTree.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"

#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>


namespace test {
  class NuEScatAna;
}


class test::NuEScatAna : public art::EDAnalyzer {
public:
  explicit NuEScatAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuEScatAna(NuEScatAna const&) = delete;
  NuEScatAna(NuEScatAna&&) = delete;
  NuEScatAna& operator=(NuEScatAna const&) = delete;
  NuEScatAna& operator=(NuEScatAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginSubRun(const art::SubRun& eSR);

  void reconfigure(fhicl::ParameterSet const& p);

  // Selected optional functions.
  void beginJob() override;
//  void beginSubRun (const art::SubRun& eSR);
  void endJob() override;

private:

  // Declare member data here.
  std::string name_HitLabel;

  TTree* fTree;
  
  int fEventID;
  int EventSubRun;
  int EventRun;

  double POT;

  int TrueIntMode;
  int TrueIntType;
  double TrueIntVerX;
  double TrueIntVerY;
  double TrueIntVerZ;

  std::vector<int> trueParts;
  std::vector<double> truePartsX;
  std::vector<double> truePartsY;
  std::vector<double> truePartsZ;
  std::vector<double> truePartsPx;
  std::vector<double> truePartsPy;
  std::vector<double> truePartsPz;
  std::vector<double> truePartsE;
  std::vector<int> truePartsStatus;

  std::vector<int> G4parts;
  std::vector<int> G4partsStatusCode;
  std::vector<std::string> G4partsGenProcess;

  std::vector<double> G4EDepos;
  std::vector<int> G4EDeposPDG;
  std::vector<double> G4EDeposT;
  std::vector<double> G4EDeposMidX;
  std::vector<double> G4EDeposMidY;
  std::vector<double> G4EDeposMidZ;
  std::vector<double> G4EDeposStartX;
  std::vector<double> G4EDeposStartY;
  std::vector<double> G4EDeposStartZ;
  std::vector<double> G4EDeposEndX;
  std::vector<double> G4EDeposEndY;
  std::vector<double> G4EDeposEndZ;

  int nTotalHits;

  std::vector<double> spX;
  std::vector<double> spY;
  std::vector<double> spZ;
  std::vector<double> spHitsPTime;
  std::vector<double> spHitsInteg;
  std::vector<int> spHitsPlane;
  std::vector<int> spHitsWire;
  std::vector<int> spHitsTPC;

};


test::NuEScatAna::NuEScatAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  this->reconfigure(p);
}

void test::NuEScatAna::reconfigure(fhicl::ParameterSet const& p)
{
//La función consta en coger el nombre indicado como "HitLabel" en el fhicl para que sea ese método el que me analice los hits.
 name_HitLabel = p.get<std::string>("HitLabelmode");
}

void test::NuEScatAna::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  fTree=tfs->make<TTree>("analysisOutputTree", "Analysis Output Tree");

  fTree->Branch("EventID", &fEventID, "EventID/I");
  fTree->Branch("EventRun", &EventRun, "EventRun/I");
  fTree->Branch("EventSubRun", &EventSubRun, "EventSubRun/I");

  fTree->Branch("POT", &POT, "POT/D");

  fTree->Branch("TrueIntMode", &TrueIntMode, "TrueIntMode/I");
  fTree->Branch("TrueIntType", &TrueIntType, "TrueIntType/I");

  fTree->Branch("TrueIntVerX", &TrueIntVerX, "TrueIntVerX/D");
  fTree->Branch("TrueIntVerY", &TrueIntVerY, "TrueIntVerY/D");
  fTree->Branch("TrueIntVerZ", &TrueIntVerZ, "TrueIntVerZ/D");

  fTree->Branch("trueParts", &trueParts);
  fTree->Branch("truePartsX", &truePartsX);
  fTree->Branch("truePartsY", &truePartsY);
  fTree->Branch("truePartsZ", &truePartsZ);
  fTree->Branch("truePartsPx", &truePartsPx);
  fTree->Branch("truePartsPy", &truePartsPy);
  fTree->Branch("truePartsPz", &truePartsPz);
  fTree->Branch("truePartsE", &truePartsE);
  fTree->Branch("truePartsStatus", &truePartsStatus);

  fTree->Branch("G4parts", &G4parts);
  fTree->Branch("G4partsStatusCode", &G4partsStatusCode);
  fTree->Branch("G4partsGenProcess", &G4partsGenProcess);

  fTree->Branch("G4EDepos", &G4EDepos);
  fTree->Branch("G4EDeposT", &G4EDeposT);
  fTree->Branch("G4EDeposPDG", &G4EDeposPDG);
  fTree->Branch("G4EDeposMidX", &G4EDeposMidX);
  fTree->Branch("G4EDeposMidY", &G4EDeposMidY);
  fTree->Branch("G4EDeposMidZ", &G4EDeposMidZ);
  fTree->Branch("G4EDeposStartX", &G4EDeposStartX);
  fTree->Branch("G4EDeposStartY", &G4EDeposStartY);
  fTree->Branch("G4EDeposStartZ", &G4EDeposStartZ);
  fTree->Branch("G4EDeposEndX", &G4EDeposEndX);
  fTree->Branch("G4EDeposEndY", &G4EDeposEndY);
  fTree->Branch("G4EDeposEndZ", &G4EDeposEndZ);

  fTree->Branch("nTotalHits", &nTotalHits, "nTotalHits/I");

  fTree->Branch("spX", &spX);
  fTree->Branch("spY", &spY);
  fTree->Branch("spZ", &spZ);
  fTree->Branch("spHitsPTime", &spHitsPTime);
  fTree->Branch("spHitsInteg", &spHitsInteg);
  fTree->Branch("spHitsPlane", &spHitsPlane);
  fTree->Branch("spHitsWire", &spHitsWire);
  fTree->Branch("spHitsTPC", &spHitsTPC);
  
}



void test::NuEScatAna::beginSubRun(const art::SubRun& eSR)
{
  //Information about the number of protons on target (POT) needed to produce the event sample

  art::Handle<sumdata::POTSummary> POTListHandle;

    eSR.getByLabel("generator",POTListHandle);
    POT=POTListHandle->totpot;
}


void test::NuEScatAna::analyze(art::Event const& e)
{
  trueParts.clear();
  truePartsX.clear();
  truePartsY.clear();
  truePartsZ.clear();
  truePartsPx.clear();
  truePartsPy.clear();
  truePartsPz.clear();
  truePartsE.clear();
  truePartsStatus.clear();

  G4EDepos.clear();
  G4EDeposT.clear();
  G4EDeposPDG.clear();
  G4EDeposMidX.clear();
  G4EDeposMidY.clear();
  G4EDeposMidZ.clear();
  G4EDeposStartX.clear();
  G4EDeposStartY.clear();
  G4EDeposStartZ.clear();
  G4EDeposEndX.clear();
  G4EDeposEndY.clear();
  G4EDeposEndZ.clear();

  spX.clear();
  spY.clear();
  spZ.clear();
  spHitsPTime.clear();
  spHitsInteg.clear();
  spHitsPlane.clear();
  spHitsWire.clear();
  spHitsTPC.clear();

//................... 

  //Event identification numbers: run, subrun and ID

  fEventID=e.id().event();
  EventSubRun=e.id().subRun();
  EventRun=e.id().run();

//...................

  //Monte-Carlo (true) information of the event. Here, information 
  //is obtained from GENIE step

  art::Handle<std::vector<simb::MCTruth>> mctruths;
  e.getByLabel("generator", mctruths);
         
  for (auto const& truth : *mctruths) {

    simb::MCNeutrino const& primInt = truth.GetNeutrino();
    simb::MCParticle const& lep = primInt.Lepton();
    
    TrueIntMode = primInt.Mode();
    TrueIntType = primInt.InteractionType();
       
    TrueIntVerX=lep.Vx();
    TrueIntVerY=lep.Vy();
    TrueIntVerZ=lep.Vz();

      for (int i=0; i<truth.NParticles(); i++){
        simb::MCParticle const& particle = truth.GetParticle(i);
        truePartsE.push_back(particle.E());
        trueParts.push_back(particle.PdgCode());
        truePartsPx.push_back(particle.Px());
        truePartsPy.push_back(particle.Py());
        truePartsPz.push_back(particle.Pz());
        truePartsX.push_back(particle.Vx());
        truePartsY.push_back(particle.Vy());
        truePartsZ.push_back(particle.Vz());
        truePartsStatus.push_back(particle.StatusCode());

      }
    

  }

//............................

  //Monte-Carlo (true) energy-depositions information of the event. Here, information 
  //is obtained from GEANT4 step

  art::Handle<std::vector<sim::SimEnergyDeposit>> SimEDeps;
  e.getByLabel("largeant", "LArG4DetectorServicevolTPCActive", SimEDeps);
  for (auto const& EDep : *SimEDeps){
    G4EDepos.push_back(EDep.Energy());
    G4EDeposT.push_back(EDep.T());
    G4EDeposPDG.push_back(EDep.PdgCode());
    G4EDeposMidX.push_back(EDep.MidPointX());
    G4EDeposMidY.push_back(EDep.MidPointY());
    G4EDeposMidZ.push_back(EDep.MidPointZ());
    G4EDeposStartX.push_back(EDep.StartX());
    G4EDeposStartY.push_back(EDep.StartY());
    G4EDeposStartZ.push_back(EDep.StartZ());
    G4EDeposEndX.push_back(EDep.EndX());
    G4EDeposEndY.push_back(EDep.EndY());
    G4EDeposEndZ.push_back(EDep.EndZ());
  }

//...............................

  //First-level reconstructed information of the event. Here, information 
  //is obtained from reco1 step

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitvector;
  e.getByLabel(name_HitLabel, hitListHandle);
  art::fill_ptr_vector(hitvector, hitListHandle);

  nTotalHits=hitvector.size();

  
//.............................

  //First-level reconstructed information of the event. Here, information 
  //is obtained from reco1 step

  art::Handle<std::vector<recob::SpacePoint>> eventSpacePoints;
  std::vector<art::Ptr<recob::SpacePoint>> eventSpacePointsVect;

  e.getByLabel("pandora", eventSpacePoints);
  art::fill_ptr_vector(eventSpacePointsVect, eventSpacePoints);

  art::FindManyP<recob::Hit> SPToHitAssoc (eventSpacePointsVect, e, "pandora");

  for (const art::Ptr<recob::SpacePoint> &SP: eventSpacePointsVect){

    spX.push_back(SP->position().X());
    spY.push_back(SP->position().Y());
    spZ.push_back(SP->position().Z());

    std::vector<art::Ptr<recob::Hit>> SPHit = SPToHitAssoc.at(SP.key());

    spHitsTPC.push_back(SPHit.at(0)->WireID().TPC);
    spHitsWire.push_back(SPHit.at(0)->WireID().Wire);
    spHitsPlane.push_back(SPHit.at(0)->WireID().Plane);
    spHitsPTime.push_back(SPHit.at(0)->PeakTime());
    spHitsInteg.push_back(SPHit.at(0)->Integral());


  }//end all space points for an event loop

//.............................


  fTree->Fill(); 
}



void test::NuEScatAna::endJob()
{


}

DEFINE_ART_MODULE(test::NuEScatAna)