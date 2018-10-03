
#include <sstream>

#include "RootIO.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

static G4String fileName = "FaserMC_Default.root";
static RootIO* instance = nullptr;

RootIO::RootIO()
{
  // initialize ROOT
  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->OpenFile(fileName);
  fAnalysisManager->CreateNtuple("hits", "Faser tracker digits"); 

  // tracker readout
  fAnalysisManager->CreateNtupleIColumn("digi_plane", fPlaneVector);
  fAnalysisManager->CreateNtupleIColumn("digi_module", fModuleVector);
  fAnalysisManager->CreateNtupleIColumn("digi_sensor", fSensorVector);
  fAnalysisManager->CreateNtupleIColumn("digi_row", fRowVector);
  fAnalysisManager->CreateNtupleIColumn("digi_strip", fStripVector);
  // not included in trackers, but useful for post hoc thresholding
  fAnalysisManager->CreateNtupleDColumn("digi_charge", fChargeVector);

  // truth values
  fAnalysisManager->CreateNtupleIColumn("truth_plane", fTruthPlaneVector);
  fAnalysisManager->CreateNtupleIColumn("truth_module", fTruthModuleVector);
  fAnalysisManager->CreateNtupleIColumn("truth_sensor", fTruthSensorVector);
  fAnalysisManager->CreateNtupleIColumn("truth_row", fTruthRowVector);
  fAnalysisManager->CreateNtupleIColumn("truth_strip", fTruthStripVector);

  fAnalysisManager->CreateNtupleIColumn("truth_track", fTruthTrackVector);
  fAnalysisManager->CreateNtupleDColumn("truth_energy", fTruthEnergyVector);

  fAnalysisManager->CreateNtupleDColumn("truth_global_x", fTruthGlobalXVector);
  fAnalysisManager->CreateNtupleDColumn("truth_global_y", fTruthGlobalYVector);
  fAnalysisManager->CreateNtupleDColumn("truth_global_z", fTruthGlobalZVector);

  fAnalysisManager->CreateNtupleDColumn("truth_local_x", fTruthLocalXVector);
  fAnalysisManager->CreateNtupleDColumn("truth_local_y", fTruthLocalYVector);
  fAnalysisManager->CreateNtupleDColumn("truth_local_z", fTruthLocalZVector);

  fAnalysisManager->CreateNtupleIColumn("truth_origin_track", fTruthOriginTrackVector);

  fAnalysisManager->FinishNtuple();
}

RootIO::~RootIO()
{
  delete fAnalysisManager;
}

RootIO* RootIO::GetInstance()
{
  if (instance == nullptr)
  {
    instance = new RootIO();
  }
  return instance;
}

void RootIO::SetFileName(G4String name)
{
  if (name == fileName) return;
  fileName = name;

  if (instance == nullptr) return;

  G4cout << "Closing previous ROOT file" << G4endl;

  instance->Close();
  delete instance;
  instance = nullptr;

  return;
}

void RootIO::AddDigits(FaserDigiCollection* dc)
{

  G4int nDigi = dc->entries();

  for (G4int i=0; i<nDigi; i++)
  {
    FaserDigi* digi = (*dc)[i];
    
    fPlaneVector.push_back(digi->Plane());
    fModuleVector.push_back(digi->Module());
    fSensorVector.push_back(digi->Sensor());
    fRowVector.push_back(digi->Row());
    fStripVector.push_back(digi->Strip());
    fChargeVector.push_back(digi->Charge()/coulomb*1e15);
  }
}

void RootIO::AddTruth(FaserSensorHitsCollection* hc)
{
  G4int nHits = hc->entries();

  for (G4int i=0; i<nHits; i++)
  {
    FaserSensorHit* hit = (*hc)[i];

    fTruthPlaneVector.push_back(hit->Plane());
    fTruthModuleVector.push_back(hit->Module());
    fTruthSensorVector.push_back(hit->Sensor());
    fTruthRowVector.push_back(hit->Row());
    fTruthStripVector.push_back(hit->Strip());

    fTruthTrackVector.push_back(hit->Track());
    fTruthEnergyVector.push_back(hit->Energy());
    
    fTruthGlobalXVector.push_back(hit->GlobalPos().x()/cm);
    fTruthGlobalYVector.push_back(hit->GlobalPos().y()/cm);
    fTruthGlobalZVector.push_back(hit->GlobalPos().z()/cm);
    
    fTruthLocalXVector.push_back(hit->LocalPos().x()/mm);
    fTruthLocalYVector.push_back(hit->LocalPos().y()/mm);
    fTruthLocalZVector.push_back(hit->LocalPos().z()/mm);

    fTruthOriginTrackVector.push_back(hit->OriginTrack());
  }
}

void RootIO::WriteEvent()
{
  fAnalysisManager->AddNtupleRow();

  fPlaneVector.clear();
  fModuleVector.clear();
  fSensorVector.clear();
  fRowVector.clear();
  fStripVector.clear();
  fChargeVector.clear();

  fTruthPlaneVector.clear();
  fTruthModuleVector.clear();
  fTruthSensorVector.clear();
  fTruthRowVector.clear();
  fTruthStripVector.clear();

  fTruthTrackVector.clear();
  fTruthEnergyVector.clear();

  fTruthGlobalXVector.clear();
  fTruthGlobalYVector.clear();
  fTruthGlobalZVector.clear();

  fTruthLocalXVector.clear();
  fTruthLocalYVector.clear();
  fTruthLocalZVector.clear();

  fTruthOriginTrackVector.clear();
}

void RootIO::Close()
{
  fAnalysisManager->Write();
  fAnalysisManager->CloseFile();  
}
