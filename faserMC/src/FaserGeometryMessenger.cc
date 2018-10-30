#include "FaserGeometryMessenger.hh"
#include "FaserDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

FaserGeometryMessenger::FaserGeometryMessenger(FaserDetectorConstruction* detectorConstruction)
  : fDetectorConstruction(detectorConstruction)
{
  geometryDirectory = new G4UIdirectory("/faser/geo/");
  geometryDirectory->SetGuidance("Detector geometry initialization parameters.");

  cmd_sensor_readoutStrips = new G4UIcmdWithAnInteger("/faser/geo/readoutStrips", this);
  cmd_sensor_readoutStrips->SetGuidance("Number of active readout strips per sensor.");
  cmd_sensor_readoutStrips->SetParameterName("readoutStrips", true, true);
  cmd_sensor_readoutStrips->AvailableForStates(G4State_PreInit);

  cmd_sensor_stripPitch = new G4UIcmdWithADoubleAndUnit("/faser/geo/stripPitch", this);
  cmd_sensor_stripPitch->SetGuidance("Distance between strips in the silicon wafers.");
  cmd_sensor_stripPitch->SetParameterName("stripPitch", true, true);
  cmd_sensor_stripPitch->SetDefaultUnit("mm");
  cmd_sensor_stripPitch->SetUnitCandidates("mm micron cm m");
  cmd_sensor_stripPitch->AvailableForStates(G4State_PreInit);

  cmd_sensor_stripLength = new G4UIcmdWithADoubleAndUnit("/faser/geo/stripLength", this);
  cmd_sensor_stripLength->SetGuidance("Length of strips in the silicon wafers.");
  cmd_sensor_stripLength->SetParameterName("stripLength", true, true);
  cmd_sensor_stripLength->SetDefaultUnit("mm");
  cmd_sensor_stripLength->SetUnitCandidates("mm micron cm m");
  cmd_sensor_stripLength->AvailableForStates(G4State_PreInit);

  cmd_sensor_gap = new G4UIcmdWithADoubleAndUnit("/faser/geo/sensorGap", this);
  cmd_sensor_gap->SetGuidance("Gap between rows of strips in the silicon wafers.");
  cmd_sensor_gap->SetParameterName("sensorGap", true, true);
  cmd_sensor_gap->SetDefaultUnit("mm");
  cmd_sensor_gap->SetUnitCandidates("mm micron cm m");
  cmd_sensor_gap->AvailableForStates(G4State_PreInit);

  cmd_sensor_sizeZ = new G4UIcmdWithADoubleAndUnit("/faser/geo/sensorSizeZ", this);
  cmd_sensor_sizeZ->SetGuidance("Set the thickness of the silicon wafers.");
  cmd_sensor_sizeZ->SetParameterName("sensorSizeZ", true, true);
  cmd_sensor_sizeZ->SetDefaultUnit("mm");
  cmd_sensor_sizeZ->SetUnitCandidates("mm micron cm m");
  cmd_sensor_sizeZ->AvailableForStates(G4State_PreInit);

  cmd_sensor_stereoAngle = new G4UIcmdWithADoubleAndUnit("/faser/geo/stereoAngle", this);
  cmd_sensor_stereoAngle->SetGuidance("Set the angle wafers on each side are rotated.");
  cmd_sensor_stereoAngle->SetParameterName("sensorStereo", true, true);
  cmd_sensor_stereoAngle->SetDefaultUnit("mrad");
  cmd_sensor_stereoAngle->SetUnitCandidates("rad mrad deg");
  cmd_sensor_stereoAngle->AvailableForStates(G4State_PreInit);

  cmd_support_sizeZ = new G4UIcmdWithADoubleAndUnit("/faser/geo/supportSizeZ", this);
  cmd_support_sizeZ->SetGuidance("Set the thickness of the support strut.");
  cmd_support_sizeZ->SetParameterName("supportSizeZ", true, true);
  cmd_support_sizeZ->SetDefaultUnit("mm");
  cmd_support_sizeZ->SetUnitCandidates("mm micron cm m");
  cmd_support_sizeZ->AvailableForStates(G4State_PreInit);

  cmd_tracker_sensorPlanes = new G4UIcmdWithAnInteger("/faser/geo/trackerPlanes", this);
  cmd_tracker_sensorPlanes->SetGuidance("Number of silicon sensor planes in the tracker.");
  cmd_tracker_sensorPlanes->SetParameterName("sensorPlanes", true, true);
  cmd_tracker_sensorPlanes->AvailableForStates(G4State_PreInit);

  cmd_tracker_planePitch = new G4UIcmdWithADoubleAndUnit("/faser/geo/planePitch", this);
  cmd_tracker_planePitch->SetGuidance("Longitudinal spacing between sensor planes.");
  cmd_tracker_planePitch->SetParameterName("planePitch", true, true);
  cmd_tracker_planePitch->SetDefaultUnit("mm");
  cmd_tracker_planePitch->SetUnitCandidates("mm micron cm m");
  cmd_tracker_planePitch->AvailableForStates(G4State_PreInit);

  cmd_sampler_sensorPlanes = new G4UIcmdWithAnInteger("/faser/geo/samplerPlanes", this);
  cmd_sampler_sensorPlanes->SetGuidance("Number of silicon sensor planes in the calorimeter.");
  cmd_sampler_sensorPlanes->SetParameterName("samplerPlanes", true, true);
  cmd_sampler_sensorPlanes->AvailableForStates(G4State_PreInit);

  cmd_sampler_absorberW = new G4UIcmdWithADouble("/faser/geo/absorberW", this);
  cmd_sampler_absorberW->SetGuidance("Calorimeter absorber tungsten thickness in radiation lengths.");
  cmd_sampler_absorberW->SetParameterName("absorberW", true, true);
  cmd_sampler_absorberW->AvailableForStates(G4State_PreInit);

  cmd_sampler_absorberCu = new G4UIcmdWithADouble("/faser/geo/absorberCu", this);
  cmd_sampler_absorberCu->SetGuidance("Calorimeter absorber copper thickness in radiation lengths.");
  cmd_sampler_absorberCu->SetParameterName("absorberCu", true, true);
  cmd_sampler_absorberCu->AvailableForStates(G4State_PreInit);

  cmd_sampler_absorberC = new G4UIcmdWithADouble("/faser/geo/absorberC", this);
  cmd_sampler_absorberC->SetGuidance("Calorimeter absorber graphite thickness in radiation lengths.");
  cmd_sampler_absorberC->SetParameterName("absorberC", true, true);
  cmd_sampler_absorberC->AvailableForStates(G4State_PreInit);

  cmd_calo_planes = new G4UIcmdWithAnInteger("/faser/geo/caloPlanes", this);
  cmd_calo_planes->SetGuidance("Number of scintillator sensor planes in a calorimeter module.");
  cmd_calo_planes->SetParameterName("caloPlanes", true, true);
  cmd_calo_planes->AvailableForStates(G4State_PreInit);

  cmd_calo_towers = new G4UIcmdWithAnInteger("/faser/geo/caloTowers", this);
  cmd_calo_towers->SetGuidance("Number of scinillator towers in a calorimeter module.");
  cmd_calo_towers->SetParameterName("caloTowers", true, true);
  cmd_calo_towers->SetRange("caloTowers == 1 || caloTowers == 4 || caloTowers == 9");
  cmd_calo_towers->AvailableForStates(G4State_PreInit);

  cmd_calo_modules = new G4UIcmdWithAnInteger("/faser/geo/caloModules", this);
  cmd_calo_modules->SetGuidance("Number of calorimeter modules in the detector.");
  cmd_calo_modules->SetParameterName("caloModules", true, true);
  cmd_calo_modules->SetRange("caloModules == 0 || caloModules == 1 || caloModules == 4 || caloModules == 9");
  cmd_calo_modules->AvailableForStates(G4State_PreInit);

  cmd_calo_scintThickness = new G4UIcmdWithADoubleAndUnit("/faser/geo/caloScintThickness", this);
  cmd_calo_scintThickness->SetGuidance("Scintillator layer thickness in the calorimeter.");
  cmd_calo_scintThickness->SetParameterName("caloScintThickness", true, true);
  cmd_calo_scintThickness->SetDefaultUnit("mm");
  cmd_calo_scintThickness->SetUnitCandidates("mm micron cm m");
  cmd_calo_scintThickness->AvailableForStates(G4State_PreInit);

  cmd_calo_absorbThickness = new G4UIcmdWithADoubleAndUnit("/faser/geo/caloAbsorbThickness", this);
  cmd_calo_absorbThickness->SetGuidance("Absorber layer thickness in the calorimeter.");
  cmd_calo_absorbThickness->SetParameterName("caloAbsorbThickness", true, true);
  cmd_calo_absorbThickness->SetDefaultUnit("mm");
  cmd_calo_absorbThickness->SetUnitCandidates("mm micron cm m");
  cmd_calo_absorbThickness->AvailableForStates(G4State_PreInit);
  
  cmd_calo_planeXY = new G4UIcmdWithADoubleAndUnit("/faser/geo/caloPlaneXY", this);
  cmd_calo_planeXY->SetGuidance("Inner width of the calorimeter.");
  cmd_calo_planeXY->SetParameterName("caloPlaneXY", true, true);
  cmd_calo_planeXY->SetDefaultUnit("cm");
  cmd_calo_planeXY->SetUnitCandidates("mm micron cm m");
  cmd_calo_planeXY->AvailableForStates(G4State_PreInit);

  cmd_calo_moduleXY = new G4UIcmdWithADoubleAndUnit("/faser/geo/caloModuleXY", this);
  cmd_calo_moduleXY->SetGuidance("Outer width of the calorimeter.");
  cmd_calo_moduleXY->SetParameterName("caloModuleXY", true, true);
  cmd_calo_moduleXY->SetDefaultUnit("cm");
  cmd_calo_moduleXY->SetUnitCandidates("mm micron cm m");
  cmd_calo_moduleXY->AvailableForStates(G4State_PreInit);
  
  cmd_detector_samplerLength = new G4UIcmdWithADoubleAndUnit("/faser/geo/samplerLength", this);
  cmd_detector_samplerLength->SetGuidance("Thickness of calorimeter sampling layer.");
  cmd_detector_samplerLength->SetParameterName("samplerSizeZ", true, true);
  cmd_detector_samplerLength->SetDefaultUnit("um");
  cmd_detector_samplerLength->SetUnitCandidates("mm micron um cm m");
  cmd_detector_samplerLength->AvailableForStates(G4State_PreInit);

  cmd_detector_decayVolumeLength = new G4UIcmdWithADoubleAndUnit("/faser/geo/decayVolumeLength", this);
  cmd_detector_decayVolumeLength->SetGuidance("Air volume in front of first sensor.");
  cmd_detector_decayVolumeLength->SetParameterName("decayVolumeLength", true, true);
  cmd_detector_decayVolumeLength->SetDefaultUnit("m");
  cmd_detector_decayVolumeLength->SetUnitCandidates("mm micron cm m");
  cmd_detector_decayVolumeLength->AvailableForStates(G4State_PreInit);

  cmd_detector_trackerLength = new G4UIcmdWithADoubleAndUnit("/faser/geo/trackerLength", this);
  cmd_detector_trackerLength->SetGuidance("Length of the tracker region of the detector");
  cmd_detector_trackerLength->SetParameterName("trackerLength", true, true);
  cmd_detector_trackerLength->SetDefaultUnit("m");
  cmd_detector_trackerLength->SetUnitCandidates("mm micron cm m");
  cmd_detector_trackerLength->AvailableForStates(G4State_PreInit);

  cmd_detector_calorimeterLength = new G4UIcmdWithADoubleAndUnit("/faser/geo/calorimeterLength", this);
  cmd_detector_calorimeterLength->SetGuidance("Length of the calorimeter region of the detector");
  cmd_detector_calorimeterLength->SetParameterName("calorimeterLength", true, true);
  cmd_detector_calorimeterLength->SetDefaultUnit("m");
  cmd_detector_calorimeterLength->SetUnitCandidates("mm micron cm m");
  cmd_detector_calorimeterLength->AvailableForStates(G4State_PreInit);

  fDetectorConstruction->setReadoutStrips( FaserDetectorConstruction::default_sensor_readoutStrips );
  fDetectorConstruction->setStripPitch( FaserDetectorConstruction::default_sensor_stripPitch );
  fDetectorConstruction->setStripLength( FaserDetectorConstruction::default_sensor_stripLength );
  fDetectorConstruction->setSensorGap( FaserDetectorConstruction::default_sensor_gap );
  fDetectorConstruction->setSensorSizeZ ( FaserDetectorConstruction::default_sensor_sizeZ );
  fDetectorConstruction->setSensorStereoAngle( FaserDetectorConstruction::default_sensor_stereoAngle );
  fDetectorConstruction->setSupportSizeZ ( FaserDetectorConstruction::default_support_sizeZ );
  fDetectorConstruction->setSensorPlanes( FaserDetectorConstruction::default_tracker_sensorPlanes );
  fDetectorConstruction->setSamplerPlanes( FaserDetectorConstruction::default_sampler_sensorPlanes );
  fDetectorConstruction->setAbsorberW( FaserDetectorConstruction::default_sampler_absorberW );
  fDetectorConstruction->setAbsorberCu( FaserDetectorConstruction::default_sampler_absorberCu );
  fDetectorConstruction->setAbsorberC( FaserDetectorConstruction::default_sampler_absorberC );
  fDetectorConstruction->setCaloPlanes( FaserDetectorConstruction::default_calo_planes );
  fDetectorConstruction->setCaloTowers( FaserDetectorConstruction::default_calo_towers );
  fDetectorConstruction->setCaloModules( FaserDetectorConstruction::default_calo_modules );
  fDetectorConstruction->setCaloScintThickness( FaserDetectorConstruction::default_calo_scintThickness );
  fDetectorConstruction->setCaloAbsorbThickness( FaserDetectorConstruction::default_calo_absorbThickness );
  fDetectorConstruction->setCaloPlaneXY( FaserDetectorConstruction::default_calo_planeXY );
  fDetectorConstruction->setCaloModuleXY( FaserDetectorConstruction::default_calo_moduleXY );
  fDetectorConstruction->setSamplerSizeZ( FaserDetectorConstruction::default_detector_samplerLength );
  fDetectorConstruction->setPlanePitch( FaserDetectorConstruction::default_detector_planePitch );
  fDetectorConstruction->setDecayVolumeLength( FaserDetectorConstruction::default_detector_decayVolumeLength );
  fDetectorConstruction->setTrackerLength( FaserDetectorConstruction::default_detector_trackerLength );
  fDetectorConstruction->setCalorimeterLength( FaserDetectorConstruction::default_detector_calorimeterLength );
}

FaserGeometryMessenger::~FaserGeometryMessenger()
{
  if (cmd_detector_calorimeterLength) delete cmd_detector_calorimeterLength;
  if (cmd_detector_trackerLength) delete cmd_detector_trackerLength;
  if (cmd_detector_decayVolumeLength) delete cmd_detector_decayVolumeLength;
  if (cmd_tracker_planePitch) delete cmd_tracker_planePitch;
  if (cmd_tracker_sensorPlanes) delete cmd_tracker_sensorPlanes;
  if (cmd_sampler_sensorPlanes) delete cmd_sampler_sensorPlanes;
  if (cmd_sampler_absorberW) delete cmd_sampler_absorberW;
  if (cmd_sampler_absorberCu) delete cmd_sampler_absorberCu;
  if (cmd_sampler_absorberC) delete cmd_sampler_absorberC;
  if (cmd_detector_samplerLength) delete cmd_detector_samplerLength;
  if (cmd_support_sizeZ) delete cmd_support_sizeZ;
  if (cmd_sensor_stereoAngle) delete cmd_sensor_stereoAngle;
  if (cmd_sensor_gap) delete cmd_sensor_gap;
  if (cmd_sensor_sizeZ) delete cmd_sensor_sizeZ;
  if (cmd_sensor_stripPitch) delete cmd_sensor_stripPitch;
  if (cmd_sensor_stripLength) delete cmd_sensor_stripLength;
  if (cmd_sensor_readoutStrips) delete cmd_sensor_readoutStrips;
  if (cmd_calo_planes) delete cmd_calo_planes;
  if (cmd_calo_towers) delete cmd_calo_towers;
  if (cmd_calo_modules) delete cmd_calo_modules;
  if (cmd_calo_scintThickness) delete cmd_calo_scintThickness;
  if (cmd_calo_absorbThickness) delete cmd_calo_absorbThickness;
  if (cmd_calo_planeXY) delete cmd_calo_planeXY;
  if (cmd_calo_moduleXY) delete cmd_calo_moduleXY;
  if (geometryDirectory) delete geometryDirectory;
}

void FaserGeometryMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if (command == cmd_sensor_readoutStrips) 
  {
    fDetectorConstruction->setReadoutStrips( cmd_sensor_readoutStrips->GetNewIntValue(newValues) );
  }
  else if (command == cmd_sensor_stripPitch)
  {
    fDetectorConstruction->setStripPitch( cmd_sensor_stripPitch->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_sensor_stripLength)
  {
    fDetectorConstruction->setStripLength( cmd_sensor_stripLength->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_sensor_gap)
  {
    fDetectorConstruction->setSensorGap( cmd_sensor_gap->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_sensor_sizeZ)
  {
    fDetectorConstruction->setSensorSizeZ( cmd_sensor_sizeZ->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_sensor_stereoAngle)
  {
    fDetectorConstruction->setSensorStereoAngle( cmd_sensor_stereoAngle->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_support_sizeZ)
  {
    fDetectorConstruction->setSupportSizeZ( cmd_support_sizeZ->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_tracker_sensorPlanes)
  {
    fDetectorConstruction->setSensorPlanes( cmd_tracker_sensorPlanes->GetNewIntValue(newValues) );
  }
  else if (command == cmd_sampler_sensorPlanes)
  {
    fDetectorConstruction->setSamplerPlanes( cmd_sampler_sensorPlanes->GetNewIntValue(newValues) );
  }
  else if (command == cmd_sampler_absorberW)
  {
    fDetectorConstruction->setAbsorberW( cmd_sampler_absorberW->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_sampler_absorberCu)
  {
    fDetectorConstruction->setAbsorberCu( cmd_sampler_absorberCu->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_sampler_absorberC)
  {
    fDetectorConstruction->setAbsorberC( cmd_sampler_absorberC->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_calo_planes)
  {
    fDetectorConstruction->setCaloPlanes( cmd_calo_planes->GetNewIntValue(newValues) );
  }
  else if (command == cmd_calo_towers)
  {
    fDetectorConstruction->setCaloTowers( cmd_calo_towers->GetNewIntValue(newValues) );
  }
  else if (command == cmd_calo_modules)
  {
    fDetectorConstruction->setCaloModules( cmd_calo_modules->GetNewIntValue(newValues) );
  }
  else if (command == cmd_calo_scintThickness)
  {
    fDetectorConstruction->setCaloScintThickness( cmd_calo_scintThickness->GetNewDoubleValue(newValues));
  }
  else if (command == cmd_calo_absorbThickness)
  {
    fDetectorConstruction->setCaloAbsorbThickness( cmd_calo_absorbThickness->GetNewDoubleValue(newValues));
  }
  else if (command == cmd_calo_planeXY)
  {
    fDetectorConstruction->setCaloPlaneXY( cmd_calo_planeXY->GetNewDoubleValue(newValues));
  }
  else if (command == cmd_calo_moduleXY)
  {
    fDetectorConstruction->setCaloModuleXY( cmd_calo_moduleXY->GetNewDoubleValue(newValues));
  }
  else if (command == cmd_detector_samplerLength)
  {
    fDetectorConstruction->setSamplerSizeZ( cmd_detector_samplerLength->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_tracker_planePitch)
  {
    fDetectorConstruction->setPlanePitch( cmd_tracker_planePitch->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_detector_decayVolumeLength)
  {
    fDetectorConstruction->setDecayVolumeLength( cmd_detector_decayVolumeLength->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_detector_trackerLength)
  {
    fDetectorConstruction->setTrackerLength( cmd_detector_trackerLength->GetNewDoubleValue(newValues) );
  }
  else if (command == cmd_detector_calorimeterLength)
  {
    fDetectorConstruction->setCalorimeterLength( cmd_detector_calorimeterLength->GetNewDoubleValue(newValues) );
  }
}

G4String FaserGeometryMessenger::GetCurrentValue(G4UIcommand* command)
{
  G4String cv;
  
  if (command == cmd_sensor_readoutStrips)
  {
    cv = cmd_sensor_readoutStrips->ConvertToString(fDetectorConstruction->getReadoutStrips());
  }
  else if (command == cmd_sensor_stripPitch)
  {
    cv = cmd_sensor_stripPitch->ConvertToString(fDetectorConstruction->getStripPitch(), "mm");
  }
  else if (command == cmd_sensor_stripLength)
  {
    cv = cmd_sensor_stripLength->ConvertToString(fDetectorConstruction->getStripLength(), "mm");
  }
  else if (command == cmd_sensor_gap)
  {
    cv = cmd_sensor_gap->ConvertToString(fDetectorConstruction->getSensorGap(), "mm");
  }
  else if (command == cmd_sensor_sizeZ)
  {
    cv = cmd_sensor_sizeZ->ConvertToString(fDetectorConstruction->getSensorSizeZ(), "mm");
  }
  else if (command == cmd_sensor_stereoAngle)
  {
    cv = cmd_sensor_stereoAngle->ConvertToString(fDetectorConstruction->getSensorStereoAngle(), "mrad");
  }
  else if (command == cmd_support_sizeZ)
  {
    cv = cmd_support_sizeZ->ConvertToString(fDetectorConstruction->getSupportSizeZ(), "mm");
  }
  else if (command == cmd_tracker_sensorPlanes)
  {
    cv = cmd_tracker_sensorPlanes->ConvertToString(fDetectorConstruction->getSensorPlanes());
  }
  else if (command == cmd_sampler_sensorPlanes)
  {
    cv = cmd_sampler_sensorPlanes->ConvertToString(fDetectorConstruction->getSamplerPlanes());
  }
  else if (command == cmd_sampler_absorberW)
  {
    cv = cmd_sampler_absorberW->ConvertToString(fDetectorConstruction->getAbsorberW());
  }
  else if (command == cmd_sampler_absorberCu)
  {
    cv = cmd_sampler_absorberCu->ConvertToString(fDetectorConstruction->getAbsorberCu());
  }
  else if (command == cmd_sampler_absorberC)
  {
    cv = cmd_sampler_absorberC->ConvertToString(fDetectorConstruction->getAbsorberC());
  }
  else if (command == cmd_calo_planes)
  {
    cv = cmd_calo_planes->ConvertToString(fDetectorConstruction->getCaloPlanes());
  }
  else if (command == cmd_calo_towers)
  {
    cv = cmd_calo_towers->ConvertToString(fDetectorConstruction->getCaloTowers());
  }
  else if (command == cmd_calo_modules)
  {
    cv = cmd_calo_modules->ConvertToString(fDetectorConstruction->getCaloModules());
  }
  else if (command == cmd_calo_scintThickness)
  {
    cv = cmd_calo_scintThickness->ConvertToString(fDetectorConstruction->getCaloScintThickness(), "mm");
  }
  else if (command == cmd_calo_absorbThickness)
  {
    cv = cmd_calo_absorbThickness->ConvertToString(fDetectorConstruction->getCaloAbsorbThickness(), "mm");
  }
  else if (command == cmd_calo_planeXY)
  {
    cv = cmd_calo_planeXY->ConvertToString(fDetectorConstruction->getCaloPlaneXY(), "cm");
  }
  else if (command == cmd_calo_moduleXY)
  {
    cv = cmd_calo_moduleXY->ConvertToString(fDetectorConstruction->getCaloModuleXY(), "cm");
  }
  else if (command == cmd_detector_samplerLength)
  {
    cv = cmd_detector_samplerLength->ConvertToString(fDetectorConstruction->getSamplerSizeZ(), "um");
  }
  else if (command == cmd_tracker_planePitch)
  {
    cv = cmd_tracker_planePitch->ConvertToString(fDetectorConstruction->getPlanePitch(), "m");
  }
  else if (command == cmd_detector_decayVolumeLength)
  {
    cv = cmd_detector_decayVolumeLength->ConvertToString(fDetectorConstruction->getDecayVolumeLength(), "m");
  }
  else if (command == cmd_detector_trackerLength)
  {
    cv = cmd_detector_trackerLength->ConvertToString(fDetectorConstruction->getTrackerLength(), "m");
  }
  else if (command == cmd_detector_calorimeterLength)
  {
    cv = cmd_detector_calorimeterLength->ConvertToString(fDetectorConstruction->getCalorimeterLength(), "m");
  }
  return cv;
}
