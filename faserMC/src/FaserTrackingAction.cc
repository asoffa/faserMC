#include "FaserTrackingAction.hh"
#include "FaserTrackInformation.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

FaserTrackingAction::FaserTrackingAction()
  : G4UserTrackingAction()
{;}

void FaserTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  FaserTrackInformation* info = dynamic_cast<FaserTrackInformation*>(aTrack->GetUserInformation());
	if (info->GetTrackingStatus() > 0)
	{
		fpTrackingManager->SetStoreTrajectory(true);
	}
	else 
	{
		fpTrackingManager->SetStoreTrajectory(false);
	}
  // if (info->GetSourceTrackID() != 0 && info->GetSourceTrackID() != aTrack->GetTrackID()) 
	// {
	// 	fpTrackingManager->SetStoreTrajectory(false);
	// }
	// else
	// {
	// 	fpTrackingManager->SetStoreTrajectory(true);
	// }
}

void FaserTrackingAction::PostUserTrackingAction(const G4Track*aTrack)
{
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if (secondaries != nullptr)
  {
      FaserTrackInformation* info = (FaserTrackInformation*) (aTrack->GetUserInformation());
      size_t nSeco = secondaries->size();
      if (nSeco >0)
      {
				for (size_t i = 0; i < nSeco; i++)
				{
						FaserTrackInformation* infoNew = nullptr;
						// kind of a kludge, but treat immediate daughters of primary particles as primaries                                        
						if (aTrack->GetParentID() == 0)
						{
							infoNew = new FaserTrackInformation((*secondaries)[i]);
							// G4cout << "Created new trackinfo for daughter of primary at z = " << aTrack->GetPosition().z() << G4endl;
							// infoNew->Print();
						}
						else
						{
							infoNew = new FaserTrackInformation(info);
							// G4cout << "Copied trackinfo to daughter of secondary at z= " << aTrack->GetPosition().z() << G4endl;
							// infoNew->Print();
						}
						(*secondaries)[i]->SetUserInformation(infoNew);
				}
      }
  }
}
