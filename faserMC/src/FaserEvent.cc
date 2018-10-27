#include "FaserEvent.hh"

#include "G4ParticleTable.hh"
#include "G4VTrajectoryPoint.hh"

#include <tuple>

////////////////////////////////////////////////////////////////////////////////


FaserEvent::FaserEvent(G4int eventNumber)
  : fEventNumber(eventNumber)
{}

//------------------------------------------------------------------------------

FaserEvent::FaserEvent()
  : FaserEvent(0)
{ }

//------------------------------------------------------------------------------

FaserEvent::~FaserEvent()
{ 
  for (auto p : fParticles) if (p) delete p;
  fParticles.clear();

  for (auto c : fClusters) if (c) delete c;
  fClusters.clear();

  for (auto s : fSpacePoints) if (s) delete s;
  fSpacePoints.clear();
}

//------------------------------------------------------------------------------

void FaserEvent::SetParticles(G4TrajectoryContainer* particles)
{
  if (particles == nullptr) return;
  for (size_t i = 0; i < particles->size(); i++)
  {
    G4VTrajectory& t = *((*particles)[i]);
    G4int pdgCode = t.GetPDGEncoding();
    G4double mass = G4ParticleTable::GetParticleTable()->FindParticle(pdgCode)->GetPDGMass();
    G4ThreeVector momentum = t.GetInitialMomentum();
    G4double energy = sqrt(momentum.mag2() + mass*mass);
    FaserTruthParticle* p = new FaserTruthParticle(t.GetTrackID(), t.GetParentID(), pdgCode,
						   (*(t.GetPoint(0))).GetPosition(), momentum, energy); 
    fParticles.push_back(p);
  }
}

//------------------------------------------------------------------------------

void FaserEvent::SetClusters()
{
  fClusters.clear();

  // If there are no digits, return empty `fClusters`
  if (fDigis.size() == 0) return;

  // Collect the digits by planes of the detector
  map<int, vector<FaserDigi*>> planeMap = mapDigitsByPlane(*this);

  // Process each plane in turn
  map<int, map<int, vector<FaserDigi*>>> rowMap;
  for (auto & plane : planeMap) 
  {
    rowMap[plane.first] = mapDigitsByRow(plane.second);
  }

  int clusterNumber = 0;
  // cluster the digits in each row of a plane
  map<int, vector<FaserCluster*>> clusterMap;
  for ( auto& plane : rowMap ) {
    for ( auto& row : plane.second ) {
      vector<vector<FaserDigi*>> clusters = clusterOneRow(row.second);
      for ( auto& c : clusters ) {
        clusterMap[plane.first].push_back(new FaserCluster(clusterNumber++, c));
      }
    }
  }

  for ( auto& plane : clusterMap ) {
    for ( auto c : plane.second ) {
      fClusters.push_back( c );
    }
  }

}

//------------------------------------------------------------------------------

void FaserEvent::SetSpacePoints()
{
  fSpacePoints.clear();

  // First find track IDs of e-/e+ or mu-/mu+ decay from dark photon
  vector<int> decayTrackIDs;
  for (FaserTruthParticle * particle : fParticles) {
    if (particle->ParentID() == 1) {
      decayTrackIDs.push_back(particle->TrackID());
    }
  }

  // Look for truth hits w/ the same track ID as truth particles having parent ID 1
  vector<FaserSensorHit*> hitsToCheck;
  for (FaserSensorHit * hit : fHits) {
    if (std::find(decayTrackIDs.begin(), decayTrackIDs.end(), hit->Track()) != decayTrackIDs.end()) {
      hitsToCheck.push_back(hit);
    }
  }

  // Group clusters on opposite sides of each module in the same row containing at least
  // one digit on the same strip as a truth hit having in the above list
  map<int, FaserCluster*> clusterMap_sensor0;
  map<int, FaserCluster*> clusterMap_sensor1;
  for (FaserCluster * cl : fClusters) {
    for (FaserDigi * digi : cl->Digis()) {
      for (FaserSensorHit * hit : hitsToCheck) {
        if (hit->Plane()  != digi->Plane() ) continue;
        if (hit->Module() != digi->Module()) continue;
        if (hit->Sensor() != digi->Sensor()) continue;
        if (hit->Row()    != digi->Row()   ) continue;
        if (hit->Strip()  != digi->Strip() ) continue;
        if (std::find(decayTrackIDs.begin(), decayTrackIDs.end(), hit->Track()) == decayTrackIDs.end()) continue;

        if (hit->Sensor() == 0) clusterMap_sensor0[hit->Track()] = cl;
        else if (hit->Sensor() == 1) clusterMap_sensor1[hit->Track()] = cl;
        goto endClusterHitCheck;
      }
    }
    endClusterHitCheck: ;
  }

  // Pair front/back clusters to form space points
  for (auto & it0 : clusterMap_sensor0) {
    for (auto & it1 : clusterMap_sensor1) {
      if (it0.first != it1.first) continue;
      FaserCluster * cl0 = it0.second;
      FaserCluster * cl1 = it1.second;
      if (cl0->Plane()  != cl1->Plane() ) continue;
      if (cl0->Module() != cl1->Module()) continue;
      if (cl0->Sensor() == cl1->Sensor()) continue; // opposite sides so skip if `==` here
      if (cl0->Row()    != cl1->Row()) continue;
      auto sp = new FaserSpacePoint;
      sp->AddCluster(cl0);
      sp->AddCluster(cl1);
      fSpacePoints.push_back(sp);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// Utility methods for cluster finding                                        //
////////////////////////////////////////////////////////////////////////////////

int FaserEvent::rowID(FaserDigi* d)
{
  return d->Row() + 2 * ( d->Sensor() + 4 * ( d->Module() + 2 * d->Plane() ) );
}

//------------------------------------------------------------------------------

vector<vector<FaserDigi*>> FaserEvent::clusterOneRow(vector<FaserDigi*> digits)
{
  vector<vector<FaserDigi*>> clusters;
  vector<FaserDigi*> cluster;
  int previousStrip = 1 << 30;
  for (auto d : digits)
  {
    if (cluster.size() > 0 && d->Strip() > previousStrip + 1)
    {
        clusters.push_back(cluster);
        cluster.clear();
    }
    cluster.push_back(d);
    previousStrip = d->Strip();
  }
  if (cluster.size() > 0) clusters.push_back(cluster);
  return clusters;
}

//------------------------------------------------------------------------------

void FaserEvent::sortDigits(vector<FaserDigi*>& v)
{
    std::sort(v.begin(), v.end(), [](FaserDigi* const& l, FaserDigi* const& r) { return l->Strip() < r->Strip(); });
}

//------------------------------------------------------------------------------

map<int, vector<FaserDigi*>> FaserEvent::mapDigitsByPlane(FaserEvent & evt)
{
  map<int, vector<FaserDigi*>> planeMap;
  for (FaserDigi* d : evt.Digis())
  {
    planeMap[d->Plane()].push_back( d );
  }
  return planeMap;
}

//------------------------------------------------------------------------------

map<int, vector<FaserDigi*>> FaserEvent::mapDigitsByRow(vector<FaserDigi*> v)
{
  map<int, vector<FaserDigi*>> rowMap;

  for (FaserDigi* d : v)
  {
    rowMap[rowID(d)].push_back( d );
  }
  for (auto& row : rowMap)
  {
    sortDigits(row.second);
  }
  return rowMap;
}
