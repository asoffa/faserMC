#ifndef FaserSamplerHit_h
#define FaserSamplerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "FaserDetectorConstruction.hh"
#include "tls.hh"

// Hit class for the Faser silicon sensor
//

class FaserSamplerHit : public G4VHit
{
   public:

    FaserSamplerHit();
    FaserSamplerHit(const FaserSamplerHit&);
    virtual ~FaserSamplerHit();

    const FaserSamplerHit& operator=(const FaserSamplerHit&);
    G4int operator==(const FaserSamplerHit&);

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // over-ridden methods
    virtual void Draw();
    virtual void Print();

    // setters
    void SetPlane(G4int plane)               		{ fPlaneID = plane; };
    void SetEdep(G4double de)           	       	{ fEdep = de; };
    void SetGlobalPos(G4ThreeVector globalXYZ) 		{ fGlobalPos = globalXYZ; };
    void SetLocalPos(G4ThreeVector localXYZ)   		{ fLocalPos = localXYZ; };
  
    void SetOriginTrack(G4int track)                    { fOriginTrackID = track; }
    void SetSourceTrack(G4int track)            { fSourceTrackID = track; }

    // getters
    G4int Plane() const    	       		        { return fPlaneID; };
    G4double Edep() const           			{ return fEdep; };
    G4ThreeVector GlobalPos() const 			{ return fGlobalPos; };
    G4ThreeVector LocalPos() const	  		{ return fLocalPos; };


    G4int OriginTrack() const			        { return fOriginTrackID; }
    G4int SourceTrack() const             { return fSourceTrackID; }

   private:

    G4int fPlaneID;
    G4double fEdep;
    G4ThreeVector fGlobalPos;
    G4ThreeVector fLocalPos;

    G4int fOriginTrackID;
    G4int fSourceTrackID;

    static const FaserDetectorConstruction* fDetectorConstruction;
};

typedef G4THitsCollection<FaserSamplerHit> FaserSamplerHitsCollection;

extern G4ThreadLocal G4Allocator<FaserSamplerHit>* FaserSamplerHitAllocator;

inline void* FaserSamplerHit::operator new(size_t)
{
  if (!FaserSamplerHitAllocator)
    FaserSamplerHitAllocator = new G4Allocator<FaserSamplerHit>;
  return (void *) FaserSamplerHitAllocator->MallocSingle();
}

inline void FaserSamplerHit::operator delete(void* hit)
{
  FaserSamplerHitAllocator->FreeSingle((FaserSamplerHit*) hit);
}

#endif
