#ifndef FASER_TRACKER_GEOMETRY_HH
#define FASER_TRACKER_GEOMETRY_HH 1

#include <vector>
#include <string>

//------------------------------------------------------------------------------

class FaserTrackerGeometry {

public:

  int     nStrips       = -10000;
  double  stripPitch    = -10000.; // mm
  double  stereoAngle   = -10000.; // rad
  double  moduleOffsetX = -10000.; // mm
  double  sensorOffsetY = -10000.; // mm
  double  rowOffsetY    = -10000.; // mm
  std::vector<double> planeZ; // mm
  std::vector<double> * p_planeZ = &planeZ;

  virtual ~FaserTrackerGeometry() {
  }

  void WriteToFile(const std::string & fileName);
};

#endif
