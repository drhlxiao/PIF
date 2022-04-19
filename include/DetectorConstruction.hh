

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

// STL //
#include <string>

// GEANT4 //
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
#include "DetectorMessenger.hh"

#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

    void SetVisAttrib(G4LogicalVolume *log, G4double red, G4double green, G4double blue, G4double alpha);
	void SetDegraderThickness(G4double v){Degrade_thickness=v;}
	void SetDegraderDistance(G4double v){Degrade_distance=v;}
	void SetCollimator(G4double v){Collimator_radius=v;}
	G4VPhysicalVolume *CreateDegrader(G4double &width, G4double &height, G4double &thickness, G4double &posX, G4double &poY, G4double &posZ);
	G4VPhysicalVolume *CreateCollimatorRing(G4double &innerR, G4double &outerR, G4double &thickness, G4double &posX, G4double &poY, G4double &posZ);


  private:
    G4VSolid * world_solid;
    G4LogicalVolume* world_logical;
    G4VPhysicalVolume* world_physical;
    
	G4double Degrade_thickness;
	G4double Degrade_distance;
	G4double Collimator_radius;
	
    
    

    DetectorMessenger *det_msg;
};

#endif

