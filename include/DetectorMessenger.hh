//
/// \file /include/DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADouble;

class DetectorMessenger: public G4UImessenger 
{
     public:
     	     DetectorMessenger(DetectorConstruction*);
     	     ~DetectorMessenger();
    
     	     void SetNewValue( G4UIcommand*, G4String );
     	     
     private:
     	     DetectorConstruction* fDetector;
     	     G4UIdirectory*             fDetectorDir;
     	   //  G4UIcmdWithAString*        fSetFileCommand;
     	     G4UIcmdWithADouble*        fSetThickness;
     	     G4UIcmdWithADouble*        fSetDistance;
     	     G4UIcmdWithADouble*        fSetCollimator;
};

#endif
