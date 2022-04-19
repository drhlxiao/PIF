/***************************************************************
 * Class to define detector messeger.
 * Author  : Hualin Xiao
 * Date    : Jan, 2015
 * Version : 1.10
 *
 ***************************************************************/

#include "globals.hh"

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction* theDet)
: fDetector(theDet) { 
	
	fDetectorDir = new G4UIdirectory( "/PIF/degrader" );
	fDetectorDir->SetGuidance("Detector control.");

	//fSetFileCommand = new G4UIcmdWithAString("/PIF/degrader", this);
//	fSetFileCommand ->AvailableForStates(G4State_PreInit);
 // 
	fSetThickness= new G4UIcmdWithADouble("/PIF/degrader/thickness", this);
	fSetDistance= new G4UIcmdWithADouble("/PIF/degrader/distance", this);
	fSetCollimator= new G4UIcmdWithADouble("/PIF/collimator/radius", this);
}

DetectorMessenger::~DetectorMessenger()
{
//	delete fSetFileCommand;
//	delete fSetWorldCommand;
	delete fDetectorDir;
	delete fSetThickness;
	delete fSetDistance;
	delete fSetCollimator;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
	if (command == fSetThickness) { 
        fDetector->SetDegraderThickness(fSetThickness->GetNewDoubleValue(newValue));
	}
	if (command == fSetDistance) { 
        fDetector->SetDegraderDistance(fSetDistance->GetNewDoubleValue(newValue));
	}
	if (command == fSetCollimator) { 
        fDetector->SetCollimator(fSetCollimator->GetNewDoubleValue(newValue));
	}
}
