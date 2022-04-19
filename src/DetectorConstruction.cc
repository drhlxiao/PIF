
// USER //
#include "DetectorConstruction.hh"


// GEANT4 //
#include "globals.hh"
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>

#include <G4Element.hh>
#include <G4Material.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Cons.hh>
#include <G4Sphere.hh>
#include <G4LogicalVolume.hh>
#include <G4Polycone.hh>
#include <G4PVReplica.hh>
#include <G4PVPlacement.hh>
#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <G4AssemblyVolume.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4Trd.hh>

#include <G4MaterialTable.hh>
#include <G4NistManager.hh>
#include <G4ElementTable.hh>
#include <G4RotationMatrix.hh>
#include <G4Transform3D.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>

#include <G4UImanager.hh>
#include <G4RunManager.hh>
#include <G4GeometryManager.hh>
#include <G4SolidStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>

#include <G4Sphere.hh>
#include <G4GDMLParser.hh>

#include "G4NistManager.hh"
#include "G4VisAttributes.hh"


#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction()
{
	Degrade_thickness=30*mm;
	Degrade_distance=10*cm;
	Collimator_radius=10*cm;

    det_msg=new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
    delete det_msg;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4bool overlap_check=true;
    G4NistManager * nist_manager = G4NistManager::Instance();
    G4Material * Air = nist_manager->FindOrBuildMaterial("G4_AIR");
    G4Material * Iron= nist_manager->FindOrBuildMaterial("G4_Fe");
    G4Material * Alum= nist_manager->FindOrBuildMaterial("G4_Al");
    G4Material * Bi= nist_manager->FindOrBuildMaterial("G4_Bi");
    G4Material * Ta= nist_manager->FindOrBuildMaterial("G4_Ta");
    G4Material * Copper= nist_manager->FindOrBuildMaterial("G4_Cu");
    G4Material * Silicon= nist_manager->FindOrBuildMaterial("G4_Si");
    G4Material * Vacuum= nist_manager->FindOrBuildMaterial("G4_Galactic");
 
    G4Element *elCu=nist_manager->FindOrBuildElement("Cu");
    G4Element *elTa=nist_manager->FindOrBuildElement("Ta");
 	G4double density=1e30*g/cm3;
    G4Material *blackHole= new G4Material("blackHole",density, 1);
    blackHole->AddElement(elCu, 1);



    //world
    world_solid = new G4Box("world_solid", 200*cm, 200*cm, 200*cm);
    world_logical = new G4LogicalVolume(world_solid, Vacuum,"world_logical",0,0,0);
    world_physical = new G4PVPlacement(0, G4ThreeVector(), world_logical, "world_physical", 0, false, 0);
    world_logical->SetVisAttributes(G4VisAttributes::Invisible);

	


	//construct degraders
	G4double degraderThickness[6]={0.1, 0.2, 0.4, 0.8, 1.6, 3.2};//3.2,1.6, 0.8,0.4,0.2,0.1};
	G4double degraderZ[6]={1.4, 1.4*2, 1.4*3, 1.4*4, 1.4*5, 1.4*7};
	G4double originToFirstDegrader=10*cm;
	// (0, 0, 0) is at 10 cm to the first degrader
	//change it if you want
	G4Box *degrader;
    G4LogicalVolume *degraderLog;
	G4VPhysicalVolume *degraderPhys[6];
	for(G4int i=0;i<6;i++){
		degrader=new G4Box("degrader",50*mm,50*mm,degraderThickness[i]*cm/2.);
		degraderLog=new G4LogicalVolume(degrader,Copper,"degrader",0,0,0);
		degraderPhys[i]=new G4PVPlacement(0, G4ThreeVector(0,0,degraderZ[i]*cm + originToFirstDegrader), degraderLog,
            "degraderPhys", world_logical, false, i);
	}

	//construct collimators
	G4double collimatorThickness[3]={0.5, 0.3, 0.1};
	G4double collimatorInnerRadius[3]={2, 3, 4};
	G4double collimatorOuterRadius[3]={3, 4, 5};
	G4double collimatorToOrigin = 10*cm + 8*1.4*cm +  originToFirstDegrader;
	//change the distance if you want

	G4Tubs *collimator;
    G4LogicalVolume *collimatorLog;
	G4VPhysicalVolume *collimatorPhys[3];

	for(G4int i=0;i<3;i++){
		collimator=new G4Tubs("collimator",collimatorInnerRadius[i]*cm,
				collimatorOuterRadius[i]*cm,collimatorThickness[i]*cm/2., 0, 360*deg);
		collimatorLog=new G4LogicalVolume(collimator,Copper,"collimator",0,0,0);
		collimatorPhys[i]=new G4PVPlacement(0, G4ThreeVector(0,0,-collimatorThickness[i]*cm/2 + collimatorToOrigin), collimatorLog,
            "collimatorPhys", world_logical, false, i);
	}
	

	//construct screen
	G4double siDetectorToLensDistance=20*cm;

    G4Box *SiDet=new G4Box("Si",50*mm,50*mm,0.03*mm/2);
    G4LogicalVolume *SiLog=new G4LogicalVolume(SiDet,Silicon,"siliconLog",0,0,0);
    G4VPhysicalVolume *siPhys=new G4PVPlacement(0, G4ThreeVector(0,0,collimatorToOrigin+siDetectorToLensDistance), SiLog,
            "siPhys", world_logical, false, 0);


    return world_physical;
}

void DetectorConstruction::SetVisAttrib(G4LogicalVolume *log, G4double red, G4double green, G4double blue, G4double alpha)
{
      G4VisAttributes *visAttrib=new G4VisAttributes(G4Colour(red,green, blue,alpha));

   	visAttrib->SetForceWireframe(true);

      visAttrib->SetForceSolid(true);

      log->SetVisAttributes(visAttrib);
}


