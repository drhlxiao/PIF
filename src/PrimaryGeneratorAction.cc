/***************************************************************
 * descriptions of partcile gun
 * Author  : Hualin Xiao
 * Date    : Feb., 2015
 * Version : 1.10
 ***************************************************************/

#include "TFile.h"
#include "TTree.h"
#include "t2sim.h"

#include "PrimaryGeneratorAction.hh"
// For Random Generator
#include "Randomize.hh"
#include <CLHEP/Random/RandFlat.h>
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction():
	sourceType(0),fTree(NULL),fFile(NULL),nEntries(0),iEntry(0)
{
	//fParticleGunMessenger = new PrimaryGeneratorMessenger(this);


	E0=0;

	//particle source
	G4int n_particle = 1;
	//fParticleGun  = new G4ParticleGun(n_particle);     
	// default particle kinematic
	particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	G4ParticleDefinition* particle
		= particleTable->FindParticle(particleName="proton");
	//fParticleGun->SetParticleDefinition(particle);
	//fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
	//fParticleGun->SetParticleEnergy(200.*MeV);
	//fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,10*cm));
	//particle source
	fParticleSource  = new G4GeneralParticleSource();

}
PrimaryGeneratorAction::~PrimaryGeneratorAction() 
{
	delete fParticleSource;
	//delete fParticleGunMessenger;
	delete fParticleGun;
}

void PrimaryGeneratorAction::GetGPS(G4ThreeVector &position, G4ThreeVector &direction, G4double &energy)
{
	position=fParticleSource->GetParticlePosition(); 
	direction=fParticleSource->GetParticleMomentumDirection();
	energy=fParticleSource->GetParticleEnergy()/keV;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) 
{
	fParticleSource->GeneratePrimaryVertex(anEvent); 

}

G4bool PrimaryGeneratorAction::InitFile()
{
}

