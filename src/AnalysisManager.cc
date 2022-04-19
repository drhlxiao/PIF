/***************************************************************
 * descriptions of AnalysisManager
 * Author  : Hualin Xiao (hualin.xiao@psi.ch)
 * Date    : Jan., 2016
 * Version : 1.10
 ***************************************************************/
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TH1F.h"
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TNamed.h"

#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"


#include "PrimaryGeneratorAction.hh"

#include "AnalysisManager.hh"
#include <fstream>
using namespace std;


void CanRebin(TH1* h){

#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
	h->SetBit(TH1::kCanRebin);
#else
	h->SetCanExtend(TH1::kXaxis);
#endif
}



AnalysisManager* AnalysisManager::fManager = 0;
AnalysisManager* AnalysisManager::GetInstance()
{
	if(!fManager) {
		fManager = new AnalysisManager();
	}
	return fManager;
}	
AnalysisManager::AnalysisManager() {

	TID=0;
	enable_tree=true;

}

void AnalysisManager::CreateTree() 
{
	if(enable_tree)
	{
		fTree=new TTree("event","event");
		fTree->Branch("edep",&fEdepSum,"edep/D");    
		fTree->Branch("event_id",&fEventID,"event_id/I");
		fTree->Branch("source_position",particle_position,"source_position[3]/D");    
		fTree->Branch("source_direction",particle_direction,"source_direction[3]/D");    
		fTree->Branch("source_energy",&particle_energy,"source_energy/D");    
		fTree->Branch("ntrack",&ntrack,"ntrack/I");    
		fTree->Branch("x",x,"x[ntrack]/D");    
		fTree->Branch("y",y,"y[ntrack]/D");    
		fTree->Branch("z",z,"z[ntrack]/D");    
		//fTree->Branch("entering",entering,"entering[ntrack]/I");    
		//fTree->Branch("escaping",escaping,"escaping[ntrack]/I");    
		fTree->Branch("px",dx,"px[ntrack]/D");    
		fTree->Branch("py",dy,"py[ntrack]/D");    
		fTree->Branch("pz",dz,"pz[ntrack]/D");    
		fTree->Branch("pdg",pdg,"pdg[ntrack]/I");    

		fTree->Branch("energy",energy,"energy[ntrack]/D");    
		fTree->Branch("parent_id",parent_id,"parent_id[ntrack]/I");    
		fTree->Branch("time",time,"time[ntrack]/D");    
		fTree->SetAutoSave(1000);
	}

	G4String particles_name[9]={"protons","neutrons","gamma-rays","e-", "e+", 
		"mu-","pi-", "other particles","unused"};
	G4String xyz[3]={"x","y","z"};

	hedep=new TH1F("hedep","Energy deposition ; Energy deposition (MeV); Counts",100,0,1);
	hin=new TH1F("hin","Energies at DUT ; Energy (MeV); Counts",100,0,1);
	//hedep->SetBit(TH1::kCanRebin);
	//hedep->SetBit(TH1::kCanRebin);
	CanRebin(hedep);
	CanRebin(hin);




}


//////////////////////////////////////////////////////////////////////////

void AnalysisManager::BeginOfRunAction(const G4Run *run)
{

	OpenFile();
	CreateTree();
}

void AnalysisManager::EndOfRunAction(const G4Run *run)
{
	WriteFile();
}

void AnalysisManager::BeginOfEventAction(const G4Event *event)
{
	fEdepSum = 0.0;
	fToRecord=false;

	ntrack=0;
	for(int i=0;i<maxTracks;i++)
	{
		x[i]=y[i]=z[i]=-1000;
		dx[i]=dy[i]=dz[i]=1000;
		pdg[i]=0;
		parent_id[i]=-1;
		time[i]=0;
	}
}
void AnalysisManager::EndOfEventAction(const G4Event *event)
{



	fEventID=event->GetEventID();

	G4RunManager* runManager = G4RunManager::GetRunManager();
	PrimaryGeneratorAction* primaryAction = (PrimaryGeneratorAction*) runManager->GetUserPrimaryGeneratorAction();
	G4ThreeVector position,direction;

	primaryAction->GetGPS(position,direction,penergy);
	particle_position[0]=position.getX();
	particle_position[1]=position.getY();
	particle_position[2]=position.getZ();
	particle_direction[0]=direction.getX();
	particle_direction[1]=direction.getY();
	particle_direction[2]=direction.getZ();
	particle_energy=penergy/1000.;
	//keV to MeV


	if(fEdepSum>0)hedep->Fill(fEdepSum/1000.);

	if(enable_tree)
	{
		if(fToRecord)fTree->Fill();
	}


}
//Stepping Action

void AnalysisManager::SteppingAction(const G4Step *aStep)
{

	const G4Track* track = aStep->GetTrack();
	G4String volName; 

	if (track->GetVolume()) volName =  track->GetVolume()->GetName(); 

	G4StepPoint* point1 = aStep->GetPreStepPoint();
	G4StepPoint* point2 = aStep->GetPostStepPoint();

	bool is_in_det=false;
	if(volName=="siPhys")
	{
		is_in_det=true;
		G4double edep;
		edep=aStep->GetTotalEnergyDeposit()/keV;
		if (edep > 0.0) 
		{
			AddEnergy(edep);
		}	
	}
	else
	{
		//G4cout<<volName<<G4endl;
		return ;
	}

	if(ntrack>=(maxTracks-1))return;


	if (point1->GetStepStatus() == fGeomBoundary) 
	{
		G4int parent=track->GetParentID(); //
		G4int particle_id=track->GetDefinition()->GetPDGEncoding();
		G4double kinetic_energy=point1->GetKineticEnergy()/MeV;

		fToRecord=true;
		x[ntrack]=track->GetPosition().x()/mm;
		y[ntrack]=track->GetPosition().y()/mm;
		z[ntrack]=track->GetPosition().z()/mm;


		pdg[ntrack]=particle_id;

		energy[ntrack]=kinetic_energy;
		hin->Fill(kinetic_energy);
		parent_id[ntrack]=parent;

		ntrack++;
	}

}






inline void AnalysisManager::AddEnergy(G4double edep) 
{
	fEdepSum += edep;
}



AnalysisManager::~AnalysisManager()
{
	//if(fTree)delete fTree;
	//if(fTFile)delete fTFile;
	//if(fManager)delete fManager;
	//fManager = 0;

	fTFile->Close();
}
void AnalysisManager::OpenFile()
{
	fTFile = new TFile(outputFilename.Data(),"recreate");
}
void AnalysisManager::WriteFile() 
{
	G4cout<<"writing histograms"<<G4endl;
	fTFile->cd();
	hedep->Write();
	hin->Write();
	if(enable_tree)fTree->Write();
	G4cout<<"histograms have been written"<<G4endl;

}


void AnalysisManager::WriteMacros(G4String macfilename)
{
}

void AnalysisManager::CloseFile()
{
	if(fTFile)fTFile->Close();
}

