//
/// \file AnalysisManager.hh
/// \brief Definition of the AnalysisManager class

#ifndef AnalysisManager_h
#define AnalysisManager_h 1

#include "globals.hh"

//#include "G4Step.hh"
//#include "G4Event.hh"
//#include "G4Run.hh"
#include "TString.h"
class G4Run;
class G4Event;
class G4Step;

class TCanvas;
class TH1F;
class TFile;
class TTree;


const G4int maxTracks=40;

class AnalysisManager
{
    public:
        static AnalysisManager* GetInstance();
        static void Dispose();

        void OpenFile();
        void CreateTree();
        void WriteFile();
        void SetEventID(G4int eventid){fEventID=eventid;}
        void SetOutputFileName(TString filen){outputFilename=filen;}

       	void AddEnergy(G4double edep);
    	//void AddNonIonizingEnergy(G4double edep);


        void BeginOfEventAction(const G4Event *event);
        void EndOfEventAction(const G4Event *event);
        void SteppingAction(const G4Step *aStep);
        void BeginOfRunAction(const G4Run *);
        void EndOfRunAction(const G4Run *);
        void WriteMacros(G4String macfilename);
        
		void CloseFile();
		void EnableTree(bool v){enable_tree=v;}

        AnalysisManager();
        ~AnalysisManager();

    private:
        static AnalysisManager* fManager;

        TString outputFilename;

        TFile* fTFile;
        TTree* fTree;
        int ngammas;
		bool enable_tree;


        long long fEventsFill;
        G4bool fToRecord;
        G4bool recordFluence;

   //     G//4double fNonIonizingEdep[3];

        G4double particle_position[3];
        G4double particle_direction[3];
        G4double particle_energy;
        TH1F *hedep;
        TH1F *hin;
        TH1F *hspec[9];
        //momentum direction at MCP0, MCP1 

        G4double penergy;

        G4int fEventID;


    	G4double fEdepSum;
    	G4double fNonIonizingEdepSum;

        //track entering MCPs 
        G4int ntrack;
        G4double x[maxTracks];
        G4double y[maxTracks];
        G4double z[maxTracks];
        G4double dx[maxTracks];
        G4double dy[maxTracks];
        G4double dz[maxTracks];
        G4int pdg[maxTracks]; //particle name
        G4double energy[maxTracks]; //particle energy
        G4int parent_id[maxTracks]; //
        G4double time[maxTracks];
        G4int entering[maxTracks];
        G4int escaping[maxTracks];

		G4double TID;

};

#endif


