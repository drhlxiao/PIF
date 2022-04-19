//
/// \file g4
/// \brief Main program for MCP
//   author: Hualin Xiao (hualin.xiao@psi.ch)
//   History:
//
//   Jan. 01st, 2016:
//   
//

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"

#include "G4UImanager.hh"
#include "G4VUserPhysicsList.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BIC.hh"
//#include "FTFP_BERT.hh"
//FTFP_BERT.hh"
#include "time.h"

#include "G4PhysListFactory.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "AnalysisManager.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif





void Help() 
{
    G4cout<<"PIF GEANT4 Simulation package"<<G4endl;
        G4cout<<"Usage:"<<G4endl
        <<"./g4PIF [OPTIONS] -m run.mac  -o OUTPUT"<<G4endl;
    G4cout<<"Options:"<<G4endl
        <<" -i                  <input.root> "<<G4endl
        <<"                     Read source information for a root file."<<G4endl
        <<"                     Its structure is defined in t2sim.h ."<<G4endl
        <<"-td                  degraderThickness in unit of mm  ."<<G4endl
        <<"-E0                  Initial Energy."<<G4endl
        <<" -h                  Print help information"<<G4endl;
}

int main(int argc,char **argv) 
{

    G4String outputFilename="run0.root";
    G4String macFilename;
    G4String inputFilename;

	G4double td=0;

    if(argc==1)Help();
    
    int s=0;
    G4String sel;
	double E0=0;

    while(s<argc-1)
    {

        sel=argv[++s];
        if(sel=="-h"||sel=="--help"||sel=="--h")
        {
            Help();
            return 0;

        }
        else if(sel=="-o")
        {
            outputFilename=argv[++s];
            if(!outputFilename.contains(".root"))
            {
                Help();
                return 0;
            }
        }
        else if(sel=="-td")
        {
			td=atof(argv[++s]);
			G4cout<<"degrader thickness:"<<td<<G4endl;
		}
        else if(sel=="-E0")
        {
			E0=atof(argv[++s]);
			G4cout<<"E0 (MeV):"<<E0<<G4endl;
		}



        else if(sel=="-m")
        {
            macFilename=argv[++s];
            if(!macFilename.contains(".mac"))
            {
                Help();
                return 0;
            }
        }
        else if(sel=="-i")
        {
            inputFilename=argv[++s];
            if(!inputFilename.contains(".root"))
            {
                Help();
                return 0;
            }
        }

    }

	if(E0>0&&td>0)
	{
		char ss[128];
		sprintf(ss,"Simout_%.1fMeV_%.1fmm.root",E0,td);
		
		outputFilename=ss;
	}

    time_t DateTime = time( NULL );
    CLHEP::HepRandom::setTheSeed( DateTime );

    AnalysisManager* analysisManager = AnalysisManager::GetInstance();
    //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);	
    analysisManager->SetOutputFileName(outputFilename);

#ifdef G4MULTITHREADED  
    G4MTRunManager * runManager = new G4MTRunManager();
    G4int nThreads = G4Threading::G4GetNumberOfCores()-1;
    runManager->SetNumberOfThreads(nThreads);
#else
    G4RunManager * runManager = new G4RunManager(); 
#endif


	DetectorConstruction *detConstruction=new DetectorConstruction();
	if(td>0)detConstruction->SetDegraderThickness(td);
    runManager->SetUserInitialization(detConstruction);




   G4String plname = "QGSP_BIC";  // set however you like ...
   G4PhysListFactory factory;
   G4VModularPhysicsList* physlist = factory.GetReferencePhysList(plname);

    runManager->SetUserInitialization(physlist);


    G4cout<<"Initializing primary generation"<<G4endl;
    PrimaryGeneratorAction *primarygen=new PrimaryGeneratorAction();
	if(E0>0)
	{
		primarygen->SetInitialEnergy(E0);
	}
    runManager->SetUserAction(primarygen);
    if(inputFilename!="")primarygen->SetParticleSourceFile(inputFilename);

    G4cout<<"Initializing run action"<<G4endl;
    runManager->SetUserAction(new RunAction());
    EventAction* evtAction = new EventAction();
    G4cout<<"Initializing event action"<<G4endl;
    runManager->SetUserAction(evtAction);
    G4cout<<"Initializing stepping action"<<G4endl;
    runManager->SetUserAction(new SteppingAction());


#ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif


    if (macFilename=="") 
    {
#ifdef G4UI_USE
        G4UIExecutive * ui = new G4UIExecutive(argc,argv,"qt");
        ui->SessionStart();
        delete ui;
#endif
    } else {      
          G4UIExecutive * ui = new G4UIExecutive(argc,argv,"tcsh");

        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        G4String command = "/control/execute "; 
        G4String fileName = macFilename; 
        G4cout<<"Applying command: "<<command+fileName<<G4endl;
        UImanager->ApplyCommand(command+fileName);
        delete ui;

    }


#ifdef G4VIS_USE
    delete visManager;
#endif
    analysisManager->WriteMacros(macFilename);
	analysisManager->CloseFile();

    delete runManager;


    return 0;
}
