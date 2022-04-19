
/*
 *  generate_diffphotBG_churazov08.C
 *
 *  Routine to produce diffuse photon background as described in the 
 *  paper from Gruber et al. ApJ 520 (1999) 124  
 *    The normalization constant is increased by 10%, as suggested by Churazov et al. A&A 467 (2007) 529
 *
 *
 *  How to run:  
 *      - generate_diffphotBG_churazov08(run_number,duration_sec, e_min_keV, e_max_keV,SphrRad_cm,CircRad_cm,"a0_output.root" )
 *
 *
 *  Created by Estela Suarez on 26-Jun-2008
 *  Copyright 2009 Universite de Geneve. All rights reserved.
 *
 *
 */

#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom.h"
#include "t0class.cc"
#include "TMath.h"
#include "TNamed.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;


Double_t fun(Double_t *x, Double_t *par)
{
	//Function from Gruber et al. ApJ 520 (1999) 124, 
	// In units [kev(cm-2 s-1 keV-1 sr-1)]
	// in fact is the spectrum multiplied by the energy
	//x[0]: energy 
	//p[0]: normalization factor (1 for Gruber, 1.1 for Churazov)
	Double_t E = x[0]; //energy in keV
	Double_t result(0);
	
	if (E<60.) 
	{
		result = par[0]*(7.877*pow(E,-0.29)*exp(-E/41.13));
	} else if (E>60.)
	{
		result = par[0]*(0.0259*pow((E/60.),-5.5) + 0.504*pow((E/60.),-1.58) + 0.0288*pow((E/60.),-1.05));
	}
	return result;
}


Double_t spectrum(Double_t *x, Double_t *par)
{
	//Function from Gruber et al. ApJ 520 (1999) 124, 
	// Divided by E to get it in normal spectral units [cm-2 s-1 keV-1 sr-1]
	//x[0]: energy 
	//p[0]: normalization factor (1 for Gruber, 1.1 for Churazov)
	Double_t E = x[0]; //energy in keV
	Double_t result(0);
	
	if (E<60.) 
	{
		result = par[0]*(7.877*pow(E,-0.29)*exp(-E/41.13))/E;
	} else if (E>60.)
	{
		result = par[0]*(0.0259*pow((E/60.),-5.5) + 0.504*pow((E/60.),-1.58) + 0.0288*pow((E/60.),-1.05))/E;
	}
	return result;
}




void generate_diffphotBG_churazov08(int run_number=1,
		     double duration=1.,
		     double e_min=5.,
		     double e_max=1000.,
		     double SphrRad=50.,
		     double CircRad=30.,
		     const char *file = "a0_DFBG_chur.root"
		     )
{

  // pol_angle is relative to x axis if possible
  // not possible if photon along x axis, in this case make it along y axis
  // energy band e_min e_max (keV)
  // duration in sec heavyside form
  // all energies in keV
  // radius in cm
  
  gRandom->SetSeed(0);	 //To produce always random output, according to computer clock
  
  TFile *a = new TFile(file, "recreate");
 
  
  //Function from fitting the data of Olivier Grimm thesis (ETHZ 2004) See reference in Produit et al. 2005 (NIMA)
  //TF1 *f1=new TF1("diffPhot",
  	//"pow(10.,0.940059-1.28089*TMath::Log10(x) - 0.2624141 *TMath::Log10(x)*TMath::Log10(x))",
	//e_min,e_max);    
 
 	
	//Function from Gruber et al. ApJ 520 (1999) 124, 
	// but with a 10% higher normailzation factor,
	// as suggested in Churazov et al. A&A 467 (2007) 529
  Double_t cte[1] = {1.1};	//Normalization constant (1.0 for Gruber)--(1.1 for Churazov)
  TF1 *f1 = new TF1("CXB",spectrum,e_min,e_max,1);
  f1->SetParameters(cte);
	
  f1->SetNpx((int)e_max*10);  
  // Calculate number of photons that should be produced 
  double norm=f1->Integral(e_min,e_max);
  double norm10 = f1->Integral(10.,1000.); 
  //double number = 3.86*2*TMath::Pi()*TMath::Pi()*CircRad*CircRad*duration*norm/norm10;
  double number = 2*TMath::Pi()*TMath::Pi()*CircRad*CircRad*duration*norm;
	
  TCanvas *c1= new TCanvas("c1","c1",400,400);	
  f1->Draw("EH");
	
	cout << norm << "   " << norm10 << "  " << 2*TMath::Pi()*TMath::Pi()*CircRad*CircRad*duration << endl;
  // // double number = 2*TMath::Pi()*TMath::Pi()*CircRad*CircRad*duration;
  // double number = 3.48468*2*TMath::Pi()*TMath::Pi()*CircRad*CircRad*duration/f1->Eval(10.);
	
  int n_phot = (long) gRandom->Poisson(number);

  double* x=new double[3];  
  double* d=new double[3];
  char name[10] = "gamma";
  double pol_frac=0.0;  // Diffuse photon background is not polarized
  double pol_angle=0.;
  
  double dt=duration/n_phot;
  double sume=0;
  
  TDatime now;
  t0class* ip=new t0class();  
  ip->fChain->SetTitle("t0");
  ip->time = 0.;
  ip->run= (Int_t) run_number;
 
  TF1 *f2 = new TF1("f2",fun,e_min,e_max,1);
  Double_t rate=f2->Integral(e_min,e_max);
  	cout << "Entering loop... " << endl;
	cout << "Spectrum Integral: " << norm << " [photons.cm-2.s-1.keV-1.sr-1] " 
	<< " Rate: " << rate << " [photons.cm-2.s-1.sr-1] " << endl;
	cout << "Number of Protons: " << n_phot << endl; 
	cout << "(Conversion constant: " << 2*TMath::Pi()*TMath::Pi()*CircRad*CircRad*duration << ")" << endl;



  for (int i=0;i<n_phot;i++)
  {
    ip->event= (Long_t) i;
    ip->time+= (Double_t) gRandom->Exp(dt);
    ip->energy= (Float_t) f1->GetRandom();
    sume+=ip->energy;
    
    //Photons come randomly above POLAR
    double ctheta = gRandom->Uniform(0,1);	//Half-sphere above POLAR -> Solid Angle = 2 Pi
    double theta = acos(ctheta);
    double phi = gRandom->Uniform()*2.*TMath::Pi();

    x[0]=SphrRad*sin(theta)*cos(phi);
    x[1]=SphrRad*sin(theta)*sin(phi);
    x[2]=SphrRad*cos(theta);
    
    d[0]=-1*sin(theta)*cos(phi);
    d[1]=-1*sin(theta)*sin(phi);
    d[2]=-1*cos(theta);
    
    double philoc=gRandom->Uniform()*2*TMath::Pi();
    double r=CircRad*sqrt(gRandom->Uniform());
    double xx=r*cos(philoc);
    double yy=r*sin(philoc);
    ip->xloc[0]= (Float_t) x[0]+cos(theta)*cos(phi)*xx-sin(phi)*yy;
    ip->xloc[1]= (Float_t) x[1]+cos(theta)*sin(phi)*xx+cos(phi)*yy;
    ip->xloc[2]= (Float_t) x[2]-sin(theta)*xx;
    ip->momentum[0]= (Float_t) d[0]*ip->energy;
    ip->momentum[1]= (Float_t) d[1]*ip->energy;
    ip->momentum[2]= (Float_t) d[2]*ip->energy;
    strcpy(ip->nameParticle,name);// photon
    double modul2= d[1]*d[1]+d[2]*d[2];
    double e_perpend[3];
    if (modul2 > 0.){
      e_perpend[0] = 0;
      e_perpend[1]=-(1./sqrt(modul2))*d[2]; 
      e_perpend[2]=(1./sqrt(modul2))*d[1];
    }
    else{
      e_perpend[0]=0;
      e_perpend[1]=1;
      e_perpend[2]=0;
    }
    double e_paralle[3];
    e_paralle[0]=e_perpend[1]*d[2]-e_perpend[2]*d[1];
    e_paralle[1]=-e_perpend[0]*d[2]+e_perpend[2]*d[0];
    e_paralle[2]=e_perpend[0]*d[1]-e_perpend[1]*d[0];
    
    double DegToRad = 3.1415/180.;
    double pangle = pol_angle*DegToRad;
    if (gRandom->Uniform()>pol_frac){
      pangle=gRandom->Uniform()*2*TMath::Pi();
    }
    ip->polarization[0] = (Float_t) cos(pangle)*e_paralle[0] + sin(pangle)*e_perpend[0];
    ip->polarization[1] = (Float_t) cos(pangle)*e_paralle[1] + sin(pangle)*e_perpend[1];
    ip->polarization[2] = (Float_t) cos(pangle)*e_paralle[2] + sin(pangle)*e_perpend[2];
    
 	ip->epr[0]=-9999.;
	ip->epr[1]=-9999.;
	ip->epr[2]=-9999.;
	ip->epr[3]=-9999.;
	ip->epr[4]=-9999.;
	ip->epr[5]=-9999.;
	ip->epr[6]=-9999.;
	ip->epr[7]=-9999.;
	ip->epr[8]=-9999.;
    ip->fChain->Fill();
  }
  
   
  Int_t InEvents =  n_phot;
  Double_t Iduration = duration;
  Char_t IpartName[10];
  strcpy(IpartName,name);
  Char_t IgeneratorName[100] = "generate_diffphotBG_churazov08";
  Int_t InumParam = 7;
  Double_t Iparam[7] = {run_number, duration, e_min, e_max, SphrRad, CircRad};
  Char_t InamesParam[200] = {"run_number, duration, e_min, e_max, SphrRad, CircRad"};


  TTree *tinput = new TTree("tinput","tinput_v0");
   
  
  tinput -> Branch("InEvents", &InEvents, "InEvents/I");   
  tinput -> Branch("Iduration", &Iduration, "Iduration/D");  
  tinput -> Branch("IpartName", IpartName, "IpartName[10]/C"); 
  tinput -> Branch("IgeneratorName", IgeneratorName, "IgeneratorName[100]/C"); 
  tinput -> Branch("InumParam", &InumParam,"InumParam/I");
  tinput -> Branch("InamesParam", &InamesParam,"InamesParam[200]/C");
  tinput -> Branch("Iparam", Iparam,"Iparam[InumParam]/D");
  
  tinput ->Fill();
  a->ls();
  tinput -> Write();  
  tinput -> Print();
    
  ip->fChain->Write();
  ip->fChain->Print();
  
  // To store command used in file
  char com[200];
  sprintf(com,"generate_diffphotBG_churazov08(%i, %g, %g, %g, %g, %g, '%s') ",
  	run_number, duration, e_min, e_max, SphrRad, CircRad, file);
 
  cout << " COMMAND USED  " << com << endl;    

  TNamed *command = new TNamed("command", com);
  command -> Write(); 
  
  
  
  // Print output to check the rate:
  int n_bins = (int) (e_max-e_min);
  TH1D *spec = new TH1D("spec","spec",n_bins,e_min,e_max);
//  TCanvas *c1= new TCanvas("c1","c1",400,400);
//  ip->fChain->Draw("energy>>spec");
	ip->fChain->Project("spec","energy");
//  delete c1;
 
  double sum10 = spec->Integral(10-(int)e_min+1,1000-(int)e_min+1);	//n_bins); 
  double sumTot = spec->Integral(1,n_bins); 
  cout << "------------------------------------------------------------" << endl;
  cout << " Number Photons [10keV, 1000keV] = " << sum10 << endl;
  cout << " Rate = " << 
  	sum10/(TMath::Pi()*CircRad*CircRad*duration*2*TMath::Pi()) <<
	" photons cm-2 s-1 sr-1 " << "   (Theory: 3.86) " << endl;
  cout << "------------------------------------------------------------" << endl;  
  cout << " Total Number Photons [" << e_min << "keV," << e_max << "keV] = " << sumTot << endl;
  cout << " Total Rate = " << 
  	sumTot/(TMath::Pi()*CircRad*CircRad*duration*2*TMath::Pi()) <<
	" photons cm-2 s-1 sr-1 " << endl;
  cout << "------------------------------------------------------------" << endl;  
   
  
  
  //delete f1;
  delete[] x;
  delete[] d;
  //printf("%f\n",sume/ergtokev/duration/CircRad/CircRad/TMath::Pi());
 
  //return ip->fChain;

  //delete a;
  
 
}


