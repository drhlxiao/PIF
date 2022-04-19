
#include <iostream>
#include <string>    
#include <fstream>    
using namespace std; 
#include <TApplication.h> 
#include <TSystem.h> 
#include <TGraph.h>
#include <TROOT.h> 
#include <TFile.h>  
#include <TTree.h>  
#include <TStyle.h>  
#include <TCanvas.h>  
#include <TF1.h>  
#include <TH1F.h> 
 #include <TH2F.h> 
int main(int argc, char *argv[])
{
	 
	double E0[]={230, 200, 72.8,70.8};
	double th=0;
	double energy=0;
	TString fname;
	ofstream outfile("result.csv"); 

	 if(!outfile.good())
	 { 
		  cout<<"can not create file..."<<endl; 
		   return 0;
	 } 


	double maxL[4]={56,45,7.5,7.5};
	  
	for(int i=0;i<4;i++)

		
	for( th=0.5; th<maxL[i];th+=0.5)
	{
		energy=E0[i];
		fname=Form("Simout_%.1fMeV_%.1fmm.root",energy,th);
		cout<<"fname:"<<fname<<endl;
		TString pngfile=Form("Simout_%.1fMeV_%.1fmm.svg",energy,th);
		TFile ff(fname);
		if(!ff.IsOpen()&&!ff.IsZombie())continue;
		TH1F *h=(TH1F*)ff.Get("hspec_CH0");
		TH1F *h1=(TH1F*)ff.Get("hspec_CH1");
		TH1F *h2=(TH1F*)ff.Get("hspec_CH2");
		TH1F *h3=(TH1F*)ff.Get("hspec_CH3");
		TH1F *hedep=(TH1F*)ff.Get("hedep");

		double mean=hedep->GetMean();


		if(!h)continue;
		if(h->GetEntries()<100)continue;

		TCanvas c1("c1","c1",1200,1000);
		c1.Divide(2,2);
		double max=		h->GetMaximum();
		double rms=		h->GetRMS();
		double mean=		h->GetMean();
		double xmin=mean-4*rms;
		double xmax=mean+4*rms;
		c1.cd(1);
		gStyle->SetOptFit(1111);
		h->Fit("gaus","R","",xmin,xmax);
		h->GetXaxis()->SetRangeUser(xmin,xmax+rms);
		h->SetLineColor(kBlack);
		h->SetTitle("protons");
		h->Draw();

		c1.cd(2);
		h1->SetLineColor(kBlack);
		h1->SetTitle("neutrons");
		h1->Draw();
		c1.cd(3);
		h2->SetLineColor(kBlack);
		h2->SetTitle("gamma-rays");
		h2->Draw();
		c1.cd(4);
		h3->SetTitle("e-");
		h3->SetLineColor(kBlack);
		h3->Draw();
		c1.Update();
		c1.Print(pngfile);




	}
	

	return 0;
}
