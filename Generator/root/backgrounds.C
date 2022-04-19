//#include <EventFormat.hh>

double cxbspc(double *x, double *par)
{
	// cxb spectrum
	return pow(10, par[0]+par[1]*TMath::Log10(x[0])+par[2]*TMath::Log10(x[0])*TMath::Log10(x[0]));

}
double elespc(double *x, double *par)
{
	return par[0]*(x[0]**(par[1]));
}
double posspc(double *x, double *par)
{
	return par[0]*(x[0]**(par[1]));
}




double cxbspec99(double *x, double *par)
{
	//Function from Gruber et al. ApJ 520 (1999) 124, 
	// Divided by E to get it in normal spectral units [cm-2 s-1 keV-1 sr-1]
	//x[0]: energy 
	//p[0]: normalization factor (1 for Gruber, 1.1 for Churazov)
	double E = x[0]; //energy in keV
	double result(0);


	if (E<60.) 
	{
		result = par[0]*(7.877*pow(E,-0.29)*exp(-E/41.13))/E;
	} else if (E>60.)
	{
		result = par[0]*(0.0259*pow((E/60.),-5.5) + 0.504*pow((E/60.),-1.58) + 0.0288*pow((E/60.),-1.05))/E;
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






void backgrounds(
		double duration=200,
		double emin=5,
		double emax=1e6,
		double SphrRad=50.,
		double CircRad=30.,
		string fileName="bkg_200s.root"

		)
{
	const double Pi=3.1415926;


	// diffuse cosmic x rays
//	TF1 *cxbb=new TF1("cxbb",cxbspc,emin,emax,3);
//	cxbb->SetNpx(500);
//	cxbb->SetParameters(0.940059, -1.28089, -0.262414);



	//diffuse cosmic x ray, Gruber
	TF1 *cxbb99=new TF1("cxbb99",cxbspec99, emin,emax,1);
	cxbb99->SetNpx(500);
	cxbb99->SetParameter(0,1.1);
	double cxbb99i=cxbb99->Integral(emin,emax);



	double eminele=40;
	emax=10000;
	// electron
	TF1 *eleb=new TF1 ("eleb",elespc,eminele,emax,2);
	eleb->SetNpx(500);
	eleb->SetParameters(9113., -2.30);
	double elebi=eleb->Integral(eminele,emax);

	//positron 
	TF1 *posb=new TF1 ("posb",posspc,eminele,emax,2);
	posb->SetNpx(500);
	posb->SetParameters(11.14., -1.70);
	double posbi=posb->Integral(eminele,emax);

	//caculate number of event
	double allIntegral=cxbb99i+posbi+elebi;
//	double allIntegral=elebi;

	double pos_frac=posbi/allIntegral;
	double cxb_frac=(cxbb99i+posbi)/allIntegral;
	double ele_frac=1;
	cout<<"background(cxb, positron, electron): "<< cxbb99i<<" "<<posbi<<" "<<elebi<<endl;

	cout<<"Sum of the backgrounds (positron, electron and CXB ) :"<<allIntegral<<endl;

	double number=2*Pi*Pi*CircRad*CircRad*duration*allIntegral;
//	double number=2*Pi*Pi*SphrRad*SphrRad*duration*allIntegral;
	cout<<"conversion factor: "<<2*Pi*Pi*CircRad*CircRad*duration<<endl;

	cout<<"Number of Photons: " << number<<endl; 
	int n_phot=(int)gRandom->Poisson(number);
	cout<<"Number of Photons to be generated : "<<n_phot<<endl;



	double dt=duration/n_phot;
	cout<<"delta time is :"<<dt<<endl;

	TDatime now;


	//structure for event output root file

	double position[3];
	double polarization[3];
	double direction[3];
	char pname[10];
	double energy;
	double time=now;

	TFile f(fileName.c_str(),"recreate");


	t2sim *ts=new t2sim();

        
	double theta,phi;

	for(int i=0;i<n_phot;i++)
	{
		//sampling 
		time+=gRandom->Exp(dt);

		double rnd=gRandom->Uniform();
		if(rnd<pos_frac)
		{//positron
			strcpy(pname,"e+");
			energy=posb->GetRandom();

		}
		else if(rnd<cxb_frac)
		{
			//cxb
			strcpy(pname,"gamma");
			energy=cxbb99->GetRandom();

		}
		else 
		{
			strcpy(pname,"e-");
			//electron
			energy=eleb->GetRandom();
		}

		theta=acos(gRandom->Uniform());
		phi=2*Pi*gRandom->Uniform();
		position[0]=SphrRad*cos(phi)*sin(theta);
		position[1]=SphrRad*sin(phi)*sin(theta);
		position[2]=SphrRad*cos(theta);

		//direction

		direction[0]=-cos(phi)*sin(theta);
		direction[1]=-sin(phi)*sin(theta);
		direction[2]=-cos(theta);

		double philoc=gRandom->Uniform()*2*TMath::Pi();
		double r=CircRad*sqrt(gRandom->Uniform());
		double xx=r*cos(philoc);
		double yy=r*sin(philoc);
		position[0]= (Float_t) position[0]+cos(theta)*cos(phi)*xx-sin(phi)*yy;
		position[1]= (Float_t) position[1]+cos(theta)*sin(phi)*xx+cos(phi)*yy;
		position[2]= (Float_t) position[2]-sin(theta)*xx;
		polarization[0]=0.;
		polarization[1]=0.;
		polarization[2]=0.;
		if(i%10000==0)cout<<"Event ID:"<<i<<endl;


		strcpy(ts->particle_name,pname);

		ts->time=time;
		ts->energy=energy;
		ts->direction[0]=direction[0];
		ts->direction[1]=direction[1];
		ts->direction[2]=direction[2];
		ts->position[0]=position[0];
		ts->position[1]=position[1];
		ts->position[2]=position[2];
		ts->polarization[0]=polarization[0];
		ts->polarization[1]=polarization[1];
		ts->polarization[2]=polarization[2];
		ts->fChain->Fill();

	}

	ts->fChain->Write();
        
	cout<<"done!"<<endl;
	cout<<"particle information has been written to :"<<fileName<<endl;


	delete cxbb;
	delete eleb;
	delete posb;
	delete t;

}
