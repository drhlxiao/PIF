
double amsd[]={
0.22,		0.31,		154e-2,
0.31,		0.44,		156e-2,
0.44,		0.62,		143e-2,
0.62,		0.85,		120e-2,
0.85,		1.15,		966e-3,
1.15,		1.54,		738e-3,
1.54,		2.02,		533e-3,
2.02,		2.62,		372e-3,
2.62,		3.38,		247e-3,
3.38,		4.31,		161e-3,
4.31,		5.45,		101e-3,
5.45,		6.86,		630e-4,
6.86,		8.6,		378e-4,
8.6,		10.7,		226e-4,
10.7,		13.3,		135e-4,
13.3,		16.5,		786e-5,
16.5,		20.5,		449e-5,
20.5,		25.3,		266e-5,
25.3,		31.2,		148e-5,
31.2,		38.4,		856e-6,
38.4,		47.3,		496e-6,
47.3,		58.2,		284e-6,
58.2,		71.5,		154e-6,
71.5,		87.8,		86.2e-6,
87.8,		108,		49.4e-6,
108,		132,		29.0e-6,
132,		162,		16.4e-6,
162,		199,		9.3e-6
};
// Physics Letter B, 2000a AMS proton vol 490 p27-33
double x[28];

double y[28];

double spc(double *xxx,double *par)
{
	double xx=xxx[0];
//	return (x+550)(1+1.8552/x)**(-1.278);

	int po=-1;
	for(int i=0;i<27;i++)
	{
		if(x[i]<xx&&xx<x[i+1])
		{
			po=i;
			break;
		}
	}
	if(po==-1)return 0;
	else
	{

	double x0=x[i];
	double y0=y[i];
	double x1=x[i+1];
	double y1=y[i+1];
	return y0+(xx-x0)*(y1-y0)/(x1-x0);
	}

}


void amsproton(
		double duration=1000,
		double emin=5,
		double emax=1e4,
		double SphrRad=15,
		double CircRad=15,
		string fileName="proton_ams.root"

		)
{
	const double Pi=3.1415926;


	for(int i=0;i<28;i++)
	{
		x[i]=(amsd[3*i]+amsd[3*i+1])/2;
		y[i]=amsd[3*i+2]/(amsd[3*i+1]-amsd[3*i]);
		//GeV to MeV
	}
//	f.Draw();

	// diffuse cosmic x rays


	//diffuse cosmic x ray, Gruber



	double eminele=40;
	// electron
//	TF1 *pspc=new TF1 ("pspc","1.1e4*x**(-2.73)",20,1000);

	TF1 *pspc=new TF1("spc",spc,0.2,200,0);
	pspc->Write();
	double allIntegral=pspc->Integral(0.2,200);
	cout<<"All Integral is: "<<allIntegral<<endl;

	//positron 

	//caculate number of event
//	double allIntegral=cxbb99i+posbi+elebi;


	double number=2*Pi*Pi*CircRad*CircRad*duration*allIntegral/1e4;
	cout<<"conversion factor: "<<2*Pi*Pi*CircRad*CircRad*duration<<endl;

	cout<<"Number of protons: " << number<<endl; 
	int n_phot=(int)gRandom->Poisson(number);
	cout<<"Number of Photons to be generated : "<<n_phot<<endl;



	double dt=duration/n_phot;
	cout<<"delta time is :"<<dt<<endl;



	//structure for event output root file

	double position[3];
	double polarization[3];
	double direction[3];
	char pname[10];
	double energy;
	double time=0;

	TFile f(fileName.c_str(),"recreate");
	//TTree *t=new TTree("eventtree","eventtree");
//	t->Branch("position",position,"position[3]/D");
//	t->Branch("direction",direction,"direction[3]/D");
//	t->Branch("polarization",polarization,"polarization[3]/D");
//	t->Branch("pname",pname,"pname[10]/C");
//	t->Branch("energy",&energy,"energy/D");
//	t->Branch("time",&time,"time/D");


	t2sim *ts=new t2sim();

	double theta,phi;
	//n_phot=100000;

	for(int i=0;i<n_phot;i++)
	{
		//sampling 
		time+=gRandom->Exp(dt);

		strcpy(pname,"proton");
			//electron
		energy=1e9*pspc->GetRandom();
		//gev to keV

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
	//	t->Fill();
	
		strcpy(ts->particle_name,"proton");

	//	direction[0]=1;
	//	direction[1]=0;
	//	direction[2]=0;
	//	position[0]=-30;
	//	position[1]=0.2;
	//	position[2]=0.2;

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
		if(i%10000==0)cout<<"Event ID:"<<i<<endl;

	}



	ts->fChain->Write();

}
