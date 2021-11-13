//Simulazione di un decadimento di un mesone B in un pione + un kaone. Il mesone è nel sistema di riferimento del laboratorio, mentre i restanti nel centro di massa.
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

#include "TRandom.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"

int main(int argc, char **argv){

	//dati del problema
	double m_B = 5.279; //GeV
	double m_pi = 0.140; //GeV
	double m_k = 0.500; //GeV
	double p_B = 0.300; //GeV (asse x)
	double p_cm = 2.614; //impulso dei prodotti nel centro di massa in GeV
	double angle, mi, mi_meas; //massa invariante
	int N=10000;
	vector<double> resol;
	resol.push_back(0.01);
	resol.push_back(0.03);
	resol.push_back(0.05);
	resol.push_back(0.10);

	//apertura file
	TString rootfname("./output.root");
	TFile rfile(rootfname, "RECREATE");
	if(!rfile.IsOpen()){
		cout <<"Errore nell'apertura del file." << endl;
		exit(1);
	}

	//dichiarazione del vettore
	


	TLorentzVector p4_B, p4_pi, p4_k, p_t, p4_pi_O, p4_k_O; //qi del mesone
	p4_B.SetPxPyPzE(p_B,0,0,sqrt(p_B*p_B+m_B*m_B)); //sistema di riferimento del lab
	TVector3 p1, p2;
	//istogramma

	//invariant mass
	TH1F h1("h1", "True Mass", 100, 5.2, 5.4);
	TH1F hm1("hm1", "Measured Mass sigma=1%", 100, 4.8, 5.8);
	TH1F hm2("hm2", "Measured Mass sigma=3%", 100, 4.8, 5.8);
	TH1F hm3("hm3", "Measured Mass sigma=5%", 100, 4.8, 5.8);
	TH1F hm4("hm4", "Measured Mass sigma=10%", 100, 4.8, 5.8);
	vector<TH1F> h1_meas;
	h1_meas.push_back(hm1);
	h1_meas.push_back(hm2);
	h1_meas.push_back(hm3);
	h1_meas.push_back(hm4);
	TH1F h2("h2", "Lab angles", 100, 2.9, 3.3);
	TH1F pi("pi", "Module", 100, 0, 6);

	//TTree
	TTree* tree = new TTree("datatree", "tree containing our data");
	int nDau=2;
	double nmass[]={m_pi, m_k};
	double p[nDau];
	double Theta[nDau];
	double Phi[nDau];
	tree->Branch("p_B", &p_B,  "p_B/D");
	tree->Branch("nDau", &nDau, "nDau/I");
	tree->Branch("nmass", nmass,  "nmass[nDau]/D");
	tree->Branch("p", p, "p[nDau]/D");
	tree->Branch("theta", Theta,  "theta[nDau]/D");
	tree->Branch("phi", Phi, "phi[nDau]/D");


	//chiamata di TRandom3

	TRandom* gen = new TRandom();
	gen->SetSeed(0);
	double x, y, z;

	//dati per la simulazione del detector
	double p_k_O, p_pi_O, p_pi_meas, p_k_meas;
	

	for(int i=0; i < N; ++i){
	  
	  //punto casuale sulla sfera			
	  gen->Sphere(x,y,z,p_cm);

	  p4_pi.SetPxPyPzE(x,y,z,sqrt(p_cm*p_cm+m_pi*m_pi)); //4-impulso pione centro di massa
	  
	  p4_k.SetPxPyPzE(-x,-y,-z,sqrt(p_cm*p_cm+m_k*m_k)); //4-impulso kaone centro di massa

	  
	  //massa invariante
	  p_t=p4_pi+p4_k; //4-impulso totale
	  mi=sqrt(p_t.Dot(p_t)); //Massa invariante 
	  h1.Fill(mi);
			
	  //boost
	  p4_pi.Boost(p4_B.BoostVector()); //si passa a sistema di riferimento del lab
	  p4_k.Boost(p4_B.BoostVector());
	  
	  //calcolo angoli
	  p1 = p4_pi.Vect(); //3-vettori
	  p2 = p4_k.Vect();

	  
			
	  angle = acos(p1.Dot(p2)/sqrt(p1.Dot(p1)*p2.Dot(p2)));
	  h2.Fill(angle);

	  //calcolo moduli (la loro distribuzione è uniforme)
	  p_pi_O = sqrt(p1.Dot(p1));
	  p_k_O = sqrt(p2.Dot(p2));

	  //si salvano gli eventi dei branch del TTree (SR lab)
	  p[0]=p_pi_O;
	  p[1]=p_k_O;
	  Theta[0]=p1.Theta();
	  Theta[1]=p2.Theta();
	  Phi[0]=p1.Phi();
	  Phi[1]=p2.Phi();
	  tree->Fill();

	  //misura di p_i e p_k a resol=0.03
	  p_pi_meas = gen->Gaus(p_pi_O,resol[1]*p_pi_O);
	  p_k_meas = gen->Gaus(p_k_O,resol[1]*p_k_O);
	  pi.Fill(p_pi_meas);

	  //misura della massa invariante a diverse resol
	  for(int j=0; j<4; j++){
	    mi_meas=gen->Gaus(mi, resol[j]*mi);
	    h1_meas[j].Fill(mi_meas);
	  }
	}
	
	h1.Write();
	for(int j=0; j<4; j++){
	  h1_meas[j].Write();
	}
	h2.Write();
	pi.Write();

	tree->Write();
	tree->Print();

	//plotting mi
	TCanvas canv1("canv1", "Canvas 1", 1280, 1024);


	//plotting  invariant masses (measured and not)	
	canv1.Divide(2,3);
	canv1.cd(1);
	h1.Draw();
	for(int j=0; j<4; j++){
	  canv1.cd(j+2);
	  h1_meas[j].Draw();
	}
	canv1.SaveAs("invariant-mass.png");

	
	//plotting angle
	canv1.Clear();
	h2.GetXaxis()->SetTitle("Distribution of angles [rad]");
	h2.Draw();
	canv1.SaveAs("true-angles.png");

	//plotting p_pi_meas
	canv1.Clear();
	pi.GetXaxis()->SetTitle("Distribution of pi [GeV]");
	pi.Draw();
	canv1.SaveAs("PI-measured-momentum.png");

	//chiusure varie
	rfile.Close();
	delete gen;


	
	return 0;
}
