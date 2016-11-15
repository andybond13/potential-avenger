//potential-avenger-launcher.c
//Andrew Stershic
//DCML, Duke University
//(c) 2013, 2014

#include <potential-avenger.h>
#include "boost/random.hpp"
#include "time.h"

using namespace std;

void writeDisplacements (const vector<double>& uIn, const double& L) {
    //write displacements to file for testing
    string dispFil = "./raw_disp.out";
    FILE * pFileP;
    pFileP = fopen(dispFil.c_str(), "w");
    fprintf( pFileP, "#raw displacements normalized by length\n");
    for (unsigned j = 0; j < uIn.size(); ++j) fprintf(pFileP, "%0.8g\n",uIn.at(j)/L);
	return;
}

template <typename T> T min(const vector<T>& in) {
    assert(in.size() > 0);
    T min = in[0];
    for (unsigned i = 1; i < in.size(); ++i) {
        if (in[i] < min) min = in[i];
    }
    return min;
}

int main(int argc, const char* argv[]) {
    assert(argc == 19);
 
	//read input from potential-avenger-launcher.py
    double strain_rate = atof(argv[1]);
    double ts_refine = atof(argv[2]);
    double end_t = atof(argv[3]);
    unsigned Nelt = atoi(argv[4]);
    double lc = atof(argv[5]);
    double alpha = atof(argv[6]);
    double CV = atof(argv[7]);
    unsigned startWithLoad = atoi(argv[8]);
    unsigned printVTK = atoi(argv[9]);
    string sm = string(argv[10]);
    int oneAtATime = atoi(argv[11]);
    double minOpenDist = atof(argv[12]);
    unsigned TLSoption = atoi(argv[13]);
    unsigned visualizeCracks = atoi(argv[14]);
    unsigned fullCompression = atoi(argv[15]);
	unsigned elemDeath = atoi(argv[16]);
	unsigned frontExtension = atoi(argv[17]);
	unsigned maxIteration = atoi(argv[18]);
    string path = ".";
 
    PotentialAvenger pa = PotentialAvenger(strain_rate, ts_refine, end_t, Nelt, lc, printVTK, oneAtATime, minOpenDist, alpha, TLSoption, visualizeCracks, fullCompression, sm, elemDeath, frontExtension, maxIteration, path);

	//"serious" hard-coded input
	unsigned Nnod = Nelt + 1;
	double E = 610.0;
	double rho = 3.9e-6;
	double A = 0.2;
	double L = 1.0;
	double Yc = 8.1969e-4;
	vector<double> pg;
	pg.push_back(0.5 - sqrt(1./12.)); pg.push_back(0.5 + sqrt(1./12.));
//	pg.push_back(0.0); pg.push_back(1.0); pg.push_back(0.5);
//	pg.push_back(0.0); pg.push_back(1.0);
	vector<double> wg;
//	wg.push_back(1./6.); wg.push_back(1./6.); wg.push_back(2./3.); 
	wg.push_back(0.5); wg.push_back(0.5);
	vector<double> phiIn = vector<double>(Nnod, 0.0);
	bool vbc = true;

	double h = 1.0/static_cast<double>(Nelt) * L;
	//double alfa = 1.0 - pow(1.0 - L * strain_rate * sqrt(rho / (2*Yc)),2);
	//double qty = sqrt(2.0 * Yc * (1.0 - alfa)/E) * static_cast<double>(startWithLoad);
	vector<double> xIn = vector<double>(Nnod,0);
	vector<double> uIn = vector<double>(Nnod,0);
	vector<double> vIn = vector<double>(Nnod,0);
	vector<double> YcvIn = vector<double>(Nnod, Yc);
/*
YcvIn.at(Nelt/2) = Yc/2.0;
*/
/*
YcvIn.at(Nelt/7) = Yc/2.0;
YcvIn.at(3*Nelt/7) = Yc/2.0;
YcvIn.at(5*Nelt/7) = Yc/2.0;

*/
/*
	double num = 6; //number of weaknesses
	double amp = 0.1; // amplitude of variation
	const double PI = atan2(0, -1);
    for (unsigned i = 0; i < YcvIn.size(); i++){
		double X = static_cast<double>(i) * h + 0.5*h;
		YcvIn[i] = Yc * (1 + amp * cos(X / L * num * (2.0*PI)));
	}
*/

	double qty = 0.999 * sqrt(2.0 * min(YcvIn) / E) * static_cast<double>(startWithLoad);
	double e0 = qty;
	vector<double> eIn = vector<double>(Nelt,e0);

	if (CV != 0.0) {
		//make weibull distribution
		double _lambda = Yc * 2.1586552 * CV;
		double _k = 2.0;
		double _loc = Yc * (1.0 - 1.9130584 * CV);
    	//Create a uniform distrib'd number
   		boost::random::uniform_real_distribution<double> distro(0., 1.);
   		boost::mt19937 rng(static_cast<unsigned> (time(0)));
	    for (unsigned i = 0; i < YcvIn.size(); i++){
        	//Create a Weibull distrib'd number from the uniform distrib'd number
        	YcvIn[i] = pow( -log(1.0 - distro(rng)), 1/_k) * _lambda + _loc;
    	}
	}

	DamageModel dm = DamageModel("Parabolic", lc);

	for (unsigned j = 0; j < Nnod; ++j) {
		xIn[j] = static_cast<double>(j) * h;
		uIn[j] = xIn[j] * e0;
		vIn[j] = static_cast<double>(vbc) * strain_rate * xIn[j];
//		phiIn[j] = max(0.0, 2.0 * h - xIn[j]);
	}	
	vIn[Nnod-1] = strain_rate*L;

	vector<Segment*> segIn; /*
	Segment* seg1 = new Segment(xIn[0], phiIn[0], -1);
	seg1->indices.push_back(0);
	for (unsigned j = 1; j < Nnod; ++j) if (phiIn[j-1] > 0) seg1->indices.push_back(j);
	segIn.push_back(seg1);
	unsigned nucleated = 1;*/
//YcvIn[0] *= 0.1;
	unsigned nucleated = 0;

    pa.run(E, rho, A, L, Yc, pg, wg, phiIn, segIn, nucleated, vbc, eIn, xIn, uIn, vIn, YcvIn, dm);
	
	//write displacements to file for testing
	//writeDisplacements(uIn,L);
}

