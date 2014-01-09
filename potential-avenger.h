//potential-avenger.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <stdio.h>
#include <vector>
#include <objects.h>
#include <damageModel.h>
#include <iostream>

#ifndef POTENTIAL_AVENGER_H
#define POTENTIAL_AVENGER_H

class PotentialAvenger{

public:
//dyndimresttls.m
void run();

PotentialAvenger(double& in0, double& in1, double& in2, unsigned& in3, double& in4, unsigned& in5, unsigned& in6, std::string& path);
~PotentialAvenger();


private:
double strain_rate;
double ts_refine;
double end_t;
unsigned Nelt;
unsigned Nnod;
double lc;
unsigned intOrder;
unsigned printVTK;

//new variables
std::string _path, _FragFile, _EnrgFile, _SThetaFile, _HistoFile;
double E, A, rho, c, L, h, dt, Yc, sigc, ec;
unsigned _numFrag, _Nt, _DtPrint;
double strain_energy, dissip_energy, kinetic_energy, max_energy, ext_energy, tot_energy;
double _fMean, _fMed, _fMax, _fMin, _fStDev, _fRange, _fSkew, _fExKurtosis;
std::vector<double> x, t, xe, d, u, v, a, s, e, phi, Y, YmYc, energy;
std::vector<unsigned> nfrags;	
DamageModel dm;
std::vector<Fragment> fragment_list;

void printRunInfo();

//nucleate.m
void nucleate(double t, const std::vector<double>& x, std::vector<double>& phi, const std::vector<double>& xnuc, const std::vector<double>& phinuc, std::vector<Segment>& newSegment);

//findFragments.m
void findFragments(DamageModel& dm, std::vector<Segment>& newSegment, const std::vector<double>& phi, unsigned& nfrags, std::vector<Fragment>& fragmentList);

//checkFailureCriteria.m
void checkFailureCriteria(double t, const std::vector<double>& x, std::vector<double>& phi, std::vector<double>& criterion, std::string elemOrNodal, const std::vector<double>& qty, bool absOrAsIs, bool phiPos, double failvalue, std::vector<Segment>& newSegment);

//analyzeDamage.m
void analyzeDamage(const std::vector<double>& x, std::vector<double>& phi, const double h, std::vector<Segment>& newSegment);

void plotEnergies ();
void plotFrags ();
void plotHisto ();
void plotSTheta ();

void display ( const unsigned timStepNumElas, const unsigned timStepNumFrac );

void printVtk ( const unsigned timStepNum ) const;
void printHeader ( const std::string& vtkFile ) const;
void printMesh ( const std::string& vtkFile ) const;
void printPointData ( const std::string& vtkFile ) const;
void printCellData ( const std::string& vtkFile ) const;
void printGlobalInfo () const;
void printFrags ();
void printSTheta ();
void printClean() const;

};



#endif

