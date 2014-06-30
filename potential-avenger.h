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

PotentialAvenger(double& in0, double& in1, double& in2, unsigned& in3, double& in4, unsigned& in5, unsigned& in6, int& in7, double& in8, double& in9, unsigned& in10, unsigned& in11, unsigned& in12, std::string& path);
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
bool oneAtATime;
double minOpenDist;
double alpha;
unsigned localOnly;
unsigned visualizeCracks;
unsigned fullCompression;

//new variables
std::string _path, _FragFile, _EnrgFile, _SThetaFile, _HistoFile;
double E, A, rho, c, L, h, dt, Yc, sigc, ec;
unsigned _numFrag, _Nt, _DtPrint;
double strain_energy, dissip_energy, dissip_energy_TLS, dissip_energy_local, kinetic_energy, max_energy, ext_energy, tot_energy;
double _fMean, _fMed, _fMax, _fMin, _fStDev, _fRange, _fSkew, _fExKurtosis;
std::vector<double> x, t, xe, d, u, v, a, s, e, phi, Y, Ycv, YmYc, energy, m, d_1, u_1, Ystat, ustat, phi_1, phi_2, phi_3, phi_4, phi_5, phi_6;
std::vector<unsigned> nfrags;	
std::vector<unsigned> d_type;	
DamageModel dm;
std::vector<unsigned> inTLS;
std::vector<unsigned> inTLSnode;
std::vector<double> d_max,d_max_alt;
std::vector<double> gradPhi;
unsigned nucleated;
std::vector<std::vector<double> > d_quad,d_quad_wt;

void printRunInfo();

void calculateEnergies(const unsigned& i, const std::vector<double>& pg, const std::vector<double>& wg);
void calculateStresses(const std::vector<double>& pg, const std::vector<double>& wg, std::vector<Segment*>& newSegment);

//nucleate.m
void nucleate(double t, const std::vector<double>& x, std::vector<double>& phi, const std::vector<double>& xnuc, const std::vector<double>& phinuc, std::vector<Segment*>& newSegment, const std::string& elemOrNodal);

//findFragments.m
std::vector<double> findFragments(DamageModel& dm, std::vector<Segment*>& newSegment, const std::vector<double>& phi, unsigned& nfrags, std::vector<Fragment*>& fragmentList);
std::vector<double> fragmentLength(std::vector<Fragment*>& fragmentList);

//checkFailureCriteria.m
void checkFailureCriteria(double t, const std::vector<double>& x, std::vector<double>& phi, std::vector<double>& criterion, std::string elemOrNodal, const std::vector<double>& qty, bool absOrAsIs, bool phiPos, double failvalue, std::vector<Segment*>& newSegment);

//analyzeDamage.m
void analyzeDamage(const std::vector<double>& x, std::vector<double>& phi, const double h, std::vector<Segment*>& newSegment);

//calculate level-set gradient
void calculateLevelSetGradient(const std::vector<double>& d, std::vector<double>& gradPhi);

//get list of which elements are in a TLS zone
void checkInTLS(const std::vector<Segment*>& segments, std::vector<unsigned>& elem, std::vector<unsigned>& node);

//update level set for nodes in TLS
void updateLevelSet(const unsigned& i, std::vector<unsigned>& nbiter, std::vector<Segment*>& segments, const std::vector<double>& pg, const std::vector<double>& wg);

double H(const unsigned, const double);
double dH(const unsigned, const double) const;

void calculateDmaxAlt(const std::vector<double>& pg, const std::vector<double>& wg);

void plotEnergies ();
void plotFrags ();
void plotHisto ();
void plotSTheta ();

void display ( const unsigned timStepNumElas, const unsigned timStepNumFrac );

void printVtk ( const unsigned timStepNum ) const;
void printHeader ( const std::string& vtkFile ) const;
unsigned printMesh ( const std::string& vtkFile ) const;
void printPointData ( const std::string& vtkFile ) const;
void printCellData ( const std::string& vtkFile, const unsigned& Ncell ) const;
void printGlobalInfo () const;
void printFrags (const std::vector<double>& fragLength);
void printSTheta ();
void printHisto (const std::vector<double>& fragLength);
void printClean() const;
void fragmentStats(const std::vector<double>& fragLength);
void killSegments(std::vector<Segment*>& seg);
void killFragments(std::vector<Fragment*>& frag);

};



#endif

