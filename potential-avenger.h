//potential-avenger.h
//Andrew Stershic
//DCML, Duke University
//(c) 2013,2014,2015

#include <objects.h>
#include <damageModel.h>
#include <math.h>

#ifndef POTENTIAL_AVENGER_H
#define POTENTIAL_AVENGER_H

class PotentialAvenger{

public:
//dyndimresttls.m
void run(const double& E, const double& rho, const double& A, const double& L, const double& Yc, const std::vector<double>& pg, const std::vector<double>& wg, const std::vector<double>& phiIn, const std::vector<Segment*> segIn, const unsigned& nucleated, bool& vbc, const std::vector<double>& eIn, const std::vector<double>& xIn, std::vector<double>& uIn, const std::vector<double>& vIn, const std::vector<double>& YcvIn, const DamageModel& dm);

PotentialAvenger(double& in0, double& in1, double& in2, unsigned& in3, double& in4, unsigned& in6, int& in7, double& in8, double& in9, unsigned& in10, unsigned& in11, unsigned& in12, std::string& sm, unsigned& in13, unsigned& in14, unsigned& in15, std::string& path);
~PotentialAvenger();


private:
double strain_rate;
double ts_refine;
double end_t;
unsigned Nelt, Nnod;
double lc;
double alpha;
unsigned startWithLoad;
unsigned printVTK;
std::string sm;
int oneAtATime;
double minOpenDist;
unsigned TLSoption;
unsigned visualizeCracks;
unsigned fullCompression;
unsigned elemDeath;
unsigned frontExtension;
unsigned maxIteration;

//new variables
std::string _path, _FragFile, _EnrgFile, _SThetaFile, _HistoFile, _FDFile;
double E, A, rho, c, L, h, dt, Yc, sigc, ec;
unsigned _numFrag, _Nt, _DtPrint;
double strain_energy, dissip_energy, dissip_energy_TLS, dissip_energy_local, kinetic_energy, max_energy, ext_energy, tot_energy;
double _fMean, _fMed, _fMax, _fMin, _fStDev, _fRange, _fSkew, _fExKurtosis, _fAltMin;
std::vector<double> x, t, xe, d, u, v, a, s, e, phiL, phiNL, Y, Ycv, energy, m, d_1, u_1, Ystat, ustat, phiNL_1, Ybar;
std::vector<unsigned> nfrags;	
DamageModel dm;
std::vector<unsigned> inTLS;
std::vector<unsigned> inTLSnode;
std::vector<double> d_max,d_max_alt;
std::vector<double> gradPhiL;
std::vector<double> gradPhiNL;
std::vector<double> gradPhiNLelem;
std::vector<double> altFragLength;
std::vector<double> phidot;
unsigned nucleated;
std::vector<std::vector<double> > d_quad,d_quad_wt,d_quad_phi;
double EPS,Fboundary;

void printRunInfo();

void calculateEnergies(const unsigned& i, const std::vector<double>& pg, const std::vector<double>& wg);
unsigned calculateStressesL(const std::vector<double>& pg, const std::vector<double>& wg, const unsigned& j);
void calculateStressesNL(const std::vector<double>& pg, const std::vector<double>& wg, std::vector<Segment*>& newSegment);
void setPeak(const std::vector<double>& phi, std::vector<Segment*>& segments, const unsigned index);

//nucleate.m
void nucleate(double t, const std::vector<double>& xnuc, const std::vector<double>& phinuc, std::vector<Segment*>& newSegment, const unsigned& elemOrNodal);

//findFragments.m
std::vector<double> findFragments(unsigned& nfrags, const std::vector<Segment*>& segments);
std::vector<double> fragmentLength(const std::vector<Segment*>& segments);
std::vector<double> localFragmentLength();

//checkFailureCriteria.m
unsigned checkFailureCriteria(unsigned ts, std::vector<double>& criterion, const unsigned& elemOrNodal, const std::vector<double>& qty, bool absOrAsIs, bool phiPos, double failvalue, std::vector<Segment*>& newSegment, const std::vector<double>& pg, const std::vector<double>& wg);

//analyzeDamage.m
void analyzeDamage(std::vector<double>& phi, const double h, std::vector<Segment*>& newSegment);

//calculate level-set gradient
void calculateLevelSetGradientL(const std::vector<double>& d, std::vector<double>& gradPhi);
void calculateLevelSetGradientNL(const std::vector<double>& d, std::vector<double>& gradPhi);

//get list of which elements are in a TLS zone
void checkInTLS(const std::vector<Segment*>& segments, std::vector<unsigned>& elem, std::vector<unsigned>& node);

//check enforcement of constraints (Ycbar < Yc), (|gradPhi|<=1)
void checkConstraints(const std::vector<double>& gradientPhiL, const std::vector<double>& gradientPhiNL, const std::vector<Segment*>& segments);

//update level set for nodes in TLS
void updateLevelSetL(std::vector<Segment*>& segments, const std::vector<double>& pg, const std::vector<double>& wg, const unsigned& j);
void updateLevelSetNL(std::vector<Segment*>& segments, const std::vector<double>& pg, const std::vector<double>& wg);

double H(const unsigned, const double);
double dH(const unsigned, const double);
double d2H(const unsigned, const double);

//calculate Ybar
unsigned calculateYbar(const std::vector<double>& pg, const std::vector<double>& wg, double& Ycavg, double& YbarmYc, double& residu_Y, double& tangent_Y, double& phimin, double& phimax, double& phiminY, double& phimaxY, unsigned& nbiter, const unsigned sbegin, const unsigned send, Segment* segment, const double& segZero);

void calculateDmaxAlt(const std::vector<double>& pg, const std::vector<double>& wg);

void plotEnergies ();
void plotForceDisp ();
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
void printForceDisp () const;
void printFrags (const std::vector<double>& fragLength, const unsigned nSegs);
void printSTheta ();
void printHisto (const std::vector<double>& fragLength);
void printClean() const;
void fragmentStats(const std::vector<double>& fragLength);
void killSegments(std::vector<Segment*>& seg);
//void killFragments(std::vector<Fragment*>& frag);
void setPeakAll(const std::vector<double>& phiin, std::vector<Segment*>& segments);
double calculateZero(Segment* segment, const std::vector<double>& phiIN);
double calculateTotal(Segment* segment, const std::vector<double>& phiIN);
};



#endif

