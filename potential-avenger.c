//potential-avenger.c
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <potential-avenger.h>
#include <damageModel.h>
#include <iostream>
#include <math.h>

using namespace std;

int main(int argc, const char* argv[]) {
    assert(argc == 7);

    double strain_rate = atof(argv[1]);
    double ts_refine = atof(argv[2]);
    double end_t = atof(argv[3]);
    unsigned Nelt = atoi(argv[4]);
    double lc = atof(argv[5]);
    unsigned intOrder = atoi(argv[6]);

    PotentialAvenger pa = PotentialAvenger(strain_rate, ts_refine, end_t, Nelt, lc, intOrder);
    pa.run();
}

PotentialAvenger::PotentialAvenger(double in0, double in1, double in2, unsigned in3, double in4, unsigned in5){
    strain_rate = in0;
    ts_refine = in1;
    end_t = in2;
    Nelt = in3;
    lc = in4;
    intOrder = in5;
};

PotentialAvenger::~PotentialAvenger(){};

void PotentialAvenger::run() {

    printRunInfo();
};

void PotentialAvenger::nucleate(const double t, const std::vector<double>& x, std::vector<double>& phi, const std::vector<double>& xnuc, const std::vector<double>& phinuc){
    //t      -time
    //x      -mesh
    //phi    -level-set calculated at this time-step
    //xnuc   -location(s) of localizations to be nucleated
    //phinuc -amount of the level-set to be set at nucleated localizations

    assert(xnuc.size() == phinuc.size());
    assert(xnuc.size() > 0);
    
    for (unsigned j = 0; j < xnuc.size(); ++j) {
        double h = 0;
        unsigned loc = 0;
        double delta = 0;
        
        for (unsigned i = 0; i < x.size()-1; ++i) {
            if ((xnuc[j] >= x[i]) && (xnuc[j] < x[i+1])) {
                loc = i;
                h = x[i+1] - x[i];                
                delta = (xnuc[j] - x[i])/h;
                break;
            }

        }
    
        assert(loc != 0);
    
        phi[loc] = phinuc[j] - delta*h;
        phi[loc+1] = phinuc[j] - (1-delta)*h;
        printf("crack nucleated, t = %f, x = %f",t,xnuc[j]);
    }
};

void PotentialAvenger::findFragments(const std::vector<double>& x, const std::vector<double>& phi, const std::vector<double>& d, unsigned& nfrags, std::vector<std::pair<unsigned,unsigned> >& fragmentList) {


    fragmentList.clear();
    unsigned sbegin = 0;
    unsigned send = 0;

    for (unsigned i = 0; i < x.size() - 1; ++i) {

        if (sbegin == 0) {
            if (i == 1) {
                sbegin = i;
            }
        }
    
        if (send == 0) {
            if (i == x.size()-2) {
                send = i+1;
            }
            if ((i > 0) && (i < x.size()-1)) {
                if (phi[i] >= phi[i-1] && phi[i] >= phi[i+1] && d[i] == 1) {
                    send = i;
                }
            }
        }
    
        if (sbegin*send != 0) {
            pair<unsigned,unsigned> list = std::pair<unsigned,unsigned>(sbegin,send);
            fragmentList.push_back(list);
            sbegin = send + 1;
            send = 0;
        }
    }

    //calculate total number of fragments (removing symmetry simplification)
    nfrags = fragmentList.size()*2;
    if (d[0] < 1.0) {
        nfrags = nfrags - 1;
    }
};

void PotentialAvenger::checkFailureCriteria(const double t, const std::vector<double>& x, std::vector<double>& phi, std::vector<double>& criterion, const std::string elemOrNodal, const std::vector<double>& qty, const bool absOrAsIs, const bool phiPos, const double failvalue){

    //x              -mesh
    //phi            -level-set calculated at this time-step
    //criterion      -criterion to compare against for failure
    //elemOrNodal    -either 'elem' or 'nodal' - is criterion elemental or nodal
    //qty            -the quantity to be compared to the criterion -e.g. s,Y
    //absOrAsIs      -whether to compare (0) qty or (1) abs(qty) to criterion
    //phiPos         -whether(1) or not(0) failure cannot occur depending on if phi>0
    //failvalue      -what to call phi at localization zone if created - e.g. h
    //failure if qty > criterion

    assert(elemOrNodal.compare("elem") || elemOrNodal.compare("nodal"));
    assert(absOrAsIs == false || absOrAsIs == true);
    if (elemOrNodal.compare("nodal")) {
        assert(x.size() == qty.size());
    } else {
        assert(x.size() == qty.size()+1);
    }

    assert((qty.size() == criterion.size()) || (criterion.size() == 1));

    if (criterion.size() == 1) {
        double val = criterion[0];
        criterion.assign(qty.size(),val);
    }
    
    vector<double> xlist;
    
    for (unsigned i = 0; i < qty.size(); ++i) {

        if (phiPos == 0) {//can't fail if phi>0
            if (elemOrNodal.compare("nodal")) {
                if (phi[i] > 0) {
                    continue;
                }
            } else {
                if (phi[i]>0 || phi[i+1]>0) {
                    continue;
                }
            }
        }
    
        double qtyc = qty[i];

        if (absOrAsIs) {
            qtyc = fabs(qtyc);
        }
        if (qtyc > criterion[i]) {
            if (elemOrNodal.compare("nodal")) {
                xlist.push_back(x[i]);
            } else {
                //assume middle of element
                xlist.push_back(0.5*(x[i]+x[i+1]));
            }
        }
    }

    //nucleate list
    if (xlist.size() > 0) {
        vector<double> failvalueList = vector<double>(xlist.size(),failvalue);
        nucleate(t,x,phi,xlist,failvalueList);
    }

};

void analyzeDamage(const vector<double>& x, vector<double>& phi, const double h, vector<vector<unsigned> >& newSegment) {

    //produce:
    //new phi based on distances - maxima

    vector<vector<unsigned> > segment;

    unsigned sum = 0;
    for (unsigned i = 0; i < phi.size(); ++i) {
        if (phi[i] > -1) sum++;
    }
    if (sum == 0) return;

    vector<unsigned> list_max;
    vector<double> value_max;

    for (unsigned i = 0; i < phi.size()-1; ++i) {
        if ((i > 0) && (i < phi.size()-2)) { //local maxima/minima
    
            if ((phi[i] >= phi[i-1]) && (phi[i+1]>=phi[i+2])) {
                double delta = (phi[i+1]-phi[i] + h)/2;
                if (phi[i]+delta < 0) {
                    continue;
                }
                list_max.push_back( x[i] + delta);
                value_max.push_back(phi[i] + delta);
            }
        }
        if ((i == 0) && (phi[0]-phi[1] >= 0)) {
            list_max.push_back(x[i]);
            value_max.push_back(phi[i]);
        }
        if ((i == phi.size()-2) && phi[i+1]>phi[i]) {
            list_max.push_back(x[i+1]);
            value_max.push_back(phi[i+1]);
        }
    
    }

    segment.resize(list_max.size());    
    assert(x.size() >= 1);
    if (value_max.size() == list_max.size() && value_max.size() == 0) return;

    vector<double>phinew = vector<double>(phi.size(),-1);
    for (unsigned i = 0; i < phi.size(); ++i) {
        
        double min = 9999999999;
        for (unsigned k = 0; k < x.size(); ++k) {
            double qty = -value_max[k] + fabs(x[i] - list_max[k]);
            if (qty < min) min = qty;
        }
        phinew[i] = -min;
       
        for (unsigned j = 0; j < list_max.size(); ++j) {
            if (phinew[i] == -(-value_max[j]+fabs(x[i] -list_max[j]))) {
                segment[j].push_back(i);
                break;
            }
        }
        phinew[i] = max(phinew[i],phi[i]);
    }

    //join segments together if same peak
    for (unsigned j = 0; j < list_max.size()-1; ++j) {
        if (fabs(list_max[j]-list_max[j+1]) < h) {
            segment[j].insert(segment[j].end(),segment[j+1].begin(),segment[j+1].end());
            segment[j+1].clear();
//            segment[j] = sort(segment[j]);//TODO turned off
        }
    }


    //new
    //split hat segments into two
    newSegment.clear();
    vector<unsigned> list_maxnew;
    vector<double> value_maxnew;
    for (unsigned i = 0; i < segment.size(); ++i) {
        vector<unsigned> indices = segment[i];
        if (indices.size() == 0) continue; //don't copy empty's
        if (indices.size() == 1) {
            //solo point - copy as is
            list_maxnew.push_back(list_max[i]);
            value_maxnew.push_back(value_max[i]);
            newSegment.push_back(segment[i]);
            continue;
        }
   
        assert(indices.size() > 1);
        if ((phinew[indices[0]] < phinew[indices[1]]) && (phinew[*indices.end()] < phinew[*(indices.end()-1)] )) {
            //hat - duplicate
            list_maxnew.push_back(list_max[i]);
            list_maxnew.push_back(list_max[i]);
            value_maxnew.push_back(value_max[i]);
            value_maxnew.push_back(value_max[i]);

            int iphimax = -1;
            double phimax = -1;
            for (vector<unsigned>::iterator j = indices.begin() ; j != indices.end(); ++j) {
                if (phi[*j] > phimax) {
                    phimax = phi[*j];
                    iphimax = *j;
                }
            }
            assert(iphimax != -1);

            vector<unsigned> vec1;
            vector<unsigned> vec2;
            for (vector<unsigned>::iterator j = indices.begin() ; j != indices.end(); ++j) {
                if (*j <= iphimax) vec1.push_back(*j);
                else vec2.push_back(*j);
            }
            newSegment.push_back(vec1);
            newSegment.push_back(vec2);

        } else if ((phinew[indices[0]] >= phinew[indices[1]]) && (phinew[*indices.end()] < phinew[*(indices.end()-1)] )) {
            //not hat - copy as is
            list_maxnew.push_back(list_max[i]);
            value_maxnew.push_back(value_max[i]);
            newSegment.push_back(segment[i]);
        } else if ((phinew[indices[0]] <= phinew[indices[1]]) && (phinew[*indices.end()] >= phinew[*(indices.end()-1)] )) {
            //not hat - copy as is
            list_maxnew.push_back(list_max[i]);
            value_maxnew.push_back(value_max[i]);
            newSegment.push_back(segment[i]);       
        } else {
            //[list_max;value_max]
            //[indices;x[indices];phinew[indices]]
            //plot(x[indices],phinew[indices])
            assert(1==0);
        }
    }
    
    for (unsigned i = 0; i < newSegment.size(); ++i) {
        vector<unsigned> indices = newSegment[i];
        for (vector<unsigned>::iterator j = indices.begin() ; j != indices.end(); ++j) {
            phinew[*j] = value_maxnew[i]-fabs(x[*j] -list_maxnew[i]);
        }
    }

    //return phinew as phi
    for (unsigned i = 0; i < phi.size(); ++i) {
        phi[i] = phinew[i];
    }
};  

void PotentialAvenger::printRunInfo() {
    cout << endl;
    cout << "Run Info:" << endl;
    cout << "Strain rate = " << strain_rate << endl;
    cout << "# Elements = " << Nelt << endl;
    cout << "End Time (s) = " << end_t << endl;
    cout << "Time-step refinement ratio = " << ts_refine << endl;
    cout << "Damage gradient length (lc) = " << lc << endl;
    cout << "Integration Order = " << intOrder << endl;
    cout << "Beginning run... " << endl;
    cout << endl;

};

