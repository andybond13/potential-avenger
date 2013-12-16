//potential-avenger.c
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <potential-avenger.h>
#include <damageModel.h>
#include <iostream>

using namespace std;

int main(int argc, const char* argv[]) {
}

void PotentialAvenger::nucleate(double t, std::vector<double>& x, std::vector<double>& phi, std::vector<double>& xnuc, std::vector<double>& phinuc){
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

void PotentialAvenger::findFragments(std::vector<double>& x, std::vector<double>& phi, std::vector<double>& d, unsigned& nfrags, std::vector<std::pair<unsigned,unsigned> >& fragmentList) {


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
