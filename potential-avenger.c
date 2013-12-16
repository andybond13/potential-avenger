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






