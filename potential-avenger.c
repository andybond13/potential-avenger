//potential-avenger.c
//Andrew Stershic
//DCML, Duke University
//(c) 2013

#include <potential-avenger.h>
#include <matrix.h>
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

PotentialAvenger::PotentialAvenger(double& in0, double& in1, double& in2, unsigned& in3, double& in4, unsigned& in5){
    strain_rate = in0;
    ts_refine = in1;
    end_t = in2;
    Nelt = in3;
    lc = in4;
    intOrder = in5;
};

PotentialAvenger::~PotentialAvenger(){};

unsigned min(const vector<unsigned> in) {
    assert(in.size() > 0);
    unsigned min = in[0];
    for (unsigned i = 1; i < in.size(); ++i) {
        if (in[i] < min) min = in[1];
    }
    return min;
}

unsigned max(const vector<unsigned> in) {
    assert(in.size() > 0);
    unsigned max = in[0];
    for (unsigned i = 1; i < in.size(); ++i) {
        if (in[i] > max) max = in[1];
    }
    return max;
}

double dotProduct(const vector<double> a, const vector<double> b) {
    double sum = 0;
    assert(a.size() == b.size() );
    for (unsigned i = 0; i < a.size(); ++i) sum+= a[i]*b[i];
    return sum;
}

unsigned median(const vector<unsigned> in) {
    return in[(in.size()-1)/2];
}

double sum(const vector<double> in) {
    return dotProduct(in,vector<double>(in.size(),1));
}


void PotentialAvenger::run() {

    printRunInfo();

    DamageModel dm = DamageModel("Cubic",lc);

    unsigned Ntim = Nelt*ts_refine*end_t;
    unsigned Nnod = Nelt+1;
    double E = 1; //(beton)
    double rho = 1; //(beton)
    double A = 1; // barre de 10cm sur 10cm
    double c = sqrt(E/rho);
    double L = 1;
    double h = 1/static_cast<double>(Nelt); //
    double cfl = 1./ts_refine;
    double dt = cfl * h/c;//
    double Yc = E/10;
    double ec = sqrt(2 * Yc / E);
    double sigc = E * ec;

    // two gauss point on the element
    vector<double> pg(2);
    pg[0] = (1-sqrt(3)/3)/2;
    pg[1] = (1+sqrt(3)/3)/2;

    Matrix d = Matrix(Ntim, Nelt);
    Matrix s = Matrix(Ntim, Nelt);
    Matrix e = Matrix(Ntim, Nelt);
    Matrix energy = Matrix(Ntim, Nelt);
    Matrix Y = Matrix(Ntim, Nelt);
    vector<double> strain_energy(Ntim,0);
    vector<double> kinetic_energy(Ntim,0);
    vector<double> dissip_energy(Ntim,0);
    vector<double> ext_energy(Ntim,0);
    vector<double> tot_energy(Ntim,0);
    vector<Fragment> fragment_list;

    //Matrix bpos = Matrix(Ntim, Nelt);
    //Matrix grad = Matrix(Ntim, Nelt);
    Matrix u = Matrix(Ntim, Nnod);
    Matrix ustat = Matrix(Ntim, Nnod);
    Matrix Ystat = Matrix(Ntim, Nnod);
    Matrix v = Matrix(Ntim, Nnod);
    Matrix a = Matrix(Ntim, Nnod);
    Matrix phi = Matrix(Ntim, Nnod);
    Matrix YmYc = Matrix(Ntim, Nelt);

    Matrix phidot = Matrix(Ntim,0);
    //vector<double> ddotbar(Ntim,0);
    vector<double> dissip(Ntim,0);
    vector<double> nbiter(Ntim,0);
    vector<unsigned> nfrags(Ntim,0);
    vector<Segment> segments;

    vector<double> m(Nnod,rho*h*A);
    m[0] = m[0]/2; m[Nnod-1] = m[Nnod-1]/2;

    //initialization
    vector<double> x(Nnod);
    for (unsigned j = 0; j < Nnod; ++j) x[j] = j*h;

    vector<double> xe(Nnod);
    for (unsigned j = 0; j < Nelt; ++j) xe[j] = 0.5*(x[j] + x[j+1]);

    vector<double> t(Ntim);
    for (unsigned j = 0; j < Ntim; ++j) t[j] = j*dt;

    // Initially the bar is loaded and all elements are at Yc
    // A tls is placed on the first element, obviously it satifies
    //   Ybar = Yc
    // This extra damage will create new stress that are less than in the next element
    // and ill thus give an unloading wave.

    //constant strain rate applied
    double csr = strain_rate;
    bool vbc = true;
    for (unsigned j = 0; j < Nnod; ++j) v(0,j) = vbc*csr*x[j];

    for (unsigned j = 0; j < Nnod; ++j) {
        u(0,j) = x[j] * ec * L * 0.999*(1-vbc);
        ustat(0,j) = x[j] * ec * L * 0.999*(1-vbc);
        phi(0,j) = (2*h-x[j])*(1-vbc)-vbc;
    }

    for (unsigned j = 0; j < Nelt; ++j) {
        e(0,j) = (u(0,j+1)-u(0,j))/h;
        vector<double> dloc(2,0.0);
        if (phi(0,j) > 0  && phi(0,j+1) > 0) {
            for (unsigned k = 0; k < 2; ++k) {
                double philoc = pg[k]*phi(0,j)+ (1-pg[k])*phi(0,j+1);
                dloc[k] = dm.dval(philoc);
                s(0,j) += 0.5 * (1-dloc[k]) * E * e(0,j);
            }
        } else if  (phi(0,j) <= 0 && phi(0,j+1) <= 0) {
            s(0,j) = E * e(0,j);
        } else if  (phi(0,j) > 0 && phi(0,j+1) <= 0) {
            double delta = fabs(phi(0,j))/(fabs(phi(0,j))+fabs(phi(0,j+1)));
            double sloc = 0;
            for (unsigned k = 0; k < 2; ++k) {
                double philoc = pg[k]*phi(0,j);
                dloc[k] = dm.dval(philoc);
                sloc += 0.5 * (1-dloc[k]) * E * e(0,j);
            }
            s(1,j) = delta *  sloc +  (1-delta) * E * e(0,j);
        } else if (phi(0,j) <= 0 && phi(0,j+1) > 0) {
            double delta = fabs(phi(0,j+1))/(fabs(phi(0,j))+fabs(phi(0,j+1)));
            double sloc = 0;
            for (unsigned k = 0; k < 2; ++k) {
                double philoc = pg[k]*phi(0,j+1);
                dloc[k] = dm.dval(philoc);
                sloc += 0.5 * (1-dloc[k]) * E * e(0,j);
            }
            s(0,j) = delta *  sloc +  (1-delta) * E * e(0,j);
        }
        d(0,j) = 0.5*(dloc[0]+dloc[1]);
    }

    //acceleration
    a(0,0) = 0;
    a(0,Nnod-1) = 0;
    for (unsigned j = 1; j < Nnod-1; ++j)  a(0,j) = A*(s(0,j) - s(0,j-1)) /m[j];

    nbiter[0] = 0;
    analyzeDamage(x,phi.row(0),h,segments);
    unsigned len = 0;
    for (unsigned l = 0; l < segments.size(); ++l) {
        if (segments[l].size() == 0) continue;
        len++;
    }
    phidot.row(0).resize(len);

    //time-integration loop
    for (unsigned i = 1; i < Ntim; ++i){
//cout << "t = " << t[i] << endl;

        //prediction
        for (unsigned j = 0; j < Nnod; ++j) {
            v(i,j)= v(i-1,j) + 0.5*dt*a(i-1,j);
            u(i,j)= u(i-1,j) + dt*v(i-1,j) + 0.5*dt*dt*a(i-1,j);
        }

        //def computation and Y update.
        for (unsigned j = 0; j < Nelt; ++j) {
            e(i,j) = (u(i,j+1)-u(i,j))/(h);
            //b=0.5*E*e(i,j)*e(i,j)-Yc;
        }

        // moving the localization front
        // we compute a = integral (Yn+1 - Yc) d' in the current non-local zone
        // then we compute b = (Yn+1-Yc) d' on the front
        // the shift in level set if the ratio of the two.
        
        for (unsigned j = 0; j < Nnod; ++j) {
            phi(i,j) = phi(i-1,j);
        }

        for (unsigned l = 0; l < segments.size(); ++l) {
            if (segments[l].size()==0) continue;//segments.erase(segments.begin()+l);
            unsigned sbegin = segments[l].begin();
            unsigned send = segments[l].end();
        
            //skip if all negative
            bool allNeg = true;
            for (unsigned k = sbegin; k <= send; ++k) {
                if (phi(i,k) != -1) {
                    allNeg = false;
                    break;
                }
            }
            if (allNeg) continue;
        
            double err_crit = 1e15;
            double dphi = 0;
            nbiter[i] = 0;
            double residu = 0;

            while (err_crit > 1.e-6) {
                nbiter[i]++;
                double residu_Y = 0; double tangent_Y = 0;
                unsigned loop_residu = 0;
                unsigned loop_tangent = 0;
                for (unsigned j = sbegin-1; j <= min(send,Nelt-2); ++j) {
                    if (j < 0) continue;
                    if (phi(i,j) > 0 && phi(i,j+1) > 0) {
                        for (unsigned k = 0; k < 2; ++k) {
                            double philoc = pg[k]*phi(i,j) + (1-pg[k])*phi(i,j+1);
                            residu_Y += h * 0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dm.dp(philoc);
                            tangent_Y += h * 0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dm.dpp(philoc);
                        }
                        loop_residu++;
                    } else if  (phi(i,j) > 0 && phi(i,j+1) <= 0) {
                        double delta = h * fabs(phi(i,j))/(fabs(phi(i,j))+fabs(phi(i,j+1))); //phi>0 portion
                        for (unsigned k = 0; k < 2; ++k) {
                            double philoc = pg[k]*phi(i,j);
                            residu_Y += delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dm.dp(philoc);
                            tangent_Y += delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dm.dpp(philoc);
                        }
                        loop_residu++;
                        if (delta < h) tangent_Y += (0.5 * E * e(i,j) * e(i,j) - Yc)* dm.dp(0.);
                        else tangent_Y += (0.5 * E * e(i,j+1) * e(i,j+1) - Yc)  * dm.dp(0.);
                        
                        loop_tangent = loop_tangent + 1;
                    } else if  (phi(i,j) <= 0 && phi(i,j+1) > 0) {
                        double delta = h * fabs(phi(i,j+1))/(fabs(phi(i,j))+fabs(phi(i,j+1))); //phi>0 portion
                        for (unsigned k = 0; k < 2; ++k) {
                            double philoc = pg[k]*phi(i,j+1);
                            residu_Y += delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dm.dp(philoc);
                            tangent_Y += delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dm.dpp(philoc);
                        }
                        loop_residu++;
                        if (delta < h) tangent_Y += (0.5 * E * e(i,j) * e(i,j) - Yc)* dm.dp(0.); //%todo-doublecheck this  
                        else tangent_Y += (0.5 * E * e(i,j-1) * e(i,j-1) - Yc)  * dm.dp(0.);   //%todo-doublecheck this
                        loop_tangent++;
                    
                    }
                }
            
                // law Ybar = Yc
                if (1==1) {
                    bool flag = true; //0 is exterior (damage centered on edge of domain); 1 is interior
                    if (l == 1 && segments[l].begin() == 0) flag = false;
                    if (l == segments.size()-1 && segments[l].end() == Nnod-1) flag = false;

                    
                    double phimax = phi(i,sbegin);
                    unsigned iphimax = sbegin;
                    for (unsigned k = sbegin; k <= send; ++k) {
                        if (phi(i,k) > phimax) {
                            phimax = phi(i,k);
                            iphimax = k;
                        }
                    }
                    phimax = segments[l].phipeak;

                    double phimaxY;
                    if (iphimax == Nnod-1) phimaxY = 0.5*E*pow(e(i,Nelt-1),2); //1/2*s(i,Nelt)*e(i,Nelt);
                    else phimaxY = 0.5*E*pow(e(i,iphimax),2);// 1/2*s(i,iphimax)*e(i,iphimax);
                    
                    double YbarmYc = residu_Y/(dm.dval(phimax));
                    double oldresidu = residu;
                    residu = YbarmYc/Yc;
                    err_crit = fabs(residu-oldresidu);
                    double tangent = (tangent_Y + flag*(phimaxY-Yc)*dm.dp(phimax)/2)/(Yc*dm.dval(phimax)) - ((1+flag/2)*dm.dp(phimax)/pow(dm.dval(phimax),2)) * (YbarmYc/Yc);
                    if (fabs(tangent) <= 1.e-10) {
                        err_crit = 0.; dphi = 0.;
                    } else {
                        dphi = - residu/tangent;
                    }
                }

                if (nbiter[i] > 50) {
                    dphi = 0;
                    //assert(1==0);
                }
            
                if (isnan(dphi)) {
                    dphi = 0;
                    assert(1==0);                
                }

                for (unsigned j = sbegin; j <=send; ++j) {
                    if (i > 6 && intOrder >= 6) {
                        vector<double> w;
                        w.push_back(60./147); w.push_back(360./147); w.push_back(-450./147); w.push_back(400./147);
                        w.push_back(-225./147); w.push_back(72./147); w.push_back(-10./147); 
                        vector<double> phihist;
                        phihist.push_back(dphi); phihist.push_back(phi(i-1,j)); phihist.push_back(phi(i-2,j));
                        phihist.push_back(phi(i-3,j)); phihist.push_back(phi(i-4,j)); phihist.push_back(phi(i-5,j));
                        phihist.push_back(phi(i-6,j));
                        phi(i,j) = dotProduct(phihist,w);
                    } else if (i > 5 && intOrder >= 5) {
                        vector<double> w;
                        w.push_back(60./137); w.push_back(300./137); w.push_back(-300./137); w.push_back(200./137);
                        w.push_back(-75./137); w.push_back(12./137);
                        vector<double> phihist;
                        phihist.push_back(dphi); phihist.push_back(phi(i-1,j)); phihist.push_back(phi(i-2,j));
                        phihist.push_back(phi(i-3,j)); phihist.push_back(phi(i-4,j)); phihist.push_back(phi(i-5,j));
                        phi(i,j) = dotProduct(phihist,w);
                    } else if (i > 4 && intOrder >= 4) {
                        vector<double> w;
                        w.push_back(12./25); w.push_back(48./25); w.push_back(-36./25); w.push_back(16./25);
                        w.push_back(-3./25);
                        vector<double> phihist;
                        phihist.push_back(dphi); phihist.push_back(phi(i-1,j)); phihist.push_back(phi(i-2,j));
                        phihist.push_back(phi(i-3,j)); phihist.push_back(phi(i-4,j));
                        phi(i,j) = dotProduct(phihist,w);
                    } else if (i > 3 && intOrder >= 3) {
                        vector<double> w;
                        w.push_back(6./11); w.push_back(18./11); w.push_back(-9./11); w.push_back(2./11);
                        vector<double> phihist;
                        phihist.push_back(dphi); phihist.push_back(phi(i-1,j)); phihist.push_back(phi(i-2,j));
                        phihist.push_back(phi(i-3,j));
                        phi(i,j) = dotProduct(phihist,w);
                    } else if (i > 2 && intOrder >= 2) {
                        vector<double> w;
                        w.push_back(2./3); w.push_back(4./3); w.push_back(-1./3);
                        vector<double> phihist;
                        phihist.push_back(dphi); phihist.push_back(phi(i-1,j)); phihist.push_back(phi(i-2,j));
                        phi(i,j) = dotProduct(phihist,w);
                    } else {
                        phi(i,j) = phi(i-1,j) + dphi;
                    }
                    phi(i,j) = max(phi(i,j),phi(i-1,j)); //constraint: dphi >= 0
                    //enforcing limit of level-set motion
                    //phi(i,j) = min(phi(i-1,j)+h,phi(i,j));
                }
            } //while
        
            err_crit = 0.0;
        } //for segments

        //check for nucleation
        vector<double> Yin;
        for (unsigned l = 0; l < Nelt; ++l)  Yin.push_back(0.5*E*e(i,l)*e(i,l));
        vector<double> YcVec(1,Yc);
        string elemOrNodal="elem";
        checkFailureCriteria(t[i],x,phi.row(i),YcVec,elemOrNodal,Yin,false,false,2*h,segments);

        //enforce phi constraints - update segments
        analyzeDamage(x,phi.row(i),h,segments);

        int index = -1;
        for (unsigned l = 0; l < segments.size(); ++l) {
            if (segments[l].size() == 0) continue;
            index++;
            //median(segments[l]);
            unsigned smid = median(segments[l].indices);

            len = 0;
            for (unsigned j = 0; j < segments.size(); ++j) {
                if (segments[j].size() == 0) continue;
                len++;
            }
            phidot.resizeRow(i,len);
            phidot(i,index) = (phi(i,smid) - phi(i-1,smid))/dt;
            if (phidot(i,index)*dt > h*1.01 ) {
                printf("level-set front advancing more than one element per time-step: t=%f, segment %u , dphi/h = %f \n",t[i],l,phidot(i,index)*dt/h);
            }
        }

        //updating the stress
        for (unsigned j = 0; j < Nelt; ++j) {
            vector<double> dlocg(2,0);
            if (phi(i,j) > 0 && phi(i,j+1) > 0) {
                for (unsigned k = 0; k < 2; ++k) {
                    double philoc = pg[k]*phi(i,j)+ (1-pg[k])*phi(i,j+1);
                    dlocg[k] = dm.dval(philoc);
                    s(i,j) = s(i,j) + 0.5 * (1-dlocg[k]) * E * e(i,j);
                }
            } else if  (phi(i,j) <= 0 && phi(i,j+1) <= 0) {
                s(i,j) = E * e(i,j);
            } else if  (phi(i,j) > 0 && phi(i,j+1) <= 0) {
                double delta = fabs(phi(i,j))/(fabs(phi(i,j))+fabs(phi(i,j+1)));
                double sloc = 0;
                for (unsigned k = 0; k < 2; ++k) {
                    double philoc = pg[k]*phi(i,j);
                    dlocg[k] = dm.dval(philoc);
                    sloc = sloc + 0.5 * (1-dlocg[k]) * E * e(i,j);
                }
                s(i,j) = delta * sloc +  (1-delta) * E * e(i,j);
            } else if  (phi(i,j) <= 0 && phi(i,j+1) > 0) {
                double delta = fabs(phi(i,j+1))/(fabs(phi(i,j))+fabs(phi(i,j+1)));
                double sloc = 0;
                for (unsigned k = 0; k < 2; ++k) {
                    double philoc = pg[k]*phi(i,j+1);
                    dlocg[k] = dm.dval(philoc);
                    sloc = sloc + 0.5 * (1-dlocg[k]) * E * e(i,j);
                }
                s(i,j) = delta * sloc +  (1-delta) * E * e(i,j);
            }
            d(i,j) = 0.5*(dlocg[0]+dlocg[1]);
        }

        //acceleration
        if (phi(i,0) <= lc) a(i,0) = 0;
        else a(i,0) =  A*s(i,0)/m[0];
        a(i,Nnod-1) = 0; //note: this is a constraint which was hidden. if free, it would be -s(i,Nnod)/m(j);
        for (unsigned j = 1; j < Nnod - 1; ++j) a(i,j) = A*(s(i,j) - s(i,j-1)) /m[j];

        //correction
        for (unsigned j = 0; j < Nnod; ++j) v(i,j)= v(i,j) + 0.5*dt*a(i,j);

        //record number of fragments and quantities per fragment
        fragment_list.clear();
        vector<double> phirow = phi.rowC(i);
        findFragments(dm, segments,phirow,nfrags[i],fragment_list);
    }

    for (unsigned i = 0; i < Ntim; ++i) {
        for (unsigned j = 0; j < Nelt; ++j) {
            Y(i,j) = 0.5*E*pow(e(i,j),2);
            YmYc(i,j) = Y(i,j)/Yc - 1;
            if (i > 0) {
                dissip[i] += h * Y(i,j) * (d(i,j)-d(i-1,j));
                kinetic_energy[i] += 0.5 * h * rho * 0.5 *
                                    ( pow(u(i,j) - u(i-1,j),2) + pow(u(i,j+1) - u(i-1,j+1),2) )/pow(dt,2);
                //ustat(i,j+1) = ustat(i,j) + h*s(0,Nelt)/(E*(1-d(i,j)));
            }
            energy(i,j)=h*Y(i,j)*(1-d(i,j));
        }
        //ustat(i,:) *= u(1,Nnod)/ustat(i,Nnod);
        for (unsigned j = 0; j < Nelt; ++j) Ystat(i,j) = 0.5*E*pow((ustat(i,j+1)-ustat(i,j))/h,2);
        if (i > 0) ext_energy[i] = (a(i,Nnod-1)*m[Nnod-1]+s(i,Nnod-2))*v(i,Nnod-1)*dt + ext_energy[i-1];
        if (i > 0) dissip_energy[i] = dissip_energy[i-1] + dissip[i];
        strain_energy[i] = sum (energy.rowC(i));
        tot_energy[i] = strain_energy[i] + kinetic_energy[i] + dissip_energy[i] - ext_energy[i];
    }

    double minfrag = L*2;
    double sumfrag = 0;
    for (unsigned i = 0; i < nfrags[Ntim-1]/2; ++i) {
        
        double fragLen = fragment_list[i].length()*h;

        if (((nfrags[Ntim-1] % 2) == 1) && (fragment_list[i].begin() == 0)) fragLen *= 2;
    
        if (fragLen < minfrag) minfrag = fragLen;

        sumfrag += fragLen;
        cout << "fragment " << i << " length = " << fragLen << endl;
    }

    if (nfrags[Ntim-1] > 0) assert(sumfrag == L);
    printf("Final number of fragments: %i \nMinimum fragment length: %f \nFinal dissipated energy: %f \n",nfrags[Ntim-1],minfrag,dissip_energy[Ntim-1]);

};

void PotentialAvenger::nucleate(const double t, const std::vector<double>& x, std::vector<double*> phi, const std::vector<double>& xnuc, const std::vector<double>& phinuc, std::vector<Segment>& newSegment){
    //t      -time
    //x      -mesh
    //phi    -level-set calculated at this time-step
    //xnuc   -location(s) of localizations to be nucleated
    //phinuc -amount of the level-set to be set at nucleated localizations

    assert(xnuc.size() == phinuc.size());
    assert(xnuc.size() > 0);

    for (unsigned j = 0; j < xnuc.size(); ++j) {
        double h = 0;
        unsigned loc = -1;
        double delta = 0;
        
        for (unsigned i = 0; i < x.size()-1; ++i) {
            if ((xnuc[j] >= x[i]) && (xnuc[j] < x[i+1])) {
                loc = i;
                h = x[i+1] - x[i];                
                delta = (xnuc[j] - x[i])/h;
                break;
            }

        }
    
        assert(loc != -1);
    
        *phi[loc] = phinuc[j] - delta*h;
        *phi[loc+1] = phinuc[j] - (1-delta)*h;
        printf("crack nucleated, t = %f, x = %f \n",t,xnuc[j]);
        
        //create two new segments
        unsigned Nnod = x.size();
        if (xnuc[j] > x[0]) {
            Segment seg1 = Segment(xnuc[j],phinuc[j]-delta*h,-1,Nnod);
            seg1.indices.push_back(loc);
            newSegment.push_back(seg1);
        }
        if (xnuc[j] < x[x.size()-1]) {
            Segment seg2 = Segment(xnuc[j],phinuc[j]-(1-delta)*h,1,Nnod);
            seg2.indices.push_back(loc+1);
            newSegment.push_back(seg2);
        }

    }
};

void PotentialAvenger::findFragments(DamageModel& dm, std::vector<Segment>& newSegment, const std::vector<double>& phi, unsigned& nfrags, std::vector<Fragment>& fragmentList) {

    fragmentList.clear();

    sort(newSegment.begin(),newSegment.end());

    Fragment f = Fragment();
    for (unsigned i = 0; i < newSegment.size(); ++i) {

        if (newSegment[i].size() == 0) continue;
        f.fragSegs.push_back(&newSegment[i]);

        double phimax = newSegment[i].phipeak;
        double slope = newSegment[i].slope;

        if (slope == 1 && dm.dval(phimax) == 1.0) {
            fragmentList.push_back(f);
            f.fragSegs.clear();
        } else {}
            //cout << slope << " , " << dm.dval(phimax) << endl;
    }
    
    //calculate total number of fragments (removing symmetry simplification)
    if (fragmentList.size() == 0) {nfrags = 0; return;};
    nfrags = fragmentList.size()*2;
    if (dm.dval(phi[0]) < 1.0) {
        nfrags = nfrags - 1;
    }
};

void PotentialAvenger::checkFailureCriteria(const double t, const std::vector<double>& x, std::vector<double*> phi, std::vector<double>& criterion, const std::string elemOrNodal, const std::vector<double>& qty, const bool absOrAsIs, const bool phiPos, const double failvalue, std::vector<Segment>& newSegment){

    //x              -mesh
    //phi            -level-set calculated at this time-step
    //criterion      -criterion to compare against for failure
    //elemOrNodal    -either 'elem' or 'nodal' - is criterion elemental or nodal
    //qty            -the quantity to be compared to the criterion -e.g. s,Y
    //absOrAsIs      -whether to compare (0) qty or (1) abs(qty) to criterion
    //phiPos         -whether(1) or not(0) failure cannot occur depending on if phi>0
    //failvalue      -what to call phi at localization zone if created - e.g. h
    //failure if qty > criterion

    assert(elemOrNodal.compare("elem") == 0 || elemOrNodal.compare("nodal") == 0);
    assert(absOrAsIs == false || absOrAsIs == true);
    if (elemOrNodal.compare("nodal") == 0) {
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
            if (elemOrNodal.compare("nodal") == 0) {
                if (*phi[i] > 0) continue;
            } else {
                if (*phi[i]>0 || *phi[i+1]>0) continue;
            }
        }

        //can't fail if already nucleated
        double h = x[i+1]-x[i];
        for (unsigned j = 0; j < newSegment.size(); ++j) {
            cout << x[i] << " vs. " << newSegment[j].xpeak << endl;
            if (fabs(newSegment[j].xpeak-x[i]) < h) break; 
        }


        double qtyc = qty[i];

        if (absOrAsIs) {
            qtyc = fabs(qtyc);
        }

        if (qtyc > criterion[i]) {

            if (elemOrNodal.compare("nodal") == 0) {
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
        nucleate(t,x,phi,xlist,failvalueList, newSegment);
    }
};

void PotentialAvenger::analyzeDamage(const vector<double>& x, vector<double*> phi, const double h, vector<Segment>& newSegment) {

    //produce:
    //new phi based on distances - maxima

    //if all negative one, skip
    unsigned sum = 0;
    for (unsigned i = 0; i < phi.size(); ++i) {
        if (*phi[i] > -1) sum++;
    }
    if (sum == 0) return;

    vector<double> list_max;
    vector<double> value_max;
    sort(newSegment.begin(), newSegment.end());
/*
    //manually look for segments
    for (unsigned i = 0; i < phi.size()-1; ++i) {
        if ((i > 0) && (i < phi.size()-2)) { //local maxima/minima
    
            if ((*phi[i] >= *phi[i-1]) && (*phi[i+1]>=*phi[i+2])) {
                double delta = (*phi[i+1]-*phi[i] + h)/2;
                if (*phi[i]+delta < 0) {
                    continue;
                }
                list_max.push_back( x[i] + delta);
                value_max.push_back(*phi[i] + delta);
            }
        }
        if ((i == 0) && (*phi[0]-*phi[1] >= 0)) {
            list_max.push_back(x[i]);
            value_max.push_back(*phi[i]);
        }
        if ((i == phi.size()-2) && *phi[i+1]>*phi[i]) {
            list_max.push_back(x[i+1]);
            value_max.push_back(*phi[i+1]);
        }
    
    }
*/
//    segment.resize(list_max.size());    
    assert(x.size() >= 1);
    //if (value_max.size() == list_max.size() && value_max.size() == 0) return;

    //simply use pre-defined segments
    if (newSegment.size() == 0) return;
    for (unsigned i = 0; i < newSegment.size(); ++i) {
        //assert(newSegment[i].size() > 0);
        list_max.push_back(newSegment[i].xpeak);
        value_max.push_back(newSegment[i].phipeak);
        newSegment[i].indices.clear();
    }

    //re-define segments.indices and phi
    vector<double> phinew = vector<double>(phi.size(),-1);
    for (unsigned i = 0; i < phi.size(); ++i) {
        
        unsigned segphimin;
        double min = 9999999999;
        for (unsigned k = 0; k < value_max.size(); ++k) {
            double qty = -value_max[k] + fabs(x[i] - list_max[k]);
            if (qty < min) {
                min = qty;
                segphimin = k;
            }
        }
        phinew[i] = -min;
        
        newSegment[segphimin].indices.push_back(i);
        /*for (unsigned j = 0; j < list_max.size(); ++j) {
            if (phinew[i] == -(-value_max[j]+fabs(x[i] -list_max[j]))) {
                newSegment[j].indices.push_back(i);
                break;
            }
        }*/
        phinew[i] = max(phinew[i],*phi[i]);
    }
/*
    //join segments together if same peak
    for (unsigned j = 0; j < list_max.size()-1; ++j) {
        if (fabs(list_max[j]-list_max[j+1]) < h) {
            segment[j].indices.insert(segment[j].indices.end(),segment[j+1].indices.begin(),segment[j+1].indices.end());
            segment.erase(segment.begin()+j+1);
        }
    }
*/
/*
    //new
    //split hat segments into two
    newSegment.clear();
    vector<unsigned> list_maxnew;
    vector<double> value_maxnew;
    for (unsigned i = 0; i < segment.size(); ++i) {
        Segment seg = segment[i];
        if (seg.size() == 0) continue; //don't copy empty's
        if (seg.size() == 1) {
            //solo point - copy as is
            list_maxnew.push_back(list_max[i]);
            value_maxnew.push_back(value_max[i]);
            newSegment.push_back(segment[i]);
            continue;
        }
   
        assert(seg.size() > 1);
        unsigned iend = seg.end();
        unsigned iend1 = seg.penult();
        if ((phinew[seg.begin()] < phinew[seg.second()]) && (phinew[iend] < phinew[iend1] )) {
            //hat - duplicate
            int iphimax = -1;
            double phimax = -1;
            for (vector<unsigned>::iterator j = seg.indices.begin() ; j != seg.indices.end(); ++j) {
                if (*phi[*j] > phimax) {
                    phimax = *phi[*j];
                    iphimax = *j;
                }
            }
            assert(iphimax != -1);

            Segment vec1;
            Segment vec2;
            for (vector<unsigned>::iterator j = seg.indices.begin() ; j != seg.indices.end(); ++j) {
                if (*j <= iphimax) vec1.indices.push_back(*j);
                else vec2.indices.push_back(*j);
            }
            vec1.xpeak = list_max[i];
            vec1.phipeak = value_max[i];
            vec1.slope = -1;
            vec2.xpeak = list_max[i];
            vec2.phipeak = value_max[i];
            vec2.slope = 1;

            newSegment.push_back(vec1);
            newSegment.push_back(vec2);

        } else if ((phinew[seg.begin()] >= phinew[seg.second()]) && (phinew[iend] <= phinew[iend1] )) { //TODO changed right to <=
            //not hat - copy as is
            list_maxnew.push_back(list_max[i]);
            value_maxnew.push_back(value_max[i]);
            newSegment.push_back(segment[i]);
        } else if ((phinew[seg.begin()] <= phinew[seg.second()]) && (phinew[iend] >= phinew[iend1] )) {
            //not hat - copy as is
            list_maxnew.push_back(list_max[i]);
            value_maxnew.push_back(value_max[i]);
            newSegment.push_back(segment[i]);       
        } else {
for (unsigned j = 0; j < seg.size(); ++j) cout << seg.indices[j] <<" , " <<  phinew[seg.indices[j]] << endl;
cout << (phinew[seg.begin()] <= phinew[seg.second()]);
cout << phinew[seg.begin()]-phinew[seg.second()] << endl;
cout << (phinew[iend1] <= phinew[iend]) << endl;
            //[list_max;value_max]
            //[indices;x[indices];phinew[indices]]
            //plot(x[indices],phinew[indices])
            assert(1==0);
        }
    }
    */
    unsigned tot_indices = 0;
    for (unsigned i = 0; i < newSegment.size(); ++i) {
        Segment seg = newSegment[i];
        for (vector<unsigned>::iterator j = seg.indices.begin() ; j != seg.indices.end(); ++j) {
            double value_maxnew = newSegment[i].phipeak;
            double list_maxnew = newSegment[i].xpeak;
            phinew[*j] = value_maxnew-fabs(x[*j] -list_maxnew);
            assert(phinew[*j] <= newSegment[i].phipeak);
            tot_indices++;
        }
    }
    assert(tot_indices == x.size());

    //return phinew as phi
    for (unsigned i = 0; i < phi.size(); ++i) {
        *phi[i] = phinew[i];
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

