#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>

#include "define.h"
#include "vector.h"

using namespace vemc2::mymath;

///Konstanten unserer Kammer:
bdt V       = 0.06252;   //Volumen der Kammer in m^3
bdt dV_pump = -.00134; //Pumpleistung in m^3/sekunde
bdt n_leck  = 0.0000345; //Leck
bdt P_const = 10;         //10 W Leistung durch Boden

///vector functions:
vector<bdt> v_add(vector<bdt> a, vector<bdt> b){
    vector<bdt> retVec;
    for (int i=0; i<a.size(); i++){
        retVec.push_back(a.at(i)+b.at(i));
    }
    return retVec;
}

vector<bdt> v_add(vector<bdt> a, vector<bdt> b, vector<bdt> c){
    vector<bdt> retVec;
    for (int i=0; i<a.size(); i++){
        retVec.push_back(a.at(i)+b.at(i)+c.at(i));
    }
    return retVec;
}

vector<bdt> v_mul(vector<bdt> a, vector<bdt> b){
    vector<bdt> retVec;
    for (int i=0; i<a.size(); i++){
        retVec.push_back(a.at(i)*b.at(i));
    }
    return retVec;
}

vector<bdt> v_mul(vector<bdt> a, bdt b){
    vector<bdt> retVec;
    for (int i=0; i<a.size(); i++){
        retVec.push_back(a.at(i)*b);
    }
    return retVec;
}

#define a_1 ((c_wN2*R_m*R_m)/(d_H*d_H))
#define a_2 ((c_wN2*R_m)/(d_H*d_H))

vector<bdt> f_strich(vector<bdt> input){
    vector<bdt> retV;

    bdt p    = input.at(0);
    bdt n_fl = input.at(1);
    bdt T    = input.at(2);

    bdt _p     = ( R_m*T*((P_const)/(d_H)+n_leck) + p*dV_pump ) / (V + T*T*T*n_fl*a_1/p);
    bdt _n_fl  = -(P_const/d_H) + n_fl * a_2 * ((T*T)/p)*_p;
    bdt _T     = (R_m * T*T*_p)/(d_H * p);

    /*bdt _p     = (p*dV_pump+R_m*T*n_leck)/((c_wN2*R_m*R_m*T*T*T)/(d_H*d_E*p*n_fl)+V); ///p'
    bdt _n_fl  = ((c_wN2*R_m)/(d_H*d_E*p))*n_fl*T*T*_p;                  ///n'_fl
    bdt _T     = T*T*(R_m/(d_H*p))*_p;                                   ///T'*/

    retV.push_back(_p);
    retV.push_back(_n_fl);
    retV.push_back(_T);
    return retV;
}

vector<bdt> RK4_schritt(vector<bdt> data, bdt dt){

    vector<bdt> k1 = f_strich(data);
    vector<bdt> k2 = f_strich(v_add(data, v_mul(k1, dt/2)));
    vector<bdt> k3 = f_strich(v_add(data, v_mul(k2, dt/2)));
    vector<bdt> k4 = f_strich(v_add(data, k3));

    data = v_add(data, v_mul(v_add(k1,v_mul(v_add(k2, k3), 2),k4), dt/6));

    ///Korrigiere Temperatur:
    data[2] = 1 / ( (((-R_m /d_H)*log(data[0]/p_t))+ (1./T_t) ) );

    return data;
}

vector<bdt> run_sim(bdt p_0, bdt n_fl, bdt dt, bdt stop, vector<bdt> p_mess, bool print, std::ostream& print_s, bdt offset){
    vector<bdt> data;
    data.push_back(p_0);
    data.push_back(n_fl);
    data.push_back(1 /( ( ((-R_m /d_H)*log(data[0]/p_t)) + (1./T_t) ) ));

    std::cout << "start values: p=" << data.at(0) << ", n_fl=" << data.at(1) << ", T=" << data.at(2) << ", P_const=" << P_const << std::endl;

    int it_pos=0;

    bdt deltaInsg=0;
    bdt t=0;
    for (; t<=stop; t+=dt){
        data = RK4_schritt(data, dt);
        deltaInsg += data.at(0)-p_mess.at(it_pos);
            if (print) print_s<<t+offset<<" "<<data.at(0)/100<<" "<<data.at(1)<<" "<<data.at(2)<<std::endl;
        if (it_pos<(p_mess.size()-1))
            it_pos++;
        else
            {}//std::cout << "p_mess out of range" << std::endl;

        if (data.at(1)<=0)
            break;
    }
    deltaInsg /= (bdt) it_pos;
    data.push_back(deltaInsg);
    data.push_back(t);
    return data;
}

bdt optimize_P(bdt stepsize, bdt p_0, bdt stop, vector<bdt> vec_p, bdt n_fl, bdt P_min, bdt P_max, bdt delta_p_epsilon){

    bdt delta_p;
    bdt p_end = vec_p.at(vec_p.size()-1);
    bdt ret_sum_dp;

    do {
        P_const = (P_min + P_max) / 2;
        vector<bdt> retData = run_sim(p_0, n_fl, stepsize, stop, vec_p, false, std::cout, 0);
        //std::cout << "run with end p=" << retData[0] << " and average sqrt(sum(delta p)^2) of " << pow(retData[3], 0.5) << ", ended at t=" << retData[4] << std::endl;
        //std::cout << "needed p_end=" << p_end << ", deltap=" << (retData[0]-p_end) << std::endl;
        delta_p = retData[0]-p_end;
        if (p_end>retData.at(0))
            P_min = (P_min + P_max) / 2;
        else
            P_max = (P_min + P_max) / 2;
        ret_sum_dp = retData[3];

    } while(delta_p_epsilon < abs((long)delta_p));

    std::cout << "found P_const=" << P_const << std::endl;

    return ret_sum_dp;
}

void optimize(bdt stepsize, vector<bdt> vec_p, bdt n_min, bdt n_max, bdt P_min, bdt P_max, bdt p_end_epsilon, bdt sum_delta_p_epsilon, bdt offset){

    bdt p_0 = vec_p.at(0);
    std::cout << "p_0:" << p_0 << std::endl;

    bdt stop = stepsize*(bdt)vec_p.size();

    bdt sum_dp;
    do {
        sum_dp = optimize_P(stepsize, p_0, stop, vec_p, (n_min+n_max)/2, P_min, P_max, sum_delta_p_epsilon);
        std::cout << "fit had sum_delta_p of " << sum_dp << std::endl;
        if (sum_dp<0)
            n_min = (n_min+n_max)/2;
        else
            n_max = (n_min+n_max)/2;


    } while(sum_delta_p_epsilon < abs((long)sum_dp));

    std::cout << "plotting latest fit ..." << std::endl;
    std::ofstream fout("fit_plot.txt");
    run_sim(p_0, (n_min+n_max)/2, stepsize, stop, vec_p, true, fout, offset);

}

int main(int argc, char* argv[]){

    ///Wir fangen da an, wo die Temperatur des Gases nahe der
    ///Dampfdruckkurvenenergie ist.
    bdt p_0   = 950 * 100; ///Anfangsdruck in mBar aka hPa *100
    bdt n_fl  = 1;  ///In Mol
    bdt T     = 77;

    bdt dt = 1;
    bdt stop = 120;

    if (argc == 1){
        vector<bdt> data;
        data.push_back(p_0);
        data.push_back(n_fl);
        data.push_back(T);
        for (bdt t=0; t<stop; t+=dt){
            std::cout<<t<<" "<<data.at(0)/100<<" "<<data.at(1)<<" "<<data.at(2)<<std::endl;
            data = RK4_schritt(data, dt);
        }
    }
    else if(argc == 2){
        std::string argument;
        argument = argv[1];
        if (argument == "plot"){
            std::cout << "Filename>";
            char inputName[256];
            std::cin.get(inputName, sizeof(inputName));

            bdt tempT;
            bdt tempP;

            vector<bdt> file_t;
            vector<bdt> file_p;
            std::cout << "reading file" << std::endl;
            std::ifstream fin(inputName);
            while (!fin.eof()){
                fin >> tempT;
                fin >> tempP;
                file_t.push_back(tempT);
                file_p.push_back(tempP*100);
            }
            std::cout << "read " << file_t.size() << " datasets." << std::endl;

            dt = file_t.back()/file_t.size();

            std::cout << "Stepsize (detected:'" << dt << "'):>"; std::cout.flush();
            std::string tempS;
            std::cin >> tempS;

            try {
                dt = stod(tempS);
            }
            catch (...){
                std::cout << "Ups." << std::endl;
                exit(1);
            }

            std::cout << "Offset>";
            bdt offset;
            std::cin >> tempS;

            try {
                offset = stod(tempS);
            }
            catch (...){
                std::cout << "Ups." << std::endl;
                exit(1);
            }

            file_p.erase(file_p.begin(), file_p.begin()+((int)(offset/dt)));

            optimize(dt, file_p, 1, 1000, 0.1, 200, 5, 5, offset);

        }
    }
    else{
        std::cout << "...>" << std::endl;
    }

}

