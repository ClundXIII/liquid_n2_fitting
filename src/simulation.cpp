#include "simulation.h"

#include <math.h>
#include <fstream>
#include <sstream>

#define a_1 ((c_wN2*R_m*R_m)/(d_H*d_H))
#define a_2 ((c_wN2*R_m)/(d_H*d_H))

simulation::simulation(bdt p_0, bdt n_fl, bdt dt, bdt stop, vector<bdt> p_mess, bool print, std::ostream& print_s, bdt offset, bdt P_const) :
    print_s(print_s)
    {
    this->p_0 = p_0;
    this->n_fl = n_fl;
    this->dt = dt;
    this->stop = stop;
    this->p_mess = p_mess;
    this->print = print;
    //(this->print_s) = &print_s;
    this->offset = offset;
    this->P_const = P_const;
}

int simulation::loadSettingsFile(char s[]){
    return loadSettingsFile((const char*)s);
}

int simulation::loadSettingsFile(const char s[]){
    std::string tempS(s);
    return loadSettingsFile(tempS);
}

int simulation::loadSettingsFile(std::string s){

    std::cout << "loading config " << s << " ..." << std::endl;

    std::ifstream f(s.c_str());

    if (!f.good())
        return -1;

    if (f.eof())
        return 1;

    try {
        while (!f.eof()){
            std::string tempS;
            f >> tempS;
            if (tempS == "V"){
                f >> V;
            }
            else if (tempS == "dV_pump"){
                f >> dV_pump;
            }
            else if (tempS == "n_leak"){
                f >> n_leak;
            }
            else {
                //Skip the wrong option
                f >> tempS;
                std::cout << "unknown option: \"" << tempS << "\"!" << std::endl;
            }
        }
    }
    catch (...){
        std::cout << "Unknown error while loading config File!" << std::endl;
        return -2;
    }

    return 0;
}

simulation::~simulation(){
    //
}

vector<bdt> simulation::run_sim(){
    vector<bdt> data;
    data.push_back(p_0);
    data.push_back(n_fl);
    data.push_back(1 /( ( ((-R_m /d_H)*log(data[0]/p_t)) + (1./T_t) ) ));

    std::cout << "start values: dt=" << dt << ", p=" << data.at(0) << ", n_fl=" << data.at(1) << ", T=" << data.at(2) << ", P_const=" << P_const << std::endl;

    unsigned int it_pos=0;

    deltaInsg   = 0;
    deltaSqInsg = 0;
    bdt       t = 0;
    for (; t<=stop; t+=dt){
        data = RK4_schritt(data, dt);
        deltaInsg += data.at(0)-p_mess.at(it_pos);
        deltaSqInsg += abs(data.at(0)-p_mess.at(it_pos));
        if (print) print_s<<t+offset<<" "<<data.at(0)/100<<" "<<data.at(1)<<" "<<data.at(2)<<std::endl;
        if (it_pos<(p_mess.size()-1))
            it_pos++;

        if (data.at(1)<=0)
            break;
    }
    deltaInsg /= (bdt) it_pos;
    data.push_back(deltaInsg);
    data.push_back(t);

    return data;
}

bdt simulation::get_avrg_sum_p(){
    return deltaInsg;
}

bdt simulation::get_avrg_sum_squared_p(){
    return deltaSqInsg;
}

vector<bdt> simulation::f_strich(vector<bdt> input){
    vector<bdt> retV;

    bdt p    = input.at(0);
    bdt n_fl = input.at(1);
    bdt T    = input.at(2);

    bdt _p     = ( R_m*T*((P_const)/(d_H)+n_leak) + p*dV_pump ) / (V + T*T*T*n_fl*a_1/p);
    bdt _n_fl  = -(P_const/d_H) + n_fl * a_2 * ((T*T)/p)*_p;
    bdt _T     = (R_m * T*T*_p)/(d_H * p);

    retV.push_back(_p);
    retV.push_back(_n_fl);
    retV.push_back(_T);

    return retV;
}

vector<bdt> simulation::RK4_schritt(vector<bdt> data, bdt dt){

    vector<bdt> k1 = f_strich(data);
    vector<bdt> k2 = f_strich(data + (k1 * dt/2));
    vector<bdt> k3 = f_strich(data + (k2 * dt/2));
    vector<bdt> k4 = f_strich(data+ k3);

    data = data + ( (k1 + (k2 + k3) * 2 + k4) * dt/6);

    ///Korrigiere Temperatur:
    if (correctTemp)
        data[2] = 1 / ( (((-R_m /d_H)*log(data[0]/p_t))+ (1./T_t) ) );

    return data;
}

//Initialize default params
bdt simulation::V       = 0.06252;   //Volume of the chamber
bdt simulation::dV_pump = -.00134;   //Pump characteristic in m^3/second
bdt simulation::n_leak  = 0.0000345; //Leakage of the chamber
