#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>

#include "define.h"
#include "vector.h"
#include "simulation.h"
#include "fit.h"

using namespace vemc2::mymath;

bdt P_const;

bdt read_bdt(){

    bool fail = true;

    bdt tempVar;

    while (fail){

        fail = false;
        try {
            std::cin >> tempVar;
        }
        catch (...){
            std::cout << "invalid input, try again>";
            fail = true;
        }
    }
    return tempVar;
}

int main(int argc, char* argv[]){

    //load Settings:
    simulation::loadSettingsFile("chamber.ini");

    if (argc == 1){
        ///We start where the Temperature of the Gas is already close to
        ///the boiling curve
        bdt p_0   = 950 * 100; ///Starting pressure in mBar aka hPa *100
        bdt n_fl  = 100;       ///In Mol
        bdt T     = 77;        ///In K

        bdt dt = 1;            ///Stepping
        bdt stop = 120;        ///Run Simulation for 120 Seconds

        vector<bdt> data;
        data.push_back(p_0);
        data.push_back(n_fl);
        data.push_back(T);

        bdt constP = 150;

        for (bdt t=0; t<stop; t+=dt){
            std::cout<<t<<" "<<data.at(0)/100<<" "<<data.at(1)<<" "<<data.at(2)<<std::endl;
            vector<bdt> k1 = simulation::f_strich(data, constP);
            vector<bdt> k2 = simulation::f_strich(data + (k1 * dt/2), constP);
            vector<bdt> k3 = simulation::f_strich(data + (k2 * dt/2), constP);
            vector<bdt> k4 = simulation::f_strich(data+ k3, constP);

            data = data + ( (k1 + (k2 + k3) * 2 + k4) * dt/6);

            ///correct Temperatur:
            data[2] = 1 / ( (((-R_m /d_H)*log(data[0]/p_t))+ (1./T_t) ) );
        }
    }
    else if(argc == 2){
        std::string argument;
        argument = argv[1];
        if (argument == "fit"){
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

            bdt dt = file_t.back()/file_t.size();

            std::cout << "Stepsize (detected:'" << dt << "'):>"; std::cout.flush();

            dt = read_bdt();

            std::string tempS;

            std::cout << "Offset>";

            bdt offset = read_bdt();

            file_p.erase(file_p.begin(), file_p.begin()+((int)(offset/dt)));

            //fit::optimize(dt, file_p, 1, 1000, 0.1, 200, 5, 5, offset);

            fit *f = new fit(dt, file_p, 1, 1000, 0.1, 200, 5, 5, offset);

            f->optimize();

        }
    }
    else{
        std::cout << "...>" << std::endl;
    }

}

