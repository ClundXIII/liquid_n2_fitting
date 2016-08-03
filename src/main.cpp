#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <fstream>

#include "define.h"
#include "vector.h"
#include "simulation.h"

using namespace vemc2::mymath;

bdt P_const = 10;         //10 W Leistung durch Boden

bdt optimize_P(bdt stepsize, bdt p_0, bdt stop, vector<bdt> vec_p, bdt n_fl, bdt P_min, bdt P_max, bdt delta_p_epsilon){

    bdt delta_p;
    bdt p_end = vec_p.at(vec_p.size()-1);
    bdt ret_sum_dp;

    do {
        P_const = (P_min + P_max) / 2;
        simulation s(p_0, n_fl, stepsize, stop, vec_p, false, std::cout, 0, P_const);
        vector<bdt> retData = s.run_sim();
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
    simulation s(p_0, (n_min+n_max)/2, stepsize, stop, vec_p, true, fout, offset, P_const);
    s.run_sim();

}

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
            //data = RK4_schritt(data, dt);
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

            dt = file_t.back()/file_t.size();

            std::cout << "Stepsize (detected:'" << dt << "'):>"; std::cout.flush();

            dt = read_bdt();

            std::string tempS;

            std::cout << "Offset>";

            bdt offset = read_bdt();

            file_p.erase(file_p.begin(), file_p.begin()+((int)(offset/dt)));

            optimize(dt, file_p, 1, 1000, 0.1, 200, 5, 5, offset);

        }
    }
    else{
        std::cout << "...>" << std::endl;
    }

}

