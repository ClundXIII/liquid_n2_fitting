#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>

#include "define.h"
#include "vector.h"

using namespace vemc2::mymath;

class simulation{
    public:
        simulation(bdt p_0, bdt n_fl, bdt dt, bdt stop, vector<bdt> p_mess, bool print, std::ostream& print_s, bdt offset, bdt P_const);
        virtual ~simulation();

        vector<bdt> run_sim(void);

        bdt get_avrg_sum_p();
        bdt get_avrg_sum_squared_p();

        bool correctTemp = true;

    protected:

        bdt p_0;
        bdt n_fl;
        bdt dt;
        bdt stop;
        vector<bdt> p_mess;
        bool print;
        std::ostream& print_s;
        bdt offset;
        bdt P_const;

        vector<bdt> f_strich(vector<bdt> input);
        vector<bdt> RK4_schritt(vector<bdt> data, bdt dt);

        ///Constants of the chamber
        bdt V       = 0.06252;   //Volume of the chamber
        bdt dV_pump = -.00134;   //Pump characteristic in m^3/second
        bdt n_leak  = 0.0000345; //Leakage of the chamber

        bdt deltaInsg=0;
        bdt deltaSqInsg=0;

    private:
};

#endif // SIMULATION_H
