#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <string>

#include "define.h"
#include "vector.h"

using namespace vemc2::mymath;

class simulation{
    public:
        simulation(bdt p_0, bdt n_fl, bdt dt, bdt stop, vector<bdt> p_mess, bool print, std::ostream& print_s, bdt offset, bdt P_const);
        static int loadSettingsFile(std::string s);
        static int loadSettingsFile(char s[]);
        static int loadSettingsFile(const char s[]);
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
        static bdt V;       //Volume of the chamber
        static bdt dV_pump; //Pump characteristic in m^3/second
        static bdt n_leak;  //Leakage of the chamber

        bdt deltaInsg=0;
        bdt deltaSqInsg=0;

    private:
};

#endif // SIMULATION_H
