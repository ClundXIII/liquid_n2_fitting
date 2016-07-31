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
        bdt get_avrg_sum_sqared_p();

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

    private:
};

#endif // SIMULATION_H
