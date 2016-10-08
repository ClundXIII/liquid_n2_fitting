#ifndef FIT_H
#define FIT_H

#include <vector>

#include "simulation.h"

extern bdt P_const;

class fit
{
    public:
        fit(bdt stepsize, vector<bdt> vec_p, bdt n_min, bdt n_max, bdt P_min, bdt P_max, bdt p_end_epsilon, bdt sum_delta_p_epsilon, bdt offset);

        void optimize();
        static void optimize(bdt stepsize, vector<bdt> vec_p, bdt n_min, bdt n_max, bdt P_min, bdt P_max, bdt p_end_epsilon, bdt sum_delta_p_epsilon, bdt offset);
        static bdt optimize_P(bdt stepsize, bdt p_0, bdt stop, vector<bdt> vec_p, bdt n_fl, bdt P_min, bdt P_max, bdt delta_p_epsilon);

        virtual ~fit();
    protected:
    private:

        bdt stepsize;
        vector<bdt> vec_p;
        bdt n_min;
        bdt n_max;
        bdt P_min;
        bdt P_max;
        bdt p_end_epsilon;
        bdt sum_delta_p_epsilon;
        bdt offset;

};

#endif // FIT_H
