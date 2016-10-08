#include "fit.h"

#include <fstream>

fit::fit(bdt stepsize, vector<bdt> vec_p, bdt n_min, bdt n_max, bdt P_min, bdt P_max, bdt p_end_epsilon, bdt sum_delta_p_epsilon, bdt offset){
    this->stepsize = stepsize;
    this->vec_p = vec_p;
    this->n_min = n_min;
    this->n_max = n_max;
    this->P_min = P_min;
    this->P_max = P_max;
    this->p_end_epsilon = p_end_epsilon;
    this->sum_delta_p_epsilon = sum_delta_p_epsilon;
    this->offset = offset;
}


bdt fit::optimize_P(bdt stepsize, bdt p_0, bdt stop, vector<bdt> vec_p, bdt n_fl, bdt P_min, bdt P_max, bdt delta_p_epsilon){


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

void fit::optimize(){
    optimize(stepsize, vec_p, n_min, n_max, P_min, P_max, p_end_epsilon, sum_delta_p_epsilon, offset);
}

void fit::optimize(bdt stepsize, vector<bdt> vec_p, bdt n_min, bdt n_max, bdt P_min, bdt P_max, bdt p_end_epsilon, bdt sum_delta_p_epsilon, bdt offset){

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


fit::~fit()
{
    //dtor
}
