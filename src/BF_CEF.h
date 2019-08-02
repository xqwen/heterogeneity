using namespace std;

#include <vector>

class BF_CEF {

    private:

        vector<double> bhat_vec;
        vector<double> sd_vec;
        double *param_list;

        double compute_log10_BF(double k, double oa2);
        int prepare_params(double k, double oa2);

    public:

        vector<double> compute_log10_BF(vector<double> &beta_vec, vector<double> &sde_vec, vector<double> &k_vec, vector<double> &omega_vec);

};

double target_f (double x, void *params);
