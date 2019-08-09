using namespace std;

#include <vector>
#include <string>
#include "GenEM_mixture.h"

class controller {

    public:
        
        void load_data(char *data_filei, int use_zval);
        void run_EM(double thresh=0.05);              
    
    private:

        int K;
        int N;
        vector<string> loc_vec;
        vector<double> k_vec;
        vector<double> omg2_vec;

        vector<vector<double> > log10_BF_matrix;

        GenEM_mixture gem;
        void make_grid(double min, double max);

};

