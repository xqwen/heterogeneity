using namespace std;

#include <vector>
#include <string>


class controller {

    public:
        
        void load_data(char *data_filei, int use_zval);
               
    
    private:

        int K;
        int N;
        vector<string> loc_vec;
        vector<double> k_vec;
        vector<double> omg2_vec;

        vector<vector<double> > log10_BF_matrix;

        void make_grid(double min, double max);

};

