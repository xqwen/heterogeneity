using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "controller.h"
#include <gsl/gsl_cdf.h>
#include "BF_CEF.h"

#define NINF -999999

void controller::load_data(char *filename, int use_zval){






    ifstream dfile(filename);
    string line;
    istringstream ins;


    string loc_id;
    double beta;
    double se_beta;

    int loc_count = 0;

    double min = 1;
    double max = 0;;

    vector<vector<double> > beta_matrix;
    vector<vector<double> > se_matrix;

    while(getline(dfile,line)){
        ins.clear();
        ins.str(line);

        ins>>loc_id;

        vector<double> beta_vec;
        vector<double> se_vec;

        while(ins >> beta){
            if(!use_zval)
                ins>>se_beta;
            else 
                se_beta = 1;

            if(use_zval == -1) // pvalue is used
                beta = gsl_cdf_ugaussian_Qinv (beta/2);

            beta_vec.push_back(beta);
            se_vec.push_back(se_beta);


            loc_vec.push_back(loc_id);

            if(se_beta<min)
                min = se_beta;

            if(beta*beta - se_beta*se_beta > max){
                max = pow(beta,2)-pow(se_beta,2);
            }


        }

        beta_matrix.push_back(beta_vec);
        se_matrix.push_back(se_vec);
    }

    dfile.close();

    // compute BF vector
    double phi_min = min/10;
    double phi_max = 2*sqrt(max);
    if(phi_max>5)
        phi_max = 5;
    if(phi_max<phi_min){
        phi_max = 8*phi_min;
    }

    //make_grid(phi_min, phi_max);
    k_vec.push_back(0.4);
    k_vec.push_back(1);
    k_vec.push_back(10);

    omg2_vec.push_back(1e-8);
    omg2_vec.push_back(1e-8);
    omg2_vec.push_back(1e-8);

    K = k_vec.size()+1;
    N = beta_matrix.size();


    fprintf(stderr, "Initializing ... \n");


    //fprintf(stderr, "N=%d\t K=%d\t L=%d\t  M=%d\n", N, grid_size, annot_size, K);

    // compute BF factors


    //    log10_BF_matrix.push_back(bf_vec);  
    BF_CEF bfc;
    for(int i=0;i<N;i++){
        vector<double> bf_vec = bfc.compute_log10_BF(beta_matrix[i], se_matrix[i], k_vec, omg2_vec);
        for(int k=0;k<bf_vec.size();k++){
            printf("%7.3f ",bf_vec[k]);
        }
        printf("\n");
    }
}


void controller::make_grid(double phi_min, double phi_max){

    vector<double> gvec;
    gvec.push_back(phi_max); //phi_max
    double phi = phi_max;
    while(phi>phi_min){
        phi = phi/sqrt(2);
        gvec.push_back(phi); //phi
    }
    std::sort(gvec.begin(),gvec.end());



    vector<double> kvalue_vec{0.01, 0.462, 0.513, 0.618, 0.820, 1.0, 2.0, 10.0};
    for(int j=0;j<kvalue_vec.size();j++){
        for(int i=0;i<gvec.size();i++){
            k_vec.push_back(kvalue_vec[j]);
            double omg2 = pow(gvec[i],2)/(1+pow(kvalue_vec[j],2));
            omg2_vec.push_back(omg2);
        }
    }

    /*
    for(int i=0;i<k_vec.size();i++){
        printf("%3d   %7.3f  %7.3f\n",i+1, k_vec[i], omg2_vec[i]);
    }
    */
    return;

}

