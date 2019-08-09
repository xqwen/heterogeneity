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
        loc_vec.push_back(loc_id);

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
    /*
       k_vec.push_back(0);
       k_vec.push_back(0);

       omg2_vec.push_back(0);
       omg2_vec.push_back(1);
       */

    make_grid(phi_min, phi_max);
    K = k_vec.size()+1;
    N = beta_matrix.size();


    //fprintf(stderr, "Initializing ... \n");


    //fprintf(stderr, "N=%d\t K=%d\t L=%d\t  M=%d\n", N, grid_size, annot_size, K);

    // compute BF factors


    //    log10_BF_matrix.push_back(bf_vec);  
    BF_CEF bfc;
    for(int i=0;i<N;i++){

        vector<double> bf_vec = bfc.compute_log10_BF(beta_matrix[i], se_matrix[i], k_vec, omg2_vec);
        bf_vec.insert(bf_vec.begin(),0.0);
        log10_BF_matrix.push_back(bf_vec);

        /*
        //printf("%4d  ", i+1);
        //for(int k=0;k<bf_vec.size();k++){
        //    printf("%7.3f ",bf_vec[k]);
        //}
        //printf("\n");
        */
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


    vector<double> ovalue_vec{0.1, 0.4, 1, 1.6, 3.2, 6.4};
    vector<double> kvalue_vec{0.01, 0.462, 0.513, 0.618, 2.0, 10.0};
    for(int j=0;j<kvalue_vec.size();j++){
        for(int i=0;i<ovalue_vec.size();i++){
            k_vec.push_back(kvalue_vec[j]);
            //double omg2 = pow(gvec[i],2)/(1+pow(kvalue_vec[j],2));
            omg2_vec.push_back(pow(ovalue_vec[i],2));
        }
    }

    /*
       for(int i=0;i<k_vec.size();i++){
       printf("%3d   %7.3f  %7.3f\n",i+1, k_vec[i], omg2_vec[i]);
       }
       */
    return;

}


void controller::run_EM(double thresh){

    vector<double> wts_vec(K, 0.05/(K-1));
    wts_vec[0] = 0.95;
    double final_loglik = gem.EM_run(log10_BF_matrix, wts_vec, thresh);


    wts_vec = gem.get_estimate();
    double null_prob = wts_vec[0];
    double rep_prob = 0;
    double irp_prob = 0;
    for(int i=1;i<wts_vec.size();i++){
        if(k_vec[i-1]<=0.618)
            rep_prob += wts_vec[i];
        else
            irp_prob += wts_vec[i];
    }


    fprintf(stderr, "\n\n%15s %7.3f\n%15s %7.3f\n%15s %7.3f\n\n", "Null:", null_prob, "Reproducible:", rep_prob, "Irreproducible:", irp_prob);;


    // get posterior probabilities
    vector<vector<double> > P_matrix = gem.get_P_matrix();
    //vector<double> rp_vec;
    for(int i=0;i<N;i++){
        double rep_prob = 0;
        double irp_prob = 0;
        for(int k=1;k<K;k++){
            if(k_vec[k-1]<=0.618)
                rep_prob += P_matrix[i][k];
            else
                irp_prob += P_matrix[i][k];
        }
    
        printf("%10s   %7.3e %7.3e %7.3e\n", loc_vec[i].c_str(), P_matrix[i][0], irp_prob, rep_prob );
    }



}





