#include "BF_CEF.h"
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cdf.h>
#include <stdio.h>

using namespace std;

double target_f (double x, void *params){
  
  double *p = (double *)params;
  double val = 1;
  int size = int(p[0]);
  double k = p[1];
  double oa2 = p[2];
  //printf("init %f %f %f\n",double(size),k,oa2);
  for(int i=0;i<size;i++){
    double ds2 = pow(p[2*i+4],2.0);
    double beta = p[2*i+3];
    
    double num2 = pow(beta-x,2);
    double dnum = ds2 + k*k*x*x;

    double fac = 1.0/(2*3.141592653*dnum);
    val *= sqrt(fac)*exp(-0.5*num2/dnum);
    //printf ("val = %e\n",val);
  }
  val *= exp(-.5*x*x/oa2);
  return val;

}
 

int BF_CEF::prepare_params(double k, double oa2){
  
  // find non-informative sets
  int info_size = 0;
  
  for(int i=0;i<sd_vec.size();i++){
    if(sd_vec[i] > 0)
      info_size++;
  }
  
  
  if(info_size==0)
    return info_size;



  param_list = new double[2*info_size+3];
  
  param_list[0] = info_size;
  param_list[1] = k;
  param_list[2] = oa2;
  
  double val = 0;
  int i=0;
  for(int j=0;j<bhat_vec.size();j++){
    
    if(sd_vec[j]<= 0)
      continue; 
    param_list[2*i+3] = bhat_vec[j];
    param_list[2*i+4] = sd_vec[j];
    i++;
  }
  
  return info_size;
}


double BF_CEF::compute_log10_BF(double k, double oa2){
  
    if(oa2<1-10)
        oa2 = 1e-10;


  int info_size = prepare_params(k,oa2);
  
  if(info_size==0)
    return 0.0;
  
  // numerical integration
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (5000);
  double result, error;
  gsl_function F;
  F.function = &target_f;
  F.params = param_list;
    
  gsl_integration_qagi (&F, 0, 1e-7, 5000,w, &result, &error); 

  gsl_integration_workspace_free (w);
  delete[] param_list;

printf("result = %7.3e\n", result);

  double logBF = log(result)-.5*log(oa2)-.5*log(2*3.1415926536);
  
  
  return logBF/log(10.0);
  
  
}


vector<double> BF_CEF::compute_log10_BF(vector<double> & bhat, vector<double> & sd, vector<double> & k_vec, vector<double> &omg2_vec){
 
    bhat_vec = bhat;
    sd_vec = sd;
    vector<double> rst_vec;
    for(int i=0;i<k_vec.size();i++){
        double rst = compute_log10_BF(k_vec[i], omg2_vec[i]);
        rst_vec.push_back(rst);
    }
    return rst_vec;

}


