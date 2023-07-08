#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <cerrno>
#include <cfenv>
#include <cstring>
#include <algorithm>
#include <fstream>

using namespace Rcpp;

// [[Rcpp::export]]
// Function for discounted
double discounted(double undiscounted, 
                  double start_year, 
                  double end_year, 
                  double cr)
{
  double dis_out = undiscounted/ (exp(cr*(end_year-start_year)));
  return(dis_out);
}


// [[Rcpp::export]]
// Function to calculate basic CVD costs
double discount(double value, 
                double A, 
                double ar){
  double dis_out = value / (pow(1+ar, A));
  return(dis_out);
  }


// [[Rcpp::export]]
double compute_c_healthy(NumericVector bgcost_coef_vec, 
                         NumericVector cvdpct_coef_vec,
                         double age_mod,
                         double cAnnualFU_afterASCVD,
                         double cNonFatalASCVD
                         ) {
  
  double intercept_bgcost = bgcost_coef_vec[0];
  double slope_bgcost = bgcost_coef_vec[1];
  double intercept_cvdpct = cvdpct_coef_vec[0];
  double slope_cvdpct = cvdpct_coef_vec[1];
  
  double cvdpct = intercept_cvdpct + slope_cvdpct * age_mod;
  if (cvdpct < 0){
    cvdpct = 0;
  }
  
  double c_healthy = (intercept_bgcost + slope_bgcost*age_mod - cvdpct*(cAnnualFU_afterASCVD + cNonFatalASCVD/10)) / (1 - cvdpct);
  
  return(c_healthy);
}

// [[Rcpp::export]]
// Function to calculate CVD costs
List calc_cost_util(NumericVector death_of_ASCVD_vec,
                    NumericVector get_CVD_vec, 
                    NumericVector time_in_model_vec, 
                    NumericVector initial_age_vec,
                    NumericVector statins_use_vec, 
                    NumericVector statins_mild_adverse_event_vec, 
                    NumericVector statins_major_adverse_event_vec, 
                    int n_pop,
                    double cr,
                    double ar,
                    NumericVector bgcost_coef_vec, 
                    NumericVector cvdpct_coef_vec,
                    double cAnnualFU_afterASCVD,
                    double cNonFatalASCVD,
                    double cFatalASCVD,
                    double cAnnualStatin,
                    double cMildStatinAdverse,
                    double cMajorStatinAdverse,
                    double rInflation2017,
                    double uHealthy,
                    double uAfterASCVD,
                    double uHealthyStatin,
                    double uPenaltyMildStatinAdverse,
                    double uPenaltyMajorStatinAdverse) {
  // Output
  NumericVector cost_vec(n_pop);
  NumericVector cost_disc_vec(n_pop);
  NumericVector util_disc_vec(n_pop);
  
  // looping over population
  for (int id = 0; id < n_pop; id++){
    // Read individual characteristics
    double death_of_ASCVD = death_of_ASCVD_vec(id);
    double get_CVD = get_CVD_vec(id);
    double time_in_model = time_in_model_vec(id);
    int initial_age = initial_age_vec(id);
    bool statins_use = statins_use_vec(id);
    double statins_mild_adverse_event = statins_mild_adverse_event_vec(id);
    double statins_major_adverse_event = statins_major_adverse_event_vec(id);
    
    // cHealthcareNonCVD
    double age_end_Healthy;
    if (NumericVector::is_na(get_CVD)){
      age_end_Healthy = time_in_model + initial_age;
    } else {
      age_end_Healthy = std::min(time_in_model, get_CVD) + initial_age;
    }
    
    for (double age_temp = initial_age; age_temp <= age_end_Healthy; age_temp++){
      cost_vec(id) += compute_c_healthy(bgcost_coef_vec, cvdpct_coef_vec, age_temp, cAnnualFU_afterASCVD,
                                      cNonFatalASCVD);
    }
    
    cost_vec(id) += (age_end_Healthy - floor(age_end_Healthy)) * compute_c_healthy(
      bgcost_coef_vec, cvdpct_coef_vec, ceil(age_end_Healthy), cAnnualFU_afterASCVD,
             cNonFatalASCVD);
    
    cost_disc_vec(id) += discounted(cost_vec(id), 0, age_end_Healthy - initial_age, cr);
    
    // FatalASCVD
    if(!NumericVector::is_na(get_CVD) && !NumericVector::is_na(death_of_ASCVD)) {
      cost_disc_vec(id) += discount(rInflation2017 * cFatalASCVD, death_of_ASCVD, ar);
    }
      
    //Non-Fatal ASCVD
    if(!NumericVector::is_na(get_CVD) && NumericVector::is_na(death_of_ASCVD)) {
      cost_disc_vec(id) += discount(rInflation2017 * cNonFatalASCVD, get_CVD, ar) + 
        discounted(rInflation2017 * cAnnualFU_afterASCVD * (time_in_model - get_CVD - 0.5), 
                   get_CVD + 0.5, time_in_model, cr);
    }
    
    // statins fee and cost for adverse events
    if (statins_use){
      double uStatins;
      
      if(!NumericVector::is_na(statins_major_adverse_event)){
        cost_disc_vec(id) += discount(cMajorStatinAdverse, 0, ar);
        uStatins = uHealthyStatin - uPenaltyMajorStatinAdverse;
      } else if(!NumericVector::is_na(statins_mild_adverse_event)){
        cost_disc_vec(id) += discount(cMildStatinAdverse, 0, ar);
        uStatins = uHealthyStatin - uPenaltyMildStatinAdverse;
      } else {
        uStatins = uHealthyStatin;
      }
        
      double yrs_end_Statins;
      if (NumericVector::is_na(get_CVD)){
        yrs_end_Statins = time_in_model;
      } else {
        yrs_end_Statins = std::min(time_in_model, get_CVD);
      }
      cost_disc_vec(id) += discounted(cAnnualStatin * yrs_end_Statins * rInflation2017, 0, yrs_end_Statins, cr);
      util_disc_vec(id) += discounted(uStatins * yrs_end_Statins, 0, yrs_end_Statins, cr) + 
        discounted(uAfterASCVD * (time_in_model - yrs_end_Statins), yrs_end_Statins, time_in_model, cr);
    } else {
      double yrs_end_Healthy;
      if (NumericVector::is_na(get_CVD)){
        yrs_end_Healthy = time_in_model;
      } else {
        yrs_end_Healthy = std::min(time_in_model, get_CVD);
      }
      util_disc_vec(id) += discounted(uHealthy * yrs_end_Healthy, 0, yrs_end_Healthy, cr) + 
        discounted(uAfterASCVD * (time_in_model - yrs_end_Healthy), yrs_end_Healthy, time_in_model, cr);
    }
  }
  
  return(List::create(
      Named("util_disc") = util_disc_vec,
      Named("cost")      = cost_vec,
      Named("cost_disc") = cost_disc_vec
  )
  );
}