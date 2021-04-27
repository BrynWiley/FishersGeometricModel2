#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <stdlib.h> 
#include <iostream>
#include <random>
#include <chrono>
using namespace Rcpp;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//[[Rcpp::export]]
IntegerVector strip_even(IntegerVector v){
  IntegerVector x = unique(v);
  int n = x.size();
  std::vector<int> stripped;
  
  for(int i=0; i<n; i++){
    int num = std::count(v.begin(),v.end(), x[i]);
    if(num % 2 != 0){
      stripped.push_back(x[i]);
    }
  }
  return(IntegerVector::import(stripped.begin(),stripped.end()));
}

double norm_vec(NumericVector x){
  int n = x.size();
  double out = 0;
  
  for(int i = 0; i < n; i++) {
    out += pow(x[i],2);
  }
  return sqrt(out);
}

double fitness(NumericVector optimum, NumericVector z, int k){
  double norm=norm_vec(optimum-z);
  if(k==10){
    if(norm>1){
      return(exp(-1*pow(norm,6)));
    } else {
      return(exp(-1*pow(norm,2)));
    }
  }
  if(k==8){
    if(norm>1){
      return(exp(-1*pow(norm,4)));
    } else {
      return(exp(-1*pow(norm,2)));
    }
  }
  return(exp(-1*pow(norm,k)));
}



//must already have order trimmed!
//[[Rcpp::export]]
List mutational_differences(List order_one, List order_two, NumericMatrix mutational_library, int l, int k,int seed){
  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_real_distribution<double> real_distribution(0.0,1.0);
  std::uniform_int_distribution<int> int_distribution(0,l-1);
  
  int n1 = order_one.size();
  int n2 = order_two.size();
  int num_comparisons_one=0;
  int num_comparisons_both=0;
  
  int num_F1_comparisons_one=0;
  int num_F1_comparisons_both=0;
  
  double mutational_diffs_one = 0.0;
  double overall_diffs_one = 0.0;
  double common_distances_one = 0.0;
  double F1_fitnesses_one = 0.0;
  double F1_distances_one = 0.0;
  

  
  double mutational_diffs_both = 0.0;
  double overall_diffs_both = 0.0;
  double common_distances_both = 0.0;
  double F1_fitnesses_both = 0.0;
  double F1_distances_both = 0.0;
  
  IntegerVector genome = seq(0,l-1);
  NumericVector optimum = NumericVector::create(2.0,0.0,0.0,0.0,0.0);
  
  //within order_one orders
  for(int i = 1; i < n1; i++){
    for(int ii = 0; ii < i; ii++){
      //calculate intersect and genomic inverses of order
      IntegerVector order1 = strip_even(order_one[i]) - 1;
      IntegerVector order2 = strip_even(order_one[ii]) - 1;
      
      
      
      IntegerVector inter = intersect(order1,order2);
      
      
      
      if(inter.size() > 0){
        double mutational_diff = (double)inter.size()/((double)order1.size()+(double)order2.size());
        
        mutational_diffs_one = mutational_diff + mutational_diffs_one;
        
        NumericVector common_z = NumericVector::create(0.0,0.0,0.0,0.0,0.0);
        for(int iii = 0; iii < inter.size(); iii++){
          common_z = common_z + mutational_library(inter[iii],_);
        }
        
        common_distances_one += norm_vec(common_z-optimum);
      } else {
        common_distances_one += 2.0;
      }
      double overall_diff_one = ((double)order1.size()+(double)order2.size()-(double)inter.size()-(double)inter.size())/(double)l;
      overall_diffs_one = overall_diffs_one + overall_diff_one;
      
      for(int iii = 0; iii < 100; iii++){
        NumericVector F1_z = NumericVector::create(0.0,0.0,0.0,0.0,0.0);
        IntegerVector indexes = IntegerVector::create(-1);
        for(int iv = 0; iv < l; iv++){
          if(real_distribution(generator) < 0.5){
            indexes.push_back(iv);
          }
        }
        indexes.erase(0);
        IntegerVector F1_order1 = intersect(order1,indexes);
        IntegerVector F1_order2 = setdiff(order2,indexes);
        
        IntegerVector F1_order = union_(F1_order1,F1_order2);
        if(F1_order.size() > 0){
          for(int iv = 0; iv < F1_order.size(); iv++ ){
            F1_z = F1_z + mutational_library(F1_order[iv],_);
          }
        }
        F1_fitnesses_one = F1_fitnesses_one + fitness(optimum, F1_z,k);
        F1_distances_one = F1_distances_one + norm_vec(optimum-F1_z);
        num_F1_comparisons_one++;
      }
      num_comparisons_one++;
    }
  }
  

  
  //between both orders
  for(int i = 0; i < n1; i++){
    for(int ii = 0; ii < n2; ii++){
      //calculate intersect and genomic inverses of order
      IntegerVector order1 = strip_even(order_one[i]) - 1;
      IntegerVector order2 = strip_even(order_two[ii]) - 1;
      
      
      
      IntegerVector inter = intersect(order1,order2);
      
      
      
      if(inter.size() > 0){
        double mutational_diff = (double)inter.size()/((double)order1.size()+(double)order2.size());
        
        mutational_diffs_both = mutational_diff + mutational_diffs_both;
        
        NumericVector common_z = NumericVector::create(0.0,0.0,0.0,0.0,0.0);
        for(int iii = 0; iii < inter.size(); iii++){
          common_z = common_z + mutational_library(inter[iii],_);
        }
        
        common_distances_both += norm_vec(common_z-optimum);
      } else {
        common_distances_both += 2.0;
      }
      double overall_diff_both = ((double)order1.size()+(double)order2.size()-(double)inter.size()-(double)inter.size())/(double)l;
      overall_diffs_both = overall_diffs_both + overall_diff_both;
      
      for(int iii = 0; iii < 100; iii++){
        NumericVector F1_z = NumericVector::create(0.0,0.0,0.0,0.0,0.0);
        IntegerVector indexes = IntegerVector::create(-1);
        for(int iv = 0; iv < l; iv++){
          if(real_distribution(generator) < 0.5){
            indexes.push_back(iv);
          }
        }
        indexes.erase(0);
        IntegerVector F1_order1 = intersect(order1,indexes);
        IntegerVector F1_order2 = setdiff(order2,indexes);
        
        IntegerVector F1_order = union_(F1_order1,F1_order2);
        if(F1_order.size() > 0){
          for(int iv = 0; iv < F1_order.size(); iv++ ){
            F1_z = F1_z + mutational_library(F1_order[iv],_);
          }
        }
        F1_fitnesses_both = F1_fitnesses_both + fitness(optimum, F1_z,k);
        F1_distances_both = F1_distances_both + norm_vec(optimum-F1_z);
        num_F1_comparisons_both++;
      }
      num_comparisons_both++;
    }
  }
  
  
  return List::create(
    _["mean_mutational_within"]=(double)mutational_diffs_one/(double)num_comparisons_one,
    _["mean_overall_within"]=(double)overall_diffs_one/(double)num_comparisons_one,
    _["mean_common_distance_within"] = (double)common_distances_one/(double)num_comparisons_one,
    _["mean_F1_fitness_within"]=(double)F1_fitnesses_one/(double)num_F1_comparisons_one,
    _["mean_F1_distance_within"]=(double)F1_distances_one/(double)num_F1_comparisons_one,
    _["mean_mutational_between"]=(double)mutational_diffs_both/(double)num_comparisons_both,
    _["mean_overall_between"]=(double)overall_diffs_both/(double)num_comparisons_both,
    _["mean_common_distance_between"] = (double)common_distances_both/(double)num_comparisons_both,
    _["mean_F1_fitness_between"]=(double)F1_fitnesses_both/(double)num_F1_comparisons_both,
    _["mean_F1_distance_between"]=(double)F1_distances_both/(double)num_F1_comparisons_both
  );
}
