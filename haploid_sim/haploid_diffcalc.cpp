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

//This file contains functions used to produce haploid F1 hybrids from populations generated
//using haploid_sim.cpp functions. Specifically, it creates these hybrids between all possible 
//population pairs present in a population, and measures various hybridization and parental 
//adaptation statistics.

//A simple function to remove any components of an integer vector 
//that occur an even number of times 
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

//A simple function to generate the 2-norm of a vector
double norm_vec(NumericVector x){
  int n = x.size();
  double out = 0;
  
  for(int i = 0; i < n; i++) {
    out += pow(x[i],2);
  }
  return sqrt(out);
}

//A function to generate the fitness of an individual 
//Inputs: 
//optimum: the phenotypic optimum, used for fitness calculation
//z: the phenotype of an individual, represented as a vector with trait values as components
//k: the k value for the fitness function
//Output: The fitness of an individual

//Note, k values of 10 and 8 represent mixed fitness functions, with k values of 2 for distances
//less than 1 from the optimum, and k values of 6 and 4 respectively for distances greater than 1
// from the optimum
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



//A function that generates hybridization and parental adaptation statistics by generating F1 hybrids from all possible
//pairs of populations in the population set listed. Also measures parallel evolution statistics between population pairs.
//Inputs:
//order: A list of fixation events in each population. Each element in the list needs to be an integer vector,
//with each element of the vector the id of the loci whose mutation fixed during the fixation event in order of 
//fixation. The id refers to the row number of mutational_library. This row specified is the mutational vector
//of the loci specified.
//mutational_library: a numeric matrix containing the mutational vectors associated with each loci in the genome of 
//the populations described by the order input. Needs to have l rows and d columns
//optimum: a numeric vector representing the phenotypic optimum for the populations
//l: the number of loci contributing to adaptation
//d: the number of trait dimensions under consideration
//k: the k value for fitness calculations
//seed: an integer used to set the random number generator seeds, for reproducibility
//Output:
//A list with 5 components:
//mean_mutational: the average proportion of common fixed mutations between all pairs of populations
//mean_overall: the average portion of the genome common between all pairs of populations
//mean_common_distance: the average phenotypic distance to the optimum only considering fixed mutations common between all pairs of populations
//mean_F1_fitness: the average F1 hybrid fitness, with hybrids generated from all possible pairs of populations
//mean_F1_distance: the average F1 hybrid phenotypic distance to the optimum, with hybrids generated from all possible
//  pairs of populations
//[[Rcpp::export]]
List mutational_differences(List order, NumericMatrix mutational_library, NumericVector optimum, int l, int d, int k,int seed){
  //set random number generators with entered seed
  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_real_distribution<double> real_distribution(0.0,1.0);
  std::uniform_int_distribution<int> int_distribution(0,l-1);
  
  int n = order.size();
  int num_comparisons=0;
  int num_F1_comparisons=0;
  
  double mutational_diffs = 0.0;
  double overall_diffs = 0.0;
  double common_distances = 0.0;
  double F1_fitnesses = 0.0;
  double F1_distances = 0.0;
  
  IntegerVector genome = seq(0,l-1);
  
  //repeat for all population pairs (number of populations is n)
  for(int i = 1; i < n; i++){
    for(int ii = 0; ii < i; ii++){
      //get the fixed mutations for each of the two selected population pairs.
      //We need to subtract 1 due to indexing differences between R and C++
      //Also, we remove fixiations that happened an even number of times, as this
      //undoes changes to the original wildtype origin position of the loci
      IntegerVector order1 = strip_even(order[i]) - 1;
      IntegerVector order2 = strip_even(order[ii]) - 1;
      
      //unique common fixed mutations between the two populations
      IntegerVector inter = intersect(order1,order2);


      //calculate proportion common mutations and phenotypic distance to the optimum only considering common mutations
      if(inter.size() > 0){
        double mutational_diff = (double)inter.size()/((double)order1.size()+(double)order2.size());

        mutational_diffs = mutational_diff + mutational_diffs;
        
        NumericVector common_z = NumericVector(d);
        for(int iii = 0; iii < inter.size(); iii++){
          common_z = common_z + mutational_library(inter[iii],_);
        }
        
        common_distances += norm_vec(common_z-optimum);
      } else {
        common_distances += 2.0;
      }
      //calculate portion of the genome that is common between the two populations
      double overall_diff = ((double)order1.size()+(double)order2.size()-(double)inter.size()-(double)inter.size())/(double)l;
      overall_diffs = overall_diff + overall_diffs;
      
      //F1 hybrid generation
      for(int iii = 0; iii < 10; iii++){
        NumericVector F1_z = NumericVector(d);
        IntegerVector indexes = IntegerVector::create(-1);
        for(int iv = 0; iv < l; iv++){
          //each loci has a 50% chance to be passed on from either parent, with free recombination between loci
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
        F1_fitnesses = F1_fitnesses + fitness(optimum, F1_z,k);
        F1_distances = F1_distances + norm_vec(optimum-F1_z);
        num_F1_comparisons++;
      }
      num_comparisons++;
    }
  }
  
  
  return List::create(
    _["mean_mutational"]=(double)mutational_diffs/(double)num_comparisons,
    _["mean_overall"]=(double)overall_diffs/(double)num_comparisons,
    _["mean_common_distance"] = (double)common_distances/(double)num_comparisons,
    _["mean_F1_fitness"]=(double)F1_fitnesses/(double)num_F1_comparisons,
    _["mean_F1_distance"]=(double)F1_distances/(double)num_F1_comparisons
  );
}
