#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <stdlib.h> 
#include <chrono>
#include <random>
using namespace Rcpp;


//This file contains functions used by haploid_monomorphic.r for
//simulating haploid population adaptation using Fisher's Geometric Model

//A simple function to generate the 2-norm of a vector
//[[Rcpp::export]]
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
//[[Rcpp::export]]
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


//A function to get the probability of fixation of a mutation
//Inputs:
//N: The population size
//s: the selection coefficient of the mutation
//Output: The probability of fixation with selection coefficient s in a population of size N
//[[Rcpp::export]]
double prob_fix(int N, double s){
  return((1-exp(-2*s))/(1-exp(-2*N*s)));
}

//A function to get the selection coefficient of a mutation
//Inputs:
//z: a phenotype vector of the individuals in the population this mutation arose in
//delta_z: a vector representing the shift in phenotype caused by the mutation
//optimum: the phenotypic optimum of the population
//k: the k value for fitness calculation
//Output: the selection coefficient of a mutation with parameters entered
//[[Rcpp::export]]
double s_coef(NumericVector z, NumericVector delta_z, NumericVector optimum, int k){
  return(fitness(optimum,z+delta_z,k)/fitness(optimum,z,k)-1);
}

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

//A function that simulates a monomorphic haploid population adapting to a given environmental optimum,
//with a pre-defined finite set of possible mutations, using Fisher's geometric model
//Inputs:
//mutation_library: a matrix representing the possible mutations available to the population.
//  Each row represents a mutational vector associated with a given mutation.
//optimum: a vector representing the phenotypic optimum
//k: the k value for fitness calculations
//N: the population size
//l: the number of loci contributing to adaptation in this simulation, also the 
//  number of available mutations
//n: the number of trait dimensions
//order: The id's of pre-existing mutations fixed in the population before simulation, in order of fixation.
//  Enter an empty vector if none have occured
//fitnesses: the previous fitness values for the population before simulation, after each previous fixation event
//  Enter an empty vector if no mutations have fixed before this simultation
//ss: the previous selection coefficients of mutations fixed in the population before simulation, in order of mutation fixation
//  Enter an empty vector if no mutations have fixed before this simultation
//distances: the previous distances from the optimum for this population before simultation, after each previous fixation event
//  Enter an empty vector if no mutations have fixed before this simulation
//endpoint: the number of fixation events needed to have occured before this simulation finishes
//seed: an integer random number generator seed
//Output:
//A list containing four components:
//fitnesses: the fitness of the individuals in the population after each fixation event in the simulation
//s: the selection coefficients of each fixed mutation, in order of fixation
//order: the id's of fixed mutations, in order of fixation. Refers to the column number
//  of mutation vectors in the mutation library associated with the fixed mutations.
//distances: the phenotypic distance to the optimum of of the individuals in the population after each fixation event
//[[Rcpp::export]]
List haploid_sim(NumericMatrix mutation_library, NumericVector optimum, int k, int N,int l, int n, IntegerVector order,NumericVector fitnesses, NumericVector ss,NumericVector distances, int endpoint,int seed){
  //establish random number generators
  std::default_random_engine generator;
  generator.seed(seed);
  std::uniform_real_distribution<double> real_distribution(0.0,1.0);
  std::uniform_int_distribution<int> int_distribution(0,l-1);
  
  //adjustment for 1 vs 0 first indexes between R and C++
  order = order-1;
  int total_fixations = order.size();
  
  NumericVector z = NumericVector(n);
  for(int i=0; i < order.size(); i++){
    z = z + mutation_library(order[i],_);
  }
  
  while(total_fixations < endpoint){
    //randomly select a loci to mutate
    int position = int_distribution(generator);
    //if the loci is already mutate, set deltaz to -deltaz
    NumericVector delta_z = mutation_library(position,_);
    if(std::count(order.begin(),order.end(),position) % 2 == 1){
      delta_z = delta_z * -1;
    }
    //extract fixation probability for mutation, use to determine if it fixes
    double s = s_coef(z,delta_z,optimum,k);
    double fixation_prob = prob_fix(N,s);
    double random_number = real_distribution(generator);
    //if it fixes, record data
    if(random_number < fixation_prob){
      z=z+delta_z;
      order.push_back(position);
      ss.push_back(s);
      fitnesses.push_back(fitness(optimum, z, k));
      distances.push_back(norm_vec(optimum-z));
      total_fixations++;
    }
  }
  
  return List::create(
    _["fitnesses"] = fitnesses,
    _["s"] = ss,
    _["order"]=order + 1,
    _["distances"]=distances
  );
  
}



