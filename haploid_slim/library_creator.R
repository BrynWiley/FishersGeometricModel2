#A simple file to create mutation libraries for use by Haploid_Finite_Poopulation_Creation.slim
#This process needs to performed before SLiM simulations, 
#otherwise overlapping SLiM simulations will create malformed mutaiton libraries
library(tidyverse)

set.seed(12345)

reference_tibble <- tibble(lambda=numeric(),
                           l=numeric(),
                           id=numeric())

library_id=1;
for(lambda in c(0.08,0.12,016)){
  for(l in c(100,500,1000)){
    for(i in 1:25){
      library <-data.frame(sphere=numeric(),
                           traitValue1=numeric(),
                           traitValue2=numeric(),
                           traitValue3=numeric(),
                           traitValue4=numeric(),
                           traitValue5=numeric())
      for(ii in 1:l){
        library <-rbind(library,c(1,rnorm(5)*rexp(1,rate=1/lambda)))
      }
      write_csv(library,paste("mutationLibrary",library_id,".csv",sep=""))
      reference_tibble <- add_row(reference_tibble,lambda=lambda,l=l,id=library_id)
      library_id=library_id+1
    }
  }
}
