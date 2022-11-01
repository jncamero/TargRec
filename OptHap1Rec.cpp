//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
NumericMatrix OR1C(NumericMatrix geno) {
// geno: chromosome, position, p11,p12,p21,p22
// NumericVector ChromEnds
 	int m=geno.nrow();
// w: length of chromosome vector, or number of chromosomes
// int w=ChromEnds.size();
// temp: holding the best haplotype genetic value and recombinaton point 
// columns: 1) Chromosome 2) recombination position 3) haplotype 1 val 
// 4) haplotype 2 val
// NumericMatrix output(w);
// holder: storing value of optimal haplotype given recombination point i
//THE FIRST RECOMBINATION POINT IS BETWEEN locus 1 and 2
    NumericMatrix holder(m,4);
   int one;    
   int two;
   NumericMatrix ll(m,3);
   NumericVector h1(4);
   NumericVector h2(4); 
   NumericVector h3(m);
// i: recombination point
    for (int i = 0; i <m; i++) {   
      h1(0)=0;
      h1(1)=0;
      h1(2)=0;
      h1(3)=0;
      
      h2(0)=0;
      h2(1)=0;
      h2(2)=0;
      h2(3)=0;
      // haplotype upstream, to recombination point
      for(int j = 0; j <= i; j++){
         h1(0)=h1(0)+geno(j,1);
         h1(1)=h1(1)+geno(j,2);
         h1(2)=h1(2)+geno(j,3);
         h1(3)=h1(3)+geno(j,4);
      }

      //haplotype downstream from recombination point, to end
      for(int k = i+1; k<m; k++){
         h2(0)=h2(0)+geno(k,1);
         h2(1)=h2(1)+geno(k,2);
         h2(2)=h2(2)+geno(k,3);
         h2(3)=h2(3)+geno(k,4);
      }    

one=which_max(h1);
two=which_max(h2);
ll(i,0)=h1(one);
ll(i,1)=h2(two);
ll(i,2)=ll(i,0)+ll(i,1);
}

return ll;
}

