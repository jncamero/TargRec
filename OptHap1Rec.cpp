//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppArmadillo;
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
NumericMatrix OR1C(NumericMatrix geno,int start, int stop) {
// geno: chromosome, position, p11, p12, p21, p22
// NumericVector ChromEnds
   IntegerVector vec = seq((start-1),(stop-1));
 	int m = vec.size(); 
// w: length of chromosome vector, or number of chromosomes
// int w=ChromEnds.size();
// temp: holding the best haplotype genetic value and recombinaton point 
// columns: 1) Chromosome 2) recombination position 3) haplotype 1 val 
// 4) haplotype 2 val
// NumericMatrix output(w);
// holder: storing value of optimal haplotype given recombination point i
//THE FIRST RECOMBINATION POINT IS BETWEEN locus 1 and 2
      int one;    
      int two;
      NumericMatrix ll(m,6);
      NumericVector h1(4);
      NumericVector h2(4); 
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
         h1(0)=h1(0)+geno(vec(j),2);
         h1(1)=h1(1)+geno(vec(j),3);

         h1(2)=h1(2)+geno(vec(j),4);
         h1(3)=h1(3)+geno(vec(j),5);
      }

      //haplotype downstream from recombination point, to end
      for(int k = i+1; k<m; k++){
         h2(0)=h2(0)+geno(vec(k),2);
         h2(1)=h2(1)+geno(vec(k),3);
         h2(2)=h2(2)+geno(vec(k),4);
         h2(3)=h2(3)+geno(vec(k),5);
      }    

one=which_max(h1);
two=which_max(h2);

Rcout << h1(one);
Rcout << h2(two);

ll(i,0)=geno(vec(0),0);
ll(i,1)=geno(i,1);
ll(i,2)=h1(one);
ll(i,3)=h2(two);
ll(i,4)=one+1;
ll(i,5)=two+1;
}

return ll;
}

