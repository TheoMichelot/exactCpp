
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Analoguous to R's rank function
arma::uvec myrank(arma::vec x) 
{
    int n = x.size();
    arma::uvec res(n);
    arma::vec ordx = x(sort_index(x));
    int j;
    
    for(int i=0 ; i<n ; i++) {
        j = 0;
        
        while(ordx(j)!=x(i))
            j = j+1;
        
        res(i) = j;
    }
    
    return res;
}

// pick integer in 1:length(probs) according to probabilities probs
int mysample(arma::vec probs) 
{
    int k = 0;
    double s = 0;
    double rand = R::runif(0,1);
    while(s<rand && k<probs.size()) {
        s = s + probs(k);
        k = k+1;
    }
    
    return k;
}