
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
    
    // only works if the probs sum to 1, unlike R's sample
    if(sum(probs)!=1)
        probs = probs/sum(probs);
    
    while(s<rand && k<probs.size()) {
        s = s + probs(k);
        k = k+1;
    }
    
    return k;
}

//' Natural to working
//' 
//' Transforms the parameters from the natural scale to the working scale.
//' 
//' @param par Vector of parameters, in the order: m, b, v.
//' @param nbState Number of states.
arma::vec n2w(arma::vec par, int nbState)
{
    arma::vec wpar(par.size());
    
    for(int i=0 ; i<2*nbState ; i++)
        wpar(i) = par(i); // m
    for(int i=2*nbState ; i<4*nbState ; i++)
        wpar(i) = log(par(i)); // log(b) and log(v)
    
    return(wpar); // wpar = c(m,log(b),log(v))
}

//' Working to natural
//' 
//' Transforms the parameters from the working scale to the natural scale.
//' 
//' @param wpar Vector of parameters, in the order: m, b, v.
//' @param nbState Number of states.
arma::vec w2n(arma::vec wpar, int nbState)
{
    arma::vec par(par.size());
    
    for(int i=0 ; i<2*nbState ; i++)
        par(i) = wpar(i); // m
    for(int i=2*nbState ; i<4*nbState ; i++)
        par(i) = exp(wpar(i)); // exp(b) and exp(v)
    
    return(par); // wpar = c(m,exp(b),exp(v))
}
