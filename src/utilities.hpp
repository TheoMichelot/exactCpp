
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// extract rows 'ind' from matrix
arma::mat rows(arma::mat matrix, arma::vec ind)
{
    arma::mat subMatrix(ind.size(),matrix.n_cols);
    
    int i = 0, k = 0;
    while(i<matrix.n_rows && k<ind.size()) {
        if(ind(k)==i) {
            subMatrix.row(k) = matrix.row(i);
            k = k+1;
        }
        i = i+1;
    }
    
    return subMatrix;
}

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
    arma::vec par(wpar.size());
    
    for(int i=0 ; i<2*nbState ; i++)
        par(i) = wpar(i); // m
    for(int i=2*nbState ; i<4*nbState ; i++)
        par(i) = exp(wpar(i)); // exp(b) and exp(v)
    
    return(par); // wpar = c(m,exp(b),exp(v))
}

//' Truncated beta
//' 
//' Random generation for the truncated Beta distribution
//' 
//' @param n Number of values to generate
//' @param shape1 Non-negative parameter
//' @param shape2 Non-negative parameter
//' @param low Lower limit of the truncated distribution
// [[Rcpp::export]]
double rtruncbeta_rcpp(double shape1, double shape2, double low)
{
    double pmin = R::pbeta(low,shape1,shape2,1,0);
    double p = pmin + R::runif(0,1)*(1-pmin);
    
    return R::qbeta(p,shape1,shape2,1,0);
}

//' Returns an array of dimensions (nbState,nbState,nbObs) filled with number 
//' of occurrences of each element in (a1,a2,a3).
//' 
//' @param a1 Vector of values for 1st dimension
//' @param a2 Vector of values for 2nd dimension
//' @param a3 Vector of values for 3rd dimension
//' @param nbState Number of states (number of rows and columns of the array)
//' @param nbObs Number of observations (number of layers of the array)
// [[Rcpp::export]]
arma::cube minitabGeneric_rcpp(arma::vec a1, arma::vec a2, arma::vec a3,
                               int nbState, int nbObs)
{
    arma::cube res(nbState,nbState,nbObs);
    res.zeros();
    
    for(int i=0 ; i<a1.size() ; i++)
        res(a1(i)-1,a2(i)-1,a3(i)-1) = res(a1(i)-1,a2(i)-1,a3(i)-1)+1;
    
    return res;
}
