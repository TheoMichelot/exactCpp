#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Find region
//' 
//' Identify the region of (x,y) on 'map'.
//' 
//' @param x 'x' coordinate
//' @param y 'y' coordinate
//' @param map A matrix of habitats
int findRegion(double x, double y, arma::mat map)
{
    // floor to have indices between 0 and (dim-1)
    int xi = floor(x);
    int yi = floor(y);

    // keep coordinates within map limits
    if(xi<0) 
        xi = 0;
    else if(xi>=map.n_cols) 
        xi = map.n_cols-1;

    if(yi<0)
        yi = 0;
    else if(yi>=map.n_rows)
        yi = map.n_rows-1;

    return map(xi,yi);
}
