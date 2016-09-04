
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "utilities.hpp"

// [[Rcpp::export]]
arma::vec updateRate_rcpp(arma::mat allData, arma::vec indSwitch, double kappa,
                          double shape1, double shape2)
{
    // enable reference by "name"
    int colX = 0, colY = 1, colTime = 2, colState = 3, colHabitat = 4, colJump = 5, colBehav = 6;
    
    arma::vec rates(6);
    
    arma::vec s1 = rows(allData,indSwitch-1).col(colState);
    arma::vec s2 = rows(allData,indSwitch).col(colState);
    arma::vec h = rows(allData,indSwitch).col(colHabitat);
    
    arma::cube counts =  minitabGeneric_rcpp(s1,s2,h,3,3);
    
    rates(0) = kappa*rtruncbeta_rcpp(shape1+counts(0,1,1),shape2+counts(0,0,1),0); // 1 -> 2
    rates(1) = kappa*rtruncbeta_rcpp(shape1+counts(0,2,2),shape2+counts(0,0,2),0); // 1 -> 3
    rates(2) = kappa*rtruncbeta_rcpp(shape1+counts(1,0,0),shape2+counts(1,1,0),0); // 2 -> 1
    rates(3) = kappa*rtruncbeta_rcpp(shape1+counts(1,2,2),shape2+counts(1,1,2),0); // 2 -> 3
    rates(4) = kappa*rtruncbeta_rcpp(shape1+counts(2,0,0),shape2+counts(2,2,0),0); // 3 -> 1
    rates(5) = kappa*rtruncbeta_rcpp(shape1+counts(2,1,1),shape2+counts(2,2,1),0); // 3 -> 2

    return rates;
}