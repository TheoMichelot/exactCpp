
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <updateMove_fisher.cpp>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat localUpdate_rcpp(arma::mat allData, arma::mat aSwitches, int jorder, arma::vec par, arma::vec lambdapar,
                           double kappa, int nbState, double SDP, arma::mat map)
{
    // enable reference by "name"
    int colX = 0, colY = 1, colTime = 2, colState = 3, colHabitat = 4, colJump = 5, colBehav = 6;
    
    double T2 = allData(jorder-2,colTime);
    double X2 = allData(jorder-2,colX);
    double Y2 = allData(jorder-2,colY);
    int S2 = allData(jorder-2,colState);
    
    double T1 = allData(jorder-1,colTime);
    double X1 = allData(jorder-1,colX);
    double Y1 = allData(jorder-1,colY);
    int S1 = allData(jorder-1,colState);
    
    int H1 = findRegion(X1,Y1,map);
    
    double Tj = allData(jorder,colTime);
    double Xj = allData(jorder,colX);
    double Yj = allData(jorder,colY);
        
    // old log-likelihood
    arma::vec oldMove2 = rawMove(par,S2,T1-T2,X2,Y2,nbState);
    double oldLogLX2 = R::dnorm(X1,X2+oldMove2(0),oldMove2(2),1);
    double oldLogLY2 = R::dnorm(Y1,Y2+oldMove2(1),oldMove2(3),1);
    
    arma::vec oldMove1 = rawMove(par,S1,Tj-T1,X1,Y1,nbState);
    double oldLogLX1 = R::dnorm(Xj,X1+oldMove1(0),oldMove1(2),1);
    double oldLogLY1 = R::dnorm(Yj,Y1+oldMove1(2),oldMove1(3),1);
    
    double oldLogL = oldLogLX1 + oldLogLY1 + oldLogLX2 + oldLogLY2;
    
    // propose new location
    double Xnew = R::rnorm(X1,SDP);
    double Ynew = R::rnorm(Y1,SDP);
    double Tnew = R::runif(T2,Tj);
    int Hnew = findRegion(Xnew,Ynew,map);
    
    // effect of habitat
    // ??
    
    // new log-likelihood
    arma::vec newMove2 = rawMove(par,S2,Tnew-T2,X2,Y2,nbState);
    double newLogLX2 = R::dnorm(Xnew,X2+newMove2(0),newMove2(2),1);
    double newLogLY2 = R::dnorm(Ynew,Y2+newMove2(1),newMove2(3),1);
    
    arma::vec newMove1 = rawMove(par,S1,Tj-Tnew,Xnew,Ynew,nbState);
    double newLogLX1 = R::dnorm(Xj,Xnew+newMove1(0),newMove1(2),1);
    double newLogLY1 = R::dnorm(Yj,Ynew+newMove1(1),newMove1(3),1);
    
    double newLogL = newLogLX1 + newLogLY1 + newLogLX2 + newLogLY2;
    
    // log Hastings ratio
    double logHR = newLogL - oldLogL;
    
    // if the local update is accepted
    if(R::runif(0,1)<exp(logHR)) {
        int prev = 1; // index of switch "jorder-1" in aSwitches
        while(aSwitches(prev,colTime)<allData(jorder-1,colTime))
            prev = prev + 1;
        
        aSwitches(prev-1,colTime) = Tnew;
        aSwitches(prev-1,colX) = Xnew;
        aSwitches(prev-1,colY) = Ynew;
    }
    
    return aSwitches;
}
