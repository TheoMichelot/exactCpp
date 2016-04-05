
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <updateMove_fisher.cpp>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec updatePar_rcpp(arma::mat allData, arma::vec par, arma::vec priorMean, arma::vec priorSD, 
                         arma::vec proposalSD, int nbState, bool mHomog, bool bHomog, bool vHomog, 
                         arma::mat obs)
{
    // enable reference by "name"
    int colX = 0, colY = 1, colTime = 2, colState = 3, colHabitat = 4, colJump = 5, colBehav = 6;
    
    // pick out points that are informative about movement
    int indMin = 0;
    while(allData(indMin,colTime)<=obs(0,colTime))
        indMin = indMin + 1;
    
    arma::mat allData2 = allData.rows(arma::span(indMin,allData.n_rows-1));
    int l = allData2.n_rows;
    
    // calculate changes in position and time
    arma::vec preX(l), preY(l), dX(l), dY(l);
    preX(0) = allData(indMin-1,colX);
    dX(0) = allData2(0,colX) - preX(0);
    preY(0) = allData(indMin-1,colY);
    dY(0) = allData2(0,colY) - preY(0);
    for(int i=1 ; i<l ; i++) {
        preX(i) = allData2(i-1,colX);
        dX(i) = allData2(i,colX) - preX(i);
        preY(i) = allData2(i-1,colY);
        dY(i) = allData2(i,colY) - preY(i);
    }
    
    // pick new movement parameters from proposal
    List moveStep = updateMove_rcpp(par,priorMean,priorSD,proposalSD,nbState,mHomog,bHomog,vHomog);
    arma::vec newPar = moveStep[0];
    double oldLogPrior = moveStep[1];
    double newLogPrior = moveStep[2];
    
    arma::vec oldMove(l), newMove(l);
    arma::vec oldLogLX(l), oldLogLY(l), newLogLX(l), newLogLY(l);
    
    for(int i=0 ; i<l ; i++) {
        int which = i-1; // index of previous switch
        double deltaT; // interval between switch and obs
        int state; // state at previous switch
        
        if(i>0) {
            deltaT = allData2(which+1,colTime) - allData2(which,colTime);
            state = allData2(which,colState);
        }
        else {
            deltaT = allData2(which+1,colTime) - allData(indMin-1,colTime);
            state = allData(indMin-1,colState);
        }
        
        oldMove = rawMove(par,state,deltaT,preX(i),preY(i),nbState);
        newMove = rawMove(newPar,state,deltaT,preX(i),preY(i),nbState);
        
        // oldMove = (emx,emy,sdx,sdy)
        oldLogLX(i) = R::dnorm(dX(i), oldMove(0), oldMove(2), 1);
        oldLogLY(i) = R::dnorm(dY(i), oldMove(1), oldMove(3), 1);
        newLogLX(i) = R::dnorm(dX(i), newMove(0), newMove(2), 1);
        newLogLY(i) = R::dnorm(dY(i), newMove(1), newMove(3), 1);
    }
    
    // log Hastings ratio
    double logHR = newLogPrior - oldLogPrior + sum(newLogLX) + sum(newLogLY) - sum(oldLogLX) - sum(oldLogLY);
    
    if(R::runif(0,1)<exp(logHR)) {
        // accept movement parameters
        return newPar;
    }
    else
        return par;
}