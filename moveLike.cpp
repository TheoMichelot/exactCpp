
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

#include <utilities.hpp>
#include <rawMove.cpp>
#include <findRegion.cpp>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double moveLike_rcpp(arma::mat subData, arma::uvec indObs, arma::uvec indSwitch, arma::vec par,
                     arma::mat aSwitches, int nbState)
{
    // enable reference by "name"
    int colX = 0, colY = 1, colTime = 2, colState = 3, colHabitat = 4, colJump = 5, colBehav = 6;
    
    arma::mat subObs = subData.rows(indObs);
    int len = subObs.n_rows;
    
    double Tbeg = subObs(0,colTime);
    double Tend = subObs(len-1,colTime);
    
    //============//
    // Likelihood //
    //============//
    
    // 1. Likelihood with potential switches
    arma::vec newLikeX(len-1), newLikeY(len-1);
    
    // compute distribution of observations, i.e. likelihood
    for(int i=0 ; i<len-1 ; i++) { 
        int which = indObs(i+1)-1; // index of previous switch
        double deltaT = subData(which+1,colTime) - subData(which,colTime); // interval between switch and obs
        int state = subData(which,colState); // state at previous switch
        
        arma::vec move = rawMove(par,state,deltaT,subData(which,colX),subData(which,colY),nbState);
        double emx = move(0);
        double emy = move(1);
        double sdx = move(2);
        double sdy = move(3);
        
        newLikeX(i) = R::dnorm(subData(indObs(i+1),colX),
                 subData(which,colX)+emx, sdx, 1);
        
        newLikeY(i) = R::dnorm(subData(indObs(i+1),colY),
                 subData(which,colY)+emy, sdy, 1);
        
    }
    
    // 2. Likelihood with actual switches
    
    // count actual switches happening between Tbeg and Tend
    int nbActual = 0;
    for(int i=0 ; i<aSwitches.n_rows ; i++) {
        if(aSwitches(i,colTime)>Tbeg && aSwitches(i,colTime)<Tend)
            nbActual = nbActual+1;
    }
    
    // indices of those switches
    arma::uvec whichActual(nbActual);
    int count = 0;
    for(int i=0 ; i<aSwitches.n_rows ; i++) {
        if(aSwitches(i,colTime)>Tbeg && aSwitches(i,colTime)<Tend) {
            whichActual(count) = i;
            count = count+1;
        }
    }
    
    // actual switches and observations between Tbeg and Tend
    arma::mat aSubData = join_cols(aSwitches.rows(whichActual),subObs);
    
    // ranks of data by time
    arma::uvec aRanks = myrank(aSubData.col(colTime));
    // indices of actual switches among observations
    if(nbActual>0)
        arma::uvec aIndSwitch = aRanks(arma::span(0,nbActual-1));
    // indices of observations among actual switches
    arma::uvec aIndObs = aRanks(arma::span(nbActual,nbActual+len-1));
    
    // order data in time
    arma::mat aSubDataCopy = aSubData;
    arma::uvec orders = sort_index(aSubData.col(colTime));
    for(int i=0 ; i<nbActual+len ; i++)
        aSubData.row(i) = aSubDataCopy.row(orders(i));
    
    arma::vec oldLikeX(len-1), oldLikeY(len-1);

    for(int i=0 ; i<len-1 ; i++) {
        int which = aIndObs(i+1)-1; // index of actual switch
        double deltaT = aSubData(which+1,colTime) - aSubData(which,colTime); // time interval between switch and obs
        int state = aSubData(which,colState); // state of actual switch
        
        arma::vec move = rawMove(par,state,deltaT,aSubData(which,colX),aSubData(which,colY),nbState);
        double emx = move(0);
        double emy = move(1);
        double sdx = move(2);
        double sdy = move(3);
        
        oldLikeX(i) = R::dnorm(aSubData(which+1,colX),
                 aSubData(which,colX)+emx, sdx, 1);
        
        oldLikeY(i) = R::dnorm(aSubData(which+1,colY),
                 aSubData(which,colY)+emy, sdy, 1);
    }

    // 3. Deduce the Hastings ratio
    arma::vec newLike = newLikeX + newLikeY;
    arma::vec oldLike = oldLikeX + oldLikeY;
    double HR = exp(sum(newLike)-sum(oldLike));

    return HR;
}