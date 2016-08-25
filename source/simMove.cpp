
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <utilities.hpp>
#include <rawMove.cpp>
#include <findRegion.cpp>

#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
List simMove_rcpp(arma::mat subObs, arma::vec par, double kappa, arma::vec lambdapar, int nbState, arma::mat map,
                  bool adapt)
{
    // enable reference by "name"
    int colX = 0, colY = 1, colTime = 2, colState = 3, colHabitat = 4, colJump = 5, colBehav = 6;
    
    int len = subObs.n_rows;
    
    double Tbeg = subObs(0,colTime);
    double Tend = subObs(len-1,colTime);
    
    //=============================//
    // Simulate potential switches //
    //=============================//
    // number of potential switches
    int nbSwitch = R::rpois((Tend-Tbeg)*kappa);
    
    arma::mat switches;
    arma::mat subData;
    if(nbSwitch>0) {
        // initialize potential switches
        switches = arma::mat(nbSwitch,7).fill(NA_REAL);
        for(int i=0 ; i<nbSwitch ; i++) {
            switches(i,colTime) = R::runif(Tbeg,Tend); // Poisson process -> uniformly distributed
            switches(i,colBehav) = 0;
        }
        
        // simulated switches are on top, observations at the bottom
        subData = join_cols(switches,subObs); // equivalent to rbind()
    }
    else
        subData = subObs; // if nbSwitch=0
    
    // ranks of data by time
    arma::uvec ranks = myrank(subData.col(colTime));
    
    arma::uvec indSwitch;
    // indices of potential switches among observations
    if(nbSwitch>0)
         indSwitch = ranks(arma::span(0,nbSwitch-1));

    // indices of observations among potential switches
    arma::uvec indObs = ranks(arma::span(nbSwitch,nbSwitch+len-1));

    // order data in time
    arma::mat subDataCopy = subData;
    arma::uvec orders = sort_index(subData.col(colTime));
    for(int i=0 ; i<nbSwitch+len ; i++)
        subData.row(i) = subDataCopy.row(orders(i));
    
    //================================//
    // Simulate movement and switches //
    //================================//
    bool bk = false;
    int t = 1;
    while(!bk && t<nbSwitch+len) {
        
        if(R_IsNA(subData(t,colHabitat))) { // if potential switch
            int state = subData(t-1,colState); // state prior to t
            double deltaT = subData(t,colTime) - subData(t-1,colTime); // time interval between t and t-1
            
            double Xfrom = subData(t-1,colX);
            double Yfrom = subData(t-1,colY);
            
            // compute distribution of location at time t
            arma::vec move = rawMove(par,state,deltaT,Xfrom,Yfrom,nbState);
            double emx = move(0);
            double emy = move(1);
            double sdx = move(2);
            double sdy = move(3);
            
            // simulate location at time t
            subData(t,colX) = R::rnorm(Xfrom+emx,sdx);
            subData(t,colY) = R::rnorm(Yfrom+emy,sdy);
            
            // map back to known habitat
            subData(t,colHabitat) = findRegion(subData(t,colX),subData(t,colY),map);
            
            // matrix of switching rates
            arma::mat A(nbState,nbState);
            int k = 0;
            for(int i=0 ; i<nbState ; i++) {
                for(int j=0 ; j<nbState ; j++) {
                    if(i==j)
                        A(i,j) = 0; // diagonal is zero
                    else {
                        A(i,j) = lambdapar(k); // switching rate i -> j
                        k = k+1;
                    }
                }
            }
            
            // probabilities of actual switch
            arma::rowvec probs(nbState);
            if(adapt) {
                probs.zeros();
                probs(subData(t,colHabitat)-1) = A(subData(t-1,colState)-1,subData(t,colHabitat)-1);
            } 
            else
                probs = A.row(subData(t-1,colState)-1);
            
            probs = probs/kappa;
            
            bool jumpNow = (R::runif(0,1)<sum(probs));
            
            if(jumpNow) {
                // if jump, pick new state
                subData(t,colState) = mysample(arma::conv_to< arma::vec >::from(probs));
                subData(t,colJump) = 1;
            }
            else {
                subData(t,colState) = subData(t-1,colState);
                subData(t,colJump) = 0;
            }
        }
        else { // if observation
            bool check1 = ((t==subData.n_rows-1) || (subData(t,colBehav)==1));
            bool check2 = (subData(t-1,colState) != subData(t,colState));
            
            if(check1 && check2) // if state at Tend cannot be changed
                bk = true;
            else
                subData(t,colState) = subData(t-1,colState);
        }
        
        t = t+1;
    }

    List res(4);
    res[0] = subData;
    res[1] = indObs+1; // +1 for compatibility with R
    res[2] = indSwitch+1; // idem
    res[3] = bk;
    
    return(res);
}
