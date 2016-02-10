
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>

#include <rawMove.cpp>
#include <findRegion.cpp>

using namespace std;
using namespace Rcpp;

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
    while(s<rand) {
        k = k+1;
        s = s + probs(k);
    }
    
    return k;
}

// [[Rcpp::export]]
List simMove_rcpp(arma::mat subObs, arma::vec par, double kappa, arma::vec lambdapar, int nbState, arma::mat map)
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
    
    cout << "1" << endl;
    
    // initialize potential switches
    arma::mat switches = arma::mat(nbSwitch,7).fill(NA_REAL);
    for(int i=0 ; i<nbSwitch ; i++) {
        switches(i,colTime) = R::runif(Tbeg,Tend); // Poisson process -> uniformly distributed
        switches(i,colBehav) = 0;
    }
    
    // simulated switches are on top, observations at the bottom
    arma::mat subData = join_cols(switches,subObs); // equivalent to rbind()
    
    // ranks of data by time
    arma::uvec ranks = myrank(subData.col(colTime));
    // indices of potential switches among observations
    arma::uvec indSwitch = ranks(arma::span(0,nbSwitch-1));
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
            arma::vec probs(nbState,0);
            probs(subData(t,colHabitat)) = A(subData(t-1,colState)-1,subData(t,colHabitat)-1);
            
            bool jumpNow = (R::runif(0,1)<sum(probs));
            
            if(jumpNow) {
                // if jump, pick new state
                subData(t,colState) = mysample(probs);
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
        }
        
        t = t+1;
    }
    
    List res(4);
    res[0] = subData;
    res[1] = indObs;
    res[2] = indSwitch;
    res[3] = bk;
    
    return(res);
}
