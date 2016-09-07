
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <utilities.hpp>
#include <rawMove.cpp>
#include <findRegion.cpp>

using namespace Rcpp;

// Multivariate normal generator
arma::mat mvrnormArma(arma::vec mu, arma::mat sigma) 
{
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(1, ncols);
    return arma::repmat(mu, 1, 1).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
List updateMove_rcpp(arma::vec par, arma::vec priorMean, arma::vec priorSD, arma::vec proposalSD,
                     int nbState, bool mHomog, bool bHomog, bool vHomog, arma::vec mty)
{
    arma::vec wpar = n2w(par,nbState);
    // unpack the vector of parameters
    arma::vec m = wpar(arma::span(0,2*nbState-1));
    arma::vec b = wpar(arma::span(2*nbState,3*nbState-1));
    arma::vec v = wpar(arma::span(3*nbState,4*nbState-1));
    
    // unpack prior
    arma::vec mPriorMean = priorMean(arma::span(0,2*nbState-1));
    arma::vec bPriorMean = priorMean(arma::span(2*nbState,3*nbState-1));
    arma::vec vPriorMean = priorMean(arma::span(3*nbState,4*nbState-1));
    arma::vec mPriorSD = priorSD(arma::span(0,2*nbState-1));
    arma::vec bPriorSD = priorSD(arma::span(2*nbState,3*nbState-1));
    arma::vec vPriorSD = priorSD(arma::span(3*nbState,4*nbState-1));
    
    // unpack proposal
    arma::vec mProposalSD = proposalSD(arma::span(0,2*nbState-1));
    arma::vec bProposalSD = proposalSD(arma::span(2*nbState,3*nbState-1));
    arma::vec vProposalSD = proposalSD(arma::span(3*nbState,4*nbState-1));
    
    arma::vec mprime = m, bprime = b, vprime = v;
    
    double oldLogPrior = 0, newLogPrior = 0;
    
    // allow for homogeneous and non-homogeneous cases
    if(mHomog) {
        double xmove = R::rnorm(0,mProposalSD(0)); // change in x
        double ymove = R::rnorm(0,mProposalSD(1)); // change in y
        for(int i=0 ; i<m.size()/2 ; i++) {
            // apply same change to all states
            mprime(2*i) = mprime(2*i) + xmove;
            mprime(2*i+1) = mprime(2*i+1) + ymove;
        }
        
        // log-prior with old parameter values
        oldLogPrior = oldLogPrior + nbState*R::dnorm(m(0),mPriorMean(0),mPriorSD(0),1) + 
            nbState*R::dnorm(m(1),mPriorMean(1),mPriorSD(1),1);
        
        // log-prior with new parameter values
        newLogPrior = newLogPrior + nbState*R::dnorm(mprime(0),mPriorMean(0),mPriorSD(0),1) + 
            nbState*R::dnorm(mprime(1),mPriorMean(1),mPriorSD(1),1);
    } else {
        for(int i=0 ; i<m.size() ; i++) {
            if(mty(floor(i/2))==2) { // only if OU
                mprime(i) = mprime(i) + R::rnorm(0,mProposalSD(i)); // different change for each state
                
                oldLogPrior = oldLogPrior + R::dnorm(m(i),mPriorMean(i),mPriorSD(i),1);
                newLogPrior = newLogPrior + R::dnorm(mprime(i),mPriorMean(i),mPriorSD(i),1);                
            }
        }
    }
    
    // if(bHomog) {
    //     double bmove = R::rnorm(0,bProposalSD(0));
    //     for(int i=0 ; i<b.size() ; i++)
    //         bprime(i) = bprime(i) + bmove;
    //     
    //     oldLogPrior = oldLogPrior + nbState*R::dnorm(b(0),bPriorMean(0),bPriorSD(0),1);
    //     newLogPrior = newLogPrior + nbState*R::dnorm(bprime(0),bPriorMean(0),bPriorSD(0),1);
    // } else {
    //     for(int i=0 ; i<b.size() ; i++) {
    //         if(mty(i)==2) { // only if OU
    //             bprime(i) = bprime(i) + R::rnorm(0,bProposalSD(i));
    //             
    //             oldLogPrior = oldLogPrior + R::dnorm(b(i),bPriorMean(i),bPriorSD(i),1);
    //             newLogPrior = newLogPrior + R::dnorm(bprime(i),bPriorMean(i),bPriorSD(i),1);
    //         }
    //     }
    // }
    // 
    // if(vHomog) {
    //     double vmove = R::rnorm(0,vProposalSD(0));
    //     for(int i=0 ; i<v.size() ; i++)
    //         vprime(i) = vprime(i) + vmove;
    //     
    //     oldLogPrior = oldLogPrior + nbState*R::dnorm(v(0),vPriorMean(0),vPriorSD(0),1);
    //     newLogPrior = newLogPrior + nbState*R::dnorm(vprime(0),vPriorMean(0),vPriorSD(0),1);
    // } else {
    //     for(int i=0 ; i<v.size() ; i++) {
    //         vprime(i) = vprime(i) + R::rnorm(0,vProposalSD(i));
    //         
    //         oldLogPrior = oldLogPrior + R::dnorm(v(i),vPriorMean(i),vPriorSD(i),1);
    //         newLogPrior = newLogPrior + R::dnorm(vprime(i),vPriorMean(i),vPriorSD(i),1);
    //     }
    // }
    
    arma::mat sigma(2,2);
    sigma(0,0) = 1;
    sigma(0,1) = -0.4;
    sigma(1,0) = -0.4;
    sigma(1,1) = 1;
    arma::vec mu(2);
    mu.zeros();
    
    for(int i=0 ; i<b.size() ; i++) {
        arma::mat update = mvrnormArma(mu, sigma);
        bprime(i) = bprime(i) + update(0);
        vprime(i) = vprime(i) + update(1);
        
        oldLogPrior = oldLogPrior + R::dnorm(b(i),bPriorMean(i),bPriorSD(i),1);
        newLogPrior = newLogPrior + R::dnorm(bprime(i),bPriorMean(i),bPriorSD(i),1);
        oldLogPrior = oldLogPrior + R::dnorm(v(i),vPriorMean(i),vPriorSD(i),1);
        newLogPrior = newLogPrior + R::dnorm(vprime(i),vPriorMean(i),vPriorSD(i),1);
    }
    
    // newPar = w2n(c(mprime,bprime,vprime),nbState)
    arma::vec newPar(m.size()+b.size()+v.size());
    for(int i=0 ; i<2*nbState ; i++)
        newPar(i) = mprime(i);
    for(int i=0 ; i<nbState ; i++)
        newPar(2*nbState+i) = bprime(i);
    for(int i=0 ; i<nbState ; i++)
        newPar(3*nbState+i) = vprime(i);
    
    newPar = w2n(newPar,nbState);
    
    // list to be returned
    List res(3);
    res[0] = newPar;
    res[1] = oldLogPrior;
    res[2] = newLogPrior;
    
    return res;
}
