#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Raw movement
//' 
//' Returns mean change in location, and standard deviation of the animal's next locations, 
//' from the movement parameters and the previous location.
//' 
//' @param par Vector of movement parameters
//' @param state State of underlying process
//' @param deltaT Time interval between current and next location
//' @param x X-coordinate of current location
//' @param y Y-coordinate of current location
//' @param nbState Number of states
//' @param mty Movement type (1 if Brownian motion, 2 if OU)
arma::vec rawMove(arma::vec par, int state, double deltaT, double x, double y, int nbState,
                  int mty)
{
    // unpack the vector of parameters
    arma::vec mpar = par(arma::span(0,2*nbState-1));
    arma::vec bpar = par(arma::span(2*nbState,3*nbState-1));
    arma::vec vpar = par(arma::span(3*nbState,4*nbState-1));

    arma::vec res(4);
    
    if(mty==1) { // Brownian motion
        double v = vpar(state-1);
        
        res(0) = 0;
        res(1) = 0;
        res(2) = v;
        res(3) = v;
    } else { // Ornstein-Uhlenbeck
        double mux = mpar(2*state-2);
        double muy = mpar(2*state-1);
        double b = -bpar(state-1);
        double v = vpar(state-1);
        double phi = exp(b*deltaT);
        double sd = sqrt(v*(1-phi*phi));
        
        res(0) = (1-phi)*(mux-x); // x mean (new mean - old location)
        res(1) = (1-phi)*(muy-y); // y mean (new mean - old location)
        res(2) = sd; // x sd
        res(3) = sd; // y sd        
    }

    return res;
}
