/*
 Created on Tue Jul 15 08:44:04 2014
 
 Chris Mirabzadeh 2014
 
 @author: ChrisM
 
 program mcaimlj.c
 
 Perform Adaptive Integration Method on two Lennard-Jones particles
 
 M. Fasnacht, R.H. Swendsen, J.M. Rosenberg, "Adaptive integration
 method for Monte Carlo simulations," Phys. Rev. E. 69, 056704 (2004).
 
 
 compile using "g++ -o mcaimlj.o mcaimlj.c"
 
 runs as "./mcaimlj -T <temperature> -N <numcycles(1e6)> -dr <delta-r> \
 -rc <max lambda> -L <boxLength> -bins <number of bins>"
 
 */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <valarray>
using namespace std;

/* Trapezoidal rule on precomputed values where x and y are arrays of the same length */
double trapz(valarray<double> x, valarray<double> y, int n)
{
    double sum = 0;
    int i=1;
    while (i < n)
    {
        sum += 0.5*(x[i] - x[i - 1])*(y[i] + y[i - 1]);
        i += 1;
    }
    cout<<sum<<endl;
    return sum;
}

/*
 Generate a random number between 0 and 1
 return a uniform number in [0,1].
 */
double unifRand()
{
    return rand() / double(RAND_MAX);
}

/* Create a function for linearly spaced arrays */
valarray<double> linspace(double min, double max, int n) {
    valarray<double> newV = valarray<double>(n);
    for (unsigned i = 0; i < newV.size(); ++i) {
        newV[i] = (double)(i)*(max - min) / (double)(n - 1) + min;
    }
    return newV;
}

/* Define lennard jones potential */
double lj(double r)
{
    if (r == 0) return 0;
    r = 1.0/r;
    double r6 = pow(r,6);
    double r12 = pow(r,12);
    return 4.0*(r12 - r6);
}

/* Derivative of the LJ potential */
double ljDeriv(double r)
{
    if (r == 0) return 0;
    r = 1.0/r;
    double r7 = pow(r,7);
    double r13 = pow(r,13);
    return 24.0*(r7 - 2.0*r13);
}

/*  The main loop */
int main ( int argc, char * argv[] ) {
    // Number of steps
    int nSteps = 100000000;
    
    //Maximum change in r
    double maxDr = 0.5;
    
    // Temperature of the system
    double T = 0.2;
    
    // For computational efficiency, use reciprocal T
    const double beta = 1.0/T;
    
    // Size of the box for periodic boundary conditions
    double L = 3.0;
    
    // Maximum interaction distance
    double rc = 3.0;
    
    // x, y position of particle one
    double One[2] = {unifRand(), unifRand()};
    
    // x, y position of particle two
    double Two[2] = {unifRand(), unifRand()};
    
    // Distance between x and y positions of the particles
    double dr[2] = {(Two[0] - One[0]), (Two[1] - One[1])};
    
    // Total distance between the particles
    double rCurrent = sqrt(dr[0]*dr[0] + dr[1]*dr[1]);
    
    // Current energy between the particles
    double uCurrent = lj(rCurrent);
    
    // Track the acceptance rate
    int acc = 0, l_accept = 0;
    
    // nbins is the number of windows/segments used in AIM
    int nbins = 200;
    
    // Range of relevent areas of the potential
    const double rng[2]={1.0, rc};
    
    // all possible values of the distance between the particles
    valarray<double> l_values = linspace(rng[0],rng[1],nbins);  //tested and working
    
    // Indexes for lambda values
    int l_current = 0, l_trial = 0;
    
    // Storage for lambda values
    double l_old = 0.0, l_new = 0.0;
    
    // Our step size
    const double lstep = 1.0/((double)(nbins));
    
    // Initialize the arrays for storage of the AIM values
    valarray<double> dfest(0.0, nbins);
    valarray<double> dfnum(0.0, nbins);
    valarray<double> dfsum(0.0, nbins);
    valarray<double> dfpmf(0.0, nbins);
    
    // Our bin size
    const double binsize = nbins;
    
    // Out counter for loops
    int i = 0;
    
    // Various variables for tracking
    valarray<double> OneTrial(0.0, 2);
    valarray<double> TwoTrial(0.0, 2);
    valarray<double> drTrial(0.0, 2);
    double xTrial = 0,yTrial=0,rTrial=0,uTrial=0,deltaU=0;
    int laccept=0;
    double uOld=0,uNew=0,du=0,df=0,mcprobe=0,flatness=0;
    
    
    /* Here we parse the command line arguments */
    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
        else if (!strcmp(argv[i],"-dr")) maxDr=atof(argv[++i]);
        else if (!strcmp(argv[i],"-rc")) rc=atof(argv[++i]);
        else if (!strcmp(argv[i],"-N")) nSteps = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-L")) L = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-bins")) nbins = atoi(argv[++i]);
        else {
            fprintf(stderr,"Error.  Argument '%s' is not recognized.\n",
                    argv[i]);
            exit(-1);
        }
    }
    
    for (i=1;i<nSteps;i++){
        /*
         Generate a random number between 0 and 1
         if it is less than 0.5 then the trial particle is particle
         one otherwise the trial particle is particle two.
         */
        if (unifRand() <= 0.5)
        {
            xTrial = One[0] + maxDr*(unifRand()-0.5);
            yTrial = One[1] + maxDr*(unifRand()-0.5);
            OneTrial[0] = xTrial;
            OneTrial[1] = yTrial;
            TwoTrial[0] = Two[0];
            TwoTrial[1] = Two[1];
        }
        else
        {
            xTrial = Two[0] + maxDr*(unifRand()-0.5);
            yTrial = Two[1] + maxDr*(unifRand()-0.5);
            TwoTrial[0] = xTrial;
            TwoTrial[1] = yTrial;
            OneTrial[0] = One[0];
            OneTrial[1] = One[1];
        }
        // Apply Zero Flux Boundary Conditions to reject a move outside of the box
        if ((xTrial > L) || (xTrial < 0)) continue;
        else if ((yTrial > L) || (yTrial < 0)) continue;
        else
        {
            // Calculate the distance between X and Y
            drTrial[0] = TwoTrial[0] - OneTrial[0];
            drTrial[1] = TwoTrial[1] - OneTrial[1];
            rTrial = sqrt(pow(drTrial[0],2) + pow(drTrial[1],2));
            // Calculate the potential energy
            uTrial = lj(rTrial);
            // Calculate the difference in potential energy
            deltaU = uTrial - uCurrent;
            // If the proposed energy is lower then we accept
            if (deltaU < 0)
            {
                One[0] = OneTrial[0];
                Two[0] = TwoTrial[0];
                One[1] = OneTrial[1];
                Two[1] = TwoTrial[1];
                rCurrent = rTrial;
                uCurrent = uTrial;
                acc += 1;
            }
            // If the proposed energy is higher accept only with probability exp()
            else if (unifRand() < exp(-deltaU*beta))
            {
                One[0] = OneTrial[0];
                Two[0] = TwoTrial[0];
                One[1] = OneTrial[1];
                Two[1] = TwoTrial[1];
                rCurrent = rTrial;
                uCurrent = uTrial;
                acc += 1;
            }
        }
        
        // AIM
        // Randomly choose a direction in the index
        if (unifRand() < 0.5) l_trial = l_current + 1;
        else l_trial = l_current - 1;
        
        laccept = -1;
        
        // Track the value of l_trial so we are within the possible values
        if (l_trial < 0) laccept = 0;
        if (l_trial > (nbins-1)) laccept = 0;
        
        // Now we take our trial move
        if (laccept < 0)
        {
            // lambda old is the old configuration of the system
            l_old = l_values[l_current];
            // lambda_new is the new configuration of the system
            l_new = l_values[l_trial];
            // Calculate the potential energy at the proposed position
            uOld = lj(l_old);
            uNew = lj(l_new);
            // Get the energy difference
            du = uNew - uOld;
            // Trapezoidal rule integrate from l_current to l_trial
            df = 0.5*(l_new-l_old)*(dfest[l_trial]+dfest[l_current]);
            mcprobe = exp(-beta*(du-df));
            // If uLamdaTrial is less than probability, then accept
            if (mcprobe >= unifRand()) laccept = 1;
            else laccept = 0;
        }
        // Update values
        if (laccept > 0)
        {
            l_accept += 1;
            l_current = l_trial;
            rCurrent = l_new;
            uCurrent = uNew;
        }
        // Update Free energy estimates
        dfnum[l_current] += 1.0;
        dfsum[l_current] += ljDeriv(l_values[l_current]);
        dfest[l_current] = dfsum[l_current]/dfnum[l_current];
    }
    
    /* Print out some results */
    
    valarray<double> boltz(0.0, nbins);
    valarray<double> boltzpmf(0.0, nbins);
    double temp01=0, temp02=0;
    int ii = 0;
    while (ii < (nbins-1))
    {
        ii += 1;
        temp02 = 0.5*(l_values[ii] - l_values[ii-1])*(dfest[ii] + dfest[ii-1]);
        dfpmf[ii] = dfpmf[ii-1] + temp02;
    }
    
    for (i = 0; i < nbins; ++i)
    {
        temp01 = exp(-lj(l_values[i])*beta);
        boltzpmf[i] = -(1.0 / beta)*log(temp01);
    }
    
    fprintf (stdout,"Acceptance rate: %.5lf percent\n", (100.0*acc/float(nSteps)));
    fprintf (stdout,"Lambda Acceptance rate: %.5lf percent\n", (100.0*l_accept/float(nSteps)));
    fprintf (stdout, "Area under the curve for the Boltzman Distribution %.5lf \n", 100*trapz(l_values,boltzpmf,nbins));
    fprintf (stdout, "Area under the curve for the AIM estimate %.5lf \n", 100*trapz(l_values,dfpmf,nbins));
    
    return 0;
}