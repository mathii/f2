/* Helper functions for inference */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/************************************************************************
* return a single sample from the integers 0:n-1 
* with probabilities given by probs (need not sum to 1)
************************************************************************/

int sample(int n, double *probs){
  int i;
  double cum_prob [n];
  double runif, p;

  p=0.0;
  for(i=0;i<n;++i){
    cum_prob[i]=p+probs[i];
    p=cum_prob[i];
  }

  runif=cum_prob[n-1]*rand()/(RAND_MAX+1.0);

  for(i=0;i<n;++i)
    if(runif<cum_prob[i])
      break;

  return i;
}

/************************************************************************
* MCMC sampler for density. Given a n x k vector of likelihoods for 
* n observations with a discrete density at k points, sample from the
* distribution of the density. alpha is a fixed hyperparameter 
* representing the confidence in the prior. 
*
* fmat: likelihood matrix, but as a column-stacked vector
* fmat_dim: dimensions of the likelihood matrix
* prior: vector of prior probabilities
* mcmc_params: number of simulations, burn in, thin, alpha
* RNG_seed: seed
* all_counts: zero vector of length k 
************************************************************************/

void mcmc_density_sampler( double *fmat, int *fmat_dim, double *prior, int *mcmc_params, int *RNG_seed, int *all_counts ){
  int total_sims,n_sims, burn_in, thin, alpha, k, n;
  int i,j,sim, ti;
  double bar=1;

  printf("Sampling............\n" );
  
  srand(*RNG_seed);
  
  n_sims=mcmc_params[0];
  burn_in=mcmc_params[1];
  thin=mcmc_params[2];
  alpha=mcmc_params[3];
  total_sims=n_sims+burn_in;

  n=fmat_dim[0];
  k=fmat_dim[1];

  double probs [k];
  double pr [k];
  int counts [k];

  for(j=0; j<k; j++){
    probs[j]=prior[j];
  }

  for(sim=0; sim<total_sims; sim++){
    if(((double) sim)/(total_sims) > bar/21){
      printf("#");
      bar+=1.0;
      fflush(stdout);
    }

    for(j=0; j<k; j++){
      counts[j]=0;
    };

    for(i=0; i<n; i++){
      for(j=0; j<k; j++){
	pr[j]=probs[j]*fmat[j*n+i];
      };
      ti=sample(k, pr);
      counts[ti]+=1;
    };
  
    for(j=0; j<k; j++){
      probs[j]=((double) counts[j])/(alpha+n)+alpha*prior[j]/(alpha+n);
    };

    if(sim%thin==0 && sim>burn_in){
      for(j=0; j<k; j++){
	all_counts[j]+=counts[j];
      } 
    }
  }
  printf("\n");
}
