#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct {
  double *s;
  double *y;
  double rho;
  double a;
} l_bfgs;

typedef struct {
  int clique;
  int rank;
  int numresets;
  int maxiter;
  int maxtime;
  int M;
  int savesoln;
  int solvemethod;
  int randstart;
} params;

int     normalize(double*);
double  function(double*, double*, double*, double*);
int     gradient(double*, double*, double*, double*);
int     projgradient(double*);
double  dot(double*, double*, int);
int     addvec(double*, double*, double, double*, int);
double  infnorm(double*, int);
int     direction(void);
int     updateLBFGS(double);
int     stableset(double*, double*);
double  f_inexact_lsearch(double*, double*, double*, double*, double,
         double*, double, double*, double*, double*, double*);
double  f_bracketing(double*, double*, double*, double*, double,
         double, double, double, double, double, double, double*,
         double, double*, double*, double*, double*);
double  f_sectioning(double*, double*, double*, double*, double*,
         double*, double*, double*, double, double, double, double, double,
         double, double, double, double, double, double, double, double*, int);
int     f_interpolation_quadratic(double, double, double, double, double,
         double, double, double*);
double  min(double, double);
double  max(double, double);
clock_t GetTime(void);
double  TimeInSec(clock_t, clock_t);
double  reset(double*, double*, double*, double*);
double  resetfull(double*, double*, double*, double*);
int     checkmaximal(void);

int    n, m, oldest, iter, outeriter, tempint, maxcliquesize, rank;
int    *edgesi, *edgesj, stablesetsize, *beststableset, beststablesize;
int    *set1, *set2, *set3, *set4;
int    size1, size2, size3, size4;
double tempval, dirgradprod;
double *lambda, sigma, *grad, *projgrad, *dir;
double  etabar, omegabar, tau, gammabar, alpha_omega;
double  beta_omega, alpha_eta, beta_eta, alpha_star, beta_star;
double  eta_star, omega_star, alpha_k, eta_k, omega_k;
FILE   *solnfile;
params *par;
l_bfgs *vecs;

int main(int argc, char *argv[])
{
  int i, j, h, lambdacount, sigmacount, update, currmethod;
  int iteronstage, resetint;
  char edge, nonedge, **graph, filesoln[256];
  double totaltime, val, aval, alpha, relnorm, prevval, check[5], check2[5];
  double *x, *ax, valu[1], avalu[1], valv[1], avalv[1], *valedges, *avaledges;
  FILE *datafile;
  clock_t begin, end;

  if(argc != 3) {
    printf("Usage: %s <graph> <parameter file>\n", argv[0]);
    return 0;
  }

  totaltime = 0.0;
  begin = GetTime();

  par = (params*)malloc(sizeof(params)); 
  datafile = fopen(argv[2], "r"); 
  fscanf(datafile, "%d\n", &par->clique);
  fscanf(datafile, "%d\n", &par->rank);
  fscanf(datafile, "%d\n", &par->numresets);
  fscanf(datafile, "%d\n", &par->maxiter);
  fscanf(datafile, "%d\n", &par->maxtime);
  fscanf(datafile, "%d\n", &par->M);
  fscanf(datafile, "%d\n", &par->savesoln);
  par->solvemethod = 3;
  par->randstart = 1;
  fclose(datafile);

  rank = par->rank;
  if(rank != 1 && rank != 2) {
    printf("Error: Specify either rank-1 or -2 in parameter file.\n");
    return 0;
  }

  strcpy(filesoln, argv[1]);
  if(par->clique == 1) strcat(filesoln, ".clq_soln");
  if(par->clique == 0) strcat(filesoln, ".ss_soln");

  printf("\n");
  printf(" ==============================================\n");
  printf("                 a heuristic for finding large\n");
  printf("                   cliques and stable sets\n");
  printf("  ** Max-AO **   \n");
  printf("                  Sam Burer, Renato Monteiro\n");
  printf("                   Georgia Tech, v0.120301\n");
  printf(" ==============================================\n\n");

  datafile = fopen(argv[1], "r");
  fscanf(datafile, "%d %d\n", &n, &m);

  printf("   graph = %s\n", argv[1]);
  printf(" V, E, %% = %d, %d, %.2e\n", n, m, (double)200*m/(n*(n-1)));
  printf("    mode = ");
  if(par->clique == 1) printf("clique (%d,%d)\n\n", rank, par->numresets);
  else printf("stable set (%d,%d)\n\n", rank, par->numresets);

  printf(" ------  --------\n");
  printf("  best     time  \n");
  printf(" ------  --------\n");

  graph = (char**)malloc((n+1)*sizeof(char*));
  for(i = 1; i <= n; i++) {
    graph[i] = (char*)malloc((i+1)*sizeof(char));
  }
  if(par->clique == 0) { edge = 'e'; nonedge = 'n'; }
  else { edge = 'n'; nonedge = 'e'; }
  for(i = 1; i <= n; i++) for(j = 1; j <= i-1; j++)
    graph[i][j] = nonedge;
  for(h = 1; h <= m; h++) {
    fscanf(datafile, "%d %d %lf\n", &i, &j, &tempval); 
/*    fscanf(datafile, "%d %d\n", &i, &j); */
      if(i > j) graph[i][j] = edge;
      else graph[j][i] = edge;
  }
  if(par->clique == 1) m = n*(n-1)/2 - m;
  edgesi = (int*)malloc((m+1)*sizeof(int));
  edgesj = (int*)malloc((m+1)*sizeof(int));
  h = 1;
  for(i = 1; i <= n; i++) for(j = 1; j <= i-1; j++)
    if(graph[i][j] == 'e') { edgesi[h] = i; edgesj[h++] = j; }
  for(i = 1; i <= n; i++) free(graph[i]);
  free(graph);

  x = (double*)malloc((rank*n + 1)*sizeof(double));
  ax = (double*)malloc((rank*n + 1)*sizeof(double));
  valedges = (double*)malloc((m + 1)*sizeof(double));
  avaledges = (double*)malloc((m + 1)*sizeof(double));
  lambda = (double*)malloc((m+1)*sizeof(double));
  for(h = 1; h <= m; h++) valedges[h] = lambda[h] = 0.0;
  sigma = 1.0e1;
  grad = (double*)malloc((rank*n + 1)*sizeof(double));
  projgrad = (double*)malloc((rank*n + 1)*sizeof(double));
  dir = (double*)malloc((rank*n + 1)*sizeof(double));
  for(i = 1; i <= rank*n; i++) x[i] = (double)i;
  srand(time(0)); 
  if(par->randstart == 1) for(i = 1; i <= rank*n; i++) x[i] = (double)rand() - 0.5*RAND_MAX;
  normalize(x);
  val = function(x, valu, valv, valedges);
  gradient(x, valu, valv, valedges);
  projgradient(x);
  relnorm = infnorm(projgrad, rank*n)/max(fabs(val), 1.0);

  if(par->solvemethod == 2 || par->solvemethod == 3) {
     vecs = (l_bfgs*)malloc((par->M + 1)*sizeof(l_bfgs));
    for(i = 1; i <= par->M; i++) {
      (vecs+i)->s = (double*)malloc((rank*n + 1)*sizeof(double));
      (vecs+i)->y = (double*)malloc((rank*n + 1)*sizeof(double));
      for(j = 1; j <= rank*n; j++) (vecs+i)->s[j] = (vecs+i)->y[j] = 0.0;
      (vecs+i)->rho = (vecs+i)->a = 0.0;
    }
    oldest = 1;
  }

  set1 = (int*)malloc((n+1)*sizeof(int));
  set2 = (int*)malloc((n+1)*sizeof(int));
  set3 = (int*)malloc((n+1)*sizeof(int));
  set4 = (int*)malloc((n+1)*sizeof(int));
  beststableset = (int*)malloc((n+1)*sizeof(int));

  iter = 1; outeriter = 1; prevval = 1.0e10; beststablesize = 1;
  lambdacount = sigmacount = 0;
  etabar = 1.0;    /* Affects choice between updating penalty or lambda. */
  omegabar = 1.0;  /* Affects subproblem stopping criteria. */
  gammabar = 0.99;  /* Used in calculation of alpha_k after penalty change. */
  alpha_omega = 1.0; /* Used only in the initial choice of omega_0. */
  alpha_eta = min(1.0, alpha_omega) - 0.01; /* Used only in the initial choice of eta_0. */
  beta_omega = 0.5; /* Used in the update of subproblem stopping criteria omega_k. */
  beta_eta = min(1.0, beta_omega) - 0.1; /* Used in the update of eta_k, which affects choice
                                            between updating penalty or lambda. */
  beta_eta = 0.1;                       /* The smaller, the more likely lambda is updated. */
  alpha_k = min(1.0, gammabar);
  omega_k = omegabar*pow(alpha_k, alpha_omega);
  eta_k = etabar*pow(alpha_k, alpha_eta);


  while(iter <= par->maxiter && outeriter <= par->numresets) {


    resetint = 0;
    while(resetint == 0) {


    iteronstage = 0;
    if(par->solvemethod == 3) currmethod = 1;
    else currmethod = par->solvemethod;

    while(relnorm >= omega_k && resetint == 0) {
      if(++iteronstage == 100 && par->solvemethod == 3) {
/*        printf("Switching methods...\n"); */
        currmethod = 2;
      }

      if(currmethod == 1) {


        for(i = 1; i <= rank*n; i++) dir[i] = -projgrad[i];
        aval = 1.0e10;
        if(iter%10 == 1) alpha = 2.0; else alpha *= 2.0;
        /* while(aval - val >= 0.005*alpha*dot(grad, projgrad, rank*n) && alpha > 1.0e-15) { */
	tempint = 1;
        while(aval >= val && tempint++ <= 30) {
           alpha *= 0.5;
          addvec(ax, x, alpha, dir, rank*n);
          normalize(ax);
          aval = function(ax, avalu, avalv, avaledges);
        }
        val = aval;
        addvec(x, ax, 0.0, ax, rank*n);
        *valu = *avalu; *valv = *avalv;
        addvec(valedges, avaledges, 0.0, avaledges, m);

        gradient(x, valu, valv, valedges);
        projgradient(x);

      }

      if(currmethod == 2) {
        direction();
        dirgradprod = dot(dir, projgrad, rank*n);
        for(i = 1; i <= rank*n; i++) (vecs+oldest)->y[i] = -projgrad[i];
        aval = f_inexact_lsearch(x, valu, valv, valedges, val, &alpha, 1.0e20, ax, avalu, avalv, avaledges);
        updateLBFGS(alpha);
        val = aval;
        addvec(x, ax, 0.0, ax, rank*n);
        *valu = *avalu; *valv = *avalv;
        addvec(valedges, avaledges, 0.0, avaledges, m);
      }

      relnorm = infnorm(projgrad, rank*n)/max(fabs(val), 1.0);
      if(iter%50 == 0) {
        if(iter > 100 && relnorm < 1.0e-3 && fabs(val - prevval)/max(fabs(val), 1.0) < 1.0e-7) {
/*          printf("Forced completion...\n"); */
          relnorm = 1.0e-99;
        }
        prevval = val;
      }

      end = GetTime();
      totaltime += (double)TimeInSec(begin, end);
      begin = end;

     tempval = (*valu)*(*valu) + (*valv)*(*valv); 
     if(iter%100 == 0) {

        /* If the objective value is close to an integer. */
        /* This is a used as a stopping criteria. */

        if( fabs(tempval - (double)floor(tempval)) < 1.0e-2 || fabs((double)ceil(tempval) - tempval) < 1.0e-2 ) {

          /* If the objective value is near or below the size of the best
             stable set found so far, then reset. */

          if(fabs(tempval - beststablesize) < 0.1 || tempval < beststablesize - 0.9) resetint = 1;
          tempint = stableset(x, valedges);
          if(tempint > beststablesize) {
            beststablesize = tempint;
            printf("  %4d    %6d\n", beststablesize, (int)totaltime);
            if(fabs(tempval - tempint) < 0.001) resetint = 1;
            if(par->savesoln == 1) {
              solnfile = fopen(filesoln, "w");
              for(i = 1; i <= n; i++) if(beststableset[i] == 1) fprintf(solnfile, "%d\n", i);
              fclose(solnfile);
            }
          }
        }
        else {
          tempint = stableset(x, valedges);
          if(tempint > beststablesize) {
            beststablesize = tempint;
            printf("  %4d    %6d\n", beststablesize, (int)totaltime);
            if(par->savesoln == 1) {
              solnfile = fopen(filesoln, "w");
              for(i = 1; i <= n; i++) if(beststableset[i] == 1) fprintf(solnfile, "%d\n", i);
              fclose(solnfile);
            }
          }
        }
      }

      tempval = (*valu)*(*valu) + (*valv)*(*valv); 
      check[iter%5] = tempval;
      if(iter > 5 && 0.2*(check[0] + check[1] + check[2] + check[3] + check[4]) < 1.0e-3) {
        resetint = 2;
      }

      check2[iter%5] = fabs((double)tempint - tempval);
      if((iter > 10 && 0.2*(check2[0] + check2[1] + check2[2] + check2[3] + check2[4]) < 1.0e-2) 
         || (iteronstage >= 500 && outeriter >= 25) ) { 
	resetint = 1;
      } 

      iter++;
      if(iter > par->maxiter || totaltime > par->maxtime || isnan(*valu) || isnan(*valv)) {
        /* printf("Got nan.\n"); */
        iter = par->maxiter + 1;
        resetint = 1;
      }

    }

    // Make sure rank=1 reset is full reset (without loss of outeriter)
    if(resetint == 1 && rank == 1) {
      resetint = 2;
      outeriter++;
    }

    if(resetint == 0) {

/*    printf("\nFinished subproblem ... "); */

    if(sqrt(dot(valedges, valedges, m)) <= eta_k) {
      if(lambdacount >= 3 || omega_k*pow(alpha_k, beta_omega) < 1.0e-6) {
/*        printf("Special situation (lambda)...\n"); */
        update = 2;
        lambdacount = 0;
        sigmacount = 1;
        if(beta_eta < 0.9) beta_eta += 0.1;
      }
      else {
        update = 1;
        sigmacount = 0;
        lambdacount++;
      }
    }
    else {
      if(sigmacount >= 2 || omegabar*pow(alpha_k, beta_omega) < 1.0e-6) {
/*        printf("Special situation (sigma)...\n"); */
        update = 1;
        sigmacount = 0;
        lambdacount = 1;
        if(beta_eta > 0.1) beta_eta -= 0.1;
        else beta_eta /= 10.0;
      }
      else {
        update = 2;
        lambdacount = 0;
        sigmacount++;
      }
    }

    if(update == 1) {
/*      printf("updating lambdas (%d) ...\n\n", lambdacount); */
      for(h = 1; h <= m; h++) lambda[h] += sigma*valedges[h];
      alpha_k = 1/sigma;
      eta_k *= pow(alpha_k, beta_eta);
      omega_k *= pow(alpha_k, beta_omega);
    }
    if(update == 2) {
/*      printf("increasing penalty (%d) ...\n\n", sigmacount); */
      sigma *= sqrt(10.0);
      alpha_k = gammabar/sigma;
      eta_k = etabar*pow(alpha_k, beta_eta);
      omega_k = omegabar*pow(alpha_k, beta_omega);
    }

    val = function(x, valu, valv, valedges);
    gradient(x, valu, valv, valedges);
    projgradient(x);
    relnorm = sqrt(dot(projgrad, projgrad, rank*n))/max(fabs(val), 1.0);

    if(par->solvemethod == 2 || par->solvemethod == 3) {
      for(i = 1; i <= par->M; i++) {
        for(j = 1; j <= rank*n; j++) (vecs+i)->s[j] = (vecs+i)->y[j] = 0.0;
        (vecs+i)->rho = (vecs+i)->a = 0.0;
      }
      oldest = 1;
    }

    }

    }

    outeriter++;
    if(resetint == 1) val = reset(x, valu, valv, valedges);
    if(resetint == 2) {
      val = resetfull(x, valu, valv, valedges);
      outeriter--;
    }
    relnorm = infnorm(projgrad, rank*n)/max(fabs(val), 1.0);
    resetint = 0;
    /* checkmaximal(); */
    prevval = 1.0e10;
    lambdacount = sigmacount = 0;
    alpha_k = min(1.0, gammabar);
    omega_k = omegabar*pow(alpha_k, alpha_omega);
    eta_k = etabar*pow(alpha_k, alpha_eta);

  }

  printf(" ------  --------\n");
  printf("  %4d    %6d\n\n", beststablesize, (int)totaltime);

  free(edgesi); free(edgesj);
  free(x); free(valedges);
  free(ax); free(avaledges);
  free(lambda);
  free(grad); free(projgrad);
  free(dir);
  if(par->solvemethod == 2 || par->solvemethod == 3) {
    for(i = 1; i <= par->M; i++) {
      free((vecs+i)->s);
      free((vecs+i)->y);
    }
    free(vecs);
  }
  free(set1); free(set2); free(set3); free(set4);
  free(beststableset);

  return 0;
 
}

double reset(double* x, double* valu, double* valv, double* valedges)
{
  int i, j;
  double val, w1, w2;

  /* printf("Resetting now...\n");  */

  w1 = (double)rand(); w2 = (double)rand();
  tempval = sqrt(w1*w1 + w2*w2); w1 /= tempval; w2 /= tempval;
/*  for(i = 1; i <= rank*n; i++) x[i] = ((double)rand() - 0.5*(double)RAND_MAX)/(double)RAND_MAX; */
  for(i = 1; i <= rank*n; i++) x[i] = (double)rand()/(double)RAND_MAX;
  tempval = 1/sqrt(beststablesize);
  for(i = 1; i <= n; i++) {
    if(beststableset[i] == 1) {
      x[0] = ((double)rand() - 0.5*(double)RAND_MAX)/((double)RAND_MAX*(double)RAND_MAX);
      x[i] = x[i+n]*tempval*w1 - (1 - x[i+n])*x[0]*w2;
      if(rank == 2) x[i+n] = x[i+n]*tempval*w2 + (1 - x[i+n])*x[0]*w1;
      /* Sam: Does this perturbation make sense for rank 1? */
    }
    else {
      x[0] = ((double)rand() - 0.5*(double)RAND_MAX)/((double)RAND_MAX*(double)RAND_MAX);
      x[i] = - (1 - x[i+n])*x[0]*w2;
      if(rank == 2) x[i+n] = (1 - x[i+n])*x[0]*w1;
      /* Sam: Does this perturbation make sense for rank 1? */
    }
  }
  normalize(x);
  val = function(x, valu, valv, valedges);
  gradient(x, valu, valv, valedges);
  projgradient(x);

  if(par->solvemethod == 2 || par->solvemethod == 3) {
    for(i = 1; i <= par->M; i++) {
      for(j = 1; j <= rank*n; j++) (vecs+i)->s[j] = (vecs+i)->y[j] = 0.0;
      (vecs+i)->rho = (vecs+i)->a = 0.0;
    }
    oldest = 1;
  }
   sigma = 1.0e2;
   /* for(i = 1; i <= m; i++) lambda[i] = 0.0; */

  return val;
}

double resetfull(double* x, double* valu, double* valv, double* valedges)
{
  int i, j;
  double val;

  for(i = 1; i <= rank*n; i++) x[i] = (double)rand() - 0.5*RAND_MAX;
  normalize(x);
  val = function(x, valu, valv, valedges);
  gradient(x, valu, valv, valedges);
  projgradient(x);

  if(par->solvemethod == 2 || par->solvemethod == 3) {
    for(i = 1; i <= par->M; i++) {
      for(j = 1; j <= rank*n; j++) (vecs+i)->s[j] = (vecs+i)->y[j] = 0.0;
      (vecs+i)->rho = (vecs+i)->a = 0.0;
    }
    oldest = 1;
  }

  sigma = 1.0e1;
  for(i = 1; i <= m; i++) lambda[i] = 0.0;

  return val;
}


int checkmaximal(void)
{
  int i, j, k, *verify;
 
  verify = (int*)malloc((n+1)*sizeof(int));
  for(i = 1; i <= n; i++) verify[i] = beststableset[i];
  for(k = 1; k <= m; k++) {
    i = edgesi[k]; j = edgesj[k];
    if(beststableset[i] == 1 && beststableset[j] == 0) verify[j] = 1;
    if(beststableset[i] == 0 && beststableset[j] == 1) verify[i] = 1;
  }
  for(i = 1; i <= n; i++) if(verify[i] == 0) {
    printf("Not maximal!\n");
    return 0;
  }
  printf("Yes, maximal.\n");
  free(verify);
  return 0;
}

int stableset(double* x, double* valedges)
{
  int i, k, h;
  double t1, t2, eps;
  int *tempervec;
  FILE *printfile;

  for(i = 1; i <= n; i++) set1[i] = set2[i] = set3[i] = set4[i] = 0;
  size1 = size2 = size3 = size4 = 0;
  for(i = 1; i <= n; i++) {
    t1 = x[i];
    if(rank == 2) t2 = x[i+n];
    else t2 = 0.0;
    /* Sam: is this correct for rank 1? */
    eps = infnorm(valedges, m);
/*     if(eps >= 0.125) return 0; */
/*     if(t1*t2 < 0 && fabs(t1) > sqrt(0.5*eps) && fabs(t2) > sqrt(0.5*eps)) { size1++; set1[i] = 1; } */
/*     if(t1*t2 > 0 && fabs(t1) > sqrt(0.5*eps) && fabs(t2) > sqrt(0.5*eps)) { size2++; set2[i] = 1; } */
/*     if(fabs(t1) > 2*sqrt(2*eps) && fabs(t2) < eps) { size3++; set3[i] = 1; } */
/*     if(fabs(t1) < eps && fabs(t2) > 2*sqrt(2*eps)) { size4++; set4[i] = 1; } */
    if(t1*t2 < 0 && fabs(t1) > sqrt(0.5*eps) && fabs(t2) > sqrt(0.5*eps)) { size1++; set1[i] = 1; }
    if(t1*t2 > 0 && fabs(t1) > sqrt(0.5*eps) && fabs(t2) > sqrt(0.5*eps)) { size2++; set2[i] = 1; }
    if(fabs(t2) < 0.5*eps && fabs(t1) > sqrt(eps*(1.0 + 0.25*eps)) && fabs(t1) + fabs(t2) > sqrt(2*eps)) { size3++; set3[i] = 1; }
    if(fabs(t1) < 0.5*eps && fabs(t2) > sqrt(eps*(1.0 + 0.25*eps)) && fabs(t1) + fabs(t2) > sqrt(2*eps)) { size4++; set4[i] = 1; }
  }

/*  printf("%d %d %d %d\n", size1, size2, size3, size4);  */

  if(iter%1000000000 == 0) {
    for(i = 1; i <= n; i++)
      printf("%d (%e, %e)\n", i, x[i], x[n+i]);
    exit(0);
  }

  if(size1*size2 < 0) {
    printfile = fopen("subgraph", "a");
    fprintf(printfile, "%d %d %d %d\n", size1, size2, size3, size4);
    tempervec = (int*)malloc((n+1)*sizeof(int));
    for(i = 1; i <= n; i++)
      tempervec[i] = set1[i] + set2[i] + set3[i] + set4[i];
    for(k = 1; k <= m; k++)
      if(tempervec[edgesi[k]]*tempervec[edgesj[k]] == 1)
        fprintf(printfile, "%d %d\n", edgesi[k], edgesj[k]);
    fprintf(printfile, "====\n");
    free(tempervec);
    fclose(printfile);
  }

  if(size1 > size2) {
    if(size3 > size4) {
      stablesetsize = size1 + size3;
      for(i = 1; i <= n; i++) {
        set1[i] += set3[i];
        /* if(set1[i] == 1) printf("%d\n", i); */
      }
      h = 1;
      for(k = 1; k <= m; k++)
        if(set1[edgesi[k]]*set1[edgesj[k]] == 1) h = 0;
      if(h == 1 && stablesetsize > beststablesize)
        for(i = 1; i <= n; i++) beststableset[i] = set1[i];
    }
    else {
      stablesetsize = size1 + size4;
      for(i = 1; i <= n; i++) {
        set1[i] += set4[i];
        /* if(set1[i] == 1) printf("%d\n", i); */
      }
      h = 1;
      for(k = 1; k <= m; k++)
        if(set1[edgesi[k]]*set1[edgesj[k]] == 1) h = 0;
      if(h == 1 && stablesetsize > beststablesize)
        for(i = 1; i <= n; i++) beststableset[i] = set1[i];
    }
  }
  else {
    if(size3 > size4) {
      stablesetsize = size2 + size3;
      for(i = 1; i <= n; i++) {
        set2[i] += set3[i];
        /* if(set2[i] == 1) printf("%d\n", i); */
      }
      h = 1;
      for(k = 1; k <= m; k++)
        if(set2[edgesi[k]]*set2[edgesj[k]] == 1) h = 0;
      if(h == 1 && stablesetsize > beststablesize)
        for(i = 1; i <= n; i++) beststableset[i] = set2[i];
    }
    else {
      stablesetsize = size2 + size4;
      for(i = 1; i <= n; i++) {
        set2[i] += set4[i];
        /* if(set2[i] == 1) printf("%d\n", i); */
      }
      h = 1;
      for(k = 1; k <= m; k++)
        if(set2[edgesi[k]]*set2[edgesj[k]] == 1) h = 0;
      if(h == 1 && stablesetsize > beststablesize)
        for(i = 1; i <= n; i++) beststableset[i] = set2[i];
    }
  }
/*  if(h == 1) {
    if(par->clique == 0) {
      printf("\nStable set of size %d found.\n\n", stablesetsize);
      fprintf(outfile, "\nStable set of size %d found.\n\n", stablesetsize);
    }
    if(par->clique == 1) {
      if(maxcliquesize > 1) {
        printf("\nClique of size %d found. (Max clique size is %d.)\n\n", stablesetsize, maxcliquesize);
        fprintf(outfile, "\nClique of size %d found. (Max clique size is %d.)\n\n", stablesetsize, maxcliquesize);
      }
      else {
        printf("\nClique of size %d found. (Max clique size is unknown.)\n\n", stablesetsize);
        fprintf(outfile, "\nClique of size %d found. (Max clique size is unknown.)\n\n", stablesetsize);
      }
    }
  } */

  if(h == 1) return stablesetsize;
  else return 0;
}


int normalize(double* x)
{
  int i;
  double norm;

  norm = sqrt(dot(x, x, rank*n));
  for(i = 1; i <= rank*n; i++) x[i] /= norm;

  return 0;
}

double function(double* x, double* valu, double* valv, double* valedges)
{
  int i, j, h;
  double *u, *v, sum;

  if(rank == 2) {
    u = x; v = x+n;
    *valu = *valv = 0.0;
    for(i = 1; i <= n; i++) { *valu += u[i]; *valv += v[i]; }
    sum = -(*valu)*(*valu) - (*valv)*(*valv);
    for(h = 1; h <= m; h++) {
      i = edgesi[h]; j = edgesj[h];
      tempval = valedges[h] = u[i]*u[j] + v[i]*v[j];
      sum += tempval*(lambda[h] + 0.5*sigma*tempval);
    }
  }
  else {
    u = x;
    *valu = *valv = 0.0;
    for(i = 1; i <= n; i++) *valu += u[i];
    sum = -(*valu)*(*valu);
    for(h = 1; h <= m; h++) {
      i = edgesi[h]; j = edgesj[h];
      tempval = valedges[h] = u[i]*u[j];
      sum += tempval*(lambda[h] + 0.5*sigma*tempval);
    }
  }

  return sum;
}

int gradient(double* x, double* valu, double* valv, double* valedges)
{
  int i, j, h;
  double *u, *v, *gradu, *gradv;
  double tempval1, tempval2;

  if(rank == 2) {
    u = x; v = x+n;
    gradu = grad; gradv = grad+n;
    tempval1 = -2*(*valu); tempval2 = -2*(*valv);
    for(i = 1; i <= n; i++) { gradu[i] = tempval1; gradv[i] = tempval2; }
    for(h = 1; h <= m; h++) {
      i = edgesi[h]; j = edgesj[h];
      tempval = lambda[h] + sigma*valedges[h];
      gradu[i] += tempval*u[j];
      gradu[j] += tempval*u[i];
      gradv[i] += tempval*v[j];
      gradv[j] += tempval*v[i];
    }
  }
  else {
    u = x;
    gradu = grad;
    tempval1 = -2*(*valu);
    for(i = 1; i <= n; i++) gradu[i] = tempval1;
    for(h = 1; h <= m; h++) {
      i = edgesi[h]; j = edgesj[h];
      tempval = lambda[h] + sigma*valedges[h];
      gradu[i] += tempval*u[j];
      gradu[j] += tempval*u[i];
    }
  }

  return 0;
}

int projgradient(double* x)
{
  int i;

  tempval = dot(grad, x, rank*n);
  for(i = 1; i <= rank*n; i++) projgrad[i] = grad[i] - tempval*x[i];

  return 0;
}

int direction(void)
{
  int i, ind;
  double beta;

  for(i = 1; i <= rank*n; i++) dir[i] = projgrad[i];

  for(i = 1; i <= par->M; i++) {
    if(oldest - i > 0) ind = oldest - i;
    else ind = oldest - i + par->M;
    (vecs+ind)->a = (vecs+ind)->rho*dot((vecs+ind)->s, dir, rank*n);
    addvec(dir, dir, -((vecs+ind)->a), (vecs+ind)->y, rank*n);
  }
  for(i = par->M; i >= 1; i--) {
    if(oldest - i > 0) ind = oldest - i;
    else ind = oldest - i + par->M;
    beta = (vecs+ind)->rho*dot((vecs+ind)->y, dir, rank*n);
    addvec(dir, dir, (vecs+ind)->a - beta, (vecs+ind)->s, rank*n);
  }

  for(i = 1; i <= rank*n; i++) dir[i] *= -1.0;

  return 0;
}

int updateLBFGS(double a)
{
  int i;

  addvec((vecs+oldest)->y, (vecs+oldest)->y, 1.0, projgrad, rank*n);
  for(i = 1; i <= rank*n; i++) (vecs+oldest)->s[i] = a*dir[i];
  (vecs+oldest)->rho = 1/dot((vecs+oldest)->y, (vecs+oldest)->s, rank*n);
  if((vecs+oldest)->rho > 1.0e10) (vecs+oldest)->rho = 1.0e10;
  /* printf("rho = %f\n", (vecs+oldest)->rho); */
  oldest = oldest%par->M + 1;

  return 0;
}


double dot(double* vec1, double* vec2, int length)
{
  int i;
  double sum;

  sum = 0.0;
  for(i = 1; i <= length; i++) sum += vec1[i]*vec2[i];

  return sum;
}

int addvec(double* vecnew, double* vec1, double scalar, double* vec2, int length)
{
  int i;

  for(i = 1; i <= length; i++)
    vecnew[i] = vec1[i] + scalar*vec2[i];

  return 0;
}

double infnorm(double* vec, int length)
{
  int i;
  double max;

  max = -1.0e15;
  for(i = 1; i <= length; i++)
    if(fabs(vec[i]) > max) max = fabs(vec[i]);

  return max;
}

double min(double a, double b)
{
  if(a < b) return a;
  else return b;
}

double max(double a, double b)
{
  if(a > b) return a;
  else return b;
}


/* 
  function to compute alpha that min. the function f(x+alpha*d)  
  
  The algorithm is as given in Practical Optimization by R. Fletcher,
  page 33 through 39.
  
  */

double f_inexact_lsearch(double* x, double* valu, double* valv, double* valedges, double f_val, double* alpha, double alpha_upbd, double* ax, double* avalu, double* avalv, double* avaledges) 
{
  /* define parameters */
  
  double sig = 0.1;                 /* Wolfe Powell parameter */
  double rho   = 0.01;                /* Goldstein parameter */
  /* the smaller, the tighter -- rho */
  double tau1  = 9.0;                 /* preset factor for jump of alpha */
  double tau2  = 0.1;                 /* preset factor for alpha in sect. */
  double tau3  = 0.5;                 /* preset factor for alpha in sect. */ 
  double f_lower_bound = -1.0e20;      /* lower bound of objective function */
  double value;

  value = f_bracketing(x,valu,valv,valedges,f_val,sig,rho,tau1,tau2,tau3,f_lower_bound,alpha,alpha_upbd,ax,avalu,avalv,avaledges);
  return value;
} 


/* function bracketing to get a bracket trial for alpha */

double f_bracketing(double* x, double* valu, double* valv, double* valedges, double f_val, double sig, double rho, double tau1, double tau2, double tau3, double f_lower_bound, double* alpha, double alpha_upbd, double* ax, double* avalu, double* avalv, double* avaledges)
{
  double a;                     /* one value of alpha bracket */
  double b;                     /* another value of alpha bracket */
  double a1;                    /* lower limit of updated alpha */
  double b1;                    /* upper limit of updated alpha */ 
  double localmu;               /* upper bound of initial guess for alpha */
  double f_zero;                /* alpha function at alpha = 0 */
  double df_zero;               /* derivative of alpha function at alpha=0 */
  double f_a;                   /* alpha function value at alpha=a */
  double df_a;                  /* derivative of alpha function at alpha=a */
  double f_b;                   /* alpha function at b */
  double df_b;                  /* derivative of alpha function at alpha=b */
  double alpha_prev = 0.0;      /* previous value of alpha, initially 0 */
  double f_alpha_prev;
  double df_alpha_prev;
  double f_alpha;
  double df_alpha;
  int k = 1;
  double tempval, value;

  f_alpha_prev = f_zero = f_val;
  df_alpha_prev = df_zero = dirgradprod;

  localmu = (f_lower_bound - f_zero)/(rho*df_zero);
  *alpha = 1.0;

  k = 1;
  while(k < 20) {
    if(k == 99) printf("Reached upper threshhold. (1)\n");

    addvec(ax, x, *alpha, dir, rank*n);
    normalize(ax);
    f_alpha = function(ax, avalu, avalv, avaledges);

    if(isnan(f_alpha)) f_alpha = 1.0e20;

    if (f_alpha <= f_lower_bound) return f_alpha; /* This probably never occurs. */

    if (f_alpha > f_zero + (*alpha)*rho*df_zero || f_alpha >= f_alpha_prev) {
      a = alpha_prev; f_a = f_alpha_prev; df_a = df_alpha_prev;
      b = *alpha;     f_b = f_alpha;
      value = f_sectioning(x, valu, valv, valedges, ax, avalu, avalv, avaledges, rho, tau2, tau3, sig, a, b, f_zero, df_zero, f_a, f_b, df_a, 0, alpha, 2);
      return value;
    }
    gradient(ax, avalu, avalv, avaledges);
    projgradient(ax);
    df_alpha = dot(dir, projgrad, rank*n);

    if (fabs(df_alpha) <= ((-sig)*df_zero)) return f_alpha;
    if (df_alpha >= 0.0) {
      a = *alpha;     f_a = f_alpha;      df_a = df_alpha;
      b = alpha_prev; f_b = f_alpha_prev; df_b = df_alpha_prev;
      value = f_sectioning(x, valu, valv, valedges, ax, avalu, avalv, avaledges, rho, tau2, tau3, sig, a, b, f_zero, df_zero, f_a, f_b, df_a, df_b, alpha, 3); 
      return value;
    }
    if(localmu <= (2*(*alpha) - alpha_prev)) {
      alpha_prev = *alpha;
      f_alpha_prev = f_alpha;
      df_alpha_prev = df_alpha;
      *alpha = localmu;
    }
    else {
      /* Store value of alpha before interpolation. */
      tempval = *alpha;
      /* Prepare for interpolation. */
      a = alpha_prev; f_a = f_alpha_prev; df_a = df_alpha_prev;
      b = *alpha;     f_b = f_alpha;      df_b = df_alpha;
      a1 = 2*(*alpha) - alpha_prev; 
      b1 = min(localmu, *alpha + tau1*(*alpha - alpha_prev));
      /* Do interpolation. */
      f_interpolation_quadratic(f_a, df_a, f_b, a, b, a1, b1, alpha);
      /* Prepare for next round. */
      alpha_prev = tempval;
      f_alpha_prev = f_alpha;
      df_alpha_prev = df_alpha;
    }
    k++;
  }

  /* printf("Line search failed (1).\n"); */
  return 0;

}


/* function to choose the section that gives the appropriate properties
   for alpha values */

double f_sectioning(double* x, double* valu, double* valv, double* valedges, double* ax, double* avalu, double* avalv, double* avaledges, double rho, double tau2, double tau3, double sig, double a, double b, double f_zero, double df_zero, double f_a, double f_b, double df_a, double df_b, double *alpha, int polytype)
{
  double a1;                 /* lower limit of alpha bracket */
  double b1;                 /* upper limit of alpha bracket */
  double f_alpha;
  double df_alpha;

  int k = 1;

  while (k < 20) {
    if(k == 99) printf("Reached upper threshhold. (2)\n");

    a1 = a + tau2*(b - a);
    b1 = b - tau3*(b - a); 

    f_interpolation_quadratic(f_a, df_a, f_b, a, b, a1, b1, alpha);

    addvec(ax, x, *alpha, dir, rank*n);
    normalize(ax);
    f_alpha = function(ax, avalu, avalv, avaledges);

    if(isnan(f_alpha)) f_alpha = 2*f_a;

    if(f_alpha > f_zero + rho*(*alpha)*df_zero || f_alpha >= f_a) {
      /* Prepare for next round. */
      a = a;      f_a = f_a;     df_a = df_a;
      b = *alpha; f_b = f_alpha;
    }
    else {
      gradient(ax, avalu, avalv, avaledges);
      projgradient(ax);
      df_alpha = dot(dir, projgrad, rank*n);

      if (fabs(df_alpha) <= ((- sig) * df_zero)) {
/*       if (df_alpha <= ((- sig) * df_zero)) { */
	return f_alpha;
      }
      if((b - a) * df_alpha >= 0.0) {
	b = a; f_b = f_a; df_b = df_a;
      }
      else {
	b = b; f_b = f_b; df_b = df_b;
      }
      a = *alpha; f_a = f_alpha; df_a = df_alpha;
    }
    k++;
  }

  /* printf("Line search failed (2).\n"); */
  return f_alpha;
}


/* function to generate next choose of alpha by quadratic interpolation */

int f_interpolation_quadratic (double f_a, double df_a, double f_b, double a, double b, double a1, double b1, double *alpha) 
{
  double za1, zb1, root;
  double endptmin;

  za1 = f_a + df_a*(a1 - a) + (f_b - f_a - (b-a)*df_a)*(a1 - a)*(a1 - a)/((b - a)*(b - a));
  zb1 = f_a + df_a*(b1 - a) + (f_b - f_a - (b-a)*df_a)*(b1 - a)*(b1 - a)/((b - a)*(b - a));
  if(za1 < zb1) endptmin = a1;
  else endptmin = b1;
  root = a - (b-a)*(b-a)*df_a/(2*(f_b - f_a - (b-a)*df_a));

  if(f_b - f_a - (b-a)*df_a < 0) {       /* Opens downward. */
    if(a1 < b1) {
      if(a1 <= root && root <= b1) *alpha = endptmin;
      if(root < a1) *alpha = b1;
      if(root > b1) *alpha = a1;
    }
    else {
      if(b1 <= root && root <= a1) *alpha = endptmin;
      if(root < b1) *alpha = a1;
      if(root > a1) *alpha = b1;
    }
  }
  else {                                 /* Opens upward. */
    if(a1 < b1) {
      if(a1 <= root && root <= b1) *alpha = root;
      if(root < a1) *alpha = a1;
      if(root > b1) *alpha = b1;
    }
    else {
      if(b1 <= root && root <= a1) *alpha = root;
      if(root < b1) *alpha = b1;
      if(root > a1) *alpha = a1;
    }
  }

  return 0;

}

clock_t GetTime(void)
{
  clock_t t;

  struct tms {
    clock_t    tms_utime;
    clock_t    tms_stime;
    clock_t    tms_cutime;
    clock_t    tms_cstime;
  } tmrec;

  t = times(&tmrec);
  t = tmrec.tms_utime+tmrec.tms_stime;

  return t;
}

double TimeInSec(clock_t head,
                 clock_t rear)
{
  double tscal;

  tscal=0.01;
  return ((double)(rear-head)*tscal);
}
