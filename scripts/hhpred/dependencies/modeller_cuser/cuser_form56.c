
// compile: gcc -shared -Wall -fPIC `/cluster/user/armin/bin/modeller9.10/bin/mod9.10 --cflags --libs` `pkg-config --cflags glib-2.0` -I/usr/include/python2.4 cuser_form56.c cuser_form56_wrap.c -o _cuser_form56.so -lm

// this is a reimplemtation of MODELLER's multiple Gaussian distance restraints
// there is an additional parameter W which sets the weight for all restraints explicity (instead of the physical_schedule in python file)

#include <glib.h>
#include "modeller.h"
#include "cuser_form56.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static const float RT = 6.0;
static const float W = 1;


float get_p(float stdev, float x, float delta)
{  
  return exp(-0.5 * delta * delta / (stdev*stdev)) / x; 
} 

///////////////////////
/* Evaluate the form */
///////////////////////
static int myform_eval(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, gboolean deriv, float *fderv, float *val)
{
  const int numComp = *modal;
  double w[numComp];
  double mean[numComp];
  double sd[numComp];
  double p[numComp];
  double v[numComp];
  int i;

  // fill in values from rsr file
  for (i=0; i<numComp; i++) {
    w[i] = pcsr[i];
    mean[i] = pcsr[numComp+i];
    sd[i] = pcsr[2*numComp+i];
  }
  
  // trick to ensure robust calculation of derivative
  double rgauss1 = 10;
  double rgauss2 = 200;

  for (i=0; i<numComp; i++) {
    // calculate normalized deviation
    v[i] = (feat[0] - mean[i])/sd[i];

    // ignore violations that are too large
    if (v[i] > rgauss2) {
      v[i] = rgauss2  - 1E-10;
    }
    if (v[i] > rgauss1 && v[i] <= rgauss2) {
      double M = 37; // maximal exponent for exp
      double A = (rgauss2 - M) / (M * (rgauss2 - rgauss1));
      double B = (rgauss2 / M) * ((M - rgauss1) / (rgauss2 - rgauss1));
      v[i] = v[i] / (A*abs(v[i]) + B);
    }     
  }
 
  double sum = 0;

  for (i=0; i<numComp; i++) {
    p[i] = 1.0 / (sd[i]*sqrt(2*M_PI)) * exp( -0.5*v[i]*v[i] );
    sum += w[i] * p[i];
  }

  *val = W * RT * (-1) * log( sum );
  
  if (deriv) {
    double dfh = 0;
    for (i=0; i<numComp; i++) {
      dfh += w[i]*p[i] * (v[i] / sd[i]);
    }
    fderv[0] = W * RT * (1.0/sum) * dfh;
  }

  return 0; 
}


// 1 ///////////////////////////////////
     /* Evaluate the minimum violation*/
     ///////////////////////////////////
static int myform_vmin(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val)
{
  *val = 100;
  return 0;
}


// 2 /////////////////////////////////
     /* Evaluate the heavy violation*/
     /////////////////////////////////
static int myform_vheavy(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val)
{
  *val = 100;
  return 0;
}

// 3 ////////////////////////////////////////////
     /* Evaluate the relative minimum violation*/
     ////////////////////////////////////////////
static int myform_rvmin(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, float *val)
{

  *val = 100;
  return 0;
}

// 4 //////////////////////////////////////////
     /* Evaluate the relative heavy violation*/
     //////////////////////////////////////////
static int myform_rvheavy(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, float *val)
{

  *val = 0; 
  return 0;
}

 //////////////////////////////////////
 /* Evaluate the minimum feature mean*/
 //////////////////////////////////////
static int myform_minmean(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val)
{ 
  *val = 100; 
  return 0;
}

 ////////////////////////////////////
 /* Evaluate the heavy feature mean*/
 ////////////////////////////////////
static int myform_heavymean(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val)
{
  *val = 100;
  return 0;
}

////////////////////////////////////////////////////
/* Create the new form, and return its identifier */
////////////////////////////////////////////////////
int myform_create(void)
{  
  return mod_user_form_new(myform_eval, NULL, myform_vmin, NULL, myform_vheavy,
                           NULL, myform_rvmin, NULL, myform_rvheavy, NULL,
                           myform_minmean, NULL, myform_heavymean, NULL);
}

