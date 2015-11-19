//go to: cd /cluster/bioprogs/modeller9v5_RSR/modlib/modeller 
//gcc -shared -Wall -fPIC `/cluster/user/armin/bin/modeller9.10/bin/mod9.10 --cflags --libs` `pkg-config --cflags glib-2.0` -I/usr/include/python2.4 cuser_form54.c cuser_form54_wrap.c -o _cuser_form54.so -lm
// new restraints with phylogenetic weight - python version 2.4, see also cuser_form59.c for python version 2.6

#include <glib.h>
#include "modeller.h"
#include "cuser_form54.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static const float RT = 6.0;

/* Decode parameters from Modeller parameter array */
inline static void get_param(const float pcsr[7], float *w, float *mean, float *stdev, float *w1, float *mean1, float *stdev1, float *phyloW)
{
  *w = pcsr[0] + 1E-5;  
  *mean =  pcsr[1]; 	
  *stdev = pcsr[2];
  *w1 = pcsr[3]+1E-5;	
  *mean1 =  pcsr[4]; 	
  *stdev1 = pcsr[5];
  *phyloW = pcsr[6];
  //*phyloW = 1.0;
  //check if sum w+w1 equals 1 
  float sum = *w + *w1;
  *w /= sum;
  *w1 /= sum;
}

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
  float w, mean, stdev, w1, mean1, stdev1, phyloW;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1, &phyloW);
  if (feat[0]<=0) {fprintf (stderr, "Error in cuser_form.c, myform_eval: negative argument of log!\n");return 1;}
  float logdist = log(feat[0]);
  float delta = logdist - mean;
  float delta1 = logdist - mean1;

  float arg = delta1*delta1/(2.*stdev1*stdev1) - delta*delta/(2.*stdev*stdev);
  float argBar = delta1/(stdev1*stdev1*feat[0]) - delta/(stdev*stdev*feat[0]);


  float com = w*stdev1/stdev * exp(arg);
  //printf("com=%f, arg=%f, argBar=%f, w=%f, feat=%f, logdist=%f, mean=%f, mean1=%f, stdev=%f, stdev1=%f\n", com, arg, argBar, w, feat[0], logdist, mean, mean1, stdev, stdev1);


  *val = RT * (-log(com + 1 - w) ) * phyloW;		
  
  
  if (deriv) {
    fderv[0] = RT * ((-1) * (com * argBar) / (com + 1 - w)) * phyloW;
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
  float w, mean, stdev, w1, mean1, stdev1, phyloW;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1, &phyloW);

  if (feat[0]<=0) {
    fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");
    return 1;
  }

  float prec = 1.0/(stdev*stdev);
  float prec1 = 1.0/(stdev1*stdev1);

  float xmin = exp((mean1*prec1 - mean*prec)/(prec1 - prec));

  *val = feat[0] - xmin;

  /* if (abs(logdist-mean) < abs(logdist-mean1)) { */
/*     *val = logdist - mean; */
/*   } */
/*   else { */
/*     *val = logdist - mean1; */
/*   } */
 
  return 0;
}


// 2 /////////////////////////////////
     /* Evaluate the heavy violation*/
     /////////////////////////////////
static int myform_vheavy(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1, phyloW;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1, &phyloW);


  if (feat[0]<=0) {
    fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");
    return 1;
  }

  float prec = 1.0/(stdev*stdev);
  float prec1 = 1.0/(stdev1*stdev1);

  float xmin = exp((mean1*prec1 - mean*prec)/(prec1-prec));

  *val = feat[0] - xmin;

  //search global minimum 

  /* if (get_p(stdev, logdist, logdist-mean) < get_p(stdev1, logdist, logdist - mean1)) { */
/*     *val = logdist - mean1; */
/*   } */
/*   else { */
/*     *val = logdist - mean; */
/*   } */
  //float max = w*get_p(stdev,x,x-mean) + w1*get_p(stdev1,x,x-mean1);

  return 0;
}

// 3 ////////////////////////////////////////////
     /* Evaluate the relative minimum violation*/
     ////////////////////////////////////////////
static int myform_rvmin(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1, phyloW;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1, &phyloW);

  if (feat[0]<=0) {
    fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");
    return 1;
  }

  float logdist = log(feat[0]);

  
  if (abs(logdist-mean) < abs(logdist-mean1)) {
    *val = (logdist - mean)/stdev;
  }
  else {
    *val = (logdist - mean1)/stdev1;
  }
  
  return 0;
}

// 4 //////////////////////////////////////////
     /* Evaluate the relative heavy violation*/
     //////////////////////////////////////////
static int myform_rvheavy(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1, phyloW;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1, &phyloW);

  if (feat[0]<=0) {
    fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");
    return 1;
  }

  float logdist = log(feat[0]); 

  if (get_p(stdev, logdist, logdist-mean) < get_p(stdev1, logdist, logdist - mean1)) {
    *val = (logdist - mean1)/stdev1;
  }
  else {
    *val = (logdist - mean)/stdev;
  }
 
  return 0;
}

 //////////////////////////////////////
 /* Evaluate the minimum feature mean*/
 //////////////////////////////////////
static int myform_minmean(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1, phyloW;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1, &phyloW);

  if (feat[0]<=0) {
    fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");
    return 1;
  }

  //search nearest minimum
  float prec = 1/(stdev*stdev);
  float prec1 = 1/(stdev1*stdev1);

  float x = exp( (mean1*prec1 - mean*prec) / (prec1 - prec) );
 
  *val = x; 

  return 0;
}

 ////////////////////////////////////
 /* Evaluate the heavy feature mean*/
 ////////////////////////////////////
static int myform_heavymean(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1, phyloW;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1, &phyloW);

  if (feat[0]<=0) {
    fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");
    return 1;
  }

  //search global minimum 
  float prec = 1/(stdev*stdev);
  float prec1 = 1/(stdev1*stdev1);

  float x = exp( (mean1*prec1 - mean*prec) / (prec1 - prec) );

  *val = x;

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
