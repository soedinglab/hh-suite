// comile
// gcc -shared -Wall -fPIC `/cluster/user/armin/bin/modeller9.10/bin/mod9.10 --cflags --libs` `pkg-config --cflags glib-2.0` -I/usr/include/python2.4 cuser_form52.c cuser_form52_wrap.c -o _cuser_form52.so -lm

#include <glib.h>
#include "modeller.h"
#include "cuser_form.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* RT at 297.15K, in kcal/mol */
//static const float RT = 0.5900991;
static const float RT = 6.0;

/* Decode parameters from Modeller parameter array */
inline static void get_param(const float pcsr[6], float *w, float *mean, float *stdev, float *w1, float *mean1, float *stdev1)
{
  *w = pcsr[0]+1E-5;  
  *mean =  pcsr[1]; 	
  *stdev = pcsr[2];
  *w1 = pcsr[3]+1E-5;	
  *mean1 =  pcsr[4]; 	
  *stdev1 = pcsr[5];
  //check if sum w+w1 equals 1 
  float sum = *w + *w1;
  *w /= sum;
  *w1 /= sum;
}

float get_p(float stdev, float x, float delta)
{  
  return exp(-0.5 * delta * delta / (stdev*stdev))/x; 
} 

///////////////////////
/* Evaluate the form */
///////////////////////
static int myform_eval(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, gboolean deriv, float *fderv, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1);   
  if (feat[0]<=0) {fprintf (stderr, "Error in cuser_form.c, myform_eval: negative argument of log!\n");return 1;}
  float logdist = log(feat[0]);
  float delta = logdist-mean;
  float delta1 = logdist-mean1;

  float arg = delta*delta/(2.*stdev*stdev);//printf("*x: %8.4f\n",arg); 
  float arg1= delta1*delta1/(2.*stdev1*stdev1);//printf("*y: %8.4f\n",arg1); 

  if (arg < arg1) {
    	*val = RT * (log(feat[0]) + arg - log(w/stdev + w1/stdev1 * exp(arg-arg1)));
  }
  else{    
    	*val = RT * (log(feat[0]) + arg1 - log(w1/stdev1 + w/stdev * exp(arg1-arg)));
  }
  
  if (deriv) {
	if (arg < arg1) {
      		fderv[0] = RT * (1./feat[0] * (1.+ ( w/(stdev*stdev*stdev)*delta + (w1)/(stdev1*stdev1*stdev1) * delta1 * exp(arg-arg1) ) / (w/stdev + (w1)/stdev1 * exp(arg-arg1)) ) );
	}
	else{
		fderv[0] = RT * (1./feat[0] * (1.+ ( w/(stdev*stdev*stdev)*delta*exp(arg1-arg) + (w1)/(stdev1*stdev1*stdev1) * delta1) / (w/stdev * exp(arg1-arg) + (w1)/stdev1) ) );
	}
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
  float w, mean, stdev, w1, mean1, stdev1;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1);
  if (feat[0]<=0) {fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");return 1;}
  //search nearest minimum
  float x = exp(mean-stdev*stdev);
  float x1 = exp(mean1-stdev1*stdev1);
  if (abs(x-feat[0]) <= abs(x1-feat[0])) {*val = feat[0]-x;} 
  else{*val=feat[0]-x1;} 
  return 0;
}


// 2 /////////////////////////////////
     /* Evaluate the heavy violation*/
     /////////////////////////////////
static int myform_vheavy(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1);
  if (feat[0]<=0) {fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");return 1;}
  //search global minimum 
  float x = exp(mean-stdev*stdev);
  float x1 = exp(mean1-stdev1*stdev1);
  float max = w*get_p(stdev,x,x-mean) + w1*get_p(stdev1,x,x-mean1);
  float max1 = w*get_p(stdev,x1,x1-mean) + w1*get_p(stdev1,x1,x1-mean1);
  if (max>=max1){*val = feat[0]-x;}
  else {*val = feat[0]-x1;}
  return 0;
}

// 3 ////////////////////////////////////////////
     /* Evaluate the relative minimum violation*/
     ////////////////////////////////////////////
static int myform_rvmin(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1);
  if (feat[0]<=0) {fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");return 1;}
  float logdist = log(feat[0]); 
  //search nearest minimum  & according stdev
  float x = exp(mean-stdev*stdev);
  float x1 = exp(mean1-stdev1*stdev1);
  if (abs(x-feat[0]) <= abs(x1-feat[0])) {*val = (logdist-mean)/stdev;} 
  else{*val = (logdist-mean1)/stdev1;} 
  return 0;
}

// 4 //////////////////////////////////////////
     /* Evaluate the relative heavy violation*/
     //////////////////////////////////////////
static int myform_rvheavy(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1);
  if (feat[0]<=0) {fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");return 1;}
  float logdist = log(feat[0]); 
  //search global minimum & according stdev
  float x = exp(mean-stdev*stdev);
  float x1 = exp(mean1-stdev1*stdev1);
  float max = w*get_p(stdev,x,x-mean) + w1*get_p(stdev1,x,x-mean1);
  float max1 = w*get_p(stdev,x1,x1-mean) + w1*get_p(stdev1,x1,x1-mean1);
  if (max>=max1){*val = (logdist-mean)/stdev;} 
  else{*val = (logdist-mean1)/stdev1;}
  return 0;
}

 //////////////////////////////////////
 /* Evaluate the minimum feature mean*/
 //////////////////////////////////////
static int myform_minmean(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1);
  if (feat[0]<=0) {fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");return 1;}
  //search nearest minimum
  float x = exp(mean-stdev*stdev);
  float x1 = exp(mean1-stdev1*stdev1);
  if (abs(x-feat[0]) <= abs(x1-feat[0])) {*val = x;} 
  else{*val=x1;} 
  return 0;
}

 ////////////////////////////////////
 /* Evaluate the heavy feature mean*/
 ////////////////////////////////////
static int myform_heavymean(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val)
{
  float w, mean, stdev, w1, mean1, stdev1;
  get_param(pcsr, &w, &mean, &stdev, &w1, &mean1, &stdev1);
  if (feat[0]<=0) {fprintf (stderr, "Error in cuser_form.c, negative argument of log!\n");return 1;}
  //search global minimum 
  float x = exp(mean-stdev*stdev);
  float x1 = exp(mean1-stdev1*stdev1);
  float max = w*get_p(stdev,x,x-mean) + w1*get_p(stdev1,x,x-mean1);
  float max1 = w*get_p(stdev,x1,x1-mean) + w1*get_p(stdev1,x1,x1-mean1);
  if (max>=max1){*val = x;}
  else {*val = x1;} 
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
/*
int myform_create(void)
{  
  return mod_user_form_new(myform_eval, NULL, myform_vheavy, NULL, myform_vheavy,
                           NULL, myform_rvheavy, NULL, myform_rvheavy, NULL,
                           myform_heavymean, NULL, myform_heavymean, NULL);
}
*/
