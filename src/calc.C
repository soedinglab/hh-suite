#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <cstdio>     // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;

float query_len = 200.0;
float template_len = 200.0;
float query_neff = 2.0;
float template_neff = 2.0;

float q_len_norm;
float t_len_norm;
float q_neff_norm;
float t_neff_norm;


/////////////////////////////////////////////////////////////////////////////////////
// Neural network parameters
/////////////////////////////////////////////////////////////////////////////////////
const int mu_hidden = 6;
const int lamda_hidden = 4;

// Bias for all hidden units
const float mu_bias[] = {-4.25264, -3.63484, -5.86653, -4.78472, -2.76356, -2.21580};
const float lamda_bias[] = {-0.73195, -1.43792, -1.18839, -3.01141};

// Weights for the neural networks (column = start unit, row = end unit)
const float mu_weights[] = {
  1.96172, 1.07181, -7.41256, 0.26471, 
  0.84643, 1.46777, -1.04800, -0.51425,
  1.42697, 1.99927, 0.64647, 0.27834,
  1.34216, 1.64064, 0.35538, -8.08311,
  2.30046, 1.31700, -0.46435, -0.46803,
  0.90090, -3.53067, 0.59212, 1.47503,
  -1.26036, 1.52812, 1.58413, -1.90409, 0.92803, -0.66871
};
const float lamda_weights[] = {
  -0.52356, -3.37650, 1.12984, -0.46796,
  -4.71361, 0.14166, 1.66807, 0.16383,
  -0.94895, -1.24358, -1.20293, 0.95434,
  -0.00318, 0.53022, -0.04914, -0.77046,
  2.45630, 3.02905, 2.53803, 2.64379
};
/////////////////////////////////////////////////////////////////////////////////////


void help()
{
  printf("\n");
  printf("Calculate lamda and mu values for given length and diversities\n\n");
  printf("Usage: calc [options]\n");
  printf(" -q_l  <float>   Query length       (default: %7.2f)\n", query_len);
  printf(" -q_n  <float>   Query diversity    (default: %7.2f)\n", query_neff);
  printf(" -t_l  <float>   Template length    (default: %7.2f)\n", template_len);
  printf(" -t_n  <float>   Template diversity (default: %7.2f)\n", template_neff);
  printf("\n");
	 
}


/////////////////////////////////////////////////////////////////////////////////////
// calculate output of hidden units
/////////////////////////////////////////////////////////////////////////////////////
float calc_hidden_output(int unit, const float* weights, const float* bias)
{
  float res;

  // calculate input of hidden unit (sum of all inputs * weights)
  res = q_len_norm * weights[0 + unit * 4] + t_len_norm * weights[1 + unit * 4] + q_neff_norm * weights[2 + unit * 4] + t_neff_norm * weights[3 + unit * 4];
  
  // calculate output of hidden unit
  res = 1 / (1 + exp(-(res + bias[unit])));

  return res;
}
/////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
  float lamda = 0;
  float mu = 0;

  // Read input
  for (int i=1; i<argc; i++) {
    
    if (!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help")) {help(); exit(0);} 
    else if (!strcmp(argv[i],"-q_l") && (i<argc-1))    query_len=atof(argv[++i]);
    else if (!strcmp(argv[i],"-q_n") && (i<argc-1))    query_neff=atof(argv[++i]);
    else if (!strcmp(argv[i],"-t_l") && (i<argc-1))    template_len=atof(argv[++i]);
    else if (!strcmp(argv[i],"-t_n") && (i<argc-1))    template_neff=atof(argv[++i]);

  }

  // normalize input values
  q_len_norm = log(query_len) / log(1000);
  t_len_norm = log(template_len) / log(1000);
  q_neff_norm = query_neff / 10;
  t_neff_norm = template_neff / 10;
  

  // calculate lamda
  for (int i = 0; i < lamda_hidden; i++) {
    lamda += calc_hidden_output(i, lamda_weights, lamda_bias) * lamda_weights[lamda_hidden * 4 + i]; 
  }
  
  // calculate mu
  for (int i = 0; i < mu_hidden; i++) {
    mu += calc_hidden_output(i, mu_weights, mu_bias) * mu_weights[mu_hidden * 4 + i]; 
  }
  // correct normalization of mu
  mu *= 20;

  ////////////
  // OUTPUT
  ////////////

  printf("\nResults:\n");
  printf("Query-length    = %7.2f  Query-diversity    = %7.2f\n", query_len, query_neff);
  printf("Template-length = %7.2f  Template-diversity = %7.2f\n", template_len, template_neff);
  printf("\n");
  printf("Lamda = %12.8f\n", lamda);
  printf("Mu    = %12.8f\n", mu);
  printf("\n");

  return 0;
}
