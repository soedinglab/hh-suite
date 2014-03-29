/*
 * hhhitlist-inl.h
 *
 *  Created on: Mar 28, 2014
 *      Author: meiermark
 */

#ifndef HHHITLIST_INL_H_
#define HHHITLIST_INL_H_

/////////////////////////////////////////////////////////////////////////////////////
// Calculate output of hidden neural network units
/////////////////////////////////////////////////////////////////////////////////////
inline float calc_hidden_output(const float* weights, const float* bias,
    float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm) {
  float res;
  // Calculate activation of hidden unit = sum of all inputs * weights + bias
  res = Lqnorm * weights[0] + Ltnorm * weights[1] + Nqnorm * weights[2]
      + Ntnorm * weights[3] + *bias;
  res = 1.0 / (1.0 + exp(-(res))); // logistic function
  return res;
}

////////////////////////////////////////////////////////////////////////////////////
//// Neural network regressions of lamda for EVD
/////////////////////////////////////////////////////////////////////////////////////
inline float lamda_NN(float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm) {
  const int inputs = 4;
  const int hidden = 4;
  const float biases[] = { -0.73195, -1.43792, -1.18839, -3.01141 }; // bias for all hidden units
  const float weights[] = { // Weights for the neural networks (column = start unit, row = end unit)
      -0.52356, -3.37650, 1.12984, -0.46796, -4.71361, 0.14166, 1.66807,
          0.16383, -0.94895, -1.24358, -1.20293, 0.95434, -0.00318, 0.53022,
          -0.04914, -0.77046, 2.45630, 3.02905, 2.53803, 2.64379 };
  float lamda = 0.0;
  for (int h = 0; h < hidden; h++) {
    lamda += calc_hidden_output(weights + inputs * h, biases + h, Lqnorm,
        Ltnorm, Nqnorm, Ntnorm) * weights[hidden * inputs + h];
  }
  return lamda;
}

////////////////////////////////////////////////////////////////////////////////////
//// Neural network regressions of mu for EVD
/////////////////////////////////////////////////////////////////////////////////////
inline float mu_NN(float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm) {
  const int inputs = 4;
  const int hidden = 6;
  const float biases[] = { -4.25264, -3.63484, -5.86653, -4.78472, -2.76356,
      -2.21580 };  // bias for all hidden units
  const float weights[] = { // Weights for the neural networks (column = start unit, row = end unit)
      1.96172, 1.07181, -7.41256, 0.26471, 0.84643, 1.46777, -1.04800, -0.51425,
          1.42697, 1.99927, 0.64647, 0.27834, 1.34216, 1.64064, 0.35538,
          -8.08311, 2.30046, 1.31700, -0.46435, -0.46803, 0.90090, -3.53067,
          0.59212, 1.47503, -1.26036, 1.52812, 1.58413, -1.90409, 0.92803,
          -0.66871 };
  float mu = 0.0;
  for (int h = 0; h < hidden; h++) {
    mu += calc_hidden_output(weights + inputs * h, biases + h, Lqnorm, Ltnorm,
        Nqnorm, Ntnorm) * weights[hidden * inputs + h];
  }
  return 20.0 * mu;
}

////////////////////////////////////////////////////////////////////////////////////
//// Neural network regressions of alpha for Evalue correction factor
/////////////////////////////////////////////////////////////////////////////////////
inline float alpha_NN(float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm) {
  const int inputs = 4;
  const int hidden = 4;
  const float biases[] = { 7.89636, 3.68944, 2.05448, 3.69149 }; // bias for all hidden units
  const float alpha_bias = 1.33439;
  const float weights[] = { // Weights for the neural networks (column = start unit, row = end unit)
      -6.72336, -4.73393, -2.15446, -4.75140, -14.54957, 4.05462, 0.57951,
          3.55780, 2.08289, -1.81976, -1.19936, -17.35097, 1.53268, -8.13514,
          -2.50677, 1.51106, 6.37397, -0.36254, 0.16279, -1.32174 };
  float alpha = 0.0;
  for (int h = 0; h < hidden; h++) {
    alpha += calc_hidden_output(weights + inputs * h, biases + h, Lqnorm,
        Ltnorm, Nqnorm, Ntnorm) * weights[hidden * inputs + h];
  }
  alpha = 1.0 / (1.0 + exp(-(alpha + alpha_bias))); // logistic function
  return alpha;
}

////////////////////////////////////////////////////////////////////////////////////
//// Neural network regressions of beta for Evalue correction factor
/////////////////////////////////////////////////////////////////////////////////////
inline float beta_NN(float Lqnorm, float Ltnorm, float Nqnorm, float Ntnorm) {
  const int inputs = 4;
  const int hidden = 4;
  const float biases[] = { 7.89636, 3.68944, 2.05448, 3.69149 }; // bias for all hidden units
  const float beta_bias = 5.43347;
  const float weights[] = { // Weights for the neural networks (column = start unit, row = end unit)
      -6.72336, -4.73393, -2.15446, -4.75140, -14.54957, 4.05462, 0.57951,
          3.55780, 2.08289, -1.81976, -1.19936, -17.35097, 1.53268, -8.13514,
          -2.50677, 1.51106, -2.27841, -7.79426, -9.53092, 3.65717 };
  float beta = 0.0;
  for (int h = 0; h < hidden; h++) {
    beta += calc_hidden_output(weights + inputs * h, biases + h, Lqnorm, Ltnorm,
        Nqnorm, Ntnorm) * weights[hidden * inputs + h];
  }
  beta = 1.0 / (1.0 + exp(-(beta + beta_bias))); // logistic function
  return beta;
}

#endif /* HHHITLIST_INL_H_ */
