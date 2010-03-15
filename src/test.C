#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream 
#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>
inline float log2(float x)  {return (x<=0? -100000.0:1.442695041*log(x));}


inline float fast_log2(float x)
{
  static float lg2[1025];         // lg2[i] = log2[1+x/1024]
  static float diff[1025];        // diff[i]= (lg2[i+1]-lg2[i])/8096 (for interpolation)
  static char initialized=0;
  if (x<=0) return -100000; 
  if (!initialized)   //First fill in the arrays lg2[i] and diff[i]
    {
      float prev = 0.0f;
      lg2[0] = 0.0f;
      for (int i=1; i<=1024; ++i)
        {
          lg2[i] = log(float(1024+i))*1.442695041-10.0f;
          diff[i-1] = (lg2[i]-prev)*1.2352E-4;
          prev = lg2[i];
        }
      initialized=1;
    }
  int a = (((*((int *)&x)) & 0x7F800000) >>23 )-0x7f; // exponent
  int b =  ((*((int *)&x)) & 0x007FE000) >>13; // first 10 bits of mantisse
  int c =  ((*((int *)&x)) & 0x00001FFF);      // further 13 bits of mantisse
  return a + lg2[b] + diff[b]*(float)(c);
}

inline float flog2(float x)
{
  if (x<=0) return -128;
  int *px = (int*)(&x);                 // store address of float as pointer to long int
  float e = (float) (((*px & 0x7F800000) >>23 )-0x7f); // shift right by 23 bits and subtract 127 = 0x7f => exponent
  *px =  ((*px & 0x007FFFFF) | 0x3f800000);  // set exponent to 127 (i.e., 0) 
  x -= 1.0;         // and calculate x-1.0 
  x *= (1.441740 + x*(-0.7077702 +x*(0.4123442 +x*(-0.1903190+x*0.0440047)))); // polynomial approximation of log(1+x)
  return x+e;
}


/////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  float sum=0.0;
  for (int i=1; i<=100000000; i++)
    {
      sum += fast_log2((float)i);
    } 
  printf("sum=%8.3f\n",sum);

  float x = 1.5;
  printf("flog2(%8.3e)=%8.4f log2(%8.3e)=%8.4f diff=%9.6f\n",x,flog2(x),x,log2(x),flog2(x)-log2(x));
  exit(0);
}
