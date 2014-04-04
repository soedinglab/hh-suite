// Substitution matrices and their background frequencies

#include "hhmatrices.h"


void SetBlosumMatrix(const float BlosumXX[])
{
  int a,b,n=0;
  if (v>=3) printf("Using the BLOSUM%2i matrix",par.matrix);
  for (a=0; a<20; ++a)
    for (pb[a]=0.0f, b=0; b<=a; ++b,++n)
      P[a][b] = BlosumXX[n];
  for (a=0; a<19; a++)
    for (b=a+1; b<20; ++b)
      P[a][b] = P[b][a];
  return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Set (global variable) substitution matrix with derived matrices and background frequencies
/////////////////////////////////////////////////////////////////////////////////////
void SetSubstitutionMatrix()
{
  int a,b;
  switch (par.matrix)
    {
    default:
    case 0:  //Gonnet matrix
      if (v>=3) std::cout<<"Using the Gonnet matrix ";
      for (a=0; a<20; ++a)
	for (pb[a]=0.0f, b=0; b<20; ++b)
	  P[a][b] = 0.000001f*Gonnet[a*20+b];
      break;

    case 30:  //BLOSUM30
      SetBlosumMatrix(Blosum30);
      break;
    case 40:  //BLOSUM40
      SetBlosumMatrix(Blosum40);
      break;
    case 50:  //BLOSUM50
      SetBlosumMatrix(Blosum50);
      break;
    case 62:  //BLOSUM62
      SetBlosumMatrix(Blosum62);
      break;
    case 65:  //BLOSUM65
      SetBlosumMatrix(Blosum65);
      break;
    case 80:  //BLOSUM80
      SetBlosumMatrix(Blosum80);
      break;
   }
  
  // Check transition probability matrix, renormalize P and calculate pb[a]
  float sumab=0.0f;
  for (a=0; a<20; a++)
    for (b=0; b<20; ++b) sumab+=P[a][b];
  for (a=0; a<20; a++)
    for (b=0; b<20; ++b) P[a][b]/=sumab;
  for (a=0; a<20; a++)
    for (pb[a]=0.0f, b=0; b<20; ++b) pb[a]+=P[a][b];

  //Compute similarity matrix for amino acid pairs (for calculating consensus sequence)
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)
      Sim[a][b] = P[a][b]*P[a][b]/P[a][a]/P[b][b];

  //Precompute matrix R for amino acid pseudocounts:
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)   
      R[a][b] = P[a][b]/pb[b]; //R[a][b]=P(a|b)
  
  //Precompute matrix R for amino acid pseudocounts:
  for (a=0; a<20; ++a)
    for (b=0; b<20; ++b)   
      S[a][b] = log2(R[a][b]/pb[a]); // S[a][b] = log2(P(a,b)/P(a)/P(b))
  
  // Evaluate sequence identity underlying substitution matrix
  if (v>=3)
    {
      float id=0.0f;
      float entropy=0.0f; 
      float entropy_pb=0.0f;
      float mut_info=0.0f;
      for (a=0; a<20; ++a) id+=P[a][a];
      for (a=0; a<20; ++a) entropy_pb-=pb[a]*log2(pb[a]);
      for (a=0; a<20; ++a) 
	  for (b=0; b<20; ++b) 
	    {
	      entropy-=P[a][b]*log2(R[a][b]);
	      mut_info += P[a][b]*S[a][b];
	    }
      
      printf(": sequence identity = %2.0f%%; entropy per column = %4.2f bits (out of %4.2f); mutual information = %4.2f bits\n",100*id,entropy,entropy_pb,mut_info);
    }

  if (v>=4) //Debugging: probability matrix and dissimilarity matrix 
    {
      std::cout<<"Check matrix: before renormalization sum P(a,b)= "<<sumab<<"...\n";//PRINT
      std::cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      std::cout<<"p[] ";
      for (a=0; a<20; a++)  printf("%4.1f ",100*pb[a]);
      std::cout<<std::endl<<"\nSubstitution matrix log2( P(a,b)/p(a)/p(b) ) (in bits):\n";
      std::cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
	  std::cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",S[a][b]);
	  std::cout<<std::endl;
	}
      std::cout<<std::endl<<"\nOdds matrix P(a,b)/p(a)/p(b):\n";
      std::cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
	  std::cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",P[b][a]/pb[a]/pb[b]);
	  std::cout<<std::endl;
	}
      std::cout<<std::endl<<"\nMatrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
      std::cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
	  std::cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.1f ",100*R[b][a]);
	  std::cout<<std::endl;
	}
      std::cout<<std::endl<<"\nProbability matrix P(a,b) (in 10E-6):\n";
      std::cout<<"      A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V\n";
      for (b=0; b<20; b++)
	{
	  std::cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%5.0f ",1000000*P[b][a]);
	  std::cout<<std::endl;
	}
      std::cout<<std::endl<<"Similarity matrix P(a,b)^2/P(a,a)/P(b,b) (in %):\n";
      std::cout<<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
	  std::cout<<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  printf("%4.0f ",100*Sim[b][a]);
	  std::cout<<std::endl;
	}
      std::cout<<std::endl;


    }
}
 

/////////////////////////////////////////////////////////////////////////////////////
// Set secondary structure substitution matrix
/////////////////////////////////////////////////////////////////////////////////////
void SetSecStrucSubstitutionMatrix()
{
  int A;        //observed ss state (determined dssp)
  int B,BB;     //predicted ss states (by psipred)
  int cf,ccf;   //confidence value of prediction
  float P73[NDSSP][NSSPRED][MAXCF];  //P73[cf][B][A] = P(A,B,cf)/P(A)/P(B,cf) = P(A|B,cf)/P(A)
  float sum;

  // S73[A][B][cf][b] = score for matching observed ss state A in query with state B in template
  // predicted with confidence cf, when query and template columns are diverged by b units 
  for (cf=0; cf<MAXCF; cf++)
    for (A=0; A<NDSSP; A++)
      for (B=0; B<NSSPRED; B++)
	{

	  P73[A][B][cf] = 1.-par.ssa + par.ssa*Ppred[cf*NSSPRED*NDSSP + B*NDSSP + A];
	  S73[A][B][cf] = log2(P73[A][B][cf]);
	}

  for (B=0; B<NSSPRED; B++)
    for (cf=0; cf<MAXCF; cf++)
      for (BB=0; BB<NSSPRED; BB++)
	for (ccf=0; ccf<MAXCF; ccf++)
	  {
	    sum=0;
	    for (A=1; A<NDSSP; A++)
	      sum += P73[A][B][cf] * P73[A][BB][ccf] * Pobs[A];
	    S33[B][cf][BB][ccf] = log2(sum);
	  }  
}
