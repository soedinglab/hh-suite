// Substitution matrices and their background frequencies

#include "hhmatrices.h"


void SetBlosumMatrix(const char matrix, const float BlosumXX[], float* pb, float P[20][20])
{
  int a,b,n=0;
  HH_LOG(LogLevel::DEBUG) << "Using the BLOSUM " << matrix << " matrix" << std::endl;
  for (a=0; a<20; ++a)
    for (pb[a]=0.0f, b=0; b<=a; ++b,++n)
      P[a][b] = BlosumXX[n];
  for (a=0; a<19; a++)
    for (b=a+1; b<20; ++b)
      P[a][b] = P[b][a];
}

/////////////////////////////////////////////////////////////////////////////////////
// Set (global variable) substitution matrix with derived matrices and background frequencies
/////////////////////////////////////////////////////////////////////////////////////
void SetSubstitutionMatrix(const char matrix, float* pb, float P[20][20], float R[20][20], float S[20][20], float Sim[20][20])
{
  int a,b;
  switch (matrix)
    {
    default:
    case 0:  //Gonnet matrix
      HH_LOG(LogLevel::DEBUG) << "Using the Gonnet matrix" << std::endl;
      for (a=0; a<20; ++a)
	for (pb[a]=0.0f, b=0; b<20; ++b)
	  P[a][b] = 0.000001f*Gonnet[a*20+b];
      break;

    case 30:  //BLOSUM30
      SetBlosumMatrix(matrix, Blosum30, pb, P);
      break;
    case 40:  //BLOSUM40
      SetBlosumMatrix(matrix, Blosum40, pb, P);
      break;
    case 50:  //BLOSUM50
      SetBlosumMatrix(matrix, Blosum50, pb, P);
      break;
    case 62:  //BLOSUM62
      SetBlosumMatrix(matrix, Blosum62, pb, P);
      break;
    case 65:  //BLOSUM65
      SetBlosumMatrix(matrix, Blosum65, pb, P);
      break;
    case 80:  //BLOSUM80
      SetBlosumMatrix(matrix, Blosum80, pb, P);
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
  if (Log::reporting_level() >= LogLevel::DEBUG) {
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
      
      HH_LOG(LogLevel::DEBUG) << "sequence identity = " << 100*id << " %; entropy per column = " << entropy << " bits (out of " << entropy_pb << "); mutual information = " << mut_info << " bits" << std::endl;
  }

  //Debugging: probability matrix and dissimilarity matrix
  if (Log::reporting_level() >= LogLevel::DEBUG1) {
      HH_LOG(LogLevel::DEBUG) << "Check matrix: before renormalization sum P(a,b)= "<<sumab<<"...\n";//PRINT
      HH_LOG(LogLevel::DEBUG) <<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      HH_LOG(LogLevel::DEBUG) <<"p[] ";
      for (a=0; a<20; a++)  HH_LOG(LogLevel::DEBUG) << 100*pb[a] << " ";
      HH_LOG(LogLevel::DEBUG) <<std::endl<<"\nSubstitution matrix log2( P(a,b)/p(a)/p(b) ) (in bits):\n";
      HH_LOG(LogLevel::DEBUG) <<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
      HH_LOG(LogLevel::DEBUG) << i2aa(b) << "   ";
	  for (a=0; a<20; a++)  HH_LOG(LogLevel::DEBUG) << S[a][b] << " ";
	  HH_LOG(LogLevel::DEBUG) << std::endl;
	}
      HH_LOG(LogLevel::DEBUG) << std::endl << "\nOdds matrix P(a,b)/p(a)/p(b):\n";
      HH_LOG(LogLevel::DEBUG) <<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
    	  HH_LOG(LogLevel::DEBUG) <<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  HH_LOG(LogLevel::DEBUG) << P[b][a]/pb[a]/pb[b] << " ";
	  HH_LOG(LogLevel::DEBUG) <<std::endl;
	}
      HH_LOG(LogLevel::DEBUG) <<std::endl<<"\nMatrix of conditional probabilities P(a|b) = P(a,b)/p(b) (in %):\n";
      HH_LOG(LogLevel::DEBUG) <<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
      HH_LOG(LogLevel::DEBUG) <<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  HH_LOG(LogLevel::DEBUG) << 100*R[b][a] << " ";
	  HH_LOG(LogLevel::DEBUG) <<std::endl;
	}
      HH_LOG(LogLevel::DEBUG) <<std::endl<<"\nProbability matrix P(a,b) (in 10E-6):\n";
      HH_LOG(LogLevel::DEBUG) <<"      A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V\n";
      for (b=0; b<20; b++)
	{
      HH_LOG(LogLevel::DEBUG) <<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  HH_LOG(LogLevel::DEBUG) << 1000000*P[b][a] << " ";
	  HH_LOG(LogLevel::DEBUG) <<std::endl;
	}
      HH_LOG(LogLevel::DEBUG) <<std::endl<<"Similarity matrix P(a,b)^2/P(a,a)/P(b,b) (in %):\n";
      HH_LOG(LogLevel::DEBUG) <<"      A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V\n";
      for (b=0; b<20; b++)
	{
      HH_LOG(LogLevel::DEBUG) <<i2aa(b)<<"   ";
	  for (a=0; a<20; a++)  HH_LOG(LogLevel::DEBUG) << 100*Sim[b][a] << " ";
	  HH_LOG(LogLevel::DEBUG) <<std::endl;
	}
      HH_LOG(LogLevel::DEBUG) <<std::endl;
  }
}
 

/////////////////////////////////////////////////////////////////////////////////////
// Set secondary structure substitution matrix
/////////////////////////////////////////////////////////////////////////////////////
void SetSecStrucSubstitutionMatrix(const float ssa, float S73[NDSSP][NSSPRED][MAXCF], float S33[NSSPRED][MAXCF][NSSPRED][MAXCF])
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

	  P73[A][B][cf] = 1.-ssa + ssa*Ppred[cf*NSSPRED*NDSSP + B*NDSSP + A];
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
