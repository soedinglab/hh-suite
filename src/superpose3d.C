// superpose3d.c
// Superpose two pdb chains by RMS minimization over a set of aligned residues
// Usage: superpose3d ali.txt in1.pdb in2.pdb out.pdb
//
// The program uses singular value decomposition (SVD) of the correlation matrix 
//        C = (1/n) Sum_k (s_k-sbar)(r_k-rbar)^t 
// Here, r_k, s_k are the coordinate vectors of the aligned pairs and rbar, sbar are their center-of-mass vectors.
// The SVD of C is 
//        C = U W V^t, 
// where U and V are orthogonal and W is a diagonal matrix. Then 
//        R = U V^t  
// is the rotation matrix that superposes the r_k with y_k with minimum RMS deviation:  
//        argmin_R' sum_k(s_k - R' r_k)^2 = R
// Reference: JH Challis (1995) A procedure for determining rigid body transformation parameters, J. Biomechanics 28, 733.

#include <stdio.h>    // printf
#include <stdlib.h>   // exit
#include <string>     // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
// #include <new>
// #include "efencepp.h"


const char VERSION[]="version 1.0  (August 2005)";
const int MINRES=-100;     // min residue index  
const int par.maxres=65536;    // max number of superposed atoms in pdb files and max residue index 
const int MAXLIN=64768;    // max number of lines when reading from standard input 
const int LINELEN=128;     // max length of line read in from input files; must be >= 82
const float blue=5.0;      // temperature factor corresponding to blue
const float cyan=30.0;     // temperature factor corresponding to cyan
const float orange=90.0;   // temperature factor corresponding to orange
const float red=99.99;     // temperature factor corresponding to red
char* program;             // set to argv[0]
int v=3;             // verbose mode

void Help()
{
  printf("\n");
  printf("%s - %s\n",program,VERSION);
  printf("Superpose two structures by RMS minimization over a set of aligned residues\n");
  printf("Usage: superpose3d [options] <alignment file> <in1.pdb> <in2.pdb> <out.pdb>\n");
  printf("\n");
  printf("The program reads in a file with pairs of aligned indices and two pdb files.\n");
  printf("It writes a pdb file with the coordinates of the first pdb file superposed\n");
  printf("onto the coordinates of the second file.\n");
  printf("The format for the alignment file (written as regex) is '(\\d+\\s+\\d+\\n)+'.\n");
  printf("The indices in this file refer to the residue numbers in the pdb files.\n");
  printf("A '-' in place of the alignment file means the program reads from stdin and a '-'\n");
  printf("in place of the output file will make the program write the pdb file to stdout.\n");
  printf("\n");
  printf("Options:\n");
  printf("-col        The superposed structures are color-coded through their temperature factor\n");
  printf("            cyan/blue and orange/red, for residues within / outside of the aligned set.\n");
  printf("-rms <rms>  use the maximum alignable subset of residue pairs for the superposition:\n");
  printf("            The program iterates and throws out the pair with the largest RMS until no\n");
  printf("            pairs are left with an RMS above <rms>.\n");
  printf("-v <int>    verbose mode\n");
  printf("\n");
  printf("(C) Johannes Soeding. This software is freely available under the GPL.\n");
}
 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg) 
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2)) 
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2)) 
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a)) 

/* Given a matrix a[0..m-1][0..n-1], this routine computes its singular value decomposition, A = U · W · V^t . 
   The matrix U replaces a on output. The diagonal matrix of singular values W is out-put as a vector w[0..n-1]. 
   The matrix V (not the transpose V^t) is output as v[0..n-1][0..n-1] . */
void SVD(float **a, int m, int n, float w[], float **v)
{
  float pythag(float a, float b);
  int flag,i,its,j,jj,k,l=1,nm=1;
  float anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  rv1=new(float[n]);
  g=scale=anorm=0.0;    
  
  // Householder reduction to bidiagonal form.
  for (i=0;i<n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i < m && i != n-1) {
      for (k=l;k<n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<m;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<n;k++) a[i][k] *= scale;
      }
    }
    anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  // Accumulation of right-hand transformations.
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g) {
	// Double division to avoid possible underflow.
	for (j=l;j<n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  // Accumulation of left-hand transformations.
  for (i=IMIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<m;j++) a[j][i] *= g;
    } else for (j=i;j<m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  // Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations.
  for (k=n-1;k>=0;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      // Test for splitting. Note that rv1[1] is always zero.
      for (l=k;l>=0;l--) {
	nm=l-1;
	if ((float)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((float)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	// Cancellation of rv1[l], if l > 1.
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((float)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      // Convergence.
      if (l == k) {
	// Singular value is made nonnegative.
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30) {fprintf(stderr,"Error in SVD: no convergence in 30 iterations"); exit(4);}
      // Shift from bottom 2-by-2 minor.
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      // Next QR transformation:
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	// Rotation can be arbitrary if z = 0.
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	
	for (jj=0;jj<m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  delete[](rv1);
}

// Computes (a2 + b2 )^1/2 without destructive underflow or overflow.
float pythag(float a, float b)
{
  float absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Removes the newline and other control characters at the end of a string (if present)
inline char* chomp(char* str) 
{
  if (str) 
    {
      int l;
      for (l=strlen(str)-1; l>=0 && str[l]<32; l--);
      str[++l]='\0';
    }
  return str;
}

// Returns leftmost integer in ptr and sets the pointer to first char after 
// the integer. If no integer is found, returns INT_MIN and sets pt to NULL
int strint(char*& ptr)
{
  int i;
  char* ptr0=ptr;
  if (!ptr) return INT_MIN;
  while (*ptr!='\0' && !(*ptr>='0' && *ptr<='9')) ptr++;
  if (*ptr=='\0') 
    {
      ptr=0;
      return INT_MIN;
    }
  if (*(ptr-1)=='-' && ptr>ptr0) i=-atoi(ptr); else i=atoi(ptr);
  while (*ptr!='\0' && (*ptr>='0' && *ptr<='9')) ptr++;
  return i;
}

// Count number of newline characters in file f
int CountNewlines(FILE* f)
{
  char c;
  unsigned int nl=0;
  rewind(f);
  while ((c = getc(f)) != EOF) { if (c=='\n') nl++; }
  rewind(f);
  return nl;
}

// Calculate the determinant of a 3x3 matrix M
inline float Determinant(float** M)
{
  return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) 
       + M[1][0]*(M[2][1]*M[0][2]-M[0][1]*M[2][2]) 
       + M[2][0]*(M[0][1]*M[1][2]-M[1][1]*M[0][2]); 
}

// Print values of rotation matrix and translation vector etc
void PrintTransformation(float R[3][3], float x, float y, float z)
{
  printf("Rotation matrix:          Translation vector:\n");
  printf("%6.2f  %6.2f  %6.2f    %6.2f\n",R[0][0],R[0][1],R[0][2],x);
  printf("%6.2f  %6.2f  %6.2f    %6.2f\n",R[1][0],R[1][1],R[1][2],y);
  printf("%6.2f  %6.2f  %6.2f    %6.2f\n",R[2][0],R[2][1],R[2][2],z);  
}

/////////////////////////////////////////////////////////////////////////////////////
//// Class PdbFile 
/////////////////////////////////////////////////////////////////////////////////////
class PdbFile
{
 public:
  char filename[LINELEN];   // name of corresponding file
  char rootname[LINELEN];   // name of corresponding file
  float xbar, ybar, zbar;   // mean value of x[], y[], and z[]
                            //           where l is the record index of the C_alpha for the k'th aligned pair
  float* rmsd;              // rms deviation for residue pair k
 
  PdbFile();
  ~PdbFile();
  void Read(char* infile);
  void DebugPrint();
  void Recolor(int *i, int num_pairs, float T0, float T1);
  int CleanPairs(PdbFile& pdb1, PdbFile& pdb2, int* i, int* j, int num_pairs);
  float MinimizeRMSD(PdbFile& pdb1, PdbFile& pdb2, int* i, int* j, int num_pairs, float R[3][3], float b[3]);
  float MaxsubRMSD  (PdbFile& pdb1, PdbFile& pdb2, int* i, int* j, int& num_pairs, float R[3][3], float b[3], float maxrms);
  void Superpose(PdbFile& pdb1, PdbFile& pdb2, float R[3][3], float b[3]);
  void Print(char* outfile);
  void Center(int* i, int n); // Calculate xbar, ybar, and zbar for residues indexed by i[]


 private:
  char header[LINELEN];     // pdb file header line
  char title[LINELEN];      // pdb file title line (not used at the moment)
  int num_records;          // total number of ATOM, TER, and HETATM record lines in pdb file
  char** record;            // lines containing ATOM, TER and HETATM records
  float *x, *y, *z;         // vectors of x, y, and z coordinates of atoms
  int* nres;                // residue numbers for ATOM/HETATM records in record[]
  int num_res;              // total number of residues found in pdbfile
  int* l_record;            // l_record[i] is the record index for the C_alpha of the i'th residue
  float* r[3];              // r[][k] = (x[l]-xbar, y[l]-ybar, z[l]-zbar),
};


PdbFile::PdbFile() 
{
  record=NULL;
  x=y=z=NULL;
  nres=NULL;
  l_record=NULL;
  rmsd=NULL;
  num_records=0;
  num_res = -1;
  r[0]=r[1]=r[2]=NULL;
}


PdbFile::~PdbFile() 
{
  if (v>=4) fprintf(stderr,"Executing ~PdbFile() for %s\n",filename);
  for (int l=0; l<num_records; l++) delete[](record[l]);
  delete[](record);
  delete[](x);
  delete[](y);
  delete[](z);
  delete[](nres);
  delete[](l_record);
  delete[](rmsd);
  delete[](r[0]);
  delete[](r[1]);
  delete[](r[2]);
}




// Read a pdb file into an instance of PdbFile
void PdbFile::Read(char* infile)
{
  char line[LINELEN]="";     // line read in from a file
  int l;                     // counts ATOM, HETATM, and TER lines
  char warn_maxres=0;
  FILE* inf=NULL; 

  inf = fopen(infile, "r");
  if (!inf) {fprintf(stderr,"\nError in %s: could not open %s\n\n",program,infile); exit(2);}
  l= CountNewlines(inf);
  if (v>=3) printf("Found %i newlines in %s\n",l,infile);
 
  // Assigne filename and rootname (filename without path)
  strcpy(filename,infile);   // copy name of file into object
  char* ptr=strrchr(filename,'/');
  if (ptr==NULL) ptr=filename; else ptr++;
  strcpy(rootname,ptr);
  
  // Allocate memory
  record = new(char*[l+1]);     
  nres = new(int[l+1]);        
  l_record = new(int[par.maxres-MINRES])-MINRES;        
  x = new(float[l+1]);
  y = new(float[l+1]);
  z = new(float[l+1]);
  rmsd = new(float[l+1]);
  r[0] = new(float[l+1]);
  r[1] = new(float[l+1]);
  r[2] = new(float[l+1]);
 
  // Initialize 
  for (int h=MINRES; h<par.maxres; h++) l_record[h]=-1; 
  num_res=0; 
  l=0;
  
  // Read line by line
  while (fgets(line,LINELEN-1,inf))
    {
      chomp(line); // remove newline at the end
      if (!strncmp(line,"ATOM  ",6) || !strncmp(line,"HETATM",6)) 
	{
	  // Store new line
	  record[l] = new(char[LINELEN]);
	  strcpy(record[l],line);
	  // Read nres[l], x[l],y[l],z[l]
	  if ( sscanf(line+22,"%i",&nres[l]) != 1 || sscanf(line+30,"%f %f %f",&x[l],&y[l],&z[l]) != 3 )
	    {fprintf(stderr,"\nError in %s: wrong format in ATOM/HETATM record of %s:\nline=%s\n",program,infile,line); exit(1);}
	  // Found a backbone atom to superpose?
	  if (!strncmp(line+13,"CA ",3))  // superpose only C_alpha atoms
	    {
	      if (nres[l]>=par.maxres) 
		{
		  if (v>=1 && warn_maxres==0) 
		    {
		      warn_maxres=1;
		      fprintf(stderr,"WARNING: %s contains residue indices above %i that will be ignored for the superposition\n",
			      program,infile,par.maxres);
		    }
		} else {
		  // Found an atom to superpose (C_alpha or backbone)?
		  l_record[nres[l]]=l;
		  if (nres[l]>num_res) num_res=nres[l];
		}
	    }
	  l++;
	} 
 
      else if (!strncmp(line,"TER   ",6)) 
	{
	  record[l] = new(char[LINELEN]);
	  strcpy(record[l],line);
	  l++;
	} 
      else if (!strncmp(line,"HEADER",6)) strcpy(header,line);
      else if (!strncmp(line,"TITLE ",6)) strcpy(title,line);
    }
  fclose(inf);
  num_records=l;
  return;
}



// Debugging printout
void PdbFile::DebugPrint()
{
  printf("\nLines read in from %s:\n",filename);
  for (int l=0; l<num_records; l++) printf("%3i: i=%-3i x=%-7.2f y=%-7.2f z=%-7.2f ->%s\n",l,nres[l],x[l],y[l],z[l],record[l]);
  printf("\nC_alpha coordinates read in from %s:\n",filename);
  for (int i=0; i<=num_res; i++) 
    if (l_record[i]>=0) 
      printf("%3i: %7.2f  %7.2f  %7.2f\n",i,x[l_record[i]],y[l_record[i]],z[l_record[i]]);
    else 
      printf("%3i: not defined\n",i);
  return;
}



// Recolor residues of pdbfile. 
// Residues whose index is contained in array *i get a temperature factor T1, all others get T0
void PdbFile::Recolor(int *i, int n, float T0, float T1)
{
  char* aligned=new(char[par.maxres]);
  int h,k,l;
  for (h=0; h<par.maxres; h++) aligned[h]=0;
  for (k=0; k<n; k++) aligned[i[k]]=1;
  for (l=0; l<num_records-1; l++)       // last line is TER -> no temperature factor
    {
      if (aligned[nres[l]])                // is residue in this line in aligned set i[]?
	sprintf(record[l]+60,"%6.2f",T1);
      else 
	sprintf(record[l]+60,"%6.2f",T0);
      record[l][66]=' ';  //re-replace the end-of-string symbol '\0' from sprintf with a blank
    }
  delete[](aligned);
  return;
}


// Eliminate all pairs from (i[l],j[l]) for which at least one residue is not present in its structure file
int PdbFile::CleanPairs(PdbFile& pdb1, PdbFile& pdb2, int* i, int* j, int num_pairs)
{
  int k,kk;
  for (kk=k=0; k<num_pairs; k++)
    {
      if (pdb1.l_record[i[k]]>=0 && pdb2.l_record[j[k]]>=0) 
	{
	  if (kk<k) {i[kk]=i[k]; j[kk]=j[k];}
	  kk++;
	}
    }
  if(kk<=2) {fprintf(stderr,"\nError in %s: found only %i pair%s of aligned residues which are represented in both structures %s and %s\n",program,kk,kk==1?"":"s",pdb1.filename,pdb2.filename); exit(1);}
  if (v>=3) printf("Found %i pairs of aligned residues which are represented in both structures\n",kk);
  return kk;
}



// Calculate 3D vector r[][k] by subtracting center-of-mass vector rbar = (1/n) Sum_k r_k from x[],y[],z[]
void PdbFile::Center(int* i, int num_pairs)
{
  // Calculate xbar, ybar, and zbar for residues indexed by i[] 
  xbar=ybar=zbar=0.0;
  for (int k=0; k<num_pairs; k++)
    {
      if (l_record[i[k]]<0) continue;
      int l=l_record[i[k]];
      xbar+=x[l];
      ybar+=y[l];
      zbar+=z[l];
   }
  xbar /= num_pairs;
  ybar /= num_pairs;
  zbar /= num_pairs;
  if (v>=3) printf("Center of mass vector for %s is (%.2f, %.2f, %.2f)\n",rootname,xbar,ybar,zbar);

  // Subtract center-of-mass vector from x[] y[] z[] (needs to be done only for aligned pairs)
  for (int k=0; k<num_pairs; k++)
    {
      if (l_record[i[k]]<0) continue;
      int l=l_record[i[k]];
      r[0][k] = x[l]-xbar;
      r[1][k] = y[l]-ybar;
      r[2][k] = z[l]-zbar;
    }
  return;
}



// Minimize RMSD by calculating SVD between aligned C_alpha atoms
float PdbFile::MinimizeRMSD(PdbFile& pdb1, PdbFile& pdb2, int* i, int* j, int num_pairs, float R[3][3], float b[3])
{
  float** U=new(float*[3]);
  float** V=new(float*[3]);
  for (int a=0; a<3; a++) { U[a]=new(float[3]); V[a]=new(float[3]);}
  float W[3];
  float RMSD=0;
  int k;

  // Center coordinates x[] y[] z[] of all aligned pairs with residues in both structures => xal[] yal[] zal[]
  pdb1.Center(i,num_pairs);
  pdb2.Center(j,num_pairs);
  
  // Calculate correlation matrix (U =) C = Sum_k (s_k-sbar)(r_k-rbar)^t
  for (int a=0; a<3; a++) // sum over x,y,z dimensions
    for (int b=0; b<3; b++) // sum over x,y,z dimensions
      {
	U[a][b] = 0.0;
	for (k=0; k<num_pairs; k++) 
	  U[a][b] += pdb2.r[a][k]*pdb1.r[b][k];
	U[a][b] /= num_res;
     }

  if (v>=3) 
    {
      printf("\nCorrelation matrix:\n");
      for (int a=0; a<3; a++) // sum over x,y,z dimensions
	printf("%6.2f  %6.2f  %6.2f  \n",U[a][0],U[a][1],U[a][2]);
    }
  // Do SVD of C: C = U W V^t
  /* Given a matrix a[0..m-1][0..n-1], this routine computes its singular value decomposition, A = U · W · V^t . 
     The matrix U replaces a on output. The diagonal matrix of singular values W is output as a vector w[0..n-1]. 
     The matrix V (not the transpose V^t) is output as v[0..n-1][0..n-1] . */
  SVD(U,3,3,W,V);
  
  // Calculate determinants of U and V
  W[0]=W[1]=1;                              // comment these two lines out to test SVD: R should be equal to C then
  W[2] = Determinant(U) * Determinant(V);   // comment these two lines out to test SVD
  if (v>=3) printf("\ndet(U) = %-6.3f   det(V) = %-6.3f   (should be +/- 1)\n",Determinant(U),Determinant(V));
  
  // Calculate R = U W V^t 
  for (int a=0; a<3; a++) // sum over x,y,z dimensions
    for (int b=0; b<3; b++) // sum over x,y,z dimensions
      {
	R[a][b]=0.0;
	for (int c=0; c<3; c++) // sum over x,y,z dimensions
	  R[a][b] += U[a][c]*W[c]*V[b][c];
      }

  // Calculate the RMSD deviation
  for (k=0; k<num_pairs; k++) // sum over aligned pairs
    {
      float rmsd2=0.0; // rmsd^2 deviation per C_alpha atom
      for (int a=0; a<3; a++) // sum over x,y,z dimensions
	{
	  float s = 0.0; // one of the coordinates of the transformed C_alpha vector of pdb1
	  for (int b=0; b<3; b++) 
	    s += R[a][b]*pdb1.r[b][k]; // transform s = R pdb1.r 
	  rmsd2 += (s-pdb2.r[a][k])*(s-pdb2.r[a][k]);
	}
      rmsd[k]=sqrt(rmsd2);
      RMSD += rmsd2;
      if (v>=3) printf("i[%i]=%i  j[%i]=%i  rmsd=%-6.2f\n",k,i[k],k,j[k],sqrt(rmsd2));
    }
  RMSD /= num_pairs;
  RMSD = sqrt(RMSD);

  // Screen display
  if (v>=4) PrintTransformation(R,pdb1.xbar-pdb2.xbar,pdb1.ybar-pdb2.ybar,pdb1.zbar-pdb2.zbar);
  if (v>=3) 
    {
      printf("\nNumber of aligned residues = %i\n",num_pairs); 
      printf("RMSD deviation =  %6.3f A \n",RMSD); 
    }

  for (int a=0; a<3; a++) { delete[](U[a]); delete[](V[a]);}
  delete[](U); delete[](V);
  return RMSD;
}



// Minimize RMSD for maximum subset of alignable reidues
float PdbFile::MaxsubRMSD(PdbFile& pdb1, PdbFile& pdb2, int* i, int* j, int& num_pairs, float R[3][3], float b[3], float maxrms)
{
  float RMSD=0;
  int kmax;
  float rmsd_max;
  int k,k1;

  while(1) 
    {
      // Minimize RMS by transforming pdb1 onto pdb2
      RMSD = MinimizeRMSD(pdb1,pdb2,i,j,num_pairs,R,b);
      
      // Determine atom pair kmax with largest rmsd
      rmsd_max=-1.0; kmax=0;
      for (k=0; k<num_pairs; k++) // sum over aligned pairs
	if (rmsd[k]>rmsd_max) {kmax=k; rmsd_max=rmsd[k];}
      
      // Found no more pair with too large rmsd?
      if (rmsd_max<=maxrms || num_pairs<5) break;
      
      // Eliminate kmax from aligned set
      for (k=kmax, k1=kmax+1; k1<num_pairs; k++,k1++)
	{
	  i[k]=i[k1]; 
	  j[k]=j[k1];
	}
      num_pairs--;
    }
  return RMSD;
}





// Transform the coordinates in pdb1 by rotating with matrix R and shifting with vector b 
// and write result into qtfile, together with the atom records of pdb2. Residues from 
// pdb1 will be marked by chain A, residues in pdb2 will have chain B. 
void PdbFile::Superpose(PdbFile& pdb1, PdbFile& pdb2, float R[3][3], float b[3])
{
  int l=0;  // counts records of output pdb file
  int ll;   // counts records of files pdb1 and pdb2

  // Make header and title
  strcpy(filename,"Superposition");
  strcpy(header,"HEADER    3D-superposition of ");
  strcat(header,pdb1.rootname);
  strcat(header," onto ");
  strcat(header,pdb2.rootname);

  // Set center-of-mass vector of pdb2
  float shift[3];
  shift[0] = pdb2.xbar;
  shift[1] = pdb2.ybar;
  shift[2] = pdb2.zbar;

  // Allocate memory for records etc
  num_records = pdb1.num_records + pdb2.num_records;
  record = new(char*[num_records]);

  // Write transformed pdb1 records
  for (ll=0; ll<pdb1.num_records; ll++, l++) 
    {
      record[l] = new(char[LINELEN]);
      strcpy(record[l],pdb1.record[ll]);
      record[l][21]='A'; // write chain name
      if (strncmp(record[l],"TER   ",6)) 
	for (int a=0; a<3; a++) // sum over x,y,z dimensions
	  {
	    float s = R[a][0]*(pdb1.x[ll]-pdb1.xbar) + R[a][1]*(pdb1.y[ll]-pdb1.ybar) + R[a][2]*(pdb1.z[ll]-pdb1.zbar) + shift[a]; 
	    sprintf(record[l]+30+a*8,"%8.3f",s);
	    record[l][30+(a+1)*8]=' '; // overwrite end-of-string symbol from sprintf
	  }
    }

  // Write original pdb2 record
  for (ll=0; ll<pdb2.num_records; ll++, l++) 
    {
      record[l] = new(char[LINELEN]);
      strcpy(record[l],pdb2.record[ll]);
      record[l][21]='B'; // write chain name
    }  
  return;
}



// Write PdbFile object into a file, or to stdout if outfile=NULL
void PdbFile::Print(char* outfile)
{
  FILE* outf=NULL;
  if (outfile) outf = fopen(outfile, "w"); else outf = stdout;

  if (!outf) {fprintf(stderr,"\nError in %s: could not open %s\n",program,outfile); exit(2);}
  if (!strncmp(header,"HEADER",6)) fprintf(outf,"%-80.80s\n",header);
  if (!strncmp(header,"TITLE ",6)) fprintf(outf,"%-80.80s\n",title);
  for (int l=0; l<num_records; l++) fprintf(outf,"%-80.80s\n",record[l]);
  fprintf(outf,"%-80.80s\n","END");
  if (outfile) fclose(outf);
  return;
}




/////////////////////////////////////////////////////////////////////////////////////
//// MAIN PROGRAM
/////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{  
  float maxrms=0;               // no iteration to determine maximum subset as default
  char col=0;                // no coloring as default

  int arg=1;
  int* i;                    // (i[k],j[k]) are the pairs of indices of aligned residues...
  int* j;                    // ... from in1.pdb and in2.pdb
  int k;                     // index over aligned pairs
  int num_pairs_in;          // number of aligned pairs in sequence alignment
  int num_pairs_pdb;         // number of aligned pairs present in both structures
  int num_pairs_sub;         // number of aligned pairs in maximum subset
  char line[LINELEN]="";     // line read in from a file
  char* ptr;
  PdbFile pdb1, pdb2, pdbout;
  float R[3][3];             // rotation matrix
  float b[3];                // difference vector between cetenters of mass of two sets of residues
  float RMSD;
  program=argv[0];
  
  // Read input options
  while (arg<argc && argv[arg][0]=='-')
    {
      if      (!strcmp(argv[arg],"-rms") && (arg<argc-1))  maxrms = atof(argv[++arg]); 
      else if (!strcmp(argv[arg],"-v") && (arg<argc-1)) v = atoi(argv[++arg]);
      else if (!strcmp(argv[arg],"-col"))  col = 1;
      arg++;
    }

  if (argc-arg<4) {Help(); exit(0);}
  char* alifile=argv[arg];
  char* infile1=argv[arg+1];
  char* infile2=argv[arg+2];
  char* outfile=argv[arg+3];
 
  if (v>=2) printf("%s - %s\n",program,VERSION);

  // Read alignment file
  FILE* alif=NULL;
  if (!strcmp(alifile,"-")) 
    {
      alif=stdin; 
      num_pairs_in=MAXLIN;
    }
  else 
    {
      alif = fopen(alifile, "r");
      if (!alif) {fprintf(stderr,"\nError in %s: could not open %s\n",program,alifile); exit(2);}
      num_pairs_in = CountNewlines(alif);
      if (v>=4) printf("Found %i newlines in %s\n",num_pairs_in,alifile);
    }
  i = new(int[num_pairs_in+1]);
  j = new(int[num_pairs_in+1]);
  num_pairs_in = 0;
  while (fgets(line,LINELEN,alif))
    {
      ptr = line;
      i[num_pairs_in] = strint(ptr); // read next integer, starting at ptr
      j[num_pairs_in] = strint(ptr); // read next integer, starting at ptr
      if (ptr!=NULL) num_pairs_in++; else break;
   }
  fclose(alif);
  if(num_pairs_in<=2) {fprintf(stderr,"\nError in %s: found only %i pair%s of aligned residues in %s\n",program,num_pairs_in,num_pairs_in==1?"":"s",alifile); exit(1);}
  // Debug:
  if (v>=4) for (k=0; k<num_pairs_in; k++) printf("%3i: (%i,%i)\n",k,i[k],j[k]);


  // Read pdb files
  pdb1.Read(infile1);
  if (v>=4) pdb1.DebugPrint();
  pdb2.Read(infile2);
  if (v>=4) pdb2.DebugPrint();

  // Minimize RMS
  num_pairs_pdb = pdb1.CleanPairs(pdb1,pdb2,i,j,num_pairs_in);
  if (maxrms>0) 
    {
      num_pairs_sub=num_pairs_pdb;
      RMSD = pdb1.MaxsubRMSD(pdb1,pdb2,i,j,num_pairs_sub,R,b,maxrms);
    }
  else 
    {
      RMSD = pdb1.MinimizeRMSD(pdb1,pdb2,i,j,num_pairs_pdb,R,b);
      num_pairs_sub=num_pairs_pdb;
    }

  // Recolor pdb files by changing the temperatue factor, to mark aligned residues
  if (col) 
    {
      pdb1.Recolor(i,num_pairs_sub,red,orange);
      pdb2.Recolor(j,num_pairs_sub,blue,cyan);
    }
  //  pdb1.Print(NULL); // print to screen


  // Transform coordinates of r -> t = R r + b
  pdbout.Superpose(pdb1,pdb2,R,b);

  // Screen display
  if (v>=2) 
    {
      printf("Number of aligned residue pairs in sequence alignment     = %i\n",num_pairs_in); 
      printf("Number of aligned residue pairs present in both pdb files = %i\n",num_pairs_pdb); 
      if (maxrms) 
	  printf("Number of residue pairs in maximum subset (RMSD < %.2f A) = %i\n",maxrms,num_pairs_sub); 
      printf("RMSD deviation =  %6.3f A \n",RMSD); 
      printf("\n");
      PrintTransformation(R,pdb1.xbar-pdb2.xbar,pdb1.ybar-pdb2.ybar,pdb1.zbar-pdb2.zbar);
    }
      
  // Write coordinates s and transformed coordinates t of r into output file
  if (!strcmp(outfile,"-")) pdbout.Print(NULL); else pdbout.Print(outfile);

  delete[](i);
  delete[](j);
  exit(0);
}
