// hhfullalignment.C

#include "hhhalfalignment.h"

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
// Methods of class HalfAlignment
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
  


/////////////////////////////////////////////////////////////////////////////////////
// Constructor
HalfAlignment::HalfAlignment(int maxseqdis)
{
  n=0; 
  sname=seq=NULL; 
  nss_dssp = nss_pred = nss_conf = nsa_dssp = ncons= -1;
  h = new int[maxseqdis];   //h[k] = next position of sequence k to be written
  s = new char*[maxseqdis];  //s[k][h] = character in column h, sequence k of output alignment
  l = new int*[maxseqdis];   //counts non-gap residues: l[k][i] = index of last residue AT OR BEFORE match state i in seq k
  m = new int*[maxseqdis];   //counts positions:        m[k][i] = position of match state i in string seq[k]
  if (!h || !s || !l || !m) MemoryError("space for formatting HMM-HMM alignment", __FILE__, __LINE__, __func__);
}

/////////////////////////////////////////////////////////////////////////////////////
// Destructor
HalfAlignment::~HalfAlignment()
{
  Unset();
  delete[] h;
  delete[] s;
  delete[] l;
  delete[] m;
}



/////////////////////////////////////////////////////////////////////////////////////
// Free memory in HalfAlignment arrays s[][], l[][], and m[][]
void HalfAlignment::Unset()
{
   // Free memory for alignment characters and residue counts
  for (int k=0; k<n; k++) 
    {
      delete[] s[k]; 
      delete[] l[k];
      delete[] m[k];
    }
  n=0; 
  sname=seq=NULL; 
  nss_dssp = nss_pred = nss_conf = nsa_dssp = ncons= -1;
}


/////////////////////////////////////////////////////////////////////////////////////
// Prepare a2m/a3m alignment: Calculate l[k][i] (residue indices) and m[k][i] (position in seq[k])
void HalfAlignment::Set(char* name, char** seq_in, char** sname_in, int n_in, int L_in, int n1, int n2, int n3, int n4, int nc)
{
  int i;     //counts match states in seq[k]
  int ll;    //counts residues LEFT from or at current position in seq[k]
  int mm;    //counts postions in string seq[k]
  int k;     //counts sequences
  char c;
  char warned=0;

  nss_dssp=n1; nss_pred=n2; nss_conf=n3; nsa_dssp=n4; ncons=nc;
  seq=seq_in;     //flat copy of sequences
  sname=sname_in; //flat copy of sequence names
  n=n_in;
  L=L_in;    
  pos=0;

  // Allocate memory for alignment characters and residue counts
  for (k=0; k<n; k++) 
    {
      s[k]=new char[LINELEN];
      l[k]=new int[L+10]; 
      m[k]=new int[L+10];
      if (!s[k] || !l[k] || !m[k]) MemoryError("space for formatting HMM-HMM alignment", __FILE__, __LINE__, __func__);
      h[k]=0; //starting positions in alignment = 0
    }

  for (k=0; k<n; k++)
    {
      m[k][0]=0;  // 0'th match state (virtual) is begin state at 0
      //if k is consensus sequence 
      if (k==nc) {
	for (i=1; i<=L; i++) m[k][i]=l[k][i]=i; 
	m[k][L+1]=l[k][L+1]=L; 
	continue;
      }
      i=1; mm=1; ll=1;
      while ((c=seq[k][mm]))
	{
	  if (MatchChr(c)==c)    //count match/delete states
	    {
	      l[k][i]=ll;
	      m[k][i]=mm;
	      i++;
	    }
	  if (WordChr(c)) ll++;  //index of next residue
	  mm++;
	}
      l[k][i]=ll-1; //set l[k][L+1] eq number of residues in seq k (-1 since there is no residue at L+1st match state)
      m[k][i]=mm;   //set m[k][L+1]
     if ((i-1)!=L && !warned) 
	{
	  std::cerr<<"WARNING: sequence "<<sname[k]<<" in HMM "<<name<<" has "<<i<<" match states but should have "<<L<<"\n";
	  warned=1;
	}
    }
  //DEBUG
  if (v>=5)
    {
      fprintf(stderr,"  i chr   m   l\n");
      for(i=0;i<=L+1;i++) fprintf(stderr,"%3i   %1c %3i %3i\n",i,seq[0][m[0][i]],m[0][i],l[0][i]);
      printf("\n");
    }
}


/////////////////////////////////////////////////////////////////////////////////////
// Fill in insert states following match state i (without inserting '.' to fill up)
void HalfAlignment::AddInserts(int i)
{
  for (int k=0; k<n; k++)                        // for all sequences...
    for (int mm=m[k][i]+1; mm<m[k][i+1]; mm++)   // for all inserts between match state i and i+1...
      s[k][h[k]++]=seq[k][mm];                   // fill inserts into output alignment s[k]
}

/////////////////////////////////////////////////////////////////////////////////////
// Fill up alignment with gaps '.' to generate flush end (all h[k] equal)
void HalfAlignment::FillUpGaps()
{ 
  int k;      //counts sequences
  pos=0;

  // Determine max position h[k]
  for (k=0; k<n; k++) pos = imax(h[k],pos);
  
  // Fill in gaps up to pos
  for (k=0; k<n; k++) 
    {
      for (int hh=h[k]; hh<pos; hh++) s[k][hh]='.';
      h[k]=pos;
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Fill in insert states following match state i and fill up gaps with '.' 
void HalfAlignment::AddInsertsAndFillUpGaps(int i) 
{ 
  AddInserts(i); 
  FillUpGaps(); 
}

/////////////////////////////////////////////////////////////////////////////////////
// Add gap column '.'
void HalfAlignment::AddChar(char c)  
{ 
  for (int k=0; k<n; k++) s[k][h[k]++]=c;                
  pos++;
}

/////////////////////////////////////////////////////////////////////////////////////
// Add match state column i as is
void HalfAlignment::AddColumn(int i) 
{ 
  for (int k=0; k<n; k++) s[k][h[k]++]=seq[k][m[k][i]];  
  pos++;
}

/////////////////////////////////////////////////////////////////////////////////////
// Add match state column i as insert state
void HalfAlignment::AddColumnAsInsert(int i) 
{ 
  char c; 
  for (int k=0; k<n; k++) 
    if ((c=seq[k][m[k][i]])!='-' && (c<'0' || c>'9')) 
      s[k][h[k]++]=InsertChr(c); 
  pos++;
}

/////////////////////////////////////////////////////////////////////////////////////
// Remove all characters c from template sequences
void HalfAlignment::RemoveChars(char c)
{ 
  int k,h,hh;
  for (k=0; k<n; k++)
    {
      for (h=hh=0; h<pos; h++)
	if (s[k][h]!=c) s[k][hh++]=s[k][h];
      s[k][++hh]='\0';
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Transform alignment sequences from A3M to A2M (insert ".")
void HalfAlignment::BuildFASTA()
{
  AddInserts(0); 
  FillUpGaps();
  for (int i=1; i<=L; i++)
    {
      AddColumn(i); 
      AddInserts(i); 
      FillUpGaps();
    }
  ToFASTA();
}

/////////////////////////////////////////////////////////////////////////////////////
// Transform alignment sequences from A3M to A2M (insert ".")
void HalfAlignment::BuildA2M()
{
  AddInserts(0); 
  FillUpGaps();
  for (int i=1; i<=L; i++)
    {
      AddColumn(i); 
      AddInserts(i); 
      FillUpGaps();
    }
  AddChar('\0');
}

/////////////////////////////////////////////////////////////////////////////////////
// Transform alignment sequences from A3M to A2M (insert ".")
void HalfAlignment::BuildA3M()
{
  AddInserts(0); 
  for (int i=1; i<=L; i++)
    {
      AddColumn(i); 
      AddInserts(i); 
    }
  AddChar('\0');
}

/////////////////////////////////////////////////////////////////////////////////////
// Transform alignment sequences from A2M to FASTA ( lowercase to uppercase and '.' to '-')
void HalfAlignment::ToFASTA()
{
  for (int k=0; k<n; k++)
    {
      uprstr(s[k]);
      strtr(s[k],".","-");
    }
}

/////////////////////////////////////////////////////////////////////////////////////
// Align query (HalfAlignment) to template (i.e. hit) match state structure
/////////////////////////////////////////////////////////////////////////////////////
void HalfAlignment::AlignToTemplate(Hit& hit, const char outformat)
{
  int i,j;
  int step;    // column of the HMM-HMM alignment (first:nstep, last:1)
  char state;

  if(0) {  //par.loc==0) { //////////////////////////////////////////// STILL NEEDED??
    // If in global mode: Add part of alignment before first MM state
    AddInserts(0); // Fill in insert states before first match state
    for (i=1; i<hit.i[hit.nsteps]; i++)
      {
	AddColumnAsInsert(i);
	AddInserts(i);
	if (outformat<=2) FillUpGaps();
      }
  }

  // Add endgaps (First state must be an MM state!!)
  for (j=1; j<hit.j[hit.nsteps]; j++)    
    {
      AddChar('-');
    }

  // Add alignment between first and last MM state
  for (step=hit.nsteps; step>=1; step--) 
  {
    state = hit.states[step];
    i = hit.i[step];

    switch(state)
      {
      case MM:  //MM pair state (both query and template in Match state)
	AddColumn(i);
	AddInserts(i);
	break;
      case DG: //D- state
      case MI: //MI state
	AddColumnAsInsert(i);
	AddInserts(i);
	break;
      case GD: //-D state
      case IM: //IM state
	AddChar('-');
	break;
      }
    if (outformat<=2) FillUpGaps();

  }

  if(0) { //par.loc==0) { //////////////////////////////////////////// STILL NEEDED??

    // If in global mode: Add part of alignment after last MM state
    for (i=hit.i[1]+1; i<=L; i++)    
      {
	AddColumnAsInsert(i);
	AddInserts(i);
	if (outformat==2) FillUpGaps();
      }
  }

  // Add endgaps 
  for (j=hit.j[1]+1; j<=hit.L; j++)    
    {
      AddChar('-');
    }

  // Add end-of-string character
  AddChar('\0');
}


/////////////////////////////////////////////////////////////////////////////////////
// Write the a2m/a3m alignment into alnfile 
/////////////////////////////////////////////////////////////////////////////////////
void HalfAlignment::Print(char* alnfile, const char append, char* commentname, const char format[])
{
  int k;      //counts sequences
  int omitted=0; // counts number of sequences with no residues in match states
  FILE *outf;
  char* tmp_name = new char[NAMELEN];
  if (strcmp(alnfile,"stdout"))
    {
      if (append) outf=fopen(alnfile,"a"); else outf=fopen(alnfile,"w");
      if (!outf) OpenFileError(alnfile, __FILE__, __LINE__, __func__);
    } 
  else
    outf = stdout;
  if (v>=3) std::cout<<"Writing alignment to "<<alnfile<<"\n";

  if (!format || strcmp(format,"psi"))
    {
      if (commentname != NULL) fprintf(outf,"#%s\n",commentname);

      for (k=0; k<n; k++)
	if (k==nss_pred || k==nss_conf || k==nss_dssp || k==nsa_dssp)
	  {
	    fprintf(outf,">%s\n",sname[k]);
	    fprintf(outf,"%s\n",s[k]);
	  }
      for (k=0; k<n; k++)
	{
	  if (!(k==nss_pred || k==nss_conf || k==nss_dssp || k==nsa_dssp))
	    {
	      // Print sequence only if it contains at least one residue in a match state
	      if (1) //strpbrk(s[k],"ABCDEFGHIKLMNPQRSTUVWXYZ1234567890")) 
		{
		  fprintf(outf,">%s\n",sname[k]);
		  fprintf(outf,"%s\n",s[k]);
		} else {
		omitted++;
		if (v>=3) printf("%-14.14s contains no residue in match state. Omitting sequence\n",sname[k]);
	      }
	    }
	}
      if (v>=2 && omitted) printf("Omitted %i sequences in %s which contained no residue in match state\n",omitted,alnfile);
    }
  else
    {
      for (k=0; k<n; k++)
          {
	    strwrd(tmp_name,sname[k],NAMELEN);
            fprintf(outf,"%-20.20s ",tmp_name);
            char* ptr=s[k];
            for (; *ptr!='\0'; ptr++)
              if (*ptr==45 || (*ptr>=65 && *ptr<=90)) fprintf(outf,"%c",*ptr);
            fprintf(outf,"\n");
          }
    }
  fclose(outf);
  delete[] tmp_name;
}
