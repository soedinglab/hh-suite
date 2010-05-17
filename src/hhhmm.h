// hhhmm.h

class CSCounts;

class HMM
{
public:
    HMM(int maxseqdis=MAXSEQDIS, int maxres=MAXRES);
    ~HMM();
    HMM& operator=(HMM&);

    int n_display;            // number of sequences stored for display of alignment (INCLUDING >ss_ and >cf_ sequences)
    char** sname;             // names of stored sequences
    char** seq;               // residues of stored sequences (first at pos 1!)
    int ncons;                // index of consensus sequence
    int nfirst;               // index of first sequence (query sequence of HMM)
    int nss_dssp;             // index of seq[] with secondary structure by dssp
    int nsa_dssp;             // index of seq[] with solvent accessibility by dssp
    int nss_pred;             // index of seq[] with predicted secondary structure
    int nss_conf;             // index of seq[] with confidence values for secondary structure prediction

    int L;                    // length of HMM = number of match states; set in declaration of HMM object
    int N_in;                 // number of sequences in alignment
    int N_filtered;           // number of sequences after filtering
    float* Neff_M;            // Neff_M[i] = diversity of subalignment of seqs that have residue in col i
    float* Neff_I;            // Neff_I[i] = diversity of subalignment of seqs that have insert in col i
    float* Neff_D;            // Neff_D[i] = diversity of subalignment of seqs that have delete in col i
    float Neff_HMM;           // average number of Neff over total length of HMM

    char* longname;           // Full name of first sequence of original alignment (NAME field)
    char name[NAMELEN];       // HMM name = first word in longname in lower case
    char file[NAMELEN];       // Basename (with path, without extension) of alignment file that was used to construct the HMM
    char fam[NAMELEN];        // family ID (derived from name) (FAM field)
    char sfam[NAMELEN];       // superfamily ID (derived from name)
    char fold[NAMELEN];       // fold ID (derived from name)
    char cl[NAMELEN];         // class ID (derived from name)

    float lamda, mu;          // coefficients for aa score distribution of HMM using parameters in 'Parameters par'
    bool has_pseudocounts;    // set to true if HMM contains pseudocounts

   // Make a flat copy of q
    void FlatCopyTo(HMM& t);

    // Read an HMM from a HHsearch .hhm file and return 0 at end of file
    int Read(FILE* dbf, char* path=NULL);

    // Read an HMM from a HMMer .hmm file; return 0 at end of file
    int ReadHMMer(FILE* dbf, char* filestr=NULL);

    // Add transition pseudocounts to HMM
    void AddTransitionPseudocounts(float gapd=par.gapd, float gape=par.gape, float gapf=par.gapf, float gapg=par.gapg, float gaph=par.gaph, float gapi=par.gapi, float gapb=par.gapb);

    // Generate an amino acid frequency matrix g[][] with full pseudocount admixture (tau=1)
    void PreparePseudocounts();

    // Add context specific amino acid pseudocounts to HMM: t.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
    void AddContextSpecificPseudocounts(char pcm=par.pcm, float pca=par.pca, float pcb=par.pcb, float pcc=par.pcc);

    // Fill CountProfile with HMM-counts for CS pseudocount calculation
    void fillCountProfile(cs::CountProfile<cs::AA> *csProfile);

    // Add amino acid pseudocounts to HMM: t.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
    void AddAminoAcidPseudocounts(char pcm=par.pcm, float pca=par.pca, float pcb=par.pcb, float pcc=par.pcc);

    // Calculate amino acid backround frequencies for HMM
    void CalculateAminoAcidBackground();

    // Add no amino acid pseudocounts to HMM: copy  t.p[i][a] = f[i][a]
    void NoAminoAcidPseudocounts() {for(int i=1; i<=L; i++) for(int a=0; a<20; a++) p[i][a]=f[i][a];};

    // Factor Null model into HMM t
    void IncludeNullModelInHMM(HMM& q, HMM& t, int columnscore=par.columnscore);

    // Write HMM to output file
    void WriteToFile(char* outfile);

    // Insert calibration line 'EVD   lamda   mu      hashvalue' into HMM file
    void InsertCalibration(char* infile);

    // Transform log to lin transition probs
    void Log2LinTransitionProbs(float beta=1.0);

    // Set query columns in His-tags etc to Null model distribution
    void NeutralizeTags();

    // Calculate effective number of sequences using profiles INCLUDING pseudocounts
    float CalcNeff();

    // Add secondary structure prediction to HMM
    void AddSSPrediction(char seq_pred[], char seq_conf[]);

    // Initialize f[i][a] with query HMM
    void MergeQueryHMM(HMM& q, float wk[]);

    // Rescale rate matrices P[a][b], R[a][b] according to HMM av. aa composition in pav[a]
    void RescaleMatrix();

    // Needed for SSE2 prefiltering with HHblits with amino acid alphabet
    float** p;                // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts

private:
    float** f;                // f[i][a] = prob of finding amino acid a in column i WITHOUT pseudocounts
    float** g;                // g[i][a] = prob of finding amino acid a in column i WITH pseudocounts
    //float** p;                // p[i][a] = prob of finding amino acid a in column i WITH OPTIMUM pseudocounts
    float** tr;               // tr[i][X2Y] = log2 of transition probabilities M2M M2I M2D I2M I2I D2M D2D

    char* ss_dssp;            // secondary structure determined by dssp 0:-  1:H  2:E  3:C  4:S  5:T  6:G  7:B
    char* sa_dssp;            // solvent accessibility state determined by dssp 0:-  1:A (absolutely buried) 2:B  3:C  4:D  5:E (exposed)
    char* ss_pred;            // predicted secondary structure          0:-  1:H  2:E  3:C
    char* ss_conf;            // confidence value of prediction         0:-  1:0 ... 10:9
    char* Xcons;              // consensus sequence in internal representation (A=0 R=1 N=2 D=3 ...)
    float pav[NAA];           // pav[a] = average freq of amino acids in HMM (including subst matrix pseudocounts)
    float pnul[NAA];          // null model probabilities used in comparison (only set in template/db HMMs)
    int* l;                   // l[i] = pos. of j'th match state in aligment
    char dont_delete_seqs;    // set to one if flat copy of seqs and sname was made to a hit object, to avoid deletion
    char trans_lin;           // transition probs are given in log or lin space? (0: p_tr  1: log(p_tr)

    // Utility for Read()
    int Warning(FILE* dbf, char line[], char name[])
    {
        if (v) cerr<<"\nWARNING: could not read line\n\'"<<line<<"\'\nin HMM "<<name<<" in "<<file<<"\n";
        while (fgetline(line,LINELEN,dbf) && !(line[0]=='/' && line[1]=='/'));
        if (line) return 2;  //return status: skip HMM
        return 0;            //return status: end of database file
    }

    friend class Hit;
    friend class Alignment;
    friend class CSCounts;
};
