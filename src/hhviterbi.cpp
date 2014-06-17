// hhviterbifilter.C
//
// Functions for Viterbi filter in HHblits

#ifndef HHVITERBI_c
#define HHVITERBI_c

#include "hhviterbi.h"
#include "util.h"
#include "hhutil.h"

// static
void Viterbi::ExcludeAlignment(ViterbiMatrix * matrix,HMMSimd* q_four, HMMSimd* t_four,int elem,
                               int * i_steps, int * j_steps, int nsteps){
    const HMM * q = (const HMM *) q_four->GetHMM(elem);
    const HMM * t = (const HMM *) t_four->GetHMM(elem);
    // Exclude cells in direct neighbourhood from all further alignments
    for(int step=1;step < nsteps;step++){
        int i=i_steps[step];
        int j=j_steps[step];
        
        for (int ii=imax(i-2,1); ii<=imin(i+2,q->L); ++ii){
            matrix->setCellOff(ii, j, elem, true);
        }
        for (int jj=imax(j-2,1); jj<=imin(j+2,t->L); ++jj){
            matrix->setCellOff(i, jj, elem, true);
        }
    }
}


Viterbi::Viterbi(int max_seq_length,bool local,float penalty_gap_query,float penalty_gap_template, float correlation, int par_min_overlap, float shift){
    this->max_seq_length=max_seq_length;
    this->local = local;
    this->penalty_gap_query = penalty_gap_query;
    this->penalty_gap_template = penalty_gap_template;
    sMM_DG_MI_GD_IM_vec = (simd_float *) malloc_simd_float(VEC_SIZE*max_seq_length*5*sizeof(float));
    this->correlation = correlation;
    this->par_min_overlap = par_min_overlap;
//    this->exclstr = new char[strlen(exclstr)+1];
//    strcpy(this->exclstr, exclstr);
    this->shift = shift;
}

Viterbi::~Viterbi(){
    free(sMM_DG_MI_GD_IM_vec);
//    delete exclstr;
}







/////////////////////////////////////////////////////////////////////////////////////
// Trace back Viterbi alignment of two profiles based on matrices bXX[][]
/////////////////////////////////////////////////////////////////////////////////////
//TODO inline
Viterbi::BacktraceResult Viterbi::Backtrace(ViterbiMatrix * matrix,int elem,int start_i[VEC_SIZE], int start_j[VEC_SIZE])
{
    // Trace back trough the matrices bXY[i][j] until first match state is found (STOP-state)
    int step;      // counts steps in path through 5-layered dynamic programming matrix
    int i,j;       // query and template match state indices
    //    InitializeBacktrace(q,t);
    const int maxAlignmentLength = start_i[elem]+start_j[elem]+2;
    int * i_steps  = new int[maxAlignmentLength];
    int * j_steps  = new int[maxAlignmentLength];
    char * states  = new char[maxAlignmentLength];
    


    // Back-tracing loop
    int matched_cols=0;         // for each MACTH (or STOP) state matched_col is incremented by 1
    step=0;                 // steps through the matrix correspond to alignment columns (from 1 to nsteps)
    char state=ViterbiMatrix::MM;           // state with maximum score must be MM state  // already set at the end of Viterbi()
    i=start_i[elem]; j=start_j[elem];   // last aligned pair is (i2,j2)
    
    while (state!=ViterbiMatrix::STOP)     // while (state!=STOP)  because STOP=0
    {
        step++;
        states[step] = state;
        i_steps[step] = i;
        j_steps[step] = j;
        
        switch (state)
        {
            case ViterbiMatrix::MM: // current state is MM, previous state is bMM[i][j]
                matched_cols++;
                state = (i <= 1||j<=1) ? ViterbiMatrix::STOP : matrix->getMatMat(i--,j--,elem);
                break;
            case ViterbiMatrix::GD: // current state is GD
                if(j<=1)
                    state = ViterbiMatrix::STOP;
                else if (matrix->getGapDel(i, j--, elem) == true)
                    state = ViterbiMatrix::MM;    // previous state is Match state
                break;                            // default: previous state is same state (GD)
            case ViterbiMatrix::IM:
                if(j<=1)
                    state = ViterbiMatrix::STOP;
                else if (matrix->getInsMat(i, j--, elem) == true)
                    state = ViterbiMatrix::MM;    // previous state is Match state
                break;                            // default: previous state is same state (IM)
            case ViterbiMatrix::DG:
                if(i<=1)
                    state = ViterbiMatrix::STOP;
                else if (matrix->getDelGap(i--, j, elem) == true)
                    state = ViterbiMatrix::MM;    // previous state is Match state
                break;                            // default: previous state is same state (DG)
            case ViterbiMatrix::MI:
                if(i<=1)
                    state=ViterbiMatrix::STOP;
                else if (matrix->getMatIns(i--, j, elem) == true)
                    state = ViterbiMatrix::MM;    // previous state is Match state
                break;                            // default: previous state is same state (MI)
            default:
                fprintf(stderr,"Error: unallowed state value %i occurred during backtracing at step %i, (i,j)=(%i,%i)\n",state,step,i,j);
                state=ViterbiMatrix::STOP;
                //    v=4;
                break;
        } //end switch (state)
    } //end while (state)
    
    states[step] = ViterbiMatrix::MM;  // first state (STOP state) is set to MM state
    int nsteps=step;
    
    Viterbi::BacktraceResult result;
    
    result.i_steps = i_steps;
    result.j_steps = j_steps;
    result.states = states;
    result.count = nsteps;
    result.matched_cols = matched_cols;

    return result;
}

//TODO: inline
Viterbi::ViterbiResult* Viterbi::Align(HMMSimd* q, HMMSimd* t,ViterbiMatrix * viterbiMatrix, int maxres){
    Viterbi::ViterbiResult* result = new Viterbi::ViterbiResult();
    if (viterbiMatrix->hasCellOff()==true) {
        this->AlignWithCellOff(q,t,viterbiMatrix, maxres, result);
    } else {
        this->AlignWithOutCellOff(q,t,viterbiMatrix, maxres, result);
    }
    viterbiMatrix->setCellOff(false); // the ViterbiAlign set all Cell of values to false

    return result;
}


//TODO: inline
Viterbi::BacktraceScore Viterbi::ScoreForBacktrace(HMMSimd* q_four, HMMSimd* t_four,
                                                          int elem,Viterbi::BacktraceResult * backtraceResult,
                                                          float alignmentScore[VEC_SIZE],
                                                          int ssm1,int ssm2)
{
    
    // Allocate new space for alignment scores
    const HMM * q = (const HMM *) q_four->GetHMM(elem);
    const HMM * t = (const HMM *) t_four->GetHMM(elem);
    char * states=backtraceResult->states;
    int * i_steps=backtraceResult->i_steps;
    int * j_steps=backtraceResult->j_steps;
    int nsteps=backtraceResult->count;
    float * S=new float[nsteps+1];
    float * S_ss=new float[nsteps+1];
    if (!S_ss) MemoryError("space for HMM-HMM alignments", __FILE__, __LINE__, __func__);
    
    // Add contribution from secondary structure score, record score a long alignment,
    // and record template consensus sequence in master-slave-alignment to query sequence
    
    float score_ss=0.0f;
    float score=alignmentScore[elem];
    float score_sort = 0.0f;
    float score_aass = 0.0f;
    float Pvalt = 1.0f;
    float logPvalt = 0.0f;
    
    int ssm=ssm1+ssm2;
    for (int step=1; step<=nsteps; step++)
    {
        switch(states[step])
        {
            case ViterbiMatrix::MM:
                S[step]    = Score(q->p[i_steps[step]],t->p[j_steps[step]]);
                S_ss[step] = ScoreSS(q,t,i_steps[step],j_steps[step],ssm);
                score_ss += S_ss[step];
                break;
            case ViterbiMatrix::MI: //if gap in template
            case ViterbiMatrix::DG:
            default: //if gap in T or Q
                S[step]=S_ss[step]=0.0f;
                break;
        }
    }
    
    if (ssm2>=1) score-=score_ss;    // subtract SS score added during alignment!!!!
    //    printf("###New score %f",score);
    // Add contribution from correlation of neighboring columns to score
    float Scorr=0;
    if (nsteps)
    {
        for (int step=2; step<=nsteps; step++) Scorr+=S[step]*S[step-1];
        for (int step=3; step<=nsteps; step++) Scorr+=S[step]*S[step-2];
        for (int step=4; step<=nsteps; step++) Scorr+=S[step]*S[step-3];
        for (int step=5; step<=nsteps; step++) Scorr+=S[step]*S[step-4];
        score+=correlation*Scorr;
    }
    
    // Set score, P-value etc.
    score_sort = score_aass = -score;
    if (t->mu)
    {
        logPvalt=logPvalue(score,t->lamda,t->mu);
        Pvalt   =Pvalue(score,t->lamda,t->mu);
    }
    else { logPvalt=0; Pvalt=1;}
    //   printf("%-10.10s lamda=%-9f  score=%-9f  logPval=%-9g\n",name,t->lamda,score,logPvalt);
    //DEBUG: Print out Viterbi path
    
    Viterbi::BacktraceScore backtraceScore;
    backtraceScore.score_ss=score_ss;
    backtraceScore.score=score;
    backtraceScore.score_sort=score_sort;
    backtraceScore.score_aass=score_aass;
    backtraceScore.Pvalt=Pvalt;
    backtraceScore.logPvalt=logPvalt;
    backtraceScore.S=S;
    backtraceScore.S_ss=S_ss;
    
    if (Log::reporting_level() >= LogLevel::DEBUG1) {
        Viterbi::PrintDebug(q,t,&backtraceScore,backtraceResult,ssm);
    }
    
    return backtraceScore;
}


void Viterbi::PrintDebug(const HMM * q,const HMM *t,Viterbi::BacktraceScore * backtraceScore,Viterbi::BacktraceResult * backtraceResult,
                         const int ssm){
    int nfirst=0;
    char * states=backtraceResult->states;
    int * i_steps=backtraceResult->i_steps;
    int * j_steps=backtraceResult->j_steps;
    int nsteps=backtraceResult->count;
    
    printf("score=%7.3f  score_ss=%7.3f\n",backtraceScore->score,backtraceScore->score_ss);
    printf("step  Q T    i    j  state   score    T Q cf ss-score\n");
    for (int step=nsteps; step>=1; step--)
    {
        switch(states[step])
        {
            case ViterbiMatrix::MM:
                printf("%4i  %1c %1c ",step,q->seq[q->nfirst][i_steps[step]],t->seq[nfirst][j_steps[step]]);
                break;
            case ViterbiMatrix::GD:
            case ViterbiMatrix::IM:
                printf("%4i  - %1c ",step,t->seq[nfirst][j_steps[step]]);
                break;
            case ViterbiMatrix::DG:
            case ViterbiMatrix::MI:
                printf("%4i  %1c - ",step,q->seq[q->nfirst][i_steps[step]]);
                break;
        }
        printf("%4i %4i     %2i %7.2f    ",i_steps[step],j_steps[step],(int)states[step],Score(q->p[i_steps[step]],t->p[j_steps[step]]));
        printf("%c %c %1i %7.2f\n",i2ss(t->ss_dssp[j_steps[step]]),i2ss(q->ss_pred[i_steps[step]]),q->ss_conf[i_steps[step]]-1,ScoreSS(q,t,i_steps[step],j_steps[step],ssm));
    }
}


/////////////////////////////////////////////////////////////////////////////////////
//// Functions that set cell off
/////////////////////////////////////////////////////////////////////////////////////
void Viterbi::InitializeForAlignment(HMM* q, HMM* t, ViterbiMatrix * matrix, int elem, bool self, int par_min_overlap)
{
    int i,j;
    int min_overlap;
    
    if (self)
    {
        // Cross out cells in lower diagonal for self-comparison?
        for (i=1; i<=q->L; ++i)
        {
            int jmax = imin(i+SELFEXCL,t->L);
            for (j=1; j<=jmax; ++j)
                matrix->setCellOff(i, j, elem, true);   // cross out cell near diagonal
            for (j=jmax+1; j<=t->L+1; ++j)
                matrix->setCellOff(i, j, elem, false);   // no other cells crossed out yet
        }
    }
    else
    {
        
        // Compare two different HMMs Q and T
        // Activate all cells in dynamic programming matrix
        for (i=1; i<=q->L; ++i)
            for (j=1; j<=t->L; ++j)
                matrix->setCellOff(i, j, elem, false);  // no other cells crossed out yet
        
        // Cross out cells that are excluded by the minimum-overlap criterion
        if (par_min_overlap==0)
            min_overlap = imin(60, (int)(0.333f*imin(q->L,t->L))+1); // automatic minimum overlap
        else
            min_overlap = imin(par_min_overlap, (int)(0.8f*imin(q->L, t->L)));
        
        for (i=0; i<min_overlap; ++i)
            for (j=i-min_overlap+t->L+1; j<=t->L; ++j) // Lt-j+i>=Ovlap => j<=i-Ovlap+Lt => jmax=min{Lt,i-Ovlap+Lt}
                matrix->setCellOff(i, j, elem, true);
        for (i=q->L-min_overlap+1; i<=q->L; ++i)
            for (j=1; j<i+min_overlap-q->L; ++j)      // Lq-i+j>=Ovlap => j>=i+Ovlap-Lq => jmin=max{1, i+Ovlap-Lq}
                matrix->setCellOff(i, j, elem, true);
        
//        // Cross out rows which are contained in range given by exclstr ("3-57,238-314")
//        if (exclstr)
//        {
//            char* ptr=exclstr;
//            int i0, i1;
//            while (1)
//            {
//                i0 = abs(strint(ptr));
//                i1 = abs(strint(ptr));
//                if (!ptr) break;
//                for (i=i0; i<=imin(i1,q->L); ++i)
//                    for (j=1; j<=t->L; ++j)
//                        matrix->setCellOff(i, j, elem, true);
//            }
//        }
    }
}


////////////////////////////////////////////////////////////////////////
#endif
