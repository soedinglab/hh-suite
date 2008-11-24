// hhfunc.C

// Read input file (HMM, HHM, or alignment format), and add pseudocounts etc.
void ReadAndPrepare(char* infile, HMM& q, Alignment* qali=NULL)
{
    char path[NAMELEN];

    // Open query file and determine file type
    char line[LINELEN]=""; // input line
    FILE* inf=NULL;
    if (strcmp(infile,"stdin"))
    {
        inf = fopen(infile, "r");
        if (!inf) OpenFileError(infile);
        Pathname(path,infile);
    }
    else
    {
        inf = stdin;
        if (v>=2) printf("Reading HMM / multiple alignment from standard input ...\n(To get a help list instead, quit and type %s -h.)\n",program_name);
        *path='\0';
    }

    fgetline(line,LINELEN-1,inf);

    // Is infile a HMMER file?
    if (!strncmp(line,"HMMER",5))
    {
        if (v>=2) cout<<"Query file is in HMMER format\n";

        // Read 'query HMMER file
        q.ReadHMMer(inf,path);
        if (v>=2 && q.Neff_HMM>11.0)
            fprintf(stderr,"WARNING: HMM %s looks too diverse (Neff=%.1f>11). Better check the underlying alignment... \n",q.name,q.Neff_HMM);

        // Don't add transition pseudocounts to query!!

        if (!*par.clusterfile) { //compute context-specific pseudocounts?
            // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
            q.PreparePseudocounts();
        } else {
            // Generate an amino acid frequency matrix from f[i][a] with full context specific pseudocount admixture (tau=1) -> g[i][a]
            q.PrepareContextSpecificPseudocounts();
        }

        // DON'T ADD amino acid pseudocounts to query: pcm=0!  q.p[i][a] = f[i][a]
        q.AddAminoAcidPseudocounts(0, par.pca, par.pcb, par.pcc);
        q.CalculateAminoAcidBackground();
    }

    // ... or is it an hhm file?
    else if (!strncmp(line,"NAME",4) || !strncmp(line,"HH",2))
    {
        if (v>=2) cout<<"Query file is in HHM format\n";

        // Rewind to beginning of line and read query hhm file
        rewind(inf);
        q.Read(inf,path);
        if (v>=2 && q.Neff_HMM>11.0)
            fprintf(stderr,"WARNING: HMM %s looks too diverse (Neff=%.1f>11). Better check the underlying alignment... \n",q.name,q.Neff_HMM);

        // Add transition pseudocounts to query -> q.p[i][a]
        q.AddTransitionPseudocounts();

        if (!*par.clusterfile) { //compute context-specific pseudocounts?
            // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
            q.PreparePseudocounts();
        } else {
            // Generate an amino acid frequency matrix from f[i][a] with full context specific pseudocount admixture (tau=1) -> g[i][a]
            q.PrepareContextSpecificPseudocounts();
        }

        // Add amino acid pseudocounts to query:  q.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
        q.AddAminoAcidPseudocounts(!q.has_pseudocounts ? par.pcm : 0, par.pca, par.pcb, par.pcc);;
        q.CalculateAminoAcidBackground();
    }
    // ... or is it an alignment file
    else
    {
        Alignment* pali;
        if (qali==NULL) pali=new(Alignment); else pali=qali;
        if (par.calibrate) {
            printf("\nError in %s: only HHM files can be calibrated.\n",program_name);
            printf("Build an HHM file from your alignment with 'hhmake -i %s' and rerun hhsearch with the hhm file\n\n",infile);
            exit(1);
        }

        if (v>=2 && strcmp(infile,"stdin")) cout<<infile<<" is in A2M, A3M or FASTA format\n";

        // Read alignment from infile into matrix X[k][l] as ASCII (and supply first line as extra argument)
        pali->Read(inf,infile,line);

        // Convert ASCII to int (0-20),throw out all insert states, record their number in I[k][i]
        // and store marked sequences in name[k] and seq[k]
        pali->Compress(infile);

        // Sort out the nseqdis most dissimilar sequences for display in the output alignments
        pali->FilterForDisplay(par.max_seqid,par.coverage,par.qid,par.qsc,par.nseqdis);

//       // Filter alignment for min score per column with core query profile, defined by coverage_core and qsc_core
//       if (par.coresc>-10) pali->HomologyFilter(par.coverage_core, par.qsc_core, par.coresc);

        // Remove sequences with seq. identity larger than seqid percent (remove the shorter of two)
        pali->N_filtered = pali->Filter(par.max_seqid,par.coverage,par.qid,par.qsc,par.Ndiff);

        // Calculate pos-specific weights, AA frequencies and transitions -> f[i][a], tr[i][a]
        pali->FrequenciesAndTransitions(q);
        if (v>=2 && q.Neff_HMM>11.0)
            fprintf(stderr,"WARNING: alignment %s looks too diverse (Neff=%.1f>11). Better check it with an alignment viewer... \n",q.name,q.Neff_HMM);

        // Add transition pseudocounts to query -> p[i][a]
        q.AddTransitionPseudocounts();

        if (!*par.clusterfile) { //compute context-specific pseudocounts?
            // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
            q.PreparePseudocounts();
        } else {
            // Generate an amino acid frequency matrix from f[i][a] with full context specific pseudocount admixture (tau=1) -> g[i][a]
            q.PrepareContextSpecificPseudocounts();
        }

        // Add amino acid pseudocounts to query:  p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
        q.AddAminoAcidPseudocounts(!q.has_pseudocounts ? par.pcm : 0, par.pca, par.pcb, par.pcc);
        q.CalculateAminoAcidBackground();

        if (qali==NULL) delete(pali);
    }
    fclose(inf);

    if (par.forward>=1) q.Log2LinTransitionProbs(1.0);

    return;
}

/////////////////////////////////////////////////////////////////////////////////////
// Do precalculations for q and t to prepare comparison
/////////////////////////////////////////////////////////////////////////////////////
void PrepareTemplate(HMM& q, HMM& t, int format)
{
    if (format==0) // HHM format
    {
        // Add transition pseudocounts to template
        t.AddTransitionPseudocounts();

        // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
        t.PreparePseudocounts();

        // Add amino acid pseudocounts to template: t.p[i][a] = (1-tau)*f[i][a] + tau*g[i][a]
        t.AddAminoAcidPseudocounts(!t.has_pseudocounts ? par.pcm : 0, par.pca, par.pcb, par.pcc);

        t.CalculateAminoAcidBackground();
    }
    else // HHMER format
    {
        // Don't add transition pseudocounts to template
        // t.AddTransitionPseudocounts(par.gapd, par.gape, par.gapf, par.gapg, par.gaph, par.gapi, 0.0);

        // Generate an amino acid frequency matrix from f[i][a] with full pseudocount admixture (tau=1) -> g[i][a]
        t.PreparePseudocounts();

        // DON'T ADD amino acid pseudocounts to temlate: pcm=0!  t.p[i][a] = t.f[i][a]
        t.AddAminoAcidPseudocounts(0, par.pca, par.pcb, par.pcc);
        t.CalculateAminoAcidBackground();
    }

    // Modify transition probabilities to include SS-dependent penalties
    if (par.ssgap) t.UseSecStrucDependentGapPenalties();

    if (par.forward>=1) t.Log2LinTransitionProbs(1.0);

    // Factor Null model into HMM t
    // ATTENTION! t.p[i][a] is divided by pnul[a] (for reasons of efficiency) => do not reuse t.p
    t.IncludeNullModelInHMM(q,t);  // Can go BEFORE the loop if not dependent on template

    return;
}

