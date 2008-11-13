// Calculate P-values and Probabilities from transitive scoring over whole database
void HitList::TransitiveScoring2()
{
  void PrintMatrix(float** V, int N);
  void PrintMatrix(double** V, int N);

  float** Z;    // matrix of intra-db Z-scores Z_kl
  float** C;    // covariance matrix for Z_k: C_kl = sum_m=1^N (Z_km * Z_lm)
  char** fold;  // fold name of HMM k
  char** fam;   // family of HMM k
  float* Prob;  // probability of HMM k
  float* Zq;    // Zq[k] = Z-score between query and database HMM k
  float* Ztq;   // Ztq[k] = transitive Z-score from query to database HMM k: Ztq[k] = sum_l[ w_ql * Z_lk] / normalization_q
  float* Zrq;   // Zrq[k] = transitive Z-score from database HMM k to query: Zrq[k] = sum_l[ w_kl * Z_lq] / normalization_k
  float* w;     // unnormalized weight matrix; w[l] is w_ql or w_kl, respectively
  int* ll;      // ll[m] is the m'th index l for which Z_lq, Z_lk > Zmin_trans
  int N;        // dimension of weight matrix is NxN
  int M;        // number of HMMs l with Z_ql>Ztrans_min (or Z_lk>Ztrans_min, respectively)
  int k,l,m,n;  // indices for database HMMs 
  char name[NAMELEN];
  Hash<int> index(MAXPROF+7);  // index{name} = index of HMM name in {1,...,N} 
  index.Null(-1);              // Set int value to return when no data can be retrieved
  Hash<int> excluded(13);      // Hash containing names of superfamilies to be excluded from fit
  excluded.Null(0);            // Set int value to return when no data can be retrieved
  Hit hit; 
  
  // Read weights matrix W with index hash and names array
  fprintf(stderr,"Reading in weights file\n");
  FILE* wfile=NULL;
  wfile = fopen(par.wfile,"rb");
  if (v>=1 && wfile==NULL) 
    {
      fprintf(stderr,"Error: %s could not be opened: ",par.wfile,N_searched);
      perror("fopen");
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  fread(&N,sizeof(int),1,wfile);  // read matrix dimension (i.e. number of HMMs in database)
  if (v>=1 && N!=N_searched) 
    {
      fprintf(stderr,"Error: Number %i of HMMs in weight file is different from number %i of HMMs in searched databases. \n",N,N_searched);
      fprintf(stderr,"Skipping caclulation of transitive P-values\n"); 
      par.trans=0;
      return;
    }
  if (v>=2) fprintf(stderr,"Calculating transitive P-values for %i HMMs\n",N);
  // Read names of HMMs (to specify mapping of HMM to matrix indices)
  for (k=0; k<N; k++) 
    {
      fread(name,sizeof(char),IDLEN,wfile);
      index.Add(name,k);
    }
  // Read symmetric Z-scores matrix
  Z = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      Z[k] = new(float[N]);
      for (l=0; l<k; l++) Z[k][l] = Z[l][k];
      fread(Z[k]+k,sizeof(float),N-k,wfile);   
    }
  // Read symmetric covariance matrix
  C = new(float*[N]);
  for (k=0; k<N; k++) 
    {
      C[k] = new(float[N]);
      for (l=0; l<k; l++) C[k][l] = C[l][k];
      fread(C[k]+k,sizeof(float),N-k,wfile);
    }
  fclose(wfile);

  // Allocate memory
  Zq = new(float[N]);
  Ztq = new(float[N]);
  Zrq = new(float[N]);
  fold = new(char*[N]);
  fam = new(char*[N]);
  Prob = new(float[N]);
  ll = new(int[N]);
  w = new(float[N]);

  // Transform P-values to normally distributed Z-scores and store in Zq vector
  fprintf(stderr,"Transform P-values to Z-scores\n");
  float Zmax_neg   = Score2Z( -log(MINEVALEXCL) + log(N_searched) ); // calculate Z-score corresponding to E-value MINEVALEXCL
  float Zmin_trans = Score2Z( -log(par.Emax_trans) + log(N_searched) ); // calculate Z-score corresponding to E-value par.Emax_trans
  printf("Zmax = %6.2f   Zmin = %6.2f \n",Zmax_neg,Zmin_trans);

  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      if (hit.irep>1) continue;
      k = index.Show(hit.name);
      if (k<0) {fprintf(stderr,"Error: no index found in weights file for domain %s\n",hit.name); exit(1);}      
      if (hit.logPvalt<0)
	Zq[k] = 0.5*Score2Z(fabs(hit.logPval)) + 0.5*Score2Z(fabs(hit.logPvalt));  // Zq[k] = 0.5*(Zkq + Zqk)
      else 
	Zq[k] = Score2Z(fabs(hit.logPval));                           // Zq[k] = Zqk 
//      printf("%4i  %-10.10s logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPvalt,Zq[k]);
//      if (isnan(Zq[k])) {
// 	fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	printf("%4i  %-10.10s logPval=%9g  logPvalt=%9g  Zq=%9f\n",k,hit.name,hit.logPval,hit.logPvalt,Zq[k]);
//   	par.trans=0;
// 	return;
//       }
      if (Zq[k]>Zmax_neg) excluded.Add(hit.fold);
      fold[k] = new(char[IDLEN]);
      fam[k] = new(char[IDLEN]);
      strcpy(fold[k],hit.fold);
      strcpy(fam[k],hit.fam);
      weight[k] = hit.weight;
      Prob[k] = hit.Probab;
   }
  
  if (v>=3) 
    {
      excluded.Reset();
      while (!excluded.End())
	{
	  excluded.ReadNext(name);
	  printf("Excluded fold %s from fitting to Ztq\n",name);
	}
    }


  ////////////////////////////////////////////////////////////////
  // Calculate transitive score (query->l) Zt[l]
  
  // Construct vector ll of indices l for which Z_lq > Zmin_trans
  m = 0;
  for (l=0; l<N; l++)
    if (Zq[l]>=Zmin_trans) ll[m++]=l;
  M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_trans
  
//   for (m=0; m<M; m++)
//     fprintf(stderr,"m=%-4i l=%-4i  %-10.10s  Zq[l]=%7f\n",m,ll[m],fam[ll[m]],Zq[ll[m]]);

  if (M<=1) 
    for (k=0; k<N; k++) Ztq[k]=0.0;
  else
    {
      // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
      double** Csub = new(double*[M]);
      double** Cinv = new(double*[M]);
      for (m=0; m<M; m++) 
	{
	  Csub[m] = new(double[M]);
	  Cinv[m] = new(double[M]);
	  for (n=0; n<M; n++)
	    Csub[m][n] = double(C[ll[m]][ll[n]]);
	}
      
      if (v>=3) 
	{
	  fprintf(stderr,"Covariance matrix\n");
	  PrintMatrix(Csub,M);
	}
      
//       // Invert Csub
//       fprintf(stderr,"Calculate inverse of covariance submatrix\n");
//       InvertMatrix(Cinv,Csub,M);
      
//       if (v>=3) 
// 	{
// 	  fprintf(stderr,"Inverse covariance matrix\n");
// 	  PrintMatrix(Cinv,M);
// 	}
      
      // Calculate weights w[l] 
      for (m=0; m<M; m++) 
	{
	  double sum = 0.0;
	  for (n=0; n<M; n++)
	    sum += 1.0 * Csub[m][n]; 
	  printf("w[%4i] = %-8.5f\n",ll[m],1.0/sum);
	  w[m] = (sum>0? Zq[ll[m]] / sum : 0.0);
	}
      for (l=0; l<M; l++) delete[](Cinv[l]);
      delete[](Cinv);
      
      // Calculate Ztq[k] for all HMMs k
      fprintf(stderr,"Calculate Ztq vector of transitive Z-scores\n");
      float norm = NormalizationFactor(Csub,w,M);
      for (k=0; k<N; k++) 
	{
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * Z[ll[m]][k];
	  Ztq[k] = sumZ/norm;   
	}
      
      for (l=0; l<M; l++) delete[](Csub[l]);
      delete[](Csub);
    }

  ////////////////////////////////////////////////////////////////
  // Calculate reverse transitive score (l->query-) Zrq[l]

  fprintf(stderr,"Calculate Zrq vector of transitive Z-scores\n");
  for (k=0; k<N; k++) 
    {
      // Construct vector ll of indices l for which Z_lk > Zmin_tran
      m = 0;
      for (l=0; l<N; l++)
	  if (Z[l][k]+Z[k][l]>=2*Zmin_trans) ll[m++]=l;  
      int M = m;  // number of indices l for which Z_lq,Z_lk > Zmin_tran


//    fprintf(stderr,"\nfam[k]: %s\n",fam[k]);
//    for (m=0; m<M; m++)
//    printf(stderr,"m=%-4i k=%-4i  l=%-4i  %-10.10s  Zq[l]=%7f  Z_lk=%7f  \n",m,k,ll[m],fold[ll[m]],Zq[ll[m]],Z[k][ll[m]]);
           
      if (M<=1) 
	{
	  Zrq[k] = Zq[k];
	}
      else 
	{
	  // Generate submatrix of C for indices l for which Z_lq,Z_lk > Zmin_trans 
	  double** Csub = new(double*[M]);
	  for (m=0; m<M; m++) 
	    {
	      Csub[m] = new(double[M]);
	      for (n=0; n<M; n++)
		Csub[m][n] = double(C[ll[m]][ll[n]]);
	    }
//        fprintf(stderr,"Covariance matrix\n");
//        PrintMatrix(Csub,M);
	      
	  if (M<=2) 
	    {
	      for (m=0; m<M; m++) w[m] = 1.0/M;
	    }
	  else 
	    {
	      
	      double** Cinv = new(double*[M]);
	      for (m=0; m<M; m++) Cinv[m] = new(double[M]);

// 	      // Invert Csub
// 	      InvertMatrix(Cinv,Csub,M);
	      
// //	      fprintf(stderr,"Inverse covariance matrix\n");
// //	      PrintMatrix(Cinv,M);
	      
	      // Calculate weights w[l] 
	      for (m=0; m<M; m++) 
		{
		  double sum = 0.0;
		  for (n=0; n<M; n++)
		    sum += 1.0 * Csub[m][n]; 
		  w[m] = (sum>0? Z[ll[m]][k] / sum : 0.0);
		}

//            for (m=0; m<M; m++) fprintf(stderr,"w[%i]=%8.2g\n",m,w[m]);

	      for (l=0; l<M; l++) delete[](Cinv[l]);
	      delete[](Cinv);
	    }

	  // Calculate Zrq[k] and normalize
	  float norm = NormalizationFactor(Csub,w,M);
	  double sumZ = 0.0;
	  for (m=0; m<M; m++) 
	    sumZ += w[m] * Zq[ll[m]];
	  Zrq[k] = sumZ/norm;   
	  
	  for (l=0; l<M; l++) delete[](Csub[l]);
	  delete[](Csub);
	} 

//    fprintf(stderr,"\nZq[k]=%8.2g  Zq1[k]=%8.2g\n",Zq[k],Zrq[k]);
    }

  // Total Z-score = weighted sum over original Z-score, forward transitive and reverse transitive Z-score
  for (k=0; k<N; k++) 
    {
      float Zqtot =  Zq[k] + par.wtrans*(Ztq[k]+Zrq[k]);
//        if (isnan(Zqtot))
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
// 	  par.trans=0;
// 	  return;
// 	}
      if (v>=2 &&  Zq[k] + Zqtot > 2*Zmin_trans) {
	printf("%4i  %-10.10s Zq=%6.2f  Ztq=%6.2f  Zrq=%6.2f  -> Zqtot=%6.2f\n",k,fam[k],Zq[k],Ztq[k],Zrq[k],Zqtot);
      }
      Ztq[k] = Zqtot;
    }

  // Calculate mean and standard deviation of Z1q
  fprintf(stderr,"Calculate mean and standard deviation of Ztq\n");
  double sumw=0.0;
  double sumZ=0.0;
  double sumZ2=0.0;
  for (k=0; k<N; k++) 
    {  
      if (excluded.Contains(fold[k])) continue;
      sumw  += weight[k];
      sumZ  += weight[k]*Ztq[k];
      sumZ2 += weight[k]*Ztq[k]*Ztq[k];
//       if (isnan(sumZ)) 
// 	{
// 	  fprintf(stderr,"Error: a floating point exception occurred. Skipping transitive scoring\n"); 
// 	  printf("%4i  %-10.10s Zq=%9f  Zrq=%9f  Ztq=%9f\n",k,fam[k],Zq[k],Zrq[k],Ztq[k]);
// 	  par.trans=0;
// 	  return;
// 	}
    }
  float mu = sumZ/sumw;  
  float sigma = sqrt(sumZ2/sumw-mu*mu);
  if (v>=2) printf("mu(Ztq)=%6.3f  sigma(Ztq)=%6.2f\n",mu,sigma);
  sigma *= 1.01;// correct different fitting of EVD and normal variables

  // Normalize Ztq and calculate P1-values
  fprintf(stderr,"Normalize Ztq and calculate P1-values\n");
  Reset();
  while (!End()) 
    {
      hit = ReadNext();
      hit.logPval = -Z2Score((Ztq[index.Show(hit.name)]-mu)/sigma);
      hit.E1val = N_searched*(hit.logPval<-100? 0.0 : exp(hit.logPval));
      hit.Probab = Probab(hit);
      hit.score_sort = hit.logPval;
      Overwrite(hit);                                        // copy hit object into current position of hitlist
   }

  for (k=0; k<N; k++) delete[](Z[k]);
  for (k=0; k<N; k++) delete[](C[k]);
  for (k=0; k<N; k++) delete[](fold[k]);
  for (k=0; k<N; k++) delete[](fam[k]);
  delete[](C);
  delete[](Z);
  delete[](fold);
  delete[](fam);
  delete[](Prob);
  delete[](ll);
  delete[](Zq);
  delete[](Ztq);
}

