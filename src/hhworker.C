// hhworker.C

///////////////////////////////////////////////////////////////////////////////////////
//// Do the pairwise comparison of q and *(t[bin]) for the database search
//////////////////////////////////////////////////////////////////////////////////////
void AlignByWorker(int bin)
{
  // Prepare q ant t and compare
  PrepareTemplate(*q,*(t[bin]),format[bin]);

  // Do HMM-HMM comparison, store results if score>SMIN, and try next best alignment
  for (hit[bin]->irep=1; hit[bin]->irep<=par.altali; hit[bin]->irep++)
    {
      if (par.forward==0)
        {
	  hit[bin]->Viterbi(*q,*(t[bin]));
          if (hit[bin]->irep>1 && hit[bin]->score <= SMIN) break;
          hit[bin]->Backtrace(*q,*(t[bin]));
        }
      else if (par.forward==1)
        {
          hit[bin]->Forward(*q,*(t[bin]));
          hit[bin]->StochasticBacktrace(*q,*(t[bin]),1); // the 1 selects maximization instead of stochastic backtracing
        }
      else if (par.forward==2)
        {
          hit[bin]->Forward(*q,*(t[bin]));
          hit[bin]->Backward(*q,*(t[bin]));
          hit[bin]->MACAlignment(*q,*(t[bin]));
          hit[bin]->BacktraceMAC(*q,*(t[bin]));
        }
      hit[bin]->score_sort = hit[bin]->score_aass;
      if (hit[bin]->score <= SMIN) hit[bin]->lastrep=1; else hit[bin]->lastrep=0;
      //printf ("%-12.12s  %-12.12s   irep=%-2i  score=%6.2f\n",hit[bin]->name,hit[bin]->fam,hit[bin]->irep,hit[bin]->score);

#ifdef PTHREAD
      pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif
      hitlist.Push(*(hit[bin]));            // insert hit at beginning of list (last repeats first!)

      if (par.early_stopping_filter)
	{
	  // Calculate Evalue
	  float q_len = log(q->L)/LOG1000;
	  float hit_len = log(hit[bin]->L)/LOG1000;
	  float q_neff = q->Neff_HMM/10.0;
	  float hit_neff = hit[bin]->Neff_HMM/10.0;
	  float lamda = lamda_NN( q_len, hit_len, q_neff, hit_neff ); 
	  float mu    =    mu_NN( q_len, hit_len, q_neff, hit_neff ); 
	  hit[bin]->logPval = logPvalue(hit[bin]->score,lamda,mu);
	  
	  float alpha = 0;
	  float log_Pcut = log(par.prefilter_evalue_thresh / par.dbsize);
	  float log_dbsize = log(par.dbsize);

	  if (par.prefilter) 
	    alpha = par.alphaa + par.alphab * (hit_neff - 1) * (1 - par.alphac * (q_neff - 1));
      
	  hit[bin]->Eval = exp(hit[bin]->logPval + log_dbsize + (alpha * log_Pcut)); 
	  hit[bin]->logEval = hit[bin]->logPval + log_dbsize + (alpha * log_Pcut); 

	  par.filter_sum -= par.filter_evals[par.filter_counter];
	  par.filter_evals[par.filter_counter] = 1.0/(1.0+hit[bin]->Eval);
	  par.filter_sum += par.filter_evals[par.filter_counter];

	  //printf("E-val: %4.2g   1/(1+Eval): %4.2g  => new sum: %16.2f\n",hit[bin]->Eval,par.filter_evals[par.filter_counter],par.filter_sum);
	  par.filter_counter++;
	  if (par.filter_counter==par.filter_length) {par.filter_counter=0;}
	}
#ifdef PTHREAD
      pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif

      if (par.forward>0) break; // find only best alignment for forward algorithm and stochastic sampling
      if (hit[bin]->score <= SMIN) break;  // break if score for previous hit is already worse than SMIN
    }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////
//// Realign q and with *(t[bin]) in all hits from same tempate using  MAC algorithm 
//////////////////////////////////////////////////////////////////////////////////////
void RealignByWorker(int bin)
{
  // Realign all hits pointed to by list List<void*>* hit[bin]->plist_phits;
  // This list is set up in HHseach and HHblits at the beginning of perform_realign()
  Hit* hit_cur;

  // Prepare MAC comparison(s)
  PrepareTemplate(*q,*(t[bin]),format[bin]);
  t[bin]->Log2LinTransitionProbs(1.0);

  hit[bin]->irep=1; 
  hit[bin]->plist_phits->Reset();
  while (! hit[bin]->plist_phits->End())
    {
      // Set pointer hit_cur to next hit to be realigned
      hit_cur = (Hit*) hit[bin]->plist_phits->ReadNext();

      // Realign only around previous Viterbi hit
      // hit[bin] = *hit_cur; is not possible because the pointers to the DP matrices would be overwritten
      hit[bin]->i1 = hit_cur->i1;
      hit[bin]->i2 = hit_cur->i2;
      hit[bin]->j1 = hit_cur->j1;
      hit[bin]->j2 = hit_cur->j2;
      hit[bin]->nsteps = hit_cur->nsteps;
      hit[bin]->i  = hit_cur->i;
      hit[bin]->j  = hit_cur->j;
      hit[bin]->realign_around_viterbi=true;
      
      // Align q to template in *hit[bin]
      hit[bin]->Forward(*q,*(t[bin]));
      hit[bin]->Backward(*q,*(t[bin]));
      hit[bin]->MACAlignment(*q,*(t[bin]));
      hit[bin]->BacktraceMAC(*q,*(t[bin]));
      
      // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
      hit[bin]->score      = hit_cur->score;
      hit[bin]->score_ss   = hit_cur->score_ss;
      hit[bin]->score_aass = hit_cur->score_aass;
      hit[bin]->score_sort = hit_cur->score_sort;
      hit[bin]->Pval       = hit_cur->Pval;
      hit[bin]->Pvalt      = hit_cur->Pvalt;
      hit[bin]->logPval    = hit_cur->logPval;
      hit[bin]->logPvalt   = hit_cur->logPvalt;
      hit[bin]->Eval       = hit_cur->Eval;
      hit[bin]->logEval    = hit_cur->logEval;
      hit[bin]->Probab     = hit_cur->Probab;
      
      // Replace original hit in hitlist with realigned hit
      hit_cur->Delete();     // delete content of pointers etc. of hit_cur (but not DP matrices)
      *hit_cur = *hit[bin];  // copy all variables and pointers from *hit[bin] into hitlist

      hit[bin]->irep++;
    }

  return;
}



#ifdef PTHREAD

///////////////////////////////////////////////////////////////////////////////////////
//// This is the main thread loop that waits for new jobs (i.e. pairwise alignment) and executes them
//////////////////////////////////////////////////////////////////////////////////////
void* WorkerLoop(void* data)
{
  int thread_id = (*((Thread_args*)data)).thread_id;  // typecast 'data' from pointer-to-void to pointer-to-Thread_args
  void (*ExecuteJob)(int) = (*((Thread_args*)data)).function; // dito; data.function is a pointer to a function(int) that returns void
  int bin;                           // bin index

  // Lock access to bin_status
  pthread_mutex_lock(&bin_status_mutex);

  ///////////////////////////////////////////////////////////////////////////////////////
  // Take jobs from one of the SUBMITTED bins and execute them. If no submitted jobs found, wait for signal 'new_job'
  while (reading_dbs || jobs_submitted>0)
    {

     if (jobs_submitted>0 ) // submitted job in one of the bins?
        {
          bin = PickBin(SUBMITTED);
          if (DEBUG_THREADS)
            fprintf(stderr,"Thread %3i:   start job in bin %i       jobs running: %i  jobs_submitted:%i \n",thread_id,bin,jobs_running,jobs_submitted);
          jobs_running++;            // the order of the following three lines is important, since...
          jobs_submitted--;          // ... the main thread tries to find a free bin...
          bin_status[bin] = RUNNING; // ... if jobs_running+jobs_submitted<bins !

          // Execute job
          rc = pthread_mutex_unlock(&bin_status_mutex); // unlock access to bin_status
          if (DEBUG_THREADS)
            fprintf(stderr,"Thread %3i:   executing HMM %-10.10s jobs running: %i  jobs_submitted:%i \n",thread_id,t[bin]->name,jobs_running,jobs_submitted);

	  // Here the function AlignByWorker or RealignByWorker is called
          ExecuteJob(bin);

          rc = pthread_mutex_lock(&bin_status_mutex); // lock access to bin_status

          bin_status[bin] = FREE;
          jobs_running--;

          // Signal completion of job
          rc = pthread_cond_signal(&finished_job);

          if (DEBUG_THREADS)
            fprintf(stderr,"Thread %3i:   finished job in bin %i    jobs running: %i  jobs_submitted:%i \n",thread_id,bin,jobs_running,jobs_submitted);
        }
     else
       {

         if (DEBUG_THREADS)
           fprintf(stderr,"Thread %3i:   waiting for new job ...   jobs running: %i  jobs_submitted:%i \n",thread_id,jobs_running,jobs_submitted);
         rc = pthread_cond_wait(&new_job, &bin_status_mutex);
       }

      // Unlock access to bin_status
      rc = pthread_mutex_unlock(&bin_status_mutex);
      // Lock access to bin_status
      rc = pthread_mutex_lock(&bin_status_mutex);
    }
  ///////////////////////////////////////////////////////////////////////////////////////

  // Unlock access to bin_status
  rc = pthread_mutex_unlock(&bin_status_mutex);

  // Exit thread automatically
  pthread_exit(NULL);
  return NULL;
}

#endif
