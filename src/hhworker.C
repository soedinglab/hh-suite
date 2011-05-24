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
	  
	  // OLD!!!
	  // float alpha = 0;
	  // float beta = 0;
	  // if (par.prefilter) 
	  //   {
	  //     alpha = alpha_NN( q_len, hit_len, q_neff, hit_neff ); 
	  //     beta = beta_NN( q_len, hit_len, q_neff, hit_neff );
	  //   }
	  // hit[bin]->Eval = exp(hit[bin]->logPval+log(hitlist.N_searched)+(alpha*par.hhblits_prefilter_logpval - beta));
	  // hit[bin]->logEval = hit[bin]->logPval+log(hitlist.N_searched)+(alpha*par.hhblits_prefilter_logpval - beta);

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
      if (hit[bin]->score <= SMIN) break;  // break if score for first hit is already worse than SMIN
    }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////
//// Realign q and *(t[bin]) with the MAC algorithm (after the database search)
//////////////////////////////////////////////////////////////////////////////////////
// Remark: The List object is not practical for multithreading, since the "current" pointer
// needs to be locked while the list is read, to prevent other threads from moving it at the 
// same time. We therefore have to recored the current index position here, and after 
// realignment, we have to run through the list again to the previous position
// (see "pos = hitlist.GetPos();" and  "hit_cur = hitlist.Read(pos);". 
// A better solution would be to store all pointers to Hit objects in a hash with 
// the index or some other unique identifier plus the irep value as key: <index>_<irep>
void RealignByWorker(int bin)
{
  Hit hit_cur;
  int nhits=0;
  int pos;
  hit[bin]->irep=1;

  // Prepare MAC comparison(s)
  PrepareTemplate(*q,*(t[bin]),format[bin]);
  t[bin]->Log2LinTransitionProbs(1.0);

#ifdef PTHREAD
          pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif

  // Search positions in hitlist with correct index of template
  hitlist.Reset();
  while (!hitlist.End())
    {
      hit_cur = hitlist.ReadNext();

      if (nhits > (2*par.realign_max) && nhits>=imax(par.B,par.Z)) break;
      if (hit_cur.Eval > par.e)
      	{
      	  if (nhits>=imax(par.B,par.Z)) continue;
      	  if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) continue;
      	  if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
      	}

      if (hit_cur.index==hit[bin]->index) // found position with correct template
        {
	  pos = hitlist.GetPos();

	  // Realign only around previous Viterbi hit
	  hit[bin]->i1 = hit_cur.i1;
	  hit[bin]->i2 = hit_cur.i2;
	  hit[bin]->j1 = hit_cur.j1;
	  hit[bin]->j2 = hit_cur.j2;
	  hit[bin]->nsteps = hit_cur.nsteps;
	  hit[bin]->i = hit_cur.i;
	  hit[bin]->j = hit_cur.j;
	  hit[bin]->realign_around_viterbi=true;

#ifdef PTHREAD
          pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif

	  // Align q to template in *hit[bin]
          hit[bin]->Forward(*q,*(t[bin]));
          hit[bin]->Backward(*q,*(t[bin]));
          hit[bin]->MACAlignment(*q,*(t[bin]));
          hit[bin]->BacktraceMAC(*q,*(t[bin]));
#ifdef PTHREAD
          pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif
          hit_cur = hitlist.Read(pos);

          // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
          hit[bin]->score      = hit_cur.score;
          hit[bin]->score_aass = hit_cur.score_aass;
          hit[bin]->score_ss   = hit_cur.score_ss; // comment out?? => Andrea
          hit[bin]->Pval       = hit_cur.Pval;
          hit[bin]->Pvalt      = hit_cur.Pvalt;
          hit[bin]->logPval    = hit_cur.logPval;
          hit[bin]->logPvalt   = hit_cur.logPvalt;
          hit[bin]->Eval       = hit_cur.Eval;
          hit[bin]->logEval    = hit_cur.logEval;
          hit[bin]->Probab     = hit_cur.Probab;

          // Replace original hit in hitlist with realigned hit
          //hitlist.ReadCurrent().Delete();
          hitlist.Delete().Delete();                // delete list record and hit object
          hitlist.Insert(*hit[bin]);
          hit[bin]->irep++;
        }
      nhits++;
    }
#ifdef PTHREAD
  pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif


  if (hit[bin]->irep==1)
    {
      fprintf(stderr,"*************************************************\n");
      fprintf(stderr,"\nError: could not find template %s in hit list (index:%i dbfile:%s ftell:%i\n\n",hit[bin]->name, hit[bin]->index,hit[bin]->dbfile,(unsigned int)hit[bin]->ftellpos);
      fprintf(stderr,"*************************************************\n");
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
  int rc;                            // return code for threading commands
  int bin;                           // bin index

  // Lock access to bin_status
  rc = pthread_mutex_lock(&bin_status_mutex);

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
