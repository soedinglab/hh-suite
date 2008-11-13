void RealignByWorker(int bin)
{
  Hit hit_cur;
  int nhits=0;
  int pos;
  hit[bin]->irep=1;

  // Prepare MAC comparison(s)
  PrepareTemplate(q,*(t[bin]),format[bin]);
  t[bin]->Log2LinTransitionProbs(1.0);
  
  // Align q to template in *hit[bin]
  hit[bin]->Forward(q,*(t[bin])); 
  hit[bin]->Backward(q,*(t[bin])); 
  hit[bin]->MACAlignment(q,*(t[bin]));
  hit[bin]->BacktraceMAC(q,*(t[bin]));

#ifdef PTHREAD
	  pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif

  // Search positions in hitlist with correct index of template
  hitlist.Reset();
  while (!hitlist.End()) 
    {
      hit_cur = hitlist.ReadNext();
      if (nhits>=imax(par.B,par.Z) && nhits.Evalue>par.e) break;
      //    fprintf(stderr,"t->name=%s  hit_cur.name=%s  hit[bin]->irep=%i  nhits=%i  hit_cur.index=%i  hit[bin]->index=%i\n",t[bin]->name,hit_cur.name,hit_cur.irep,nhits,hit_cur.index,hit[bin]->index); 
      if (nhits>=imax(par.b,par.z) && hit_cur.Probab < par.p) break;
      if (nhits>=imax(par.b,par.z) && hit_cur.Eval > par.E) continue;
      if (hit_cur.index==hit[bin]->index) // found position with correct template
	{
	  //  	  fprintf(stderr,"  t->name=%s   hit_cur.irep=%i  hit[bin]->irep=%i  nhits=%i\n",t[bin]->name,hit_cur.irep,hit[bin]->irep,nhits); 
	  if (hit[bin]->irep > 1) 
	    {
	      pos = hitlist.GetPos();
#ifdef PTHREAD
	      pthread_mutex_unlock(&hitlist_mutex); // unlock access to hitlist
#endif
	      // Align q to template in *hit[bin]
	      hit[bin]->Forward(q,*(t[bin])); 
	      hit[bin]->Backward(q,*(t[bin])); 
	      hit[bin]->MACAlignment(q,*(t[bin]));
	      hit[bin]->BacktraceMAC(q,*(t[bin]));
#ifdef PTHREAD
	      pthread_mutex_lock(&hitlist_mutex);   // lock access to hitlist
#endif	  
	      hit_cur = hitlist.Read(pos);
	    }

	  // Overwrite *hit[bin] with Viterbi scores, Probabilities etc. of hit_cur
	  hit[bin]->score      = hit_cur.score;
	  hit[bin]->score_aass = hit_cur.score_aass;
	  hit[bin]->score_ss   = hit_cur.score_ss;
	  hit[bin]->Pval       = hit_cur.Pval;
	  hit[bin]->Pvalt      = hit_cur.Pvalt;
	  hit[bin]->logPval    = hit_cur.logPval;
	  hit[bin]->logPvalt   = hit_cur.logPvalt;
	  hit[bin]->logP1val   = hit_cur.logP1val;
	  hit[bin]->Eval       = hit_cur.Eval;
	  hit[bin]->E1val      = hit_cur.E1val;
	  hit[bin]->Probab     = hit_cur.Probab;
	  
	  // Replace original hit in hitlist with realigned hit
	  hitlist.Delete();
	  hitlist.Insert(*hit[bin]);	 	  
	  hit[bin]->irep++;
	}
      nhits++;
    }
#if
