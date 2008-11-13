#include <stdio.h>       /* standard I/O routines                     */
#include <pthread.h>     /* pthread functions and data structures     */
#include <stdlib.h>      /* rand() and srand() functions              */

const int MAXTHREADS=256;
const int MAXBINS=300;
enum bin_states {FREE=0, SUBMITTED=1, RUNNING=2};

int threads=3;
int bins=6; 

char bin_status[MAXBINS];
float* t[MAXBINS];
int jobs_running;
int jobs_submitted;
int read_from_db=1;
int v=3;

// This condition variable signals that a new job has been submitted
pthread_cond_t new_job = PTHREAD_COND_INITIALIZER;

// Mutex assures exclusive access by threads to bin_status[] and new_job
pthread_mutex_t bin_status_mutex = PTHREAD_MUTEX_INITIALIZER;



///////////////////////////////////////////////////////////////////////////////////////
//// Return a bin with the desired status, return -1 if no such bin found
//////////////////////////////////////////////////////////////////////////////////////
inline int pick_bin(char status)
{
  for (int b=0; b<bins; b++) {if (bin_status[b]==status) return b;}
  return -1;
}


///////////////////////////////////////////////////////////////////////////////////////
//// This is the main thread loop that waits for new jobs and executes them
//////////////////////////////////////////////////////////////////////////////////////
void* execute_jobs_loop(void* data)
{
  struct timespec exec_time;         // simulate job execution
  exec_time.tv_sec = 0;
 
  int thread_id = *((int*)data);  
  int rc;                            // return code for threading commands
  int bin;                           // bin index
 
  // Lock access to bin_status
  rc = pthread_mutex_lock(&bin_status_mutex);

  // If no submitted jobs are in the queue we have to wait for a new job ...
  if (jobs_submitted==0) 
    {
      printf("Thread %3i:   waiting for new job ...   jobs running: %i  jobs_submitted:%i \n",thread_id,jobs_running,jobs_submitted);
      rc = pthread_cond_wait(&new_job, &bin_status_mutex);
    }

  ///////////////////////////////////////////////////////////////////////////////////////
  // Take jobs from one of the SUBMITTED bins and execute them. If no submitted jobs found, wait for signal 'new_job'
  while (read_from_db>0 || jobs_submitted>0)
    {

      bin = pick_bin(SUBMITTED);
      if (bin>=0) // found submitted job in bin?
	{
	  if (v>=3) 
	    printf("Thread %3i:   start job in bin %i       jobs running: %i  jobs_submitted:%i \n",thread_id,bin,jobs_running,jobs_submitted);

	  bin_status[bin] = RUNNING;
	  jobs_submitted--;
	  jobs_running++;

	  // Unlock access to bin_status
	  rc = pthread_mutex_unlock(&bin_status_mutex);

	  // Do job 
	  long int time = (long int)(*(t[bin])*10000000.0);
	  if (v>=3) 
	    printf("Thread %3i:   execution time=%5.2f     jobs running: %i  jobs_submitted:%i \n",thread_id,*(t[bin]),jobs_running,jobs_submitted);
	  for (int i=0; i<time; i++) ;

	  // Lock access to bin_status
	  rc = pthread_mutex_lock(&bin_status_mutex);

	  bin_status[bin] = FREE;
	  jobs_running--;

	  if (v>=3) 
	    printf("Thread %3i:   finished job in bin %i    jobs running: %i  jobs_submitted:%i \n",thread_id,bin,jobs_running,jobs_submitted);
	}
 
      // If no submitted jobs are in the queue we have to wait for a new job ...
      if (jobs_submitted==0 && read_from_db) 
	{
	  printf("Thread %3i:   waiting for new job ...   jobs running: %i  jobs_submitted:%i \n",thread_id,jobs_running,jobs_submitted);
	  rc = pthread_cond_wait(&new_job, &bin_status_mutex);
	}
    }
  ///////////////////////////////////////////////////////////////////////////////////////
  
  // Unlock access to bin_status
  rc = pthread_mutex_unlock(&bin_status_mutex);
  
  // Exit thread automatically
  return NULL;
}
  


///////////////////////////////////////////////////////////////////////////////////////
////  Main
///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

  int bin;                           // bin index
  int thread_id[MAXTHREADS];         // so each thread knows its number (for debugging)
  pthread_t pthread[MAXTHREADS];     // info on thread's structures (needed by system)
  int rc;                            // return code for threading commands
  
  struct timespec wait_finished;     // Time delay between checking if one of the threads has finished its job
  wait_finished.tv_sec = 0;
  wait_finished.tv_nsec = 1000;      // 1us = 1000 ns delay
 
  // Initialize bins, threads
  jobs_running = jobs_submitted = 0;
  for (bin=0; bin<bins; bin++) 
    {
      t[bin]=new(float);
      bin_status[bin] = FREE;
    }
  for (int j=0; j<threads; j++) 
    {
      thread_id[j] = j;
      pthread_create(&pthread[j], NULL, execute_jobs_loop, (void*)&thread_id[j]);
    }

  ///////////////////////////////////////////////////////////////////////////////////////
  // The loop submits jobs into free bins as soon as they become available
  while(read_from_db!=0) 
    {
      
      // Submit jobs until no bin is free anymore
      while (jobs_submitted+jobs_running<bins) 
	{
	  // Submit next job (i.e. read next HMM from database, allocate free bin, and send signal to threads)
	  bin = pick_bin(FREE);
	  if (bin<0) {
	    fprintf(stderr,"Error: found no free bin! jobs running: %i  jobs_submitted:%i  threads:%i\n",jobs_running,jobs_submitted,threads); 
	    for (bin=0; bin<bins; bin++) fprintf(stderr,"bin_status[%i]=%i\n",bin,bin_status[bin]);
	    exit(1);
	  }

	  *(t[bin]) = 100.*rand()/RAND_MAX;
	  if (*(t[bin])<2.0) read_from_db=0;
	  long int time = (long int)(100.*rand()/RAND_MAX*2000000.0);
	  for (int i=0; i<time; i++) ;
	  
	  if (read_from_db==0) break;

	  // Lock access to bin_status
	  rc = pthread_mutex_lock(&bin_status_mutex);
	  
	  // Send the job in bin to a thread
	  bin_status[bin] = SUBMITTED;
	  jobs_submitted++;
	  
	  if (v>=3) 
	    printf("Main: put job into bin %i  data=%6.3f  Jobs running: %i  jobs_submitted:%i \n",bin,*(t[bin]),jobs_running,jobs_submitted);

	  // Unlock access to bin_status and send signal
	  rc = pthread_mutex_unlock(&bin_status_mutex);
	  
	  // Restart threads jobs waiting for a signal
	  while (jobs_submitted>0 && jobs_running<threads) rc = pthread_cond_signal(&new_job);
	}
      
      // Wait until job finishes and a bin becomes free
      while (jobs_submitted+jobs_running==bins) 
	{	
	  // In case queue had been empty and some jobs are now waiting for a signal
	  if (jobs_submitted>0) rc = pthread_cond_signal(&new_job); 
	  nanosleep(&wait_finished,NULL);
	}
    }
  ///////////////////////////////////////////////////////////////////////////////////////

  // No more HMMs in database => wait until all threads have finished
  if (v>=3) 
    printf("No more jobs read from database         Jobs running:%i  jobs_submitted:%i \n",jobs_running,jobs_submitted);
  while (jobs_running>0 || jobs_submitted>0) 
    {
      if (jobs_submitted>0) rc = pthread_cond_signal(&new_job); 
      nanosleep(&wait_finished,NULL);
    }

  printf("Glory,  we are done.\n");
  exit(0);
}

 

