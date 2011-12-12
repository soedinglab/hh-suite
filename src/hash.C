// hash.C
// Class for Hash data structure 
// * works in the same way as a hash in Perl
// * keys are strings of type char*
// * data elements are of type Typ 
// * objects have to be declared with maximal size, e.g. Hash hash1(10000) (num_slots should not be a power of 2)
// * works also if maximal size is exceeded, but gets slower by a factor ~num_keys/num_slots 

#ifndef HASH
#define HASH

#ifndef MAIN
#include <iostream>   // cin, cout, cerr
#include <cstdio>     // printf
#include <stdlib.h>   // exit
#include <string.h>   // strcmp, strstr
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <ctype.h>    // islower, isdigit etc
#include <time.h>     // clock
#include <errno.h>    // perror()
using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::ifstream;
using std::ofstream;
#endif

#ifndef JLIST
#define JLIST
#include "list.h"       // List<Typ>
////  #include "list.C"       ////////////////////////////////// DEBUG
#endif

#include "hash.h"



////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Methods of class Hash ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////
//                                Private Methods


////////////////////////////////////////////////////////////////////////////////////////////
//                        Constructor and Destructor of Hash

////////////////////////////////////////////////////////////////////////////////////////////
// Constructor of class Hash
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Hash<Typ>::Hash()
{
  num_keys=0; max_len=0; prev=curr=num_slots = 0; slot=NULL;
}

template<class Typ> 
Hash<Typ>::Hash(int nslots)
{
  num_keys=0; max_len=0; prev=curr=num_slots = nslots;
  slot = new Slot<Typ>*[num_slots];             //Create array of num_slots slots
  for (int i=0; i<num_slots; i++) slot[i]=NULL; //set pointers to NULL
  fail = static_cast<Typ>(0);
}

template<class Typ> 
Hash<Typ>::Hash(int nslots, Typ f)
{
  num_keys=0; max_len=0; prev=curr=num_slots = nslots;
  slot = new Slot<Typ>*[num_slots];             //Create array of num_slots slots
  for (int i=0; i<num_slots; i++) slot[i]=NULL; //set pointers to NULL
  fail=f;
}


////////////////////////////////////////////////////////////////////////////////////////////
// Destructor of class Hash
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Hash<Typ>::~Hash()
{
  RemoveAll();
  delete[] slot;
}


////////////////////////////////////////////////////////////////////////////////////////////
// Hash function
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
inline unsigned int Hash<Typ>::HashValue(char* key)  //returns the hash value for key 
    {
      // Calculate a hash value by the division method: 
      // Transform key into a natural number k = sum ( key[i]*128^(L-i) ) and calculate i= k % num_slots. 
      // Since calculating k would lead to an overflow, i is calculated iteratively 
      // and at each iteration the part divisible by num_slots is subtracted, i.e. (% num_slots is taken).
      if (key==NULL) {printf("Warning from hash.C: key=NULL\n"); return 0;}
      unsigned int i=0;     // Start of iteration: k is zero
      char* c = key;
      while(*c) i = ((i<<7) + *(c++)) % num_slots; 
      key_len = c - key;
      //cerr<<"      Hash value for \'"<<key<<"\' is "<<i<<"\n";
      return i;
    }

////////////////////////////////////////////////////////////////////////////////////////////
// Create new hash (and delete any data present)
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
void Hash<Typ>::New(int nslots, Typ f)
{
  fail=f; 
  RemoveAll(); 
  delete[] slot;
  num_keys=0; max_len=0; prev=curr=num_slots = nslots;
  slot = new Slot<Typ>*[num_slots];             //Create array of num_slots slots
  for (int i=0; i<num_slots; i++) slot[i]=NULL; //set pointers to NULL
}



////////////////////////////////////////////////////////////////////////////////////////////
//                  Methods that work with a key supplied as an argument 

////////////////////////////////////////////////////////////////////////////////////////////
// Return data element for key. Returns 'fail' if key does not exist 
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Typ Hash<Typ>::Show(char* key)
{
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot) return fail;
  pslot->Reset();
  do{
    if(!strcmp(pslot->ReadNext().key,key)) return pslot->ReadCurrent().data;
  } while(!pslot->End());
  return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Add/replace key/data pair to hash and return address of data element
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Typ* Hash<Typ>::Add(char* key, Typ data)
{
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot) { num_keys++; KeyLen(); slot[i]=new(Slot<Typ>); return slot[i]->Push(key_len,key,data);}
  pslot->Reset();
  do
    {
      pairp = pslot->ReadNextAddress();
      if(!strcmp(pairp->key,key))
        {
          pairp->data=data;
          pslot->Overwrite(*pairp);
          return &(pairp->data);
        }
    } while(!pslot->End());
  num_keys++;
  KeyLen();
  return pslot->Push(key_len,key,data);
}



////////////////////////////////////////////////////////////////////////////////////////////
// Add key to hash and return address of data element. 
// If key exists leave data element unchanged, else set it to 'fail'.
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Typ* Hash<Typ>::Add(char* key)
{
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot) { num_keys++; KeyLen(); slot[i]=new(Slot<Typ>); return slot[i]->Push(key_len,key,fail);}
  pslot->Reset();
  do
    {
      if(!strcmp(pslot->ReadNext().key,key))
        {
          return &((pslot->ReadCurrentAddress())->data);
        }
    } while(!pslot->End());
  num_keys++;
  KeyLen();
  return pslot->Push(key_len,key,fail);
}


/////////////////////////////////////////////////////////////////////////////////////////////
// Remove key from hash and return data element for key ('fail' if key does not exist)
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Typ Hash<Typ>::Remove(char* key)
{
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot) return fail;
  pslot->Reset();
  do
    {
      if(!strcmp(pslot->ReadNext().key,key)) 
	{
	  Pair<Typ> pair = pslot->ReadCurrent();
	  num_keys--; 
	  pslot->Delete();
	  // Delete key-Array
	  delete[] pair.key;
	  // if key was the only element in pslot then delete whole list
	  if (pslot->Size()==0) {delete pslot; slot[i]=0;} 
	  //	  return pslot->ReadCurrent().data;
	  return pair.data;
	} 
    } while(!pslot->End()); 
  return fail;
}


////////////////////////////////////////////////////////////////////////////////////////////
// Remove all keys from hash;
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
void Hash<Typ>::RemoveAll()
{
  for(int i=0; i<num_slots; i++) 
    if(slot[i]) {delete slot[i]; slot[i]=NULL;}
  num_keys=0;
  max_len=0;
  curr=prev=num_slots;
}





///////////////////////////////////////////////////////////////////////////////////////////
//                  Methods that work with an internal "current key":


////////////////////////////////////////////////////////////////////////////////////////////
// Return data of next key. Return 'fail' data and empty key if at end  
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Typ Hash<Typ>::ReadNext()
{
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  if (curr>=num_slots) {return fail;}
  pslot = slot[curr];  // current list is never empty, except when current=num_slots
  pairp = pslot->ReadNextAddress(); 
  if (pslot->End()) {
    prev=curr;
    do   // move on to next non-empty list
      {
	if (++curr>=num_slots) return pairp->data;
	pslot = slot[curr];
      } while (!pslot);
    pslot->Reset();
  }
  return pairp->data;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Write next key into variable key and return data. Return 'fail' data and empty key if at end  
// Attention: 'key' must have memory of at least char[MaxLen()+1] allocated!
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Typ Hash<Typ>::ReadNext(char* key)
{
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  if (curr>=num_slots) {*key='\0'; return fail;}
  pslot = slot[curr];  // current list is never empty, except when current=num_slots
  pairp = pslot->ReadNextAddress(); 
  strcpy(key,pairp->key);
  if (pslot->End()) {
    prev=curr;
    do   // move on to next non-empty list
      {
	if (++curr>=num_slots) return pairp->data;
	pslot = slot[curr];
      } while (!pslot);
    pslot->Reset();
  }
  return pairp->data;
}



////////////////////////////////////////////////////////////////////////////////////////////
// Return data of current key 
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Typ Hash<Typ>::ReadCurrent()
{
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  curr=prev;
  if (curr>=num_slots) {return fail;}
  pslot = slot[curr];  // current list is never empty, except when current=num_slots
  Pair<Typ> pair = pslot->ReadCurrent();
  pairp = &pair; 
  if (pslot->End()) {
    do   // move on to next non-empty list
      {
	if (++curr>=num_slots) return pairp->data;
	pslot = slot[curr];
      } while (!pslot);
    pslot->Reset();
  }
  return pairp->data;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Write key last read into variable key and return data
// Attention: 'key' must have memory of at least char[MaxLen()+1] allocated!
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Typ Hash<Typ>::ReadCurrent(char* key)
{
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  curr=prev;
  if (curr>=num_slots) {*key='\0'; return fail;}
  pslot = slot[curr];  // current list is never empty, except when current=num_slots
  Pair<Typ> pair = pslot->ReadCurrent();
  pairp = &pair; 
  strcpy(key,pairp->key);
  if (pslot->End()) {
    do   // move on to next non-empty list
      {
	if (++curr>=num_slots) return pairp->data;
	pslot = slot[curr];
      } while (!pslot);
    pslot->Reset();
  }
  return pairp->data;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Remove current key, return data, and advance to next key (after Reset() remove first element)
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Typ Hash<Typ>::RemoveCurrent()
{
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;
  curr=prev;

  if (curr>=num_slots) {return fail;}
  pslot = slot[curr];  // current list is never empty, except when current=num_slots
  Pair<Typ> pair = pslot->Delete();
  pairp = &pair; 
  num_keys--;
  // if key was the only element in pslot then delete whole list
  if (pslot->Size()==0) {delete pslot; slot[curr]=0;}  
  if (!pslot || pslot->End()) {
    do   // move on to next non-empty list
      {
	if (++curr>=num_slots) {prev=curr; return pairp->data;}
	pslot = slot[curr];
      } while (!pslot);
    pslot->Reset();
  }
  prev=curr;
  return pairp->data;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Remove current key, return data, copy current key into key, and advance to next key 
// (After Reset() remove first element)
// Attention: 'key' must have memory of at least char[MaxLen()+1] allocated!
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
Typ Hash<Typ>::RemoveCurrent(char* key)
{
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  curr=prev;
  if (curr>=num_slots) {*key='\0'; return fail;}
  pslot = slot[curr];  // current list is never empty, except when current=num_slots
  pairp = &(pslot->Delete()); 
  strcpy(key,pairp->key);
  num_keys--;
  // if key was the only element in pslot then delete whole list
  if (pslot->Size()==0) {delete pslot; slot[curr]=0;}  
  if (!pslot || pslot->End()) {
    do   // move on to next non-empty list
      {
	if (++curr>=num_slots) {prev=curr; return pairp->data;}
	pslot = slot[curr];
      } while (!pslot);
    pslot->Reset();
  }
  prev=curr;
  return pairp->data;
}



////////////////////////////////////////////////////////////////////////////////////////////
// Reset current position to beginning of hash
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
void Hash<Typ>::Reset()
  {
    curr=-1;
    Slot<Typ>* pslot;
    do
      {
	curr++;
	if (curr>=num_slots) {prev=curr; return;}
	pslot = slot[curr];
      } while (!pslot);
     pslot->Reset();
     prev=curr;
     return;
  }


/////////////////////////////////////////////////////////////////////////////////////////////
//            Methods that return usefull information about the data stored in Hash:


////////////////////////////////////////////////////////////////////////////////////////////
// Returns 1 if the hash contains key, 0 otherwise 
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
int Hash<Typ>::Contains(char* key)
{
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot) return 0;
  pslot->Reset();
  do{
    if(!strcmp(pslot->ReadNext().key,key)) return 1;
  } while(!pslot->End());
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//print out list of keys and data
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
void Hash<Typ>::Print()
{
  char key[MaxLen()+1];

  cout<<"\nPrint hash:\n";
  Reset();
  while(!End()) 
    cout<<key<<"->"<<ReadNext(key)<<"\n";
}
/////////////////////////////////////////////////////////////////////////////////////////////
//print out list of keys and data
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
void Hash<Typ>::PrintKeys()
{
  char key[MaxLen()+1];

  cout<<"\nPrint hash-keys:\n";
  Reset();
  while(!End()) 
    cout<<key<<"\n";
}


/////////////////////////////////////////////////////////////////////////////////////////////
//Print out hash with internal representation as array 
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ> 
void Hash<Typ>::DebugPrint()
{
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  cout<<"\n";
  cout<<"Debug-print hash:";
  for(int i=0; i<num_slots; i++)
    {
      pslot = slot[i];
      if (pslot) 
	{
	  cout<<"\nhash value "<<i;
	  pslot->Reset();
	  while(!pslot->End())
	    {
	      pairp = pslot->ReadNextAddress(); 
	      cout<<"  "<<pairp->key<<"->"<<pairp->data;
	    }
	}
    }
  cout<<"\n\n";
  return;
}



#endif /* HASH */



////////////////////////////////////////////////////////////////////////////////
// Main program: test class Hash
////////////////////////////////////////////////////////////////////////////////
//  int main()
//  {
//    Hash<int> ihash(5);
//    char* key=new char[128];
//    int data;

//    ihash.Fail(1000);
//    ihash.Add("So many monsters",36);
//    ihash.Add("So many chickens",25);

// //    cerr<<"Address of ihash(\"So many monsters\") is "<<ihash("So many monsters")<<"\n";
// //    cerr<<"Address of ihash(\"So many monsters\") is "<<ihash("So many monsters")<<"\n";
// //    *ihash("So many monsters")=2;
// //    (*ihash("So many monsters"))++;
// //    (*ihash("So many monsters"))++;
//    cout<<"Size of hash="<<ihash.Size()<<"\n";
//    ihash.DebugPrint();
//    ihash.Print();
   
//    strcpy(key,"Some more monsters");
//    ihash.Add(key,3);
//    strcpy(key,"Even more monsters");
//    ihash.Add(key,7);
//    cout<<"Size of hash="<<ihash.Size()<<"\n";
//    cout<<"Maximum key length = "<<ihash.MaxLen()<<"\n";
//    ihash.Print();
   
//    cout<<"ihash.Remove(\"Even more monsters\") returns "<<ihash.Remove("Even more monsters")<<"\n";
//    ihash.Print();


//    cout<<"ihash.Remove(\"More monsters\") returns "<<ihash.Remove("More monsters")<<"\n";
//    ihash.Add("So many chickens",999);
//    ihash.Add("So many monster",1);
//    ihash.Print();

//    cout<<"So many chickens:"<<ihash.Show("So many chickens")<<"\n";
//    cout<<"Size of hash="<<ihash.Size()<<"\n";



//    ihash.Reset();
//    while (!ihash.End())
//      {
//        data = ihash.ReadNext(key);
//        cout<<" "<<key<<"->"<<data<<endl;
//        data = ihash.ReadCurrent(key);
//        cout<<" "<<key<<"->"<<data<<endl;
//      }
//    cout<<"Removing all keys..."<<endl;
//    cout<<"Size of hash="<<ihash.Size()<<"\n";
//    ihash.Reset();
//    while (ihash.Size())
//      {
//        data = ihash.RemoveCurrent(key);
//        cout<<" "<<key<<"->"<<data<<"\n";
//        cout<<"Size of hash="<<ihash.Size()<<"\n";
//      }

//    ihash.ReadCurrent(key);

//    ihash.Print();
//    return 0;
//  }
