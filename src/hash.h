// Class for Hash data structure 
// * works in the same way as a hash in Perl
// * keys are strings of type char*
// * data elements are of type Typ 
// * objects have to be declared with maximal size, e.g. Hash hash1(10000) (num_slots should not be a power of 2)
// * works also if maximal size is exceeded, but gets slower by a factor ~num_keys/num_slots 
//
// Applications
// * fast storage and retrieval of data by a key
// * fast check if a string occurs in a list of strings (data field not used)
//
// Time complexity: (L is the length of the key string)
// * Show(key), Add(key), Remove(key): O(L) for calculating hash value & compare keys in slot list
// * ReadNext, ReadCurrent, RemoveCurrent: O(L) for copying key into returned string
// * Contains: O(L) for calculating hash value
// * Size, MaxLen, Reset: O(1)
// * RemoveAll(), Hash() and ~Hash(): O(num_slots)
//
// Memory complexity: ~3*num_keys*(datasize+12bytes) + num_slots*4bytes + total added length of keys (including \0)
//
// Implementation:
// Hash is an array of pointers to lists of size num_slots. The lists, called slots, contain key/data pairs. 
// When a key/data pair is added (e.g. with Add()) the array index i for the key (0<=i<num_slots) 
// is calculated with the HashValue() function. 
// When data is to be retrieved for key (e.g. with Show()) the hash value for key is calculated and
// the corresponding list (=slot) is searched for the occurence of the key. 
// Array elements of hash values that have no keys associated yet contain the null pointer 
// for faster rejection of undefined keys.
// 


////////////////////////////////////////////////////////////////////////////////////////////
// Declaration of Pair, consisting of a key string and a data element of type Typ 
////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HASH_H_
#define HASH_H_

#include <cstring>
#include <iostream>
#include "list.h"


template<class Typ> 
class Pair
{
public:
  char* key;             //hash key
  Typ data;              //data for key
  Pair() {}
  Pair(char* k, Typ& d) {key = new char[strlen(k)+1]; strcpy(key,k); data=d;}
  Pair(int& l, char* k, Typ& d) {key = new char[l+1]; strcpy(key,k); data=d;}
};


////////////////////////////////////////////////////////////////////////////////////////////
// Declaration of Slot, a list of key/data pairs
////////////////////////////////////////////////////////////////////////////////////////////
 
template<class Typ>
class Slot : public List< Pair<Typ> >
{
public:
  //Destructor of Slot deletes whole list TOGETHER WITH THE KEY STRINGS
  ~Slot() {this->Reset(); while (!this->End()) delete[] this->Pop().key; }
  
  // Push key/data pair onto slot list and return address of data element
  // Attention: l must be at least length of key
  inline Typ* Push(int& l, char* key, Typ& data)
  {
    Pair<Typ> pair(l,key,data);  //create a pair with key/data
    return &(List<Pair<Typ> >::Push(pair)->data);
  }
};

 
////////////////////////////////////////////////////////////////////////////////////////////
// Declaration of Hash, an array of slots, i.e. pointers to lists of key/data 
////////////////////////////////////////////////////////////////////////////////////////////

template<class Typ> 
class Hash  
{
private:
  int num_slots;         //number of slots in slot[]. num_slots is set with the constructor
  int curr;              //index of current slot
  int prev;              //index of slot from previous ReadNext() 
  size_t num_keys;          //total number of keys in hash
  int max_len;           //length of longest key in hash
  int key_len;           //length of key in argument
  Typ fail;
  
  Slot<Typ>** slot;      //each slot[i] (i<num_slots) points to a list of key/data pairs for this slot

  inline unsigned int HashValue(const char* key);  //returns the hash value for key

public:
  Hash();
  Hash(int nslots);
  Hash(int nslots, Typ n);
  ~Hash(); // note: if <Typ> data is a pointer to another data structure, that structure is not deleted!

  // Set Fail element to be returned when the current key or supplied key are not defined
  inline void Null(Typ f) {fail=f;}
  inline void Fail(Typ f) {fail=f;}

  // Set 'fail' element to be returned when the current key or supplied key are not defined (default: 0)
  void New(int nslots, Typ n=static_cast<Typ>(0)); 

  // Update maximum key length and caculate key_len;
  inline void KeyLen() {if(key_len>max_len) max_len=key_len; return;}


////////////////////////////////////////////////////////////////////////////////////////////
//                 Methods that work with a key supplied as an argument

  // Return data element for key. Returns 'fail' if key does not exist
  Typ Show(const char* key);
  inline Typ operator[](const char* key) {return Show(key);}

  // Add/replace key/data pair to hash and return address of data element for key
  Typ* Add(char* key, Typ data);

  // Add key to hash and return address of data element. If key exists leave data element unchanged, else set it to 'fail'.
  Typ* Add(char* key);
  inline Typ* operator()(char* key) {return Add(key);} 

  // Remove key from hash and return data element for key ('fail' if key does not exist)
  Typ Remove(char* key);

  // Remove all keys from hash
  void RemoveAll();


///////////////////////////////////////////////////////////////////////////////////////////
//                  Methods that work with an internal "current key":
// It is set to the first key by Reset() and moves to the next key with ReadNext or RemoveCurrent
// Note:the methods above (e.g. Store, Show, [], Add, (), etc. DO NOT CHANGE the current key

  // Return data of next key. Return 'fail' data and empty key if at end  
  Typ ReadNext();

  // Write next key into variable key and return data. Return 'fail' data and empty key if at end  
  // Attention: 'key' must have memory of at least char[MaxLen()+1] allocated!
  Typ ReadNext(char* key);

  // Return data of current key 
  Typ ReadCurrent();

  // Write key last read into variable key and return data
  // Attention: 'key' must have memory of at least char[MaxLen()+1] allocated!
  Typ ReadCurrent(char* key);

  // Remove current key, return data and advance to next key 
  Typ RemoveCurrent();

  // Remove current key, return data, copy current key into key, and advance to next key 
  // (After Reset() remove first element)
  // Attention: 'key' must have memory of at least char[MaxLen()+1] allocated!
  Typ RemoveCurrent(char* key);

  // Reset readout of keys to beginning of hash
  void Reset(); 

  // Returns 1 if the current key has arrived at the end, 0 otherwise
  int End() {return (curr>=num_slots);}



///////////////////////////////////////////////////////////////////////////////////////////
//            Methods that return usefull information about the data stored in Hash:

  // Returns 1 if the hash contains key, 0 otherwise 
  int Contains(const char* key);

  // Return number of slots
  size_t Size()  {return num_keys;}

  // Return length of longest key INCLUDING DELETED KEYS (excluding \0)
  int MaxLen()  {return max_len;}  

  //print out list of keys and data
  void Print();

  //print out list of keys
  void PrintKeys();

  //Print out hash with internal representation as array 
  void DebugPrint();
};

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
Hash<Typ>::Hash() {
  num_keys = 0;
  max_len = 0;
  prev = curr = num_slots = 0;
  slot = NULL;
}

template<class Typ>
Hash<Typ>::Hash(int nslots) {
  num_keys = 0;
  max_len = 0;
  prev = curr = num_slots = nslots;
  slot = new Slot<Typ>*[num_slots];            //Create array of num_slots slots
  for (int i = 0; i < num_slots; i++)
    slot[i] = NULL; //set pointers to NULL
  fail = static_cast<Typ>(0);
}

template<class Typ>
Hash<Typ>::Hash(int nslots, Typ f) {
  num_keys = 0;
  max_len = 0;
  prev = curr = num_slots = nslots;
  slot = new Slot<Typ>*[num_slots];            //Create array of num_slots slots
  for (int i = 0; i < num_slots; i++)
    slot[i] = NULL; //set pointers to NULL
  fail = f;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Destructor of class Hash
// Note: if <Typ> data is a pointer to another data structure, that structure is not deleted!
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Hash<Typ>::~Hash() {
  RemoveAll();
  delete[] slot;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Hash function
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
inline unsigned int Hash<Typ>::HashValue(const char* key) //returns the hash value for key
    {
  // Calculate a hash value by the division method:
  // Transform key into a natural number k = sum ( key[i]*128^(L-i) ) and calculate i= k % num_slots.
  // Since calculating k would lead to an overflow, i is calculated iteratively
  // and at each iteration the part divisible by num_slots is subtracted, i.e. (% num_slots is taken).
  if (key == NULL) {
    printf("Warning from hash.cpp: key=NULL\n");
    return 0;
  }
  unsigned int i = 0;     // Start of iteration: k is zero
  const char* c = key;
  while (*c)
    i = ((i << 7) + *(c++)) % num_slots;
  key_len = c - key;
  //cerr<<"      Hash value for \'"<<key<<"\' is "<<i<<"\n";
  return i;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Create new hash (and delete any data present)
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
void Hash<Typ>::New(int nslots, Typ f) {
  fail = f;
  RemoveAll();
  delete[] slot;
  num_keys = 0;
  max_len = 0;
  prev = curr = num_slots = nslots;
  slot = new Slot<Typ>*[num_slots];            //Create array of num_slots slots
  for (int i = 0; i < num_slots; i++)
    slot[i] = NULL; //set pointers to NULL
}

////////////////////////////////////////////////////////////////////////////////////////////
//                  Methods that work with a key supplied as an argument

////////////////////////////////////////////////////////////////////////////////////////////
// Return data element for key. Returns 'fail' if key does not exist
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Typ Hash<Typ>::Show(const char* key) {
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot)
    return fail;
  pslot->Reset();
  while (!pslot->End()) {
    if (!strcmp(pslot->ReadNext().key, key))
      return pslot->ReadCurrent().data;
  }
  return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Add/replace key/data pair to hash and return address of data element
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Typ* Hash<Typ>::Add(char* key, Typ data) {
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot) {
    num_keys++;
    KeyLen();
    slot[i] = new (Slot<Typ> );
    return slot[i]->Push(key_len, key, data);
  }
  pslot->Reset();
  while (!pslot->End()) {
    pairp = pslot->ReadNextAddress();
    if (!strcmp(pairp->key, key)) {
      pairp->data = data;
      pslot->Overwrite(*pairp);
      return &(pairp->data);
    }
  }
  num_keys++;
  KeyLen();
  return pslot->Push(key_len, key, data);
}

////////////////////////////////////////////////////////////////////////////////////////////
// Add key to hash and return address of data element.
// If key exists leave data element unchanged, else set it to 'fail'.
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Typ* Hash<Typ>::Add(char* key) {
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot) {
    num_keys++;
    KeyLen();
    slot[i] = new (Slot<Typ> );
    return slot[i]->Push(key_len, key, fail);
  }
  pslot->Reset();
  while (!pslot->End()) {
    if (!strcmp(pslot->ReadNext().key, key)) {
      return &((pslot->ReadCurrentAddress())->data);
    }
  }
  num_keys++;
  KeyLen();
  return pslot->Push(key_len, key, fail);
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Remove key from hash and return data element for key ('fail' if key does not exist)
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Typ Hash<Typ>::Remove(char* key) {
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot)
    return fail;
  pslot->Reset();
  while (!pslot->End()) {
    if (!strcmp(pslot->ReadNext().key, key)) {
      Pair<Typ> pair = pslot->ReadCurrent();
      num_keys--;
      pslot->Delete();
      // Delete key-Array
      delete[] pair.key;
      // if key was the only element in pslot then delete whole list
      if (pslot->Size() == 0) {
        delete pslot;
        slot[i] = 0;
      }
      //	  return pslot->ReadCurrent().data;
      return pair.data;
    }
  }
  return fail;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Remove all keys from hash;
// Note: if <Typ> data is a pointer to another data structure, that structure is not deleted!
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
void Hash<Typ>::RemoveAll() {
  for (int i = 0; i < num_slots; i++)
    if (slot[i]) {
      delete slot[i];
      slot[i] = NULL;
    }
  num_keys = 0;
  max_len = 0;
  curr = prev = num_slots;
}

///////////////////////////////////////////////////////////////////////////////////////////
//                  Methods that work with an internal "current key":

////////////////////////////////////////////////////////////////////////////////////////////
// Return data of next key. Return 'fail' data and empty key if at end
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Typ Hash<Typ>::ReadNext() {
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  if (curr >= num_slots) {
    return fail;
  }
  pslot = slot[curr]; // current list is never empty, except when current=num_slots
  pairp = pslot->ReadNextAddress();
  if (pslot->End()) {
    prev = curr;
    do   // move on to next non-empty list
    {
      if (++curr >= num_slots)
        return pairp->data;
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
Typ Hash<Typ>::ReadNext(char* key) {
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  if (curr >= num_slots) {
    *key = '\0';
    return fail;
  }
  pslot = slot[curr]; // current list is never empty, except when current=num_slots
  pairp = pslot->ReadNextAddress();
  strcpy(key, pairp->key);
  if (pslot->End()) {
    prev = curr;
    do   // move on to next non-empty list
    {
      if (++curr >= num_slots)
        return pairp->data;
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
Typ Hash<Typ>::ReadCurrent() {
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  curr = prev;
  if (curr >= num_slots) {
    return fail;
  }
  pslot = slot[curr]; // current list is never empty, except when current=num_slots
  Pair<Typ> pair = pslot->ReadCurrent();
  pairp = &pair;
  if (pslot->End()) {
    do   // move on to next non-empty list
    {
      if (++curr >= num_slots)
        return pairp->data;
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
Typ Hash<Typ>::ReadCurrent(char* key) {
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  curr = prev;
  if (curr >= num_slots) {
    *key = '\0';
    return fail;
  }
  pslot = slot[curr]; // current list is never empty, except when current=num_slots
  Pair<Typ> pair = pslot->ReadCurrent();
  pairp = &pair;
  strcpy(key, pairp->key);
  if (pslot->End()) {
    do   // move on to next non-empty list
    {
      if (++curr >= num_slots)
        return pairp->data;
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
Typ Hash<Typ>::RemoveCurrent() {
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;
  curr = prev;

  if (curr >= num_slots) {
    return fail;
  }
  pslot = slot[curr]; // current list is never empty, except when current=num_slots
  Pair<Typ> pair = pslot->Delete();
  pairp = &pair;
  num_keys--;
  // if key was the only element in pslot then delete whole list
  if (pslot->Size() == 0) {
    delete pslot;
    slot[curr] = 0;
  }
  if (!pslot || pslot->End()) {
    do   // move on to next non-empty list
    {
      if (++curr >= num_slots) {
        prev = curr;
        return pairp->data;
      }
      pslot = slot[curr];
    } while (!pslot);
    pslot->Reset();
  }
  prev = curr;
  return pairp->data;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Remove current key, return data, copy current key into key, and advance to next key
// (After Reset() remove first element)
// Attention: 'key' must have memory of at least char[MaxLen()+1] allocated!
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
Typ Hash<Typ>::RemoveCurrent(char* key) {
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  curr = prev;
  if (curr >= num_slots) {
    *key = '\0';
    return fail;
  }
  pslot = slot[curr]; // current list is never empty, except when current=num_slots
  pairp = &(pslot->Delete());
  strcpy(key, pairp->key);
  num_keys--;
  // if key was the only element in pslot then delete whole list
  if (pslot->Size() == 0) {
    delete pslot;
    slot[curr] = 0;
  }
  if (!pslot || pslot->End()) {
    do   // move on to next non-empty list
    {
      if (++curr >= num_slots) {
        prev = curr;
        return pairp->data;
      }
      pslot = slot[curr];
    } while (!pslot);
    pslot->Reset();
  }
  prev = curr;
  return pairp->data;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Reset current position to beginning of hash
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
void Hash<Typ>::Reset() {
  curr = -1;
  Slot<Typ>* pslot;
  do {
    curr++;
    if (curr >= num_slots) {
      prev = curr;
      return;
    }
    pslot = slot[curr];
  } while (!pslot);
  pslot->Reset();
  prev = curr;
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//            Methods that return usefull information about the data stored in Hash:

////////////////////////////////////////////////////////////////////////////////////////////
// Returns 1 if the hash contains key, 0 otherwise
////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
int Hash<Typ>::Contains(const char* key) {
  Slot<Typ>* pslot;
  int i = HashValue(key);

  pslot = slot[i];
  if (!pslot)
    return 0;
  pslot->Reset();
  while (!pslot->End()) {
    if (!strcmp(pslot->ReadNext().key, key))
      return 1;
  }
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//print out list of keys and data
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
void Hash<Typ>::Print() {
  char key[MaxLen() + 1];

  std::cout << "\nPrint hash:\n";
  Reset();
  while (!End())
    std::cout << key << "->" << ReadNext(key) << "\n";
}
/////////////////////////////////////////////////////////////////////////////////////////////
//print out list of keys and data
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
void Hash<Typ>::PrintKeys() {
  char key[MaxLen() + 1];

  std::cout << "\nPrint hash-keys:\n";
  Reset();
  while (!End())
    std::cout << key << "\n";
}

/////////////////////////////////////////////////////////////////////////////////////////////
//Print out hash with internal representation as array
/////////////////////////////////////////////////////////////////////////////////////////////
template<class Typ>
void Hash<Typ>::DebugPrint() {
  Pair<Typ>* pairp;
  Slot<Typ>* pslot;

  std::cout << "\n";
  std::cout << "Debug-print hash:";
  for (int i = 0; i < num_slots; i++) {
    pslot = slot[i];
    if (pslot) {
      std::cout << "\nhash value " << i;
      pslot->Reset();
      while (!pslot->End()) {
        pairp = pslot->ReadNextAddress();
        std::cout << "  " << pairp->key << "->" << pairp->data;
      }
    }
  }
  std::cout << "\n\n";
  return;
}





#endif
