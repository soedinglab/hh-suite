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
  int num_keys;          //total number of keys in hash
  int max_len;           //length of longest key in hash
  int key_len;           //length of key in argument
  Typ fail;
  
  Slot<Typ>** slot;      //each slot[i] (i<num_slots) points to a list of key/data pairs for this slot

  inline unsigned int HashValue(char* key);  //returns the hash value for key 

public:
  Hash();
  Hash(int nslots);
  Hash(int nslots, Typ n);
  ~Hash();

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
  Typ Show(char* key);
  inline Typ operator[](char* key) {return Show(key);}

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
  int Contains(char* key);

  // Return number of slots
  int Size()  {return num_keys;}  

  // Return length of longest key INCLUDING DELETED KEYS (excluding \0)
  int MaxLen()  {return max_len;}  

  //print out list of keys and data
  void Print();

  //Print out hash with internal representation as array 
  void DebugPrint();
};

