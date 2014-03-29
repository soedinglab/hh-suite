// hash.cpp
// Class for Hash data structure 
// * works in the same way as a hash in Perl
// * keys are strings of type char*
// * data elements are of type Typ 
// * objects have to be declared with maximal size, e.g. Hash hash1(10000) (num_slots should not be a power of 2)
// * works also if maximal size is exceeded, but gets slower by a factor ~num_keys/num_slots 

#include "list.h"
#include "hash.h"

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
