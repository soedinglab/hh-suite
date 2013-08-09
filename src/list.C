// list.C
// Class for double-linked list 
#ifndef JLIST
#define JLIST

#ifndef MAIN
#include <iostream>   // cin, cout
#include <stdlib.h>   // 
#include <stdio.h>    // 
using std::cout;
using std::cerr;
#endif

#include "list.h"

////////////////////////////////////////////////////////////////////////////////
// Double-linked list implementation with head and tail dummy elements
// We set head->prev=head and tail->next=tail.
// This makes sure that repeated current=current->next; ends up in tail
// and repeated current=current->prev; ends up in head.
// head and tail optionally contain a NULL element of Typ defined by method Null(Typ)
////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////
// Constructors and destructor

////////////////////////////////////////////////////////////////////////////
// Creates empty List with two dummy elements
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
List<Typ>::List()
{
  head=new ListEl<Typ>();
  if (!head) { cerr<<"Could not create new element\n"; return; }
  tail=new ListEl<Typ>(head,NULL);
  if (!tail) { cerr<<"Could not create new element\n"; return; }
  tail->next = tail;
  head->prev = head;
  head->next = tail;
  current = head;
  size=0;
}

////////////////////////////////////////////////////////////////////////////
// Creates List with one element 
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
List<Typ>::List(Typ d)
{
  head=new ListEl<Typ>();
  if (!head) { cerr<<"Could not create new element\n"; return; }
  tail=new ListEl<Typ>();
  if (!tail) { cerr<<"Could not create new element\n"; return; }
  ListEl<Typ>* el = new ListEl<Typ>(d,head,tail);
  if (!el)   { cerr<<"Could not create new element\n"; return; }
  head->prev = head;
  head->next = el;
  tail->prev = el;
  tail->next = tail;
  current = head;
  size=1;
}

////////////////////////////////////////////////////////////////////////////
// Destructor deletes List object
// Note: If <Typ> data is a pointer to another data structure, that structure is not deleted!
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
List<Typ>::~List()
{
  ListEl<Typ>* n=head->next;
  while(head!=n)
    {
      delete(head);
      head=n;
      n=head->next;
  }
  delete(head);
  size=0;
}

////////////////////////////////////////////////////////////////////////////
// Flat copy
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
inline List<Typ>& List<Typ>::operator=(List<Typ>& l)
{
  head = l.head;
  tail = l.tail;
  current = l.current;
  size = l.size;
}


////////////////////////////////////////////////////////////////////////////
// Reverse order of list
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
void List<Typ>::Reverse()
{
  ListEl<Typ> *n; // next list element; also for swapping
  ListEl<Typ> *c; // current element to be sorted in. Everything to the left is already sorted
  if (Size()<=1) return;
  for (c=head; c!=tail; c=n)
    {
      // Swap prev and next pointers of all list elements
      n = c->next;
      c->next = c->prev;
      c->prev = n;
    }
  
  // Swap prev and next pointers of tail
  tail->next = tail->prev;
  tail->prev = tail;

  // Swap head an tail
  n = head;
  head = tail;
  tail = n;
}




////////////////////////////////////////////////////////////////////////////////////////////
// Methods that act at the end of the list 

////////////////////////////////////////////////////////////////////////////
// Insert Element after LAST of list  (and return address of data element)
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
Typ* List<Typ>::Push(Typ d)
{ 
  ListEl<Typ>* t=new ListEl<Typ>(d,tail->prev,tail); 
  if (!t) { cerr<<"Could not create new element\n"; return 0; }
  tail->prev->next=t;
  tail->prev = t;
  size++;
  return &(t->data);
} 

////////////////////////////////////////////////////////////////////////////
// Remove and return LAST element of list. Returns head->data if empty
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
Typ List<Typ>::Pop()
{
  if (!size) return head->data;
  ListEl<Typ>* t=tail->prev;
  if (current==t) current=tail;
  Typ d=t->data;
  t->prev->next=tail;
  tail->prev=t->prev;
  delete t;
  size--;
  return d;
}



////////////////////////////////////////////////////////////////////////////////////////////
// Methods that act at the beginning of the list

////////////////////////////////////////////////////////////////////////////
// Insert element as FIRST element of list (and return address of data element)
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
Typ* List<Typ>::Enqueue(Typ d)
{
  ListEl<Typ>* h = new ListEl<Typ>(d,head,head->next); 
  if (!h) { cerr<<"Could not create new element\n"; return 0; }
  h->next->prev = h;
  head->next=h;
  size++;
  return &(h->data);
} 

////////////////////////////////////////////////////////////////////////////
// Remove element at BEGINNING of list
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
Typ List<Typ>::Dequeue()
{
  if (!size) return head->data;
  ListEl<Typ>* h=head->next;
  if (current==h) current=head;
  Typ d=h->data;
  h->next->prev=head;
  head->next=h->next;
  delete h;
  size--;
  return d;
} 

////////////////////////////////////////////////////////////////////////////////////////////
// Methods that work with 'current' position in the list

////////////////////////////////////////////////////////////////////////////
// Reads next element; advances current position by 1
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
inline Typ List<Typ>::ReadNext()
{
  current = current->next;
#ifdef HHDEBUG
  if (current == tail) cerr<<"WARNING in list.C function ReadNext(): Attempting to read tail element. Urgh!"<<endl;
#endif
  return current->data;
}

////////////////////////////////////////////////////////////////////////////
// Reads current element again (NULL if nothing read yet)
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
inline Typ List<Typ>::ReadCurrent()
{
#ifdef HHDEBUG
  if (current == head) cerr<<"WARNING in list.C function ReadCurrent(): Attempting to read head element. Urgh!"<<endl;
  if (current == tail) cerr<<"WARNING in list.C function ReadCurrent(): Attempting to read tail element. Urgh!"<<endl;
#endif
  return current->data;
}

////////////////////////////////////////////////////////////////////////////
// Reads previous element; moves current position back by 1
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
inline Typ List<Typ>::ReadPrevious()
{
  current = current->prev;
#ifdef HHDEBUG
  if (current == head) cerr<<"WARNING in list.C function ReadPrevious(): Attempting to read head element. Urgh!"<<endl;
#endif
  return current->data;
}

////////////////////////////////////////////////////////////////////////////
// Reads address next data element; advances current position by 1
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
inline Typ* List<Typ>::ReadNextAddress()
{
  current = current->next;
  if (current==tail) return NULL;
  return &(current->data);
}

////////////////////////////////////////////////////////////////////////////
// Reads address of data element if current element again, returns NULL if at end of list
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
inline Typ* List<Typ>::ReadCurrentAddress()
{
  if (current==tail) return NULL;
  return &(current->data);
}

////////////////////////////////////////////////////////////////////////////
// Sets current position to k and reads k'th element (first=1)
////////////////////////////////////////////////////////////////////////////
template <class Typ>
Typ List<Typ>::Read(int pos)
{
  if (pos>size) {current = tail; return tail->data;}
  if (pos<=0)   {current = head; return head->data;}
  current = head->next;
  for (; pos>1; pos--) current = current->next; //If pos==2 do 1 iteration
  return current->data;
}

////////////////////////////////////////////////////////////////////////////
// Inserts element d AFTER current element and sets current element to inserted
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
void List<Typ>::Insert(Typ d)
{
  ListEl<Typ>* el = new ListEl<Typ>(d,current,current->next); 
  if (!el) { cerr<<"Could not create new element\n"; return; }
  (current->next)->prev = el;
  current->next = el;
  current=el;
  size++;
}

////////////////////////////////////////////////////////////////////////////
// Deletes current element and returns content of deleted element. Current element 
// will be previous one after Delete(). After Reset() delete first element (not 0'th)
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
Typ List<Typ>::Delete()
{
  Typ d;
  ListEl<Typ>* p;
  if (!size || current==tail) return tail->data;
  if (current==head) current = head->next; // After Reset() delete first element (not 0'th)     
  (current->prev)->next = current->next;
  (current->next)->prev = current->prev;
  d = current->data;
  p = current->prev;
  delete current;
  current = p;
  size--;
  return d;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Methods that return useful information about the list

////////////////////////////////////////////////////////////////////////////
// Get current position within list (0 <= pos <= Size+1) 
////////////////////////////////////////////////////////////////////////////
template <class Typ>
int List<Typ>::GetPos()
{
  int pos=0;
  ListEl<Typ>* el;
  for (el = head; el!=current; el=el->next) pos++; 
  return pos;
}

////////////////////////////////////////////////////////////////////////////
//print out list
////////////////////////////////////////////////////////////////////////////
template <class Typ>
void List<Typ>::PrintList()
{
  int j=0;
  ListEl<Typ>* c=current;
  Reset();
  printf("List: ");
  while (!End()) 
    {
      j++;
      cout<<j<<" "<<ReadNext()<<"  ";
      if (!(j%10)) cout<<"\n      ";
    }
  cout<<"\n";
  current=c;
}

////////////////////////////////////////////////////////////////////////////
// Get largest data element
////////////////////////////////////////////////////////////////////////////
template <class Typ>
Typ List<Typ>::Largest()
{
  Typ* result= &((tail->prev)->data);
  Reset();
  while (!End()) 
    {
      if (*result<ReadNext()) result=ReadCurrentAddress();
    }
  return *result;
}

////////////////////////////////////////////////////////////////////////////
// Get smallest data element
////////////////////////////////////////////////////////////////////////////
template <class Typ>
Typ List<Typ>::Smallest()
{
  Typ* result= &((tail->prev)->data);
  Reset();
  while (!End()) 
    {
      if (ReadNext()<*result) result=ReadCurrentAddress();
    }
  return *result;
}

////////////////////////////////////////////////////////////////////////////////////////////
// Methods that manipulate the list as a whole

////////////////////////////////////////////////////////////////////////////
// Copies list 0 into list object
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
void List<Typ>::Copy(List<Typ>* list)
{
  if (list==this) return;
  while (!End()) Pop(); //empty list
  list->Reset();
  while (!list->End()) Push(list->ReadNext());
}

////////////////////////////////////////////////////////////////////////////
// Appends a copy of list2 to class object
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
void List<Typ>::AppendCopy(List<Typ>* list2)
{
  List<Typ>* cpy=new List<Typ>;
  cpy->Copy(list2);
  Append(cpy);
  delete cpy;
}

////////////////////////////////////////////////////////////////////////////
// Appends list2 to class object
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
void List<Typ>::Append(List<Typ>* list)
{
  if (this==list) { AppendCopy(list); return;}
  (tail->prev)->next = list->head->next;
  (list->head->next)->prev = tail->prev;
  if (current==tail) current=tail->prev; 
  ListEl<Typ>* t=tail;
  tail = list->tail;
  size += list->size;

// Reuse old tail as new tail t for list2 and initialize pointers for empty list 
  list->tail=t;
  list->head->next=t;
  t->prev=list->head;
  t->next=t;
  list->head->prev=list->head;
  t->data=list->head->data;
  list->current=list->head;
  list->size = 0;
}

////////////////////////////////////////////////////////////////////////////
// Use QUICKSORT to sort list in ascending order between two list elements
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
void List<Typ>::SortList(ListEl<Typ>* left, ListEl<Typ>* right, int sz) 
{
  if (sz<=1) return; // when SortList() is called, left=head->next, right=tail->prev
  ListEl<Typ> *l=left->prev, *r=right->next;

  // Choose *random* pivot element!! 
  // (Otherwise, complexity for an already sorted list is N^2 => recursive calls may lead to stack overflow)
  ListEl<Typ> *c=left;
  for (int i=1; i<(int)(float(rand())*sz/(RAND_MAX+0.999)); i++) c = c->next; 
  SwapContent(left,c);

  Typ pivot = left->data;
//   Typ* pivot= &(left->data);
  int sz0=sz+1;
  //  cout<<"Sorting between "<<left->data<<" and "<<right->data<<". Pivot="<<pivot<<endl;
  while(1)
    {
//        PrintList();
      do {r=r->prev; sz0--;} while (pivot < r->data); 
      do l=l->next; while (l->data < pivot);
      if (l==r || l->prev==r) break;
      SwapContent(l,r);
    }
  SortList(left,r,sz0);
  SortList(r->next,right,sz-sz0);
  pivot = tail->data; // to avoid calling the destructor of Typ on some real data element
 }

////////////////////////////////////////////////////////////////////////////
// Use QUICKSORT to sort list of POINTERS by comparing the objects the pointers point to
////////////////////////////////////////////////////////////////////////////
template <class Typ> 
void List<Typ>::SortPointerList(ListEl<Typ>* left, ListEl<Typ>* right) 
{
  if (right==left || right->next==left) return;
  ListEl<Typ> *l=left->prev, *r=right->next;
  Typ pivot=left->data;
//    cout<<"Sorting between "<<left->data<<" and "<<right->data<<". Pivot="<<pivot<<endl;
  while(1)
    {
//        PrintList();
      do
	{
	  r=r->prev;
//  	  cout<<"r=r->prev. r->data="<<r->data<<endl;
	} while(*pivot < *(r->data)); 
      do
	{
	  l=l->next;
//  	  cout<<"l=l->next l->data="<<l->data<<endl;
	} while (*(l->data) < *pivot);
      if (l==r || l->prev==r) break;
      SwapContent(l,r);
    }
  SortPointerList(left,r);
  SortPointerList(r->next,right);
}

// Use INSERTSORT to sort list in asscending order between two list elements. Use only for presorted lists, otherwise time O(N^2)!
template <class Typ> 
void List<Typ>::ResortList() 
{
  ListEl<Typ> *c; // current element to be sorted in. Everything to the left is already sorted
  ListEl<Typ> *n; // next element to be sorted in
  ListEl<Typ> *p; // pointer for looping through sorted part of list
  ListEl<Typ> *pnext; // for swapping
  if (Size()<=1) return;
  c=head->next->next;
  while (c!=tail) 
    {
      p=c->prev;
      n=c->next;
      if (c->data < p->data)
	{
	  do {p=p->prev;} while (p!=head && c->data < p->data);
	  // Connect c->prev to c->next ...
	  c->next->prev=c->prev;
	  c->prev->next=c->next;
	  // ... and insert c between p and p->next ...
	  pnext=p->next;
	  p->next=c;
	  c->next=pnext;
	  pnext->prev=c;
	  c->prev=p;
	}
      c=n;
    }
}


#endif /* JLIST */


// //Main program: test class List

//  int main()
//  {
//    int p;
//    List<int>* plist=new List<int>(11);
//    List<int> list(22);

//    plist->Push(24);
//    plist->Push(18);
//    plist->Push(3);
//    plist->Enqueue(17);
//    plist->Enqueue(29);
//    printf("List 1 with pushed and enqueued elements:\n");
//    plist->PrintList();

//    list.Push(222);
//    printf("List 1 with list 2 appended:\n");
//    plist->Append(&list);
//    plist->PrintList();

//    printf("Pushing one element three times into list 2:\n");
//    list.Push(333);
//    list.Push(444);
//    list.Push(555);
//    printf("Printing plist and list with three elements:\n");
//    list.PrintList();
//    plist->PrintList();
  
//    printf("list.Copy(plist). Printing list 1 and 2:\n");
//    list.Copy(plist);
//    plist->PrintList(); 
//    list.PrintList();

//    printf("Appending list 1 to itself:\n");
//    plist->Append(plist);
//    plist->PrintList();

//    cout<<"Popping "<<plist->Pop()<<"\n";
//    cout<<"Popping "<<plist->Pop()<<"\n";
//    plist->PrintList();

//    cout<<"Dequeing "<<plist->Dequeue()<<"\n";
//    cout<<"Dequeing "<<plist->Dequeue()<<"\n";
//    plist->PrintList();

//    cout<<"Reversing list\n";
//    plist->Reverse();
//    plist->PrintList();

//    cout<<"Reversing to original list\n";
//    plist->Reverse();
//    plist->PrintList();

//    for (p=plist->Reset(); p>=5;p--)
//      {cout<<plist->GetPos()<<": "<<plist->Read(p)<<"\n";}
  
//    cout<<"Deleting "<<plist->Delete()<<"\n";
//    cout<<"Deleting "<<plist->Delete()<<"\n";
//    plist->PrintList();

//    plist->Append(plist);
//    plist->PrintList();
//    cout<<"List 1 sorted:\n";
//    plist->SortList();
//    plist->PrintList();

//  }
