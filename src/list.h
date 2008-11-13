// list.h

////////////////////////////////////////////////////////////////////////////////
// Double-linked list implementation with head and tail dummy elements
// We set head->prev=head and tail->next=tail.
// This makes sure that repeated current=current->next; ends up in tail
// and repeated current=current->prev; ends up in head.
// head and tail optionally contain a NULL element of Typ defined by method Null(Typ)
////////////////////////////////////////////////////////////////////////////////

template<class Typ> 
class List
{
protected:
template<class Typ1> 
class ListEl             //elements of List; essentially a data structure
  {
  public:
    Typ1 data;           //Typ is type of data to be stored in list
    ListEl* prev;        //points to previous list element
    ListEl* next;        //points to next list element
    ListEl() : prev(0), next(0) {}
    ListEl(Typ1 d) : data(d), prev(0), next(0) {}
    ListEl(ListEl* p, ListEl* n) : prev(p), next(n) {}
    ListEl(Typ1 d, ListEl* p, ListEl* n) : data(d), prev(p), next(n) {}
  };
  
  ListEl<Typ>* head;     //points to dummy element at beginning of list
  ListEl<Typ>* tail;     //points to dummy element at end of list    
  ListEl<Typ>* current;  //current element position within list
  int size;              //Number of elements in list

  // Use QUICKSORT to sort list in asscending order between two list elements
  void SortList(ListEl<Typ>*, ListEl<Typ>*, int); 
  // Use QUICKSORT to sort list of pointers by comparing elements they point to
  void SortPointerList(ListEl<Typ>*, ListEl<Typ>*);

  // Swap two list elements by making a flat copy (don't need two copies of data)
  // Warning: Gets slow if Typ is composite type with many variables (>=5)
  void SwapContent(ListEl<Typ>* e1, ListEl<Typ>* e2)
  { Typ d; if (e1!=e2) {d=e1->data; e1->data=e2->data; e2->data=d;} }

public:
////////////////////////////////////////////////////////////////////////////////////////////
// General methods
  List();
  List(Typ d);
  ~List();
  List<Typ>& operator=(List<Typ>&);
    
  // Set Null element that will be returned when trying to read from an empty list
  void Null(Typ null) {head->data = tail->data = null;}


////////////////////////////////////////////////////////////////////////////////////////////
// Methods that act at the end of the list 

  // Insert Element after LAST element of list (and return address of data element)
  Typ* Push(Typ);

  // Remove and return LAST element of list. Returns head->data if list empty
  Typ Pop();

  // return LAST element of list. Returns null element in head->data if list empty
  Typ ReadLast() {return tail->prev->data;}


////////////////////////////////////////////////////////////////////////////////////////////
// Methods that act at the beginning of the list

  // Insert element as FIRST element of list (and return address of data element)
  Typ* Enqueue(Typ);

  // Remove and return element at BEGINNING of list. Returns head->data if list empty
  Typ Dequeue();

  // return FIRST element of list. Returns null element in head->data if list empty
  Typ ReadFirst() {if (size) return head->next->data; else return head->data;}


////////////////////////////////////////////////////////////////////////////////////////////
// Methods that work with 'current' position in the list

  // Advances current position by 1 and reads next element; returns head->data if at end of list.
  Typ ReadNext(); 

  // Reads current element again
  Typ ReadCurrent(); 

  // Moves current position back by 1 and reads previous element; returns head->data if at beginning of list.
  Typ ReadPrevious(); 

  // Advances current position by 1 and reads address of next element; returns NULL if at end of list.
  Typ* ReadNextAddress(); 

  // Reads address of current element again, returns NULL if at end of list
  Typ* ReadCurrentAddress(); 

  // Sets current position to k and reads k'th element (first=1). Returns head->data if current points to no data element
  Typ Read(int);

  // Inserts element AFTER CURRENT element; current element will be set to inserted element
  void Insert(Typ);

  // Removes and returns element at CURRENT position. New position is one BEFORE current position. 
  // Returns head->data if current points to no data element. After Reset() delete first element (not 0'th)
  Typ Delete();

  // Overwrites data at current position with new data 
  void Overwrite(Typ d) {current->data=d;} 

  // Reset current position to 0 (one BEFORE the first)
  int Reset() {current = head; return size;} 

  // Reset current position to End (one AFTER the last)
  int SetToEnd() {current = tail; return size;} 


////////////////////////////////////////////////////////////////////////////////////////////
// Methods that return information about the list

  // Return number of list elements (size>=0)
  int Size()  {return size;}  

  // return true if end of list, i.e. ReadNext would give tail->data (i.e. current position >= Size)
  char End()  {return (current==tail || current==tail->prev);}
  char End(void* curr)  {return ( curr == tail || curr == tail->prev);}

  // return true if start of list, i.e. ReadPrevious would give head->data (i.e. current position <=1)
  char Start()  {return (current==head || current==head->next);}

  // Get current position within list (0 <= pos <= Size+1) 
  int GetPos();

  //print out list (elements assumed int)
  void PrintList();

  // Get largest data element (Null element for empty list)
  Typ Largest();

  // Get smallest data element (Null element for empty list)
  Typ Smallest();

////////////////////////////////////////////////////////////////////////////////////////////
// Methods that manipulate the list as a whole

  // Reverse list 
  void Reverse();

  // Copies list into list object
  void Copy(List<Typ>* list);

  // Appends a copy of list to class object
  void AppendCopy(List<Typ>* list);

  // Appends list to class object list
  void Append(List<Typ>* list); 

  // Use QUICKSORT to sort list in ascending order. Use only for UNSORTED lists, otherwise time O(N^2) instead of O(N*log(N))
/*   void SortList() {if (size>1) SortList(head->next, tail->prev);}  */
  void SortList() {if (size>1) SortList(head->next, tail->prev, size);} 
  void QuickSort() {if (size>1) SortList(head->next, tail->prev, size);} 

  // Use QUICKSORT to sort list of pointers in ascending order. Use only for UNSORTED lists, otherwwise time O(N^2)!
  void SortPointerList() {if (size>1) SortPointerList(head->next, tail->prev);} 
  void QuickSortPointer() {if (size>1) SortPointerList(head->next, tail->prev);} 

  // Use INSERTSORT to sort list in asscending order. Use only for PRESORTED lists, otherwise time O(N^2)!
  void ResortList(); 
};


