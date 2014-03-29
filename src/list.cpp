// list.C
// Class for double-linked list 
#include "list.h"


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
