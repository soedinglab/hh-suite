// Utility subroutines

#include "util.h"




/////////////////////////////////////////////////////////////////////////////////////
// String functions
/////////////////////////////////////////////////////////////////////////////////////

void split(const std::string& s, char c, std::vector<std::string>& v) {
  std::string::size_type i = 0;
  std::string::size_type j = s.find(c);

  while (j != std::string::npos) {
    v.push_back(s.substr(i, j - i));
    i = ++j;
    j = s.find(c, j);
  }

  if (j == std::string::npos) {
    v.push_back(s.substr(i, s.length()));
  }
}

// Copy substring str[a,b] into substr and returns substr
char *substr(char* substr, char* str, int a, int b) {
  if (b < a) {
    int i = b;
    b = a;
    a = i;
  }
  if (b - a > 1000) {
    printf("Function substr: >1000 chars to copy. Exiting.\n");
    exit(6);
  }
  char* dest = substr;
  char* source = str + a;
  char* send = str + b;
  while (*source != '\0' && source <= send)
    *(dest++) = *(source++);
  *dest = '\0';
  return substr;
}

// Similar to Perl's tr/abc/ABC/: Replaces all chars in str found in one list with characters from the second list
// Returns the number of replaced characters
int strtr(char* str, const char oldchars[], const char newchars[]) {
  char* ptr;
  const char *plist;
  int ntr = 0;
  for (ptr = str; *ptr != '\0'; ptr++)
    for (plist = oldchars; *plist != '\0'; plist++)
      if (*ptr == *plist) {
        *ptr = newchars[plist - oldchars];
        ntr++;
        break;
      }
  return ntr;
}

// Similar to Perl's tr/abc//d: deletes all chars in str found in the list
// Returns number of removed characters
int strtrd(char* str, const char chars[]) {
  char* ptr0 = str;
  char* ptr1 = str;
  const char *plist;
  while (*ptr1 != '\0') {
    for (plist = chars; *plist != '\0'; plist++)
      if (*ptr1 == *plist)
        break;
    if (*plist == '\0') {
      *ptr0 = *ptr1;
      ptr0++;
    }
    ptr1++;
  }
  *ptr0 = *ptr1;
  return ptr1 - ptr0;
}

// Similar to Perl's tr/a-z//d: deletes all chars in str found in the list
// Returns number of removed characters
int strtrd(char* str, char char1, char char2) {
  char* ptr0 = str;
  char* ptr1 = str;
  while (*ptr1 != '\0') {
    if (*ptr1 >= char1 && *ptr1 <= char2) {
      *ptr0 = *ptr1;
      ptr0++;
    }
    ptr1++;
  }
  *ptr0 = *ptr1;
  return ptr1 - ptr0;
}

// Counts the number of characters in str that are in range between char1 and char2
int strcount(char* str, char char1, char char2) {
  char* ptr = str;
  int count = 0;
  while (*ptr != '\0') {
    if (*ptr >= char1 && *ptr <= char2)
      count++;
    ptr++;
  }
  return count;
}

// transforms str into an all uppercase string
char* uprstr(char* str) {
  char* s = str;
  while (*s != '\0') {
    if (*s >= 'a' && *s <= 'z')
      *s += 'A' - 'a';
    s++;
  }
  return (str);
}

// transforms str into an all uppercase string
char* lwrstr(char* str) {
  char* s = str;
  while (*s != '\0') {
    if (*s >= 'A' && *s <= 'Z')
      *s += 'a' - 'A';
    s++;
  }
  return (str);
}

// Returns leftmost integer in ptr and sets the pointer to first char after
// the integer. If no integer is found, returns INT_MIN and sets pt to NULL

int strint(char*& ptr) {
  const char* ptr0 = ptr;
  if (!ptr)
    return INT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9'))
    ptr++;
  if (*ptr == '\0') {
    ptr = 0;
    return INT_MIN;
  }
  int i;
  if (ptr > ptr0 && *(ptr - 1) == '-')
    i = -atoi(ptr);
  else
    i = atoi(ptr);
  while (*ptr >= '0' && *ptr <= '9')
    ptr++;
  return i;
}

//
//int strint(const char* ptr) {
//  const char* ptr0 = ptr;
//  if (!ptr)
//    return INT_MIN;
//  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9'))
//    ptr++;
//  if (*ptr == '\0') {
//    ptr = 0;
//    return INT_MIN;
//  }
//  int i;
//  if (ptr > ptr0 && *(ptr - 1) == '-')
//    i = -atoi(ptr);
//  else
//    i = atoi(ptr);
//  while (*ptr >= '0' && *ptr <= '9')
//    ptr++;
//  return i;
//}

// Same as strint, but interpretes '*' as default
int strinta(char*& ptr, int deflt) {
  int i;
  if (!ptr)
    return INT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9') && *ptr != '*')
    ptr++;
  if (*ptr == '\0') {
    ptr = 0;
    return INT_MIN;
  }
  if (*ptr == '*') {
    ptr++;
    return deflt;
  }
  if (*(ptr - 1) == '-')
    i = atoi(ptr - 1);
  else
    i = atoi(ptr);
  while (*ptr >= '0' && *ptr <= '9')
    ptr++;
  return i;
}

// Returns leftmost float in ptr and sets the pointer to first char after
// the float. If no float is found, returns FLT_MIN and sets pt to NULL
float strflt(char*& ptr) {
  float i;
  char* ptr0 = ptr;
  if (!ptr)
    return FLT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9'))
    ptr++;
  if (*ptr == '\0') {
    ptr = 0;
    return FLT_MIN;
  }
  if (ptr > ptr0 && *(ptr - 1) == '-')
    i = -atof(ptr);
  else
    i = atof(ptr);
  while ((*ptr >= '0' && *ptr <= '9') || *ptr == '.')
    ptr++;
  return i;
}

// Same as strint, but interpretes '*' as default
float strflta(char*& ptr, float deflt) {
  float i;
  if (!ptr)
    return FLT_MIN;
  while (*ptr != '\0' && !(*ptr >= '0' && *ptr <= '9') && *ptr != '*')
    ptr++;
  if (*ptr == '\0') {
    ptr = 0;
    return FLT_MIN;
  }
  if (*ptr == '*') {
    ptr++;
    return deflt;
  }
  if (*(ptr - 1) == '-')
    i = -atof(ptr);
  else
    i = atof(ptr);
  while ((*ptr >= '0' && *ptr <= '9') || *ptr == '.')
    ptr++;
  return i;
}

// QSort sorting routine. time complexity of O(N ln(N)) on average
// Sorts the index array k between elements i='left' and i='right' in such a way that afterwards
// v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
void QSortInt(int v[], int k[], int left, int right, int up) {
  int i;
  int last;   // last element to have been swapped

  if (left >= right)
    return;        // do nothing if less then 2 elements to sort
  // Put pivot element in the middle of the sort range to the side (to position 'left') ...
  swapi(k, left, (left + right) / 2);
  last = left;
  // ... and swap all elements i SMALLER than the pivot
  // with an element that is LARGER than the pivot (element last+1):
  if (up == 1) {
    for (i = left + 1; i <= right; i++)
      if (v[k[i]] < v[k[left]])
        swapi(k, ++last, i);
  }
  else
    for (i = left + 1; i <= right; i++)
      if (v[k[i]] > v[k[left]])
        swapi(k, ++last, i);

  // Put the pivot to the right of the elements which are SMALLER, left to elements which are LARGER
  swapi(k, left, last);

  // Sort the elements left from the pivot and right from the pivot
  QSortInt(v, k, left, last - 1, up);
  QSortInt(v, k, last + 1, right, up);
}

// QSort sorting routine. time complexity of O(N ln(N)) on average
// Sorts the index array k between elements i='left' and i='right' in such a way that afterwards
// v[k[i]] is sorted downwards (up=-1) or upwards (up=+1)
void QSortFloat(float v[], int k[], int left, int right, int up) {
  int i;
  int last;   // last element to have been swapped
  void swapi(int k[], int i, int j);

  if (left >= right)
    return;        // do nothing if less then 2 elements to sort
  // Put pivot element in the middle of the sort range to the side (to position 'left') ...
  swapi(k, left, (left + right) / 2);
  last = left;
  // ... and swap all elements i SMALLER than the pivot
  // with an element that is LARGER than the pivot (element last+1):
  if (up == 1) {
    for (i = left + 1; i <= right; i++)
      if (v[k[i]] < v[k[left]])
        swapi(k, ++last, i);
  }
  else
    for (i = left + 1; i <= right; i++)
      if (v[k[i]] > v[k[left]])
        swapi(k, ++last, i);

  // Put the pivot to the right of the elements which are SMALLER, left to elements which are LARGER
  swapi(k, left, last);

  // Sort the elements left from the pivot and right from the pivot
  QSortFloat(v, k, left, last - 1, up);
  QSortFloat(v, k, last + 1, right, up);
}

void readU16(char** ptr, uint16_t &result) {
  unsigned char array[2];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[0] | (array[1] << 8);
}

void readU32(char** ptr, uint32_t &result) {
  unsigned char array[4];

  array[0] = (unsigned char) (**ptr);
  (*ptr)++;
  array[1] = (unsigned char) (**ptr);
  (*ptr)++;
  array[2] = (unsigned char) (**ptr);
  (*ptr)++;
  array[3] = (unsigned char) (**ptr);
  (*ptr)++;

  result = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);
}
