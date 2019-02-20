#ifndef __FMEMOPEN_H__
#define __FMEMOPEN_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct fmem {
    size_t pos;
    size_t size;
    char *buffer;
};
typedef struct fmem fmem_t;

FILE *fmemopen(void *, size_t, const char *);

#endif /* __FMEMOPEN_H__ */
