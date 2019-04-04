#ifndef __FMEMOPEN_H__
#define __FMEMOPEN_H__

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* __FMEMOPEN_H__ */

