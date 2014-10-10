extern "C" {
#include <stdio.h>
}

#ifdef DEBUG

#define DEBUG_FUNCNAME fprintf(stderr, "%s\n", __PRETTY_FUNCTION__);

#else 

#define DEBUG_FUNCNAME

#endif
