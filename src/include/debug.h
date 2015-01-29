extern "C" {
#include <stdio.h>
}

#ifdef DEBUG

#define DEBUG_FUNCNAME do {\
    fprintf(stderr, "%s\n", __PRETTY_FUNCTION__); \
} while(0)

#else 

#define DEBUG_FUNCNAME do {\
} while(0)

#endif
