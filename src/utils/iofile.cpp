#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include <string.h>
#include <stdlib.h>

bool isHDF5Name(const char *fileName) {
    char *name = strdup(fileName);
    char *ext = name;
    char *ptr = name;
    if (fileName == NULL) {
        free(name);
        return false;
    }
    while ((ptr = strstr(ext, ".")) != NULL) {
        ext = ptr + 1;
    }
    if (strcasecmp(ext, "hdf5") == 0 || strcasecmp(ext, "h5") == 0) {
        free(name);
        return true;
    }
    free(name);
    return false;
}
