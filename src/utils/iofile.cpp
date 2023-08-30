#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include <string.h>
#include <stdlib.h>

bool isHDF5Name(const char *fileName) {
    if (fileName == NULL) {
        return false;
    }

    const char *ext;
    if((ext = strstr(fileName, ".")) == NULL) {
        return false;
    }
    if (strcasecmp(ext, ".hdf5") == 0 || strcasecmp(ext, ".h5") == 0) {
        return true;
    }
    return false;
}
