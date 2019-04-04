// To keep compatibility with configure
#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "debug.h"
#include "utils.h"
#include "ester.h"
#include <stdlib.h>

int main(int argc,char *argv[]) {

    star2d A2d;
    if(A2d.read(argv[1]) == 0) {
        A2d.dump_info();
        return 0;
    }

    star1d A1d;
    if(A1d.read(argv[1]) == 0) {
        A1d.dump_info();
        return 0;
    }

    ester_err("Unknown file format");
    return 1;
}
