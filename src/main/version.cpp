// To keep compatibility with configure
#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "stdio.h"

int main() {
	printf("ESTER %s\n", PACKAGE_VERSION);
    return 0;
}
