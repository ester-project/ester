#include "ester-config.h"
#include "stdio.h"

int main() {
	printf("ESTER version %s\n", PACKAGE_VERSION);
#ifdef USE_HDF5
    const char *h5supp = "yes";
#else
    const char *h5supp = "no";
#endif
#ifdef USE_PGPLOT
    const char *pgsupp = "yes";
#else
    const char *pgsupp = "no";
#endif
	printf("HDF5 support: %s\n", h5supp);
	printf("PGPLOT support: %s\n", pgsupp);
    return 0;
}
