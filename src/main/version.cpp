#include "version.h"
#include"stdio.h"

int main() {

	printf("Version %d.%d (rev %d",ESTER_VERSION_MAJOR,ESTER_VERSION_MINOR,ESTER_REVISION);
	if(ESTER_VERSION_SVN) printf(".svn");
	printf(")\n");

}
