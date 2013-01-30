#include"ester.h"
#include<stdlib.h>

int main(int argc,char *argv[]) {

	star2d *A;
	
	A=new star2d();
	if(A->read(argv[1])) {
		A->dump_info();
		delete A;
		return 0;
	}
	delete A;
		
	A=new star1d();
	if(A->read(argv[1])) {
		A->dump_info();
		delete A;
		return 0;
	}
	delete A;
	
	printf("Unknown file format\n");
	exit(1);
	
}
