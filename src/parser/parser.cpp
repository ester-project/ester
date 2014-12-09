#include "config.h"
#include "parser.h"

extern "C" {
#include <string.h>
}

void cmdline_parser::open(int argc_in, char *argv_in[]) {

	argc=argc_in;
	argv=argv_in;
	i=1;
	
}

int cmdline_parser::get(char *&arg,char *&val) {

	if(i>=argc) return 0;
	while(argv[i][0]=='\0') {
		i++;
		if(i==argc) break;
	}
	if(i==argc) return 0;
	if(argv[i][0]!='-') {
		printf("Invalid argument %s\n",argv[i]);
		arg=argv[i]+1;
		val=NULL;
		i++;
		return -1;
	}
	arg=argv[i]+1;
	val=NULL;
	i++;
	if(i<argc)
		if(argv[i][0]!='\0') {
			bool cond;
			char c;
			c=argv[i][0];
			cond=(c!='-');
			c=argv[i][1];
			cond=cond||c=='0'||c=='1'||c=='2'||c=='3'||c=='4'||
					c=='5'||c=='6'||c=='7'||c=='8'||c=='9'||c=='.'; 
			if(cond)
				val=argv[i++];
		}
	return 1;

}

void cmdline_parser::ack(char *arg,char *val) {

	*(arg-1)='\0';
	if(val!=NULL) *val='\0';

}

void cmdline_parser::close(){}


int file_parser::open(const char *file) {

	iline=0;
	if(!(fp=fopen(file,"rt"))) return 0;
	return 1;
}
	
int file_parser::get(char *&arg,char *&val){

	int i;

	if(fgets(line,1024,fp)) {
		iline++;
		line[strlen(line)-1]='\0';
		arg=line;
		while((*arg==' '||*arg=='\t')&&*arg!='\0'&&*arg!='#') arg++;
		if(*arg=='\0'||*arg=='#') return get(arg,val);
		val=arg;
		while(*val!='='&&*val!='\0'&&*val!='#') {
			if(*val==' '||*arg=='\t') *val='\0';
			val++;
		}
		if(*val=='\0'||*val=='#') {
			*val='\0';
			val=NULL;
		}
		else {
			*val='\0';val++;
			while((*val==' '||*val=='\t')&&*val!='\0'&&*val!='#') val++;
			if(*val=='\0'||*val=='#') val=NULL;
			else {
				i=0;
				while(val[i]!='\0'&&val[i]!='#') i++;
				val[i]='\0';
				i=strlen(val)-1;
				while(val[i]==' '||val[i]=='\t') val[i--]='\0';
			}
		}
		return iline;	
	} else return 0;
	
}

void file_parser::close() {
	
	fclose(fp);
	
}








