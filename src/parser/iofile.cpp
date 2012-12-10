#include"parser.h"
#include<string.h>
#include<stdlib.h>


int OUTFILE::open(const char *name,char mode_set) {

	char m[3];

	if(mode_set=='t') mode='t';
	else mode='b';
	m[0]='w';m[1]=mode;m[2]='\0';
	fp=fopen(name,m);
	
	if(fp==NULL) return 0;
	
	if(mode=='t') {
		char tag[33]="ESTERdata_t";
		fprintf(fp,"%s\n",tag);
	} else {
		char tag[33]="ESTERdata_b";
		fwrite(tag,sizeof(char),32,fp);
	}
	
	return 1;

}

void OUTFILE::close() {

	fclose(fp);

}

void OUTFILE::write_tag(const char *tag,unsigned long n) {

	if(mode=='t') fprintf(fp,"%s\n%lu\n",tag,n);
	else {
		int l=strlen(tag);
		fwrite(&l,sizeof(int),1,fp);
		fwrite(tag,sizeof(char),l,fp);
		fwrite(&n,sizeof(unsigned long),1,fp);
	}

}

void OUTFILE::write(const char *tag,const matrix *a) {

	if(mode=='t') {
		write_tag(tag,a->ncols()+1);
		fprintf(fp,"%d %d\n",a->nrows(),a->ncols());
	} else {
		write_tag(tag,a->ncols()*a->nrows()*sizeof(double)+2*sizeof(int));
		int n=a->nrows();
		fwrite(&n,sizeof(int),1,fp);
		n=a->ncols();
		fwrite(&n,sizeof(int),1,fp);
	}
	
	a->write(fp,mode);
	
}

void OUTFILE::write(const char *tag,void *x,unsigned long n,size_t size) {

	write_tag(tag,n*size);
	fwrite(x,size,n,fp);

}

int INFILE::open(const char *name,char mode_set) {

	char m[3];

	if(mode_set=='t') mode='t';
	else mode='b';
	m[0]='r';m[1]=mode;m[2]='\0';
	fp=fopen(name,m);
	
	if(fp==NULL) return 0;
	
	char tag[33];
	tag[0]='\0';
	
	if(mode=='t') {
		getline(tag,33);
		if(strcmp(tag,"ESTERdata_t")) {
			fclose(fp);
			return 0;
		}
	} else {
		fread(tag,sizeof(char),32,fp);
		tag[32]='\0';
		if(strcmp(tag,"ESTERdata_b")) {
			fclose(fp);
			return 0;
		}
	}
	
	return 1;

}

void INFILE::close() {

	fclose(fp);

}

char *INFILE::getline(char *line,int n) {

	if(!fgets(line,n,fp)) return NULL;
	if(line[strlen(line)-1]=='\n') {
		line[strlen(line)-1]='\0';
		return line;
	}
	char buf[1024];
	while(fgets(buf,1024,fp)) {
		if(buf[strlen(buf)-1]=='\n') break;
	}

	return line;

}

unsigned long INFILE::seek(const char *tag) {

	unsigned long n=0;
	char *tag2;

	if(mode=='t') {
		bool found=false;
		tag2=new char[512];
		rewind(fp);
		getline(tag2,512);
		while(!found && getline(tag2,512)) {
			if(!strcmp(tag,tag2)) {
				found=true;
				getline(tag2,512);
				n=strtoul(tag2,NULL,0);
			} else {
				getline(tag2,512);
				n=strtoul(tag2,NULL,0);
				for(unsigned long i=0;i<n;i++) getline(tag2,512);
				n=0;
			}
		}
	
	} else {
		bool found=false;
		int tag_len;
		fseek(fp,32,SEEK_SET);
		while(!found && fread(&tag_len,sizeof(int),1,fp)) {
			tag2=new char[tag_len+1];
			fread(tag2,sizeof(char),tag_len,fp);
			tag2[tag_len]='\0';
			fread(&n,sizeof(unsigned long),1,fp);
			if(!strcmp(tag,tag2)) found=true;
			else {fseek(fp,n,SEEK_CUR);n=0;}
			delete [] tag2;
		}
	}
	return n;
}

unsigned long INFILE::len(const char *tag) {

	return seek(tag);

}

int INFILE::read(const char *tag,matrix *a) {

	if(!seek(tag)) return 0;
	int nr,nc;
	if(mode=='t') {
		char line[512];
		getline(line,512);
		sscanf(line,"%d %d",&nr,&nc);
	} else {
		fread(&nr,sizeof(int),1,fp);
		fread(&nc,sizeof(int),1,fp);
	}

	if(!a->read(nr,nc,fp,mode)) return 0;
	else return 1;

}

int INFILE::read(const char *tag,void *x) {

	unsigned long n=seek(tag);

	if(!n) return 0;
	
	if(n!=fread(x,1,n,fp)) return 0;
	
	return 1;
	
}


