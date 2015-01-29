#include "ester-config.h"
#include "utils.h"
#include "parser.h"

extern "C" {
#include <string.h>
#include <stdlib.h>
}

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
		fwrite(tag,sizeof(tag[0]),32,fp);
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
		fwrite(&l,sizeof(l),1,fp);
		fwrite(tag,sizeof(tag[0]),l,fp);
		fwrite(&n,sizeof(n),1,fp);
	}

}

void OUTFILE::write(const char *tag,const matrix *a) {

	if(mode=='t') {
		write_tag(tag,a->ncols()+1);
		fprintf(fp,"%d %d\n",a->nrows(),a->ncols());
	} else {
		write_tag(tag, a->ncols()*a->nrows()*sizeof(double)+2*sizeof(int));
		int n=a->nrows();
		fwrite(&n,sizeof(n),1,fp);
		n=a->ncols();
		fwrite(&n,sizeof(n),1,fp);
	}

	a->write(fp,mode);

}

void OUTFILE::write(const char *tag,const matrix_map *a) {

	if(mode=='t') {
		int nrows=1,nitems=0;
		matrix_map::const_iterator it;
		for(it=a->begin();it!=a->end();it++) {
			nrows+=(it->second).ncols()+2;
			nitems++;
		}
		write_tag(tag,nrows);
		fprintf(fp,"%d\n",nitems);
		for(it=a->begin();it!=a->end();it++) {
			fprintf(fp,"%s\n",(it->first).c_str());
			fprintf(fp,"%d %d\n",(it->second).nrows(),(it->second).ncols());
			(it->second).write(fp,mode);
		}
	} else {
		int nbytes=sizeof(int),nitems=0;
		matrix_map::const_iterator it;
		for(it=a->begin();it!=a->end();it++) {
			nitems++;
			nbytes+=(it->second).ncols()*(it->second).nrows()*sizeof(double)
					+3*sizeof(int)+(it->first).length()*sizeof(char);
		}
		write_tag(tag,nbytes);
		fwrite(&nitems,sizeof(int),1,fp);
		for(it=a->begin();it!=a->end();it++) {
			int n;
			n=(it->first).length();
			fwrite(&n,sizeof(int),1,fp);
			fwrite((it->first).c_str(),sizeof(char),n,fp);
			n=(it->second).nrows();
			fwrite(&n,sizeof(int),1,fp);
			n=(it->second).ncols();
			fwrite(&n,sizeof(int),1,fp);
			(it->second).write(fp,mode);
		}
	}

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
		if (fread(tag,sizeof(tag[0]),32,fp) < 32)
            ester_warn("could not read tag");
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
		while(!found && fread(&tag_len,sizeof(tag_len),1,fp)) {
			tag2=new char[tag_len+1];
			if (fread(tag2,sizeof(tag2[0]),tag_len,fp) < (size_t) tag_len)
                ester_warn("could not read tag2");
			tag2[tag_len]='\0';
			if (fread(&n,sizeof(n),1,fp) < 1)
                ester_warn("could not read n");
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

int INFILE::read(const char *tag, matrix *a) {
	if(!seek(tag)) return 1;
	int nr,nc;
	if(mode=='t') {
		char line[512];
		getline(line,512);
		sscanf(line,"%d %d",&nr,&nc);
	} else {
		if (fread(&nr,sizeof(nr),1,fp) < 1)
            ester_warn("could not read nr");
		if (fread(&nc,sizeof(nc),1,fp) < 1)
            ester_warn("could not read nc");
	}

	if(!a->read(nr,nc,fp,mode)) return 1;
	return 0;
}

int INFILE::read(const char *tag,matrix_map *a) {

	a->clear();
	if(!seek(tag)) return 1;
	int nr,nc,nitems;
	matrix m;
	if(mode=='t') {
		char line[512];
		getline(line,512);
		sscanf(line,"%d",&nitems);
		for(int n=0;n<nitems;n++) {
			getline(line,512);
			std::string item(line);
			getline(line,512);
			sscanf(line,"%d %d",&nr,&nc);
			if(!m.read(nr,nc,fp,mode)) return 1;
			(*a)[item]=m;
			getline(line,512);
		}
	} else {
		if (fread(&nitems,sizeof(nitems),1,fp) < 1)
            ester_warn("could not read nitems");
		char item_s[512];
		for(int n=0;n<nitems;n++) {
			int l;
			if (fread(&l,sizeof(l),1,fp) < 1)
                ester_warn("Could not read '%s' value", tag);
			if (fread(item_s,sizeof(item_s[0]),l,fp) < (size_t) l)
                ester_warn("Could not read '%s' value", tag);
			item_s[l]='\0';
			std::string item(item_s);
			if (fread(&nr,sizeof(nr),1,fp) < 1)
                ester_warn("Could not read '%s' value", tag);
			if (fread(&nc,sizeof(nc),1,fp) < 1)
                ester_warn("Could not read '%s' value", tag);
			if(!m.read(nr,nc,fp,mode)) return 1;
			(*a)[item]=m;
		}
	}

	return 0;
}


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
