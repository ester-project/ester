#ifndef _PARSER_H
#define _PARSER_H

#include "ester-config.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <string.h>
#include "matrix.h"

class cmdline_parser {
	int argc,i;
	char **argv;
public:
	void open(int argc_in,char *argv_in[]);
	int get(char *&arg,char *&val);
	void ack(char *arg,char *val);
	void close();
};

class file_parser {
	FILE *fp;
	int iline;
	char line[1024];
public:
	int open(const char *file);
	int get(char *&arg,char *&val);
	void close();
};

class OUTFILE {
	FILE *fp;
	char mode;
	void write_tag(const char *tag,unsigned long n);
public:
	int open(const char *name, char mode='b');
	void write(const char *tag,const matrix *);
	void write(const char *tag,const matrix_map *);
	void write(const char *tag,matrix *m) {write(tag,(const matrix *) m);};
	void write(const char *tag,matrix_map *m) {write(tag,(const matrix_map *) m);};
	template <class T>
	void write(const char *tag,T *x,unsigned long n=1) {
		size_t size=sizeof(T);
		if(mode=='t') {
			if(typeid(x[0])==typeid(char)) write_tag(tag,1);
			else write_tag(tag,n);
		} else write_tag(tag,n*size);
		if(mode=='t') {
			std::stringstream temp;
			std::string str;
			temp.setf(std::ios::scientific, std::ios::floatfield);
			for(unsigned long i=0;i<n;i++) {
				temp.str("");
				temp.precision(16);
				temp<<x[i]<<std::endl;
				temp>>str;
				fprintf(fp,"%s",str.c_str());
				if(typeid(x[i])!=typeid(char)||i==n-1) fprintf(fp,"\n");
			}
		} else fwrite(x,size,n,fp);
	};
	template <class T>
	void write_fmt(const char *tag,const char *fmt,T *x,unsigned long n=1) {
		write_tag(tag,n);
		for(unsigned long i=0;i<n;i++) {
			fprintf(fp,fmt,x[i]);
			fprintf(fp,"\n");
		}
	}; 
	void close();
};

class INFILE {
	FILE *fp;
	char mode;
	unsigned long seek(const char *tag);
	char *getline(char *,int n);
public:
	int open(const char *name, char mode='b');
	unsigned long len(const char *tag);
	int read(const char *tag,matrix *);
	int read(const char *tag,matrix_map *);
	template <class T>
	int read(const char *tag,T *x) {
		unsigned long n=seek(tag);
		if(!n) return 1;
		if(mode=='t') {
			char str[512];
			if(typeid(x[0])==typeid(char)) {
				getline(str,512);
				strcpy((char *)x,str);
			} else {
				for(unsigned long i=0;i<n;i++) {
					getline(str,512);
					std::stringstream temp(str);
					temp>>x[i];
				}
			}
		}
		else if(n!=fread(x,1,n,fp)) return 1;
		return 0;
	};
	template <class T>
	int read_fmt(const char *tag,const char *fmt,T *x) {
		unsigned long n=seek(tag);
		if(!n) return 0;
		for(unsigned long i=0;i<n;i++) 
			if(!fscanf(fp,fmt,&x[i])) return 0;
		return 1;
	}; 
	void close();
};




#endif

