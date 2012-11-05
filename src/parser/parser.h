#ifndef _PARSER_H
#define _PARSER_H

#include<stdio.h>
#include"matrix.h"

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
	void write(const char *tag,void *,unsigned long n,size_t size);
public:
	int open(const char *name, char mode='b');
	void write(const char *tag,const matrix *);
	void write(const char *tag,matrix *m) {write(tag,(const matrix *) m);};
	template <class T>
	void write(const char *tag,T *x,unsigned long n=1) {write(tag,(void *) x,n,sizeof(T));};
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
	int read(const char *tag,void *);
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

