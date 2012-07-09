#ifndef _PARSER_H
#define _PARSER_H

#include<stdio.h>

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

#endif

