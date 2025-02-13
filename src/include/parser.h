#ifndef _PARSER_H
#define _PARSER_H

#include <string.h>
#include "matrix.h"

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <typeinfo>
#include <clocale>

inline void force_locale() {
    setlocale(LC_NUMERIC, "C");
}

class cmdline_parser {
    int argc, i;
    char **argv;
public:
    cmdline_parser() { force_locale(); }
    void open(int argc_in, char *argv_in[]);
    int get(char *&arg, char *&val);
    void ack(char *arg, char *val);
    void close();
};

class file_parser {
    FILE *fp;
    int iline;
    char line[1024];
public:
    file_parser() { force_locale(); }
    int open(const char *file);
    int get(char *&arg,char *&val);
    void close();
};

#endif

