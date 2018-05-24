#ifndef DEBUG_H
#define DEBUG_H

#include "utils.h"
#include "matrix.h"

#include <string>
#include <iostream>
#include <vector>
#include <exception>


void sigfpe_handler(int);
void enable_sigfpe();
void disable_sigfpe();

class runtime_exception : public std::exception {
};

class star1d;

class repl {
    public:
        repl(const std::string& name);
        std::string read();
        void print();
        virtual int eval(const std::string&) = 0;

        static std::vector<std::string> get_args(const std::string&);
};

class debugger : public repl {
    private:
        star1d& star;
        void help();
        matrix get_var(const std::string& var);

    public:
        debugger(int argc, char *argv[], star1d&);
        virtual int eval(const std::string& cmd);
        int exec();
};

#endif
