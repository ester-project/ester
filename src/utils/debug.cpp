#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "debug.h"
#include "star.h"
#include "matplotlib.h"

#include <cfenv>
#include <csignal>
#include <cstdlib>
#include <iostream>


void sigfpe_handler(int) {
    ester_critical("SIGFPE");
}

void enable_sigfpe() {
    // FE_DIVBYZERO | FE_INEXACT | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW 
    // FE_ALL_EXCEPT
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    signal(SIGFPE, sigfpe_handler);
}

void disable_sigfpe() {
    feclearexcept(FE_ALL_EXCEPT);
}

repl::repl(const std::string& name) {
    std::cout << name << "\n";
    print();
}

void repl::print() {
    std::flush(std::cout);
    std::cout << "> ";
    std::flush(std::cout);
}

debugger::debugger(int argc, char * argv[], star1d& s)
    : repl("ESTER Debugger:"), star(s) {
        plt::close();
    }

std::string repl::read() {
    char line[256];
    std::cin.getline(line, 256);
    return std::string(line);
}

std::vector<std::string> repl::get_args(const std::string& cmd) {
    std::vector<std::string> args;
    std::string arg{""};

    for (auto c: cmd) {
        if (c == ' ' || c == '\0' || c == '\n') {
            args.push_back(arg);
            arg = "";
        }
        else {
            arg += c;
        }
    }

    if (arg != "") args.push_back(arg);
    return args;
}

void debugger::help() {
    std::cout << "Commands availiable:\n";
    std::cout << "  help     # prints this help\n";
    std::cout << "  plot var # plots varaible `var'\n";
    std::cout << "  show     # shows plots\n";
    std::cout << "  save     # saves model to file \"debug.h5\"\n";
    std::cout << "  exit     # exits debugger\n";
}

int debugger::eval(const std::string& _cmd) {
    auto args = get_args(_cmd);
    std::string cmd = args[0];
    if (cmd == "clear") {
        plt::clf();
    }
    else if (cmd == "help") {
        help();
    }
    else if (cmd == "show") {
        plt::show(true);
    }
    else if (cmd == "plot") {
        if (args.size() < 2) {
            LOGE("Missong argument to command `%s\'\n", args[0].c_str());
            return 1;
        }
        std::string var = args[1];
        matrix m = get_var(var);
        plt::plot(star.r, m, "$" + var + "$");
        plt::legend();
    }
    else if (cmd == "exit" || cmd == "\0") {
        return 0;
    }
    else if (cmd == "save") {
        star.hdf5_write("debug.h5");
    }
    else if (cmd == "spectrum") {
        if (args.size() < 2) {
            LOGE("Missing argument to command `%s\'\n", args[0].c_str());
            return 1;
        }
        std::string var = args[1];
        matrix m = get_var(var);

        plt::semilogy(star.spectrum(m), "spectrum " + var);
        plt::legend();
    }
    else if (cmd == "") {
    }
    else {
        LOGE("Unknown command `%s\'\n", cmd.c_str());
    }
    return 1;
}

matrix debugger::get_var(const std::string& var) {
    if (var == "rho") return star.rho;
    if (var == "T") return star.T;
    if (var == "p") return star.p;
    if (var == "D,rho") return (star.D,star.rho);
    if (var == "D,T") return (star.D,star.T);
    if (var == "D,p") return (star.D,star.p);
    LOGE("Unknown varaible `%s\'\n", var.c_str());
    return zeros(1, 1);
}

int debugger::exec() {
    int r = 1;
    while (r) {
        auto cmd = read();
        r = eval(cmd);
        print();
    }
    return r;
}
