#include "debug.h"

#include <cfenv>
#include <csignal>
#include <cstdlib>
#include <iostream>

void sigfpe_handler(int) {
    ester_err("SIGFPE");
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
