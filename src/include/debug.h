#ifndef DEBUG_H
#define DEBUG_H

#include "utils.h"

void sigfpe_handler(int);
void enable_sigfpe();
void disable_sigfpe();

#endif
