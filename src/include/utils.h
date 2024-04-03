#ifndef ESTER_UTILS_H
#define ESTER_UTILS_H

#include "stack.h"

#include <stdlib.h>
#include <stdio.h>
#include <execinfo.h>
#include <unistd.h>

#ifdef DEBUG
#include "debug.h"
#endif

#define RED "\033[0;31m"
#define GREEN "\033[0;32m"
#define YELLOW "\033[0;33m"
#define RESET "\033[0;0m"

/*
Define utility log functions:
ester_debug|info|warn|err|critical

`ester_debug`:
print to stderr in DEBUG
none else

`ester_info`:
prefixed with green "Info"
print to stdout in both modes

`ester_warn`:
prefixed with yellow "Warning:"
print to stderr
+ indicate line and file in DEBUG

`ester_err`:
same as `ester_warn` with red "Error:"
+ print stack in DEBUG

`ester_critical`:
same as `ester_err` with red "CRITICAL:"
+ exit program or throw eception
*/

#ifdef DEBUG
#define ester_debug(...) do { \
    fprintf(stderr, __VA_ARGS__); \
} while (0)

#define ester_info(...) do { \
    printf(GREEN); \
    printf("Info: "); \
    printf(RESET); \
    printf(__VA_ARGS__); \
    fflush(stdout); \
} while (0)

#define ester_warn(...) do { \
    fprintf(stderr, YELLOW); \
    fprintf(stderr, "Warning at %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, RESET); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
} while(0)

#define ester_err(...) do { \
    fprintf(stderr, RED); \
    fprintf(stderr, "Error at %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, RESET); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    print_stack();
} while(0)

#define ester_critical(...) do { \
    fprintf(stderr, RED); \
    fprintf(stderr, "CRITICAL at %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, RESET); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    print_stack();
    throw runtime_exception(); \
} while(0)

#else // DEBUG
#define ester_debug(...) do {} while (0)

#define ester_info(...) do { \
    printf(GREEN); \
    printf("Info: "); \
    printf(RESET); \
    printf(__VA_ARGS__); \
    fflush(stdout); \
} while (0)

#define ester_warn(...) do { \
    fprintf(stderr, YELLOW); \
    fprintf(stderr, "Warning: "); \
    fprintf(stderr, RESET); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
} while(0)

#define ester_err(...) do { \
    fprintf(stderr, RED); \
    fprintf(stderr, "Error: "); \
    fprintf(stderr, RESET); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
} while(0)

#define ester_critical(...) do { \
    fprintf(stderr, RED); \
    fprintf(stderr, "CRITICAL: "); \
    fprintf(stderr, RESET); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    exit(EXIT_FAILURE); \
} while(0)

#endif // DEBUG

bool isHDF5Name(const char *fileName);

#endif // ESTER_UTILS_H
