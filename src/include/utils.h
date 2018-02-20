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

#ifdef __GNUG__
#define _print_stack()   { \
    print_stack(); \
}
#else
#define _print_stack()   { }
#endif

#ifdef DEBUG
#define ester_err(...) do { \
    LOGE("Error at %s:%d: ", __FILE__, __LINE__); \
    LOGE(__VA_ARGS__); \
    LOGE("\n"); \
    print_stack(); \
    throw runtime_exception(); \
} while(0)
#else
#define ester_err(...) do { \
    LOGE("Error: "); \
    LOGE(__VA_ARGS__); \
    LOGE("\n"); \
    print_stack(); \
    exit(EXIT_FAILURE); \
} while(0)
#endif

#ifdef DEBUG
#define ester_warn(...) do { \
    LOGW("Warning at %s:%d: ", __FILE__, __LINE__); \
    LOGW(__VA_ARGS__); \
    LOGW("\n"); \
} while(0)
#else
#define ester_warn(...) do { \
    LOGW("Warning: "); \
    LOGW(__VA_ARGS__); \
    LOGW("\n"); \
} while(0)
#endif

#ifdef DEBUG
#define ester_debug(...) do { \
    fprintf(stderr, __VA_ARGS__); \
} while (0)
#else
#define ester_debug(...) do {} while (0)
#endif

bool isHDF5Name(const char *fileName);

#define RED "\033[0;31m"
#define GREEN "\033[0;32m"
#define YELLOW "\033[0;33m"
#define RESET "\033[0;0m"

#define LOGI(...)   { printf(GREEN); printf(__VA_ARGS__); printf(RESET); fflush(stdout); }
#define LOGW(...)   { fprintf(stderr, YELLOW); fprintf(stderr, __VA_ARGS__); fprintf(stderr, RESET); }
#define LOGE(...)   { fprintf(stderr, RED); fprintf(stderr, __VA_ARGS__); fprintf(stderr, RESET); }

#endif // ESTER_UTILS_H
