#ifndef WITH_CMAKE
#include "config.h"
#endif
#include "stack.h"

#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <cstdint>

inline std::string addr2str(uintptr_t addr) {
    std::ostringstream stream;
    stream << "0x" << std::hex << addr;
    return stream.str();
}

void print_stack() {
    void *stack[STACK_SIZE];
    int size = backtrace(stack, STACK_SIZE);
    std::cerr << std::setfill('-') << std::setw(20) << '['
        << "Error Stack Trace]";
    std::cerr << std::setw(20) << '\n';
    std::cerr << std::setfill(' ');
#if HAVE_ADDR2LINE
    Dl_info *info = new Dl_info;

    // start at index 1 to skip print_stack...
    for (int i=1; i<size; i++) {
        int st = dladdr(stack[i], info);
        if (st != 0 && info->dli_fname != NULL && info->dli_fname[0] != '\0') {
            std::string filename(info->dli_fname);
            uintptr_t addr = uintptr_t(stack[i]);
            uintptr_t offset = uintptr_t(info->dli_fbase);
            addr = addr - offset;
            std::string cmd = "addr2line --functions --demangle=auto --exe ";
            cmd += filename + " " + addr2str(addr);

            FILE *pipe = popen(cmd.c_str(), "r");
            char *buf = NULL;

            if (pipe) {
                size_t linelen = 0;
                int nline = 0;
                while (getline(&buf, &linelen, pipe) > 0) {
                    if (nline%2 == 0) {
                        std::cerr << "#" << std::setw(3) << i << " in " << buf;
                    }
                    else {
                        std::string loc(buf);
                        size_t end = loc.find_first_of(' ');
                        loc = loc.substr(0, end);
                        std::cerr << "     at: "
                            <<  loc;
                        if (end != std::string::npos) {
                            std::cerr << '\n';
                        }
                    }
                    linelen = 0;
                    free (buf);
                    nline++;
                }
            }
            else {
                std::cerr << "popen failed\n";
            }

            pclose(pipe);
        }
    }
    delete info;
#else
    backtrace_symbols_fd(stack, size, STDERR_FILENO);
#endif
}
