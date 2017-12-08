#ifndef MATPLOTLIB_H
#define MATPLOTLIB_H

#include "ester-config.h"
#include "utils.h"
#include "matrix.h"

#include <Python.h>
#include <vector>

class matplotlib {
    private:
        static std::map<std::string, PyObject *> py;
        static std::vector<std::string> functions;
        static void call(const std::string& name,
                PyObject *args = NULL,
                PyObject *kwargs = NULL);
        static void block();
        static bool noplot;

    public:
        static void init(bool noplot = false);
        matplotlib();
        static void subplot(int, bool clear_axis = false);
        static void plot(const matrix&, std::string label = "");
        static void plot(const matrix&, const matrix&, std::string label = "");
        static void semilogx(const matrix&, std::string label = "");
        static void semilogy(const matrix&, std::string label = "");
        static void loglog(const matrix&, const matrix&, std::string label = "");
        static void pcolormesh(const matrix&, const matrix&, const matrix&);
        static void axvline(double);
        static void text(double, double, std::string);
        static void show(bool block = false);
        static void ion();
        static void ioff();
        static void clf();
        static void draw();
        static void legend();
        static void colorbar();
        static void close();
        static void savefig(const std::string&);
        static void title(const std::string&);
        static void pause(double = 1e-4);
};

extern matplotlib plt;

#endif
