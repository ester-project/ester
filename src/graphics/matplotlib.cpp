#include "matplotlib.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>


std::map<std::string, PyObject *> matplotlib::py;
std::vector<std::string> matplotlib::functions = {
    "plot",
    "semilogy",
    "semilogy",
    "loglog",
    "subplot",
    "text",
    "show",
    "legend",
    "colorbar",
    "title",
    "clf",
    "ion",
    "ioff",
    "draw",
    "axvline",
    "savefig",
    "close",
    "pcolormesh",
    "pause"};
bool matplotlib::noplot = false;

matplotlib plt;

void matplotlib::init(bool noplot) {
    matplotlib::noplot = noplot;
    if (noplot) return;

    static bool init = false;
    if (init == false) {

        // Should the plot window be close at exit (bug with TkAgg backend)
        bool close = false;

        Py_SetProgramName((char *) std::string("ester").c_str());
        Py_Initialize();

        import_array();

        PyObject *res;
        PyObject *matplotlib_name = PyString_FromString("matplotlib");
        if (!matplotlib_name) {
            ester_warn("error initializing matplotlib");
            matplotlib::noplot = true;
            return;
        }
        PyObject *matplotlib = PyImport_Import(matplotlib_name);
        if (!matplotlib) {
            ester_warn("import matplotlib failed\n");
            matplotlib::noplot = true;
            return;
        }
        Py_DECREF(matplotlib_name);

        res = PyObject_CallMethod(matplotlib,
                const_cast<char *>(std::string("use").c_str()),
                const_cast<char *>(std::string("s").c_str()),
                "GTKAgg");
        if (res) Py_DECREF(res);
        else {
            close = true;
        }

        PyObject *pyplot_name = PyString_FromString("matplotlib.pyplot");
        if (!pyplot_name) {
            ester_warn("error initializing matplotlib.pyplot");
            matplotlib::noplot = true;
            return;
        }
        PyObject *pyplot = PyImport_Import(pyplot_name);
        if (!pyplot) {
            ester_warn("import matplotlib.pyplot failed\n");
            matplotlib::noplot = true;
            return;
        }
        Py_DECREF(pyplot_name);

        for (auto function: matplotlib::functions) {
            py[function] = PyObject_GetAttrString(pyplot, function.c_str());
            if (!py[function]) {
                ester_warn("Failed loading python function `%s'", function.c_str());
                matplotlib::noplot = true;
                return;
            }
            if (!PyFunction_Check(py[function])) {
                ester_warn("`%s' is not a function...", function.c_str());
                matplotlib::noplot = true;
                return;
            }
        }

        PyObject *gcf = PyObject_GetAttrString(pyplot, std::string("gcf").c_str());
        PyObject *args = PyTuple_New(0);
        res = PyObject_CallObject(gcf, args);
        Py_DECREF(args);
        if (res) {
            PyObject *set_size = PyObject_GetAttrString(res, std::string("set_size_inches").c_str());
            if (set_size) {
                PyObject *args = PyTuple_New(2);
                PyTuple_SetItem(args, 0, PyFloat_FromDouble(12.8));
                PyTuple_SetItem(args, 1, PyFloat_FromDouble(4.8));
                PyObject *r = PyObject_CallObject(set_size, args);
                Py_DECREF(args);
                if (r) Py_DECREF(r);
                Py_DECREF(res);
                Py_DECREF(set_size);
            }
        }

        atexit(Py_Finalize);
        if (close == false) {
            atexit(matplotlib::block);
        }
        init = true;
    }
}

matplotlib::matplotlib() {
    matplotlib::init();
}

PyObject *matrix_to_py(const matrix& m) {
    npy_intp dims[2];
    dims[0] = (size_t) m.nrows();
    dims[1] = (size_t) m.ncols();


    return PyArray_New(&PyArray_Type, 2, dims,
            NPY_DOUBLE, NULL, m.data(),
            sizeof(double), NPY_ARRAY_FARRAY_RO, NULL);

}

void matplotlib::plot(const matrix& x, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, pyx);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyString_FromString(label.c_str()));

    call("plot", args, kwargs);
}

void matplotlib::plot(const matrix& x, const matrix& y, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyString_FromString(label.c_str()));

    call("plot", args, kwargs);
}

void matplotlib::semilogx(const matrix& x, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, pyx);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyString_FromString(label.c_str()));

    call("semilogx", args, kwargs);
}

void matplotlib::semilogy(const matrix& x, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, pyx);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyString_FromString(label.c_str()));

    call("semilogy", args, kwargs);
}

void matplotlib::legend() {
    if (noplot) return;

    call("legend");
}

void matplotlib::colorbar() {
    if (noplot) return;

    call("colorbar");
}

void matplotlib::title(const std::string& title) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyString_FromString(title.c_str()));
    call("title", args);
}

void matplotlib::clf() {
    if (noplot) return;

    call("clf");
}

void matplotlib::ion() {
    if (noplot) return;

    call("ion");
}

void matplotlib::ioff() {
    if (noplot) return;

    call("ioff");
}

void matplotlib::close() {
    if (noplot) return;

    call("close");
}

void matplotlib::show(bool block) {
    if (noplot) return;

    PyObject *args = PyTuple_New(0);
    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "block", PyBool_FromLong((long) block));
    call("show", args, kwargs);
}

void matplotlib::block() {
    if (noplot) return;

    show(true);
}

void matplotlib::pause(double t) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(t));
    call("pause", args);
}

void matplotlib::subplot(int subplot, bool clear_axis) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyInt_FromLong((long) subplot));
    PyObject *res = PyObject_CallObject(py["subplot"], args);
    Py_DECREF(args);
    if(res) {
        if (clear_axis) {
            PyObject *set_axis_off = PyObject_GetAttrString(res, "set_axis_off");
            PyObject *args = PyTuple_New(0);
            PyObject *r = PyObject_CallObject(set_axis_off, args);
            if (r) Py_DECREF(r);
            Py_DECREF(args);
        }
        Py_DECREF(res);
    }
    else ester_warn("matplotlib.pyplot.pause failed");
}

void matplotlib::draw() {
    if (noplot) return;

    call("draw");
}

void matplotlib::axvline(double x) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));
    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "color", PyString_FromString("gray"));
    PyDict_SetItemString(kwargs, "linestyle", PyString_FromString("--"));
    PyDict_SetItemString(kwargs, "alpha", PyFloat_FromDouble(0.5));
    call("axvline", args, kwargs);
}

void matplotlib::text(double x, double y, std::string text) {
    if (noplot) return;

    PyObject *args = PyTuple_New(3);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));
    PyTuple_SetItem(args, 1, PyFloat_FromDouble(y));
    PyTuple_SetItem(args, 2, PyString_FromString(text.c_str()));
    call("text", args);
}

void matplotlib::savefig(const std::string& filename) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyString_FromString(filename.c_str()));
    call("savefig", args);
}

void matplotlib::call(const std::string& name, PyObject *args, PyObject *kwargs) {
    if (noplot) return;

    PyObject *res;

    if (args == NULL)
        args = PyTuple_New(0);

    if (kwargs == NULL)
        res = PyObject_CallObject(py[name], args);
    else
        res = PyObject_Call(py[name], args, kwargs);

    if (args) Py_DECREF(args);
    if (kwargs) Py_DECREF(kwargs);
    if(res) Py_DECREF(res);
    else ester_warn("call to matplotlib.pyplot.%s failed", name.c_str());
}

void matplotlib::pcolormesh(const matrix& x, const matrix& y, const matrix& c) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);
    PyObject *pyc = matrix_to_py(c);

    PyObject *args = PyTuple_New(3);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);
    PyTuple_SetItem(args, 2, pyc);

    call("pcolormesh", args);
}
