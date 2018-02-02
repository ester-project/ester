#include "matplotlib.h"

#ifndef USE_DEPRECATED_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif
#include <numpy/arrayobject.h>

#include <sstream>

std::map<std::string, PyObject *> plt::py;
std::vector<std::string> plt::functions = {
    "plot",
    "semilogx",
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
    "figure",
    "axis",
    "pause"};
bool plt::noplot = false;

static PyObject *import_module(const std::string& name) {
    PyObject *module_name = PyUnicode_FromString(name.c_str());
    if (!module_name) {
        ester_warn("error initializing %s", name.c_str());
        PyErr_PrintEx(0);
        return NULL;
    }
    PyObject *module = PyImport_Import(module_name);
    if (!module) {
        ester_warn("import %s failed", name.c_str());
        PyErr_PrintEx(0);
        return NULL;
    }
    Py_DECREF(module_name);
    return module;
}

#if PY_MAJOR_VERSION >= 3
void *
#else
void
#endif
import_array_wrapper() {
    import_array();
#if PY_MAJOR_VERSION >= 3
    return NULL;
#endif
}

void plt::init(bool noplot) {
    plt::noplot = noplot;

    if (noplot) return;

    static bool init = false;
    if (init == false) {

        bool close = false;

#if PY_MAJOR_VERSION >= 3
        Py_SetProgramName((wchar_t *) std::wstring(L"ester").c_str());
#else
        Py_SetProgramName((char *) std::string("ester").c_str());
#endif
        Py_Initialize();

        import_array_wrapper();

        PyObject *res;

        // try to use system's python libs first
        std::ostringstream str_stream;
        str_stream << "import sys; sys.path.insert(0, '/usr/lib/python"
            << PYTHON_VERSION
            << "/dist-packages')\n";
        PyRun_SimpleString(str_stream.str().c_str());

        PyObject *matplotlib = import_module("matplotlib");
        if (!matplotlib) {
            ester_warn("import matplotlib failed");
            PyErr_PrintEx(0);
            plt::noplot = true;
            return;
        }

#if defined(__linux__)
#if USE_GTK
        res = PyObject_CallMethod(matplotlib,
                const_cast<char *>(std::string("use").c_str()),
                const_cast<char *>(std::string("s").c_str()),
                "GTKAgg");
        if (res) Py_DECREF(res);
        else {
            close = true;
        }
#endif
#endif

        PyObject *pyplot = import_module("matplotlib.pyplot");
        if (!pyplot) {
            ester_warn("import matplotlib.pyplot failed");
            PyErr_PrintEx(0);
            plt::noplot = true;
            return;
        }

        for (auto function: plt::functions) {
            py[function] = PyObject_GetAttrString(pyplot, function.c_str());
            if (!py[function]) {
                ester_warn("Failed loading python function `%s'", function.c_str());
                PyErr_PrintEx(0);
                plt::noplot = true;
                return;
            }
            if (!PyFunction_Check(py[function])) {
                ester_warn("`%s' is not a function...", function.c_str());
                PyErr_PrintEx(0);
                plt::noplot = true;
                return;
            }
        }

        // PyObject *gcf = PyObject_GetAttrString(pyplot, std::string("gcf").c_str());
        // PyObject *args = PyTuple_New(0);
        // res = PyObject_CallObject(gcf, args);
        // Py_DECREF(args);
        // if (res) {
        //     PyObject *set_size = PyObject_GetAttrString(res, std::string("set_size_inches").c_str());
        //     if (set_size) {
        //         PyObject *args = PyTuple_New(2);
        //         PyTuple_SetItem(args, 0, PyFloat_FromDouble(12.8));
        //         PyTuple_SetItem(args, 1, PyFloat_FromDouble(4.8));
        //         PyObject *r = PyObject_CallObject(set_size, args);
        //         Py_DECREF(args);
        //         if (r) Py_DECREF(r);
        //         Py_DECREF(res);
        //         Py_DECREF(set_size);
        //     }
        // }

        atexit(Py_Finalize);
        if (close == false) {
            atexit(plt::block);
        }
        init = true;
    }
}

PyObject *matrix_to_py(const matrix& m) {
    npy_intp dims[2];
    dims[0] = (size_t) m.nrows();
    dims[1] = (size_t) m.ncols();


#ifndef USE_DEPRECATED_NUMPY
#define NUMPY_ARRAY NPY_ARRAY_FARRAY_RO
#else
#define NUMPY_ARRAY NPY_FARRAY_RO
#endif
    return PyArray_New(&PyArray_Type, 2, dims,
            NPY_DOUBLE, NULL, m.data(),
            sizeof(double), NUMPY_ARRAY, NULL);
#undef NUMPY_ARRAY
}

void plt::plot(const matrix& x, std::string label, std::string style) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, PyUnicode_FromString(style.c_str()));

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("plot", args, kwargs);
}

void plt::plot(const matrix& x, const matrix& y, std::string label, std::string style) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(3);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);
    PyTuple_SetItem(args, 2, PyUnicode_FromString(style.c_str()));

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("plot", args, kwargs);
}

void plt::semilogx(const matrix& x, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, pyx);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("semilogx", args, kwargs);
}

void plt::semilogx(const matrix& x, const matrix& y, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("semilogx", args, kwargs);
}

void plt::semilogy(const matrix& x, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, pyx);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("semilogy", args, kwargs);
}

void plt::semilogy(const matrix& x, const matrix& y, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("semilogy", args, kwargs);
}

void plt::loglog(const matrix& x, const matrix& y, std::string label) {
    if (noplot) return;

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("loglog", args, kwargs);
}

void plt::legend(std::string loc) {
    if (noplot) return;

    if (loc == "") {
        call("legend");
    }
    else {
        PyObject *args = PyTuple_New(0);

        PyObject *kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "loc", PyUnicode_FromString(loc.c_str()));
        call("legend", args, kwargs);
    }
}

void plt::colorbar() {
    if (noplot) return;

    call("colorbar");
}

void plt::title(const std::string& title) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyUnicode_FromString(title.c_str()));
    call("title", args);
}

void plt::clf() {
    if (noplot) return;

    call("clf");
}

void plt::ion() {
    if (noplot) return;

    call("ion");
}

void plt::ioff() {
    if (noplot) return;

    call("ioff");
}

void plt::close() {
    if (noplot) return;

    call("close");
}

void plt::show(bool block) {
    if (noplot) return;

    PyObject *args = PyTuple_New(0);
    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "block", PyBool_FromLong((long) block));
    call("show", args, kwargs);
}

void plt::block() {
    if (noplot) return;

    show(true);
}

void plt::pause(double t) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(t));
    call("pause", args);
}

void plt::subplot(int subplot, bool clear_axis) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyLong_FromLong((long) subplot));
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
    else {
        ester_warn("matplotlib.pyplot.pause failed");
        PyErr_PrintEx(0);
    }
}

void plt::draw() {
    if (noplot) return;

    call("draw");
}

void plt::axvline(double x) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));
    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "color", PyUnicode_FromString("gray"));
    PyDict_SetItemString(kwargs, "linestyle", PyUnicode_FromString("--"));
    PyDict_SetItemString(kwargs, "alpha", PyFloat_FromDouble(0.5));
    call("axvline", args, kwargs);
}

void plt::text(double x, double y, std::string text) {
    if (noplot) return;

    PyObject *args = PyTuple_New(3);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));
    PyTuple_SetItem(args, 1, PyFloat_FromDouble(y));
    PyTuple_SetItem(args, 2, PyUnicode_FromString(text.c_str()));
    call("text", args);
}

void plt::savefig(const std::string& filename) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyUnicode_FromString(filename.c_str()));
    call("savefig", args);
}

void plt::call(const std::string& name, PyObject *args, PyObject *kwargs) {
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

void plt::pcolormesh(const matrix& x, const matrix& y, const matrix& c) {
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

void plt::figure(const int& id, int width, int height) {
    if (noplot) return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyInt_FromLong(id));

    if (width != -1 || height != -1) {
        PyObject *kwargs = PyDict_New();
        PyObject *size = PyTuple_New(2);
        PyTuple_SetItem(size, 0, PyInt_FromLong(width));
        PyTuple_SetItem(size, 1, PyInt_FromLong(height));
        PyDict_SetItemString(kwargs, "figsize", size);
        call("figure", args, kwargs);
    }
    else {
        call("figure", args);
    }
}

void plt::axis(const double& x0, const double& x1, const double& y0, const double& y1) {
    PyObject *args = PyTuple_New(1);
    PyObject *list = PyList_New(4);

    PyList_SetItem(list, 0, PyFloat_FromDouble(x0));
    PyList_SetItem(list, 1, PyFloat_FromDouble(x1));
    PyList_SetItem(list, 2, PyFloat_FromDouble(y0));
    PyList_SetItem(list, 3, PyFloat_FromDouble(y1));

    PyTuple_SetItem(args, 0, list);
    call("axis", args);
}

void plt::axis(const std::string& a) {
    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyUnicode_FromString(a.c_str()));
    call("axis", args);
}
