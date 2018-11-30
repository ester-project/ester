#ifndef WITH_CMAKE
#include "ester-config.h"
#endif
#include "matplotlib.h"

#if ENABLE_PLT

#ifndef USE_DEPRECATED_NUMPY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif
#include <numpy/arrayobject.h>

#include <sstream>
#include <iostream>
#endif

#if USE_VTK
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTable.h>
#include <vtkChartXY.h>
#include <vtkContextActor.h>
#include <vtkContextScene.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkAxis.h>
#include <vtkRendererCollection.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

double plt::viewport[4];
vtkSmartPointer<vtkRenderWindow> plt::renderWindow;
vtkSmartPointer<vtkRenderWindowInteractor> plt::renderWindowInteractor;
vtkSmartPointer<vtkChartXY> plt::chart;
vtkSmartPointer<vtkRenderer> plt::textRenderer;
int plt::ncolor;

double colors[12] = {1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    0, 0, 0,
};
const int maxColor = 4;

#endif


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
    "axhline",
    "savefig",
    "close",
    "pcolormesh",
    "figure",
    "axis",
    "pause"};

bool plt::noplot = false;
PyObject *plt::matplotlib = nullptr;

#if ENABLE_PLT
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
#endif

#if ENABLE_PLT
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
#endif



void plt::init(bool noplot) {
#if ENABLE_PLT
    plt::noplot = noplot;

    if (noplot) return;

    static bool init = false;
    if (init == false) {

#if PY_MAJOR_VERSION >= 3
        Py_SetProgramName((wchar_t *) std::wstring(L"ester").c_str());
#else
        Py_SetProgramName((char *) std::string("ester").c_str());
#endif
        Py_Initialize();

        import_array_wrapper();

        matplotlib = import_module("matplotlib");
        if (!matplotlib) {
            ester_warn("import matplotlib failed");
            PyErr_PrintEx(0);
            plt::noplot = true;
            return;
        }

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

        init = true;
    }
#endif
#ifdef USE_VTK
    if (renderWindow != nullptr || noplot) return;

    viewport[0] = 0.;
    viewport[1] = 1.;
    viewport[2] = 0.;
    viewport[3] = 1.;

    renderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(768, 575);

    renderWindowInteractor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();

    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindow->SetWindowName("ESTER");

    chart = nullptr;
    textRenderer = nullptr;
#endif
}

PyObject *matrix_to_py(const matrix& m) {
#if ENABLE_PLT
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
#else // ENABLE_PLT
    return nullptr;
#endif
}

void plt::plot(const matrix& x, std::string label, std::string style) {
    if (noplot) return;
    if (x.nrows() < 2) return;
#if ENABLE_PLT

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, PyUnicode_FromString(style.c_str()));

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("plot", args, kwargs);
#endif
#ifdef USE_VTK
    for (auto j=0; j<x.ncols(); j++) {
        vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

        vtkSmartPointer<vtkFloatArray> arrX =
            vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("");
        table->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY =
            vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName(label.c_str());
        table->AddColumn(arrY);

        table->SetNumberOfRows(x.nrows());
        for (auto i=0; i<x.nrows(); i++) {
            table->SetValue(i, 0, i);
            table->SetValue(i, 1, x(i, j));
        }


        if (chart == nullptr) {
            vtkSmartPointer<vtkRenderer> renderer =
                vtkSmartPointer<vtkRenderer>::New();

            renderWindow->AddRenderer(renderer);
            renderer->SetViewport(viewport);
            renderer->SetBackground(1, 1, 1);

            chart = vtkSmartPointer<vtkChartXY>::New();

            chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("");
            chart->GetAxis(vtkAxis::LEFT)->SetTitle("");

            vtkSmartPointer<vtkContextScene> chartScene = vtkSmartPointer<vtkContextScene>::New();
            vtkSmartPointer<vtkContextActor> chartActor = vtkSmartPointer<vtkContextActor>::New();

            chartScene->AddItem(chart);
            chartActor->SetScene(chartScene);

            renderer->AddActor(chartActor);
            chartScene->SetRenderer(renderer);
        }

        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetColor(colors[ncolor], colors[ncolor+1], colors[ncolor+2]);
        line->SetInputData(table, 0, 1);

    }


    ncolor = (ncolor + 3)%(3*maxColor);
#endif
}

void plt::plot(const matrix& x, const matrix& y, std::string label, std::string style) {
    assert(x.nrows() == y.nrows());
    assert(x.ncols() == y.ncols());
    if (noplot) return;
    if (x.nrows() < 2) return;
#if ENABLE_PLT

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(3);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);
    PyTuple_SetItem(args, 2, PyUnicode_FromString(style.c_str()));

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("plot", args, kwargs);
#endif

#ifdef USE_VTK
    init();

    for (auto j=0; j<x.ncols(); j++) {
        vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

        vtkSmartPointer<vtkFloatArray> arrX =
            vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("");
        table->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY =
            vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName(label.c_str());
        table->AddColumn(arrY);

        table->SetNumberOfRows(x.nrows());
        for (auto i=0; i<x.nrows(); i++) {
            table->SetValue(i, 0, x(i, j));
            table->SetValue(i, 1, y(i, j));
        }


        if (chart == nullptr) {
            vtkSmartPointer<vtkRenderer> renderer =
                vtkSmartPointer<vtkRenderer>::New();

            renderWindow->AddRenderer(renderer);
            renderer->SetViewport(viewport);
            renderer->SetBackground(1, 1, 1);

            chart = vtkSmartPointer<vtkChartXY>::New();

            chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("");
            chart->GetAxis(vtkAxis::LEFT)->SetTitle("");

            vtkSmartPointer<vtkContextScene> chartScene = vtkSmartPointer<vtkContextScene>::New();
            vtkSmartPointer<vtkContextActor> chartActor = vtkSmartPointer<vtkContextActor>::New();

            chartScene->AddItem(chart);
            chartActor->SetScene(chartScene);

            renderer->AddActor(chartActor);
            chartScene->SetRenderer(renderer);
        }

        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetColor(colors[ncolor], colors[ncolor+1], colors[ncolor+2]);
        line->SetInputData(table, 0, 1);

    }


    ncolor = (ncolor + 3)%(3*maxColor);
#endif
}

void plt::semilogx(const matrix& x, std::string label) {
    if (noplot) return;
    if (x.nrows() < 2) return;
#if ENABLE_PLT

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, pyx);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("semilogx", args, kwargs);
#endif
#ifdef USE_VTK
    for (auto j=0; j<x.ncols(); j++) {
        vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

        vtkSmartPointer<vtkFloatArray> arrX =
            vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("");
        table->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY =
            vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName(label.c_str());
        table->AddColumn(arrY);

        table->SetNumberOfRows(x.nrows());
        for (auto i=0; i<x.nrows(); i++) {
            table->SetValue(i, 0, i);
            table->SetValue(i, 1, x(i, j));
        }


        if (chart == nullptr) {
            vtkSmartPointer<vtkRenderer> renderer =
                vtkSmartPointer<vtkRenderer>::New();

            renderWindow->AddRenderer(renderer);
            renderer->SetViewport(viewport);
            renderer->SetBackground(1, 1, 1);

            chart = vtkSmartPointer<vtkChartXY>::New();

            chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("");
            chart->GetAxis(vtkAxis::LEFT)->SetTitle("");

            chart->GetAxis(vtkAxis::BOTTOM)->SetLogScale(true);

            vtkSmartPointer<vtkContextScene> chartScene = vtkSmartPointer<vtkContextScene>::New();
            vtkSmartPointer<vtkContextActor> chartActor = vtkSmartPointer<vtkContextActor>::New();

            chartScene->AddItem(chart);
            chartActor->SetScene(chartScene);

            renderer->AddActor(chartActor);
            chartScene->SetRenderer(renderer);
        }

        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetColor(colors[ncolor], colors[ncolor+1], colors[ncolor+2]);
        line->SetInputData(table, 0, 1);

    }


    ncolor = (ncolor + 3)%(3*maxColor);
#endif
}

void plt::semilogx(const matrix& x, const matrix& y, std::string label) {
    assert(x.nrows() == y.nrows());
    assert(x.ncols() == y.ncols());
    if (noplot) return;
    if (x.nrows() < 2) return;
#if ENABLE_PLT

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("semilogx", args, kwargs);
#endif
#ifdef USE_VTK
    for (auto j=0; j<x.ncols(); j++) {
        vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

        vtkSmartPointer<vtkFloatArray> arrX =
            vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("");
        table->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY =
            vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName(label.c_str());
        table->AddColumn(arrY);

        table->SetNumberOfRows(x.nrows());
        for (auto i=0; i<x.nrows(); i++) {
            table->SetValue(i, 0, x(i, j));
            table->SetValue(i, 1, y(i, j));
        }


        if (chart == nullptr) {
            vtkSmartPointer<vtkRenderer> renderer =
                vtkSmartPointer<vtkRenderer>::New();

            renderWindow->AddRenderer(renderer);
            renderer->SetViewport(viewport);
            renderer->SetBackground(1, 1, 1);

            chart = vtkSmartPointer<vtkChartXY>::New();

            chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("");
            chart->GetAxis(vtkAxis::LEFT)->SetTitle("");

            chart->GetAxis(vtkAxis::BOTTOM)->SetLogScale(true);

            vtkSmartPointer<vtkContextScene> chartScene = vtkSmartPointer<vtkContextScene>::New();
            vtkSmartPointer<vtkContextActor> chartActor = vtkSmartPointer<vtkContextActor>::New();

            chartScene->AddItem(chart);
            chartActor->SetScene(chartScene);

            renderer->AddActor(chartActor);
            chartScene->SetRenderer(renderer);
        }

        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetColor(colors[ncolor], colors[ncolor+1], colors[ncolor+2]);
        line->SetInputData(table, 0, 1);

    }


    ncolor = (ncolor + 3)%(3*maxColor);
#endif
}

void plt::semilogy(const matrix& x, std::string label) {
    if (noplot) return;
    if (x.nrows() < 2) return;
#if ENABLE_PLT

    PyObject *pyx = matrix_to_py(x);

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, pyx);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("semilogy", args, kwargs);
#endif
#ifdef USE_VTK
    for (auto j=0; j<x.ncols(); j++) {
        vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

        vtkSmartPointer<vtkFloatArray> arrX =
            vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("");
        table->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY =
            vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName(label.c_str());
        table->AddColumn(arrY);

        table->SetNumberOfRows(x.nrows());
        for (auto i=0; i<x.nrows(); i++) {
            table->SetValue(i, 0, i);
            table->SetValue(i, 1, x(i, j));
        }


        if (chart == nullptr) {
            vtkSmartPointer<vtkRenderer> renderer =
                vtkSmartPointer<vtkRenderer>::New();

            renderWindow->AddRenderer(renderer);
            renderer->SetViewport(viewport);
            renderer->SetBackground(1, 1, 1);

            chart = vtkSmartPointer<vtkChartXY>::New();

            chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("");
            chart->GetAxis(vtkAxis::LEFT)->SetTitle("");

            chart->GetAxis(vtkAxis::LEFT)->SetLogScale(true);

            vtkSmartPointer<vtkContextScene> chartScene = vtkSmartPointer<vtkContextScene>::New();
            vtkSmartPointer<vtkContextActor> chartActor = vtkSmartPointer<vtkContextActor>::New();

            chartScene->AddItem(chart);
            chartActor->SetScene(chartScene);

            renderer->AddActor(chartActor);
            chartScene->SetRenderer(renderer);
        }

        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetColor(colors[ncolor], colors[ncolor+1], colors[ncolor+2]);
        line->SetInputData(table, 0, 1);

    }


    ncolor = (ncolor + 3)%(3*maxColor);
#endif
}

void plt::semilogy(const matrix& x, const matrix& y, std::string label) {
    assert(x.nrows() == y.nrows());
    assert(x.ncols() == y.ncols());
    if (noplot) return;
    if (x.nrows() < 2) return;
#if ENABLE_PLT

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("semilogy", args, kwargs);
#endif
#ifdef USE_VTK
    ester_warn("semilogy not yet implemented");
#endif
}

void plt::loglog(const matrix& x, const matrix& y, std::string label) {
    assert(x.nrows() == y.nrows());
    assert(x.ncols() == y.ncols());
    if (noplot) return;
    if (x.nrows() < 2) return;
#if ENABLE_PLT

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);

    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "label", PyUnicode_FromString(label.c_str()));

    call("loglog", args, kwargs);
#endif
#ifdef USE_VTK
    ester_warn("loglog not yet implemented");
#endif
}

void plt::legend(std::string loc) {
    if (noplot) return;
#if ENABLE_PLT

    if (loc == "") {
        call("legend");
    }
    else {
        PyObject *args = PyTuple_New(0);

        PyObject *kwargs = PyDict_New();
        PyDict_SetItemString(kwargs, "loc", PyUnicode_FromString(loc.c_str()));
        call("legend", args, kwargs);
    }
#endif
#ifdef USE_VTK
    if (chart != nullptr) {
        chart->SetShowLegend(true);
    }
#endif
}

void plt::colorbar() {
    if (noplot) return;
#if ENABLE_PLT

    call("colorbar");
#endif
#ifdef USE_VTK
    ester_warn("colorbar not yet implemented");
#endif
}

void plt::title(const std::string& title) {
    if (noplot) return;
#if ENABLE_PLT

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyUnicode_FromString(title.c_str()));
    call("title", args);
#endif
#ifdef USE_VTK
    if (chart != nullptr) {
        chart->SetTitle(title);
    }
    else {
        ester_warn("Cannot set chart title (%s) before calls to \"plot\", \"semilogx\"...",
                title.c_str());
    }
#endif
}

void plt::clf() {
    if (noplot) return;
#if ENABLE_PLT
    call("clf");
#endif
#ifdef USE_VTK
    if (renderWindow != nullptr) {
        while (auto renderer = renderWindow->GetRenderers()->GetNextItem()) {
            renderWindow->RemoveRenderer(renderer);
        }
    }
#endif
}

void plt::ion() {
    if (noplot) return;
#if ENABLE_PLT

    call("ion");
#endif
#ifdef USE_VTK
    ester_warn("ion not yet implemented");
#endif
}

void plt::ioff() {
    if (noplot) return;
#if ENABLE_PLT

    call("ioff");
#endif
#ifdef USE_VTK
    ester_warn("ioff not yet implemented");
#endif
}

void plt::close() {
    if (noplot) return;
#if ENABLE_PLT

    call("close");
#endif
#ifdef USE_VTK
    ester_warn("close not yet implemented");
#endif
}

void plt::show(bool block) {
    if (noplot) return;
#if ENABLE_PLT

    PyObject *args = PyTuple_New(0);
    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "block", PyBool_FromLong((long) block));
    call("show", args, kwargs);
#endif
#ifdef USE_VTK
    if (block == true)
        renderWindowInteractor->Start();
    else
        renderWindowInteractor->Render();
#endif
}

void plt::block() {
    if (noplot) return;
    show(true);
}

void plt::pause(double t) {
    if (noplot) return;
#if ENABLE_PLT

    return;

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(t));
    call("pause", args);
#endif
#ifdef USE_VTK
     renderWindow->Render();
#endif
}

void plt::subplot(int subplot, bool clear_axis) {
    if (noplot) return;
#if ENABLE_PLT

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
#endif
#ifdef USE_VTK
    init();

    int nr =  subplot / 100;
    int nc =  (subplot - nr*100) / 10;
    int i =  subplot%10 - 1;

    int ix = i%nc;
    int iy = (nr-1) - i/nc;

    viewport[0] = ix * 1.0/nc;
    viewport[1] = iy * 1.0/nr;
    viewport[2] = (ix+1) * 1.0/nc;
    viewport[3] = (iy+1) * 1.0/nr;

    chart = nullptr;
    textRenderer = nullptr;
    ncolor = 0;
#endif
}

void plt::draw() {
#if ENABLE_PLT

    if (noplot || matplotlib == nullptr) return;

     PyObject* rc = PyObject_GetAttrString(matplotlib, "rcParams");
     if (!rc) {
         ester_warn("call to draw() failed\n");
         return;
     }

     PyObject *backendString = PyUnicode_FromString("backend");
     assert(backendString);

     PyObject *backend = PyObject_GetItem(rc, backendString);
     if (!backend) {
         ester_warn("call to draw() failed\n");
         return;
     }

     // skip the backend in _matplotlib.rcsetup.interactive_bk: test

     PyObject *pylabHelpers = PyObject_GetAttrString(matplotlib, "_pylab_helpers");
     assert(pylabHelpers);
     PyObject *gcf = PyObject_GetAttrString(pylabHelpers, "Gcf");
     assert(gcf);
     PyObject *figManager = PyObject_CallMethod(gcf, "get_active", nullptr);

     if (figManager) {
         PyObject *canvas = PyObject_GetAttrString(figManager, "canvas");
         assert(canvas);

         PyObject *figure = PyObject_GetAttrString(canvas, "figure");
         assert(figure);
         PyObject *stale = PyObject_GetAttrString(figure, "stale");

         if (stale) {
             PyObject_CallMethod(canvas, "draw_idle", nullptr);
         }
         PyObject_CallMethod(canvas, "start_event_loop", "d", 1e-2);
     }

    // backend = _plt.rcParams['backend']
    // if backend in _matplotlib.rcsetup.interactive_bk:
    //     figManager = _matplotlib._pylab_helpers.Gcf.get_active()
    //     if figManager is not None:
    //         canvas = figManager.canvas
    //         if canvas.figure.stale:
    //             canvas.draw_idle()
    //         canvas.start_event_loop(1e-2)
    //         return
#endif
#ifdef USE_VTK
     renderWindow->Render();
#endif
}

void plt::axvline(double x) {
    if (noplot) return;
#if ENABLE_PLT

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));
    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "color", PyUnicode_FromString("gray"));
    PyDict_SetItemString(kwargs, "linestyle", PyUnicode_FromString("--"));
    PyDict_SetItemString(kwargs, "alpha", PyFloat_FromDouble(0.5));
    call("axvline", args, kwargs);
#endif
#ifdef USE_VTK
    if (chart != nullptr && false) {
        vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

        vtkSmartPointer<vtkFloatArray> arrX =
            vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("");
        table->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY =
            vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName("-");
        table->AddColumn(arrY);

        table->SetNumberOfRows(2);

        table->SetValue(0, 0, x);
        table->SetValue(0, 1, chart->GetAxis(vtkAxis::LEFT)->GetMinimum());
        table->SetValue(1, 0, x);
        table->SetValue(1, 1, chart->GetAxis(vtkAxis::LEFT)->GetMaximum());
        vtkPlot *line = chart->AddPlot(vtkChart::LINE);
        line->SetColor(0.5, 0.5, 0.5);
        line->SetInputData(table, 0, 1);
    }
#endif
}

void plt::axhline(double y) {
    if (noplot) return;
#if ENABLE_PLT

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(y));
    PyObject *kwargs = PyDict_New();
    PyDict_SetItemString(kwargs, "color", PyUnicode_FromString("gray"));
    PyDict_SetItemString(kwargs, "linestyle", PyUnicode_FromString("--"));
    PyDict_SetItemString(kwargs, "alpha", PyFloat_FromDouble(0.5));
    call("axhline", args, kwargs);
#endif
#ifdef USE_VTK
    ester_warn("axhline not yet implemented");
#endif
}

void plt::text(double x, double y, std::string text) {
    if (noplot) return;
#if ENABLE_PLT
    PyObject *args = PyTuple_New(3);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));
    PyTuple_SetItem(args, 1, PyFloat_FromDouble(y));
    PyTuple_SetItem(args, 2, PyUnicode_FromString(text.c_str()));
    call("text", args);
#endif
#ifdef USE_VTK
      vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
      std::string line = text;
      for (auto i=0.0; i<y; i+=0.1)
          line.append("\n");
      textActor->SetInput(line.c_str());
      textActor->GetTextProperty()->SetFontSize(12);
      textActor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);

      if (textRenderer == nullptr) {
          textRenderer =
              vtkSmartPointer<vtkRenderer>::New();
          renderWindow->AddRenderer(textRenderer);
          textRenderer->SetViewport(viewport);
          textRenderer->SetBackground(1, 1, 1);
      }
      textRenderer->AddActor2D(textActor);
#endif
}

void plt::savefig(const std::string& filename) {
    if (noplot) return;
#if ENABLE_PLT

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyUnicode_FromString(filename.c_str()));
    call("savefig", args);
#endif
#ifdef USE_VTK
    ester_warn("savefig not yet implemented");
#endif
}

void plt::call(const std::string& name, PyObject *args, PyObject *kwargs) {
#if ENABLE_PLT
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
#endif
}

void plt::pcolormesh(const matrix& x, const matrix& y, const matrix& c) {
    if (noplot) return;
#if ENABLE_PLT

    PyObject *pyx = matrix_to_py(x);
    PyObject *pyy = matrix_to_py(y);
    PyObject *pyc = matrix_to_py(c);

    PyObject *args = PyTuple_New(3);
    PyTuple_SetItem(args, 0, pyx);
    PyTuple_SetItem(args, 1, pyy);
    PyTuple_SetItem(args, 2, pyc);

    call("pcolormesh", args);
#endif
#ifdef USE_VTK
    ester_warn("pcolormesh not yet implemented");
#endif
}

void plt::figure(const int& id, int width, int height) {
    if (noplot) return;
#if ENABLE_PLT

    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyLong_FromLong(id));

    if (width != -1 || height != -1) {
        PyObject *kwargs = PyDict_New();
        PyObject *size = PyTuple_New(2);
        PyTuple_SetItem(size, 0, PyLong_FromLong(width));
        PyTuple_SetItem(size, 1, PyLong_FromLong(height));
        PyDict_SetItemString(kwargs, "figsize", size);
        call("figure", args, kwargs);
    }
    else {
        call("figure", args);
    }
#endif
#ifdef USE_VTK
    renderWindow =
        vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetSize(width*100, height*100);

    renderWindowInteractor =
        vtkSmartPointer<vtkRenderWindowInteractor>::New();

    renderWindowInteractor->SetRenderWindow(renderWindow);
    renderWindow->SetWindowName("ESTER");

    chart = nullptr;
    textRenderer = nullptr;
#endif
}

void plt::axis(const double& x0, const double& x1, const double& y0, const double& y1) {
    if (noplot) return;
#if ENABLE_PLT
    PyObject *args = PyTuple_New(1);
    PyObject *list = PyList_New(4);

    PyList_SetItem(list, 0, PyFloat_FromDouble(x0));
    PyList_SetItem(list, 1, PyFloat_FromDouble(x1));
    PyList_SetItem(list, 2, PyFloat_FromDouble(y0));
    PyList_SetItem(list, 3, PyFloat_FromDouble(y1));

    PyTuple_SetItem(args, 0, list);
    call("axis", args);
#endif
#ifdef USE_VTK
    if (chart != nullptr) {
        chart->GetAxis(vtkAxis::BOTTOM)->SetMinimum(x0);
        chart->GetAxis(vtkAxis::BOTTOM)->SetMaximum(x1);
        chart->GetAxis(vtkAxis::LEFT)->SetMinimum(y0);
        chart->GetAxis(vtkAxis::LEFT)->SetMaximum(y1);
    }
#endif
}

void plt::axis(const std::string& a) {
    if (noplot) return;
#if ENABLE_PLT
    PyObject *args = PyTuple_New(1);
    PyTuple_SetItem(args, 0, PyUnicode_FromString(a.c_str()));
    call("axis", args);
#endif
#ifdef USE_VTK
    ester_warn("axis not yet implemented");
#endif
}
