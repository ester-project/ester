do_subst = sed -e 's,[@]PYTHON[@],$(PYTHON),g'

if BUILD_ESTER_PY
EXTRA_DIST = ester.in ester_wrap.py ester_wrap.cxx ester.i num.in cesam.in

BUILT_SOURCES = __init__.py cesam.py num.py

if HAVE_SWIG
SWIG_SRC = $(srcdir)/ester.i
BUILT_SOURCES += ester_wrap.cpp
endif

pkgpython_PYTHON = __init__.py ester_wrap.py num.py cesam.py
pkgpyexec_LTLIBRARIES = _ester_wrap.la

_ester_wrap_la_SOURCES = ester_wrap.cpp
_ester_wrap_la_CPPFLAGS = -I$(top_srcdir)/src/include -I$(top_srcdir)/src/graphics
_ester_wrap_la_LDFLAGS = -module
_ester_wrap_la_LIBADD = $(top_builddir)/src/libester.la

if HAVE_SWIG
ester_wrap.cpp: $(SWIG_SRC)
	$(SWIG) -c++ -python $(SWIG_PYTHON_OPT) -I$(top_srcdir) -o $@ $<
else
ester_wrap.cpp: $(srcdir)/ester_wrap.cxx
	cp $< $@
ester_wrap.py: $(srcdir)/ester_wrap.py
	cp $< $@
endif

__init__.py: $(srcdir)/ester.in Makefile
	$(do_subst) < $< > $@

num.py: $(srcdir)/num.in Makefile
	$(do_subst) < $< > $@

cesam.py: $(srcdir)/cesam.in Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
endif

install-exec-hook:
	@echo
	@echo
	@echo
	@echo "----------------------------------------------------------------------"
	@echo "Binaries are installed in:"
	@echo "    $(prefix)/bin"
	@echo "    (add it to your PATH)"
	@echo
	@echo "Libraries are installed in:"
	@echo "    $(libdir)"
	@echo
if BUILD_ESTER_PY
	@echo "Python modules are in:"
	@echo "    $(pkgpythondir)"
	@echo "    $(pkgpyexecdir)"
	@echo "    (add it to your PYTHONPATH)"
else
	@echo "Python modules were not built"
endif
	@echo
	@echo "Examples have been built in:"
	@echo "    $(abs_top_builddir)/lib-example/bin"
	@echo "----------------------------------------------------------------------"
	@echo
	@echo
	@echo
