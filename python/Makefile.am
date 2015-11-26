do_subst = sed -e 's,[@]PYTHON[@],$(PYTHON),g'

EXTRA_DIST = ester.in

BUILT_SOURCES = ester_wrap.cpp ester.py
SWIG_SOURCES = ester_wrap.i

pkgpython_PYTHON = ester_wrap.py ester.py
pkgpyexec_LTLIBRARIES = _ester_wrap.la

_ester_wrap_la_SOURCES = ester_wrap.cpp $(SWIG_SOURCES)
_ester_wrap_la_CPPFLAGS = $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)
_ester_wrap_la_LDFLAGS = -module
_ester_wrap_la_LIBADD = $(top_builddir)/src/libester.la

ester_wrap.cpp: $(SWIG_SOURCES)
	$(SWIG) -c++ -python $(SWIG_PYTHON_OPT) -I$(top_srcdir) -o $@ $<

ester.py: $(srcdir)/ester.in Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)
