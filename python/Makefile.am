do_subst = sed -e 's,[@]PYTHON[@],$(PYTHON),g'

if BUILD_ESTER_PY
EXTRA_DIST = ester.in ester_wrap.py ester_wrap.cxx ester.i

BUILT_SOURCES = ester.py

if HAVE_SWIG
SWIG_SRC = $(srcdir)/ester.i
BUILT_SOURCES += ester_wrap.cpp
endif

pkgpython_PYTHON = ester.py ester_wrap.py
pkgpyexec_LTLIBRARIES = _ester_wrap.la

_ester_wrap_la_SOURCES = ester_wrap.cpp
_ester_wrap_la_CPPFLAGS = -I$(top_srcdir)
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

ester.py: $(srcdir)/ester.in Makefile
	$(do_subst) < $< > $@

CLEANFILES = $(BUILT_SOURCES)

install-exec-hook:
	@echo
	@echo
	@echo
	@echo "Binaries are installed in:"
	@echo "    $(prefix)/bin"
	@echo
	@echo "Libraries are installed in:"
	@echo "    $(libdir)"
	@echo
	@echo "Python modules are in:"
	@echo "    $(pkgpythondir)"
	@echo "    $(pkgpyexecdir)"
	@echo
	@echo
	@echo

else

install-exec-hook:
	@echo
	@echo
	@echo
	@echo "Binaries are installed in:"
	@echo "    $(prefix)/bin"
	@echo
	@echo "Libraries are installed in:"
	@echo "    $(libdir)"
	@echo
	@echo
	@echo

endif
