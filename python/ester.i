%module ester_wrap

#pragma SWIG nowarn=362

%ignore mapping::ext_map;
%ignore double_map;
%ignore create_double_map;
%ignore matrix_map_elem;

%ignore operator*;
%ignore operator+;
%ignore operator-;
%ignore operator/;
%ignore operator,;
%ignore operator=;
%ignore operator*=;
%ignore operator+=;
%ignore operator&&;
%ignore operator||;
%ignore operator[];
%ignore operator();
%ignore operator==;
%ignore operator!=;
%ignore operator>=;
%ignore operator<=;
%ignore operator>;
%ignore operator<;
%ignore operator matrix_map;
%ignore operator double_map;
%ignore operator matrix;
%ignore star2d::version_struct;
%ignore star2d::units_struct;
%ignore star2d::config_struct;

%{
#include "matrix.h"
#include "star.h"
#include "physics.h"
#include "mapping.h"
#include "numdiff.h"
#include "matplotlib.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

%include "matrix.h"
%include "physics.h"
%include "mapping.h"
%include "numdiff.h"
%include "star.h"
%include "matplotlib.h"
