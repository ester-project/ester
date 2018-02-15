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
#include "../src/include/matrix.h"
#include "../src/include/star.h"
#include "../src/include/physics.h"
#include "../src/include/mapping.h"
#include "../src/include/numdiff.h"
#include "../src/graphics/matplotlib.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

%include "../src/include/matrix.h"
%include "../src/include/physics.h"
%include "../src/include/mapping.h"
%include "../src/include/numdiff.h"
%include "../src/include/star.h"
%include "../src/graphics/matplotlib.h"
