%module ester_wrap

%{
#include "../src/include/star.h"
%}

%include "carrays.i"
%array_class(double, doubleArray);

%include "../src/include/physics.h"
%include "../src/include/matrix.h"
%include "../src/include/mapping.h"
%include "../src/include/star.h"
