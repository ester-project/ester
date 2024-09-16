#ifndef WITH_CMAKE
#include "ester-config.h"
#endif

#include "physics.h"

matrix composition_map::X() const {

	return (*this)["H"];

}

matrix composition_map::Y() const {

	return (*this)["He3"]+(*this)["He4"];

}

matrix composition_map::Z() const {

	return 1.-X()-Y();

}

