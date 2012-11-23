ESTER_DIR=../..

all:
	sed -e '1a#include"starR.h"' -e 's/star1d/star1dR/g' $(ESTER_DIR)/src/main/star1d.cpp >star1dR.cpp
	sed -e '1a#include"starR.h"' -e 's/star2d/star2dR/g' $(ESTER_DIR)/src/main/star2d.cpp >star2dR.cpp
	$(ESTER_DIR)/utils/ester_build star1dR.cpp star1dR_class.cpp -o $(ESTER_DIR)/bin/star1dR
	$(ESTER_DIR)/utils/ester_build star2dR.cpp star2dR_class.cpp -o $(ESTER_DIR)/bin/star2dR
