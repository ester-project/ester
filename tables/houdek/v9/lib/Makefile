#
# Makefile for the lib-dir of the 
# opacity-interpolation package
#
# G. Houdek
# 16.9.1992
#
#----------------------------------------
LIB  =  libopint.a
#----------------------------------------
$(LIB): $(LIB_OBJ)
	( \
         ar ruv $(LIB) $? ;\
         $(RANLIB) $(LIB) ;\
)
