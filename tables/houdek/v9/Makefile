#
# Makefile for opacity-interpolation-package
#
# G. Houdek
# 16.9.1992
#
#--------------------------------------------
#
RANLIB	= ranlib	#use ranlib
#RANLIB	= echo		# no ranlib
#

FC     = gfortran
#FFLAGS = -cxdb -fi -cs
#FFLAGS = -native -g -C
#FFLAGS = -fast
FFLAGS =  -O2 -fd-lines-as-comments
#FFLAGS = -O1 -32
#FFLAGS = -fast -xld -C

.KEEP_STATE:

all: dopxext_s dopxext_s_95 doptesf dopdof dopdos exakop95\
     da2b dexak2 kur2l2g a2b95 exop95 alex94bext kur2l2g_95
     

#--------------------------------------------

TAB_DIR    = opal92
TAB_DIR_95 = opal95
ALX_DIR    = alex94
KUR_DIR    = kuru91
RKS_DIR    = rksuite
LIB_DIR	   = lib
LIB_LIB	   = libopint.a
RKS_LIB    = librksuite.a

#RKS_DEX = -l$(RKS_LIB:lib%.a=%)
#LIB_DEX = -l$(LIB_LIB:lib%.a=%)

RKS_DEX = -lrksuite
LIB_DEX = -lopint

RKS_DIL = -L$(RKS_DIR)
LIB_DIL = -L$(LIB_DIR)

LIB_OBJ = condux_ad.o     dopinit.o       dreadkx.o       ratmat.o\
          condux_cf.o     dopintc.o       dshep.o         rbival.o\
          dadjac.o        dopintf.o       duvip3d.o       rbivpd.o\
          dctang.o        dopints.o       duvip3p.o       rde95.o \
          dextrp.o        dordgr.o        gi.o            rdi95.o \
          dhunt.o         dordgri.o       gid.o           readal.o\
          dintrp.o        dpderix.o       inttwo.o        rtperm.o\
          dintrpi.o       dpderiy.o       itoh_coeff_cf.o sbilin.o\
          dmaceps.o       dpdmin.o        itoh_ec_ad.o    sfa_ad.o\
          dmasub.o        dpdste.o        itoh_ec_cf.o    sfa_cf.o\
          dmasube.o       dreade.o        itoh_init_ad.o  trdisa.o\
          dmasubi.o       dreadi.o        itoh_init_cf.o  trdisb.o\
          dmaxmn.o        dreadk.o        rat2d.o         dcspline.o

RKS_OBJ = chkfl.o         glberr.o        softfl.o        stiffb.o\
          const.o         intpo.o         stat.o          stiffc.o\
          ct.o            mconst.o        step.o          stiffd.o\
          dotprd.o        reset.o         stepa.o         truerr.o\
          envirn.o        rkmsg.o         stepb.o         ut.o    \
          evali.o         rksit.o         stiff.o                 \
          formi.o         setup.o         stiffa.o                

LIB_SRC = $(LIB_OBJ:%.o=$(LIB_DIR)/%.f)
RKS_SRC = $(RKS_OBJ:%.o=$(RKS_DIR)/%.f)

lib:    $(LIB_DIR)/$(LIB_LIB)
rks:    $(RKS_DIR)/$(RKS_LIB)

dopdof: $(LIB_DIR)/$(LIB_LIB) dopdof.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)

dopdos: $(LIB_DIR)/$(LIB_LIB) dopdos.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)

doptesf:$(LIB_DIR)/$(LIB_LIB) doptesf.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)

dpderiv:$(LIB_DIR)/$(LIB_LIB) dpderiv.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)

dopxext_s:$(LIB_DIR)/$(LIB_LIB) $(RKS_DIR)/$(RKS_LIB) dopxext_s.o
	$(FC) $(FFLAGS) -o $@ $@.o $(RKS_DIL) $(LIB_DIL) $(RKS_DEX) $(LIB_DEX)

dopxext_s_95:$(LIB_DIR)/$(LIB_LIB) $(RKS_DIR)/$(RKS_LIB) dopxext_s_95.o
	$(FC) $(FFLAGS) -o $@ $@.o $(RKS_DIL) $(LIB_DIL) $(RKS_DEX) $(LIB_DEX)

dexak2: $(LIB_DIR)/$(LIB_LIB) dexak2.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)

exop95: $(LIB_DIR)/$(LIB_LIB) exop95.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)
	
exakop95: $(LIB_DIR)/$(LIB_LIB) exakop95.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)

da2b:   da2b.o
	$(FC) $(FFLAGS) -o $@ $@.o

a2b95:  a2b95.o
	$(FC) $(FFLAGS) -o $@ $@.o

alex94bext: $(LIB_DIR)/$(LIB_LIB) alex94bext.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)

dopdalex94: $(LIB_DIR)/$(LIB_LIB) dopdalex94.o
	$(FC) $(FFLAGS) -o $@ $@.o $(LIB_DIL) $(LIB_DEX)

kur2l2g:$(LIB_DIR)/$(LIB_LIB) kur2l2g.o kurset.o
	$(FC) $(FFLAGS) -o $@ $@.o kurset.o $(LIB_DIL) $(LIB_DEX)

kur2l2g_95:$(LIB_DIR)/$(LIB_LIB) kur2l2g_95.o kurset.o
	$(FC) $(FFLAGS) -o $@ $@.o kurset.o $(LIB_DIL) $(LIB_DEX)

$(LIB_DIR)/$(LIB_LIB): FORCE
	@(cd $(@D); make FC="$(FC)" LIB="$(LIB_LIB)" LIB_OBJ="$(LIB_OBJ)" \
                             FFLAGS="$(FFLAGS)" RANLIB="$(RANLIB)" \
                             F77FLAGS="$(FFLAGS)" \
        )

$(RKS_DIR)/$(RKS_LIB): FORCE
	@(cd $(@D); make FC="$(FC)" LIB="$(RKS_LIB)" LIB_OBJ="$(RKS_OBJ)" \
                             FFLAGS="$(FFLAGS)" RANLIB="$(RANLIB)" \
                             F77FLAGS="$(FFLAGS)" \
        )

tab92:  da2b dexak2 kur2l2g dopxext_s
	(cd $(TAB_DIR); ../da2b; ../dexak2)
	(cd $(KUR_DIR); ../kur2l2g; ../dopxext_s)

tab95:  a2b95 exakop95 kur2l2g_95 dopxext_s_95
	(cd $(TAB_DIR_95); ../a2b95; ../exakop95)
	(cd $(KUR_DIR); ../kur2l2g_95; ../dopxext_s_95)

tab95a: alex94bext dopdalex94
	(cd $(ALX_DIR); ../alex94bext)
	(./dopdalex94)

FORCE:

clean:
	touch dummy.trace
	touch $(LIB_DIR)/dummy.o
	touch $(RKS_DIR)/dummy.o
	touch $(TAB_DIR)/dummy.o
	touch $(KUR_DIR)/dummy.o
	rm -f  dopdof dopdos dopxext_s dopxext_s_95 dopdalex94 \
               doptesf da2b dexak2 kur2l2g kur2l2g_95 exakop95 \
               a2b95 exop95 alex94bext                         \
               *.o fort.* core p*.dat *.trace *.ps *.er *.lst *.fln
	rm -f $(LIB_DIR)/*.o $(LIB_DIR)/*.a
	rm -f $(RKS_DIR)/*.o $(RKS_DIR)/*.a
	rm -f $(TAB_DIR)/*.o $(TAB_DIR)/fort.* $(TAB_DIR)/core
	rm -f $(KUR_DIR)/*.o $(KUR_DIR)/fort.* $(KUR_DIR)/core
