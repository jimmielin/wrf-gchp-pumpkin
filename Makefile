#
#  Changing this file requires mods to 
#	WRFV3_top/chem/KPP/util/wkc/change_chem_Makefile.c
#	Right now it looks at the line in MODULES with module_data_sorgam
#	AND it cares about the word DEPENDENCIES.  DO NOT EVEN CHANGE
#	SPACING ON THESE TWO LINES.

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f

MODULES =                                 \
        module_data_radm2.o \
        module_data_soa_vbs.o \
        module_data_sorgam.o \
        module_chem_utilities.o \
        module_interpolate.o \
        module_vertmx_wrf.o \
        module_aer_opt_out.o \
        module_tropopause.o \
        module_upper_bc_driver.o \
        module_data_mosaic_asect.o \
        module_mosaic_driver.o \
        module_radm.o \
        module_dep_simple.o \
        module_mosaic_wetscav.o \
        module_aerosols_sorgam.o \
        module_gocart_aerosols.o \
        module_aerosols_soa_vbs.o \
        module_input_chem_data.o \
        module_input_tracer.o

OBJS    =                           \
        chemics_init.o              \
        chem_driver.o

LIBTARGET    =  chemics
TARGETDIR    =  ./
$(LIBTARGET) :  MODULE DRIVERS
		$(AR) $(ARFLAGS) ../main/$(LIBWRFLIB) $(MODULES) $(OBJS)

MODULE  : $(MODULES)

DRIVERS : $(OBJS)

include ../configure.wrf

clean:
	@ echo 'use the clean script'

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

include depend.chem

