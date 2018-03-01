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
        module_data_sorgam.o \
        module_data_mosaic_asect.o \
        module_chem_utilities.o \
        module_interpolate.o \
        module_vertmx_wrf.o \
        module_aer_opt_out.o \
        module_tropopause.o \
        module_upper_bc_driver.o \
        module_mosaic_driver.o \
        module_aerosols_sorgam.o \
        module_gocart_aerosols.o \
        module_aerosols_soa_vbs.o \
        module_input_chem_data.o \
        module_input_tracer.o \
        module_input_chem_bioemiss.o

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
	rm -f *.o
	rm -f *.f90
	rm -f *.mod
	@echo "Cleaning chem may not be enough - check subdirectories."

devclean:
	@echo "This is for development purposes only.\n"
	git pull origin master
	rm -f *.f90
	rm -f *.o
	rm -f *.mod
	git checkout -- .
	@echo "Done"

about:
	@echo "*****************************************************************"
	@echo "  __          _______  ______       _____  _____ _    _ _____    "
	@echo "  \ \        / /  __ \|  ____|     / ____|/ ____| |  | |  __ \   "
	@echo "   \ \  /\  / /| |__) | |__ ______| |  __| |    | |__| | |__) |  "
	@echo "    \ \/  \/ / |  _  /|  __|______| | |_ | |    |  __  |  ___/   "
	@echo "     \  /\  /  | | \ \| |         | |__| | |____| |  | | |       "
	@echo "      \/  \/   |_|  \_\_|          \_____|\_____|_|  |_|_|       "
	@echo "*****************************************************************"
	@echo "   THIS IS THE WRF-GCHP 'PUMPKIN' CHEMISTRY ABSTRACTION LAYER    "
	@echo "                    FOR THE WRF MODEL VERSION 3                  "
	@echo "*****************************************************************"
	@echo " FOR ERRORS, SUGGESTIONS AND FEEDBACK, CONTACT HAIPENG LIN AT    "
	@echo "           LINHAIPENG@PKU.EDU.CN | JIMMIE.LIN@GMAIL.COM          "
	@echo "*****************************************************************"
	@echo " (c) 2018 Haipeng Lin                                            "
	@echo " Peking University, Atmospheric Chemistry and Climate Group      "
	@echo "*****************************************************************"

include depend.chem

