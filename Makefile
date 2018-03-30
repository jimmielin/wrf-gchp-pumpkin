################################################################################
#
#   WRF-GCHP
#   GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
#
#   WRF & GCHP are (c) their original authors.
#   WRF-GCHP coupling layer (WGCL) is (c) Atmospheric Chemistry and Climate Group, Peking University
#
#   Developed by Haipeng Lin <linhaipeng@pku.edu.cn>, Xu Feng, 2018-01
#   Peking University, School of Physics
#
################################################################################
#
#   Codename Pumpkin: Abstracted Bindings for Chemistry-to-WRF
#
#   This Chemical Interface (chem/) is written after comprehensive study of
#   the original chem_driver.f from WRF-Chem v3.6.1
#   which is (c) their respective authors.
#
################################################################################
#
#  Makefile
#
#  A generic makefile that compiles a hollow "chem" shell.
#  You will want to customize this so it compiles your chemistry along with this
#  component.
#
################################################################################

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f

MODULES =                                 \
        module_data_radm2.o \
        module_data_sorgam.o \
        module_data_mosaic_asect.o \
        module_chem_utilities.o \
        module_interpolate.o \
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

.PHONY: devclean install_registry compile_chem

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

install_registry:
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
	@echo " (c) 2018 Haipeng Lin                                            "
	@echo " Peking University, Atmospheric Chemistry and Climate Group      "
	@echo "*****************************************************************"
	@echo "THIS WILL INSTALL AND REPLACE THE WRF-CHEM STANDARD REGISTRY.    "
	mv ../Registry/registry.chem ../Registry/registry.chem.bak
	cp ./registry.chem ../Registry/registry.chem
	@echo "Your original registry is now in registry.chem.bak.              "
	@echo "To compile, return to root directory and run ./compile em_real   "


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
	@echo "Commands:                                                        "
	@echo "    make about - Show this about (help) screen                   "
	@echo "    make clean - Clean chemistry (might not be thorough)         "
	@echo "    make install_registry - Install chemistry species into WRF   "
	@echo "       (replacing existing Registry.chem)                        "
	@echo "    make compile_chem - Compile target chemistry                 "
	@echo "    make devclean - Dev purposes only, get from git origin/master"
	@echo "                                                                 "
	@echo "  * Not all commands, especially install_registry, compile_chem  "
	@echo "    are always available. Check with your chemistry maintainer.  "
	@echo "    These commands are mostly developed for WRF-GCHP's GIGC part."
	@echo "  * If you compiled WRF without install_registry, you need to    "
	@echo "    clean all (./clean -a) to rebuild module_state_description.  "
	@echo "    A non-clean compile cannot refresh the registry.             "
	@echo "*****************************************************************"

# DEPENDENCIES : only dependencies after this line
module_data_radm2.o:

module_data_sorgam.o: module_data_radm2.o

# Required by module_cu_kfcup
module_data_mosaic_asect.o:

module_chem_utilities.o:

module_interpolate.o:

module_vertmx_wrf.o:

module_aer_opt_out.o:

module_tropopause.o: module_interpolate.o

module_upper_bc_driver.o: module_tropopause.o

module_input_chem_bioemiss.o:

module_mosaic_driver.o:

module_aerosols_sorgam.o:

module_aerosols_soa_vbs.o:

module_input_chem_data.o: module_data_sorgam.o

chemics_init.o: module_input_chem_data.o module_tropopause.o module_upper_bc_driver.o

chem_driver.o: ../dyn_em/module_convtrans_prep.o module_input_chem_data.o module_chem_utilities.o module_tropopause.o module_upper_bc_driver.o
