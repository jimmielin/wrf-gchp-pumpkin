!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   WRF-GCHP
!   GEOS-Chem High Performance-powered Chemistry Add-On for WRF Model
!
!   WRF & GCHP are (c) their original authors.
!   WRF-GCHP coupling layer (WGCL) is (c) Atmospheric Chemistry and Climate Group, Peking University
!
!   Developed by Haipeng Lin <linhaipeng@pku.edu.cn>, Xu Feng, 2018-01
!   Peking University, School of Physics
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Codename Pumpkin: Abstracted Bindings for Chemistry-to-WRF
!
!   This Chemical Interface (chem/) is written after comprehensive study of
!   the original chem_driver.f from WRF-Chem v3.6.1
!   which is (c) their respective authors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  MODULE: module_input_chem_data
!  DESCRIPTION: Input Chemistry Data Module for "Pumpkin" Abstraction Layer
!               Satisfying all WRF Calls through stubbing or redirection to external
!               parameters, as of WRF v3.9.1.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_input_chem_data
   use module_io_domain
   use module_domain
   use module_get_file_names, only:eligible_file_name, number_of_eligible_files, unix_ls

   implicit none

   ! last_chem_time
   ! req. by chem/chem_driver for timestep calculations.
   type(WRFU_Time), dimension(max_domains) :: last_chem_time

contains
   ! get_last_gas
   ! Get the index of the last gas species depending on mechanism. req. by chem/chem_driver
   !!!! CHEMISTRY DEVELOPERS: YOU MUST UPDATE THIS FOR THE RIGHT INDEXES DEP. ON YOUR REGISTRY.
   integer function get_last_gas(chem_opt)
      implicit none
      integer, intent(in) :: chem_opt

      select case (chem_opt)
      case (0)
         get_last_gas = 0
      case (1)
         get_last_gas = p_ho2

      case default
         call wrf_error_fatal("Pumpkin module_input_chem_data::get_last_gas: could not decipher chem_opt value")

      end select

   end function get_last_gas

   ! setup_gasprofile_maps
   ! Sets up the cross reference mapping indices and fractional
   ! apportionment of the default species profiles for use with
   ! ICs and BCs.
   ! req. by chem/chemics_init
   !!!! CHEMISTRY DEVELOPERS: YOU MUST UPDATE THIS FOR THE RIGHT INDEXES DEP. ON YOUR REGISTRY.
   subroutine setup_gasprofile_maps(chem_opt, numgas)
      integer, intent(in) :: chem_opt, numgas
      select case (chem_opt)
      case (1)
         !! CHEMISTRY DEVELOPERS: CHEMISTRY NEEDS TO BE ADDED HERE.
         !! Either by including a new module file (good practice) or hardcode it...
         ! call setup_gasprofile_map_geoschem
      end select
   end subroutine setup_gasprofile_maps

   ! an example for setup_gasprofile_maps is set below.
   ! you should preferably add this to a module file!
   ! this is setup_gasprofile_map_radm_racm from WRF-Chem 3.6.1
   ! (c) original authors
   ! subroutine setup_gasprofile_map_radm_racm

   !    iref(:) = 7 !default value
   !    iref(1:41) = (/12, 19, 2, 2, 1, 3, 4, 9, 8, 5, 5, 32, 6, 6, 6, 30, 30, 10, 26, 13, 11, 6, 6, &
   !                   14, 15, 15, 23, 23, 32, 16, 23, 31, 17, 23, 23, 23, 23, 23, 7, 28, 29/)

   !    fracref(:) = 1. !default value
   !    fracref(1:41) = (/1., 1., .75, .25, 1., 1., 1., 1., 1., 1., &
   !                      .5, .5, 6.25E-4, 7.5E-4, 6.25E-5, .1, &
   !                      .9, 1., 1., 1., 1., 8.E-3, 1., 1., 1., .5, &
   !                      1., 1., .5, 1., 1., 1., 1., 1., 1., 1., 1., &
   !                      1., 1., 1., 1./)

   !    ggnam(:) = 'JUNK' !default value
   !    ggnam(1:41) = (/'SO2 ', 'SULF', 'NO2 ', 'NO  ', 'O3  ', 'HNO3', &
   !                    'H2O2', 'ALD ', 'HCHO', 'OP1 ', 'OP2 ', 'PAA ', &
   !                    'ORA1', 'ORA2', 'NH3 ', 'N2O5', 'NO3 ', 'PAN ', &
   !                    'HC3 ', 'HC5 ', 'HC8 ', 'ETH ', 'CO  ', 'OL2 ', &
   !                    'OLT ', 'OLI ', 'TOL ', 'XYL ', 'ACO3', 'TPAN', &
   !                    'HONO', 'HNO4', 'KET ', 'GLY ', 'MGLY', 'DCB ', &
   !                    'ONIT', 'CSL ', 'ISO ', 'HO  ', 'HO2 '/)

   ! end subroutine setup_gasprofile_map_radm_racm

   ! chem_dbg
   ! Chemistry Debug function req. by chem/chem_driver
   !!!! CHEMISTRY DEVELOPERS: YOU WILL LIKELY NEED TO UPDATE THIS SPECIES LIST FOR DEBUGGING.
   !!!! THIS IS SEPARATE FROM REGISTRY.CHEM
#ifdef CHEM_DBG_I
   subroutine chem_dbg(i, j, k, dtstep, itimestep, &
                       dz8w, t_phy, p_phy, rho_phy, chem, &
                       emis_ant, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       kemit, &
                       ph_macr, ph_o31d, ph_o33p, ph_no2, ph_no3o2, ph_no3o, ph_hno2, &
                       ph_hno3, ph_hno4, ph_h2o2, ph_ch2or, ph_ch2om, ph_ch3cho, &
                       ph_ch3coch3, ph_ch3coc2h5, ph_hcocho, ph_ch3cocho, &
                       ph_hcochest, ph_ch3o2h, ph_ch3coo2h, ph_ch3ono2, ph_hcochob, ph_n2o5, &
                       ph_o2)

      implicit none
      INTEGER, INTENT(IN) :: i, j, k, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte, &
                             kemit
      real, intent(in) :: dtstep
      integer, intent(in) :: itimestep
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(INOUT) :: chem
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: dz8w, t_phy, p_phy, rho_phy
      REAL, DIMENSION(ims:ime, kms:kemit, jms:jme, num_emis_ant), &
         INTENT(IN) :: emis_ant
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN), OPTIONAL :: &
         ph_macr, ph_o31d, ph_o33p, ph_no2, ph_no3o2, ph_no3o, ph_hno2, &
         ph_hno3, ph_hno4, ph_h2o2, ph_ch2or, ph_ch2om, ph_ch3cho, &
         ph_ch3coch3, ph_ch3coc2h5, ph_hcocho, ph_ch3cocho, &
         ph_hcochest, ph_ch3o2h, ph_ch3coo2h, ph_ch3ono2, ph_hcochob, ph_n2o5, &
         ph_o2

      integer :: n
      real :: conva, convg

      print *, "itimestep =", itimestep

      print *, "MET DATA AT (i,k,j):", i, k, j
      print *, "t_phy,p_phy,rho_phy=", t_phy(i, k, j), p_phy(i, k, j), rho_phy(i, k, j)

      if (dz8w(i, k, j) /= 0.) then
         conva = dtstep/(dz8w(i, k, j)*60.)
         convg = 4.828e-4/rho_phy(i, k, j)*dtstep/(dz8w(i, k, j)*60.)
         print *, "ADJUSTED EMISSIONS (PPM) AT (i,k,j):", i, k, j
         print *, "dtstep,dz8w(i,k,j):", dtstep, dz8w(i, k, j)
         print *, "e_pm25 i,j:", emis_ant(i, k, j, p_e_pm25i)*conva, &
            emis_ant(i, k, j, p_e_pm25j)*conva
         print *, "e_ec i,j:", emis_ant(i, k, j, p_e_eci)*conva, &
            emis_ant(i, k, j, p_e_ecj)*conva
         print *, "e_org i,j:", emis_ant(i, k, j, p_e_orgi)*conva, &
            emis_ant(i, k, j, p_e_orgj)*conva
         print *, "e_so2:", emis_ant(i, k, j, p_e_so2)*convg
         print *, "e_no:", emis_ant(i, k, j, p_e_no)*convg
         print *, "e_co:", emis_ant(i, k, j, p_e_co)*convg
         print *, "e_eth:", emis_ant(i, k, j, p_e_eth)*convg
         print *, "e_hc3:", emis_ant(i, k, j, p_e_hc3)*convg
         print *, "e_hc5:", emis_ant(i, k, j, p_e_hc5)*convg
         print *, "e_hc8:", emis_ant(i, k, j, p_e_hc8)*convg
         print *, "e_xyl:", emis_ant(i, k, j, p_e_xyl)*convg
         print *, "e_ol2:", emis_ant(i, k, j, p_e_ol2)*convg
         print *, "e_olt:", emis_ant(i, k, j, p_e_olt)*convg
         print *, "e_oli:", emis_ant(i, k, j, p_e_oli)*convg
         print *, "e_tol:", emis_ant(i, k, j, p_e_tol)*convg
         print *, "e_csl:", emis_ant(i, k, j, p_e_csl)*convg
         print *, "e_hcho:", emis_ant(i, k, j, p_e_hcho)*convg
         print *, "e_ald:", emis_ant(i, k, j, p_e_ald)*convg
         print *, "e_ket:", emis_ant(i, k, j, p_e_ket)*convg
         print *, "e_ora2:", emis_ant(i, k, j, p_e_ora2)*convg
         print *, "e_pm25:", emis_ant(i, k, j, p_e_pm_25)*conva
         print *, "e_pm10:", emis_ant(i, k, j, p_e_pm_10)*conva
         print *, "e_nh3:", emis_ant(i, k, j, p_e_nh3)*convg
         print *, "e_no2:", emis_ant(i, k, j, p_e_no2)*convg
         print *, "e_ch3oh:", emis_ant(i, k, j, p_e_ch3oh)*convg
         print *, "e_c2h5oh:", emis_ant(i, k, j, p_e_c2h5oh)*convg
         print *, "e_iso:", emis_ant(i, k, j, p_e_iso)*convg
         print *, "e_so4 f,c:", emis_ant(i, k, j, p_e_so4j)*conva
         print *, "e_no3 f,c:", emis_ant(i, k, j, p_e_no3j)*conva
         print *, "e_orgc:", emis_ant(i, k, j, p_e_orgc)*conva
         print *, "e_ecc:", emis_ant(i, k, j, p_e_ecc)*conva
         print*
      else
         print *, "dz8w=0 so cannot show adjusted emissions"
      end if
      print *, "CHEM_DBG PRINT (PPM or ug/m^3) AT (i,k,j):", i, k, j
      do n = 1, num_chem
         print *, n, chem(i, k, j, n)
      end do
      if (present(ph_macr)) then
         print *, "PHOTOLYSIS DATA:"
         print *, "ph_macr:", ph_macr(i, :, j)
         print *, "ph_o31d:", ph_o31d(i, :, j)
         print *, "ph_o33p:", ph_o33p(i, :, j)
         print *, "ph_no2:", ph_no2(i, :, j)
         print *, "ph_no3o2:", ph_no3o2(i, :, j)
         print *, "ph_no3o:", ph_no3o(i, :, j)
         print *, "ph_hno2:", ph_hno2(i, :, j)
         print *, "ph_hno3:", ph_hno3(i, :, j)
         print *, "ph_hno4:", ph_hno4(i, :, j)
         print *, "ph_h2o2:", ph_h2o2(i, :, j)
         print *, "ph_ch2or:", ph_ch2or(i, :, j)
         print *, "ph_ch2om:", ph_ch2om(i, :, j)
         print *, "ph_ch3cho:", ph_ch3cho(i, :, j)
         print *, "ph_ch3coch3:", ph_ch3coch3(i, :, j)
         print *, "ph_ch3coc2h5:", ph_ch3coc2h5(i, :, j)
         print *, "ph_hcocho:", ph_hcocho(i, :, j)
         print *, "ph_ch3cocho:", ph_ch3cocho(i, :, j)
         print *, "ph_hcochest:", ph_hcochest(i, :, j)
         print *, "ph_ch3o2h:", ph_ch3o2h(i, :, j)
         print *, "ph_ch3coo2h:", ph_ch3coo2h(i, :, j)
         print *, "ph_ch3ono2:", ph_ch3ono2(i, :, j)
         print *, "ph_hcochob:", ph_hcochob(i, :, j)
         print *, "ph_n2o5:", ph_n2o5(i, :, j)
         print *, "ph_o2:", ph_o2(i, :, j)
      end if
      print*
   end subroutine chem_dbg
#endif

end module module_input_chem_data
