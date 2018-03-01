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

   ! input_chem_profile
   ! req. by real_em
   ! Based on original WRF-Chem module_input_chem_data::input_chem_profile, (c) original authors
   subroutine input_chem_profile(si_grid)
      implicit none
      type(domain) :: si_grid

      INTEGER :: i, j, k, &
                 ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 ips, ipe, jps, jpe, kps, kpe
      INTEGER :: fid, ierr, numgas
      INTEGER :: debug_level

      REAL, ALLOCATABLE, DIMENSION(:, :, :) :: si_zsigf, si_zsig

      ! Get grid dimensions
      call get_ijk_from_grid(si_grid, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             ips, ipe, jps, jpe, kps, kpe)

      ! Get scalar grid point heights
      allocate (si_zsigf(ims:ime, kms:kme, jms:jme))
      allocate (si_zsig(ims:ime, kms:kme, jms:jme))

      si_zsigf = (si_grid%ph_1 + si_grid%phb)/grav

      do k = 1, kde - 1
         si_zsig(:, k, :) = 0.5*(si_zsigf(:, k, :) + si_zsigf(:, k + 1, :))
      enddo
      si_zsig(:, kde, :) = 0.5*(3.*si_zsigf(:, kde, :) - si_zsigf(:, kde - 1, :))

      ! Determine the index of the last gas species
      numgas = get_last_gas(si_grid%chem_opt)

      ! Setup the cross reference mappings between the default profiles and
      ! the gas mechanism species (wig, 2-May-2007)
      call setup_gasprofile_maps(si_grid%chem_opt, numgas)

      ! Interpolate the chemistry data to the SI grid. These values should typically
      ! be set to match the values in bdy_chem_value_tracer so that the boundaries
      ! and interior match each other.
      if (si_grid%chem_opt == 0) then
         ! this is set as an example from the original module_input_chem_data
         ! si_grid%chem(ims:ime, kms:kme, jms:jme, 1:numgas) = 0.
      else
         call make_chem_profile(ims, ime, jms, jme, kms, kme, num_chem, numgas, &
                                si_grid%chem_opt, si_zsig, si_grid%chem)
      end if

      call wrf_debug(100, ' input_chem_profile: exit subroutine ')

      deallocate (si_zsigf); deallocate (si_zsig)
      return
   end subroutine input_chem_profile

   ! make_chem_profile
   ! req. by input_chem_profile, in turn required by real_em
   ! Based on original WRF-Chem module_input_chem_data::input_chem_profile, (c) original authors
   subroutine make_chem_profile(nx1, nx2, ny1, ny2, nz1, nz2, nch, numgas, &
                                chem_opt, zgrid, chem)
      implicit none

      integer, intent(in) :: nx1, ny1, nz1, nx2, ny2, nz2
      integer, intent(in) :: nch, numgas, chem_opt
      real, dimension(nx1:nx2, nz1:nz2, ny1:ny2) :: zgrid

      integer :: i, j, k, l, is
      real, dimension(nx1:nx2, nz1:kx, ny1:ny2, lo + 1) :: chprof
      real, dimension(nx1:nx2, nz1:kx, ny1:ny2) :: zprof

      real, dimension(nx1:nx2, nz1:nz2, ny1:ny2, nch) :: chem
      real, dimension(nx1:nx2, nz1:nz2, ny1:ny2, lo) :: stor

      if (nch .NE. num_chem) then
         call wrf_error_fatal("Pumpkin module_input_chem_data: wrong number of chemical species.")
      endif

      ! Vertically flip the chemistry data as it is given top down and
      ! heights are bottom up. Fill temp 3D chemical and profile array,
      ! keep chem slot 1 open as vinterp_chem assumes there is no data.
      do j = ny1, ny2
         do k = 1, kx
            do i = nx1, nx2
               chprof(i, k, j, 2:lo + 1) = xl(1:lo, kx - k + 1)
               zprof(i, k, j) = 0.5*(zfa(k) + zfa(k + 1))
            enddo
         enddo
      enddo

      ! Change number concentrations to mixing ratios for short-lived NALROM species
      do k = 1, kx
         chprof(:, k, :, lo - 5:lo + 1) = chprof(:, k, :, lo - 5:lo + 1)/dens(k)
      end do

      ! Interpolate temp 3D chemical and profile array to WRF grid
      call vinterp_chem(nx1, nx2, ny1, ny2, nz1, kx, nz2, lo, zprof, zgrid, &
                        chprof, chem, .false.)

      ! place interpolated data into temp storage array
      stor(nx1:nx2, nz1:nz2, ny1:ny2, 1:lo) = chem(nx1:nx2, nz1:nz2, ny1:ny2, 2:lo + 1)

      ! Here is where the chemistry profile is constructed
      chem(nx1:nx2, nz1:nz2, ny1:ny2, 1) = -999.

      do l = 2, numgas
         is = iref(l - 1)
         do j = ny1, ny2
            do k = nz1, nz2
               do i = nx1, nx2
                  chem(i, k, j, l) = fracref(l - 1) * stor(i, k, j, is) * 1.E6
               enddo
            enddo
         enddo
      enddo

      return
   end subroutine make_chem_profile

   ! vinterp_chem
   ! Interpolates columns of chemistry data from one set of height surfaces to another.
   ! req. by make_chem_profile
   ! WRF-Chem v3.6.1, (c) original authors
   subroutine vinterp_chem(nx1, nx2, ny1, ny2, nz1, nz_in, nz_out, nch, z_in, z_out, &
                           data_in, data_out, extrapolate)
      INTEGER, INTENT(IN)                :: nx1, nx2
      INTEGER, INTENT(IN)                :: ny1, ny2
      INTEGER, INTENT(IN)                :: nz1
      INTEGER, INTENT(IN)                :: nz_in
      INTEGER, INTENT(IN)                :: nz_out
      INTEGER, INTENT(IN)                :: nch
      REAL, INTENT(IN)                   :: z_in(nx1:nx2, nz1:nz_in, ny1:ny2)
      REAL, INTENT(IN)                   :: z_out(nx1:nx2, nz1:nz_out, ny1:ny2)
      REAL, INTENT(IN)                   :: data_in(nx1:nx2, nz1:nz_in, ny1:ny2, nch)
      REAL, INTENT(OUT)                  :: data_out(nx1:nx2, nz1:nz_out, ny1:ny2, nch)
      LOGICAL, INTENT(IN)                :: extrapolate

      INTEGER                            :: i, j, l
      INTEGER                            :: k, kk
      REAL                               :: desired_z
      REAL                               :: dvaldz
      REAL                               :: wgt0

      ! Loop over the number of chemical species
      chem_loop: DO l = 2, nch
         data_out(:, :, :, l) = -99999.9

         DO j = ny1, ny2
            DO i = nx1, nx2
               output_loop: DO k = nz1, nz_out
                  desired_z = z_out(i, k, j)
                  IF (desired_z .LT. z_in(i, 1, j)) THEN
                     IF ((desired_z - z_in(i, 1, j)) .LT. 0.0001) THEN
                        data_out(i, k, j, l) = data_in(i, 1, j, l)
                     ELSE
                        IF (extrapolate) THEN
                           ! Extrapolate downward because desired height level is below
                           ! the lowest level in our input data.  Extrapolate using simple
                           ! 1st derivative of value with respect to height for the bottom 2
                           ! input layers.

                           ! Add a check to make sure we are not using the gradient of
                           ! a very thin layer
                           IF ((z_in(i, 1, j) - z_in(i, 2, j)) .GT. 0.001) THEN
                              dvaldz = (data_in(i, 1, j, l) - data_in(i, 2, j, l))/ &
                                       (z_in(i, 1, j) - z_in(i, 2, j))
                           ELSE
                              dvaldz = (data_in(i, 1, j, l) - data_in(i, 3, j, l))/ &
                                       (z_in(i, 1, j) - z_in(i, 3, j))
                           ENDIF
                           data_out(i, k, j, l) = MAX(data_in(i, 1, j, l) + &
                                                      dvaldz*(desired_z - z_in(i, 1, j)), 0.)
                        ELSE
                           data_out(i, k, j, l) = data_in(i, 1, j, l)
                        ENDIF
                     ENDIF
                  ELSE IF (desired_z .GT. z_in(i, nz_in, j)) THEN
                     IF ((z_in(i, nz_in, j) - desired_z) .LT. 0.0001) THEN
                        data_out(i, k, j, l) = data_in(i, nz_in, j, l)
                     ELSE
                        IF (extrapolate) THEN
                           ! Extrapolate upward
                           IF ((z_in(i, nz_in - 1, j) - z_in(i, nz_in, j)) .GT. 0.0005) THEN
                              dvaldz = (data_in(i, nz_in, j, l) - data_in(i, nz_in - 1, j, l))/ &
                                       (z_in(i, nz_in, j) - z_in(i, nz_in - 1, j))
                           ELSE
                              dvaldz = (data_in(i, nz_in, j, l) - data_in(i, nz_in - 2, j, l))/ &
                                       (z_in(i, nz_in, j) - z_in(i, nz_in - 2, j))
                           ENDIF
                           data_out(i, k, j, l) = MAX(data_in(i, nz_in, j, l) + &
                                                      dvaldz*(desired_z - z_in(i, nz_in, j)), 0.)
                        ELSE
                           data_out(i, k, j, l) = data_in(i, nz_in, j, l)
                        ENDIF
                     ENDIF
                  ELSE
                     ! We can trap between two levels and linearly interpolate
                     input_loop: DO kk = 1, nz_in - 1
                        IF (desired_z .EQ. z_in(i, kk, j)) THEN
                           data_out(i, k, j, l) = data_in(i, kk, j, l)
                           EXIT input_loop
                        ELSE IF (desired_z .EQ. z_in(i, kk + 1, j)) THEN
                           data_out(i, k, j, l) = data_in(i, kk + 1, j, l)
                           EXIT input_loop
                        ELSE IF ((desired_z .GT. z_in(i, kk, j)) .AND. &
                                 (desired_z .LT. z_in(i, kk + 1, j))) THEN
                           wgt0 = (desired_z - z_in(i, kk + 1, j))/ &
                                  (z_in(i, kk, j) - z_in(i, kk + 1, j))
                           data_out(i, k, j, l) = MAX(wgt0*data_in(i, kk, j, l) + &
                                                      (1.-wgt0)*data_in(i, kk + 1, j, l), 0.)
                           EXIT input_loop
                        ENDIF
                     ENDDO input_loop
                  ENDIF
               ENDDO output_loop
            ENDDO
         ENDDO
      ENDDO chem_loop

      return
   end subroutine vinterp_chem

   ! chem_dbg
   ! Chemistry Debug function req. by chem/chem_driver
   !!!! CHEMISTRY DEVELOPERS: YOU WILL LIKELY NEED TO UPDATE THIS SPECIES LIST FOR DEBUGGING.
   !!!! THIS IS SEPARATE FROM REGISTRY.CHEM
#ifdef CHEM_DBG_I
   subroutine chem_dbg(i, j, k, dtstep, itimestep, &
                       dz8w, t_phy, p_phy, rho_phy, chem, emis_ant, &
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
