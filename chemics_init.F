!WRF:MODEL_LAYER:INITIALIZATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Codename Pumpkin: Abstracted Bindings for Chemistry-to-WRF
!   (c) Haipeng Lin <linhaipeng@pku.edu.cn>, 2018-01
!   Developed at Peking University, School of Physics
!
!   This Chemical Interface (chem/) is written after comprehensive study of
!   the original chem_driver.f from WRF-Chem v3.6.1
!   which is (c) their respective authors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chem_init(id, chem, emis_ant, scalar, dt, bioemdt, photdt, chemdt, stepbioe, &
                    stepphot, stepchem, stepfirepl, plumerisefire_frq, z_at_w, xlat, xlong, &
                    g, aerwrf, config_flags, grid, alt, t, p, CONVFAC, ttday, tcosz, julday, gmt, &
                    tauaer1, tauaer2, tauaer3, tauaer4, &
                    gaer1, gaer2, gaer3, gaer4, &
                    waer1, waer2, waer3, waer4, &
                    l2aer, l3aer, l4aer, l5aer, l6aer, l7aer, &
                    extaerlw1, extaerlw2, extaerlw3, extaerlw4, &
                    extaerlw5, extaerlw6, extaerlw7, extaerlw8, &
                    extaerlw9, extaerlw10, extaerlw11, extaerlw12, &
                    extaerlw13, extaerlw14, extaerlw15, extaerlw16, &
                    tauaerlw1, tauaerlw2, tauaerlw3, tauaerlw4, &
                    tauaerlw5, tauaerlw6, tauaerlw7, tauaerlw8, &
                    tauaerlw9, tauaerlw10, tauaerlw11, tauaerlw12, &
                    tauaerlw13, tauaerlw14, tauaerlw15, tauaerlw16, &
                    dgnum4d, dgnumwet4d, dgnum_a1, dgnum_a2, dgnum_a3, &
                    dgnumwet_a1, dgnumwet_a2, dgnumwet_a3, &
                    pm2_5_dry, pm2_5_water, pm2_5_dry_ec, &
                    last_chem_time_year, last_chem_time_month, &
                    last_chem_time_day, last_chem_time_hour, &
                    last_chem_time_minute, last_chem_time_second, &
                    chem_in_opt, kemit, num_vert_mix, &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte)

      ! WRF Core Dependencies
      USE module_domain
      USE module_configure
      USE module_state_description
      USE module_dm

      ! WRF Physics Dependencies
      ! Just maintaining old code to not break CAM Physics Support.
      USE module_cam_support, only: gas_pcnst_modal_aero, gas_pcnst_modal_aero_pos

      USE module_input_chem_data, only: get_last_gas, last_chem_time, setup_gasprofile_maps

      implicit none

      ! provides information on WRF DMPAR parallel (MPI) process information
      ! namely whether this is the master process, using wrf_dm_on_monitor()
      logical, external :: wrf_dm_on_monitor

      real, intent(in) :: bioemdt, photdt, chemdt, dt, gmt
      integer, intent(IN) :: plumerisefire_frq
      integer, intent(IN) :: chem_in_opt
      integer, intent(INOUT) :: num_vert_mix
      integer, intent(IN) :: id, julday, kemit, &
                             last_chem_time_year, &
                             last_chem_time_month, &
                             last_chem_time_day, &
                             last_chem_time_hour, &
                             last_chem_time_minute, &
                             last_chem_time_second, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(INOUT) :: &
         pm2_5_dry, pm2_5_water, pm2_5_dry_ec, &
         tauaer1, tauaer2, tauaer3, tauaer4, &
         extaerlw1, extaerlw2, extaerlw3, extaerlw4, &
         extaerlw5, extaerlw6, extaerlw7, extaerlw8, &
         extaerlw9, extaerlw10, extaerlw11, extaerlw12, &
         extaerlw13, extaerlw14, extaerlw15, extaerlw16, &
         tauaerlw1, tauaerlw2, tauaerlw3, tauaerlw4, &
         tauaerlw5, tauaerlw6, tauaerlw7, tauaerlw8, &
         tauaerlw9, tauaerlw10, tauaerlw11, tauaerlw12, &
         tauaerlw13, tauaerlw14, tauaerlw15, tauaerlw16, &
         gaer1, gaer2, gaer3, gaer4, &
         waer1, waer2, waer3, waer4

      !-- Arrays needed to output dgnum when restarting
      real, dimension(ims:ime, kms:kme, jms:jme, 3), &
         intent(INOUT) :: dgnum4d, dgnumwet4d
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(INOUT) :: &
         dgnum_a1, dgnum_a2, dgnum_a3, &
         dgnumwet_a1, dgnumwet_a2, dgnumwet_a3
      !-- end dgnum restart arrays

      real, dimension(ims:ime, kms:kme, jms:jme, 1:4), &
         intent(INOUT) :: l2aer, l3aer, l4aer, l5aer, l6aer, l7aer
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(IN) :: z_at_w, t, p, alt, convfac
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), &
         intent(INOUT) :: chem
      real, dimension(ims:ime, 1:kemit, jms:jme, num_emis_ant), &
         intent(INOUT) :: emis_ant
      real, dimension(ims:ime, kms:kme, jms:jme, num_scalar), &
         intent(INOUT) :: scalar
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(INOUT) :: aerwrf
      real, dimension(ims:ime, jms:jme), &
         intent(INOUT) :: ttday, tcosz, xlat, xlong
      real, intent(IN) :: g
      integer, intent(out) :: stepbioe, stepphot, stepchem, stepfirepl
      type(grid_config_rec_type), intent(in) :: config_flags
      type(domain), intent(inout) :: grid

      ! Local Variables
      character*256 :: mminlu_loc
      character*256 :: message_txt
      type(WRFU_TimeInterval) :: tmpTimeInterval
      integer :: i, j, k, l, numgas, n, kk, nv

      ! From module_aerosols_sorgam / dep_simple, hplin 2-25-2018
      ! num_chem is from module_state_description
      logical, allocatable :: is_aerosol(:) ! true if field is aerosol (any phase)

      ! From module ctrans_grell, hplin 2-25-2018
      real :: epsilc

      ! WRF DM (MPI) Parallel Information - is master process?
      logical :: Am_I_Root
      integer :: WRF_DM_MyProc, WRF_DM_NProc

      if(wrf_dm_on_monitor()) then
         Am_I_Root = .true.
      else
         Am_I_Root = .false.
      endif

      call wrf_get_nproc(WRF_DM_NProc)
      call wrf_get_myproc(WRF_DM_MyProc)

#ifdef CHEM_DBG_I
      call print_chem_species_index(config_flags%chem_opt)
#endif

      call wrf_message("*****************************************************************")
      call wrf_message("  __          _______  ______       _____  _____ _    _ _____    ")
      call wrf_message("  \ \        / /  __ \|  ____|     / ____|/ ____| |  | |  __ \   ")
      call wrf_message("   \ \  /\  / /| |__) | |__ ______| |  __| |    | |__| | |__) |  ")
      call wrf_message("    \ \/  \/ / |  _  /|  __|______| | |_ | |    |  __  |  ___/   ")
      call wrf_message("     \  /\  /  | | \ \| |         | |__| | |____| |  | | |       ")
      call wrf_message("      \/  \/   |_|  \_\_|          \_____|\_____|_|  |_|_|       ")
      call wrf_message("*****************************************************************")
      call wrf_message("       WRF-GCHP MODEL, VERSION 1802.01 - WRF BINDINGS            ")
      call wrf_message("*****************************************************************")

      numgas = get_last_gas(config_flags%chem_opt)

      chem_select: SELECT CASE(config_flags%chem_opt)
      CASE (1)
        CALL wrf_debug(15, 'Pumpkin chemics_init: Initializing Chemistry Driver')
      END SELECT chem_select

      ! WRF-Chem used to perform configuration flags validation here.
      ! This can be done too in the future to prevent "gotchas"

      IF (config_flags%chem_opt == 0 .AND. config_flags%aer_ra_feedback .NE. 0) THEN
         call wrf_error_fatal("chemics_init: If no chemistry, aer_ra_feedback must be 0")
      ENDIF

      CALL nl_get_mminlu(1, mminlu_loc)

      IF (trim(mminlu_loc) /= 'USGS' .and. trim(mminlu_loc) /= 'MODIFIED_IGBP_MODIS_NOAH') THEN
         print *, mminlu_loc
         message_txt = "chemics_init: Chemistry routines require USGS or MODIS_NOAH land use maps."
         call wrf_error_fatal(trim(message_txt))
      ELSE
         IF (trim(mminlu_loc) == 'USGS' .and. grid%num_land_cat <= 23) THEN
            message_txt = "chemics_init: USGS land use map should have 24 or more categories."
            call wrf_error_fatal(trim(message_txt))
         ELSEIF (trim(mminlu_loc) == 'MODIFIED_IGBP_MODIS_NOAH' .and. grid%num_land_cat <= 19) THEN
            message_txt = "chemics_init: MODIS_NOAH land use map should have 20 or more categories."
            call wrf_error_fatal(trim(message_txt))
         ENDIF
      ENDIF

      !-- Load dgnum arrays when restart is active
      if (config_flags%restart) then
         do j = jts, jte
            do k = kts, kte
               do i = its, ite
                  dgnum4d(i, k, j, 1) = dgnum_a1(i, k, j)
                  dgnum4d(i, k, j, 2) = dgnum_a2(i, k, j)
                  dgnum4d(i, k, j, 3) = dgnum_a3(i, k, j)

                  dgnumwet4d(i, k, j, 1) = dgnumwet_a1(i, k, j)
                  dgnumwet4d(i, k, j, 2) = dgnumwet_a2(i, k, j)
                  dgnumwet4d(i, k, j, 3) = dgnumwet_a3(i, k, j)
               end do
            end do
         end do
      endif
      !-- end load dgnum arrays

      if (.NOT. config_flags%restart) then
         do j = jts, jte
            do k = kts, kte
               do i = its, ite
                  tauaer1(i, k, j) = 0.
                  tauaer2(i, k, j) = 0.
                  tauaer3(i, k, j) = 0.
                  tauaer4(i, k, j) = 0.
                  gaer1(i, k, j) = 0.
                  gaer2(i, k, j) = 0.
                  gaer3(i, k, j) = 0.
                  gaer4(i, k, j) = 0.
                  waer1(i, k, j) = 0.
                  waer2(i, k, j) = 0.
                  waer3(i, k, j) = 0.
                  waer4(i, k, j) = 0.
                  l2aer(i, k, j, 1) = 0.
                  l2aer(i, k, j, 2) = 0.
                  l2aer(i, k, j, 3) = 0.
                  l2aer(i, k, j, 4) = 0.
                  l3aer(i, k, j, 1) = 0.
                  l3aer(i, k, j, 2) = 0.
                  l3aer(i, k, j, 3) = 0.
                  l3aer(i, k, j, 4) = 0.
                  l4aer(i, k, j, 1) = 0.
                  l4aer(i, k, j, 2) = 0.
                  l4aer(i, k, j, 3) = 0.
                  l4aer(i, k, j, 4) = 0.
                  l5aer(i, k, j, 1) = 0.
                  l5aer(i, k, j, 2) = 0.
                  l5aer(i, k, j, 3) = 0.
                  l5aer(i, k, j, 4) = 0.
                  l6aer(i, k, j, 1) = 0.
                  l6aer(i, k, j, 2) = 0.
                  l6aer(i, k, j, 3) = 0.
                  l6aer(i, k, j, 4) = 0.
                  l7aer(i, k, j, 1) = 0.
                  l7aer(i, k, j, 2) = 0.
                  l7aer(i, k, j, 3) = 0.
                  l7aer(i, k, j, 4) = 0.
                  extaerlw1(i, k, j) = 0.
                  extaerlw2(i, k, j) = 0.
                  extaerlw3(i, k, j) = 0.
                  extaerlw4(i, k, j) = 0.
                  extaerlw5(i, k, j) = 0.
                  extaerlw6(i, k, j) = 0.
                  extaerlw7(i, k, j) = 0.
                  extaerlw8(i, k, j) = 0.
                  extaerlw9(i, k, j) = 0.
                  extaerlw10(i, k, j) = 0.
                  extaerlw11(i, k, j) = 0.
                  extaerlw12(i, k, j) = 0.
                  extaerlw13(i, k, j) = 0.
                  extaerlw14(i, k, j) = 0.
                  extaerlw15(i, k, j) = 0.
                  extaerlw16(i, k, j) = 0.
                  tauaerlw1(i, k, j) = 0.
                  tauaerlw2(i, k, j) = 0.
                  tauaerlw3(i, k, j) = 0.
                  tauaerlw4(i, k, j) = 0.
                  tauaerlw5(i, k, j) = 0.
                  tauaerlw6(i, k, j) = 0.
                  tauaerlw7(i, k, j) = 0.
                  tauaerlw8(i, k, j) = 0.
                  tauaerlw9(i, k, j) = 0.
                  tauaerlw10(i, k, j) = 0.
                  tauaerlw11(i, k, j) = 0.
                  tauaerlw12(i, k, j) = 0.
                  tauaerlw13(i, k, j) = 0.
                  tauaerlw14(i, k, j) = 0.
                  tauaerlw15(i, k, j) = 0.
                  tauaerlw16(i, k, j) = 0.
               end do
            end do
         end do

         ! Zero out anthro. emission arrays.
         do l = 1, num_emis_ant
           do j = jts, jte
              do k = 1, kemit
                 do i = its, ite
                    emis_ant(i, k, j, l) = 0.
                 end do
              end do
           end do
         end do
      end if

      ! Catch some unsupported options.
      ! In early development, 
      if (config_flags%progn > 0) then
         call wrf_error_fatal("Pumpkin chemics_init: progn > 0 option is not supported by this Chemistry abstraction layer.")
      endif

      if (config_flags%have_bcs_upper) then
         call wrf_error_fatal("Pumpkin chemics_init: have_bcs_upper = .true. is not supported by this Chemistry abstraction layer.")
      endif

      ! Don't run if there isn't anything to do...
      if(config_flags%chem_opt == 0) then
        CALL wrf_debug(15, 'Pumpkin chemics_init: **** NO CHEMISTRY IS USED IN THIS WRF RUN ****')
        return
      endif

      ! Set the num_vert_mix variable if using ACM
      num_vert_mix = 0
      IF (config_flags%bl_pbl_physics == ACMPBLSCHEME) THEN
         mix_select: select case(config_flags%chem_opt)
            case (radm2sorg_aq, radm2sorg_aqchem, racmsorg_aq)
               ! if you have aqueous chemistry, then the num_vert_mix needs to be # of gas species & not # of chemistry_species.
               num_vert_mix = numgas
            case default
               num_vert_mix = num_chem
         end select mix_select

         if (num_vert_mix .gt. config_flags%ndepvel) then
            write (message_txt, '(A30,2(I8,2x))') 'chem_init: num_vert_mix and ndepvel ', num_vert_mix, config_flags%ndepvel
            call wrf_message(trim(message_txt))
            call wrf_error_fatal("Pumpkin chemics_init: FATAL - num_vert_mix > ndepvel ")
         endif
      endif

      ! Calculate "timesteps" for out intent.
      stepbioe = nint(bioemdt*60./dt)
      stepphot = nint(photdt*60./dt)
      stepchem = nint(chemdt*60./dt)
      stepfirepl = nint(plumerisefire_frq*60/dt)
      stepbioe = max(stepbioe, 1)
      stepphot = max(stepphot, 1)
      stepchem = max(stepchem, 1)
      stepfirepl = max(stepfirepl, 1)
      call wrf_debug(15, 'Pumpkin chemics_init: Starting chemistry init')

      ! Setup the cross reference mappings between the default profiles and
      ! the gas mechanism species (wig, 2-May-2007)
      ! Sets up the cross reference mapping indices and fractional
      ! apportionment of the default species profiles for use with
      ! ICs and BCs. (module_input_chem_data)
      call setup_gasprofile_maps(config_flags%chem_opt, numgas)

      ! Aerosol Information (is_aerosol array) allocation.
      if (.not. allocated(is_aerosol)) then
         allocate (is_aerosol(num_chem))
      else
         if (size(is_aerosol) /= num_chem) &
            call wrf_error_fatal("The number of chemistry species has changed between nests. Use same chem_opt for all domains.")
      end if

      aer_select: select case(config_flags%chem_opt)
         case(1)
            call wrf_debug(15, 'chemics_init: stubbed out: aerosols initialization routine')
            if (config_flags%chem_in_opt == 0) then
               if (.not. config_flags%restart) then
                  do j = jts, jte
                     do k = kts, kte
                        do i = its, ite
                        !    chem(i, k, j, p_dms) = ...
                        enddo
                     enddo
                  enddo
               endif
            endif

            !...convert aerosols to mixing ratio
            if (.not. config_flags%restart) then
              if (config_flags%chem_in_opt == 0 .and. num_chem .gt. numgas) then
              do l = numgas + 1, num_chem
                 do j = jts, jte
                    do k = kts, kte
                       kk = min(k, kde - 1)
                       do i = its, ite
                          chem(i, k, j, l) = chem(i, kk, j, l)*alt(i, kk, j)
                       enddo
                    enddo
                 enddo
              enddo
              endif
            endif

            ! ensure chem has a min value >= epsilc
            epsilc=1.e-30
            chem(its:ite, kts:min(kte, kde - 1), jts:jte, :) = max(chem(its:ite, kts:min(kte, kde - 1), jts:jte, :), epsilc)
      end select aer_select

      do nv = 1, num_chem
         do j = jts, jte
            do i = its, ite
               chem(i, kde, j, nv) = chem(i, kde - 1, j, nv)
            enddo
         enddo
      enddo

      ! Setup time info for chemistry
      if (config_flags%restart) then
         call wrf_debug(15, "Setting last_chem_time from restart file")
         call WRFU_TimeSet(last_chem_time(id), &
                           YY=last_chem_time_year, &
                           MM=last_chem_time_month, &
                           DD=last_chem_time_day, &
                           H=last_chem_time_hour, &
                           M=last_chem_time_minute, &
                           S=last_chem_time_second)
      else
         call wrf_debug(15, "Setting last_chem_time to model start time-dt")
         call WRFU_TimeIntervalSet(tmpTimeInterval, s_=real(dt, 8))
         last_chem_time(id) = domain_get_current_time(grid) - tmpTimeInterval
      end if

      ! Keep the following for compatibility with WRF CAM/MAM Physics Scheme.
      gas_pcnst_modal_aero_pos = max(1, gas_pcnst_modal_aero) ! BSINGH - Fix for non-MAM simulations
END SUBROUTINE chem_init

!-----------------------------------------------------------------------
!  CHEMISTRY DEBUG: print_chem_species_index(chem_opt)
!  Outputs Index Numbers for Chemistry species
!-----------------------------------------------------------------------

#ifdef CHEM_DBG_I
   subroutine print_chem_species_index(chem_opt)
      use module_state_description
      implicit none

      integer, intent(in) :: chem_opt

      print *, 'Index numbers for chemistry species:'
!
! Gas species...
!
      select case (chem_opt)
      case (1)
         print *, p_so2, "so2"
         print *, p_sulf, "sulf"
         print *, p_no2, "no2"
         print *, p_no, "no"
         print *, p_o3, "o3"
         print *, p_hno3, "hno3"
         print *, p_h2o2, "h2o2"
         print *, p_ald, "ald"
         print *, p_hcho, "hcho"
         print *, p_op1, "op1"
         print *, p_op2, "op2"
         print *, p_paa, "paa"
         print *, p_ora1, "ora1"
         print *, p_ora2, "ora2"
         print *, p_nh3, "nh3"
         print *, p_n2o5, "n2o5"
         print *, p_no3, "no3"
         print *, p_pan, "pan"
         print *, p_hc3, "hc3"
         print *, p_hc5, "hc5"
         print *, p_hc8, "hc8"
         print *, p_eth, "eth"
         print *, p_co, "co"
         print *, p_ol2, "ol2"
         print *, p_olt, "olt"
         print *, p_oli, "oli"
         print *, p_tol, "tol"
         print *, p_xyl, "xyl"
         print *, p_aco3, "aco3"
         print *, p_tpan, "tpan"
         print *, p_hono, "hono"
         print *, p_hno4, "hno4"
         print *, p_ket, "ket"
         print *, p_gly, "gly"
         print *, p_mgly, "mgly"
         print *, p_dcb, "dcb"
         print *, p_onit, "onit"
         print *, p_csl, "csl"
         print *, p_iso, "iso"
         print *, p_ho, "ho"
         print *, p_ho2, "ho2"
      end select

!
! Aerosol species...
!
      select case (chem_opt)
      case (1)
         print *, p_so4aj, "so4aj"
         print *, p_so4ai, "so4ai"
         print *, p_nh4aj, "nh4aj"
         print *, p_nh4ai, "nh4ai"
         print *, p_no3aj, "no3aj"
         print *, p_no3ai, "no3ai"
         print *, p_orgaro1j, "orgaro1j"
         print *, p_orgaro1i, "orgaro1i"
         print *, p_orgaro2j, "orgaro2j"
         print *, p_orgaro2i, "orgaro2i"
         print *, p_orgalk1j, "orgalk1j"
         print *, p_orgalk1i, "orgalk1i"
         print *, p_orgole1j, "orgole1j"
         print *, p_orgole1i, "orgole1i"
         print *, p_orgba1j, "orgba1j"
         print *, p_orgba1i, "orgba1i"
         print *, p_orgba2j, "orgba2j"
         print *, p_orgba2i, "orgba2i"
         print *, p_orgba3j, "orgba3j"
         print *, p_orgba3i, "orgba3i"
         print *, p_orgba4j, "orgba4j"
         print *, p_orgba4i, "orgba4i"
         print *, p_orgpaj, "orgpaj"
         print *, p_orgpai, "orgpai"
         print *, p_ecj, "ecj"
         print *, p_eci, "eci"
         print *, p_p25j, "p25j"
         print *, p_p25i, "p25i"
         print *, p_antha, "antha"
         print *, p_seas, "seas"
         print *, p_soila, "soila"
         print *, p_nu0, "nu0"
         print *, p_ac0, "ac0"
         print *, p_corn, "corn"
      end select
   end subroutine print_chem_species_index
#endif
