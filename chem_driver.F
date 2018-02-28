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

subroutine chem_driver(grid, config_flags &
#include "dummy_new_args.inc"
)

   !! WRF Interfaces
   USE module_domain, only:domain
   USE module_configure
   USE module_driver_constants
   USE module_machine
   USE module_tiles
   USE module_dm
   USE module_model_constants
   USE module_state_description

   !! Chemistry Stub Interfaces
   USE module_chem_utilities
   USE module_aer_opt_out, only:aer_opt_out
   USE module_input_chem_data, only:last_chem_time, &
#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
                                    chem_dbg, &
#endif
      get_last_gas

   USE module_upper_bc_driver, only:upper_bc_driver
   USE module_tropopause, only:tropopause_driver

   IMPLICIT NONE

   ! BSINGH(PNNL)- Lahey compiler forces to declare the following interface
   interface
      SUBROUTINE sum_pm_driver(config_flags, &
                               alt, chem, h2oaj, h2oai, &
                               pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, &
                               hoa_a01, hoa_a02, hoa_a03, hoa_a04, &
                               bboa_a01, bboa_a02, bboa_a03, bboa_a04, &
                               soa_a01, soa_a02, soa_a03, soa_a04, &
                               bbsoa_a01, bbsoa_a02, bbsoa_a03, bbsoa_a04, &
                               hsoa_a01, hsoa_a02, hsoa_a03, hsoa_a04, &
                               biog_a01, biog_a02, biog_a03, biog_a04, &
                               asmpsoa_a01, asmpsoa_a02, asmpsoa_a03, asmpsoa_a04, &
                               arosoa_a01, arosoa_a02, arosoa_a03, arosoa_a04, &
                               totoa_a01, totoa_a02, totoa_a03, totoa_a04, &
                               hsoa_c, hsoa_o, bbsoa_c, bbsoa_o, &
                               biog_v1, biog_v2, biog_v3, biog_v4, &
                               ant_v1, ant_v2, ant_v3, ant_v4, &
                               smpa_v1, smpbb_v1, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte)

         USE module_configure

         IMPLICIT NONE

         integer, INTENT(IN)    :: &
            ids, ide, jds, jde, kds, kde, &
            ims, ime, jms, jme, kms, kme, &
            its, ite, jts, jte, kts, kte

         real, dimension(ims:ime, kms:kme, jms:jme, num_chem), &
            INTENT(IN) :: chem

         real, dimension(ims:ime, kms:kme, jms:jme), &
            INTENT(IN) :: alt
         real, dimension(ims:ime, kms:kme, jms:jme), &
            OPTIONAL, &
            INTENT(IN) :: h2oaj, h2oai

         real, dimension(ims:ime, kms:kme, jms:jme), &
            OPTIONAL, &
            INTENT(OUT) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, &
                           hoa_a01, hoa_a02, hoa_a03, hoa_a04, &
                           bboa_a01, bboa_a02, bboa_a03, bboa_a04, &
                           soa_a01, soa_a02, soa_a03, soa_a04, &
                           bbsoa_a01, bbsoa_a02, bbsoa_a03, bbsoa_a04, &
                           hsoa_a01, hsoa_a02, hsoa_a03, hsoa_a04, &
                           biog_a01, biog_a02, biog_a03, biog_a04, &
                           arosoa_a01, arosoa_a02, arosoa_a03, arosoa_a04, &
                           totoa_a01, totoa_a02, totoa_a03, totoa_a04, &
                           hsoa_c, hsoa_o, bbsoa_c, bbsoa_o, &
                           biog_v1, biog_v2, biog_v3, biog_v4, &
                           ant_v1, ant_v2, ant_v3, ant_v4, &
                           smpa_v1, &
                           smpbb_v1, &
                           asmpsoa_a01, asmpsoa_a02, asmpsoa_a03, asmpsoa_a04

         TYPE(grid_config_rec_type), INTENT(IN)    :: config_flags

      end SUBROUTINE sum_pm_driver
   end interface

   !  Input data.
   TYPE(domain), TARGET :: grid

   !
   !  Definitions of dummy arguments to solve
# include <dummy_new_decl.inc>
# define NO_I1_OLD

   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags

   integer :: ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, &
              ips, ipe, jps, jpe, kps, kpe, &
              its, ite, jts, jte, kts, kte

! ..
! .. Local Scalars ..
   integer :: stepave, i, j, k, l, numgas, nv, n, nr, ktau, k_start, k_end, idf, jdf, kdf
   integer :: ijulian
! UoC dust scheme option
   integer :: imod

   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) ::vcsulf_old, vcso2_old, vch2o2_old

   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) ::n2o5_het
   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: &
      p_phy, u_phy, v_phy &
      , t_phy, dz8w, t8w, p8w &
      , rho, rri, z_at_w, vvel, zmid, rh
   real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: pbl_h
   real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33, 5) :: seasin, dustin
   real, dimension(grid%sm32:grid%em32-1) :: QL, TL
   real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: REXNSFC, FACTRS, TOT, TSFC

   ! temporary arrays for old chemistry values
   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33, num_chem_ct) :: chem_old
   integer, dimension(num_chem_ct) :: chem_ct_indices

! Variables for adaptive time steps...
   TYPE(WRFU_TimeInterval) :: tmpTimeInterval
   real(KIND=8) :: curr_secs
   real(KIND=8) :: real_time_r8 !ext. function in adapt_timestep_em.F
   logical      :: adapt_step_flag, do_chemstep, do_photstep

   real :: DAYI, DPL, FICE, FRAIN, HOUR, PLYR, QI, QR, QW, RADT, TIMES, WC, TDUM, WMSK, RWMSK
   integer                         :: ij
   integer                         :: im, num_3d_m, ic, num_3d_c, num_3d_s
   integer                         :: ijds, ijde
   integer                         :: ksubt

   real :: epsilc
   real :: chem_minval, dtstepc

   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: &
      del_h2so4_gasprod ! change to h2so4 during gas-phase chemistry, cam_mam packages

   logical :: haveaer
   CHARACTER(LEN=1000) :: msg
   CHARACTER(LEN=256) :: current_date_char
   integer :: current_month

   INTRINSIC max, min

!
! Setup the adaptive timestep for the chem routines. Most of this follows
! what is in solve_em, except for the call to adjust time_step.
!
   ! The necesssary variables exist for the EM core and using the adaptive
   ! techniques will work even with a constant time step. In fact, they
   ! prevent issues with restarts and changed time steps. So, we will
   ! always use them with the EM core.
   adapt_step_flag = .TRUE.
   ktau = grid%itimestep
   tmpTimeInterval = domain_get_time_since_sim_start(grid)
   curr_secs = real_time_r8(tmpTimeInterval)
   ijulian = ifix(grid%julian)

   ! initializing diagnostics and macros

   if (ktau == 1) then
      grid%conv_ct(:, :, :, :) = 0.
      grid%chem_ct(:, :, :, :) = 0.
      grid%vmix_ct(:, :, :, :) = 0.
   endif

   if (config_flags%chemdiag == USECHEMDIAG) then
      ! modify tendency list here
      chem_ct_indices(p_chem_co) = p_co
      chem_ct_indices(p_chem_o3) = p_o3
      chem_ct_indices(p_chem_no) = p_no
      chem_ct_indices(p_chem_no2) = p_no2
      chem_ct_indices(p_chem_hno3) = p_hno3
      chem_ct_indices(p_chem_iso) = p_iso
      chem_ct_indices(p_chem_ho) = p_ho
      chem_ct_indices(p_chem_ho2) = p_ho2
   endif

   !! -- Check whether to do chemistry steps (do_chemstep bool) -- !!
   do_chemstep = .false.
   IF (ktau == 1) then
      do_chemstep = .true.
      grid%ktauc = 1
   ELSE IF (adapt_step_flag) THEN
      IF ((grid%chemdt <= 0) .or. &
          (curr_secs + real(grid%dt, 8) + 0.01 >= &
           (INT(curr_secs/real(grid%chemdt*60., 8) + 1, 8)*real(grid%chemdt*60., 8))) &
          ) then
         do_chemstep = .true.
         grid%ktauc = grid%ktauc + 1
         last_chem_time(grid%id) = domain_get_current_time(grid)
         call WRFU_TimeGet(last_chem_time(grid%id), &
                           YY=grid%last_chem_time_year, &
                           MM=grid%last_chem_time_month, &
                           DD=grid%last_chem_time_day, &
                           H=grid%last_chem_time_hour, &
                           M=grid%last_chem_time_minute, &
                           S=grid%last_chem_time_second)
      ENDIF
   ELSE IF ((MOD(ktau, grid%stepchem) == 0) .or. (grid%stepchem == 1)) THEN
      do_chemstep = .true.
      grid%ktauc = max(ktau/grid%stepchem, 1)
   ENDIF

   CALL get_ijk_from_grid(grid, &
                          ids, ide, jds, jde, kds, kde, &
                          ims, ime, jms, jme, kms, kme, &
                          ips, ipe, jps, jpe, kps, kpe)

! following two lines needed for MEGAN
   CALL domain_clock_get(grid, current_timestr=current_date_char)
   read (current_date_char(6:7), FMT='(I2)') current_month

!
! initialize
!
   grid%raincv_b(:, :) = grid%raincv(:, :)

   num_3d_m = num_moist
   num_3d_c = num_chem
   num_3d_s = num_scalar
   numgas = get_last_gas(config_flags%chem_opt)

! Compute these starting and stopping locations for each tile and number of tiles.
   CALL set_tiles(grid, ids, ide, jds, jde, ips, ipe, jps, jpe)
   k_start = kps
   k_end = kpe

   ijds = min(ids, jds)
   ijde = max(ide, jde)

   ! From module_ctrans_grell (175)
   epsilc=1.e-30

   chem_minval = epsilc !chem_minval can be case dependant and set below...
   chem_select: SELECT CASE(config_flags%chem_opt)
CASE (1)
   CALL wrf_debug(15, 'Pumpkin chem_driver: chem_opt = 1, initializing chemistry option')
   haveaer = .false.
CASE DEFAULT
   call wrf_error_fatal("Pumpkin chem_driver: Unrecognized chem_opt")
   END SELECT chem_select

   do nv = 1, num_chem
      do j = jps, jpe
         do k = kps, kpe
            do i = ips, ipe
               chem(i, k, j, nv) = max(chem(i, k, j, nv), chem_minval) ! prevent underflow in chemistry
            enddo
         enddo
      enddo
   enddo

   do j = jps, min(jde - 1, jpe)
      do k = kps, kpe
         do i = ips, min(ide - 1, ipe)
            vvel(i, k, j) = grid%w_2(i, k, j)
            zmid(i, k, j) = grid%z(i, k, j)
         enddo
      enddo
   enddo
   do j = jps, min(jde - 1, jpe)
      do k = kps, min(kde - 1, kpe)
         do i = ips, min(ide - 1, ipe)
            rri(i, k, j) = grid%alt(i, k, j)
         enddo
      enddo
   enddo
   do j = jps, min(jde - 1, jpe)
      do i = ips, min(ide - 1, ipe)
         pbl_h(i, j) = grid%pblh(i, j)
      enddo
   enddo

!------------------------------------------------------------------------
! Main chemistry tile loop
!------------------------------------------------------------------------

!$OMP PARALLEL DO   &
!$OMP PRIVATE ( ij, its, ite, jts, jte )

!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
   chem_tile_loop_1: DO ij = 1, grid%num_tiles
      its = grid%i_start(ij)
      ite = min(grid%i_end(ij), ide - 1)
      jts = grid%j_start(ij)
      jte = min(grid%j_end(ij), jde - 1)

      kts = k_start
      kte = min(k_end, kde - 1)

      CALL wrf_debug(15, ' call chem_prep')

      !!! chem_prep @ module_chem_utilities.f !!!
      CALL chem_prep(config_flags, &
                     grid%u_2, grid%v_2, grid%p, grid%pb, &
                     grid%alt, grid%ph_2, grid%phb, grid%t_2, &
                     moist, num_3d_m, rho, &
                     p_phy, u_phy, v_phy, &
                     p8w, t_phy, t8w, grid%z, z_at_w, &
                     dz8w, rh, grid%fnm, grid%fnp, &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, &
                     k_start, k_end)

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
      if ((its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and. &
          (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and. &
          (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)) then
         call wrf_debug(15, "calling chem_dbg at top of chem_driver")
         call chem_dbg(CHEM_DBG_I, CHEM_DBG_J, CHEM_DBG_K, grid%dt, ktau, &
                       dz8w, t_phy, p_phy, rho, chem, emis_ant, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       config_flags%kemit)
      end if
#endif

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
      if ((its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and. &
          (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and. &
          (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)) then
         call wrf_debug(15, 'calling chem_dbg after emissions_driver')
         call chem_dbg(CHEM_DBG_I, CHEM_DBG_J, CHEM_DBG_K, grid%dt, ktau, &
                       dz8w, t_phy, p_phy, rho, chem, emis_ant, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       config_flags%kemit)
      end if
#endif

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
      if ((its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and. &
          (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and. &
          (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)) then
         call wrf_debug(15, 'calling chem_dbg after photolysis_driver')
         call chem_dbg(CHEM_DBG_I, CHEM_DBG_J, CHEM_DBG_K, grid%dt, ktau, &
                       dz8w, t_phy, p_phy, rho, chem, emis_ant, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       config_flags%kemit, &
                       grid%ph_macr, grid%ph_o31d, grid%ph_o33p, grid%ph_no2, grid%ph_no3o2, &
                       grid%ph_no3o, grid%ph_hno2, grid%ph_hno3, grid%ph_hno4, grid%ph_h2o2, &
                       grid%ph_ch2or, grid%ph_ch2om, grid%ph_ch3cho, grid%ph_ch3coch3, &
                       grid%ph_ch3coc2h5, grid%ph_hcocho, grid%ph_ch3cocho, &
                       grid%ph_hcochest, grid%ph_ch3o2h, grid%ph_ch3coo2h, &
                       grid%ph_ch3ono2, grid%ph_hcochob, grid%ph_n2o5, grid%ph_o2)
      end if
#endif

! save old concentrations for vertical mixing tendencies
      DO nv = PARAM_FIRST_SCALAR, num_chem_ct
         chem_old(its:ite, kts:kte, jts:jte, nv) = chem(its:ite, kts:kte, jts:jte, chem_ct_indices(nv))
      ENDDO

! accumulate vertical mixing tendencies
!   DO nv = PARAM_FIRST_SCALAR,num_chem_ct
!      grid%vmix_ct(its:ite,kts:kte,jts:jte,nv) = grid%vmix_ct(its:ite,kts:kte,jts:jte,nv) + &
!                                                        (chem(its:ite,kts:kte,jts:jte,chem_ct_indices(nv)) - &
!                                                     chem_old(its:ite,kts:kte,jts:jte,nv))
!  ENDDO

#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
      if ((its <= CHEM_DBG_I .and. ite >= CHEM_DBG_I) .and. &
          (jts <= CHEM_DBG_J .and. jte >= CHEM_DBG_J) .and. &
          (k_start <= CHEM_DBG_K .and. k_end >= CHEM_DBG_K)) then
         call wrf_debug(15, 'calling chem_dbg after dry_deposition_driver')
         call chem_dbg(CHEM_DBG_I, CHEM_DBG_J, CHEM_DBG_K, grid%dt, ktau, &
                       dz8w, t_phy, p_phy, rho, chem, emis_ant, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte, &
                       config_flags%kemit, &
                       grid%ph_macr, grid%ph_o31d, grid%ph_o33p, grid%ph_no2, grid%ph_no3o2, &
                       grid%ph_no3o, grid%ph_hno2, grid%ph_hno3, grid%ph_hno4, grid%ph_h2o2, &
                       grid%ph_ch2or, grid%ph_ch2om, grid%ph_ch3cho, grid%ph_ch3coch3, &
                       grid%ph_ch3coc2h5, grid%ph_hcocho, grid%ph_ch3cocho, &
                       grid%ph_hcochest, grid%ph_ch3o2h, grid%ph_ch3coo2h, &
                       grid%ph_ch3ono2, grid%ph_hcochob, grid%ph_n2o5, grid%ph_o2 &
                       )
      end if
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CHEMISTRY DRIVER/ MECHANISM (do_timestep)

   if (do_chemstep) then

      !! -- Stubbed out: do_chemstep (used to call mechanism_driver) -- !!

   end if ! Chemistry time step check

! Fill top level to prevent spurious interpolation results (no extrapolation)
   do nv = 1, num_chem
      do j = jts, jte
         do i = its, ite
            chem(i, k_end, j, nv) = chem(i, kte, j, nv)
         enddo
      enddo
   enddo

   call wrf_debug(15, 'done tileloop in chem_driver')

   if (grid%OPT_PARS_OUT == 1) then
      call wrf_debug(15, 'calculate optical output stuff')
      call aer_opt_out(TAUAER300=grid%tauaer1, TAUAER400=grid%tauaer2                            &
     &        , TAUAER600=grid%tauaer3, TAUAER999=grid%tauaer4                                    &
     &        , GAER300=grid%gaer1, GAER400=grid%gaer2, GAER600=grid%gaer3, GAER999=grid%gaer4    &
     &        , WAER300=grid%waer1, WAER400=grid%waer2, WAER600=grid%waer3, WAER999=grid%waer4 &
              , ext_coeff=grid%ext_coef, bscat_coeff=grid%bscat_coef, asym_par=grid%asym_par &
              , num_ext_coef=num_ext_coef, num_bscat_coef=num_bscat_coef, num_asym_par=num_asym_par &
     &        , dz8w=dz8w                                                                         &
     &        , IDS=ids, IDE=ide, JDS=jds, JDE=jde, KDS=kds, KDE=kde                                 &
     &        , IMS=ims, IME=ime, JMS=jms, JME=jme, KMS=kms, KME=kme                                 &
     &        , its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

   endif

!-- set upper boundary condition
   if (config_flags%have_bcs_upper) then
      call wrf_debug(15, 'set upper boundary condition')
      call tropopause_driver(grid%id, grid%dt, current_date_char, &
                             t_phy, p_phy, p8w, zmid, z_at_w, &
                             grid%tropo_lev, grid%tropo_p, grid%tropo_z, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte)
      call upper_bc_driver(grid%id, grid%dt, current_date_char, &
                           chem, p_phy, p8w, grid%tropo_lev, &
                           ids, ide, jds, jde, kds, kde, &
                           ims, ime, jms, jme, kms, kme, &
                           its, ite, jts, jte, kts, kte)
   endif

END DO chem_tile_loop_1

!-- Work around for dgnum and dgnumwet not being written to restart files.
!-- NOTE hplin 2/20/2018: Might need to stub this out for WRF to operate correctly (3.6.1)
! grid%dgnum_a1(its:ite, kts:kte, jts:jte) = grid%dgnum4d(its:ite, kts:kte, jts:jte, 1)
! grid%dgnum_a2(its:ite, kts:kte, jts:jte) = grid%dgnum4d(its:ite, kts:kte, jts:jte, 2)
! grid%dgnum_a3(its:ite, kts:kte, jts:jte) = grid%dgnum4d(its:ite, kts:kte, jts:jte, 3)

! grid%dgnumwet_a1(its:ite, kts:kte, jts:jte) = grid%dgnumwet4d(its:ite, kts:kte, jts:jte, 1)
! grid%dgnumwet_a2(its:ite, kts:kte, jts:jte) = grid%dgnumwet4d(its:ite, kts:kte, jts:jte, 2)
! grid%dgnumwet_a3(its:ite, kts:kte, jts:jte) = grid%dgnumwet4d(its:ite, kts:kte, jts:jte, 3)

END subroutine chem_driver
