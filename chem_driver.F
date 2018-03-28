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
   use module_domain, only:domain
   use module_configure
   use module_driver_constants
   use module_machine
   use module_tiles
   use module_dm
   use module_model_constants
   use module_state_description

   !! chemistry stub interfaces
   use module_chem_utilities
   use module_aer_opt_out, only:aer_opt_out
   use module_input_chem_data, only:last_chem_time, &
#if (defined(CHEM_DBG_I) && defined(CHEM_DBG_J) && defined(CHEM_DBG_K))
                                    chem_dbg, &
#endif
      get_last_gas

   use module_upper_bc_driver, only: upper_bc_driver
   use module_tropopause, only: tropopause_driver

   implicit none

   !  Input data.
   TYPE(domain), TARGET :: grid

   !  Definitions of dummy arguments to solve
# include <dummy_new_decl.inc>
# define NO_I1_OLD

   TYPE(grid_config_rec_type), INTENT(IN) :: config_flags

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

   !-- ips           start index for i in patch
   !-- ipe           end index for i in patch
   !-- jps           start index for j in patch
   !-- jpe           end index for j in patch

   integer :: ids, ide, jds, jde, kds, kde, &
              ims, ime, jms, jme, kms, kme, &
              ips, ipe, jps, jpe, kps, kpe, &
              its, ite, jts, jte, kts, kte


   ! Local Scalars
   integer :: stepave, i, j, k, l, numgas, nv, n, nr, ktau, k_start, k_end, idf, jdf, kdf
   integer :: ijulian

   ! Memorder XZY: %s-em31 is I, %s-em32 is K, %s-em33 is J
   ! Local Variables computed for chemistry usage.
   !
   ! To pass physics (meteo) variable to CHEMISTRY, edit this to create new local variables
   ! (sized as ims:ime, kms:kme, jms:jme)
   !
   ! Modify the call to module_chem_utilities::chem_prep,
   ! and modify the calls to external_driver as necessary.
   real, dimension(grid%sm31:grid%em31, grid%sm32:grid%em32, grid%sm33:grid%em33) :: &
      p_phy, u_phy, v_phy, t_phy, dz8w, t8w, p8w, rho, rri, z_at_w, vvel, zmid, rh

   ! PBL Height, dimension is (I x J) obviously...
   real, dimension(grid%sm31:grid%em31, grid%sm33:grid%em33) :: pbl_h

   ! Variables for calculating adaptive timesteps.
   TYPE(WRFU_TimeInterval) :: tmpTimeInterval
   real(KIND=8) :: curr_secs
   real(KIND=8) :: real_time_r8 ! ext. function in adapt_timestep_em.F
   logical      :: do_chemstep

   integer :: ij
   integer :: im, num_3d_m, ic, num_3d_c, num_3d_s
   integer :: ijds, ijde
   integer :: ksubt

   ! From module_ctrans_grell - used for setting a minimum range on chemistry species.
   real :: epsilc
   parameter(epsilc=1.e-30)

   real :: chem_minval, dtstepc

   character(len=256) :: current_date_char
   integer :: current_month

   intrinsic max, min

   !
   ! Setup the adaptive timestep for the chem routines. Most of this follows
   ! what is in solve_em, except for the call to adjust time_step.
   !
   ! The necesssary variables exist for the EM core and using the adaptive
   ! techniques will work even with a constant time step. In fact, they
   ! prevent issues with restarts and changed time steps. So, we will
   ! always use them with the EM core.
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

   ! Check whether to do chemistry steps (do_chemstep logical)
   do_chemstep = .false.
   if (ktau == 1) then
      do_chemstep = .true.
      grid%ktauc = 1
   else
      if ((grid%chemdt <= 0) .or. &
          (curr_secs + real(grid%dt, 8) + 0.01 >= &
           (int(curr_secs/real(grid%chemdt*60., 8) + 1, 8)*real(grid%chemdt*60., 8))) &
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
      endif
   endif

   call get_ijk_from_grid(grid, &
                          ids, ide, jds, jde, kds, kde, &
                          ims, ime, jms, jme, kms, kme, &
                          ips, ipe, jps, jpe, kps, kpe)

   call domain_clock_get(grid, current_timestr=current_date_char)
   read (current_date_char(6:7), FMT='(I2)') current_month

   grid%raincv_b(:, :) = grid%raincv(:, :)

   ! Number of moist variables
   num_3d_m = num_moist
   num_3d_c = num_chem
   num_3d_s = num_scalar
   numgas = get_last_gas(config_flags%chem_opt)

   ! Compute these starting and stopping locations for each tile and number of tiles.
   call set_tiles(grid, ids, ide, jds, jde, ips, ipe, jps, jpe)
   k_start = kps
   k_end = kpe

   ijds = min(ids, jds)
   ijde = max(ide, jde)

   chem_minval = epsilc ! chem_minval can be case dependant and set below...
   chem_select: select case(config_flags%chem_opt)
      case (1)
         call wrf_debug(15, 'Pumpkin chem_driver: chem_opt = 1. Running chemistry option.')
      case default
         call wrf_error_fatal("Pumpkin chem_driver: Unrecognized chem_opt")
   end select chem_select

   do nv = 1, num_chem
      do j = jps, jpe
         do k = kps, kpe
            do i = ips, ipe
               chem(i, k, j, nv) = max(chem(i, k, j, nv), chem_minval) ! prevent underflow in chemistry
            enddo
         enddo
      enddo
   enddo

   ! Some leftovers from original WRF chem_select. Don't touch unless you know what you are doing (I don't)
   ! $ hplin 20180327 $
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
   chem_tile_loop_1: DO ij = 1, grid%num_tiles
      its = grid%i_start(ij)
      ite = min(grid%i_end(ij), ide - 1)
      jts = grid%j_start(ij)
      jte = min(grid%j_end(ij), jde - 1)

      kts = k_start
      kte = min(k_end, kde - 1)

      CALL wrf_debug(15, 'Pumpkin chem_driver: Calling module_chem_utilities::chem_prep')

      ! chem_prep @ module_chem_utilities.f
      ! The chem_prep routine reads in the grid & configuration flags (grid%...)
      ! And does conversions onto memorder XZY local variables, so they can be 
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (do_chemstep) then
      !! -- Stubbed out: do_chemstep (used to call mechanism_driver) -- !!
      ! An example call would be below (check what you have available):
      ! call mechanism_driver(grid%id, curr_secs, ktau, grid%dt, grid%ktauc, dtstepc, config_flags, &
      !                         grid%gmt, ijulian, t_phy, moist, p8w, t8w, grid%gd_cldfr, &
      !                         p_phy, chem, rho, dz8w, grid%dx, g, &
      !                         zmid, z_at_w, grid%xlat, grid%xlong, &
      !                         grid%ttday, grid%tcosz, &
      !                         grid%ph_macr, grid%ph_o31d, grid%ph_o33p, grid%ph_no2, grid%ph_no3o2, &
      !                         grid%ph_no3o, grid%ph_hno2, grid%ph_hno3, grid%ph_hno4, grid%ph_h2o2, &
      !                         grid%ph_ch2or, grid%ph_ch2om, grid%ph_ch3cho, grid%ph_ch3coch3, &
      !                         grid%ph_ch3coc2h5, grid%ph_hcocho, grid%ph_ch3cocho, grid%ph_hcochest, &
      !                         grid%ph_ch3o2h, grid%ph_ch3coo2h, grid%ph_ch3ono2, grid%ph_hcochob, &
      !                         grid%ph_n2o5, grid%ph_o2, grid%backg_oh, grid%backg_h2o2, grid%backg_no3, &
      !                         grid%addt, grid%addx, grid%addc, grid%etep, &
      !                         grid%oltp, grid%olip, grid%cslp, grid%limp, grid%hc5p, grid%hc8p, grid%tolp, &
      !                         grid%xylp, grid%apip, grid%isop, grid%hc3p, grid%ethp, grid%o3p, grid%tco3, &
      !                         grid%mo2, grid%o1d, grid%olnn, grid%rpho, grid%xo2, &
      !                         grid%ketp, grid%olnd, &
      !                         ids, ide, jds, jde, kds, kde, &
      !                         ims, ime, jms, jme, kms, kme, &
      !                         its, ite, jts, jte, kts, kte)

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
end do chem_tile_loop_1

!-- Work around for dgnum and dgnumwet not being written to restart files.
!-- NOTE hplin 2/20/2018: Might need to stub this out for WRF to operate correctly (3.6.1)
! grid%dgnum_a1(its:ite, kts:kte, jts:jte) = grid%dgnum4d(its:ite, kts:kte, jts:jte, 1)
! grid%dgnum_a2(its:ite, kts:kte, jts:jte) = grid%dgnum4d(its:ite, kts:kte, jts:jte, 2)
! grid%dgnum_a3(its:ite, kts:kte, jts:jte) = grid%dgnum4d(its:ite, kts:kte, jts:jte, 3)

! grid%dgnumwet_a1(its:ite, kts:kte, jts:jte) = grid%dgnumwet4d(its:ite, kts:kte, jts:jte, 1)
! grid%dgnumwet_a2(its:ite, kts:kte, jts:jte) = grid%dgnumwet4d(its:ite, kts:kte, jts:jte, 2)
! grid%dgnumwet_a3(its:ite, kts:kte, jts:jte) = grid%dgnumwet4d(its:ite, kts:kte, jts:jte, 3)

END subroutine chem_driver
