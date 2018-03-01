! This is stubbed for the purpose of satisfying a dependency in solve_em.F

module module_input_tracer
   use module_state_description, only:tracer_smoke,tracer_test1,tracer_test2,param_first_scalar,p_tr17_1,p_tr17_2,p_tr17_3,p_tr17_4,p_tr17_5,p_tr17_6,p_tr17_7,p_tr17_8
contains
   subroutine initialize_tracer(chem, chem_in_opt, &
                                tracer_opt, num_chem, &
                                ids, ide, jds, jde, kds, kde, & ! domain dims
                                ims, ime, jms, jme, kms, kme, & ! memory dims
                                ips, ipe, jps, jpe, kps, kpe, & ! patch  dims
                                its, ite, jts, jte, kts, kte)
      integer, intent(in)    :: chem_in_opt, tracer_opt, num_chem
      integer, intent(in)    :: ids, ide, jds, jde, kds, kde
      integer, intent(in)    :: ims, ime, jms, jme, kms, kme
      integer, intent(in)    :: ips, ipe, jps, jpe, kps, kpe
      integer, intent(in)    :: its, ite, jts, jte, kts, kte
      real, dimension(ims:ime, kms:kme, jms:jme, num_chem), intent(inout) :: chem
      if (chem_in_opt == 1) return
   end subroutine initialize_tracer

#if (EM_CORE == 1 )
   subroutine flow_dep_bdy_tracer(chem, &
                                  chem_bxs, chem_btxs, &
                                  chem_bxe, chem_btxe, &
                                  chem_bys, chem_btys, &
                                  chem_bye, chem_btye, &
                                  dt, &
                                  spec_bdy_width, z, &
                                  have_bcs_chem, &
                                  u, v, tracer_opt, alt, &
                                  t, pb, p, t0, p1000mb, rcp, ph, phb, g, &
                                  spec_zone, ic, &
                                  ids, ide, jds, jde, kds, kde, & ! domain dims
                                  ims, ime, jms, jme, kms, kme, & ! memory dims
                                  ips, ipe, jps, jpe, kps, kpe, & ! patch  dims
                                  its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent(in)    :: tracer_opt
      integer, intent(in)    :: ids, ide, jds, jde, kds, kde
      integer, intent(in)    :: ims, ime, jms, jme, kms, kme
      integer, intent(in)    :: ips, ipe, jps, jpe, kps, kpe
      integer, intent(in)    :: its, ite, jts, jte, kts, kte
      integer, intent(in)    :: spec_zone, spec_bdy_width, ic
      real, intent(in)    :: dt

      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: chem
      real, dimension(jms:jme, kds:kde, spec_bdy_width), intent(in) :: chem_bxs, chem_bxe, chem_btxs, chem_btxe
      real, dimension(ims:ime, kds:kde, spec_bdy_width), intent(in) :: chem_bys, chem_bye, chem_btys, chem_btye
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: z
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: alt
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: u
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: v
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(in) :: ph, phb, t, pb, p
      real, intent(in) :: g, rcp, t0, p1000mb

      integer    :: i, j, k, numgas
      integer    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      integer    :: i_inner, j_inner
      integer    :: b_dist
      integer    :: i_bdy_method
      real tempfac, convfac
      logical, optional    :: have_bcs_chem

   end subroutine flow_dep_bdy_tracer
#else
#if ( WRF_CHEM == 1 )
   subroutine flow_dep_bdy_tracer(chem, chem_b, chem_bt, dt, &
                                  spec_bdy_width, z, &
                                  ijds, ijde, have_bcs_chem, &
                                  u, v, tracer_opt, alt, &
                                  t, pb, p, t0, p1000mb, rcp, ph, phb, g, &
                                  spec_zone, ic, &
                                  ids, ide, jds, jde, kds, kde, & ! domain dims
                                  ims, ime, jms, jme, kms, kme, & ! memory dims
                                  ips, ipe, jps, jpe, kps, kpe, & ! patch  dims
                                  its, ite, jts, jte, kts, kte)

      implicit none

      integer, intent(in)    :: tracer_opt
      integer, intent(in)    :: ids, ide, jds, jde, kds, kde
      integer, intent(in)    :: ims, ime, jms, jme, kms, kme
      integer, intent(in)    :: ips, ipe, jps, jpe, kps, kpe
      integer, intent(in)    :: its, ite, jts, jte, kts, kte
      integer, intent(in)    :: ijds, ijde
      integer, intent(in)    :: spec_zone, spec_bdy_width, ic
      real, intent(in)    :: dt

      real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: chem
      real, dimension(ijds:ijde, kds:kde, spec_bdy_width, 4), intent(in) :: chem_b
      real, dimension(ijds:ijde, kds:kde, spec_bdy_width, 4), intent(in) :: chem_bt
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: z
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: alt
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: u
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: v
      real, dimension(ims:ime, kms:kme, jms:jme), &
         intent(in) :: &
         ph, phb, t, pb, p
      real, intent(in) :: g, rcp, t0, p1000mb

      integer    :: i, j, k, numgas
      integer    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      integer    :: i_inner, j_inner
      integer    :: b_dist
      integer    :: i_bdy_method
      real tempfac, convfac
      real       :: tracer_bv_def
      logical    :: have_bcs_chem

   end subroutine flow_dep_bdy_tracer
#endif
#endif
   subroutine set_tracer(dtstep, ktau, pbl_h, tracer, t, tracer_opt, num_tracer, &
                         z, ht, ids, ide, jds, jde, kds, kde, &
                         ims, ime, jms, jme, kms, kme, &
                         its, ite, jts, jte, kts, kte)
      integer, intent(in)    :: ktau, tracer_opt, num_tracer
      integer, intent(in)    :: ids, ide, jds, jde, kds, kde
      integer, intent(in)    :: ims, ime, jms, jme, kms, kme
      integer, intent(in)    :: its, ite, jts, jte, kts, kte
      real, dimension(ims:ime, kms:kme, jms:jme, num_tracer), intent(inout) :: tracer
      real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: t, z
      real, dimension(ims:ime, jms:jme), intent(in) :: PBL_H, HT
      real, intent(in) :: dtstep
      integer :: count_trop, count_pbl

   end subroutine set_tracer

   subroutine bdy_tracer_value(trac, trac_b, trac_bt, dt, ic)

      implicit none

      real, intent(out)  :: trac
      real, intent(in)   :: trac_b
      real, intent(in)   :: trac_bt
      real, intent(in)   :: dt
      integer, intent(in)   :: ic

      real :: epsilc = 1.e-12

      trac = max(epsilc, trac_b + trac_bt*dt)

      return
   end subroutine bdy_tracer_value

END MODULE module_input_tracer
