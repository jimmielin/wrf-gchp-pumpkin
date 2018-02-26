! This is stubbed for the purpose of satisfying a dependency in solve_em.F

MODULE module_input_tracer
   USE module_state_description, only:tracer_smoke,tracer_test1,tracer_test2,param_first_scalar,p_tr17_1,p_tr17_2,p_tr17_3,p_tr17_4,p_tr17_5,p_tr17_6,p_tr17_7,p_tr17_8
CONTAINS
   SUBROUTINE initialize_tracer(chem, chem_in_opt, &
                                tracer_opt, num_chem, &
                                ids, ide, jds, jde, kds, kde, & ! domain dims
                                ims, ime, jms, jme, kms, kme, & ! memory dims
                                ips, ipe, jps, jpe, kps, kpe, & ! patch  dims
                                its, ite, jts, jte, kts, kte)
      INTEGER, INTENT(IN)    :: chem_in_opt, tracer_opt, num_chem
      INTEGER, INTENT(IN)    :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN)    :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN)    :: ips, ipe, jps, jpe, kps, kpe
      INTEGER, INTENT(IN)    :: its, ite, jts, jte, kts, kte
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), INTENT(INOUT) :: chem
      if (chem_in_opt == 1) return
   END SUBROUTINE initialize_tracer

#if (EM_CORE == 1 )
   SUBROUTINE flow_dep_bdy_tracer(chem, &
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

      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: tracer_opt
      INTEGER, INTENT(IN)    :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN)    :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN)    :: ips, ipe, jps, jpe, kps, kpe
      INTEGER, INTENT(IN)    :: its, ite, jts, jte, kts, kte
      INTEGER, INTENT(IN)    :: spec_zone, spec_bdy_width, ic
      REAL, INTENT(IN)    :: dt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: chem
      REAL, DIMENSION(jms:jme, kds:kde, spec_bdy_width), INTENT(IN) :: chem_bxs, chem_bxe, chem_btxs, chem_btxe
      REAL, DIMENSION(ims:ime, kds:kde, spec_bdy_width), INTENT(IN) :: chem_bys, chem_bye, chem_btys, chem_btye
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: z
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: alt
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: u
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: v
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: &
         ph, phb, t, pb, p
      real, INTENT(IN) :: g, rcp, t0, p1000mb

      INTEGER    :: i, j, k, numgas
      INTEGER    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      INTEGER    :: i_inner, j_inner
      INTEGER    :: b_dist
      integer    :: i_bdy_method
      real tempfac, convfac
      logical, optional    :: have_bcs_chem

   END SUBROUTINE flow_dep_bdy_tracer
#else
#if ( WRF_CHEM == 1 )
   SUBROUTINE flow_dep_bdy_tracer(chem, chem_b, chem_bt, dt, &
                                  spec_bdy_width, z, &
                                  ijds, ijde, have_bcs_chem, &
                                  u, v, tracer_opt, alt, &
                                  t, pb, p, t0, p1000mb, rcp, ph, phb, g, &
                                  spec_zone, ic, &
                                  ids, ide, jds, jde, kds, kde, & ! domain dims
                                  ims, ime, jms, jme, kms, kme, & ! memory dims
                                  ips, ipe, jps, jpe, kps, kpe, & ! patch  dims
                                  its, ite, jts, jte, kts, kte)

      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: tracer_opt
      INTEGER, INTENT(IN)    :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN)    :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN)    :: ips, ipe, jps, jpe, kps, kpe
      INTEGER, INTENT(IN)    :: its, ite, jts, jte, kts, kte
      INTEGER, INTENT(IN)    :: ijds, ijde
      INTEGER, INTENT(IN)    :: spec_zone, spec_bdy_width, ic
      REAL, INTENT(IN)    :: dt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: chem
      REAL, DIMENSION(ijds:ijde, kds:kde, spec_bdy_width, 4), INTENT(IN) :: chem_b
      REAL, DIMENSION(ijds:ijde, kds:kde, spec_bdy_width, 4), INTENT(IN) :: chem_bt
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: z
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: alt
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: u
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: v
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: &
         ph, phb, t, pb, p
      real, INTENT(IN) :: g, rcp, t0, p1000mb

      INTEGER    :: i, j, k, numgas
      INTEGER    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      INTEGER    :: i_inner, j_inner
      INTEGER    :: b_dist
      integer    :: i_bdy_method
      real tempfac, convfac
      real       :: tracer_bv_def
      logical    :: have_bcs_chem

   END SUBROUTINE flow_dep_bdy_tracer
#endif
#endif
   SUBROUTINE set_tracer(dtstep, ktau, pbl_h, tracer, t, tracer_opt, num_tracer, &
                         z, ht, ids, ide, jds, jde, kds, kde, &
                         ims, ime, jms, jme, kms, kme, &
                         its, ite, jts, jte, kts, kte)
      INTEGER, INTENT(IN)    :: ktau, tracer_opt, num_tracer
      INTEGER, INTENT(IN)    :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN)    :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN)    :: its, ite, jts, jte, kts, kte
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_tracer), INTENT(INOUT) :: tracer
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: t, z
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: PBL_H, HT
      REAL, INTENT(IN) :: dtstep
      INTEGER:: count_trop, count_pbl

   END SUBROUTINE set_tracer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE bdy_tracer_value(trac, trac_b, trac_bt, dt, ic)

      IMPLICIT NONE

      REAL, intent(OUT)  :: trac
      REAL, intent(IN)   :: trac_b
      REAL, intent(IN)   :: trac_bt
      REAL, intent(IN)   :: dt
      INTEGER, intent(IN)   :: ic

      REAL                  :: epsilc = 1.E-12

      trac = max(epsilc, trac_b + trac_bt*dt)

      RETURN
   END SUBROUTINE bdy_tracer_value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE module_input_tracer
