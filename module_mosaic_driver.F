module module_mosaic_driver

#define cw_species_are_in_registry

contains

   subroutine mosaic_aerchem_driver( &
      id, curr_secs, ktau, dtstep, ktauc, dtstepc, config_flags, &
      t_phy, rho_phy, p_phy, &
      moist, chem, vbs_nbin, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      use module_configure, only:grid_config_rec_type, &
         p_qv, &
         p_so2, p_ho2, p_so4aj, p_corn, p_hcl, p_mtf, &
         p_so4_a01, p_water_a01, p_num_a01, &
         p_so4_a04, p_water_a04, p_num_a04

      use module_state_description, only:num_moist, num_chem

      implicit none

      integer, intent(in) :: &
         id, ktau, ktauc, &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte, &
         vbs_nbin(1)
!   id - domain index
!   ktau - time step number
!   ktauc - gas and aerosol chemistry time step number

!   [ids:ide, kds:kde, jds:jde] - spatial (x,z,y) indices for "domain"
!   [ims:ime, kms:kme, jms:jme] - spatial (x,z,y) indices for "memory"
!   Most arrays that are arguments to chem_driver
!   are dimensioned with these spatial indices.
!   [its:ite, kts:kte, jts:jte] - spatial (x,z,y) indices for "tile"
!   chem_driver and routines under it do calculations
!   over these spatial indices.

      real(kind=8), intent(in) :: curr_secs
      real, intent(in) :: dtstep, dtstepc
!   dtstep - main model time step (s)
!   dtstepc - time step for gas and aerosol chemistry(s)

      real, intent(in), &
         dimension(ims:ime, kms:kme, jms:jme) :: &
         t_phy, rho_phy, p_phy
!   t_phy - temperature (K)
!   rho_phy - air density (kg/m^3)
!   p_phy - air pressure (Pa)

      real, intent(in), &
         dimension(ims:ime, kms:kme, jms:jme, 1:num_moist) :: &
         moist
!   moist - mixing ratios of moisture species (water vapor,
!   cloud water, ...) (kg/kg for mass species, #/kg for number species)

      real, intent(inout), &
         dimension(ims:ime, kms:kme, jms:jme, 1:num_chem) :: &
         chem
!   chem - mixing ratios of trace gase (ppm) and aerosol species
!   (ug/kg for mass species, #/kg for number species)

      type(grid_config_rec_type), intent(in) :: config_flags
!   config_flags - configuration and control parameters

      return
   end subroutine mosaic_aerchem_driver

!-----------------------------------------------------------------------
   subroutine sum_pm_mosaic( &
      alt, chem, &
      pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      USE module_state_description, only:num_chem
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10

      REAL :: mass

      INTEGER :: i, imax, j, jmax, k, kmax, n, itype, iphase

   end subroutine sum_pm_mosaic

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine sum_pm_mosaic_vbs2( &
      alt, chem, &
      pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      USE module_state_description, only:num_chem
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10

      REAL :: mass

      INTEGER :: i, imax, j, jmax, k, kmax, n, itype, iphase

   end subroutine sum_pm_mosaic_vbs2

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine sum_pm_mosaic_vbs0( &
      alt, chem, &
      pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      USE module_state_description, only:num_chem
      USE module_data_mosaic_asect
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10

      REAL :: mass

      INTEGER :: i, imax, j, jmax, k, kmax, n, itype, iphase
   end subroutine sum_pm_mosaic_vbs0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine sum_pm_mosaic_vbs4( &
      alt, chem, &
      pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      USE module_state_description, only:num_chem
      USE module_data_mosaic_asect
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10

      REAL :: mass

      INTEGER :: i, imax, j, jmax, k, kmax, n, itype, iphase

   end subroutine sum_pm_mosaic_vbs4

   subroutine sum_vbs0(alt, chem, &
                       hoa_a01, hoa_a02, hoa_a03, hoa_a04, &
                       bboa_a01, bboa_a02, bboa_a03, bboa_a04, &
                       soa_a01, soa_a02, soa_a03, soa_a04, &
                       bbsoa_a01, bbsoa_a02, bbsoa_a03, bbsoa_a04, &
                       biog_a01, biog_a02, biog_a03, biog_a04, &
                       asmpsoa_a01, asmpsoa_a02, asmpsoa_a03, asmpsoa_a04, &
                       arosoa_a01, arosoa_a02, arosoa_a03, arosoa_a04, &
                       totoa_a01, totoa_a02, totoa_a03, totoa_a04, &
                       biog_v1, biog_v2, biog_v3, biog_v4, &
                       ant_v1, ant_v2, ant_v3, ant_v4, &
                       smpa_v1, smpbb_v1, &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte)

      USE module_state_description, only:num_chem
      USE module_data_mosaic_asect
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: hoa_a01, hoa_a02, hoa_a03, hoa_a04, &
                        bboa_a01, bboa_a02, bboa_a03, bboa_a04, &
                        soa_a01, soa_a02, soa_a03, soa_a04, &
                        bbsoa_a01, bbsoa_a02, bbsoa_a03, bbsoa_a04, &
                        biog_a01, biog_a02, biog_a03, biog_a04, &
                        asmpsoa_a01, asmpsoa_a02, asmpsoa_a03, asmpsoa_a04, &
                        arosoa_a01, arosoa_a02, arosoa_a03, arosoa_a04, &
                        totoa_a01, totoa_a02, totoa_a03, totoa_a04, &
                        biog_v1, biog_v2, biog_v3, biog_v4, &
                        ant_v1, ant_v2, ant_v3, ant_v4, smpa_v1, smpbb_v1

      INTEGER :: i, imax, j, jmax, k, kmax, n, itype, iphase

   end subroutine sum_vbs0

!-----------------------------------------------------------------------

   subroutine sum_vbs2( &
      alt, chem, &
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
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      USE module_state_description, only:num_chem
      USE module_data_mosaic_asect
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: hoa_a01, hoa_a02, hoa_a03, hoa_a04, &
                        bboa_a01, bboa_a02, bboa_a03, bboa_a04, &
                        soa_a01, soa_a02, soa_a03, soa_a04, &
                        bbsoa_a01, bbsoa_a02, bbsoa_a03, bbsoa_a04, &
                        biog_a01, biog_a02, biog_a03, biog_a04, &
                        hsoa_a01, hsoa_a02, hsoa_a03, hsoa_a04, &
                        arosoa_a01, arosoa_a02, arosoa_a03, arosoa_a04, &
                        totoa_a01, totoa_a02, totoa_a03, totoa_a04, &
                        hsoa_c, hsoa_o, bbsoa_c, bbsoa_o, &
                        biog_v1, biog_v2, biog_v3, biog_v4, &
                        ant_v1, ant_v2, ant_v3, ant_v4

      INTEGER :: i, imax, j, jmax, k, kmax, n, itype, iphase

   end subroutine sum_vbs2

!----------------------------------------------------------------------

!-----------------------------------------------------------------------

   subroutine sum_vbs4( &
      alt, chem, &
      hoa_a01, hoa_a02, hoa_a03, hoa_a04, &
      soa_a01, soa_a02, soa_a03, soa_a04, &
      biog_a01, biog_a02, biog_a03, biog_a04, &
      totoa_a01, totoa_a02, totoa_a03, totoa_a04, &
      biog_v1, biog_v2, biog_v3, biog_v4, &
      ant_v1, ant_v2, ant_v3, ant_v4, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      USE module_state_description, only:num_chem
      USE module_data_mosaic_asect
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: hoa_a01, hoa_a02, hoa_a03, hoa_a04, &
                        soa_a01, soa_a02, soa_a03, soa_a04, &
                        biog_a01, biog_a02, biog_a03, biog_a04, &
                        totoa_a01, totoa_a02, totoa_a03, totoa_a04, &
                        biog_v1, biog_v2, biog_v3, biog_v4, &
                        ant_v1, ant_v2, ant_v3, ant_v4

      INTEGER :: i, imax, j, jmax, k, kmax, n, itype, iphase

   end subroutine sum_vbs4

!----------------------------------------------------------------------
   subroutine sum_vbs9(alt, chem, &
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
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      USE module_state_description, only: num_chem
      USE module_data_mosaic_asect
      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: hoa_a01, hoa_a02, hoa_a03, hoa_a04, &
                        bboa_a01, bboa_a02, bboa_a03, bboa_a04, &
                        soa_a01, soa_a02, soa_a03, soa_a04, &
                        bbsoa_a01, bbsoa_a02, bbsoa_a03, bbsoa_a04, &
                        biog_a01, biog_a02, biog_a03, biog_a04, &
                        hsoa_a01, hsoa_a02, hsoa_a03, hsoa_a04, &
                        arosoa_a01, arosoa_a02, arosoa_a03, arosoa_a04, &
                        totoa_a01, totoa_a02, totoa_a03, totoa_a04, &
                        hsoa_c, hsoa_o, bbsoa_c, bbsoa_o, &
                        biog_v1, biog_v2, biog_v3, biog_v4, &
                        ant_v1, ant_v2, ant_v3, ant_v4

      INTEGER :: i, imax, j, jmax, k, kmax, n, itype, iphase

   end subroutine sum_vbs9

   subroutine mapaer_tofrom_host(imap, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte, &
                                 it, jt, ktmaps, ktmape, &
                                 num_moist, num_chem, moist, chem, &
                                 t_phy, p_phy, rho_phy)
      implicit none

!   imap determines mapping direction (chem-->rsub if <=0, rsub-->chem if >0)
      integer, intent(in) :: imap
!   wrf array dimensions
      integer, intent(in) :: num_moist, num_chem
      integer, intent(in) :: ims, ime, jms, jme, kms, kme
      integer, intent(in) :: its, ite, jts, jte, kts, kte
!   do mapping for wrf i,k,j grid points = [it,ktmaps:ktmape,jt]
      integer, intent(in) :: it, jt, ktmaps, ktmape
!
      real, intent(in), dimension(ims:ime, kms:kme, jms:jme) :: &
         t_phy, rho_phy, p_phy

      real, intent(in), &
         dimension(ims:ime, kms:kme, jms:jme, 1:num_moist) :: &
         moist

      real, intent(inout), &
         dimension(ims:ime, kms:kme, jms:jme, 1:num_chem) :: &
         chem

      return

   end subroutine mapaer_tofrom_host

   subroutine init_data_mosaic_asect(n2o5_hetchem, is_aerosol)
      implicit none
      use module_state_description, only: num_chem
      logical, intent(out) :: is_aerosol(num_chem)
      integer, intent(in) :: n2o5_hetchem

   end subroutine init_data_mosaic_asect

   subroutine init_data_mosaic_ptr(is_aerosol)
      implicit none
      use module_state_description, only: num_chem
      logical, intent(out) :: is_aerosol(num_chem)

      !! Stubbed out completely

      return
   end subroutine init_data_mosaic_ptr

   subroutine aerchem_debug_dump(iflag, iclm, jclm, dtchem)
      implicit none
      integer iflag, iclm, jclm
      real dtchem

      return
   end subroutine aerchem_debug_dump

end module module_mosaic_driver
