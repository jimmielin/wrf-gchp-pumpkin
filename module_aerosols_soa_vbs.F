MODULE module_aerosols_soa_vbs

   USE module_state_description
   USE module_data_soa_vbs

   IMPLICIT NONE
#define cw_species_are_in_registry

CONTAINS

   SUBROUTINE soa_vbs_driver(id, ktau, dtstep, t_phy, moist, aerwrf, p8w, &
                             t8w, alt, p_phy, chem, rho_phy, dz8w, rh, z, z_at_w, &
                             h2oaj, h2oai, nu3, ac3, cor3, asulf, ahno3, anh3, &
                             vcsulf_old, &
                             vdrog3, &
                             kemit, brch_ratio, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte)

!   USE module_configure, only: grid_config_rec_type
!   TYPE (grid_config_rec_type), INTENT (in) :: config_flags

      INTEGER, INTENT(IN)  :: ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              its, ite, jts, jte, kts, kte, &
                              kemit, id, ktau

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_moist), &
         INTENT(IN) ::                                      moist

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(INOUT) ::                                   chem
!
! following are aerosol arrays that are not advected
!
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT) :: &
         h2oaj, h2oai, nu3, ac3, cor3, asulf, ahno3, anh3

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT) ::    brch_ratio

!           cvasoa1,cvasoa2,    &
!           cvasoa3,cvasoa4,cvbsoa1,cvbsoa2,cvbsoa3,cvbsoa4

      REAL, DIMENSION(ims:ime, kms:kme - 0, jms:jme, ldrog_vbs), &
         INTENT(IN) :: VDROG3
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) ::                             t_phy, &
                       alt, &
                       p_phy, &
                       dz8w, &
                       rh, & ! fractional relative humidity
                       z, &
                       t8w, p8w, z_at_w, &
                       aerwrf, &
                       rho_phy
      REAL, DIMENSION(ims:ime, kms:kme - 0, jms:jme), &
         INTENT(IN) ::   vcsulf_old
      REAL, INTENT(IN) ::   dtstep

      REAL drog_in(ldrog_vbs) ! anthropogenic AND biogenic organic aerosol precursors [ug m**-3 s**-1]

!      REAL condvap_in(lspcv) ! condensable vapors [ug m**-3]
      REAL, PARAMETER :: rgas = 8.314510
      REAL convfac, convfac2

!...BLKSIZE set to one in column model ciarev02
      INTEGER, PARAMETER :: blksize = 1

!...number of aerosol species
!  number of species (gas + aerosol)
      INTEGER nspcsda
      PARAMETER(nspcsda=l1ae) !bs
! (internal aerosol dynamics)
!bs # of anth. cond. vapors in SOA_VBS
      INTEGER nacv
      PARAMETER(nacv=lcva) !bs # of anth. cond. vapors in CTM
!bs total # of cond. vapors in SOA_VBS
      INTEGER ncv
      PARAMETER(ncv=lspcv) !bs
!bs total # of cond. vapors in CTM
      REAL cblk(blksize, nspcsda) ! main array of variables
      ! particles [ug/m^3/s]
      REAL soilrat_in
      ! emission rate of soil derived coars
      ! input HNO3 to CBLK [ug/m^3]
      REAL nitrate_in
      ! input NH3 to CBLK  [ug/m^3]
      REAL nh3_in
      ! input SO4 vapor    [ug/m^3]
      REAL hcl_in

      REAL vsulf_in

      REAL so4rat_in
      ! input SO4 formation[ug/m^3/sec]
      REAL epm25i(blksize), epm25j(blksize), epmcoarse(blksize)
      ! Emission rate of i-mode EC [ug m**-3 s**-1]
      REAL eeci_in
      ! Emission rate of j-mode EC [ug m**-3 s**-1]
      REAL eecj_in
      ! Emission rate of j-mode org. aerosol [ug m**-
      REAL eorgi_in

      REAL eorgj_in ! Emission rate of j-mode org. aerosol [ug m**-
      REAL pres ! pressure in cb
      REAL temp ! temperature in K
      !     REAL relhum   ! rel. humidity (0,1)
      REAL brrto

      REAL :: p(kts:kte), t(kts:kte), rh0(kts:kte)

!...molecular weights                   ciarev02
! these molecular weights aren't used at all

! molecular weight for SO4
      REAL mwso4
      PARAMETER(mwso4=96.0576)

! molecular weight for HNO3
      REAL mwhno3
      PARAMETER(mwhno3=63.01287)

! molecular weight for NH3
      REAL mwnh3
      PARAMETER(mwnh3=17.03061)

! molecular weight for HCL
      REAL mwhcl
      PARAMETER(mwhcl=36.46100)

!bs molecular weight for Elemental Carbon
      REAL mwec
      PARAMETER(mwec=12.0)

      INTEGER :: i, j, k, l, debug_level

   END SUBROUTINE soa_vbs_driver

   SUBROUTINE sum_pm_soa_vbs( &
      alt, chem, h2oaj, h2oai, &
      pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10, dust_opt, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      INTEGER, INTENT(IN) ::     dust_opt, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(IN) :: chem

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: alt, h2oaj, h2oai

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(OUT) :: pm2_5_dry, pm2_5_water, pm2_5_dry_ec, pm10

      INTEGER :: i, ii, j, jj, k, n

   END SUBROUTINE sum_pm_soa_vbs
! ///////////////////////////////////////////////////

   SUBROUTINE soa_vbs_depdriver(id, config_flags, ktau, dtstep, &
                                ust, t_phy, moist, p8w, t8w, rmol, znt, pbl, &
                                alt, p_phy, chem, rho_phy, dz8w, rh, z, z_at_w, &
                                h2oaj, h2oai, nu3, ac3, cor3, asulf, ahno3, anh3, &
                                aer_res, vgsa, &
                                numaer, &
                                ids, ide, jds, jde, kds, kde, &
                                ims, ime, jms, jme, kms, kme, &
                                its, ite, jts, jte, kts, kte)

      USE module_configure, only:grid_config_rec_type
      TYPE(grid_config_rec_type), INTENT(IN) :: config_flags

      INTEGER, INTENT(IN)    ::       numaer, &
                                ids, ide, jds, jde, kds, kde, &
                                ims, ime, jms, jme, kms, kme, &
                                its, ite, jts, jte, kts, kte, &
                                id, ktau

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_moist), &
         INTENT(IN) ::                                   moist
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(INOUT) ::                                   chem
!
! following are aerosol arrays that are not advected
!
      REAL, DIMENSION(its:ite, jts:jte, numaer), &
         INTENT(INOUT) :: &
         vgsa
      REAL, DIMENSION(its:ite, jts:jte), &
         INTENT(INOUT) :: &
         aer_res

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT) :: &
         h2oaj, h2oai, nu3, ac3, cor3, asulf, ahno3, anh3

! no vapors
!cvaro1,cvaro2,cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) ::                            t_phy, &
                       alt, &
                       p_phy, &
                       dz8w, &
                       rh, &
                       z, &
                       t8w, p8w, z_at_w, &
                       rho_phy
      REAL, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN) ::                     ust, rmol, pbl, znt
      REAL, INTENT(IN) ::                                 dtstep

      REAL, PARAMETER   ::   rgas = 8.314510
      REAL convfac, convfac2
!...BLKSIZE set to one in column model ciarev02

      INTEGER, PARAMETER   :: blksize = 1

!...number of aerosol species
!  number of species (gas + aerosol)
      INTEGER nspcsda
      PARAMETER(nspcsda=l1ae) !bs
! (internal aerosol dynamics)
!bs # of anth. cond. vapors in SOA_VBS
      INTEGER nacv
      PARAMETER(nacv=lcva) !bs # of anth. cond. vapors in CTM
!bs total # of cond. vapors in SOA_VBS
      INTEGER, PARAMETER :: ncv = lspcv ! number of bins=8
!bs total # of cond. vapors in CTM
      REAL cblk(blksize, nspcsda) ! main array of variables
      ! particles [ug/m^3/s]
      REAL soilrat_in
      ! emission rate of soil derived coars
      ! input HNO3 to CBLK [ug/m^3]
      REAL nitrate_in
      ! input NH3 to CBLK  [ug/m^3]
      REAL nh3_in
      ! input SO4 vapor    [ug/m^3]
      REAL vsulf_in

      REAL so4rat_in
      ! input SO4 formation[ug/m^3/sec]
      ! pressure in cb
      REAL pres
      ! temperature in K
      REAL temp
      !bs
      REAL relhum
      ! rel. humidity (0,1)
      REAL ::  p(kts:kte), t(kts:kte), rh0(kts:kte)

!...molecular weights                   ciarev02

! molecular weight for SO4
      REAL mwso4
      PARAMETER(mwso4=96.0576)

! molecular weight for HNO3
      REAL mwhno3
      PARAMETER(mwhno3=63.01287)

! molecular weight for NH3
      REAL mwnh3
      PARAMETER(mwnh3=17.03061)

!bs molecular weight for Organic Spec
!     REAL mworg
!     PARAMETER (mworg=175.0)

!bs molecular weight for Elemental Ca
      REAL mwec
      PARAMETER(mwec=12.0)

      INTEGER NUMCELLS ! actual number of cells in arrays ( default is 1 in box model)
!ia                       kept to 1 in current version of column model
      PARAMETER(NUMCELLS=1)

      REAL RA(BLKSIZE) ! aerodynamic resistance [ s m**-1 ]
      REAL USTAR(BLKSIZE) ! surface friction velocity [ m s**-1 ]
      REAL WSTAR(BLKSIZE) ! convective velocity scale [ m s**-1 ]
      REAL PBLH(BLKSIZE) ! PBL height (m)
      REAL ZNTT(BLKSIZE) ! Surface roughness length (m)
      REAL RMOLM(BLKSIZE) ! Inverse of Monin-Obukhov length (1/m)

      REAL BLKPRS(BLKSIZE) ! pressure in cb
      REAL BLKTA(BLKSIZE) ! temperature in K
      REAL BLKDENS(BLKSIZE) ! Air density in kg/m3
!
! *** OUTPUT:
!
! *** atmospheric properties

      REAL XLM(BLKSIZE) ! atmospheric mean free path [ m ]
      REAL AMU(BLKSIZE) ! atmospheric dynamic viscosity [ kg/m s ]

! *** followng is for future version
      REAL VSED(BLKSIZE, NASPCSSED) ! sedimentation velocity [ m s**-1 ]
      REAL VDEP(BLKSIZE, NASPCSDEP) ! deposition velocity [ m s**-1 ]

! *** modal diameters: [ m ]
      REAL DGNUC(BLKSIZE) ! nuclei mode geometric mean diameter  [ m ]
      REAL DGACC(BLKSIZE) ! accumulation geometric mean diameter [ m ]
      REAL DGCOR(BLKSIZE) ! coarse mode geometric mean diameter  [ m ]

! *** aerosol properties:
! *** Modal mass concentrations [ ug m**3 ]
      REAL PMASSN(BLKSIZE) ! mass concentration in Aitken mode
      REAL PMASSA(BLKSIZE) ! mass concentration in accumulation mode
      REAL PMASSC(BLKSIZE) ! mass concentration in coarse mode

! *** average modal particle densities  [ kg/m**3 ]
      REAL PDENSN(BLKSIZE) ! average particle density in nuclei mode
      REAL PDENSA(BLKSIZE) ! average particle density in accumulation mode
      REAL PDENSC(BLKSIZE) ! average particle density in coarse mode

! *** average modal Knudsen numbers
      REAL KNNUC(BLKSIZE) ! nuclei mode  Knudsen number
      REAL KNACC(BLKSIZE) ! accumulation Knudsen number
      REAL KNCOR(BLKSIZE) ! coarse mode  Knudsen number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      INTEGER :: i, j, k, l

   END SUBROUTINE soa_vbs_depdriver
! ///////////////////////////////////////////////////

   SUBROUTINE actcof(cat, an, gama, molnu, phimult)

! Normal, successful completion
      INTEGER xstat0
      PARAMETER(xstat0=0)
! File I/O error
      INTEGER xstat1
      PARAMETER(xstat1=1)
! Execution error
      INTEGER xstat2
      PARAMETER(xstat2=2)
! Special  error
      INTEGER xstat3
      PARAMETER(xstat3=3)
      CHARACTER*120 xmsg

!...........PARAMETERS and their descriptions:
! number of cations
      INTEGER ncat
      PARAMETER(ncat=2)

! number of anions
      INTEGER nan
      PARAMETER(nan=3)

!...........ARGUMENTS and their descriptions
! tot # moles of all ions
      REAL molnu
! multicomponent paractical osmo
      REAL phimult
      REAL cat(ncat) ! cation conc in moles/kg (input
      REAL an(nan) ! anion conc in moles/kg (input)
      REAL gama(ncat, nan)

   END SUBROUTINE actcof

   SUBROUTINE aeroproc(blksize, nspcsda, numcells, layer, cblk, dt, blkta, blkprs, &
                       blkdens, blkrh, so4rat, organt1rat, organt2rat, organt3rat, organt4rat, &
                       orgbio1rat, orgbio2rat, orgbio3rat, orgbio4rat, drog, ldrog_vbs, ncv, nacv, epm25i, &
                       epm25j, eorgi, eorgj, eeci, eecj, epmcoarse, esoil, eseas, xlm, amu, dgnuc, &
                       dgacc, dgcor, pmassn, pmassa, pmassc, pdensn, pdensa, pdensc, knnuc, knacc, &
                       kncor, fconcn, fconca, fconcn_org, fconca_org, dmdt, dndt, cgrn3, cgra3, urn00, &
                       ura00, brna01, c30, deltaso4a, igrid, jgrid, kgrid, brrto)

!USE module_configure, only: grid_config_rec_type
!TYPE (grid_config_rec_type), INTENT (in) :: config_flags

!     IMPLICIT NONE
! dimension of arrays
      INTEGER blksize
! number of species in CBLK
      INTEGER nspcsda
! actual number of cells in arrays
      INTEGER numcells
! number of k-level
      INTEGER layer
! of organic aerosol precursor
      INTEGER ldrog_vbs
      REAL cblk(blksize, nspcsda) ! main array of variables (INPUT a

      REAL dt
! *** Meteorological information:

! synchronization time  [s]
      REAL blkta(blksize) ! Air temperature [ K ]
      REAL blkprs(blksize) ! Air pressure in [ Pa ]
      REAL blkdens(blksize) ! Air density  [ kg/ m**3 ]
      REAL blkrh(blksize)
! *** Chemical production rates: [ ug / m**3 s ]

! Fractional relative humidity
      REAL so4rat(blksize)
! sulfate gas-phase production rate
! total # of cond. vapors & SOA species
      INTEGER ncv
      INTEGER nacv
!bs * organic condensable vapor production rate
! # of anthrop. cond. vapors & SOA speci
      REAL drog(blksize, ldrog_vbs) !bs
! *** anthropogenic organic aerosol mass production rates from aromatics
! Delta ROG conc. [ppm]
      REAL organt1rat(blksize)

! *** anthropogenic organic aerosol mass production rates from aromatics
      REAL organt2rat(blksize)

! *** anthropogenic organic aerosol mass production rates from alkanes &
      REAL organt3rat(blksize)

! *** anthropogenic organic aerosol mass production rates from alkenes &
      REAL organt4rat(blksize)

! *** biogenic organic aerosol production rates
      REAL orgbio1rat(blksize)

! *** biogenic organic aerosol production rates
      REAL orgbio2rat(blksize)

! *** biogenic organic aerosol production rates
      REAL orgbio3rat(blksize)

! *** biogenic organic aerosol production rates
      REAL orgbio4rat(blksize)

! *** Primary emissions rates: [ ug / m**3 s ]
! *** emissions rates for unidentified PM2.5 mass
      REAL epm25i(blksize) ! Aitken mode
      REAL epm25j(blksize)
! *** emissions rates for primary organic aerosol
! Accumululaton mode
      REAL eorgi(blksize) ! Aitken mode
      REAL eorgj(blksize)
! *** emissions rates for elemental carbon
! Accumululaton mode
      REAL eeci(blksize) ! Aitken mode
      REAL eecj(blksize)
! *** emissions rates for coarse mode particles
! Accumululaton mode
      REAL esoil(blksize) ! soil derived coarse aerosols
      REAL eseas(blksize) ! marine coarse aerosols
      REAL epmcoarse(blksize)

! *** OUTPUT:
! *** atmospheric properties
! anthropogenic coarse aerosols
      REAL xlm(blksize) ! atmospheric mean free path [ m ]
      REAL amu(blksize)
! *** modal diameters: [ m ]

! atmospheric dynamic viscosity [ kg
      REAL dgnuc(blksize) ! nuclei mode geometric mean diamete
      REAL dgacc(blksize) ! accumulation geometric mean diamet
      REAL dgcor(blksize)

! *** aerosol properties:
! *** Modal mass concentrations [ ug m**3 ]
! coarse mode geometric mean diamete
      REAL pmassn(blksize) ! mass concentration in Aitken mode
      REAL pmassa(blksize) ! mass concentration in accumulation
      REAL pmassc(blksize)
! *** average modal particle densities  [ kg/m**3 ]

! mass concentration in coarse mode
      REAL pdensn(blksize) ! average particle density in nuclei
      REAL pdensa(blksize) ! average particle density in accumu
      REAL pdensc(blksize)
! *** average modal Knudsen numbers

! average particle density in coarse
      REAL knnuc(blksize) ! nuclei mode  Knudsen number
      REAL knacc(blksize) ! accumulation Knudsen number
      REAL kncor(blksize)
! ***  modal condensation factors ( see comments in NUCLCOND )

! coarse mode  Knudsen number
      REAL fconcn(blksize)
      REAL fconca(blksize)
!bs
      REAL fconcn_org(blksize)
      REAL fconca_org(blksize)
!bs

! *** Rates for secondary particle formation:

! *** production of new mass concentration [ ug/m**3 s ]
      REAL dmdt(blksize) !                                 by particle formation

! *** production of new number concentration [ number/m**3 s ]

! rate of production of new mass concen
      REAL dndt(blksize) !                                 by particle formation

! *** growth rate for third moment by condensation of precursor
!      vapor on existing particles [ 3rd mom/m**3 s ]

! rate of producton of new particle num
      REAL cgrn3(blksize) !  Aitken mode
      REAL cgra3(blksize)
! *** Rates for coaglulation: [ m**3/s ]

! *** Unimodal Rates:

!  Accumulation mode
      REAL urn00(blksize) ! Aitken mode 0th moment self-coagulation ra
      REAL ura00(blksize)

! *** Bimodal Rates:  Aitken mode with accumulation mode ( d( Aitken mod

! accumulation mode 0th moment self-coagulat
      REAL brna01(blksize)
! *** 3rd moment intermodal transfer rate replaces coagulation rate ( FS
! rate for 0th moment
      REAL c30(blksize) ! by intermodal c
      REAL brrto

! *** other processes

! intermodal 3rd moment transfer r
      REAL deltaso4a(blksize) !                                  sulfate aerosol by condensation   [ u

!      INTEGER NN, VV ! loop indICES
! increment of concentration added to

      RETURN
   END SUBROUTINE aeroproc
!//////////////////////////////////////////////////////////////////

! *** Time stepping code advances the aerosol moments one timestep;
   SUBROUTINE aerostep(layer, blksize, nspcsda, numcells, cblk, dt, so4rat &
                       , organt1rat, organt2rat, organt3rat, organt4rat, orgbio1rat, orgbio2rat &
                       , orgbio3rat, orgbio4rat, epm25i, epm25j, eorgi, eorgj, eeci, eecj, esoil, eseas &
                       , epmcoarse, dgnuc, dgacc, fconcn, fconca, fconcn_org, fconca_org, pmassn &
                       , pmassa, pmassc, dmdt, dndt, deltaso4a, urn00, ura00, brna01, c30, cgrn3, cgra3, &
                       igrid, jgrid, kgrid)

! dimension of arrays
      INTEGER blksize
! actual number of cells in arrays
      INTEGER numcells
! nmber of species in CBLK
      INTEGER nspcsda
! model layer
      INTEGER layer
      REAL cblk(blksize, nspcsda) ! main array of variables
      INTEGER igrid, jgrid, kgrid
      REAL dt
! *** Chemical production rates: [ ug / m**3 s ]

! time step [sec]
      REAL so4rat(blksize) ! sulfate gas-phase production rate

! anthropogenic organic aerosol mass production rates
      REAL organt1rat(blksize)
      REAL organt2rat(blksize)
      REAL organt3rat(blksize)
      REAL organt4rat(blksize)

! biogenic organic aerosol production rates
      REAL orgbio1rat(blksize)
      REAL orgbio2rat(blksize)
      REAL orgbio3rat(blksize)
      REAL orgbio4rat(blksize)

! *** Primary emissions rates: [ ug / m**3 s ]
! *** emissions rates for unidentified PM2.5 mass
      REAL epm25i(blksize) ! Aitken mode
      REAL epm25j(blksize)
! *** emissions rates for primary organic aerosol
! Accumululaton mode
      REAL eorgi(blksize) ! Aitken mode
      REAL eorgj(blksize)
! *** emissions rates for elemental carbon
! Accumululaton mode
      REAL eeci(blksize) ! Aitken mode
      REAL eecj(blksize)
! *** emissions rates for coarse mode particles
! Accumululaton mode
      REAL esoil(blksize) ! soil derived coarse aerosols
      REAL eseas(blksize) ! marine coarse aerosols
      REAL epmcoarse(blksize)
! anthropogenic coarse aerosols
      REAL dgnuc(blksize) ! nuclei mode mean diameter [ m ]
      REAL dgacc(blksize)
! accumulation
      REAL fconcn(blksize) ! Aitken mode  [ 1 / s ]
! reciprocal condensation rate
      REAL fconca(blksize) ! acclumulation mode [ 1 / s ]
! reciprocal condensation rate
      REAL fconcn_org(blksize) ! Aitken mode  [ 1 / s ]
! reciprocal condensation rate for organ
      REAL fconca_org(blksize) ! acclumulation mode [ 1 / s ]
! reciprocal condensation rate for organ
      REAL dmdt(blksize) ! by particle formation [ ug/m**3 /s ]
! rate of production of new mass concent
      REAL dndt(blksize) ! by particle formation [ number/m**3 /s
! rate of producton of new particle numb
      REAL deltaso4a(blksize) ! sulfate aerosol by condensation [ ug/m
! increment of concentration added to
      REAL urn00(blksize) ! Aitken intramodal coagulation rate
      REAL ura00(blksize) ! Accumulation mode intramodal coagulati
      REAL brna01(blksize) ! bimodal coagulation rate for number
      REAL c30(blksize) ! by intermodal coagulation
! intermodal 3rd moment transfer rate by
      REAL cgrn3(blksize) ! growth rate for 3rd moment for Aitken
      REAL cgra3(blksize)
! *** Modal mass concentrations [ ug m**3 ]

! growth rate for 3rd moment for Accumul
      REAL pmassn(blksize) ! mass concentration in Aitken mode
      REAL pmassa(blksize) ! mass concentration in accumulation
      REAL pmassc(blksize)

! *** Local Variables

! mass concentration in coarse mode
      INTEGER l, lcell, spc
! ** following scratch variables are used for solvers

! *** variables needed for modal dynamics solvers:
! Loop indices
      REAL*8 a, b, c
      REAL*8 m1, m2, y0, y
      REAL*8 dhat, p, pexpdt, expdt
      REAL*8 loss, prod, pol, lossinv
! mass intermodal transfer by coagulation
      REAL mstrnsfr

      REAL factrans

! *** CODE additions for renaming
      REAL getaf2
      REAL aaa, xnum, xm3, fnum, fm3, phnum, phm3 ! Defined below
      REAL erf, & ! Error and complementary error function
         erfc

      REAL xx
! dummy argument for ERF and ERFC
! a numerical value for a minimum concentration

! *** This value is smaller than any reported tropospheric concentration

! *** Statement function given for error function. Source is
!     Meng, Z., and J.H.Seinfeld (1994) On the source of the submicromet
!      droplet mode of urban and regional aerosols. Aerosol Sci. and Tec
!      20:253-265. They cite Reasearch & Education Asociation (REA), (19
!      Handbook of Mathematical, Scientific, and Engineering Formulas,
!      Tables, Functions, Graphs, Transforms: REA, Piscataway, NJ. p. 49

      erf(xx) = sqrt(1.0-exp(-4.0*xx*xx/pirs))
      erfc(xx) = 1.0-erf(xx)
!     ::::::::::::::::::::::::::::::::::::::::

      RETURN
   END SUBROUTINE aerostep
!#######################################################################

   SUBROUTINE awater(irhx, mso4, mnh4, mno3, wh2o)

!     IMPLICIT NONE
      INTEGER irhx, irh
      REAL mso4, mnh4, mno3
      REAL tso4, tnh4, tno3, wh2o, x
      REAL aw, awc
!     REAL poly4, poly6
      REAL mfs0, mfs1, mfs15, mfs2
      REAL c0(4), c1(4), c15(4), c2(4)
      REAL y, y0, y1, y15, y2, y3, y40, y140, y1540, yc
      REAL kso4(6), kno3(6), mfsso4, mfsno3
      REAL mwso4, mwnh4, mwno3, mw2, mwano3

! *** molecular weights:
      PARAMETER(mwso4=96.0636, mwnh4=18.0985, mwno3=62.0049, &
                mw2=mwso4 + 2.0*mwnh4, mwano3=mwno3 + mwnh4)

!     The polynomials use data for aw as a function of mfs from Tang and
!     Munkelwitz, JGR 99: 18801-18808, 1994.
!     The polynomials were fit to Tang's values of water activity as a
!     function of mfs.

! *** coefficients of polynomials fit to Tang and Munkelwitz data
!     now give mfs as a function of water activity.

      DATA c1/0.9995178, -0.7952896, 0.99683673, -1.143874/
      DATA c15/1.697092, -4.045936, 5.833688, -3.463783/
      DATA c2/2.085067, -6.024139, 8.967967, -5.002934/

! *** the following coefficients are a fit to the data in Table 1 of
!     Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
!      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
! *** New data fit to data from
!       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
!       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
!       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
      DATA c0/0.798079, -1.574367, 2.536686, -1.735297/

! *** polynomials for ammonium nitrate and ammonium sulfate are from:
!     Chan et al.1992, Atmospheric Environment (26A): 1661-1673.

      DATA kno3/0.2906, 6.83665, -26.9093, 46.6983, -38.803, 11.8837/
      DATA kso4/2.27515, -11.147, 36.3369, -64.2134, 56.8341, -20.0953/

      RETURN
   END SUBROUTINE awater
!//////////////////////////////////////////////////////////////////////

   SUBROUTINE coagrate(blksize, nspcsda, numcells, cblk, blkta, pdensn, pdensa, amu, &
                       dgnuc, dgacc, knnuc, knacc, urn00, ura00, brna01, c30)
!***********************************************************************
!**    DESCRIPTION:  calculates aerosol coagulation rates for unimodal
!       and bimodal coagulation using E. Whitby 1990's prescription.

!.......   Rates for coaglulation:
!.......   Unimodal Rates:
!.......   URN00:  nuclei       mode 0th moment self-coagulation rate
!.......   URA00:  accumulation mode 0th moment self-coagulation rate

!.......   Bimodal Rates:  (only 1st order coeffs appear)
!.......   NA-- nuclei  with accumulation coagulation rates,
!.......   AN-- accumulation with nuclei coagulation rates
!.......   BRNA01:  rate for 0th moment ( d(nuclei mode 0) / dt  term)
!.......   BRNA31:           3rd        ( d(nuclei mode 3) / dt  term)
!**    Revision history:
!       prototype 1/95 by Uma and Carlie
!       Revised   8/95 by US for calculation of density from stmt func
!                 and collect met variable stmt funcs in one include fil
!      REVISED 7/25/96 by FSB to use block structure
!      REVISED 9/13/96 BY FSB for Uma's FIXEDBOTH case only.
!      REVISED 11/08/96 BY FSB the Whitby Shankar convention on signs
!                              changed. All coagulation coefficients
!                              returned with positive signs. Their
!                              linearization is also abandoned.
!                              Fixed values are used for the corrections
!                              to the free-molecular coagulation integra
!                              The code forces the harmonic means to be
!                              evaluated in 64 bit arithmetic on 32 bit
!     REVISED 11/14/96 BY FSB  Internal units are now MKS, moment / unit

!      REVISED 1/12/98 by FSB   C30 replaces BRNA31 as an array. This wa
!                              because BRNA31 can become zero on a works
!                              because of limited precision. With the ch
!                              aerostep to omit update of the 3rd moment
!                              C30 is the only variable now needed.
!                              the logic using ONE88 to force REAL*8 ari
!                              has been removed and all intermediates ar
!                              REAL*8.
!     IMPLICIT NONE

! dimension of arrays
      INTEGER blksize
! actual number of cells in arrays
      INTEGER numcells

      INTEGER nspcsda

! nmber of species in CBLK
      REAL cblk(blksize, nspcsda) ! main array of variables
      REAL blkta(blksize) ! Air temperature [ K ]
      REAL pdensn(blksize) ! average particel density in Aitk
      REAL pdensa(blksize) ! average particel density in accu
      REAL amu(blksize) ! atmospheric dynamic viscosity [
      REAL dgnuc(blksize) ! Aitken mode mean diameter [ m ]
      REAL dgacc(blksize) ! accumulation mode mean diameter
      REAL knnuc(blksize) ! Aitken mode Knudsen number
      REAL knacc(blksize)
! *** output:

! accumulation mode Knudsen number
      REAL urn00(blksize) ! intramodal coagulation rate (Ait
      REAL ura00(blksize)
! intramodal coagulation rate (acc
      REAL brna01(blksize) ! intermodal coagulaton rate (numb
      REAL c30(blksize) ! by inter

! *** Local variables:
! intermodal 3rd moment transfer r
      REAL*8 kncnuc, & ! coeffs for unimodal NC coag rate
         kncacc
      REAL*8 kfmnuc, & ! coeffs for unimodal FM coag rate
         kfmacc
      REAL*8 knc, & ! coeffs for bimodal NC, FM coag rate
         kfm
      REAL*8 bencnn, & ! NC 0th moment coag rate (both modes)
         bencna
      REAL*8 & ! NC 3rd moment coag rate (nuc mode)
         bencm3n
      REAL*8 befmnn, & ! FM 0th moment coag rate (both modes)
         befmna
      REAL*8 & ! FM 3rd moment coag rate (nuc mode)
         befm3n
      REAL*8 betann, & ! composite coag rates, mom 0 (both mode
         betana
      REAL*8 & ! intermodal coagulation rate for 3rd mo
         brna31
      REAL*8 & ! scratch subexpression
         s1
      REAL*8 t1, & ! scratch subexpressions
         t2
      REAL*8 t16, & ! T1**6, T2**6
         t26
      REAL*8 rat, & ! ratio of acc to nuc size and its inver
         rin
      REAL*8 rsqt, & ! sqrt( rat ), rsqt**4
         rsq4
      REAL*8 rsqti, & ! sqrt( 1/rat ), sqrt( 1/rat**3 )
         rsqi3
      REAL*8 & ! dgnuc**3
         dgn3
      REAL*8 & !                                 in 64 bit arithmetic
         dga3
! dgacc**3

      INTEGER lcell
! *** Fixed values for correctionss to coagulation
!      integrals for free-molecular case.
! loop counter
      REAL*8 bm0
      PARAMETER(bm0=0.8D0)
      REAL*8 bm0i
      PARAMETER(bm0i=0.9D0)
      REAL*8 bm3i
      PARAMETER(bm3i=0.9D0)
      REAL*8 & ! approx Cunningham corr. factor
         a
      PARAMETER(a=1.246D0)
!.......................................................................
!   begin body of subroutine  COAGRATE

!...........   Main computational grid-traversal loops
!...........   for computing coagulation rates.

! *** Both modes have fixed std devs.
      DO lcell = 1, &
                 numcells
! *** moment independent factors

!  loop on LCELL
         s1 = two3*boltz*blkta(lcell)/amu(lcell)

! For unimodal coagualtion:

         kncnuc = s1
         kncacc = s1

         kfmnuc = sqrt(3.0*boltz*blkta(lcell)/pdensn(lcell))
         kfmacc = sqrt(3.0*boltz*blkta(lcell)/pdensa(lcell))

! For bimodal coagulation:

         knc = s1
         kfm = sqrt(6.0*boltz*blkta(lcell)/(pdensn(lcell) + pdensa(lcell)))

!...........   Begin unimodal coagulation rate calculations:
!...........   Near-continuum regime.

         dgn3 = dgnuc(lcell)**3
         dga3 = dgacc(lcell)**3

         t1 = sqrt(dgnuc(lcell))
         t2 = sqrt(dgacc(lcell))
         t16 = & ! = T1**6
            dgn3
         t26 = &
            dga3
!.......   Note rationalization of fractions and subsequent cancellation
!.......   from the formulation in  Whitby et al. (1990)

! = T2**6
         bencnn = kncnuc*(1.0+esn08+a*knnuc(lcell)*(esn04+esn20))

         bencna = kncacc*(1.0+esa08+a*knacc(lcell)*(esa04+esa20))

!...........   Free molecular regime. Uses fixed value for correction
!               factor BM0

         befmnn = kfmnuc*t1*(en1 + esn25+2.0*esn05)*bm0
         befmna = kfmacc*t2*(ea1 + esa25+2.0*esa05)*bm0

!...........   Calculate half the harmonic mean between unimodal rates
!...........   free molecular and near-continuum regimes

! FSB       64 bit evaluation

         betann = bencnn*befmnn/(bencnn + befmnn)
         betana = bencna*befmna/(bencna + befmna)

         urn00(lcell) = betann
         ura00(lcell) = betana

! *** End of unimodal coagulation calculations.

!...........   Begin bimodal coagulation rate calculations:

         rat = dgacc(lcell)/dgnuc(lcell)
         rin = 1.0D0/rat
         rsqt = sqrt(rat)
         rsq4 = rat**2

         rsqti = 1.0D0/rsqt
         rsqi3 = rin*rsqti

!...........   Near-continuum coeffs:
!...........   0th moment nuc mode bimodal coag coefficient

         bencnn = knc*(2.0+a*knnuc(lcell)*(esn04+rat*esn16*esa04) + a*knacc(lcell &
                                                                            )*(esa04+rin*esa16*esn04) + (rat + rin)*esn04*esa04)

!...........   3rd moment nuc mode bimodal coag coefficient

         bencm3n = knc*dgn3*(2.0*esn36+a*knnuc(lcell)*(esn16+rat*esn04*esa04) + a &
                             *knacc(lcell)*(esn36*esa04+rin*esn64*esa16) + rat*esn16*esa04+ &
                             rin*esn64*esa04)

!...........   Free molecular regime coefficients:
!...........   Uses fixed value for correction
!               factor BM0I, BM3I
!...........   0th moment nuc mode coeff

         befmnn = kfm*bm0i*t1*(en1 + rsqt*ea1 + 2.0*rat*en1*esa04+rsq4*esn09*esa16+ &
                               rsqi3*esn16*esa09+2.0*rsqti*esn04*ea1)

!...........   3rd moment nuc mode coeff

         befm3n = kfm*bm3i*t1*t16*(esn49+rsqt*esn36*ea1 + 2.0*rat*esn25*esa04+ &
                                   rsq4*esn09*esa16+rsqi3*esn100*esa09+2.0*rsqti*esn64*ea1)

!...........   Calculate half the harmonic mean between bimodal rates
!...........   free molecular and near-continuum regimes

! FSB       Force 64 bit evaluation

         brna01(lcell) = bencnn*befmnn/(bencnn + befmnn)

         brna31 = bencm3n* & ! BRNA31 now is a scala
                  befm3n/(bencm3n+befm3n)
         c30(lcell) = brna31*cblk(lcell, vac0)*cblk(lcell, vnu0)
!       print *,c30(lcell),brna31,cblk(lcell,vac0),cblk(lcell,vnu0)
         ! 3d moment transfer by intermodal coagula
!         End bimodal coagulation rate.

      END DO
! end of main lop over cells
      RETURN
   END SUBROUTINE coagrate
!------------------------------------------------------------------

! subroutine  to find the roots of a cubic equation / 3rd order polynomi
! formulae can be found in numer. recip.  on page 145
!   kiran  developed  this version on 25/4/1990
!   dr. francis binkowski modified the routine on 6/24/91, 8/7/97
! ***
!234567
! coagrate
   SUBROUTINE cubic(a2, a1, a0, nr, crutes)
!     IMPLICIT NONE
      INTEGER nr
      REAL*8 a2, a1, a0
      REAL crutes(3)
      REAL*8 qq, rr, a2sq, theta, sqrt3, one3rd
      REAL*8 dum1, dum2, part1, part2, part3, rrsq, phi, yy1, yy2, yy3
      REAL*8 costh, sinth
      DATA sqrt3/1.732050808/, one3rd/0.333333333/
!bs
      REAL*8 onebs
      PARAMETER(onebs=1.0)
!bs
      a2sq = a2*a2
      qq = (a2sq - 3.*a1)/9.
      rr = (a2*(2.*a2sq - 9.*a1) + 27.*a0)/54.
! CASE 1 THREE REAL ROOTS or  CASE 2 ONLY ONE REAL ROOT
      dum1 = qq*qq*qq
      rrsq = rr*rr
      dum2 = dum1 - rrsq
      IF (dum2 >= 0.) THEN
! NOW WE HAVE THREE REAL ROOTS
         phi = sqrt(dum1)
         IF (abs(phi) < 1.E-20) THEN
            print *, ' cubic phi small, phi = ', phi
            crutes(1) = 0.0
            crutes(2) = 0.0
            crutes(3) = 0.0
            nr = 0
            CALL wrf_error_fatal('PHI < CRITICAL VALUE')
         END IF
         theta = acos(rr/phi)/3.0
         costh = cos(theta)
         sinth = sin(theta)
! *** use trig identities to simplify the expressions
! *** binkowski's modification
         part1 = sqrt(qq)
         yy1 = part1*costh
         yy2 = yy1 - a2/3.0
         yy3 = sqrt3*part1*sinth
         crutes(3) = -2.0*yy1 - a2/3.0
         crutes(2) = yy2 + yy3
         crutes(1) = yy2 - yy3
! *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
         IF (crutes(1) < 0.0) crutes(1) = 1.0E9
         IF (crutes(2) < 0.0) crutes(2) = 1.0E9
         IF (crutes(3) < 0.0) crutes(3) = 1.0E9
! *** put smallest positive root in crutes(1)
         crutes(1) = min(crutes(1), crutes(2), crutes(3))
         nr = 3
!     NOW HERE WE HAVE ONLY ONE REAL ROOT
      ELSE
! dum IS NEGATIVE
         part1 = sqrt(rrsq - dum1)
         part2 = abs(rr)
         part3 = (part1 + part2)**one3rd
         crutes(1) = -sign(onebs, rr)*(part3 + (qq/part3)) - a2/3.
!bs     &        -sign(1.0,rr) * ( part3 + (qq/part3) ) - a2/3.
         crutes(2) = 0.
         crutes(3) = 0.
!IAREV02...ADDITIONAL CHECK on NEGATIVE ROOTS
! *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
!     if(crutes(1) .lt. 0.0) crutes(1) = 1.0e9
         nr = 1
      END IF
      RETURN
   END SUBROUTINE cubic
!///////////////////////////////////////////////////////////////////////

!    Calculate the aerosol chemical speciation and water content.
! cubic
   SUBROUTINE eql3(blksize, nspcsda, numcells, cblk, blkta, blkrh, igrid, jgrid, kgrid)
!***********************************************************************
!**    DESCRIPTION:
!        Calculates the distribution of ammonia/ammonium, nitric acid/nitrate,
!        and water between the gas and aerosol phases as the total sulfate,
!        ammonia, and nitrate concentrations, relative humidity and
!        temperature change.  The evolution of the aerosol mass concentration
!        due to the change in aerosol chemical composition is calculated.
!**    REVISION HISTORY:
!       prototype 1/95 by Uma and Carlie
!       Revised   8/95 by US to calculate air density in stmt func
!                 and collect met variable stmt funcs in one include fil
!       Revised 7/26/96 by FSB to use block concept.
!       Revise 12/1896 to do do i-mode calculation.
!**********************************************************************
!     IMPLICIT NONE

! dimension of arrays
      INTEGER blksize
! actual number of cells in arrays
      INTEGER numcells
! nmber of species in CBLK
      INTEGER nspcsda, igrid, jgrid, kgrid
      REAL cblk(blksize, nspcsda)
! *** Meteorological information in blocked arays:

! main array of variables
      REAL blkta(blksize) ! Air temperature [ K ]
      REAL blkrh(blksize)

! Fractional relative humidity

      INTEGER lcell
! loop counter
! air temperature
      REAL temp
!iamodels3
      REAL rh
! relative humidity
      REAL so4, no3, nh3, nh4, hno3
      REAL aso4, ano3, ah2o, anh4, gnh3, gno3
! Fraction of dry sulfate mass in i-mode
      REAL fraci
!.......................................................................
      REAL fracj
!
! ISOROPIA variables double precision
!
      real(kind=8) wi(5), wt(5), wt_save(5)
      real(kind=8) rhi, tempi, cntrl(2)
      real(kind=8) gas(3), aerliq(12), aersld(9), other(6)
      character*15 scasi

!      WRITE(20,*) ' IN EQL 3 '

! Fraction of dry sulfate mass in j-mode
      DO lcell = 1, &
                 numcells
! *** Fetch temperature, fractional relative humidity, and
!     air density

!  loop on cells
         temp = blkta(lcell)
         rh = blkrh(lcell)

         rhi = amin1(rh, 0.995)
         tempi = temp
         cntrl(1) = 0.d0 ! 0 = forward problem
         cntrl(2) = 0.d0 ! 0 = solids and liquid allowed

         wi(1) = (cblk(lcell, vnaaj) + cblk(lcell, vnaai))/mw_na_aer*1.e-6 ! sodium

         wi(2) = (cblk(lcell, vsulf)/(mw_so4_aer + 2.) + &
                  (cblk(lcell, vso4aj) + cblk(lcell, vso4ai))/mw_so4_aer)*1.e-6 ! sulfate

         wi(3) = (cblk(lcell, vnh3)/(mw_nh4_aer - 1.) + &
                  (cblk(lcell, vnh4aj) + cblk(lcell, vnh4ai))/mw_nh4_aer)*1.e-6 ! ammoinum

         wi(4) = (cblk(lcell, vhno3)/(mw_no3_aer + 1.) + &
                  (cblk(lcell, vno3aj) + cblk(lcell, vno3ai))/mw_no3_aer)*1.e-6 ! nitrate

         wi(5) = (cblk(lcell, vhcl)/(mw_cl_aer - 1.) + &
                  (cblk(lcell, vclaj) + cblk(lcell, vclai))/mw_cl_aer)*1.e-6 ! chloride

         wt_save(1) = wi(1) ! sodium
         wt_save(2) = wi(2) ! sulfate
         wt_save(3) = wi(3) ! ammoinum
         wt_save(4) = wi(4) ! nitrate
         wt_save(5) = wi(5) ! chloride
         if (igrid .eq. 28 .and. jgrid .eq. 24 .and. kgrid .eq. 1) then
            print *, vhcl, vclai
            print *, wi(1), wi(2), wi(3), wi(4), wi(5)
         endif

!        call isoropia(wi,rhi,tempi,cntrl,wt,gas,aerliq,aersld,scasi,other)

! *** the following is an interim procedure. Assume the i-mode has the
!     same relative mass concentrations as the total mass. Use SO4 as
!     the surrogate.

! *** update gas / vapor phase

         cblk(lcell, vnh3) = gas(1)*1.e6*17.
         cblk(lcell, vhno3) = gas(2)*1.e6*63.
         cblk(lcell, vhcl) = gas(3)*1.e6*36.
!        if(igrid.eq.28.and.jgrid.eq.24.and.kgrid.eq.1)then
!         print *,vhcl,vnh3,vhno3
!         print *,cblk(lcell,vnh3),cblk(lcell,vhno3),cblk(lcell,vhcl)
!        endif

! *** get modal fraction
         fraci = cblk(lcell, vso4ai)/(cblk(lcell, vso4aj) + cblk(lcell, vso4ai))
         fracj = 1.0-fraci

! *** update do i-mode
         cblk(lcell, vh2oai) = fraci*aerliq(8)*18.*1.e6
         cblk(lcell, vnh4ai) = fraci*(wt_save(3) - gas(1))
         cblk(lcell, vno3ai) = fraci*(wt_save(4) - gas(2))
         cblk(lcell, vclai) = fraci*(wt_save(5) - gas(3))
         cblk(lcell, vnaai) = fraci*wi(1)

! *** update accumulation mode:
         cblk(lcell, vh2oaj) = fracj*aerliq(8)*18.*1.e6
         cblk(lcell, vnh4aj) = fracj*(wt_save(3) - gas(1))
         cblk(lcell, vno3aj) = fracj*(wt_save(4) - gas(2))
         cblk(lcell, vclaj) = fracj*(wt_save(5) - gas(3))
         cblk(lcell, vnaaj) = fracj*wi(1)
!        if(igrid.eq.28.and.jgrid.eq.24.and.kgrid.eq.1)then
!         print *,vh2oaj,vnh4aj,vno3aj,vclaj,vnaaj
!         print *,cblk(lcell,vnh4aj),cblk(lcell,vno3aj),cblk(lcell,vclaj),aerliq(8)
!        endif

      END DO
!  end loop on cells

      RETURN

!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
   END SUBROUTINE eql3
! eql3
!    Calculate the aerosol chemical speciation and water content.

! cubic
   SUBROUTINE eql4(blksize, nspcsda, numcells, cblk, blkta, blkrh)
!***********************************************************************
!**    DESCRIPTION:
!       Calculates the distribution of ammonia/ammonium, nitric acid/nitrate,
!       and water between the gas and aerosol phases as the total sulfate,
!       ammonia, and nitrate concentrations, relative humidity and
!       temperature change.  The evolution of the aerosol mass concentration
!       due to the change in aerosol chemical composition is calculated.
!**    REVISION HISTORY:
!       prototype 1/95 by Uma and Carlie
!       Revised   8/95 by US to calculate air density in stmt func
!                 and collect met variable stmt funcs in one include fil
!       Revised 7/26/96 by FSB to use block concept.
!       Revise 12/1896 to do do i-mode calculation.
!**********************************************************************
!     IMPLICIT NONE

! dimension of arrays
      INTEGER blksize
! actual number of cells in arrays
      INTEGER numcells
! nmber of species in CBLK
      INTEGER nspcsda
      REAL cblk(blksize, nspcsda)
! *** Meteorological information in blocked arays:

! main array of variables
      REAL blkta(blksize) ! Air temperature [ K ]
      REAL blkrh(blksize)

! Fractional relative humidity

      INTEGER lcell
! loop counter
! air temperature
      REAL temp
!iamodels3
      REAL rh
! relative humidity
      REAL so4, no3, nh3, nh4, hno3
      REAL aso4, ano3, ah2o, anh4, gnh3, gno3
! Fraction of dry sulfate mass in i-mode
      REAL fraci
!.......................................................................
      REAL fracj
! Fraction of dry sulfate mass in j-mode
      DO lcell = 1, &
                 numcells
! *** Fetch temperature, fractional relative humidity, and
!     air density

!  loop on cells
         temp = blkta(lcell)
         rh = blkrh(lcell)

! *** the following is an interim procedure. Assume the i-mode has the
!     same relative mass concentrations as the total mass. Use SO4 as
!     the surrogate. The results of this should be the same as those
!     from the original RPM.

! *** do total aerosol
         so4 = cblk(lcell, vso4aj) + cblk(lcell, vso4ai)

!iamodels3
         no3 = cblk(lcell, vno3aj) + cblk(lcell, vno3ai)
!    &                        + CBLK(LCELL, VHNO3)

         hno3 = cblk(lcell, vhno3)

!iamodels3

         nh3 = cblk(lcell, vnh3)

         nh4 = cblk(lcell, vnh4aj) + cblk(lcell, vnh4ai)
!    &                        + CBLK(LCELL, VNH3)

!bs           CALL rpmares(SO4,HNO3,NO3,NH3,NH4,RH,TEMP,
!bs     &             ASO4,ANO3,AH2O,ANH4,GNH3,GNO3)
!bs
!bs * call old version of rpmares
!bs
         CALL rpmares_old(so4, hno3, no3, nh3, nh4, rh, temp, aso4, ano3, ah2o, anh4, &
                          gnh3, gno3)
!bs

! *** get modal fraction
         fraci = cblk(lcell, vso4ai)/(cblk(lcell, vso4aj) + cblk(lcell, vso4ai))
         fracj = 1.0-fraci

! *** update do i-mode

         cblk(lcell, vh2oai) = fraci*ah2o
         cblk(lcell, vnh4ai) = fraci*anh4
         cblk(lcell, vno3ai) = fraci*ano3

! *** update accumulation mode:

         cblk(lcell, vh2oaj) = fracj*ah2o
         cblk(lcell, vnh4aj) = fracj*anh4
         cblk(lcell, vno3aj) = fracj*ano3

! *** update gas / vapor phase
         cblk(lcell, vnh3) = gnh3
         cblk(lcell, vhno3) = gno3

      END DO
!  end loop on cells
      RETURN

!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
   END SUBROUTINE eql4
! eql4

   SUBROUTINE fdjac(n, x, fjac, ct, cs, imw)
!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
!bs                                                                    !
!bs  Description:                                                      !
!bs                                                                    !
!bs  Get the Jacobian of the function                                  !
!bs                                                                    !
!bs         ( a1 * X1^2 + b1 * X1 + c1 )                               !
!bs         ( a2 * X2^2 + b2 * X1 + c2 )                               !
!bs         ( a3 * X3^2 + b3 * X1 + c3 )                               !
!bs  F(X) = ( a4 * X4^2 + b4 * X1 + c4 ) = 0.                          !
!bs         ( a5 * X5^2 + b5 * X1 + c5 )                               !
!bs         ( a6 * X6^2 + b6 * X1 + c6 )                               !
!bs                                                                    !
!bs   a_i = IMW_i                                                      !
!bs   b_i = SUM(X_j * IMW_j)_j.NE.i + CSAT_i * IMX_i - CTOT_i * IMW_i  !
!bs   c_i = - CTOT_i * [ SUM(X_j * IMW_j)_j.NE.i + M ]                 !
!bs                                                                    !
!bs          delta F_i    ( 2. * a_i * X_i + b_i           if i .EQ. j !
!bs  J_ij = ----------- = (                                            !
!bs          delta X_j    ( X_i * IMW_j - CTOT_i * IMW_j   if i .NE. j !
!bs                                                                    !
!bs                                                                    !
!bs  Called by:       NEWT                                             !
!bs                                                                    !
!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
!bs
!     IMPLICIT NONE
!bs
!bs
!dimension of problem
      INTEGER n
      REAL x(n) !bs
!     INTEGER NP                !bs maximum expected value of N
!     PARAMETER (NP = 6)
!bs initial guess of CAER
      REAL ct(np)
      REAL cs(np)
      REAL imw(np)
!bs
      REAL fjac(n, n)
!bs
      INTEGER i, & !bs loop index
         j
      REAL a(np)
      REAL b(np)
      REAL b1(np)
      REAL b2(np)
      REAL sum_jnei
!bs
      DO i = 1, n
         a(i) = imw(i)
         sum_jnei = 0.
         DO j = 1, n
            sum_jnei = sum_jnei + x(j)*imw(j)
         END DO
         b1(i) = sum_jnei - (x(i)*imw(i))
         b2(i) = cs(i)*imw(i) - ct(i)*imw(i)
         b(i) = b1(i) + b2(i)
      END DO
      DO j = 1, n
         DO i = 1, n
            IF (i == j) THEN
               fjac(i, j) = 2.*a(i)*x(i) + b(i)
            ELSE
               fjac(i, j) = x(i)*imw(j) - ct(i)*imw(j)
            END IF
         END DO
      END DO
!bs
      RETURN
   END SUBROUTINE fdjac
!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
   FUNCTION fmin(x, fvec, n, ct, cs, imw, m)
!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
!bs                                                                    !
!bs  Description:                                                      !
!bs                                                                    !
!bs  Adopted from Numerical Recipes in FORTRAN, Chapter 9.7, 2nd ed.   !
!bs                                                                    !
!bs  Returns f = 0.5 * F*F at X. SR FUNCV(N,X,F) is a fixed name,      !
!bs  user-supplied routine that returns the vector of functions at X.  !
!bs  The common block NEWTV communicates the function values back to   !
!bs  NEWT.                                                             !
!bs                                                                    !
!bs  Called by:       NEWT                                             !
!bs                                                                    !
!bs  Calls:           FUNCV                                            !
!bs                                                                    !
!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!

!     IMPLICIT NONE

!bs
!bs
      INTEGER n
!     INTEGER NP
!     PARAMETER (NP = 6)
      REAL ct(np)
      REAL cs(np)
      REAL imw(np)
      REAL m, fmin
      REAL x(*), fvec(np)

      INTEGER i
      REAL sum

      CALL funcv(n, x, fvec, ct, cs, imw, m)
      sum = 0.
      DO i = 1, n
         sum = sum + fvec(i)**2
      END DO
      fmin = 0.5*sum
      RETURN
   END FUNCTION fmin
!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
   SUBROUTINE funcv(n, x, fvec, ct, cs, imw, m)
!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
!bs                                                                    !
!bs  Description:                                                      !
!bs                                                                    !
!bs  Called by:       FMIN                                             !
!bs                                                                    !
!bs  Calls:           None                                             !
!bs                                                                    !
!bs ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS ** ** BS ** BS *!
!bs
!     IMPLICIT NONE
!bs
!bs
      INTEGER n
      REAL x(*)
      REAL fvec(n)
!bs
!     INTEGER NP
!     PARAMETER (NP = 6)
      REAL ct(np)
      REAL cs(np)
      REAL imw(np)
      REAL m
!bs
      INTEGER i, j
      REAL sum_jnei
      REAL a(np)
      REAL b(np)
      REAL c(np)
!bs
      DO i = 1, n
         a(i) = imw(i)
         sum_jnei = 0.
         DO j = 1, n
            sum_jnei = sum_jnei + x(j)*imw(j)
         END DO
         sum_jnei = sum_jnei - (x(i)*imw(i))
         b(i) = sum_jnei + cs(i)*imw(i) - ct(i)*imw(i)
         c(i) = -ct(i)*(sum_jnei + m)
         fvec(i) = a(i)*x(i)**2 + b(i)*x(i) + c(i)
      END DO
!bs
      RETURN
   END SUBROUTINE funcv
   REAL FUNCTION getaf(ni, nj, dgni, dgnj, xlsgi, xlsgj, sqrt2)
! *** set up new processor for renaming of particles from i to j modes
!     IMPLICIT NONE
      REAL aa, bb, cc, disc, qq, alfa, l, yji
      REAL ni, nj, dgni, dgnj, xlsgi, xlsgj, sqrt2

      alfa = xlsgi/xlsgj
      yji = log(dgnj/dgni)/(sqrt2*xlsgi)
      aa = 1.0-alfa*alfa
      l = log(alfa*nj/ni)
      bb = 2.0*yji*alfa*alfa
      cc = l - yji*yji*alfa*alfa
      disc = bb*bb - 4.0*aa*cc
      IF (disc < 0.0) THEN
         getaf = - & ! error in intersection
                 5.0
         RETURN
      END IF
      qq = -0.5*(bb + sign(1.0, bb)*sqrt(disc))
      getaf = cc/qq
      RETURN
! *** subroutine to implement Kulmala, Laaksonen, Pirjola
   END FUNCTION getaf
!     Parameterization for sulfuric acid/water
!     nucleation rates, J. Geophys. Research (103), pp 8301-8307,
!     April 20, 1998.

!ia rev01 27.04.99 changes made to calculation of MDOT see RBiV p.2f
!ia rev02 27.04.99 security check on MDOT > SO4RAT

!ia      subroutine klpnuc( Temp, RH, H2SO4,NDOT, MDOT, M2DOT)
! getaf
   SUBROUTINE klpnuc(temp, rh, h2so4, ndot1, mdot1, so4rat)
!     IMPLICIT NONE

! *** Input:

! ambient temperature [ K ]
      REAL temp
! fractional relative humidity
      REAL rh
! sulfuric acid concentration [ ug / m**3 ]
      REAL h2so4

      REAL so4rat
! *** Output:

!sulfuric acid production rate [ ug / ( m**3 s )]
! particle number production rate [ # / ( m**3 s )]
      REAL ndot1
! particle mass production rate [ ug / ( m**3 s )]
      REAL mdot1
      ! [ m**2 / ( m**3 s )]
      REAL m2dot

! *** Internal:

! *** NOTE, all units are cgs internally.
! particle second moment production rate

      REAL ra
! fractional relative acidity
! sulfuric acid vaper concentration [ cm ** -3 ]
      REAL nav
! water vapor concentration   [ cm ** -3 ]
      REAL nwv
! equilibrium sulfuric acid vapor conc. [ cm ** -3 ]
      REAL nav0
      ! to produce a nucleation rate of 1 [ cm ** -3  s ** -1
      REAL nac
! critical sulfuric acid vapor concentration [ cm ** -3
! mole fractio of the critical nucleus
      REAL xal
      REAL nsulf, & ! see usage
         delta
      REAL*8 & ! factor to calculate Jnuc
         chi
      REAL*8 &
         jnuc
! nucleation rate [ cm ** -3  s ** -1 ]
      REAL tt, & ! dummy variables for statement functions
         rr
      REAL pi
      PARAMETER(pi=3.14159265)

      REAL pid6
      PARAMETER(pid6=pi/6.0)

! avogadro's constant [ 1/mol ]
      REAL avo
      PARAMETER(avo=6.0221367E23)

! universal gas constant [ j/mol-k ]
      REAL rgasuniv
      PARAMETER(rgasuniv=8.314510)

! 1 atmosphere in pascals
      REAL atm
      PARAMETER(atm=1013.25E+02)

! formula weight for h2so4 [ g mole **-1 ]
      REAL mwh2so4
      PARAMETER(mwh2so4=98.07948)

! diameter of a 3.5 nm particle in cm
      REAL d35
      PARAMETER(d35=3.5E-07)
      REAL d35sq
      PARAMETER(d35sq=d35*d35)
! volume of a 3.5 nm particle in cm**3
      REAL v35
      PARAMETER(v35=pid6*d35*d35sq)
!ia rev01

      REAL mp
! ***  conversion factors:
! mass of sulfate in a 3.5 nm particle
      ! number per cubic cm.
      REAL ugm3_ncm3
! micrograms per cubic meter to
      PARAMETER(ugm3_ncm3=(avo/mwh2so4)*1.0E-12)
!ia rev01
! molecules to micrograms
      REAL nc_ug
      PARAMETER(nc_ug=(1.0E6)*mwh2so4/avo)

! *** statement functions **************

      REAL pdens, &
         rho_p
! particle density [ g / cm**3]
      REAL ad0, ad1, ad2, &
         ad3
! coefficients for density expression
      PARAMETER(ad0=1.738984, ad1=-1.882301, ad2=2.951849, ad3=-1.810427)
! *** Nair and Vohra, Growth of aqueous sulphuric acid droplets
!     as a function of relative humidity,
!     J. Aerosol Science, 6, pp 265-271, 1975.

!ia rev01

! fit to Nair & Vohra data
      ! the mass of sulfate in a 3.5 nm particle
      REAL mp35
! arithmetic statement function to compute
      REAL a0, a1, a2, & ! coefficients for cubic in mp35
         a3
      PARAMETER(a0=1.961385E2, a1=-5.564447E2, a2=8.828801E2, a3=-5.231409E2)

      REAL ph2so4, & ! for h2so4 and h2o vapor pressures [ Pa ]
         ph2o

! arithmetic statement functions
      pdens(rr) = ad0 + rr*(ad1 + rr*(ad2 + rr*ad3))

      ph2o(tt) = exp(77.34491296-7235.4246512/tt - 8.2*log(tt) + tt*5.7113E-03)

      ph2so4(tt) = exp(27.78492066-10156.0/tt)

! *** both ph2o and ph2so4 are  as in Kulmala et al.  paper

!ia rev01

! *** function for the mass of sulfate in   a 3.5 nm sphere
! *** obtained from a fit to the number of sulfate monomers in
!     a 3.5 nm particle. Uses data from Nair & Vohra
      mp35(rr) = nc_ug*(a0 + rr*(a1 + rr*(a2 + rr*a3)))

! *** begin code:

!     The 1.0e-6 factor in the following converts from MKS to cgs units

! *** get water vapor concentration [ molecles / cm **3 ]

      nwv = rh*ph2o(temp)/(rgasuniv*temp)*avo*1.0E-6

! *** calculate the equilibrium h2so4 vapor concentration.

! *** use Kulmala corrections:

! ***
      nav0 = ph2so4(temp)/(rgasuniv*temp)*avo*1.0E-6

! *** convert sulfuric acid vapor concentration from micrograms
!     per cubic meter to molecules per cubic centimeter.

      nav = ugm3_ncm3*h2so4

! *** calculate critical concentration of sulfuric acid vapor

      nac = exp(-14.5125+0.1335*temp - 10.5462*rh + 1958.4*rh/temp)

! *** calculate relative acidity

      ra = nav/nav0

! *** calculate temperature correction

      delta = 1.0+(temp - 273.15)/273.14

! *** calculate molar fraction

      xal = 1.2233-0.0154*ra/(ra + rh) + 0.0102*log(nav) - 0.0415*log(nwv) + &
            0.0016*temp

! *** calculate Nsulf
      nsulf = log(nav/nac)

! *** calculate particle produtcion rate [ # / cm**3 ]

      chi = 25.1289*nsulf - 4890.8*nsulf/temp - 1743.3/temp - &
            2.2479*delta*nsulf*rh + 7643.4*xal/temp - 1.9712*xal*delta/rh

      jnuc = exp(chi)
! [ # / cm**3 ]
      ndot1 = (1.0E06)*jnuc
!      write(91,*) ' inside klpnuc '
!     write(91,*) ' Jnuc = ', Jnuc
!     write(91,*) ' NDOT = ', NDOT1

! *** calculate particle density

      rho_p = pdens(rh)

!     write(91,*) ' rho_p =', rho_p

! *** get the mass of sulfate in a 3.5 nm particle

      mp = mp35(rh) ! in a 3.5 nm particle at ambient RH

! *** calculate mass production rate [ ug / m**3]
!     assume that the particles are 3.5 nm in diameter.

!     MDOT1 =  (1.0E12) * rho_p * v35 * Jnuc

!ia rev01

! number of micrograms of sulfate
      mdot1 = mp*ndot1

!ia rev02

      IF (mdot1 > so4rat) THEN

         mdot1 = &
            so4rat
! limit nucleated mass by available ma
         ndot1 = mdot1/ &
                 mp
! adjust DNDT to this
      END IF

      IF (mdot1 == 0.) ndot1 = 0.

! *** calculate M2 production rate [ m**2 / (m**3 s)]

      m2dot = 1.0E-04*d35sq*ndot1

      RETURN

   END SUBROUTINE klpnuc
!------------------------------------------------------------------------------

   SUBROUTINE modpar(blksize, nspcsda, numcells, cblk, blkta, blkprs, pmassn, &
                     pmassa, pmassc, pdensn, pdensa, pdensc, xlm, amu, dgnuc, dgacc, dgcor, knnuc, &
                     knacc, kncor)

!**    DESCRIPTION:
!       Calculates modal parameters and derived variables,
!       log-squared of std deviation, mode mean size, Knudsen number)
!       based on current values of moments for the modes.
! FSB   Now calculates the 3rd moment, mass, and density in all 3 modes.
!**
!**    Revision history:
!       Adapted 3/95 by US and CJC from EAM2's MODPAR and INIT3
!       Revised  7/23/96 by FSB to use COMMON blocks and small blocks
!        instead of large 3-d arrays, and to assume a fixed std.
!       Revised 12/06/96 by FSB to include coarse mode
!       Revised 1/10/97 by FSB to have arrays passed in call vector
!**********************************************************************

!     IMPLICIT NONE

!     Includes:

! *** input:

! dimension of arrays
      INTEGER blksize
! actual number of cells in arrays
      INTEGER numcells

      INTEGER nspcsda

! nmber of species in CBLK
      REAL cblk(blksize, nspcsda) ! main array of variables
      REAL blkta(blksize) ! Air temperature [ K ]
      REAL blkprs(blksize)
! *** output:

! Air pressure in [ Pa ]
! concentration lower limit [ ug/m*
! lowest particle diameter ( m )
      REAL dgmin
      PARAMETER(dgmin=1.0E-09)

! lowest particle density ( Kg/m**3
      REAL densmin
      PARAMETER(densmin=1.0E03)

      REAL pmassn(blksize) ! mass concentration in nuclei mode
      REAL pmassa(blksize) ! mass concentration in accumulation
      REAL pmassc(blksize) ! mass concentration in coarse mode
      REAL pdensn(blksize) ! average particel density in Aitken
      REAL pdensa(blksize) ! average particel density in accumu
      REAL pdensc(blksize) ! average particel density in coarse
      REAL xlm(blksize) ! atmospheric mean free path [ m]
      REAL amu(blksize) ! atmospheric dynamic viscosity [ kg
      REAL dgnuc(blksize) ! Aitken mode mean diameter [ m ]
      REAL dgacc(blksize) ! accumulation
      REAL dgcor(blksize) ! coarse mode
      REAL knnuc(blksize) ! Aitken mode Knudsen number
      REAL knacc(blksize) ! accumulation
      REAL kncor(blksize)

! coarse mode

      INTEGER lcell
!      WRITE(20,*) ' IN MODPAR '

! *** set up  aerosol  3rd moment, mass, density

! loop counter
      DO lcell = 1, numcells

! *** Aitken-mode
!        cblk(lcell,vnu3) = max(conmin,(so4fac*cblk(lcell, & ! ghan
         cblk(lcell, vnu3) = so4fac*cblk(lcell, &
                                         vso4ai) + nh4fac*cblk(lcell, vnh4ai) + h2ofac*cblk(lcell, &
                                                                                            vh2oai) + no3fac*cblk(lcell, vno3ai) + &
                             nafac*cblk(lcell, vnaai) + clfac*cblk(lcell, vclai) + &
                             orgfac*cblk(lcell, &
                                         vasoa1i) + orgfac*cblk(lcell, vasoa2i) + orgfac*cblk(lcell, &
                                                                       vasoa3i) + orgfac*cblk(lcell, vasoa4i) + orgfac*cblk(lcell, &
                                                                       vbsoa1i) + orgfac*cblk(lcell, vbsoa2i) + orgfac*cblk(lcell, &
                                                                       vbsoa3i) + orgfac*cblk(lcell, vbsoa4i) + orgfac*cblk(lcell, &
                                                                  vorgpai) + anthfac*cblk(lcell, vp25ai) + anthfac*cblk(lcell, veci)
!          vorgpai)+anthfac*cblk(lcell,vp25ai)+anthfac*cblk(lcell,veci))) ! ghan

! *** Accumulation-mode
!        cblk(lcell,vac3) = max(conmin,(so4fac*cblk(lcell, & ! ghan
         cblk(lcell, vac3) = so4fac*cblk(lcell, &
                                         vso4aj) + nh4fac*cblk(lcell, vnh4aj) + h2ofac*cblk(lcell, &
                                                                                            vh2oaj) + no3fac*cblk(lcell, vno3aj) + &
                             nafac*cblk(lcell, vnaaj) + clfac*cblk(lcell, vclaj) + &
                             orgfac*cblk(lcell, &
                                         vasoa1j) + orgfac*cblk(lcell, vasoa2j) + orgfac*cblk(lcell, &
                                                                       vasoa3j) + orgfac*cblk(lcell, vasoa4j) + orgfac*cblk(lcell, &
                                                                       vbsoa1j) + orgfac*cblk(lcell, vbsoa2j) + orgfac*cblk(lcell, &
                                                                       vbsoa3j) + orgfac*cblk(lcell, vbsoa4j) + orgfac*cblk(lcell, &
                                                                  vorgpaj) + anthfac*cblk(lcell, vp25aj) + anthfac*cblk(lcell, vecj)
!          vorgpaj)+anthfac*cblk(lcell,vp25aj)+anthfac*cblk(lcell,vecj))) ! ghan

! *** coarse mode
!        cblk(lcell,vcor3) = max(conmin,(soilfac*cblk(lcell, & ! ghan rely on conmin applied to mass, not moment
!          vsoila)+seasfac*cblk(lcell,vseas)+anthfac*cblk(lcell,vantha)))

         cblk(lcell, vcor3) = soilfac*cblk(lcell, &
                                           vsoila) + seasfac*cblk(lcell, vseas) + anthfac*cblk(lcell, vantha)

! *** now get particle mass and density

! *** Aitken-mode:
         pmassn(lcell) = max(conmin, (cblk(lcell, vso4ai) + cblk(lcell, &
                                                                 vnh4ai) + cblk(lcell, vh2oai) + cblk(lcell, vno3ai) + cblk(lcell, &
                                                              vasoa1i) + cblk(lcell, vasoa2i) + cblk(lcell, vasoa3i) + cblk(lcell, &
                                                              vasoa4i) + cblk(lcell, vbsoa1i) + cblk(lcell, vbsoa2i) + cblk(lcell, &
                                                              vbsoa3i) + cblk(lcell, vbsoa4i) + cblk(lcell, vorgpai) + cblk(lcell, &
                                                                                                       vp25ai) + cblk(lcell, veci)))

! *** Accumulation-mode:
         pmassa(lcell) = max(conmin, (cblk(lcell, vso4aj) + cblk(lcell, &
                                                                 vnh4aj) + cblk(lcell, vh2oaj) + cblk(lcell, vno3aj) + cblk(lcell, &
                                                              vasoa1j) + cblk(lcell, vasoa2j) + cblk(lcell, vasoa3j) + cblk(lcell, &
                                                              vasoa4j) + cblk(lcell, vbsoa1j) + cblk(lcell, vbsoa2j) + cblk(lcell, &
                                                              vbsoa3j) + cblk(lcell, vbsoa4j) + cblk(lcell, vorgpaj) + cblk(lcell, &
                                                                                                       vp25aj) + cblk(lcell, vecj)))

! *** coarse mode:
         pmassc(lcell) = max(conmin, cblk(lcell, vsoila) + cblk(lcell, vseas) + cblk( &
                             lcell, vantha))

      END DO
! *** now get particle density, mean free path, and dynamic viscosity

! aerosol  3rd moment and  mass
      DO lcell = 1, &
                 numcells
! *** density in [ kg m**-3 ]

! Density and mean free path
         pdensn(lcell) = max(densmin, (f6dpim9*pmassn(lcell)/cblk(lcell, vnu3)))
         pdensa(lcell) = max(densmin, (f6dpim9*pmassa(lcell)/cblk(lcell, vac3)))
         pdensc(lcell) = max(densmin, (f6dpim9*pmassc(lcell)/cblk(lcell, vcor3)))

! *** Calculate mean free path [ m ]:
         xlm(lcell) = 6.6328E-8*pss0*blkta(lcell)/(tss0*blkprs(lcell))

! *** 6.6328E-8 is the sea level values given in Table I.2.8
! *** on page 10 of U.S. Standard Atmosphere 1962

! ***         Calculate dynamic viscosity [ kg m**-1 s**-1 ]:

! *** U.S. Standard Atmosphere 1962 page 14 expression
!     for dynamic viscosity is:
!     dynamic viscosity =  beta * T * sqrt(T) / ( T + S)
!     where beta = 1.458e-6 [ kg sec^-1 K**-0.5 ], s = 110.4 [ K ].

         amu(lcell) = 1.458E-6*blkta(lcell)*sqrt(blkta(lcell))/ &
                      (blkta(lcell) + 110.4)

      END DO
!...............   Standard deviation fixed in both modes, so
!...............   diagnose diameter from 3rd moment and number concentr

!  density and mean free path
      DO lcell = 1, &
                 numcells

! calculate diameters
         dgnuc(lcell) = max(dgmin, (cblk(lcell, vnu3)/(cblk(lcell, vnu0)*esn36))** &
                            one3)

         dgacc(lcell) = max(dgmin, (cblk(lcell, vac3)/(cblk(lcell, vac0)*esa36))** &
                            one3)

         dgcor(lcell) = max(dgmin, (cblk(lcell, vcor3)/(cblk(lcell, vcorn)*esc36)) &
                            **one3)

! when running with cloudborne aerosol, apply some very mild bounding
! to avoid unrealistic dg values
         if (cw_phase > 0) then
            dgnuc(lcell) = max(dgnuc(lcell), dginin*0.2) !  > 0.002 um
            dgnuc(lcell) = min(dgnuc(lcell), dginin*10.0) !  < 0.10  um
            dgacc(lcell) = max(dgacc(lcell), dginia*0.2) !  > 0.014 um
            dgacc(lcell) = min(dgacc(lcell), dginia*10.0) !  < 0.7 um
            dgcor(lcell) = max(dgcor(lcell), dginic*0.2) !  > 0.2 um
            dgcor(lcell) = min(dgcor(lcell), dginic*10.0) ! < 10.0 um
         end if

      END DO
! end loop on diameters
      DO lcell = 1, &
                 numcells
! Calculate Knudsen numbers
         knnuc(lcell) = 2.0*xlm(lcell)/dgnuc(lcell)

         knacc(lcell) = 2.0*xlm(lcell)/dgacc(lcell)

         kncor(lcell) = 2.0*xlm(lcell)/dgcor(lcell)

      END DO

! end loop for  Knudsen numbers
      RETURN

   END SUBROUTINE modpar
!------------------------------------------------------------------------------

   SUBROUTINE nuclcond(blksize, nspcsda, numcells, cblk, dt, layer, blkta, blkprs, &
                       blkrh, so4rat, organt1rat, organt2rat, organt3rat, organt4rat, orgbio1rat, &
                       orgbio2rat, orgbio3rat, orgbio4rat, drog, ldrog_vbs, ncv, nacv, dgnuc, dgacc, &
                       fconcn, fconca, fconcn_org, fconca_org, dmdt, dndt, deltaso4a, cgrn3, cgra3, igrid, jgrid, kgrid, brrto)

!***********************************************************************
!**    DESCRIPTION:  calculates aerosol nucleation and condensational
!**    growth rates using Binkowski and Shankar (1995) method.

! *** In this version, the method od RPM is followed where
!     the diffusivity, the average molecular ve3locity, and
!     the accomodation coefficient for sulfuric acid are used for
!     the organics. This is for consistency.
!       Future versions will use the correct values.  FSB 12/12/96

!**
!**    Revision history:
!       prototype 1/95 by Uma and Carlie
!       Corrected 7/95 by Uma for condensation of mass not nucleated
!       and mass conservation check
!       Revised   8/95 by US to calculate air density in stmt function
!                 and collect met variable stmt funcs in one include fil
!       Revised 7/25/96 by FSB to use block structure.
!       Revised 9/17/96 by FSB to use Y&K or K&W Nucleation mechanism
!       Revised 11/15/96 by FSB to use MKS,  and mom m^-3 units.
!       Revised 1/13/97 by FSB to pass arrays and simplify code.
!       Added   23/03/99 by BS growth factors for organics
!**********************************************************************
!     IMPLICIT NONE

!     Includes:
! *** arguments

! *** input;
!USE module_configure, only: grid_config_rec_type
!TYPE (grid_config_rec_type), INTENT (in) :: config_flags

! dimension of arrays
      INTEGER blksize
      INTEGER layer
! number of species in CBLK
      INTEGER nspcsda
! actual number of cells in arrays
      INTEGER numcells
      INTEGER igrid, jgrid, kgrid

      INTEGER ldrog_vbs
! # of organic aerosol precursor
      REAL cblk(blksize, nspcsda) ! main array of variables
! model time step in  SECONDS
      REAL dt
      REAL blkta(blksize) ! Air temperature [ K ]
      REAL blkprs(blksize) ! Air pressure in [ Pa ]
      REAL blkrh(blksize) ! Fractional relative humidity
      REAL so4rat(blksize) ! rate [  ug/m**3 /s ]
      REAL brrto
!bs
! sulfate gas-phase production
! total # of cond. vapors & SOA spe
      INTEGER ncv
!bs
      INTEGER nacv
!bs * anthropogenic organic condensable vapor production rate
! # of anthrop. cond. vapors & SOA
      REAL drog(blksize, ldrog_vbs) !bs
! Delta ROG conc. [ppm]

! anthropogenic vapor production rates
      REAL organt1rat(blksize)
      REAL organt2rat(blksize)
      REAL organt3rat(blksize)
      REAL organt4rat(blksize)

! biogenic vapor production rates
      REAL orgbio1rat(blksize)
      REAL orgbio2rat(blksize)
      REAL orgbio3rat(blksize)
      REAL orgbio4rat(blksize)

! biogenic organic aerosol production
      REAL dgnuc(blksize) ! accumulation
      REAL dgacc(blksize)
! *** output:

! coarse mode
      REAL fconcn(blksize) ! Aitken mode  [ 1 / s ]
! reciprocal condensation rate
      REAL fconca(blksize) ! acclumulation mode [ 1 / s ]
! reciprocal condensation rate
      REAL fconcn_org(blksize) ! Aitken mode  [ 1 / s ]
! reciprocal condensation rate
      REAL fconca_org(blksize) ! acclumulation mode [ 1 / s ]
! reciprocal condensation rate
      REAL dmdt(blksize) ! by particle formation [ ug/m**3 /s ]
! rate of production of new mass concent
      REAL dndt(blksize) ! concentration by particle formation [#
! rate of producton of new particle numb
      REAL deltaso4a(blksize) ! sulfate aerosol by condensation [ ug/m
! increment of concentration added to
      REAL cgrn3(blksize) ! Aitken mode [ 3rd mom/m **3 s ]
! growth rate for 3rd moment for
      REAL cgra3(blksize) ! Accumulation mode

!...........    SCRATCH local variables and their descriptions:

! growth rate for 3rd moment for

      INTEGER lcell

! LOOP INDEX
! conv rate so2 --> so4 [mom-3/g/s]
      REAL chemrat
! conv rate for organics [mom-3/g/s]
      REAL chemrat_org
      REAL am1n, & ! 1st mom density (nuc, acc modes) [mom_
         am1a
      REAL am2n, & ! 2nd mom density (nuc, acc modes) [mom_
         am2a
      REAL gnc3n, & ! near-cont fns (nuc, acc) for mom-3 den
         gnc3a
      REAL gfm3n, & ! free-mol  fns (nuc, acc) for mom-3 den
         gfm3a
! total reciprocal condensation rate
      REAL fconc

      REAL td
! d * tinf (cgs)
      REAL*8 & ! Cnstant to force 64 bit evaluation of
         one88
      PARAMETER(one88=1.0D0)
!  *** variables to set up sulfate and organic condensation rates

! sulfuric acid vapor at current time step
      REAL vapor1
!                                    chemistry and emissions
      REAL vapor2
! Sulfuric acid vapor prior to addition from
!bs
      REAL deltavap
!bs * start update
!bs
! change to vapor at previous time step
      REAL diffcorr

!bs *
      REAL csqt_org
!bs * end update
!bs
      REAL csqt
!.......................................................................
!   begin body of subroutine  NUCLCOND

!...........   Main computational grid-traversal loop nest
!...........   for computing condensation and nucleation:

      DO lcell = 1, &
                 numcells
! *** First moment:

!  1st loop over NUMCELLS
         am1n = cblk(lcell, vnu0)*dgnuc(lcell)*esn04
         am1a = cblk(lcell, vac0)*dgacc(lcell)*esa04

!..............   near-continuum factors [ 1 / sec ]
!bs
!bs * adopted from code of FSB
!bs * correction to DIFFSULF and DIFFORG for temperature and pressure
!bs
         diffcorr = (pss0/blkprs(lcell))*(blkta(lcell)/273.16)**1.
!bs
         gnc3n = cconc*am1n*diffcorr
         gnc3a = cconc*am1a*diffcorr

! *** Second moment:

         am2n = cblk(lcell, vnu0)*dgnuc(lcell)*dgnuc(lcell)*esn16
         am2a = cblk(lcell, vac0)*dgacc(lcell)*dgacc(lcell)*esa16

         csqt = ccofm*sqrt(blkta(lcell))
!...............   free molecular factors [ 1 / sec ]

! put in temperature fac
         gfm3n = csqt*am2n
         gfm3a = csqt*am2a

! *** Condensation factors in [ s**-1] for h2so4
! *** In the future, separate factors for condensing organics will
!      be included. In this version, the h2so4 values are used.

!...............   Twice the harmonic mean of fm, nc functions:
! *** Force 64 bit evaluation:

         fconcn(lcell) = one88*gnc3n*gfm3n/(gnc3n+gfm3n)
         fconca(lcell) = one88*gnc3a*gfm3a/(gnc3a+gfm3a)
         fconc = fconcn(lcell) + fconca(lcell)

! *** NOTE: FCONCN and FCONCA will be redefined below <<<<<<
!bs
!bs * start modifications for organcis
!bs
         gnc3n = cconc_org*am1n*diffcorr
         gnc3a = cconc_org*am1a*diffcorr
!bs
         csqt_org = ccofm_org*sqrt(blkta(lcell))
         gfm3n = csqt_org*am2n
         gfm3a = csqt_org*am2a
!bs
         fconcn_org(lcell) = one88*gnc3n*gfm3n/(gnc3n+gfm3n)
         fconca_org(lcell) = one88*gnc3a*gfm3a/(gnc3a+gfm3a)
!bs
!bs * end modifications for organics
!bs
! *** calculate the total change to sulfuric acid vapor from production
!                      and condensation

         vapor1 = cblk(lcell, vsulf) ! curent sulfuric acid vapor
         vapor2 = cblk(lcell, vsulf) - so4rat(lcell)* &
                  dt
! vapor at prev
         vapor2 = max(0.0, vapor2)
         deltavap = max(0.0, (so4rat(lcell)/fconc - vapor2)*(1.0-exp(-fconc*dt)))

! *** Calculate increment in total sufate aerosol mass concentration

! *** This follows the method of Youngblood & Kreidenweis.!bs
!bs        DELTASO4A( LCELL) = MAX( 0.0, SO4RAT(LCELL) * DT - DELTAVAP)
!bs
!bs * allow DELTASO4A to be negative, but the change must not be larger
!bs * than the amount of vapor available.
!bs
         deltaso4a(lcell) = max(-1.*cblk(lcell, vsulf), &
                                so4rat(lcell)*dt - deltavap)

! *** zero out growth coefficients
         cgrn3(lcell) = 0.0
         cgra3(lcell) = 0.0

      END DO

! *** Select method of nucleation
! End 1st loop over NUMCELLS
      IF (inucl == 1) THEN

! *** Do Youngblood & Kreidenweis  Nucleation

!         CALL BCSUINTF(DT,SO4RAT,FCONCN,FCONCA,BLKTA,BLKRH,
!     &        DNDT,DMDT,NUMCELLS,BLKSIZE,
!     &        VAPOR1)
!       IF (firstime) THEN
!         WRITE (6,*)
!         WRITE (6,'(a,i2)') 'INUCL =', inucl
!         WRITE (90,'(a,i2)') 'INUCL =', inucl
!         firstime = .FALSE.
!       END IF

      ELSE IF (inucl == 0) THEN

! *** Do Kerminen & Wexler Nucleation

!         CALL nuclKW(DT,SO4RAT,FCONCN,FCONCA,BLKTA,BLKRH,
!     &        DNDT,DMDT,NUMCELLS,BLKSIZE)
!       IF (firstime) THEN
!         WRITE (6,*)
!         WRITE (6,'(a,i2)') 'INUCL =', inucl
!         WRITE (90,'(a,i2)') 'INUCL =', inucl
!         firstime = .FALSE.
!       END IF

      ELSE IF (inucl == 2) THEN

!bs ** Do Kulmala et al. Nucleation
!       if(dndt(1).lt.-10.)print *,'before klpnuc',blkta(1),blkrh(1),vapor1,dndt(1),dmdt(1),so4rat(1)

         if (blkta(1) .ge. 233.15 .and. blkrh(1) .ge. 0.1) then
            CALL klpnuc(blkta(1), blkrh(1), vapor1, dndt(1), dmdt(1), so4rat(1))
         else
            dndt(1) = 0.
            dmdt(1) = 0.
         endif

!       CALL klpnuc(blkta(1),blkrh(1),vapor1,dndt(1),dmdt(1),so4rat(1))
!       if(dndt(1).lt.-10.)print *,'after klpnuc',dndt(1),dmdt(1)
         IF (dndt(1) == 0.) dmdt(1) = 0.
         IF (dmdt(1) == 0.) dndt(1) = 0.
!       IF (firstime) THEN
!         WRITE (6,*)
!         WRITE (6,'(a,i2)') 'INUCL =', inucl
!         WRITE (90,'(a,i2)') 'INUCL =', inucl
!         firstime = .FALSE.
!       END IF
!     ELSE
!       WRITE (6,'(a)') '*************************************'
!       WRITE (6,'(a,i2,a)') '  INUCL =', inucl, ',  PLEASE CHECK !!'
!       WRITE (6,'(a)') '        PROGRAM TERMINATED !!'
!       WRITE (6,'(a)') '*************************************'
!       STOP

      END IF
!bs
!bs * Secondary organic aerosol module (SOA_VBS)
!bs
! end of selection of nucleation method

      CALL soa_vbs(layer, blkta, blkprs, organt1rat, organt2rat, organt3rat, &
                   organt4rat, orgbio1rat, orgbio2rat, orgbio3rat, orgbio4rat, drog, ldrog_vbs, ncv, &
                   nacv, cblk, blksize, nspcsda, numcells, dt, igrid, jgrid, kgrid, brrto)
!bs
!bs *  Secondary organic aerosol module (SOA_VBS)
!bs
      DO lcell = 1, numcells

! *** redefine FCONCN & FCONCA to be the nondimensional fractionaL
!     condensation factors

         td = 1.0/(fconcn(lcell) + fconca(lcell))
         fconcn(lcell) = td*fconcn(lcell)
         fconca(lcell) = td*fconca(lcell)
!bs
         td = 1.0/(fconcn_org(lcell) + fconca_org(lcell))
         fconcn_org(lcell) = td*fconcn_org(lcell)
         fconca_org(lcell) = td*fconca_org(lcell)
!bs
      END DO

! ***  Begin second loop over cells

      DO lcell = 1, numcells
! *** note CHEMRAT includes  species other than sulfate.

! 3rd loop on NUMCELLS
         chemrat = so4fac*so4rat(lcell) ! [mom3 m**-3 s-
         chemrat_org = orgfac*(organt1rat(lcell) + organt2rat(lcell) + organt3rat( &
                               lcell) + organt4rat(lcell) + orgbio1rat(lcell) + orgbio2rat(lcell) + &
                               orgbio3rat(lcell) + orgbio4rat(lcell))

! *** Calculate the production rates for new particle
! [mom3 m**-3 s-
         cgrn3(lcell) = so4fac*dmdt(lcell)
! Rate of increase of 3rd
         chemrat = chemrat - cgrn3(lcell) !bs 3rd moment production fro

!bs Remove the rate of new pa
         chemrat = max(chemrat, 0.0)
! *** Now calculate the rate of condensation on existing particles.

! Prevent CHEMRAT from being negativ
         cgrn3(lcell) = cgrn3(lcell) + chemrat*fconcn(lcell) + &
                        chemrat_org*fconcn_org(lcell)
         cgra3(lcell) = chemrat*fconca(lcell) + chemrat_org*fconca_org(lcell)
! ***
      END DO
!  end 2nd loop over NUMCELLS
      RETURN

   END SUBROUTINE nuclcond
!------------------------------------------------------------------------------

! nuclcond
   REAL FUNCTION poly4(a, x)
      REAL a(4), x

      poly4 = a(1) + x*(a(2) + x*(a(3) + x*(a(4))))
      RETURN
   END FUNCTION poly4
   REAL FUNCTION poly6(a, x)
      REAL a(6), x

      poly6 = a(1) + x*(a(2) + x*(a(3) + x*(a(4) + x*(a(5) + x*(a(6))))))
      RETURN
   END FUNCTION poly6
!-----------------------------------------------------------------------

   SUBROUTINE rpmares_old(so4, hno3, no3, nh3, nh4, rh, temp, aso4, ano3, ah2o, anh4, &
                          gnh3, gno3)
! Description:

!   ARES calculates the chemical composition of a sulfate/nitrate/
!   ammonium/water aerosol based on equilibrium thermodynamics.

!   This code considers two regimes depending upon the molar ratio
!   of ammonium to sulfate.

!   For values of this ratio less than 2,the code solves a cubic for
!   hydrogen ion molality, HPLUS,  and if enough ammonium and liquid
!   water are present calculates the dissolved nitric acid. For molal
!   ionic strengths greater than 50, nitrate is assumed not to be present

!   For values of the molar ratio of 2 or greater, all sulfate is assumed
!   to be ammonium sulfate and a calculation is made for the presence of
!   ammonium nitrate.

!   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!   obtain the activity coefficients. Abandoned -7/30/97 FSB

!   The Bromley method of calculating the activity coefficients is used in this version

!   The calculation of liquid water is done in subroutine water. Details for both calculations are given
!   in the respective subroutines.

!   Based upon MARS due to
!   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld,
!   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.

!   and SCAPE due to
!   Kim, Seinfeld, and Saxeena, Aerosol Ceience and Technology,
!   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.

! NOTE: All concentrations supplied to this subroutine are TOTAL
!       over gas and aerosol phases

! Parameters:

!  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate (IN)
!  HNO3  : Nitric Acid in MICROGRAMS/M**3 as nitric acid (IN)
!  NO3   : Total nitrate in MICROGRAMS/M**3 as nitric acid (IN)
!  NH3   : Total ammonia in MICROGRAMS/M**3 as ammonia (IN)
!  NH4   : Ammonium in MICROGRAMS/M**3 as ammonium (IN)
!  RH    : Fractional relative humidity (IN)
!  TEMP  : Temperature in Kelvin (IN)
!  GNO3  : Gas phase nitric acid in MICROGRAMS/M**3 (OUT)
!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 (OUT)
!  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 (OUT)
!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 (OUT)
!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 (OUT)
!  AH2O  : Aerosol phase water in MICROGRAMS/M**3 (OUT)
!  NITR  : Number of iterations for obtaining activity coefficients  (OU
!  NR    : Number of real roots to the cubic in the low ammonia case (OU

! Revision History:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   11/10/87  Received the first version of the MARS code
!   S.Roselle   12/30/87  Restructured code
!   S.Roselle   2/12/88   Made correction to compute liquid-phase
!                         concentration of H2O2.
!   S.Roselle   5/26/88   Made correction as advised by SAI, for
!                         computing H+ concentration.
!   S.Roselle   3/1/89    Modified to operate with EM2
!   S.Roselle   5/19/89   Changed the maximum ionic strength from
!                         100 to 20, for numerical stability.
!   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!                         using equations for nitrate budget.
!   F.Binkowski 6/18/91   New ammonia poor case which
!                         omits letovicite.
!   F.Binkowski 7/25/91   Rearranged entire code, restructured
!                         ammonia poor case.
!   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!                         as SO4--
!   F.Binkowski 12/6/91   Changed the ammonia defficient case so that
!                         there is only neutralized sulfate (ammonium
!                         sulfate) and sulfuric acid.
!   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreemen
!                          with the Cohen et al. (1987)  maximum molalit
!                          of 36.2 in Table III.( J. Phys Chem (91) page
!                          4569, and Table IV p 4587.)
!   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to rem
!                         possibility for denomenator becoming zero;
!                         this involved solving for HPLUS first.
!                         Note that for a relative humidity
!                          less than 50%, the model assumes that there i
!                          aerosol nitrate.
!   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System
!                          Redid logic as follows
!                         1. Water algorithm now follows Spann & Richard
!                         2. Pitzer Multicomponent method used
!                         3. Multicomponent practical osmotic coefficien
!                            use to close iterations.
!                         4. The model now assumes that for a water
!                            mass fraction WFRAC less than 50% there is
!                            no aerosol nitrate.
!   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia p
!                         case, and changed the WFRAC criterion to 40%.
!                         For ammonium to sulfate ratio less than 1.0
!                         all ammonium is aerosol and no nitrate aerosol
!                         exists.
!   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case
!                         allow gas-phase ammonia to exist.
!   F.Binkowski 7/26/95   Changed equilibrium constants to values from
!                         Kim et al. (1993)
!   F.Binkowski 6/27/96   Changed to new water format
!   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent
!                         activity coefficients. The binary activity coe
!                         are the same as the previous version
!   F.Binkowski 8/1/97    Chenged minimum sulfate from 0.0 to 1.0e-6 i.e
!                         1 picogram per cubic meter

!-----------------------------------------------------------------------
!     IMPLICIT NONE
!...........INCLUDES and their descriptions
!cc      INCLUDE SUBST_CONST          ! constants
!...........PARAMETERS and their descriptions:

! molecular weight for NaCl
      REAL mwnacl
      PARAMETER(mwnacl=58.44277)

! molecular weight for NO3
      REAL mwno3
      PARAMETER(mwno3=62.0049)

! molecular weight for HNO3
      REAL mwhno3
      PARAMETER(mwhno3=63.01287)

! molecular weight for SO4
      REAL mwso4
      PARAMETER(mwso4=96.0576)

! molecular weight for HSO4
      REAL mwhso4
      PARAMETER(mwhso4=mwso4 + 1.0080)

! molecular weight for H2SO4
      REAL mh2so4
      PARAMETER(mh2so4=98.07354)

! molecular weight for NH3
      REAL mwnh3
      PARAMETER(mwnh3=17.03061)

! molecular weight for NH4
      REAL mwnh4
      PARAMETER(mwnh4=18.03858)

! molecular weight for Organic Species
      REAL mworg
      PARAMETER(mworg=16.0)

! molecular weight for Chloride
      REAL mwcl
      PARAMETER(mwcl=35.453)

! molecular weight for AIR
      REAL mwair
      PARAMETER(mwair=28.964)

! molecular weight for Letovicite
      REAL mwlct
      PARAMETER(mwlct=3.0*mwnh4 + 2.0*mwso4 + 1.0080)

! molecular weight for Ammonium Sulfa
      REAL mwas
      PARAMETER(mwas=2.0*mwnh4 + mwso4)

! molecular weight for Ammonium Bisul
      REAL mwabs
      PARAMETER(mwabs=mwnh4 + mwso4 + 1.0080)

!...........ARGUMENTS and their descriptions

!iamodels3
      REAL so4
! Total sulfate in micrograms / m**3
! Total nitric acid in micrograms / m
      REAL hno3
! Total nitrate in micrograms / m**3
      REAL no3
! Total ammonia in micrograms / m**3
      REAL nh3
! Total ammonium in micrograms / m**3
      REAL nh4
! Fractional relative humidity
      REAL rh
! Temperature in Kelvin
      REAL temp
! Aerosol sulfate in micrograms / m**
      REAL aso4
! Aerosol nitrate in micrograms / m**
      REAL ano3
! Aerosol liquid water content water
      REAL ah2o
! Aerosol ammonium in micrograms / m*
      REAL anh4
! Gas-phase nitric acid in micrograms
      REAL gno3
      REAL gnh3
!...........SCRATCH LOCAL VARIABLES and their descriptions:

! Gas-phase ammonia in micrograms / m
! Index set to percent relative humid
      INTEGER irh
! Number of iterations for activity c
      INTEGER nitr
! Loop index for iterations
      INTEGER nnn
      INTEGER nr
! Number of roots to cubic equation f
      REAL*8 & ! Coefficients and roots of
         a0
      REAL*8 & ! Coefficients and roots of
         a1
      REAL*8 & ! Coefficients and roots of
         a2
! Coefficients and discriminant for q
      REAL aa
! internal variables ( high ammonia c
      REAL bal
! Coefficients and discriminant for q
      REAL bb
! Variables used for ammonia solubili
      REAL bhat
! Coefficients and discriminant for q
      REAL cc
! Factor for conversion of units
      REAL convt
! Coefficients and discriminant for q
      REAL dd
! Coefficients and discriminant for q
      REAL disc
! Relative error used for convergence
      REAL eror
!  Free ammonia concentration , that
      REAL fnh3
! Activity Coefficient for (NH4+, HSO
      REAL gamaab
! Activity coefficient for (NH4+, NO3
      REAL gamaan
! Variables used for ammonia solubili
      REAL gamahat
! Activity coefficient for (H+ ,NO3-)
      REAL gamana
! Activity coefficient for (2H+, SO4-
      REAL gamas1
! Activity coefficient for (H+, HSO4-
      REAL gamas2
! used for convergence of iteration
      REAL gamold
! internal variables ( high ammonia c
      REAL gasqd
! Hydrogen ion (low ammonia case) (mo
      REAL hplus
! Equilibrium constant for ammoniua t
      REAL k1a
! Equilibrium constant for sulfate-bi
      REAL k2sa
! Dissociation constant for ammonium
      REAL k3
! Equilibrium constant for ammonium n
      REAL kan
! Variables used for ammonia solubili
      REAL khat
! Equilibrium constant for nitric aci
      REAL kna
! Henry's Law Constant for ammonia
      REAL kph
! Equilibrium constant for water diss
      REAL kw
! Internal variable using KAN
      REAL kw2
! Nitrate (high ammonia case) (moles
      REAL man
! Sulfate (high ammonia case) (moles
      REAL mas
! Bisulfate (low ammonia case) (moles
      REAL mhso4
! Nitrate (low ammonia case) (moles /
      REAL mna
! Ammonium (moles / kg water)
      REAL mnh4
! Total number of moles of all ions
      REAL molnu
! Sulfate (low ammonia case) (moles /
      REAL mso4
! Practical osmotic coefficient
      REAL phibar
! Previous value of practical osmotic
      REAL phiold
! Molar ratio of ammonium to sulfate
      REAL ratio
! Internal variable using K2SA
      REAL rk2sa
! Internal variables using KNA
      REAL rkna
! Internal variables using KNA
      REAL rknwet
      REAL rr1
      REAL rr2
! Ionic strength
      REAL stion
! Internal variables for temperature
      REAL t1
! Internal variables for temperature
      REAL t2
! Internal variables of convenience (
      REAL t21
! Internal variables of convenience (
      REAL t221
! Internal variables for temperature
      REAL t3
! Internal variables for temperature
      REAL t4
! Internal variables for temperature
      REAL t6
! Total ammonia and ammonium in micro
      REAL tnh4
! Total nitrate in micromoles / meter
      REAL tno3
! Tolerances for convergence test
      REAL toler1
! Tolerances for convergence test
      REAL toler2
! Total sulfate in micromoles / meter
      REAL tso4
! 2.0 * TSO4  (high ammonia case) (mo
      REAL twoso4
! Water mass fraction
      REAL wfrac
      ! micrograms / meter **3 on output
      REAL wh2o
      ! internally it is 10 ** (-6) kg (wat
      ! the conversion factor (1000 g = 1 k
      ! for AH2O output
! Aerosol liquid water content (inter
! internal variables ( high ammonia c
      REAL wsqd
! Nitrate aerosol concentration in mi
      REAL xno3
! Variable used in quadratic solution
      REAL xxq
! Ammonium aerosol concentration in m
      REAL ynh4
! Water variable saved in case ionic
      REAL zh2o

      REAL zso4
! Total sulfate molality - mso4 + mhs
      REAL cat(2) ! Array for cations (1, H+); (2, NH4+
      REAL an(3) ! Array for anions (1, SO4--); (2, NO
      REAL crutes(3) ! Coefficients and roots of
      REAL gams(2, 3) ! Array of activity coefficients
! Minimum value of sulfate laerosol c
      REAL minso4
      PARAMETER(minso4=1.0E-6/mwso4)
      REAL floor
      PARAMETER(floor=1.0E-30)
!-----------------------------------------------------------------------
!  begin body of subroutine RPMARES

!...convert into micromoles/m**3
!cc      WRITE( 10, * ) 'SO4, NO3, NH3 ', SO4, NO3, NH3
!iamodels3 merge NH3/NH4 , HNO3,NO3 here
! minimum concentration
      tso4 = max(0.0, so4/mwso4)
      tno3 = max(0.0, (no3/mwno3 + hno3/mwhno3))
      tnh4 = max(0.0, (nh3/mwnh3 + nh4/mwnh4))
!cc      WRITE( 10, * ) 'TSO4, TNO3, TNH4, RH ', TSO4, TNO3, TNH4, RH

!...now set humidity index IRH as a percent

      irh = nint(100.0*rh)

!...Check for valid IRH

      irh = max(1, irh)
      irh = min(99, irh)
!cc      WRITE(10,*)'RH,IRH ',RH,IRH

!...Specify the equilibrium constants at  correct
!...  temperature.  Also change units from ATM to MICROMOLE/M**3 (for KA
!...  KPH, and K3 )
!...  Values from Kim et al. (1993) except as noted.

      convt = 1.0/(0.082*temp)
      t6 = 0.082E-9*temp
      t1 = 298.0/temp
      t2 = alog(t1)
      t3 = t1 - 1.0
      t4 = 1.0+t2 - t1
      kna = 2.511E+06*exp(29.17*t3 + 16.83*t4)*t6
      k1a = 1.805E-05*exp(-1.50*t3 + 26.92*t4)
      k2sa = 1.015E-02*exp(8.85*t3 + 25.14*t4)
      kw = 1.010E-14*exp(-22.52*t3 + 26.92*t4)
      kph = 57.639*exp(13.79*t3 - 5.39*t4)*t6
!cc      K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6
      khat = kph*k1a/kw
      kan = kna*khat

!...Compute temperature dependent equilibrium constant for NH4NO3
!...  ( from Mozurkewich, 1993)
      k3 = exp(118.87-24084.0/temp - 6.025*alog(temp))

!...Convert to (micromoles/m**3) **2
      k3 = k3*convt*convt
      wh2o = 0.0
      stion = 0.0
      ah2o = 0.0
      mas = 0.0
      man = 0.0
      hplus = 0.0
      toler1 = 0.00001
      toler2 = 0.001
      nitr = 0
      nr = 0
      ratio = 0.0
      gamaan = 1.0
      gamold = 1.0

!...set the ratio according to the amount of sulfate and nitrate
      IF (tso4 > minso4) THEN
         ratio = tnh4/tso4

!...If there is no sulfate and no nitrate, there can be no ammonium
!...  under the current paradigm. Organics are ignored in this version.

      ELSE

         IF (tno3 == 0.0) THEN

! *** If there is very little sulfate and no nitrate set concentrations
!      to a very small value and return.
            aso4 = max(floor, aso4)
            ano3 = max(floor, ano3)
            wh2o = 0.0
            ah2o = 0.0
            gnh3 = max(floor, gnh3)
            gno3 = max(floor, gno3)
            RETURN
         END IF

!...For the case of no sulfate and nonzero nitrate, set ratio to 5
!...  to send the code to the high ammonia case

         ratio = 5.0
      END IF

!....................................
!......... High Ammonia Case ........
!....................................

      IF (ratio > 2.0) THEN

         gamaan = 0.1

!...Set up twice the sulfate for future use.

         twoso4 = 2.0*tso4
         xno3 = 0.0
         ynh4 = twoso4

!...Treat different regimes of relative humidity

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  start with ammomium sulfate solution without nitrate

         CALL awater(irh, tso4, ynh4, tno3, ah2o) !**** note TNO3
         wh2o = 1.0E-3*ah2o
         aso4 = tso4*mwso4
         ano3 = 0.0
         anh4 = ynh4*mwnh4
         wfrac = ah2o/(aso4 + anh4 + ah2o)
!cc        IF ( WFRAC .EQ. 0.0 )  RETURN   ! No water
         IF (wfrac < 0.2) THEN

!... dry  ammonium sulfate and ammonium nitrate
!...  compute free ammonia

            fnh3 = tnh4 - twoso4
            cc = tno3*fnh3 - k3

!...check for not enough to support aerosol

            IF (cc <= 0.0) THEN
               xno3 = 0.0
            ELSE
               aa = 1.0
               bb = -(tno3 + fnh3)
               disc = bb*bb - 4.0*cc

!...Check for complex roots of the quadratic
!...  set nitrate to zero and RETURN if complex roots are found

               IF (disc < 0.0) THEN
                  xno3 = 0.0
                  ah2o = 1000.0*wh2o
                  ynh4 = twoso4
                  gno3 = tno3*mwhno3
                  gnh3 = (tnh4 - ynh4)*mwnh3
                  aso4 = tso4*mwso4
                  ano3 = 0.0
                  anh4 = ynh4*mwnh4
                  RETURN
               END IF

!...to get here, BB .lt. 0.0, CC .gt. 0.0 always

               dd = sqrt(disc)
               xxq = -0.5*(bb + sign(1.0, bb)*dd)

!...Since both roots are positive, select smaller root.

               xno3 = min(xxq/aa, cc/xxq)

            END IF
            ah2o = 1000.0*wh2o
            ynh4 = 2.0*tso4 + xno3
            gno3 = (tno3 - xno3)*mwhno3
            gnh3 = (tnh4 - ynh4)*mwnh3
            aso4 = tso4*mwso4
            ano3 = xno3*mwno3
            anh4 = ynh4*mwnh4
            RETURN

         END IF

!...liquid phase containing completely neutralized sulfate and
!...  some nitrate.  Solve for composition and quantity.

         mas = tso4/wh2o
         man = 0.0
         xno3 = 0.0
         ynh4 = twoso4
         phiold = 1.0

!...Start loop for iteration

!...The assumption here is that all sulfate is ammonium sulfate,
!...  and is supersaturated at lower relative humidities.

         DO nnn = 1, 150
            nitr = nnn
            gasqd = gamaan*gamaan
            wsqd = wh2o*wh2o
            kw2 = kan*wsqd/gasqd
            aa = 1.0-kw2
            bb = twoso4 + kw2*(tno3 + tnh4 - twoso4)
            cc = -kw2*tno3*(tnh4 - twoso4)

!...This is a quadratic for XNO3 [MICROMOLES / M**3] of nitrate in solut

            disc = bb*bb - 4.0*aa*cc

!...Check for complex roots, if so set nitrate to zero and RETURN

            IF (disc < 0.0) THEN
               xno3 = 0.0
               ah2o = 1000.0*wh2o
               ynh4 = twoso4
               gno3 = tno3*mwhno3
               gnh3 = (tnh4 - ynh4)*mwnh3
               aso4 = tso4*mwso4
               ano3 = 0.0
               anh4 = ynh4*mwnh4
!cc            WRITE( 10, * ) ' COMPLEX ROOTS '
               RETURN
            END IF

            dd = sqrt(disc)
            xxq = -0.5*(bb + sign(1.0, bb)*dd)
            rr1 = xxq/aa
            rr2 = cc/xxq

!...Check for two non-positive roots, if so set nitrate to zero and RETURN
            IF (rr1 <= 0.0 .AND. rr2 <= 0.0) THEN
               xno3 = 0.0
               ah2o = 1000.0*wh2o
               ynh4 = twoso4
               gno3 = tno3*mwhno3
               gnh3 = (tnh4 - ynh4)*mwnh3
               aso4 = tso4*mwso4
               ano3 = 0.0
               anh4 = ynh4*mwnh4
!            WRITE(*,*) 'TWO NON-POSITIVE ROOTS!!! '
               RETURN
            END IF

!...choose minimum positve root

            IF ((rr1*rr2) < 0.0) THEN
               xno3 = max(rr1, rr2)
            ELSE
               xno3 = min(rr1, rr2)
            END IF
            xno3 = min(xno3, tno3)

!...This version assumes no solid sulfate forms (supersaturated )
!...  Now update water

            CALL awater(irh, tso4, ynh4, xno3, ah2o)

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  The conversion from micromoles to moles is done by the units of WH

            wh2o = 1.0E-3*ah2o

!...Ionic balance determines the ammonium in solution.

            man = xno3/wh2o
            mas = tso4/wh2o
            mnh4 = 2.0*mas + man
            ynh4 = mnh4*wh2o

!...MAS, MAN and MNH4 are the aqueous concentrations of sulfate, nitrate
!...  and ammonium in molal units (moles/(kg water) ).

            stion = 3.0*mas + man
            cat(1) = 0.0
            cat(2) = mnh4
            an(1) = mas
            an(2) = man
            an(3) = 0.0
            CALL actcof(cat, an, gams, molnu, phibar)
            gamaan = gams(2, 2)

!...Use GAMAAN for convergence control

            eror = abs(gamold - gamaan)/gamold
            gamold = gamaan

!...Check to see if we have a solution

            IF (eror <= toler1) THEN
!cc            WRITE( 11, * ) RH, STION, GAMS( 1, 1 ),GAMS( 1, 2 ), GAMS
!cc     &      GAMS( 2, 1 ), GAMS( 2, 2 ), GAMS( 2, 3 ), PHIBAR

               aso4 = tso4*mwso4
               ano3 = xno3*mwno3
               anh4 = ynh4*mwnh4
               gno3 = (tno3 - xno3)*mwhno3
               gnh3 = (tnh4 - ynh4)*mwnh3
               ah2o = 1000.0*wh2o
               RETURN
            END IF

         END DO

!...If after NITR iterations no solution is found, then:

         aso4 = tso4*mwso4
         ano3 = 0.0
         ynh4 = twoso4
         anh4 = ynh4*mwnh4
         CALL awater(irh, tso4, ynh4, xno3, ah2o)
         gno3 = tno3*mwhno3
         gnh3 = (tnh4 - ynh4)*mwnh3
         RETURN

      ELSE
!......................................
!......... Low Ammonia Case ...........
!......................................

!...coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)

!...All cases covered by this logic
         wh2o = 0.0
         CALL awater(irh, tso4, tnh4, tno3, ah2o)
         wh2o = 1.0E-3*ah2o
         zh2o = ah2o
!...convert 10**(-6) kg water/(cubic meter of air) to micrograms of wate
!...  per cubic meter of air (1000 g = 1 kg)

         aso4 = tso4*mwso4
         anh4 = tnh4*mwnh4
         ano3 = 0.0
         gno3 = tno3*mwhno3
         gnh3 = 0.0

!...Check for zero water.
         IF (wh2o == 0.0) RETURN
         zso4 = tso4/wh2o

!...ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4

!cc         IF ( ZSO4 .GT. 11.0 ) THEN

!...do not solve for aerosol nitrate for total sulfate molality
!...  greater than 11.0 because the model parameters break down
!...  greater than  9.0 because the model parameters break down

         IF (zso4 > 9.0) & ! 18 June 97
            THEN
            RETURN
         END IF

!...First solve with activity coeffs of 1.0, then iterate.
         phiold = 1.0
         gamana = 1.0
         gamas1 = 1.0
         gamas2 = 1.0
         gamaab = 1.0
         gamold = 1.0

!...All ammonia is considered to be aerosol ammonium.
         mnh4 = tnh4/wh2o

!...MNH4 is the molality of ammonium ion.
         ynh4 = tnh4

!...loop for iteration
         DO nnn = 1, 150
            nitr = nnn

!...set up equilibrium constants including activities
!...  solve the system for hplus first then sulfate & nitrate
!          print*,'gamas,gamana',gamas1,gamas2,gamana
            rk2sa = k2sa*gamas2*gamas2/(gamas1*gamas1*gamas1)
            rkna = kna/(gamana*gamana)
            rknwet = rkna*wh2o
            t21 = zso4 - mnh4
            t221 = zso4 + t21

!...set up coefficients for cubic

            a2 = rk2sa + rknwet - t21
            a1 = rk2sa*rknwet - t21*(rk2sa + rknwet) - rk2sa*zso4 - rkna*tno3
            a0 = -(t21*rk2sa*rknwet + rk2sa*rknwet*zso4 + rk2sa*rkna*tno3)

            CALL cubic(a2, a1, a0, nr, crutes)

!...Code assumes the smallest positive root is in CRUTES(1)

            hplus = crutes(1)
            bal = hplus**3 + a2*hplus**2 + a1*hplus + a0
            mso4 = rk2sa*zso4/(hplus + rk2sa) ! molality of sulfat
            mhso4 = zso4 - & ! molality of bisulf
                    mso4
            mna = rkna*tno3/(hplus + rknwet) ! molality of nitrat
            mna = max(0.0, mna)
            mna = min(mna, tno3/wh2o)
            xno3 = mna*wh2o
            ano3 = mna*wh2o*mwno3
            gno3 = (tno3 - xno3)*mwhno3

!...Calculate ionic strength
            stion = 0.5*(hplus + mna + mnh4 + mhso4 + 4.0*mso4)

!...Update water
            CALL awater(irh, tso4, ynh4, xno3, ah2o)

!...Convert 10**(-6) kg water/(cubic meter of air) to micrograms of wate
!...  per cubic meter of air (1000 g = 1 kg)

            wh2o = 1.0E-3*ah2o
            cat(1) = hplus
            cat(2) = mnh4
            an(1) = mso4
            an(2) = mna
            an(3) = mhso4
!          print*,'actcof',cat(1),cat(2),an(1),an(2),an(3),gams,molnu,phibar
            CALL actcof(cat, an, gams, molnu, phibar)

            gamana = gams(1, 2)
            gamas1 = gams(1, 1)
            gamas2 = gams(1, 3)
            gamaan = gams(2, 2)

            gamahat = (gamas2*gamas2/(gamaab*gamaab))
            bhat = khat*gamahat
!cc          EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
!cc          PHIOLD = PHIBAR
            eror = abs(gamold - gamahat)/gamold
            gamold = gamahat

!...write out molalities and activity coefficient
!...  and return with good solution

            IF (eror <= toler2) THEN
!cc            WRITE(12,*) RH, STION,HPLUS,ZSO4,MSO4,MHSO4,MNH4,MNA
!cc            WRITE(11,*) RH, STION, GAMS(1,1),GAMS(1,2),GAMS(1,3),
!cc     &                  GAMS(2,1),GAMS(2,2),GAMS(2,3), PHIBAR
               RETURN
            END IF

         END DO

!...after NITR iterations, failure to solve the system, no ANO3

         gno3 = tno3*mwhno3
         ano3 = 0.0
         CALL awater(irh, tso4, tnh4, tno3, ah2o)
         RETURN

      END IF
! ratio .gt. 2.0
   END SUBROUTINE rpmares_old

!ia*********************************************************
!ia                                                        *
!ia BEGIN OF AEROSOL ROUTINE                                   *
!ia                                                           *
!ia*********************************************************

!***********************************************************************
!           BEGIN OF AEROSOL CALCULATIONS
!***********************************************************************
!ia                                                                          *
!ia        MAIN AEROSOL DYNAMICS ROUTINE                                        *
!ia        based on MODELS3 formulation by FZB                                *
!ia        Modified by IA in May 97                                        *
!ia     THIS PROGRAMME IS THE LINK BETWEEN GAS PHASE AND AEROSOL PHASE
!ia     CALCULATIONS IN THE COLUMN MODEL. IT CONVERTS ALL DATA AND
!ia     VARIABLES BETWEEN THE TWO PARTS AND DRIVES THE AEROSOL
!ia     CALCULATIONS.
!ia     INPUT DATA REQUIRED FOR AEROSOL DYNAMICS ARE SET UP HERE FOR
!ia     ONE GRID CELL!!!!
!ia     and passed to dynamics calcs. subroutines.
!ia                                                                        *
!ia        Revision history                                                *
!ia        When        WHO        WHAT                                                *
!ia        ----        ----        ----                                                *
!ia        ????        FZB        BEGIN                                                *
!ia        05/97        IA        Adapted for use in CTM2-S                        *
!ia                        Modified renaming/bug fixing                        *
!ia     11/97   IA      Modified for new model version
!ia                     see comments under iarev02
!ia     03/98   IA      corrected error on pressure units
!ia                                                                        *
!ia        Called BY:        CHEM                                                *
!ia                                                                        *
!ia        Calls to:        OUTPUT1,AEROPRC                                        *
!ia                                                                        *
!ia*********************************************************************

! end RPMares
! convapr_in is removed, it wasn't used indeed
   SUBROUTINE rpmmod3(nspcsda, blksize, layer, dtsec, pres, temp, relhum, &
                      nitrate_in, nh3_in, vsulf_in, so4rat_in, drog_in, ldrog_vbs, ncv, &
                      nacv, eeci_in, eecj_in, eorgi_in, eorgj_in, epm25i, epm25j, epmcoarse, &
                      soilrat_in, cblk, igrid, jgrid, kgrid, brrto)

!USE module_configure, only: grid_config_rec_type
!TYPE (grid_config_rec_type), INTENT (in) :: config_flags

!     IMPLICIT NONE

!     Includes:
!iarev02       INCLUDE  AEROINCL.EXT
! block size, set to 1 in column model  ciarev0
      INTEGER blksize
!ia                             kept to 1 in current version of column model
! actual number of cells in arrays ( default is
      INTEGER, PARAMETER  :: numcells = 1

      INTEGER layer
! number of layer (default is 1 in

! index for cell in blocked array (default is 1 in
      INTEGER, PARAMETER :: ncell = 1
! *** inputs
! Input temperature [ K ]
      REAL temp
! Input relative humidity  [ fraction ]
      REAL relhum
! Input pressure [ hPa ]
      REAL pres
! Input number for Aitken mode [ m**-3 ]
      REAL numnuc_in
! Input number for accumulation mode [ m**-3 ]
      REAL numacc_in
! Input number for coarse mode  [ m**-3 ]
      REAL numcor_in
      ! sulfuric acid [ ug m**-3 ]
      REAL vsulf_in
! total sulfate vapor as sulfuric acid as
      ! sulfuric acid [ ug m**-3 ]
      REAL asulf_in
! total sulfate aerosol as sulfuric acid as
! i-mode sulfate input as sulfuric acid [ ug m*
      REAL asulfi_in
! ammonia gas [  ug m**-3 ]
      REAL nh3_in
! input value of nitric acid vapor [ ug m**-3 ]
      REAL nitrate_in
! Production rate of sulfuric acid   [ ug m**-3
      REAL so4rat_in
      ! aerosol [ ug m**-3 s**-1 ]
      REAL soilrat_in
! Production rate of soil derived coarse
! Emission rate of i-mode EC [ug m**-3 s**-1]
      REAL eeci_in
! Emission rate of j-mode EC [ug m**-3 s**-1]
      REAL eecj_in
! Emission rate of j-mode org. aerosol [ug m**-
      REAL eorgi_in
      REAL eorgj_in
! Emission rate of j-mode org. aerosol [ug m**-
! total # of cond. vapors & SOA species
      INTEGER ncv
! # of anthrop. cond. vapors & SOA speci
      INTEGER nacv
! # of organic aerosol precursor
      INTEGER ldrog_vbs
      REAL drog_in(ldrog_vbs) ! organic aerosol precursor [ppm]
! Input delta ROG concentration of
      REAL condvap_in(ncv) ! cond. vapor input [ug m^-3]
      REAL drog(blksize, ldrog_vbs) ! organic aerosol precursor [ppm]

      REAL brrto
!bs
! *** Primary emissions rates: [ ug / m**3 s ]

! *** emissions rates for unidentified PM2.5 mass
! Delta ROG concentration of
      REAL epm25i(blksize) ! Aitken mode
      REAL epm25j(blksize)
! *** emissions rates for primary organic aerosol
! Accumululaton mode
      REAL eorgi(blksize) ! Aitken mode
      REAL eorgj(blksize)
! *** emissions rates for elemental carbon
! Accumululaton mode
      REAL eeci(blksize) ! Aitken mode
      REAL eecj(blksize)
! *** Primary emissions rates [ ug m**-3 s -1 ] :

! Accumululaton mode
      REAL epm25(blksize) ! emissions rate for PM2.5 mass
      REAL esoil(blksize) ! emissions rate for soil derived coarse a
      REAL eseas(blksize) ! emissions rate for marine coarse aerosol
      REAL epmcoarse(blksize)
! emissions rate for anthropogenic coarse

      REAL dtsec
! time step [ s ], PASSED FROM MAIN COLUMN MODE

      REAL newm3
      REAL totaersulf
! total aerosol sulfate
! loop index for time steps
      INTEGER numsteps
      REAL step

! *** arrays for aerosol model codes:

! synchronization time  [s]

      INTEGER nspcsda

! number of species in CBLK ciarev02
      REAL cblk(blksize, nspcsda)

! *** Meteorological information in blocked arays:

! *** Thermodynamic variables:

! main array of variables
      REAL blkta(blksize) ! Air temperature [ K ]
      REAL blkprs(blksize) ! Air pressure in [ Pa ]
      REAL blkdens(blksize) ! Air density  [ kg m^-3 ]
      REAL blkrh(blksize)

! *** Chemical production rates [ ug m**-3 s -1 ] :

! Fractional relative humidity
      REAL so4rat(blksize) ! rate [ug/m^3/s]
! sulfuric acid vapor-phase production
      REAL organt1rat(blksize) ! production rate from aromatics [ ug /
! anthropogenic organic aerosol mass
      REAL organt2rat(blksize) ! production rate from aromatics [ ug /
! anthropogenic organic aerosol mass
      REAL organt3rat(blksize) ! rate from alkanes & others [ ug / m^3
! anthropogenic organic aerosol mass pro
      REAL organt4rat(blksize) ! rate from alkanes & others [ ug / m^3
! anthropogenic organic aerosol mass pro
      REAL orgbio1rat(blksize) ! rate [ ug / m^3 s ]
! biogenic organic aerosol production
      REAL orgbio2rat(blksize) ! rate [ ug / m^3 s ]
! biogenic organic aerosol production
      REAL orgbio3rat(blksize) ! rate [ ug / m^3 s ]
! biogenic organic aerosol production
      REAL orgbio4rat(blksize) ! rate [ ug / m^3 s ]
!bs
! *** atmospheric properties

! biogenic organic aerosol production
      REAL xlm(blksize) ! atmospheric mean free path [ m ]
      REAL amu(blksize)
! *** aerosol properties:

! *** modal diameters:

! atmospheric dynamic viscosity [ kg
      REAL dgnuc(blksize) ! nuclei mode geometric mean diamete
      REAL dgacc(blksize) ! accumulation geometric mean diamet
      REAL dgcor(blksize)

! *** Modal mass concentrations [ ug m**3 ]

! coarse mode geometric mean diamete
      REAL pmassn(blksize) ! mass concentration in Aitken mode
      REAL pmassa(blksize) ! mass concentration in accumulation
      REAL pmassc(blksize)
! *** average modal particle densities  [ kg/m**3 ]

! mass concentration in coarse mode
      REAL pdensn(blksize) ! average particle density in nuclei
      REAL pdensa(blksize) ! average particle density in accumu
      REAL pdensc(blksize)
! *** average modal Knudsen numbers

! average particle density in coarse
      REAL knnuc(blksize) ! nuclei mode  Knudsen number
      REAL knacc(blksize) ! accumulation Knudsen number
      REAL kncor(blksize)
! *** reciprocal modal condensation rates for sulfuric acid [ 1/s ]

! coarse mode  Knudsen number
      REAL fconcn(blksize)
! reciprocal condensation rate Aitke
      REAL fconca(blksize) !bs
! reciprocal condensation rate acclu
      REAL fconcn_org(blksize)
      REAL fconca_org(blksize)

! *** Rates for secondary particle formation:

! *** production of new mass concentration [ ug/m**3 s ]
      REAL dmdt(blksize) !                                 by particle formation

! *** production of new number concentration [ number/m**3 s ]

! rate of production of new mass concen
      REAL dndt(blksize) !                                 by particle formation
! *** growth rate for third moment by condensation of precursor
!      vapor on existing particles [ 3rd mom/m**3 s ]

! rate of producton of new particle num
      REAL cgrn3(blksize) !  Aitken mode
      REAL cgra3(blksize)
! *** Rates for coaglulation: [ m**3/s ]

! *** Unimodal Rates:

!  Accumulation mode
      REAL urn00(blksize) ! Aitken mode 0th moment self-coagulation ra
      REAL ura00(blksize)

! *** Bimodal Rates:  Aitken mode with accumulation mode ( Aitken mode)
! accumulation mode 0th moment self-coagulat
      REAL brna01(blksize) ! rate for 0th moment
      REAL brna31(blksize)
! *** other processes

! rate for 3rd moment
      REAL deltaso4a(blksize) !                                  sulfate aerosol by condensation   [ u

! *** housekeeping variables:
! increment of concentration added to
      INTEGER unit
      PARAMETER(unit=30)
      CHARACTER*16 pname
      PARAMETER(pname=' BOX            ')
      INTEGER isp, igrid, jgrid, kgrid

! loop index for species.
      INTEGER ii, iimap(8)
      DATA iimap/1, 2, 18, 19, 21, 22, 23, 24/

!   begin body  of program box

! *** Set up files and other info
! *** set up experimental conditions
! *** initialize model variables
!ia *** not required any more

!ia       DO ISP = 1, NSPCSDA
!ia       CBLK(BLKSIZE,ISP) = 1.0e-19 ! set CBLK to a very small number
!ia       END DO

      step = dtsec ! set time step

      blkta(blksize) = temp ! T in Kelvin

      blkprs(blksize) = pres*100. ! P in  Pa (pres is given in

      blkrh(blksize) = relhum ! fractional RH

      blkdens(blksize) = blkprs(blksize)/(rdgas*blkta(blksize)) !rs      CBLK(BLKSIZE,VSULF) = vsulf_in

!rs      CBLK(BLKSIZE,VHNO3) = nitrate_in
!rs      CBLK(BLKSIZE,VNH3) =  nh3_in
!bs
!rs      CBLK(BLKSIZE,VCVARO1) = condvap_in(PSOAARO1)
!rs      CBLK(BLKSIZE,VCVARO2) = condvap_in(PSOAARO2)
!rs      CBLK(BLKSIZE,VCVALK1) = condvap_in(PSOAALK1)
!rs      CBLK(BLKSIZE,VCVOLE1) = condvap_in(PSOAOLE1)
!rs      CBLK(BLKSIZE,VCVAPI1) = condvap_in(PSOAAPI1)
!rs      CBLK(BLKSIZE,VCVAPI2) = condvap_in(PSOAAPI2)
!rs      CBLK(BLKSIZE,VCVLIM1) = condvap_in(PSOALIM1)
!rs      CBLK(BLKSIZE,VCVLIM2) = condvap_in(PSOALIM2)

      DO isp = 1, ldrog_vbs
         drog(blksize, isp) = drog_in(isp)
      END DO

!      print*,'drog in rpm',drog
!bs
!ia *** 27/05/97 the following variables are transported quantities
!ia *** of the column-model now and thuse do not need this init.
!ia *** step.

!     CBLK(BLKSIZE,VNU0) = numnuc_in
!     CBLK(BLKSIZE,VAC0) = numacc_in
!     CBLK(BLKSIZE,VSO4A) =  asulf_in
!     CBLK(BLKSIZE,VSO4AI) = asulfi_in
!     CBLK(BLKSIZE, VCORN) = numcor_in

      so4rat(blksize) = so4rat_in

!...INITIALISE EMISSION RATES

!     epm25i(blksize) = & ! unidentified PM2.5 mass
!       0.0
!     epm25j(blksize) = &
!       0.0
! unidentified PM2.5 m
      eorgi(blksize) = & ! primary organic
         eorgi_in
      eorgj(blksize) = &
         eorgj_in
! primary organic
      eeci(blksize) = & ! elemental carbon
         eeci_in
      eecj(blksize) = &
         eecj_in
! elemental carbon
      epm25(blksize) = & !currently from input file ACTIONIA
         0.0
      esoil(blksize) = & ! ACTIONIA
         soilrat_in
      eseas(blksize) = & !currently from input file ACTIONIA
         0.0
!     epmcoarse(blksize) = & !currently from input file ACTIONIA
!       0.0
      dgnuc(blksize) = dginin
      dgacc(blksize) = dginia
      dgcor(blksize) = dginic
      newm3 = 0.0

! *** Set up initial total 3rd moment factors

      totaersulf = 0.0
      newm3 = 0.0
! ***  time loop
! write(50,*) ' numsteps dgnuc dgacc ', ' aso4 aso4i Ni Nj ah2o ah2oi M3i m3j'

! *** Call aerosol routines
      CALL aeroproc(blksize, nspcsda, numcells, layer, cblk, step, blkta, blkprs, &
                    blkdens, blkrh, so4rat, organt1rat, organt2rat, organt3rat, &
                    organt4rat, orgbio1rat, orgbio2rat, orgbio3rat, orgbio4rat, drog, ldrog_vbs, ncv, &
                    nacv, epm25i, epm25j, eorgi, eorgj, eeci, eecj, epmcoarse, esoil, eseas, xlm, &
                    amu, dgnuc, dgacc, dgcor, pmassn, pmassa, pmassc, pdensn, pdensa, pdensc, knnuc, &
                    knacc, kncor, fconcn, fconca, fconcn_org, fconca_org, dmdt, dndt, cgrn3, cgra3, &
                    urn00, ura00, brna01, brna31, deltaso4a, igrid, jgrid, kgrid, brrto)

! *** write output
!      WRITE(UNIT,*) ' AFTER AEROPROC '
!      WRITE(UNIT,*) ' NUMSTEPS = ', NUMSTEPS

! *** Write out file for graphing.

!     write(50,*) NUMSTEPS, DGNUC,DGACC,(CBLK(1,iimap(ii)),ii=1,8)

! *** update sulfuric acid vapor
!ia 21.04.98 this update is not required here
!ia artefact from box model
!       CBLK(BLKSIZE,VSULF) = CBLK(BLKSIZE,VSULF) +
!    &            SO4RAT(BLKSIZE) * STEP

      RETURN
   END SUBROUTINE rpmmod3
!---------------------------------------------------------------------------
   SUBROUTINE soa_vbs(layer, blkta, blkprs, organt1rat, organt2rat, organt3rat, &
                      organt4rat, orgbio1rat, orgbio2rat, orgbio3rat, orgbio4rat, drog, ldrog_vbs, ncv, &
                      nacv, cblk, blksize, nspcsda, numcells, dt, igrid, jgrid, kgrid, brrto)

      USE module_configure, only:grid_config_rec_type

      ! model layer
      INTEGER layer
      ! dimension of arrays
      INTEGER blksize
      ! number of species in CBLK
      INTEGER nspcsda ! actual number of cells in arrays
      INTEGER numcells ! # of organic aerosol precursor
      INTEGER ldrog_vbs ! total # of cond. vapors & SOA sp
      INTEGER ncv ! # of anthrop. cond. vapors & SOA
      INTEGER nacv
      INTEGER igrid, jgrid, kgrid

      REAL cblk(blksize, nspcsda) ! main array of variables
      REAL dt ! model time step in  SECONDS
      REAL blkta(blksize) ! Air temperature [ K ]
      REAL blkprs(blksize) ! Air pressure in [ Pa ]

      REAL, INTENT(OUT) :: brrto ! branching ratio for NOx conditions

      ! anthropogenic organic vapor production rates

      REAL organt1rat(blksize) ! rates from
      REAL organt2rat(blksize) ! rates from
      REAL organt3rat(blksize) ! rates from
      REAL organt4rat(blksize) ! rates from

      ! biogenic organic vapor production rates
      REAL orgbio1rat(blksize)
      REAL orgbio2rat(blksize)
      REAL orgbio3rat(blksize)
      REAL orgbio4rat(blksize)
      REAL drog(blksize, ldrog_vbs) !blksize=1, ldrog_vbs=9+1, the last ldrog_vbs is actually is the branching ratio

      !bs * local variable declaration
      ! Delta ROG conc. [ppm]
      !bs numerical value for a minimum thresh
      REAL, PARAMETER :: thrsmin = 1.E-19
      !bs numerical value for a minimum thresh
      !bs
      !bs universal gas constant [J/mol-K]
      REAL, PARAMETER :: rgas = 8.314510

      !sam reference temperature T0 = 300 K, a change from original 298K
      REAL, PARAMETER :: tnull = 300.

      !bs molecular weight for C
      REAL, PARAMETER :: mwc = 12.0
      !bs molecular weight for organic species
      REAL, PARAMETER :: mworg = 175.0
      !bs molecular weight for SO4
      REAL, PARAMETER :: mwso4 = 96.0576
      !bs molecular weight for NH4
      REAL, PARAMETER :: mwnh4 = 18.03858
      !bs molecular weight for NO3
      REAL, PARAMETER :: mwno3 = 62.01287
      ! molecular weight for AIR

!     REAL mwair
!     PARAMETER (mwair=28.964)
      !bs relative tolerance for mass check
      REAL, PARAMETER :: CABSMIN = .00001 ! Minimum amount of absorbing material - needed in iteration method
      !sm number of basis set variables in CMU partitioning scheme
      INTEGER, PARAMETER :: nbin = 4 ! we use 4 bin volatility according to Robinson A. et al.

      ! we have 2 type of SOA - anthropogenic and biogenic
      !sm number of SAPRC species variables in CMU lumped partitioning table
      !sm 1=ALK4(hc5),2=ALK5(hc8),3=OLE1(ol2),4=OLE2(oli),5=ARO1(tol)
      !sm 6=AOR2(xyl),7=ISOP(iso),8=SESQ(?),9=TERP(alp)
      INTEGER, PARAMETER :: nsaprc = 9 ! number of precursor classes

      !bs loop indices
      INTEGER lcell, n, l, ll, bn, cls
      !bs conversion factor ppm --> ug/m^3
      REAL convfac
      !bs difference of inverse temperatures
      REAL ttinv
      !bs initial organic absorbing mass [ug/m^3]
      REAL minit
      !bs inorganic mass [ug/m^3]
      REAL mnono
      !bs total organic mass [ug/m^3]
      REAL mtot

!     REAL msum(ncv)  !bs input total mass [ug/m^3]
      REAL mwcv(ncv) !bs molecular weight of cond. vapors [g/
      REAL pnull(ncv) !bs vapor pres. of pure cond. vapor [Pa]
      REAL dhvap(ncv) !bs heat of vaporisation of compound i [
      REAL pvap(ncv) !bs vapor pressure cond. vapor [Pa]
      REAL ctot(ncv) !bs total conc. of cond. vapor aerosol +
      REAL cgas(ncv) !bs gasphase concentration of cond. vapors
      REAL caer(ncv) !bs aerosolphase concentration of cond.
      REAL asav(ncv) !bs saved CAER for iteration
      REAL aold(ncv) !bs saved CAER for rate determination
      REAL csat(ncv) !bs saturation conc. of cond. vapor ug/,

      ! in basis set approach we need only 4 csat
      REAL ccsat(nbin)
      REAL ccaer(nbin)
      REAL cctot(nbin)
      REAL w1(nbin), w2(nbin)

      REAL prod(ncv) !bs production of condensable vapor ug/
      REAL p(ncv) !bs PROD(L) * TIMEFAC [ug/m^3]
      REAL f(ldrog_vbs) !bs scaling factor for ind. oxidant

      REAL alphlowN(nbin, nsaprc) ! sm normalized (1 g/m3 density) yield for condensable vapors low NOx condition
      REAL alphhiN(nbin, nsaprc) ! sm normalized (1 g/m3 density) yield for condensable vapors high NOx condition
      REAL alphai(nbin, nsaprc) ! mass-based stoichometric yield for product i and csti is the effective saturation
      ! concentration in ug m^-3
      REAL mwvoc(nsaprc) ! molecular weight of the SOA precusors

      REAL PnGtotal, DUM, FMTOT, FMTOT2, DUM2 ! Real constants used in Newton iteration
      integer, save :: icall


      RETURN
   END SUBROUTINE soa_vbs
!
! *** this routine calculates the dry deposition and sedimentation
!     velocities for the three modes.
!     coded 1/23/97 by Dr. Francis S. Binkowski. Follows
!     FSB's original method, i.e. uses Jon Pleim's expression for deposition
!     velocity but includes Marv Wesely's wstar contribution.
!ia eliminated Stokes term for coarse mode deposition calcs.,
!ia see comments below

   SUBROUTINE VDVG(BLKSIZE, NSPCSDA, NUMCELLS, &
                   LAYER, &
                   CBLK, &
                   BLKTA, BLKDENS, RA, USTAR, WSTAR, AMU, &
                   DGNUC, DGACC, DGCOR, &
                   KNNUC, KNACC, KNCOR, &
                   PDENSN, PDENSA, PDENSC, &
                   VSED, VDEP)

! *** calculate size-averaged particle dry deposition and
!     size-averaged sedimentation velocities.

!     IMPLICIT NONE

      INTEGER BLKSIZE ! dimension of arrays
      INTEGER NSPCSDA ! number of species in CBLK
      INTEGER NUMCELLS ! actual number of cells in arrays
      INTEGER LAYER ! number of layer

      REAL CBLK(BLKSIZE, NSPCSDA) ! main array of variables
      REAL BLKTA(BLKSIZE) ! Air temperature [ K ]
      REAL BLKDENS(BLKSIZE) ! Air density  [ kg m^-3 ]
      REAL RA(BLKSIZE) ! aerodynamic resistance [ s m**-1 ]
      REAL USTAR(BLKSIZE) ! surface friction velocity [ m s**-1 ]
      REAL WSTAR(BLKSIZE) ! convective velocity scale [ m s**-1 ]
      REAL AMU(BLKSIZE) ! atmospheric dynamic viscosity [ kg m**-1 s**-1 ]
      REAL DGNUC(BLKSIZE) ! nuclei mode mean diameter [ m ]
      REAL DGACC(BLKSIZE) ! accumulation
      REAL DGCOR(BLKSIZE) ! coarse mode
      REAL KNNUC(BLKSIZE) ! nuclei mode Knudsen number
      REAL KNACC(BLKSIZE) ! accumulation
      REAL KNCOR(BLKSIZE) ! coarse mode
      REAL PDENSN(BLKSIZE) ! average particel density in nuclei mode [ kg / m**3 ]
      REAL PDENSA(BLKSIZE) ! average particel density in accumulation mode [ kg / m**3 ]
      REAL PDENSC(BLKSIZE) ! average particel density in coarse mode [ kg / m**3 ]

! *** modal particle diffusivities for number and 3rd moment, or mass:

      REAL DCHAT0N(BLKSIZE), DCHAT0A(BLKSIZE), DCHAT0C(BLKSIZE)
      REAL DCHAT3N(BLKSIZE), DCHAT3A(BLKSIZE), DCHAT3C(BLKSIZE)

! *** modal sedimentation velocities for number and 3rd moment, or mass:

      REAL VGHAT0N(BLKSIZE), VGHAT0A(BLKSIZE), VGHAT0C(BLKSIZE)
      REAL VGHAT3N(BLKSIZE), VGHAT3A(BLKSIZE), VGHAT3C(BLKSIZE)

! *** deposition and sedimentation velocities

      REAL VDEP(BLKSIZE, NASPCSDEP) ! sedimantation velocity [ m s**-1 ]
      REAL VSED(BLKSIZE, NASPCSSED) ! deposition  velocity [ m s**-1 ]

      INTEGER LCELL
      REAL DCONST1, DCONST1N, DCONST1A, DCONST1C
      REAL DCONST2, DCONST3N, DCONST3A, DCONST3C
      REAL SC0N, SC0A, SC0C ! Schmidt numbers for number
      REAL SC3N, SC3A, SC3C ! Schmidt numbers for 3rd moment
      REAL ST0N, ST0A, ST0C ! Stokes numbers for number
      REAL ST3N, ST3A, ST3C ! Stokes numbers for 3rd moment
      REAL RD0N, RD0A, RD0C ! canopy resistance for number
      REAL RD3N, RD3A, RD3C ! canopy resisteance for 3rd moment
      REAL UTSCALE ! scratch function of USTAR and WSTAR.
      REAL NU !kinematic viscosity [ m**2 s**-1 ]
      REAL USTFAC ! scratch function of USTAR, NU, and GRAV
      REAL BHAT
      PARAMETER(BHAT=1.246) ! Constant from Cunningham slip correction.

! *** check layer value.

      IF (LAYER .EQ. 1) THEN ! calculate diffusitities and
!                                    sedimentation velocities

         DO LCELL = 1, NUMCELLS

            DCONST1 = BOLTZ*BLKTA(LCELL)/ &
                      (THREEPI*AMU(LCELL))
            DCONST1N = DCONST1/DGNUC(LCELL)
            DCONST1A = DCONST1/DGACC(LCELL)
            DCONST1C = DCONST1/DGCOR(LCELL)
            DCONST2 = GRAV/(18.0*AMU(LCELL))
            DCONST3N = DCONST2*PDENSN(LCELL)*DGNUC(LCELL)**2
            DCONST3A = DCONST2*PDENSA(LCELL)*DGACC(LCELL)**2
            DCONST3C = DCONST2*PDENSC(LCELL)*DGCOR(LCELL)**2

! *** i-mode

            DCHAT0N(LCELL) = DCONST1N &
                             *(ESN04+BHAT*KNNUC(LCELL)*ESN16)

            DCHAT3N(LCELL) = DCONST1N &
                             *(ESNM20+BHAT*KNNUC(LCELL)*ESNM32)

            VGHAT0N(LCELL) = DCONST3N &
                             *(ESN16+BHAT*KNNUC(LCELL)*ESN04)

            VGHAT3N(LCELL) = DCONST3N &
                             *(ESN64+BHAT*KNNUC(LCELL)*ESN28)

! *** j-mode

            DCHAT0A(LCELL) = DCONST1A &
                             *(ESA04+BHAT*KNACC(LCELL)*ESA16)

            DCHAT3A(LCELL) = DCONST1A &
                             *(ESAM20+BHAT*KNACC(LCELL)*ESAM32)

            VGHAT0A(LCELL) = DCONST3A &
                             *(ESA16+BHAT*KNACC(LCELL)*ESA04)

            VGHAT3A(LCELL) = DCONST3A &
                             *(ESA64+BHAT*KNACC(LCELL)*ESA28)

! *** coarse mode

            DCHAT0C(LCELL) = DCONST1C &
                             *(ESC04+BHAT*KNCOR(LCELL)*ESC16)

            DCHAT3C(LCELL) = DCONST1C &
                             *(ESCM20+BHAT*KNCOR(LCELL)*ESCM32)

            VGHAT0C(LCELL) = DCONST3C &
                             *(ESC16+BHAT*KNCOR(LCELL)*ESC04)

            VGHAT3C(LCELL) = DCONST3C &
                             *(ESC64+BHAT*KNCOR(LCELL)*ESC28)

         END DO

! *** now calculate the deposition and sedmentation velocities

!ia  07.05.98
! *** NOTE In the deposition velocity for coarse mode,
!     the impaction term  10.0 ** (-3.0 / st) is eliminated because
!     coarse particles are likely to bounce on impact and the current
!     formulation does not account for this.

         DO LCELL = 1, NUMCELLS

            NU = AMU(LCELL)/BLKDENS(LCELL)
            USTFAC = USTAR(LCELL)*USTAR(LCELL)/(GRAV*NU)
            UTSCALE = USTAR(LCELL) + &
                      0.24*WSTAR(LCELL)*WSTAR(LCELL)/USTAR(LCELL)

! *** first do number

! *** nuclei or Aitken mode  ( no sedimentation velocity )

            SC0N = NU/DCHAT0N(LCELL)
            ST0N = MAX(VGHAT0N(LCELL)*USTFAC, 0.01)
            RD0N = 1.0/(UTSCALE* &
                        (SC0N**(-TWO3) + 10.0**(-3.0/ST0N)))

            VDEP(LCELL, VDNNUC) = VGHAT0N(LCELL) + &
                                  1.0/( &
                                  RA(LCELL) + RD0N+RD0N*RA(LCELL)*VGHAT0N(LCELL))

            VSED(LCELL, VSNNUC) = VGHAT0N(LCELL)

! *** accumulation mode

            SC0A = NU/DCHAT0A(LCELL)
            ST0A = MAX(VGHAT0A(LCELL)*USTFAC, 0.01)
            RD0A = 1.0/(UTSCALE* &
                        (SC0A**(-TWO3) + 10.0**(-3.0/ST0A)))

            VDEP(LCELL, VDNACC) = VGHAT0A(LCELL) + &
                                  1.0/( &
                                  RA(LCELL) + RD0A+RD0A*RA(LCELL)*VGHAT0A(LCELL))

            VSED(LCELL, VSNACC) = VGHAT0A(LCELL)

! *** coarse mode

            SC0C = NU/DCHAT0C(LCELL)
!ia        ST0C = MAX( VGHAT0C(LCELL) * USTFAC, 0.01 )
!ia        RD0C = 1.0 / ( UTSCALE *
!ia     &            ( SC0C**(-TWO3) + 10.0**(-3.0 / ST0C) ) )

            RD0C = 1.0/(UTSCALE* &
                        (SC0C**(-TWO3))) ! eliminate impaction term

            VDEP(LCELL, VDNCOR) = VGHAT0C(LCELL) + &
                                  1.0/( &
                                  RA(LCELL) + RD0C+RD0C*RA(LCELL)*VGHAT0C(LCELL))

            VSED(LCELL, VSNCOR) = VGHAT0C(LCELL)

! *** now do m3 for the deposition of mass

! *** nuclei or Aitken mode

            SC3N = NU/DCHAT3N(LCELL)
            ST3N = MAX(VGHAT3N(LCELL)*USTFAC, 0.01)
            RD3N = 1.0/(UTSCALE* &
                        (SC3N**(-TWO3) + 10.0**(-3.0/ST3N)))

            VDEP(LCELL, VDMNUC) = VGHAT3N(LCELL) + &
                                  1.0/( &
                                  RA(LCELL) + RD3N+RD3N*RA(LCELL)*VGHAT3N(LCELL))

            VSED(LCELL, VSMNUC) = VGHAT3N(LCELL)

! *** accumulation mode

            SC3A = NU/DCHAT3A(LCELL)
            ST3A = MAX(VGHAT3A(LCELL)*USTFAC, 0.01)
            RD3A = 1.0/(UTSCALE* &
                        (SC3A**(-TWO3) + 10.0**(-3.0/ST3A)))

            VDEP(LCELL, VDMACC) = VGHAT3A(LCELL) + &
                                  1.0/( &
                                  RA(LCELL) + RD3A+RD3A*RA(LCELL)*VGHAT3A(LCELL))

! *** fine mass deposition velocity: combine Aitken and accumulation
!     mode deposition velocities. Assume density is the same
!     for both modes.

!       VDEP(LCELL,VDMFINE) = (
!    &    CBLK(LCELL,VNU3) * VDEP(LCELL, VDMNUC) +
!    &    CBLK(LCELL,VAC3) * VDEP(LCELL, VDMACC) ) /
!    &    ( CBLK(LCELL,VAC3) + CBLK(LCELL,VNU3) )

! *** fine mass sedimentation velocity

!       VSED( LCELL, VSMFINE) = (
!    &    CBLK(LCELL, VNU3) * VGHAT3N(LCELL) +
!    &     CBLK(LCELL, VAC3) * VGHAT3A(LCELL) ) /
!    &    ( CBLK(LCELL, VNU3) + CBLK(LCELL, VAC3)  )

            VSED(LCELL, VSMACC) = VGHAT3A(LCELL)

! *** coarse mode

            SC3C = NU/DCHAT3C(LCELL)
!ia        ST3C = MAX( VGHAT3C(LCELL) * USTFAC, 0.01 )
!ia        RD3C = 1.0 / ( UTSCALE *
!ia     &            ( SC3C**(-TWO3) + 10.0**(-3.0 / ST3C) ) )

            RD3C = 1.0/(UTSCALE* &
                        (SC3C**(-TWO3))) ! eliminate impaction term
            VDEP(LCELL, VDMCOR) = VGHAT3C(LCELL) + &
                                  1.0/( &
                                  RA(LCELL) + RD3C+RD3C*RA(LCELL)*VGHAT3C(LCELL))

! *** coarse mode sedmentation velocity

            VSED(LCELL, VSMCOR) = VGHAT3C(LCELL)

         END DO

      ELSE ! LAYER greater than 1

! *** for layer greater than 1 calculate  sedimentation velocities only

         DO LCELL = 1, NUMCELLS

            DCONST2 = GRAV/(18.0*AMU(LCELL))

            DCONST3N = DCONST2*PDENSN(LCELL)*DGNUC(LCELL)**2
            DCONST3A = DCONST2*PDENSA(LCELL)*DGACC(LCELL)**2
            DCONST3C = DCONST2*PDENSC(LCELL)*DGCOR(LCELL)**2

            VGHAT0N(LCELL) = DCONST3N &
                             *(ESN16+BHAT*KNNUC(LCELL)*ESN04)

! *** nucleation mode number sedimentation velocity

            VSED(LCELL, VSNNUC) = VGHAT0N(LCELL)

            VGHAT3N(LCELL) = DCONST3N &
                             *(ESN64+BHAT*KNNUC(LCELL)*ESN28)

! *** nucleation mode volume sedimentation velocity

            VSED(LCELL, VSMNUC) = VGHAT3N(LCELL)

            VGHAT0A(LCELL) = DCONST3A &
                             *(ESA16+BHAT*KNACC(LCELL)*ESA04)

! *** accumulation mode number sedimentation velocity

            VSED(LCELL, VSNACC) = VGHAT0A(LCELL)

            VGHAT3A(LCELL) = DCONST3A &
                             *(ESA64+BHAT*KNACC(LCELL)*ESA28)

! *** fine mass sedimentation velocity

!           VSED( LCELL, VSMFINE) = (
!    &       CBLK(LCELL, VNU3) * VGHAT3N(LCELL) +
!    &        CBLK(LCELL, VAC3) * VGHAT3A(LCELL) ) /
!    &       ( CBLK(LCELL, VNU3) + CBLK(LCELL, VAC3)  )

            VSED(LCELL, VSMACC) = VGHAT3A(LCELL)

            VGHAT0C(LCELL) = DCONST3C &
                             *(ESC16+BHAT*KNCOR(LCELL)*ESC04)

! *** coarse mode sedimentation velocity

            VSED(LCELL, VSNCOR) = VGHAT0C(LCELL)

            VGHAT3C(LCELL) = DCONST3C &
                             *(ESC64+BHAT*KNCOR(LCELL)*ESC28)

! *** coarse mode mass sedimentation velocity

            VSED(LCELL, VSMCOR) = VGHAT3C(LCELL)

         END DO

      END IF ! check on layer

   END SUBROUTINE VDVG
!
!---------------------------------------------------------------------------
!
! *** this routine calculates the dry deposition and sedimentation
!     velocities for the three modes.
!   Stu McKeen 10/13/08
!   Gaussian Quadrature numerical integration over diameter range for each mode.
! Quadrature taken from Abramowitz and Stegun (1974), equation 25.4.46 and Table 25.10
! Quadrature points are the zeros of Hermite polynomials of order NGAUSdv
!   Numerical Integration allows more complete discription of the
!   Cunningham Slip correction factor, Interception Term (not included previously),
!   and the correction due to rebound for higher diameter particles.
!   Sedimentation velocities the same as original Binkowski code, also the
!   Schmidt number and Brownian diffusion efficiency dependence on Schmidt number the
!   same as Binkowski.
!   Stokes number, and efficiency dependence on Stokes number now according to
!   Peters and Eiden (1992).  Interception term taken from Slinn (1982) with
!   efficiency at .2 micron diam. (0.3%) tuned to yield .2 cm/s deposition velocitiy
!   for needleaf evergreen trees (Pryor et al., Tellus, 2008). Rebound correction
!   term is that of Slinn (1982)
!
!     Original code 1/23/97 by Dr. Francis S. Binkowski. Follows
!     FSB's original method, i.e. uses Jon Pleim's expression for deposition
!     velocity but includes Marv Wesely's wstar contribution.
!ia eliminated Stokes term for coarse mode deposition calcs.,
!ia see comments below

! CBLK is eliminated since the subroutine doesn't use it!
   SUBROUTINE VDVG_2(BLKSIZE, NSPCSDA, NUMCELLS, &
                     LAYER, &
                     BLKTA, BLKDENS, &
                     RA, USTAR, PBLH, ZNTT, RMOLM, AMU, &
                     DGNUC, DGACC, DGCOR, XLM, &
                     KNNUC, KNACC, KNCOR, &
                     PDENSN, PDENSA, PDENSC, &
                     VSED, VDEP)

! *** calculate size-averaged particle dry deposition and
!     size-averaged sedimentation velocities.
!     IMPLICIT NONE

      INTEGER BLKSIZE ! dimension of arrays
      INTEGER NSPCSDA ! number of species in CBLK
      INTEGER NUMCELLS ! actual number of cells in arrays
      INTEGER LAYER ! number of layer
      INTEGER, PARAMETER :: iprnt = 0

!     REAL CBLK( BLKSIZE, NSPCSDA ) ! main array of variables
      REAL BLKTA(BLKSIZE) ! Air temperature [ K ]
      REAL BLKDENS(BLKSIZE) ! Air density  [ kg m^-3 ]
      REAL RA(BLKSIZE) ! aerodynamic resistance [ s m**-1 ]
      REAL USTAR(BLKSIZE) ! surface friction velocity [ m s**-1 ]
      REAL PBLH(BLKSIZE) ! PBL height (m)
      REAL ZNTT(BLKSIZE) ! Surface roughness length (m)
      REAL RMOLM(BLKSIZE) ! Inverse of Monin-Obukhov length (1/m)
      REAL AMU(BLKSIZE) ! atmospheric dynamic viscosity [ kg m**-1 s**-1 ]
      REAL XLM(BLKSIZE) ! mean free path of dry air [ m ]
      REAL DGNUC(BLKSIZE) ! nuclei mode mean diameter [ m ]
      REAL DGACC(BLKSIZE) ! accumulation
      REAL DGCOR(BLKSIZE) ! coarse mode
      REAL KNNUC(BLKSIZE) ! nuclei mode Knudsen number
      REAL KNACC(BLKSIZE) ! accumulation
      REAL KNCOR(BLKSIZE) ! coarse mode
      REAL PDENSN(BLKSIZE) ! average particle density in nuclei mode [ kg / m**3 ]
      REAL PDENSA(BLKSIZE) ! average particle density in accumulation mode [ kg / m**3 ]
      REAL PDENSC(BLKSIZE) ! average particle density in coarse mode [ kg / m**3 ]

! *** deposition and sedimentation velocities

      REAL VDEP(BLKSIZE, NASPCSDEP) ! sedimentation velocity [ m s**-1 ]
      REAL VSED(BLKSIZE, NASPCSSED) ! deposition  velocity [ m s**-1 ]

      INTEGER LCELL, N
      REAL DCONST1, DCONST2, DCONST3, DCONST3N, DCONST3A, DCONST3C
      REAL UTSCALE, CZH ! scratch functions of USTAR and WSTAR.
      REAL NU !kinematic viscosity [ m**2 s**-1 ]
      REAL BHAT
      PARAMETER(BHAT=1.246) ! Constant from Binkowski-Shankar approx to Cunningham slip correction.
      REAL COLCTR_BIGD, COLCTR_SMALD
      PARAMETER(COLCTR_BIGD=2.E-3, COLCTR_SMALD=20.E-6) ! Collector diameters in Stokes number and Interception Efficiency (Needleleaf Forest)
      REAL SUM0, SUM3, DQ, KNQ, CUNQ, VSEDQ, SCQ, STQ, RSURFQ, vdplim
      REAL Eff_dif, Eff_imp, Eff_int, RBcor
      INTEGER ISTOPvd0, IdoWesCor
      PARAMETER(ISTOPvd0=0) ! ISTOPvd0 = 1 means dont call VDVG, particle dep. velocities are set = 0; ISTOPvd0 = 0 means do depvel calcs.

      ! no Wesley deposition, otherwise EC is too low
      PARAMETER(IdoWesCor=0) ! IdoWesCor = 1 means do Wesley (85) convective correction to PM dry dep velocities; 0 means don't do correction
      IF (ISTOPvd0 .EQ. 1) THEN
         RETURN
      ENDIF
! *** check layer value.

      IF (iprnt .eq. 1) print *, 'In VDVG, Layer=', LAYER
      IF (LAYER .EQ. 1) THEN ! calculate diffusitities and sedimentation velocities

         DO LCELL = 1, NUMCELLS
            DCONST1 = BOLTZ*BLKTA(LCELL)/ &
                      (THREEPI*AMU(LCELL))
            DCONST2 = GRAV/(18.0*AMU(LCELL))
            DCONST3 = USTAR(LCELL)/(9.*AMU(LCELL)*COLCTR_BIGD)

! *** now calculate the deposition velocities at layer 1

            NU = AMU(LCELL)/BLKDENS(LCELL)

            UTSCALE = 1.
            IF (IdoWesCor .EQ. 1) THEN
! Wesley (1985) Monin-Obukov dependence for convective conditions (SAM 10/08)
               IF (RMOLM(LCELL) .LT. 0.) THEN
                  CZH = -1.*PBLH(LCELL)*RMOLM(LCELL)
                  IF (CZH .GT. 30.0) THEN
                     UTSCALE = 0.45*CZH**0.6667
                  ELSE
                     UTSCALE = 1.+(-300.*RMOLM(LCELL))**0.6667
                  ENDIF
               ENDIF
            ENDIF ! end of (IdoWesCor.EQ.1) test

            UTSCALE = USTAR(LCELL)*UTSCALE
            IF (iprnt .eq. 1) THEN
               print *, 'NGAUSdv,xxlsga,USTAR,UTSCALE'
               print *, NGAUSdv, xxlsga, USTAR(LCELL), UTSCALE
               print *, 'DCONST2,PDENSA,DGACC,GRAV,AMU'
               print *, DCONST2, PDENSA(LCELL), DGACC(LCELL), GRAV, AMU(LCELL)
            endif

! *** nuclei mode

            SUM0 = 0.
            SUM3 = 0.
            DO N = 1, NGAUSdv
               DQ = DGNUC(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsgn) ! Diameter (m) at quadrature point
               KNQ = 2.*XLM(LCELL)/DQ ! Knudsen number at quadrature point
               CUNQ = 1.+KNQ*(1.257+.4*exp(-1.1/KNQ)) ! Cunningham correction factor; Pruppacher and Klett (1980) Eq (12-16)
               VSEDQ = PDENSN(LCELL)*DCONST2*CUNQ*DQ*DQ ! Gravitational sedimentation velocity m/s
               SCQ = NU*DQ/DCONST1/CUNQ ! Schmidt number, Brownian diffusion parameter - Same as Binkowski and Shankar
               Eff_dif = SCQ**(-TWO3) ! Efficiency term for diffusion - Same as Binkowski and Shankar
               STQ = DCONST3*PDENSN(LCELL)*DQ**2 ! Stokes number, Peters and Eiden (1992)
               Eff_imp = (STQ/(0.8+STQ))**2 ! Efficiency term for impaction - Peters and Eiden (1992)
               !       Eff_int=0.3*DQ/(COLCTR_SMALD+DQ) ! Slinn (1982) Interception term, 0.3 prefac insures .2 cm/s at .2 micron diam.
               Eff_int = (0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 ! McKeen(2008) Intercptn trm, val of .00421 @ ustr=0.475, diam=.1414 micrn, stable, needleleaf evergreen
               RBcor = 1. ! Rebound correction factor
               vdplim = UTSCALE*(Eff_dif + Eff_imp + Eff_int)*RBcor
               !       vdplim=.002*UTSCALE
               vdplim = min(vdplim, .02)
               RSURFQ = RA(LCELL) + 1./vdplim
               !       RSURFQ=RA(LCELL)+1./(UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor) ! Total surface resistence
               !
!   limit this here to be consisten with the gocart routine, which bases this on Walcek et al. 1986
!
               !       RSURFQ=max(RSURFQ,50.)
               SUM0 = SUM0 + WGAUS(N)*(VSEDQ + 1./RSURFQ) ! Quadrature sum for 0 moment
               SUM3 = SUM3 + WGAUS(N)*(VSEDQ + 1./RSURFQ)*DQ**3 ! Quadrature sum for 3rd moment
            ENDDO
            VDEP(LCELL, VDNNUC) = SUM0/sqrtpi ! normalize 0 moment vdep quadrature sum to sqrt(pi) (and number =1 per unit volume)
            VDEP(LCELL, VDMNUC) = SUM3/(sqrtpi*EXP((1.5*sqrt2*xxlsgn)**2)*DGNUC(LCELL)**3) !normalize 3 moment quad. sum to sqrt(pi) and 3rd moment analytic sum

! *** accumulation mode

            SUM0 = 0.
            SUM3 = 0.
            DO N = 1, NGAUSdv
               DQ = DGACC(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsga) ! Diameter (m) at quadrature point
               KNQ = 2.*XLM(LCELL)/DQ ! Knudsen number at quadrature point
               CUNQ = 1.+KNQ*(1.257+.4*exp(-1.1/KNQ)) ! Cunningham correction factor; Pruppacher and Klett (1980) Eq (12-16)
               VSEDQ = PDENSA(LCELL)*DCONST2*CUNQ*DQ*DQ ! Gravitational sedimentation velocity m/s
               SCQ = NU*DQ/DCONST1/CUNQ ! Schmidt number, Brownian diffusion parameter - Same as Binkowski and Shankar
               Eff_dif = SCQ**(-TWO3) ! Efficiency term for diffusion - Same as Binkowski and Shankar
               STQ = DCONST3*PDENSA(LCELL)*DQ**2 ! Stokes number, Peters and Eiden (1992)
               Eff_imp = (STQ/(0.8+STQ))**2 ! Efficiency term for impaction - Peters and Eiden (1992)
               !       Eff_int=0.3*DQ/(COLCTR_SMALD+DQ) ! Slinn (1982) Interception term, 0.3 prefac insures .2 cm/s at .2 micron diam.
               Eff_int = (0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 ! McKeen(2008) Intercptn term, val of .00421 @ ustr=0.475, diam=.1414 micrn, stable, needleleaf evergreen
               RBcor = 1. ! Rebound correction factor
               vdplim = UTSCALE*(Eff_dif + Eff_imp + Eff_int)*RBcor
               vdplim = min(vdplim, .02)
               RSURFQ = RA(LCELL) + 1./vdplim
!       RSURFQ=RA(LCELL)+1./(UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor) ! Total surface resistence
!
!   limit this here to be consisten with the gocart routine, which bases this on Walcek et al. 1986
!
!       RSURFQ=max(RSURFQ,50.)
               SUM0 = SUM0 + WGAUS(N)*(VSEDQ + 1./RSURFQ) ! Quadrature sum for 0 moment
               SUM3 = SUM3 + WGAUS(N)*(VSEDQ + 1./RSURFQ)*DQ**3 ! Quadrature sum for 3rd moment
               IF (iprnt .eq. 1) THEN
                  print *, 'N,Y_GQ,WGAUS,DQ,KNQ,CUNQ,VSEDQ,SCQ,STQ,RSURFQ'
                  print *, N, Y_GQ(N), WGAUS(N), DQ, KNQ, CUNQ, VSEDQ, SCQ, STQ, RSURFQ
                  print *, 'N,Eff_dif,imp,int,SUM0,SUM3'
                  print *, N, Eff_dif, Eff_imp, Eff_int, SUM0, SUM3
               endif
            ENDDO
            VDEP(LCELL, VDNACC) = SUM0/sqrtpi ! normalize 0 moment vdep quadrature sum to sqrt(pi) (and number =1 per unit volume)
            VDEP(LCELL, VDMACC) = SUM3/(sqrtpi*EXP((1.5*sqrt2*xxlsga)**2)*DGACC(LCELL)**3) !normalize 3 moment quad. sum to sqrt(pi) and 3rd moment analytic sum

! *** coarse mode

            SUM0 = 0.
            SUM3 = 0.
            DO N = 1, NGAUSdv
               DQ = DGCOR(LCELL)*EXP(Y_GQ(N)*sqrt2*xxlsgc) ! Diameter (m) at quadrature point
               KNQ = 2.*XLM(LCELL)/DQ ! Knudsen number at quadrature point
               CUNQ = 1.+KNQ*(1.257+.4*exp(-1.1/KNQ)) ! Cunningham correction factor; Pruppacher and Klett (1980) Eq (12-16)
               VSEDQ = PDENSC(LCELL)*DCONST2*CUNQ*DQ*DQ ! Gravitational sedimentation velocity m/s
               SCQ = NU*DQ/DCONST1/CUNQ ! Schmidt number, Brownian diffusion parameter - Same as Binkowski and Shankar
               Eff_dif = SCQ**(-TWO3) ! Efficiency term for diffusion - Same as Binkowski and Shankar
               STQ = DCONST3*PDENSC(LCELL)*DQ**2 ! Stokes number, Peters and Eiden (1992)
               Eff_imp = (STQ/(0.8+STQ))**2 ! Efficiency term for impaction - Peters and Eiden (1992)
!          Eff_int=0.3*DQ/(COLCTR_SMALD+DQ) ! Slinn (1982) Interception term, 0.3 prefac insures .2 cm/s at .2 micron diam.
               Eff_int = (0.00116+0.0061*ZNTT(LCELL))*DQ/1.414E-7 ! McKeen(2008) Interception term, val of .00421 @ ustr=0.475, diam=.1414 micrn, stable, needleleaf evergreen
               EFF_int = min(1., EFF_int)
               RBcor = exp(-2.0*(STQ**0.5)) ! Rebound correction factor used in Slinn (1982)
               vdplim = UTSCALE*(Eff_dif + Eff_imp + Eff_int)*RBcor
               vdplim = min(vdplim, .02)
               RSURFQ = RA(LCELL) + 1./vdplim
!       RSURFQ=RA(LCELL)+1./(UTSCALE*(Eff_dif+Eff_imp+Eff_int)*RBcor) ! Total surface resistence
!
!   limit this here to be consisten with the gocart routine, which bases this on Walcek et al. 1986
!
!       RSURFQ=max(RSURFQ,50.)
               SUM0 = SUM0 + WGAUS(N)*(VSEDQ + 1./RSURFQ) ! Quadrature sum for 0 moment
               SUM3 = SUM3 + WGAUS(N)*(VSEDQ + 1./RSURFQ)*DQ**3 ! Quadrature sum for 3rd moment
            ENDDO
            VDEP(LCELL, VDNCOR) = SUM0/sqrtpi ! normalize 0 moment vdep quadrature sum to sqrt(pi) (and number =1 per unit volume)
            VDEP(LCELL, VDMCOR) = SUM3/(sqrtpi*EXP((1.5*sqrt2*xxlsgc)**2)*DGCOR(LCELL)**3) !normalize 3 moment quad. sum to sqrt(pi) and 3rd moment analytic sum
         END DO

      ENDIF ! ENDOF LAYER = 1 test

! *** Calculate gravitational sedimentation velocities for all layers - as in Binkowski and Shankar (1995)

      DO LCELL = 1, NUMCELLS

         DCONST2 = GRAV/(18.0*AMU(LCELL))
         DCONST3N = DCONST2*PDENSN(LCELL)*DGNUC(LCELL)**2
         DCONST3A = DCONST2*PDENSA(LCELL)*DGACC(LCELL)**2
         DCONST3C = DCONST2*PDENSC(LCELL)*DGCOR(LCELL)**2

! *** nucleation mode number and mass sedimentation velociticies
         VSED(LCELL, VSNNUC) = DCONST3N &
                               *(ESN16+BHAT*KNNUC(LCELL)*ESN04)
         VSED(LCELL, VSMNUC) = DCONST3N &
                               *(ESN64+BHAT*KNNUC(LCELL)*ESN28)

! *** accumulation mode number and mass sedimentation velociticies
         VSED(LCELL, VSNACC) = DCONST3A &
                               *(ESA16+BHAT*KNACC(LCELL)*ESA04)
         VSED(LCELL, VSMACC) = DCONST3A &
                               *(ESA64+BHAT*KNACC(LCELL)*ESA28)

! *** coarse mode number and mass sedimentation velociticies
         VSED(LCELL, VSNCOR) = DCONST3C &
                               *(ESC16+BHAT*KNCOR(LCELL)*ESC04)
         VSED(LCELL, VSMCOR) = DCONST3C &
                               *(ESC64+BHAT*KNCOR(LCELL)*ESC28)
      END DO
   END SUBROUTINE VDVG_2
!------------------------------------------------------------------------------

   SUBROUTINE aerosols_soa_vbs_init(chem, convfac, z_at_w, &
                                    pm2_5_dry, pm2_5_water, pm2_5_dry_ec, &
                                    chem_in_opt, aer_ic_opt, is_aerosol, &
                                    ids, ide, jds, jde, kds, kde, &
                                    ims, ime, jms, jme, kms, kme, &
                                    its, ite, jts, jte, kts, kte, config_flags)

      USE module_configure, only:grid_config_rec_type

      implicit none
      INTEGER, INTENT(IN) ::  chem_in_opt, aer_ic_opt
      INTEGER, INTENT(IN) ::  ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      LOGICAL, INTENT(OUT) :: is_aerosol(num_chem)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(INOUT) :: &
         chem
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(INOUT) :: &
         pm2_5_dry, pm2_5_water, pm2_5_dry_ec
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: &
         convfac
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: &
         z_at_w
      TYPE(grid_config_rec_type), INTENT(in) :: config_flags

      integer i, j, k, l, ii, jj, kk
      real tempfac, mwso4, zz
!    real,dimension(its:ite,kts:kte,jts:jte) :: convfac
      REAL splitfac
      !between gas and aerosol phase
      REAL so4vaptoaer
!factor for splitting initial conc. of SO4
!3rd moment i-mode [3rd moment/m^3]
      REAL m3nuc
!3rd MOMENT j-mode [3rd moment/m^3]
      REAL m3acc
!       REAL ESN36
      REAL m3cor
      DATA splitfac/.98/
      DATA so4vaptoaer/.999/


      return
   END SUBROUTINE aerosols_soa_vbs_init

!
   SUBROUTINE soa_vbs_addemiss(id, dtstep, u10, v10, alt, dz8w, xland, chem, &
                               ebu, &
                               slai, ust, smois, ivgtyp, isltyp, &
                               emis_ant, dust_emiss_active, &
                               seasalt_emiss_active, kemit, biom, num_soil_layers, emissopt, &
                               dust_opt, ktau, p8w, u_phy, v_phy, rho_phy, g, dx, erod, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte)

      USE module_state_description, only:num_chem

      INTEGER, INTENT(IN)   ::    seasalt_emiss_active, kemit, emissopt, &
                               dust_emiss_active, num_soil_layers, id, &
                               ktau, dust_opt, biom, &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte

      REAL, INTENT(IN) ::    dtstep

! trace species mixing ratios (aerosol mass = ug/kg-air; number = #/kg-air)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(INOUT) ::   chem
!
! aerosol emissions arrays ((ug/m3)*m/s)
!
      REAL, DIMENSION(ims:ime, kms:kemit, jms:jme, num_emis_ant), &
         INTENT(IN) ::    emis_ant

! biomass burning aerosol emissions arrays ((ug/m3)*m/s)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_ebu), &
         INTENT(IN) ::    ebu

! 1/(dry air density) and layer thickness (m)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: &
         alt, dz8w

      ! add for gocart dust
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: p8w, u_phy, v_phy, rho_phy
      REAL, INTENT(IN) :: dx, g
      REAL, DIMENSION(ims:ime, jms:jme, 3), &
         INTENT(IN) :: erod

      REAL, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN) :: &
         u10, v10, xland, slai, ust
      INTEGER, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN) ::   ivgtyp, isltyp
      REAL, DIMENSION(ims:ime, num_soil_layers, jms:jme), &
         INTENT(INOUT) ::   smois


   END SUBROUTINE soa_vbs_addemiss

!------------------------------------------------------------------------
   SUBROUTINE soa_vbs_seasalt_emiss( &
      dtstep, u10, v10, alt, dz8w, xland, chem, &
      ids, ide, jds, jde, kds, kde, &
      ims, ime, jms, jme, kms, kme, &
      its, ite, jts, jte, kts, kte)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      REAL, INTENT(IN) ::    dtstep

! 10-m wind speed components (m/s)
      REAL, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN) ::   u10, v10, xland

! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(INOUT) ::   chem

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) ::   alt, dz8w

   END SUBROUTINE soa_vbs_seasalt_emiss
!----------------------------------------------------------------------

   subroutine soa_vbs_dust_emiss(slai, ust, smois, ivgtyp, isltyp, &
                                 id, dtstep, u10, v10, alt, dz8w, xland, num_soil_layers, &
                                 chem, &
                                 ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte)
!
! adds dust emissions for mosaic aerosol species (i.e. emission tendencies
! over time dtstep are applied to the aerosol mixing ratios)
!
! This is a simple dust scheme based on Shaw et al. (2008) to appear in
! Atmospheric Environment, recoded by Jerome Fast
!
! NOTE:
! 1) This version only works with the 8-bin version of MOSAIC.
! 2) Dust added to MOSAIC's other inorganic specie, OIN.  If Ca and CO3 are
!    activated in the Registry, a small fraction also added to Ca and CO3.
! 3) The main departure from Shaw et al., is now alphamask is computed since
!    the land-use categories in that paper and in WRF differ.  WRF currently
!    does not have that many land-use categories and adhoc assumptions had to
!    be made. This version was tested for Mexico in the dry season.  The main
!    land-use categories in WRF that are likely dust sources are grass, shrub,
!    and savannna (that WRF has in the desert regions of NW Mexico).  Having
!    dust emitted from these types for other locations and other times of the
!    year is not likely to be valid.
! 4) An upper bound on ustar was placed because the surface parameterizations
!    in WRF can produce unrealistically high values that lead to very high
!    dust emission rates.
! 5) Other departures' from Shaw et al. noted below, but are probably not as
!    important as 2) and 3).
!
      USE module_configure, only:grid_config_rec_type
      USE module_state_description, only:num_chem, param_first_scalar
      USE module_data_mosaic_asect

      IMPLICIT NONE

!  TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

      INTEGER, INTENT(IN) :: id, num_soil_layers, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte

      REAL, INTENT(IN) ::    dtstep

! 10-m wind speed components (m/s)
      REAL, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN) ::   u10, v10, xland, slai, ust
      INTEGER, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN) ::   ivgtyp, isltyp

! trace species mixing ratios (aerosol mass = ug/kg; number = #/kg)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(INOUT) ::   chem

! alt  = 1.0/(dry air density) in (m3/kg)
! dz8w = layer thickness in (m)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) ::   alt, dz8w

      REAL, DIMENSION(ims:ime, num_soil_layers, jms:jme), &
         INTENT(INOUT) ::   smois

! local variables
      integer i, j, k, l, l_oin, l_ca, l_co3, n, ii
      integer iphase, itype, izob
      integer p1st

      real dum, dumdlo, dumdhi, dumlandfrac, dumspd10
      real factaa, factbb, fracoin, fracca, fracco3, fractot
      real ustart, ustar1, ustart0
      real alphamask, f8, f50, f51, f52, wetfactor, sumdelta, ftot
      real smois_grav, wp, pclay
      real :: beta(4, 7)
      real :: gamma(4), delta(4)
      real :: sz(8)
      real :: dustflux, densdust, mass1part
      real :: dp_meanvol_tmp

      return

   END subroutine soa_vbs_dust_emiss

   subroutine soa_vbs_dust_gocartemis(ktau, dt, num_soil_layers, alt, u_phy, &
                                      v_phy, chem, rho_phy, dz8w, smois, u10, v10, p8w, erod, &
                                      ivgtyp, isltyp, xland, dx, g, &
                                      ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte)
      USE module_configure
      USE module_state_description
      USE module_model_constants, ONLY:mwdry
      USE module_data_mosaic_asect
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ktau, num_soil_layers, &
                             ids, ide, jds, jde, kds, kde, &
                             ims, ime, jms, jme, kms, kme, &
                             its, ite, jts, jte, kts, kte
      INTEGER, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN) :: &
         ivgtyp, &
         isltyp
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme, num_chem), &
         INTENT(INOUT) ::                                   chem
      REAL, DIMENSION(ims:ime, num_soil_layers, jms:jme), &
         INTENT(INOUT) ::                               smois
      REAL, DIMENSION(ims:ime, jms:jme, 3), &
         INTENT(IN) ::    erod
      REAL, DIMENSION(ims:ime, jms:jme), &
         INTENT(IN) :: &
         u10, &
         v10, &
         xland
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), &
         INTENT(IN) :: &
         alt, &
         dz8w, p8w, &
         u_phy, v_phy, rho_phy

      REAL, INTENT(IN) :: dt, dx, g

   end subroutine soa_vbs_dust_gocartemis

   SUBROUTINE soa_vbs_source_du(nmx, dt1, i, j, &
                                erod, ilwi, dxy, w10m, gwet, airden, airmas, &
                                bems, month, g0)

      INTEGER, INTENT(IN)    :: nmx
      REAL*8, INTENT(IN)  :: erod(ndcls)
      INTEGER, INTENT(IN)    :: ilwi, month

      REAL*8, INTENT(IN)    :: w10m, gwet
      REAL*8, INTENT(IN)    :: dxy
      REAL*8, INTENT(IN)    :: airden, airmas
      REAL*8, INTENT(OUT)   :: bems(nmx)

      REAL*8    :: den(nmx), diam(nmx)
      REAL*8    :: tsrc, u_ts0, cw, u_ts, dsrc, srce
      REAL, intent(in)    :: g0
      REAL    :: rhoa, g, dt1
      INTEGER :: i, j, n, m, k

   END SUBROUTINE soa_vbs_source_du

!===========================================================================

!===========================================================================

END Module module_aerosols_soa_vbs
