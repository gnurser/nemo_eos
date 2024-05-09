MODULE eosbn2
   !!==============================================================================
   !!                       ***  MODULE  eosbn2  ***
   !! Equation Of Seawater : in situ density - Brunt-Vaisala frequency
   !!==============================================================================
   !! History :  OPA  ! 1989-03  (O. Marti)  Original code
   !!            6.0  ! 1994-07  (G. Madec, M. Imbard)  add bn2
   !!            6.0  ! 1994-08  (G. Madec)  Add Jackett & McDougall eos
   !!            7.0  ! 1996-01  (G. Madec)  statement function for e3
   !!            8.1  ! 1997-07  (G. Madec)  density instead of volumic mass
   !!             -   ! 1999-02  (G. Madec, N. Grima) semi-implicit pressure gradient
   !!            8.2  ! 2001-09  (M. Ben Jelloul)  bugfix on linear eos
   !!   NEMO     1.0  ! 2002-10  (G. Madec)  add eos_init
   !!             -   ! 2002-11  (G. Madec, A. Bozec)  partial step, eos_insitu_2d
   !!             -   ! 2003-08  (G. Madec)  F90, free form
   !!            3.0  ! 2006-08  (G. Madec)  add tfreez function (now eos_fzp function)
   !!            3.3  ! 2010-05  (C. Ethe, G. Madec)  merge TRC-TRA
   !!             -   ! 2010-10  (G. Nurser, G. Madec)  add alpha/beta used in ldfslp
   !!            3.7  ! 2012-03  (F. Roquet, G. Madec)  add primitive of alpha and beta used in PE computation
   !!             -   ! 2012-05  (F. Roquet)  add Vallis and original JM95 equation of state
   !!             -   ! 2013-04  (F. Roquet, G. Madec)  add eos_rab, change bn2 computation and reorganize the module
   !!             -   ! 2014-09  (F. Roquet)  add TEOS-10, S-EOS, and modify EOS-80
   !!             -   ! 2015-06  (P.A. Bouttier) eos_fzp functions changed to subroutines for AGRIF
   !!           4.2.x ! 2024-04  (G. Nurser) modified to be wrapped by f2py
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   eos           : generic interface of the equation of state
   !!   eos_insitu    : Compute the in situ density
   !!   eos_insitu_pot: Compute the insitu and surface referenced potential volumic mass
   !!   eos_insitu_2d : Compute the in situ density for 2d fields
   !!   bn2           : compute the Brunt-Vaisala frequency
   !!   eos_pt_from_ct: compute the potential temperature from the Conservative Temperature
   !!   eos_rab       : generic interface of in situ thermal/haline expansion ratio
   !!   eos_rab_3d    : compute in situ thermal/haline expansion ratio
   !!   eos_rab_2d    : compute in situ thermal/haline expansion ratio for 2d fields
   !!   eos_fzp_2d    : freezing temperature for 2d fields
   !!   eos_fzp_0d    : freezing temperature for scalar
   !!   eos_init      : set eos parameters (namelist)
   !!----------------------------------------------------------------------

   USE OMP_LIB, ONLY : omp_set_num_threads, omp_get_thread_num, omp_get_max_threads
   IMPLICIT NONE

   INTEGER     ::   neos            ! Identifier for equation of state used

   INTEGER, PRIVATE , PARAMETER ::   np_teos10    = -1 ! parameter for using TEOS10
   INTEGER, PRIVATE , PARAMETER ::   np_eos80     =  0 ! parameter for using EOS80
   INTEGER, PRIVATE , PARAMETER ::   np_old_eos80 =  2 ! parameter for using Macdougall and Jackett EOS80
   INTEGER, PRIVATE , PARAMETER ::   np_seos      =  1 ! parameter for using Simplified Equation of state

   REAL*8 ::  rho0        = 1026.d0          !: volumic mass of reference     [kg/m3]
   REAL*8, PRIVATE ::  r1_rho0                    ! reciprocal of volumic mass of reference     [kg/m3]
   ! REAL*8, PRIVATE ::  grav     = 9.80665d0       !: gravity                            [m/s2]

   !                               !!!  simplified eos coefficients (default value: Vallis 2006)
   REAL*8 ::   rn_a0      = 1.6550d-1     ! thermal expansion coeff.
   REAL*8 ::   rn_b0      = 7.6554d-1     ! saline  expansion coeff.
   REAL*8 ::   rn_lambda1 = 5.9520d-2     ! cabbeling coeff. in T^2
   REAL*8 ::   rn_lambda2 = 5.4914d-4     ! cabbeling coeff. in S^2
   REAL*8 ::   rn_mu1     = 1.4970d-4     ! thermobaric coeff. in T
   REAL*8 ::   rn_mu2     = 1.1090d-5     ! thermobaric coeff. in S
   REAL*8 ::   rn_nu      = 2.4341d-3     ! cabbeling coeff. in theta*salt

   ! TEOS10/EOS80 parameters
   REAL*8, PRIVATE ::   r1_S0, r1_T0, r1_Z0, rdeltaS

   REAL*8, PRIVATE, PARAMETER ::     R00 = 4.6494977072d+01
   REAL*8, PRIVATE, PARAMETER ::     R01 = -5.2099962525d+00
   REAL*8, PRIVATE, PARAMETER ::     R02 = 2.2601900708d-01
   REAL*8, PRIVATE, PARAMETER ::     R03 = 6.4326772569d-02
   REAL*8, PRIVATE, PARAMETER ::     R04 = 1.5616995503d-02
   REAL*8, PRIVATE, PARAMETER ::     R05 = -1.7243708991d-03

   ! EOS parameters
   REAL*8, PRIVATE ::   EOS000 , EOS100 , EOS200 , EOS300 , EOS400 , EOS500 , EOS600
   REAL*8, PRIVATE ::   EOS010 , EOS110 , EOS210 , EOS310 , EOS410 , EOS510
   REAL*8, PRIVATE ::   EOS020 , EOS120 , EOS220 , EOS320 , EOS420
   REAL*8, PRIVATE ::   EOS030 , EOS130 , EOS230 , EOS330
   REAL*8, PRIVATE ::   EOS040 , EOS140 , EOS240
   REAL*8, PRIVATE ::   EOS050 , EOS150
   REAL*8, PRIVATE ::   EOS060
   REAL*8, PRIVATE ::   EOS001 , EOS101 , EOS201 , EOS301 , EOS401
   REAL*8, PRIVATE ::   EOS011 , EOS111 , EOS211 , EOS311
   REAL*8, PRIVATE ::   EOS021 , EOS121 , EOS221
   REAL*8, PRIVATE ::   EOS031 , EOS131
   REAL*8, PRIVATE ::   EOS041
   REAL*8, PRIVATE ::   EOS002 , EOS102 , EOS202
   REAL*8, PRIVATE ::   EOS012 , EOS112
   REAL*8, PRIVATE ::   EOS022
   REAL*8, PRIVATE ::   EOS003 , EOS103
   REAL*8, PRIVATE ::   EOS013

   ! ALPHA parameters
   REAL*8, PRIVATE ::   ALP000 , ALP100 , ALP200 , ALP300 , ALP400 , ALP500
   REAL*8, PRIVATE ::   ALP010 , ALP110 , ALP210 , ALP310 , ALP410
   REAL*8, PRIVATE ::   ALP020 , ALP120 , ALP220 , ALP320
   REAL*8, PRIVATE ::   ALP030 , ALP130 , ALP230
   REAL*8, PRIVATE ::   ALP040 , ALP140
   REAL*8, PRIVATE ::   ALP050
   REAL*8, PRIVATE ::   ALP001 , ALP101 , ALP201 , ALP301
   REAL*8, PRIVATE ::   ALP011 , ALP111 , ALP211
   REAL*8, PRIVATE ::   ALP021 , ALP121
   REAL*8, PRIVATE ::   ALP031
   REAL*8, PRIVATE ::   ALP002 , ALP102
   REAL*8, PRIVATE ::   ALP012
   REAL*8, PRIVATE ::   ALP003

   ! BETA parameters
   REAL*8, PRIVATE ::   BET000 , BET100 , BET200 , BET300 , BET400 , BET500
   REAL*8, PRIVATE ::   BET010 , BET110 , BET210 , BET310 , BET410
   REAL*8, PRIVATE ::   BET020 , BET120 , BET220 , BET320
   REAL*8, PRIVATE ::   BET030 , BET130 , BET230
   REAL*8, PRIVATE ::   BET040 , BET140
   REAL*8, PRIVATE ::   BET050
   REAL*8, PRIVATE ::   BET001 , BET101 , BET201 , BET301
   REAL*8, PRIVATE ::   BET011 , BET111 , BET211
   REAL*8, PRIVATE ::   BET021 , BET121
   REAL*8, PRIVATE ::   BET031
   REAL*8, PRIVATE ::   BET002 , BET102
   REAL*8, PRIVATE ::   BET012
   REAL*8, PRIVATE ::   BET003

   ! PEN parameters
   REAL*8, PRIVATE ::   PEN000 , PEN100 , PEN200 , PEN300 , PEN400
   REAL*8, PRIVATE ::   PEN010 , PEN110 , PEN210 , PEN310
   REAL*8, PRIVATE ::   PEN020 , PEN120 , PEN220
   REAL*8, PRIVATE ::   PEN030 , PEN130
   REAL*8, PRIVATE ::   PEN040
   REAL*8, PRIVATE ::   PEN001 , PEN101 , PEN201
   REAL*8, PRIVATE ::   PEN011 , PEN111
   REAL*8, PRIVATE ::   PEN021
   REAL*8, PRIVATE ::   PEN002 , PEN102
   REAL*8, PRIVATE ::   PEN012

   ! ALPHA_PEN parameters
   REAL*8, PRIVATE ::   APE000 , APE100 , APE200 , APE300
   REAL*8, PRIVATE ::   APE010 , APE110 , APE210
   REAL*8, PRIVATE ::   APE020 , APE120
   REAL*8, PRIVATE ::   APE030
   REAL*8, PRIVATE ::   APE001 , APE101
   REAL*8, PRIVATE ::   APE011
   REAL*8, PRIVATE ::   APE002

   ! BETA_PEN parameters
   REAL*8, PRIVATE ::   BPE000 , BPE100 , BPE200 , BPE300
   REAL*8, PRIVATE ::   BPE010 , BPE110 , BPE210
   REAL*8, PRIVATE ::   BPE020 , BPE120
   REAL*8, PRIVATE ::   BPE030
   REAL*8, PRIVATE ::   BPE001 , BPE101
   REAL*8, PRIVATE ::   BPE011
   REAL*8, PRIVATE ::   BPE002

CONTAINS

     SUBROUTINE set_eos_threads(nthreads)
       INTEGER*4, INTENT(IN)  :: nthreads
       CALL omp_set_num_threads(nthreads)
     END SUBROUTINE set_eos_threads

     SUBROUTINE get_eos_threads(nthreads)
       INTEGER*4, INTENT(OUT)  :: nthreads
       nthreads = omp_get_max_threads()
     END SUBROUTINE get_eos_threads

     SUBROUTINE set_eos(neos_in)
       INTEGER*4, INTENT(IN)  :: neos_in
       neos = neos_in
     END SUBROUTINE set_eos

     SUBROUTINE get_r0( depth_km, r0)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sigma_n4  ***
      !!
      !! ** Purpose :   Get vertical reference profile r0
      !!
      !!    where [Eq. (13) of Roquet et al. (2015)] the full potential density
      !!   rho(S_A,CT,z) = r0(z) + r(S_A,CT,z)
      !!  with
      !!     r0(z) = rho(35.16504g/kg,4 deg,z) - rho(35.16504g/kg,4 deg,0 m)
      !!       as defined in A.1 of Roquet et al. (2015)
      !! Check value for Z = −1000 m: r0 = 4.59763035 kg /m^3
      !!
      !!----------------------------------------------------------------------
      REAL*8, INTENT(IN) ::   depth_km
      REAL*8, INTENT(OUT) ::  r0
      !f2py intent (in) depth_kmx
      !f2py intent (out) r0
      !
      REAL*8 ::   zh              ! local scalars
      !!----------------------------------------------------------------------
      !
      zh = depth_km * 1.d3 * r1_Z0
      r0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
    END SUBROUTINE get_r0

    SUBROUTINE eos_insitu4(T, S, depth, rho, n)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE eos_insitu  ***
      !!
      !! ** Purpose :   Compute the in unmasked situ (density deviation from 1000kg/m^3
      !!  rho(t,s,z)) from
      !!       potential temperature salinity and depth using an equation of state
      !!       selected in the nameos namelist
      !!
      !! ** Method  :   rho(t,s,z) = r0(z) + r(z)=[rho(t,s,z) - r0(z)] - 1000.
      !!                t      TEOS10: CT or EOS80: PT      Celsius
      !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
      !!                z      depth                        meters
      !!                rho    in situ density anomaly      kg/m^3
      !!
      !!     np_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
      !!               Note global mean r0(z)
      !!         Check value: rho - r0 = 28.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
      !!
      !!     np_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
      !!         Check value: rho-r? = 28.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
      !!
      !!     np_seos : simplified equation of state
      !!              rho(t,s,z) = rho0 -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0) - 1000.
      !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
      !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
      !!              Vallis like equation: use default values of coefficients
      !!
      !!   np_old_eos80: original Jackett and Macdougall (1994) EOS80
      !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
      !!          t = 40 deg celcius, s=40 psu
      !!
      !! ** Action  :   compute rho , the in situ density anomaly (kg/m^3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::  n
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(in)             ::  depth   ! depth                                 [m]
      REAL*4, DIMENSION(n), INTENT(out)            ::  rho     ! deviation of rho from 1000          [kg/m^3]
      !
      INTEGER  ::  i
      REAL*8  :: zt, zh ! local scalars
      REAL*8  :: zs! local scalars
      REAL*8  :: zn1, zn2!   -      -
      REAL*8  :: zn, zn0, zn3!   -      -
      REAL*8  :: zr0
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
      REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0        !    -         -
      !!----------------------------------------------------------------------
      !
      SELECT CASE( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         !$omp parallel do private(zt,zs,zh,zn,zn0,zn1,zn2,zn3,zr0)
         DO i=1,n
            !
            zh  = depth(i) * r1_Z0                                  ! depth
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            zn3 = EOS013*zt   &
               &   + EOS103*zs+EOS003
               !
            zn2 = (EOS022*zt   &
               &   + EOS112*zs+EOS012)*zt   &
               &   + (EOS202*zs+EOS102)*zs+EOS002
               !
            zn1 = (((EOS041*zt   &
               &   + EOS131*zs+EOS031)*zt   &
               &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
               &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
               &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
               !
            zn0 = (((((EOS060*zt   &
               &   + EOS150*zs+EOS050)*zt   &
               &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
               &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
               &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
               &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
               &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
    !!----------------------------------------------------------------------
            !!  Subtract 1000. to give density anomaly
            IF (neos == np_teos10) THEN
            !!  Add reference profile zr0 to anomaly
               zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
               rho(i) = REAL(zn - 1000.d0 + zr0, KIND=4)
            ELSE
               rho(i) = REAL(zn - 1000.d0, KIND=4)
            END IF
    !!----------------------------------------------------------------------
            !
            ! prd(i) = (  zn * r1_rho0 - 1.d0  ) * ztm  ! density anomaly (masked)
            !
         END DO
         !$omp  end parallel do
         !
      CASE( np_seos )                !==  simplified EOS  ==!
         !
         DO i=1,n
            zt  = T(i) - 10.d0
            zs  = S(i) - 35.d0
            zh  = depth (i)
            !
            zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
               &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
               &  - rn_nu * zt * zs
               !
            ! prd(i) = zn * r1_rho0 * ztm                ! density anomaly (masked)
            rho(i) = REAL(zn + rho0 - 1000.d0, KIND=4)
         END DO
         !
      CASE(np_old_eos80)
         !$omp parallel do private(zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze, zbw,zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
         DO i=1,n
            zt = T(i)
            zs = S(i)
            zh = depth(i)                 ! depth
            zsr= SQRT( ABS( zs) )        ! square root salinity
            !
            ! compute volumic mass pure water at atm pressure
            zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4 )*zt   &
                 &                          -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594d0
            ! seawater volumic mass atm pressure
            zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
                 &                                         -4.0899e-3 ) *zt+0.824493d0
            zr3= ( -1.6546e-6*zt+1.0227e-4 )    *zt-5.72466e-3
            zr4= 4.8314e-4
            !
            ! potential volumic mass (reference to the surface)
            zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
            !
            ! add the compression terms
            ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
            zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
            zb = zbw + ze * zs
            !
            zd = -2.042967e-2
            zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
            zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859d0 ) *zt - 4.721788d0
            za = ( zd*zsr + zc ) *zs + zaw
            !
            zb1=   (  -0.1909078d0  *zt+7.390729d0    ) *zt-55.87545d0
            za1= ( (   2.326469e-3*zt+1.553190d0    ) *zt-65.00517d0 ) *zt + 1044.077d0
            zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638d0 ) *zt + 2098.925d0 ) *zt+190925.6d0
            zk0= ( zb1*zsr + za1 )*zs + zkw
            !
            ! masked in situ density anomaly
            rho(i) = REAL(zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0, KIND=4)
         END DO
         !$omp  end parallel do
      END SELECT
      !
    END SUBROUTINE eos_insitu4


    SUBROUTINE eos_insitu4_m(fillvalue, mask, T, S, depth, rho, n)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE eos_insitu  ***
      !!
      !! ** Purpose :   Compute the masked in situ density deviation from 1000kg/m^3 ( rho(t,s,z)) from
      !!       potential temperature salinity and depth using an equation of state
      !!       selected in the nameos namelist
      !!
      !! ** Method  :   rho(t,s,z)
      !!                t      TEOS10: CT or EOS80: PT      Celsius
      !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
      !!                z      depth                        meters
      !!                rho    in situ density              kg/m^3
      !!
      !!     ln_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
      !!               Note global mean r0(z) is added to perturbation r(t,s,z) = rho(t,s,z) - r0(z)
      !!         Check value: rho = 1028.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
      !!
      !!     ln_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
      !!         Check value: rho = 1028.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
      !!
      !!     ln_seos : simplified equation of state
      !!              rho(t,s,z) = -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0)
      !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
      !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
      !!              Vallis like equation: use default values of coefficients
      !!
      !!   np_old_eos80: original Macdougall and Jackett EOS80
      !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
      !!          t = 40 deg celcius, s=40 psu
      !!
      !! ** Action  :   compute rho , the in situ density (kg/m^3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      LOGICAL(KIND=1), DIMENSION(n), INTENT(IN)    ::  mask   !logical mask of land points
      REAL*4, INTENT(IN)                           :: fillvalue
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(in)             ::  depth   ! depth                                 [m]
      REAL*4, DIMENSION(n), INTENT(out)            ::  rho     ! deviation of rho from 1000          [kg/m^3]
      !
      INTEGER  ::  i
      REAL*8  :: zt, zh ! local scalars
      REAL*8  :: zs! local scalars
      REAL*8  :: zn1, zn2!   -      -
      REAL*8  :: zn, zn0, zn3!   -      -
      REAL*8  :: zr0
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
      REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0        !    -         -
      !!----------------------------------------------------------------------
      !
      SELECT CASE( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         !$omp parallel do private(zt,zs,zh,zn,zn0,zn1,zn2,zn3,zr0)
         DO i=1,n
            IF (mask(i)) THEN
               rho(i) = fillvalue
               CYCLE
            END IF
            !
            zh  = depth(i) * r1_Z0                                  ! depth
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            zn3 = EOS013*zt   &
               &   + EOS103*zs+EOS003
               !
            zn2 = (EOS022*zt   &
               &   + EOS112*zs+EOS012)*zt   &
               &   + (EOS202*zs+EOS102)*zs+EOS002
               !
            zn1 = (((EOS041*zt   &
               &   + EOS131*zs+EOS031)*zt   &
               &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
               &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
               &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
               !
            zn0 = (((((EOS060*zt   &
               &   + EOS150*zs+EOS050)*zt   &
               &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
               &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
               &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
               &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
               &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
    !!----------------------------------------------------------------------
            !!  Subtract 1000. to give density anomaly
            IF (neos == np_teos10) THEN
            !!  Add reference profile zr0 to anomaly
               zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
               rho(i) = REAL(zn - 1000.d0 + zr0, KIND=4)
            ELSE
               rho(i) = REAL(zn - 1000.d0, KIND=4)
            END IF
    !!----------------------------------------------------------------------
            !
            ! prd(i) = (  zn * r1_rho0 - 1.d0  ) * ztm  ! density anomaly (masked)
            !
         END DO
         !$omp  end parallel do
         !
      CASE( np_seos )                !==  simplified EOS  ==!
         !
         DO i=1,n
            IF (mask(i)) THEN
               rho(i) = fillvalue
               CYCLE
            END IF
            zt  = T(i) - 10.d0
            zs  = S(i) - 35.d0
            zh  = depth (i)
            !
            zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
               &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
               &  - rn_nu * zt * zs
               !
            ! prd(i) = zn * r1_rho0 * ztm                ! density anomaly (masked)
            rho(i) = REAL(zn + rho0 - 1000.d0, KIND=4)
         END DO
         !
      CASE(np_old_eos80)
         !$omp parallel do private(zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze, zbw,zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
         DO i=1,n
            IF (mask(i)) THEN
               rho(i) = fillvalue
               CYCLE
            END IF

            zt = T(i)
            zs = S(i)
            zh = depth(i)                 ! depth
            zsr= SQRT( ABS( zs) )        ! square root salinity
            !
            ! compute volumic mass pure water at atm pressure
            zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4 )*zt   &
                 &                          -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594d0
            ! seawater volumic mass atm pressure
            zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
                 &                                         -4.0899e-3 ) *zt+0.824493d0
            zr3= ( -1.6546e-6*zt+1.0227e-4 )    *zt-5.72466e-3
            zr4= 4.8314e-4
            !
            ! potential volumic mass (reference to the surface)
            zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
            !
            ! add the compression terms
            ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
            zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
            zb = zbw + ze * zs
            !
            zd = -2.042967e-2
            zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
            zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859d0 ) *zt - 4.721788d0
            za = ( zd*zsr + zc ) *zs + zaw
            !
            zb1=   (  -0.1909078d0  *zt+7.390729d0    ) *zt-55.87545d0
            za1= ( (   2.326469e-3*zt+1.553190d0    ) *zt-65.00517d0 ) *zt + 1044.077d0
            zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638d0 ) *zt + 2098.925d0 ) *zt+190925.6d0
            zk0= ( zb1*zsr + za1 )*zs + zkw
            !
            ! masked in situ density anomaly
            rho(i) = REAL(zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0, KIND=4)
         END DO
         !$omp  end parallel do
       END SELECT
      !
    END SUBROUTINE eos_insitu4_m

    SUBROUTINE eos_sigman4(T, S, depth_km, rho, n)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE eos_insitu  ***
      !!
      !! ** Purpose :   Compute the in unmasked situ density deviation from 1000kg/m^3 ( rho(t,s,z)) from
      !!       potential temperature salinity and depth using an equation of state
      !!       selected in the nameos namelist
      !!
      !! ** Method  :   rho(t,s,z)
      !!                t      TEOS10: CT or EOS80: PT      Celsius
      !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
      !!                depth  reference depth                    km
      !!                rho    potential density              kg/m^3
      !!
      !!     ln_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
      !!               Note global mean r0(z) is added to perturbation r(t,s,z) = rho(t,s,z) - r0(z)
      !!         Check value: rho = 1028.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
      !!
      !!     ln_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
      !!         Check value: rho = 1028.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
      !!
      !!     ln_seos : simplified equation of state
      !!              rho(t,s,z) = -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0)
      !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
      !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
      !!              Vallis like equation: use default values of coefficients
      !!
      !! ** Action  :   compute rho , the in situ density (kg/m^3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::  n
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, INTENT(in)                           ::  depth_km! reference depth                    [km]
      REAL*4, DIMENSION(n), INTENT(out)            ::  rho     ! deviation of rho from 1000          [kg/m^3]
      !
      INTEGER  ::  i
      REAL*8  :: zt, zh ! local scalars
      REAL*8  :: zs! local scalars
      REAL*8  :: zn1, zn2!   -      -
      REAL*8  :: zn, zn0, zn3!   -      -
      REAL*8  :: zr0
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
      REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0        !    -         -
      !!----------------------------------------------------------------------
      !
      SELECT CASE( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         zh = depth_km * 1.d3 * r1_Z0
         IF (neos == np_teos10) THEN
            ! Define reference profile zr0 to be added to anomaly
             zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
         END IF
         DO i=1,n
            !
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            zn3 = EOS013*zt   &
               &   + EOS103*zs+EOS003
               !
            zn2 = (EOS022*zt   &
               &   + EOS112*zs+EOS012)*zt   &
               &   + (EOS202*zs+EOS102)*zs+EOS002
               !
            zn1 = (((EOS041*zt   &
               &   + EOS131*zs+EOS031)*zt   &
               &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
               &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
               &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
               !
            zn0 = (((((EOS060*zt   &
               &   + EOS150*zs+EOS050)*zt   &
               &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
               &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
               &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
               &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
               &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
    !!----------------------------------------------------------------------
            !!  Subtract 1000. to give density anomaly
            IF (neos == np_teos10) THEN
            !!  Add reference profile zr0 to anomaly
               rho(i) = REAL(zn - 1000.d0 + zr0, KIND=4)
            ELSE
               rho(i) = REAL(zn - 1000.d0, KIND=4)
            END IF
    !!----------------------------------------------------------------------
            !
            ! prd(i) = (  zn * r1_rho0 - 1.d0  ) * ztm  ! density anomaly (masked)
            !
         END DO
         !
      CASE( np_seos )                !==  simplified EOS  ==!
         !
         zh = depth_km * 1.d3
         DO i=1,n
            zt  = T(i) - 10.d0
            zs  = S(i) - 35.d0
            !
            zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
               &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
               &  - rn_nu * zt * zs
               !
            ! prd(i) = zn * r1_rho0 * ztm                ! density anomaly (masked)
            rho(i) = REAL(zn + rho0 - 1000.d0, KIND=4)
         END DO
         !
      CASE(np_old_eos80)
          zh = depth_km*1000.d0                  ! depth
          DO i=1,n
             zt = T(i)
             zs = S(i)
             zsr= SQRT( ABS( zs) )        ! square root salinity
             !
             ! compute volumic mass pure water at atm pressure
             zr1= ( ( ( ( 6.536332d-9*zt-1.120083d-6 )*zt+1.001685d-4 )*zt   &
                  &                          -9.095290d-3 )*zt+6.793952d-2 )*zt+999.842594d0
             ! seawater volumic mass atm pressure
             zr2= ( ( ( 5.3875d-9*zt-8.2467d-7 ) *zt+7.6438d-5 ) *zt   &
                  &                                         -4.0899d-3 ) *zt+0.824493d0
             zr3= ( -1.6546d-6*zt+1.0227d-4 )    *zt-5.72466d-3
             zr4= 4.8314d-4
             !
             ! potential volumic mass (reference to the surface)
             zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
             !
             ! add the compression terms
             ze = ( -3.508914d-8*zt-1.248266d-8 ) *zt-2.595994d-6
             zbw= (  1.296821d-6*zt-5.782165d-9 ) *zt+1.045941d-4
             zb = zbw + ze * zs
             !
             zd = -2.042967d-2
             zc =   (-7.267926d-5*zt+2.598241d-3 ) *zt+0.1571896
             zaw= ( ( 5.939910d-6*zt+2.512549d-3 ) *zt-0.1028859d0 ) *zt - 4.721788d0
             za = ( zd*zsr + zc ) *zs + zaw
             !
             zb1=   (  -0.1909078d0  *zt+7.390729d0    ) *zt-55.87545d0
             za1= ( (   2.326469d-3*zt+1.553190d0    ) *zt-65.00517d0 ) *zt + 1044.077d0
             zkw= ( ( (-1.361629d-4*zt-1.852732d-2 ) *zt-30.41638d0 ) *zt + 2098.925d0 ) *zt+190925.6d0
             zk0= ( zb1*zsr + za1 )*zs + zkw
             !
             ! unmasked in situ density anomaly
             rho(i) = REAL(zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0, KIND=4)
          end do
      END SELECT
      !
   END SUBROUTINE eos_sigman4



   SUBROUTINE eos_sigman4_m(fillvalue, mask, T, S, depth_km, rho, n)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE eos_insitu  ***
      !!
      !! ** Purpose :   Compute the masked situ density deviation refrerenced to depth depth_km*1000. m
      !!   from 1000kg/m^3 ( rho(t,s,z)) from
      !!       potential temperature salinity and depth using an equation of state
      !!       selected in the nameos namelist
      !!
      !! ** Method  :   rho(t,s,z)
      !!                t      TEOS10: CT or EOS80: PT      Celsius
      !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
      !!                depth_km reference depth                    km
      !!                rho    potential density              kg/m^3
      !!
      !!     ln_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
      !!               Note global mean r0(z) is added to perturbation r(t,s,z) = rho(t,s,z) - r0(z)
      !!         Check value: rho = 1028.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
      !!
      !!     ln_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
      !!         Check value: rho = 1028.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
      !!
      !!     ln_seos : simplified equation of state
      !!              rho(t,s,z) = -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0)
      !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
      !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
      !!              Vallis like equation: use default values of coefficients
      !!
      !! ** Action  :   compute rho , the in situ density (kg/m^3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      LOGICAL(KIND=1), DIMENSION(n), INTENT(IN)    :: mask   !logical mask of land points
      REAL*4, INTENT(IN)                           ::   fillvalue
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, INTENT(in)                           ::  depth_km!                    [km]
      REAL*4, DIMENSION(n), INTENT(out)            ::  rho     ! deviation of rho from 1000          [kg/m^3]
      !
      INTEGER  ::  i
      REAL*8  :: zt, zh! local scalars
      REAL*8  :: zs! local scalars
      REAL*8  :: zn1, zn2!   -      -
      REAL*8  :: zn, zn0, zn3!   -      -
      REAL*8  :: zr0
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
      REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0        !    -         -
      !!----------------------------------------------------------------------
      !
      SELECT CASE( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         zh = depth_km * 1.d3 * r1_Z0
         IF (neos == np_teos10) THEN
            ! Define reference profile zr0 to be added to anomaly
              zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
         END IF
         DO i=1,n
            IF (mask(i)) THEN
               rho(i) = fillvalue
               CYCLE
            END IF
            !
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            zn3 = EOS013*zt   &
               &   + EOS103*zs+EOS003
               !
            zn2 = (EOS022*zt   &
               &   + EOS112*zs+EOS012)*zt   &
               &   + (EOS202*zs+EOS102)*zs+EOS002
               !
            zn1 = (((EOS041*zt   &
               &   + EOS131*zs+EOS031)*zt   &
               &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
               &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
               &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
               !
            zn0 = (((((EOS060*zt   &
               &   + EOS150*zs+EOS050)*zt   &
               &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
               &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
               &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
               &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
               &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
    !!----------------------------------------------------------------------
            !!  Subtract 1000. to give density anomaly
            IF (neos == np_teos10) THEN
            !!  Add reference profile zr0 to anomaly
               rho(i) = REAL(zn - 1000.d0 + zr0, KIND=4)
            ELSE
               rho(i) = REAL(zn - 1000.d0, KIND=4)
            END IF
    !!----------------------------------------------------------------------
            !
            ! prd(i) = (  zn * r1_rho0 - 1.d0  ) * ztm  ! density anomaly (masked)
            !
         END DO
         !
      CASE( np_seos )                !==  simplified EOS  ==!
         !
         zh = depth_km * 1.d3
         DO i=1,n
            IF (mask(i)) THEN
               rho(i) = fillvalue
               CYCLE
            END IF
            zt  = T(i) - 10.d0
            zs  = S(i) - 35.d0
            !
            zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
               &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
               &  - rn_nu * zt * zs
               !
            ! prd(i) = zn * r1_rho0 * ztm                ! density anomaly (masked)
            rho(i) = REAL(zn + rho0 - 1000.d0, KIND=4)
         END DO
         !
      CASE(np_old_eos80)
         zh = depth_km*1000.d0                  ! depth
         DO i=1,n
            IF (mask(i)) THEN
               rho(i) = fillvalue
               CYCLE
            END IF
            zt = T(i)
            zs = S(i)
            zsr= SQRT( ABS( zs) )        ! square root salinity
            !
            ! compute volumic mass pure water at atm pressure
            zr1= ( ( ( ( 6.536332d-9*zt-1.120083d-6 )*zt+1.001685d-4 )*zt   &
                 &                          -9.095290d-3 )*zt+6.793952d-2 )*zt+999.842594d0
            ! seawater volumic mass atm pressure
            zr2= ( ( ( 5.3875d-9*zt-8.2467d-7 ) *zt+7.6438d-5 ) *zt   &
                 &                                         -4.0899d-3 ) *zt+0.824493d0
            zr3= ( -1.6546d-6*zt+1.0227d-4 )    *zt-5.72466d-3
            zr4= 4.8314d-4
            !
            ! potential volumic mass (reference to the surface)
            zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
            !
            ! add the compression terms
            ze = ( -3.508914d-8*zt-1.248266d-8 ) *zt-2.595994d-6
            zbw= (  1.296821d-6*zt-5.782165d-9 ) *zt+1.045941d-4
            zb = zbw + ze * zs
            !
            zd = -2.042967d-2
            zc =   (-7.267926d-5*zt+2.598241d-3 ) *zt+0.1571896
            zaw= ( ( 5.939910d-6*zt+2.512549d-3 ) *zt-0.1028859d0 ) *zt - 4.721788d0
            za = ( zd*zsr + zc ) *zs + zaw
            !
            zb1=   (  -0.1909078d0  *zt+7.390729d0    ) *zt-55.87545d0
            za1= ( (   2.326469d-3*zt+1.553190d0    ) *zt-65.00517d0 ) *zt + 1044.077d0
            zkw= ( ( (-1.361629d-4*zt-1.852732d-2 ) *zt-30.41638d0 ) *zt + 2098.925d0 ) *zt+190925.6d0
            zk0= ( zb1*zsr + za1 )*zs + zkw
            !
            ! unmasked in situ density anomaly
            rho(i) = REAL(zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0, KIND=4)
          END DO
      END SELECT
      !
   END SUBROUTINE eos_sigman4_m

   SUBROUTINE eos_sigma04(T, S, sigma0, n)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE eos_insitu_pot  ***
      !!
      !! ** Purpose :   Compute unmasked sigma0, potential density referenced to the surface (Kg/m3)
      !!      from potential temperature and
      !!      salinity fields using an equation of state selected in the
      !!     namelist.
      !!
      !! ** Action  : - sigma0, the deviation of the potential volumic mass fro 1000 (Kg/m3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]

      REAL*4, DIMENSION(n), INTENT(out)            ::  sigma0     ! deviation of rho from 1000          [kg/m^3]
      !
      INTEGER  ::  i

      REAL*8  :: zt ! local scalars
      REAL*8  :: zs! local scalars
      REAL*8  :: zn, zn0!   -      -
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop  ! temporary scalars
      !!----------------------------------------------------------------------
      SELECT CASE ( neos )
         !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         DO i=1,n
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            zn0 = (((((EOS060*zt   &
                 &   + EOS150*zs+EOS050)*zt   &
                 &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
                 &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
                 &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
                 &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
                 &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
            !
            sigma0(i) = REAL(zn0 - 1000.d0, KIND=4)                          ! potential density referenced at the surface
         END DO

      CASE( np_seos )                !==  simplified EOS  ==!
         !
         DO i=1,n
            zt  = T(i) - 10.d0
            zs  = S(i) - 35.d0
            !                                                     ! potential density referenced at the surface
            zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt ) * zt   &
                 &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs ) * zs   &
                 &  - rn_nu * zt * zs
            sigma0(i) =  REAL(rho0 + zn - 1000.d0, KIND=4)
            !
         END DO
         !
      CASE(np_old_eos80)
         DO i=1,n
             zt = T(i)
             zs = S(i)
             zsr= SQRT( ABS( zs) )        ! square root salinity
             !
             ! compute volumic mass pure water at atm pressure
             zr1= ( ( ( ( 6.536332d-9*zt-1.120083d-6 )*zt+1.001685d-4 )*zt   &
                  &                          -9.095290d-3 )*zt+6.793952d-2 )*zt+999.842594d0
             ! seawater volumic mass atm pressure
             zr2= ( ( ( 5.3875d-9*zt-8.2467d-7 ) *zt+7.6438d-5 ) *zt   &
                  &                                         -4.0899d-3 ) *zt+0.824493d0
             zr3= ( -1.6546d-6*zt+1.0227d-4 )    *zt-5.72466d-3
             zr4= 4.8314d-4
             !
             ! potential volumic mass (reference to the surface)
             zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
             ! masked in situ density anomaly
             sigma0(i) = REAL(zrhop - 1000.d0, KIND=4)
          END DO
      END SELECT
      !
   END SUBROUTINE eos_sigma04


   SUBROUTINE eos_sigma04_m(fillvalue, mask, T, S, sigma0, n)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE eos_insitu_pot  ***
      !!
      !! ** Purpose :   Compute masked sigma0, potential density referenced to the surface (Kg/m3)
      !!      from potential temperature and
      !!      salinity fields using an equation of state selected in the
      !!     namelist.
      !!
      !! ** Action  : - sigma0, the deviation of the potential volumic mass fro 1000 (Kg/m3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      LOGICAL(KIND=1), DIMENSION(n), INTENT(IN)    :: mask   !logical mask of land points
      REAL*4, INTENT(IN)                           ::  fillvalue
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(out)            ::  sigma0     ! deviation of rho from 1000          [kg/m^3]
      !
      INTEGER  ::  i

      REAL*8  :: zt ! local scalars
      REAL*8  :: zs! local scalars
      REAL*8  :: zn, zn0!   -      -
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop  ! temporary scalars
        !!----------------------------------------------------------------------
      SELECT CASE ( neos )
         !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         DO i=1,n
            IF (mask(i)) THEN
               sigma0(i) = fillvalue
               CYCLE
            END IF
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            zn0 = (((((EOS060*zt   &
                 &   + EOS150*zs+EOS050)*zt   &
                 &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
                 &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
                 &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
                 &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
                 &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
            !
            sigma0(i) = REAL(zn0 - 1000.d0, KIND=4)                          ! potential density referenced at the surface
         END DO

      CASE( np_seos )                !==  simplified EOS  ==!
         !
         DO i=1,n
            IF (mask(i)) THEN
               sigma0(i) = fillvalue
               CYCLE
            END IF
            zt  = T(i) - 10.d0
            zs  = S(i) - 35.d0
            !                                                     ! potential density referenced at the surface
            zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt ) * zt   &
                 &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs ) * zs   &
                 &  - rn_nu * zt * zs
            sigma0(i) =  REAL(rho0 + zn - 1000.d0, KIND=4)
            !
         END DO
         !
      CASE(np_old_eos80)
         DO i=1,n
             if (mask(i)) then
                sigma0(i) = fillvalue
                cycle
             end if
             zt = T(i)
             zs = S(i)
             zsr= SQRT( ABS( zs) )        ! square root salinity
             !
             ! compute volumic mass pure water at atm pressure
             zr1= ( ( ( ( 6.536332d-9*zt-1.120083d-6 )*zt+1.001685d-4 )*zt   &
                  &                          -9.095290d-3 )*zt+6.793952d-2 )*zt+999.842594d0
             ! seawater volumic mass atm pressure
             zr2= ( ( ( 5.3875d-9*zt-8.2467d-7 ) *zt+7.6438d-5 ) *zt   &
                  &                                         -4.0899d-3 ) *zt+0.824493d0
             zr3= ( -1.6546d-6*zt+1.0227d-4 )    *zt-5.72466d-3
             zr4= 4.8314d-4
             !
             ! potential volumic mass (reference to the surface)
             zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
             ! masked in situ density anomaly
             sigma0(i) = REAL(zrhop - 1000.d0, KIND=4)
          END DO
      END SELECT
      !
   END SUBROUTINE eos_sigma04_m


   SUBROUTINE eos_rab4( T, S, depth, alpha, beta, n )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE rab_3d  ***
      !!
      !! ** Purpose :   Calculates thermal/haline expansion ratio at T-points
      !!
      !! ** Method  :   calculates (unmasked) alpha / beta at T-points
      !!
      !!                T      TEOS10: CT or EOS80: PT      [Celsius]
      !!                S      TEOS10: SA or EOS80: SP      TEOS10: [g/kg] or EOS80: [psu]
      !!                depth      depth                        [meters]
      !!                alpha     thermal expansion coeficient [K^{-1}]
      !!                beta     saline contraction coeficient TEOS-10 [kg/g] or EOS80 [psu^{-1}]
      !!
      !!     np_teos10 : polynomial TEOS-10
      !!          Check values for SA = 30 g/kg, CT = 10◦C, Z = −1000 m:
      !!          a = alpha*rho0 = 0.179646281 kg m−3 K−1 , b = beta*rho0 = 0.765555368 kg m−3 (g/kg)−1 .
      !!
      !!
      !!     np_eos80 : polynomial EOS-80 equation of state is used for alpha, beta.
      !!         Check value: ??
      !!
      !!     np_seos : simplified equation of state
      !!              rho(t,s,z) = rho0 -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0) - 1000.
      !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
      !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
      !!              Vallis like equation: use default values of coefficients
      !!
      !!   np_old_eos80: 
      !!         References :   McDougall, J. Phys. Oceanogr., 17, 1950-1964, 1987.
      !!         Check value: alpha/beta = 0.34763 psu/K at T=10 Celsius, psu = 40 psu, depth=4000m
      !!         Check value: beta = 0.72088 x 1.e-3 [psu^{-1}]   at T=10 Celsius, psu = 40 psu, depth=4000m
      !!
      !! ** Action  : - outputs alpha, beta     : thermal/haline expansion ratio at T-points
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
       !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(in)             ::  depth   ! depth                                 [m]
      REAL*4, DIMENSION(n), INTENT(out)            ::  alpha   ! thermal expansion coefficent          [Celsius^{-1}]
      REAL*4, DIMENSION(n), INTENT(out)            ::  beta    ! saline contraction coefficent          [psu^{-1}/kg/g]
      !
      INTEGER  ::  i
      !
      REAL*8 ::   zt , zh , zs         ! local scalars
      REAL*8 ::   zn , zn0, zn1, zn2, zn3   !   -      -
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zalbet  ! temporary scalar
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         do i=1,n
            !
            zh  = depth (i) * r1_Z0                                ! depth
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            ! alpha
            zn3 = ALP003
            !
            zn2 = ALP012*zt + ALP102*zs+ALP002
            !
            zn1 = ((ALP031*zt   &
               &   + ALP121*zs+ALP021)*zt   &
               &   + (ALP211*zs+ALP111)*zs+ALP011)*zt   &
               &   + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
               !
            zn0 = ((((ALP050*zt   &
               &   + ALP140*zs+ALP040)*zt   &
               &   + (ALP230*zs+ALP130)*zs+ALP030)*zt   &
               &   + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   &
               &   + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   &
               &   + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            !
            alpha(i) = REAL(zn * r1_rho0, KIND=4)
            !
            ! beta
            zn3 = BET003
            !
            zn2 = BET012*zt + BET102*zs+BET002
            !
            zn1 = ((BET031*zt   &
               &   + BET121*zs+BET021)*zt   &
               &   + (BET211*zs+BET111)*zs+BET011)*zt   &
               &   + ((BET301*zs+BET201)*zs+BET101)*zs+BET001
               !
            zn0 = ((((BET050*zt   &
               &   + BET140*zs+BET040)*zt   &
               &   + (BET230*zs+BET130)*zs+BET030)*zt   &
               &   + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   &
               &   + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   &
               &   + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            !
            beta(i) =  REAL(zn / zs * r1_rho0, KIND=4)
            !
         END DO
         !
      CASE( np_seos )                  !==  simplified EOS  ==!
         !
         do i=1,n
            zt  = T (i) - 10.d0   ! pot. temperature anomaly (t-T0)
            zs  = S (i) - 35.d0   ! abs. salinity anomaly (s-S0)
            zh  = depth(i)                ! depth in meters at t-point

            !
            zn  = rn_a0 * ( 1.d0 + rn_lambda1*zt + rn_mu1*zh ) + rn_nu*zs
            alpha(i) = REAL(zn * r1_rho0, KIND=4)   ! alpha
            !
            zn  = rn_b0 * ( 1.d0 - rn_lambda2*zs - rn_mu2*zh ) - rn_nu*zt
            beta(i) = REAL(zn * r1_rho0, KIND=4)   ! beta
            !
         END DO
         !
      CASE( np_old_eos80 )                  !==  simplified EOS  ==!
         DO i=1,n
            zt = T(i)
            zs = S(i) - 35.d0
            zh = depth(i)                 ! depth
            zalbet = ( ( ( - 0.255019d-07 * zt + 0.298357d-05 ) * zt   &   ! ratio alpha/bdta
                 &                                  - 0.203814d-03 ) * zt   &
                 &                                  + 0.170907d-01 ) * zt   &
                 &   +         0.665157d-01                                 &
                 &   +     ( - 0.678662d-05 * zs                            &
                 &           - 0.846960d-04 * zt + 0.378110d-02 ) * zs   &
                 &   +   ( ( - 0.302285d-13 * zh                            &
                 &           - 0.251520d-11 * zs                            &
                 &           + 0.512857d-12 * zt * zt              ) * zh   &
                 &           - 0.164759d-06 * zs                            &
                 &        +(   0.791325d-08 * zt - 0.933746d-06 ) * zt   &
                 &                                  + 0.380374d-04 ) * zh
            !
            beta(i) = REAL( &
                 &   ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
                 &                               - 0.301985d-05 ) * zt      &
                 &   +       0.785567d-03                                   &
                 &   + (     0.515032d-08 * zs                              &
                 &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
                 &   + ( (   0.121551d-17 * zh                              &
                 &         - 0.602281d-15 * zs                              &
                 &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
                 &                                + 0.408195d-10   * zs     &
                 &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
                 &                                - 0.121555d-07 ) * zh  &
                 &    , KIND=4)
            alpha(i) = REAL(zalbet*beta(i), KIND=4)
         END DO
      END SELECT
      !
   END SUBROUTINE eos_rab4


   SUBROUTINE eos_rab4_m( fillvalue, mask, T, S, depth, alpha, beta, n )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE rab_3d  ***
      !!
      !! ** Purpose :   Calculates thermal/haline expansion ratio at T-points
      !!
      !! ** Method  :   calculates (unmasked) alpha / beta at T-points
      !!
      !! ** Action  : - outputs ajpha, beta     : thermal/haline expansion ratio at T-points
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      LOGICAL(KIND=1), DIMENSION(n), INTENT(IN)    :: mask   !logical mask of land points
      REAL*4, INTENT(IN)                       ::  fillvalue
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(in)             ::  depth   ! depth                                 [m]
      REAL*4, DIMENSION(n), INTENT(out)            ::  alpha   ! thermal expansion coefficent          [Celsius^{-1}]
      REAL*4, DIMENSION(n), INTENT(out)            ::  beta    ! saline contraction coefficent          [psu^{-1}/kg/g]
      !
      INTEGER  ::  i
      !
      REAL*8 ::   zt , zh , zs         ! local scalars
      REAL*8 ::   zn , zn0, zn1, zn2, zn3   !   -      -
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zalbet  ! temporary scalar
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         DO i=1,n
            IF (mask(i)) THEN
               alpha(i) = fillvalue
               beta(i) = fillvalue
               CYCLE
            END IF
            !
            zh  = depth (i) * r1_Z0                                ! depth
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            ! alpha
            zn3 = ALP003
            !
            zn2 = ALP012*zt + ALP102*zs+ALP002
            !
            zn1 = ((ALP031*zt   &
               &   + ALP121*zs+ALP021)*zt   &
               &   + (ALP211*zs+ALP111)*zs+ALP011)*zt   &
               &   + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
               !
            zn0 = ((((ALP050*zt   &
               &   + ALP140*zs+ALP040)*zt   &
               &   + (ALP230*zs+ALP130)*zs+ALP030)*zt   &
               &   + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   &
               &   + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   &
               &   + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            !
            alpha(i) = REAL(zn * r1_rho0, KIND=4)
            !
            ! beta
            zn3 = BET003
            !
            zn2 = BET012*zt + BET102*zs+BET002
            !
            zn1 = ((BET031*zt   &
               &   + BET121*zs+BET021)*zt   &
               &   + (BET211*zs+BET111)*zs+BET011)*zt   &
               &   + ((BET301*zs+BET201)*zs+BET101)*zs+BET001
               !
            zn0 = ((((BET050*zt   &
               &   + BET140*zs+BET040)*zt   &
               &   + (BET230*zs+BET130)*zs+BET030)*zt   &
               &   + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   &
               &   + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   &
               &   + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            !
            beta(i) =  REAL(zn / zs * r1_rho0, KIND=4)
            !
         END DO
         !
      CASE( np_seos )                  !==  simplified EOS  ==!
         !
         DO i=1,n
            IF (mask(i)) THEN
               alpha(i) = fillvalue
               beta(i) = fillvalue
               CYCLE
            END IF
            zt  = T (i) - 10.d0   ! pot. temperature anomaly (t-T0)
            zs  = S (i) - 35.d0   ! abs. salinity anomaly (s-S0)
            zh  = depth(i)                ! depth in meters at t-point

            !
            zn  = rn_a0 * ( 1.d0 + rn_lambda1*zt + rn_mu1*zh ) + rn_nu*zs
            alpha(i) = REAL(zn * r1_rho0, KIND=4)   ! alpha
            !
            zn  = rn_b0 * ( 1.d0 - rn_lambda2*zs - rn_mu2*zh ) - rn_nu*zt
            beta(i) = REAL(zn * r1_rho0, KIND=4)   ! beta
         END DO
         !
      CASE( np_old_eos80 )                  !==  old EOS  ==!
         !
         DO i=1,n
            IF (mask(i)) THEN
               alpha(i) = fillvalue
               beta(i) = fillvalue
               CYCLE
            END IF
            zt = T(i)
            zs = S(i) - 35.d0
            zh = depth(i)                 ! depth
            zalbet = ( ( ( - 0.255019d-07 * zt + 0.298357d-05 ) * zt   &   ! ratio alpha/bdta
                 &                                  - 0.203814d-03 ) * zt   &
                 &                                  + 0.170907d-01 ) * zt   &
                 &   +         0.665157d-01                                 &
                 &   +     ( - 0.678662d-05 * zs                            &
                 &           - 0.846960d-04 * zt + 0.378110d-02 ) * zs   &
                 &   +   ( ( - 0.302285d-13 * zh                            &
                 &           - 0.251520d-11 * zs                            &
                 &           + 0.512857d-12 * zt * zt              ) * zh   &
                 &           - 0.164759d-06 * zs                            &
                 &        +(   0.791325d-08 * zt - 0.933746d-06 ) * zt   &
                 &                                  + 0.380374d-04 ) * zh 
            !
            beta(i) = REAL( &
                 &   ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
                 &                               - 0.301985d-05 ) * zt      &
                 &   +       0.785567d-03                                   &
                 &   + (     0.515032d-08 * zs                              &
                 &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
                 &   + ( (   0.121551d-17 * zh                              &
                 &         - 0.602281d-15 * zs                              &
                 &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
                 &                                + 0.408195d-10   * zs     &
                 &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
                 &                                - 0.121555d-07 ) * zh &
                 &   , KIND=4)
            alpha(i) = REAL( zalbet*beta(i), KIND=4)
         END DO
         !
      END SELECT
      !
   END SUBROUTINE eos_rab4_m

   SUBROUTINE eos_rab_ref4( T, S, depth_km, alpha, beta, n )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE rab_3d  ***
      !!
      !! ** Purpose :   Calculates thermal/haline expansion ratio at T-points
      !!
      !! ** Method  :   calculates (unmasked) alpha / beta at T-points
      !!
      !! ** Action  : - outputs ajpha, beta     : thermal/haline expansion ratio at T-points
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, INTENT(in)                           ::  depth_km! reference depth                    [km]
      REAL*4, DIMENSION(n), INTENT(out)            ::  alpha   ! thermal expansion coefficent          [Celsius^{-1}]
      REAL*4, DIMENSION(n), INTENT(out)            ::  beta    ! saline contraction coefficent          [psu^{-1}/kg/g]
      !
      INTEGER  ::  i
      !
      REAL*8 ::   zt , zh , zs         ! local scalars
      REAL*8 ::   zn , zn0, zn1, zn2, zn3   !   -      -
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zalbet  ! temporary scalar
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         zh  = depth_km*1.d3 * r1_Z0                                ! depth
         !
         do i=1,n
            !
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            ! alpha
            zn3 = ALP003
            !
            zn2 = ALP012*zt + ALP102*zs+ALP002
            !
            zn1 = ((ALP031*zt   &
               &   + ALP121*zs+ALP021)*zt   &
               &   + (ALP211*zs+ALP111)*zs+ALP011)*zt   &
               &   + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
               !
            zn0 = ((((ALP050*zt   &
               &   + ALP140*zs+ALP040)*zt   &
               &   + (ALP230*zs+ALP130)*zs+ALP030)*zt   &
               &   + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   &
               &   + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   &
               &   + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            !
            alpha(i) = REAL(zn * r1_rho0, KIND=4)
            !
            ! beta
            zn3 = BET003
            !
            zn2 = BET012*zt + BET102*zs+BET002
            !
            zn1 = ((BET031*zt   &
               &   + BET121*zs+BET021)*zt   &
               &   + (BET211*zs+BET111)*zs+BET011)*zt   &
               &   + ((BET301*zs+BET201)*zs+BET101)*zs+BET001
               !
            zn0 = ((((BET050*zt   &
               &   + BET140*zs+BET040)*zt   &
               &   + (BET230*zs+BET130)*zs+BET030)*zt   &
               &   + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   &
               &   + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   &
               &   + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            !
            beta(i) =  REAL(zn / zs * r1_rho0, KIND=4)
            !
         END DO
         !
      CASE( np_seos )                  !==  simplified EOS  ==!
         !
         zh  = depth_km*1.d3                                ! depth
         do i=1,n
            zt  = T (i) - 10.d0   ! pot. temperature anomaly (t-T0)
            zs  = S (i) - 35.d0   ! abs. salinity anomaly (s-S0)
            !
            zn  = rn_a0 * ( 1.d0 + rn_lambda1*zt + rn_mu1*zh ) + rn_nu*zs
            alpha(i) = REAL(zn * r1_rho0, KIND=4)   ! alpha
            !
            zn  = rn_b0 * ( 1.d0 - rn_lambda2*zs - rn_mu2*zh ) - rn_nu*zt
            beta(i) = REAL(zn * r1_rho0, KIND=4)   ! beta
            !
         END DO
         !
      CASE( np_old_eos80 )                  !==  simplified EOS  ==!
         zh  = depth_km*1.d3
         DO i=1,n
            zt = T(i)
            zs = S(i) - 35.d0
            zalbet = ( ( ( - 0.255019d-07 * zt + 0.298357d-05 ) * zt   &   ! ratio alpha/bdta
                 &                                  - 0.203814d-03 ) * zt   &
                 &                                  + 0.170907d-01 ) * zt   &
                 &   +         0.665157d-01                                 &
                 &   +     ( - 0.678662d-05 * zs                            &
                 &           - 0.846960d-04 * zt + 0.378110d-02 ) * zs   &
                 &   +   ( ( - 0.302285d-13 * zh                            &
                 &           - 0.251520d-11 * zs                            &
                 &           + 0.512857d-12 * zt * zt              ) * zh   &
                 &           - 0.164759d-06 * zs                            &
                 &        +(   0.791325d-08 * zt - 0.933746d-06 ) * zt   &
                 &                                  + 0.380374d-04 ) * zh
            !
            beta(i) = REAL( &
                 &   ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
                 &                               - 0.301985d-05 ) * zt      &
                 &   +       0.785567d-03                                   &
                 &   + (     0.515032d-08 * zs                              &
                 &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
                 &   + ( (   0.121551d-17 * zh                              &
                 &         - 0.602281d-15 * zs                              &
                 &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
                 &                                + 0.408195d-10   * zs     &
                 &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
                 &                                - 0.121555d-07 ) * zh  &
                 &    , KIND=4)
            alpha(i) = REAL(zalbet*beta(i), KIND=4)
         END DO
      END SELECT
      !
   END SUBROUTINE eos_rab_ref4


   SUBROUTINE eos_rab_ref4_m( fillvalue, mask, T, S, depth_km, alpha, beta, n )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE rab_3d  ***
      !!
      !! ** Purpose :   Calculates thermal/haline expansion ratio at T-points
      !!
      !! ** Method  :   calculates (unmasked) alpha / beta at T-points
      !!
      !! ** Action  : - outputs ajpha, beta     : thermal/haline expansion ratio at T-points
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      LOGICAL(KIND=1), DIMENSION(n), INTENT(IN)    :: mask   !logical mask of land points
      REAL*4, INTENT(IN)                           ::  fillvalue
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, INTENT(in)                           ::  depth_km! reference depth                    [km]
      REAL*4, DIMENSION(n), INTENT(out)            ::  alpha   ! thermal expansion coefficent          [Celsius^{-1}]
      REAL*4, DIMENSION(n), INTENT(out)            ::  beta    ! saline contraction coefficent          [psu^{-1}/kg/g]
      !
      INTEGER  ::  i
      !
      REAL*8 ::   zt , zh , zs         ! local scalars
      REAL*8 ::   zn , zn0, zn1, zn2, zn3   !   -      -
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zalbet  ! temporary scalar
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         zh  = depth_km*1.d3 * r1_Z0                                ! depth
         !
         DO i=1,n
            IF (mask(i)) THEN
               alpha(i) = fillvalue
               beta(i) = fillvalue
               CYCLE
            END IF
            !
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            ! alpha
            zn3 = ALP003
            !
            zn2 = ALP012*zt + ALP102*zs+ALP002
            !
            zn1 = ((ALP031*zt   &
               &   + ALP121*zs+ALP021)*zt   &
               &   + (ALP211*zs+ALP111)*zs+ALP011)*zt   &
               &   + ((ALP301*zs+ALP201)*zs+ALP101)*zs+ALP001
               !
            zn0 = ((((ALP050*zt   &
               &   + ALP140*zs+ALP040)*zt   &
               &   + (ALP230*zs+ALP130)*zs+ALP030)*zt   &
               &   + ((ALP320*zs+ALP220)*zs+ALP120)*zs+ALP020)*zt   &
               &   + (((ALP410*zs+ALP310)*zs+ALP210)*zs+ALP110)*zs+ALP010)*zt   &
               &   + ((((ALP500*zs+ALP400)*zs+ALP300)*zs+ALP200)*zs+ALP100)*zs+ALP000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            !
            alpha(i) = REAL(zn * r1_rho0, KIND=4)
            !
            ! beta
            zn3 = BET003
            !
            zn2 = BET012*zt + BET102*zs+BET002
            !
            zn1 = ((BET031*zt   &
               &   + BET121*zs+BET021)*zt   &
               &   + (BET211*zs+BET111)*zs+BET011)*zt   &
               &   + ((BET301*zs+BET201)*zs+BET101)*zs+BET001
               !
            zn0 = ((((BET050*zt   &
               &   + BET140*zs+BET040)*zt   &
               &   + (BET230*zs+BET130)*zs+BET030)*zt   &
               &   + ((BET320*zs+BET220)*zs+BET120)*zs+BET020)*zt   &
               &   + (((BET410*zs+BET310)*zs+BET210)*zs+BET110)*zs+BET010)*zt   &
               &   + ((((BET500*zs+BET400)*zs+BET300)*zs+BET200)*zs+BET100)*zs+BET000
               !
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            !
            beta(i) =  REAL(zn / zs * r1_rho0, KIND=4)
            !
         END DO
         !
      CASE( np_seos )                  !==  simplified EOS  ==!
         !
         zh  = depth_km*1.d3                                ! depth
         do i=1,n
            IF (mask(i)) THEN
               alpha(i) = fillvalue
               beta(i) = fillvalue
               CYCLE
            END IF
            zt  = T (i) - 10.d0   ! pot. temperature anomaly (t-T0)
            zs  = S (i) - 35.d0   ! abs. salinity anomaly (s-S0)

            !
            zn  = rn_a0 * ( 1.d0 + rn_lambda1*zt + rn_mu1*zh ) + rn_nu*zs
            alpha(i) = REAL(zn * r1_rho0, KIND=4)   ! alpha
            !
            zn  = rn_b0 * ( 1.d0 - rn_lambda2*zs - rn_mu2*zh ) - rn_nu*zt
            beta(i) = REAL(zn * r1_rho0, KIND=4)   ! beta
            !
         END DO
         !
      CASE( np_old_eos80 )                  !==  simplified EOS  ==!
         zh  = depth_km*1.d3
         DO i=1,n
            zt = T(i)
            zs = S(i) - 35.d0
            zalbet = ( ( ( - 0.255019d-07 * zt + 0.298357d-05 ) * zt   &   ! ratio alpha/bdta
                 &                                  - 0.203814d-03 ) * zt   &
                 &                                  + 0.170907d-01 ) * zt   &
                 &   +         0.665157d-01                                 &
                 &   +     ( - 0.678662d-05 * zs                            &
                 &           - 0.846960d-04 * zt + 0.378110d-02 ) * zs   &
                 &   +   ( ( - 0.302285d-13 * zh                            &
                 &           - 0.251520d-11 * zs                            &
                 &           + 0.512857d-12 * zt * zt              ) * zh   &
                 &           - 0.164759d-06 * zs                            &
                 &        +(   0.791325d-08 * zt - 0.933746d-06 ) * zt   &
                 &                                  + 0.380374d-04 ) * zh
            !
            beta(i) = REAL( &
                 &   ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
                 &                               - 0.301985d-05 ) * zt      &
                 &   +       0.785567d-03                                   &
                 &   + (     0.515032d-08 * zs                              &
                 &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
                 &   + ( (   0.121551d-17 * zh                              &
                 &         - 0.602281d-15 * zs                              &
                 &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
                 &                                + 0.408195d-10   * zs     &
                 &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
                 &                                - 0.121555d-07 ) * zh &
                 & , KIND=4)
            alpha(i) = REAL(zalbet*beta(i), KIND=4)
         END DO
      END SELECT
      !
    END SUBROUTINE eos_rab_ref4_m

!    ! SUBROUTINE bn2_t( pts, pab, ktab, pn2, ktn2, Kmm )
!    !    !!----------------------------------------------------------------------
!    !    !!                  ***  ROUTINE bn2  ***
!    !    !!
!    !    !! ** Purpose :   Compute the local Brunt-Vaisala frequency at the
!    !    !!                time-step of the input arguments
!    !    !!
!    !    !! ** Method  :   pn2 = grav * (alpha dk[T] + beta dk[S] ) / e3w
!    !    !!      where alpha and beta are given in pab, and computed on T-points.
!    !    !!      N.B. N^2 is set one for all to zero at jk=1 in istate module.
!    !    !!
!    !    !! ** Action  :   pn2 : square of the brunt-vaisala frequency at w-point
!    !    !!
!    !    !!----------------------------------------------------------------------
!    !    INTEGER                                , INTENT(in   ) ::  Kmm   ! time level index
!    !    INTEGER                                , INTENT(in   ) ::  ktab, ktn2
!    !    REAL*8, DIMENSION(jpi,jpj,  jpk,jpts), INTENT(in   ) ::  pts   ! pot. temperature and salinity   [Celsius,psu]
!    !    REAL*8, DIMENSION(A2D_T(ktab),JPK,JPTS), INTENT(in   ) ::  pab   ! thermal/haline expansion coef.  [Celsius-1,psu-1]
!    !    REAL*8, DIMENSION(A2D_T(ktn2),JPK     ), INTENT(  out) ::  pn2   ! Brunt-Vaisala frequency squared [1/s^2]
!    !    !
!    !    INTEGER  ::   ji, jj, jk      ! dummy loop indices
!    !    REAL*8 ::   zaw, zbw, zrw   ! local scalars
!    !    !!----------------------------------------------------------------------
!    !    !
!    !    IF( ln_timing )   CALL timing_start('bn2')
!    !    !
!    !    DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 2, jpkm1 )      ! interior points only (2=< jk =< jpkm1 ); surface and bottom value set to zero one for all in istate.F90
!    !       zrw =   ( gdepw(i  ,Kmm) - gdept(i,Kmm) )   &
!    !          &  / ( gdept(i-1,Kmm) - gdept(i,Kmm) )
!    !          !
!    !       zaw = pab(i,jp_tem) * (1. - zrw) + pab(i-1,jp_tem) * zrw
!    !       zbw = pab(i,jp_sal) * (1. - zrw) + pab(i-1,jp_sal) * zrw
!    !       !
!    !       pn2(i) = grav * (  zaw * ( pts(i-1,jp_tem) - pts(i,jp_tem) )     &
!    !          &                    - zbw * ( pts(i-1,jp_sal) - pts(i,jp_sal) )  )  &
!    !          &            / e3w(i,Kmm) * wmask(i)
!    !    END_3D
!    !    !
!    !    IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=CASTDP(pn2), clinfo1=' bn2  : ' )
!    !    !
!    !    IF( ln_timing )   CALL timing_stop('bn2')
!    !    !
!    ! END SUBROUTINE bn2_t


   SUBROUTINE eos_pot_from_CT_SA4( CT, SA, pot, n )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE eos_pot_from_CT_SA  ***
      !!
      !! ** Purpose :   Compute pot.temp. from cons. temp. [Celsius] and absolute salinity [g/kg]
      !!
      !! ** Method  :   rational approximation (5/3th order) of TEOS-10 algorithm
      !!       checkvalue: pt=20.02391895 Celsius for sa=35.7g/kg, ct=20degC
      !!
      !! Reference  :   TEOS-10, UNESCO
      !!                Rational approximation to TEOS10 algorithm (rms error on WOA13 values: 4.0d-5 degC)
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      REAL*4, DIMENSION(n), INTENT(in)             ::  CT       ! conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  SA       ! absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(out)            ::  pot     ! potential temperature [Celsius]

      INTEGER  ::   i                    ! dummy loop indices
      REAL*8 ::   zt , zs              ! local scalars
      REAL*8 ::   zn , zd              ! local scalars
      REAL*8 ::   zdeltaS , z1_S0 , z1_T0
      !!----------------------------------------------------------------------
      !
      zdeltaS = 5.d0
      z1_S0   = 0.875d0/35.16504d0
      z1_T0   = 1.d0/40.d0
      !
      DO i=1,n
         !
         zt  = CT(i) * z1_T0
         zs  = SQRT( ABS( SA(i) + zdeltaS ) * z1_S0 )
         !
         zn = ((((-2.1385727895d-01*zt   &
            &   - 2.7674419971d-01*zs+1.0728094330d0)*zt   &
            &   + (2.6366564313d0*zs+3.3546960647d0)*zs-7.8012209473d0)*zt   &
            &   + ((1.8835586562d0*zs+7.3949191679d0)*zs-3.3937395875d0)*zs-5.6414948432d0)*zt   &
            &   + (((3.5737370589d0*zs-1.5512427389d+01)*zs+2.4625741105d+01)*zs   &
            &      +1.9912291000d+01)*zs-3.2191146312d+01)*zt   &
            &   + ((((5.7153204649d-01*zs-3.0943149543d0)*zs+9.3052495181d0)*zs   &
            &      -9.4528934807d0)*zs+3.1066408996d0)*zs-4.3504021262d-01
            !
         zd = (2.0035003456d0*zt   &
            &   -3.4570358592d-01*zs+5.6471810638d0)*zt   &
            &   + (1.5393993508d0*zs-6.9394762624d0)*zs+1.2750522650d+01
            !
         pot(i) = REAL( zt / z1_T0 + zn / zd, KIND=4)
            !
      END DO
      !
   END SUBROUTINE eos_pot_from_CT_SA4


  !  SUBROUTINE  eos_fzp_2d_t( psal, ptf, kttf, pdep )
  !     !!----------------------------------------------------------------------
  !     !!                 ***  ROUTINE eos_fzp  ***
  !     !!
  !     !! ** Purpose :   Compute the freezing point temperature [Celsius]
  !     !!
  !     !! ** Method  :   UNESCO freezing point (ptf) in Celsius is given by
  !     !!       ptf(t,z) = (-.0575+1.710523d-3*sqrt(abs(s))-2.154996d-4*s)*s - 7.53d-4*z
  !     !!       checkvalue: tf=-2.588567 Celsius for s=40psu, z=500m
  !     !!
  !     !! Reference  :   UNESCO tech. papers in the marine science no. 28. 1978
  !     !!----------------------------------------------------------------------
  !     INTEGER                       , INTENT(in   )           ::   kttf
  !     REAL*8, DIMENSION(jpi,jpj)  , INTENT(in   )           ::   psal   ! salinity   [psu]
  !     REAL*8, DIMENSION(jpi,jpj)  , INTENT(in   ), OPTIONAL ::   pdep   ! depth      [m]
  !     REAL*8, DIMENSION(A2D_T(kttf)), INTENT(out  )           ::   ptf    ! freezing temperature [Celsius]
  !     !
  !     INTEGER  ::   ji, jj          ! dummy loop indices
  !     REAL*8 ::   zt, zs, z1_S0   ! local scalars
  !     !!----------------------------------------------------------------------
  !     !
  !     SELECT CASE ( neos )
  !     !
  !     CASE ( np_teos10, np_seos )      !==  CT,SA (TEOS-10 and S-EOS formulations) ==!
  !        !
  !        z1_S0 = 1.d0 / 35.16504d0
  !        DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
  !           zs= SQRT( ABS( psal(ji,jj) ) * z1_S0 )           ! square root salinity
  !           ptf(ji,jj) = ((((1.46873d-03*zs-9.64972d-03)*zs+2.28348d-02)*zs &
  !              &          - 3.12775d-02)*zs+2.07679d-02)*zs-5.87701d-02
  !        END_2D
  !        ptf(:,:) = ptf(:,:) * psal(:,:)
  !        !
  !        IF( PRESENT( pdep ) )   ptf(:,:) = ptf(:,:) - 7.53d-4 * pdep(:,:)
  !        !
  !     CASE ( np_eos80 )                !==  PT,SP (UNESCO formulation)  ==!
  !        !
  !        ptf(:,:) = ( - 0.0575d0 + 1.710523d-3 * SQRT( psal(:,:) )   &
  !           &                     - 2.154996d-4 *       psal(:,:)   ) * psal(:,:)
  !           !
  !        IF( PRESENT( pdep ) )   ptf(:,:) = ptf(:,:) - 7.53d-4 * pdep(:,:)
  !        !
  !     CASE DEFAULT
  !        WRITE(ctmp1,*) '          bad flag value for neos = ', neos
  !        CALL ctl_stop( 'eos_fzp_2d:', ctmp1 )
  !        !
  !     END SELECT
  !     !
  ! END SUBROUTINE eos_fzp_2d_t


  ! SUBROUTINE eos_fzp_0d( psal, ptf, pdep )
  !     !!----------------------------------------------------------------------
  !     !!                 ***  ROUTINE eos_fzp  ***
  !     !!
  !     !! ** Purpose :   Compute the freezing point temperature [Celsius]
  !     !!
  !     !! ** Method  :   UNESCO freezing point (ptf) in Celsius is given by
  !     !!       ptf(t,z) = (-.0575+1.710523d-3*sqrt(abs(s))-2.154996d-4*s)*s - 7.53d-4*z
  !     !!       checkvalue: tf=-2.588567 Celsius for s=40psu, z=500m
  !     !!
  !     !! Reference  :   UNESCO tech. papers in the marine science no. 28. 1978
  !     !!----------------------------------------------------------------------
  !     REAL*8, INTENT(in )           ::   psal         ! salinity   [psu]
  !     REAL*8, INTENT(in ), OPTIONAL ::   pdep         ! depth      [m]
  !     REAL*8, INTENT(out)           ::   ptf          ! freezing temperature [Celsius]
  !     !
  !     REAL*8 :: zs   ! local scalars
  !     !!----------------------------------------------------------------------
  !     !
  !     SELECT CASE ( neos )
  !     !
  !     CASE ( np_teos10, np_seos )      !==  CT,SA (TEOS-10 and S-EOS formulations) ==!
  !        !
  !        zs  = SQRT( ABS( psal ) / 35.16504d0 )           ! square root salinity
  !        ptf = ((((1.46873d-03*zs-9.64972d-03)*zs+2.28348d-02)*zs &
  !                 &          - 3.12775d-02)*zs+2.07679d-02)*zs-5.87701d-02
  !        ptf = ptf * psal
  !        !
  !        IF( PRESENT( pdep ) )   ptf = ptf - 7.53d-4 * pdep
  !        !
  !     CASE ( np_eos80 )                !==  PT,SP (UNESCO formulation)  ==!
  !        !
  !        ptf = ( - 0.0575d0 + 1.710523d-3 * SQRT( psal )   &
  !           &                - 2.154996d-4 *       psal   ) * psal
  !           !
  !        IF( PRESENT( pdep ) )   ptf = ptf - 7.53d-4 * pdep
  !        !
  !     CASE DEFAULT
  !        WRITE(ctmp1,*) '          bad flag value for neos = ', neos
  !        CALL ctl_stop( 'eos_fzp_0d:', ctmp1 )
  !        !
  !     END SELECT
  !     !
  !  END SUBROUTINE eos_fzp_0d


   SUBROUTINE eos_pen4( T, S, depth, alpha_pe, beta_pe, ppen, n )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE eos_pen  ***
      !!
      !! ** Purpose :   Calculates nonlinear anomalies of alpha_PE, beta_PE and PE at T-points
      !!
      !! ** Method  :   PE is defined analytically as the vertical
      !!                   primitive of EOS times -g integrated between 0 and z>0.
      !!                pen is the nonlinear bsq-PE anomaly: pen = ( PE - rho0 gz ) / rho0 gz - rd
      !!                                                      = 1/z * /int_0^z rd dz - rd
      !!                                where rd is the density anomaly (see eos_rhd function)
      !!                ab_pe are partial derivatives of PE anomaly with respect to T and S:
      !!                    ab_pe(1) = - 1/(rho0 gz) * dPE/dT + drd/dT = - d(pen)/dT
      !!                    ab_pe(2) =   1/(rho0 gz) * dPE/dS + drd/dS =   d(pen)/dS
      !!
      !! ** Action  : - pen         : PE anomaly given at T-points
      !!            : - pab_pe  : given at T-points
      !!                    pab_pe(:,:,:,jp_tem) is alpha_pe
      !!                    pab_pe(:,:,:,jp_sal) is beta_pe
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      REAL*4, DIMENSION(n), INTENT(in)             ::  T       ! potential/conservative temperature  [Celsius]
      REAL*4, DIMENSION(n), INTENT(in)             ::  S       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(in)             ::  depth   ! depth                                [m]
      REAL*4, DIMENSION(n), INTENT(  out)        ::  alpha_pe
      REAL*4, DIMENSION(n), INTENT(  out)        ::  beta_pe
      REAL*4, DIMENSION(n), INTENT(  out)        ::  ppen     ! potential energy anomaly
      !
      INTEGER  ::   i                         ! dummy loop indices
      REAL*8 ::   zt , zh , zs              ! local scalars
      REAL*8 ::   zn , zn0, zn1, zn2        !   -      -
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( neos )
       !
       CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         DO i=1,n
            !
            zh  = depth(i) * r1_Z0                                ! depth
            zt  = T (i) * r1_T0                           ! temperature
            zs  = SQRT( ABS( S(i) + rdeltaS ) * r1_S0 )   ! square root salinity
            !
            ! potential energy non-linear anomaly
            zn2 = (PEN012)*zt   &
               &   + PEN102*zs+PEN002
               !
            zn1 = ((PEN021)*zt   &
               &   + PEN111*zs+PEN011)*zt   &
               &   + (PEN201*zs+PEN101)*zs+PEN001
               !
            zn0 = ((((PEN040)*zt   &
               &   + PEN130*zs+PEN030)*zt   &
               &   + (PEN220*zs+PEN120)*zs+PEN020)*zt   &
               &   + ((PEN310*zs+PEN210)*zs+PEN110)*zs+PEN010)*zt   &
               &   + (((PEN400*zs+PEN300)*zs+PEN200)*zs+PEN100)*zs+PEN000
               !
            zn  = ( zn2 * zh + zn1 ) * zh + zn0
            !
            ppen(i)  = REAL(zn * zh * r1_rho0, KIND=4)
            !
            ! alphaPE non-linear anomaly
            zn2 = APE002
            !
            zn1 = (APE011)*zt   &
               &   + APE101*zs+APE001
               !
            zn0 = (((APE030)*zt   &
               &   + APE120*zs+APE020)*zt   &
               &   + (APE210*zs+APE110)*zs+APE010)*zt   &
               &   + ((APE300*zs+APE200)*zs+APE100)*zs+APE000
               !
            zn  = ( zn2 * zh + zn1 ) * zh + zn0
            !
            alpha_pe(i) = REAL(zn * zh * r1_rho0, KIND=4)
            !
            ! betaPE non-linear anomaly
            zn2 = BPE002
            !
            zn1 = (BPE011)*zt   &
               &   + BPE101*zs+BPE001
               !
            zn0 = (((BPE030)*zt   &
               &   + BPE120*zs+BPE020)*zt   &
               &   + (BPE210*zs+BPE110)*zs+BPE010)*zt   &
               &   + ((BPE300*zs+BPE200)*zs+BPE100)*zs+BPE000
               !
            zn  = ( zn2 * zh + zn1 ) * zh + zn0
            !
            beta_pe(i) = REAL(zn / zs * zh * r1_rho0, KIND=4)
            !
         END DO
         !
      CASE( np_seos )                !==  Vallis (2006) simplified EOS  ==!
         !
         Do i=1,n
            zt  = T(i) - 10.d0  ! temperature anomaly (t-T0)
            zs =  S(i) - 35.d0  ! abs. salinity anomaly (s-S0)
            zh  = depth(i)              ! depth in meters  at t-point

            zn  = 0.5d0 * zh * r1_rho0
            !                                    ! Potential Energy
            ppen(i) = REAL(( rn_a0 * rn_mu1 * zt + rn_b0 * rn_mu2 * zs ) * zn, KIND=4)
            !                                    ! alphaPE
            alpha_pe(i) = - REAL(rn_a0 * rn_mu1 * zn, KIND=4)
            beta_pe(i) =   REAL(rn_b0 * rn_mu2 * zn, KIND=4)
            !
         END DO
         !
      END SELECT
   END SUBROUTINE eos_pen4

   SUBROUTINE eos_insitu04(T0, S0, depth_km, depth, drho0, n)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE eos_insitu04  ***
      !!
      !! ** Purpose :   Compute the unmasked difference in density between the in-situ density
      !!                  rho(T0,S0,z) and the reference density rho(T0, S0, 1.e3*depth_km).
      !!
      !! ** Method  :   drho(T0,S,z) = rho(T0,S0,z) - rho(T0, S0, 1.e3*depth_km)
      !!                t      TEOS10: CT or EOS80: PT      Celsius
      !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
      !!                z      depth                        meters
      !!                drho    in situ density anomaly      kg/m^3
      !!
      !!     ln_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
      !!               Note global mean r0(z)
      !!         Check value: rho = 28.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
      !!
      !!     ln_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
      !!         Check value: rho = 28.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
      !!
      !!     ln_seos : simplified equation of state
      !!              rho(t,s,z) = rho0 -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0) - 1000.
      !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
      !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
      !!              Vallis like equation: use default values of coefficients
      !!
      !! ** Action  :   compute rho , the in situ density anomaly (kg/m^3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      REAL*4, INTENT(in)                       ::  T0       ! potential/conservative temperature  [Celsius]
      REAL*4, INTENT(in)                       ::  S0       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(in)         ::  depth    ! depth                                 [m]
      REAL*4, INTENT (IN)                      ::  depth_km ! reference depth
      REAL*4, DIMENSION(n), INTENT(out)        ::  drho0    ! deviation of rho0 from value at reference depth  [kg/m^3]
      !
      INTEGER  ::  i
      REAL*8  :: zt, zh ! local scalars
      REAL*8  :: zs! local scalars
      REAL*8  :: zn1, zn2!   -      -
      REAL*8  :: zn, zn0, zn3!   -      -
      REAL*8  :: rho00
      REAL*8  :: zr0
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
      REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0        !    -         -
      !!----------------------------------------------------------------------
      !
      SELECT CASE( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         zt  = T0 * r1_T0                           ! temperature
         zs  = SQRT( ABS( S0 + rdeltaS ) * r1_S0 )   ! square root salinity
         !
         zn3 = EOS013*zt   &
            &   + EOS103*zs+EOS003
            !
         zn2 = (EOS022*zt   &
            &   + EOS112*zs+EOS012)*zt   &
            &   + (EOS202*zs+EOS102)*zs+EOS002
            !
         zn1 = (((EOS041*zt   &
            &   + EOS131*zs+EOS031)*zt   &
            &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
            &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
            &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
            !
         zn0 = (((((EOS060*zt   &
            &   + EOS150*zs+EOS050)*zt   &
            &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
            &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
            &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
            &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
            &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
            !
         IF (depth_km < 1.e-4) THEN
            rho00 = zn0 - 1000.d0
         ELSE
            zh  = depth_km * 1.d3 * r1_Z0 
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            IF (neos == np_teos10) THEN
            !!  Add reference profile zr0 to anomaly
               zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
               rho00 = zn - 1000.d0 + zr0
            ELSE
               rho00 = zn - 1000.d0
            END IF
         END IF
         DO i=1,n
            !
            zh  = depth(i) * r1_Z0                                  ! depth
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
    !!----------------------------------------------------------------------
            !!  Subtract 1000. to give density anomaly
            IF (neos == np_teos10) THEN
            !!  Add reference profile zr0 to anomaly
               zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
               drho0(i) = REAL(zn - 1000.d0 + zr0 - rho00, KIND=4)
            ELSE
               drho0(i) = REAL(zn - 1000.d0 - rho00, KIND=4)
            END IF
    !!----------------------------------------------------------------------
         END DO

      CASE( np_seos )                !==  simplified EOS  ==!
         !
         zt  = T0
         zs  = S0
         zh  = 1.d3*depth_km
         zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
            &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
            &  - rn_nu * zt * zs
            !
         rho00 = zn - 1000.d0
         !
         DO i=1,n
            zh = depth(i)
            zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
               &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
               &  - rn_nu * zt * zs
               !
            ! prd(i) = zn * r1_rho0 * ztm                ! density anomaly (masked)
            drho0(i) = REAL(zn + rho0 - 1000.d0, KIND=4)
          END DO
         !
      CASE(np_old_eos80)
         zt = T0
         zs = S0
         zsr= SQRT( ABS( zs) )        ! square root salinity
         !
         ! compute volumic mass pure water at atm pressure
         zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4 )*zt   &
              &                          -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594d0
         ! seawater volumic mass atm pressure
         zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
              &                                         -4.0899e-3 ) *zt+0.824493d0
         zr3= ( -1.6546e-6*zt+1.0227e-4 )    *zt-5.72466e-3
         zr4= 4.8314e-4
         !
         ! potential volumic mass (reference to the surface)
         zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
         !
         ! add the compression terms
         ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
         zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
         zb = zbw + ze * zs
         !
         zd = -2.042967e-2
         zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
         zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859d0 ) *zt - 4.721788d0
         za = ( zd*zsr + zc ) *zs + zaw
         !
         zb1=   (  -0.1909078d0  *zt+7.390729d0    ) *zt-55.87545d0
         za1= ( (   2.326469e-3*zt+1.553190d0    ) *zt-65.00517d0 ) *zt + 1044.077d0
         zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638d0 ) *zt + 2098.925d0 ) *zt+190925.6d0
         zk0= ( zb1*zsr + za1 )*zs + zkw

         zh = 1000.d0*depth_km
         rho00 = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )

         DO i=1,n
            zh = depth(i)
            ! unmasked in situ density anomaly
            drho0(i) = REAL(zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rho00, KIND=4)
            END DO
      END SELECT
      !
   END SUBROUTINE eos_insitu04

   SUBROUTINE eos_insitu04_m(fillvalue, mask, T0, S0, depth_km, depth, drho0, n)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE eos_insitu04  ***
      !!
      !! ** Purpose :   Compute the masked difference in density between the in-situ density
      !!                  rho(T0,S0,z) and the reference density rho(T0, S0, 1.e3*depth_km).
      !!
      !! ** Method  :   drho(T0,S,z) = rho(T0,S0,z) - rho(T0, S0, 1.e3*depth_km)
      !!                t      TEOS10: CT or EOS80: PT      Celsius
      !!                s      TEOS10: SA or EOS80: SP      TEOS10: g/kg or EOS80: psu
      !!                z      depth                        meters
      !!                drho    in situ density anomaly      kg/m^3
      !!
      !!     ln_teos10 : polynomial TEOS-10 equation of state is used for rho(t,s,z).
      !!               Note global mean r0(z)
      !!         Check value: rho = 28.21993233072 kg/m^3 for z=3000 dbar, ct=3 Celsius, sa=35.5 g/kg
      !!
      !!     ln_eos80 : polynomial EOS-80 equation of state is used for rho(t,s,z).
      !!         Check value: rho = 28.35011066567 kg/m^3 for z=3000 dbar, pt=3 Celsius, sp=35.5 psu
      !!
      !!     ln_seos : simplified equation of state
      !!              rho(t,s,z) = rho0 -a0*(1+lambda/2*(T-T0)+mu*z+nu*(S-S0))*(T-T0) + b0*(S-S0) - 1000.
      !!              linear case function of T only: rn_alpha<>0, other coefficients = 0
      !!              linear eos function of T and S: rn_alpha and rn_beta<>0, other coefficients=0
      !!              Vallis like equation: use default values of coefficients
      !!
      !! ** Action  :   compute rho , the in situ density anomaly (kg/m^3)
      !!
      !! References :   Roquet et al, Ocean Modelling (2015)
      !!                Vallis, Atmospheric and Oceanic Fluid Dynamics, 2006
      !!                TEOS-10 Manual, 2010
      !!----------------------------------------------------------------------
      INTEGER*4, INTENT(IN)                        ::   n
      LOGICAL(KIND=1), DIMENSION(n), INTENT(IN)    ::  mask   !logical mask of land points
      REAL*4, INTENT(IN)                       ::  fillvalue
      REAL*4, INTENT(in)                       ::  T0       ! potential/conservative temperature  [Celsius]
      REAL*4, INTENT(in)                       ::  S0       ! practical/absolute salinity        [psu/g/kg]
      REAL*4, DIMENSION(n), INTENT(in)         ::  depth    ! depth                                 [m]
      REAL*4, INTENT (IN)                      ::  depth_km ! reference depth
      REAL*4, DIMENSION(n), INTENT(out)        ::  drho0    ! deviation of rho0 from value at reference depth  [kg/m^3]
      !
      INTEGER  ::  i
      REAL*8  :: zt, zh ! local scalars
      REAL*8  :: zs! local scalars
      REAL*8  :: zn1, zn2!   -      -
      REAL*8  :: zn, zn0, zn3!   -      -
      REAL*8  :: rho00
      REAL*8  :: zr0
      ! For old_eos80 (Jackett and McDougall, J. Atmos. Ocean. Tech., 1994)
      REAL*8 ::   zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
      REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0        !    -         -
      !!----------------------------------------------------------------------
      !
      SELECT CASE( neos )
      !
      CASE( np_teos10, np_eos80 )                !==  polynomial TEOS-10 / EOS-80 ==!
         !
         zt  = T0 * r1_T0                           ! temperature
         zs  = SQRT( ABS( S0 + rdeltaS ) * r1_S0 )   ! square root salinity
         !
         zn3 = EOS013*zt   &
            &   + EOS103*zs+EOS003
            !
         zn2 = (EOS022*zt   &
            &   + EOS112*zs+EOS012)*zt   &
            &   + (EOS202*zs+EOS102)*zs+EOS002
            !
         zn1 = (((EOS041*zt   &
            &   + EOS131*zs+EOS031)*zt   &
            &   + (EOS221*zs+EOS121)*zs+EOS021)*zt   &
            &   + ((EOS311*zs+EOS211)*zs+EOS111)*zs+EOS011)*zt   &
            &   + (((EOS401*zs+EOS301)*zs+EOS201)*zs+EOS101)*zs+EOS001
            !
         zn0 = (((((EOS060*zt   &
            &   + EOS150*zs+EOS050)*zt   &
            &   + (EOS240*zs+EOS140)*zs+EOS040)*zt   &
            &   + ((EOS330*zs+EOS230)*zs+EOS130)*zs+EOS030)*zt   &
            &   + (((EOS420*zs+EOS320)*zs+EOS220)*zs+EOS120)*zs+EOS020)*zt   &
            &   + ((((EOS510*zs+EOS410)*zs+EOS310)*zs+EOS210)*zs+EOS110)*zs+EOS010)*zt   &
            &   + (((((EOS600*zs+EOS500)*zs+EOS400)*zs+EOS300)*zs+EOS200)*zs+EOS100)*zs+EOS000
         
         IF (depth_km < 1.e-4) THEN
            rho00 = zn0 - 1000.d0
         ELSE
            zh  = depth_km * 1.d3 * r1_Z0 
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
            IF (neos == np_teos10) THEN
            !!  Add reference profile zr0 to anomaly
               zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
               rho00 = zn - 1000.d0 + zr0
            ELSE
               rho00 = zn - 1000.d0
            END IF
         END IF
         DO i=1,n
            IF (mask(i)) THEN
               drho0(i) = fillvalue
               CYCLE
            END IF
            !
            zh  = depth(i) * r1_Z0                                  ! depth
            zn  = ( ( zn3 * zh + zn2 ) * zh + zn1 ) * zh + zn0
    !!----------------------------------------------------------------------
            !!  Subtract 1000. to give density anomaly
            IF (neos == np_teos10) THEN
            !!  Add reference profile zr0 to anomaly
               zr0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
               drho0(i) = REAL(zn - 1000.d0 + zr0 - rho00, KIND=4)
            ELSE
               drho0(i) = REAL(zn - 1000.d0 - rho00, KIND=4)
            END IF
    !!----------------------------------------------------------------------
            !
            ! prd(i) = (  zn * r1_rho0 - 1.d0  ) * ztm  ! density anomaly (masked)
            !
         END DO
         !
      CASE( np_seos )                !==  simplified EOS  ==!
         !
         zt  = T0
         zs  = S0
         zh  = 1.d3*depth_km
         zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
            &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
            &  - rn_nu * zt * zs
            !
         rho00 = zn - 1000.d0
         !
         DO i=1,n
            IF (mask(i)) THEN
               drho0(i) = fillvalue
               CYCLE
            END IF
            zh = depth(i)
            zn =  - rn_a0 * ( 1.d0 + 0.5d0*rn_lambda1*zt + rn_mu1*zh ) * zt   &
               &  + rn_b0 * ( 1.d0 - 0.5d0*rn_lambda2*zs - rn_mu2*zh ) * zs   &
               &  - rn_nu * zt * zs
               !
            ! prd(i) = zn * r1_rho0 * ztm                ! density anomaly (masked)
            drho0(i) = REAL(zn + rho0 - 1000.d0, KIND=4)
         END DO
         !
      CASE(np_old_eos80)
         zt = T0
         zs = S0
         zsr= SQRT( ABS( zs) )        ! square root salinity
         !
         ! compute volumic mass pure water at atm pressure
         zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4 )*zt   &
              &                          -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594d0
         ! seawater volumic mass atm pressure
         zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
              &                                         -4.0899e-3 ) *zt+0.824493d0
         zr3= ( -1.6546e-6*zt+1.0227e-4 )    *zt-5.72466e-3
         zr4= 4.8314e-4
         !
         ! potential volumic mass (reference to the surface)
         zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
         !
         ! add the compression terms
         ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
         zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
         zb = zbw + ze * zs
         !
         zd = -2.042967e-2
         zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
         zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859d0 ) *zt - 4.721788d0
         za = ( zd*zsr + zc ) *zs + zaw
         !
         zb1=   (  -0.1909078d0  *zt+7.390729d0    ) *zt-55.87545d0
         za1= ( (   2.326469e-3*zt+1.553190d0    ) *zt-65.00517d0 ) *zt + 1044.077d0
         zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638d0 ) *zt + 2098.925d0 ) *zt+190925.6d0
         zk0= ( zb1*zsr + za1 )*zs + zkw

         zh = 1000.d0*depth_km
         rho00 = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )

         DO i=1,n
            IF (mask(i)) THEN
               drho0(i) = fillvalue
               CYCLE
            END IF
            zh = depth(i)
            ! masked in situ density anomaly
            drho0(i) = REAL(zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rho00, KIND=4)
            END DO
      END SELECT
      !
   END SUBROUTINE eos_insitu04_m

   SUBROUTINE eos_init(neos_in)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE eos_init  ***
      !!
      !! ** Purpose :   initializations for the equation of state
      !!
      !! ** Method  :   Read the namelist nameos and control the parameters
      !!----------------------------------------------------------------------
      ! NAMELIST/nameos/ rn_a0, rn_b0, rn_lambda1, rn_mu1,   &
      !    &                                             rn_lambda2, rn_mu2, rn_nu
      !!----------------------------------------------------------------------
      !
      !
      INTEGER*4 neos_in


      neos = neos_in
      SELECT CASE( neos )         ! check option
      !
      CASE( np_teos10 )                       !==  polynomial TEOS-10  ==!
         rdeltaS = 32.d0
         r1_S0  = 0.875d0/35.16504d0
         r1_T0  = 1.d0/40.d0
         r1_Z0  = 1.d-4
         !
         EOS000 = 8.0189615746d+02
         EOS100 = 8.6672408165d+02
         EOS200 = -1.7864682637d+03
         EOS300 = 2.0375295546d+03
         EOS400 = -1.2849161071d+03
         EOS500 = 4.3227585684d+02
         EOS600 = -6.0579916612d+01
         EOS010 = 2.6010145068d+01
         EOS110 = -6.5281885265d+01
         EOS210 = 8.1770425108d+01
         EOS310 = -5.6888046321d+01
         EOS410 = 1.7681814114d+01
         EOS510 = -1.9193502195d0
         EOS020 = -3.7074170417d+01
         EOS120 = 6.1548258127d+01
         EOS220 = -6.0362551501d+01
         EOS320 = 2.9130021253d+01
         EOS420 = -5.4723692739d0
         EOS030 = 2.1661789529d+01
         EOS130 = -3.3449108469d+01
         EOS230 = 1.9717078466d+01
         EOS330 = -3.1742946532d0
         EOS040 = -8.3627885467d0
         EOS140 = 1.1311538584d+01
         EOS240 = -5.3563304045d0
         EOS050 = 5.4048723791d-01
         EOS150 = 4.8169980163d-01
         EOS060 = -1.9083568888d-01
         EOS001 = 1.9681925209d+01
         EOS101 = -4.2549998214d+01
         EOS201 = 5.0774768218d+01
         EOS301 = -3.0938076334d+01
         EOS401 = 6.6051753097d0
         EOS011 = -1.3336301113d+01
         EOS111 = -4.4870114575d0
         EOS211 = 5.0042598061d0
         EOS311 = -6.5399043664d-01
         EOS021 = 6.7080479603d0
         EOS121 = 3.5063081279d0
         EOS221 = -1.8795372996d0
         EOS031 = -2.4649669534d0
         EOS131 = -5.5077101279d-01
         EOS041 = 5.5927935970d-01
         EOS002 = 2.0660924175d0
         EOS102 = -4.9527603989d0
         EOS202 = 2.5019633244d0
         EOS012 = 2.0564311499d0
         EOS112 = -2.1311365518d-01
         EOS022 = -1.2419983026d0
         EOS003 = -2.3342758797d-02
         EOS103 = -1.8507636718d-02
         EOS013 = 3.7969820455d-01
         !
         ALP000 = -6.5025362670d-01
         ALP100 = 1.6320471316d0
         ALP200 = -2.0442606277d0
         ALP300 = 1.4222011580d0
         ALP400 = -4.4204535284d-01
         ALP500 = 4.7983755487d-02
         ALP010 = 1.8537085209d0
         ALP110 = -3.0774129064d0
         ALP210 = 3.0181275751d0
         ALP310 = -1.4565010626d0
         ALP410 = 2.7361846370d-01
         ALP020 = -1.6246342147d0
         ALP120 = 2.5086831352d0
         ALP220 = -1.4787808849d0
         ALP320 = 2.3807209899d-01
         ALP030 = 8.3627885467d-01
         ALP130 = -1.1311538584d0
         ALP230 = 5.3563304045d-01
         ALP040 = -6.7560904739d-02
         ALP140 = -6.0212475204d-02
         ALP050 = 2.8625353333d-02
         ALP001 = 3.3340752782d-01
         ALP101 = 1.1217528644d-01
         ALP201 = -1.2510649515d-01
         ALP301 = 1.6349760916d-02
         ALP011 = -3.3540239802d-01
         ALP111 = -1.7531540640d-01
         ALP211 = 9.3976864981d-02
         ALP021 = 1.8487252150d-01
         ALP121 = 4.1307825959d-02
         ALP031 = -5.5927935970d-02
         ALP002 = -5.1410778748d-02
         ALP102 = 5.3278413794d-03
         ALP012 = 6.2099915132d-02
         ALP003 = -9.4924551138d-03
         !
         BET000 = 1.0783203594d+01
         BET100 = -4.4452095908d+01
         BET200 = 7.6048755820d+01
         BET300 = -6.3944280668d+01
         BET400 = 2.6890441098d+01
         BET500 = -4.5221697773d0
         BET010 = -8.1219372432d-01
         BET110 = 2.0346663041d0
         BET210 = -2.1232895170d0
         BET310 = 8.7994140485d-01
         BET410 = -1.1939638360d-01
         BET020 = 7.6574242289d-01
         BET120 = -1.5019813020d0
         BET220 = 1.0872489522d0
         BET320 = -2.7233429080d-01
         BET030 = -4.1615152308d-01
         BET130 = 4.9061350869d-01
         BET230 = -1.1847737788d-01
         BET040 = 1.4073062708d-01
         BET140 = -1.3327978879d-01
         BET050 = 5.9929880134d-03
         BET001 = -5.2937873009d-01
         BET101 = 1.2634116779d0
         BET201 = -1.1547328025d0
         BET301 = 3.2870876279d-01
         BET011 = -5.5824407214d-02
         BET111 = 1.2451933313d-01
         BET211 = -2.4409539932d-02
         BET021 = 4.3623149752d-02
         BET121 = -4.6767901790d-02
         BET031 = -6.8523260060d-03
         BET002 = -6.1618945251d-02
         BET102 = 6.2255521644d-02
         BET012 = -2.6514181169d-03
         BET003 = -2.3025968587d-04
         !
         PEN000 = -9.8409626043d0
         PEN100 = 2.1274999107d+01
         PEN200 = -2.5387384109d+01
         PEN300 = 1.5469038167d+01
         PEN400 = -3.3025876549d0
         PEN010 = 6.6681505563d0
         PEN110 = 2.2435057288d0
         PEN210 = -2.5021299030d0
         PEN310 = 3.2699521832d-01
         PEN020 = -3.3540239802d0
         PEN120 = -1.7531540640d0
         PEN220 = 9.3976864981d-01
         PEN030 = 1.2324834767d0
         PEN130 = 2.7538550639d-01
         PEN040 = -2.7963967985d-01
         PEN001 = -1.3773949450d0
         PEN101 = 3.3018402659d0
         PEN201 = -1.6679755496d0
         PEN011 = -1.3709540999d0
         PEN111 = 1.4207577012d-01
         PEN021 = 8.2799886843d-01
         PEN002 = 1.7507069098d-02
         PEN102 = 1.3880727538d-02
         PEN012 = -2.8477365341d-01
         !
         APE000 = -1.6670376391d-01
         APE100 = -5.6087643219d-02
         APE200 = 6.2553247576d-02
         APE300 = -8.1748804580d-03
         APE010 = 1.6770119901d-01
         APE110 = 8.7657703198d-02
         APE210 = -4.6988432490d-02
         APE020 = -9.2436260751d-02
         APE120 = -2.0653912979d-02
         APE030 = 2.7963967985d-02
         APE001 = 3.4273852498d-02
         APE101 = -3.5518942529d-03
         APE011 = -4.1399943421d-02
         APE002 = 7.1193413354d-03
         !
         BPE000 = 2.6468936504d-01
         BPE100 = -6.3170583896d-01
         BPE200 = 5.7736640125d-01
         BPE300 = -1.6435438140d-01
         BPE010 = 2.7912203607d-02
         BPE110 = -6.2259666565d-02
         BPE210 = 1.2204769966d-02
         BPE020 = -2.1811574876d-02
         BPE120 = 2.3383950895d-02
         BPE030 = 3.4261630030d-03
         BPE001 = 4.1079296834d-02
         BPE101 = -4.1503681096d-02
         BPE011 = 1.7676120780d-03
         BPE002 = 1.7269476440d-04
         !
      CASE( np_eos80 )                        !==  polynomial EOS-80 formulation  ==!
         !
         rdeltaS = 20.d0
         r1_S0  = 1.d0/40.d0
         r1_T0  = 1.d0/40.d0
         r1_Z0  = 1.d-4
         !
         EOS000 = 9.5356891948d+02
         EOS100 = 1.7136499189d+02
         EOS200 = -3.7501039454d+02
         EOS300 = 5.1856810420d+02
         EOS400 = -3.7264470465d+02
         EOS500 = 1.4302533998d+02
         EOS600 = -2.2856621162d+01
         EOS010 = 1.0087518651d+01
         EOS110 = -1.3647741861d+01
         EOS210 = 8.8478359933d0
         EOS310 = -7.2329388377d0
         EOS410 = 1.4774410611d0
         EOS510 = 2.0036720553d-01
         EOS020 = -2.5579830599d+01
         EOS120 = 2.4043512327d+01
         EOS220 = -1.6807503990d+01
         EOS320 = 8.3811577084d0
         EOS420 = -1.9771060192d0
         EOS030 = 1.6846451198d+01
         EOS130 = -2.1482926901d+01
         EOS230 = 1.0108954054d+01
         EOS330 = -6.2675951440d-01
         EOS040 = -8.0812310102d0
         EOS140 = 1.0102374985d+01
         EOS240 = -4.8340368631d0
         EOS050 = 1.2079167803d0
         EOS150 = 1.1515380987d-01
         EOS060 = -2.4520288837d-01
         EOS001 = 1.0748601068d+01
         EOS101 = -1.7817043500d+01
         EOS201 = 2.2181366768d+01
         EOS301 = -1.6750916338d+01
         EOS401 = 4.1202230403d0
         EOS011 = -1.5852644587d+01
         EOS111 = -7.6639383522d-01
         EOS211 = 4.1144627302d0
         EOS311 = -6.6955877448d-01
         EOS021 = 9.9994861860d0
         EOS121 = -1.9467067787d-01
         EOS221 = -1.2177554330d0
         EOS031 = -3.4866102017d0
         EOS131 = 2.2229155620d-01
         EOS041 = 5.9503008642d-01
         EOS002 = 1.0375676547d0
         EOS102 = -3.4249470629d0
         EOS202 = 2.0542026429d0
         EOS012 = 2.1836324814d0
         EOS112 = -3.4453674320d-01
         EOS022 = -1.2548163097d0
         EOS003 = 1.8729078427d-02
         EOS103 = -5.7238495240d-02
         EOS013 = 3.8306136687d-01
         !
         ALP000 = -2.5218796628d-01
         ALP100 = 3.4119354654d-01
         ALP200 = -2.2119589983d-01
         ALP300 = 1.8082347094d-01
         ALP400 = -3.6936026529d-02
         ALP500 = -5.0091801383d-03
         ALP010 = 1.2789915300d0
         ALP110 = -1.2021756164d0
         ALP210 = 8.4037519952d-01
         ALP310 = -4.1905788542d-01
         ALP410 = 9.8855300959d-02
         ALP020 = -1.2634838399d0
         ALP120 = 1.6112195176d0
         ALP220 = -7.5817155402d-01
         ALP320 = 4.7006963580d-02
         ALP030 = 8.0812310102d-01
         ALP130 = -1.0102374985d0
         ALP230 = 4.8340368631d-01
         ALP040 = -1.5098959754d-01
         ALP140 = -1.4394226233d-02
         ALP050 = 3.6780433255d-02
         ALP001 = 3.9631611467d-01
         ALP101 = 1.9159845880d-02
         ALP201 = -1.0286156825d-01
         ALP301 = 1.6738969362d-02
         ALP011 = -4.9997430930d-01
         ALP111 = 9.7335338937d-03
         ALP211 = 6.0887771651d-02
         ALP021 = 2.6149576513d-01
         ALP121 = -1.6671866715d-02
         ALP031 = -5.9503008642d-02
         ALP002 = -5.4590812035d-02
         ALP102 = 8.6134185799d-03
         ALP012 = 6.2740815484d-02
         ALP003 = -9.5765341718d-03
         !
         BET000 = 2.1420623987d0
         BET100 = -9.3752598635d0
         BET200 = 1.9446303907d+01
         BET300 = -1.8632235232d+01
         BET400 = 8.9390837485d0
         BET500 = -1.7142465871d0
         BET010 = -1.7059677327d-01
         BET110 = 2.2119589983d-01
         BET210 = -2.7123520642d-01
         BET310 = 7.3872053057d-02
         BET410 = 1.2522950346d-02
         BET020 = 3.0054390409d-01
         BET120 = -4.2018759976d-01
         BET220 = 3.1429341406d-01
         BET320 = -9.8855300959d-02
         BET030 = -2.6853658626d-01
         BET130 = 2.5272385134d-01
         BET230 = -2.3503481790d-02
         BET040 = 1.2627968731d-01
         BET140 = -1.2085092158d-01
         BET050 = 1.4394226233d-03
         BET001 = -2.2271304375d-01
         BET101 = 5.5453416919d-01
         BET201 = -6.2815936268d-01
         BET301 = 2.0601115202d-01
         BET011 = -9.5799229402d-03
         BET111 = 1.0286156825d-01
         BET211 = -2.5108454043d-02
         BET021 = -2.4333834734d-03
         BET121 = -3.0443885826d-02
         BET031 = 2.7786444526d-03
         BET002 = -4.2811838287d-02
         BET102 = 5.1355066072d-02
         BET012 = -4.3067092900d-03
         BET003 = -7.1548119050d-04
         !
         PEN000 = -5.3743005340d0
         PEN100 = 8.9085217499d0
         PEN200 = -1.1090683384d+01
         PEN300 = 8.3754581690d0
         PEN400 = -2.0601115202d0
         PEN010 = 7.9263222935d0
         PEN110 = 3.8319691761d-01
         PEN210 = -2.0572313651d0
         PEN310 = 3.3477938724d-01
         PEN020 = -4.9997430930d0
         PEN120 = 9.7335338937d-02
         PEN220 = 6.0887771651d-01
         PEN030 = 1.7433051009d0
         PEN130 = -1.1114577810d-01
         PEN040 = -2.9751504321d-01
         PEN001 = -6.9171176978d-01
         PEN101 = 2.2832980419d0
         PEN201 = -1.3694684286d0
         PEN011 = -1.4557549876d0
         PEN111 = 2.2969116213d-01
         PEN021 = 8.3654420645d-01
         PEN002 = -1.4046808820d-02
         PEN102 = 4.2928871430d-02
         PEN012 = -2.8729602515d-01
         !
         APE000 = -1.9815805734d-01
         APE100 = -9.5799229402d-03
         APE200 = 5.1430784127d-02
         APE300 = -8.3694846809d-03
         APE010 = 2.4998715465d-01
         APE110 = -4.8667669469d-03
         APE210 = -3.0443885826d-02
         APE020 = -1.3074788257d-01
         APE120 = 8.3359333577d-03
         APE030 = 2.9751504321d-02
         APE001 = 3.6393874690d-02
         APE101 = -5.7422790533d-03
         APE011 = -4.1827210323d-02
         APE002 = 7.1824006288d-03
         !
         BPE000 = 1.1135652187d-01
         BPE100 = -2.7726708459d-01
         BPE200 = 3.1407968134d-01
         BPE300 = -1.0300557601d-01
         BPE010 = 4.7899614701d-03
         BPE110 = -5.1430784127d-02
         BPE210 = 1.2554227021d-02
         BPE020 = 1.2166917367d-03
         BPE120 = 1.5221942913d-02
         BPE030 = -1.3893222263d-03
         BPE001 = 2.8541225524d-02
         BPE101 = -3.4236710714d-02
         BPE011 = 2.8711395266d-03
         BPE002 = 5.3661089288d-04
         !
         IF(neos == np_teos10) r1_S0  = 0.875d0/35.16504d0   ! Used to convert CT to potential temperature when using bulk formulae
                                         !   (eos_pot_from_CT)
      CASE( np_seos )                        !==  Simplified EOS     ==!

      END SELECT
      !
      r1_rho0     = 1.d0 / rho0
      ! rcp         = 3991.86795711963d0      !: heat capacity     [J/K]
      ! rho0_rcp    = rho0 * rcp
      ! r1_rcp      = 1.d0 / rcp
      ! r1_rho0_rcp = 1.d0 / rho0_rcp
      ! !
      !
   END SUBROUTINE eos_init

!    !!======================================================================
END MODULE eosbn2
