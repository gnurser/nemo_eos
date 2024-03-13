module eos
  ! to compile for python on linux with ifort, do
  ! dof2py2 -x '-openmp -D__OPENMP -fast ' --open_mp nemo_rho.F90
  USE OMP_LIB, ONLY : omp_set_num_threads
  IMPLICIT NONE
  REAL*8 :: grav  = 9.80665
  REAL*8 :: rn_alpha = 2.e-4
  REAL*8 ::rau0 = 1035.d0
  REAL*8 :: rn_beta = 7.7e-4
  INTEGER*4 :: neos = 0
  REAL*8 :: fillvalue = 0.d0
  ! real*8, allocatable :: re3w(:,:,:)
  ! LOGICAL(KIND=1), allocatable :: mask3d(:,:,:)
  REAL*8 :: S00 = 35.
  REAL*8 :: theta00 = 10.
  REAL*8 :: rho000 = 26.9524116
contains
  SUBROUTINE set_eos_threads(nthreads)
    INTEGER*4, INTENT(IN)  :: nthreads
    CALL omp_set_num_threads(nthreads)
  END SUBROUTINE set_eos_threads

  SUBROUTINE set_eos(neos_in)
    INTEGER*4, INTENT(IN)  :: neos_in
    neos = neos_in
  END SUBROUTINE set_eos

  SUBROUTINE rho_mn8( fillvalue,mask,theta,S,depth,rho,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_insitu_pot  ***
    !!
    !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
    !!      potential volumic mass (Kg/m3) from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter nn_eos.
    !!
    !! ** Method  :
    !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!
    !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
    !!          t = 40 deg celcius, s=40 psu
    !!
    !!
    !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
    !!                Brown and Campana, Mon. Weather Rev., 1978
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*8, INTENT(IN) ::   fillvalue,theta(n),S(n)
    REAL*4, INTENT(IN) ::   depth(n)
    REAL*8, INTENT(OUT) ::  rho(n)
    !f2py intent (in) theta,S,depth,n
    !f2py intent (out) rho

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zrau0r       !    -         -
    !!----------------------------------------------------------------------

    !$omp parallel do private(zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze, zbw,zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
    do i=1,n
       if (mask(i)) then
          rho(i) = fillvalue
          cycle
       end if

       zt = theta(i)
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
       rho(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0
    end do
    !$omp  end parallel do
    !
  END SUBROUTINE rho_mn8

  SUBROUTINE rho_mn4( fillvalue,mask,theta,S,depth,rho,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_insitu_pot  ***
    !!
    !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
    !!      potential volumic mass (Kg/m3) from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter nn_eos.
    !!
    !! ** Method  :
    !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!
    !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
    !!          t = 40 deg celcius, s=40 psu
    !!
    !!
    !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
    !!                Brown and Campana, Mon. Weather Rev., 1978
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue,theta(n),S(n)
    REAL*4, INTENT(IN) ::   depth(n)
    REAL*4, INTENT(OUT) ::  rho(n)
    !f2py intent (in) theta,S,depth,n
    !f2py intent (out) rho

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zrau0r       !    -         -
    !!----------------------------------------------------------------------

    !$omp parallel do private(zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze, zbw,zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
    do i=1,n
       if (mask(i)) then
          rho(i) = fillvalue
          cycle
       end if

       zt = theta(i)
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
       rho(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0
    end do
    !$omp  end parallel do
    !
  END SUBROUTINE rho_mn4

  SUBROUTINE drho0_mn4( fillvalue,mask,theta0,S0,depth_km,depth,drho0,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_insitu_pot  ***
    !!
    !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
    !!      potential volumic mass (Kg/m3) from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter nn_eos.
    !!
    !! ** Method  :
    !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!
    !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
    !!          t = 40 deg celcius, s=40 psu
    !!
    !!
    !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
    !!                Brown and Campana, Mon. Weather Rev., 1978
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue,theta0,S0,depth(n)
    REAL*8, INTENT (IN) :: depth_km
    REAL*4, INTENT(OUT) ::  drho0(n)
    !f2py intent (in) theta0,S0,depth_km,depth,n
    !f2py intent (out) drho0

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, rho00       !    -         -
    !!----------------------------------------------------------------------

    zt = theta0
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

   !$omp parallel do private(zh)
    do i=1,n
       if (mask(i)) then
          drho0(i) = fillvalue
          cycle
       end if
       zh = depth(i)
       ! masked in situ density anomaly
       drho0(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rho00
    end do
    !$omp  end parallel do
    !
  END SUBROUTINE drho0_mn4

  SUBROUTINE drho0_mn8( fillvalue,mask,theta0,S0,depth_km,depth,drho0,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_insitu_pot  ***
    !!
    !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
    !!      potential volumic mass (Kg/m3) from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter nn_eos.
    !!
    !! ** Method  :
    !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!
    !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
    !!          t = 40 deg celcius, s=40 psu
    !!
    !!
    !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
    !!                Brown and Campana, Mon. Weather Rev., 1978
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*8, INTENT(IN) ::   fillvalue,theta0,S0
    REAL*4, INTENT(IN) ::  depth(n)
    REAL*8, INTENT (IN) :: depth_km
    REAL*8, INTENT(OUT) ::  drho0(n)
    !f2py intent (in) theta0,S0,depth_km,depth,n
    !f2py intent (out) drho0

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, rho00       !    -         -
    !!----------------------------------------------------------------------

    zt = theta0
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

    !$omp parallel do private(zh)
    do i=1,n
       if (mask(i)) then
          drho0(i) = fillvalue
          cycle
       end if
       zh = depth(i)
       ! masked in situ density anomaly
       drho0(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rho00
    end do
    !$omp  end parallel do
    !
  END SUBROUTINE drho0_mn8

  SUBROUTINE sigma_n( theta,S,n,depth_km,rho)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_insitu_pot  ***
    !!
    !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
    !!      potential volumic mass (Kg/m3) from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter nn_eos.
    !!
    !! ** Method  :
    !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!
    !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
    !!          t = 40 deg celcius, s=40 psu
    !!
    !!
    !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
    !!                Brown and Campana, Mon. Weather Rev., 1978
    !!----------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n
    REAL*8, INTENT(IN) ::   theta(n),S(n)
    REAL*8, INTENT(IN) ::   depth_km
    REAL*8, INTENT(OUT) ::   rho(n)
    !f2py intent (in) n
    !f2py intent (in) theta,S
    !f2py intent (in) depth_km
    !f2py intent (out) rho
    INTEGER ::   i  ! loop counter
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zrau0r       !    -         -
    !!----------------------------------------------------------------------

    if (neos==0) then   !==  Jackett and McDougall (1994) formulation  ==!
       if (depth_km > 1.e-4) then
          !$omp parallel do private(zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze, zbw,zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
          do i=1,n
             zt = theta(i)
             zs = S(i)
             zh = depth_km*1000.d0                  ! depth
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
             ! masked in situ density anomaly
             rho(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0
          !$omp  end parallel do
          end do
       else   ! sigma0
          !$omp parallel do private(zt,zs,zsr,zr1,zr2,zr3,zr4,zrhop)
          do i=1,n
             zt = theta(i)
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
             rho(i) = zrhop - 1000.d0
          !$omp  end parallel do
          end do
       end if
       ! ************************* old NEMO way *************
       ! nemo uses prd = (rho-rau0)/rau0
       ! else if (neos==1) then  !==  Linear formulation = F( temperature )  ==!
       !    zs = 1.d0 + 0.0285d0
       !    ! prd =  0.0285 - rn_alpha * theta
       !    do i=1,n
       !       rho(i) =  rau0*(zs - rn_alpha *theta(i) ) - 1000.d0
       !    end do
       ! else if (neos==2) then  !==  Linear formulation = F( temperature , salinity )  ==!
       !    ! prd =  rn_beta * S - rn_alpha * theta
       !    do i=1,n
       !       rho(i) =  rau0*(1.d0 + rn_beta*S(i) - rn_alpha *theta(i) ) - 1000.d0
       !    end do
       ! **************************************************8

       ! here we want sigma_0, so reference wrt sigma(theta00,S00,0)= rho00 - 1000.
    else if (neos==1) then  !==  Linear formulation = F( temperature )  ==!
       ! prd =  0.0285 - rn_alpha * theta
       do i=1,n
          rho(i) =  rau0*rn_alpha*(theta00 - theta(i)) + rho000
       end do
    else if (neos==2) then  !==  Linear formulation = F( temperature , salinity )  ==!
       ! prd =  rn_beta * S - rn_alpha * theta
       do i=1,n
          rho(i) =  rau0*(rn_alpha*(theta00 - theta(i)) + rn_beta*(S(i) - S00)) + rho000
       end do
    else
       print *,'unknown eos'
       stop
    end if
  END SUBROUTINE sigma_n

  SUBROUTINE sigma_n8(fillvalue,mask,theta,S,n,depth_km,rho)
    !!----------------------------------------------------------------------
    !!                  ***  sigma_n8  ***
    !!
    !! ** Purpose :   Compute the  density (ratio rho/rau0) and the
    !!      potential volumic mass (Kg/m3)potential density anomaly at depth depth_km
    !!     from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter neos. A boolean mask (mask=.T. over land) must be supplied
    !!      and a fill value for those masked points.
    !!
    !!    All fields are assumed to be real*8
    !!
    !! ** Method  :
    !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!
    !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
    !!          t = 40 deg celcius, s=40 psu
    !!
    !!
    !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
    !!                Brown and Campana, Mon. Weather Rev., 1978
    !!----------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*8, INTENT(IN) ::   fillvalue
    REAL*8, INTENT(IN) ::   theta(n),S(n)
    REAL*8, INTENT(IN) ::   depth_km
    REAL*8, INTENT(OUT) ::   rho(n)
    !f2py intent (in) n
    !f2py intent (in) theta,S
    !f2py intent (in) depth_km
    !f2py intent (out) rho
    INTEGER ::   i  ! loop counter
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zrau0r       !    -         -
    !!----------------------------------------------------------------------

    if (neos==0) then   !==  Jackett and McDougall (1994) formulation  ==!
       if (depth_km > 1.e-4) then
          !$omp parallel do private(zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze, zbw,zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
          do i=1,n
             if (mask(i)) then
                rho(i) = fillvalue
                cycle
             end if
             zt = theta(i)
             zs = S(i)
             zh = depth_km*1000.d0                  ! depth
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
             ! masked in situ density anomaly
             rho(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0
          end do
          !$omp  end parallel do
       else   ! sigma0
          !$omp parallel do private(zt,zs,zsr,zr1,zr2,zr3,zr4,zrhop)
          do i=1,n
             if (mask(i)) then
                rho(i) = fillvalue
                cycle
             end if
             zt = theta(i)
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
             rho(i) = zrhop - 1000.d0
          end do
          !$omp  end parallel do
       end if
       ! nemo uses prd = (rho-rau0)/rau0
    else if (neos==1) then  !==  Linear formulation = F( temperature )  ==!
       zs = 1.d0 + 0.0285d0
       ! prd =  0.0285 - rn_alpha * theta
       do i=1,n
          rho(i) =  rau0*(zs - rn_alpha *theta(i) ) - 1000.d0
       end do
    else if (neos==2) then  !==  Linear formulation = F( temperature , salinity )  ==!
       ! prd =  rn_beta * S - rn_alpha * theta
       do i=1,n
          rho(i) =  rau0*(1.d0 + rn_beta*S(i) - rn_alpha *theta(i) ) - 1000.d0
       end do
    else
       print *,'unknown eos'
       stop
    end if
  END SUBROUTINE sigma_n8

  SUBROUTINE sigma_n4(fillvalue, mask, theta, S, n, depth_km, rho)
    !!----------------------------------------------------------------------
    !!                  ***  sigma_n4  ***
    !!
    !! ** Purpose :   Compute the  density (ratio rho/rau0) and the
    !!      potential volumic mass (Kg/m3)potential density anomaly at depth depth_km
    !!     from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter neos. A boolean mask (mask=.T. over land) must be supplied
    !!      and a fill value for those masked points.
    !!
    !!    All fields are assumed to be real*4
    !
    !! ** Method  :
    !!      nn_eos = 0 : Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!
    !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
    !!          t = 40 deg celcius, s=40 psu
    !!
    !!
    !! References :   Jackett and McDougall, J. Atmos. Ocean. Tech., 1994
    !!                Brown and Campana, Mon. Weather Rev., 1978
    !!----------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue
    REAL*4, INTENT(IN) ::   theta(n),S(n)
    REAL*8, INTENT(IN) ::   depth_km
    REAL*4, INTENT(OUT) ::   rho(n)
    !f2py intent (in) n
    !f2py intent (in) theta,S
    !f2py intent (in) depth_km
    !f2py intent (out) rho
    INTEGER ::   i  ! loop counter
    REAL*8 ::   zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw   ! temporary scalars
    REAL*8 ::   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0, zrau0r       !    -         -
    !!----------------------------------------------------------------------


    if (neos==0) then   !==  Jackett and McDougall (1994) formulation  ==!
       if (depth_km > 1.e-4) then
          !$omp parallel do private(zt,zs,zh,zsr,zr1,zr2,zr3,zr4,zrhop,ze, zbw,zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
          do i=1,n
             if (mask(i)) then
                rho(i) = fillvalue
                cycle
             end if
             zt = theta(i)
             zs = S(i)
             zh = depth_km*1000.d0                  ! depth
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
             ! masked in situ density anomaly
             rho(i) = zrhop / (  1.0d0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - 1000.d0
          end do
          !$omp  end parallel do
       else   ! sigma0
          !$omp parallel do private(zt,zs,zsr,zr1,zr2,zr3,zr4,zrhop)
          do i=1,n
             if (mask(i)) then
                rho(i) = fillvalue
                cycle
             end if
             zt = theta(i)
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
             rho(i) = zrhop - 1000.d0
          end do
          !$omp  end parallel do
       end if
       ! nemo uses prd = (rho-rau0)/rau0
    else if (neos==1) then  !==  Linear formulation = F( temperature )  ==!
       zs = 1.d0 + 0.0285d0
       ! prd =  0.0285 - rn_alpha * theta
       do i=1,n
          rho(i) =  rau0*(zs - rn_alpha *theta(i) ) - 1000.d0
       end do
    else if (neos==2) then  !==  Linear formulation = F( temperature , salinity )  ==!
       ! prd =  rn_beta * S - rn_alpha * theta
       do i=1,n
          rho(i) =  rau0*(1.d0 + rn_beta*S(i) - rn_alpha *theta(i) ) - 1000.d0
       end do
    else
       print *,'unknown eos'
       stop
    end if
  END SUBROUTINE sigma_n4


  SUBROUTINE alpha_beta4( fillvalue,mask,theta,S,depth,alpha,beta,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_bn2  ***
    !!
    !! ** Purpose :   Compute the in situ alpha and beta for masked real*4 inputs
    !!
    !! ** Method :
    !! References :   McDougall, J. Phys. Oceanogr., 17, 1950-1964, 1987.
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue,theta(n),S(n),depth(n)
    REAL*4, INTENT(OUT) ::  alpha(n), beta(n)
    !f2py intent (in) theta,S,depth,n
    !f2py intent (out) alpha,beta

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zalbet

    !$omp parallel do private (zt,zs,zh, zalbet)
    do i=1,n
       if (mask(i)) then
          alpha(i) = fillvalue
          beta(i) = fillvalue
          cycle
       end if

       zt = theta(i)
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
       beta(i) = ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
            &                               - 0.301985d-05 ) * zt      &
            &   +       0.785567d-03                                   &
            &   + (     0.515032d-08 * zs                              &
            &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
            &   + ( (   0.121551d-17 * zh                              &
            &         - 0.602281d-15 * zs                              &
            &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
            &                                + 0.408195d-10   * zs     &
            &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
            &                                - 0.121555d-07 ) * zh
       alpha(i) = zalbet*beta(i)
    END DO
    !$omp  end parallel do
  END SUBROUTINE alpha_beta4

  SUBROUTINE alpha_beta_n4( fillvalue,mask,theta,S,depth_km,alpha,beta,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_bn2  ***
    !!
    !! ** Purpose :   Compute the alpha and beta at depth depth_km  for masked real*4 inputs
    !!
    !! ** Method :
    !! References :   McDougall, J. Phys. Oceanogr., 17, 1950-1964, 1987.
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue,theta(n),S(n),depth_km
    REAL*4, INTENT(OUT) ::  alpha(n), beta(n)
    !f2py intent (in) theta,S,depth_km,n
    !f2py intent (out) alpha,beta

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zalbet

    !$omp parallel do private (zt,zs,zh, zalbet)
    do i=1,n
       if (mask(i)) then
          alpha(i) = fillvalue
          beta(i) = fillvalue
          cycle
       end if

       zt = theta(i)
       zs = S(i) - 35.d0
       zh = depth_km*1000.0d0                 ! depth
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
       beta(i) = ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
            &                               - 0.301985d-05 ) * zt      &
            &   +       0.785567d-03                                   &
            &   + (     0.515032d-08 * zs                              &
            &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
            &   + ( (   0.121551d-17 * zh                              &
            &         - 0.602281d-15 * zs                              &
            &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
            &                                + 0.408195d-10   * zs     &
            &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
            &                                - 0.121555d-07 ) * zh
       alpha(i) = zalbet*beta(i)
    END DO
    !$omp  end parallel do
  END SUBROUTINE alpha_beta_n4

  SUBROUTINE alpha_beta_04( fillvalue,mask,theta,S,alpha,beta,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_bn2  ***
    !!
    !! ** Purpose :   Compute surface alpha and beta for masked real*4 S & T
    !!
    !! ** Method :
    !! References :   McDougall, J. Phys. Oceanogr., 17, 1950-1964, 1987.
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue,theta(n),S(n)
    REAL*4, INTENT(OUT) ::  alpha(n), beta(n)
    !f2py intent (in) theta,S,n
    !f2py intent (out) alpha,beta

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zsr, zalbet

    do i=1,n
       if (mask(i)) then
          alpha(i) = fillvalue
          beta(i) = fillvalue
          cycle
       end if

       zt = theta(i)
       zs = S(i) - 35.d0
       zh = 0.d0                 ! depth
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
       beta(i) = ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
            &                               - 0.301985d-05 ) * zt      &
            &   +       0.785567d-03                                   &
            &   + (     0.515032d-08 * zs                              &
            &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
            &   + ( (   0.121551d-17 * zh                              &
            &         - 0.602281d-15 * zs                              &
            &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
            &                                + 0.408195d-10   * zs     &
            &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
            &                                - 0.121555d-07 ) * zh
       alpha(i) = zalbet*beta(i)
    END DO
  END SUBROUTINE alpha_beta_04

  SUBROUTINE alpha_beta8( fillvalue,mask,theta,S,depth,alpha,beta,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_bn2  ***
    !!
    !! ** Purpose :   Compute the in situ alpha and beta for masked real*8 inputs
    !!
    !! ** Method :
    !! References :   McDougall, J. Phys. Oceanogr., 17, 1950-1964, 1987.
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*8, INTENT(IN) ::   fillvalue,theta(n),S(n),depth(n)
    REAL*8, INTENT(OUT) ::  alpha(n), beta(n)
    !f2py intent (in) theta,S,depth,n
    !f2py intent (out) alpha,beta

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zsr, zalbet

    do i=1,n
       if (mask(i)) then
          alpha(i) = fillvalue
          beta(i) = fillvalue
          cycle
       end if

       zt = theta(i)
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
       beta(i) = ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
            &                               - 0.301985d-05 ) * zt      &
            &   +       0.785567d-03                                   &
            &   + (     0.515032d-08 * zs                              &
            &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
            &   + ( (   0.121551d-17 * zh                              &
            &         - 0.602281d-15 * zs                              &
            &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
            &                                + 0.408195d-10   * zs     &
            &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
            &                                - 0.121555d-07 ) * zh
       alpha(i) = zalbet*beta(i)
    END DO
  END SUBROUTINE alpha_beta8

  SUBROUTINE alpha_beta8_nomask(theta,S,depth_km,alpha,beta,n)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE eos_bn2  ***
    !!
    !! ** Purpose :   Compute the alpha and beta at depth depth_km for unmasked real*8 T & S
    !!
    !! ** Method :
    !! References :   McDougall, J. Phys. Oceanogr., 17, 1950-1964, 1987.
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    REAL*8, INTENT(IN) ::   theta(n),S(n),depth_km
    REAL*8, INTENT(OUT) ::  alpha(n), beta(n)
    !f2py intent (in) theta,S,depth,n
    !f2py intent (out) alpha,beta

    INTEGER*4 :: i
    REAL*8 ::   zt, zs, zh, zsr, zalbet

    do i=1,n
       zt = theta(i)
       zs = S(i) - 35.d0
       zh = depth_km*1000.0d0                 ! depth
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
       beta(i) = ( ( -0.415613d-09 * zt + 0.555579d-07 ) * zt    &   ! beta
            &                               - 0.301985d-05 ) * zt      &
            &   +       0.785567d-03                                   &
            &   + (     0.515032d-08 * zs                              &
            &         + 0.788212d-08 * zt - 0.356603d-06 ) * zs     &
            &   + ( (   0.121551d-17 * zh                              &
            &         - 0.602281d-15 * zs                              &
            &         - 0.175379d-14 * zt + 0.176621d-12 ) * zh     &
            &                                + 0.408195d-10   * zs     &
            &     + ( - 0.213127d-11 * zt + 0.192867d-09 ) * zt     &
            &                                - 0.121555d-07 ) * zh
       alpha(i) = zalbet*beta(i)
    END DO
  END SUBROUTINE alpha_beta8_nomask

end module eos
