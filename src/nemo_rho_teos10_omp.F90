MODULE eos_teos10
  !!==============================================================================
  !!                       ***  MODULE  eos_teos  ***
  !! Equation Of Seawater for TEOS-10 : potential density, expansion coefficients
  !!==============================================================================
  !!             -   ! 2013-04  (F. Roquet, G. Madec)  add eos_rab, change bn2 computation and reorganize the module
  !!             -   ! 2014-09  (F. Roquet)  add TEOS-10, S-EOS, and modify EOS-80
  !!----------------------------------------------------------------------
  IMPLICIT NONE
  REAL*8, PARAMETER ::     rdeltaS = 32.d0
  REAL*8, PARAMETER ::     r1_S0  = 0.875d0/35.16504d0
  REAL*8, PARAMETER ::     r1_T0  = 1.d0/40.d0
  REAL*8, PARAMETER ::     r1_Z0  = 1.d-4
  REAL*8, PARAMETER ::     r1_rau0 = 1.d0/1026.d0

  REAL*8, PARAMETER ::     R00 = 4.6494977072d+01
  REAL*8, PARAMETER ::     R01 = -5.2099962525d+00
  REAL*8, PARAMETER ::     R02 = 2.2601900708d-01
  REAL*8, PARAMETER ::     R03 = 6.4326772569d-02
  REAL*8, PARAMETER ::     R04 = 1.5616995503d-02
  REAL*8, PARAMETER ::     R05 = -1.7243708991d-03
  
  REAL*8, PARAMETER ::     EOS000 = 8.0189615746d+02
  REAL*8, PARAMETER ::     EOS100 = 8.6672408165d+02
  REAL*8, PARAMETER ::     EOS200 = -1.7864682637d+03
  REAL*8, PARAMETER ::     EOS300 = 2.0375295546d+03
  REAL*8, PARAMETER ::     EOS400 = -1.2849161071d+03
  REAL*8, PARAMETER ::     EOS500 = 4.3227585684d+02
  REAL*8, PARAMETER ::     EOS600 = -6.0579916612d+01
  REAL*8, PARAMETER ::     EOS010 = 2.6010145068d+01
  REAL*8, PARAMETER ::     EOS110 = -6.5281885265d+01
  REAL*8, PARAMETER ::     EOS210 = 8.1770425108d+01
  REAL*8, PARAMETER ::     EOS310 = -5.6888046321d+01
  REAL*8, PARAMETER ::     EOS410 = 1.7681814114d+01
  REAL*8, PARAMETER ::     EOS510 = -1.9193502195d0
  REAL*8, PARAMETER ::     EOS020 = -3.7074170417d+01
  REAL*8, PARAMETER ::     EOS120 = 6.1548258127d+01
  REAL*8, PARAMETER ::     EOS220 = -6.0362551501d+01
  REAL*8, PARAMETER ::     EOS320 = 2.9130021253d+01
  REAL*8, PARAMETER ::     EOS420 = -5.4723692739d0
  REAL*8, PARAMETER ::     EOS030 = 2.1661789529d+01
  REAL*8, PARAMETER ::     EOS130 = -3.3449108469d+01
  REAL*8, PARAMETER ::     EOS230 = 1.9717078466d+01
  REAL*8, PARAMETER ::     EOS330 = -3.1742946532d0
  REAL*8, PARAMETER ::     EOS040 = -8.3627885467d0
  REAL*8, PARAMETER ::     EOS140 = 1.1311538584d+01
  REAL*8, PARAMETER ::     EOS240 = -5.3563304045d0
  REAL*8, PARAMETER ::     EOS050 = 5.4048723791d-01
  REAL*8, PARAMETER ::     EOS150 = 4.8169980163d-01
  REAL*8, PARAMETER ::     EOS060 = -1.9083568888d-01
  REAL*8, PARAMETER ::     EOS001 = 1.9681925209d+01
  REAL*8, PARAMETER ::     EOS101 = -4.2549998214d+01
  REAL*8, PARAMETER ::     EOS201 = 5.0774768218d+01
  REAL*8, PARAMETER ::     EOS301 = -3.0938076334d+01
  REAL*8, PARAMETER ::     EOS401 = 6.6051753097d0
  REAL*8, PARAMETER ::     EOS011 = -1.3336301113d+01
  REAL*8, PARAMETER ::     EOS111 = -4.4870114575d0
  REAL*8, PARAMETER ::     EOS211 = 5.0042598061d0
  REAL*8, PARAMETER ::     EOS311 = -6.5399043664d-01
  REAL*8, PARAMETER ::     EOS021 = 6.7080479603d0
  REAL*8, PARAMETER ::     EOS121 = 3.5063081279d0
  REAL*8, PARAMETER ::     EOS221 = -1.8795372996d0
  REAL*8, PARAMETER ::     EOS031 = -2.4649669534d0
  REAL*8, PARAMETER ::     EOS131 = -5.5077101279d-01
  REAL*8, PARAMETER ::     EOS041 = 5.5927935970d-01
  REAL*8, PARAMETER ::     EOS002 = 2.0660924175d0
  REAL*8, PARAMETER ::     EOS102 = -4.9527603989d0
  REAL*8, PARAMETER ::     EOS202 = 2.5019633244d0
  REAL*8, PARAMETER ::     EOS012 = 2.0564311499d0
  REAL*8, PARAMETER ::     EOS112 = -2.1311365518d-01
  REAL*8, PARAMETER ::     EOS022 = -1.2419983026d0
  REAL*8, PARAMETER ::     EOS003 = -2.3342758797d-02
  REAL*8, PARAMETER ::     EOS103 = -1.8507636718d-02
  REAL*8, PARAMETER ::     EOS013 = 3.7969820455d-01

  REAL*8, PARAMETER :: ALP000 = -6.5025362670d-01
  REAL*8, PARAMETER :: ALP100 = 1.6320471316d0
  REAL*8, PARAMETER :: ALP200 = -2.0442606277d0
  REAL*8, PARAMETER :: ALP300 = 1.4222011580d0
  REAL*8, PARAMETER :: ALP400 = -4.4204535284d-01
  REAL*8, PARAMETER :: ALP500 = 4.7983755487d-02
  REAL*8, PARAMETER :: ALP010 = 1.8537085209d0
  REAL*8, PARAMETER :: ALP110 = -3.0774129064d0
  REAL*8, PARAMETER :: ALP210 = 3.0181275751d0
  REAL*8, PARAMETER :: ALP310 = -1.4565010626d0
  REAL*8, PARAMETER :: ALP410 = 2.7361846370d-01
  REAL*8, PARAMETER :: ALP020 = -1.6246342147d0
  REAL*8, PARAMETER :: ALP120 = 2.5086831352d0
  REAL*8, PARAMETER :: ALP220 = -1.4787808849d0
  REAL*8, PARAMETER :: ALP320 = 2.3807209899d-01
  REAL*8, PARAMETER :: ALP030 = 8.3627885467d-01
  REAL*8, PARAMETER :: ALP130 = -1.1311538584d0
  REAL*8, PARAMETER :: ALP230 = 5.3563304045d-01
  REAL*8, PARAMETER :: ALP040 = -6.7560904739d-02
  REAL*8, PARAMETER :: ALP140 = -6.0212475204d-02
  REAL*8, PARAMETER :: ALP050 = 2.8625353333d-02
  REAL*8, PARAMETER :: ALP001 = 3.3340752782d-01
  REAL*8, PARAMETER :: ALP101 = 1.1217528644d-01
  REAL*8, PARAMETER :: ALP201 = -1.2510649515d-01
  REAL*8, PARAMETER :: ALP301 = 1.6349760916d-02
  REAL*8, PARAMETER :: ALP011 = -3.3540239802d-01
  REAL*8, PARAMETER :: ALP111 = -1.7531540640d-01
  REAL*8, PARAMETER :: ALP211 = 9.3976864981d-02
  REAL*8, PARAMETER :: ALP021 = 1.8487252150d-01
  REAL*8, PARAMETER :: ALP121 = 4.1307825959d-02
  REAL*8, PARAMETER :: ALP031 = -5.5927935970d-02
  REAL*8, PARAMETER :: ALP002 = -5.1410778748d-02
  REAL*8, PARAMETER :: ALP102 = 5.3278413794d-03
  REAL*8, PARAMETER :: ALP012 = 6.2099915132d-02
  REAL*8, PARAMETER :: ALP003 = -9.4924551138d-03
  !
  REAL*8, PARAMETER :: BET000 = 1.0783203594d+01
  REAL*8, PARAMETER :: BET100 = -4.4452095908d+01
  REAL*8, PARAMETER :: BET200 = 7.6048755820d+01
  REAL*8, PARAMETER :: BET300 = -6.3944280668d+01
  REAL*8, PARAMETER :: BET400 = 2.6890441098d+01
  REAL*8, PARAMETER :: BET500 = -4.5221697773d0
  REAL*8, PARAMETER :: BET010 = -8.1219372432d-01
  REAL*8, PARAMETER :: BET110 = 2.0346663041d0
  REAL*8, PARAMETER :: BET210 = -2.1232895170d0
  REAL*8, PARAMETER :: BET310 = 8.7994140485d-01
  REAL*8, PARAMETER :: BET410 = -1.1939638360d-01
  REAL*8, PARAMETER :: BET020 = 7.6574242289d-01
  REAL*8, PARAMETER :: BET120 = -1.5019813020d0
  REAL*8, PARAMETER :: BET220 = 1.0872489522d0
  REAL*8, PARAMETER :: BET320 = -2.7233429080d-01
  REAL*8, PARAMETER :: BET030 = -4.1615152308d-01
  REAL*8, PARAMETER :: BET130 = 4.9061350869d-01
  REAL*8, PARAMETER :: BET230 = -1.1847737788d-01
  REAL*8, PARAMETER :: BET040 = 1.4073062708d-01
  REAL*8, PARAMETER :: BET140 = -1.3327978879d-01
  REAL*8, PARAMETER :: BET050 = 5.9929880134d-03
  REAL*8, PARAMETER :: BET001 = -5.2937873009d-01
  REAL*8, PARAMETER :: BET101 = 1.2634116779d0
  REAL*8, PARAMETER :: BET201 = -1.1547328025d0
  REAL*8, PARAMETER :: BET301 = 3.2870876279d-01
  REAL*8, PARAMETER :: BET011 = -5.5824407214d-02
  REAL*8, PARAMETER :: BET111 = 1.2451933313d-01
  REAL*8, PARAMETER :: BET211 = -2.4409539932d-02
  REAL*8, PARAMETER :: BET021 = 4.3623149752d-02
  REAL*8, PARAMETER :: BET121 = -4.6767901790d-02
  REAL*8, PARAMETER :: BET031 = -6.8523260060d-03
  REAL*8, PARAMETER :: BET002 = -6.1618945251d-02
  REAL*8, PARAMETER :: BET102 = 6.2255521644d-02
  REAL*8, PARAMETER :: BET012 = -2.6514181169d-03
  REAL*8, PARAMETER :: BET003 = -2.3025968587d-04

CONTAINS

  SUBROUTINE set_eos_threads(nthreads)
    INTEGER*4, INTENT(IN)  :: nthreads
    CALL omp_set_num_threads(nthreads)
  END SUBROUTINE set_eos_threads

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
    !!
    !!----------------------------------------------------------------------
    REAL*8, INTENT(IN) ::   depth_km
    REAL*8, INTENT(OUT) ::  r0
    !f2py intent (in) depth_km
    !f2py intent (out) r0
    !
    REAL*8 ::   zh              ! local scalars
    !!----------------------------------------------------------------------
    !
    zh = depth_km * 1.d3 * r1_Z0
    r0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
  END SUBROUTINE get_r0
  
  SUBROUTINE sigma_n4( fillvalue, mask, CT, SA, n, depth_km, rho)
    !!----------------------------------------------------------------------
    !!                  ***  ROUTINE sigma_n4  ***
    !!
    !! ** Purpose :   Compute the potential density - 1000 kg/m^3 at depth_km depth
    !!      from conservative temperature and absolute salinity using TEOS-10
    !!
    !!   Note that this routine returns
    !!   rho(S_A,CT,z) = r0(z) + r(S_A,CT,z)
    !!  with
    !!     r0(z) = rho(35.16504g/kg,4 deg,z) - rho(35.16504g/kg,4 deg,0 m)
    !!
    !!----------------------------------------------------------------------
    INTEGER, INTENT(IN) :: n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue
    REAL*4, INTENT(IN) ::   CT(n), SA(n)
    REAL*8, INTENT(IN) ::   depth_km
    REAL*4, INTENT(OUT) ::   rho(n)
    !f2py intent (in) fillvalue,mask
    !f2py intent (in) CT,SA
    !f2py intent (in) n
    !f2py intent (in) depth_km
    !f2py intent (out) rho
    !
    !
    INTEGER  ::   i                ! dummy loop indices
    REAL*8 ::   zt , zh , zs, r0              ! local scalars
    REAL*8 ::   zn , zn0, zn1, zn2, zn3   !   -      -
    !!----------------------------------------------------------------------
    !
    zh = depth_km * 1.d3 * r1_Z0
    !!----------------------------------------------------------------------
    !!    get reference profile r0
    r0 = (((((R05 * zh + R04) * zh + R03 ) * zh + R02 ) * zh + R01) * zh + R00) * zh
    !!----------------------------------------------------------------------

    !!----------------------------------------------------------------------
    !!    now get anomaly part r(S_A,CT,z)
    !!----------------------------------------------------------------------
    !$omp parallel do private(zt, zs, zn, zn0, zn1, zn2, zn3)
    DO i=1,n   ! vector opt.
       if (mask(i)) then
          rho(i) = fillvalue
          cycle
       end if
       !
       zt  = CT (i) * r1_T0                           ! temperature
       zs  = SQRT( ABS( SA(i) + rdeltaS ) * r1_S0 )   ! square root salinity
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
       !
    !!----------------------------------------------------------------------
    !!  Add reference profile r0 to anomaly part & subtract 1000.
      rho(i) = REAL(zn - 1000.d0 + r0, KIND=4)
    !!----------------------------------------------------------------------
       ! prd(ji,jj) = zn * r1_rau0 - 1.               ! unmasked in situ density anomaly
    END DO
    !$omp  end parallel do
  END SUBROUTINE sigma_n4



  SUBROUTINE alpha_beta_n4( fillvalue,mask,CT,SA,depth_km,alpha,beta,n )
    !!----------------------------------------------------------------------
    !!                 ***  ROUTINE rab_2d  ***
    !!
    !! ** Purpose :  Compute the (boussinesq, assuming rho0=1026.) alpha and beta
    !!
    !!----------------------------------------------------------------------
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue, CT(n), SA(n), depth_km
    REAL*4, INTENT(OUT) ::  alpha(n), beta(n)
    !f2py intent (in) fillvalue,mask
    !f2py intent (in) CT,SA,depth_km,n
    !f2py intent (out) alpha,beta
    !
    INTEGER  ::   i                       ! dummy loop indices
    REAL*8 ::   zt , zh , zs              ! local scalars
    REAL*8 ::   zn , zn0, zn1, zn2, zn3   !   -      -
    !!----------------------------------------------------------------------
    !
    !
    !$omp parallel do private (zt,zs,zh, zn, zn0, zn1, zn2, zn3)
    do i=1,n
       if (mask(i)) then
          alpha(i) = fillvalue
          beta(i) = fillvalue
          cycle
       end if
       !
       zh  = depth_km * r1_Z0                                  ! depth
       zt  = CT(i) * r1_T0                           ! temperature
       zs  = SQRT( ABS( SA(i) + rdeltaS ) * r1_S0 )   ! square root salinity
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
       alpha(i) = REAL( zn * r1_rau0, KIND=4)
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
       beta(i) = REAL((zn / zs * r1_rau0), KIND=4)
       !
    END DO
    !$omp  end parallel do
  END SUBROUTINE alpha_beta_n4




  SUBROUTINE pt_from_ct(fillvalue, mask, CT, SA, pt, n)
    !!----------------------------------------------------------------------
    !!                 ***  ROUTINE eos_pt_from_ct  ***
    !!
    !! ** Purpose :   Compute pot.temp. from cons. temp. [Celsius]
    !!
    !! ** Method  :   rational approximation (5/3th order) of TEOS-10 algorithm
    !!       checkvalue: pt=20.02391895 Celsius for sa=35.7g/kg, ct=20degC
    !!
    !! Reference  :   TEOS-10, UNESCO
    !!                Rational approximation to TEOS10 algorithm (rms error on WOA13 values: 4.0e-5 degC)
    !!----------------------------------------------------------------------
    REAL*8, PARAMETER ::     zdeltaS = 5.d0
    INTEGER*4, INTENT(IN) ::   n
    LOGICAL(KIND=1), INTENT(IN) :: mask(n)
    REAL*4, INTENT(IN) ::   fillvalue
    REAL*4, INTENT(in) ::   CT(n)  ! Cons. Temp   [Celsius]
    REAL*4, INTENT(in) ::   SA(n)  ! absolute salinity [g/kg]
    REAL*4, INTENT(OUT) ::  pt(n)   ! potential temperature [Celsius]
    !f2py intent (in) fillvalue,mask
    !f2py intent (in) CT,SA,n
    !f2py intent (out) pt
    !
    INTEGER  ::   i                    ! dummy loop indices
    REAL*8 ::   zt , zs              ! local scalars
    REAL*8 ::   zn , zd              ! local scalars
    !!----------------------------------------------------------------------

    !$omp parallel do private (zt,zs,zn,zd)
    do i=1,n
       if (mask(i)) then
          pt(i) = fillvalue
          cycle
       end if
       !
       zt  = CT(i) * r1_T0
       zs  = SQRT( ABS( SA(i) + zdeltaS ) * r1_S0 )
       !
       zn = ((((-2.1385727895e-01*zt   &
            &   - 2.7674419971e-01*zs+1.0728094330)*zt   &
            &   + (2.6366564313*zs+3.3546960647)*zs-7.8012209473)*zt   &
            &   + ((1.8835586562*zs+7.3949191679)*zs-3.3937395875)*zs-5.6414948432)*zt   &
            &   + (((3.5737370589*zs-1.5512427389e+01)*zs+2.4625741105e+01)*zs   &
            &      +1.9912291000e+01)*zs-3.2191146312e+01)*zt   &
            &   + ((((5.7153204649e-01*zs-3.0943149543)*zs+9.3052495181)*zs   &
            &      -9.4528934807)*zs+3.1066408996)*zs-4.3504021262e-01
       !
       zd = (2.0035003456*zt   &
            &   -3.4570358592e-01*zs+5.6471810638)*zt   &
            &   + (1.5393993508*zs-6.9394762624)*zs+1.2750522650e+01
       !
       pt(i) = REAL(zt / r1_T0 + zn / zd, KIND=4)
       !
    END DO
    !$omp  end parallel do
  END SUBROUTINE pt_from_ct
END MODULE eos_teos10
