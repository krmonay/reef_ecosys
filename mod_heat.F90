
!!!=== ver 2016/03/16   Copyright (c) 2012-2016 Takashi NAKAMURA  =====

#include "cppdefs.h"


!!!**** Heat MODULE ************************************

  MODULE mod_heat

    implicit none

    TYPE T_HEAT

      real(8), pointer :: tg(:,:) 

    END TYPE T_HEAT

    TYPE (T_HEAT), allocatable :: HEAT(:)

  CONTAINS

!!! **********************************************************************
!!!  set initial conditions for heat
!!! **********************************************************************

    SUBROUTINE initialize_heat(ng, Ngrids, LBi, UBi, LBj, UBj)

      implicit none
! input parameters
      integer, intent(in) :: ng, Ngrids, LBi, UBi, LBj, UBj
      integer i,j,n

      IF (ng.eq.1) allocate ( HEAT(Ngrids) )

      allocate( HEAT(ng)%tg(LBi:UBi,LBj:UBj)  )

!  Set initial conditions
      do j=LBj,UBj
        do i=LBi,UBi
          HEAT(ng)%tg(i,j)=29.1d0
        enddo
      enddo

      RETURN
    END SUBROUTINE initialize_heat


! **********************************************************************
!  Heat and mass balance of water column
! **********************************************************************

    SUBROUTINE heat_mass_balance  &
!          input parameters
     &            (ng, i, j       &   ! ng: nested grid number; i,j: position
     &            ,N              &   ! Number of vertical grid (following ROMS vertical grid)
     &            ,dt             &   ! Time step (sec)
     &            ,dz             &   ! dz(N): vertical grid size (m)
     &            ,ssradi         &   ! Surface shortwave radiation (W m-2)
     &            ,tair           &   ! air temperature (oC)
     &            ,pair           &   ! atm pressure (hPa)
     &            ,eair           &   ! vapor pressur (hPa)
     &            ,u10            &   ! wind speed (m s-1)
     &            ,fvol_pre       &   ! Precipitation volume flux (m s-1)
#ifdef LONGWAVE_OUT
     &            ,dw_lwradi      &   ! Downward longwave radiation (W m-2)
#endif
     &            ,Tmp            &   ! Tmp(N): Temperature (oC)
     &            ,Sal            &   ! Sal(N): Salinity (PSU)

!          output parameters
     &            ,dTmp_dt        &   ! dTmp_dt(N): dTmp/dt (K s-1)
     &             )
!
!-----------------------------------------------------------------------
!                                                                       
!                     rho point    Face                                 
!                       (i,j)                                           
!                    _____|______  _N    Surface                        
!                   /     |      /|                                     
!      z     Layer /___________ / |                                     
!                  |           |  |_N-1                                 
!     dz(N) {   N  |           | /|                                     
!                  |___________|/ : :                                   
!                  |           |  :      Water column                   
!               :  :           :  |_2                                   
!               :  :           : /|                                     
!               :  |___________|/ |                                     
!                  |           |  |_1                                   
!     dz(2) {   2  |           | /|                                     
!                  |___________|/ |                                     
!                  |           |  |_0    Bottom                         
!     dz(1) {   1  |           | /                                      
!                  |___________|/                                       
!                                                                       
!                                                                       
!      A vertical section of the ecosys grid showing water column.      
!-----------------------------------------------------------------------
!
      implicit none

! input parameters
      integer, intent(in) :: ng, i, j        
      integer, intent(in) :: N          
      real(8), intent(in) :: dt         
      real(8), intent(in) :: dz(N)      
      real(8), intent(in) :: ssradi     
      real(8), intent(in) :: tair    ! air temperature (oC)
      real(8), intent(in) :: pair    ! atm pressure (hPa)
      real(8), intent(in) :: eair    ! vapor pressur (hPa)
      real(8), intent(in) :: u10        
      real(8), intent(in) :: fvol_pre  ! Precipitation volume flux (m s-1)
#ifdef LONGWAVE_OUT
      real(8), intent(in) :: dw_lwradi ! Downward longwave radiation (W m-2)
#endif

      real(8), intent(in) :: Tmp(N)        
      real(8), intent(in) :: Sal(N)        
! output parameters
      real(8), intent(out) :: dTmp_dt(N)     

      real(8), parameter :: Cpw   = 3.99d6   ! Seawater heat capacity 3.99d6 J m-3 K-1
      real(8), parameter :: Eevap = 2.44d6   ! Heat of evaporation 2.44d6 J m-3

      real(8) :: lradi, sens
      real(8) :: fvol_evap
      real(8) :: radi, radir
      real(8) :: tr, bsed,cbsed
        !Tr: transmission, 
        !bsed: turbidity parameter; offshere water = 1, 
        !                           coast water    = 2-3, 
        !                           cloudy water   = 5-9
        !cbsed:adjusting parameter for bsed
        
      real(8) :: Ftbot    ! Bottom heat flux (K m s-1)
      real(8) :: Diff     ! Vertical eddy diffusion coefficient (m2 s-1)
            
      integer :: k


      cbsed=10. !!!“K“–!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Check after Sed changed

!-----------------------------------------------------------------------
! Initialize all time difference values (d_XXX)
!-----------------------------------------------------------------------
      
      DO k=1,N
        dTmp_dt(k) = 0.0d0
!        dSal_dt(k) = 0.0d0
      END DO
      
!-----------------------------------------------------------------------
! Sea suface heat flux
!-----------------------------------------------------------------------

! Long wave radiation (W m-2)
#ifdef LONGWAVE_OUT
      lradi = up_long_wave_radi(Tmp(N))-dw_lwradi
#else
      lradi = net_long_wave_radi(Tmp(N))
#endif

! Sensible heat (W m-2)
      sens = sensheat(Tmp(N), tair, pair, eair, u10)

! Evaporation volume flux (m s-1)
      fvol_evap=evaporation(Tmp(N), Sal(N), tair, pair, eair, u10) 
      


      dTmp_dt(N) = dTmp_dt(N) + ( (-lradi -sens -fvol_evap*Eevap)/Cpw     &
     &                             +fvol_pre*(tair-Tmp(N))                &
     &                          )  /dz(N)

!-----------------------------------------------------------------------
! Water column heat flux
!-----------------------------------------------------------------------
      
! Short wave radiation

! Ref. Kondo et al.(1979) J. Phys. Oceanogr. 9, 360-372.
!      Kondo (1994) book

      radi=(1.0d0-0.07)*ssradi
      
      do k=N,1,-1
        bsed=2.0e0!+cbsed*Ci(nSed,i,j,1)
        tr= 0.14*exp(-500.0d0 *dz(k))         &
     &     +0.23*exp( -12.0d0 *dz(k))         &
     &     +0.14*exp( -25.0d0 *dz(k))         &
     &     +0.13*exp(  -0.4d0 *dz(k))         &
     &     +0.20*exp(  -0.1d0 *dz(k)*bsed)    &
     &     +0.16*exp(  -0.04d0*dz(k)*bsed)

        dTmp_dt(k) = dTmp_dt(k) + radi*(1.0d0-tr)/Cpw/dz(k)

        radi=radi*tr

      end do

!     Refrection from bottom
      radir=(1.0d0-0.5)*radi !TN:Albedo of bottom sediment = 0.7            
      
!-----------------------------------------------------------------------
! Bottom heat flux
!-----------------------------------------------------------------------
      
      Diff = 1.0d-3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      Ftbot= Diff *(HEAT(ng)%tg(i,j)-Tmp(1))/(dz(1)/2)
      dTmp_dt(1) = dTmp_dt(1) + Ftbot/dz(1)

      HEAT(ng)%tg(i,j)=HEAT(ng)%tg(i,j)                                &
     &         +(-1.0d0/(0.1*86400.)*(HEAT(ng)%tg(i,j)-29.0)  &  ! 0.1 day relaxation time
     &            -Ftbot                             &
     &            +(radi-radir)/3.0d6/0.1            &  ! ‰·“x‚Ìã‚ª‚éŠâ”Õ‚ÌŒú‚³‚ð0.2m‚Æ‰¼’è
     &           )*dt                                   ! Šâ”Õ‚Ì”M—e—Ê3.0e6 J m-3 K-1‚Æ‰¼’è

!-----------------------------------------------------------------------
! Water column heat flux
!-----------------------------------------------------------------------

      do k=1,N
        bsed=2.0d0!+cbsed*Ci(nSed,i,j,1)
        tr= 0.14*exp(-500.0d0 *dz(k))         &
     &     +0.23*exp( -12.0d0 *dz(k))         &
     &     +0.14*exp( -25.0d0 *dz(k))         &
     &     +0.13*exp(  -0.4d0 *dz(k))         &
     &     +0.20*exp(  -0.1d0 *dz(k)*bsed)    &
     &     +0.16*exp(  -0.04d0*dz(k)*bsed)
     
        dTmp_dt(k) = dTmp_dt(k) + radir*(1.0d0-tr)/Cpw/dz(k)

        radir=radir*tr
        
      end do
          
      RETURN

    END SUBROUTINE heat_mass_balance


!! **********************************************************************
!!  Caliculate heat flux between water column and soil. 
!! **********************************************************************
!
!    subroutine calgsoil
!!
!      implicit none
!
!
!      dzzzz = 0.03                        !  [m]
!      dzzzz2 = dzzzz*dzzzz
!      ag = 5.9d-7                      !  [m2/s]
!      cgrg =3.0d+6                     !  [J/m3/K]
!
!!$OMP parallel
!!$OMP do
!      do j=1,jm
!        do i=1,im
!
!          Tg(i,j,1)=t(i,j,kbm1)
!
!          do loop=1,20
!
!            do k=2,kmaxG-1
!              Tgn(i,j,k)=(Tgb(i,j,k)
!     &            +dti*ag*(Tg(i,j,k+1)+Tg(i,j,k-1))/dzzzz2)         &
!     &                                   /(1.0+2.*ag*dti/dzzzz2)
!            enddo
!
!            error = 0.e0
!            do k=2,kmaxG-1
!              error = max( error,abs(Tgn(i,j,k)-Tg(i,j,k)) )
!            enddo
!
!            if(error.le.1.0d-6) goto 1000
!
!            do k=2,kmaxG-1
!              Tg(i,j,k) = Tgn(i,j,k)
!            enddo
!
!!:      b.c for Tg
!
!            Tg(i,j,1)    = t(i,j,kbm1)
!            Tg(i,j,kmaxG)=Tg(i,j,kmaxG-1)
!
!          enddo
!
! 1000   continue
!
!          Gsoil(i,j) = 0.e0
!          Gsoil(i,j) = Gsoil(i,j)         &
!     &                 + cgrg*(Tg(i,j,2)-Tgb(i,j,2))/dti*dzzzz
!          do k=2,kmaxG
!            Gsoil(i,j) = Gsoil(i,j)         &
!     &                 + cgrg*(Tg(i,j,k)-Tgb(i,j,k))/dti*dzzzz   ! [J/m2/s]
!          enddo
!
!          Gsoil(i,j)=Gsoil(i,j)/rhoref         &                           ! [m*K/s]
!     &             /(3.958-0.0523*s(i,j,kbm1)+0.000837*t(i,j,kbm1))/1000  !Bromley, L.A. etal,(1967)
!
!          if(iint.gt.int(3600.e0/dti)) then
!            Gcoef = 1.0
!          else
!            Gcoef = 0.5*(1.-cos(pi*(iint*dti)/3600.0))
!          endif
!
!          Gsoil(i,j) = Gcoef * Gsoil(i,j)
!
!        enddo
!      enddo
!!$OMP end do
!
!!$OMP do
!      do j=1,jm
!        do i=1,im
!          do k=1,kmaxG
!            Tgb(i,j,k)=Tg(i,j,k)
!         enddo
!        enddo
!      enddo
!!$OMP end do
!!$OMP end parallel
!
!      return
!    end subroutine calgsoil


! **********************************************************************
!  Evaporation volume flux (m s-1)
! **********************************************************************

    real(8) function evaporation(tsurf, Ssurf, tair, pair, eair, u10) ! Evaporation volume flux (m s-1)

      implicit none
!                          input parameters
      real(8), intent(in) :: tsurf   ! Sea surface temperature (oC)
      real(8), intent(in) :: Ssurf   ! Sea surface salinity (psu)
      real(8), intent(in) :: tair    ! air temperature (oC)
      real(8), intent(in) :: pair    ! atm pressure (hPa)
      real(8), intent(in) :: eair    ! vapor pressur (hPa)
      real(8), intent(in) :: u10     ! wind speed (m s-1)
!      real(8), intent(in) :: Hum     !relative humidity (%)

      real(8) rhoatm, esw, spHumw, spHuma   

!      tair=27.0
!      pair=1013.0
!      Hum=50.d0   ! Humidity

!      eair=sat_vapor_press(tair)*Hum/100.d0    ! vapor pressur (hPa)
      rhoatm=dens_air(tair,pair,eair)                   !density of air (kg m-3)
      spHuma=specific_humidity(eair,pair)               !air specific humidity (kg kg-1)
      esw=(1.d0-5.37d-4*Ssurf)*sat_vapor_press(tsurf)   !vapor pressur on sea surface(hPa)
      spHumw=specific_humidity(esw,pair)                !specific humidity at sea surface (kg kg-1)

      evaporation=rhoatm*1.2d-3*u10*(spHumw-spHuma)     ! Evaporation volume flux (m s-1)

      return
    end function evaporation


! **********************************************************************
!  Sensible heat (W m-2)  [sea surface to air is postive]
! **********************************************************************

    real(8) function sensheat(tsurf, tair, pair, eair, u10) !sensible heat (W m-2)

      implicit none
      
      real(8), intent(in) :: tsurf   ! Sea surface temperature (oC)
      real(8), intent(in) :: tair    ! air temperature (oC)
      real(8), intent(in) :: pair    ! atm pressure (hPa)
      real(8), intent(in) :: eair    ! vapor pressur (hPa)
      real(8), intent(in) :: u10     ! wind speed (m s-1)
!      real(8), intent(in) :: Hum     !relative humidity (%)

      real(8) rhoatm   !density of air (kg m-3)

!      tair=27.0
!      pair=1013.0
!      Hum=50.e0   ! Humidity
!      eair=sat_vapor_press(tair)*Hum/100.e0    ! vapor pressur (hPa)

      rhoatm=dens_air(tair,pair,eair)                 !density of air (kg m-3)
      sensheat=1005.d0*rhoatm*1.2d-3*u10*(tsurf-tair) !sensible heat (W m-2)
              !1005 J kg-1 K-1 -> ‹ó‹C‚Ì’èˆ³”ä”M
      return
    end function sensheat

! **********************************************************************
!  Long wave radiation (W m-2) from sea to atmosphere
! **********************************************************************

    real(8) function up_long_wave_radi(tsurf) ! Long wave radiation (W m-2)

      implicit none
!                          input parameters
      real(8), intent(in) :: tsurf   ! Sea surface temperature (oC)
      
      up_long_wave_radi = 0.96*(5.670d-8*(tsurf+273.15)**4)
      
      return
    end function up_long_wave_radi

! **********************************************************************
!  Net long wave radiation (W m-2)
! **********************************************************************

    real(8) function net_long_wave_radi(tsurf,tair,eair) ! Long wave radiation (W m-2)

      implicit none
!                          input parameters
      real(8), intent(in) :: tsurf   ! Sea surface temperature (oC)
      real(8), intent(in) :: tair    ! air temperature (oC)
      real(8), intent(in) :: eair    ! vapor pressur (hPa)
      real(8) :: Ta,Ts
      
      Ta = tair + 273.15d0
      Ts = tsurf+ 273.15d0
      
      net_long_wave_radi = 0.97 * 5.670d-8 * Ta**3.0d0 * (              &
     &     Ta*(0.39d0-0.05d0*SQRT(eair))*(1.0d0-x*n*n)                  &
     &   + 4.0d0*Ta*(Ts-Ta) )
      
      return
    end function net_long_wave_radi


! **********************************************************************
!  density of air (kg m-3)
! **********************************************************************

    real(8) function dens_air(tair,pair,eair)!  density of air (kg m-3)

      implicit none

      real(8), intent(in) :: tair    ! air temperature (oC)
      real(8), intent(in) :: pair    ! atm pressure (hPa)
      real(8), intent(in) :: eair    ! vapor pressur (hPa)


      dens_air=1.293*273.15/(273.15+tair)*(pair/1013.25)  &
     &          *(1-0.378*eair/pair)

      return
    end function dens_air


! **********************************************************************
!  Saturation vapor pressure (hPa)
! **********************************************************************

    real(8) function sat_vapor_press(tair)

      implicit none

      real(8), intent(in) :: tair    ! air temperature (oC)

      sat_vapor_press=6.1078*10.**(7.5*tair/(237.3+tair))

      return
    end function sat_vapor_press
      
      
! **********************************************************************
!  Specific humidity (kg kg-1)
! **********************************************************************

    real(8) function specific_humidity(eair,pair)

      implicit none
      
      real(8), intent(in) :: pair    ! atm pressure (hPa)
      real(8), intent(in) :: eair    ! vapor pressur (hPa)
      
      real(8) e_per_p

      e_per_p=eair/pair
      specific_humidity=0.622*e_per_p/(1.-0.378*e_per_p)

      return
    end function specific_humidity
    
    
  END MODULE mod_heat


