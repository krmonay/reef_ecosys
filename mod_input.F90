
!!!=== ver 2017/03/07   Copyright (c) 2012-2017 Takashi NAKAMURA  =====

!--------------------------------------------------------------------------------
!
!              Input module
!
!--------------------------------------------------------------------------------

#include "cppdefs.h"

  module mod_input
  
  contains


! **********************************************************************
!  Read timeseries data file
! **********************************************************************

    subroutine read_timeseies
    
      USE mod_param
      
      implicit none
      
      integer, parameter :: datamax = 1000
      real(8) :: fdata(4,datamax)
      
      integer i,j
!
! ----- READ Tide data (m) ---------------------------------------------

      open(77,file='./input/tide2009summer.dat')
      dt_tide=1.d0  !1.0 hour interval
      
      nm_tide=1      
      do i=1,datamax
        read(77,*,end=1000) fdata(1,i) !tide(i)
        nm_tide=nm_tide+1
      enddo
 1000 close(77)

      allocate( tide(nm_tide) )

      do i=1,nm_tide
        tide(i)=fdata(1,i)
      enddo

      etide=0.0e0

! --------- READ Wind data (m/s)------------------------------------------

      open(77,file='./input/wind2009summer.dat')
      dt_wind=1.d0  !1.0 hour interval
      
      nm_wind=1
      do i=1,datamax
        read(77,*,end=1100) fdata(1,i),fdata(2,i) !windu(i), windv(i)
        nm_wind=nm_wind+1
      enddo
 1100 close(77)
 
      allocate( windu(nm_wind), windv(nm_wind) )

      do i=1,nm_wind
        windu(i)=fdata(1,i)
        windv(i)=fdata(2,i)
      enddo


! ----- READ meterological dataset --------------------------------

      open(77,file='./input/air2009summer.dat')
      dt_air=1.d0  !1.0 hour interval
      
      nm_air=1
      do i=1,datamax
        read(77,*,end=1200) fdata(1,i),fdata(2,i),fdata(3,i),fdata(4,i) !sradi(i),eTair(i),eEair(i),ePsea(i)
        nm_air=nm_air+1
      enddo
 1200 close(77)

      allocate( sradi(nm_air),eTair(nm_air),eEair(nm_air),ePsea(nm_air))

      do i=1,nm_air
        sradi(i)=fdata(1,i)
        eTair(i)=fdata(2,i)
        eEair(i)=fdata(3,i)
        ePsea(i)=fdata(4,i)
      enddo
      
#ifdef LONGWAVE_OUT
      
! ----- READ downward longwave radiation (W m-2) ---------------------------------------------

      open(77,file='./input/dwlwradi2009summer.dat')
      dt_dlwrad=1.d0  !1.0 hour interval
      
      nm_dlwrad=1      
      do i=1,datamax
        read(77,*,end=1000) fdata(1,i) !tide(i)
        nm_dlwrad=nm_dlwrad+1
      enddo
 1000 close(77)

      allocate( dlwrad(nm_dlwrad) )

      do i=1,nm_tide
        dlwrad(i)=fdata(1,i)
      enddo

#endif

! --------- READ Wave data ---------------------------------------------

      open(77,file='./input/wave2009summer.dat')
      dt_wave=2.d0  !2.0 hour interval

      nm_wave=1
      do i=1,datamax
        read(77,*,end=1300) fdata(1,i),fdata(2,i) ! Ht(i),Tmt(i)
        nm_wave=nm_wave+1
      end do
 1300 close(77)
 
      allocate( Hsin(nm_wave),Tpin(nm_wave) )
 
      do i=1,nm_wave
        Hsin(i) =fdata(1,i)
        Tpin(i)=fdata(2,i)
      enddo

!-----------READ river data---------------------------------------------




      return

    end subroutine read_timeseies


! **********************************************************************
!  Read chamber condition data file
! **********************************************************************

    subroutine read_chambercondition
    
      USE mod_param
      
      implicit none
      
      integer, parameter :: N_rep = 5    ! Repeat number of same light condition
      integer :: N_data
      integer :: ios
      
      integer i,j
!
      
! ===== READ PPFD data ======================================================
!         time (hour), PPFD (umol m-2 s-1)

      open(77,file='./input/pfd_04.txt')
! ----- Count data -----
      N_PFD = 0
      do
        read(77,*,iostat=ios)
        if(ios==-1) exit
        N_PFD = N_PFD + 1
      end do
      allocate( PFD_time(N_PFD*N_rep), PFD_data(N_PFD*N_rep) )   ! Repeat same light condition for 5 times
      rewind(77) 
! ----- Read data -----
      do j=1, N_rep
        do i=1, N_PFD
          read(77,*) PFD_time(i+N_PFD*(j-1)), PFD_data(i+N_PFD*(j-1))
          if(j >= 2) then
            PFD_time(i+N_PFD*(j-1)) = PFD_time(i+N_PFD*(j-1)) + PFD_time(N_PFD)*dble(j-1)
          end if
        end do
        rewind(77) 
      end do
      N_PFD = N_PFD*N_rep
      close(77)
      
!      do i=1, N_PFD   !!!!!! for DEBUG
!        write(99,*) PFD_time(i)
!      end do


! ===== READ TA and DIC data =================================================
!         time (hour), TA & DIC (umol kg-1)

      open(77,file='./input/site04.txt')
! ----- Count data -----
      N_WQ = 0
      do
        read(77,*,iostat=ios)
        if(ios==-1) exit
        N_WQ = N_WQ + 1
      end do
      allocate( WQ_time(N_WQ), TA_data(N_WQ), DIC_data(N_WQ) )
      rewind(77) 
! ----- Read data -----
      do i=1, N_WQ
        read(77,*) WQ_time(i), TA_data(i), DIC_data(i)
        WQ_time(i) = WQ_time(i) + 4.0d0*24.0d0
      end do
      close(77)
      
!      do i=1, N_WQ   !!!!!! for DEBUG
!        write(99,*) WQ_time(i)
!      end do


      return

    end subroutine read_chambercondition


! **********************************************************************
!  Read 2D mapping data file
! **********************************************************************
!
!    subroutine read_map(im,jm)
!
!      use mod_ecosys
!      
!      implicit none
!      
!      integer im,jm
!      
!      integer i,j
!
!! --------- READ Spring points data ---------------------------------------------
!      
!      open(77,file='input/SpringPoints.txt')
!      do j=jm,1,-1
!        read(77,*) (sprpoint(i,j),i=1,im)
!      enddo
!      close(77)
!
!! --------- READ Coral coverage data ---------------------------------------------
!      
!      open(77,file='./input/coralcover2007-3.txt')
!      do j=jm,1,-1
!        read(77,*) (p_coral(1,i,j),i=1,im)
!      enddo
!      close(77)
!
!! --------- READ seagrass coverage data ---------------------------------------------
!
!      open(77,file='./input/seagrasscover2007-4.txt')
!      do j=jm,1,-1
!        read(77,*) (p_sgrass(1,i,j),i=1,im)
!      enddo
!      close(77)
!
!      return
!
!    end subroutine read_map
!
!
! **********************************************************************
!  Set environmental condition
! **********************************************************************

    subroutine setdata(nSetting)
!     Setting of condition (nSetting)
!
!     nSetting = 1: Stable condition
!                2: Closed chamber condition
!                3: Constant flow condition
!                4: Reef simulation condition
!                5: Incubation chamber condition
    
      USE mod_param
      
      implicit none
      
      integer, intent(in) :: nSetting

!  -- Set Tide data ---------------------------------------------

!      etide=lin_interpol(time,tide,dt_tide,nm_tide)
!      etide=-lin_interpol(time,tide,dt_tide)!!!!!Debug test
      etide=0.!!!!!Debug test

!   -- Set wind data -----------------------------------------------

      Uwind=0.
      Vwind=0.
!      Uwind=lin_interpol(time,windu,dt_wind,nm_wind)
!      Vwind=lin_interpol(time,windv,dt_wind,nm_wind)


      if (nSetting .eq. 5) then
      
        CALL lin_interpol2(time*24.0d0,PFD_time,PFD_data,N_PFD, PFDsurf, i_PFD)
      
      else

!   -- Set solar radiation data -----------------------------------------------

!        ssradi=lin_interpol(time,sradi,dt_air,nm_air)   !data from file

!        ssradi=solar_radi(time,1400.d0,0.d0)!9.d0/24.d0)                !Artificial solar radiation

!        ssradi=light_and_dark(time, 140.0d0, 10./60./24., 0.5/60./24.)   !light and dark method 10 min interval
!        ssradi=light_and_dark(time, 350.0d0, 30./60./24., 0.5/60./24.)   !light and dark method 1 hour interval

!        ssradi=short_radi(time, 0.0d0, 1.0d0, 24.0d0, 25.0d0, 50.0d0, 1) !shortwave radiation by Zillman equation
        ssradi=short_radi(time, 0.0d0, 77.0d0, 24.0d0, 25.0d0, 50.0d0, 2) !shortwave radiation around 3/21 by Zillman equation
      
!     Convert solar radiation (W m-2) to photosynthetic photon flux 
!     density: 2.1 umol m-2 s-1 per W m-2 (Britton and Dodd 1976)
!     Sea surface albedo: 0.07

        PFDsurf = 2.1d0 * ssradi*(1.0d0-0.07d0)
      end if

!----- Temperature controle -----------------------------------------

!      if(time > 5.0d0) then
!
!        do j=1,Jm
!          do i=1,Im
!            do k=1,N
!              C(i,j,k,1,iTemp) = 27.0d0   !27.0d0 32.0d0
!            enddo
!          enddo
!        enddo
!      
!      end if
      
!----- Flux calculation -----------------------------------------

      tau = 1024*0.01*0.0**2. *0.5 !densSW*Cd*Ub**2    (0 cm s-1)
!      tau = 1024*0.01*0.02**2. *0.5  !densSW*Cd*Ub**2  (2 cm s-1)
!      tau = 1024*0.01*0.2**2. *0.5 !densSW*Cd*Ub**2   (20 cm s-1)
      
      pCO2air = 370.0d0 !(uatm)
      
      U10 = SQRT( Uwind*Uwind + Vwind*Vwind )
      
      fvol_pre =0.0d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
#if defined USE_HEAT
!  --- Set heat data ------------------------------------------------
      Tair=lin_interpol(time,eTair,dt_air,nm_air)   !data from file
      Eair=lin_interpol(time,eEair,dt_air,nm_air)   !data from file
      Psea=lin_interpol(time,ePsea,dt_air,nm_air)   !data from file
#endif

      return

    end subroutine setdata

    real(8) function lin_interpol(time,dataset,dt_data,datamax)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  linear interpolation between input data.            *
! *              time:     progress time (day)                         *
! *              dataset(i): data set                                  *
! *              dt_data: data interval (hour)                         *
! *                                                                    *
! **********************************************************************
!
      implicit none

      integer i
!      real lin_interpol
      real(8) dataset(datamax)
      real(8) time,dt_data
      integer datamax

      i=int(time*24.e0/dt_data)+1
      
      if(i.ge.datamax) i=datamax-1

      lin_interpol=dataset(i)+(dataset(i+1)-dataset(i))      &
     &       *(time*24.e0/dt_data-real(i-1))

      return

    end function lin_interpol
    
! **********************************************************************

    subroutine lin_interpol2(time,t_dat,d_dat,N_dat, out_dat, i_dat)
!
      implicit none

!      real lin_interpol
      real(8), intent(in ) :: time
      real(8), intent(in ) :: t_dat(N_dat)
      real(8), intent(in ) :: d_dat(N_dat)
      integer, intent(in ) :: N_dat
      real(8), intent(out) :: out_dat
      integer, intent(inout) :: i_dat

      real(8) :: dt
      integer ::  i,j


      do i=i_dat, N_dat-1
        if( t_dat(i+1)>=time .and. t_dat(i+1)>t_dat(i) ) exit
      end do
      if(i+1.gt.N_dat) then
        out_dat = d_dat(i)
      else
        dt = t_dat(i+1)-t_dat(i)
        out_dat = d_dat(i)+(d_dat(i+1)-d_dat(i)) * (time-t_dat(i))/dt
      end if
      i_dat = i

      return

    end subroutine lin_interpol2

    real(8) function solar_radi(time,PFD,start_time)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  caliculate solar radiation (W/m2).                  *
! *              time:       progress time (day)                       *
! *              start_time: start time (day)                          *
! *                 ex.) 0:00-> start_time=0; 12:00-> start_time=0.5   *
! *                                                                    *
! **********************************************************************

      implicit none
      
      real(8), parameter :: pi = 3.141592654d0

      real(8) time, PFD,start_time

!      solar_radi=PFD/4.57*sin(2.*pi*(time+start_time-0.25))
      solar_radi=PFD/1.82d0*sin(2.*pi*(time+start_time-0.25))
                      !convert photon flux density (umol m-2 s-1) to solar radiation (W m-2)

      if(solar_radi.lt.0.) then
        solar_radi=0.
      endif

      return

    end function solar_radi


    real(8) function short_radi(time,start_time,start_day, lat   &
     &                         ,temp,Hum, mode )
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  Caliculate short wave radiation (W/m2).             *
! *                 by Zillman equation                                *
! *                                                                    *
! *              time:       progress time (day)                       *
! *              start_time: start time (day)                          *
! *                 ex.) 0:00-> start_time=0; 12:00-> start_time=0.5   *
! *              start_day:  start day (day)                           *
! *                 ex.) Jan. 1-> start_day=1; Jun 22-> start_day=174  *
! *                                                                    *
! *              lat: latitude (degree)                                *
! *              Hum: rerative humidity (%)                            *
! *              temp: temperature (degree C)                          *
! *                                                                    *
! *              Hum: rerative humidity (%)                            *
! *                                                                    *
! *              mode: mode =1 -> day time progress                    *
! *                    mode =2 -> loop same day condition              *
! *                                                                    *
! **********************************************************************
!
      implicit none
      
      real(8), parameter :: pi = 3.141592654d0
      real(8), parameter :: Sc = 1353.0d0        !Solar constant (W m-2)

      real(8) :: time, start_time,lat,temp, Hum
      integer :: mode
      real(8) :: start_day
      real(8) :: yDay !day of year
      real(8) :: cosZ
      real(8) :: DA    !DA=2.*pi*yDay/365.
      real(8) :: e_vep !vapor pressure (Pa)
      real(8) :: dec   !declination (radian)
      real(8) :: HA    !hour angle (radian)
      real(8) :: lati
      integer :: i

      e_vep=611.*10.**(7.5*(temp/(temp+273.15-35.86))) * Hum/100.

      HA=(0.5-time)*2.*pi

      if(mode .eq. 1) then
        yDay=start_day+Int(time) ! day time progress
      else
        yDay=start_day            ! light condition releated same day.
      endif
      
      yDay=mod(yDay-1.,365.)+1.

      dec=23.44*cos((172.-yDay)*2.*pi/365.) *pi/180.
!      DA=2.*pi*yDay/365.
!      dec=0.006918-0.399912*cos(DA)   +0.070257*sin(DA)
!     &            -0.006758*cos(2.*DA)+0.000907*sin(2.*DA)
!     &            -0.002697*cos(3.*DA)+0.001480*sin(3.*DA)
!
      lati=lat *pi/180. !degree -> radian

      cosZ=sin(lati)*sin(dec)+cos(lati)*cos(dec)*cos(HA)



      if(cosZ.lt.0.) then
        short_radi=0.
      else
        short_radi=Sc*cosZ**2./((cosZ+2.7)*e_vep*1.e-5 + 1.085*cosZ +0.10)
      endif

!----for debug---------------------------------------------------------c
!        if(i.eq.7200) then
!!        if(iprint.eq.1000) then
!!        if(iprint.eq.1) then
!         write(52,*) time,HA,yDay, dec,lati,cosZ,short_radi,e_vep
!        endif
!------------------------------------------------------------------c        

      return

    end function short_radi


    real(8) function light_and_dark(time,PFD,interval,transperi)
! **********************************************************************
! *                                                                    *
! * FUNCTION    :  caliculate solar radiation (W/m2).                  *
! *              time:       progress time (day)                       *
! *              interval: interval of light & dark (day)              *
! *              transperi: transition period between light & dark (day)*
! *                                                                    *
! **********************************************************************
!
      implicit none
      
      real(8), parameter :: pi = 3.141592654d0

      real(8) :: time, PFD,interval, transperi
      real(8) :: Imax   !light intensity (W/m2)

!            !convert photon flux density (umol m-2 s-1) to solar radiation (W m-2)
      Imax = PFD/4.57 ! for Al-Horani et al. (2003) experiment
!      Imax = 350./4.57 ! for Kuhl et al. (1995) experiment                      

!  Dark->light->Dark->light

!      if(mod(time,interval*2.).lt.interval) then
!        if(mod(time,interval*2.).lt.transperi) then
!          light_and_dark=
!     &           Imax/2.*(Cos(pi*mod(time,interval*2.)/transperi)+1.)
!        else
!          light_and_dark=0.
!        endif
!      else
!        if(mod(time,interval*2.).lt.interval+transperi) then
!          light_and_dark=
!     &          Imax/2.*(-Cos(pi*(mod(time,interval*2.)-interval)
!     &                             /transperi)+1.)
!        else
!          light_and_dark= Imax
!        endif    
!      endif


!  light->Dark->light->Dark

      if(mod(time,interval*2.).lt.interval) then
        if(mod(time,interval*2.).lt.transperi) then
          light_and_dark=                                         &
     &          Imax/2.*(-Cos(pi*(mod(time,interval*2.)-interval) &
     &                             /transperi)+1.)
        else
          light_and_dark= Imax
        endif
      else
        if(mod(time,interval*2.).lt.interval+transperi) then
          light_and_dark=                                            &
     &           Imax/2.*(Cos(pi*mod(time,interval*2.)/transperi)+1.)
        else
          light_and_dark=0.

        endif    
      endif

      return

    end function light_and_dark
      
  end module mod_input

