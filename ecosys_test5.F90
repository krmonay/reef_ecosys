
!!!=== ver 2017/03/06   Copyright (c) 2012-2017 Takashi NAKAMURA  =====

#include "cppdefs.h"


      PROGRAM ecosys_test
! **********************************************************************
! *                                                                    *
! *   Test program of mod_reef_ecosys3                                 *
! *                                                                    *
! **********************************************************************
!
      USE mod_param
      USE mod_reef_ecosys
#if defined USE_HEAT
      USE mod_heat
#endif
      USE mod_input
      USE mod_geochem

      implicit none
      
      integer, parameter :: Im = 2
      integer, parameter :: Jm = 1
      integer, parameter :: N  = 1

      real(8), parameter :: dt = 1.0d0
      integer, parameter :: isplitc  = 20
      integer, parameter :: isplitsed  = 1
      integer, parameter :: Ngrids = 1


      integer :: i,j,k,id, Nid
      integer :: istep, iprint
      integer :: nSetting, nheat
      real(8) :: time    !(day)
      real(8) :: R13C

      real(8) :: PFDsurf    
      real(8) :: tau        
      real(8) :: pCO2air    
      real(8) :: U10        
      
      real(8) :: sspH      
      real(8) :: ssfCO2    
      real(8) :: ssWarg    
      real(8) :: ssCO2flux 
      real(8) :: ssO2flux  
      real(8) :: PFDbott   
      real(8) :: coral_Pg  
      real(8) :: coral_R   
      real(8) :: coral_Gn  
      real(8) :: sgrass_Pg 
      real(8) :: sgrass_R  
      
      real(8) :: d13C_DIC
      
      real(8) :: fvol_cre      ! volume flux through the reef crest (m3 m-2 s-1)
      real(8) :: fvol_cha      ! volume flux through the channel(m3 m-2 s-1)
      real(8) :: fvol_pre      ! Precipitation volume flux (m s-1)
!      real(8) :: dw_lwradi     ! Downward longwave radiation (W m-2)
      real(8), save :: ereef = 0.0d0   ! sea surface elevation on the reef flat (m)
      real(8), parameter :: z_crest = -0.15d0 ! reef crest position (m)
      real(8), parameter :: kvol_cre   = 5.0d-2 ! reef crest conductivity
      real(8), parameter :: kvol_cha = 1.0d-1 ! Channel conductivety
      real(8) :: wave_setup   ! wave setup (m)
      
!  For Output      
      real(8), save :: dsec = 0.d0 !sec


      open(52,file='./output/eco4-crl_his.txt')
      open(53,file='./output/eco4-crl_ave.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(54,file='./output/eco4-zoo_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(55,file='./output/eco4-box_his.txt')!!!!!!!!!!!!!!!!!!!for debug

      open(56,file='./output/eco4-sedDIC_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(57,file='./output/eco4-sedTA_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(58,file='./output/eco4-sedDO_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(59,file='./output/eco4-sedpH_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(60,file='./output/eco4-sedWarg_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(61,file='./output/eco4-sedNH4_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(62,file='./output/eco4-sedNO2_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(63,file='./output/eco4-sedNO3_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(64,file='./output/eco4-sedPO4_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(65,file='./output/eco4-sedDOC_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(66,file='./output/eco4-sedPOC_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(67,file='./output/eco4-sedDON_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(68,file='./output/eco4-sedPON_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(69,file='./output/eco4-sedDOP_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(70,file='./output/eco4-sedPOP_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(71,file='./output/eco4-sedPg_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(72,file='./output/eco4-sedRdoc_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(73,file='./output/eco4-sedRpoc_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(74,file='./output/eco4-sedGn_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(75,file='./output/eco4-sedNit1_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(76,file='./output/eco4-sedNit2_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(78,file='./output/eco4-sedDNd_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(79,file='./output/eco4-sedDNp_his.txt')!!!!!!!!!!!!!!!!!!!for debug


      istep=0
      iprint=0
      
      
      
! Setting of condition (nsetting)
!
!     nSetting = 1: Stable condition
!                2: Closed chamber condition
!                3: Constant flow condition
!                4: Reef simulation condition
!                5: Incubation chamber condition
      
      nSetting = 5  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
!----- Set initial conditions -------------------------

      CALL initialize_params(1, Im, 1, Jm, N, Nid)

!----- Import data -----------------------------------

      CALL read_timeseies

!----- Set initial conditions -------------------------

#if defined USE_HEAT
      CALL initialize_heat(1, Ngrids, 1, Im, 1, Jm)
#endif
      CALL initialize_reef_ecosys(1, Ngrids, 1, Im, 1, Jm)

!      call Coral_iniSizeDis
!      call Coral_Size2Cover

      time=0.

!----- Main loop -------------------------------------------

!      do istep=1, int(60.*60./dt)+1       ! 1 hour
      do istep=1, int(24.*60.*60./dt) * 5 +1      ! 5 days
!      do istep=1, int(24.*60.*60./dt) * 7 +1      ! 7 days
!      do istep=1, int(24.*60.*60./dt) * 30 +1      ! 30 days
!      do istep=1, int(24.*60.*60./dt) *365*3 +1      ! 3 year

!------ Set environmental parameters ----------------------------

        CALL setdata(time)

!      if(time > 5.0d0) then

        do j=1,Jm
          do i=1,Im
            do k=1,N
!              C(i,j,k,1,iTemp) = 27.0d0   !27.0d0 32.0d0
            enddo
          enddo
        enddo
      
!      end if
      
!----- Flux calculation -----------------------------------------

        
!       Fraction of shortwave radiation that is photosynthetically active
!       (nondimensional), {0.43d0}.          
!       convert radiation (W m-2) to photon flux density (umol m-2 s-1)  1 W m-2 = 4.24 umol m-2 s-1
!       Thus convert solar radiation (W m-2) to photosynthetic photon flux density: 0.43*4.24 = 1.82

        PFDsurf = 1.82d0 * ssradi
        
        
        tau = 1024*0.01*0.0**2. *0.5 !densSW*Cd*Ub**2    (0 cm s-1)
!        tau = 1024*0.01*0.02**2. *0.5  !densSW*Cd*Ub**2  (2 cm s-1)
!        tau = 1024*0.01*0.2**2. *0.5 !densSW*Cd*Ub**2   (20 cm s-1)
        
        pCO2air = 370.0d0 !(uatm)
        
        U10 = SQRT( Uwind*Uwind + Vwind*Vwind )
        
        fvol_pre =0.0d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined USE_HEAT

!----- Heat & mass balance model ----------------------------------------

        CALL heat_mass_balance    &
!          input parameters
     &            (1,1,1          &   ! ng: nested grid number; i,j: position
     &            ,N              &   ! Number of vertical grid (following ROMS vertical grid)
     &            ,dt             &   ! Time step (sec)
     &            ,dz(1,1,:)      &   ! dz(N): vertical grid size (m)
     &            ,ssradi         &   ! Surface shortwave radiation (W m-2)
     &            ,Tair           &   ! air temperature (oC)
     &            ,Psea           &   ! atm pressure (hPa)
     &            ,Eair           &   ! vapor pressur (hPa)
     &            ,U10            &   ! wind speed (m s-1)
     &            ,fvol_pre       &   ! Precipitation volume flux (m s-1)
#ifdef LONGWAVE_OUT
     &            ,dw_lwradi      &   ! Downward longwave radiation (W m-2)
#endif

     &            ,C(1,1,:,1,iTemp)     &   ! Tmp(N): Temperature (oC)
     &            ,C(1,1,:,1,iSalt)     &   ! Sal(N): Salinity (PSU)

!          output parameters
     &            ,dC_dt(1,1,:,iTemp)   &   ! dDIC_dt(N): dDIC/dt (umol kg-1 s-1)
     &             )
#endif

!----- Ecosystem model ----------------------------------------

        CALL reef_ecosys          &
!          input parameters
     &            (1, 1, 1        &   ! ng: nested grid number; i,j: position
     &            ,N              &   ! Number of vertical grid (following ROMS vertical grid)
     &            ,isplitc        &   ! Internal loop counts of coral polyp model
     &            ,isplitsed      &   ! Internal loop counts of sediment ecosystem model
     &            ,dt             &   ! Time step (sec)
     &            ,dz(1,1,:)      &   ! dz(N): vertical grid size (m)
     &            ,PFDsurf        &   ! Sea surface photon flux density (umol m-2 s-1)
     &            ,tau            &   ! bottom shear stress (N m-2)
     &            ,pCO2air        &   ! Air CO2 pertial pressure (uatm)
     &            ,U10            &   ! wind speed (m s-1)
#ifdef CORAL_POLYP
     &            ,p_coral(:,1,1) &   ! Coral coverage (0-1)
#endif
#ifdef SEAGRASS
     &            ,p_sgrass(1,1)  &   ! seagrass coverage (0-1)
#endif
#ifdef MACROALGAE
     &            ,p_algae(1,1)   &   ! algal coverage (0-1)
#endif
#ifdef SEDIMENT_ECOSYS
     &            ,p_sand(1,1)    &   ! sediment coverage (0-1)
#endif

     &            ,C(1,1,:,1,iTemp)     &   ! Tmp(N): Temperature (oC)
     &            ,C(1,1,:,1,iSalt)     &   ! Sal(N): Salinity (PSU)
     &            ,C(1,1,:,1,iTIC_)     &   ! DIC(N): Total dissolved inorganic carbon (DIC: umol kg-1)
     &            ,C(1,1,:,1,iTAlk)     &   ! TA (N): Total alkalinity (TA: umol kg-1)
     &            ,C(1,1,:,1,iOxyg)     &   ! DOx(N): Dissolved oxygen (umol L-1)
#if defined ORGANIC_MATTER
     &            ,C(1,1,:,1,iDOC_)     &   ! DOC(N): Dissolved organic carbon (DOC: umol L-1)
     &            ,C(1,1,:,1,iPOC_)     &   ! POC(N): Particulate organic carbon (POC: umol L-1)
     &            ,C(1,1,:,1,iPhyt1)     &   ! PHY(N): phytoplankton (umol C L-1)
     &            ,C(1,1,:,1,iPhyt2)     &   ! PHY(N): phytoplankton (umol C L-1)
     &            ,C(1,1,:,1,iZoop)     &   ! ZOO(N): zooplankton (umol C L-1)
#endif
#if defined CARBON_ISOTOPE
     &            ,C(1,1,:,1,iT13C)     &   ! DI13C(N): 13C of DIC (umol kg-1)
#endif
#if defined NUTRIENTS         
     &            ,C(1,1,:,1,iNO3_)     &   ! NO3(N): NO3 (umol L-1)
     &            ,C(1,1,:,1,iNO2_)     &   ! NO2(N): NO2 (umol L-1)
     &            ,C(1,1,:,1,iNH4_)     &   ! NH4(N): NH4 (umol L-1)
     &            ,C(1,1,:,1,iPO4_)     &   ! PO4(N): PO4 (umol L-1)
# if defined ORGANIC_MATTER
     &            ,C(1,1,:,1,iDON_)     &   ! DON(N): Dissolved organic nitrogen (DON: umol L-1)
     &            ,C(1,1,:,1,iPON_)     &   ! PON(N): Particulate organic nitrogen (PON: umol L-1)
     &            ,C(1,1,:,1,iDOP_)     &   ! DOP(N): Dissolved organic phosporius (DOP: umol L-1)
     &            ,C(1,1,:,1,iPOP_)     &   ! POP(N): Particulate organic phosporius (POP: umol L-1)
# endif
#endif
#if defined COT_STARFISH         
     &            ,C(1,1,:,1,iCOTe)     &   ! COTe(N): COT starfish egg (umol L-1)
     &            ,C(1,1,:,1,iCOTl)     &   ! COTl(N): COT starfish larvae (umol L-1)
#endif
!          output parameters
     &            ,dC_dt(1,1,:,iTIC_)   &   ! dDIC_dt(N): dDIC/dt (umol kg-1 s-1)
     &            ,dC_dt(1,1,:,iTAlk)   &   ! dTA_dt (N): dTA/dt (umol kg-1 s-1)
     &            ,dC_dt(1,1,:,iOxyg)   &   ! dDOx_dt(N): dDO/dt (umol L-1 s-1)
#if defined ORGANIC_MATTER
     &            ,dC_dt(1,1,:,iDOC_)   &   ! dDOC_dt(N): dDOC/dt (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iPOC_)   &   ! dPOC_dt(N): dPOC/dt (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iPhyt1)   &   ! dPHY_dt(N): dPHY/dt (umol C L-1 s-1)
     &            ,dC_dt(1,1,:,iPhyt2)   &   ! dPHY_dt(N): dPHY/dt (umol C L-1 s-1)
     &            ,dC_dt(1,1,:,iZoop)   &   ! dZOO_dt(N): dZOO/dt (umol C L-1 s-1)
#endif
#if defined CARBON_ISOTOPE
     &            ,dC_dt(1,1,:,iT13C)   &   ! dDI13C_dt(N): dDI13C/dt (umol kg-1 s-1)
#endif
#if defined NUTRIENTS 
     &            ,dC_dt(1,1,:,iNO3_)   &   ! dNO3_dt(N): dNO3/dt (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iNO2_)   &   ! dNO2_dt(N): dNO2/dt (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iNH4_)   &   ! dNH4_dt(N): dNH4/dt (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iPO4_)   &   ! dPO4_dt(N): dPO4/dt (umol L-1 s-1)
# if defined ORGANIC_MATTER
     &            ,dC_dt(1,1,:,iDON_)   &   ! dDON_dt(N): dDON/dt (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iPON_)   &   ! dPON_dt(N): dPON/dt (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iDOP_)   &   ! dDOP_dt(N): dDOP/dt (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iPOP_)   &   ! dPOP_dt(N): dPOP/dt (umol L-1 s-1)
# endif
#endif
#if defined COT_STARFISH         
     &            ,dC_dt(1,1,:,iCOTe)   &   ! dCOTe/dt(N): (umol L-1 s-1)
     &            ,dC_dt(1,1,:,iCOTl)   &   ! dCOTl/dt(N): (umol L-1 s-1)
#endif
     &            ,sspH           &   ! sea surface pH
     &            ,ssfCO2         &   ! sea surface fCO2 (uatm)
     &            ,ssWarg         &   ! sea surface aragonite saturation state
     &            ,ssCO2flux      &   ! sea surface CO2 flux (mmol m-2 s-1)
     &            ,ssO2flux       &   ! sea surface O2 flux (mmol m-2 s-1)
     &            ,PFDbott        &   ! Bottom photon flux density (umol m-2 s-1)
     &             )                    
!

        do k=1,N
        

!---------- for Stable condition ---------------------------------
              
          if (nSetting .eq. 1) then
              
            ! nothing to calculate

!---------- for Closed Chamber condition -------------------------

          else if (nSetting .eq. 2) then
            
            C(1,1,k,1,iTemp)=C(1,1,k,1,iTemp)+0.
            C(1,1,k,1,iSalt)=C(1,1,k,1,iSalt)+0.
            C(1,1,k,1,iSedi)=C(1,1,k,1,iSedi)+0.
            
            do id=4,Nid
              C(1,1,k,1,id)=C(1,1,k,1,id) + dC_dt(1,1,k,id)*dt
            end do

!---------- Constant Flow condition ------------------------------------

          else if (nSetting .eq. 3) then









!---------- Reef condition ------------------------------------

          else if (nSetting .eq. 4) then
!            wave_setup=0.20d0 !(m)
            wave_setup=0.12*Hs !(m)
            wave_setup=0.12*Hs +0.1d0 !(m)

            if (max(etide+wave_setup,ereef) .gt. z_crest) then
              fvol_cre = kvol_cre*(etide+wave_setup - ereef) * ABS(etide+wave_setup - ereef)  ! volume flux (m3 m-2 s-1)
!              if (etide+wave_setup .gt. ereef) then
!                fvol_cre =  kvol_cre*(ereef-z_crest)* (etide+wave_setup - ereef)**0.5d0  ! volume flux (m3 m-2 s-1)
!              else
!                fvol_cre = -kvol_cre*(etide+wave_setup-z_crest)* (ereef - (etide+wave_setup))**0.5d0  ! volume flux (m3 m-2 s-1)
!              end if
            else
              fvol_cre = 0.0d0
            end if
            
            fvol_cha = kvol_cha*(etide - ereef) * ABS(etide - ereef)  ! volume flux (m3 m-2 s-1)
!            if (etide .gt. ereef) then
!              fvol_cha =  kvol_cha *(ereef-(z_crest-1.0d0))* (etide - ereef)**0.5d0  ! volume flux (m3 m-2 s-1)
!            else
!              fvol_cha = -kvol_cha *(etide-(z_crest-1.0d0))* (ereef - etide)**0.5d0  ! volume flux (m3 m-2 s-1)
!            end if
            
            ereef = ereef + (fvol_cre + fvol_cha) * dt
            dz(1,1,k) = dz(1,1,k) + (fvol_cre + fvol_cha) * dt/N 
            
            do id=1,Nid
              C(1,1,k,1,id)=C(1,1,k,1,id) + dC_dt(1,1,k,id)*dt               &
     &                      +(0.5d0*(ABS(fvol_cre)+fvol_cre)* (C(2,1,k,1,id)-C(1,1,k,1,id))  &  !  (t unit) m s-1
!     &                       -0.5d0*(ABS(fvol_cre)-fvol_cre)* C(1,1,k,1,id)  &
     &                       +0.5d0*(ABS(fvol_cha)+fvol_cha)* (C(2,1,k,1,id)-C(1,1,k,1,id))   &
!     &                       -0.5d0*(ABS(fvol_cha)-fvol_cha)* C(1,1,k,1,id)  &
     &                       )*dt/dz(1,1,k)
            end do

!---------- Incubation chamber condition ------------------------------------

          else if (nSetting .eq. 5) then

            C(1,1,k,1,iTemp)=C(1,1,k,1,iTemp)+0.
            C(1,1,k,1,iSalt)=C(1,1,k,1,iSalt)+0.
            C(1,1,k,1,iSedi)=C(1,1,k,1,iSedi)+0.
            
            do id=4,Nid
              C(1,1,k,1,id)=C(1,1,k,1,id) + dC_dt(1,1,k,id)*dt
            end do

          end if

!              depsed(i,j)=0.
!              radi(i,j,k)=-ssradi

              
         enddo


!        p_coral(1,i,j)=p_coral(1,i,j)
!     &          +(g_coral(1,i,j)-m_coral(1,i,j))*p_coral(1,i,j)
!
!        if (p_coral(1,i,j) .lt. 1.e-5) then
!          p_coral(1,i,j)=0.
!        end if


!        time=dt*istep/86400.
        time=time+dt/86400.

#if defined TESTMODE

        if(time .ge. dsec/86400.) then
        
          write(*,*) 'Time (day): ', time  ! Output for standard out
          
# if defined CARBON_ISOTOPE
          d13C_DIC=d13C_fromR13C(C(1,1,1,1,iT13C)/C(1,1,1,1,iTIC_))
# endif
        
          write(55,*) time, PFDsurf                              &
     &       ,coral_Pg, coral_Pg-coral_R, coral_Gn, sgrass_Pg, sgrass_R  &
     &       ,dz(1,1,1),C(1,1,1,1,iTemp),C(1,1,1,1,iSalt)        &
     &       ,sspH, ssfCO2, ssWarg, U10, ssCO2flux, ssO2flux     &
     &       ,etide, ereef, dz(1,1,1)                            &
     &       ,C(1,1,1,1,iTAlk),C(1,1,1,1,iTIC_),C(1,1,1,1,iOxyg) &
#if defined ORGANIC_MATTER
     &       ,C(1,1,1,1,iDOC_),C(1,1,1,1,iPOC_),C(1,1,1,1,iPhyt1) &
     &       ,C(1,1,1,1,iPhyt2),C(1,1,1,1,iZoop)                  &
#else
     &       ,0.0d0, 0.0d0, 0.0d0 &
     &       ,0.0d0        &
#endif
#if defined NUTRIENTS            
     &       ,C(1,1,1,1,iNO3_),C(1,1,1,1,iNO2_),C(1,1,1,1,iNH4_) &
     &       ,C(1,1,1,1,iPO4_)                                   &
# if defined ORGANIC_MATTER
     &       ,C(1,1,1,1,iDON_),C(1,1,1,1,iPON_),C(1,1,1,1,iDOP_) &
     &       ,C(1,1,1,1,iPOP_)                                   &
# else
     &       ,0.0d0, 0.0d0, 0.0d0 &
     &       ,0.0d0        &
# endif
#else
     &       ,0.0d0, 0.0d0, 0.0d0 &
     &       ,0.0d0        &
#endif
#if defined COT_STARFISH
     &       ,C(1,1,1,1,iCOTe),C(1,1,1,1,iCOTl)
#else
     &       ,0.0d0, 0.0d0
# endif
     
     
!          dsec=dsec+10.*60.  !!!!!!!!!!!!!!! print 10 min interval 
        dsec=dsec+60.*60.  !!!!!!!!!!!!!!! print 1 hour interval 

        endif
#endif

      enddo
      
!----- End loop --------------------------------------

      end
!----------------------------------------------------------------------!

!     End of main program

!-----------------------------------------------------------------------

