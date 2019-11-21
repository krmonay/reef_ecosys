
!!!=== Copyright (c) 2012-2019 Takashi NAKAMURA  =====

!--------------------------------------------------------------------------------
!
!  Parameter module
!
!--------------------------------------------------------------------------------

#include "cppdefs.h"

MODULE mod_param
  
  implicit none
!
!-----------------------------------------------------------------------
!  Grid nesting parameters.
!-----------------------------------------------------------------------

  real(8) :: time    !(day)
  
  real(8), allocatable :: tide_data(:)
  real(8) :: tide
  real(8) :: dt_tide
  integer :: nm_tide
  
  real(8), allocatable :: Hs_data(:)
  real(8), allocatable :: Tp_data(:)
  real(8) :: Hs
  real(8) :: Tp
  real(8) :: dt_wave
  integer :: nm_wave
  
  real(8), allocatable :: Pair_data (:)
  real(8), allocatable :: Uwind_data(:)
  real(8), allocatable :: Vwind_data(:)
  real(8), allocatable :: Tair_data (:)
  real(8), allocatable :: Qair_data (:)
  real(8), allocatable :: rain_data (:)
  real(8), allocatable :: swrad_data(:)
  real(8), allocatable :: cloud_data(:)
  real(8) :: Pair
  real(8) :: Uwind
  real(8) :: Vwind
  real(8) :: Tair
  real(8) :: Qair
  real(8) :: rain
  real(8) :: swrad
  real(8) :: cloud
  real(8) :: dt_air
  integer :: nm_air
  
  real(8), allocatable :: dlwrad_data(:)
  real(8) :: dlwrad
  real(8) :: dt_dlwrad
  integer :: nm_dlwrad
  
  real(8), allocatable :: offtemp_data(:)
  real(8) :: offtemp
  real(8) :: dt_offtemp
  integer :: nm_offtemp
  
  real(8), allocatable :: PFD_time(:)
  real(8), allocatable :: PFD_data(:)
  integer :: N_PFD
  integer :: i_PFD = 1
  real(8), allocatable :: WQ_time(:)
  real(8), allocatable :: TA_data(:)
  real(8), allocatable :: DIC_data(:)
  real(8), allocatable :: DO_data(:)
  integer :: N_WQ
  integer :: i_WQ = 1


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

  real(8) :: fvol_cre      ! volume flux through the reef crest (m3 m-2 s-1)
  real(8) :: fvol_cha      ! volume flux through the channel(m3 m-2 s-1)
  real(8) :: fvol_pre      ! Precipitation volume flux (m s-1)
!  real(8) :: dw_lwradi     ! Downward longwave radiation (W m-2)
  real(8) :: ereef         ! sea surface elevation on the reef flat (m)
  real(8), parameter :: z_crest = -0.15d0 ! reef crest position (m)
  real(8), parameter :: kvol_cre   = 5.0d-2 ! reef crest conductivity
  real(8), parameter :: kvol_cha = 1.0d-1 ! Channel conductivety
  real(8) :: wave_setup   ! wave setup (m)


  integer :: iTemp                  ! Temperature
  integer :: iSalt                  ! Salinity
  integer :: iSedi                  ! Sediment
  integer :: iTIC_                  ! Total inorganic carbon
  integer :: iTAlk                  ! Total alkalinity
  integer :: iOxyg                  ! Dissolved oxygen concentration
#if defined ORGANIC_MATTER
  integer :: iDOC_                  ! Dissolved organic C-concentration
  integer :: iPOC_                  ! Particulate organic C-concentration
  integer :: iPhyt1                  ! Phytoplankton concentration
  integer :: iPhyt2                  ! Phytoplankton concentration
  integer :: iZoop                  ! Zooplankton concentration
#endif
#if defined CARBON_ISOTOPE
  integer :: iT13C                  ! Corbon 13 of total inorganic carbon
# if defined ORGANIC_MATTER
  integer :: iDO13                  ! Dissolved organic 13C-concentration
  integer :: iPO13                  ! Particulate organic 13C-concentration
  integer :: iPh13                  ! Phytoplankton 13C-concentration
  integer :: iZo13                  ! Zooplankton 13C-concentration
# endif
#endif
#if defined NUTRIENTS
  integer :: iNO3_                  ! Nitrate concentration
  integer :: iNO2_                  ! Nitrite concentration
  integer :: iNH4_                  ! Ammonium concentration
  integer :: iPO4_                  ! Ammonium concentration
# if defined ORGANIC_MATTER
  integer :: iDON_                  ! Dissolved organic N-concentration
  integer :: iPON_                  ! Particulate organic N-concentration
  integer :: iDOP_                  ! Dissolved organic P-concentration
  integer :: iPOP_                  ! Particulate organic P-concentration
# endif
#endif
#if defined COT_STARFISH
  integer :: iCOTe                  ! Eggs of crown-of-thorns starfish
  integer :: iCOTl                  ! Larvae of crown-of-thorns starfish
#endif

  real(8), allocatable :: dz(:,:,:)    
  real(8), allocatable :: C(:,:,:,:,:)
  real(8), allocatable :: dC_dt(:,:,:,:)
  real(8), allocatable :: Co(:,:,:)
  
  real(8), allocatable :: p_coral(:,:,:)
  real(8), allocatable :: p_sgrass(:,:)
  real(8), allocatable :: p_algae(:,:)
  real(8), allocatable :: p_sand(:,:)

contains

!!! **********************************************************************
!!!  Set initial conditions
!!! **********************************************************************

  subroutine initialize_params(LBi, UBi, LBj, UBj, N, Nid)

    use mod_geochem
    
    implicit none
! input parameters
    integer, intent(in) :: LBi, UBi, LBj, UBj, N
    integer, intent(out) :: Nid
    real(8)  R13C
    integer i,j,k
!  Local variable declarations
!
    integer :: ic

!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
    ic=0

    i=1
    iTemp=ic+i
    i=i+1
    iSalt=ic+i

    i=i+1
    iSedi=ic+i

    i=i+1
    iTIC_=ic+i
    i=i+1
    iTAlk=ic+i
    i=i+1
    iOxyg=ic+i  !  4
#if defined ORGANIC_MATTER
    i=i+1
    iDOC_=ic+i
    i=i+1
    iPOC_=ic+i
    i=i+1
    iPhyt1=ic+i
    i=i+1
    iPhyt2=ic+i
    i=i+1
    iZoop=ic+i
#endif
#if defined CARBON_ISOTOPE
    i=i+1
    iT13C=ic+i  ! +1
#endif
#if defined NUTRIENTS
    i=i+1
    iNO3_=ic+i
    i=i+1
    iNO2_=ic+i
    i=i+1
    iNH4_=ic+i
    i=i+1
    iPO4_=ic+i
# if defined ORGANIC_MATTER
    i=i+1
    iDON_=ic+i
    i=i+1
    iPON_=ic+i
    i=i+1
    iDOP_=ic+i
    i=i+1
    iPOP_=ic+i
# endif
#endif
#if defined COT_STARFISH
    i=i+1
    iCOTe=ic+i
    i=i+1
    iCOTl=ic+i
#endif

!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------

    Nid=i
!---------------------------------------------------------------------


   allocate( dz(LBi:UBi,LBj:UBj,N) )

   allocate( C(LBi:UBi, LBj:UBj, N, 1, Nid) , &
             dC_dt(LBi:UBi, LBj:UBj, N, Nid)   )
   allocate( Co(LBi:UBi, LBj:UBj, Nid) )
   
   allocate( p_coral (2,LBi:UBi, LBj:UBj) , &
             p_sgrass(LBi:UBi, LBj:UBj) , &
             p_sand  (LBi:UBi, LBj:UBj) , &
             p_algae (LBi:UBi, LBj:UBj)   )   


!-----------------------------------------------------------------------
!  Set initial conditions
!-----------------------------------------------------------------------

    do j=LBj,UBj
      do i=LBi,UBi
        do k=1,N
        
!          dz(i,j,k)=0.3d0 !(m)
          dz(i,j,k)=0.01d0 !(m)
        
!          C(i,j,k,1,iTemp) = 28.0d0   !27.0d0 32.0d0
          C(i,j,k,1,iTemp) = 32.0d0   !32.0d0 32.0d0
          
          C(i,j,k,1,iSalt) = 34.0d0
!          C(i,j,k,1,iSalt) = 35.9d0

          C(i,j,k,1,iSedi) = 0.0d0    !Sediment concentration (g m-3) 0.e0, 1.e0

          C(i,j,k,1,iTIC_) = 1915.0d0     !DIC  (umol kg-1) 1915.0d0 1930.0d0, 1700.0d0, 2100.0d0, 2300.0d0, 2500.0d0, 2700.0d0, 3000.0d0 
          C(i,j,k,1,iTAlk) = 2232.0d0     !TA  (umol kg-1)  2232.0d0 2275.0d0, 2500.0d0, 2150.0d0, 2000.0d0, 1800.0d0, 1700.0d0, 1600.0d0 

!          C(i,j,k,1,iOxyg) = 200.0d0      !DO  (umol L-1)
          C(i,j,k,1,iOxyg) = O2satu(C(i,j,k,1,iTemp)+273.15d0, C(i,j,k,1,iSalt))

! Offshore values
          Co(i,j,iTemp) = 33.0d0   !27.0d0 32.0d0
          Co(i,j,iSalt) = 34.0d0

          Co(i,j,iSedi) = 0.0d0    !Sediment concentration (g m-3) 0.e0, 1.e0

          Co(i,j,iTIC_) = 1915.0d0     !DIC  (umol kg-1) 1915.0d0 1930.0d0, 1700.0d0, 2100.0d0, 2300.0d0, 2500.0d0, 2700.0d0, 3000.0d0 
          Co(i,j,iTAlk) = 2232.0d0     !TA  (umol kg-1)  2232.0d0 2275.0d0, 2500.0d0, 2150.0d0, 2000.0d0, 1800.0d0, 1700.0d0, 1600.0d0 

          Co(i,j,iOxyg) = O2satu(C(i,j,k,1,iTemp)+273.15d0, C(i,j,k,1,iSalt))
#if defined ORGANIC_MATTER
          C(i,j,k,1,iDOC_) = 65.0d0       !DOC  (umol L-1) 
          C(i,j,k,1,iPOC_) = 4.2d0        !POC  (umol L-1) 
          C(i,j,k,1,iPhyt1) =  2.23d0       !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          C(i,j,k,1,iPhyt2) =  2.23d0       !Phytoplankton2 0.3(umolC L-1) 							diatom
          C(i,j,k,1,iZoop) =  1.3d0       !Zooplankton (umol L-1)1.3
! Offshore values
          Co(i,j,iDOC_) = 65.0d0       !DOC  (umol L-1) 
          Co(i,j,iPOC_) = 4.2d0        !POC  (umol L-1) 
          Co(i,j,iPhyt1) =  2.23d0       !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
          Co(i,j,iPhyt2) =  2.23d0       !Phytoplankton2 0.3(umolC L-1) 							diatom
          Co(i,j,iZoop) =  1.3d0       !Zooplankton (umol L-1)1.3
#endif
#if defined CARBON_ISOTOPE
          R13C=R13C_fromd13C(0.7d0)
!          Ci(nDI13C,i,j,k)=R13C/(1.+R13C)*Ci(nDIC,i,j,k) !DI13C (umol kg-1)
          C(i,j,k,1,iT13C) =R13C*C(i,j,k,1,iTIC_) !DI13C (umol kg-1) 
! Offshore values
          Co(i,j,iT13C) =R13C*C(i,j,k,1,iTIC_) !DI13C (umol kg-1) 
#endif
#if defined NUTRIENTS            
!!!!!!!!!!!!
          C(i,j,k,1,iNO3_) =  0.2d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
!          C(i,j,k,1,iNO3_) =  1.8d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
!					 C(i,j,k,1,iNO3_) =  3.5d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0

          C(i,j,k,1,iPO4_) =  0.04d0      !PO4  (umol L-1) 0.03d0 0.05d0, 2.0d0 control0.04d0 P0.2d0
!          C(i,j,k,1,iPO4_) =  0.2d0      !PO4  (umol L-1) 0.03d0 0.05d0, 2.0d0 control0.04d0 P0.2d0
!!!!!!!!!!!!
          C(i,j,k,1,iNO2_) =  0.02d0      !NO2  (umol L-1) 
          C(i,j,k,1,iNH4_) =  0.26d0       !NH4  (umol L-1) 
! Offshore values
          Co(i,j,iNO3_) =  0.2d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
          Co(i,j,iPO4_) =  0.04d0      !PO4  (umol L-1) 0.03d0 0.05d0, 2.0d0 control0.04d0 P0.2d0
          Co(i,j,iNO2_) =  0.02d0      !NO2  (umol L-1) 
          Co(i,j,iNH4_) =  0.26d0       !NH4  (umol L-1) 

# if defined ORGANIC_MATTER
          C(i,j,k,1,iDON_) =  10.0d0 !3.0d0       !DON  (umol L-1) 
          C(i,j,k,1,iPON_) =  0.6d0 !0.05d0      !PON  (umol L-1) 
          C(i,j,k,1,iDOP_) =  0.6d0 !0.1d0       !DOP  (umol L-1) 
          C(i,j,k,1,iPOP_) =  0.04d0 !0.015d0     !POP  (umol L-1) 
! Offshore values
          Co(i,j,iDON_) =  10.0d0 !3.0d0       !DON  (umol L-1) 
          Co(i,j,iPON_) =  0.6d0 !0.05d0      !PON  (umol L-1) 
          Co(i,j,iDOP_) =  0.6d0 !0.1d0       !DOP  (umol L-1) 
          Co(i,j,iPOP_) =  0.04d0 !0.015d0     !POP  (umol L-1) 
# endif
#endif
#if defined COT_STARFISH
          C(i,j,k,1,iCOTe) =  0.0d0       !COTS eggs (umol L-1) 
          C(i,j,k,1,iCOTl) =  0.0d0       !COTS larvae (umol L-1) 
#endif
        enddo
          
#if defined CHAMBER_SITE4
        p_coral(1,i,j) = 0.208d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.0d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 1.0d0-p_coral(1,i,j)
#elif defined CHAMBER_SITE5
        p_coral(1,i,j) = 0.265d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.0d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 1.0d0-p_coral(1,i,j)
#elif defined CHAMBER_SITE6
        p_coral(1,i,j) = 0.37d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.0d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 1.0d0-p_coral(1,i,j)
#elif defined CHAMBER_SITE7
        p_coral(1,i,j) = 0.231d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.0d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 1.0d0-p_coral(1,i,j)
#elif defined CHAMBER_SITE9
        p_coral(1,i,j) = 0.0d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.412d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 0.0d0
#elif defined CHAMBER_SITE10
        p_coral(1,i,j) = 0.0d0  ! Site4: 0.208d0, Site5: 0.265d0, Site6: 0.37d0, Site7: 0.231d0
        p_coral(2,i,j) = 0.486d0    ! Site9: 0.412d0, Site10: 0.486d0
        p_sand(i,j)  = 0.0d0

#else
        p_coral(1,i,j) = 0.2d0
        p_coral(2,i,j) = 0.0d0
        p_sand(i,j)  = 0.0d0
#endif
        p_algae(i,j) = 0.0d0
        p_sgrass(i,j)= 0.0d0
        
        
      enddo
    enddo
    
    ereef = 0.0d0

    return
  end subroutine initialize_params

END MODULE mod_param

