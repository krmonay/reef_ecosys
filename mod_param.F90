
!!!=== ver 2017/03/06   Copyright (c) 2012-2017 Takashi NAKAMURA  =====

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
!
!  Number of nested and/or connected grids to solve.
!

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

    real(8), allocatable, save :: dz(:,:,:)    
    real(8), allocatable, save :: C(:,:,:,:,:)
    real(8), allocatable, save :: dC_dt(:,:,:,:)
    
    real(8), allocatable, save :: p_coral(:,:,:)
    real(8), allocatable, save :: p_sgrass(:,:)
    real(8), allocatable, save :: p_algae(:,:)
    real(8), allocatable, save :: p_sand(:,:)

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
     &          dC_dt(LBi:UBi, LBj:UBj, N, Nid)   )
     
      allocate( p_coral (2,LBi:UBi, LBj:UBj) , &
     &          p_sgrass(LBi:UBi, LBj:UBj) , &
     &          p_sand  (LBi:UBi, LBj:UBj) , &
     &          p_algae (LBi:UBi, LBj:UBj)   )   


!-----------------------------------------------------------------------
!  Set initial conditions
!-----------------------------------------------------------------------

      do j=LBj,UBj
        do i=LBi,UBi
          do k=1,N
          
            dz(i,j,k)=1.5d0 !(m)
          
            C(i,j,k,1,iTemp) = 27.0d0   !27.0d0 32.0d0
            C(i,j,k,1,iSalt) = 34.0d0

            C(i,j,k,1,iSedi) = 0.0d0    !Sediment concentration (g m-3) 0.e0, 1.e0

            C(i,j,k,1,iTIC_) = 1915.0d0     !DIC  (umol kg-1) 1915.0d0 1930.0d0, 1700.0d0, 2100.0d0, 2300.0d0, 2500.0d0, 2700.0d0, 3000.0d0 
            C(i,j,k,1,iTAlk) = 2232.0d0     !TA  (umol kg-1)  2232.0d0 2275.0d0, 2500.0d0, 2150.0d0, 2000.0d0, 1800.0d0, 1700.0d0, 1600.0d0 
!            C(i,j,k,1,iOxyg) = 200.0d0      !DO  (umol L-1)
            C(i,j,k,1,iOxyg) = O2satu(C(i,j,k,1,iTemp)+273.15d0, C(i,j,k,1,iSalt))
#if defined ORGANIC_MATTER
            C(i,j,k,1,iDOC_) = 65.0d0       !DOC  (umol L-1) 
            C(i,j,k,1,iPOC_) = 4.2d0        !POC  (umol L-1) 
            C(i,j,k,1,iPhyt1) =  2.23d0       !Phytoplankton1 0.3(umolC L-1) all0.561 0.746-gC/gchla-1h 2.23-gC/gchla-30 4.47-gC/gchla-60	others
            C(i,j,k,1,iPhyt2) =  2.23d0       !Phytoplankton2 0.3(umolC L-1) 							diatom
            C(i,j,k,1,iZoop) =  1.3d0       !Zooplankton (umol L-1)1.3
#endif
#if defined CARBON_ISOTOPE
            R13C=R13C_fromd13C(0.7d0)
!            Ci(nDI13C,i,j,k)=R13C/(1.+R13C)*Ci(nDIC,i,j,k) !DI13C (umol kg-1)
            C(i,j,k,1,iT13C) =R13C*C(i,j,k,1,iTIC_) !DI13C (umol kg-1) 
#endif
#if defined NUTRIENTS            
!!!!!!!!!!!!
            C(i,j,k,1,iNO3_) =  0.2d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
!            C(i,j,k,1,iNO3_) =  1.8d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0
!						 C(i,j,k,1,iNO3_) =  3.5d0       !NO3  (umol L-1)  0.5d0, 10.0d0 control0.2d0 N1.8d0 N*2 3.5d0

            C(i,j,k,1,iPO4_) =  0.04d0      !PO4  (umol L-1) 0.03d0 0.05d0, 2.0d0 control0.04d0 P0.2d0
!            C(i,j,k,1,iPO4_) =  0.2d0      !PO4  (umol L-1) 0.03d0 0.05d0, 2.0d0 control0.04d0 P0.2d0
!!!!!!!!!!!!
            C(i,j,k,1,iNO2_) =  0.02d0      !NO2  (umol L-1) 
            C(i,j,k,1,iNH4_) =  0.26d0       !NH4  (umol L-1) 

# if defined ORGANIC_MATTER
            C(i,j,k,1,iDON_) =  10.0d0 !3.0d0       !DON  (umol L-1) 
            C(i,j,k,1,iPON_) =  0.6d0 !0.05d0      !PON  (umol L-1) 
            C(i,j,k,1,iDOP_) =  0.6d0 !0.1d0       !DOP  (umol L-1) 
            C(i,j,k,1,iPOP_) =  0.04d0 !0.015d0     !POP  (umol L-1) 
# endif
#endif
#if defined COT_STARFISH
            C(i,j,k,1,iCOTe) =  0.0d0       !COTS eggs (umol L-1) 
            C(i,j,k,1,iCOTl) =  0.0d0       !COTS larvae (umol L-1) 
#endif
          enddo
          
          p_coral(1,i,j) = 1.0d0
          p_coral(2,i,j) = 0.0d0
          p_sgrass(i,j)= 1.0d0
          p_sand(i,j)  = 1.0d0
          p_algae(i,j) = 0.0d0
          
          
        enddo
      enddo

      return
    end subroutine initialize_params

  END MODULE mod_param

