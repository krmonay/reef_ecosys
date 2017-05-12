
!!!=== ver 2017/05/12   Copyright (c) 2017 Takashi NAKAMURA  =====

#include "cppdefs.h"


!!!**** box model MODULE ************************************

  MODULE mod_hydro

    implicit none

    TYPE T_REEF

      real(8), pointer :: Wir(:,:)
      real(8), pointer :: hir(:,:)
      real(8), pointer :: Lir(:,:)
      real(8), pointer :: Wrc(:,:)
      real(8), pointer :: hrc(:,:)
      real(8), pointer :: Lrc(:,:)
      real(8), pointer :: hch(:,:)

      real(8), pointer :: el(:,:) 

    END TYPE T_REEF

    TYPE (T_REEF), allocatable :: REEF(:)

  CONTAINS

!!! **********************************************************************
!!!  set initial conditions for reef hydrodynamics
!!! **********************************************************************

    SUBROUTINE initialize_hydro(ng, Ngrids, LBi, UBi, LBj, UBj)

      implicit none
! input parameters
      integer, intent(in) :: ng, Ngrids, LBi, UBi, LBj, UBj
      integer i,j,n

      IF (ng.eq.1) allocate ( REEF(Ngrids) )

      allocate( REEF(ng)%Wir(LBi:UBi,LBj:UBj)  )
      allocate( REEF(ng)%hir(LBi:UBi,LBj:UBj)  )
      allocate( REEF(ng)%Lir(LBi:UBi,LBj:UBj)  )
      allocate( REEF(ng)%Wrc(LBi:UBi,LBj:UBj)  )
      allocate( REEF(ng)%hrc(LBi:UBi,LBj:UBj)  )
      allocate( REEF(ng)%Lrc(LBi:UBi,LBj:UBj)  )
      allocate( REEF(ng)%hch(LBi:UBi,LBj:UBj)  )
      allocate( REEF(ng)%el (LBi:UBi,LBj:UBj)  )

!  Set initial conditions
      do j=LBj,UBj
        do i=LBi,UBi
          REEF(ng)%Wir(i,j)=1.0d3   ! Inner reef width (m)
          REEF(ng)%hir(i,j)=3.0d0   ! Inner reef depth (m)
          REEF(ng)%Lir(i,j)=1.0d3   ! Inner reef length (m)
          REEF(ng)%Wrc(i,j)=0.9d0   ! Reef crest width (m)
          REEF(ng)%hrc(i,j)=0.5d0   ! Reef crest depth (m)
          REEF(ng)%Lrc(i,j)=0.0d0   ! Reef crest length (m)
          REEF(ng)%hrc(i,j)=2.0d0   ! Channel depth (m)
          
          REEF(ng)%el (i,j)=0.0d0
        enddo
      enddo

      RETURN
    END SUBROUTINE initialize_hydro


! **********************************************************************
!  Heat and mass balance of water column
! **********************************************************************

    SUBROUTINE reef_hydro         &
!          input parameters
     &            (ng, i, j       &   ! ng: nested grid number; i,j: position
     &            ,N              &   ! Number of vertical grid (following ROMS vertical grid)
     &            ,dt             &   ! Time step (sec)
     &            ,Hs_o           &   ! Significant wave hight at offshore (m)
     &            ,k              &   ! Wavenumber
     &            ,el_o           &   ! Sea surface elevation at offshore (m)

!          output parameters
     &            ,Qrc            &   ! volume flux through the reef crest (m3 s-1)
     &            ,Qch            &   ! volume flux through the channel (m3 s-1)
     &             )
!
!-----------------------------------------------------------------------
!                                                                       
!                   ___________ _____                                   
!                  |           :     :                                  
!                  |           :Chan-:                                  
!                  |           : nel :Wch = Wir-Wch                     
!                  |           : hch :                                  
!        Shoreline |           :_____:                                  
!                  |   inner   |     |     Offshere                     
!               Wir|   reef    |     |                                  
!                  |           |Reef |                                  
!                  |    hir    |crest|Wrc                               
!                  |           |     |                                  
!                  |           | hcr |                                  
!                  |___________|____ |                                  
!                        Lir     Lrc                                    
!                                                                       
!      Plan view of an idealized reef-channel system.                   
!-----------------------------------------------------------------------
!
      implicit none

! input parameters
      integer, intent(in) :: ng, i, j        
      integer, intent(in) :: N          
      real(8), intent(in) :: dt         
      real(8), intent(in) :: Hs_o
      real(8), intent(in) :: k
      real(8), intent(in) :: el_o
! output parameters
      real(8), intent(out) :: Qin
      real(8), intent(out) :: Qout

!-----------------------------------------------------------------------
      real(8), parameter :: rho = 1024.0d0   ! Seawater density (kg m-3)
      real(8), parameter :: g   = 9.80665d0  ! Gravitational acceleration (m s-2)
      real(8), parameter :: pi  = 3.14159265359d0  ! Circle ratio

      real(8) :: Hs_i          ! Hignificant wave hight in the reef (m)
      real(8) :: Ew_i          ! Wave energy in the reef (J)
      real(8) :: Ew_o          ! Wave energy in offshore (J)
      real(8) :: Sxx_i         ! Radiation stress in the reef (N m-1)
      real(8) :: Sxx_o         ! Radiation stress in offshore (N m-1)
      real(8) :: d             ! Water depth on the reef crest (m)
      real(8) :: del           ! difference of elevation between inner and outer reef (m)
      real(8) :: dSxx          ! difference of radiation stress between inner and outer reef (N m-2)
            
!-----------------------------------------------------------------------

     d    = REEF(ng)%el(i,j) + REEF(ng)%hrc(i,j) 
     del  = REEF(ng)%el(i,j) - el_o
     
     Hs_i = 0.8d0*d
     
     Ew_o = 0.125d0*rho*g*Hs_o*Hs_o
     Ew_i = 0.125d0*rho*g*Hs_i*Hs_i
     
     Sxx_o = Ew_o * (2.0d0*k*d/(sinh(2*k*d)) + 0.5)
     Sxx_i = Ew_i * (2.0d0*k*d/(sinh(2*k*d)) + 0.5)
     dSxx  = Sxx_i - Sxx_o
     
     if ( el_o < -REEF(ng)%hrc(i,j) ) then
       Qrc = 0.0d0
     else
     ! Momentum equation
       Qrc = sqrt(-g*d*del - d*dSxx/rho)
     endif
     

    END SUBROUTINE reef_hydro

  END MODULE mod_hydro


