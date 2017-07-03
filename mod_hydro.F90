
!!!=== ver 2017/05/13   Copyright (c) 2017 Takashi NAKAMURA  =====

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

      real(8), pointer :: Urc(:,:) 
      real(8), pointer :: Uch(:,:) 
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
      
      allocate( REEF(ng)%Urc(LBi:UBi,LBj:UBj)  )
      allocate( REEF(ng)%Uch(LBi:UBi,LBj:UBj)  )
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
          
          REEF(ng)%Urc(i,j)=0.0d0   ! Velocity through the reef crest (m s-1)
          REEF(ng)%Uch(i,j)=0.0d0   ! Velocity through the channel (m s-1)
          REEF(ng)%el (i,j)=0.0d0   ! Seawater elevation inside the reef (m)
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
     &            ,Qrc            &   ! volume flux through the reef crest (m3 s-1 m-2)
     &            ,Qch            &   ! volume flux through the channel (m3 s-1 m-2)
     &             )
!
!-----------------------------------------------------------------------
!                                                                       
!                   ___________ _____                                   
!                  |           :     :                                  
!                  |           :CHAN-:Wch                               
!                  |           : NEL<-- Uch,Qch                         
!                  |           : hch :                                  
!        SHORELINE |           :_____:                                  
!                  |   INNER   |     |     OFFSHORE                     
!               Wir|   REEF    |REEF |                                  
!                  |           |CREST|Wrc                               
!                  |    hir,   |    <-- Urc,Qrc                         
!                  |    el     | hcr |                                  
!                  |           |     |                                  
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
      real(8), intent(out) :: Qrc
      real(8), intent(out) :: Qch

!-----------------------------------------------------------------------
      real(8), parameter :: rho = 1024.0d0   ! Seawater density (kg m-3)
      real(8), parameter :: g   = 9.80665d0  ! Gravitational acceleration (m s-2)
      real(8), parameter :: pi  = 3.14159265359d0  ! Circle ratio

      real(8) :: Hs_rc         ! Significant wave hight on the reef crest (m)
      real(8) :: Hs_ch         ! Significant wave hight on the reef crest (m)
      real(8) :: Ew_rc         ! Wave energy on the reef crest (J)
      real(8) :: Ew_ch         ! Wave energy on the reef crest (J)
      real(8) :: Ew_o          ! Wave energy in offshore (J)
      real(8) :: Sxx_rc        ! Radiation stress in the reef (N m-1)
      real(8) :: Sxx_ch        ! Radiation stress in the reef (N m-1)
      real(8) :: Sxx_o         ! Radiation stress in offshore (N m-1)
      real(8) :: d             ! Water depth on the reef crest (m)
      real(8) :: del           ! difference of elevation between inner and outer reef (m)
      real(8) :: dSxx_rc       ! difference of radiation stress between inner and outer reef (N m-2)
      real(8) :: dSxx_ch       ! difference of radiation stress between inner and outer reef (N m-2)
            
!-----------------------------------------------------------------------

     d    = REEF(ng)%el(i,j) + REEF(ng)%hrc(i,j) 
     del  = REEF(ng)%el(i,j) - el_o
     
     Ew_o = 0.125d0*rho*g*Hs_o*Hs_o
     Sxx_o = Ew_o * (2.0d0*k*d/(sinh(2*k*d)) + 0.5)
     
     if ( el_o < -REEF(ng)%hrc(i,j) ) then
       REEF(ng)%Urc = 0.0d0
     else
       Hs_i = 0.8d0*d
       Sxx_i = Ew_i * (2.0d0*k*d/(sinh(2*k*d)) + 0.5)
       Ew_i = 0.125d0*rho*g*Hs_i*Hs_i
       
       dSxx  = Sxx_i - Sxx_o
     ! Momentum equation
       REEF(ng)%Urc = sqrt(-g*d*del - d*dSxx/rho)
     endif
     

    END SUBROUTINE reef_hydro
    
! **********************************************************************
!  Momentum equation
! **********************************************************************

    SUBROUTINE momentum_eq             &
!          input parameters
     &            (dt             &   ! Time step (sec)
     &            ,Hs_o           &   ! Significant wave hight at offshore (m)
     &            ,k              &   ! Wavenumber
     &            ,el_o           &   ! Sea surface elevation at offshore (m)
     &            ,el_i           &   ! Sea surface elevation of inside (m)
     &            ,h_i            &   ! Depth of inside (m)

!          output parameters
     &            ,U              &   ! Velocity (o->i is positive, m s-1)
     &             )
!-----------------------------------------------------------------------
!
      implicit none

! input parameters
      real(8), intent(in ) :: dt         
      real(8), intent(in ) :: Hs_o
      real(8), intent(in ) :: k
      real(8), intent(in ) :: el_o
      real(8), intent(in ) :: el_i
      real(8), intent(in ) :: h_i
! output parameters
      real(8), intent(out) :: U

!-----------------------------------------------------------------------
      real(8), parameter :: rho = 1024.0d0   ! Seawater density (kg m-3)
      real(8), parameter :: g   = 9.80665d0  ! Gravitational acceleration (m s-2)
      real(8), parameter :: pi  = 3.14159265359d0  ! Circle ratio

      real(8) :: Hs_i          ! Significant wave hight on the reef crest (m)
      real(8) :: Ew_i          ! Wave energy on the reef crest (J)
      real(8) :: Ew_o          ! Wave energy in offshore (J)
      real(8) :: Sxx_i         ! Radiation stress in the reef (N m-1)
      real(8) :: Sxx_o         ! Radiation stress in offshore (N m-1)
      real(8) :: d             ! Water depth on the reef crest (m)
      real(8) :: del           ! difference of elevation between inner and outer reef (m)
      real(8) :: dSxx          ! difference of radiation stress between inner and outer reef (N m-2)
            
!-----------------------------------------------------------------------

     d    = REEF(ng)%el(i,j) + REEF(ng)%hrc(i,j) 
     del  = REEF(ng)%el(i,j) - el_o
     
     Ew_o = 0.125d0*rho*g*Hs_o*Hs_o
     Sxx_o = Ew_o * (2.0d0*k*d/(sinh(2*k*d)) + 0.5)
     
     if ( el_o < -REEF(ng)%hrc(i,j) ) then
       REEF(ng)%Urc = 0.0d0
     else
       Hs_i = 0.8d0*d
       Sxx_i = Ew_i * (2.0d0*k*d/(sinh(2*k*d)) + 0.5)
       Ew_i = 0.125d0*rho*g*Hs_i*Hs_i
       
       dSxx  = Sxx_i - Sxx_o
     ! Momentum equation
       REEF(ng)%Urc = sqrt(-g*d*del - d*dSxx/rho)
     endif
     

    END SUBROUTINE moment

  END MODULE mod_hydro


