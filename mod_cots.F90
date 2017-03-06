
!!!=== ver 2013/11/01   Copyright (c) 2012-2013 Takashi NAKAMURA  =====

#include "cppdefs.h"


!!!**** MODULE OF COTS (Crown-of-thorns starfish) MODEL ****************

  module mod_cots

    implicit none

    integer, parameter :: Nstage = 5    !! Number of COTS life stage

    real(8), allocatable, save :: COTS(:,:,:)       ! COTS density (Individual m-2)
    real(8), allocatable, save :: p_coral(:,:,:)       ! COTS density (Individual m-2)

  contains


!!! **********************************************************************
!!!  set initial conditions for COTS (Crown-of-thorns starfish) model
!!! **********************************************************************

    subroutine initialize_cots(LBi, UBi, LBj, UBj, N)

      implicit none
! input parameters
      integer, intent(in) :: LBi, UBi, LBj, UBj, N

      integer i,j,k


      allocate( COTS(Nstage,LBi:UBi,LBj:UBj)      &
     &          p_coral(LBi:UBi,LBj:UBj)    ) !!!!!!!!!!!!!!!!!!

!----------set data -----------------------

!  Set initial conditions

      do j=LBj,UBj
        do i=LBi,UBi
          COTS (1,i,j)=0.d0 ! Life stage 1: 0.5-5.5 month (corallin algae eater)
          COTS (2,i,j)=0.d0 ! Life stage 2: 5.5 month - 1 year (coral eater)
          COTS (3,i,j)=0.d0 ! Life stage 3: 1-2 years (coral eater)
          COTS (4,i,j)=0.d0 ! Life stage 4: 2-3 years (coral eater)
          COTS (5,i,j)=0.d0 ! Life stage 5: 3- years (coral eater, adult stage)
          p_coral(i,j)=0.d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
      enddo

      return
    end subroutine initialize_cots
      

!!! **********************************************************************
!!!  Main program of COTS (Crown-of-thorns starfish) model 
!!! **********************************************************************

    subroutine cots               &
!          input parameters
     &            (N              &   ! Number of vertical grid (following ROMS vertical grid)
     &            ,i,j            &   ! i,j: position
     &            ,dt             &   ! Time step (sec)
     &            ,dx             &   ! dx: x grid size (m)
     &            ,dy             &   ! dy: y grid size (m)
     &            ,dz             &   ! dz(N): vertical grid size (m)
     &            ,LARcots        &   ! COTS larvae (individual m-3)
!     &            ,Tamb           &   ! Temperature (oC)
!     &            ,Samb           &   ! Salinity (PSU)
!     &            ,PHYamb         &   ! phytoplankton (umol C L-1)
!     &            ,tau_amb        &   ! bottom shear stress (N m-2)
!          output parameters
     &            ,EGGuptake      &   ! Egg uptake rate (individual m-2 s-1)  * direction of water column to coral is positive
!     &            ,DICuptake      &   ! DIC uptake rate (mmol m-2 s-1)  * direction of water column to coral is positive
!     &            ,DOuptake       &   ! DO  uptake rate (mmol m-2 s-1)  * direction of water column to coral is positive
!     &            ,PHYuptake      &   ! Phytoplankton ingestion rate (nmol cm-2 s-1)  * direction of water column to coral is positive
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
      integer, intent(in) :: N          
      integer, intent(in) :: i,j        
      real(8), intent(in) :: dt
      real(8), intent(in) :: dx      
      real(8), intent(in) :: dy      
      real(8), intent(in) :: dz(N)      
      real(8), intent(in) :: LARcots   
!      real(8), intent(in) :: Tamb     
!      real(8), intent(in) :: Samb     
!      real(8), intent(in) :: PHYamb   
!      real(8), intent(in) :: tau_amb  
! output parameters
      real(8), intent(out) :: EGGuptake
!      real(8), intent(out) :: DICuptake
!      real(8), intent(out) :: DOuptake 
!      real(8), intent(out) :: PHYuptake

!!!------------Set parameters  ----------------------------------

!----- Physical constants ------------------------
      real(8), parameter :: Rgas=8.314d0 ! Gas constant (J mol-1)


!--- Other variables ----------------------------------------------
      real(8) ::  Umax(Nstage)
      real(8) ::  Diff(Nstage)
      real(8) ::  Mort(Nstage)
      
! COTS maximum moving (m s-1)
      Umax(1)=1.d0
      Umax(2)=2.d0
      Umax(3)=3.d0
      Umax(4)=4.d0
      Umax(5)=5.d0

! COTS random movement coefficients (m2 s-1)
      Diff(1)=1.d0
      Diff(2)=2.d0
      Diff(3)=3.d0
      Diff(4)=4.d0
      Diff(5)=5.d0

! COTS mortality rate (Individual s-1)
      Mort(1)=1.d0
      Mort(2)=2.d0 * (1.d0 - p_coral(i,j))
      Mort(3)=3.d0 * (1.d0 - p_coral(i,j))
      Mort(4)=4.d0 * (1.d0 - p_coral(i,j))
      Mort(5)=5.d0 * (1.d0 - p_coral(i,j))

!  Output
      integer, save :: t_day = 0.d0  !day after spawning

! ---- COTS dynamics -----------------------------------------------

      do n=1,Nstage
      
        COTS(n,i,j)=COTS(n,i,j)+(                               &
          ! Advection term (upstream difference scheme)
          !   x axis
     &       -( 0.5d0*Umax(n)*ABS(p_coral(i+1,j)-p_coral(i,j))  &
     &         +0.5d0*Umax(n)*   (p_coral(i+1,j)-p_coral(i,j))) &
     &         *(COTS(n,i+1,j)-COTS(n,i,j))/dx                  &
     &       +( 0.5d0*Umax(n)*ABS(p_coral(i,j)-p_coral(i-1,j))  &
     &         +0.5d0*Umax(n)*   (p_coral(i,j)-p_coral(i-1,j))) &
     &         *(COTS(n,i,j)-COTS(n,i-1,j))/dx                  &
          !   y axis
     &       -( 0.5d0*Umax(n)*ABS(p_coral(i,j+1)-p_coral(i,j))  &
     &         +0.5d0*Umax(n)*   (p_coral(i,j+1)-p_coral(i,j))) &
     &         *(COTS(n,i,j+1)-COTS(n,i,j))/dy                  &
     &       +( 0.5d0*Umax(n)*ABS(p_coral(i,j)-p_coral(i,j-1))  &
     &         +0.5d0*Umax(n)*   (p_coral(i,j)-p_coral(i,j-1))) &
     &         *(COTS(n,i,j)-COTS(n,i,j-1))/dy                  &
          ! Diffusion term
     &       +Diff(n)*(COTS(n,i+1,j)-2.d0*COTS(n,i,j)+COTS(n,i-1,j))/dx/dx  &
     &       +Diff(n)*(COTS(n,i,j+1)-2.d0*COTS(n,i,j)+COTS(n,i,j-1))/dy/dy  &
          ! Reaction term
     &       -Mort(n)                                      &
     &       ) *dt
     
      end do

! ---- Coral dynamics -----------------------------------------------

      p_coral(i,j)=p_coral(i,j)+(                               &
     &       -Mort(n)                                           &
     &       ) *dt



!--- Change each life stage ------ --------------------------------
      if (t_day=15) then
        COTS(1,i,j)=LARcots
      end if
      if (t_day=150) then
        COTS(2,i,j)=COTS(1,i,j)
        COTS(1,i,j)=0.d0
      end if
      if (t_day=365) then
        COTS(5,i,j)=COTS(5,i,j)+COTS(4,i,j)
        COTS(4,i,j)=COTS(3,i,j)
        COTS(3,i,j)=COTS(2,i,j)
        COTS(2,i,j)=0.d0
        t_day=0.d0
        EGGuptake=-COTS(5,i,j)*1.0d6
      end if
      
!------------------------------------------------------------------------
! Print section (for debug)


       return
     end subroutine cots

      
  end module mod_cots

