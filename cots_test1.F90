
!!!=== ver 2013/11/01   Copyright (c) 2012-2013 Takashi NAKAMURA  =====

#include "cppdefs.h"


      PROGRAM cots_test
! **********************************************************************
! *                                                                    *
! *   Test program of cots.F                                           *
! *                                                                    *
! **********************************************************************
!
      USE mod_cots

      implicit none
      
      integer, parameter :: Im = 200
      integer, parameter :: Jm = 200
      integer, parameter :: N  = 1

      real(8), parameter :: dt = 1.0d0
      integer, parameter :: Nburst  = 20
      
      integer, parameter :: iEgg = 1      ! ID for COTS egg
      integer, parameter :: iLar = 2      ! ID for COTS larvae

      real(8), parameter :: Diff = 1.0d0  !!! Eddy diffusion coefficient

      real(8) :: h(0:Im,0:Jm)      ! Depth (m)   
      real(8) :: mask_rho(0:Im,0:Jm)      ! mask for rho points (m)   
      real(8) :: mask_u(1:Im,0:Jm)      ! mask for u points (m)   
      real(8) :: mask_v(0:Im,1:Jm)      ! mask for v points (m)   

      real(8) :: dx  = 100.d0  ! X Grid size (m)
      real(8) :: dy  = 100.d0  ! Y Grid size (m)   
      real(8) :: dz(Im,Jm,N)  = 1.d0   
      real(8) :: C(2,0:Im,0:Jm) = 0.d0
      real(8) :: dC_dt(2,0:Im,0:Jm) = 0.d0
      real(8) :: p_coral(0:Im,0:Jm) = 0.d0

      integer :: i,j,k,id, Nid
      integer :: istep, iburst, iprint
      real(8) :: time    !(day)
!  For Output      
      real(8), save :: dsec = 0.d0 !sec


!      open(52,file='./output/cots_his.txt')
!      open(53,file='./output/eco4-crl_ave.txt')!!!!!!!!!!!!!!!!!!!for debug
!      open(54,file='./output/eco4-zoo_his.txt')!!!!!!!!!!!!!!!!!!!for debug
!      open(55,file='./output/eco4-box_his.txt')!!!!!!!!!!!!!!!!!!!for debug

      istep=0
      iprint=0
      

!----- Set parameters -------------------------

      do j=0,Jm
        do i=0, Im
        
          if(h(i,j)<=0.d0) then
            mask_rho(i,j)=0
          end if
        
        end do
      end do
      do j=0,Jm
        do i=1, Im
          mask_u(i,j)=mask_rho(i-1,j)*mask(i,j)
        end do
      end do
      do j=1,Jm
        do i=0, Im
          mask_v(i,j)=mask_rho(i,j-1)*mask(i,j)
        end do
      end do
     
!----- Set initial conditions -------------------------

      CALL initialize_cots(1, Im, 1, Jm, N)

      time=0.

!----- Main loop -------------------------------------------

      do istep=1, int(24.*60.*60./dt)*365 * 5 +1      ! 5 years

!-----------------------------------------------------------
!    COTS model
!-----------------------------------------------------------

        CALL cots                 &
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

!-----------------------------------------------------------


!-----Reaction diffusion of egg and larvae -----------------

        if (t_day<=15d0) then      ! 15 days

          do iburst=1, Nburst

            if (t_day<=0.5d0) then
              do j=1,Jm-1
                do i=1, Im-1
                  dC_dt(iEGG,i,j)=-EGGuptake/h(i,j)  ! (individual m-3 s-1)
                enddo
              enddo
            end if
            
            cff=1.d-5
            if (t_day<=2.0d0) then
              do j=1,Jm-1
                do i=1, Im-1
                  dC_dt(iLAR,i,j)=cff*C(iEGG,i,j)*C(iEGG,i,j)   & ! (individual m-3 s-1)
    &                             -1.d0 *C(iLAR,i,j)   !!! mortality (individual m-3 s-1)
                  dC_dt(iEGG,i,j)=-2.d0 * dC_dt(iLAR,i,j)       &
    &                             -1.d0 *C(iEGG,i,j)   !!! mortality (individual m-3 s-1)
                enddo
              enddo
            end if


            do id=1,2
      
      !     Cycle boundary condition
              do i=1,Im-1
                C(id,i,0 )=C(id,i,Jm-1)
                C(id,i,Jm)=C(id,i,1)
              enddo
              do j=1,Jm-1
                C(id,0, j)=C(id,Im-1,j)
                C(id,Im,j)=C(id,1,j)
              enddo
      
      !     Diffusion equation
              do j=1,Jm-1
                do i=1, Im-1

                  C(id,i,j)=C(id,i,j)+(                                     &
    &                Diff*( 0.5d0*(h(i+1,j)+h(i,j))*(C(id,i+1,j)-C(id,i,j))*mask_u(i+1,j)  &
    &                      -0.5d0*(h(i,j)+h(i-1,j))*(C(id,i,j)-C(id,i-1,j))*mask_u(i,j)    &
    &                     )/dx/dx                                            &
    &               +Diff*( 0.5d0*(h(i,j+1)+h(i,j))*(C(id,i,j+1)-C(id,i,j))*mask_v(i,j+1)  &
    &                      -0.5d0*(h(i,j)+h(i,j-1))*(C(id,i,j)-C(id,i,j-1))*mask_v(i,j)    &
    &                     )/dy/dy                                            &
    &               )/h(i,j)*dt                                              &
    &              +dC_dt(id,i,j) *dt

                enddo
              enddo
            
            enddo ! End of id
            
          enddo ! END of iburst
        endif

!-----------------------------------------------------------



      enddo
      
!----- End loop --------------------------------------

      end PROGRAM cots_test
!----------------------------------------------------------------------!

!     End of main program

!-----------------------------------------------------------------------

