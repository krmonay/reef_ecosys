
!!!=== ver 2013/12/03   Copyright (c) 2012-2013 Takashi NAKAMURA  =====

#include "cppdefs.h"


!!!**** MODULE OF FOOD WEB *******************************************

  MODULE mod_foodweb

    implicit none

  CONTAINS
  
!!! **********************************************************************
!!!  Main program of foodweb model (Modified from Yamamoto et al. under review)
!!! **********************************************************************

    SUBROUTINE foodweb            &
!          input parameters
     &            (ng, n, i, j    &   ! ng: nested grid number; n: coral compartment; i,j: position
     &            ,dt             &   ! Time step (sec)
     &            ,PFD            &   ! Photon flux density (umol m-2 s-1)
     &            ,rho_sw         &   ! Density of seawater (g cm-3)
     &            ,Tmp            &   ! Temperature (oC)
     &            ,Sal            &   ! Salinity (PSU)
     &            ,DIC            &   ! Total dissolved inorganic carbon (DIC: umol kg-1)
     &            ,TA             &   ! Total alkalinity (TA: umol kg-1)
     &            ,DOx            &   ! Dissolved oxygen (umol L-1)
     &            ,DOC            &   ! Dissolved organic carbon (DOC: umol L-1)
     &            ,POC            &   ! Particulate organic carbon (DOC: umol L-1)
     &            ,PHY            &   ! phytoplankton (umol C L-1)
     &            ,ZOO            &   ! zooplankton (umol C L-1)
#if defined CARBON_ISOTOPE
     &            ,DI13C          &   !13C of DIC (umol kg-1)
#endif
#if defined NUTRIENTS         
     &            ,NO3            &   ! NO3 (umol L-1)
     &            ,NO2            &   ! NO2 (umol L-1)
     &            ,NH4            &   ! NH4 (umol L-1)
     &            ,PO4            &   ! PO4 (umol L-1)
     &            ,DON            &   ! Dissolved organic nitrogen (DON: umol L-1)
     &            ,PON            &   ! Particulate organic nitrogen (PON: umol L-1)
     &            ,DOP            &   ! Dissolved organic phosporius (DOP: umol L-1)
     &            ,POP            &   ! Particulate organic phosporius (POP: umol L-1)
#endif
#if defined COT_STARFISH         
     &            ,COTe           &   ! COT starfish egg (umol L-1)
     &            ,COTl           &   ! COT starfish larvae (umol L-1)
#endif
!          output parameters
     &            ,dDIC_dt        &   ! dDIC/dt  (umol kg-1 s-1)  1 mmol m-3 = 1 umol L-1 = 1/1.024 umol kg-1
     &            ,dTA_dt         &   ! dTA/dt   (umol kg-1 s-1) 
     &            ,dDOx_dt         &  ! dDOx/dt  (umol L-1 s-1) 
     &            ,dDOC_dt        &   ! dDOC/dt  (umol L-1 s-1) 
     &            ,dPOC_dt        &   ! dPOC/dt  (umol L-1 s-1) 
     &            ,dPHY_dt        &   ! dPHY/dt  (umol L-1 s-1)  
     &            ,dZOO_dt        &   ! dZOO/dt  (umol L-1 s-1)  
#if defined CARBON_ISOTOPE
     &            ,dDI13C_dt      &   ! dDI13C/dt (umol kg-1 s-1)
#endif
#if defined NUTRIENTS
     &            ,dNO3_dt        &   ! dNO3/dt (umol L-1 s-1)
     &            ,dNO2_dt        &   ! dNO2/dt (umol L-1 s-1)
     &            ,dNH4_dt        &   ! dNH4/dt (umol L-1 s-1)
     &            ,dPO4_dt        &   ! dPO4/dt (umol L-1 s-1)
     &            ,dDON_dt        &   ! dDON/dt (umol L-1 s-1)
     &            ,dPON_dt        &   ! dPON/dt (umol L-1 s-1)
     &            ,dDOP_dt        &   ! dDOP/dt (umol L-1 s-1)
     &            ,dPOP_dt        &   ! dPOP/dt (umol L-1 s-1)
#endif                                
#if defined COT_STARFISH         
     &            ,dCOTe_dt       &   ! dCOTe/dt (umol L-1 s-1)
     &            ,dCOTl_dt       &   ! dCOTl/dt (umol L-1 s-1)
#endif
     &             )

!-----------------------------------------------------------------------
!
      implicit none

!          input parameters
      integer, intent(in) :: ng, n, i, j    ! ng: nested grid number; n: coral compartment; i,j: position
      real(8), intent(in) :: dt             ! Time step (sec)
      real(8), intent(in) :: PFD            ! Photon flux density (umol m-2 s-1)
      real(8), intent(in) :: rho_sw         ! Density of seawater (g cm-3)
      real(8), intent(in) :: Tmp            ! Temperature (oC)
      real(8), intent(in) :: Sal            ! Salinity (PSU)
      real(8), intent(in) :: DIC            ! Total dissolved inorganic carbon (DIC: umol kg-1)
      real(8), intent(in) :: TA             ! Total alkalinity (TA: umol kg-1)
      real(8), intent(in) :: DOx            ! Dissolved oxygen (umol L-1)
      real(8), intent(in) :: DOC            ! Dissolved organic carbon (DOC: umol L-1)
      real(8), intent(in) :: POC            ! Particulate organic carbon (DOC: umol L-1)
      real(8), intent(in) :: PHY            ! phytoplankton (umol C L-1)
      real(8), intent(in) :: ZOO            ! zooplankton (umol C L-1)
#if defined CARBON_ISOTOPE
      real(8), intent(in) :: DI13C          !13C of DIC (umol kg-1)
#endif
#if defined NUTRIENTS         
      real(8), intent(in) :: NO3            ! NO3 (umol L-1)
      real(8), intent(in) :: NO2            ! NO2 (umol L-1)
      real(8), intent(in) :: NH4            ! NH4 (umol L-1)
      real(8), intent(in) :: PO4            ! PO4 (umol L-1)
      real(8), intent(in) :: DON            ! Dissolved organic nitrogen (DON: umol L-1)
      real(8), intent(in) :: PON            ! Particulate organic nitrogen (PON: umol L-1)
      real(8), intent(in) :: DOP            ! Dissolved organic phosporius (DOP: umol L-1)
      real(8), intent(in) :: POP            ! Particulate organic phosporius (POP: umol L-1)
#endif
#if defined COT_STARFISH         
      real(8), intent(in) :: COTe           ! COT starfish egg (umol L-1)
      real(8), intent(in) :: COTl           ! COT starfish larvae (umol L-1)
#endif
!          output parameters
      real(8), intent(out) :: dDIC_dt       ! dDIC/dt  (umol kg-1 s-1)  1 mmol m-3 = 1 umol L-1 = 1/1.024 umol kg-1
      real(8), intent(out) :: dTA_dt        ! dTA/dt   (umol kg-1 s-1) 
      real(8), intent(out) :: dDOx_dt       ! dDO/dt   (umol L-1 s-1) 
      real(8), intent(out) :: dDOC_dt       ! dDOC/dt  (umol L-1 s-1) 
      real(8), intent(out) :: dPOC_dt       ! dPOC/dt  (umol L-1 s-1) 
      real(8), intent(out) :: dPHY_dt       ! dPHY/dt  (umol L-1 s-1)  
      real(8), intent(out) :: dZOO_dt       ! dZOO/dt  (umol L-1 s-1)  
#if defined CARBON_ISOTOPE
      real(8), intent(out) :: dDI13C_dt     ! dDI13C/dt (umol kg-1 s-1)
#endif
#if defined NUTRIENTS
      real(8), intent(out) :: dNO3_dt       ! dNO3/dt (umol L-1 s-1)
      real(8), intent(out) :: dNO2_dt       ! dNO2/dt (umol L-1 s-1)
      real(8), intent(out) :: dNH4_dt       ! dNH4/dt (umol L-1 s-1)
      real(8), intent(out) :: dPO4_dt       ! dPO4/dt (umol L-1 s-1)
      real(8), intent(out) :: dDON_dt       ! dDON/dt (umol L-1 s-1)
      real(8), intent(out) :: dPON_dt       ! dPON/dt (umol L-1 s-1)
      real(8), intent(out) :: dDOP_dt       ! dDOP/dt (umol L-1 s-1)
      real(8), intent(out) :: dPOP_dt       ! dPOP/dt (umol L-1 s-1)
#endif
#if defined COT_STARFISH
      real(8), intent(out) :: dCOTe_dt      ! dCOTe/dt (umol L-1 s-1)
      real(8), intent(out) :: dCOTl_dt      ! dCOTl/dt (umol L-1 s-1)
#endif

!!!------------Set parameters  ----------------------------------

!------- Phytoplankton parameters ------------------------
      real(8), parameter :: k_Pphy =  1.0d0/86400.0d0     ! (s-1)          PHY maximum photosynthesis rate at 0 oC (1.0d0 d-1; Kawamiya et al., 1995)
      real(8), parameter :: b_Pphy =  0.063d0             ! (degC-1)       Temperature coefficient for PHY photosynthesis (Kawamiya et al., 1995)
      real(8), parameter :: Iphy   = 48.83d0*1.82d0       ! (umol m-2 s-1) PHY optimum light intensity (48.83d0 J m2 s-1; Kawamiya et al., 1995)
      real(8), parameter :: k_Rphy =  0.03d0/86400.0d0    ! (s-1)          PHY respiration rate at 0 oC (Kawamiya et al., 1995)
      real(8), parameter :: b_Rphy =  0.0519d0            ! (degC-1)       Temperature coefficient for PHY respiration rate (0.03d0 d-1; Kawamiya et al., 1995)
      real(8), parameter :: k_Mphy =  0.00281d0/86400.0d0 ! (umol-1 s-1)   PHY mortality rate at 0 oC (0.00281d0 umol-1 d-1; Kishi et al., 2001)
      real(8), parameter :: b_Mphy =  0.069d0             ! (degC-1)       Temperature coefficient for PHY mortality (Kawamiya et al., 1995)
      real(8), parameter :: k_Ephy =  0.135d0             ! (no dim.)      PHY ratio of extracellular excretion to production (Kawamiya et al., 1995)
#if defined NUTRIENTS         
      real(8), parameter :: K_NH4 = 0.2d0   ! (s-1)         
      real(8), parameter :: K_NO3 = 0.1d0   ! (s-1)         
      real(8), parameter :: K_PO4 = 0.01d0   ! (s-1)         
      real(8), parameter :: psi =  0.01d0    !(s-1)          
#endif
!------- Zooplankton parameters ------------------------
      real(8), parameter :: k_Gphy2zoo = 0.3d0/86400.0d0 ! (s-1)          Maximum grazing rate of PHY by ZOO (0.3d0 d-1; Kawamiya et al., 1995)
      real(8), parameter :: b_Gphy2zoo = 0.063d0         ! (degC-1)       Temperature coefficient for ZOO grazing (Kawamiya et al., 1995)
      real(8), parameter :: e_Gphy2zoo = 0.7d0           ! (no dim.)      Assimilation efficiency of ZOO (Kawamiya et al., 1995)
      real(8), parameter :: k_Rzoo = 0.005d0/86400.0d0   !!! (s-1)          ZOO respiration rate at 0  oC  !!!(Tuning)
      real(8), parameter :: b_Rzoo = 0.0693d0            ! (degC-1)       Temperature coefficient for ZOO respiration rate (Kawamiya et al., 1995)
      real(8), parameter :: k_Mzoo = 3.0d-1!!!/86400.0d0  !!! (umol-1 s-1)   ZOO mortality rate at 0 oC (0.0088d0 umol-1 d-1; Kawamiya et al., 1995)
      real(8), parameter :: b_Mzoo = 0.0693d0            ! (degC-1)       Temperature coefficient for ZOO mortality (Kawamiya et al., 1995)
!------- Microbial loop parameters --------------------
      real(8), parameter :: k_Gpoc2zoo = 0.0d0          ! (s-1)          PHY mortality rate at 0 oC (0.3d0 d-1; Kawamiya et al., 1995)
      real(8), parameter :: b_Gdoc2zoo = 0.0d0          ! (degC-1)       Temperature coefficient for ZOO grazing (Kawamiya et al., 1995)
      real(8), parameter :: k_Gdoc2zoo = 0.0d0          ! (s-1)          PHY mortality rate at 0 oC (0.3d0 d-1; Kawamiya et al., 1995)
      real(8), parameter :: b_Gpoc2zoo = 0.0d0          ! (degC-1)       Temperature coefficient for ZOO grazing (Kawamiya et al., 1995)
!------- Decomposition parameters --------------------
      real(8), parameter :: k_Ddoc = 0.03d0/86400.0d0   ! (umol L-1 s-1)          PHY mortality rate at 0 oC (0.3d0 d-1; Kishi et al., 2001)
      real(8), parameter :: b_Ddoc = 0.0693d0           ! (degC-1)       Temperature coefficient for ZOO grazing (Kishi et al., 2001)
      real(8), parameter :: k_Dpoc = 0.01d0/86400.0d0   ! (s-1)          PHY mortality rate at 0 oC (0.3d0 d-1; Kishi et al., 2001)
      real(8), parameter :: b_Dpoc = 0.0693d0           ! (degC-1)       Temperature coefficient for ZOO grazing (Kishi et al., 2001)
!------- Physical parameters --------------------
#if defined NUTRIENTS         
      real(8), parameter :: rCNphy = 106.0d0/16.0d0   ! (no dim.) PHY C:N ratio (Redfield ratio)
      real(8), parameter :: rCPphy = 106.0d0/1.0d0    ! (no dim.) PHY C:P ratio (Redfield ratio)
      real(8), parameter :: rCNzoo = 106.0d0/16.0d0   ! (no dim.) ZOO C:N ratio (Redfield ratio)
      real(8), parameter :: rCPzoo = 106.0d0/1.0d0    ! (no dim.) ZOO C:P ratio (Redfield ratio)
#endif                                

!------- Local variables --------------------
      real(8) :: Pphy, Rphy, Mphy, Ephy, Aphy
      real(8) :: Gphy2zoo, Rzoo, Mzoo
      real(8) :: Gdoc2zoo, Gpoc2zoo
      real(8) :: Ddoc, Dpoc, Dpoc2doc
#if defined NUTRIENTS
      real(8) :: Gdon2zoo, Gpon2zoo
      real(8) :: Gdop2zoo, Gpop2zoo
      real(8) :: Ddon, Dpon, Dpon2don
      real(8) :: Ddop, Dpop, Dpop2dop
      real(8) :: V_NH4, V_NO3
      real(8) :: V_PO4
#endif

!!!------- Phytoplankton reaction ------------------------

!----- Gross photosynthetic rate (umolC L-1 s-1) -----------------
!      Pphy = k_Pphy * exp(b_Pphy*Tmp) * PFD/Iphy*exp(1.0d0-PFD/Iphy) * PHY
      Pphy = k_Pphy * exp(b_Pphy*Tmp) * tanh(PFD/Iphy) * PHY
      IF(DIC <= 0.d0) THEN !-----For Error handling
        Pphy = 0.d0
      ENDIF

!----- Assimilation rate (umolC L-1 s-1) -----------------

#if defined NUTRIENTS         
      V_NH4 = NH4/(NH4+K_NH4)
      V_NO3 = NO3/(NO3+K_NO3) * exp(-psi * NH4)
      V_PO4 = PO4/(PO4+K_PO4)
      
      Aphy = MIN( V_NH4+V_NO3, V_PO4 ) ! Assimilation rate
      Aphy = Aphy * Pphy
#else
      Aphy = (1.0d0-k_Ephy) * Pphy
#endif

!----- Excretion rate (umolC L-1 s-1) -----------------
#if defined NUTRIENTS
      Ephy = Pphy - Aphy
#else
      Ephy = k_Ephy * Pphy
#endif

!----- Respiration rate (umolC L-1 s-1) -----------------
      Rphy = k_Rphy * exp(b_Rphy*Tmp) * PHY
      IF(DOx <= 0.d0) THEN !-----For Error handling
        Rphy = 0.d0
      ENDIF

!----- Mortality (umolC L-1 s-1) -----------------
      Mphy = k_Mphy * exp(b_Mphy*Tmp) * PHY*PHY
!      Mphy = k_Mphy * exp(b_Mphy*Tmp) * PHY

!!!------- Zooplankton reaction ------------------------

!----- Grazing rate of PHY by ZOO (umolC L-1 s-1) -----------------
      Gphy2zoo = k_Gphy2zoo * exp(b_Gphy2zoo*Tmp) * PHY * ZOO

!----- Respiration rate (umolC L-1 s-1) -----------------
      Rzoo = k_Rzoo * exp(b_Rzoo*Tmp) * ZOO
      IF(DOx <= 0.d0) THEN !-----For Error handling
        Rzoo = 0.d0
      ENDIF

!----- Mortality (umolC L-1 s-1) -----------------
!      Mzoo = k_Mzoo * exp(b_Mzoo*Tmp) * ZOO*ZOO
!      Mzoo = k_Mzoo * exp(b_Mzoo*Tmp) * ZOO
      Mzoo = k_Mzoo * exp(b_Mzoo*Tmp) * Gphy2zoo

!!!------- Microbial loop (inplicitly assumed) ---------------

!----- Grazing rate of DOC by ZOO (umolC L-1 s-1) -----------------
      Gdoc2zoo = k_Gpoc2zoo * exp(b_Gdoc2zoo*Tmp) * DOC * ZOO

!----- Grazing rate of POC by ZOO (umolC L-1 s-1) -----------------
      Gpoc2zoo = k_Gdoc2zoo * exp(b_Gpoc2zoo*Tmp) * POC * ZOO

!!!------- Decomposition ------------------------

!----- Decomposition rate of DOC (umolC L-1 s-1) -----------------
!      Ddoc = k_Ddoc * exp(b_Ddoc*Tmp) * DOC
      Ddoc = k_Ddoc * exp(b_Ddoc*Tmp) * DOC**2.0d0/(60.0d0**2.0d0+DOC**2.0d0)
      IF(DOx <= 0.d0) THEN !-----For Error handling
        Ddoc = 0.d0
      ENDIF
#if defined NUTRIENTS
      Ddon = k_Ddoc * exp(b_Ddoc*Tmp) * DON**2.0d0/(60.0d0**2.0d0+DON**2.0d0)
      Ddop = k_Ddoc * exp(b_Ddoc*Tmp) * DOP**2.0d0/(60.0d0**2.0d0+DOP**2.0d0)
#endif

!----- Decomposition rate of POC (umolC L-1 s-1) -----------------
      Dpoc = k_Dpoc * exp(b_Dpoc*Tmp) * POC**2.0d0/(60.0d0**2.0d0+POC**2.0d0)
      IF(DOx <= 0.d0) THEN !-----For Error handling
        Dpoc = 0.d0
      ENDIF
#if defined NUTRIENTS
      Dpon = k_Dpoc * exp(b_Dpoc*Tmp) * PON**2.0d0/(60.0d0**2.0d0+DOC**2.0d0)
      Dpop = k_Dpoc * exp(b_Dpoc*Tmp) * POP**2.0d0/(60.0d0**2.0d0+DOC**2.0d0)
#endif
!----- Decomposition rate from POC to DOC (umolC L-1 s-1) -----------------
      Dpoc2doc = 0.d0 !!!
#if defined NUTRIENTS
      Dpon2don = 0.d0 !!!
      Dpop2dop = 0.d0 !!!
#endif

!!!------- Mass barance equations ------------------------

      dPHY_dt = Aphy - Rphy - Mphy - Ephy - Gphy2zoo
      dZOO_dt = Gphy2zoo*e_Gphy2zoo + Gdoc2zoo + Gpoc2zoo - Rzoo -Mzoo
#if defined COT_STARFISH
      dCOTe_dt = 0.0d0
      dCOTl_dt = 0.0d0
#endif

      dDIC_dt = -Pphy + Rphy + Rzoo + Ddoc + Dpoc
      dDIC_dt = dDIC_dt/rho_sw !(umol L-1 s-1) -> (umol kg-1 s-1)
      dTA_dt  = 0.0d0
      dTA_dt  = dTA_dt/rho_sw !(umol L-1 s-1) -> (umol kg-1 s-1)
      dDOC_dt = Ephy - Gdoc2zoo - Ddoc + Dpoc2doc
      dPOC_dt = Mphy + Gphy2zoo*(1.0d0-e_Gphy2zoo) + Mzoo - Gpoc2zoo - Dpoc - Dpoc2doc

#if defined NUTRIENTS         
      dNO3_dt = -V_NO3/(V_NO3+V_NH4)*Aphy/rCNphy
      dNO2_dt = 0.0d0
      dNH4_dt = (-V_NH4/(V_NO3+V_NH4)*Aphy + Rphy)/rCNphy + Rzoo/rCNzoo + Ddon + Dpon
      dPO4_dt = (-Aphy + Rphy)/rCPphy + Rzoo/rCPzoo + Ddop + Dpop

      dDON_dt = - Gdoc2zoo/rCNzoo - Ddon + Dpon2don
      dPON_dt = (Mphy + Gphy2zoo*(1.0d0-e_Gphy2zoo))/rCNphy            &
     &         +(Mzoo - Gpon2zoo)/rCNzoo                               &
     &         - Dpon -Dpon2don

      dDOP_dt = - Gdop2zoo/rCPzoo - Ddop + Dpop2dop
      dPOP_dt = (Mphy + Gphy2zoo*(1.0d0-e_Gphy2zoo))/rCPphy            &
     &         +(Mzoo - Gpop2zoo)/rCPzoo                               &
     &         - Dpop -Dpop2dop
#endif

      dDOx_dt  = Pphy - Rphy - Rzoo - Ddoc - Dpoc

      RETURN

    END SUBROUTINE foodweb

  END MODULE mod_foodweb


