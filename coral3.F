c
!*************** CORALS ************************************************
c
!#define TESTMODE
c
      subroutine coral_setini
C **********************************************************************
C *                                                                    *
C * FUNCTION    : Calculate Flux (umol m-2 s-1).                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k_v9.blk'
c
      double precision  R13C_fromd13C
      double precision  R13C

!----------set data -----------------------

c   Coral minimum quota
      data QC0coral /286.  /   ! Minimum quota (umol cm-2)   (Tanaka et al. 2010)  ＊minimumかどうかは分からない。
      data QN0coral / 28.  /   ! Minimum quota (umol cm-2)  
      data QP0coral /  0.85/   ! Minimum quota (umol cm-2) 
c   Coral cell capacitycapasity of nutrients
c      data NO3coral_Max /4.e-3  /   ! (umol cm-2) !  70 umol L-1 Approxmatry estimated from Agostini et al. 2012
c      data NO2coral_Max /5.e-4  /   ! (umol cm-2) !   9 umol L-1 Approxmatry estimated from Agostini et al. 2012
c      data NH4coral_Max /5.5e-2  /  ! (umol cm-2) ! 540 umol L-1 Approxmatry estimated from Agostini et al. 2012
c
c      data NO3coral_Max / 70./   ! (umol kg-1) !  70 umol L-1 Approxmatry estimated from Agostini et al. 2012
c      data NO2coral_Max /  9./   ! (umol kg-1) !   9 umol L-1 Approxmatry estimated from Agostini et al. 2012
c      data NH4coral_Max /540./   ! (umol kg-1) ! 540 umol L-1 Approxmatry estimated from Agostini et al. 2012
c      data PO4coral_Max /920./   ! (umol kg-1) ! 920 umol L-1 Approxmatry estimated from Agostini et al. 2012
c
c   Zooxanthellae minimum cell quota      
      data QC0zoox   /8. / !(pmolC/cell)  100 pg cell-1 ~ 10 pmol cell-1 = (Verde & McCloskey 1998)
      data QN0zoox   /1.1/ !(pmolN/cell)  calculated from C:N:P = 116:16:1
      data QP0zoox   /0.07/ !(pmolP/cell) 
      
      data V_zoox    / 5.e-10/ ! Volume of zooxanthellae (cm3/cells)
                               ! ~500 um3 =5.e-10 cm3 (e.g. Jone 1997; Jones and Yellowlees 1997)
c   Zooxanthellae cell capacitycapasity of nutrients
c      data NO3zoox_Max /0.02  /   ! (pmolN/cell) !! Flynn and Flynn, 1998
c      data NO2zoox_Max /0.005  /  ! (pmolN/cell) !! Flynn and Flynn, 1998
c      data NH4zoox_Max /0.01  /   ! (pmolN/cell) !! Flynn and Flynn, 1998
c      data PO4zoox_Max /0.01  /   ! (pmolP/cell) !! Flynn and Flynn, 1998
c
!------------------------------------------

      integer i,j,k,n,m
c

c  Set initial conditions
      do j=1,jm
        do i=1,im
          do n=1,ccom

c        coral internal conditions

            TAcal (n,i,j)=2550.
            TAcoe (n,i,j)=2100.
            DICcal (n,i,j)=2100.
            DICcoe (n,i,j)=2000.
            DOcoe (n,i,j)=10. !(umol/kg)
            CH2O (n,i,j)=3.e5 !450 ug/cm2 = 15 umol cm-2 -> 5.e6 umol kg-1 (腔腸内平均濃度とした）(Achituv et al. 1994)
            QCcoral(n,i,j)=QC0coral(n)*1.1  !(umol cm-2)
            QNcoral(n,i,j)=QN0coral(n)*1.1  !(umol cm-2)
            QPcoral(n,i,j)=QP0coral(n)*1.1  !(umol cm-2)
            ROS(n,i,j)=0.  !(umol kg-1)
c            NO3coral(n,i,j)=0.2*NO3coral_Max(n)
c            NO2coral(n,i,j)=0.2*NO2coral_Max(n)
c            NH4coral(n,i,j)=0.2*NH4coral_Max(n)
c            PO4coral(n,i,j)=0.2*PO4coral_Max(n)
            NO3coral(n,i,j)= 25.   !(umol kg-1) Agostini et al. 2012
            NO2coral(n,i,j)=  5.2  !(umol kg-1) Agostini et al. 2012
            NH4coral(n,i,j)=275.   !(umol kg-1) Agostini et al. 2012
            PO4coral(n,i,j)=275.   !(umol kg-1) Agostini et al. 2012
c
            R13C=R13C_fromd13C(-15.d0)
            c13CH2O (n,i,j)=R13C/(1.+R13C)*CH2O(n,i,j)
c
            R13C=R13C_fromd13C(0.d0)
            DI13Ccoe (n,i,j)=R13C/(1.+R13C)*DICcoe (n,i,j)
c
            R13C=R13C_fromd13C(0.d0)
            DI13Ccal (n,i,j)=R13C/(1.+R13C)*DICcal (n,i,j)

c        zooxanthellae conditions

            QCzoox(n,i,j)=QC0zoox*1.1    !(pmol cell-1)
            QNzoox(n,i,j)=QN0zoox*1.1    !(pmol cell-1)
            QPzoox(n,i,j)=QP0zoox*1.1    !(pmol cell-1)
            DENSzoox(n,i,j) = 1.e6        !(cells cm-2)
            
            dsed_coral(n,i,j)=0.        !(cm)
            
c            NO3zoox(n,i,j)=0.2*NO3zoox_Max
c            NO2zoox(n,i,j)=0.2*NO2zoox_Max
c            NH4zoox(n,i,j)=0.2*NH4zoox_Max
c            PO4zoox(n,i,j)=0.2*PO4zoox_Max
            NO3zoox(n,i,j)=0.01     ! (pmolN/cell) !! Flynn and Flynn, 1998
            NO2zoox(n,i,j)=0.002    ! (pmolN/cell) !! Flynn and Flynn, 1998
            NH4zoox(n,i,j)=0.1     ! (pmolN/cell) !! Flynn and Flynn, 1998
            PO4zoox(n,i,j)=0.01     ! (pmolP/cell) !! Flynn and Flynn, 1998
c
          enddo
c
c          do n=1,scom
c            s_QCb(n,i,j)=p_sgrass(n,i,j)!!coverage to biomass (molC/m2)
c            s_QNb(n,i,j)=s_QCb(n,i,j)*s_CNP(nN)/s_CNP(nC)
c            s_QPb(n,i,j)=s_QCb(n,i,j)*s_CNP(nP)/s_CNP(nC)
cc
c            s_QC(n,i,j)=s_QCb(n,i,j)
c            s_QN(n,i,j)=s_QNb(n,i,j)
c            s_QP(n,i,j)=s_QPb(n,i,j)
c          enddo
c
c          do n=1,acom
c            a_QCb(n,i,j)=p_algae(n,i,j)!!coverage to biomass (molC/m2)
c            a_QNb(n,i,j)=a_QCb(n,i,j)*a_CNP(nN)/a_CNP(nC)
c            a_QPb(n,i,j)=a_QCb(n,i,j)*a_CNP(nP)/a_CNP(nC)
cc
c            a_QC(n,i,j)=a_QCb(n,i,j)
c            a_QN(n,i,j)=a_QNb(n,i,j)
c            a_QP(n,i,j)=a_QPb(n,i,j)
c          enddo
        enddo
      enddo
c
! --------------- Isotope ----------------------
c
      a_phot=exp(-20./1000.)!!!!!!!!!!!!!!!!!!適当
      a_resp=exp(0./1000.)
      a_calc=exp(-7./1000.)
      a_co2 =exp(-8./1000.)!!!!!!!!!!!!!!!!!!適当-8-4.4
c
      S_PFD_dt=0.
      S_Gn_dt=0.
      S_d13CargxGn_dt=0.
      S_d13CCH2O_dt=0.
      S_Pg_dt=0.
      S_R_dt=0.
      S_CH2O_dt=0.
      iday=1
c
      return
      end
      
!______________________________________________________________________
c      
c      subroutine coral_simple
cc          input parameters
c     &            (n,i,j
c     &             ,PFD
c     &             )
cC **********************************************************************
cC *                                                                    *
cC * FUNCTION    : Coral Metabolism.                                    *
cC *                                                                    *
cC **********************************************************************
cC
c      implicit none
cC
c      include 'pom2k_v9.blk'
cc input parameters
c      integer n,i,j
c      double precision  PFD
c      dimension nc(nb),pc(nb)
cc
cc --- C:N:P ratio of benthic organisms(1:coral, 2:seagrass) -
cc *** Animal ***
cc      nc(1)=1./6.37 !G.Muller-Parker et al.(1994)
cc      pc(1)=1./172.
cc *** Zooxanthellae ***
c      nc(1)=1./19.66 !G.Muller-Parker et al.(1994)
c      pc(1)=1./365.
c
cc --- Photosynthesis and Calcification Parameters (coral) ---
cc      pmax(1) = 85.9 !referenced from Hata et al(2002)
cc      pIk(1) = 526.
cc      p0(1) = 0. !-24.8
cc      cmax(1) = 14.5
cc      cIk(1) = 35.
cc      c0(1) = -1.7
c      pmax(1) = 58.1 !1999 metabolism referenced from Kayanne et al(2005)
c      pIk(1)  = 213.
c      p0(1)   = -16.5
c      cmax(1) = 10.4
c      cIk(1)  = 30.
c      c0(1)   = -0.6
c      
c      porg(n,i,j)=(pmax(n)*tanh(ephq/pIk(n))+p0(n))/3600
c     &            *p_coral(n,i,j)!*aruft(i,j)  !!TN:remove       !Light response curve [mmolC/m2/s]
c      ginorg(n,i,j)=(cmax(n)*tanh(ephq/cIk(n))+c0(n))/3600
c     &            *p_coral(n,i,j)                                !Light response curve [mmolC/m2/s]
c
c      if(porg(n,i,j).gt.0.) then
c        fno3(n,i,j)=porg(n,i,j)*nc(n)*npref(n,i,j)         !Flux caliculated in Carbon assimilated by Photosynthesis times C:N:P ratio [mmol/m2/s]
c        fnh4(n,i,j)=porg(n,i,j)*nc(n)*(1.-npref(n,i,j))
c        fpo4(n,i,j)=porg(n,i,j)*pc(n)
c      else
c        fno3(n,i,j)=0.d0                   !!TN:change
c        fnh4(n,i,j)=porg(n,i,j)*nc(n)    !!TN:change
c        fpo4(n,i,j)=porg(n,i,j)*pc(n)    !!TN:change
c      endif
c
c      return
c      end
c
!______________________________________________________________________

      subroutine coral_polyp
c          input parameters
     &            (n,i,j
     &            ,PFD
     &            ,Tamb
     &            ,Samb
     &            ,DICamb
     &            ,TAamb
     &            ,DOamb
     &            ,NO3amb
     &            ,NO2amb
     &            ,NH4amb
     &            ,PO4amb
     &            ,TOCamb
     &            ,TONamb
     &            ,TOPamb
     &            ,PP1amb
     &            ,ZP1amb
     &            ,DI13Camb
     &            ,tau_amb
     &            ,Fsed
c          output parameters
     &            ,DICuptake
     &            ,TAuptake
     &            ,DOuptake
     &            ,NO3uptake
     &            ,NO2uptake
     &            ,NH4uptake
     &            ,PO4uptake
     &            ,TOCuptake
     &            ,TONuptake
     &            ,TOPuptake
     &            ,PP1uptake
     &            ,ZP1uptake
     &            ,DI13Cuptake
     &             )
c
C **********************************************************************
C *                                                                    *
C * FUNCTION    : Coral Metabolism.                                    *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k_v9.blk'
c input parameters
      integer n,i,j
      double precision
     &             PFD
     &            ,Tamb
     &            ,Samb
     &            ,DICamb
     &            ,TAamb
     &            ,DOamb
     &            ,NO3amb
     &            ,NO2amb
     &            ,NH4amb
     &            ,PO4amb
     &            ,TOCamb
     &            ,TONamb
     &            ,TOPamb
     &            ,PP1amb
     &            ,ZP1amb
     &            ,DI13Camb
     &            ,tau_amb
     &            ,Fsed
c          output parameters
     &            ,DICuptake
     &            ,TAuptake
     &            ,DOuptake
     &            ,NO3uptake
     &            ,NO2uptake
     &            ,NH4uptake
     &            ,PO4uptake
     &            ,TOCuptake
     &            ,TONuptake
     &            ,TOPuptake
     &            ,PP1uptake
     &            ,ZP1uptake
     &            ,DI13Cuptake

      double precision
     &             Pg
     &            ,R
     &            ,Pn
     &            ,Gn
c
      double precision  pH_fromATCT,cHCO3_fromATpH,cCO3_fromATpH
     &    ,cCO2aq_fromATpH,fCO2_fromATpH,Warg_fromcCO3
     &    ,fCO2_fromcCO2aq
c
     &    ,d13C_fromR13C,R13C_fromd13C

ccc----- Physical parameters -------------------
      double precision  TKamb, Rgas   !Temperature (K)
      double precision  Si(ncom)

c----- for zooxanthellae interaction -------------------
      double precision
     &             Expul                    ! Zooxanthellae expulsion rate by host coral (cell cm-2 s-1)
     &            ,C_excre,N_excre,P_excre  ! Organic matter excretion controled by host coral (nmol cm-2 s-1)
     &            ,C_waste,N_waste,P_waste  ! Waste organic matter (nmol cm-2 s-1)
     &            ,F_ROS                    ! Flux of Reactive Oxygen Species  (nmol cm-2 s-1)
     &            ,F_ONOO                   ! Flux of peroxynitrite  (nmol cm-2 s-1)
     &            ,PFDsurf
     
     &            ,k_excre(ccom)
     &            ,Rel_max(ccom)
     &            ,Fk_damage(ccom)
     &            ,lamb
     &            ,NO3_trans, NO2_trans, NH4_trans ! Nutrients transport rate from host to Zooxanthellae (nmol cm-2 s-1)
     &            ,PO4_trans

c----- for respiration --------------------------
      double precision  Rc           !Coral respiration rate (nmol cm-2 s-1)
      double precision  Rz           !Zooxanthellae respiration rate (nmol cm-2 s-1)
      double precision  Rmax(ccom)
      double precision  K_CH2O(ccom)
      double precision  K_QC(ccom)
      double precision  K_DO(ccom)
c       Temperature related factor
      double precision  Eact(ccom)
      double precision  f_temp
      double precision  SQC
c----- for calcification ------------------------
      double precision  E_tot,E_ca,E_other,F_H,F_CO2
      double precision  pHcal,fCO2cal,cCO3cal,cCO2aqcal,Wargcal
      double precision  pHcoe,fCO2coe,cCO2aqcoe,cHCO3coe
      double precision  pHamb,fCO2amb,cCO3amb,Wargamb  !!!!!不要：アウトプット用
      double precision  k_Gn(ccom)
      double precision  eeff,cal_tick,coe_tick
      double precision  k_CO2i
      double precision  leak_TA,leak_DIC
      double precision  k_TA, k_DIC
c----- for coral metabolism -------------------
      double precision  C_mucus, N_mucus, P_mucus
      double precision  C_growth,N_growth,P_growth
      double precision  g_max(ccom)

      double precision  C_ingest,N_ingest,P_ingest
      double precision  assim
      
      double precision  RMsed    ! sediment removing rate (g cm-2 s-1)

c----- Nutrients related parameters -------------------
      double precision  
     &   N_dissim       ! Nitrogen dissimilation rate (pmol cell-1 s-1)
     &  ,Vmax_NO3
     &  ,Vmax_NO2
     &  ,Vmax_NH4
     &  ,Vmax_PO4
     &  ,k_NH4rel
     &  ,k_NO3rel
     &  ,k_NO2rel
     &  ,k_PO4rel

c----- ROS related parameters -------------------
      double precision  V_detox(ccom)        ! Maximum detox rate (nmol cm-2 s-1)
      double precision  K_ROS(ccom)          ! half saturation constant (umol kg-1)
      double precision  k_damage(ccom)       ! Damage rate
      double precision  k_detox(ccom)        ! Detox constant
      double precision  F_detox              ! Detox rate  (nmol cm-2 s-1)
      double precision  F_damage             ! Damage rate (nmol cm-2 s-1)
      double precision  Damage               ! Damage function

c----- for cabon isotope calculation -------------------
      double precision  R13Ccoe,R13Ccal,R13CH2O,R13Camb,R13Carg
      double precision  d13C_DICcoe,d13C_DICcal,d13C_CH2O,d13C_DICamb
     &                 ,d13C_arg
      double precision  R13Cleak,R13Cfco2,R13Cuptake
      double precision  F_13CO2,k_13CO2i,c13CO2aqcoe,c13CO2aqcal
      double precision  leak_DI13C

c----- Size dynamics parameters -------------------
      double precision  km_CNP, am_CNP, km_dam, km_min
c
!------------Set Constants  ----------------------------------
c
c----- Physical constants ------------------------
c
      data Rgas  /8.314/ ! Gas constant (J mol-1)
c
c
!------------Coral parameters ---------------------------------
c
c----- Coral basic parameters ----------------------------------------
c
      data coe_tick /0.3/ !Tickness of coral tissue (cm)
      data cal_tick /0.0005/ !Tickness of coral calcifying fluid (cm)
      
c----- for zooxanthellae interaction -------------------

      data lamb /1.e3/   !         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*適当
      
c      data k_excre  /1.e-14/    ! (cm s-1)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*適当
      data k_excre   /1.e-5/      ! k_excre (umol cm-2 s-1 cell-1)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*適当
      data Rel_max   /1.e1/  !         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*適当
      data Fk_damage /1.e-3/   !         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*適当

c----- Nutrients related parameters -------------------

      data Vmax_NH4  /1.e-2/  !(nmol cm-2 s-1) Godinot et al. (2011)
      data Vmax_NO3  /8.e-3/  !(nmol cm-2 s-1) Godinot et al. (2011) 
      data Vmax_NO2  /8.e-3/  !(nmol cm-2 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data Vmax_PO4  /5.e-3/  !(nmol cm-2 s-1) Godinot et al. (2011)
      data k_NH4rel  /5.e-5/  !(cm s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data k_NO3rel  /2.42e-5/  !(cm s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data k_NO2rel  /3.e-6/  !(cm s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data k_PO4rel  /1.e-6/  !(cm s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      
c----- for respiration --------------------------
c
c      R= k_R(n)*CH2O(n,i,j)/(K_CH2O(n)+CH2O(n,i,j))
c     &            *DOcoe(n,i,j)/(K_DO(n)+DOcoe(n,i,j))
c
      data Rmax    /0.37/!/0.37/ !Maximum respiration rate (nmol O2 cm-2 s-1) *Tuned using Kuhl et al (1995), Al-Hprani et al. (2003) data
      data K_CH2O  /1.e5/ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*適当
      data K_QC    /280./ !/10./ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*適当
      data K_DO    /150./ ! ca. 10 (umol kg-1) Shick et al. (1990) *Tuned value
c
c     Temperature dependency parameters (Hikosaka et al. 2006) 
c        - Parameters estimated from data by Fujimura et al. (2008)
c
      data Eact   /5.117e4/ !Activation energy (J mol-1)) *Tuned
c
c----- for calcification ------------------------
c
      data k_Gn /1.e-2/

      data eeff /0.5/ !Energy efficency of coral
      data E_other /1.e3/ ! Energy flux of other metabolisms (nJ cm-2 s-1)
c
      data k_CO2i /1.5e-3/  !permeability coefficient (cm s-1): Gutknecht et al. 1977: (0.35 cm s-1)
      data k_TA   /3.0e-3/ !conductivity of TA through the leak pass (cm s-1)
      data k_DIC  /3.0e-3/ !conductivity of DIC through the leak pass (cm s-1)
c
c---- for mucus relrase --------------------------



c----- for coral metabolism -------------------
      data g_max    /3.e-2/       ! Maximum growth rate (nmol cm-2 s-1)  !!!!!!!!!!!!!!!!!!適当!!!!!!!!!!!!!!
      data assim    /0.5 /       ! assimilation efficiency              !!!!!!!!!!!!!!!!!!適当!!!!!!!!!!!!!!
      
c----- for carbon isotope  ------------------------

c----- ROS related parameters -------------------

      data  K_ROS    /1./          ! half saturation constant (umol kg-1)!  適当!!!!!!!!!!!!!!
      data  k_detox  /6.e-4/           ! Detox constant                    !  適当!!!!!!!!!!!!!!
      data  k_damage /1.e-4/           ! Damage rate                       !  適当!!!!!!!!!!!!!!
 
c----- Size dynamics parameters -------------------
      
      data km_CNP    /  1.e-10/ !  適当!!!!!!!!!!!!!!
      data am_CNP    /  1.e1 /
      data km_dam    /  1.e-6/ !  適当!!!!!!!!!!!!!!
      data km_min    /  1.e-10/ !  適当!!!!!!!!!!!!!!

!-----Set Physical parameters --------------------------------------------
c
c----------- Temperature (K)  --------------------------------
c
      TKamb=Tamb+273.15 !(oC->K)
c
c----------- calculate mass transfer velocity (cm s-1) -------
c
      Si(nDIC)=(65.7*tau_amb**0.4+4.7)*1.e-4
      Si(nTA )=(65.7*tau_amb**0.4+4.7)*1.e-4
      Si(nLDOC)=(65.7*tau_amb**0.4+4.7)*1.e-4
c      Si(nNO3)=(65.7*tau_amb**0.4+4.7)*1.e-4
c      Si(nNO2)=(65.7*tau_amb**0.4+4.7)*1.e-4
c      Si(nNH4)=(65.7*tau_amb**0.4+4.7)*1.e-4
c      Si(nPO4)=(65.7*tau_amb**0.4+4.7)*1.e-4
      Si(nDO )=(65.7*tau_amb**0.4+4.7)*1.e-4
      Si(nPP1)=(65.7*tau_amb**0.4+4.7)*1.e-4
      Si(nZP1)=(65.7*tau_amb**0.4+4.7)*1.e-4

      Si(nNO3)=(65.7*tau_amb**0.4 +  27.)*1.e-4   ! 27.e-4 cm s-1 @ zero flow was estimated from Vmax/K=Si (Data from Godinot et al. 2011) 
      Si(nNO2)=(65.7*tau_amb**0.4 +  27.)*1.e-4   ! 27: てきとう
      Si(nNH4)=(65.7*tau_amb**0.4 + 213.)*1.e-4   ! (Data from Godinot et al. 2011) 
      Si(nPO4)=(65.7*tau_amb**0.4 +  52.)*1.e-4   ! (Data from Godinot et al. 2011) 

c----------- CO2 system in calcified fulid -------------------
c
      pHcal=pH_fromATCT(TAcal(n,i,j),DICcal(n,i,j),TKamb, Samb)
      cCO2aqcal=cCO2aq_fromATpH(TAcal(n,i,j),pHcal,TKamb, Samb)
      fCO2cal=fCO2_fromATpH(TAcal(n,i,j),pHcal,TKamb, Samb)
      cCO3cal=cCO3_fromATpH(TAcal(n,i,j),pHcal,TKamb, Samb)
      Wargcal=Warg_fromcCO3(cCO3cal,TKamb, Samb)
c
c----------- CO2 system in coelenteron -------------------
c
      pHcoe=pH_fromATCT(TAcoe(n,i,j),DICcoe(n,i,j),TKamb, Samb)
      cCO2aqcoe=cCO2aq_fromATpH(TAcoe(n,i,j),pHcoe,TKamb, Samb)
      fCO2coe=fCO2_fromcCO2aq(cCO2aqcoe,TKamb, Samb)
      cHCO3coe=cHCO3_fromATpH(TAcoe(n,i,j),pHcoe,TKamb, Samb)
c
c----------- CO2 system in ambient seawater ------------------- !!不要：アウトプット用
c
#if defined TESTMODE
      pHamb=pH_fromATCT(TAamb,DICamb,TKamb, Samb)
      fCO2amb=fCO2_fromATpH(TAamb,pHamb,TKamb, Samb)
      cCO3amb=cCO3_fromATpH(TAamb,pHamb,TKamb, Samb)
      Wargamb=Warg_fromcCO3(cCO3amb,TKamb, Samb)
#endif
c
c---------- ROS related fluxes (nmol cm-2 s-1) ----------------------------

!     V_detoxは余剰の栄養素の関数
c
c      V_detox=k_detox(n)*exp(
c     &                     min(  1.-QP0coral(n)/QPcoral(n,i,j)  ,
c     &                       min(1.-QN0coral(n)/QNcoral(n,i,j),
c     &                           1.-QC0coral(n)/QCcoral(n,i,j) ) )
c     &                       )
c
c
c      F_detox=V_detox(n)*ROS(n,i,j)/(K_ROS(n)+ROS(n,i,j))
      F_detox=0.15*ROS(n,i,j)*1.024*coe_tick  !(s-1)*(umol kg-1)*1024 (kg cm-3)*(cm)*0.001(n/u)=(nmol cm-2 s-1)
             !0.15
c      F_detox=0.02*ROS(n,i,j)/(K_ROS(n)+ROS(n,i,j))

c      F_detox=0.035*ROS(n,i,j)/(0.1+ROS(n,i,j))
      
c      F_damage=k_damage(n)*ROS(n,i,j)
      F_damage=0.016*ROS(n,i,j)*1.024*coe_tick


c----- peroxynitrite (ONOO-) flux (nmol cm-2 s-1)------------------------------------
c
c      F_ONOO=1.e0 * NO2zoox(n,i,j)*ROS(n,i,j)   ! (pmol cell-1 s-1)!!!!!!!!!!!!!!!!!!てきとう
c      
c      F_ONOO=F_ONOO*DENSzoox(n,i,j)*1.e-3   ! (nmol cm-2 s-1)

      F_ONOO=1.e-6* NO2zoox(n,i,j)*ROS(n,i,j)*DENSzoox(n,i,j)   ! (pmol cell-1 s-1)!!!!!!!!!!!!!!!!!!てきとう



!------ Zooxanthellae interction section ----------------------------

c------- Nutrients flux (nmol cm-2 s-1) ------------------------
c
c NH4 uptake
      NH4uptake=Vmax_NH4 * NH4amb /(Vmax_NH4/Si(nNH4) + NH4amb)
c     &          -2.* Vmax_NH4
c     &             *(NH4coral(n,j,i)/NH4coral_Max(n))**4.
c     &               /((NH4coral(n,j,i)/NH4coral_Max(n))**4. + 1.)

     &          -k_NH4rel * NH4coral(n,j,i)

c NO3 uptake
      NO3uptake=Vmax_NO3 * NO3amb /(Vmax_NO3/Si(nNO3) + NO3amb)
c     &          -2.* Vmax_NO3
c     &             *(NO3coral(n,j,i)/NO3coral_Max(n))**4.
c     &               /((NO3coral(n,j,i)/NO3coral_Max(n))**4. + 1.)

     &          -k_NO3rel * NO3coral(n,j,i)

c NO2 uptake
      NO2uptake=Vmax_NO2 * NO2amb /(Vmax_NO2/Si(nNO2) + NO2amb)
c     &          -2.* Vmax_NO2
c     &             *(NO2coral(n,j,i)/NO2coral_Max(n))**4.
c     &               /((NO2coral(n,j,i)/NO2coral_Max(n))**4. + 1.)

     &          -k_NO2rel * NO2coral(n,j,i)


c PO4 uptake
      PO4uptake=Vmax_PO4 * PO4amb /(Vmax_PO4/Si(nPO4) + PO4amb)
c     &          -2.* Vmax_PO4
c     &             *(PO4coral(n,j,i)/PO4coral_Max(n))**4.
c     &               /((PO4coral(n,j,i)/PO4coral_Max(n))**4. + 1.)

     &          -k_PO4rel * PO4coral(n,j,i)

c------ Organic matter excretion controled by host coral (nmol cm-2 s-1) -----

c      C_excre=k_excre(n)*(                                               ! k_excre (cm s-1)
c     &                  (QCzoox(n,i,j)-QC0zoox)/V_zoox *1.e-3         ! (pmol cell-1)/(cm3 cell-1)*1.e-3= (nmol cm-3)
c     &                 -(QCcoral(n,i,j)-QC0coral(n))/coe_tick *1.e3   ! (umol cm-2)/(cm)*1.e3 = (nmol cm-3)
c     &                 )*DENSzoox(n,i,j)
c      
c      N_excre=k_excre(n)*(
c     &                  (QNzoox(n,i,j)-QN0zoox)/V_zoox *1.e-3
c     &                 -(QNcoral(n,i,j)-QN0coral(n))/coe_tick *1.e3
c     &                 )*DENSzoox(n,i,j)
c      
c      P_excre=k_excre(n)*(
c     &                  (QPzoox(n,i,j)-QP0zoox)/V_zoox *1.e-3
c     &                 -(QPcoral(n,i,j)-QP0coral(n))/coe_tick *1.e3
c     &                 )*DENSzoox(n,i,j)
c     
      
      C_excre=k_excre(n)*(                            ! k_excre (umol cm-2 s-1 cell-1)
     &                  QCzoox(n,i,j)/QC0zoox         ! 
     &                 -QCcoral(n,i,j)/QC0coral(n)   !
     &                 )*DENSzoox(n,i,j)
      
      N_excre=k_excre(n)*(
     &                  QNzoox(n,i,j)/QN0zoox
     &                 -QNcoral(n,i,j)/QN0coral(n)
     &                 )*DENSzoox(n,i,j)
      
      P_excre=k_excre(n)*(
     &                  QPzoox(n,i,j)/QP0zoox
     &                 -QPcoral(n,i,j)/QP0coral(n)
     &                 )*DENSzoox(n,i,j)
           
c      N_excre=0.
c      P_excre=0.


c------ Zooxanthellae expulsion rate by host coral (cell cm-2 s-1) -----

c      Expul=0.!Rel_max(n)*tanh(F_damage/Fk_damage(n))
c      Expul=200.*F_damage
c     &            *DENSzoox(n,i,j)/(5.e-5+DENSzoox(n,i,j))

     
c      Damage=exp(F_damage + 3.*F_ONOO)-1.   !Damage by ROS (F_damage) and peroxynitrite (F_ONOO)

      Damage=exp(F_damage + F_ONOO)-1.   !Damage by ROS (F_damage) and peroxynitrite (F_ONOO)
     
      Expul=0.6e2*Damage
c     &            *DENSzoox(n,i,j)/(5.e-5+DENSzoox(n,i,j))     
      
c------ Zooxanthellae population dynamics model ----------------------------- 

c   Sediment dependent term
      PFDsurf=PFD*exp(-lamb*dsed_coral(n,i,j))  !!!!!!!!!!!!!!!!!!!!!赤土の影響　要検討!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call zooxanthellae
c          input parameters
     &            (n,i,j
     &            ,PFDsurf
     &            ,TKamb
     &            ,cCO2aqcoe, cHCO3coe
c     &            ,PO4uptake                       ! PO4 uptake (nmol cm-2 s-1)
     &            ,Expul                           ! Zooxanthellae expulsion rate by host coral (cell cm-2 s-1)
     &            ,C_excre,N_excre,P_excre         ! Organic matter excretion controled by host coral (nmol cm-2 s-1)
c          output parameters
     &            ,Pg                              ! Gross photosynthesis rate (nmol cm-2 s-1)
     &            ,Rz                              ! Zooxanthellae respiration rate (nmol cm-2 s-1)
     &            ,NO3_trans, NO2_trans, NH4_trans ! Nutrients transport rate from host to Zooxanthellae (nmol cm-2 s-1)
     &            ,PO4_trans
     &            ,C_waste,N_waste,P_waste         ! Waste organic matter (nmol cm-2 s-1)
     &            ,F_ROS                           ! Flux of Reactive Oxygen Species  (nmol cm-2 s-1)
     &            ,F_ONOO                          ! Flux of peroxynitrite  (nmol cm-2 s-1)
     &             )


!---- Coral metabolism section ------------------------------------------
c
c----- Respilation rate (nmolO2 cm-2 s-1) ---------------------------

c   Temperature dependent term
      f_temp=1.
!      f_temp=exp(Eact(n)*(TKamb-298.)/298./Rgas/TKamb)
c
cc   Respilation rate (nmolO2 cm-2 s-1) 
c      R= Rmax(n)*f_temp*CH2O(n,i,j)/(K_CH2O(n)+CH2O(n,i,j))
c     &            *DOcoe(n,i,j)/(K_DO(n)+DOcoe(n,i,j))
c
      Rc= Rmax(n)*f_temp*QCcoral(n,i,j)/(K_QC(n)+QCcoral(n,i,j))
     &            *DOcoe(n,i,j)/(K_DO(n)+DOcoe(n,i,j))

c      SQC=QCcoral(n,i,j)-QC0coral(n)
c
c      Rc= Rmax(n)*f_temp*SQC/(K_QC(n)+SQC)
c     &            *DOcoe(n,i,j)/(K_DO(n)+DOcoe(n,i,j))
      R = Rc+Rz

c ---Net photosynthesis rate (not used for calculation, output only)
c
      Pn=(Pg-R)!(nmol cm-2 s-1)
c

c----- Nitrogen dissimilation rate (nmol cm-2 s-1) ---------------------

      N_dissim = Rc * QNcoral(n,i,j)/QCcoral(n,i,j)   ! Dissimilation rate (pmol cell-1 s-1)

c----- Calcification rate (nmol cm-2 s-1) ------------------------------------

      Gn=k_Gn(n)*(Wargcal-1.)


c----- H+ flux through Ca2+ ATPase  (nmol cm-2 s-1) --------------------------

c   energy Flux (nJ cm-2 s-1) partition

c  C6H12O6 + 6 O2 + 38 ADP + 38 Pi -> 6 CO2 + 6 H2O + 38 ATP  dG=-2873.4 kJ/mol
c  ATP + H2O -> ADP + Pi  dG=-7.3kcal/mol=-30.5kJ/mol

      E_tot=30.5e3*38./6.*R  !nJ cm-2 s-1
c      E_other=1000.           !constant -> set at the parameter section
      E_ca=E_tot-E_other
      E_ca=max(E_ca,0.) !Error handring

c   H+ flux (nmol cm-2 s-1)

      F_H= E_ca * eeff/
     &  (2.3*8.31*(Tamb+273.15)*max((pHcal-pHcoe)+0.01,0.01d0))


c---------- Coral mucus (nmol cm-2 s-1) ------------
c
c      m_sed=k_m_sed(n)*dsed
c
      if(dsed_coral(n,i,j).gt.0.) then
        C_mucus=1.e6 * dsed_coral(n,i,j)
     &               * max(0.,
     &                     min(1.-QN0coral(n)/QNcoral(n,i,j),
     &                         1.-QC0coral(n)/QCcoral(n,i,j) )  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!適当
     &                    )
      else
        C_mucus=0.
      endif
      
      Rmsed=1.e-7 * C_mucus  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!適当
      
c      C_mucus=25.e-4       !
      N_mucus= C_mucus*1./25.
      P_mucus= 0.     !!!!!!!!!!!!!!!!適当

c      C_mucus=25.e-4       ! 9   (nmol cm-2 h-1) ~ 25.e-4 (nmol cm-2 s-1) (Tanaka et al. 2010)
c      N_mucus= 1.e-4       ! 0.4 (nmol cm-2 h-1) ~ 0.1e-3 (nmol cm-2 s-1) (Tanaka et al. 2010)
c      P_mucus= 0.06e-4     !!!!!!!!!!!!!!!!適当
c
c      C_mucus=0.       ! 9   (nmol cm-2 h-1) ~ 25.e-4 (nmol cm-2 s-1) (Tanaka et al. 2010)
c      N_mucus=0.       ! 0.4 (nmol cm-2 h-1) ~ 0.1e-3 (nmol cm-2 s-1) (Tanaka et al. 2010)
c      P_mucus=0.     !!!!!!!!!!!!!!!!適当
c
      C_mucus=C_mucus
     &        +1.e-6* max(0., 1.-QC0coral(n)/QCcoral(n,i,j))       ! 9   (nmol cm-2 h-1) ~ 25.e-4 (nmol cm-2 s-1) (Tanaka et al. 2010)
      N_mucus=N_mucus
     &        +1.e-6* max(0., 1.-QN0coral(n)/QNcoral(n,i,j))       ! 0.4 (nmol cm-2 h-1) ~ 0.1e-3 (nmol cm-2 s-1) (Tanaka et al. 2010)
      P_mucus=P_mucus
     &        +1.e-6* max(0., 1.-QP0coral(n)/QPcoral(n,i,j))     !!!!!!!!!!!!!!!!適当

c---------- Coral growth (nmol cm-2 s-1) ------------
c
c    Droop model (1973) + Liebig's minimum law

      C_growth=g_max(n)*min(    1.-QP0coral(n)/QPcoral(n,i,j)  ,
     &                      min(1.-QN0coral(n)/QNcoral(n,i,j),
     &                          1.-QC0coral(n)/QCcoral(n,i,j) ) )
      
      if(C_growth .lt. 0.) then
        C_growth =0.
      end if
      
      N_growth=C_growth*QNcoral(n,i,j)/QCcoral(n,i,j)
      
      P_growth=C_growth*QPcoral(n,i,j)/QCcoral(n,i,j)


!------ Output parameters: g_coral m_coral -------------------------------
c----- Convert internal conditions to size dynamics parameters

      g_coral(n,i,j)=C_growth/QCcoral(n,i,j)*1.e-3  ! growth rate (cm2 cm-2 s-1)
c     (cm2 cm-2 s-1)=(nmol cm-2 s-1)/(umol cm-2)*1.e-3

      m_coral(n,i,j)= km_CNP                          ! mortality (cm2 cm-2 s-1) 
     &                 *exp(-am_CNP*                                      !* 足りない栄養素に依存（要検討!!!）
     &                     min(  1.-QP0coral(n)/QPcoral(n,i,j)  ,
     &                       min(1.-QN0coral(n)/QNcoral(n,i,j),
     &                           1.-QC0coral(n)/QCcoral(n,i,j) ) )
     &                       )
     &               +km_dam * Damage                                !* ダメージ依存（要検討!!!）
c     &               +km_min                                    !* 定数

!--------------------------------------------------------------------------

c---------- Other fluxes (nmol cm-2 s-1) ----------------------------

c   CO2* flux (nmol cm-2 s-1)
c      F_CO2=k_CO2i*(fCO2coe-fCO2cal)
      F_CO2=k_CO2i*(cCO2aqcoe-cCO2aqcal)*1.024

c   Leak flux (nmol cm-2 s-1) 
      leak_TA =k_TA *(TAcoe(n,i,j)-TAcal(n,i,j))*1.024
      leak_DIC=k_DIC*(DICcoe(n,i,j)-DICcal(n,i,j))*1.024

c   DIC and TA uptake (nmol cm-2 s-1)
      DICuptake=Si(nDIC)*(DICamb-DICcoe(n,i,j))*1.024
      TAuptake =Si(nTA)*(TAamb-TAcoe(n,i,j))*1.024

c   DO uptake (nmol cm-2 s-1)
      DOuptake =Si(nDO)*(DOamb-DOcoe(n,i,j))*1.024


c----- Net TOC, TON, TOP uptake (nmol cm-2 s-1)  ------------------

c      TOCuptake=Si(nTOC)*TOCamb*1.024-mucus
c      TONuptake=Si(nTOC)*TONamb*1.024
c     &           -mucus*cMucus_CNP(nN)/cMucus_CNP(nC)
c      TOPuptake=Si(nTOC)*TOPamb*1.024
c     &           -mucus*cMucus_CNP(nP)/cMucus_CNP(nC)

      TOCuptake=-C_mucus-C_waste
     &          -m_coral(n,i,j)/QCcoral(n,i,j)*p_coral(n,i,j)*1.e3  !umol cm-2 s-1 -> nmol cm-2 s-1

      TONuptake=-N_mucus-N_waste
     &          -m_coral(n,i,j)/QNcoral(n,i,j)*p_coral(n,i,j)*1.e3  !umol cm-2 s-1 -> nmol cm-2 s-1

      TOPuptake=-P_mucus-P_waste
     &          -m_coral(n,i,j)/QPcoral(n,i,j)*p_coral(n,i,j)*1.e3  !umol cm-2 s-1 -> nmol cm-2 s-1


c---- Plankton predation (nmol cm-2 s-1)  ------------------
c      PP1uptake=Si(nPP1)*PP1amb*1.024
      PP1uptake=0.
      ZP1uptake=Si(nZP1)*ZP1amb*1.024

c   Temporal ingestion
!      C_ingest=ZP1uptake
!      N_ingest=ZP1uptake*zPlkt_CNP(nN)/zPlkt_CNP(nC)
!      P_ingest=ZP1uptake*zPlkt_CNP(nP)/zPlkt_CNP(nC)

c   Not assimilated materials      
      TOCuptake=TOCuptake - C_ingest*(1.-assim)
      TONuptake=TONuptake - N_ingest*(1.-assim)
      TOPuptake=TOPuptake - P_ingest*(1.-assim)

c   Assimilated materials
      C_ingest=C_ingest*assim
      N_ingest=N_ingest*assim
      P_ingest=P_ingest*assim



!------------ internal condition change -----------------------------

c---- Nutrients (umol cm-2) ------------------

      NO3coral(n,i,j)=NO3coral(n,i,j) 
     &              +(
     &                 NO3uptake     ! (nmol cm-2 s-1)
     &                -NO3_trans     ! (nmol cm-2 s-1)
c     &               )*1.e-3 *dtc     !umol cm-2
     &          )/coe_tick/1.024 *dtc !umol kg-1
      NO3coral(n,i,j)=max(NO3coral(n,i,j),0.d0) !Error handring

      NO2coral(n,i,j)=NO2coral(n,i,j)
     &              +(
     &                 NO2uptake     ! (nmol cm-2 s-1)
     &                -NO2_trans     ! (nmol cm-2 s-1)
c     &               )*1.e-3 *dtc     !umol cm-2
     &          )/coe_tick/1.024 *dtc !umol kg-1
      NO2coral(n,i,j)=max(NO2coral(n,i,j),0.d0) !Error handring

      NH4coral(n,i,j)=NH4coral(n,i,j)
     &              +(
     &                 NH4uptake     ! (nmol cm-2 s-1)
     &                -NH4_trans     ! (nmol cm-2 s-1)
     &                +N_dissim      ! (nmol cm-2 s-1)
c     &               )*1.e-3 *dtc     !umol cm-2
     &          )/coe_tick/1.024 *dtc !umol kg-1
      NH4coral(n,i,j)=max(NH4coral(n,i,j),0.d0) !Error handring


      PO4coral(n,i,j)=PO4coral(n,i,j)
     &              +(
     &                 PO4uptake     ! (nmol cm-2 s-1)
     &                -PO4_trans     ! (nmol cm-2 s-1)
c     &               )*1.e-3 *dtc     !umol cm-2
     &          )/coe_tick/1.024 *dtc !umol kg-1
      PO4coral(n,i,j)=max(PO4coral(n,i,j),0.d0) !Error handring


c     CH2O (umol kg-1)
c     Pg,R: nmol cm-2 s-1
      CH2O (n,i,j)=CH2O (n,i,j)
c     &         +(Pg-Rcarb-F_c2l)*1.e-3 *dtc
     &         +(Pg-R)/coe_tick/1.024 *dtc    !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1

c    Quota (umol cm-2)
      QCcoral(n,i,j)=QCcoral(n,i,j)
     &         +( C_excre      ! (nmol cm-2 s-1)
     &           +C_ingest     ! (nmol cm-2 s-1)
     &           -R            ! (nmol cm-2 s-1)
     &           -C_growth
     &           -C_mucus
     &          )*1.e-3 *dtc    !umol cm-2

      QNcoral(n,i,j)=QNcoral(n,i,j)
     &         +( N_excre
     &           +N_ingest
     &           -N_growth
     &           -N_mucus
     &          )*1.e-3 *dtc    !umol cm-2

      QPcoral(n,i,j)=QPcoral(n,i,j)
     &         +( P_excre
     &           +P_ingest
     &           -P_growth
     &           -P_mucus
     &          )*1.e-3 *dtc    !umol cm-2

c    ROS concentration (umol kg-1)
      ROS(n,i,j)=ROS(n,i,j)
     &         +( F_ROS
     &           -F_detox
     &           -F_damage
     &           -F_ONOO
     &           -Si(nDIC)*ROS(n,i,j)*1.024
     &          )/coe_tick/1.024 *dtc            !umol kg-1
      ROS(n,i,j)=max(ROS(n,i,j),0.d0) !Error handring
      
c    DIC & TA (umol kg-1) in calcified fulid
      DICcal (n,i,j)=DICcal (n,i,j)
     &            +(-Gn
     &              +F_CO2
     &              +leak_DIC
     &             )/cal_tick/1.024 *dtc    !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1
      DICcal(n,i,j)=max(DICcal(n,i,j),0.d0) !Error handring
      DICcal(n,i,j)=min(DICcal(n,i,j),100000.d0) !Error handring

      TAcal (n,i,j)=TAcal (n,i,j)
     &     +(-2.*Gn
     &       +F_H
     &       +leak_TA
     &      )/cal_tick/1.024 *dtc    !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1
      TAcal(n,i,j)=max(TAcal(n,i,j),0.d0) !Error handring
      TAcal(n,i,j)=min(TAcal(n,i,j),100000.d0) !Error handring


c    DIC & TA (umol kg-1) in coelectron
      DICcoe (n,i,j)=DICcoe (n,i,j)
     &            +(-Pg+R
     &              +DICuptake
     &              -F_CO2
     &              -leak_DIC
     &             )/coe_tick/1.024 *dtc   !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1
      DICcoe(n,i,j)=max(DICcoe(n,i,j),0.d0) !Error handring
      DICcoe(n,i,j)=min(DICcoe(n,i,j),100000.d0) !Error handring

      TAcoe (n,i,j)=TAcoe (n,i,j)
     &     +(-F_H
     &        +TAuptake
     &        -leak_TA
     &      )/coe_tick/1.024 *dtc   !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1
      TAcoe(n,i,j)=max(TAcoe(n,i,j),0.d0) !Error handring
      TAcoe(n,i,j)=min(TAcoe(n,i,j),100000.d0) !Error handring


c    DO (umol kg-1) in coelectron  !!TN: Add 20120718
      DOcoe (n,i,j)=DOcoe (n,i,j)
     &            +(Pg-R
     &              +DOuptake
     &             )/coe_tick/1.024 *dtc   !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1
     
c    Sediment depth on coral (cm)
      dsed_coral(n,i,j)=dsed_coral(n,i,j)  ! (cm)
     &            +(Fsed
     &             -Rmsed 
     &              )/1.4 *dtc      !(g cm-2 s-1)/1.4(g cm-3)*(s) = (cm)
      dsed_coral(n,i,j)=max(dsed_coral(n,i,j),0.d0) !Error handring                                          !       1.4 (g cm-3) bulk density of sediment

! ----- Isotope calculation ----------------------------------------------------
c
c Carbon isotope ratio R13C=13C/12C
c
      R13CH2O=c13CH2O(n,i,j)/(CH2O(n,i,j)-c13CH2O(n,i,j))   !coral organism
      R13Ccoe=DI13Ccoe(n,i,j)/(DICcoe(n,i,j)-DI13Ccoe(n,i,j))   !DIC in coelenteron
      R13Ccal=DI13Ccal(n,i,j)/(DICcal(n,i,j)-DI13Ccal(n,i,j))   !DIC in calcifying fluid
      R13Camb=DI13Camb/(DICamb-DI13Camb)!!!!!  DIC in ambient seawater
      
      c13CO2aqcoe=cCO2aqcoe*R13Ccoe*a_co2/(1.+R13Ccoe*a_co2)
      c13CO2aqcal=cCO2aqcal*R13Ccal*a_co2/(1.+R13Ccal*a_co2)
c
c Carbon isotope fluxes
c
      F_13CO2=k_CO2i/1.0044*(c13CO2aqcoe-c13CO2aqcal)*1.024
c      F_13CO2= 0.5*( dabs(F_CO2)*R13Ccoe*a_co2/(1.+R13Ccoe*a_co2)
c     &                   +F_CO2 *R13Ccoe*a_co2/(1.+R13Ccoe*a_co2)
c     &              -dabs(F_CO2)*R13Ccal*a_co2/(1.+R13Ccal*a_co2)
c     &                   +F_CO2 *R13Ccal*a_co2/(1.+R13Ccal*a_co2))
c
      DI13Cuptake=Si(nDIC)*(DI13Camb-DI13Ccoe(n,i,j))*1.024
c      DI13Cuptake=0.5*( dabs(DICuptake)*R13Camb/(1.+R13Camb)
c     &                      +DICuptake *R13Camb/(1.+R13Camb)
c     &                 -dabs(DICuptake)*R13Ccoe/(1.+R13Ccoe)
c     &                      +DICuptake *R13Ccoe/(1.+R13Ccoe))
c
      leak_DI13C=k_DIC*(DI13Ccoe(n,i,j)-DI13Ccal(n,i,j))*1.024
c      leak_DI13C= 0.5*( dabs(leak_DIC)*R13Ccoe/(1.+R13Ccoe)
c     &                      +leak_DIC *R13Ccoe/(1.+R13Ccoe)
c     &                 -dabs(leak_DIC)*R13Ccal/(1.+R13Ccal)
c     &                      +leak_DIC *R13Ccal/(1.+R13Ccal))
c
c      a_calc=exp((3.-0.140*Gn/0.011098)/1000.)    ! kinetic effect calcualtion: 1 mm y-1 
c                    -> (10cm*1.4g/cm) /(365*24*60*60s*40g/mol)* 10**6 = 0.011098 umol cm-2 s-1
c
      c13CH2O (n,i,j)=c13CH2O (n,i,j)
     &         +(Pg*R13Ccoe*a_phot/(1.+R13Ccoe*a_phot)
     &           -R*R13CH2O*a_resp/(1.+R13CH2O*a_resp)
     &          )/coe_tick/1.024 *dtc    !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1
c
      DI13Ccoe (n,i,j)=DI13Ccoe (n,i,j)
     &            +(-Pg*R13Ccoe*a_phot/(1.+R13Ccoe*a_phot)
     &              +R*R13CH2O*a_resp/(1.+R13CH2O*a_resp)
     &              +DI13Cuptake
     &              -F_13CO2
     &              -leak_DI13C
     &             )/coe_tick/1.024 *dtc   !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1
c
      DI13Ccal (n,i,j)=DI13Ccal (n,i,j)
     &            +(-Gn*R13Ccal*a_calc/(1.+R13Ccal*a_calc)
     &              +F_13CO2
     &              +leak_DI13C
     &             )/cal_tick/1.024 *dtc    !nmol cm-3 s-1 = 1./1.024 umol kg-1 s-1
c
c------ Calculation of d13C -------------------------
c
      R13CH2O=c13CH2O(n,i,j)/(CH2O(n,i,j)-c13CH2O(n,i,j))   !coral organism
      R13Ccoe=DI13Ccoe(n,i,j)/(DICcoe(n,i,j)-DI13Ccoe(n,i,j))   !DIC in coelenteron
      R13Ccal=DI13Ccal(n,i,j)/(DICcal(n,i,j)-DI13Ccal(n,i,j))   !DIC in calcifying fluid
c
      d13C_DICamb=d13C_fromR13C(R13Camb)
      d13C_DICcoe=d13C_fromR13C(R13Ccoe)
      d13C_CH2O=d13C_fromR13C(R13CH2O)
      d13C_DICcal=d13C_fromR13C(R13Ccal)
c
      R13Carg=R13Ccal*a_calc
c
      d13C_arg=d13C_fromR13C(R13Carg)
c



!------------------------------------------------------------------------
! Print section (for debug)
c
#if defined TESTMODE

      if(iprint.eq.0) then

        write(52,*) time,PFD
     &   ,Pg,R,Pn,Gn,E_ca
     &   ,TAcal(n,i,j),TAcoe(n,i,j),TAamb
     &   ,DICcal(n,i,j),DICcoe(n,i,j),DICamb
     &   ,DOcoe(n,i,j),DOamb
     &   ,pHcal,pHcoe,pHamb,Wargcal,Wargamb,fCO2cal,fCO2coe,fCO2amb
     &   ,d13C_DICamb,d13C_DICcoe,d13C_CH2O,d13C_DICcal,d13C_arg
     &   ,d13C_arg*Gn
     &   ,QCcoral(n,i,j),QNcoral(n,i,j),QPcoral(n,i,j)
     &   ,ROS(n,i,j)
     &   ,g_coral(n,i,j),m_coral(n,i,j)
     &   ,g_coral(n,i,j)-m_coral(n,i,j)
     
     &   ,F_ROS,F_detox,F_damage
     &   ,Damage
     &   ,TOCuptake,TONuptake,TOPuptake
     &   ,C_excre,C_growth,C_mucus,C_ingest
c     &   ,N_excre,N_growth,N_mucus,N_ingest
c     &   ,P_excre,P_growth,P_mucus,P_ingest
     &   ,p_coral(n,i,j)
     &   ,dsed_coral(n,i,j),Fsed, Rmsed
     
     &   ,NO3coral(n,i,j),NO2coral(n,i,j),NH4coral(n,i,j)
     &   ,NO3uptake,NO2uptake,NH4uptake
     &   ,F_ONOO
     &   ,PO4coral(n,i,j),PO4uptake
     &   ,N_dissim, DICuptake

c
      endif

#endif
c
c
!------ Output parameters (2) ------------------------------------------------

c------ Total flux (mmol m-2 s-1) --------------------------------------

      DICuptake=DICuptake*1.d-2 *p_coral(n,i,j)!(nmol cm-2 s-1)->(mmol m-2 s-1)
      TAuptake =TAuptake *1.d-2 *p_coral(n,i,j)!(nmol cm-2 s-1)->(mmol m-2 s-1)
      NH4uptake=NH4uptake*1.d-2 *p_coral(n,i,j)
      NO3uptake=NO3uptake*1.d-2 *p_coral(n,i,j)
      NO2uptake=NO2uptake*1.d-2 *p_coral(n,i,j)
      PO4uptake=PO4uptake*1.d-2 *p_coral(n,i,j)
      TOCuptake=TOCuptake*1.d-2 *p_coral(n,i,j)
      TONuptake=TONuptake*1.d-2 *p_coral(n,i,j)
      TOPuptake=TOPuptake*1.d-2 *p_coral(n,i,j)
      DOuptake =DOuptake *1.d-2 *p_coral(n,i,j)
      PP1uptake=PP1uptake*1.d-2 *p_coral(n,i,j)
      ZP1uptake=ZP1uptake*1.d-2 *p_coral(n,i,j)
c
      DI13Cuptake=DI13Cuptake*1.d-2 *p_coral(n,i,j)!(nmol cm-2 s-1)->(mmol m-2 s-1)


c
      
c Coral record calculation section
c
      S_PFD_dt       =S_PFD_dt+PFD*dtc
      S_Gn_dt         =S_Gn_dt+Gn*dtc
      S_d13CargxGn_dt=S_d13CargxGn_dt+d13C_arg*Gn*dtc
      S_d13CCH2O_dt  =S_d13CCH2O_dt+d13C_CH2O*dtc
      S_Pg_dt        =S_Pg_dt+Pg*dtc
      S_R_dt         =S_R_dt +R *dtc
      S_CH2O_dt      =S_CH2O_dt+CH2O(n,i,j)*dtc

#if defined TESTMODE

      if(time.ge.real(iday)) then

        write(53,*) iday
     &   ,S_PFD_dt*1.e-6  !Photon flux density (mol m-2 d-1)
     &   ,S_Gn_dt*1.e-3   !Calcification rate (mmol cm-2 d-1)
     &   ,S_Pg_dt*1.e-3   !Gross photosynthesis rate (mmol cm-2 d-1)
     &   ,S_R_dt*1.e-3    !Respiration rate (mmol cm-2 d-1)
     &   ,S_CH2O_dt/24./60./60. ! 1 day avaraged value of CH2O
     &   ,(S_Pg_dt-S_R_dt)*1.e-3  !Net photosynthesis rate (mmol cm-2 d-1)
     &   ,S_d13CargxGn_dt
     &   ,S_d13CargxGn_dt/S_Gn_dt   !d13C
     &   ,S_d13CCH2O_dt/24./60./60. ! 1 day avaraged value of d13C_CH2O
c
        iday=iday+1
        S_PFD_dt=0.
        S_Gn_dt=0.
        S_d13CargxGn_dt=0.
        S_d13CCH2O_dt=0.
        S_Pg_dt=0.
        S_R_dt=0.
        S_CH2O_dt=0.
c
      endif
#endif
ccc
      return
      end




!*************** Zooxanthellae ************************************************
c
      subroutine zooxanthellae
c          input parameters
     &            (n,i,j
     &            ,PFDsurf
     &            ,TKamb                           !Temperature (K)
     &            ,cCO2aqcoe, cHCO3coe
c     &            ,PO4uptake                       ! PO4 uptake (nmol cm-2 s-1)
     &            ,Expul                           ! Zooxanthellae release rate by host coral (individual cm-2 s-1)
     &            ,C_excre,N_excre,P_excre         ! Organic matter excretion controled by host coral (nmol cm-2 s-1)
c          output parameters
     &            ,Pg                              ! Gross photosynthesis rate (nmol cm-2 s-1)
     &            ,Rz                              ! Zooxanthellae respiration rate (nmol cm-2 s-1)
     &            ,NO3_trans, NO2_trans, NH4_trans ! Nutrients transport rate from host to Zooxanthellae (nmol cm-2 s-1)
     &            ,PO4_trans                       ! Nutrients transport rate from host to Zooxanthellae (nmol cm-2 s-1)
     &            ,C_waste,N_waste,P_waste         ! Waste organic matter (nmol cm-2 s-1)
     &            ,F_ROS                           ! Flux of Reactive Oxygen Species  (nmol cm-2 s-1)
     &            ,F_ONOO                          ! Flux of peroxynitrite  (nmol cm-2 s-1)
     &             )
c
C **********************************************************************
C *                                                                    *
C * FUNCTION    : Zooxanthellae population dynamics.                   *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom2k_v9.blk'

c          input parameters
      integer n,i,j
      double precision
     &             PFDsurf
     &            ,TKamb                          !Temperature (K)
     &            ,cCO2aqcoe, cHCO3coe
c     &            ,PO4uptake                      ! PO4 uptake (nmol cm-2 s-1)
     &            ,Expul                      ! Zooxanthellae release rate by host coral
     &            ,C_excre,N_excre,P_excre        ! Organic matter excretion controled by host coral
c          output parameters
     &            ,Pg                             ! Gross photosynthesis rate (nmol cm-2 s-1)
     &            ,Rz                             ! Zooxanthellae respiration rate (nmol cm-2 s-1)
     &            ,NO3_trans, NO2_trans, NH4_trans ! Nutrients transport rate from host to Zooxanthellae (nmol cm-2 s-1)
     &            ,PO4_trans                       ! Nutrients transport rate from host to Zooxanthellae (nmol cm-2 s-1)
     &            ,C_waste,N_waste,P_waste        ! Waste organic matter (nmol cm-2 s-1)
     &            ,F_ROS                          ! Flux of Reactive Oxygen Species  (nmol cm-2 s-1)
     &            ,F_ONOO                          ! Flux of peroxynitrite  (nmol cm-2 s-1)

c
      double precision  Rgas   ! Gas constant
c
c----- for gross photosynthesi -------------------
      double precision  lamb,Pgmax(ccom),Ik(ccom)
      double precision  K_HCO3(ccom),K_CO2aq(ccom)
      double precision  V_HCO3(ccom),V_CO2aq(ccom)
c      Temperature related factor
      double precision  f_temp
      double precision  Hact(ccom),Hdeact(ccom),Entr(ccom)
c
c----- Respiration parameters -------------------
      double precision  Rmax                     ! Respiration rate
      double precision  K_QC
      double precision  K_DO
      double precision  SQC

c----- Nutrients related parameters -------------------
      double precision  
     &    N_dissim       ! Nitrogen dissimilation rate (pmol cell-1 s-1)
     &   ,N_assim        ! Nitrogen assimilation rate (pmol cell-1 s-1)
c     &   ,P_dissim       ! P dissimilation rate (pmol cell-1 s-1)
     &   ,P_assim        ! P assimilation rate (pmol cell-1 s-1)
     &   ,NO3_reduc      ! NO3 reductance rate (pmol cell-1 s-1)
     &   ,NO2_reduc      ! NO2 reductance rate (pmol cell-1 s-1)

     &   ,Vmax_NO3tra
     &   ,Vmax_NO2tra
     &   ,Vmax_NH4tra
     &   ,Vmax_PO4tra
     
     &   ,K_NO3tra
     &   ,K_NO2tra
     &   ,K_NH4tra
     &   ,K_PO4tra

     &  ,k_NH4rel
     &  ,k_NO3rel
     &  ,k_NO2rel
     &  ,k_PO4rel

     &   ,Vmax_NO3red
     &   ,Vmax_NO2red
     &   ,K_NO3red
     &   ,K_NO2red

     &   ,Vmax_Nass
     &   ,K_Nass
     &   ,Vmax_Pass
     &   ,K_Pass
                                            
c----- Metabolic parameters -------------------
      double precision  C_repro,N_repro,P_repro
      double precision  u_max 

c----- Population dynamics parameters -------------------
      double precision  mortality, km_CNP, am_CNP, km_ROS, km_min
      double precision  Repro, Death


!------------Set Constants  ----------------------------------
c----- Physical constants ------------------------
c
      data Rgas  /8.314/ ! Gas constant (J mol-1)

c----- for gross photosynthesi -------------------
c
c  Kuhl et al. 1995, MEPS 
c     Pg=Pgmax*tanh(I/Ik) for Acropora sp.
c
c      Pg=Pgmax(n)*tanh(PFDsurf/Ik(n))*(1-Br(i,j))
c     &       *cHCO3coe/(K_HCO3+cHCO3coe)
c     &       *coe_tick
c
      data lamb    /1000./!*適当(未使用）
c
c      data Pgmax   /1.0/ ! 11.*coe_tick (nmolO2 cm-3 s-1); 0.47 (nmol O2 cm-2 s-1) Kuhl et al. (1995) *Tuned
      data Pgmax   /1.e-3/ !1.0 (nmol cm-2 -s)=1.0/DENSzoox=1.0e-6 (nmol cell-1 s-1) =1.0e-3 (pmol cell-1 s-1)
      data Ik      /275./    !(275 uEin m-2 s-1) * refference value
c
!------ Goreau et al. (1996)------------------------------------------   
!                    Coral    FIZ    CZ
!          K_HCO3   408.    71.     178.   uM
!          V_HCO3     1.45   1.09    13.   nmol O2 (ug Chl a)-1 min-1
!          K_CO2aq    5.     0.67     2.7  uM
!          V_CO2aq    2.     1.15    12.5  nmol O2 (ug Chl a)-1 min-1
!---------------------------------------------------------------------
      data K_HCO3  /408./   !uM
      data V_HCO3  /1.45/  !nmol O2 (ug Chl a)-1 min-1
      data K_CO2aq /5./!/5./  !uM
      data V_CO2aq /2./!/2./  !nmol O2 (ug Chl a)-1 min-1
c
c     Temperature dependency parameters (Hikosaka et al. 2006) 
c        - Parameters estimated from data by Al-Horani (2005)
c
      data Hact   /1.855e5/ !Activation energy (J mol-1)) *tuned
      data Hdeact /3.055e5/ !Deactivation energy (J mol-1) *tuned
      data Entr   /1028./   !Entropy term (J K-1) *tuned
c
c----- Respitration constants ---------------------------
c
      data Rmax    /3.e-5/!/0.37/ !Maximum respiration rate (pmol cell-1 s-1) *Tuned using Kuhl et al (1995), Al-Horani et al. (2003) data
      data K_QC    /8./ ! (pmol cell-1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*適当
      data K_DO    /5./ ! ca. 10 (umol kg-1) Shick et al. (1990) *Tuned value


c----- Nutrients related parameters -------------------

      data Vmax_NO3tra  /1.67e-6/  ! (pmol cell-1 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう 
      data Vmax_NO2tra  /1.5e-7/  ! (pmol cell-1 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data Vmax_NH4tra  /1.e-5/  ! (pmol cell-1 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data Vmax_PO4tra  /3.e-7/  ! (pmol cell-1 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう

      data K_NO3tra  /50./  ! (umol kg-1)  !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data K_NO2tra  /10./  ! (umol kg-1)  !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data K_NH4tra  /300./ ! (umol kg-1)  !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう      
      data K_PO4tra  /300./ ! (umol kg-1)  !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう      

      data k_NH4rel  /1.e-5/  !(s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data k_NO3rel  /1.9e-8/  !(s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data k_NO2rel  /1.e-5/  !(s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data k_PO4rel  /1.e-5/  !(s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう



      data Vmax_NO3red  /1.3e-6/  ! (pmol cell-1 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう 
      data Vmax_NO2red  /1.3e-6/  ! (pmol cell-1 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう

      data K_NO3red  /0.014/  !(pmol cell-1)  !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう
      data K_NO2red  /0.0015 /  !(pmol cell-1)  !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう

      data Vmax_Nass  /2.e-5/  ! (pmol cell-1 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう 
      data K_Nass  /1.e-1/  !(pmol cell-1)  !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう

      data Vmax_Pass  /5.e-7/  ! (pmol cell-1 s-1) !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう 
      data K_Pass  /1.e-2/  !(pmol cell-1)  !!!!!!!!!!!!!!!!!!!!!!!!!!てきとう

c
c----- Population dynamics parameters -------------------

      data u_max    / 2.e-5/ !Maximum reproduction rate (pmol cell-1 s-1)  !!!!!!!!!!!!!!!!!!適当!!!!!!!!!!!!!!
      
      data km_CNP   /  1.e-8/ ! (cell cm-2 s-1)  適当!!!!!!!!!!!!!!
      data am_CNP   /  1.e0/
      data km_ROS   / 1.e-11/ !                  適当!!!!!!!!!!!!!!
      data km_min   /  1.e-14/   ! (cell cm-2 s-1)  適当!!!!!!!!!!!!!!
      
      

!----- Zooxanthellae metabolism  ----------------------------------

c----- Gross photosynthesis rate (pmol cell-1 s-1) -----------------
c

c
c   Temperature dependent term
      f_temp=1.
!      f_temp=exp(Hact(n)*(TKamb-298.)/298./Rgas/TKamb)
     &       *(1.+exp((298.*Entr(n)-Hdeact(n))/298./Rgas))
     &       /(1.+exp((TKamb*Entr(n)-Hdeact(n))/TKamb/Rgas))

c   Gross photosynthesis rate (pmol cell-1 s-1)
      Pg=Pgmax(n)*f_temp*tanh(PFDsurf/Ik(n))
     &     *1./(V_HCO3(n)+V_CO2aq(n))
     &       *( V_HCO3(n)*cHCO3coe/(K_HCO3(n)+cHCO3coe)
     &         +V_CO2aq(n)*cCO2aqcoe/(K_CO2aq(n)+cCO2aqcoe)
     &        )

c
c----- Respilation rate (pmol cell-1 s-1) ---------------------------

c   Temperature dependent term
c      f_temp=1.
c      f_temp=exp(Eact(n)*(TKamb-298.)/298./Rgas/TKamb)
c
cc   Respilation rate (nmolO2 cm-2 s-1) 
c      R= Rmax(n)*f_temp*CH2O(n,i,j)/(K_CH2O(n)+CH2O(n,i,j))
c     &            *DOcoe(n,i,j)/(K_DO(n)+DOcoe(n,i,j))
c
      Rz= Rmax*QCzoox(n,i,j)/(K_QC+QCzoox(n,i,j)) ! Respilation rate (pmol cell-1 s-1)
     &         *DOcoe(n,i,j)/(K_DO+DOcoe(n,i,j))

c      SQC=QCzoox(n,i,j)-QC0zoox
c
c      Rz= Rmax*SQC/(K_QC+SQC)
c     &         *DOcoe(n,i,j)/(K_DO+DOcoe(n,i,j))

!----- Nutrient internal dynamics -----------------------------------

c----- Nutrients transport rate (pmol cell-1 s-1) ---------------------

c  These formulations are basically followed to Flynn and Flynn (1998)

      NH4_trans =Vmax_NH4tra*NH4coral(n,i,j)/(K_NH4tra+NH4coral(n,i,j))

c     &          -2.* Vmax_NH4tra
c     &             *(NH4zoox(n,j,i)/NH4zoox_Max)**4.
c     &               /((NH4zoox(n,j,i)/NH4zoox_Max)**4. + 1.)

     &          -k_NH4rel * NH4zoox(n,j,i)

      NO3_trans =Vmax_NO3tra*NO3coral(n,i,j)/(K_NO3tra+NO3coral(n,i,j))
c     &               *(1.-(NH4zoox(n,i,j)/NH4zoox_Max)**1.)
     &               *(1.-NH4zoox(n,i,j)/0.1)

c     &          -2.* Vmax_NO3tra
c     &             *(NO3zoox(n,j,i)/NO3zoox_Max)**4.
c     &               /((NO3zoox(n,j,i)/NO3zoox_Max)**4. + 1.)

     &          -k_NO3rel * NO3zoox(n,j,i)

      NO2_trans =Vmax_NO2tra*NO2coral(n,i,j)/(K_NO2tra+NO2coral(n,i,j))
c     &               *(1.-(NH4zoox(n,i,j)/NH4zoox_Max)**1.)
     &               *(1.-NH4zoox(n,i,j)/0.1)
c
c     &          -2.* Vmax_NO2tra
c     &             *(NO2zoox(n,j,i)/NO2zoox_Max)**4.
c     &               /((NO2zoox(n,j,i)/NO2zoox_Max)**4. + 1.)

     &          -k_NO2rel * NO2zoox(n,j,i)

      PO4_trans =Vmax_PO4tra*PO4coral(n,i,j)/(K_PO4tra+PO4coral(n,i,j))

c     &          -2.* Vmax_PO4tra
c     &             *(PO4zoox(n,j,i)/PO4zoox_Max)**4.
c     &               /((PO4zoox(n,j,i)/PO4zoox_Max)**4. + 1.)

c     &          -k_PO4rel * PO4zoox(n,j,i)

c----- Nutrients reductance rate (pmol cell-1 s-1) ---------------------

      NO3_reduc = Vmax_NO3red*NO3zoox(n,i,j)/(K_NO3red+NO3zoox(n,i,j))
      NO2_reduc = Vmax_NO2red*NO2zoox(n,i,j)/(K_NO2red+NO2zoox(n,i,j))

c----- Nitrogen assimilation rate (pmol cell-1 s-1) ---------------------

      N_assim = Vmax_Nass* NH4zoox(n,i,j)/(K_Nass+NH4zoox(n,i,j)) 


c----- P assimilation rate (pmol cell-1 s-1) ---------------------

      P_assim = Vmax_Pass* PO4zoox(n,i,j)/(K_Pass+PO4zoox(n,i,j)) 


c----- Nitrogen dissimilation rate (pmol cell-1 s-1) ---------------------

      N_dissim = Rz * QNzoox(n,i,j)/QCzoox(n,i,j)   ! Dissimilation rate (pmol cell-1 s-1)


!----- interaction between zooxanthellae and coral -----------------

c----- Reproduction rate (pmol cell-1 s-1) ----------------------------

      C_repro=u_max*min(    1.-QP0zoox/QPzoox(n,i,j)  ,
     &                  min(1.-QN0zoox/QNzoox(n,i,j),
     &                      1.-QC0zoox/QCzoox(n,i,j) ) )
     
      if(C_repro .lt. 0.) then
        C_repro =0.
      end if
      
      N_repro=C_repro*QNzoox(n,i,j)/QCzoox(n,i,j)
      
      P_repro=C_repro*QPzoox(n,i,j)/QCzoox(n,i,j)


!----- Convert cell conditions to population dynamics parameters

      Repro=C_repro/QCzoox(n,i,j)*DENSzoox(n,i,j)  ! Reproduction rate (cell cm-2 s-1)

      mortality= km_CNP                          ! mortality (s-1) 
     &            *exp(-am_CNP*                                      !* 足りない栄養素に依存（要検討!!!）
     &                min(  1.-QP0zoox/QPzoox(n,i,j)  ,
     &                  min(1.-QN0zoox/QNzoox(n,i,j),
     &                      1.-QC0zoox/QCzoox(n,i,j) ) )
     &             )
c     &          +km_ROS * ROS(n,i,j)                                !* ROS濃度依存（要検討!!!）  今のところゼロにしている
c     &          +km_min                                    !* 定数

      Death=mortality*DENSzoox(n,i,j)              ! Death rate (cell cm-2 s-1)


!--------Time step progressing ------------------------

c---- Nutrients (pmol/cell) ------------------

      NO3zoox(n,i,j)=NO3zoox(n,i,j)
     &              +(
     &                 NO3_trans
     &                -NO3_reduc
     &               )*dtc
      NO3zoox(n,i,j)=max(NO3zoox(n,i,j),0.d0) !Error handring

      NO2zoox(n,i,j)=NO2zoox(n,i,j)
     &              +(
     &                 NO2_trans
     &                +NO3_reduc
     &                -NO2_reduc
     &               )*dtc
      NO2zoox(n,i,j)=max(NO2zoox(n,i,j),0.d0) !Error handring

      NH4zoox(n,i,j)=NH4zoox(n,i,j)
     &              +(
     &                 NH4_trans
     &                +NO2_reduc
     &                -N_assim
     &                +N_dissim
     &               )*dtc
      NH4zoox(n,i,j)=max(NH4zoox(n,i,j),0.d0) !Error handring
      

      PO4zoox(n,i,j)=PO4zoox(n,i,j)
     &              +(
     &                 PO4_trans
     &                -P_assim
     &               )*dtc
      PO4zoox(n,i,j)=max(PO4zoox(n,i,j),0.d0) !Error handring

c-----Cell quota (pmol/cell) -----------------!!!!!単位変換チェック!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      QCzoox(n,i,j)=QCzoox(n,i,j)
     &              +(
     &                 Pg-Rz
     &                -C_excre/DENSzoox(n,i,j)*1.e3      !(nmol cm-2 s-1)/(cell cm-2)*1.e3=(pmol cell-1 s-1)
     &                -C_repro
     &               )*dtc
      QCzoox(n,i,j)=max(QCzoox(n,i,j),0.d0) !Error handring
           
      QNzoox(n,i,j)=QNzoox(n,i,j)
     &              +(
     &                 N_assim
     &                -N_dissim
c     &                 NO3uptake/DENSzoox(n,i,j)*1.e3    !(nmol cm-2 s-1)/(cell cm-2)*1.e3=(pmol cell-1 s-1)
c     &                +NH4uptake/DENSzoox(n,i,j)*1.e3    !(nmol cm-2 s-1)/(cell cm-2)*1.e3=(pmol cell-1 s-1)
     &                -N_excre/DENSzoox(n,i,j)*1.e3      !(nmol cm-2 s-1)/(cell cm-2)*1.e3=(pmol cell-1 s-1)
     &                -N_repro
     &               )*dtc
      QNzoox(n,i,j)=max(QNzoox(n,i,j),0.d0) !Error handring     
      
      QPzoox(n,i,j)=QPzoox(n,i,j)
     &              +(
c     &                 PO4uptake/DENSzoox(n,i,j)*1.e3    !(nmol cm-2 s-1)/(cell cm-2)*1.e3=(pmol cell-1 s-1)
     &                 P_assim
     &                -P_excre/DENSzoox(n,i,j)*1.e3      !(nmol cm-2 s-1)/(cell cm-2)*1.e3=(pmol cell-1 s-1)
     &                -P_repro
     &               )*dtc
      QPzoox(n,i,j)=max(QPzoox(n,i,j),0.d0) !Error handring

c-----Zooxanthellae density (cell cm-2) -----------------
c
c     DENSzoox ~ 1.e6 cells cm-2  (e.g. Higuchi et al., 2012)
c
      DENSzoox(n,i,j)=DENSzoox(n,i,j)
     &             +(
     &                Repro   ! Reproduction
     &               -Death   ! Death
     &               -Expul   ! Expulsion by host coral
     &              )*dtc
      DENSzoox(n,i,j)=max(DENSzoox(n,i,j),0.d0) !Error handring

!-----Output parameters  (nmol cm-2 s-1) --------------------------------

c----- ROS flux (nmol cm-2 s-1)------------------------------------
c
c  Estimated from Saragosti et al. (2010) data
c      
c      F_ROS=0.002382*PFDsurf+0.006891*exp(0.1395*(TKamb-273.15))   ! (fmol min-1 cell-1)
      
c      F_ROS=(0.000354*PFDsurf+0.0672)*exp(0.0643*(TKamb-273.15))   ! (fmol min-1 cell-1)


c      F_ROS=(0.000354*PFDsurf+0.0672)*exp(0.0643*(TKamb-273.15))   ! (fmol min-1 cell-1)

      F_ROS=(0.000354*PFDsurf+0.0672)*exp(0.4*(TKamb-273.15))*2.e-4   ! (fmol min-1 cell-1)

      
      F_ROS=F_ROS*DENSzoox(n,i,j)/60.*1.e-6   ! (nmol cm-2 s-1)
      
cc----- peroxynitrite (ONOO-) flux (nmol cm-2 s-1)------------------------------------
cc
c      F_ONOO=1.e0 * NO2zoox(n,i,j)*ROS(n,i,j)   ! (pmol cell-1 s-1)!!!!!!!!!!!!!!!!!!てきとう
c      
c      F_ONOO=F_ONOO*DENSzoox(n,i,j)*1.e-3   ! (nmol cm-2 s-1)
c


!------------------------------------------------------------------------
! Print section (for debug)
c
#if defined TESTMODE

      if(iprint.eq.0) then

        write(54,*) 
     &    time,PFDsurf,Pg,Rz
     &   ,DENSzoox(n,i,j),QCzoox(n,i,j),QNzoox(n,i,j),QPzoox(n,i,j)
     &   ,ROS(n,i,j)
     &   ,C_repro,N_repro,P_repro
     &   ,Repro,Death,Expul, F_ROS
     &   ,NO3zoox(n,i,j),NO2zoox(n,i,j),NH4zoox(n,i,j),F_ONOO
     &   ,NO3_trans,NO2_trans,NH4_trans
     &   ,NO3_reduc,NO2_reduc
     &   ,N_assim,N_dissim
     &   ,PO4zoox(n,i,j),PO4_trans,P_assim
c
      endif

#endif
c      
c
!-----------------------------------------------------------------------

c
c
c----- Outoput for coral internal system model (nmol cm-2 s-1)------------------------------------

      Pg=Pg*DENSzoox(n,i,j)*1.e-3  !=(pmol cell-1 s-1)*(cell cm-1)*1.e-3=(nmol cm-2 s-1)
      Rz=Rz*DENSzoox(n,i,j)*1.e-3  !=(pmol cell-1 s-1)*(cell cm-1)*1.e-3=(nmol cm-2 s-1)

      NO3_trans=NO3_trans*DENSzoox(n,i,j)*1.e-3  !=(pmol cell-1 s-1)*(cell cm-1)*1.e-3=(nmol cm-2 s-1)
      NO2_trans=NO2_trans*DENSzoox(n,i,j)*1.e-3  !=(pmol cell-1 s-1)*(cell cm-1)*1.e-3=(nmol cm-2 s-1)
      NH4_trans=NH4_trans*DENSzoox(n,i,j)*1.e-3  !=(pmol cell-1 s-1)*(cell cm-1)*1.e-3=(nmol cm-2 s-1)
      PO4_trans=PO4_trans*DENSzoox(n,i,j)*1.e-3  !=(pmol cell-1 s-1)*(cell cm-1)*1.e-3=(nmol cm-2 s-1)

      C_waste = (Death+Expul)*QCzoox(n,i,j)*1.e-3  !=(cell cm-2 s-1)*(pmol cell-1)*1.e-3=(nmol cm-2 s-1)
      N_waste = (Death+Expul)*QNzoox(n,i,j)*1.e-3
      P_waste = (Death+Expul)*QPzoox(n,i,j)*1.e-3
      

c---------------------------------------------------------------------

      return
      end
cc


