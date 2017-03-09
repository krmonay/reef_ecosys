
!!!=== ver 2017/03/09   Copyright (c) 2017 Takashi NAKAMURA  =====

!--------------------------------------------------------------------------------
!
!              Output module
!
!--------------------------------------------------------------------------------

#include "cppdefs.h"

  MODULE mod_output
  
  contains
  
! **********************************************************************
!  Open files
! **********************************************************************

    SUBROUTINE files_open
    
      implicit none

#if defined CHAMBER_SITE4
      open(10,file='./output/site04-env_his.csv')
# if defined CORAL_TESTMODE
      open(11,file='./output/site04-crl1_his.csv')
      open(12,file='./output/site04-crl2_his.csv')
      open(21,file='./output/site04-crl1_ave.csv')
      open(22,file='./output/site04-crl2_ave.csv')
      open(31,file='./output/site04-zoo1_his.csv')
# endif
# if defined ECOSYS_TESTMODE
      open(40,file='./output/site04-ecosys_his.csv')
# endif
#elif defined CHAMBER_SITE5
      open(10,file='./output/site05-env_his.csv')
# if defined CORAL_TESTMODE
      open(11,file='./output/site05-crl1_his.csv')
      open(12,file='./output/site05-crl2_his.csv')
      open(21,file='./output/site05-crl1_ave.csv')
      open(22,file='./output/site05-crl2_ave.csv')
      open(31,file='./output/site05-zoo1_his.csv')
# endif
# if defined ECOSYS_TESTMODE
      open(40,file='./output/site05-ecosys_his.csv')
# endif
#elif defined CHAMBER_SITE6
      open(10,file='./output/site06-env_his.csv')
# if defined CORAL_TESTMODE
      open(11,file='./output/site06-crl1_his.csv')
      open(12,file='./output/site06-crl2_his.csv')
      open(21,file='./output/site06-crl1_ave.csv')
      open(22,file='./output/site06-crl2_ave.csv')
      open(31,file='./output/site06-zoo1_his.csv')
# endif
# if defined ECOSYS_TESTMODE
      open(40,file='./output/site06-ecosys_his.csv')
# endif
#elif defined CHAMBER_SITE7
      open(10,file='./output/site07-env_his.csv')
# if defined CORAL_TESTMODE
      open(11,file='./output/site07-crl1_his.csv')
      open(12,file='./output/site07-crl2_his.csv')
      open(21,file='./output/site07-crl1_ave.csv')
      open(22,file='./output/site07-crl2_ave.csv')
      open(31,file='./output/site07-zoo1_his.csv')
# endif
# if defined ECOSYS_TESTMODE
      open(40,file='./output/site07-ecosys_his.csv')
# endif
#elif defined CHAMBER_SITE9
      open(10,file='./output/site09-env_his.csv')
# if defined CORAL_TESTMODE
      open(11,file='./output/site09-crl1_his.csv')
      open(12,file='./output/site09-crl2_his.csv')
      open(21,file='./output/site09-crl1_ave.csv')
      open(22,file='./output/site09-crl2_ave.csv')
      open(31,file='./output/site09-zoo1_his.csv')
# endif
# if defined ECOSYS_TESTMODE
      open(40,file='./output/site09-ecosys_his.csv')
# endif
#elif defined CHAMBER_SITE10
      open(10,file='./output/site10-env_his.csv')
# if defined CORAL_TESTMODE
      open(11,file='./output/site10-crl1_his.csv')
      open(12,file='./output/site10-crl2_his.csv')
      open(21,file='./output/site10-crl1_ave.csv')
      open(22,file='./output/site10-crl2_ave.csv')
      open(31,file='./output/site10-zoo1_his.csv')
# endif
# if defined ECOSYS_TESTMODE
      open(40,file='./output/site10-ecosys_his.csv')
# endif


#else
      open(10,file='./output/eco5-env_his.csv')
# if defined CORAL_TESTMODE
      open(11,file='./output/eco5-crl1_his.csv')
      open(12,file='./output/eco5-crl2_his.csv')
      open(21,file='./output/eco5-crl1_ave.csv')
      open(22,file='./output/eco5-crl2_ave.csv')
      open(31,file='./output/eco5-zoo1_his.csv')
# endif
# if defined ECOSYS_TESTMODE
      open(40,file='./output/eco5-ecosys_his.csv')
# endif

#endif
      
#if defined SEDIMENT_TESTMODE
      open(56,file='./output/eco5-sedDIC_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(57,file='./output/eco5-sedTA_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(58,file='./output/eco5-sedDO_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(59,file='./output/eco5-sedpH_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(60,file='./output/eco5-sedWarg_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(61,file='./output/eco5-sedNH4_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(62,file='./output/eco5-sedNO2_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(63,file='./output/eco5-sedNO3_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(64,file='./output/eco5-sedPO4_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(65,file='./output/eco5-sedDOC_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(66,file='./output/eco5-sedPOC_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(67,file='./output/eco5-sedDON_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(68,file='./output/eco5-sedPON_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(69,file='./output/eco5-sedDOP_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(70,file='./output/eco5-sedPOP_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(71,file='./output/eco5-sedPg_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(72,file='./output/eco5-sedRdoc_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(73,file='./output/eco5-sedRpoc_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(74,file='./output/eco5-sedGn_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(75,file='./output/eco5-sedNit1_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(76,file='./output/eco5-sedNit2_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(78,file='./output/eco5-sedDNd_his.txt')!!!!!!!!!!!!!!!!!!!for debug
      open(79,file='./output/eco5-sedDNp_his.txt')!!!!!!!!!!!!!!!!!!!for debug
#endif

      RETURN

    END SUBROUTINE files_open

! **********************************************************************
!  Close files
! **********************************************************************

    SUBROUTINE files_close
    
      implicit none
      
      close(10)
#if defined CORAL_TESTMODE
      close(11)
      close(12)
      close(21)
      close(22)
      close(31)
#endif
#if defined ECOSYS_TESTMODE
      close(40)
#endif
#if defined SEDIMENT_TESTMODE
      close(56)
      close(57)
      close(58)
      close(59)
      close(60)
      close(61)
      close(62)
      close(63)
      close(64)
      close(65)
      close(66)
      close(67)
      close(68)
      close(69)
      close(70)
      close(71)
      close(72)
      close(73)
      close(74)
      close(75)
      close(76)
      close(78)
      close(79)
#endif

      RETURN

    END SUBROUTINE files_close

! **********************************************************************
!  Write lavel of environmental data
! **********************************************************************

    SUBROUTINE write_env_lavel(fid)
    
      USE mod_param
      
      implicit none
      
      integer, intent(in) :: fid
      
      write(fid,*) 'time, ', 'PFDsurf, '                               &
     &    ,'Temp, ','Salt, ','TA, ','DIC, ','DO, '                     &
#if defined ORGANIC_MATTER
     &    ,'DOC, ','POC, '     &
     &    ,'Phyt1, ','Phyt2, ','Zoop, '                                &
#endif
#if defined NUTRIENTS            
     &    ,'NO3, ','NO2, ','NH4, ','PO4, '                             &
# if defined ORGANIC_MATTER
     &    ,'DON, ','PON, ','DOP, ','POP, '                             &
# endif
#endif
#if defined COT_STARFISH
     &    ,'COT, ','COTl '                                             &
# endif
     &    ,'dz, ', 'etide, ', 'ereef'
      

      RETURN

    END SUBROUTINE write_env_lavel

! **********************************************************************
!  Write timeseries environmental data
! **********************************************************************

    SUBROUTINE write_env_data(fid)
    
      USE mod_param
#ifdef CORAL_POLYP
      USE mod_coral
#endif
#ifdef SEAGRASS
      USE mod_seagrass
#endif
#ifdef MACROALGAE
      USE mod_macroalgae
#endif
#ifdef SEDIMENT_ECOSYS
      USE mod_sedecosys
#endif
#if defined CARBON_ISOTOPE
      USE mod_geochem
#endif
      
      implicit none
      
      integer, intent(in) :: fid
#if defined CARBON_ISOTOPE
      real(8) :: d13C_DIC

      d13C_DIC=d13C_fromR13C(C(1,1,1,1,iT13C)/C(1,1,1,1,iTIC_))
#endif
        
      write(fid,*) time,',', PFDsurf,','                              &
     &       ,C(1,1,1,1,iTemp),',',C(1,1,1,1,iSalt),','              &
     &       ,C(1,1,1,1,iTAlk),',',C(1,1,1,1,iTIC_),','              &
     &       ,C(1,1,1,1,iOxyg),','                                   &
#if defined ORGANIC_MATTER
     &       ,C(1,1,1,1,iDOC_), ',',C(1,1,1,1,iPOC_),','             &
     &       ,C(1,1,1,1,iPhyt1),',',C(1,1,1,1,iPhyt2),','            &
     &       ,C(1,1,1,1,iZoop), ','                                  &
#endif
#if defined NUTRIENTS            
     &       ,C(1,1,1,1,iNO3_),',',C(1,1,1,1,iNO2_),','              &
     &       ,C(1,1,1,1,iNH4_),',',C(1,1,1,1,iPO4_),','              &
# if defined ORGANIC_MATTER
     &       ,C(1,1,1,1,iDON_),',',C(1,1,1,1,iPON_),','              &
     &       ,C(1,1,1,1,iDOP_),',',C(1,1,1,1,iPOP_),','              &
# endif
#endif
#if defined CARBON_ISOTOPE
     &       ,d13C_DIC,','                                           &
#endif
#if defined COT_STARFISH
     &       ,C(1,1,1,1,iCOTe),',',C(1,1,1,1,iCOTl),','              &
# endif
     &       ,dz(1,1,1),',',etide,',', ereef


      RETURN

    END SUBROUTINE write_env_data
    
#if defined CORAL_TESTMODE
! **********************************************************************
!  Write lavel of coral internal conditions
! **********************************************************************

    SUBROUTINE write_crl_his_lavel(fid)
    
      USE mod_param
      
      implicit none
      
      integer, intent(in) :: fid
      
      write(fid,*) 'time,', 'PFD,'                                &
     &   ,'Pg,', 'R,', 'Pn,', 'G,','QC,'                          &
     &   ,'TAcal,',  'TAcoe,',  'TAamb,'                          &
     &   ,'DICcal,', 'DICcoe,', 'DICamb,'                         &
     &   ,'DOcoe,', 'DOamb,'                                      &
     &   ,'pHcal,','pHcoe, ','pHamb,','Wacal,','Waamb,'           &
     &   ,'fCO2cal,','fCO2coe,','fCO2amb,'                        &
     &   ,'CO2aqcal,','HCO3cal,','CO3cal,'                        &
     &   ,'CO2aqcoe,','HCO3coe,','CO3coe,'                        &
# if defined CORAL_CARBON_ISOTOPE
     &   ,'d13C_DICamb,','d13C_DICcoe,','d13C_QC,d13C_DICcal,'    &
     &   ,'d13C_arg','d13C_arg*Gn,'                               &
     &   ,'d13C_CO2aqcal,','d13C_HCO3cal,','d13C_CO3cal,'         &
     &   ,'d13C_CO2aqcoe,','d13C_HCO3coe,','d13C_CO3coe,'         &
     &   ,'13CO2aqcal,','H13CO3cal,','13CO3cal,'                  &
     &   ,'13CO2aqcoe,','H13CO3coe,','13CO3coe,'                  &
# endif
# if defined CORAL_BORON_ISOTOPE
     &   ,'d11Barg,'                                              &
# endif
# if defined CORAL_MUCUS
     &   ,'DOCuptake,'                                            &
# endif
# if defined CORAL_INGESTION
     &   ,'ZOOuptake,'                                            &
# endif
# if defined CORAL_SIZE_DYNAMICS
     &   ,'growth,','mort,','Damage,','F_Cgrowth,'                &
     &   ,'F_damage,','F_detox,'                                  &
# endif
     &   ,'E_ca'

      RETURN

    END SUBROUTINE write_crl_his_lavel
    
    
    SUBROUTINE write_crl_ave_lavel(fid)
    
      USE mod_param
      
      implicit none
      
      integer, intent(in) :: fid
      
        write(fid,*) 'day,'         &
     &   ,'S_PFD_dt,'               &   !Photon flux density (mol m-2 d-1)
     &   ,'S_Gn_dt,'                &   !Calcification rate (umol cm-2 d-1)
# if defined CORAL_CARBON_ISOTOPE
     &   ,'S_d13CargxGn_dt,'        &
     &   ,'d13Carg,'                &   !d13C
     &   ,'S_d13C_QC_dt,'           &   ! 1 day avaraged value of d13C_QC
# endif
# if defined CORAL_BORON_ISOTOPE
     &   ,'d11Barg,'   & 
# endif
     &   ,'S_Pg_dt,'                &   !Gross photosynthesis rate (umol cm-2 d-1)
     &   ,'S_R_dt,'                 &   !Respiration rate (umol cm-2 d-1)
     &   ,'S_QC_dt,'                &   ! 1 day avaraged value of QC
     &   ,'S_Pn_dt'                     !Net photosynthesis rate (umol cm-2 d-1)

      RETURN

    END SUBROUTINE write_crl_ave_lavel
#endif
#if defined ECOSYS_TESTMODE
! **********************************************************************
!  Write lavel of coral internal conditions
! **********************************************************************

    SUBROUTINE write_ecosys_his_lavel(fid)
    
      USE mod_param
      
      implicit none
      
      integer, intent(in) :: fid
      
      write(fid,*) 'time,', 'PFDbott,'                                 &
#ifdef CORAL_POLYP
     &    ,'coral1_Pg, ', 'coral1_R, ', 'coral1_Pn, ', 'coral1_G, '    &
     &    ,'coral2_Pg, ', 'coral2_R, ', 'coral2_Pn, ', 'coral2_G, '    &
#endif
#ifdef SEAGRASS
     &    ,'sgrass_Pg, ', 'sgrass_R, ', 'sgrass_Pn, '                  &
#endif
#ifdef MACROALGAE
     &    ,'algae_Pg, ' , 'algae_R, ' , 'algae_Pn, '                   &
#endif
#ifdef SEDIMENT_ECOSYS
     &    ,'sedeco_Pg, ', 'sedeco_R, ', 'sedeco_Pn, ', 'sedeco_G, '    &
#endif
     &    ,'dDIC_dt,','dTA_dt,','dDOx_dt,'                             &
# if defined ORGANIC_MATTER
     &    ,'dDOC_dt,','dPOC_dt(1),'                                    &
# endif
# if defined CARBON_ISOTOPE
     &    ,'dDI13C_dt(1),'                                             &
# endif
# if defined NUTRIENTS
     &    ,'dNO3_dt,','dNO2_dt,','dNH4_dt,'                            &
     &    ,dPO4_dt(1),','                                              &
#  if defined ORGANIC_MATTER
     &    ,'dDON_dt,','dPON_dt,','dDOP_dt,','dPOP_dt,'                 &
#  endif
# endif
     &    ,'pH, ', 'fCO2, ', 'Warg, '                                  &
     &    ,'U10, ', 'CO2flux, ' , 'O2flux'

      RETURN

    END SUBROUTINE write_ecosys_his_lavel
#endif
      
  END MODULE mod_output

