
!!!=== ver 2017/03/07   Copyright (c) 2017 Takashi NAKAMURA  =====

!--------------------------------------------------------------------------------
!
!              Output module
!
!--------------------------------------------------------------------------------

#include "cppdefs.h"

  MODULE mod_output
  
  contains

! **********************************************************************
!  Write lavel of environmental data
! **********************************************************************

    SUBROUTINE write_env_lavel(fid)
    
      USE mod_param
      
      implicit none
      
      integer, intent(in) :: fid
      
      write(fid,*) 'time, ', 'PFDsurf, '                               &
#ifdef CORAL_POLYP
     &    ,'coral1_Pg, ', 'coral1_R, ', 'coral1_Pn, ', 'coral1_G, '    &
     &    ,'coral1_QC, '                                               &
     &    ,'coral2_Pg, ', 'coral2_R, ', 'coral2_Pn, ', 'coral2_G, '    &
     &    ,'coral2_QC, '                                               &
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
     &    ,'pH, ', 'fCO2, ', 'Warg, '                                  &
     &    ,'U10, ', 'CO2flux, ' , 'O2flux, '                           &
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
#ifdef CORAL_POLYP
     &       ,CORAL(1)%Pg(1,1,1),',', CORAL(1)%R (1,1,1),','         &
     &       ,CORAL(1)%Pg(1,1,1)-CORAL(1)%R (1,1,1),','              &
     &       ,CORAL(1)%G (1,1,1),',',CORAL(1)%QC(1,1,1),','          &
     &       ,CORAL(1)%Pg(2,1,1),',', CORAL(1)%R (2,1,1),','         &
     &       ,CORAL(1)%Pg(2,1,1)-CORAL(1)%R (2,1,1),','              &
     &       ,CORAL(1)%G (2,1,1),',',CORAL(1)%QC(2,1,1),','          &
#endif
#ifdef SEAGRASS
     &       ,SGRASS(1)%Pg(1,1,1),',', SGRASS(1)%R (1,1,1),','       &
     &       ,SGRASS(1)%Pg(1,1,1)-SGRASS(1)%R (1,1,1),','            &
#endif
#ifdef MACROALGAE
     &       ,ALGAE(1)%Pg(1,1,1),',', ALGAE(1)%R (1,1,1),','         &
     &       ,ALGAE(1)%Pg(1,1,1)-ALGAE(1)%R (1,1,1),','              &
#endif
#ifdef SEDIMENT_ECOSYS
     &       ,SEDECO(1)%Pg(1,1),',', SEDECO(1)%R (1,1),','           &
     &       ,SEDECO(1)%Pg(1,1)-SEDECO(1)%R (1,1),','                &
     &       ,SEDECO(1)%G (1,1),','                                  &
#endif
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
     &       ,sspH,',', ssfCO2,',', ssWarg,','                       &
     &       ,U10,',',ssCO2flux,',', ssO2flux,','                    &
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
     &   ,'S_Gn_dt(n),'             &   !Calcification rate (umol cm-2 d-1)
# if defined CORAL_CARBON_ISOTOPE
     &   ,'S_d13CargxGn_dt,'        &
     &   ,'d13Carg,'                &   !d13C
     &   ,'S_d13C_QC_dt(n),'        &   ! 1 day avaraged value of d13C_QC
# endif
# if defined CORAL_BORON_ISOTOPE
     &   ,'d11Barg,'   & 
# endif
     &   ,'S_Pg_dt,'                &   !Gross photosynthesis rate (umol cm-2 d-1)
     &   ,'S_R_dt,'                 &   !Respiration rate (umol cm-2 d-1)
     &   ,'S_QC_dt(n),'             &   ! 1 day avaraged value of QC
     &   ,'S_Pn_dt'                     !Net photosynthesis rate (umol cm-2 d-1)

      RETURN

    END SUBROUTINE write_crl_ave_lavel
#endif
      
  END MODULE mod_output

