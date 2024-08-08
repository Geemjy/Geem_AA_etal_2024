
#==============================================
# BEFORE RUNNING 
#==============================================

'''
** THIS IS THE SAME FILE THAT UPLOADED AT https://github.com/Geemjy/Geem_etal_MNRAS_2022**
This is the code for deriving the $q$ and $u$ values of the asteroids taken by the wide field grism spectrograph 2 
(WFGS2; Uehara et al. 2004; Kawakami et al. 2022) on the 2.0-m Nayuta telescope at the Nishi–Harima Astronomical Observatory. 


1. 
 - Input file:  
   Phot_{DATE}_{Object_name}.csv        Photometric result of each images generated by 2.WFGS2_aperture_photometry.py.
 
 - Outout file:
   Pol_{DATE}_{Object_name}.csv         Polarimetric result of each set
     
'''     

#==============================================
# INPUT VALUE
#==============================================
Observatory = {'lon': 134.3356,
               'lat': 35.0253,
               'elevation': 0.449} #NHAO
Target_name = 3200
subpath = 'Directory path where Phot_{DATE}_{Object_name}.csv is saved.' 
filename = 'Phot_{DATE}_{Object_name}.csv'





#==============================================
# IMPORT PACKAGES AND DEFINE THE FUNCTION
#==============================================



import os
import pandas as pd
import numpy as np
from astropy.time import Time
from astroquery.jplhorizons import Horizons


def NHAO_qu(kappa,err_kappa,INST_PA,INSROT):

    k_0 = kappa[0]
    k_45 = kappa[1]
    k_22 = kappa[2]
    k_67 = kappa[3]

    ek_0 = err_kappa[0]
    ek_45 = err_kappa[1]
    ek_22 = err_kappa[2]
    ek_67 = err_kappa[3]

    aQ = np.sqrt(k_0/k_45)
    aU = np.sqrt(k_22/k_67)

    q = (1-aQ)/(1+aQ) #Q/I
    u = (1-aU)/(1+aU) #U/I

    q_ran = aQ/((aQ + 1)**2)  *  np.sqrt(ek_0 + ek_45)
    u_ran = aU/((aU + 1)**2)  *  np.sqrt(ek_22 + ek_67)    


    ###==================== 
    ## Correct Efficiency
    ###==================== 
    eff = 1 #assume
    efferr = 0 #assume

    qq = q/eff
    uu = u/eff

    #random error of corrected q,u
    qq_ran = q_ran/eff
    uu_ran = u_ran/eff

    #the systematic errors
    qq_sys = np.abs(q)*efferr/eff**2
    uu_sys = np.abs(u)*efferr/eff**2


    ###==================== 
    ## Correc Instrumental polarization
    ###====================     
    q_inst = -0.00043
    u_inst =  0.00273
    eq_inst =  0.00012
    eu_inst =  0.00012 #for 2021 November
    
#     q_inst = -0.00042
#     u_inst =  0.00178
#     eq_inst =  0.00016
#     eu_inst =  0.00011 # For 2021 October

    INSROT = np.deg2rad(INSROT)
    qqq = qq - (q_inst*np.cos(2*INSROT) - u_inst*np.sin(2*INSROT))
    uuu = uu - (q_inst*np.sin(2*INSROT) + u_inst*np.cos(2*INSROT))

    #random error of corrected q,u
    qqq_ran = qq_ran
    uuu_ran = uu_ran    

    #the systematic errors    
    qqq_sys = np.sqrt( qq_sys**2 + (eq_inst*np.cos(2*INSROT))**2 + (eu_inst*np.sin(2*INSROT))**2)
    uuu_sys = np.sqrt( uu_sys**2 + (eq_inst*np.sin(2*INSROT))**2 + (eu_inst*np.cos(2*INSROT))**2)

    ###==================== 
    ## Transform_CelestialCoord
    ###====================    
    the = -5.19 - INST_PA
    the_err = 0.15
    
    theta = the
    theta = np.deg2rad(theta)
    the_err = np.deg2rad(the_err)


    qqqq = qqq * np.cos(2*theta) + uuu*np.sin(2*theta)
    uuuu = -qqq * np.sin(2*theta) + uuu*np.cos(2*theta)

    qqqq_ran = np.sqrt( (qqq_ran*np.cos(2*theta))**2 + (uuu_ran*np.sin(2*theta))**2 )
    uuuu_ran = np.sqrt( (qqq_ran*np.sin(2*theta))**2 + (uuu_ran*np.cos(2*theta))**2 )

    qqqq_sys = np.sqrt( (qqq_sys*np.cos(2*theta))**2 +                         (uuu_sys*np.sin(2*theta))**2 +                         (np.pi/180*2*uuuu*the_err)**2 )
    uuuu_sys = np.sqrt( (qqq_sys*np.sin(2*theta))**2 +                         (uuu_sys*np.cos(2*theta))**2 +                         (np.pi/180*2*qqqq*the_err)**2 ) 
    return(qqqq, qqqq_ran, qqqq_sys, uuuu, uuuu_ran, uuuu_sys)

def weight(x,err):
    x = np.array(x)
    err = np.array(err)
    
    w = 1/err**2
    sumW = np.sum(w)
    weight = w/sumW
    
    xav = np.sum(weight*x)
    Err = 1/np.sqrt(sumW)
    
    return(xav,Err)








#==============================================
# DERIVING q AND u
#==============================================

#======================================#
#             Polarimetry              #
#======================================#    
filename = os.path.join(subpath,filename)
Phot = pd.read_csv(filename)

Pol_log = pd.DataFrame({})
order = np.arange(0,len(Phot),8)     

for i in order:
    OBJECT = Phot['Object'][i]
    
    filenn = Phot['filename'].values[i].split('/')[-1].split('.')[0] +' ~ ' + Phot['filename'].values[i+7].split('/')[-1].split('.')[0]
    level = max(Phot['level'].values[i:i+4])
    object_name = Target_name
    Flux_0_e = Phot['Flux'][i]
    Flux_0_o = Phot['Flux'][i+1]
    eFlux_0_e = Phot['eFlux'][i]
    eFlux_0_o = Phot['eFlux'][i+1]
    err_0 = (Flux_0_e/Flux_0_o**2 * eFlux_0_o)**2 + (1/Flux_0_o * eFlux_0_e)**2

    Flux_45_e = Phot['Flux'][i+2]
    Flux_45_o = Phot['Flux'][i+3]
    eFlux_45_e = Phot['eFlux'][i+2]
    eFlux_45_o = Phot['eFlux'][i+3]
    err_45 = (Flux_45_e/Flux_45_o**2 * eFlux_45_o)**2 + (1/Flux_45_o * eFlux_45_e)**2

    Flux_22_e = Phot['Flux'][i+4]
    Flux_22_o = Phot['Flux'][i+5]
    eFlux_22_e = Phot['eFlux'][i+4]
    eFlux_22_o = Phot['eFlux'][i+5]
    err_22 = (Flux_22_e/Flux_22_o**2 * eFlux_22_o)**2 + (1/Flux_22_o * eFlux_22_e)**2

    Flux_67_e = Phot['Flux'][i+6]
    Flux_67_o = Phot['Flux'][i+7]
    eFlux_67_e = Phot['eFlux'][i+6]
    eFlux_67_o = Phot['eFlux'][i+7]
    err_67 = (Flux_67_e/Flux_67_o**2 * eFlux_67_o)**2 + (1/Flux_67_o * eFlux_67_e)**2


    kappa = [Flux_0_e/Flux_0_o, Flux_45_e/Flux_45_o, Flux_22_e/Flux_22_o, Flux_67_e/Flux_67_o]
    ekappa = [err_0, err_45, err_22, err_67]
    Inst_PA = np.mean(Phot['INST-PA'].values[i:i+7])
    INSROT = np.mean(Phot['INSROT'].values[i:i+7])
    
    q, ran_q, sys_q, u, ran_u, sys_u = NHAO_qu(kappa, ekappa,Inst_PA,INSROT)
   
    P = np.sqrt(q**2 + u**2)
    P_ran = np.sqrt( (q*ran_q)**2 + (u*ran_u)**2 )/P
    P_sys = np.sqrt( (q*sys_q)**2 + (u*sys_u)**2 )/P
    P_error = np.sqrt(P_ran**2 + P_sys**2) #Polarization error  

    if P**2 - P_ran**2 < 0:
        Pcor = 0
    else:
        Pcor = np.sqrt(P**2 - P_ran**2)


    theta_pol = np.rad2deg(1/2* np.arctan2(u,q)) 
    ran_PolAng = 1/2 * 180/3.14 * P_ran/Pcor
    sys_PolAng = 1/2 * 180/3.14 * P_sys/Pcor
    etheta = np.sqrt(ran_PolAng**2 + sys_PolAng**2)  

    PA_av = np.mean(Phot['PA'].values[i:i+7])
    PsANG_av = np.mean(Phot['PsANG'].values[i:i+7])
    SNR_av = np.mean(Phot['SNR'].values[i:i+7])
    Aper_av = np.mean(Phot['Aper_pix'].values[i:i+7])

    if PsANG_av + 90 < 180:
        pi = PsANG_av + 90
    else:
        pi = PsANG_av - 90

    Theta_r = theta_pol - pi
    Pr = Pcor * np.cos(2*np.deg2rad(Theta_r))
    Pol_log = Pol_log.append({'filename':filenn,
                              'Object':OBJECT,
                              'DATE':Phot['DATE'].values[0],
                              'set':int(i/8),
                              'Theta':theta_pol,
                              'eTheta':etheta,
                              'P':Pcor,
                              'eP':P_error,
                              'q':q,
                              'u':u,
                              'ran_q':ran_q,
                              'sys_q':sys_q,
                              'ran_u':ran_u,
                              'sys_u':sys_u,
                              'eq':np.sqrt(ran_q**2 + sys_q**2),
                              'eu':np.sqrt(ran_u**2 + sys_u**2),
                              'Pr':Pr,
                              'Theta_r':Theta_r,
                              'PA':PA_av,
                              'PsANG':PsANG_av,
                              'SNR':SNR_av,
                              'Aper_Dia_pix':Aper_av,
                              'level':level},
                              ignore_index = True)
    
q_av, ranq_av = weight(Pol_log['q'],Pol_log['ran_q'])
u_av, ranu_av = weight(Pol_log['u'],Pol_log['ran_u'])
sysq_av = np.mean(Pol_log['sys_q'])
sysu_av = np.mean(Pol_log['sys_u'])
errq_av = (ranq_av**2 + sysq_av**2)**0.5
erru_av = (ranu_av**2 + sysu_av**2)**0.5

P = np.sqrt(q_av**2+u_av**2)
ran_P = np.sqrt((q_av*ranq_av)**2 + (u_av*ranu_av)**2)/P
sys_P = np.sqrt((q_av*sysq_av)**2 + (u_av*sysu_av)**2)/P
eP = np.sqrt(ran_P**2 + sys_P**2)

Theta = 1/2*np.rad2deg(np.arctan2(u_av,q_av))
if P**2 - ran_P**2 < 0:
    Pcor = 0
    ran_PolAng = 51.96
    sys_PolAng = 51.96
    PolAng_error = 51.96
    # Naghizadeh-Khouei & Clarke 1993
else:
    Pcor = np.sqrt(P**2 - ran_P**2)
    ran_PolAng = 1/2 * 180/3.14 * P_ran/Pcor
    sys_PolAng = 1/2 * 180/3.14 * P_sys/Pcor
    PolAng_error = np.sqrt(ran_PolAng**2 + sys_PolAng**2)  

PsANG_av = np.mean(Pol_log['PsANG'])
PA_av = np.mean(Pol_log['PA'])
Aper_av = np.mean(Pol_log['Aper_Dia_pix'])

#Scattering plane
if PsANG_av + 90 < 180:
    pi = PsANG_av + 90
else:
    pi = PsANG_av - 90

theta_r = Theta - pi
Pr = Pcor * np.cos(2*np.deg2rad(theta_r))
SNR_av = np.mean(Pol_log['SNR'])

file_name = 'Weighted_average'
Pol_log = Pol_log.append({'filename':file_name,
                          'Object':OBJECT,
                          'DATE':Phot['DATE'].values[0],
                          'Theta':Theta,
                          'eTheta':PolAng_error,
                          'P':Pcor,
                          'eP':eP,
                          'q':q_av,
                          'u':u_av, 
                          'ran_q':ranq_av,
                          'sys_q':sysq_av,
                          'ran_u':ranu_av,
                          'sys_u':sysq_av,
                          'eq':np.sqrt(ranq_av**2 + sysq_av**2),
                          'eu':np.sqrt(ranu_av**2 + sysq_av**2),
                          'Pr':Pr,
                          'Theta_r':theta_r,
                          'PA':PA_av,
                          'PsANG':PsANG_av,
                          'SNR':SNR_av,
                          'Aper_Dia_pix':Aper_av},
                          ignore_index = True)    


new_index = ['filename','DATE','Object','set','PA','Pr','eP',
             'Theta_r','eTheta','q','u','ran_q','sys_q',
             'ran_u','sys_u','eq','eu','P','Theta','PsANG','SNR','Aper_Dia_pix','level']
Pol_log = Pol_log.reindex(columns = new_index)  
Pol_log = Pol_log.round({'PA':1,'Pr':4,'eP':4,'Theta_r':4,
                         'eTheta':4,'q':4,'u':4,'ran_q':4,'sys_q':4,
                         'ran_u':4,'sys_u':4,'eq':4,'eu':4,'P':4,
                         'Theta':2,'PsANG':2,'SNR':1,'Aper_Dia_pix':1})


filename = os.path.join(path,'Phot_{0}_{1}.csv'.format(Phot['DATE'].values[0].replace('-','_'),OBJECT))
Pol_log.to_csv(filename)    
print(filename + ' is saved.')
