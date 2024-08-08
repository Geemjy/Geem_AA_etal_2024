class MSI():
    import numpy as np
    import pandas as pd
    
    
    RESULT = [] #For __str__
    
    
    def __init__(self, PsANG, Filter, calib_file):
        import pandas as pd
        
        self.RESULT = []
        self.PSANG = PsANG
        self.Filter = Filter
        self.calib_file =calib_file
    
    def Get_qu(self,I,S):
        import numpy as np
        import pandas as pd
    
        '''
        Estimate q,u value. 
        Every value is in [electron] unit
        I = [I0,I45,I22,I67]
        S = [S0,S45,S22,S67]
        Input vaue is I_ex/I_o at each RET-ANG2.
        i.e., 
        I0  = I_ex/I_o at RET-ANG2 = 0deg
        I45 = I_ex/I_o at RET-ANG2 = 45deg
        I22 = I_ex/I_o at RET-ANG2 = 22.5deg
        I67 = I_ex/I_o at RET-ANG2 = 67.5deg
        
        
        S0 = (Sigma_{I_ex,0}/I_ex,0)**2 + (Sigma_{I_o,0}/I_o,0)**2  
        
        Return q,u,q_ran,q_sys,u_ran,u_sys
        '''
        I0 = I[0]
        I45 = I[1]
        I22 = I[2]
        I67 = I[3]
        
        S0 = S[0]
        S45 = S[1]
        S22 = S[2]
        S67 = S[3]
        
            
        Rq = np.sqrt(I0/I45)
        Ru = np.sqrt(I22/I67)
                
        q = (Rq - 1)/(Rq + 1)
        u = (Ru - 1)/(Ru + 1)
        
        
        q_ran = Rq/((Rq + 1)**2)  *  np.sqrt(S0 + S45)
        u_ran = Ru/((Ru + 1)**2)  *  np.sqrt(S22 + S67)
        
        q_sys = 0
        u_sys = 0

        q_err = np.sqrt(q_ran**2 + q_sys**2)
        u_err = np.sqrt(u_ran**2 + u_sys**2)
        
        self.RESULT = [q,u,q_ran,q_sys,u_ran,u_sys,self.PSANG,'Get_qu']

        return q,u,q_ran,q_sys,u_ran,u_sys
    
    
    def Correct_Efficiency(self,q,u,q_ran,q_sys,u_ran,u_sys,verbose=False):
        import numpy as np
        import pandas as pd
        #q,u obtained fr===om the instrument system by the polarization efficiency, p
        df = pd.read_csv(self.calib_file)
        if self.Filter == 'Rc':
            eff = df[df['param']=='eff_rc']['value'].values[0]
            efferr = df[df['param']=='err_eff_rc']['value'].values[0]
        elif self.Filter == 'V':
            eff = df[df['param']=='eff_v']['value'].values[0]
            efferr = df[df['param']=='err_eff_v']['value'].values[0]
        elif self.Filter == 'B':
            eff = df[df['param']=='eff_b']['value'].values[0]
            efferr = df[df['param']=='err_eff_b']['value'].values[0]
           # print('B is choosen.')      
        elif self.Filter == 'Ic':
            eff = df[df['param']=='eff_ic']['value'].values[0]
            efferr = df[df['param']=='err_eff_ic']['value'].values[0]
           # print('Ic is choosen.')      

        
        
        qq = q/eff
        uu = u/eff
        
        #random error of corrected q,u
        qq_ran = q_ran/eff
        uu_ran = u_ran/eff
        
        #the systematic errors
        qq_sys = np.abs(q)*efferr/eff
        uu_sys = np.abs(u)*efferr/eff
        
        self.RESULT = [qq,uu,qq_ran,qq_sys,uu_ran,uu_sys,self.PSANG,'Effici']
        
        if verbose==False:
            return qq,uu,qq_ran,qq_sys,uu_ran,uu_sys
        elif verbose==True:
            return qq,uu,qq_ran,qq_sys,uu_ran,uu_sys, eff, efferr
    


    def Correct_Inst_Pol(self,qq,uu,qq_ran,qq_sys,uu_ran,uu_sys,
                         INR_STR,INR_END,tilted_ang,verbose=False):
        import numpy as np
        import pandas as pd
        '''
        INR_STR = [INR_STR(0),INR_STR(45),INR_STR(22),INR_STR(67)]
        INR_END = [INR_END(0),INR_END(45),INR_END(22),INR_END(67)]
        Return qqq,uuu,qqq_ran,qqq_sys,uuu_ran,uuu_sys
        '''
        STR0,STR45,STR22,STR67 = INR_STR[0],INR_STR[1],INR_STR[2],INR_STR[3]
        END0,END45,END22,END67 = INR_END[0],INR_END[1],INR_END[2],INR_END[3]
        rq = (STR0 + END0 + STR45 + END45)/4.+tilted_ang
        ru = (STR22 + END22 + STR67 + END67)/4.+tilted_ang
        

        
        
        if rq < 0:
                rq = rq+180
        if rq >180:
            rq = rq-180
 
         #averaged value of frame-reader value
        rq = np.deg2rad(rq)
         # Instrument star, end value ( 0,0,45,45 ) in rad
        ru = np.deg2rad(ru)
         # Instrument star, end value ( 22.5,22.5,67.5,67.5) in rad  
            
        df = pd.read_csv(self.calib_file)    
        if self.Filter == 'Rc':
            q_inst = df[df['param']=='q_inst_rc']['value'].values[0]
            u_inst =df[df['param']=='u_inst_rc']['value'].values[0]
            err_qin = df[df['param']=='err_q_inst_rc']['value'].values[0]
            err_uin = df[df['param']=='err_u_inst_rc']['value'].values[0]
            
        elif self.Filter == 'V':
            q_inst = df[df['param']=='q_inst_v']['value'].values[0]
            u_inst =df[df['param']=='u_inst_v']['value'].values[0]
            err_qin = df[df['param']=='err_q_inst_v']['value'].values[0]
            err_uin = df[df['param']=='err_u_inst_v']['value'].values[0]
            
        elif self.Filter == 'B':
            q_inst = df[df['param']=='q_inst_b']['value'].values[0]
            u_inst =df[df['param']=='u_inst_b']['value'].values[0]
            err_qin = df[df['param']=='err_q_inst_b']['value'].values[0]
            err_uin = df[df['param']=='err_u_inst_b']['value'].values[0]
            
        elif self.Filter == 'Ic':
            q_inst = df[df['param']=='q_inst_ic']['value'].values[0]
            u_inst =df[df['param']=='u_inst_ic']['value'].values[0]
            err_qin = df[df['param']=='err_q_inst_ic']['value'].values[0]
            err_uin = df[df['param']=='err_u_inst_ic']['value'].values[0]
            
        qqq = qq - ((q_inst * np.cos(2*rq)) - (u_inst * np.sin(2*rq)))
        uuu = uu - ((q_inst * np.sin(2*ru)) + (u_inst * np.cos(2*ru)))
        
        qqq_sys = np.sqrt( qq_sys**2 + (err_qin * np.cos(2*rq))**2 +\
                          (err_uin*np.sin(2*rq))**2 )
        uuu_sys = np.sqrt( uu_sys**2 + (err_qin * np.sin(2*ru))**2 +\
                          (err_uin*np.cos(2*ru))**2 )
        
        qqq_ran = qq_ran
        uuu_ran = uu_ran
        
        self.RESULT = [qqq,uuu,qqq_ran,qqq_sys,uuu_ran,uuu_sys,self.PSANG,'Inst']

    
        if verbose==False:
            return qqq,uuu,qqq_ran,qqq_sys,uuu_ran,uuu_sys
        elif verbose==True:
            return (qqq,uuu,qqq_ran,qqq_sys,uuu_ran,uuu_sys
                    ,q_inst, err_qin, u_inst, err_uin)
    
    
    def Transform_CelestialCoord(self,qqq,uuu,
                                 qqq_ran,qqq_sys,
                                 uuu_ran,uuu_sys,PosANG,tilted_ang,
                                 verbose=False):
        import numpy as np
        import pandas as pd
        df = pd.read_csv(self.calib_file)    
        if self.Filter == 'Rc':
            the = df[df['param']=='the_rc']['value'].values[0]
            the_err =df[df['param']=='err_the_rc']['value'].values[0]
        elif self.Filter == 'V':
            the = df[df['param']=='the_v']['value'].values[0]
            the_err =df[df['param']=='err_the_v']['value'].values[0]
        elif self.Filter == 'B':
            the = df[df['param']=='the_b']['value'].values[0]
            the_err =df[df['param']=='err_the_b']['value'].values[0]
        elif self.Filter == 'Ic':
            the = df[df['param']=='the_ic']['value'].values[0]
            the_err =df[df['param']=='err_the_ic']['value'].values[0]

         
        PosANG = PosANG
        theta = np.deg2rad(the) - np.deg2rad(PosANG) + np.deg2rad(tilted_ang)
        qqqq = qqq * np.cos(2*theta) + uuu*np.sin(2*theta)
        uuuu = -qqq * np.sin(2*theta) + uuu*np.cos(2*theta)
        
        qqqq_ran = np.sqrt( (qqq_ran*np.cos(2*theta))**2 + (uuu_ran*np.sin(2*theta))**2 )
        uuuu_ran = np.sqrt( (qqq_ran*np.sin(2*theta))**2 + (uuu_ran*np.cos(2*theta))**2 )
        
        qqqq_sys = np.sqrt( (qqq_sys*np.cos(2*theta))**2 + \
                            (uuu_sys*np.sin(2*theta))**2 + \
                            (np.pi/180*2*uuuu*the_err)**2 )
        uuuu_sys = np.sqrt( (qqq_sys*np.sin(2*theta))**2 + \
                            (uuu_sys*np.cos(2*theta))**2 + \
                            (np.pi/180*2*qqqq*the_err)**2 )  
        
        self.RESULT = [qqqq,uuuu,qqqq_ran,qqqq_sys,uuuu_ran,uuuu_sys,self.PSANG,\
                       'Celestial']
        
        if verbose==False:
            return qqqq,uuuu,qqqq_ran,qqqq_sys,uuuu_ran,uuuu_sys
        elif verbose==True:
            return (qqqq,uuuu,qqqq_ran,qqqq_sys,uuuu_ran,uuuu_sys, 
                    the, the_err)
    
    
    
    def Scattering_plane(self,q,u,q_ran,q_sys,u_ran,u_sys):
        import numpy as np
        '''
        return P_r, P_error, theta_l, PolAng_error    
        '''
        
        P = np.sqrt(q**2 + u**2)
        P_ran = np.sqrt( (q*q_ran)**2 + (u*u_ran)**2 )/P
        P_sys = np.sqrt( (q*q_sys)**2 + (u*u_sys)**2 )/P
        theta_pol = np.rad2deg(1/2* np.arctan2(u,q))

        PsANG = self.PSANG
        
        if PsANG <= 90:
            pi = PsANG + 90
        elif PsANG > 90:
            pi = PsANG - 90
            
            
        theta_l = theta_pol - pi
        
        
        def set_angle(theta_l):
            if -45 < theta_l < 45:
                return theta_l 
            elif 45 < theta_l < 135:
                return theta_l
            #Or =========================
            elif theta_l < -45:
                theta_l = theta_l + 180
                return set_angle(theta_l)
            elif theta_l > 135:
                theta_l = theta_l - 180
                return set_angle(theta_l)
        
        if P**2 > P_ran**2:
            print('Random error bias correction is done.')
            P_cor = np.sqrt(P**2 - P_ran**2)

        elif P**2 < P_ran**2 :
            print('Due to P < randome error, random error bias correction is NOT done.')
            P_cor = P
            
            
        P_error = np.sqrt(P_ran**2 + P_sys**2) #Polarization error
        theta_l = set_angle(theta_l) # Position angle

        ran_PolAng = 1/2 * 180/3.14 * P_ran/P_cor
        sys_PolAng = 1/2 * 180/3.14 * P_sys/P_cor
        PolAng_error = np.sqrt(ran_PolAng**2 + sys_PolAng**2)
        
        P_r = P_cor * np.cos(2*np.deg2rad(theta_l))

        
        return P_r, P_error, theta_l, PolAng_error    
    
    def __str__(self):
        import numpy as np
        
        q = self.RESULT[0]
        u = self.RESULT[1]
        q_ran = self.RESULT[2]
        q_sys = self.RESULT[3]
        u_ran = self.RESULT[4]
        u_sys = self.RESULT[5]
        PsANG = self.RESULT[6]
        mode = self.RESULT[7]
        
        if mode == 'Get_qu':
            text = 'BEFORE ANY CORRECTION ==================================='
        elif mode == 'Effici':
            text = 'AFTER EFFICIENCY CORRECTION=============================='
        elif mode == 'Inst':
            text = 'AFTER INSTRUMENT POLATIZATION CORRECTION================='
        elif mode == 'Celestial':
            text = 'AFTER TRANSFORMING TO CELESTIAL COORDINATE==============='
            
            
        P = np.sqrt(q**2 + u**2)
        P_ran = np.sqrt( (q*q_ran)**2 + (u*u_ran)**2 )/P
        P_sys = np.sqrt( (q*q_sys)**2 + (u*u_sys)**2 )/P   
        theta_pol = np.rad2deg(1/2* np.arctan2(u,q))


        
        if PsANG <= 90:
            pi = PsANG + 90
        elif PsANG > 90:
            pi = PsANG - 90 
            
        theta_l = theta_pol - pi
        
        def set_angle(theta_l):
            if -45 < theta_l < 45:
                return theta_l 
            elif 45 < theta_l < 135:
                return theta_l
            #Or =========================
            elif theta_l < -45:
                theta_l = theta_l + 180
                return set_angle(theta_l)
            elif theta_l > 135:
                theta_l = theta_l - 180
                return set_angle(theta_l)
        
        if P**2 > P_ran**2:
            text2 = 'Random error bias correction is done.'
            P_cor = np.sqrt(P**2 - P_ran**2)

        elif P**2 < P_ran**2 :
            text2 = 'Due to P < randome error, random error bias correction is NOT done.'
            P_cor = P
            
            
        P_error = np.sqrt(P_ran**2 + P_sys**2) #Polarization error
        theta_l = set_angle(theta_l) # Position angle
        
        ran_PolAng = 1/2 * 180/3.14 * P_ran/P_cor
        sys_PolAng = 1/2 * 180/3.14 * P_sys/P_cor
        PolAng_error = np.sqrt(ran_PolAng**2 + sys_PolAng**2)
        
        P_r = P_cor * np.cos(2*np.deg2rad(theta_l))
        
        printed_text = text + '\n'\
        + '  Pr = {0:.2f} +/- {1:.2f}'.format(P_r*100,P_error*100) \
        + '  (' + text2 + ')\n'\
        + '  Theta_r = {0:.2f} +/- {1:.2f}'.format(theta_l,PolAng_error)

        return printed_text
