# Implementataion of glioblastoma radiation response model from Leder et al. 2014 Cell

from __future__ import division
import numpy as np

def a(t, t0, parameters):
    '''Proliferation of differentiated cells'''

    r_d = parameters[10]
    r_s = parameters[11]
    a_s = parameters[13]
    nu = parameters[14]

    a = gamma(t0, parameters)*np.exp(-nu*t) \
        + (1 - gamma(t0, parameters))*J1_d(t, r_d, parameters) \
        + a_s*gamma(t0, parameters)*nu/(r_s + nu)*(1/(r_s - r_d)*(J2(t, r_s, parameters) - J2(t, r_d, parameters)) \
        + 1/r_s*(J3(t, r_s, parameters) - J2(t, r_s, parameters)) \
        + 1/r_s*(J4(t, r_s, parameters) - J4(t, 0, parameters))) \
        + a_s*gamma(t0, parameters)*nu/(r_s + nu)*(1/(nu + r_d)*(J2(t, -nu, parameters) - J2(t, r_d, parameters)) \
        + 1/nu*(J3(t, -nu, parameters) - J2(t, -nu, parameters)) \
        + 1/nu*(J4(t, -nu, parameters) - J4(t, 0, parameters)))
    
    return a


def b(t, parameters):
    '''Proliferation of cells that differentiated from stem-like to
    differentiated cells since previous radiation dose'''
    
    r_d = parameters[10]
    r_s = parameters[11]
    R = parameters[12]
    a_s = parameters[13]
    
    b = a_s/(R*(r_s - r_d))*(J2(t, r_s, parameters) - J2(t, r_d, parameters)) \
        + a_s/(R*r_s)*(J3(t, r_s, parameters) - J2(t, r_s, parameters)) \
        + a_s/(R*r_s)*(J4(t, r_s, parameters) - J4(t, 0, parameters))
    
    return b

                
def c(t, t0, parameters):
    '''Proliferation of cells that dedifferentiated from differentiated to
    stem-like cells since previous radiation dose'''
    
    r_s = parameters[11]
    R = parameters[12]
    nu = parameters[14]
    
    c = gamma(t0, parameters)*R/(nu + r_s)*(nu*J1_s(t, r_s, parameters) \
        + r_s*J1_s(t, -nu, parameters)) \
        - gamma(t0, parameters)*R*np.exp(-nu*t)

    return c


def d(t, parameters):
    '''Proliferation of stem-like cells'''

    r_s = parameters[11]

    d = J1_s(t, r_s, parameters)

    return d
   
   
def gamma(t0, parameters):
    '''Time-dependent dedifferentiation from stem-like to differentiated cells'''
    
    gamma0 = parameters[15]
    mu = parameters[16]
    sigma2 = parameters[17]
    
    if (t0 == 'firstDose'):
        
        gamma = 0
        
    else:
    
        gamma = gamma0*np.exp(-(t0 - mu)**2/sigma2)
    
    
    return gamma
             
             
def J1_d(t, z, parameters):

    L_d = parameters[6]
    lambda_d = parameters[8]

    if (z + lambda_d != 0):
        
        if (t == L_d):
            
            J1_d = 0
            
        # If time is less than minimum quiescence time cells do not proliferate
        elif (t < L_d):
            
            J1_d = 1
            
        elif (t > L_d):
            
            J1_d = lambda_d/(z + lambda_d)*np.exp(z*(t - L_d)) \
                + z/(z + lambda_d)*np.exp(-lambda_d*(t - L_d))
        
    elif (z + lambda_d == 0):
        
        if (t == L_d):
            
            J1_d = 0
           
        # If time is less than minimum quiescence time cells do not proliferate 
        elif (t < L_d):
            
            J1_d = 1
            
        elif (t > L_d):
            
            J1_d = lambda_d*np.exp(z*t + lambda_d*L_d)*(t - L_d) \
                + np.exp(-lambda_d*(t - L_d))

    return J1_d
    
      
def J1_s(t, z, parameters):

    L_s = parameters[7]
    lambda_s = parameters[9]
    
    if (t == L_s):
        
        J1_s = 0
    
    elif (t < L_s):
        
        J1_s = 1
        
    elif (t > L_s):
        
        J1_s = lambda_s/(z + lambda_s)*np.exp(z*(t - L_s)) \
            + z/(z + lambda_s)*np.exp(-lambda_s*(t - L_s))

    return J1_s
    
    
def J2(t, z, parameters):

    eta_d = parameters[4]
    M_d = parameters[5]
    L_s = parameters[7]
    lambda_s = parameters[9]    

    if (lambda_s + z != 0):
        
        if (t <= L_s + M_d):
            
            J2 = 0
            
        elif (t > L_s + M_d):
            
            J2 = lambda_s*eta_d/(lambda_s + z)*(
                np.exp(z*(t - L_s))*np.exp(eta_d*M_d)*(
                np.exp(-M_d*(z + eta_d)) - np.exp(-(t - L_s)*(z + eta_d)))/(eta_d + z) \
                - np.exp(-lambda_s*(t - L_s))*np.exp(eta_d*M_d)*(
                np.exp((t - L_s)*(lambda_s - eta_d)) - np.exp(M_d*(lambda_s - eta_d)))/(lambda_s - eta_d))

    elif (lambda_s + z == 0):

        if (t <= L_s + M_d):

            J2 = 0

        elif (t > L_s + M_d):

            J2 = lambda_s*eta_d*np.exp(-lambda_s*t)*np.exp(lambda_s*L_s + eta_d*M_d)/(z + eta_d)*(
                (t - L_s)*(np.exp(-M_d*(z + eta_d)) - np.exp(-(t - L_s)*(z + eta_d))) \
                - M_d*np.exp(-M_d*(z + eta_d)) \
                + (t - L_s)*np.exp*(-(t - L_s)*(z + eta_d)) \
                - (np.exp(-M_d*(z + eta_d)) - np.exp(-(t - L_s)*(z + eta_d)))/(z + eta_d))
    
    return J2
    
    
def J3(t, z, parameters):

    eta_d = parameters[4]
    M_d = parameters[5]
    L_s = parameters[7]
    lambda_s = parameters[9]

    if (t <= L_s + M_d):
        
        J3 = 0
        
    elif (t > L_s + M_d):
        
        J3 = lambda_s*np.exp(z*t + lambda_s*L_s)/(lambda_s + z)*(
            np.exp(-L_s*(z + lambda_s)) - np.exp(-(t - M_d)*(z + lambda_s))) \
            - lambda_s*np.exp(z*t + lambda_s*L_s)*np.exp(-eta_d*(t - M_d))/(z + lambda_s - eta_d)*(
            np.exp(-L_s*(z + lambda_s - eta_d)) - np.exp(-(t - M_d)*(z + lambda_s - eta_d)))
            
    if (np.isnan(J3) == True):
        
        J3 = 0
    
    return J3
  
      
def J4(t, z, parameters):

    eta_d = parameters[4]
    M_d = parameters[5]
    L_s = parameters[7]
    lambda_s = parameters[9]

    if (t == L_s + M_d):
         
         J4 = 0
    
    elif (t < L_s + M_d):
        
        if (t <= L_s):
            
            J4 = 0
            
        elif (t > L_s):
            
            J4 = lambda_s*np.exp(z*t)*np.exp(lambda_s*L_s)/(lambda_s + z)*(
                np.exp(-L_s*(lambda_s + z)) - np.exp(-t*(lambda_s + z)))
        
    elif (t > L_s + M_d):
        
        J4 = lambda_s*np.exp(z*t)*np.exp(lambda_s*L_s + eta_d*M_d)*(
            np.exp(-eta_d*t)/(z + lambda_s - eta_d)*(
            np.exp(-L_s*(z + lambda_s - eta_d)) - np.exp(-(t - M_d)*(z + lambda_s - eta_d))) \
            + np.exp(-eta_d*M_d)/(lambda_s + z)*(
            np.exp(-(lambda_s + z)*(t - M_d)) - np.exp(-t*(lambda_s + z))))
            
    if (np.isnan(J4) == True):
        
        J4 = 0
    
    return J4
    

def cell_fractions(schedule, parameters, evaluationTime):
    
    alpha_d = parameters[0]
    beta_d = parameters[1]
    alpha_s = parameters[2]
    beta_s = parameters[3]
    R = parameters[12]
    
    dose = np.append(schedule[0], np.zeros(int((evaluationTime - schedule[1][-1])/24.0)))
    time = np.append(schedule[1], np.arange(schedule[1][-1] + 24, evaluationTime, 24))
    
    F_d = np.ones(np.size(time))
    F_s = np.ones(np.size(time))
    relativeTumorSize = np.ones(np.size(time))
    
    for i in range(0, np.size(time) - 1):

        # During treatment
        if (dose[i] > 0):
                
            t = time[i + 1] - time[i]
            
            if (i >= 1):
            
                t0 = time[i] - time[i - 1]
         
            elif (i == 0):
            
                t0 = 'firstDose'
                
            F_d[i + 1] = F_d[i]*np.exp(-alpha_d*dose[i] - beta_d*dose[i]**2)*a(t, t0, parameters) \
                + F_s[i]*np.exp(-alpha_s*dose[i] - beta_s*dose[i]**2)*b(t, parameters)
        
            F_s[i + 1] = F_d[i]*np.exp(-alpha_d*dose[i] - beta_d*dose[i]**2)*c(t, t0, parameters) \
                + F_s[i]*np.exp(-alpha_s*dose[i] - beta_s*dose[i]**2)*d(t, parameters)
                  
        # After treatment
        elif (dose[i] == 0):
                
            t = time[i + 1] - schedule[1][-1]
            t0 = time[i] - schedule[1][-1]
            
            F_d[i + 1] = F_d[np.size(schedule[1])]*np.exp(-alpha_d*dose[i] - beta_d*dose[i]**2)*a(t, t0, parameters) \
                + F_s[np.size(schedule[1])]*np.exp(-alpha_s*dose[i] - beta_s*dose[i]**2)*b(t, parameters)
        
            F_s[i + 1] = F_d[np.size(schedule[1])]*np.exp(-alpha_d*dose[i] - beta_d*dose[i]**2)*c(t, t0, parameters) \
                + F_s[np.size(schedule[1])]*np.exp(-alpha_s*dose[i] - beta_s*dose[i]**2)*d(t, parameters)
                        
    relativeTumorSize = (1/R*F_s + F_d)/(1/R + 1)
    fractionalVolumeChange = (relativeTumorSize - relativeTumorSize[0])/relativeTumorSize[0]
    
    return (time, F_d, F_s, fractionalVolumeChange)

