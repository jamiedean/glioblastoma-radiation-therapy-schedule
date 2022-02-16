import os
import numpy as np
import matplotlib.pyplot as plt
import gbm_radiation_response_model as model
import warnings

# Suppress warnings due to overflow. These do not affect results.
warnings.filterwarnings('ignore')

os.chdir('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule/code')

def compareSchedules(schedules, schedule_labels, schedule_colors, parameters, evaluationTime, filename):
    
    plt.figure()
    
    for n in range(0, np.size(schedule_labels)):
        
        schedule = schedules[n]
        time, F_d, F_s, fractionalVolumeChange = model.cell_fractions(schedule, parameters, evaluationTime)

        plt.plot(time/24.0, fractionalVolumeChange, label = schedule_labels[n],
                 linewidth = 4, color = schedule_colors[n])
        plt.xlabel('Time since start of treatment (days)', fontsize = 16)
        plt.ylabel('Normalized volume change', fontsize = 16)
        plt.xlim(0, 100)
        plt.ylim(-1.25, 2)
        plt.yticks(np.arange(-1, 3, 1))
        plt.legend(loc = 'upper left', frameon = False, fontsize = 14)
        plt.tick_params(direction = 'in', labelsize = 14)
    
    plt.savefig('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule//figures/' + filename + '_volume.pdf')
    
    plt.figure()
    
    for n in range(0, np.size(schedule_labels)):
        
        schedule = schedules[n]
        time, F_d, F_s, fractionalVolumeChange = model.cell_fractions(schedule, parameters, evaluationTime)

        plt.plot(time/24.0, F_s, label = schedule_labels[n],
                 linewidth = 4, color = schedule_colors[n])
        plt.xlabel('Time since start of treatment (days)', fontsize = 16)
        plt.ylabel('Normalized stem-like cell number', fontsize = 16)
        plt.xlim(0, 100)
        plt.legend(loc = 'upper left', frameon = False, fontsize = 14)
        plt.tick_params(direction = 'in', labelsize = 14)
    
    plt.savefig('/Users/jamiedean/Documents/Papers/In Revision/Novel Glioblastoma Radiation Schedule/Github Repository/glioblastoma-radiation-therapy-schedule//figures/' + filename + '_stem_cell_fraction.pdf')
    
rtog1205 = \
    np.array([[3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5],
              [0, 24, 48, 72, 96, 168, 192, 216, 240, 264]])

tid_10d = np.array([np.repeat(1.72, 30),
    [0, 3.25, 6.5, 24, 27.25, 30.5, 48, 51.25, 54.5, 72, 75.25, 78.5, 96, 99.25, 102.5, 168, 171.25, 174.5,
     192, 195.25, 198.5, 216, 219.25, 222.5, 240, 243.25, 246.5, 264, 267.25, 270.5]])
    
schedules_tid_10d = [rtog1205, tid_10d]
schedule_labels_tid_10d = ['3.5 Gy x 10 (QD)', '1.72 Gy x 30 (TID)']
schedule_colors_tid_10d = ['#f89821', '#02679a']
    
tid_5d = \
    np.array([[4.52, 4.52, 4.52, 4.52, 4.52,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 48, 72, 96,
     168, 171.25, 174.5, 192, 195.25, 198.5, 216, 219.25, 222.5, 240, 243.25, 246.5, 264, 267.25, 270.5]])
    
tid_3d = \
    np.array([[3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 48, 72, 96, 168, 192,
     216, 219.25, 222.5, 240, 243.25, 246.5, 264, 267.25, 270.5]])
    
schedules = [rtog1205, tid_5d, tid_3d]
schedule_labels = ['3.5 Gy x 10 (QD)', '4.52 Gy x 5 (QD) + \n 1.0 Gy x 15 (TID)',
                   '3.96 Gy x 7 (QD) + \n 1.0 Gy x 9 (TID)']
schedule_colors = ['#f89821', '#02679a', '#3ac6f3']

tid_3d_2_25h = np.array([[
        3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 48, 72, 96, 168, 192,
     216, 218.25, 220.5, 240, 242.25, 244.5, 264, 266.25, 268.5]])
    
tid_3d_4_25h = np.array([[
        3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 48, 72, 96, 168, 192,
     216, 220.25, 224.5, 240, 244.25, 248.5, 264, 268.25, 272.5]])

tid_3d_1_75h = np.array([[
        3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 48, 72, 96, 168, 192,
     216, 217.75, 219.5, 240, 241.75, 243.5, 264, 265.75, 267.5]])
    
tid_3d_4_75h = np.array([[
        3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 48, 72, 96, 168, 192,
     216, 220.75, 225.5, 240, 244.75, 249.5, 264, 268.75, 273.5]])

tid_interval = [rtog1205, tid_3d, tid_3d_2_25h, tid_3d_4_25h, tid_3d_1_75h, tid_3d_4_75h]
tid_interval_labels = ['3.5 Gy x 10 (QD)', '3.25 h', '2.25 h', '4.25 h', '1.75 h', '4.75 h']
tid_interval_colors = ['#f89821', '#02679a', '#3ac6f3', '#4d4d4f', 'pink', 'darkgreen']

qd_4h = \
    np.array([[3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 28, 48, 76, 96, 172, 192,
     216, 219.25, 222.5, 240, 243.25, 246.5, 264, 267.25, 270.5]])
    
qd_8h = \
    np.array([[3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 32, 48, 80, 96, 176, 192,
     216, 219.25, 222.5, 240, 243.25, 246.5, 264, 267.25, 270.5]])

qd_interval = [rtog1205, tid_3d, qd_4h, qd_8h]
qd_interval_labels = ['3.5 Gy x 10 (QD)', '0 h', '4 h', '8 h']
qd_interval_colors = ['#f89821', '#02679a', '#3ac6f3', '#4d4d4f']

tid_tue = \
    np.array([[3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 48, 72, 144, 168, 192,
     216, 219.25, 222.5, 240, 243.25, 246.5, 312, 315.25, 318.5]])
    
tid_wed = \
    np.array([[3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 48, 120, 144, 168, 192,
     216, 219.25, 222.5, 288, 291.25, 294.5, 312, 315.25, 318.5]])
    
tid_thu = \
    np.array([[3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 24, 96, 120, 144, 168, 192,
     264, 267.25, 270.5, 288, 291.25, 294.5, 312, 315.25, 318.5]])
    
tid_fri = \
    np.array([[3.96, 3.96, 3.96, 3.96, 3.96, 3.96, 3.96,
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0, 72, 96, 120, 144, 168, 240,
     288, 291.25, 294.5, 312, 315.25, 318.5, 336, 339.25, 342.5]])

start_day = [rtog1205, tid_3d, tid_tue, tid_wed, tid_thu, tid_fri]
start_day_labels = ['3.5 Gy x 10 (QD)', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri']
start_day_colors = ['#f89821', '#02679a', '#3ac6f3', '#4d4d4f', 'pink', 'darkgreen']

# Parameter values from Leder et al. Cell 2014 table 1 (final iteration)
alpha_d = 0.0987
beta_d = 1.14E-7
rho = 0.1
alpha_s = rho*alpha_d
beta_s = rho*beta_d
eta_d = 0.054
M_d = 366.3
L_d = 193.32
L_s = 477.02
lambda_d = 0.1
lambda_s = 0.0328
r_d = 0.0038
r_s = 0.0008
R = 20.0
a_s = 0.0019
nu = 0.45
gamma0 = 0.4
mu = 3.25
sigma2 = 1.46
parameters = np.array([alpha_d, beta_d, alpha_s, beta_s, eta_d, M_d, L_d, L_s,
    lambda_d, lambda_s, r_d, r_s, R, a_s, nu, gamma0, mu, sigma2])

evaluationTime = 100*24

compareSchedules(schedules_tid_10d, schedule_labels_tid_10d, schedule_colors_tid_10d, parameters, evaluationTime,
                 'compare_schedules_tid_10d')
compareSchedules(schedules, schedule_labels, schedule_colors, parameters, evaluationTime, 
                 'compare_schedules_hypo_tid')
compareSchedules(tid_interval, tid_interval_labels, tid_interval_colors, parameters, evaluationTime, 
                 'compare_schedules_tid_interval')
compareSchedules(qd_interval, qd_interval_labels, qd_interval_colors, parameters, evaluationTime, 
                 'compare_schedules_qd_interval')
compareSchedules(start_day, start_day_labels, start_day_colors, parameters, evaluationTime, 
                 'compare_schedules_start_day')
