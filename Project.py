#Project

import numpy as np
import matplotlib.pyplot as plt
import f2py_projectmodule

#radius of layers
R_base = 0.0047 # units = m, layer 1 (base)
R_upper = 0.0020 # units = m, layer 2, 3, and 4

#thickness of layers - layer 1, 2, 3, 4
Z = [0.0015, 5e-6, 5e-7, 3e-6] # units = m

#Thermal Conductivity - layer 1, 2, 3, 4
k = [35.0, 73.0, 21.0, 5.3] # units = W/ m*K

#Specific Heat - layer 1, 2, 3, 4
Cp = [880.0, 125.6, 544.0, 250.0] # units = J / Kg*K

#Density - layer 1, 2, 3, 4
Rho = [4000.0, 21400.0, 4500.0, 11000.0] # units Kg / m^3

#Initial uniform temperature - units = K
initial_temp = 300.0

#Uniform temperature boundary condition - units = K
#boundary_temp = 0.0

#Peak heat generation - units = W/m^3
heat_gen_peak = 300000000000.0

#Time step
dt = 0.2 #Units = seconds

#Total Time
tMax = 60.0 #units = seconds


dr = 0.0001 #units = m
#dz - layer 1, 2, 3, 4
dz = [7.5e-5, 2.5e-7, 2.5e-8, 1.5e-7] #units = m

Nr_init = int(R_base / dr)
remainder = R_upper % (R_base/Nr_init)
while (remainder > 0.00000001 and abs(remainder - (R_base/Nr_init)) > 0.00000001):
    Nr_init = Nr_init + 1

dr = R_base / Nr_init
Nr_upper = int(Nr_init * (R_upper / R_base))
Nr = [Nr_init, Nr_upper, Nr_upper, Nr_upper]

#time step = 1 s
time_step = 1.0 #seconds
Nt = int(time_step / dt)
#dt = tMax / Nt

Nz = [0, 0, 0, 0]

for i in range(4):
    Nz[i] = int(Z[i] / dz[i])
    dz[i] = Z[i] / Nz[i]

total = (Nr[0]*Nz[0]) + (Nr[1]*Nz[1]) + (Nr[2]*Nz[2]) + (Nr[3]*Nz[3]) 

def initialize2():
    
    temp = np.ones([total]) * initial_temp  
    np.savetxt("temp.txt", temp, delimiter = ' ')
   
t = 0
def process():                   
    global t
    while t < tMax:
        f2py_projectmodule.temp(dr, dz, dt, k, Rho, Cp, heat_gen_peak, Nr, Nz, Nt)
        t = t + (dt * Nt)
        

def plot_middle(a):
    f = open("temp.txt")

    temp_sol = f.read()

    temp_sol = temp_sol.split()

    temp_sol = map(float, temp_sol)
    
    temp_layer1 = np.ones([Nr[0], Nz[0]])

    for c in range(Nr[0]*Nz[0]):
        j = c / Nr[0]
        i = c % Nr[0]
        temp_layer1[i][j] = temp_sol[c]
   
    total_1 = Nr[0]*Nz[0]

    temp_layer2 = np.ones([Nr[1], Nz[1]])    

    for c in range(Nr[1]*Nz[1]):
        j = c / Nr[1]
        i = c % Nr[1]
        temp_layer2[i][j] = temp_sol[c + total_1]

    total_2 = Nr[1]*Nz[1]

    temp_layer3 = np.ones([Nr[2], Nz[2]])    

    for c in range(Nr[2]*Nz[2]):
        j = c / Nr[2]
        i = c % Nr[2]
        temp_layer3[i][j] = temp_sol[c + total_1 + total_2]

    total_3 = Nr[2]*Nz[2]

    temp_layer4 = np.ones([Nr[3], Nz[3]])    

    for c in range(Nr[3]*Nz[3]):
        j = c / Nr[3]
        i = c % Nr[3]
        temp_layer4[i][j] = temp_sol[c + total_1 + total_2 + total_3]
    
    middle1 = np.ones([Nz[0]])
    middle2 = np.ones([Nz[1]])
    middle3 = np.ones([Nz[2]])
    middle4 = np.ones([Nz[3]])
    
    for c in range (Nz[0]):
        middle1[c] = temp_layer1[a][c]

    for c in range (Nz[1]):
        middle2[c] = temp_layer2[a][c]

    for c in range (Nz[2]):
        middle3[c] = temp_layer3[a][c]

    for c in range (Nz[3]):
        middle4[c] = temp_layer3[a][c]

    y = np.append(middle1, [middle2, middle3])
    y = np.append(y, [middle4])

    return y
    

    
def generate_plot():
    
    title = 'Temperature vs. Vertical Position'   
    plt.xlabel('Finite Difference Element in Vertical Direction')
    plt.ylabel('Temperature (K)')
    plt.title(title)
    
    x = np.linspace(1, 80, 80)
    plt.figure(1)
   
    y1 = plot_middle(0)
    plt.plot(x, y1, label = "r = 0.0 cm")
    y1 = plot_middle(4)
    plt.plot(x, y1, label = "r = 0.05 cm")
    y1 = plot_middle(9)
    plt.plot(x, y1, label = "r = 0.10 cm")
    y1 = plot_middle(14)
    plt.plot(x, y1, label = "r = 0.15 cm")
    y1 = plot_middle(19)
    plt.plot(x, y1, label = "r = 0.20 cm")
    plt.show()

    plt.legend(loc='upper left')

def plot_contour1():
       
    f = open("temp.txt")

    temp_sol = f.read()

    temp_sol = temp_sol.split()

    temp_sol = map(float, temp_sol)
    
    temp_layer1 = np.ones([Nr[0], Nz[0]])

    for c in range(Nr[0]*Nz[0]):
        j = c / Nr[0]
        i = c % Nr[0]
        temp_layer1[i][j] = temp_sol[c]
    
    #n1 = 100
    R = np.linspace(1, 47, 47)
    Z = np.linspace(1, 20, 20)
    
    R, Z = np.meshgrid(R, Z, sparse=False, indexing='ij')
    
    plt.figure(1, figsize=(8,8))
    
    #adding x and y labels
    #plt.title('Temperature Contour of Layer 1 (60 s, 300 GW/m^3)')
    plt.xlabel('Finite Difference Element in Radial Direction')
    plt.ylabel('Finite Difference Element in Vertical Direction')
    
    contour_levels = [320.85, 320.90, 321.0, 321.10, 321.20, 321.3]
    
    plt.text(2.9, 20.0 + 0.06, '321.3 K', fontsize=10)
    plt.text(12.6, 20.0 + 0.06, '321.2 K', fontsize=10)
    plt.text(17.6, 20.0 + 0.06, '321.1 K', fontsize=10)
    plt.text(22.7, 20.0 + 0.06, '321.0 K', fontsize=10)
    plt.text(29.8, 20.0 + 0.06, '320.9 K', fontsize=10)
    plt.text(43.7, 20.0 + 0.06, '320.85 K', fontsize=10)

    #defined contour levels in plot with 4th argument
    plt.contour(R, Z, temp_layer1, contour_levels, colors='k')
    plt.show()
    
def plot_contour2():
       
    f = open("temp.txt")

    temp_sol = f.read()

    temp_sol = temp_sol.split()

    temp_sol = map(float, temp_sol)
    
    temp_layer2 = np.ones([Nr[1], Nz[1]])

    total_1 = Nr[0]*Nz[0]
    
    for c in range(Nr[1]*Nz[1]):
        j = c / Nr[1]
        i = c % Nr[1]
        temp_layer2[i][j] = temp_sol[c + total_1]
    
    #n1 = 100
    R = np.linspace(1, 20, 20)
    Z = np.linspace(1, 20, 20)
    
    R, Z = np.meshgrid(R, Z, sparse=False, indexing='ij')
    
    plt.figure(1, figsize=(8,8))
    
    #adding x and y labels
    #plt.title('Temperature Contour of Layer 1 (60 s, 300 GW/m^3)')
    plt.xlabel('Finite Difference Element in Radial Direction')
    plt.ylabel('Finite Difference Element in Vertical Direction')
    
    contour_levels = [321.20, 321.3, 321.4, 321.5, 321.6, 321.7]
    
    plt.text(9.5, 20.0 + 0.06, '321.7 K', fontsize=10)
    plt.text(16.5, 20.0 + 0.06, '321.6 K', fontsize=10)
    plt.text(20.10, 18.6, '321.5 K', fontsize=10)
    plt.text(20.10, 13.9, '321.4 K', fontsize=10)
    plt.text(20.10, 9.3, '321.3 K', fontsize=10)
    plt.text(20.10, 4.7, '321.2 K', fontsize=10)

    #defined contour levels in plot with 4th argument
    plt.contour(R, Z, temp_layer2, contour_levels, colors='k')
    plt.show()    
    
def plot_contour3():
       
    f = open("temp.txt")

    temp_sol = f.read()

    temp_sol = temp_sol.split()

    temp_sol = map(float, temp_sol)
    
    temp_layer3 = np.ones([Nr[2], Nz[2]])

    total_1 = Nr[0]*Nz[0]

    total_2 = Nr[1]*Nz[1]
    
    for c in range(Nr[2]*Nz[2]):
        j = c / Nr[2]
        i = c % Nr[2]
        temp_layer3[i][j] = temp_sol[c + total_1 + total_2]
    
    #n1 = 100
    R = np.linspace(1, 20, 20)
    Z = np.linspace(1, 20, 20)
    
    R, Z = np.meshgrid(R, Z, sparse=False, indexing='ij')
    
    plt.figure(1, figsize=(8,8))
    
    #adding x and y labels
    #plt.title('Temperature Contour of Layer 1 (60 s, 300 GW/m^3)')
    plt.xlabel('Finite Difference Element in Radial Direction')
    plt.ylabel('Finite Difference Element in Vertical Direction')
    
    contour_levels = [321.6, 321.7, 321.8, 321.9, 322.0, 322.1]
    
    plt.text(10.6, 20.0 + 0.06, '322.1 K', fontsize=10)
    plt.text(17.1, 20.0 + 0.06, '322.0 K', fontsize=10)
    plt.text(20.10, 18.1, '321.9 K', fontsize=10)
    plt.text(20.10, 13.1, '321.8 K', fontsize=10)
    plt.text(20.10, 8.3, '321.7 K', fontsize=10)
    plt.text(20.10, 3.5, '321.6 K', fontsize=10)

    #defined contour levels in plot with 4th argument
    plt.contour(R, Z, temp_layer3, contour_levels, colors='k')
    plt.show()   
    
def plot_contour4():
       
    f = open("temp.txt")

    temp_sol = f.read()

    temp_sol = temp_sol.split()

    temp_sol = map(float, temp_sol)
    
    temp_layer4 = np.ones([Nr[3], Nz[3]])

    total_1 = Nr[0]*Nz[0]

    total_2 = Nr[1]*Nz[1]

    total_3 = Nr[2]*Nz[2]
    
    for c in range(Nr[3]*Nz[3]):
        j = c / Nr[3]
        i = c % Nr[3]
        temp_layer4[i][j] = temp_sol[c + total_1 + total_2 + total_3]
    
    #n1 = 100
    R = np.linspace(1, 20, 20)
    Z = np.linspace(1, 20, 20)
    
    R, Z = np.meshgrid(R, Z, sparse=False, indexing='ij')
    
    plt.figure(1, figsize=(8,8))
    
    #adding x and y labels
    #plt.title('Temperature Contour of Layer 1 (60 s, 300 GW/m^3)')
    plt.xlabel('Finite Difference Element in Radial Direction')
    plt.ylabel('Finite Difference Element in Vertical Direction')
    
    contour_levels = [322.0, 322.1, 322.2, 322.3, 322.35]
    
    plt.text(7.6, 20.0 + 0.06, '322.35 K', fontsize=10)
    plt.text(12.8, 20.0 + 0.06, '322.3 K', fontsize=10)
    plt.text(18.1, 20.0 + 0.06, '322.2 K', fontsize=10)
    plt.text(20.10, 9.5, '322.1 K', fontsize=10)
    plt.text(20.10, 3.0, '322.0 K', fontsize=10)
 
    #defined contour levels in plot with 4th argument
    plt.contour(R, Z, temp_layer4, contour_levels, colors='k')
    plt.show()   
    
