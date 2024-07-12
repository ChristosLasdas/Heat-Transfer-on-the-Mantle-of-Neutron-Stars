import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scienceplots
import time

start_time = time.time()

plt.style.use(['science', 'notebook'])
plt.rcParams.update({"text.usetex" : True})


n = 10000 #number of iterations. iterations here are a measure of time

nx = 21 #x resolution
ny = 21 #y resolution


Temp = np.zeros((nx,ny))
Eq = np.zeros((nx,ny))
Result = np.zeros((3,nx*ny))


x = np.zeros(nx)
y = np.zeros(ny)


x_min = 0.
x_max = 1.


y_min = 0.
y_max = 1.


t = 0.

λ = 0.
k_e = 1.
C_v = 1.

dx = (x_max - x_min) / (nx - 1)
dy = (y_max - y_min) / (ny - 1)

dt = dx*dy*0.1


#Setting the values of x and y

for i in range(0,nx):
    x[i] = i*dx
    
for i in range(0,ny):
    y[i] = i*dy


#Boundary condition Temp(x,y_min) = 0
for i in range(0,nx):
    Temp[i, 0] = 0.



#Boundary conditions Temp(x_min,y) = 0
#Boundary conditions Temp(x_max,y) = 0
for i in range(0,ny):
    Temp[0, i] = 0.
    Temp[nx-1, i] = 0.
    

#Solution of Heat Equation without Magnetic Field

for k in range(0, n):
    
    #In this loop we evaluate the derivatives
    
    for i in range (1, nx-1):
        for j in range(1, ny-1):
            dxxTemp = (Temp[i+1, j] - 2*Temp[i, j] + Temp[i-1, j]) / dx**2
            dyyTemp = (Temp[i, j+1] - 2*Temp[i, j] + Temp[i, j-1]) / dy**2
                
            
            Eq[i,j] = (k_e*(dxxTemp + dyyTemp) + \
                1e3*np.exp(-((x[i] - 0.5)**2 + (y[j] - 0.5)**2) / (0.1)**2)*np.exp(-λ*t)) / C_v
        
       # At the top border, Luminosity = dT/dy = a0 * T^4
    for i in range(1, nx-1):
        Temp[i, ny-1] = Temp[i,ny-2] - ((Temp[i, ny-2])**4)*dy    
    
    #In this loop we evaluate the new value of Temp 
    for i in range (1, nx-1):
        for j in range(1, ny-1): 
            Temp[i,j] = Temp[i,j] + Eq[i,j]*dt

    t = t + dt

#print(f"t = {t}")

#In this loop we save the results in a single array           
for i in range(0,nx):
    for j in range(0, ny):
        k = i*nx + j
        Result[0, k] = x[i]
        Result[1, k] = y[j]
        Result[2, k] = Temp[i,j]

#Here we plot the results
        
X = Result[0, :].reshape(nx,ny)
Y = Result[1, :].reshape(nx,ny)
Z = Result[2, :].reshape(nx,ny)

fig = plt.figure(dpi=1000)
ax1 = fig.add_subplot(1,1,1, aspect = 1, xlim = [x_min, x_max], ylim = [y_min, y_max])

plt.contourf(X, Y, Z, levels = 50)
ax1 = ax1.contourf(X,Y,Z)
plt.colorbar()
plt.title("Thermal Evolution without Magnetic Field" + "\n" + f"$\lambda$ = {λ}, $t$ = {round(t)}")
plt.xlabel("x", fontsize = 16)
plt.ylabel("y", fontsize = 16)
 
plt.show()

end_time = time.time()

print("Time taken: {:.2f} seconds".format(end_time - start_time))