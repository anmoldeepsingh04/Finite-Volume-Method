"""
Question2 of Assignment2
Completed by: Anmoldeep Singh 180030002
Dated: 01-04-2021
"""
# Importing relevant libraries and modules
import numpy as np
import matplotlib.pyplot as plt


# defining the TDMA Algorithm
def tdma(num, low, up, dia, func, sol):
    d1 = np.zeros(num, float)
    rhs1 = np.zeros(num, float)
    d1[0] = dia[0]
    rhs1[0] = func[0]
    for j in range(1, len(dia)):
        d1[j] = dia[j] - (low[j] * up[j - 1] / d1[j - 1])
        rhs1[j] = func[j] - (low[j] * rhs1[j - 1] / d1[j - 1])
    sol[num - 1] = rhs1[num - 1] / d1[num - 1]
    for p in range(len(dia) - 2, -1, -1):
        sol[p] = (rhs1[p] - up[p] * sol[p + 1]) / d1[p]


# Defining required quantities
num_mesh = int(input("Enter the total number of internal mesh points: "))
L = 0.02
del_x = L / num_mesh
k = 0.5
q = 1000000.0
alpha = q*del_x
beta = k/del_x
T_A = 100.0
T_B = 200.0

# Defining required arrays
x_mesh = np.zeros(num_mesh, float)
exact_sol = np.zeros(num_mesh, float)
numerical_sol = np.zeros(num_mesh, float)
rhs = np.zeros(num_mesh, float)
main_dia = np.zeros(num_mesh, float)
lower_dia = np.zeros(num_mesh, float)
upper_dia = np.zeros(num_mesh, float)

# Calculating grid points and exact solution
for i in range(0, num_mesh):
    x_mesh[i] = (i + 0.5) * del_x
    x = x_mesh[i]
    exact_sol[i] = T_A + ((T_B - T_A) * x/L) + (q*L*x*(1 - (x/L))/(2*k))
    rhs[i] = alpha

# Inserting values in diagonals
for i in range(0, len(numerical_sol)):
    lower_dia[i] = -beta
    main_dia[i] = 2.0 * beta
    upper_dia[i] = -beta


# Defining the boundary conditions
def apply_bc_dia():
    lower_dia[0] = 0.0
    upper_dia[0] = -beta
    main_dia[0] = 3.0*beta
    lower_dia[num_mesh - 1] = -beta
    upper_dia[num_mesh - 1] = 0.0
    main_dia[num_mesh - 1] = 3.0*beta
    rhs[0] = alpha + 2.0 * beta * T_A
    rhs[num_mesh - 1] = alpha + 2.0 * beta * T_B


# Applying the boundary conditions
apply_bc_dia()

# Calling the function
tdma(num_mesh, lower_dia, upper_dia, main_dia, rhs, numerical_sol)

# Extending the arrays to plot boundary points
for i in range(0, 1):
    exact_sol = np.insert(exact_sol, 0, T_A)
    exact_sol = np.insert(exact_sol, len(exact_sol), T_B)
    x_mesh = np.insert(x_mesh, 0, 0.0)
    x_mesh = np.insert(x_mesh, len(x_mesh), L)
    numerical_sol = np.insert(numerical_sol, 0, T_A)
    numerical_sol = np.insert(numerical_sol, len(numerical_sol), T_B)

# Plotting the required results
plt.plot(x_mesh, numerical_sol, 'b*', label='Numerical Solution')
plt.plot(x_mesh, exact_sol, 'r', label='Exact Solution')
plt.title('Computational solution and Exact solution for %d grid points using TDMA' % num_mesh)
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.legend(fontsize = 'x-large', shadow = True, loc = 'upper left')
plt.grid()
plt.show()
