"""
Question1 of Assignment2
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
del_x = 0.5 / num_mesh
k = 1000.0
A = 1 / 100.0
beta = k * A / del_x
T_A = 200.0
T_B = 600.0

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
    exact_sol[i] = 800.0 * x_mesh[i] + 200.0


# Defining the boundary conditions
def apply_bc_dia():
    lower_dia[0] = 0.0
    upper_dia[0] = -beta
    main_dia[0] = 3.0 * beta
    lower_dia[num_mesh - 1] = -beta
    upper_dia[num_mesh - 1] = 0.0
    main_dia[num_mesh - 1] = 3.0 * beta
    rhs[0] = 2.0 * beta * T_A
    rhs[num_mesh - 1] = 2.0 * beta * T_B


# Applying the boundary conditions
apply_bc_dia()

# Inserting values in diagonals
for i in range(1, len(numerical_sol) - 1):
    lower_dia[i] = -beta

    if i == len(numerical_sol) - 1:
        main_dia[i] = 3.0 * beta
        upper_dia[i] = 0.0
    else:
        main_dia[i] = 2.0 * beta
        upper_dia[i] = -beta

# Calling the function
tdma(num_mesh, lower_dia, upper_dia, main_dia, rhs, numerical_sol)

# Extending the arrays to plot boundary points
for i in range(0, 1):
    exact_sol = np.insert(exact_sol, 0, 200.0)
    exact_sol = np.insert(exact_sol, len(exact_sol), 600.0)
    x_mesh = np.insert(x_mesh, 0, 0.0)
    x_mesh = np.insert(x_mesh, len(x_mesh), 0.5)
    numerical_sol = np.insert(numerical_sol, 0, 200.0)
    numerical_sol = np.insert(numerical_sol, len(numerical_sol), 600.0)

# Plotting the required results
plt.plot(x_mesh, numerical_sol, 'g^', label='Numerical Solution')
plt.plot(x_mesh, exact_sol, 'r', label='Exact Solution')
plt.title('Computational solution and Exact solution for %d grid points using TDMA' % num_mesh)
plt.xlim(-0.01, 0.51)
plt.ylim(190, 610)
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.legend(fontsize='x-large', shadow=True)
plt.grid()
plt.show()
