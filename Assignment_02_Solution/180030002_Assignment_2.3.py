"""
Question3 of Assignment2
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
    for l in range(1, len(dia)):
        d1[l] = dia[l] - (low[l] * up[l - 1] / d1[l - 1])
        rhs1[l] = func[l] - (low[l] * rhs1[l - 1] / d1[l - 1])
    sol[num - 1] = rhs1[num - 1] / d1[num - 1]
    for k in range(len(dia) - 2, -1, -1):
        sol[k] = (rhs1[k] - up[k] * sol[k + 1]) / d1[k]


# Defining required quantities
num_mesh = 20
L = 0.6
del_x = L / num_mesh
h = 25.0
k1 = 20.0
k2 = 1.5
k3 = 50.0
R_total = (1 / h) + (0.3 / k1) + (0.15 / k2) + (0.15 / k3)
alpha = h * del_x
beta1 = k1 / del_x
beta2 = k2 / del_x
beta3 = k3 / del_x
gamma1_Harmonic = (2 * k1 * k2) / (del_x * (k1 + k2))
gamma2_Harmonic = (2 * k2 * k3) / (del_x * (k2 + k3))
gamma1_Arithmetic = (k1 + k2) / (2.0 * del_x)
gamma2_Arithmetic = (k2 + k3) / (2.0 * del_x)
gamma = (gamma1_Harmonic, gamma2_Harmonic, gamma1_Arithmetic, gamma2_Arithmetic)
T_inf = 1073.0
T_B = 293.0
T_R = (T_inf - T_B) / R_total
T_A = T_inf - (T_R * (1 / h))
T_1 = (T_inf - T_A) * h

# Defining required arrays
x_mesh = np.zeros(num_mesh, float)
numerical_sol = np.zeros(num_mesh, float)
rhs = np.zeros(num_mesh, float)
main_dia = np.zeros(num_mesh, float)
lower_dia = np.zeros(num_mesh, float)
upper_dia = np.zeros(num_mesh, float)
exact_sol = np.zeros(num_mesh, float)

# Calculating solution for different formulations of k
for j in range(0, 3, 2):
    gamma1 = gamma[j]
    gamma2 = gamma[j + 1]

    # Calculating grid points
    for i in range(0, num_mesh):
        x_mesh[i] = (i + 0.5) * del_x

    # Inserting values in upper and lower diagonals
    for i in range(1, num_mesh):
        if 1 <= i <= 9:
            lower_dia[i] = -beta1
            upper_dia[i - 1] = -beta1
        if i == 10:
            lower_dia[i] = -gamma1
            upper_dia[i - 1] = -gamma1
        if 11 <= i <= 14:
            lower_dia[i] = -beta2
            upper_dia[i - 1] = -beta2
        if i == 15:
            lower_dia[i] = -gamma2
            upper_dia[i - 1] = -gamma2
        if 16 <= i < num_mesh:
            lower_dia[i] = -beta3
            upper_dia[i - 1] = -beta3

    # Inserting values in main diagonal
    for i in range(1, num_mesh):
        if 1 <= i <= 8:
            main_dia[i] = 2 * beta1
        if i == 9:
            main_dia[i] = gamma1 + beta1
        if i == 10:
            main_dia[i] = gamma1 + beta2
        if 11 <= i <= 13:
            main_dia[i] = 2 * beta2
        if i == 14:
            main_dia[i] = gamma2 + beta2
        if i == 15:
            main_dia[i] = gamma2 + beta3
        if 16 <= i <= num_mesh:
            main_dia[i] = 2 * beta3

    # Defining the boundary conditions
    def apply_bc_dia():
        lower_dia[0] = 0.0
        upper_dia[0] = -beta1
        main_dia[0] = 3.0 * beta1
        lower_dia[num_mesh - 1] = -beta3
        upper_dia[num_mesh - 1] = 0.0
        main_dia[num_mesh - 1] = 3.0 * beta3
        rhs[0] = (2.0 * beta1 + alpha) * T_A - alpha * T_inf
        rhs[num_mesh - 1] = 2.0 * beta3 * T_B


    # Applying the boundary conditions
    apply_bc_dia()

    # Calling the function
    tdma(num_mesh, lower_dia, upper_dia, main_dia, rhs, numerical_sol)

    # Extending the arrays to plot boundary points
    for i in range(0, 1):
        x_mesh = np.insert(x_mesh, 0, 0.0)
        x_mesh = np.insert(x_mesh, len(x_mesh), L)
        numerical_sol = np.insert(numerical_sol, 0, T_A)
        numerical_sol = np.insert(numerical_sol, len(numerical_sol), T_B)

        # Plotting the numerical for harmonic and arithmetic formulations results respectively
        if j == 0:
            plt.plot(x_mesh, numerical_sol, 'bs--', label='Harmonic Formulation')
        if j != 0:
            plt.plot(x_mesh, numerical_sol, 'r*--', label='Arithmetic Formulation')

# Defining required arrays for exact solution
x_mesh_exact = np.zeros(num_mesh, float)
exact_sol = np.zeros(num_mesh, float)
exact_sol[0] = T_A

# Calculating grid points and solution for exact solution
for i in range(0, num_mesh):
    x_mesh_exact[i] = i * del_x
for i in range(0, num_mesh):
    if 0 <= i <= 9:
        exact_sol[i + 1] = exact_sol[i] - T_1 * del_x / k1
    if 10 <= i <= 14:
        exact_sol[i + 1] = exact_sol[i] - T_1 * del_x / k2
    if 15 <= i < num_mesh - 1:
        exact_sol[i + 1] = exact_sol[i] - T_1 * del_x / k3

# Extending the arrays to plot boundary points
x_mesh_exact = np.insert(x_mesh_exact, len(x_mesh_exact), L)
exact_sol = np.insert(exact_sol, len(exact_sol), T_B)

# Plotting the required results
plt.plot(x_mesh_exact, exact_sol, 'gd-', label='Exact Solution')
plt.title('Computational solution and Exact solution for %d grid points using TDMA' % num_mesh)
plt.xlim(-0.01, 0.61)
plt.ylim(200, 900)
plt.xlabel('$x$', fontsize=14)
plt.ylabel('$Temperature$', fontsize=14)
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(fontsize='x-large', shadow=True)
plt.minorticks_on()
plt.show()
