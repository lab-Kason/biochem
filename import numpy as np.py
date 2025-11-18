import numpy as np
import matplotlib.pyplot as plt

# Transformed integer data (units: [S] in µM, v in nmol/min)
S = np.array([200, 250, 333, 500, 1000], dtype=float)   # µM
v_no = np.array([167, 200, 250, 333, 500], dtype=float)  # nmol/min
v_in = np.array([143, 167, 200, 250, 333], dtype=float)  # nmol/min

# reciprocals for Lineweaver-Burk
x = 1.0 / S            # µM^-1
y_no = 1.0 / v_no      # (nmol/min)^-1
y_in = 1.0 / v_in

# linear fits (slope m, intercept b)
m_no, b_no = np.polyfit(x, y_no, 1)
m_in, b_in = np.polyfit(x, y_in, 1)

# derived parameters: Vmax and Km (return to original convenient units)
Vmax_no_nm = 1.0 / b_no           # nmol/min
Km_no_uM = m_no * Vmax_no_nm      # µM

Vmax_in_nm = 1.0 / b_in
Km_in_uM = m_in * Vmax_in_nm

# print results (rounded)
print("Without inhibitor: Vmax ≈ {:.3g} nmol/min ({:.3g} µmol/min), Km ≈ {:.3g} µM ({:.3g} mM)".format(
    Vmax_no_nm, Vmax_no_nm/1000, Km_no_uM, Km_no_uM/1000))
print("With inhibitor:    Vmax ≈ {:.3g} nmol/min ({:.3g} µmol/min), Km ≈ {:.3g} µM ({:.3g} mM)".format(
    Vmax_in_nm, Vmax_in_nm/1000, Km_in_uM, Km_in_uM/1000))

# plotting
xs = np.linspace(0, x.max()*1.1, 200)
ys_no_fit = m_no * xs + b_no
ys_in_fit = m_in * xs + b_in

plt.figure(figsize=(6,5))
plt.scatter(x, y_no, color='C0', label='No inhibitor (data)')
plt.plot(xs, ys_no_fit, color='C0', linestyle='--', label='No inhibitor (fit)')
plt.scatter(x, y_in, color='C1', label='With inhibitor (data)')
plt.plot(xs, ys_in_fit, color='C1', linestyle='--', label='With inhibitor (fit)')
plt.xlabel('1 / [S] (µM⁻¹)')
plt.ylabel('1 / v (min / nmol)')
plt.title('Lineweaver–Burk plot')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()