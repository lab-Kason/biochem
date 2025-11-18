#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Data
S_mM = np.array([0.2, 0.25, 0.333, 0.5, 1.0])
v_no_inhib = np.array([0.167, 0.2, 0.25, 0.333, 0.5])
v_with_inhib = np.array([0.143, 0.167, 0.2, 0.25, 0.333])

# Scale to integers (requirement 1)
scale = 1000
S_int = (S_mM * scale).astype(int)
v_no_int = np.rint(v_no_inhib * scale).astype(int)
v_with_int = np.rint(v_with_inhib * scale).astype(int)

print('Transformed integers (scale = 1000):')
print('Substrate_mM_original, Substrate_scaled_int, v_no_orig, v_no_int, v_with_orig, v_with_int')
for s_o, s_i, vno_o, vno_i, vwi_o, vwi_i in zip(S_mM, S_int, v_no_inhib, v_no_int, v_with_inhib, v_with_int):
    print(f'{s_o:.3f}, {s_i}, {vno_o:.3f}, {vno_i}, {vwi_o:.3f}, {vwi_i}')

# Lineweaver-Burk reciprocals
invS = 1.0 / S_mM
invv_no = 1.0 / v_no_inhib
invv_with = 1.0 / v_with_inhib

# Linear fits
m_no, b_no = np.polyfit(invS, invv_no, 1)
m_with, b_with = np.polyfit(invS, invv_with, 1)

Vmax_no = 1.0 / b_no
Km_no = m_no / b_no
Vmax_with = 1.0 / b_with
Km_with = m_with / b_with

print('\nCalculated parameters:')
print(f'Without inhibitor: Vmax = {Vmax_no:.6f} micromol/min, Km = {Km_no:.6f} mM')
print(f'With inhibitor:    Vmax = {Vmax_with:.6f} micromol/min, Km = {Km_with:.6f} mM')

# Differences
dV = Vmax_with - Vmax_no
dK = Km_with - Km_no
print('\nDifferences (with - without):')
print(f'  Delta Vmax = {dV:.6f} micromol/min')
print(f'  Delta Km   = {dK:.6f} mM')

# Interpretation
eps = 1e-9
if abs(dV) < 0.05 * abs(Vmax_no + eps):
    # Vmax approx unchanged
    if dK > 0.05 * abs(Km_no + eps):
        interpretation = 'Competitive inhibition (Vmax ~ unchanged, Km increased)'
    elif dK < -0.05 * abs(Km_no + eps):
        interpretation = 'Uncompetitive-like (Km decreased while Vmax unchanged)'
    else:
        interpretation = 'No clear change or very small differences'
else:
    if abs(dK) < 0.05 * abs(Km_no + eps):
        interpretation = 'Noncompetitive inhibition-like (Vmax decreased, Km ~ unchanged)'
    else:
        interpretation = 'Mixed or uncompetitive inhibition (both Vmax and Km changed)'

print('\nInterpretation:')
print(interpretation)

# Plot
x_line = np.linspace(invS.min()*0.9, invS.max()*1.1, 200)
y_no_line = m_no * x_line + b_no
y_with_line = m_with * x_line + b_with

plt.figure(figsize=(7,6))
plt.scatter(invS, invv_no, label='No inhibitor', color='C0', zorder=5)
plt.plot(x_line, y_no_line, color='C0', linestyle='--', label=f'Fit no inh: slope={m_no:.3f}, intercept={b_no:.3f}')
plt.scatter(invS, invv_with, label='With inhibitor', color='C1', zorder=5)
plt.plot(x_line, y_with_line, color='C1', linestyle='--', label=f'Fit with inh: slope={m_with:.3f}, intercept={b_with:.3f}')
plt.xlabel('1 / [S] (1/mM)')
plt.ylabel('1 / v (min / micromol)')
plt.title('Lineweaver–Burk plot')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)
text_x = invS.max()*0.6
plt.text(text_x, (invv_no.max()+invv_with.max())*0.5, f'No inh: Vmax={Vmax_no:.3f}\nKm={Km_no:.3f} mM', color='C0')
plt.text(text_x, (invv_no.max()+invv_with.max())*0.35, f'With inh: Vmax={Vmax_with:.3f}\nKm={Km_with:.3f} mM', color='C1')

out_path = '/Users/kasonchiu/Desktop/jupyter/lineweaver_plot_Kason.png'
plt.savefig(out_path, dpi=150)
print(f'\nPlot saved to: {out_path}')
