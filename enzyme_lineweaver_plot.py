import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Data from the user
substrate_mM = np.array([0.200, 0.250, 0.333, 0.500, 1.000])
# Columns: Condition1..Condition5
v_data = np.array([
    [0.222, 0.200, 0.182, 0.167, 0.154],
    [0.250, 0.222, 0.200, 0.182, 0.167],
    [0.286, 0.250, 0.222, 0.200, 0.182],
    [0.333, 0.286, 0.250, 0.222, 0.200],
    [0.400, 0.333, 0.286, 0.250, 0.222]
])  # shape (5 substrate points, 5 conditions)

# Transpose to shape (5 conditions, 5 substrate points)
v_data = v_data.T

# 1) Transform the data to integers (scale by 1000 and round)
# We'll present integer-transformed V0 (units: original umol/s * 1000 -> integer units)
v_data_int = np.rint(v_data * 1000).astype(int)

# For calculations (Km and Vmax) use the original float v_data

# Prepare Lineweaver-Burk variables: x = 1/[S], y = 1/v
inv_S = 1.0 / substrate_mM
inv_v = 1.0 / v_data  # shape (5 conditions, 5 points)

conditions = [f"Condition {i+1}" for i in range(v_data.shape[0])]
colors = ["C0","C1","C2","C3","C4"]

# Function to compute linear fit using closed-form least squares and print steps
def compute_linear_fit_steps(x, y, cond_name, verbose=True):
    # x and y are numpy arrays
    n = x.size
    sum_x = x.sum()
    sum_y = y.sum()
    sum_xy = (x * y).sum()
    sum_x2 = (x * x).sum()

    # slope = (n*sum(xy) - sum(x)*sum(y)) / (n*sum(x^2) - (sum(x))^2)
    slope_num = n * sum_xy - sum_x * sum_y
    slope_den = n * sum_x2 - sum_x * sum_x
    slope = slope_num / slope_den

    # intercept = (sum(y) - slope * sum(x)) / n
    intercept = (sum_y - slope * sum_x) / n

    # Vmax = 1 / intercept ; Km = slope * Vmax
    Vmax = 1.0 / intercept
    Km = slope * Vmax

    if verbose:
        print(f"\nDetailed linear fit steps for {cond_name}:")
        print(f"n = {n}")
        print("1/[S] values (x):", ", ".join(f"{v:.6f}" for v in x))
        print("1/V0 values (y):", ", ".join(f"{v:.6f}" for v in y))
        print(f"sum(x) = {sum_x:.6f}")
        print(f"sum(y) = {sum_y:.6f}")
        print(f"sum(x*y) = {sum_xy:.6f}")
        print(f"sum(x^2) = {sum_x2:.6f}")
        print(f"slope numerator = n*sum(xy) - sum(x)*sum(y) = {slope_num:.6f}")
        print(f"slope denominator = n*sum(x^2) - (sum(x))^2 = {slope_den:.6f}")
        print(f"slope = {slope_num:.6f} / {slope_den:.6f} = {slope:.6f}")
        print(f"intercept = (sum(y) - slope*sum(x)) / n = ({sum_y:.6f} - {slope:.6f}*{sum_x:.6f}) / {n} = {intercept:.6f}")
        print(f"Vmax = 1 / intercept = 1 / {intercept:.6f} = {Vmax:.6f} umol/s")
        print(f"Km = slope * Vmax = {slope:.6f} * {Vmax:.6f} = {Km:.6f} mM")

    return float(slope), float(intercept), float(Vmax), float(Km)

# Fit linear regression (1/v = slope*(1/[S]) + intercept) for each condition using closed-form
fits = {}
for i, cond in enumerate(conditions):
    s, itc, Vmax, Km = compute_linear_fit_steps(inv_S, inv_v[i], cond, verbose=False)
    fits[cond] = {"slope": s, "intercept": itc, "Vmax": Vmax, "Km": Km}

# (Detailed steps for chosen conditions will be printed later after the conditions are defined.)

# Choose the two conditions for Km and Vmax calculation/display.
# Assumption: "the two conditions" = Condition 1 (control) and Condition 5 (highest inhibitor).
cond_a = "Condition 1"
cond_b = "Condition 5"

# Print detailed steps for all conditions (1 through 5)
print("\n=== Detailed calculations for all conditions ===")
for i, cond in enumerate(conditions):
    compute_linear_fit_steps(inv_S, inv_v[i], cond, verbose=True)

# Plot Lineweaver-Burk: 1/[S] vs 1/v for all conditions and draw fitted lines
fig, ax = plt.subplots(figsize=(8,6))
for i, cond in enumerate(conditions):
    ax.scatter(inv_S, inv_v[i], label=cond, color=colors[i], s=60)
    # fitted line
    xs = np.linspace(inv_S.min()*0.9, inv_S.max()*1.1, 100)
    ys = fits[cond]["slope"] * xs + fits[cond]["intercept"]
    ax.plot(xs, ys, color=colors[i], alpha=0.6)

# Annotate Vmax and Km for all conditions on the plot (concise summary)
lines = []
for i, cond in enumerate(conditions):
    f = fits[cond]
    lines.append(f"{cond}: Vmax={f['Vmax']:.4f}, Km={f['Km']:.4f}")

# Put the text block in the upper-left corner of the axes
text_block = "\n".join(lines)
ax.text(0.02, 0.98, text_block, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

# Move legend if it overlaps; put it at lower right
ax.legend(loc='lower right')

ax.set_xlabel('1 / [S] (mM^-1)')
ax.set_ylabel('1 / V0 (s / umol)')
ax.set_title('Lineweaver-Burk Plot')
ax.legend()
ax.grid(True)

# Save plot
out_path = Path(__file__).resolve().parent / 'lineweaver_plot.png'
fig.tight_layout()
fig.savefig(out_path, dpi=200)

# Print integer-transformed table and calculated Km/Vmax for chosen conditions
print("Integer-transformed V0 (original umol/s scaled by 1000 and rounded):")
print("Substrate_mM:\t", '\t'.join([f"{s:.3f}" for s in substrate_mM]))
for i, cond in enumerate(conditions):
    ints = '\t'.join(str(x) for x in v_data_int[i])
    print(f"{cond}:\t{ints}")

print('\nCalculated Lineweaver-Burk fits and derived parameters:')
for cond in conditions:
    f = fits[cond]
    print(f"{cond}: slope = {f['slope']:.6f}, intercept = {f['intercept']:.6f}, Vmax = {f['Vmax']:.6f} umol/s, Km = {f['Km']:.6f} mM")

print(f"\nPlot saved to: {out_path}")

# Also print a short determination of inhibition type based on Km and Vmax comparison
Vmax_a = fits[cond_a]['Vmax']
Km_a = fits[cond_a]['Km']
Vmax_b = fits[cond_b]['Vmax']
Km_b = fits[cond_b]['Km']

if np.isclose(Vmax_a, Vmax_b, rtol=0.05) and (not np.isclose(Km_a, Km_b, rtol=0.05)):
    inh_type = 'Competitive (Vmax ~ same; Km changed)'
elif (not np.isclose(Vmax_a, Vmax_b, rtol=0.05)) and np.isclose(Km_a, Km_b, rtol=0.05):
    inh_type = 'Noncompetitive (Vmax changed; Km ~ same)'
elif (not np.isclose(Vmax_a, Vmax_b, rtol=0.05)) and (not np.isclose(Km_a, Km_b, rtol=0.05)):
    inh_type = 'Mixed inhibition (both Vmax and Km changed)'
else:
    inh_type = 'No clear change (values similar)'

print(f"\nInhibition determination (Condition 1 vs Condition 5): {inh_type}")
