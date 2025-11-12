import numpy as np
import matplotlib.pyplot as plt

# Parameters (example values, adapt to match your test case)
ka = 11.4
kb = 0.15
kc = 1.8
kd = 0.65
fa = 0.9         # radians
fb = 1.02        # radians
fc = -1.55        # radians
fd = -2.5        # radians
epsd = 2
eps0 = 0.27
eps1 = 0.14
eps2 = 0.26

phi_vals = np.linspace(-np.pi, np.pi, 500)
energy_vals = []

for phi in phi_vals:
    dphia = phi - fa
    dphib = phi - fb
    dphib2 = dphib + 2.0*np.pi
    dphic = phi - fc
    dphic2 = dphic - 2.0*np.pi
    dphid = phi - fd
    dphid2 = dphid - 2.0*np.pi

    pa = -ka * dphia
    pb = -kb * dphib**3
    pb2 = -kb * dphib2**3
    pc = -kc * dphic
    pc2 = -kc * dphic2
    pd = -kd * dphid**3
    pd2 = -kd * dphid2**3

    ea = pa * dphia - epsd
    eb = pb * dphib + eps0
    eb2 = pb2 * dphib2 + eps0
    ec = pc * dphic + epsd + eps1
    ec2 = pc2 * dphic2 + epsd + eps1
    ed = pd * dphid + eps1 + eps2
    ed2 = pd2 * dphid2 + eps1 + eps2

    fea = np.exp(ea)
    feb = np.exp(eb)
    feb2 = np.exp(eb2)
    fec = np.exp(ec)
    fec2 = np.exp(ec2)
    fed = np.exp(ed)
    fed2 = np.exp(ed2)

    Z = fea + feb + feb2 + fec + fec2 + fed + fed2
    energy = -np.log(Z)

    energy_vals.append(energy)

# Plot
plt.figure(figsize=(8, 5))
plt.plot(phi_vals, energy_vals, label="Dihedral Energy", color='navy')
plt.xlabel("ϕ (radians)")
plt.ylabel("Energy (kJ/mol)")
plt.title("Double-Gaussian Dihedral Potential")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()