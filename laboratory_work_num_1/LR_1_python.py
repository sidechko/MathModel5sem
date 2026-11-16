import numpy as np
import matplotlib.pyplot as plt

# --- clear() / clf() ---
plt.close('all')

t = np.linspace(1.0, 3.0, 2001)

W1 = 1.2 * 2 * np.pi 
W2 = 0.6 * 2 * np.pi
W0 = 20.0 * np.pi
K  = 10.0

A = np.sin(t * W1) + np.cos(t * W2)
Famp = A * np.cos(t * W0)
Ffm = np.cos((W0 + K * A) * t)

# --- Построение графиков (3 подграфика, как subplot(3,1,N)) ---
fig, axes = plt.subplots(3, 1, figsize=(10, 9))

# --- Подграфик 1: A, -A, Famp ---
ax = axes[0]
ax.plot(t,  A,    label='A')
ax.plot(t, -A,    label='-A')
ax.plot(t, Famp,  label='F AM')
ax.legend()
ax.set_title('F amp')

# --- Подграфик 2: A, Ffm ---
ax = axes[1]
ax.plot(t, A,    label='A')
ax.plot(t, Ffm,  label='F FM')
ax.legend()
ax.set_title('F fm')

# --- Подграфик 3: A, -A, Famp, Ffm ---
ax = axes[2]
ax.plot(t, A,    label='A')
ax.plot(t, -A,   label='-A')
ax.plot(t, Famp, label='F AM')
ax.plot(t, Ffm,  label='F FM')
ax.set_title('Mixed')

fig.tight_layout()
plt.show()