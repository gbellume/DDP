import numpy as np
import matplotlib.pyplot as plt

sigma = np.array([0, 1, 2, 3, 5, 10])
short_term = np.array([8.2755E-10, 1.6917E-09, 3.4814E-09, 5.3440E-09, 9.2152E-09, 1.9564E-08])
long_term = np.array([8.3065e-13, 2.1277e-12, 5.5478e-12, 9.8094e-12, 2.0280e-11, 5.5009e-11])

plt.figure(figsize=(8,6))
plt.semilogy(sigma, long_term, marker='s', color='b', label='Long Term (Worst Fluence)')
plt.xlabel(r'Cross Section [$cm^2/bit$]')
plt.ylabel(r'SEU Rate [$bit^{-1}/s$]')

plt.semilogy(sigma, short_term, marker='v', color='r', label='Short Term (Worst Flux)')
plt.xlabel(r'Saturation Cross Section [$cm^2/bit$]')
plt.ylabel(r'SEU Rate [$bit^{-1}/s$]')
plt.title(r'Estimation of SEU rate for different $\sigma_{sat}$')
plt.grid(which='major', linewidth=0.3)
plt.legend()
plt.show()