import pyfits
import matplotlib.pyplot as plt
import numpy as np
from scipy import (radians, polyfit, poly1d)
from scipy.optimize import curve_fit


#Antenna Diping

w_a = 0.0898

w_s = [0.0886, 0.0884, 0.0883, 0.088, 0.0873, 0.0867, 0.0858, 0.0841,
        0.0812, 0.0463]

z = [4.99, 5.47, 6.04, 6.76, 7.66, 8.85, 10.48, 12.84, 16.6, 23.5]

ln = []
for i in range(len(w_s)):
    W = w_a - w_s[i]
    ln.append(np.log(W))

sec= []

for i in range(len(z)):
    S =np.sin(radians(z[i]))
    sec.append(1. / S)

coef = polyfit(sec, ln, 1)

recta = poly1d(coef)

fig1, ax1 = plt.subplots()
ax1.plot(sec, recta(sec))
ax1.plot(sec, ln, '*')
ax1.set_xlabel('sec(z)')
ax1.set_ylabel('log(W_a - W_sky)')
ax1.set_title('Polyfit datos para calcular opacidad atmosfera')
