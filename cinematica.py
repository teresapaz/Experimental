import pyfits
import matplotlib.pyplot as plt
import numpy as np
from scipy import radians
from scipy.optimize import curve_fit

#ctes
R_sun = 2.6228*(10**17)
V_sun = 220.
G = 6.674*(10**-11)
M_sun = 1*10**12
radii = 15. #kpc
thick = 0.6 #kpc

#esta funcion les entrega un arreglo con los valores reales del eje correspondiente
#1 es velocidad, 2 es longirud galactica, 3 es latitud galactica
def values(h,j):
	N=h['NAXIS'+str(j)];
	val=np.zeros(N);
	for i in range(0,N):
		val[i] = (i+1-float(h['CRPIX'+str(j)]))*float(h['CDELT'+str(j)]) + float(h['CRVAL'+str(j)]);
	return val;

cubo	 = pyfits.open("Cubo_de_datos.fits") #abrir objeto cubo de datos
data 	 = cubo[0].data #extraer matriz de datos
header = cubo[0].header #extraer el header del archivo fits

# data_cube, header_data_cube = pyfits.getdata("data_cube.fits", 0, header=True)

#Estos seran los tres arreglos con los valores reales de los tres ejes del cubo
vel = values(header,1)
lon = values(header,2)
lat = values(header,3)


# Velocidad Terminal
Vmax = []
bmax = []
for l in range(len(lon)):
	Vb = []
	b = []
	V = []
	for b in range(len(lat)):
		A = data[b, l, :]
		if max(A) >= 0.9:
			for v in range(len(vel)):
				if A[v] >= 0.9:
					Vb.append(abs(vel[v]))
			V.append(max(Vb))
		else:
			V.append(0)
	Vmax.append(max(V))
	bmax.append(V.index(max(V)))

b_ang = []
for b in bmax:
	b_ang.append(lat[b])

#curva rotacion
def W(v, l):
	w_sun = V_sun / R_sun
	w1 = v / (R_sun*np.sin(radians(l)))
	return -w1 + w_sun

def V(v,l):
	return abs(W(v,l)*R_sun*np.sin(radians(l)))


#corrugacion:
def z(l,b):
	p = 3.0857 * (10**13)
	z = (R_sun / p) * np.cos(radians(l)) * np.tan(radians(b))
	return z

#modelos de masa:

def m_puntual(r, m):
	v = np.sqrt(G * m / r)
	return v

def m_disco(r, s):
	M = np.pi * (r**2) * s
	v = np.sqrt(G * M / r)
	return v

def m_esfera(r, p):
	M = np.pi * (r**3) * p * (4./3)
	v = np.sqrt(G * M / r)
	return v

def m_puntualesfera(r, m, p):
	M = (np.pi * (r**3) * p * (4./3)) + m
	v = np.sqrt(G * M / r)
	return v

def m_puntualdisco(r, m, s):
	M = (np.pi * (r**2) * s) + m
	v = np.sqrt(G * M / r)
	return v

opt_puntual, cov_puntual = curve_fit(m_puntual, -np.sin(radians(lon)), V(Vmax,lon), M_sun)
opt_disco, cov_disco = curve_fit(m_disco, -np.sin(radians(lon)), V(Vmax,lon), (M_sun / (np.pi * radii**2)))
opt_esfera, cov_esfera = curve_fit(m_esfera, -np.sin(radians(lon)), V(Vmax,lon), (M_sun / (np.pi * radii**3 * (4./3))))
opt_puntualdisco, cov_puntualdisco = curve_fit(m_puntualdisco, -np.sin(radians(lon)), V(Vmax,lon), [M_sun, (M_sun / (np.pi * radii**2))])
opt_puntualesfera, cov_puntualesfera = curve_fit(m_puntualesfera, -np.sin(radians(lon)), V(Vmax,lon), [M_sun, (M_sun / (np.pi * radii**3 * (4./3)))])

mdisco_opt, sdisco_opt = opt_puntualdisco
mesfera_opt, pesfera_opt = opt_puntualesfera

fig1, ax1 = plt.subplots()
ax1.plot(-np.sin(radians(lon)), V(Vmax,lon), '*')
ax1.plot(-np.sin(radians(lon)), m_puntual(-np.sin(radians(lon)), opt_puntual),
			'-', label='Masa Puntual')
ax1.plot(-np.sin(radians(lon)), m_disco(-np.sin(radians(lon)), opt_disco), '-',
 			label='Disco Uniforme')
ax1.plot(-np.sin(radians(lon)), m_esfera(-np.sin(radians(lon)), opt_esfera),
			'-', label='Esfera Uniforme')
ax1.plot(-np.sin(radians(lon)), m_puntualdisco(-np.sin(radians(lon)),
			mdisco_opt, sdisco_opt), '-', label='Masa Puntual + Disco Uniforme')
ax1.plot(-np.sin(radians(lon)), m_puntualesfera(-np.sin(radians(lon)),
			mesfera_opt, pesfera_opt), '-', label='Masa Puntual + Esfera Uniforme')
ax1.set_xlabel('-R/R_0')
ax1.set_ylabel('Velocidad [km/sg]')
ax1.set_title('Velocidad asociada a Perfiles de Masa')
ax1.legend(loc=4)

fig2, ax2 = plt.subplots()
ax2.plot(-np.sin(radians(lon)), z(lon,b_ang), '-')
ax2.set_xlabel('-R/R_0')
ax2.set_ylabel('z [pc]')
ax2.set_title('Corrugacion Plano Galactico')

fig3,ax3 = plt.subplots()
ax3.plot(lon, Vmax, '-')
ax3.set_title('Velocidad Terminal segun longitud')
ax3.set_ylabel('Velocidad [km/sg]')
ax3.set_xlabel('Longitud [grados]')

fig4, ax4 = plt.subplots()
ax4.plot(-np.sin(radians(lon)), V(Vmax, lon), '*')
ax4.set_ylim(100,350)
ax4.set_ylabel('Velocidad [km/sg]')
ax4.set_xlabel('-R/R_0')
ax4.set_title('Curva de Rotacion')
