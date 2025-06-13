import numpy as np
import os
import re
import matplotlib.pyplot as plt

def leer_configuracion(archivo, indices_interes):
    """Lee las posiciones x, y de las partículas de interés en un archivo."""
    with open(archivo, 'r') as f:
        lineas = f.readlines()
    
    posiciones = []
    for i in indices_interes:
        columnas = lineas[i].split()
        x, y = float(columnas[0]), float(columnas[1])  # solo x, y
        posiciones.append([x, y])
    
    return np.array(posiciones)

def obtener_archivos_configuracion(directorio, inicio=10500):
    """Obtiene y ordena los archivos confXXXXX.dat por su número."""
    archivos = [f for f in os.listdir(directorio)
         if re.match(r"conf\d+\.dat", f) and inicio <= int(re.findall(r'\d+', f)[0])]
    archivos.sort(key=lambda x: int(re.findall(r'\d+', x)[0]))
    return [os.path.join(directorio, f) for f in archivos]

def reconstruir_trayectorias(posiciones_por_tiempo, Lx, Ly):
    """Reconstruye la trayectoria desenrollada de cada partícula."""
    unwrap = [posiciones_por_tiempo[0]]  # primera configuración sin cambio

    for i in range(1, len(posiciones_por_tiempo)):
        r_prev = unwrap[-1]
        r_wrapped_prev = posiciones_por_tiempo[i - 1]
        r_wrapped_now = posiciones_por_tiempo[i]
        
        delta = r_wrapped_now - r_wrapped_prev
        # aplicar mínima imagen
        delta[:, 0] = (delta[:, 0] + Lx / 2) % Lx - Lx / 2
        delta[:, 1] = (delta[:, 1] + Ly / 2) % Ly - Ly / 2

        r_unwrapped = r_prev + delta
        unwrap.append(r_unwrapped)
    
    return unwrap

def calcular_msd(trayectorias_unwrapped):
    """Calcula el MSD a partir de trayectorias ya desenrolladas."""
    n_pasos = len(trayectorias_unwrapped)
    msd_por_delta_t = []

    for delta_t in range(1, n_pasos):
        msds = []
        for t0 in range(n_pasos - delta_t):
            r0 = trayectorias_unwrapped[t0]
            rt = trayectorias_unwrapped[t0 + delta_t]
            dr = rt - r0
            dr2 = np.sum(dr**2, axis=1)
            msd = np.mean(dr2)
            msds.append(msd)
        
        msd_promedio = np.mean(msds)
        msd_por_delta_t.append(msd_promedio)
    
    return msd_por_delta_t

def main():
    directorio = "/home/pgrau/h1/Lambda10"
    indices_particulas_interes = [195,196,197,198,199,200,201,202,203,204]  # cambiar según necesidad

    Lx = np.sqrt(400/0.1)  # tamaño de la caja en x
    Ly = Lx  # tamaño de la caja en y
    
    # Leer posición inicial desde pos10000.dat
    archivo_inicial = os.path.join(directorio, "pos10000.dat")
    posiciones_iniciales = leer_configuracion(archivo_inicial, indices_particulas_interes)

    # Leer configuraciones posteriores
    archivos = obtener_archivos_configuracion(directorio)
    posiciones_envueltas = [leer_configuracion(f, indices_particulas_interes) for f in archivos]

    # Insertar la configuración inicial al principio
    posiciones_envueltas.insert(0, posiciones_iniciales)

    # Reconstruir trayectoria desenrollada
    posiciones_desenrolladas = reconstruir_trayectorias(posiciones_envueltas, Lx, Ly)

    # Calcular MSD
    msd = calcular_msd(posiciones_desenrolladas)
    deltas = list(range(1, len(msd) + 1))
    # Mostrar resultados
    # for delta_t, valor in enumerate(msd, start=1):
        # print(f"Δt = {delta_t}: MSD = {valor:.4f}")
    
    return deltas, msd
        
dt, msd = main()

x = np.array(dt)
y = np.array(msd)

filtre = x <= 150
x_reg = x[filtre]
y_reg = y[filtre]
coef, cov = np.polyfit(x_reg,y_reg,deg=1,cov = True)
m, n = coef
dm = np.sqrt(cov[0,0])
dn = np.sqrt(cov[1,1])

D = m/4
print(f"Coeficient de difusio D = {D:.5f}")

y_pred = m*x_reg+n
y_mean = np.mean(y_reg)

ss_res = np.sum((y_reg - y_pred)**2)
ss_tot = np.sum((y_reg - y_mean)**2)

r2 = 1 - ss_res/ss_tot

plt.scatter(x_reg, y_reg)
plt.plot(x_reg,y_pred,color = 'red',label = f"$y = {m:.2f}x + {n:.2f}$, $R^2 = {r2:.4f}$") 
plt.title('MSD')
plt.grid()
#plt.legend()
plt.xlabel('$\Delta$t')
plt.ylabel('MSD')
plt.savefig('MSD_h1_L10.png')
print("dm = ",dm)
print("dn = ",dn)
