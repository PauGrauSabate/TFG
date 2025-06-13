import numpy as np
import matplotlib.pyplot as plt
import os

def leer_posiciones_xy(archivo):
    datos = np.loadtxt(archivo, usecols=(0, 1))  # solo x, y
    return datos

def distancia_minima_imagen(dr, Lx, Ly):
    dr[:, 0] = dr[:, 0] - Lx * np.round(dr[:, 0] / Lx)
    dr[:, 1] = dr[:, 1] - Ly * np.round(dr[:, 1] / Ly)
    return dr

def calcular_rdf(posiciones, rho = 0.05, dr=0.1):
    N = len(posiciones)
    Lx = np.sqrt(N/rho)
    Ly = Lx
    r_max = Lx / 2

    bins = np.arange(0, r_max + dr, dr)
    histograma = np.zeros(len(bins) - 1)

    for i in range(N):
        dr_vectores = posiciones - posiciones[i]
        dr_vectores = distancia_minima_imagen(dr_vectores, Lx, Ly)
        distancias = np.sqrt(np.sum(dr_vectores**2, axis=1))
        distancias = distancias[distancias > 0]  # excluir i == j
        hist, _ = np.histogram(distancias, bins=bins)
        histograma += hist

    # Normalización
    r = 0.5 * (bins[1:] + bins[:-1])
    area_anillos = np.pi * ((bins[1:]**2) - (bins[:-1]**2))
    norm = rho * N * area_anillos
    g_r = histograma / norm

    return r, g_r

def main():
    directorio = "/home/pgrau/h1/Lambda10"

    inicio = 10900
    fin = 11000
    step = 1

    archivos = [os.path.join(directorio, f"conf{t}.dat") for t in range(inicio, fin + 1, step)]
    rdf_total = []

    for archivo in archivos:
        posiciones = leer_posiciones_xy(archivo)
        r, g_r = calcular_rdf(posiciones)
        rdf_total.append(g_r)

    g_r_promedio = np.mean(rdf_total, axis = 0)

    np.savetxt("g_r_promedio_h1_L10.dat", np.column_stack((r, g_r_promedio)))

    plt.plot(r, g_r_promedio)
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.title("Funció de distribució radial")
    plt.grid()
    plt.show()
    plt.savefig("FDR_h1_L10.png")

main()
