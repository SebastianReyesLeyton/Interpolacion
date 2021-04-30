import matplotlib.pyplot as plt

times = [0.00010180473327636719, 8.153915405273438e-05, 0.0001583099365234375, 0.0001914501190185547]
labels = ['Polinomial', 'Lagrange', 'Newton', 'Lineal a trozos']

plt.figure("Pronostico del tiempo")
plt.title("Pronostico del tiempo")
plt.grid()
plt.bar(labels, times, color=['blue', 'green', 'gray', 'red'])
plt.ylabel('Tiempo')
plt.xlabel('Interpolación')
plt.show()