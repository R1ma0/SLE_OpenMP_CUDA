import matplotlib.pyplot as plt

time = [0.013060, 0.126935, 2.241650, 27.705708]
size = [1000, 3000, 10000, 20000]

plt.plot(size, time)
plt.xlabel("N")
plt.ylabel("sec")
plt.show()
