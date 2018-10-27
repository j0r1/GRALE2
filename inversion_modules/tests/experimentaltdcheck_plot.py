import matplotlib.pyplot as plt
import numpy as np

data = np.array([ list(map(float,l.split())) for l in open("dat") if l and not l.startswith("#")])
print(data)

plt.figure(figsize=(20,5))
plt.subplot(1,4,1)
plt.plot(data[:,0], data[:,1], "+")
plt.gca().set_title("Current time delay fitness measure")
plt.gca().set_xlabel("Scale factor")

plt.subplot(1,4,2)
plt.gca().set_yscale("log")
plt.plot(data[:,0], data[:,1], "+")
plt.gca().set_title("Current time delay fitness measure (log)")
plt.gca().set_xlabel("Scale factor")

plt.subplot(1,4,3)
plt.plot(data[:,0], data[:,2], "+")
plt.gca().set_title("Experimental time delay fitness measure")
plt.gca().set_xlabel("Scale factor")

plt.subplot(1,4,4)
plt.gca().set_yscale("log")
plt.plot(data[:,0], data[:,2], "+")
plt.gca().set_title("Experimental time delay fitness measure (log)")
plt.gca().set_xlabel("Scale factor")

plt.savefig("td_fitness_comparison.png", bbox_inches="tight")

