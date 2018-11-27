import matplotlib.pyplot as plt
import numpy as np

data = np.array([ list(map(float,l.split())) for l in open("reldat") if l and not l.startswith("#")])
print(data)

plt.figure(figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(data[:,0], data[:,1], "+")
plt.gca().set_title("Current time delay fitness measure")
plt.gca().set_xlabel("Scale factor")

plt.subplot(2,2,2)
plt.plot(data[:,0], data[:,2], "+")
plt.gca().set_title("Experimental time delay fitness measure II ")
plt.gca().set_xlabel("Scale factor")

plt.subplot(2,2,3)
plt.gca().set_yscale("log")
plt.plot(data[:,0], data[:,3], "+")
plt.gca().set_title("Current time delay fitness measure rel (log)")
plt.gca().set_xlabel("Scale factor")

plt.subplot(2,2,4)
plt.plot(data[:,0], data[:,4], "+")
plt.gca().set_title("Experimental time delay fitness measure II rel")
plt.gca().set_xlabel("Scale factor")

plt.savefig("rel_fitness_comparison.png", bbox_inches="tight")

