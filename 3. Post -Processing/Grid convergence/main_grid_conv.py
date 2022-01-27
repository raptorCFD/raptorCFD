import matplotlib.pyplot as plt

name = 'Grid Convergence'

elements = [131925, 242866, 354788, 459479, 564604, 652517]
cd = [0.014174, 0.013063, 0.012853, 0.012809, 0.012774, 0.012755]

plt.plot(elements, cd, 'b', linewidth=2)
plt.suptitle(name, size = 20)
plt.ylabel('Drag Coefficient []', size=10)
plt.xlabel('Number of elements []', size= 10)
labels = [cd[i] for i in range(0, len(elements))]
for i in range(0, len(elements)):
    plt.text(elements[i], cd[i], labels[i], va='bottom')
plt.grid(True)
plt.ylim((0.0125, 0.0145))

plt.savefig('grid_convergence.png')
plt.savefig('grid_convergence.pdf')
plt.close()

