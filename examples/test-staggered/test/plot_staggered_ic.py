from boutdata.collect import collect
import matplotlib.pyplot as plt

n  = collect("n")
v  = collect("v")

plt.plot(n[0,2,:,0])
plt.plot(v[0,2,:,0])
plt.show()

