import matplotlib.pyplot as plt
import numpy as np

vx = 41.922409349774263
vy = 27.653668647301799
wx = 40.397759903758043
wy = 20.000000000000004

p0x = 40.406241130654529
p0y = 20.03233309801665
p1x = 38.844130511042778
p1y = 27.811769362379923

plt.plot([p0x], [p0y], marker=".", color="k")
plt.plot([p0x, p1x], [p0y, p1y], color="k")
plt.plot([wx], [wy], marker=".", color="r")
plt.plot([wx, vx], [wy, vy], color="r")


