import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
import numpy as np
import matplotlib

x, y = np.random.rand(2, 100)

fig, ax = plt.subplots()

p, = plt.plot(x, y, 'o')

cursor = Cursor(ax,
                horizOn=True,
                vertOn=True,
                color='green',
                linewidth=2.0)

def onclick(event):
    x1, y1 = event.xdata, event.ydata
    print(x1, y1)

fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()