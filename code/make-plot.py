# python plotting of diffusion data to avoid julia issues

#want to make a contour plot 

import matplotlib.pyplot as plt
import numpy as np
import os
# plt.style.use('_mpl-gallery-nogrid')
dst = os.path.join(os.getcwd(),"pipeline/tmp/")
# load data

import csv
agg_data = np.loadtxt("pipeline/imgs/agg-data-rewired-vs-er-modmexico-city-seir-0.1-0.05-uniform.txt",delimiter=',')

agg_data.shape


levels = np.linspace(0, 1, 10)

# plot
fig, ax = plt.subplots()
cf = ax.contourf(np.linspace(1, 26, 25),np.linspace(0, 16, 16),agg_data, levels=11,
            vmin=0,vmax=1,cmap = "hot")
# ax.contourf(agg_data, levels=levels)
ax.set_xlabel("Powerlaw     <-    Original   ->    ER \nRewiring Percent")
ax.set_ylabel("Quarantine Percent")
cbar = fig.colorbar(cf,ticks=np.linspace(0, 1, 6),label="Fraction of Infected Nodes")


fig.savefig(dst+'line_plot.png', bbox_inches='tight')


