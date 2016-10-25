import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from generate_region_mask import generate_region_mask

lat,lon,latv,lonv,LandFrac,RegionNames,RegionMask=generate_region_mask()

#Make a hover plot in plotly, over AIS, that shows breakdown of contributions to point
#Color plot by region of most significance

fig = plt.figure()
ax = plt.subplot(111)
plt.subplots_adjust(left=0.25, bottom=0.25)
s=np.flipud(RegionMask[:,:,1])
l=ax.imshow(s)
rn=RegionNames[1]
plt.title(rn[1])
c=ax.contour(np.flipud(LandFrac),levels=[0,0],colors='white')

axcolor = 'lightgoldenrodyellow'
sliderax = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
slider = Slider(sliderax, 'Mask', 1, len(RegionNames), valinit=0, valfmt='%i')

def update(val):
    k=int(slider.val)
    rn=RegionNames[k]
    l.set_data(np.flipud(RegionMask[:,:,k]))
    ax.set_title(rn[1])
    fig.canvas.draw()
slider.on_changed(update)

#resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
#button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
#def reset(event):
#    slider.reset()
#button.on_clicked(reset)

plt.show()
