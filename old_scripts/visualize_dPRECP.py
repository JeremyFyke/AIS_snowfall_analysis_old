import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.widgets import Slider
from generate_region_mask import generate_region_mask

#Make monthly spatial plots of dPRECP/dSIC

###User Input Here###
Variable="PRECT"
CaseNamePrefix="composite_ICE_wtag_"
AtmSrcDirPrefix="/glade/scratch/hailong/amwg/climo/composite_ICE_wtag_"
BasinFile="/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/data/AIS_Full_basins_Zwally_CESMgrid.nc"
BasinVariable="Zwallybasins"
###User Input Ends###
lat,lon,latv,lonv,Area,LandFrac,RegionNames,RegionMask=generate_region_mask()
SecPYr=3.154e7
KgPM3=1000.
MmPM=1000.
GtPKg=1.e12

FIELD_TOT=np.zeros(( np.size(Area,0),np.size(Area,1), 12, 2 ))

for nSIC,SIC in enumerate(["low","high"]):

    print 'Processing '+SIC+' simulation...'

    #Load original data
    for nMonth in np.arange(0,12):
	Month='%02d'%(nMonth+1)
	CaseName=CaseNamePrefix+SIC
	AtmSrcDir=AtmSrcDirPrefix+SIC
	f=nc.Dataset(AtmSrcDir+'/'+CaseName+"_"+Month+"_climo.nc")
	FIELD_TOT[:,:,nMonth,nSIC]=np.squeeze(f.variables[Variable+"_"+"H2O"][:,:,:])*KgPM3*SecPYr/GtPKg #Load total monthly mean flux field

dFIELD_TOT=(FIELD_TOT[:,:,:,1]-FIELD_TOT[:,:,:,0])

vmin=-1.e-10
vmax=1.e-10
fig = plt.figure()
ax = plt.subplot(111)
plt.subplots_adjust(left=0.25, bottom=0.25)
m = Basemap(projection='spstere',boundinglat=-60,lon_0=0,resolution='l')
x, y = m(*np.meshgrid(lonv,latv))
l=m.pcolor(x,y,dFIELD_TOT[:,:,0],ax=ax,vmin=vmin,vmax=vmax)

cb=m.colorbar(l,"bottom")
m.drawcoastlines()
# draw parallels.
parallels = np.arange(-90,90,10.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
meridians = np.arange(0.,360.,10.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

axcolor = 'lightgoldenrodyellow'
sliderax = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
slider = Slider(sliderax, 'Month', 1, 12, valinit=1, valfmt='%i')

def update(val):
    k=int(slider.val)-1
    ax.clear()
    l=m.pcolor(x,y,dFIELD_TOT[:,:,k],ax=ax,vmin=vmin,vmax=vmax)
    m.drawcoastlines(ax=ax)
    # draw parallels.
    parallels = np.arange(-90,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    meridians = np.arange(0.,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    fig.canvas.draw_idle()
slider.on_changed(update)

plt.show()
