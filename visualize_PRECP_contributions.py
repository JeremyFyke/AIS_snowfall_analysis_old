import numpy as np
from scipy import spatial
from netCDF4 import MFDataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
from generate_region_mask import generate_region_mask

Variable="PRECT"
Regions=["LND","SOCN","SIO","SPO","SAO"]
Month="01"

class viewer_2d(object):
    def __init__(self,Total,Contributions,latv,lonv,lat,lon,RegionLongName,Variable,Month):

        self.Total=Total
        self.Contributions=Contributions
	self.lonv=lonv+0.625
	self.latv=latv
	self.lon=lon
	self.lat=lat	
	self.RegionLongName=RegionLongName
	self.Variable=Variable
	self.Month=Month
	
        self.fig=plt.figure()
        #Doing some layout with subplots:
        self.fig.subplots_adjust(0.05,0.05,0.98,0.98,0.1)
        self.overview=plt.subplot2grid((8,4),(0,0),rowspan=7,colspan=2) #main plot
        self.m = Basemap(projection='spstere',boundinglat=-60,lon_0=0,resolution='l',ax=self.overview)
	self.x, self.y = self.m(*np.meshgrid(self.lonv,self.latv))

	self.cs=self.m.pcolor(self.x,self.y,self.Total,vmin=0,vmax=2.5e-8)
        #self.cb=self.m.colorbar(self.cs,"bottom")
	self.overview
        self.m.drawcoastlines()
	# draw parallels.
        self.parallels = np.arange(-90,90,10.)
        self.m.drawparallels(self.parallels,labels=[1,0,0,0],fontsize=10)
        # draw meridians
        self.meridians = np.arange(0.,360.,10.)
        self.m.drawmeridians(self.meridians,labels=[0,0,0,1],fontsize=10)	

        self.pie_subplot=plt.subplot2grid((8,4),(0,2),rowspan=4,colspan=2) #evolving plot

        #Adding widgets, to not be gc'ed, they are put in a list:
        cursor=Cursor(self.overview, useblit=True, color='black', linewidth=2 )
        self._widgets=[cursor]
        #connect events
        self.fig.canvas.mpl_connect('button_press_event',self.click)

    def click(self,event):
        if event.inaxes==self.overview:
            #Get nearest data
	    i1=event.xdata
	    j1=event.ydata    
            plon,plat=self.m(i1,j1,inverse=True)
	    if plon < 0.:
	       plon=plon+360.
	    
	    i=np.argmin(np.absolute(self.lonv-plon))
	    j=np.argmin(np.absolute(self.latv-plat))
            print 'Lat/lon: ',plat,'/',plon

	    self.pie_subplot.clear()
	    patches, texts = self.pie_subplot.pie(self.Contributions[j,i,:],shadow=True,)
	    self.pie_subplot.legend(patches,self.RegionLongName,bbox_to_anchor=(1,0))
	    self.pie_subplot.set(ylabel=self.Variable+" contributors, Month "+self.Month,aspect='equal')
	    
	    plt.draw()

if __name__=='__main__':

    AtmSrcDir="/glade/scratch/hailong/archive/composite_ICE_wtag_mean/atm/hist/"

    #Get lon/lat grids
    lat,lon,latv,lonv,LandFrac,AllRegionNames,AllRegionMask=generate_region_mask()    

    f=MFDataset(AtmSrcDir+"composite_ICE_wtag_mean.cam.h0.*"+Month+".nc")
    PRECT_TOT=np.mean(np.squeeze(f.variables[Variable+"_"+"H2O"][:,:]),0)

    RegionLongName=[]
    for r in Regions: #Get long names for regions of interest (from list above)
         for k,(RgName,LongName) in enumerate(AllRegionNames):
	     if r==RgName:
	         RegionLongName.append(LongName)
    Regions.append("Res")
    RegionLongName.append("Residual")
    
    PRECT=np.zeros(( np.shape(lat)[0], np.shape(lat)[1], len(Regions) ))  
    for k in np.arange(0,len(Regions)-1): #omit residual field, which is calculated next
        vname=Variable+"_"+Regions[k]
	PRECT[:,:,k]=np.mean(np.squeeze(f.variables[vname][:,:]),0)
    f.close()
    #calculate residual
    PRECT[:,:,-1]=PRECT_TOT-np.sum(PRECT[:,:,0:-2],2)
    for k in np.arange(0,len(Regions)):
        PRECT[:,:,k]=PRECT[:,:,k]/PRECT_TOT
    fig_v=viewer_2d(PRECT_TOT,PRECT,latv,lonv,lat,lon,RegionLongName,Variable,Month)

    plt.show()
