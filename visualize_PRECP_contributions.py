import numpy as np
from scipy import spatial
from netCDF4 import MFDataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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
        self.m.drawcoastlines()
	# draw parallels.
        self.parallels = np.arange(-90,90,10.)
        self.m.drawparallels(self.parallels,labels=[1,0,0,0],fontsize=10)
        self.meridians = np.arange(0.,360.,10.)
        self.m.drawmeridians(self.meridians,labels=[0,0,0,1],fontsize=10)	

        self.interactive_subplot=plt.subplot2grid((8,4),(0,2),rowspan=4,colspan=2) #evolving plot

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

	    self.interactive_subplot.clear()
	    #patches, texts = self.interactive_subplot.pie(self.Contributions[j,i,:],shadow=True,)
	    x=np.arange(0,12)
	    y=np.transpose(np.squeeze(self.Contributions[j,i,:,:]))
	    stack_coll =self.interactive_subplot.stackplot(x,y)
	    proxy_rects = [Rectangle((0, 0), 1, 1, fc=pc.get_facecolor()[0]) for pc in stack_coll]
	    self.interactive_subplot.legend(proxy_rects,RegionLongName,bbox_to_anchor=(0.9,-.25))
	    self.interactive_subplot.get_yaxis().set_visible(False)

	    #self.interactive_subplot.legend(r,self.RegionLongName,bbox_to_anchor=(1,0))
	    self.interactive_subplot.axis('tight')
	    #self.interactive_subplot.legend(patches,self.RegionLongName,bbox_to_anchor=(1,0))
	    self.interactive_subplot.set(xlabel=self.Variable+" contributions by climo. month ")
	    
	    plt.draw()

if __name__=='__main__':
    
    #Get lon/lat grids
    lat,lon,latv,lonv,LandFrac,AllRegionNames,AllRegionMask=generate_region_mask()    

    AtmSrcDir="/glade/scratch/hailong/archive/composite_ICE_wtag_mean/atm/hist/"
    RegionLongName=[]
    for r in Regions: #Get long names for regions of interest (from list above)
         for k,(RgName,LongName) in enumerate(AllRegionNames):
	     if r==RgName:
	         RegionLongName.append(LongName)
    Regions.append("Res")
    RegionLongName.append("Residual")
    
    PRECT_TOT=np.zeros(( np.shape(lat)[0], np.shape(lat)[1], 12 )) 
    PRECT=np.zeros(( np.shape(lat)[0], np.shape(lat)[1], 12, len(Regions) )) 
    for m in np.arange(0,12):
	Month='%02d'%(m+1)
	print 'Building climatology for month: '+Month+'...'
	f=MFDataset(AtmSrcDir+"composite_ICE_wtag_mean.cam.h0.*"+Month+".nc")
	PRECT_TOT[:,:,m]=np.mean(np.squeeze(f.variables[Variable+"_"+"H2O"][:,:,:]),0)
	for k in np.arange(0,len(Regions)-1): #omit residual field, which is calculated next
            vname=Variable+"_"+Regions[k]
	    PRECT[:,:,m,k]=np.mean(np.squeeze(f.variables[vname][:,:,:]),0)
	f.close()
	#calculate residual
	PRECT[:,:,m,-1]=PRECT_TOT[:,:,m]-np.sum(np.squeeze(PRECT[:,:,m,::-2]),axis=2) #Note, including residual field here since it's zero
	#for k in np.arange(0,len(Regions)):
        #    PRECT[:,:,m,k]=PRECT[:,:,m,k]/PRECT_TOT[:,:,m]
    
    print 'Constructing viewer...'    
    fig_v=viewer_2d(np.mean(PRECT_TOT,2),PRECT,latv,lonv,lat,lon,RegionLongName,Variable,Month)

    plt.show()
