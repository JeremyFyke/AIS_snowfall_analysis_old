import numpy as np
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import Cursor, Button
from generate_region_mask import generate_region_mask
from math import radians, sin, cos, asin, sqrt, pi, atan2

###User Input Here###
Variable="PRECT"
Regions=["LND","SOCN","SIO","SPO","SAO"]
CaseNamePrefix="composite_ICE_wtag_"
AtmSrcDirPrefix="/glade/scratch/hailong/amwg/climo/composite_ICE_wtag_"
MaskFile="/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/data/AIS_Full_basins_Zwally_CESMgrid.nc"
###User Input Ends###

class viewer_2d(object):
    def __init__(self,Mask,Contributions,dContributions_normalized,latv,lonv,lat,lon,RegionLongName,Variable,Month):

        self.Mask=Mask
        self.Contributions=Contributions
	self.dContributions_normalized=dContributions_normalized
	self.lonv=lonv+0.625
	self.latv=latv
	self.lon=lon
	self.lat=lat	
	self.RegionLongName=RegionLongName
	self.Variable=Variable
	self.Month=Month

        self.fig=plt.figure(figsize=(10,10))
	self.fig.set_size_inches(16,8)
        self.fig.subplots_adjust(0.05,0.05,0.98,0.98,0.2)
	
        #Build main plot	
        self.overview=plt.subplot2grid((15,4),(7,0),rowspan=5,colspan=2) #main plot
        self.m = Basemap(projection='spstere',boundinglat=-60,lon_0=0,resolution='l',ax=self.overview)
	self.x, self.y = self.m(*np.meshgrid(self.lonv,self.latv))
	self.cs=self.m.pcolor(self.x,self.y,self.Mask)
        #self.cb=self.m.colorbar(self.cs,"bottom")
        self.m.drawcoastlines()
	# draw parallels.
        self.parallels = np.arange(-90,90,10.)
        self.m.drawparallels(self.parallels,labels=[1,0,0,0],fontsize=10)
        self.meridians = np.arange(0.,360.,10.)
        self.m.drawmeridians(self.meridians,labels=[0,0,0,1],fontsize=10)	

        #Initialize axis for evolving plot      
        self.interactive_subplot1=plt.subplot2grid((15,4),(2,2),rowspan=5,colspan=2) #evolving plot 1
        self.interactive_subplot2=plt.subplot2grid((15,4),(8,2),rowspan=5,colspan=2) #evolving plot 1	

        #Initialize interactivity of evolving plot via an active cursor over the overview plot
        cursor=Cursor(self.overview, useblit=True, color='black', linewidth=1 )
        self._widgets=[cursor]
        #connect events
        self.fig.canvas.mpl_connect('button_press_event',self.click)

    #Define procedure that occurs when a mouse click occurs on the overview plot.
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
            print '*********'
	    print 'Lat/lon: ',plat,'/',plon
	    print 'Plotting basin #: '+str(self.Mask[j,i])
	    
	    MaskValue=np.round(self.Mask[j,i]) 
	    MaskValue=MaskValue.astype(int)-1#convert to integer and reduce by one since this is used to index into the Contributions array

            if MaskValue < 0.:
	        print 'No region defined here, not replotting.'
	    else: 
	        #Replot contributions plot for mean-SIC case, for clicked basin
		self.interactive_subplot1.clear()
		x=np.arange(0,12) #vector of values representing months
		y=np.transpose(np.squeeze(self.Contributions[MaskValue,:,:,1]))
		stack_coll =self.interactive_subplot1.stackplot(x,y)
		proxy_rects = [Rectangle((0, 0), 1, 1, fc=pc.get_facecolor()[0]) for pc in stack_coll]
		self.interactive_subplot1.legend(proxy_rects,RegionLongName,bbox_to_anchor=(-0.25,.7))	
		
		self.interactive_subplot1.axis('tight')
		self.interactive_subplot1.set(xlabel="Month",
		                             ylabel='Gt/yr flux')
					     
		#Replot fractional change in contributions plot, for clicked basin
		self.interactive_subplot2.clear()		     
		x=np.arange(0,12) #vector of values representing months
		y=np.squeeze(self.dContributions_normalized[MaskValue,:,:])			     
                self.interactive_subplot2.plot(x,y)
		self.interactive_subplot2.axis('tight')
		self.interactive_subplot2.set(xlabel="Month",
		                             ylabel=self.Variable+' change in flux (change in %, high-low SIC)')
		
		plt.draw()

if __name__=='__main__':
    
    #Get lon/lat grids
    lat,lon,latv,lonv,Area,LandFrac,AllRegionNames,AllRegionMask=generate_region_mask()    
    
    SecPYr=3.154e7
    KgPM3=1000.
    MmPM=1000.
    GtPKg=1.e12
    
    #Get long names for regions of interest (from list above)
    RegionLongName=[]
    for r in Regions:
         for k,(RgName,LongName) in enumerate(AllRegionNames):
	     if r==RgName:
	         RegionLongName.append(LongName)
    Regions.append("Res")
    RegionLongName.append("Residual")

    f=nc.Dataset(MaskFile)
    IMBIEBasinMask=f.variables["Zwallybasins"][:,:]
    IMBIEBasinMask=np.round(IMBIEBasinMask)
    nIMBIEBasins=np.amax(IMBIEBasinMask)
    nIMBIEBasins=nIMBIEBasins.astype(int)
    f.close()
      
    def get_shortest_in(needleLat,needleLon, haystackLat, haystackLon):
	"""
            haystack is a numpy array to find the point in
            that has the shortest distance to needle
	"""
	radius_km=6371
	haystackLat=np.radians(haystackLat)
	needleLat=np.radians(needleLat)
	haystackLon=np.radians(haystackLon)
	needleLon=np.radians(needleLon)
	dlat = haystackLat - needleLat
	dlon = haystackLon - needleLon
	a = np.square(np.sin(dlat/2.0)) + np.cos(needleLat) * np.cos(haystackLat) * np.square(np.sin(dlon/2.0))
	great_circle_distance = 2 * np.arcsin(np.minimum(np.sqrt(a), np.repeat(1, len(a))))
	d = radius_km * great_circle_distance
	return np.min(d)    
    
    #Generate regions that are a function of distance from coast.
    CoastMask=np.zeros((np.size(Area,0),np.size(Area,1)))
    ContinentalMask=np.zeros((np.size(Area,0),np.size(Area,1)))
    CoastMask[IMBIEBasinMask>0.]=1
    CoastLon=[]
    CoastLat=[]
    #Make lat/lon list coastal ocean points around AIS
    for nlon in np.arange(0,np.size(Area,1)):
	nlat=0
	CoastVal=1
	while CoastVal==1:
	    nlat=nlat+1
	    CoastVal=CoastMask[nlat,nlon]
	CoastLon.append(lon[nlat,nlon])
        CoastLat.append(lat[nlat,nlon])

    BinInterval=200.;
    for nlat in np.arange(0,np.size(Area,0)):
        for nlon in np.arange(0,np.size(Area,1)):
	    if CoastMask[nlat,nlon]==1:
	        ProximityToCoast=get_shortest_in( lat[nlat,nlon],lon[nlat,nlon] , CoastLat,CoastLon  )
		nbin=np.ceil(ProximityToCoast/BinInterval)
		ContinentalMask[nlat,nlon]=nbin
    nContinentalBasins=np.amax(ContinentalMask)
    nContinentalBasins=nContinentalBasins.astype(int)
    
    #Initialize total and contribution fields
    
    FIELD=np.zeros(( np.size(Area,0),np.size(Area,1), 12, len(Regions), 3 ))
    FIELD_TOT=np.zeros(( np.size(Area,0),np.size(Area,1), 12, 3 )) 
    FIELD_normalized=np.zeros(( nIMBIEBasins, 12, len(Regions), 3 ))
    
    BASIN_FIELD=np.zeros(( nIMBIEBasins, 12, len(Regions), 3 ))
    BASIN_FIELD_TOT=np.zeros(( nIMBIEBasins, 12, 3 ))    
    BASIN_FIELD_normalized=np.zeros(( nIMBIEBasins, 12, len(Regions), 3 ))    
    
    #TOT=np.zeros(( 12, 3 )) 

    MaskedArea=Area*LandFrac
    
    for nSIC,SIC in enumerate(["low","mean","high"]):
    
        print 'Processing '+SIC+' simulation...'
	
	TMP_FIELD_TOT_ANN=np.zeros(np.shape(Area))
	
	#Load original data
	for nMonth in np.arange(0,12):
	    Month='%02d'%(nMonth+1)
	    CaseName=CaseNamePrefix+SIC
	    AtmSrcDir=AtmSrcDirPrefix+SIC
	    f=nc.Dataset(AtmSrcDir+'/'+CaseName+"_"+Month+"_climo.nc")
	    FIELD_TOT[:,:,nMonth,nSIC]=np.squeeze(f.variables[Variable+"_"+"H2O"][:,:,:])*KgPM3*SecPYr/GtPKg #Load total monthly mean flux field
            for nTag in np.arange(0,len(Regions)-1): #Note: omit residual field, which is calculated next
        	vname=Variable+"_"+Regions[nTag]
		FIELD[:,:,nMonth,nTag,nSIC]=np.squeeze(f.variables[vname][:,:,:])*KgPM3*SecPYr/GtPKg #Load tagged fields
	    f.close()
		
	    #Calculate residual
            FIELD[:,:,nMonth,-1,nSIC]=FIELD_TOT[:,:,nMonth,nSIC]-np.sum(FIELD[:,:,nMonth,:,nSIC],axis=2)
            
	    #Calculate basin-specific values
            for nBasin in np.arange(0,nContinentalBasins):
	        iMask=np.where(ContinentalMask==(nBasin+1))
		TMP=np.squeeze(FIELD_TOT[:,:,nMonth,nSIC])
	        BASIN_FIELD_TOT[nBasin,nMonth,nSIC]=np.sum(TMP[iMask]*MaskedArea[iMask])
		for nTag in np.arange(0,len(Regions)): #calculate monthly integrated flux into each basin, from each tagged source region
		    TMP=np.squeeze(FIELD[:,:,nMonth,nTag,nSIC])
		    BASIN_FIELD[nBasin,nMonth,nTag,nSIC]=np.sum(TMP[iMask]*MaskedArea[iMask])
		    BASIN_FIELD_normalized[nBasin,nMonth,nTag,nSIC]=BASIN_FIELD[nBasin,nMonth,nTag,nSIC]/BASIN_FIELD_TOT[nBasin,nMonth,nSIC]		

        TMP_FIELD_TOT=np.mean(BASIN_FIELD_TOT[:,:,nSIC],axis=1)
	n=np.sum(TMP_FIELD_TOT)
        print 'total SMB for '+SIC+' simulation='+str(n)
	
    #Plot percent change in fractional contribution, for each basin, for each month, for high minus low SIC

    dFIELD_normalized=(BASIN_FIELD_normalized[:,:,:,2]-BASIN_FIELD_normalized[:,:,:,0])
	
    #Call viewer class to make plot.
    #print 'Constructing viewer...' 
    fig_v=viewer_2d(ContinentalMask,BASIN_FIELD,dFIELD_normalized,latv,lonv,lat,lon,RegionLongName,Variable,Month)

    plt.show()
