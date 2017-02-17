import numpy as np
from scipy import spatial
import netCDF4 as nc
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import Cursor, Button
from generate_region_mask import generate_region_mask

###User Input Here###
Variable="PRECT"
Regions=["LND","SOCN","SIO","SPO","SAO"]
CaseNamePrefix="composite_ICE_wtag_"
AtmSrcDirPrefix="/glade/scratch/hailong/amwg/climo/composite_ICE_wtag_"
BasinFile="/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/data/AIS_Full_basins_Zwally_CESMgrid.nc"
BasinVariable="Zwallybasins"
###User Input Ends###

#Get lon/lat grids
lat,lon,latv,lonv,Area,LandFrac,AllRegionNames,AllRegionMask=generate_region_mask()
 
#Get long names for regions of interest (from list above)
RegionLongName=[]
for r in Regions:
     for k,(RgName,LongName) in enumerate(AllRegionNames):
	 if r==RgName:
	     RegionLongName.append(LongName)
Regions.append("Res")
RegionLongName.append("Residual")
    
SecPYr=3.154e7
KgPM3=1000.
MmPM=1000.
GtPKg=1.e12

FIELD=np.zeros(( np.size(Area,0),np.size(Area,1), 12, len(Regions), 3 ))
FIELD_TOT=np.zeros(( np.size(Area,0),np.size(Area,1), 12, 3 )) 

for nSIC,SIC in enumerate(["low","mean","high"]):
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
        
	#Normalize values
	for nTag in np.arange(0,len(Regions)):
	    FIELD[:,:,nMonth,nTag,nSIC]=FIELD[:,:,nMonth,nTag,nSIC]/FIELD_TOT[:,:,nMonth,nSIC]

FIELD_ANN_AVG=np.mean(FIELD,axis=2)


for nTag,Tag in enumerate(Regions):
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0.25, bottom=0.25)
    m = Basemap(projection='spstere',boundinglat=-60,lon_0=0,resolution='l')
    x, y = m(*np.meshgrid(lonv,latv))
    TMP=FIELD_ANN_AVG[:,:,nTag,1]
    TMP[LandFrac==0.]=0.
    l=m.pcolor(x,y,FIELD_ANN_AVG[:,:,nTag,1],ax=ax,vmin=0,vmax=0.5)

    cb=m.colorbar(l,"right")
    m.drawcoastlines()
    # draw parallels.
    parallels = np.arange(-90,90,10.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    meridians = np.arange(0.,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

    plt.title(Tag+'_PRECT_contributions')

    plt.savefig('figs/'+Tag+'_PRECT_contributions.png')







