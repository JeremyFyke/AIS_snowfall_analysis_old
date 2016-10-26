def generate_region_mask(AtmSrcDir):

    import numpy as np
    import netCDF4 as nc

    print 'Generating region masks...'

    #Define vapor tagging regions.  Fields in each array are:
    #1) short name
    #2) long name
    #3) lat/lon bounds
    #4) mask
    #%
    RegionNames=[
       ("H2O","the globe"),
       ("ARO","Arctic Ocean"),
       ("ARL","Arctic Land"),
       ("NPO","N. Pacific"),
       ("NAO","N. Atlantic"),
       ("SNP","Subtropical N.P."),
       ("GOM","Gulf of Mexico"),
       ("SNA","Subtropical N.A."),
       ("NIO","N. Indian Ocean"),
       ("BAB","a subset of the Southern Ocean (East)"),
       ("DAS","Davis Strait and Baffin Bay"),
       ("LAS","Labrador Sea"),
       ("GLS","Greenland Sea"),
       ("ILS","Iceland Sea"),
       ("NGS","a subset of the Southern Ocean (West)"),
       ("LND","Land"),
       ("PWP","Pacific Warm Pool"),
       ("EPO","Equatorial Pacific"),
       ("EAO","Equatorial Atlantic"),
       ("SIO","S. Indian Ocean"),
       ("SPO","S. Pacific"),
       ("SAO","S. Atlantic"),
       ("SOCN","S. Ocean"),
       ("AMR","Amundsen"),
       ("DML","DML Region"),
       ("WLR","Wilkes Land"),
       ]

    #Grab landmask from arbitrary file
    f=nc.Dataset(AtmSrcDir+"composite_ICE_wtag_mean.cam.h0.0001-01.nc")
    LandFrac=np.squeeze(f.variables["LANDFRAC"][:,:])
    f.close()

    #Build lat/lon grids
    f=nc.Dataset(AtmSrcDir+"composite_ICE_wtag_mean.cam.h1.0011-12-07-00000.nc")
    latv=f.variables["lat"][:]
    lonv=f.variables["lon"][:]
    lat=np.transpose(np.tile(latv,(len(lonv),1)))
    lon=np.tile(lonv,(len(latv),1))
    f.close()

    RegionMask=np.zeros(( np.shape(lat)[0], np.shape(lat)[1], len(RegionNames) ))
    for k,(RgName,LongName) in enumerate(RegionNames):

	if RgName == "H2O":
            RegionMask[:,:,k] = 1.

	if RgName == "LND":
            RegionMask[:,:,k] = LandFrac[:,:]

	for i in np.arange(0,np.shape(lat)[0]):

	    for j in np.arange(0,np.shape(lat)[1]):

        		if RgName == "ARO": 
        		  if lat[i,j] >= 65.: 
		              RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "ARL": 
        		  if lat[i,j] >= 65.: 
		              RegionMask[i,j,k] = LandFrac[i,j]

        		if RgName == "NPO": 
        		  if lat[i,j] > 30. and lat[i,j] < 65.: 
        		    if lon[i,j] > 100. and lon[i,j] <= 260.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "BAB": 
        		  if lat[i,j] > -90. and lat[i,j] <= -55.: 
        		    if lon[i,j] > 285. and lon[i,j] <= 360.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "DAS": 
        		  if lat[i,j] > 60. and lat[i,j] <= 82.: 
        		    if lon[i,j] > 285. and lon[i,j] <= 315.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "LAS": 
        		  if lat[i,j] > 47. and lat[i,j] <= 60.: 
        		    if lon[i,j] > 300. and lon[i,j] <= 315.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "GLS": 
        		  if lat[i,j] > 60. and lat[i,j] <= 65.: 
        		    if lon[i,j] > 315. and lon[i,j] <= 340.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "ILS": 
        		  if lat[i,j] > 65. and lat[i,j] <= 72.: 
        		    if lon[i,j] > 320. and lon[i,j] <= 353.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "NGS": 
        		  if lat[i,j] > -90. and lat[i,j] <= -55.: 
        		    if lon[i,j] > 120. and lon[i,j] <= 210.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "NAO": 
        		  if lat[i,j] > 30. and lat[i,j] < 65.: 
        		    if lon[i,j] > 250. and lon[i,j] <= 360. or lon[i,j] <= 35.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "SNP": 
        		  if lat[i,j] > 10. and lat[i,j] <= 30.: 
        		    if lon[i,j] > 105. and lon[i,j] <= 260.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "GOM": 
        		  if lat[i,j] > 10. and lat[i,j] <= 30.: 
        		    if lon[i,j] > 260. and lon[i,j] <= 300.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "SNA": 
        		  if lat[i,j] > 10. and lat[i,j] <= 30.: 
        		    if lon[i,j] > 300. and lon[i,j] <= 360.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "NIO": 
        		  if lat[i,j] > 10. and lat[i,j] <= 30.: 
        		    if lon[i,j] > 35. and lon[i,j] <= 105.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "PWP": 
        		  if lat[i,j] > -10. and lat[i,j] <= 10.: 
        		    if lon[i,j] > 25. and lon[i,j] <= 190.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "EPO": 
        		  if lat[i,j] > -10. and lat[i,j] <= 10.: 
        		    if lon[i,j] > 190. and lon[i,j] <= 285.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "EAO": 
        		  if lat[i,j] > -10. and lat[i,j] <= 10.: 
        		    if lon[i,j] > 290. and lon[i,j] <= 360. or lon[i,j] <= 25.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "SIO": 
        		  if lat[i,j] > -50. and lat[i,j] <= -10.: 
        		    if lon[i,j] > 25. and lon[i,j] <= 130.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "SPO": 
        		  if lat[i,j] > -50. and lat[i,j] <= -10.: 
        		    if lon[i,j] > 130. and lon[i,j] <= 290.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "SAO": 
        		  if lat[i,j] > -50. and lat[i,j] <= -10.: 
        		    if lon[i,j] > 290. and lon[i,j] <= 360. or lon[i,j] <= 25.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "SOCN": 
        		  if lat[i,j] <= -50.: RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "AMR": 
        		  if lat[i,j] > -80. and lat[i,j] <= -60.: 
        		    if lon[i,j] > 210. and lon[i,j] <= 285.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "DML": 
        		  if lat[i,j] > -70. and lat[i,j] <= -53.: 
        		    if lon[i,j] > 30. and lon[i,j] <= 60.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]

        		if RgName == "WLR": 
        		  if lat[i,j] > -68. and lat[i,j] <= -55.: 
        		    if lon[i,j] > 90. and lon[i,j] <= 120.:
        		      RegionMask[i,j,k] = 1.-LandFrac[i,j]
			      
    return lat,lon,latv,lonv,LandFrac,RegionNames,RegionMask
