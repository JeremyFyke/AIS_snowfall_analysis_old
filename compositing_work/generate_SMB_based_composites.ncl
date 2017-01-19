load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 

begin
    print("Loading masks, lons, lats, areas, etc.")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/masks_cice_sh.nc","r")
    cicemask_sh = in->tmask
    area_sh = in->uarea
    lat_sh = in->TLAT
    lon_sh = in->TLON
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/masks_cice_nh.nc","r")
    cicemask_nh = in->tmask
    area_nh = in->uarea
    lat_nh_gx1v6 = in->TLAT
    lon_nh_gx1v6 = in->TLON
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/masks_pop.nc","r")
    poparea = in->TAREA
    lat_gx1v6 = in->TLAT
    lon_gx1v6 = in->TLONG
    delete_VarAtts(lat_gx1v6,"_CoordinateAxisType")
    delete_VarAtts(lon_gx1v6,"_CoordinateAxisType")
    nlat=dimsizes(poparea(:,0))
    nlon=dimsizes(poparea(0,:))
    in = addfile("/glade/scratch/hailong/archive/composite_ICE_wtag_mean/atm/hist/composite_ICE_wtag_mean.cam.h0.0001-01.nc","r")
    LandFrac = rm_single_dims(in->LANDFRAC)
    in = addfile("/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/data/AIS_Full_basins_Zwally_CESMgrid.nc","r")
    AISMask = in->Zwallybasins 
    LandFrac = where(AISMask.gt.0, LandFrac, 0.)
    in = addfile("/glade/p/cesm/cseg/inputdata/lnd/clm2/griddata/griddata_0.9x1.25_070212.nc","r")
    LandArea = in->AREA    

    function replace_mask_values(field,ReplaceValue)
    local field1d
    begin
	field1d=ndtooned(field)
	field1d (ind(ismissing( field1d ))) = ReplaceValue;
	field=onedtond(field1d,dimsizes(field))
	delete_VarAtts(field,(/"missing_value","_FillValue","coordinates"/))
	return(field)
    end
    
    print("Loading, calculating integrated annual PRECIP for compositing.")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/PRECIP_B1850C5CN_400-2200.nc","r") 
    time = in->time
    PRECIP = in->SNOW + in->RAIN
    PRECIP=replace_mask_values(PRECIP,0.)
    
    ; calculate yearly totals
    Composite_year = new(dimsizes(time),"double")
    do i = 0, dimsizes(time)-1,1
	    Composite_year(i) = sum(PRECIP(i,:,:)*LandFrac*LandArea*AISMask) ; integrated AIS PRECIP
    end do
    print("Average AIS PRECIP")
    print(avg(Composite_year))
    print("Standard Deviation AIS PRECIP")    
    print(stddev(Composite_year))
    print("Min AIS PRECIP")    
    print(min(Composite_year))    
    print("Max AIS PRECIP")    
    print(max(Composite_year))        

    ;load monthly SST data
    print("Loading global SST")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/SST_monthly_B1850C5CN_400-2200.nc","r")
    SST = in->SST(:,0,:,:) ; highest ocean level     
    SST = replace_mask_values(SST,-1.8)

    ;load monthly CICE data: is split up hemispheres - so need to glue them together.
    print("loading SH sea ice")
    print("   ...SH...")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/aice_monthly_B1850C5CN_400-2200-sh.nc","r")
    tmp = in->aice
    nt=dimsizes(tmp(:,0,0))
    ICEFRACm = new((/nt,nlat,nlon/),"double")
    ICEFRACm(:,:,:)=0.0
    nlats = dimsizes(tmp(0,:,0))
    ICEFRACm(:,0:nlats-1,:) = tmp
    print("   ...NH...")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/aice_monthly_B1850C5CN_400-2200-nh.nc","r")
    delete(tmp)
    tmp = in->aice
    nlatn = dimsizes(tmp(0,:,0))
    ICEFRACm(:,nlat-nlatn:nlat-1,:) = tmp
    delete(tmp)
    ICEFRACm=replace_mask_values(ICEFRACm,100.)
    ; change SIC to fractional representation
    ICEFRACm=ICEFRACm/100.0     

    ;load monthly PRECIP data
    print("Loading global PRECIP")
    print("   ...RAIN...")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/RAIN_monthly_B1850C5CN_400-2200.nc","r")
    PRECIPm = in->RAIN     
    print("   ...SNOW...")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/SNOW_monthly_B1850C5CN_400-2200.nc","r")
    PRECIPm = PRECIPm + in->SNOW  
    PRECIPm=replace_mask_values(PRECIPm,0.)

    print("Compositing on annual AIS integrated SMB...")

    Composite_year2 = Composite_year
    qsort(Composite_year2)
    composite_avg=avg(Composite_year2)
    composite_std=stddev(Composite_year2)
    
    
    function generate_composites(indt,SST,ICEFRACm,PRECIP)
    local indt2,i,j
    begin
        indt2 = new(dimsizes(indt)*12,"integer")
        do i = 0, dimsizes(indt)-1,1	
	    do j = 0,11
	    indt2(i*12+j) = indt(i)*12+j
	    end do
        end do
	SIC_composite = ICEFRACm(indt2,:,:)
	SST_composite = SST(indt2,:,:)
	PRECIP_composite = PRECIP(indt2,:,:)
        return([/SIC_composite,SST_composite,PRECIP_composite/])
    end
    
    print("Generating high composite...")
    ;threshold = Composite_year2(dimsizes(time)-31) ; 30 highest years
    threshold = composite_avg + 2*composite_std
    indt = ind(Composite_year.gt.threshold)
    ans=generate_composites(indt,SST,ICEFRACm,PRECIPm)
    SIC_high_gx1v6=ans[0]
    SST_high_gx1v6=ans[1]
    PRECIP_high_fv09x125=ans[2]

    print("Generating low composite...")    
    ;threshold = Composite_year2(30) ; 30 lowest years
    threshold = composite_avg - 2*composite_std
    indt = ind(Composite_year.lt.threshold)
    ans=generate_composites(indt,SST,ICEFRACm,PRECIPm)
    SIC_low_gx1v6=ans[0]
    SST_low_gx1v6=ans[1]
    PRECIP_low_fv09x125=ans[2]    
    
    print("Generating mean composite...")
    ;threshold = Composite_year2(dimsizes(time)/2-15) ; 30 mean years
    ;threshold2 = Composite_year2(dimsizes(time)/2+16)
    threshold = composite_avg - composite_std
    threshold2 = composite_avg + composite_std
    indt = ind(Composite_year.gt.threshold.and.Composite_year.lt.threshold2) 
    ans=generate_composites(indt,SST,ICEFRACm,PRECIPm)
    SIC_mean_gx1v6=ans[0]
    SST_mean_gx1v6=ans[1]
    PRECIP_mean_fv09x125=ans[2]
    
    print("Done compositing.") 
    delete(ICEFRACm)
    delete(SST)
    delete(PRECIP)    
       
    procedure write_composite_data(CompositePeriod,SST,SIC,PRECIP,lat,lon,OutPutFileName)
    local fname,fout,dimNames,dimSizes,dimUnlim
    begin
    
        delete(time)
	time    = new(360,"double")
	date    = new(360,"integer")
	datesec = new(360,"integer")
	it=0
	totaltime=0.
	do yr = 0000,0029
        	year=yr
        	year@calendar="365_day"
        	do month = 01,12
                	totaltime=totaltime+int2flt(days_in_month(year,month))
			t=totaltime-int2flt(days_in_month(year,month))/2.
			time(it)=t
			t@calendar="365_day"
			t@units="days since 0000-1-1"
			date(it)=cd_calendar(t,-2)-1
			datesec(it)=floattointeger((totaltime-t)*86400)
			it=it+1
		end do
	end do       

	fname=str_concat((/"output/",CompositePeriod,"/",OutPutFileName/))
	print("Writing "+fname)
	system("rm -f "+fname)
	fout = addfile(fname,"c")
        
	dimNames = (/"time", "lat", "lon"/)
	dsizes_lat=dimsizes(lat)
	dsizes_lon=dimsizes(lon)
	dimSizes = (/dimsizes(year), dsizes_lat, dsizes_lon/)
	dimUnlim = (/ True, False, False/)   
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

	filevardef(fout, "time"     ,typeof(time)     , (/"time"/))
	filevardef(fout, "date"     ,typeof(date)     , (/"time"/))
	filevardef(fout, "datesec"  ,typeof(datesec)  , (/"time"/))	
	filevardef(fout, "lat"      ,typeof(lat)      , (/"lat"/) )                          
	filevardef(fout, "lon"      ,typeof(lon)      , (/"lon"/) )	                
	filevardef(fout, "ICEFRAC"  ,typeof(SIC)      , (/"time", "lat", "lon" /) )
	filevardef(fout, "SST"      ,typeof(SST)      , (/"time", "lat", "lon" /) )
	filevardef(fout, "PRECIP"   ,typeof(PRECIP)   , (/"time", "lat", "lon" /) )	
              
	filevarattdef(fout,"lat"    ,lat)                  
	filevarattdef(fout,"lon"    ,lon)      
	filevarattdef(fout,"ICEFRAC",SIC)    
	filevarattdef(fout,"SST",SST)
	filevarattdef(fout,"PRECIP",PRECIP)	    

	setfileoption(fout,"DefineMode",False)

	fout->time    = (/time/)
	fout->date    = (/date/)
	fout->datesec = (/datesec/)
	fout->lat     = (/lat/)
	fout->lon     = (/lon/) 
	fout->ICEFRAC = (/SIC/) 
	fout->SST     = (/SST/)
	fout->PRECIP  = (/PRECIP/)	
    end
    
    print("Remapping SST/SIC to CAM grid, writing CAM-gridded BC files.")

    in = addfile("CAM_lat_lon.nc","r")
    lat_fv09x125 = in->lat
    lon_fv09x125 = in->lon
      
    print("Processing low composite...")
    SST = PopLatLon(SST_low_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")    
    Arr1D=ndtooned(SST)
    Arr1D (ind(ismissing( Arr1D ))) = -1.8
    SST=onedtond(Arr1D,dimsizes(SST))
    delete_VarAtts(SST,(/"missing_value","_FillValue"/))
    
    SIC = PopLatLon(SIC_low_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    delete(Arr1D)
    Arr1D=ndtooned(SIC)
    Arr1D (ind(ismissing( Arr1D ))) = 1.
    SIC=onedtond(Arr1D,dimsizes(SIC))
    delete_VarAtts(SIC,(/"missing_value","_FillValue"/))       
    write_composite_data("low",SST,SIC,PRECIP_low_fv09x125,lat_fv09x125,lon_fv09x125,"ssticetemp_fv0.9x1.25.nc")

    print("Processing mean composite...")
    delete(SST)
    delete(SIC)
    delete(Arr1D)
    SST = PopLatLon(SST_mean_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")    
    Arr1D=ndtooned(SST)
    Arr1D (ind(ismissing( Arr1D ))) = -1.8
    SST=onedtond(Arr1D,dimsizes(SST))
    delete_VarAtts(SST,(/"missing_value","_FillValue"/))
    
    SIC = PopLatLon(SIC_mean_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    delete(Arr1D)
    Arr1D=ndtooned(SIC)
    Arr1D (ind(ismissing( Arr1D ))) = 1.
    SIC=onedtond(Arr1D,dimsizes(SIC))
    delete_VarAtts(SIC,(/"missing_value","_FillValue"/))       
    write_composite_data("mean",SST,SIC,PRECIP_mean_fv09x125,lat_fv09x125,lon_fv09x125,"ssticetemp_fv0.9x1.25.nc")    

    print("Processing high composite...") 
    delete(SST)
    delete(SIC)
    delete(Arr1D)
    SST = PopLatLon(SST_high_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")    
    Arr1D=ndtooned(SST)
    Arr1D (ind(ismissing( Arr1D ))) = -1.8
    SST=onedtond(Arr1D,dimsizes(SST))
    delete_VarAtts(SST,(/"missing_value","_FillValue"/))
    
    SIC = PopLatLon(SIC_high_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    delete(Arr1D)
    Arr1D=ndtooned(SIC)
    Arr1D (ind(ismissing( Arr1D ))) = 1.
    SIC=onedtond(Arr1D,dimsizes(SIC))
    delete_VarAtts(SIC,(/"missing_value","_FillValue"/))       
    write_composite_data("high",SST,SIC,PRECIP_high_fv09x125,lat_fv09x125,lon_fv09x125,"ssticetemp_fv0.9x1.25.nc")
    
end
