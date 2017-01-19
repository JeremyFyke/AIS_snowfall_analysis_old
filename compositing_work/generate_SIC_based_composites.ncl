load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 

begin
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

    ; annual data used for composites
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/aice_yearly_B1850C5CN_400-2200-sh.nc","r") 
    time = in->time
    ICEFRAC = in->aice

    ; calculate yearly totals
    SIC_year = new(dimsizes(time),"double")
    print("Calculating ice concentration")
    do i = 0, dimsizes(time)-1,1
	    SIC_year(i) = sum(ICEFRAC(i,:,:)*area_sh*cicemask_sh) ; southern hemisphere totals.
    end do

    ;load monthly SST data
    print("Loading global SST")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/SST_monthly_B1850C5CN_400-2200.nc","r")
    SST = in->SST(:,0,:,:) ; highest ocean level     
    
    ;load monthly CICE data: is split up hemispheres - so need to glue them together.
    print("loading SH sea ice")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/aice_monthly_B1850C5CN_400-2200-sh.nc","r")
    tmp = in->aice
    nt=dimsizes(tmp(:,0,0))
    ICEFRACm = new((/nt,nlat,nlon/),"double")
    ICEFRACm(:,:,:)=0.0
    nlats = dimsizes(tmp(0,:,0))
    ICEFRACm(:,0:nlats-1,:) = tmp
    print("loading NH sea ice")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/aice_monthly_B1850C5CN_400-2200-nh.nc","r")
    delete(tmp)
    tmp = in->aice
    nlatn = dimsizes(tmp(0,:,0))
    ICEFRACm(:,nlat-nlatn:nlat-1,:) = tmp
    delete(tmp)
    
    ;fill in missing values
    SST1D=ndtooned(SST)
    SST1D (ind(ismissing( SST1D ))) = -1.8;
    SST=onedtond(SST1D,dimsizes(SST))
    delete_VarAtts(SST,(/"missing_value","_FillValue","coordinates"/))
    SIC1D=ndtooned(ICEFRACm)
    SIC1D (ind(ismissing( SIC1D ))) = 100.
    ICEFRACm=onedtond(SIC1D,dimsizes(ICEFRACm))
    delete_VarAtts(ICEFRACm,(/"missing_value","_FillValue","coordinates"/))
    
    ; change SIC to fractional representation
    ICEFRACm=ICEFRACm/100.0

    delete(SST1D)
    delete(SIC1D)

    print("Compositing monthly data ice years in sea ice extent")
    ; find the 30 most extreme years in sea ice extent
    SIC_year2 = SIC_year
    qsort(SIC_year2)
    threshold = SIC_year2(dimsizes(time)-31) ; 30 highest years
    indt = ind(SIC_year.gt.threshold)

    indt2 = new(dimsizes(indt)*12,"integer")

    do i = 0, dimsizes(indt)-1,1	
	    do j = 0,11
	    indt2(i*12+j) = indt(i)*12+j
	    end do
    end do

    SIC_year_high_gx1v6 = ICEFRACm(indt2,:,:) ; 
    SST_year_high_gx1v6 = SST(indt2,:,:) ; 

    threshold = SIC_year2(30) ; 30 highest years
    indt = ind(SIC_year.lt.threshold) 

    do i = 0, dimsizes(indt)-1,1	
	    do j = 0,11
	    indt2(i*12+j) = indt(i)*12+j
	    end do
    end do

    SIC_year_low_gx1v6 = ICEFRACm(indt2,:,:) ; 
    SST_year_low_gx1v6 = SST(indt2,:,:) ; 

    threshold = SIC_year2(dimsizes(time)/2-15) ; 30 highest years
    threshold2 = SIC_year2(dimsizes(time)/2+16) ; 30 highest years

    indt = ind(SIC_year.gt.threshold.and.SIC_year.lt.threshold2) 

    do i = 0, dimsizes(indt)-1,1	
	    do j = 0,11
	    indt2(i*12+j) = indt(i)*12+j
	    end do
    end do

    SIC_year_mean_gx1v6 = ICEFRACm(indt2,:,:)
    SST_year_mean_gx1v6 = SST(indt2,:,:)
    
    delete(ICEFRACm)
    delete(SST)
    
    procedure write_composite_data(CompositePeriod,SST,SIC,lat,lon,OutPutFileName)
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
	if (num(dsizes_lat) .eq. 1) then ; if on a regular CAM grid
	         dimSizes = (/dimsizes(year), dsizes_lat, dsizes_lon/) 
	else                           ; if on a tripole POP grid
	         dimSizes = (/dimsizes(year), dimsizes(lat(:,0)), dimsizes(lat(0,:))/) 
	end if
	dimUnlim = (/ True, False, False/)   
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

	filevardef(fout, "time"     ,typeof(time)     , (/"time"/))
	filevardef(fout, "date"     ,typeof(date)     , (/"time"/))
	filevardef(fout, "datesec"  ,typeof(datesec)  , (/"time"/))
	if (num(dsizes_lat) .eq. 1) then ; if on a regular CAM grid	
		 filevardef(fout, "lat"      ,typeof(lat)      , (/"lat"/) )                          
		 filevardef(fout, "lon"      ,typeof(lon)      , (/"lon"/) )	
	else                           ; if on a tripole POP grid  	
		 filevardef(fout, "lat"      ,typeof(lat)      , (/"lat","lon"/) )                          
		 filevardef(fout, "lon"      ,typeof(lon)      , (/"lat","lon"/) )
	end if	                       
	filevardef(fout, "ICEFRAC"  ,typeof(SIC)      , (/"time", "lat", "lon" /) )
	filevardef(fout, "SST"      ,typeof(SST)      , (/"time", "lat", "lon" /) ) 
              
	filevarattdef(fout,"lat"    ,lat)                  
	filevarattdef(fout,"lon"    ,lon)      
	filevarattdef(fout,"ICEFRAC",SIC)    
	filevarattdef(fout,"SST",SST)    

	setfileoption(fout,"DefineMode",False)

	fout->time    = (/time/)
	fout->date    = (/date/)
	fout->datesec = (/datesec/)
	fout->lat     = (/lat/)
	fout->lon     = (/lon/) 
	fout->ICEFRAC = (/SIC/) 
	fout->SST     = (/SST/) 
    end

    write_composite_data("low",SST_year_low_gx1v6,SIC_year_low_gx1v6,lat_gx1v6,lon_gx1v6,"ssticetemp_gx1v6.nc")
    write_composite_data("mean",SST_year_mean_gx1v6,SIC_year_mean_gx1v6,lat_gx1v6,lon_gx1v6,"ssticetemp_gx1v6.nc")
    write_composite_data("high",SST_year_high_gx1v6,SIC_year_high_gx1v6,lat_gx1v6,lon_gx1v6,"ssticetemp_gx1v6.nc")
    
    print("Remapping SST/SIC to CAM grid, writing CAM-gridded BC files.")
    
    ;low-SIC composite
    in = addfile("CAM_lat_lon.nc","r")
    lat_fv09x125 = in->lat
    lon_fv09x125 = in->lon
    SST = PopLatLon(SST_year_low_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")    
    Arr1D=ndtooned(SST)
    Arr1D (ind(ismissing( Arr1D ))) = -1.8
    SST=onedtond(Arr1D,dimsizes(SST))
    delete_VarAtts(SST,(/"missing_value","_FillValue"/))
    
    SIC = PopLatLon(SIC_year_low_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    delete(Arr1D)
    Arr1D=ndtooned(SIC)
    Arr1D (ind(ismissing( Arr1D ))) = 1.
    SIC=onedtond(Arr1D,dimsizes(SIC))
    delete_VarAtts(SIC,(/"missing_value","_FillValue"/))       
    write_composite_data("low",SST,SIC,lat_fv09x125,lon_fv09x125,"ssticetemp_fv0.9x1.25.nc")
    
    ;mean-SIC composite
    delete(SST)
    delete(SIC)
    delete(Arr1D)
    SST = PopLatLon(SST_year_mean_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")    
    Arr1D=ndtooned(SST)
    Arr1D (ind(ismissing( Arr1D ))) = -1.8
    SST=onedtond(Arr1D,dimsizes(SST))
    delete_VarAtts(SST,(/"missing_value","_FillValue"/))
    
    SIC = PopLatLon(SIC_year_mean_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    delete(Arr1D)
    Arr1D=ndtooned(SIC)
    Arr1D (ind(ismissing( Arr1D ))) = 1.
    SIC=onedtond(Arr1D,dimsizes(SIC))
    delete_VarAtts(SIC,(/"missing_value","_FillValue"/))       
    write_composite_data("mean",SST,SIC,lat_fv09x125,lon_fv09x125,"ssticetemp_fv0.9x1.25.nc")    
   
    ;high-SIC composite
    delete(SST)
    delete(SIC)
    delete(Arr1D)
    SST = PopLatLon(SST_year_high_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")    
    Arr1D=ndtooned(SST)
    Arr1D (ind(ismissing( Arr1D ))) = -1.8
    SST=onedtond(Arr1D,dimsizes(SST))
    delete_VarAtts(SST,(/"missing_value","_FillValue"/))
    
    SIC = PopLatLon(SIC_year_high_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    delete(Arr1D)
    Arr1D=ndtooned(SIC)
    Arr1D (ind(ismissing( Arr1D ))) = 1.
    SIC=onedtond(Arr1D,dimsizes(SIC))
    delete_VarAtts(SIC,(/"missing_value","_FillValue"/))       
    write_composite_data("high",SST,SIC,lat_fv09x125,lon_fv09x125,"ssticetemp_fv0.9x1.25.nc")
    
end
