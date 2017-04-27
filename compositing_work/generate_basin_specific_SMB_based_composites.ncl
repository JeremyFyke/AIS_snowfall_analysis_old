load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
load "/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/local_utilies.ncl"
load "/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/compositing_work/composite_utilities.ncl"

begin

    year2sec=31557600.
    actually_generate_composites=True
    generate_basin_and_longitude_monthly_timeseries=False
    
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
    LandFrac = get_LandFrac()
    AISMask = get_AIS_mask()
    nAISMask = doubletoint(max(AISMask))
    LandFrac = where(AISMask.gt.0, LandFrac, 0.)
    LandArea = get_LandArea()   
    in = addfile("CAM_lat_lon.nc","r")
    lat_fv09x125 = in->lat
    lon_fv09x125 = in->lon
    lev_fv09x125 = in->lev
    ;load monthly SLP data
    delete(in)
    BasinMask=new(dimsizes(AISMask),"double")
    
    if (actually_generate_composites) then

       ;load monthly PRECIP data
       print("Loading monthly PRECIP")
       print("   ...RAIN...")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/RAIN/*.B1850C5CN.*nc")
       in = addfiles(all_files,"r")
       PRECIPm = in[:]->RAIN
       delete(in)  
       print("   ...SNOW...")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/SNOW/*.B1850C5CN.*nc")
       in = addfiles(all_files,"r")
       PRECIPm = PRECIPm + in[:]->SNOW  
       PRECIPm=replace_mask_values(PRECIPm,0.)
       delete(in)
       ds=dimsizes(PRECIPm)
       nt=ds(0)
       
       if (generate_basin_and_longitude_monthly_timeseries) then
	  ;Print monthly basin-integrated precip to file
	  Basin_integrated_TS_mon = new((/nt,nAISMask/),"double")
	  do b = 1,nAISMask
             print(b)
             bm1=b-1
             BasinMask=where(AISMask.eq.b,1,0)
             Basin_integrated_TS_mon(:,bm1)=calculate_integrated_timeseries(PRECIPm,LandFrac,LandArea,BasinMask)
             Basin_integrated_TS_mon(:,bm1)=Basin_integrated_TS_mon(:,bm1) / 1.e9 ; integrated AIS PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr  
	  end do
	  save_timeseries(Basin_integrated_TS_mon,"BasinAccumulationTimeSeries_monthly.nc")
	  ;lonstride=2
	  ;nl=toint(ds(2)/lonstride)
	  ;Longitude_integrated_TS_mon = new((/nt,nl/),"double") 
	  ;LongMask=new(dimsizes(AISMask),"integer")
	  ;print(dimsizes(LongMask))
	  ;do l = 0,nl-1
          ;   ls=l*2
          ;   print(l)
	  ;   LongMask(:,:)=0
	  ;   LongMask(:,ls:ls+1)=1
	  ;   BasinMask=where((AISMask.gt.0).and.(LongMask.eq.1),1,0)
          ;   wks = gsn_open_wks ("png", "BasinMask")
          ;   res=True
          ;   plotv  = gsn_csm_contour(wks, BasinMask, res)	  
          ;   Longitude_integrated_TS_mon(:,l)=calculate_integrated_timeseries(PRECIPm,LandFrac,LandArea,BasinMask)
          ;   Longitude_integrated_TS_mon(:,l)=Longitude_integrated_TS_mon(:,l) / 1.e9 ; integrated AIS PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr
	  ;end do
	  ;save_timeseries(Longitude_integrated_TS_mon,"LongitudeAccumulationTimeSeries_monthly.nc")       
       end if

       print("Loading monthly SLP")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PSL/*.B1850C5CN.*.nc")
       in = addfiles(all_files,"r")
       ListSetType (in, "cat")
       SLPm=in[:]->PSL
       delete(in)

       print("Loading monthly TREFHT")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/TREFHT/*.B1850C5CN.*.nc")
       in = addfiles(all_files,"r")
       ListSetType (in, "cat")
       TREFHTm=in[:]->TREFHT
       delete(in)

       print("Loading monthly SST")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/SST/*.B1850C5CN.*nc")
       in = addfiles(all_files,"r")
       ListSetType (in, "cat")       
       SSTm = in[:]->SST(:,0,:,:) ; highest ocean level     
       SSTm = replace_mask_values(SSTm,-1.8)
       delete(in)
       
       ;print("Loading monthly EVAP")
       ;all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/EVAP_F/*.B1850C5CN.*nc")
       ;in = addfiles(all_files,"r")
       ;ListSetType (in, "cat")        
       ;EVAPm = in[:]->EVAP_F ; highest ocean level     
       ;EVAPm = replace_mask_values(EVAPm,0.0)
       ;delete(in)
       
       print("Loading monthly EVAP")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/QFLX/*.B1850C5CN.*nc")
       in = addfiles(all_files,"r")
       ListSetType (in, "cat")        
       EVAPm = in[:]->QFLX    
       delete(in)       

       ;load monthly CICE data: is split up hemispheres - so need to glue them together.
       print("loading monthly SH sea ice")
       print("   ...SH...")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/aice/*.B1850C5CN.*sh*nc")
       in = addfiles(all_files,"r")
       ListSetType (in, "cat")        
       tmp = in[:]->aice
       delete(in)
       nt=dimsizes(tmp(:,0,0))
       ICEFRACm = new((/nt,nlat,nlon/),"double")
       ICEFRACm(:,:,:)=0.0
       nlats = dimsizes(tmp(0,:,0))
       ICEFRACm(:,0:nlats-1,:) = tmp
       delete(tmp)
       print("   ...NH...")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/aice/*.B1850C5CN.*nh*nc")	 
       in = addfiles(all_files,"r")
       tmp = in[:]->aice
       nlatn = dimsizes(tmp(0,:,0))
       ICEFRACm(:,nlat-nlatn:nlat-1,:) = tmp
       delete(tmp)
       delete(in)

       ICEFRACm=replace_mask_values(ICEFRACm,100.)
       ; change SIC to fractional representation
       ICEFRACm=ICEFRACm/100.0     

       print("Loading VQ")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/VQ/*.B1850C5CN.*nc")
       in = addfiles(all_files,"r")
       ListSetType (in, "cat")
       VQm=in[:]->VQ
       delete(in)
       
       print("Loading UQ")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/UQ/*.B1850C5CN.*.nc")
       in = addfiles(all_files,"r")
       ListSetType (in, "cat")
       UQm=in[:]->UQ
       delete(in)       

       ;load monthly UQ_vert_integral
       print("Loading monthly UQ_vert_integral")
       in = addfile("/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/jer_data_dir/Vertical_Int_Moisture_Transport_UQ.nc","r")
       UQ_VERT_INTEGRALm = in->UQ_vert_integral
       delete(in)  

       ;load monthly VQ_vert_integral
       print("Loading monthly VQ_vert_integral")
       in= addfile("/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/jer_data_dir/Vertical_Int_Moisture_Transport_VQ.nc","r")
       VQ_VERT_INTEGRALm = in->VQ_vert_integral
       delete(in)
       
       ;load monthly Z500
       print("Loading monthly z500")
       in= addfile("/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/jer_data_dir/z500.nc","r")
       Z500m = in->z500
       delete(in)
       
       ;load monthly Z850
       print("Loading monthly z850")
       in= addfile("/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/jer_data_dir/z850.nc","r")
       Z850m = in->z850
       delete(in)       
    
    end if
    
    print("Calculating integrated annual PRECIP for compositing.")
    in = addfile("/glade/u/home/lenaerts/work/CESM_LE/PRECIP_B1850C5CN_400-2200.nc","r")
    time = in->time
    PRECIP_annual = (in->SNOW + in->RAIN) * year2sec / 1000.    ;mm/s -> m/yr
    PRECIP_annual=replace_mask_values(PRECIP_annual,0.)
    nt=dimsizes(time)
    BasinMask=new(dimsizes(AISMask),"double")
    Basin_integrated_TS_avg = new((/nAISMask/),"double")
    Basin_integrated_TS_std = new((/nAISMask/),"double")
    Basin_integrated_TS_5th_percentile = new((/nAISMask/),"double")
    Basin_integrated_TS_95th_percentile = new((/nAISMask/),"double")    
    Basin_integrated_TS = new((/nt,nAISMask/),"double")
    AIS_integrated_TS = new((/nt/),"double")
    print("Compositing on annual integrated Precip for individual basins,")
    print("Generating CAM-grid composite output files...")

    AIS_integrated_TS(:)=0.
    ; calculate basin-specific annual mass flux time series.
    do b = 1,nAISMask
       print("  ")
       print(str_concat((/"Building annual integrated precip timeseries for basin ",tostring(b),":"/)))
       bm1=b-1
       BasinMask=where(AISMask.eq.b,1,0)
       ; calculate yearly totals
       Basin_integrated_TS(:,bm1)=calculate_integrated_timeseries(PRECIP_annual,LandFrac,LandArea,BasinMask)
       Basin_integrated_TS(:,bm1)=Basin_integrated_TS(:,bm1) / 1.e9 ; integrated AIS PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr
       statx=stat_dispersion(Basin_integrated_TS(:,bm1),False)
       Basin_integrated_TS_avg(bm1)=statx(0)
       Basin_integrated_TS_std(bm1)=statx(1)
       Basin_integrated_TS_5th_percentile(bm1)=statx(24)
       Basin_integrated_TS_95th_percentile(bm1)=statx(25)
       AIS_integrated_TS(:)=AIS_integrated_TS(:)+dble2flt(Basin_integrated_TS(:,bm1))
       print("Average basin PRECIP="+Basin_integrated_TS_avg(bm1))
       print("Standard Deviation basin PRECIP="+Basin_integrated_TS_std(bm1))
       print("5/95 percentile basin PRECIP="+Basin_integrated_TS_5th_percentile(bm1)+"/"+Basin_integrated_TS_95th_percentile(bm1))
       print("Min basin PRECIP="+statx(2))    
       print("Max basin PRECIP="+statx(14))             
     end do
     print("Min/Max basin PRECIP="+min(Basin_integrated_TS_avg)+"/"+max(Basin_integrated_TS_avg))
     AIS_integrated_TS_avg=avg(AIS_integrated_TS)
     AIS_integrated_TS_std=stddev(AIS_integrated_TS)
     
     print("Average AIS PRECIP="+AIS_integrated_TS_avg)
     print("Standard Deviation AIS PRECIP="+AIS_integrated_TS_std)         

     AIS_integrated_CV=AIS_integrated_TS_std/AIS_integrated_TS_avg
     Basin_integrated_CV=avg(Basin_integrated_TS_std/Basin_integrated_TS_avg)
     print("Average AIS CV="+AIS_integrated_CV)
     print("Basin-average CV="+Basin_integrated_CV)
     
     ;plot_AIS_and_basin_accum_PDFs(AIS_integrated_TS,AIS_integrated_TS_avg,\
     ;                              Basin_integrated_TS,Basin_integrated_TS_avg,\
     ;				   nAISMask)			   
     save_timeseries(Basin_integrated_TS,"BasinAccumulationTimeSeries.nc")
     save_timeseries(AIS_integrated_TS,"AISAccumulationTimeSeries.nc")
     save_timeseries(time,"TimeSeriestime.nc")
     

     if (actually_generate_composites) then	
	do b = 1,nAISMask
	   bm1=b-1
	   TS=Basin_integrated_TS(:,bm1)
	   fname=str_concat((/"Basin",tostring(b),"_ssticetemp_fv0.9x1.25.nc"/))
	   print(str_concat((/"Compositing for basin file ",fname,":"/)))

	   print("Generating high composite...")
	   threshold = Basin_integrated_TS_95th_percentile(bm1)
	   ;threshold = Basin_integrated_TS_avg(bm1) + 2*Basin_integrated_TS_std(bm1)
	   indt = ind(TS.gt.threshold)
	   generate_composites("high",indt,\
	                       SSTm,ICEFRACm,PRECIPm,SLPm,\
			       UQm,VQm,UQ_VERT_INTEGRALm,VQ_VERT_INTEGRALm,EVAPm,\
			       TREFHTm,Z850m,Z500m,\
			       lev_fv09x125,lat_fv09x125,lon_fv09x125,fname)     
           delete(indt)

	   print("Generating low composite...")
	   threshold = Basin_integrated_TS_5th_percentile(bm1)   
	   ;threshold = Basin_integrated_TS_avg(bm1) - 2*Basin_integrated_TS_std(bm1)
	   indt = ind(TS.lt.threshold)
	   generate_composites("low",indt,\
			       SSTm,ICEFRACm,PRECIPm,SLPm,\
			       UQm,VQm,UQ_VERT_INTEGRALm,VQ_VERT_INTEGRALm,EVAPm,\
			       TREFHTm,Z850m,Z500m,\
			       lev_fv09x125,lat_fv09x125,lon_fv09x125,fname)
           delete(indt)

	   ;print("Generating mean composite...")
	   ;threshold = Basin_integrated_TS_avg(bm1) - 0.5*Basin_integrated_TS_std(bm1)
	   ;threshold2 = Basin_integrated_TS_avg(bm1) + 0.5*Basin_integrated_TS_std(bm1)
	   ;indt = ind(TS.gt.threshold.and.TS.lt.threshold2) 
	   ;generate_composites("mean",indt,\
	   ;                    SSTm,ICEFRACm,PRECIPm,SLPm,\
	;		       UQm,VQm,UQ_VERT_INTEGRALm,VQ_VERT_INTEGRALm,EVAPm,\
	;		       TREFHTm,Z850m,Z500m,\
	;		       lev_fv09x125,lat_fv09x125,lon_fv09x125,fname)
         ;  delete(indt)	   
	end do
     end if
        
end
