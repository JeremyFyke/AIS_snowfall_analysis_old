load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
load "/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/local_utilies.ncl"
load "/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/compositing_work/composite_utilities.ncl"

begin

    year2sec=31557600.
    actually_generate_composites=False
    
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
    ;load monthly SLP data
    delete(in)

    if (actually_generate_composites) then

       print("Loading monthly SLP")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PSL/*.B1850C5CN.*.nc")
       in = addfiles(all_files,"r")
       ListSetType (in, "cat")
       SLPm=in[:]->PSL
       delete(in)

       print("Loading monthly SST")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/SST/*.B1850C5CN.*nc")
       in = addfiles(all_files,"r")
       SSTm = in[:]->SST(:,0,:,:) ; highest ocean level     
       SSTm = replace_mask_values(SSTm,-1.8)
       delete(in)
       
       print("Loading monthly EVAP")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/ocn/proc/tseries/monthly/EVAP_F/*.B1850C5CN.*nc")
       in = addfiles(all_files,"r")
       EVAPm = in[:]->EVAP_F ; highest ocean level     
       EVAPm = replace_mask_values(EVAPm,0.0)
       delete(in)       

       ;load monthly CICE data: is split up hemispheres - so need to glue them together.
       print("loading monthly SH sea ice")
       print("   ...SH...")
       all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/ice/proc/tseries/monthly/aice/*.B1850C5CN.*sh*nc")
       in = addfiles(all_files,"r")
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
    Basin_integrated_TS = new((/nt,nAISMask/),"double")
    AIS_integrated_TS = new((/nt/),"double")
    print("Compositing on annual integrated Precip for individual basins,")
    print("Generating CAM-grid composite output files...")

    AIS_integrated_TS(:)=0.
    ; calculate basin-specific mass flux time series.
    do b = 1,nAISMask
       print(str_concat((/"Building annual integrated precip timeseries for basin ",tostring(b),":"/)))
       bm1=b-1
       BasinMask=where(AISMask.eq.b,1,0)
       ; calculate yearly totals
       Basin_integrated_TS(:,bm1)=calculate_integrated_timeseries(PRECIP_annual,LandFrac,LandArea,BasinMask)
       Basin_integrated_TS(:,bm1)=Basin_integrated_TS(:,bm1) / 1.e9 ; integrated AIS PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr
       Basin_integrated_TS_avg(bm1)=avg(Basin_integrated_TS(:,bm1))
       Basin_integrated_TS_std(bm1)=stddev(Basin_integrated_TS(:,bm1))
       print("Average basin PRECIP")
       print(Basin_integrated_TS_avg(bm1))
       print("Standard Deviation basin PRECIP")
       print(Basin_integrated_TS_std(bm1))
       print("Min basin PRECIP")    
       print(min(Basin_integrated_TS(:,bm1)))    
       print("Max basin PRECIP")    
       print(max(Basin_integrated_TS(:,bm1)))             
       AIS_integrated_TS(:)=AIS_integrated_TS(:)+dble2flt(Basin_integrated_TS(:,bm1))
       
     end do
     AIS_integrated_TS_avg=avg(AIS_integrated_TS)
     AIS_integrated_TS_std=stddev(AIS_integrated_TS)
     
     print("Average AIS PRECIP")
     print(AIS_integrated_TS_avg)
     print("Standard Deviation AIS PRECIP")
     print(AIS_integrated_TS_std)         

     AIS_integrated_CV=AIS_integrated_TS_std/AIS_integrated_TS_avg
     Basin_integrated_CV=avg(Basin_integrated_TS_std/Basin_integrated_TS_avg)
     print("Average AIS CV")
     print(AIS_integrated_CV)
     print("Basin-average CV")
     print(Basin_integrated_CV)

     wks = gsn_open_wks ("png","../figs/Precip_normalized_distribution")
     Vector=AIS_integrated_TS(:)/AIS_integrated_TS_avg
     bin_min=0.2
     bin_max=2.0
     opt=True
     opt@bin_min=bin_min
     opt@bin_max=bin_max
     nbins=75
     probability=pdfx( Vector , nbins , opt)
     res=True
     res@xyLineThicknessF       = 4
     res@tiYAxisString          = "PDF (%)"
     res@tiXAxisString          = "Normalized integrated accumulation" 
     res@trXMinF = bin_min
     res@trXMaxF = bin_max
     res@gsnDraw          = False
     res@gsnFrame         = False
     
     plot0=gsn_csm_xy (wks,probability@bin_center,probability,res) 
     
     res@xyLineThicknessF       = 2
     res@xyLineColor            = "grey"         
     do b = 1,nAISMask
        bm1=b-1
        Vector=Basin_integrated_TS(:,bm1)/Basin_integrated_TS_avg(bm1)
        probability=pdfx( Vector , nbins , opt)     
	plot1=gsn_csm_xy (wks,probability@bin_center,probability,res)
	overlay(plot0,plot1)
     end do
     draw(plot0)
     frame(wks)
     save_timeseries(Basin_integrated_TS,"BasinAccumulationTimeSeries.nc")
     save_timeseries(AIS_integrated_TS,"AISAccumulationTimeSeries.nc")
     save_timeseries(time,"TimeSeriestime.nc")         

     if (actually_generate_composites) then
	do b = 1,nAISMask
	   bm1=b-1
	   fname=str_concat((/"Basin",tostring(b),"_ssticetemp_fv0.9x1.25.nc"/))
	   print(str_concat((/"Compositing for basin file ",fname,":"/)))
	   print("Generating high composite...")
	   threshold = Basin_integrated_TS_avg(bm1) + 2*Basin_integrated_TS_std(bm1)
	   TS=Basin_integrated_TS(:,bm1)
	   indt = ind(TS.gt.threshold)
	   ans=generate_composites(indt,SSTm,ICEFRACm,PRECIPm,SLPm,UQ_VERT_INTEGRALm,VQ_VERT_INTEGRALm,EVAPm)
	   SIC_gx1v6=ans[0]
	   SST_gx1v6=ans[1]
	   PRECIP=ans[2]
	   SLP=ans[3]
	   UQ_VERT_INTEGRAL=ans[4]
	   VQ_VERT_INTEGRAL=ans[5]
	   EVAP_gx1v6=ans[6]
	   SST = PopLatLon(SST_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   SST = replace_mask_values(SST,-1.8)
	   SIC = PopLatLon(SIC_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   SIC = replace_mask_values(SIC,1.)
	   EVAP = PopLatLon(EVAP_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   EVAP = replace_mask_values(EVAP,0.)
	   write_composite_data("high",SST,SIC,PRECIP,SLP,UQ_VERT_INTEGRAL,VQ_VERT_INTEGRAL,EVAP,lat_fv09x125,lon_fv09x125,fname)
	   delete([/TS,indt,SST,SIC,SIC_gx1v6,SST_gx1v6,EVAP_gx1v6,PRECIP,SLP,UQ_VERT_INTEGRAL,VQ_VERT_INTEGRAL,EVAP/])       

	   print("Generating low composite...")    
	   threshold = Basin_integrated_TS_avg(bm1) - 2*Basin_integrated_TS_std(bm1)
	   TS=Basin_integrated_TS(:,bm1)
	   indt = ind(TS.lt.threshold)
	   ans=generate_composites(indt,SSTm,ICEFRACm,PRECIPm,SLPm,UQ_VERT_INTEGRALm,VQ_VERT_INTEGRALm,EVAPm)
	   SIC_gx1v6=ans[0]
	   SST_gx1v6=ans[1]
	   PRECIP=ans[2]
	   SLP=ans[3]
	   UQ_VERT_INTEGRAL=ans[4]
	   VQ_VERT_INTEGRAL=ans[5]
	   EVAP_gx1v6=ans[6]
	   SST = PopLatLon(SST_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   SST = replace_mask_values(SST,-1.8)
	   SIC = PopLatLon(SIC_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   SIC = replace_mask_values(SIC,1.)
	   EVAP = PopLatLon(EVAP_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   EVAP = replace_mask_values(EVAP,0.)
	   write_composite_data("low",SST,SIC,PRECIP,SLP,UQ_VERT_INTEGRAL,VQ_VERT_INTEGRAL,EVAP,lat_fv09x125,lon_fv09x125,fname)
	   delete([/TS,indt,SST,SIC,SIC_gx1v6,SST_gx1v6,EVAP_gx1v6,PRECIP,SLP,UQ_VERT_INTEGRAL,VQ_VERT_INTEGRAL,EVAP/])       

	   print("Generating mean composite...")
	   threshold = Basin_integrated_TS_avg(bm1) - Basin_integrated_TS_std(bm1)
	   threshold2 = Basin_integrated_TS_avg(bm1) + Basin_integrated_TS_std(bm1) 
	   TS=Basin_integrated_TS(:,bm1)     
	   indt = ind(TS.gt.threshold.and.TS.lt.threshold2) 
	   ans=generate_composites(indt,SSTm,ICEFRACm,PRECIPm,SLPm,UQ_VERT_INTEGRALm,VQ_VERT_INTEGRALm,EVAPm)
	   SIC_gx1v6=ans[0]
	   SST_gx1v6=ans[1]
	   PRECIP=ans[2]
	   SLP=ans[3]
	   UQ_VERT_INTEGRAL=ans[4]
	   VQ_VERT_INTEGRAL=ans[5]
	   EVAP_gx1v6=ans[6]
	   SST = PopLatLon(SST_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   SST = replace_mask_values(SST,-1.8)
	   SIC = PopLatLon(SIC_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   SIC = replace_mask_values(SIC,1.)
	   EVAP = PopLatLon(EVAP_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
	   EVAP = replace_mask_values(EVAP,0.)
	   write_composite_data("mean",SST,SIC,PRECIP,SLP,UQ_VERT_INTEGRAL,VQ_VERT_INTEGRAL,EVAP,lat_fv09x125,lon_fv09x125,fname)
	   delete([/TS,indt,SST,SIC,SIC_gx1v6,SST_gx1v6,EVAP_gx1v6,PRECIP,SLP,UQ_VERT_INTEGRAL,VQ_VERT_INTEGRAL,EVAP/])
	end do
     end if
        
end
