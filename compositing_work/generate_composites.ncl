load "/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/local_utilies.ncl"
load "/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/compositing_work/composite_utilities.ncl"

year2sec=31557600.

CompositePrecipCLM=True
CompositeSLP      =True

print("Loading masks, lons, lats, areas, etc.")
   in = addfile("/glade/u/home/lenaerts/work/CESM_LE/masks_pop.nc","r")
   nlat=dimsizes(in->TLAT)
   nlon=dimsizes(in->TLONG)
   LandFrac = get_LandFrac()
   AISMask = get_AIS_mask()
   nAISMask = doubletoint(max(AISMask))
   LandFrac = where(AISMask.gt.0, LandFrac, 0.)
   LandArea = get_LandArea()   
   in := addfile("CAM_lat_lon.nc","r")
   lat_fv09x125 = in->lat
   lon_fv09x125 = in->lon
   lev_fv09x125 = in->lev

print("Calculating integrated annual PRECIP for time series analysis and compositing.")
   in := addfile("/glade/u/home/lenaerts/work/CESM_LE/PRECIP_B1850C5CN_400-2200.nc","r")
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
   AIS_integrated_TS(:)=0.
   ; calculate basin-specific annual mass flux time series.
   do b = 1,nAISMask
      print("  ")
      print(str_concat((/"Building annual integrated precip timeseries for basin ",tostring(b),":"/)))
      bm1=b-1
      BasinMask=where(AISMask.eq.b,1,0)
     ; calculate yearly totals
      Basin_integrated_TS(:,bm1)=calculate_integrated_timeseries(PRECIP_annual,LandFrac,LandArea,BasinMask)
      Basin_integrated_TS(:,bm1)=Basin_integrated_TS(:,bm1) / 1.e9
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
  ;				 Basin_integrated_TS,Basin_integrated_TS_avg,\
  ;				  nAISMask)
   save_timeseries(Basin_integrated_TS,"BasinAccumulationTimeSeries.nc")
   save_timeseries(AIS_integrated_TS,"AISAccumulationTimeSeries.nc")
   save_timeseries(time,"TimeSeriestime.nc")

;   indt_high= NewList("lifo")
;   indt_low= NewList("lifo")    
;   do b = 1,nAISMask
;      bm1=b-1
;      TS=Basin_integrated_TS(:,bm1)
;      ListAppend(indt_high,ind(TS.gt.Basin_integrated_TS_95th_percentile(bm1)))
;      ListAppend(indt_low,ind(TS.lt.Basin_integrated_TS_5th_percentile(bm1)))
;   end do

;   if (CompositePrecipCLM) then
;      print("Compositing CLM Precipitation")
;      DataNameRAIN="RAIN"
;      DataPathRAIN="ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/RAIN/*.B1850C5CN.*nc"
;      DataNameSNOW="SNOW"
;      DataPathSNOW="ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/lnd/proc/tseries/monthly/SNOW/*.B1850C5CN.*nc"  	
;      MonthlyData=LoadMonthlyData(DataNameRAIN,DataPathRAIN)+LoadMonthlyData(DataNameSNOW,DataPathSNOW)
;      MonthlyData=replace_mask_values(MonthlyData,0.)
;      DataName="PRECIP"
;      IsOceanGrid=False
;      do b = 1,nAISMask
;	GenerateComposite("high",DataName,b,MonthlyData,indt_high(b-1),IsOceanGrid) 
;	GenerateComposite("low",DataName,b,MonthlyData,indt_low(b-1),IsOceanGrid)   
;      end do
;   end if

;   if (CompositeSLP) then
;      print("Compositing CAM SLP")
;      DataName="PSL"
;      DataPath="ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PSL/*.B1850C5CN.*.nc"
;      MonthlyData=LoadMonthlyData(DataName,DataPath)
;      IsOceanGrid=False
;      do b = 1,nAISMask
;	GenerateComposite("high",DataName,b,MonthlyData,indt_high(b-1),IsOceanGrid) 
;	GenerateComposite("low",DataName,b,MonthlyData,indt_low(b-1),IsOceanGrid)   
;      end do	
;   end if




