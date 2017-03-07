load "local_utilies.ncl"

begin
   year2sec=31557600.
   AISMask=get_AIS_mask()
   LandFrac = get_LandFrac()
   LandArea = get_LandArea()
   AISMask=where(AISMask.gt.0,1,0)
   in = addfile("/glade/u/home/lenaerts/work/CESM_LE/PRECIP_B1850C5CN_400-2200.nc","r")
   lat=in->lat
   lon=in->lon
   SNOW_annual = (in->SNOW) * year2sec / 1000.    ;mm/s -> m/yr
   SNOW_annual=replace_mask_values(SNOW_annual,0.)

   SNOW_AIS_integrated=calculate_integrated_timeseries(SNOW_annual,LandFrac,LandArea,AISMask)
   SNOW_AIS_integrated=SNOW_AIS_integrated / 1.e9 ; integrated AIS PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr
   print("Climatological integrated mean snow:")
   print(avg(SNOW_AIS_integrated))
   
   RAIN_annual = (in->RAIN) * year2sec / 1000.    ;mm/s -> m/yr
   RAIN_annual=replace_mask_values(RAIN_annual,0.)
   
   RAIN_AIS_integrated=calculate_integrated_timeseries(RAIN_annual,LandFrac,LandArea,AISMask)
   RAIN_AIS_integrated=RAIN_AIS_integrated / 1.e9 ; integrated AIS PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr
   print("Climatological integrated mean rain:")
   print(avg(RAIN_AIS_integrated))
 
   PRECIP_annual=SNOW_annual+RAIN_annual
   PRECIP_climatology=mask(dim_avg_n(PRECIP_annual,0),AISMask,1)
   PRECIP_climatology=add_lat_lon_arrays(PRECIP_climatology,lat,lon)
   PRECIP_variability=mask(dim_stddev_n(PRECIP_annual,0),AISMask,1)
   PRECIP_variability=add_lat_lon_arrays(PRECIP_variability,lat,lon)
      
   PRECIP_AIS_integrated=RAIN_AIS_integrated+SNOW_AIS_integrated
   PRECIP_AIS_integrated_avg=avg(PRECIP_AIS_integrated)
   PRECIP_AIS_integrated_std=stddev(PRECIP_AIS_integrated)
   print("Total integrated prectip (average/std:")
   print(PRECIP_AIS_integrated_avg)
   print(PRECIP_AIS_integrated_std)   
   
   wks=gsn_open_wks("png","figs/AIS_mean_precip_map")
   res=get_default_res()
   res@gsnDraw  	     = True	      ; Do not draw plot
   res@gsnFrame 	     = True	      ; Do not advance frame
   res@gsnPolar             = "SH"
   res@mpMaxLatF            = -60 ; specify min lat
   res@cnFillOn             = True
   res@cnLinesOn           = False
   res@lbTitleString       = "Precipitation (m/yr w.e.)"
   AIS_climatology_plot = gsn_csm_contour_map_polar(wks,PRECIP_climatology,res)

   delete(res)
   wks=gsn_open_wks("png","figs/AIS_std_precip_map")
   res=get_default_res()
   res@gsnDraw  	    = True	      ; Do not draw plot
   res@gsnFrame 	    = True	      ; Do not advance frame
   res@gsnPolar             = "SH"
   res@mpMaxLatF            = -60 ; specify min lat
   res@cnFillOn             = True
   res@cnLinesOn            = False
   res@lbTitleString        = "Precipitation SD (m/yr w.e.)"
   AIS_climatology_plot = gsn_csm_contour_map_polar(wks,PRECIP_variability,res)   
   
end
