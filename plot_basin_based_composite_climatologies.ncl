load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "local_utilies.ncl"
load "plotting_utilities.ncl"

begin   

   plot_AIS_SMB_spatial_fields=0
   plot_AISMask=0
   plot_dSIC=0
   plot_dSIC_plus_dZ500=0
   plot_dSST_global=0
   plot_dSLP_global=0
   plot_dMoisture_plus_SIC=0
   plot_dMoisture_dZ500_transport=1
   plot_non_annular_difference=0
   plot_dMoisture_plus_basin_correlation=0
   plot_evap_to_transport_ratio=0
   plot_spatial_correlation_CESM=0
   plot_matrix_of_CESM_interbasin_correlations=0
   plot_spatial_correlation_RACMO=0
   plot_matrix_of_RACMO_interbasin_correlations=0
   plot_spatial_correlation_PAGES=0
   plot_matrix_of_PAGES_interbasin_correlations=0 
   plot_seasonal_cycle_difference=0
   plot_non_annular_difference=0
   plot_ice_core_correlations=0
   plot_ice_core_variability=0  
   plot_basin_overview=0
   
   year2sec=31557600.
   wks_type=get_default_wks()   
   default_res = get_default_res()
   AISMask = get_AIS_mask()
   nAISMask = max(AISMask)
   Area=get_LandArea() ;m^2
   LandFrac = get_LandFrac()
   LandArea = get_LandArea() 
   SD="./"
   
   in = addfile("/glade/u/home/lenaerts/work/CESM_LE/PRECIP_B1850C5CN_400-2200.nc","r")
   lat=in->lat
   lon=in->lon    
   ZwallyBasinPolygonData=readAsciiTable("AIS_basin_masks/Ant_Full_DrainageSystem_Polygons.txt",3,"float",7)  
   
   ;;;Load RACMO AIS data
   infile = addfile("AIS_basin_masks/AIS_Full_basins_Zwally_ANT27grid.nc","r")
   AISMask_RACMO=infile->Zwallybasins
   infile = addfile("jer_data_dir/smb_monthly_ANT27_ERAIN_r490_1979_2015.nc","r")
   PRECIP_annual_RACMO=month_to_annual(infile->smb,1)
   PRECIP_annual_RACMO=dtrend_n(PRECIP_annual_RACMO,False,0)
   PRECIP_climatology_RACMO=dim_avg_n(PRECIP_annual_RACMO,0)
   PRECIP_variability_RACMO=dim_stddev_n(PRECIP_annual_RACMO,0)
   PRECIP_CV_RACMO=PRECIP_variability_RACMO/where(PRECIP_climatology_RACMO.ne.0.0,PRECIP_climatology_RACMO,PRECIP_climatology_RACMO@_FillValue); avoid divide by zer
   
   ds_RACMO=dimsizes(PRECIP_annual_RACMO)
   nt_RACMO=ds_RACMO(0)
   infile = addfile("jer_data_dir/maskers_RACMO2.nc","r")
   LandFrac_RACMO=infile->mask2d
   infile = addfile("jer_data_dir/gridarea_RACMO2.nc","r")
   LandArea_RACMO=infile->cell_area ;m^2  
   time_RACMO=fspan(1,nt_RACMO,nt_RACMO)
   Basin_integrated_TS_RACMO = new((/nt_RACMO,toint(nAISMask)/),"double")
   do b = 1,nAISMask
      bm1=b-1
      BasinMask_RACMO=where(AISMask_RACMO.eq.b,1,0)
      Basin_integrated_TS_RACMO(:,bm1)=calculate_integrated_timeseries(PRECIP_annual_RACMO,LandFrac_RACMO,LandArea_RACMO,BasinMask_RACMO)
   end do

   ;Load ice core data, stack by basin
   out:=composite_PAGES_data(lat,lon,AISMask,nAISMask)
   PAGES_core_lat=out[0]
   PAGES_core_lon=out[1]
   time_PAGES=out[3]
   Basin_integrated_TS_PAGES=out[5]
      
   if (plot_AISMask.eq.1) then
      wks = gsn_open_wks (wks_type,"AISMask")	 	 
      plot  = plot_AISMask_overview(AISMask,wks,1,-60)
      plot  = plot_core_locations(plot,wks,PAGES_core_lat,PAGES_core_lon) 	 
      draw(plot)
      frame(wks)
   end if
   
   if (plot_AIS_SMB_spatial_fields.eq.1)
      AISMaskBinary=where(AISMask.gt.0,1,0)
      in = addfile("/glade/u/home/lenaerts/work/CESM_LE/PRECIP_B1850C5CN_400-2200.nc","r")
      lat=in->lat
      lon=in->lon
      SNOW_annual = (in->SNOW) * year2sec / 1000.    ;mm/s -> m/yr
      SNOW_annual=replace_mask_values(SNOW_annual,0.)
      SNOW_AIS_integrated=calculate_integrated_timeseries(SNOW_annual,LandFrac,LandArea,AISMaskBinary)
      SNOW_AIS_integrated=SNOW_AIS_integrated / 1.e9 ; integrated AIS PRECIP - m^3/yr -> kg/yr -> Gt/yr
      RAIN_annual = (in->RAIN) * year2sec / 1000.    ;mm/s -> m/yr
      RAIN_annual=replace_mask_values(RAIN_annual,0.)
      RAIN_AIS_integrated=calculate_integrated_timeseries(RAIN_annual,LandFrac,LandArea,AISMaskBinary)
      RAIN_AIS_integrated=RAIN_AIS_integrated / 1.e9 ; integrated AIS PRECIP - m^3/yr -> kg/yr -> Gt/yr

      PRECIP_annual=SNOW_annual+RAIN_annual
      PRECIP_climatology=mask(dim_avg_n(PRECIP_annual,0),AISMaskBinary,1)
      PRECIP_climatology=add_lat_lon_arrays(PRECIP_climatology,lat,lon)
      PRECIP_variability=mask(dim_stddev_n(PRECIP_annual,0),AISMaskBinary,1)
      PRECIP_variability=add_lat_lon_arrays(PRECIP_variability,lat,lon)
      PRECIP_AIS_integrated=RAIN_AIS_integrated+SNOW_AIS_integrated
      PRECIP_AIS_integrated_avg=avg(PRECIP_AIS_integrated)
      PRECIP_AIS_integrated_std=stddev(PRECIP_AIS_integrated)
      PRECIP_CV=PRECIP_variability/PRECIP_climatology
      PRECIP_CV=add_lat_lon_arrays(PRECIP_CV,lat,lon)
      
      print("CV averaged from point-wise CVs="+wgt_areaave2(PRECIP_CV,Area,0))
      print("Climatological AIS integrated mean snow:")
      print(avg(SNOW_AIS_integrated))
      print("Climatological AIS integrated mean rain:")
      print(avg(RAIN_AIS_integrated))      
      print("Total integrated AIS precip (average/std:")
      print(PRECIP_AIS_integrated_avg)
      print(PRECIP_AIS_integrated_std)   

      wks=gsn_open_wks(wks_type,"figs/AIS_mean_precip_map")
      cnLevels = array_append_record(fspan(0.0,0.3,31),fspan(0.4,1.2,9),0)
      AIS_climatology_plot = plot_precip_climatology(PRECIP_climatology,cnLevels,wks,1)
      AIS_climatology_plot = plot_basin_polygons(ZwallyBasinPolygonData,AIS_climatology_plot,wks,nAISMask)
      draw (AIS_climatology_plot)
      frame(wks)      

      wks=gsn_open_wks(wks_type,"figs/AIS_std_precip_map")
      cnLevels := array_append_record(fspan(0.0,0.04,21),fspan(0.05,0.22,18),0)
      AIS_variability_plot = plot_precip_climatology(PRECIP_variability,cnLevels,wks,1)
      AIS_variability_plot = plot_basin_polygons(ZwallyBasinPolygonData,AIS_variability_plot,wks,nAISMask)      
      draw (AIS_variability_plot)
      frame(wks)
      
      wks=gsn_open_wks(wks_type,"figs/AIS_CV_map")
      cnLevels := fspan(0.1,0.37,37)      
      AIS_CV_plot = plot_precip_climatology(PRECIP_CV,cnLevels,wks,1)
      AIS_CV_plot = plot_basin_polygons(ZwallyBasinPolygonData,AIS_CV_plot,wks,nAISMask)      
      draw (AIS_CV_plot)
      frame(wks)
      
      wks=gsn_open_wks(wks_type,"figs/AIS_CV_map_RACMO")
      cnLevels := fspan(0.1,0.37,37)      
      AIS_CV_plot = plot_precip_climatology(PRECIP_CV_RACMO,cnLevels,wks,0)     
      draw (AIS_CV_plot)
      frame(wks)            
   end if
   
   ;Load CVDP diagnostics for CESM LE run
   in=addfile("/glade/p/cesmLE/CESM-CAM5-BGC-LE/CVDP/Controls/CESM1_LENS_Coupled_Control.cvdp_data.401-2200.nc","r")
   nino34=month_to_annual(in->nino34,1)
   sam=in->sam_pc_ann
   psa1=in->psa1_pc_ann
   psa2=in->psa2_pc_ann
   sam_eof=-(in->sam_ann)*100. ; convert to Pa from hPa, flip sign
   psa1_eof=-(in->psa1_sst_regression_ann) ;
   lat_CVDP=in->lat
   lon_CVDP=in->lon
   sam_eof=add_lat_lon_arrays(sam_eof,lat_CVDP,lon_CVDP)
   psa1_eof=add_lat_lon_arrays(psa1_eof,lat_CVDP,lon_CVDP)

   if (plot_dSLP_global.eq.1) then
      ;plot projection of SAM on SLP
      wks = gsn_open_wks (wks_type,"figs/SAM_eof")
      SAM_eof_plot=plot_SLP_for_comparison(sam_eof,wks)
      draw (SAM_eof_plot)
      frame(wks)       
   end if 
   if (plot_dSST_global.eq.1) then
      ;plot projection of PSA1 on SST
      wks = gsn_open_wks (wks_type,"figs/PSA1_eof")
      PSA1_eof_plot=plot_SST_for_comparison(psa1_eof,wks)
      draw (PSA1_eof_plot)
      frame(wks)      
   end if

   ;Load basin accumulation time series
   in = addfile("compositing_work/output/AISAccumulationTimeSeries.nc","r")
   AISAccum=in->time_series
   in := addfile("compositing_work/output/BasinAccumulationTimeSeries.nc","r")
   BasinAccum=in->time_series
   in := addfile("compositing_work/output/TimeSeriestime.nc","r")
   time=in->time_series
   BasinAccum!0="time"
   BasinAccum!1="basin"
   ds=dimsizes(BasinAccum)
   cc=escorc(AISAccum,BasinAccum(basin|:,time|:))
   
   ;Figure out which basins tend to contribute most to AIS-wide integrated accumulation anomaly
   BasinImportance=new(ds(1),"integer")
   BasinImportance(:)=0 
   ;get anomaly in total AIS snowfall
   AISAnom=dim_rmvmean(AISAccum)
   BasinAnoms=dim_rmvmean_n(BasinAccum,0)
   BasinAnoms!0="time"
   BasinAnoms!1="basin"   
   cc=escorc(AISAnom,BasinAnoms(basin|:,time|:))
   do t=1,ds(0)
      tm1=t-1
      AbsDif=abs(BasinAnoms(tm1,:)-AISAnom(tm1))
      ib=minind(AbsDif)
      BasinImportance(ib)=BasinImportance(ib)+1
   end do
   do b=1,ds(1)
     print("***Basin "+b+" Contribution to total AIS variability***")
     print("    #of years where basin leads in contribution to AIS variability: "+BasinImportance(b-1))
     print("    Correlation between basin variability and AIS variability: "+cc(b-1))
   end do

;;;Calculate basin-specific annual area-integrated precipitation correlations
   out:=calculate_interbasin_correlations(BasinAccum,nAISMask)
   Xcor_basin_by_basin_CESM=out[4]
   Rtest_basin_by_basin_CESM=out[5]
   print("***CESM: Average number of total positive correlations per basin="+avg(out[0]))   
   print("***CESM: Average number of significant positive correlations per basin="+avg(out[2]))
   print("***CESM: Average number of total negative correlations per basin="+avg(out[1]))   
   print("***CESM: Average number of significant negative correlations per basin="+avg(out[3]))   

   out:=calculate_interbasin_correlations(Basin_integrated_TS_RACMO,nAISMask)
   Xcor_basin_by_basin_RACMO=out[4]
   Rtest_basin_by_basin_RACMO=out[5]
   print("***RACMO: Average number of total positive correlations per basin="+avg(out[0]))   
   print("***RACMO: Average number of significant positive correlations per basin="+avg(out[2]))
   print("***RACMO: Average number of total negative correlations per basin="+avg(out[1]))   
   print("***RACMO: Average number of significant negative correlations per basin="+avg(out[3]))  

   out:=calculate_interbasin_correlations(Basin_integrated_TS_PAGES,nAISMask)
   Xcor_basin_by_basin_PAGES=out[4]
   Rtest_basin_by_basin_PAGES=out[5]
   print("***PAGES: Average number of total positive correlations per basin="+avg(out[0]))   
   print("***PAGES: Average number of significant positive correlations per basin="+avg(out[2]))
   print("***PAGES: Average number of total negative correlations per basin="+avg(out[1]))   
   print("***PAGES: Average number of significant negative correlations per basin="+avg(out[3]))
   
   if (plot_matrix_of_CESM_interbasin_correlations.eq.1) then
       Xcor_basin_by_basin_CESM_copy=Xcor_basin_by_basin_CESM
       Rtest_basin_by_basin_CESM_copy=Rtest_basin_by_basin_CESM
       do b=1,toint(nAISMask)
          endm1=toint(nAISMask)-1
          bm1=b-1
	  Xcor_basin_by_basin_CESM_copy(bm1,bm1:endm1)=Xcor_basin_by_basin_CESM@_FillValue
	  Rtest_basin_by_basin_CESM_copy(bm1,bm1:endm1)=Rtest_basin_by_basin_CESM@_FillValue	  
       end do
       wks = gsn_open_wks (wks_type,"figs/CESM_inter_basin_correlation_matrix")            
       corr_mat_plot=plot_basin_correlation_matrix(Xcor_basin_by_basin_CESM_copy,wks)
       significance_plot=plot_significance(Rtest_basin_by_basin_CESM_copy,wks)
       overlay(corr_mat_plot,significance_plot)
       draw(corr_mat_plot)
       frame(wks)       
   end if
   
   if (plot_matrix_of_RACMO_interbasin_correlations.eq.1) then
       Xcor_basin_by_basin_RACMO_copy=Xcor_basin_by_basin_RACMO
       Rtest_basin_by_basin_RACMO_copy=Rtest_basin_by_basin_RACMO
       do b=1,toint(nAISMask)
          endm1=toint(nAISMask)-1
          bm1=b-1
	  Xcor_basin_by_basin_RACMO_copy(bm1,bm1:endm1)=Xcor_basin_by_basin_RACMO@_FillValue
	  Rtest_basin_by_basin_RACMO_copy(bm1,bm1:endm1)=Rtest_basin_by_basin_RACMO@_FillValue	  
       end do
       wks = gsn_open_wks (wks_type,"figs/RACMO_inter_basin_correlation_matrix")            
       corr_mat_plot=plot_basin_correlation_matrix(Xcor_basin_by_basin_RACMO_copy,wks)
       significance_plot=plot_significance(Rtest_basin_by_basin_RACMO_copy,wks)
       overlay(corr_mat_plot,significance_plot)
       draw(corr_mat_plot)
       frame(wks)       
   end if
   
   if (plot_matrix_of_PAGES_interbasin_correlations.eq.1) then
       Xcor_basin_by_basin_PAGES_copy=Xcor_basin_by_basin_PAGES
       Rtest_basin_by_basin_PAGES_copy=Rtest_basin_by_basin_PAGES
       do b=1,toint(nAISMask)
          endm1=toint(nAISMask)-1
          bm1=b-1
	  Xcor_basin_by_basin_PAGES_copy(bm1,bm1:endm1)=Xcor_basin_by_basin_PAGES@_FillValue
	  Rtest_basin_by_basin_PAGES_copy(bm1,bm1:endm1)=Rtest_basin_by_basin_PAGES@_FillValue	  
       end do
       wks = gsn_open_wks (wks_type,"figs/PAGES_inter_basin_correlation_matrix")            
       corr_mat_plot=plot_basin_correlation_matrix(Xcor_basin_by_basin_PAGES_copy,wks)
       significance_plot=plot_significance(Rtest_basin_by_basin_PAGES_copy,wks)
       overlay(corr_mat_plot,significance_plot)
       draw(corr_mat_plot)
       frame(wks)       
   end if   
        
   do b=1,nAISMask
      bm1=b-1
      
      ;;;Load basin-specific composite data
      flow = addfile(str_concat((/SD,"compositing_work/output/low/Basin",tostring(b),"_ssticetemp_fv0.9x1.25.nc"/)),"r")
      fmean = addfile(str_concat((/SD,"compositing_work/output/mean/Basin",tostring(b),"_ssticetemp_fv0.9x1.25.nc"/)),"r") 
      fhigh = addfile(str_concat((/SD,"compositing_work/output/high/Basin",tostring(b),"_ssticetemp_fv0.9x1.25.nc"/)),"r")
      SIC_mean_low =  flow->annual_climo_ICEFRAC
      SIC_mean_high = fhigh->annual_climo_ICEFRAC
      SIC_var_low = flow->annual_var_ICEFRAC
      SIC_var_high = fhigh->annual_var_ICEFRAC
      SLP_mean_low =  flow->annual_climo_SLP
      SLP_mean_high = fhigh->annual_climo_SLP
      SLP_var_low = flow->annual_var_SLP
      SLP_var_high = fhigh->annual_var_SLP
      nSamples_low = flow->nSamples
      nSamples_high = fhigh->nSamples
      lat=flow->lat
      lon=flow->lon      
      
      ;Calculate SIC difference
      SIC_diff=SIC_mean_high-SIC_mean_low  
      SIC_diff=add_lat_lon_arrays(SIC_diff,lat,lon)     
      ;Calculate significance of difference
      SIC_alpha = ttest(SIC_mean_low,SIC_var_low,nSamples_low,SIC_mean_high,SIC_var_high,nSamples_high, False, False)
      SIC_alpha = add_lat_lon_arrays(SIC_alpha,lat,lon)
      ;Calculate SLP diff
      SLP_diff= SLP_mean_high-SLP_mean_low 
      SLP_diff= add_lat_lon_arrays(SLP_diff,lat,lon)
      ;Calculate significance of difference
      SLP_alpha = ttest(SLP_mean_low,SLP_var_low,nSamples_low,SLP_mean_high,SLP_var_high,nSamples_high, False, False)
      SLP_alpha = add_lat_lon_arrays(SLP_alpha,lat,lon)   
      ;Calculate moisture transport diff      
      UQ_VERT_INTEGRAL_mean_low =  flow->annual_climo_UQ_VERT_INTEGRAL
      UQ_VERT_INTEGRAL_mean_high = fhigh->annual_climo_UQ_VERT_INTEGRAL
      UQ_VERT_INTEGRAL_var_low = flow->annual_var_UQ_VERT_INTEGRAL
      UQ_VERT_INTEGRAL_var_high = fhigh->annual_var_UQ_VERT_INTEGRAL
      VQ_VERT_INTEGRAL_mean_low =  flow->annual_climo_VQ_VERT_INTEGRAL
      VQ_VERT_INTEGRAL_mean_high = fhigh->annual_climo_VQ_VERT_INTEGRAL
      VQ_VERT_INTEGRAL_var_low = flow->annual_var_VQ_VERT_INTEGRAL
      VQ_VERT_INTEGRAL_var_high = fhigh->annual_var_VQ_VERT_INTEGRAL      
      UQ_VERT_INTEGRAL_diff=UQ_VERT_INTEGRAL_mean_high - UQ_VERT_INTEGRAL_mean_low
      VQ_VERT_INTEGRAL_diff=VQ_VERT_INTEGRAL_mean_high - VQ_VERT_INTEGRAL_mean_low
      UQ_VERT_INTEGRAL_diff=UQ_VERT_INTEGRAL_diff/doubletofloat(Area)
      VQ_VERT_INTEGRAL_diff=VQ_VERT_INTEGRAL_diff/doubletofloat(Area)
      UQ_VERT_INTEGRAL_diff=add_lat_lon_arrays(UQ_VERT_INTEGRAL_diff,lat,lon) 
      VQ_VERT_INTEGRAL_diff=add_lat_lon_arrays(VQ_VERT_INTEGRAL_diff,lat,lon)
      UQ_VERT_INTEGRAL_alpha=ttest(UQ_VERT_INTEGRAL_mean_low,UQ_VERT_INTEGRAL_var_low,nSamples_low,\
                       UQ_VERT_INTEGRAL_mean_high,UQ_VERT_INTEGRAL_var_high,nSamples_high,\
		       False, False)
      UQ_VERT_INTEGRAL_alpha=add_lat_lon_arrays(UQ_VERT_INTEGRAL_alpha,lat,lon)           
      VQ_VERT_INTEGRAL_alpha=ttest(VQ_VERT_INTEGRAL_mean_low,VQ_VERT_INTEGRAL_var_low,nSamples_low,\
                       VQ_VERT_INTEGRAL_mean_high,VQ_VERT_INTEGRAL_var_high,nSamples_high,\
		       False, False)		       
      VQ_VERT_INTEGRAL_alpha=add_lat_lon_arrays(VQ_VERT_INTEGRAL_alpha,lat,lon) 
           
      ;Calculate integrated evaporation diff
      EVAP_diff=fhigh->annual_climo_EVAP - flow->annual_climo_EVAP
      EVAP_diff=-EVAP_diff ; so increased evaporation is positive
      EVAP_diff=EVAP_diff*doubletofloat(Area)
      EVAP_diff= add_lat_lon_arrays(EVAP_diff,lat,lon)
      SST_diff=fhigh->annual_climo_SST - flow->annual_climo_SST
      SST_diff= add_lat_lon_arrays(SST_diff,lat,lon)
      Z500_mean_low =  flow->annual_climo_Z500
      Z500_mean_high = fhigh->annual_climo_Z500
      Z500Z3_low = dim_rmvmean_Wrap(Z500_mean_low)
      Z500Z3_high = dim_rmvmean_Wrap(Z500_mean_high)
      dZ500Z3 = Z500Z3_high - Z500Z3_low
      dZ500Z3= add_lat_lon_arrays(dZ500Z3,lat,lon)      
      Z500_var_low = flow->annual_var_Z500
      Z500_var_high = fhigh->annual_var_Z500      
      Z500_diff = Z500_mean_high-Z500_mean_low
      Z500_diff= add_lat_lon_arrays(Z500_diff,lat,lon)
      Z500_alpha=ttest(Z500_mean_low,Z500_var_low,nSamples_low,\
                       Z500_mean_high,Z500_var_high,nSamples_high,\
		       False, False)
      Z500_alpha=add_lat_lon_arrays(Z500_alpha,lat,lon)		       
      
      rc=regline_stats(ndtooned(Z500_diff),ndtooned(SLP_diff))
      print("***Relationship between Z500 and SLP variability for basin "+b+"***")
          print("regression coefficient: "+rc@r)

      ;Get monthly precipitation climatologies
      PRECIPm_high=(fhigh->monthly_climo_PRECIP) * year2sec / 1000.    ;mm/s -> m/yr
      PRECIPm_low =(flow->monthly_climo_PRECIP)  * year2sec / 1000.    ;mm/s -> m/yr   

      BasinMask=new(dimsizes(AISMask),"integer")      
      BasinMask=where(AISMask.eq.b,1,0);set basinmask to one where basin equals 1
      BasinMask=add_lat_lon_arrays(BasinMask,lat,lon)
      out:=make_correlation_mask(AISMask,nAISMask,Xcor_basin_by_basin_CESM,time,bm1)
      CorrMask_CESM=add_lat_lon_arrays(out[0],lat,lon)
      RtestMask_CESM=add_lat_lon_arrays(out[1],lat,lon)
      out:=make_correlation_mask(AISMask_RACMO,nAISMask,Xcor_basin_by_basin_RACMO,time_RACMO,bm1)
      CorrMask_RACMO=out[0]
      RtestMask_RACMO=out[1]
      out:=make_correlation_mask(AISMask,nAISMask,Xcor_basin_by_basin_PAGES,time_PAGES,bm1)
      CorrMask_PAGES=add_lat_lon_arrays(out[0],lat,lon)
      RtestMask_PAGES=add_lat_lon_arrays(out[1],lat,lon)    
      
       ;;;;Various basin-specific plots!;;;;

      if(plot_dMoisture_dZ500_transport.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport_dZ500",sprinti("%0.2i", b)/)))
	 moisture_plot=plot_vector_scalar_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,Z500_diff,wks,0,-30)
	 sigfield=where(UQ_VERT_INTEGRAL_alpha.gt.VQ_VERT_INTEGRAL_alpha ,\
	                UQ_VERT_INTEGRAL_alpha,VQ_VERT_INTEGRAL_alpha)
	 sigfield=add_lat_lon_arrays(Z500_alpha,lat,lon)
	 significance_plot=plot_significance(sigfield,wks)
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,1,-30)
	 basin_plot=plot_basin(BasinMask,wks)
	 overlay (correlation_plot,basin_plot)
	 overlay(correlation_plot,moisture_plot)
	 overlay(correlation_plot,significance_plot)	 
	 draw(correlation_plot)
	 frame(wks)	    
      end if
      
      if (plot_non_annular_difference.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dz500Z3",sprinti("%0.2i", b)/)))      
         dZ500Z3_plot=plot_SLP_for_comparison(dZ500Z3,wks)
      	 overlay (dZ500Z3_plot,plot_basin(BasinMask,wks)) 
         draw (dZ500Z3_plot)
	 frame(wks)
	 
      end if
      if(plot_dMoisture_plus_basin_correlation.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport_plus_basin_correlation",sprinti("%0.2i", b)/)))
	 moisture_plot=plot_vector_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,wks)
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,1,-45)
	 basin_plot=plot_basin(BasinMask,wks)
	 overlay(correlation_plot,moisture_plot)
	 overlay (correlation_plot,basin_plot)
	 draw(correlation_plot)
	 frame(wks)	    
      end if
      
      if(plot_dSIC_plus_dZ500.eq.1) then    
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSIC_plus_dZ500",sprinti("%0.2i", b)/)))
         dSIC_plot=plot_SIC(SIC_diff,wks)
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,0,0)
	 basin_plot=plot_basin(BasinMask,wks)
	 dZ500_plot=plot_SLP(Z500_diff,wks)
         overlay(dSIC_plot,correlation_plot)
         overlay(dSIC_plot,basin_plot) 
         overlay(dSIC_plot,dZ500_plot)	 
	 draw(dSIC_plot)
	 frame(wks)	      
      end if  
      
      if(plot_dSIC.eq.1) then    
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSIC",sprinti("%0.2i", b)/)))
         dSIC_plot=plot_SIC(SIC_diff,wks)
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,0,0)
	 basin_plot=plot_basin(BasinMask,wks)
         overlay(dSIC_plot,correlation_plot)
         overlay(dSIC_plot,basin_plot)   
	 draw(dSIC_plot)
	 frame(wks)	      
      end if        
      
      if(plot_dMoisture_plus_SIC.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport_plus_SIC",sprinti("%0.2i", b)/)))
         dSIC_plot=plot_SIC(SIC_diff,wks)
         moisture_plot=plot_vector_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,wks)
	 correlation_plot=plot_basin_correlation(CorrMask_CESM,wks,0,0)
	 basin_plot=plot_basin(BasinMask,wks)
         overlay(dSIC_plot,correlation_plot)
         overlay(dSIC_plot,basin_plot)
         overlay(dSIC_plot,moisture_plot)	   
	 draw(dSIC_plot)
	 frame(wks)	 
      end if
      
      if (plot_dSST_global.eq.1) then      
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSST_global",sprinti("%0.2i", b)/)))
         dSST_plot = plot_SST_for_comparison(SST_diff,wks)
         overlay (dSST_plot,plot_basin(BasinMask,wks)) 
         draw (dSST_plot)
	 frame(wks)	 	  
      end if
      
      if (plot_dSLP_global.eq.1) then 
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSLP_global",sprinti("%0.2i", b)/)))
	 dSLP_plot = plot_SLP_for_comparison(SLP_diff,wks)
	 overlay (dSLP_plot,plot_basin(BasinMask,wks)) 
         draw (dSLP_plot)
	 frame(wks)
      end if 
      
      if (plot_evap_to_transport_ratio.eq.1) then
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/dEVAP",sprinti("%0.2i", b)/)))
	 res = default_res   
	 res@gsnPolar             = "SH"
	 res@mpMaxLatF            = -50 
	 res@tiMainString         = "EVAP Difference (high - low)"
	 res@cnFillOn             = False       
	 dEVAP_plot = gsn_csm_contour_map_polar(wks,EVAP_diff,res)  
	 draw (dEVAP_plot)
	 frame(wks)      
         ;TO BE CONTINUED!  See Burrows email, and Tsukernik 2013 for approach
      end if

      if (plot_spatial_correlation_CESM.eq.1) then
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Correlation",sprinti("%0.2i", b)/)))
	 correlation_plot  = plot_basin_correlation(CorrMask_CESM,wks,1,-60) 
	 significance_plot = plot_significance(RtestMask_CESM,wks)
	 basin_plot        = plot_basin(BasinMask,wks)
	 overlay (correlation_plot,significance_plot) 
	 overlay (correlation_plot,basin_plot)	 	 
	 draw(correlation_plot)
	 frame(wks)
      end if
      
      if (plot_spatial_correlation_RACMO.eq.1) then
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Correlation_RACMO",sprinti("%0.2i", b)/)))	 	 
	 correlation_plot  = plot_basin_correlation(CorrMask_RACMO,wks,0,0)	 
	 significance_plot = plot_significance(RtestMask_RACMO,wks)	 
	 ;basin_plot        = plot_basin(BasinMask_RACMO,wks)		  
	 overlay (correlation_plot,significance_plot) 
	 ;overlay (correlation_plot,basin_plot)	 	 
	 draw(correlation_plot)
	 frame(wks)
      end if      
      
      if (plot_spatial_correlation_PAGES.eq.1) then
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Correlation_PAGES",sprinti("%0.2i", b)/)))	 	 
	 correlation_plot  = plot_basin_correlation(CorrMask_PAGES,wks,1,-60)
	 correlation_plot  = plot_core_locations(correlation_plot,wks,PAGES_core_lat,PAGES_core_lon)
	 significance_plot = plot_significance(RtestMask_PAGES,wks)	 
	 basin_plot        = plot_basin(BasinMask,wks)		  
	 overlay (correlation_plot,significance_plot) 
	 overlay (correlation_plot,basin_plot)	 	 
	 draw(correlation_plot)
	 frame(wks)
      end if            
           
      if (plot_seasonal_cycle_difference.eq.1) then
      
         Monthly_integrated_cycle_low=calculate_integrated_timeseries(PRECIPm_low,LandFrac,LandArea,BasinMask)
	 Monthly_integrated_cycle_low=Monthly_integrated_cycle_low / 1.e9 ; integrated PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr
         Monthly_integrated_cycle_high=calculate_integrated_timeseries(PRECIPm_high,LandFrac,LandArea,BasinMask)
	 Monthly_integrated_cycle_high=Monthly_integrated_cycle_high / 1.e9 ; integrated PRECIP for basin - m^3/yr -> kg/yr -> Gt/yr

         dMonthly_integrated_cycle=Monthly_integrated_cycle_high-Monthly_integrated_cycle_low
	 
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Seasonal_cycle" ,sprinti("%0.2i", b)/)))
	 months=ispan(0,11,1)
	 res=default_res
	 res@xyLineThicknessF = 4      
	 res@trYMinF = 0.
	 res@trXMaxF = 11	 
	 res@xyLineColor      = "red"    
	 res@tmXBMode             = "Explicit"      ; label independently
	 res@tmXBValues           = months
         res@tmXBLabels           = (/"J","F","M","A","M","J",\
                                      "J","A","S","O","N","D"/)
	 res@tiYAxisString = "Integrated accumulation (Gt/yr)"
	 res@tiXAxisString = "Month"
         plot0=gsn_csm_xy (wks,months,Monthly_integrated_cycle_high,res) 
	 res@xyLineColor      = "blue" 	 
	 plot1=gsn_csm_xy (wks,months,Monthly_integrated_cycle_low,res)
	 res@xyLineColor      = "purple" 
	 res@xyLineThicknessF = 12	 	 
	 plot2=gsn_csm_xy (wks,months,dMonthly_integrated_cycle,res)	 
	 
	 overlay(plot0,plot1)
	 overlay(plot0,plot2)	 
	 draw(plot0)
	 frame(wks)

      end if

   end do
 
end
