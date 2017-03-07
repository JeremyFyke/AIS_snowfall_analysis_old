load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "local_utilies.ncl"

begin   

   plot_AIS_SMB_spatial_fields=0
      plot_basin_outline=0
   plot_dSIC=0
   plot_dSIC_plus_SLP =0  
   plot_dSST_global=1
   plot_dSLP_global=0
   plot_dMoisture_plus_SIC=0
   plot_dMoisture_dSLP_transport=0
   plot_dMoisture_plus_basin_correlation=0 
   plot_evap_to_transport_ratio=0
   plot_spatial_correlation=0
   plot_seasonal_cycle_difference=0

   year2sec=31557600.
   
   function plot_basin_correlation(CorrMask,wks,isPolar)
      begin
      	 res = get_default_res()  
	 res@gsnPolar             = "SH"
	 res@mpMaxLatF            = -60 
	 res@cnFillOn             = True
	 res@cnLinesOn  	 = False       ; do not draw contour lines
	 res@cnFillMode	    = "RasterFill"     ; Turn on raster color 
	 res@cnFillPalette     = "MPL_bwr"     ; set color map  
	 symMinMaxPlt (CorrMask,100,False,res)
	 if (isPolar.eq.1) then
	    correlation_plot = gsn_csm_contour_map_polar(wks,CorrMask,res)
	 else
	    correlation_plot = gsn_csm_contour(wks,CorrMask,res)	 
	 end if
         return correlation_plot
      end
   
   function plot_significance(SigField,wks)
      begin
         res=get_default_res()
	 res@cnMinLevelValF      = 0.00	     ; set min contour level
	 res@cnMaxLevelValF      = 1.00	     ; set max contour level
	 res@cnLevelSpacingF     = 0.05	     ; set contour spacing
	 res@cnLinesOn	         = False       ; do not draw contour lines
	 res@cnLineLabelsOn      = False     ; do not draw contour labels	 
	 res@cnFillScaleF        = 0.4	     ; add extra density 
	 significance_plot = gsn_csm_contour(wks,SigField,res)
	 opt		          = True	      
	 opt@gsnShadeFillType     = "pattern"
	 opt@gsnShadeLow	  = 6 
	 significance_plot = gsn_contour_shade(significance_plot, 0.05, 30., opt)  ; shade all areas >95%
         return significance_plot
      end
      
      function plot_basin(BasinMask,wks)
         begin
	    res = get_default_res() 
	    res@cnMinLevelValF      = 0.00
	    res@cnMaxLevelValF      = 1.05
	    res@cnLevelSpacingF     = 0.05
	    res@cnFillOn	    = False	
	    res@cnLineColor	    = "black"
	    res@cnMonoLineThickness = True
	    res@cnLineThicknessF = 12
	    res@cnFillMode	    = "RasterFill"    
	    res@cnLinesOn	    = True
	    res@cnLevelSelectionMode = "ManualLevels"  ; manually set the contour levels with the following 3 resources
            res@cnMinLevelValF  = 0.4		       ; set the minimum contour level
            res@cnMaxLevelValF  = 0.5		       ; set the maximum contour level
	    res@cnMaxLevelCount      = 1
	    res@cnLineLabelsOn      = False	  ; do not draw contour labels
	    res@lbLabelBarOn	    = False
	    basin_plot=gsn_csm_contour(wks,BasinMask,res)
	    return basin_plot    
         end
	 
      function plot_vector_scalar_field(u,v,s,wks,isPolar)
	 begin
	    u(0,:)=0.0 ; remove long pole vector
	    v(0,:)=0.0 ; remove long pole vector	          
	    res = get_default_res()
	    res@gsnPolar             = "SH"
	    res@mpMaxLatF            = -45
	    res@vcGlyphStyle     = "CurlyVector" 
	    res@vcRefMagnitudeF  = 1.e-8
	    res@vcRefLengthF     = 1.
	    res@vcRefAnnoOn      = False
	    res@vcMinDistanceF   = 0.01
	    res@vcLineArrowThicknessF = 3	        	
	    res@gsnScalarContour=True
	    res@cnLineThicknessF = 6
	    if (isPolar.eq.1) then
	      vector_scalar_plot=gsn_csm_vector_scalar_map_polar(wks,u(::3,::3),v(::3,::3),s,res)
	    else 
	      vector_scalar_plot=gsn_csm_vector_scalar(wks,u(::3,::3),v(::3,::3),s,res)
	    end if
	    vector_scalar_plot=ColorNegDashZeroPosContour(vector_scalar_plot,"blue","grey","red") 
	    return vector_scalar_plot
	 end

      function plot_vector_field(u,v,wks)
	 begin
	    u(0,:)=0.0 ; remove long pole vector
	    v(0,:)=0.0 ; remove long pole vector
	    res = get_default_res()
	    res@vcGlyphStyle     = "CurlyVector" 
	    res@vcRefMagnitudeF  = 1.e-8
	    res@vcRefLengthF     = 1.
	    res@vcMinDistanceF   = 0.01
	    res@vcRefAnnoOn      = False
	    res@vcLineArrowThicknessF = 3	
	    res@gsnScalarContour=True
	    vector_plot=gsn_csm_vector(wks,u(::3,::3),v(::3,::3),res) 
	    return vector_plot
	 end	    

      function plot_precip_climatology(PRECIP_climatology,wks)
	 begin
	    res=get_default_res()
	    res@gsnPolar            = "SH"
	    res@mpMaxLatF           = -60 ; specify min lat
	    res@cnFillOn            = True
	    res@cnLinesOn           = False
	    res@cnFillPalette       =  "GMT_drywet"
	    accum_plot = gsn_csm_contour_map_polar(wks,PRECIP_climatology,res)
	    return accum_plot
	 end
	 
      function plot_SIC(SIC_diff,wks)
         begin
	    res = get_default_res()   
	    res@gsnPolar             = "SH"
	    res@mpMaxLatF            = -60
	    res@cnFillOn             = True 
	    res@cnLinesOn           = True  
	    res@lbLabelAutoStride    = True
	    res@lbOrientation = "vertical"
	    res@cnMinLevelValF       = -0.2               ; set min contour level
	    res@cnMaxLevelValF       =  0.2               ; set max contour level
	    res@cnLevelSpacingF      =  0.025             ; set contour spacing  
	    res@cnFillPalette        = "cmp_b2r"          ; set color map
	    dSIC_plot = gsn_csm_contour_map_polar(wks,SIC_diff,res)
	    return dSIC_plot
	 end
	 
      function plot_SLP(SLP_diff,wks)
	 begin
	    res = get_default_res()
	    res@txFontHeightF     = 0.015
	    res@mpOutlineOn = True
	    res@mpGeophysicalLineThicknessF = 6
	    res@cnLineThicknessF = 12
	    dSLP_plot = gsn_csm_contour(wks,SLP_diff,res)  ; create the plot
	    dSLP_plot = ColorNegDashZeroPosContour(dSLP_plot,"blue","grey","red")	    
	    return dSLP_plot
	 end
	 
      function plot_SLP_for_comparison(SLP_diff,wks)
         begin
            res = get_default_res()  
	    res@cnLineLabelsOn    = False 
	    res@txFontHeightF     = 0.015
	    res@mpOutlineOn = True
	    res@mpGeophysicalLineThicknessF = 6
	    res@cnLineThicknessF = 6
	    res@cnMinLevelValF       = -1000               ; set min contour level
	    res@cnMaxLevelValF       = 1000               ; set max contour level
	    res@cnLevelSpacingF = 50. 	    
	    res@gsnPolar            = "SH"
	    res@mpMaxLatF           = -60 ; specify min lat        
	    res@cnFillOn = False
	    res@mpMinLonF              =   0
            res@mpMaxLonF              = 360
            ;res@mpCenterLonF           = 180     ; This is necessary to get the correct map
	    res@mpMinLatF              = -90
            res@mpMaxLatF              = -20 
	    dSLP_plot=gsn_csm_contour_map(wks,SLP_diff,res)
	    dSLP_plot = ColorNegDashZeroPosContour(dSLP_plot,"blue","grey","red")
	    return dSLP_plot
	  end
	  
       function plot_SST_for_comparison(SST_diff,wks)
          begin	  
	    res = get_default_res() 
	    res@mpPerimOn         = False 
	    res@cnLineLabelsOn    = False 
	    res@txFontHeightF     = 0.015
	    res@mpOutlineOn = True
	    res@mpGeophysicalLineThicknessF = 6
	    res@cnLineThicknessF = 6
	    res@mpMinLonF              =   0
            res@mpMaxLonF              = 360
            res@mpCenterLonF           = 180     ; This is necessary to get the correct map
	    res@mpMinLatF              =   -90
            res@mpMaxLatF              = 90
	    res@mpLandFillColor = "white"
	    res@mpProjection = "Robinson"	 	 
	    res@cnFillDrawOrder  = "PreDraw"  ; draw contours first
            res@cnFillOn = True
            res@cnFillPalette = "WhBlGrYeRe"
	    res@cnLinesOn  	 = False
	    res@cnMinLevelValF       = -0.4               ; set min contour level
	    res@cnMaxLevelValF       =  0.4               ; set max contour level
	    res@cnLevelSpacingF = 0.05
	    dSST_plot = gsn_csm_contour_map(wks,SST_diff,res)
	    return dSST_plot	  
	  end


   year2sec=31557600.
   wks_type=get_default_wks()   
   default_res = get_default_res()
   AISMask = get_AIS_mask()
   nAISMask = max(AISMask)
   Area=get_LandArea() ;m^2
   LandFrac = get_LandFrac()
   LandArea = get_LandArea() 
   SD="./"
   
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
      
      print("Climatological integrated mean snow:")
      print(avg(SNOW_AIS_integrated))
      print("Climatological integrated mean rain:")
      print(avg(RAIN_AIS_integrated))      
      print("Total integrated prectip (average/std:")
      print(PRECIP_AIS_integrated_avg)
      print(PRECIP_AIS_integrated_std)   

      wks=gsn_open_wks(wks_type,"figs/AIS_mean_precip_map")
      AIS_climatology_plot = plot_precip_climatology(PRECIP_climatology,wks)
      draw (AIS_climatology_plot)
      frame(wks)      

      wks=gsn_open_wks(wks_type,"figs/AIS_std_precip_map")
      AIS_variability_plot = plot_precip_climatology(PRECIP_variability,wks) 
      draw (AIS_variability_plot)
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
   lat=in->lat
   lon=in->lon
   sam_eof=add_lat_lon_arrays(sam_eof,lat,lon)
   psa1_eof=add_lat_lon_arrays(psa1_eof,lat,lon)

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
   
   do b=1,nAISMask
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
      UQ_VERT_INTEGRAL_diff=fhigh->annual_climo_UQ_VERT_INTEGRAL - flow->annual_climo_UQ_VERT_INTEGRAL
      VQ_VERT_INTEGRAL_diff=fhigh->annual_climo_VQ_VERT_INTEGRAL - flow->annual_climo_VQ_VERT_INTEGRAL
      UQ_VERT_INTEGRAL_diff=UQ_VERT_INTEGRAL_diff/doubletofloat(Area)
      VQ_VERT_INTEGRAL_diff=VQ_VERT_INTEGRAL_diff/doubletofloat(Area)
      UQ_VERT_INTEGRAL_diff=add_lat_lon_arrays(UQ_VERT_INTEGRAL_diff,lat,lon) 
      VQ_VERT_INTEGRAL_diff=add_lat_lon_arrays(VQ_VERT_INTEGRAL_diff,lat,lon)
      
      ;Calculate integrated evaporation diff
      EVAP_diff=fhigh->annual_climo_EVAP - flow->annual_climo_EVAP
      EVAP_diff=-EVAP_diff ; so increased evaporation is positive
      EVAP_diff=EVAP_diff*doubletofloat(Area)
      EVAP_diff= add_lat_lon_arrays(EVAP_diff,lat,lon)

      SST_diff=fhigh->annual_climo_SST - flow->annual_climo_SST
      SST_diff= add_lat_lon_arrays(SST_diff,lat,lon)

      ;Get monthly precipitation climatologies
      PRECIPm_high=(fhigh->monthly_climo_PRECIP) * year2sec / 1000.    ;mm/s -> m/yr
      PRECIPm_low =(flow->monthly_climo_PRECIP)  * year2sec / 1000.    ;mm/s -> m/yr     
      
      BasinMask=new(dimsizes(AISMask),"integer")
      AISBinaryMask=new(dimsizes(AISMask),"integer")
      BasinMask=where(AISMask.eq.b,1,0)
      AISBinaryMask=where(AISMask.eq.0,0,1)
      BasinMask=mask(BasinMask,AISBinaryMask,1)
      BasinMask=add_lat_lon_arrays(BasinMask,lat,lon)	

      ;;;Calculate basin-integrated annual area-integrated precipitation correlations
      in = addfile("compositing_work/output/AISAccumulationTimeSeries.nc","r")
      AISAccum=in->time_series
      delete(in)
      in = addfile("compositing_work/output/BasinAccumulationTimeSeries.nc","r")
      BasinAccum=in->time_series
      delete(in)
      in = addfile("compositing_work/output/TimeSeriestime.nc","r")
      time=in->time_series
      delete(in)
      CorrMask=new(dimsizes(AISMask),"double")
      CorrMask(:,:)=0.0
      do bc=1,nAISMask ;loop over all basins, color in based on correlation at that basin
	 Xcor=escorc(BasinAccum(:,b-1),BasinAccum(:,bc-1))
	 CorrMask=where(AISMask.eq.bc,Xcor,CorrMask)
      end do
      Rtest=rtest(CorrMask,dimsizes(time),0)
      CorrMask=mask(CorrMask,AISBinaryMask,1)
      CorrMask=add_lat_lon_arrays(CorrMask,lat,lon)
      Rtest=	mask(Rtest,AISBinaryMask,1)
      Rtest=add_lat_lon_arrays(Rtest,lat,lon)
      siglvl = 0.05
      Rtest=	mask(Rtest,AISBinaryMask,1)
      Rtest=add_lat_lon_arrays(Rtest,lat,lon)
      
      if(plot_dMoisture_dSLP_transport.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport_dSLP",sprinti("%0.2i", b)/)))
	 moisture_plot=plot_vector_scalar_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,SLP_diff,wks,0)
	 correlation_plot=plot_basin_correlation(CorrMask,wks,1)
	 basin_plot=plot_basin(BasinMask,wks)
	 overlay (correlation_plot,basin_plot)
	 overlay(correlation_plot,moisture_plot)
	 draw(correlation_plot)
	 frame(wks)	    
      end if
      
      if(plot_dMoisture_plus_basin_correlation.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport_plus_basin_correlation",sprinti("%0.2i", b)/)))
	 moisture_plot=plot_vector_field(UQ_VERT_INTEGRAL_diff,VQ_VERT_INTEGRAL_diff,wks)
	 correlation_plot=plot_basin_correlation(CorrMask,wks,1)
	 basin_plot=plot_basin(BasinMask,wks)
	 overlay(correlation_plot,moisture_plot)
	 overlay (correlation_plot,basin_plot)
	 draw(correlation_plot)
	 frame(wks)	    
      end if
      
      if(plot_dSIC_plus_SLP.eq.1) then    
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSIC_plus_dSLP",sprinti("%0.2i", b)/)))
         dSIC_plot=plot_SIC(SIC_diff,wks)
	 correlation_plot=plot_basin_correlation(CorrMask,wks,0)
	 basin_plot=plot_basin(BasinMask,wks)
	 dSLP_plot=plot_SLP(SLP_diff,wks)
         overlay(dSIC_plot,correlation_plot)
         overlay(dSIC_plot,basin_plot) 
         overlay(dSIC_plot,dSLP_plot)	 
	 draw(dSIC_plot)
	 frame(wks)	      
      end if  
      
      if(plot_dSIC.eq.1) then    
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSIC",sprinti("%0.2i", b)/)))
         dSIC_plot=plot_SIC(SIC_diff,wks)
	 correlation_plot=plot_basin_correlation(CorrMask,wks,0)
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
	 correlation_plot=plot_basin_correlation(CorrMask,wks,0)
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
	 delete(res)  
	 draw (dEVAP_plot)
	 frame(wks)      
         ;TO BE CONTINUED!
      end if

      if (plot_spatial_correlation.eq.1) then
	 
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/Correlation",sprinti("%0.2i", b)/)))
	 correlation_plot  = plot_basin_correlation(CorrMask,wks,1)
	 significance_plot = plot_significance(Rtest,wks)
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
      
      if (plot_basin_outline.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/Basin_outline",sprinti("%0.2i", b)/)))
	 climatology_plot=plot_precip_climatology(PRECIP_climatology,wks)
	 basin_plot=plot_basin(BasinMask,wks)
	 overlay(climatology_plot,basin_plot)
	 draw(climatology_plot)
	 frame(wks)	 
      end if

   end do
   
end
