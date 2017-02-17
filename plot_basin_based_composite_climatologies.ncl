load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "local_utilies.ncl"

begin

   plot_dSIC=0
   plot_dSLP_global=0
   plot_dMoisture_transport=0
   plot_evap_to_transport_ratio=0
   plot_timeseries_statistics=1
   
   default_res = get_default_res()
   AISMask = get_AIS_mask()
   nAISMask = max(AISMask)
   Area=get_LandArea() ;m^2

   SD="./"

   wks_type=get_default_wks()

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
      
      BasinMask=new(dimsizes(AISMask),"integer")
      AISBinaryMask=new(dimsizes(AISMask),"integer")
      BasinMask=where(AISMask.eq.b,1,0)
      AISBinaryMask=where(AISMask.eq.0,0,1)
      BasinMask=mask(BasinMask,AISBinaryMask,1)
      BasinMask=add_lat_lon_arrays(BasinMask,lat,lon)	
          
      ;Plots!
      if(plot_dSIC.eq.1) then      
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/dSIC",sprinti("%0.2i", b)/)))      
	 ;plot SIC difference
	 res = default_res   
	 res@gsnPolar             = "SH"
	 res@mpMaxLatF            = -50 
	 res@tiMainString         = "SIC Difference (high - low)"
	 res@cnFillOn             = True 
	 res@cnLinesOn           = True  
	 res@lbLabelAutoStride    = True
	 res@lbOrientation = "vertical"
	 res@cnMinLevelValF       = -0.2               ; set min contour level
	 res@cnMaxLevelValF       =  0.2               ; set max contour level
	 res@cnLevelSpacingF      =  0.025             ; set contour spacing  
	 res@cnFillPalette        = "BlWhRe"          ; set color map
	 dSIC_plot = gsn_csm_contour_map_polar(wks,SIC_diff,res) 
	 delete(res)

	 ;plot significance
	 res = default_res
	 res@cnMinLevelValF      = 0.00	     ; set min contour level
	 res@cnMaxLevelValF      = 1.05	     ; set max contour level
	 res@cnLevelSpacingF     = 0.05	     ; set contour spacing
	 res@cnLinesOn	       = False       ; do not draw contour lines
	 res@cnLineLabelsOn      = False       ; do not draw contour labels
	 res@cnFillScaleF        = 0.4	     ; add extra density 
	 significance_plot  = gsn_csm_contour(wks,SIC_alpha, res)
         delete(res)

	 res		          = True	      
	 res@gsnShadeFillType     = "pattern"
	 res@gsnShadeLow	       = 17
	 significance_plot   = gsn_contour_shade(significance_plot, 0.05, 30, res)  ; shade all areas >95%
         delete(res)

	 res = default_res    
	 res@cnMinLevelValF	 = 0.00        ; set min contour level
	 res@cnMaxLevelValF	 = 1.05        ; set max contour level
	 res@cnLevelSpacingF	 = 0.05        ; set contour spacing   
	 res@cnFillOn		 = True 	    ; Turn on color fill
	 res@cnFillMode 	 = "RasterFill"     ; Turn on raster color	
	 res@cnLinesOn  	 = False       ; do not draw contour lines
	 res@cnLineLabelsOn	 = False       ; do not draw contour labels
	 res@lbLabelBarOn	 = False      
	 res@cnFillPalette	 = "GMT_cool"	       ; set color map
	 basin_plot   = gsn_csm_contour(wks,BasinMask,res)
	 delete(res)
	 
	 res		       = True	      
	 res@gsnShadeFillType     = "pattern"
	 res@gsnShadeLow	       = 17
	 basin_plot   = gsn_contour_shade(basin_plot, 0.05, 30, res)  ; shade all areas >95% 
	 delete(res)
	      
	 res = default_res
	 res@cnFillOn	       = False
	 res@lbLabelBarOn      = False
	 SLP_plot = gsn_csm_contour(wks,SLP_diff,res)
	 SLP_plot = ColorNegDashZeroPosContour(SLP_plot,"red","black","blue")
	 delete(res)
	 
	 overlay (dSIC_plot,significance_plot)
	 overlay (dSIC_plot,basin_plot)
	 overlay (dSIC_plot,SLP_plot)      
	 draw (dSIC_plot)
	 frame(wks)
      end if
      
      if(plot_dMoisture_transport.eq.1) then
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dMoisture_transport",sprinti("%0.2i", b)/)))
	    print(max(UQ_VERT_INTEGRAL_diff))
	    res = default_res   
	    res@gsnPolar             = "SH"
	    res@mpMaxLatF            = -50
	    res@tiMainString         = "Moisture transport Difference (high - low)"
	    res@vcGlyphStyle     = "CurlyVector" 
	    res@vcRefMagnitudeF  = 1.e-8
	    res@vcRefLengthF     = 1.
            ;res@vcRefMagnitudeF  = 1	      
            ;res@vcRefLengthF     = 0.01
	    res@vcMinDistanceF   = 0.01	    	
	    res@gsnScalarContour=True
	    dMOisture_plot=gsn_csm_vector_scalar_map_polar(wks,UQ_VERT_INTEGRAL_diff(::3,::3),VQ_VERT_INTEGRAL_diff(::3,::3),SLP_diff,res)
	    dMOisture_plot = ColorNegDashZeroPosContour(dMOisture_plot,"red","grey","blue") 
	    delete(res)
	    
	    res = default_res	 
	    res@cnMinLevelValF      = 0.00	  ; set min contour level
	    res@cnMaxLevelValF      = 1.05	  ; set max contour level
	    res@cnLevelSpacingF     = 0.05	  ; set contour spacing   
	    res@cnFillOn	    = False		; Turn on color fill
	    res@cnLineColor	    = "dodgerblue"
	    res@cnMonoLineThickness = True
	    res@cnLineThicknessF = 5
	    res@cnFillMode	    = "RasterFill"     ; Turn on raster color	   
	    res@cnLinesOn	    = True	 ; do not draw contour lines
	    res@cnLevelSelectionMode = "ManualLevels"  ; manually set the contour levels with the following 3 resources
            res@cnMinLevelValF  = 0.4		       ; set the minimum contour level
            res@cnMaxLevelValF  = 0.5		       ; set the maximum contour level
	    res@cnMaxLevelCount      = 1
	    res@cnLineLabelsOn      = False	  ; do not draw contour labels
	    res@lbLabelBarOn	    = False
	    basin_plot=gsn_csm_contour(wks,BasinMask,res)
	    delete(res)
	    
	    overlay (dMOisture_plot,basin_plot)
	    draw(dMOisture_plot)
	    frame(wks)	    
      end if
      
      if (plot_dSLP_global.eq.1) then
      
         wks = gsn_open_wks (wks_type,str_concat((/"figs/dSLP_global",sprinti("%0.2i", b)/)))
      
         res = default_res 
         res@tiMainString         = "SLP Difference (high - low)"
	 res@mpPerimOn         = False             ; turn off box around plot
	 res@mpFillOn          = False
	 res@cnFillPalette     = "BlueDarkRed18"     ; set color map
	 res@cnLineLabelsOn    = False             ; turn off contour lines
	 res@txFontHeightF     = 0.015
	 res@mpOutlineOn = True
	 res@mpGeophysicalLineThicknessF = 3
	 dSLP_plot = gsn_csm_contour_map(wks,SLP_diff,res)  ; create the plot
	 dSLP_plot = ColorNegDashZeroPosContour(dSLP_plot,"red","grey","blue")
	 delete(res)
	 
	 ;plot significance
	 res = default_res
	 res@cnMinLevelValF      = 0.00	     ; set min contour level
	 res@cnMaxLevelValF      = 1.05	     ; set max contour level
	 res@cnLevelSpacingF     = 0.05	     ; set contour spacing
	 res@cnLinesOn	       = False       ; do not draw contour lines
	 res@cnLineLabelsOn      = False       ; do not draw contour labels
	 res@cnFillScaleF        = 0.4	     ; add extra density
	 significance_plot  = gsn_csm_contour(wks,SLP_alpha, res)
         delete(res)
	 
	 res		          = True	      
	 res@gsnShadeFillType     = "pattern"
	 res@gsnShadeLow	       = 17
	 significance_plot   = gsn_contour_shade(significance_plot, 0.05, 30, res)  ; shade all areas >95%	 
         delete(res)
	 
	 res = default_res	 
	 res@cnMinLevelValF      = 0.00	  ; set min contour level
	 res@cnMaxLevelValF      = 1.05	  ; set max contour level
	 res@cnLevelSpacingF     = 0.05	  ; set contour spacing   
	 res@cnFillOn	    = False		; Turn on color fill
	 res@cnLineColor	    = "dodgerblue"
	 res@cnMonoLineThickness = True
	 res@cnLineThicknessF = 5
	 res@cnFillMode	    = "RasterFill"     ; Turn on raster color	   
	 res@cnLinesOn	    = True	 ; do not draw contour lines
	 res@cnLevelSelectionMode = "ManualLevels"  ; manually set the contour levels with the following 3 resources
         res@cnMinLevelValF  = 0.4		       ; set the minimum contour level
         res@cnMaxLevelValF  = 0.5		       ; set the maximum contour level
	 res@cnMaxLevelCount      = 1
	 res@cnLineLabelsOn      = False	  ; do not draw contour labels
	 res@lbLabelBarOn	    = False
	 basin_plot=gsn_csm_contour(wks,BasinMask,res)
	 delete(res)	 
	    
	 overlay (dSLP_plot,significance_plot)
	 overlay (dSLP_plot,basin_plot)  
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

      end if

      if (plot_timeseries_statistics.eq.1) then
      
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
	 
	 siglvl = 0.05
	 Rtest=	mask(Rtest,AISBinaryMask,1)
	 Rtest=add_lat_lon_arrays(Rtest,lat,lon)
	 wks = gsn_open_wks (wks_type,str_concat((/"figs/CorrMap",sprinti("%0.2i", b)/)))
	 res = default_res   
	 res@gsnPolar             = "SH"
	 res@mpMaxLatF            = -50 
	 res@tiMainString         = "Basin correlations"
	 res@cnFillOn             = True
	 res@cnLinesOn  	 = False       ; do not draw contour lines
	 res@cnFillMode	    = "RasterFill"     ; Turn on raster color 
	 res@cnFillPalette     = "MPL_bwr"     ; set color map  
	 symMinMaxPlt (CorrMask,100,False,res)
	 Corrplot = gsn_csm_contour_map_polar(wks,CorrMask,res) 
	 delete(res)  

	 ;plot significance
	 res = default_res
	 res@cnMinLevelValF      = 0.00	     ; set min contour level
	 res@cnMaxLevelValF      = 1.00	     ; set max contour level
	 res@cnLevelSpacingF     = 0.01	     ; set contour spacing
	 res@cnLinesOn	       = False       ; do not draw contour lines
	 res@cnLineLabelsOn      = False     ; do not draw contour labels	 
	 res@cnFillScaleF        = 0.4	     ; add extra density 
	 significance_plot  = gsn_csm_contour(wks,Rtest,res)
         delete(res)	 
	 res		          = True	      
	 res@gsnShadeFillType     = "pattern"
	 res@gsnShadeLow	       = 6 
	 significance_plot   = gsn_contour_shade(significance_plot, 0.01, 0.01, res)  ; shade all areas >95%	 
         delete(res)	 
	 
	 overlay (Corrplot,significance_plot)
	 draw(Corrplot)
	 frame(wks)
	 
      end if

   end do
   
end
