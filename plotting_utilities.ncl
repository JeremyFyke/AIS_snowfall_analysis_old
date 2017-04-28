;Plotting functions (useful for overlaying information)
   function plot_basin_correlation(CorrMask,wks,isPolar,mpMaxLatF)
      begin
      	 res = get_default_res()  
	 res@mpMaxLatF            = mpMaxLatF
	 res@gsnPolar             = "SH"
	 res@mpMaxLatF            = mpMaxLatF
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
   
   function plot_AISMask_overview(AISMask,wks,isPolar,mpMaxLatF)
      begin
      	 res = get_default_res()  
	 res@cnFillOn             = True
	 res@cnLinesOn  	 = False       ; do not draw contour lines
	 res@cnFillMode	    = "RasterFill"     ; Turn on raster color 
	 res@cnFillPalette     = "MPL_bwr"     ; set color map
	 if (isPolar.eq.1) then
	    res@mpMaxLatF            = mpMaxLatF
	    res@gsnPolar             = "SH"
	    res@mpMaxLatF            = mpMaxLatF	 
	    plot = gsn_csm_contour_map_polar(wks,AISMask,res)
	 else
	    plot = gsn_csm_contour_map(wks,AISMask,res)	 
	 end if
         return plot
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
	 opt@gsnShadeFillMode     = "RasterFill"  ; Turn on raster fill     
	 opt@gsnShadeFillType     = "pattern"
	 opt@gsnShadeFillScales   = 1.0
	 opt@gsnShadeLow	  = 17 
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
	
      function plot_basin_polygons(ZwallyBasinPolygonData,plot,wks,nAISMask)
	 begin
	    res = get_default_res()
	    ZwallyBasinPolygonLats=ZwallyBasinPolygonData(:,0)
	    ZwallyBasinPolygonLons=ZwallyBasinPolygonData(:,1)
	    ZwallyBasinPolygonInds=toint(ZwallyBasinPolygonData(:,2))
	    res@gsnPolar             = "SH"
	    res@mpMaxLatF            = -60
	    res@mpFillOn = False
	    ;basin_polygon_plot=gsn_csm_map_polar(wks,res)
	    res_poly             = True
	    res_poly@gsEdgesOn   = True       ; draw border around polygons
	    res_poly@gsEdgeColor = "black"
	    res_poly@gsFillColor = "white"
	    res_poly@gsFillOpacityF = 0.0
	    res_poly@gsEdgeThicknessF = 10
	    a=new(toint(nAISMask),graphic)
	    do b=1,nAISMask
	       PolyInds=ind(ZwallyBasinPolygonInds.eq.b)
	       a(b-1)=gsn_add_polygon(wks,plot,\
	        	   ZwallyBasinPolygonLons(PolyInds),\
	        	   ZwallyBasinPolygonLats(PolyInds),\
			   res_poly)
	       delete(PolyInds)
            end do
	    plot@a=a ;even though this isn't explicitly used, it needs to be retained.  Weird.         
	    return plot
	 end
 
      function plot_core_locations(plot,wks,ice_lat,ice_lon)
	 begin
	    ds=dimsizes(ice_lat)
	    putmark = new(ds,"graphic")
	    res=True
	    res@gsMarkerIndex = 4  ; overdraw with empty dot
	    res@gsMarkerColor = "black"
	    res@gsMarkerIndex = 16
	    do c = 1,ds(0)
	       cm1=c-1
	       putmark(cm1) = gsn_add_polymarker(wks,plot,ice_lon(cm1),ice_lat(cm1),res)
	    end do
	    plot@putmark=putmark
	    return plot
	 end 
 
      function plot_vector_scalar_field(u,v,s,wks,isPolar,mpMaxLatF)
	 begin
	    u(0,:)=0.0 ; remove long pole vector
	    v(0,:)=0.0 ; remove long pole vector	          
	    res = get_default_res()
	    res@gsnPolar             = "SH"
	    res@mpMaxLatF            = mpMaxLatF
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

      function plot_precip_climatology(PRECIP_climatology,cnLevels,wks,isPolar)
	 begin
	    res=get_default_res()
	    res@gsnPolar            = "SH"
	    res@mpMaxLatF           = -60 ; specify min lat
	    res@cnFillOn            = True
	    res@cnLinesOn           = False
	    res@cnFillPalette       =  "GMT_drywet"
	    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
	    res@cnLevels=cnLevels
	    if (isPolar.eq.1) then
	      accum_plot = gsn_csm_contour_map_polar(wks,PRECIP_climatology,res)
	    else
	      accum_plot = gsn_csm_contour(wks,PRECIP_climatology,res)
	    end if
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
	    ;res@cnMinLevelValF       = -1000               ; set min contour level
	    ;res@cnMaxLevelValF       = 1000               ; set max contour level
	    ;res@cnLevelSpacingF = 50. 	    
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

       function plot_basin_correlation_matrix(Xcor_basin_by_basin,wks) 
          begin
	    res = get_default_res()
	    res@gsnMaximize          = True          ; Maximize plot in frame.
	    res@cnFillOn             = True          ; Turn on contour fill
	    res@cnFillMode           = "RasterFill"  ; Turn on raster fill
	    res@cnLinesOn            = False         ; Turn off contour
	    res@cnFillPalette     = "BlWhRe"     ; set color map  
	    res@cnMinLevelValF  = -.9			; set the minimum contour level
	    res@cnMaxLevelValF  = .9			; set the maximum contour level
	    res@cnLevelSpacingF = .1			; set the interval between contours       
	    ;res@sfXArray = fspan(0,27,28)    ; One more element in both
	    ;res@sfYArray = fspan(0,27,28)    ; X and Y
	    ;res@tmXMajorGrid = True
	    ;res@tmYMajorGrid = True    
	    ;res@tmXMinorGrid = True
	    ;res@tmYMinorGrid = True
	    ;res@trXMinF     =  -0.5
	    ;res@trYMinF     =  -0.5        
	    res@trXMaxF     =  25.5
	    res@trYMaxF     =  25.5                  
	    corr_mat_plot = gsn_csm_contour(wks,Xcor_basin_by_basin,res)	    
	    return corr_mat_plot
	  end  
