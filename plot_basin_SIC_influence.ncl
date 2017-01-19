begin

   SD="/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/"
   fAISMask = addfile(str_concat((/SD,"data/AIS_Full_basins_Zwally_CESMgrid.nc"/)),"r")
   AISMask = fAISMask->Zwallybasins
   nAISMask = max(AISMask)

   do b=1,nAISMask
      ;load composites for low and high (10 years of monthly data each)

      flow = addfile(str_concat((/SD,"compositing_work/output/low/Basin",tostring(b),"_ssticetemp_fv0.9x1.25.nc"/)),"r")
      fhigh = addfile(str_concat((/SD,"compositing_work/output/high/Basin",tostring(b),"_ssticetemp_fv0.9x1.25.nc"/)),"r")
      low = month_to_annual(flow->ICEFRAC,1)
      high = month_to_annual(fhigh->ICEFRAC,1)

      ;for each point, calculate significance of SIC difference
      dimXY = dimsizes(low)
      printVarSummary(low)
      ntim  = dimXY(0)
      nlat  = dimXY(1)
      mlon  = dimXY(2)

      aveX = dim_avg_n (low,0)		 ; calculate means at each grid point 
      aveY = dim_avg_n (high,0)
      varX = dim_variance_n (low,0)	 ; calculate variances
      varY = dim_variance_n (high,0)
      sX = ntim
      sY = ntim
      alpha = ttest(aveX,varX,sX, aveY,varY,sY, False, False)
      ;plot SIC difference

      ;plot significance
      wks = gsn_open_wks ("png",str_concat((/"figs/",tostring(b)/)))
      res = True                                 ; plot mods desired
      res@gsnDraw              = False           ; Do not draw plot
      res@gsnFrame             = False           ; Do not advance frame  
      res@tiMainString         = "SIC Difference (high - low)"
      res@gsnCenterString      = "5% stippled"
      res@cnFillOn             = True               ; turn on color fill
      res@lbLabelAutoStride    = True
      res@lbOrientation = "vertical"
      diff= aveY-aveX
      
      
      plot = gsn_csm_contour(wks,diff, res) 

      res2 = True                            ; res2 probability plots

      res2@gsnDraw             = False       ; Do not draw plot
      res2@gsnFrame            = False       ; Do not advance frome

      res2@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
      res2@cnMinLevelValF      = 0.00        ; set min contour level
      res2@cnMaxLevelValF      = 1.05        ; set max contour level
      res2@cnLevelSpacingF     = 0.05        ; set contour spacing

      res2@cnInfoLabelOn       = False       ; turn off info label

      res2@cnLinesOn           = False       ; do not draw contour lines
      res2@cnLineLabelsOn      = False       ; do not draw contour labels

      res2@cnFillScaleF        = 0.6         ; add extra density   

      plot2   = gsn_csm_contour(wks,gsn_add_cyclic_point(alpha), res2)
      opt     = True
      opt@gsnShadeFillType = "pattern"
      opt@gsnShadeLow = 17 
      plot2   = gsn_contour_shade(plot2, 0.05, 30, opt)  ; shade all areas less than the
                                        	 ; 0.05 contour level 

      BasinMask=new(dimsizes(AISMask),"integer")
      AISBinaryMask=new(dimsizes(AISMask),"integer")

      BasinMask=where(AISMask.eq.b,1,0)
      AISBinaryMask=where(AISMask.eq.0,0,1)
      BasinMask=mask(BasinMask,AISBinaryMask,1)

      res3 = True    
      res3@gsnDraw              = False           ; Do not draw plot
      res3@gsnFrame             = False           ; Do not advance frame ;
      res3@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
      res3@cnMinLevelValF      = 0.00        ; set min contour level
      res3@cnMaxLevelValF      = 1.05        ; set max contour level
      res3@cnLevelSpacingF     = 0.05        ; set contour spacing   
      res3@cnFillOn             = True               ; turn on color fill   
      res3@cnLinesOn           = False       ; do not draw contour lines
      res3@cnLineLabelsOn      = False       ; do not draw contour labels
      res3@lbLabelBarOn = False
      plot3   = gsn_csm_contour(wks,BasinMask,res3)

      overlay (plot, plot2)
      overlay (plot,plot3)
      draw (plot)
      frame(wks)
   end do    
end
