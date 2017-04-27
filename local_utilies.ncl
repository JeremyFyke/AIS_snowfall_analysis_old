;Local code snippets that are re-used among analysis routines
function get_default_wks()
  begin
  
     wks_type = "png"
     wks_type@wkWidth = 2500
     wks_type@wkHeight = 2500  
     return wks_type
     
  end
;;;

function get_default_res()
  begin
     
     default_res = True
     default_res@gsnDraw               = False           ; Do not draw plot
     default_res@gsnFrame              = False           ; Do not advance frame
     default_res@cnLevelSelectionMode  = "ManualLevels"     ; set manual contour levels
     default_res@mpOutlineOn          = False
     default_res@gsnAddCyclic          = True       
     default_res@cnInfoLabelOn         = False       ; turn off info label
     tt=5
     default_res@tiMainFontThicknessF = tt
     default_res@tiXAxisFontThicknessF = tt
     default_res@tiYAxisFontThicknessF = tt
     default_res@lbLabelFontThicknessF = tt  
     default_res@tmXBLabelFontThicknessF = tt
     default_res@tmXTLabelFontThicknessF = tt  
     default_res@tmYLLabelFontThicknessF = tt      
     default_res@tmYRLabelFontThicknessF = tt
     default_res@lbTitleOn = False
     default_res@lbTitlePosition = "Bottom"
     default_res@lbTitleFontThicknessF = tt
     return default_res
     
  end
;;;     
function get_AIS_mask()
   begin
      fAISMask = addfile("/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/AIS_basin_masks/AIS_Full_basins_Zwally_CESMgrid.nc","r")
      AISMask = fAISMask->Zwallybasins
      return AISMask
   end
;;;      
function get_LandFrac()
   begin
      in = addfile("/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/LANDFRAC/b.e11.B1850C5CN.f09_g16.005.cam.h0.LANDFRAC.040001-049912.nc","r")
      LandFrac = dim_avg_n(rm_single_dims(in->LANDFRAC),0)
      return LandFrac
   end
;;;
function get_LandArea()
   begin
      in = addfile("/glade/p/cesm/cseg/inputdata/lnd/clm2/griddata/griddata_0.9x1.25_070212.nc","r")
      LandArea = in->AREA   *1.e6;km^2 -> m^2 
      return LandArea
   end
;;;
function add_lat_lon_arrays(field,lat,lon)
begin
   field!0="lat"
   field!1="lon"
   field&lat=lat
   field&lon=lon
   field&lat@units="degrees_north"
   field&lon@units="degrees_east" 
   return field
end
;;;
function replace_mask_values(field,ReplaceValue)
local field1d
begin
   field1d=ndtooned(field)
   field1d (ind(ismissing( field1d ))) = ReplaceValue;
   field=onedtond(field1d,dimsizes(field))
   delete_VarAtts(field,(/"missing_value","_FillValue","coordinates"/))
   return field
end
;;;
function calculate_integrated_timeseries(Field,LandFrac,LandArea,Mask)
local i,nt,dimsize
begin
   ; calculate integrated time series
   dimsize=dimsizes(Field)
   nt=dimsize(0)
   TS=new(nt,"double")
   do i = 0, nt-1,1 
      TS(i) = sum(Field(i,:,:)*LandFrac*LandArea*Mask)
   end do
   return TS
end

procedure plot_pdf(Vector,FileName,Xlabel,Ylabel,MainString)
local wks,res,plot,probability
begin
   wks = gsn_open_wks ("png",FileName)
   plot = new(2,graphic) ; create a plot array
   probability=pdfx( Vector , 50 , False)
   res=True
   res@xyLineThicknessF       = 2
   res@tiYAxisString          = "PDF (%)"
   res@tiXAxisString          = "Integrated precipitation (Gt/yr)"
   res@tiMainString           = MainString 
   res@gsnDraw  = False                          ; don't draw
   res@gsnFrame = False                          ; don't advance frame   
   plot(0)=gsn_csm_xy (wks,probability@bin_center,probability,res)     
   delete(res)
   res=True
   res@gsnDraw  = False                          ; don't draw
   res@gsnFrame = False                          ; don't advance frame   
   plot(1) = gsn_csm_y(wks,Vector,res)
   gsn_panel(wks,plot,(/2,1/),False)
end

procedure plot_AIS_and_basin_accum_PDFs(AIS_integrated_TS,AIS_integrated_TS_avg,\
                                        Basin_integrated_TS,Basin_integrated_TS_avg,\
					nAISMASK)
local wks,res,plot,Vector,bin_min,bin_max,opt,nbins,probability,plot0,b,bm1,plot1
begin					
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
end

function composite_PAGES_data(climate_lat_vec,climate_lon_vec,AISMask,nAISMask)
local names,tmp
begin
   ;;;Read ice core data;;;
   SD="/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/jer_data_dir/PAGES_data/"
   names = asciiread(str_concat((/SD,"names.txt"/)),-1,"string")
   ; geo data
   tmp = asciiread(str_concat((/SD,"lat_lon_elev.txt"/)),(/80,3/),"double")
   ice_lat = tmp(:,0)
   ice_lon = tmp(:,1)
   ice_lon=where(ice_lon<0.,ice_lon+360.,ice_lon)
   ice_elev = tmp(:,2)
   ; actual data 
   ; first column is year, other 80 columns are individual cores
   tmp := asciiread(str_concat((/SD,"data.txt"/)),(/796,81/),"double")
   ice_year = tmp(:,0)
   ice_coresmb = tmp(:,1::)
   ds=dimsizes(ice_coresmb)
   nrecs=ds(0)   
   ncores=ds(1)
   ice_coresmb@_FillValue = -9999.
   ice_coresmb = where(ice_coresmb.lt.0.,ice_coresmb@_FillValue,ice_coresmb)
   ;ice_coresmb=dtrend_quadratic_msg_n(ice_coresmb,False,False,0)
   ;ice_coresmb=dim_standardize_n(ice_coresmb,0,1)

   ice_cv=dim_stddev_n_Wrap(ice_coresmb,0)/dim_avg_n_Wrap(ice_coresmb,0)
   
   climate_lat=conform(AISMask,climate_lat_vec,0)
   climate_lon=conform(AISMask,climate_lon_vec,1)

   print("###Cliamte lat/lon range###")   
   print(min(climate_lat))
   print(max(climate_lat))
   print(min(climate_lon))
   print(max(climate_lon))
   
   print("###Ice core lat/lon range###")
   print(min(ice_lat))
   print(max(ice_lat))
   print(min(ice_lon))
   print(max(ice_lon))   

   count=new((/nrecs,toint(nAISMask)/),"integer")
   count(:,:)=0
   basin_composite_core=new((/nrecs,toint(nAISMask)/),typeof(ice_coresmb))
   basin_composite_core(:,:)=0.0
   basin_composite_core@_FillValue = -9999.
   AISMask!0="lat"
   AISMask!1="lon"
   AISMask&lat=climate_lat_vec
   AISMask&lon=climate_lon_vec

   ;;;;Stack ice core data in basin-specific composites;;;
   do c=1,ncores
      cm1=c-1   
      core_grid_loc_indices=getind_latlon2d(climate_lat,climate_lon,ice_lat(cm1),ice_lon(cm1))
      i=core_grid_loc_indices(0,0)
      j=core_grid_loc_indices(0,1)
      lat=ice_lat(cm1)
      lon=ice_lon(cm1)
      print("Core name="+names(cm1))
      print("lat="+lat)
      print("lon="+lon)
      b=round(AISMask(i,j),3)     
      ;b=round(AISMask({lat},{lon}),3)
      print("basin="+b)

      if (b.gt.0) then ;If core falls in a defined basin (and not the ocean), stack core data
         print("basin="+b)
         bm1=b-1
         do t = 1,nrecs
	    tm1=t-1
	    if (.not. ismissing(ice_coresmb(tm1,cm1))) then
               basin_composite_core(tm1,bm1)=basin_composite_core(tm1,bm1)+ice_coresmb(tm1,cm1) ;accumulate basin-falling core
	       count(tm1,bm1)=count(tm1,bm1)+1 ; bump basin counter
	    end if
	  end do
      end if
   end do
   do b = 1,nAISMask
      do t=1,nrecs
	 bm1=b-1
	 tm1=t-1
	 if (count(tm1,bm1).gt.0) then
                basin_composite_core(tm1,bm1)=basin_composite_core(tm1,bm1)/count(tm1,bm1)
	 else
	        basin_composite_core(tm1,bm1)=basin_composite_core@_FillValue
	 end if
      end do
   end do
   
   plot_core_correlations=0
   if (plot_core_correlations.eq.1) then
   
      ;;; Jan code follows
      smbcorr = new((/80,80/),"double")
      do i = 0,79
	      do j = 0,79
	      smbcorr(i,j) = escorc(ice_coresmb(:,i),ice_coresmb(:,j))
	      end do
      end do   
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      wks = gsn_open_wks("pdf","Ant_cores_correlation")


      res   = True
      res@mpGridAndLimbOn        = False
      res@mpFillOn               = False
      res@gsnPolar               = "SH"
      res@mpMinLatF  = -90.
      res@mpMaxLatF  = -62.
      res@cnFillOn          = True              ; color plot desired
      res@cnLineLabelsOn    = False             ; turn off contour lines
      res@txFontHeightF     = 0.015 
      res@mpDataBaseVersion = "MediumRes"            ; choose higher resolution
      ;res@mpOutlineDrawOrder = "PreDraw"

      res@gsnDraw = False
      res@gsnFrame = False

      gsn_define_colormap(wks,"temp_diff_18lev")      ; -1: blue // 0: white // 1: red

      do i = 0,79,1

      map = gsn_csm_map(wks,res)


      putmark = new(80,"graphic")
      putmark2= new(80,"graphic")

      resp = True 
      resp@gsMarkerIndex = 16          ; Use filled dots for markers.

      do j = 0, 79,1
	  resp@gsMarkerIndex = 4          ; overdraw with empty dot
	  resp@gsMarkerColor      = "black"

	  putmark2(j) = gsn_add_polymarker(wks,map,ice_lon(j),ice_lat(j),resp)
	  tmp := tointeger(smbcorr(i,j)+2.6+(smbcorr(i,j)+0.6)*14)
	  if (ismissing(tmp)) then
		  resp@gsMarkerColor      = "white"
	  else
		  resp@gsMarkerColor      = tmp
	  end if
	  resp@gsMarkerIndex = 16          ; Use filled dots for markers.
	  putmark(j) = gsn_add_polymarker(wks,map,ice_lon(j),ice_lat(j),resp)
	  end do


	  ; now plot location of core on top
	  pres = True
	  pres@gsMarkerSizeF = 0.02
	  pres@gsMarkerIndex = 16
	  actualcore = gsn_add_polymarker(wks,map,ice_lon(i),ice_lat(i),pres)


	  ; plot name on top
	  rest = True
	  rest@txFontHeightF = 0.03
	  gsn_text_ndc(wks,names(i),0.5,0.1,rest)

	  draw(map)
	  frame(wks)

	  delete(map)
	  delete(putmark)
	  delete(putmark2)
	  delete(actualcore)
      end do
   end if

   return([/ice_lat,ice_lon,ice_elev,ice_year,ice_coresmb,basin_composite_core,ice_cv/])
   end

function calculate_interbasin_correlations(BasinAccum,nAISMask)
   begin
   Xcor_basin_by_basin=new((/toint(nAISMask),toint(nAISMask)/),"double")
   Rtest_basin_by_basin=new((/toint(nAISMask),toint(nAISMask)/),"double")
   nPosCorrs_total=new(toint(nAISMask),"double")
   nNegCorrs_total=new(toint(nAISMask),"double")
   nPosCorrs_significant=new(toint(nAISMask),"double")
   nNegCorrs_significant=new(toint(nAISMask),"double")   
   nPosCorrs_total(:)=0.0
   nNegCorrs_total(:)=0.0   
   nPosCorrs_significant(:)=0.0
   nNegCorrs_significant(:)=0.0
   siglvl=0.05
   ds=dimsizes(BasinAccum)
   nt=ds(0)
   nb=ds(1)
   do b=1,nAISMask
      bm1=b-1      
      do bc=1,nAISMask ;loop over all basins, color in based on correlation at that basin
	 bcm1=bc-1
	 Xcor_basin_by_basin(bm1,bcm1)=escorc(BasinAccum(:,bm1),BasinAccum(:,bcm1))
	 Rtest_basin_by_basin(bm1,bcm1)=rtest(Xcor_basin_by_basin(bm1,bcm1),nt,0)
	 if (.not. ismissing(Xcor_basin_by_basin(bm1,bcm1))) then
	    if (Xcor_basin_by_basin(bm1,bcm1).gt.0) then
	       nPosCorrs_total(bm1)=nPosCorrs_total(bm1)+1
	       if (Rtest_basin_by_basin(bm1,bcm1).lt.siglvl) then
		  nPosCorrs_significant(bm1)=nPosCorrs_significant(bm1)+1
	       end if
	    end if
	    if (Xcor_basin_by_basin(bm1,bcm1).lt.0) then
	       nNegCorrs_total(bm1)=nNegCorrs_total(bm1)+1
	       if (Rtest_basin_by_basin(bm1,bcm1).lt.siglvl) then
		  nNegCorrs_significant(bm1)=nNegCorrs_significant(bm1)+1
	       end if
	    end if
	 end if 
      end do
   end do   
   return([/nPosCorrs_total,nNegCorrs_total,nPosCorrs_significant,nNegCorrs_significant,Xcor_basin_by_basin,Rtest_basin_by_basin/])
   end

function make_correlation_mask(AISMask,nAISMask,Xcor,time,bm1)
   begin
   AISBinaryMask=where(AISMask.eq.0,0,1)   
   CorrMask=new(dimsizes(AISMask),"double")
   CorrMask(:,:)=0.0 
   do bc=1,nAISMask ;loop over all basins, color in based on correlation at that basin
      bcm1=bc-1
      CorrMask=where(AISMask.eq.bc,Xcor(bm1,bcm1),CorrMask)
   end do
   CorrMask=mask(CorrMask,AISBinaryMask,1)     
   RtestMask=rtest(CorrMask,dimsizes(time),0)
   RtestMask=mask(RtestMask,AISBinaryMask,1)
   return([/CorrMask,RtestMask/])
   end
   
