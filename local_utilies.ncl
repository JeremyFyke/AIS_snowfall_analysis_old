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
     return default_res
     
  end
;;;     
function get_AIS_mask()
   begin
      SD="/glade/u/home/jfyke/liwg/AIS_snowfall_analysis/fyke_analysis/"
      fAISMask = addfile(str_concat((/SD,"AIS_basin_masks/AIS_Full_basins_Zwally_CESMgrid.nc"/)),"r")
      AISMask = fAISMask->Zwallybasins
      return AISMask
   end
;;;      
function get_LandFrac()
   begin
      in = addfile("/glade/scratch/hailong/archive/composite_ICE_wtag_mean/atm/hist/composite_ICE_wtag_mean.cam.h0.0001-01.nc","r")
      LandFrac = rm_single_dims(in->LANDFRAC)
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
