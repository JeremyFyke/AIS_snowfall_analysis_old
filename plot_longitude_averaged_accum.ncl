load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "local_utilies.ncl"
in = addfile("compositing_work/output/LongitudeAccumulationTimeSeries_monthly.nc","r")
Accum=in->time_series    
ds=dimsizes(Accum)
do t=0,ds(0)-1
   Accum(t,:)=rmAnnCycle1D(Accum(t,:))   
end do
do l=0,ds(1)-1
   Accum(:,l)=dim_standardize(Accum(:,l),1)
end do
opt = 0
wgt  = (/0.0270, 0.05856, 0.09030, 0.11742, 0.13567, \ 
	 0.1421, 0.13567, 0.11742, 0.09030, 0.05856, 0.027/)
do t=0,ds(0)-1
   Accum(t,:)=wgt_runave(Accum(t,:), wgt, opt)   
end do
do l=0,ds(1)-1
   Accum(:,l)=wgt_runave(Accum(:,l),1)
end do


print(max(Accum))
print(min(Accum))
wks = gsn_open_wks ("png", "Accum")
res=get_default_res()

res@cnMinLevelValF       = -4              ; set min contour level
res@cnMaxLevelValF       =  4              ; set max contour level
res@cnLevelSpacingF      =  0.2
res@cnFillOn = True            
res@cnFillPalette = "MPL_bwr"
res@cnLinesOn  	 = False
plot=gsn_csm_contour(wks, Accum(1001:1120,:), res )
draw(plot)
frame(wks)
