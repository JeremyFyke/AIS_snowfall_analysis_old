load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


begin

; core names

names = asciiread("names.txt",-1,"string")
print(names)

; geo data

tmp = asciiread("lat_lon_elev.txt",(/80,3/),"double")
lat = tmp(:,0)
lon = tmp(:,1)
elev = tmp(:,2)
delete(tmp)

; actual data 
; first column is year, other 80 columns are individual cores
tmp = asciiread("data.txt",(/796,81/),"double")
year = tmp(:,0)
coresmb = tmp(:,1::)
coresmb@_FillValue = -9999.
coresmb = where(coresmb.lt.0.,coresmb@_FillValue,coresmb)
delete(tmp)
;**********************************************************

; calculate cross correlation
smbcorr = new((/80,80/),"double")

do i = 0,79
	do j = 0,79
	smbcorr(i,j) = escorc(coresmb(:,i),coresmb(:,j))
	end do
end do

print(coresmb(:,0))
print(smbcorr(0,:))

;*************************************************************

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

print(names(i))
printMinMax(smbcorr(i,:),False)
printMinMax(coresmb(:,i),False)

do j = 0, 79,1
resp@gsMarkerIndex = 4          ; overdraw with empty dot
resp@gsMarkerColor      = "black"

putmark2(j) = gsn_add_polymarker(wks,map,lon(j),lat(j),resp)


; assign colors to markers
; -0.6 should be rgb color 2 and 0.6 rgb color 20 (see https://www.ncl.ucar.edu/Document/Graphics/ColorTables/temp_diff_18lev.shtml)
; looks like something like this:
;print(fspan(-0.6,0.6,20)+2.6+fspan(0,1.2,20)*14)
;exit
tmp = tointeger(smbcorr(i,j)+2.6+(smbcorr(i,j)+0.6)*14)
if (ismissing(tmp)) then
	resp@gsMarkerColor      = "white"
else
	resp@gsMarkerColor      = tmp
end if
resp@gsMarkerIndex = 16          ; Use filled dots for markers.
putmark(j) = gsn_add_polymarker(wks,map,lon(j),lat(j),resp)
end do


; now plot location of core on top
pres = True
pres@gsMarkerSizeF = 0.02
pres@gsMarkerIndex = 16
actualcore = gsn_add_polymarker(wks,map,lon(i),lat(i),pres)


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

end