
print("Loading lev,hyam,hybm,p0mb,PSFC")
in=addfile("/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PS/b.e11.B1850C5CN.f09_g16.005.cam.h0.PS.040001-049912.nc","r")
lev=in->lev
hyam = in->hyam
hybm = in->hybm
P0mb=in->P0
P0mb=P0mb/100. ;convert to mb from Pa
PSFC=in->PS
delete(in)

print("Loading Z3")
in = addfile("/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/Z3/b.e11.B1850C5CN.f09_g16.005.cam.h0.Z3.040001-049912.nc","r")
Z=in->Z3
delete(in)

print("Calculating z500 and z850")
;vinth2p documentation:
;http://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p.shtml
Znew = vinth2p (Z,hyam,hybm,(/ 850.0,500.0 /),PSFC,1,P0mb,1,True)
Z850=Znew(:,0,:,:)
Z500=Znew(:,1,:,:)

;make some quick plots of first timeslice of data
print("Making some test plots")
wks = gsn_open_wks ("png", "z500")
res=True
plotv  = gsn_csm_contour(wks, Z500(0,:,:), res )

wks = gsn_open_wks ("png", "z850")
res=True
plotv  = gsn_csm_contour(wks, Z850(0,:,:), res )

