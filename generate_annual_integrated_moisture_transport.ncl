load "local_utilies.ncl"

print("Loading lev")
in=addfile("/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/VQ/b.e11.BRCP85C5CNBDRD.f09_g16.105.cam.h0.VQ.200601-210012.nc","r")
lev=in->lev
delete(in)

print("Loading PS")
all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PS/*.B1850C5CN.*.nc")
in = addfiles(all_files,"r")
ListSetType (in, "cat")
PSFC=in[:]->PS
delete(in)

print("Loading VQ")
all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/VQ/*.B1850C5CN.*nc")
in = addfiles(all_files,"r")
ListSetType (in, "cat")
VQ=in[:]->VQ
delete(in)
VQ_vert_int=wgt_vert_avg_beta(lev,VQ,PSFC,1,0)

print("Loading UQ")
all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/UQ/*.B1850C5CN.*.nc")
in = addfiles(all_files,"r")
ListSetType (in, "cat")
UQ=in[:]->UQ
delete(in)
UQ_vert_int=wgt_vert_avg_beta(lev,UQ,PSFC,1,0)

print("Writing files")
system("/bin/rm -f jer_data_dir/Vertical_Int_Moisture_Transport_UQ.nc")
setfileoption("nc","Format","LargeFile")
ncdf = addfile("jer_data_dir/Vertical_Int_Moisture_Transport_UQ.nc","c")


;Write out.

;===================================================================
; create global attributes of the file (optional)
;===================================================================
fAtt		   = True	     ; assign file attributes
fAtt@title	   = ""
fAtt@source_file   = "Vertical_Int_Moisture_Transport_UQ.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )	    ; copy file attributes
filedimdef(ncdf,"time",-1,True)    
ncdf->UQ_vert_integral = UQ_vert_int

system("/bin/rm -f jer_data_dir/Vertical_Int_Moisture_Transport_VQ.nc")
setfileoption("nc","Format","LargeFile")
ncdf = addfile("jer_data_dir/Vertical_Int_Moisture_Transport_VQ.nc","c")
;===================================================================
; create global attributes of the file (optional)
;===================================================================
fAtt               = True            ; assign file attributes
fAtt@title         = ""
fAtt@source_file   = "Vertical_Int_Moisture_Transport_VQ.nc"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes
filedimdef(ncdf,"time",-1,True)   
ncdf->VQ_vert_integral = VQ_vert_int
