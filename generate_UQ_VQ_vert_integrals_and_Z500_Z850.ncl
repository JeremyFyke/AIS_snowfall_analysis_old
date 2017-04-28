load "local_utilies.ncl"

calculate_vertically_integrated_moisture_transport=1
calculate_z500_and_z850=0

print("Loading lev,hyam,hybm,p0mb")
in=addfile("/glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/VQ/b.e11.BRCP85C5CNBDRD.f09_g16.105.cam.h0.VQ.200601-210012.nc","r")
lev=in->lev ;hPa, used lines 760-770 of amwg_diagostics/code/taylor_utils.ncl for reference to determine hPA units...
lev=lev*100. ; to Pa
hyam = in->hyam
hybm = in->hybm
P0mb=(in->P0)/100. ; convert to hPa/mb
delete(in)

print("Loading PS")
all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/PS/*.B1850C5CN.*.nc")
in = addfiles(all_files,"r")
ListSetType (in, "cat")
PSFC=in[:]->PS ;Pa
delete(in)

;wgt_vert_avg_beta
;p=pressure levels
;datai=input data to vertically integrate
;psfc=surface pressures
;punits: units of p and psfc.  0: hPa/mb.  1: Pa.
;sum or average: 0=average, 1=sum

g=9.80665

if (calculate_vertically_integrated_moisture_transport.eq.1) then
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
   print("Calculating vertically integrated moisture transport")
   UQ_vert_int=wgt_vert_avg_beta(lev,UQ,PSFC,1,0)

   print("Writing UQ_vert_int to file")
   system("/bin/rm -f jer_data_dir/Vertical_Int_Moisture_Transport_UQ.nc")
   setfileoption("nc","Format","LargeFile")
   ncdf = addfile("jer_data_dir/Vertical_Int_Moisture_Transport_UQ.nc","c")
   fAtt		   = True	     ; assign file attributes
   fAtt@title	   = ""
   fAtt@source_file   = "Vertical_Int_Moisture_Transport_UQ.nc"
   fAtt@Conventions   = "None"
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )	    ; copy file attributes
   filedimdef(ncdf,"time",-1,True)    
   ncdf->UQ_vert_integral = UQ_vert_int

   print("Writing VQ_vert_int to file")  
   system("/bin/rm -f jer_data_dir/Vertical_Int_Moisture_Transport_VQ.nc")
   setfileoption("nc","Format","LargeFile")
   ncdf = addfile("jer_data_dir/Vertical_Int_Moisture_Transport_VQ.nc","c")
   fAtt		   = True	     ; assign file attributes
   fAtt@title	   = ""
   fAtt@source_file   = "Vertical_Int_Moisture_Transport_VQ.nc"
   fAtt@Conventions   = "None"
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )	    ; copy file attributes
   filedimdef(ncdf,"time",-1,True)   
   ncdf->VQ_vert_integral = VQ_vert_int
end if

if (calculate_z500_and_z850.eq.1) then
   print("Loading Z3")
   all_files = systemfunc("ls /glade/p/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/monthly/Z3/*.B1850C5CN.*.nc")
   in = addfiles(all_files,"r")
   ListSetType (in, "cat")
   Z=in[:]->Z3
   delete(in)

   print("Calculating z500 and z850")
   Znew = vinth2p (Z,hyam,hybm,(/ 850.0,500.0 /),PSFC,1,P0mb,1,True)
   Z850=Znew(:,0,:,:)
   Z500=Znew(:,1,:,:)

   print("Writing z500 to file")   
   system("/bin/rm -f jer_data_dir/z500.nc")
   setfileoption("nc","Format","LargeFile")
   ncdf = addfile("jer_data_dir/z500.nc","c")
   fAtt               = True            ; assign file attributes
   fAtt@title         = ""
   fAtt@source_file   = "z500.nc"
   fAtt@Conventions   = "None"
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )            ; copy file attributes
   filedimdef(ncdf,"time",-1,True)   
   ncdf->z500 = Z500

   print("Writing z850 to file") 
   system("/bin/rm -f jer_data_dir/z850.nc")
   setfileoption("nc","Format","LargeFile")
   ncdf = addfile("jer_data_dir/z850.nc","c")
   fAtt               = True            ; assign file attributes
   fAtt@title         = ""
   fAtt@source_file   = "z850.nc"
   fAtt@Conventions   = "None"
   fAtt@creation_date = systemfunc ("date")
   fileattdef( ncdf, fAtt )            ; copy file attributes
   filedimdef(ncdf,"time",-1,True)   
   ncdf->z850 = Z850
end if
