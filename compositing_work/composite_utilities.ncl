procedure GenerateComposite(CompositeLevel,DataName,BasinNumber,MonthlyData,indt,IsOceanGrid)
   local indt2,i,j
   begin
  
      indt2 = new(dimsizes(indt)*12,"integer")
      do i = 0, dimsizes(indt)-1,1	
	  do j = 0,11
	      indt2(i*12+j) = indt(i)*12+j
	  end do
      end do
      rank=dimsizes(dimsizes(MonthlyData))
      
      if (rank.eq.3) then
	CompMonData=MonthlyData(indt2,:,:)
	CompMonMean=clmMonTLL(CompMonData)
      else if (rank.eq.4) then
	CompMonData=MonthlyData(indt2,:,:,:)
	CompMonMean=clmMonTLLL(CompMonData)      
      else
	print("Error: input monthly data to compositing procedure not rank [time,lat,lon] or [time,lev,lat,lon]")
      end if
      end if
      
      CompAnn=month_to_annual(CompMonData,1)
      CompAnnMean=dim_avg_n(CompAnn,0)
      CompAnnVar =dim_variance_n(CompAnn,0)  

      if (IsOceanGrid) then
      	CompMonMean:=PopLatLon(CompMonMean, "gx1v6","fv0.9x1.25","bilin","da","090309")
	CompAnn:=PopLatLon(CompAnn, "gx1v6","fv0.9x1.25","bilin","da","090309")
	CompAnnMean:=PopLatLon(CompAnnMean, "gx1v6","fv0.9x1.25","bilin","da","090309")
	CompAnnVar:=PopLatLon(CompAnnVar, "gx1v6","fv0.9x1.25","bilin","da","090309")
      end if      

      fname=str_concat((/"output/",CompositeLevel,"/",sprinti("%0.2i", BasinNumber),"/",DataName,".nc"/))
      print("Writing "+fname)
      system("rm -f "+fname)
      ncdf=addfile(fname ,"c")
      fAtt=True
      fAtt@creation_date = systemfunc ("date")
      fileattdef( ncdf, fAtt ) ; copy file attributes
      filedimdef(ncdf,"time",-1,True) 
      ncdf->CompMonData=CompMonData
      ncdf->CompMonMean=CompMonMean
      ncdf->CompAnnMean=CompAnnMean
      ncdf->CompAnnVar=CompAnnVar
      
   end
    
function LoadMonthlyData(Files,Var)
   begin
      in = addfiles(Files,"r")
      return in[:]->$Var$
   end

procedure save_timeseries(TS,fname)  
   local ncdf,fAtt
   begin
      system(str_concat((/"/bin/rm -f output/",fname/)))   ; remove any pre-existing file
      ncdf = addfile(str_concat((/"output/",fname/)),"c")  ; open output netCDF file
      fAtt               = True            ; assign file attributes
      fAtt@title         = "Time series"
      fAtt@source_file   =  fname
      fAtt@Conventions   = "None"
      fAtt@creation_date = systemfunc ("date")
      fileattdef( ncdf, fAtt )            ; copy file attributes
      filedimdef(ncdf,"time",-1,True)
      ncdf->time_series=TS
   end    
