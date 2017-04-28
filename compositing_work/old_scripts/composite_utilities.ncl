procedure generate_composites(CompositePeriod,indt,SSTm,ICEFRACm,PRECIPm,SLPm,UQm,VQm,UQ_VERT_INTEGRALm,VQ_VERT_INTEGRALm,EVAPm,TREFHTm,Z850m,Z500m,lev,lat,lon,OutPutFileName)
local indt2,i,j
begin
    indt2 = new(dimsizes(indt)*12,"integer")
    do i = 0, dimsizes(indt)-1,1	
	do j = 0,11
	    indt2(i*12+j) = indt(i)*12+j
	end do
    end do
    ;Calculate monthly mean climatologies of composites
    SIC_gx1v6        = ICEFRACm(indt2,:,:)
    SST_gx1v6        = SSTm(indt2,:,:)
    PRECIP           = PRECIPm(indt2,:,:)
    SLP              = SLPm(indt2,:,:)
    UQ               = UQm(indt2,:,:,:)
    VQ               = VQm(indt2,:,:,:)
    UQ_VERT_INTEGRAL = UQ_VERT_INTEGRALm(indt2,:,:)
    VQ_VERT_INTEGRAL = VQ_VERT_INTEGRALm(indt2,:,:)
    EVAP             = EVAPm(indt2,:,:)    
    ;EVAP_gx1v6       = EVAPm(indt2,:,:)
    TREFHT           = TREFHTm(indt2,:,:)
    Z850             = Z850m(indt2,:,:)
    Z500             = Z500m(indt2,:,:)
    
    SST = PopLatLon(SST_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    SST = replace_mask_values(SST,-1.8)
    SIC = PopLatLon(SIC_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    SIC = replace_mask_values(SIC,1.)
    ;EVAP = PopLatLon(EVAP_gx1v6, "gx1v6","fv0.9x1.25","bilin","da","090309")
    ;EVAP = replace_mask_values(EVAP,0.)

    ;Calculate monthly climatologies, annual averages/variances of composites for output
    SICm_mean=clmMonTLL(SIC)
    SIC_ann=month_to_annual(SIC,1)
    SIC_ann_mean=dim_avg_n (SIC_ann,0)
    SIC_ann_var =dim_variance_n (SIC_ann,0)
    
    SSTm_mean=clmMonTLL(SST)
    SST_ann=month_to_annual(SST,1)
    SST_ann_mean=dim_avg_n (SST_ann,0)
    SST_ann_var =dim_variance_n (SST_ann,0)

    PRECIPm_mean=clmMonTLL(PRECIP)
    PRECIP_ann=month_to_annual(PRECIP,1)
    PRECIP_ann_mean=dim_avg_n (PRECIP_ann,0)
    PRECIP_ann_var =dim_variance_n (PRECIP_ann,0)

    SLPm_mean=clmMonTLL(SLP)
    SLP_ann=month_to_annual(SLP,1)
    SLP_ann_mean=dim_avg_n (SLP_ann,0)
    SLP_ann_var =dim_variance_n (SLP_ann,0)
    
    UQm_mean=clmMonTLLL(UQ)    
    UQ_ann=month_to_annual(UQ,1)
    UQ_ann_mean=dim_avg_n (UQ_ann,0)
    UQ_ann_var =dim_variance_n (UQ_ann,0) 
       
    VQm_mean=clmMonTLLL(VQ)
    VQ_ann=month_to_annual(VQ,1)
    VQ_ann_mean=dim_avg_n (VQ_ann,0)
    VQ_ann_var =dim_variance_n (VQ_ann,0)    
    
    UQ_VERT_INTEGRALm_mean=clmMonTLL(UQ_VERT_INTEGRAL)    
    UQ_VERT_INTEGRAL_ann=month_to_annual(UQ_VERT_INTEGRAL,1)
    UQ_VERT_INTEGRAL_ann_mean=dim_avg_n (UQ_VERT_INTEGRAL_ann,0)
    UQ_VERT_INTEGRAL_ann_var =dim_variance_n (UQ_VERT_INTEGRAL_ann,0) 
       
    VQ_VERT_INTEGRALm_mean=clmMonTLL(VQ_VERT_INTEGRAL)
    VQ_VERT_INTEGRAL_ann=month_to_annual(VQ_VERT_INTEGRAL,1)
    VQ_VERT_INTEGRAL_ann_mean=dim_avg_n (VQ_VERT_INTEGRAL_ann,0)
    VQ_VERT_INTEGRAL_ann_var =dim_variance_n (VQ_VERT_INTEGRAL_ann,0)
    
    EVAPm_mean=clmMonTLL(EVAP)
    EVAP_ann=month_to_annual(EVAP,1)
    EVAP_ann_mean=dim_avg_n (EVAP_ann,0)
    EVAP_ann_var =dim_variance_n (EVAP_ann,0)

    TREFHTm_mean=clmMonTLL(TREFHT)
    TREFHT_ann=month_to_annual(TREFHT,1)
    TREFHT_ann_mean=dim_avg_n (TREFHT_ann,0)
    TREFHT_ann_var =dim_variance_n (TREFHT_ann,0)    

    Z500m_mean=clmMonTLL(Z500)
    Z500_ann=month_to_annual(Z500,1)
    Z500_ann_mean=dim_avg_n (Z500_ann,0)
    Z500_ann_var =dim_variance_n (Z500_ann,0)
    
    Z850m_mean=clmMonTLL(Z850)
    Z850_ann=month_to_annual(Z850,1)
    Z850_ann_mean=dim_avg_n (Z850_ann,0)
    Z850_ann_var =dim_variance_n (Z850_ann,0)      

    arrSize=dimsizes(SIC_ann)
    nSamples=arrSize(0)

    time    = new(12,"double")
    date    = new(12,"integer")
    datesec = new(12,"integer")
    it=0
    totaltime=0.
    year=0000
    year@calendar="365_day"
    do month = 01,12
       totaltime=totaltime+int2flt(days_in_month(year,month))
       t=totaltime-int2flt(days_in_month(year,month))/2.
       time(it)=t
       t@calendar="365_day"
       t@units="days since 0000-1-1"
       date(it)=cd_calendar(t,-2)-1
       datesec(it)=floattointeger((totaltime-t)*86400)
       it=it+1
    end do

    fname=str_concat((/"output/",CompositePeriod,"/",OutPutFileName/))
    print("Writing "+fname)
    system("rm -f "+fname)
    fout = addfile(fname,"c")

    dimNames = (/"time", "lev", "lat", "lon"/)
    dsizes_lev=dimsizes(lev)
    dsizes_lat=dimsizes(lat)
    dsizes_lon=dimsizes(lon)
    dimSizes = (/dimsizes(year), dimsizes(lev), dsizes_lat, dsizes_lon/)
    dimUnlim = (/ True, False, False, False/)   
    filedimdef(fout,dimNames,dimSizes,dimUnlim)

    filevardef(fout, "nSamples",typeof(nSamples)  , "ncl_scalar")
    filevardef(fout, "time"     ,typeof(time)     , (/"time"/))
    filevardef(fout, "date"     ,typeof(date)     , (/"time"/))
    filevardef(fout, "datesec"  ,typeof(datesec)  , (/"time"/))	
    filevardef(fout, "lev"      ,typeof(lev)      , (/"lat"/) )     
    filevardef(fout, "lat"      ,typeof(lat)      , (/"lat"/) )                          
    filevardef(fout, "lon"      ,typeof(lon)      , (/"lon"/) )	 

    filevardef(fout, "monthly_climo_ICEFRAC"          , typeof(SICm_mean)	        , (/"time", "lat", "lon" /) )
    filevardef(fout, "monthly_climo_SST"              , typeof(SSTm_mean)	        , (/"time", "lat", "lon" /) )
    filevardef(fout, "monthly_climo_PRECIP"           , typeof(PRECIPm_mean)	        , (/"time", "lat", "lon" /) )
    filevardef(fout, "monthly_climo_SLP"              , typeof(SLPm_mean)	        , (/"time", "lat", "lon" /) )
    filevardef(fout, "monthly_climo_UQ"               , typeof(UQm_mean)  ,               (/"time", "lev", "lat", "lon" /) )	 
    filevardef(fout, "monthly_climo_VQ"               , typeof(VQm_mean)  ,               (/"time", "lev", "lat", "lon" /) )    
    filevardef(fout, "monthly_climo_UQ_VERT_INTEGRAL" , typeof(UQ_VERT_INTEGRALm_mean)  , (/"time", "lat", "lon" /) )	 
    filevardef(fout, "monthly_climo_VQ_VERT_INTEGRAL" , typeof(VQ_VERT_INTEGRALm_mean)  , (/"time", "lat", "lon" /) )
    filevardef(fout, "monthly_climo_EVAP"             , typeof(EVAPm_mean)	        , (/"time", "lat", "lon" /) )
    filevardef(fout, "monthly_climo_TREFHT"           , typeof(TREFHTm_mean)	        , (/"time", "lat", "lon" /) )    
    filevardef(fout, "monthly_climo_Z500"             , typeof(Z500m_mean)	        , (/"time", "lat", "lon" /) )    
    filevardef(fout, "monthly_climo_Z850"             , typeof(Z850m_mean)	        , (/"time", "lat", "lon" /) )    

    filevardef(fout, "annual_climo_ICEFRAC"          , typeof(SIC_ann_mean)		 , (/"lat", "lon" /) )
    filevardef(fout, "annual_climo_SST"              , typeof(SST_ann_mean)		 , (/"lat", "lon" /) )
    filevardef(fout, "annual_climo_PRECIP"           , typeof(PRECIP_ann_mean)		 , (/"lat", "lon" /) )
    filevardef(fout, "annual_climo_SLP"              , typeof(SLP_ann_mean)		 , (/"lat", "lon" /) )
    filevardef(fout, "annual_climo_UQ" , typeof(UQ_ann_mean) , (/ "lev", "lat", "lon" /) )	      
    filevardef(fout, "annual_climo_VQ" , typeof(VQ_ann_mean) , (/ "lev", "lat", "lon" /) )    
    filevardef(fout, "annual_climo_UQ_VERT_INTEGRAL" , typeof(UQ_VERT_INTEGRAL_ann_mean) , (/"lat", "lon" /) )	      
    filevardef(fout, "annual_climo_VQ_VERT_INTEGRAL" , typeof(VQ_VERT_INTEGRAL_ann_mean) , (/"lat", "lon" /) )
    filevardef(fout, "annual_climo_EVAP"             , typeof(EVAP_ann_mean)		 , (/"lat", "lon" /) )
    filevardef(fout, "annual_climo_TREFHT"           , typeof(TREFHT_ann_mean)	         , (/"lat", "lon" /) )    
    filevardef(fout, "annual_climo_Z500"             , typeof(Z500_ann_mean)		 , (/"lat", "lon" /) )    
    filevardef(fout, "annual_climo_Z850"             , typeof(Z850_ann_mean)		 , (/"lat", "lon" /) )     
    
    filevardef(fout, "annual_var_ICEFRAC"           , typeof(SIC_ann_var)	       , (/"lat", "lon" /) )
    filevardef(fout, "annual_var_SST"               , typeof(SST_ann_var)	       , (/"lat", "lon" /) )
    filevardef(fout, "annual_var_PRECIP"            , typeof(PRECIP_ann_var)	       , (/"lat", "lon" /) )
    filevardef(fout, "annual_var_SLP"               , typeof(SLP_ann_var)	       , (/"lat", "lon" /) )
    filevardef(fout, "annual_var_UQ"  , typeof(UQ_ann_var) , (/"lev", "lat", "lon" /) )	       
    filevardef(fout, "annual_var_VQ"  , typeof(VQ_ann_var) , (/"lev", "lat", "lon" /) )    
    filevardef(fout, "annual_var_UQ_VERT_INTEGRAL"  , typeof(UQ_VERT_INTEGRAL_ann_var) , (/"lat", "lon" /) )	       
    filevardef(fout, "annual_var_VQ_VERT_INTEGRAL"  , typeof(VQ_VERT_INTEGRAL_ann_var) , (/"lat", "lon" /) )
    filevardef(fout, "annual_var_EVAP"              , typeof(EVAP_ann_var)	       , (/"lat", "lon" /) )	
    filevardef(fout, "annual_var_TREFHT"            , typeof(TREFHT_ann_var)	       , (/"lat", "lon" /) )    
    filevardef(fout, "annual_var_Z500"              , typeof(Z500_ann_var)	       , (/"lat", "lon" /) )    
    filevardef(fout, "annual_var_Z850"              , typeof(Z850_ann_var)	       , (/"lat", "lon" /) )    
    
    ;assign variablte attributes to nc variables.  This can be done explicitly, or by using filevarattdef to copy variable attributes over to nc variable definitions		
    lat@units="degreeN"
    lat@longname="latitude"
    filevarattdef(fout,"lat"    ,lat)
    lon@units="degreeE"
    lon@longname="longitude"	                  
    filevarattdef(fout,"lon"    ,lon)
    lev@units="hybrid_levels"
    lev@longname="levels"	                  
    filevarattdef(fout,"lev"    ,lev)   

    setfileoption(fout,"DefineMode",False)

    fout->nSamples = nSamples
    fout->time    = (/time/)
    fout->date    = (/date/)
    fout->datesec = (/datesec/)
    fout->lat     = (/lat/)
    fout->lon     = (/lon/)
    fout->lev     = (/lev/)    
    fout->monthly_climo_ICEFRAC              = (/SICm_mean/) 
    fout->monthly_climo_SST                  = (/SSTm_mean/)
    fout->monthly_climo_PRECIP               = (/PRECIPm_mean/)
    fout->monthly_climo_SLP                  = (/SLPm_mean/)
    fout->monthly_climo_UQ                   = (/UQm_mean/)    
    fout->monthly_climo_VQ                   = (/VQm_mean/)    
    fout->monthly_climo_UQ_VERT_INTEGRAL     = (/UQ_VERT_INTEGRALm_mean/)    
    fout->monthly_climo_VQ_VERT_INTEGRAL     = (/VQ_VERT_INTEGRALm_mean/)
    fout->monthly_climo_EVAP                 = (/EVAPm_mean/)
    fout->monthly_climo_TREFHT               = (/TREFHTm_mean/)    
    fout->monthly_climo_Z500                 = (/Z500m_mean/)    
    fout->monthly_climo_Z850                 = (/Z850m_mean/)    
    
    fout->annual_climo_ICEFRAC               = (/SIC_ann_mean/) 
    fout->annual_climo_SST                   = (/SST_ann_mean/)
    fout->annual_climo_PRECIP                = (/PRECIP_ann_mean/)
    fout->annual_climo_SLP                   = (/SLP_ann_mean/)
    fout->annual_climo_UQ                    = (/UQ_ann_mean/)	  
    fout->annual_climo_VQ                    = (/VQ_ann_mean/)	    
    fout->annual_climo_UQ_VERT_INTEGRAL      = (/UQ_VERT_INTEGRAL_ann_mean/)	  
    fout->annual_climo_VQ_VERT_INTEGRAL      = (/VQ_VERT_INTEGRAL_ann_mean/)	  
    fout->annual_climo_EVAP                  = (/EVAP_ann_mean/)
    fout->annual_climo_TREFHT                = (/TREFHT_ann_mean/)    
    fout->annual_climo_Z500                  = (/Z500_ann_mean/)    
    fout->annual_climo_Z850                  = (/Z850_ann_mean/)    
    
    fout->annual_var_ICEFRAC                 = (/SIC_ann_var/) 
    fout->annual_var_SST                     = (/SST_ann_var/)
    fout->annual_var_PRECIP                  = (/PRECIP_ann_var/)
    fout->annual_var_SLP                     = (/SLP_ann_var/)	
    fout->annual_var_UQ	                     = (/UQ_ann_var/)    
    fout->annual_var_VQ	                     = (/VQ_ann_var/)     
    fout->annual_var_UQ_VERT_INTEGRAL	     = (/UQ_VERT_INTEGRAL_ann_var/)    
    fout->annual_var_VQ_VERT_INTEGRAL	     = (/VQ_VERT_INTEGRAL_ann_var/)    
    fout->annual_var_EVAP	             = (/EVAP_ann_var/) 
    fout->annual_var_TREFHT	             = (/TREFHT_ann_var/)     
    fout->annual_var_Z500	             = (/Z500_ann_var/)     
    fout->annual_var_Z850	             = (/Z850_ann_var/)
    delete(fout)
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

