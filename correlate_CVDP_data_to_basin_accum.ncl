;Load CVDP data
NYearsTotal=1800; ;Calculate number of years based on years CVDP analyzed.
NMonthsTotal=NYearsTotal*12;
in=addfile("/glade/p/cesmLE/CESM-CAM5-BGC-LE/CVDP/Controls/CESM1_LENS_Coupled_Control.cvdp_data.401-2200.nc","r")
vNames=getfilevarnames(in)
nNames=dimsizes(vNames)
YearlyCVDPTimeSeries=NewList("fifo")
do n=0,nNames-1
   ds=getfilevardimsizes(in,vNames(n))
   if (dimsizes(ds).eq.1) then
      
      if (ds.eq.NYearsTotal) then
         StringTest=str_match_ic_regex(in->$vNames(n)$@long_name,"area-average")
	 print(StringTest)
         if  (ismissing(StringTest)) then
            ListPush(YearlyCVDPTimeSeries,in->$vNames(n)$)
	 end if
      end if
    end if
   delete(ds)
end do
delete(in)
;Load Basin-based accumulation timeseries
in=addfile("compositing_work/output/BasinAccumulationTimeSeries.nc","r")
BasinTS=tofloat(in->time_series(1::,:))

;For each basin, do cross correlation to climate indices
ds=dimsizes(BasinTS)
nBasins=ds(1)
nCVDPYearTS=ListCount(YearlyCVDPTimeSeries)
mxlag=2
do b=0,nBasins-1
   ccr=new((/nCVDPYearTS,mxlag+1/),float)
   dims=dimsizes(ccr)
   do nCVDP=0,nCVDPYearTS-1
      ccr(nCVDP,:)=esccr(YearlyCVDPTimeSeries[nCVDP],BasinTS(:,b),mxlag); build up cross correlations at all positive lags (ie BasinTS lags CVDP time series)
   end do
   ccr1d=ndtooned(ccr)
   inds = ind_resolve(maxind(ccr1d),dims)
   print("   ")
   print("***BASIN "+(b+1)+"***")
   print("Most correlated CDVP Timeseries="+YearlyCVDPTimeSeries[inds(0,0)]@long_name)
   print("Max correlation="+ccr(inds(0,0),inds(0,1)))
   print("Lag of max correlation="+inds(0,1))
   delete(ccr)
end do
