import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy import stats

f=nc.Dataset("/glade/p/cesmLE/CESM-CAM5-BGC-LE/CVDP/Controls/CESM1_LENS_Coupled_Control.cvdp_data.401-2200.nc")
#Load CVDP diagnostics for CESM LE run
index=np.zeros((27,1800))
nino34=f.variables["nino3"][:].reshape((1800,12))
index[0,:]=np.mean(nino34,axis=1)
index[1,:]=f.variables["sam_pc_ann"][:]
index[2,:]=f.variables["psa1_pc_ann"][:]
index[3,:]=f.variables["psa2_pc_ann"][:]

f=nc.Dataset("compositing_work/output/BasinAccumulationTimeSeries.nc")
BasinTS=f.variables["time_series"][1::,:]
BasinTSmean=np.mean(BasinTS,axis=1)
for k in range(BasinTS.shape[1]):
   MeanVal=np.mean(BasinTS[:,k])
   BasinTS[:,k]=BasinTS[:,k]-MeanVal

f=nc.Dataset("compositing_work/output/AISAccumulationTimeSeries.nc")
AISTS=f.variables["time_series"][1::]

time=np.arange(100)
y=BasinTS[0:100,:]
#y=np.swapaxes(BasinTS[0:50,:],0,1)
plt.plot(time,y)
plt.savefig("figs/accumulation_contributions")

#Index=np.array(['Nino3','SAM','PSA1','PSA2'])
#
#for b in range(0,27):
   #print b
   #c=BasinTS[:,b]
   #plt.close("all") 
   #for i in range(0,4):
   #   ax=plt.subplot(2,2,i+1)
   #   plt.scatter(index[i,:],c,c="black",s=10,edgecolors='none',alpha=0.3)     
   #   ax.set_xlabel(Index[i])
   #   if (np.mod(i,2)<0.5):
   #      ax.set_ylabel('Accumulation (Gt/yr)')
   #   gradient, intercept, r_value, p_value, std_err =stats.linregress(index[i,:],c)
   #   ax.text(0.05, 0.9,"$r=%.2f$"%r_value,transform=ax.transAxes)
   #   ax.text(0.6, 0.9,"$p=%.2f$"%p_value,transform=ax.transAxes)      
   #   plt.axis('tight')
   #bp1=b+1
   #plt.savefig("figs/indice_regressions%02d"%bp1)
