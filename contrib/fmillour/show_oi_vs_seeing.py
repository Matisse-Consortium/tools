from mat_show_oifits import *
import os
from os.path import expanduser

home = expanduser("~")

path = os.path.join(home,"COM1B/OIFITS/OIFITS_2018-07-17/2018-07-13_OIFITS")
path = os.path.join(home,"COM1B/OIFITS/ALL_AT/compact_quad")
print(path)

list_of_dic = open_oi_dir(path)


filt_dic = filter_oi_list(list_of_dic, dates=[], bands="L", spectral_resolutions="LOW", DIT_range=[0.07,0.08], targets=[],BCD1='IN',BCD2='IN')
filt_dic = filter_oi_list(list_of_dic, dates=[], bands="L", spectral_resolutions="LOW", DIT_range=[0.07,0.08], targets=[],BCD1='OUT',BCD2='OUT')


#show_oi_vs_anything(filt_dic, [3.2,3.7], key="TF2", datatype="TF2", xaxis="HIERARCH ESO ISS AMBI FWHM START", showvis=False,plot_errorbars=True)

#show_oi_vs_anything(filt_dic, [3.2,3.7], key="TF2", datatype="TF2", xaxis="HIERARCH ESO ISS AMBI TAU0 START", showvis=False,plot_errorbars=True)

show_oi_vs_time(filt_dic, [3.2,3.7], key="TF2", datatype="TF2", showvis=False,plot_errorbars=True)
