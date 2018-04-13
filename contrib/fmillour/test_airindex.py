from airindex import computeAirIndexMathar
from airindex import computeAirIndexCiddor
import numpy as np
import matplotlib.pyplot as plt

print("Coucou c'est moi")

P = 1.e5;
T = 300.;
H = 0.;

wlen = np.linspace(1.3,2.5,30)*1e-6;
n1 = computeAirIndexMathar(wlen, P, T, H);
plt.plot(wlen,n1,color="red");

wlen = np.linspace(2.8,4.2,30)*1e-6;
n1 = computeAirIndexMathar(wlen, P, T, H);
plt.plot(wlen,n1,color="red");

wlen = np.linspace(4.35,5.2,30)*1e-6;
n1 = computeAirIndexMathar(wlen, P, T, H);
plt.plot(wlen,n1,color="red");

wlen = np.linspace(7.5,14.1,30)*1e-6;
n1 = computeAirIndexMathar(wlen, P, T, H);
plt.plot(wlen,n1,color="red");


    # CO2 content in Chile as of 31 dec. 2008. See
    # http://www.esrl.noaa.gov/gmd/ccgg/carbontracker/co2weather.php?type=glb#imagetable
    # for details
xc = 384;
#xc = 450;
wl = np.linspace(1.3,14.1,30)*1e-6;
print("Wlen")
print(wl)

n_prop = computeAirIndexCiddor(wl, P, T, H, xc)
print("N_prop2")
print(n_prop)
plt.plot(wl,n_prop,color="blue")




H = 50.;

wlen = np.linspace(1.3,2.5,30)*1e-6;
n1 = computeAirIndexMathar(wlen, P, T, H);
plt.plot(wlen,n1,color="green");

wlen = np.linspace(2.8,4.2,30)*1e-6;
n1 = computeAirIndexMathar(wlen, P, T, H);
plt.plot(wlen,n1,color="green");

wlen = np.linspace(4.35,5.2,30)*1e-6;
n1 = computeAirIndexMathar(wlen, P, T, H);
plt.plot(wlen,n1,color="green");

wlen = np.linspace(7.5,14.1,30)*1e-6;
n1 = computeAirIndexMathar(wlen, P, T, H);
plt.plot(wlen,n1,color="green");



n_prop = computeAirIndexCiddor(wl, P, T, H, xc)
plt.plot(wl,n_prop,color="cyan")




plt.show()
