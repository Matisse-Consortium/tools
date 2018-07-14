## calibOB.py
##
## PURPOSE
##    write a MATISSE "LOLO" OB for a calibrator from the MSDFCC
##
## USE
##    python3 calibOB.py star_name
##
## PARAMETERS
##    star_name   Any name by which Simbad knows this star
##
## HISTORY
##    2018-07-13   leo     created
##
import sys, os, pickle
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy import coordinates
from astropy import units as u
import numpy as np

objectname = sys.argv[1]
OB_name = objectname+".obx"

print("Generating " + OB_name)

msdfcc = "msdfcc-v5.fits"
hdu = fits.open(msdfcc)
cat = hdu[1].data
##
## generate list of coordinates from MSDFCC
## this will take a few minutes the first time this script is run and it will then be cached for this version of the MSDFCC
coord_cache_file = msdfcc.split(".fits")[0]+".pickle"
if not os.path.isfile(coord_cache_file):
	print("Cached coordinate file not found; generating it. This will take a few minutes...")
	ra = cat["RAJ2000"]
	dec = cat["DEJ2000"]
	c = ra + " " + dec
	coords = coordinates.SkyCoord(c, unit=(u.hourangle,u.deg))
	with open(coord_cache_file,"wb") as f:
		pickle.dump(coords,f)
else:
	with open(coord_cache_file,"rb") as f:
		coords = pickle.load(f)
##
## resolve object name with simbad
query = Simbad.query_object(objectname)
ra = query["RA"].data.data[0]
dec = query["DEC"].data.data[0]
simbad_id = query["MAIN_ID"].data.data[0].decode()
coords_star_string = ra + " " + dec
coords_star = coordinates.SkyCoord(coords_star_string,unit=(u.hourangle,u.deg))
##
## match by coord with msdfcc
sep = coords.separation(coords_star)
sep_threshold = coordinates.angles.Angle(1*u.arcsec)

if np.min(sep) < sep_threshold:
	rec_id = np.argmin(sep)
else:
	raise ValueError("Could not find a calibrator within 1 arcsec of " + objectname)

##
## retrieve record from msdfcc
rec =  hdu[1].data[rec_id]

##
## write to obx file
with open(OB_name,"w") as f:
	f.write("IMPEX.VERSION \"2.0\"\n")
	name_str = "{0} L {1:5.2f} Jy N {2:5.2f} Jy LOLO".format(objectname, rec["median_Lflux"], rec["median_Nflux"])
	f.write("name                             \""+name_str+"\"\n")
	userComments_str = "L={0:5.2f} Jy N = {1:5.2f}; Simbad id = {2}".format(rec["median_Lflux"], rec["median_Nflux"], simbad_id)
	f.write("userComments                     \""+userComments_str+"\"\n")
	f.write("InstrumentComments               \"\"\n")
	f.write("userPriority                     \"1\"\n")
	f.write("type                             \"O\"\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("TARGET.NAME                      \""+objectname+"\"\n")
	propRA_str = "{0:7.5f}".format(rec["pmRA"]/1000)
	f.write("propRA                           \""+propRA_str+"\"\n")
	propDEC_str = "{0:7.5f}".format(rec["pmDEC"]/1000)
	f.write("propDec                          \""+propDEC_str+"\"\n")
	f.write("diffRA                           \"0.000000\"\n")
	f.write("diffDec                          \"0.000000\"\n")
	f.write("equinox                          \"J2000\"\n")
	f.write("epoch                            \"2000.0\"\n")
	f.write("ra                               \""+rec["RAJ2000"]+"\"\n")
	f.write("dec                              \""+rec["DEJ2000"]+"\"\n")
	f.write("objectClass                      \"                        \"\n")
	f.write("\n")
	f.write("\n")
	f.write("CONSTRAINT.SET.NAME              \"No name\"\n")
	f.write("seeing                           \"0.0\"\n")
	f.write("sky_transparency                 \"Photometric\"\n")
	f.write("air_mass                         \"5.0\"\n")
	f.write("fractional_lunar_illumination    \"1.0\"\n")
	f.write("moon_angular_distance            \"30\"\n")
	f.write("strehlratio                      \"0.0\"\n")
	f.write("twilight                         \"0\"\n")
	f.write("watervapour                      \"0.0\"\n")
	f.write("atm                              \"no constraint\"\n")
	f.write("contrast                         \"0.0\"\n")
	f.write("\n")
	f.write("\n")
	f.write("OBSERVATION.DESCRIPTION.NAME \""+objectname+" LOLO\"\n")
	f.write("instrument \"MATISSE\"\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("ACQUISITION.TEMPLATE.NAME \"MATISSE_img_acq\"\n")
	f.write("DET1.READ.CURNAME      \"SCI-SLOW-SPEED\"\n")
	f.write("DET2.READ.CURNAME      \"SCI-HIGH-GAIN\"\n")
	f.write("SEQ.ACQ.DET1.DIT       \"0.005\"\n")
	f.write("SEQ.ACQ.DET2.DIT       \"0.005\"\n")
	f.write("SEQ.ACQ.INS.FIL.NAME   \"L-\"\n")
	f.write("SEQ.ACQ.INS.FIN.NAME   \"N-\"\n")
	f.write("SEQ.ACQ.SKY.DURATION   \"30\"\n")
	f.write("SEQ.ACQ.ST             \"T\"\n")
	f.write("SEQ.ACQ.TARG.DURATION  \"10\"\n")
	f.write("SEQ.CHOP.ST            \"T\"\n")
	f.write("SEQ.DIL.WL0            \"3.5\"\n")
	flux_L_str = "{0:5.2f}".format(rec["median_Lflux"])
	f.write("SEQ.FLUX.L             \""+flux_L_str+"\"\n")
	flux_N_str = "{0:5.2f}".format(rec["median_Nflux"])
	f.write("SEQ.FLUX.N             \""+flux_N_str+"\"\n")
	f.write("SEQ.FS.DET1.DIT        \"0.1\"\n")
	f.write("SEQ.FS.DET2.DIT        \"0.02\"\n")
	f.write("SEQ.FS.INS.DIL.NAME    \"LOW\"\n")
	f.write("SEQ.FS.INS.DIN.NAME    \"LOW\"\n")
	f.write("SEQ.FS.INS.FIL.NAME    \"LM\"\n")
	f.write("SEQ.FS.INS.FIN.NAME    \"OPEN\"\n")
	f.write("SEQ.FS.INS.SFL.NAME    \"1.50\"\n")
	f.write("SEQ.FS.INS.SFN.NAME    \"2.00\"\n")
	f.write("SEQ.FS.SKY.DURATION    \"30\"\n")
	f.write("SEQ.FS.ST              \"T\"\n")
	f.write("SEQ.FS.TARG.DURATION   \"60\"\n")
	f.write("SEQ.LADC.ST            \"F\"\n")
	f.write("SEQ.OPDM.L.ST          \"F\"\n")
	f.write("SEQ.OPDM.N.ST          \"T\"\n")
	f.write("SEQ.PUP.DET1.DIT       \"0.06\"\n")
	f.write("SEQ.PUP.SKY.DURATION   \"10\"\n")
	f.write("SEQ.PUP.ST             \"F\"\n")
	f.write("SEQ.PUP.TARG.DURATION  \"10\"\n")
	f.write("SEQ.SKY.OFFS.ALPHA     \"1.0\"\n")
	f.write("SEQ.SKY.OFFS.DELTA     \"15.0\"\n")
	f.write("SEQ.TADC.ST            \"F\"\n")
	f.write("SEQ.TRACK.BAND         \"L\"\n")
	f.write("SEQ.TRACK.ST           \"T\"\n")
	f.write("TEL.TARG.ADDVELALPHA   \"0\"\n")
	f.write("TEL.TARG.ADDVELDELTA   \"0\"\n")
	f.write("TEL.GS1.ALPHA          \"0.\"\n")
	f.write("TEL.GS1.DELTA          \"0.\"\n")
	f.write("TEL.GS1.MAG            \"7.6\"\n")
	f.write("TEL.CHOP.FREQ          \"0.5\"\n")
	f.write("TEL.CHOP.POSANG        \"0.0\"\n")
	f.write("TEL.CHOP.PVRATIO       \"1.0\"\n")
	f.write("TEL.CHOP.THROW         \"5\"\n")
	f.write("COU.AG.ALPHA           \"0.\"\n")
	f.write("COU.AG.DELTA           \"0.\"\n")
	f.write("COU.AG.EPOCH           \"2000\"\n")
	f.write("COU.AG.EQUINOX         \"2000\"\n")
	f.write("COU.AG.GSSOURCE        \"SCIENCE\"\n")
	f.write("COU.AG.PMA             \"0.0\"\n")
	f.write("COU.AG.PMD             \"0.0\"\n")
	f.write("COU.AG.TYPE            \"DEFAULT\"\n")
	f.write("COU.GS.MAG             \"7.6\"\n")
	f.write("DEL.REF.OPL            \"30\"\n")
	f.write("ISS.IAS.HMAG           \"5.2\"\n")
	f.write("\n")
	f.write("\n")
	f.write("TEMPLATE.NAME \"MATISSE_hyb_obs\"\n")
	f.write("DET1.DIT              \"0.075\"\n")
	f.write("DET1.READ.CURNAME     \"SCI-SLOW-SPEED\"\n")
	f.write("DET2.DIT              \"0.02\"\n")
	f.write("DET2.READ.CURNAME     \"SCI-HIGH-GAIN\"\n")
	f.write("SEQ.CHOP.ST           \"T\"\n")
	f.write("SEQ.DIL.WL0           \"3.5\"\n")
	f.write("SEQ.FRINGES.BCD.SEQ   \"2STEPS\"\n")
	f.write("SEQ.FRINGES.DURATION  \"60\"\n")
	f.write("SEQ.FRINGES.NCYCLES   \"2\"\n")
	f.write("SEQ.LADC.ST           \"F\"\n")
	f.write("SEQ.OPDM.L.ST         \"F\"\n")
	f.write("SEQ.OPDM.N.ST         \"T\"\n")
	f.write("SEQ.PHOTO.DURATION    \"0\"\n")
	f.write("SEQ.PHOTO.ST          \"T\"\n")
	f.write("SEQ.SKY.DURATION      \"60\"\n")
	f.write("SEQ.SKY.OFFS.ALPHA    \"1.0\"\n")
	f.write("SEQ.SKY.OFFS.DELTA    \"15.0\"\n")
	f.write("SEQ.TADC.ST           \"F\"\n")
	f.write("SEQ.TRACK.BAND        \"L\"\n")
	f.write("SEQ.TRACK.ST          \"T\"\n")
	f.write("TEL.CHOP.FREQ         \"1\"\n")
	f.write("TEL.CHOP.POSANG       \"0.0\"\n")
	f.write("TEL.CHOP.PVRATIO      \"1.0\"\n")
	f.write("TEL.CHOP.THROW        \"5\"\n")
	f.write("INS.DIL.NAME          \"LOW\"\n")
	f.write("INS.DIN.NAME          \"LOW\"\n")
	f.write("INS.FIL.NAME          \"L\"\n")
	f.write("INS.FIN.NAME          \"OPEN\"\n")
	f.write("INS.POL.NAME          \"OPEN\"\n")
	f.write("INS.PON.NAME          \"OPEN\"\n")
	f.write("INS.SFL.NAME          \"1.50\"\n")
	f.write("INS.SFN.NAME          \"2.00\"\n")
	f.write("DPR.CATG              \"CALIB\"\n")

print("Successfully generated OB for " + objectname)