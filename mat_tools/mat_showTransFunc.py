#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file is part of the Matisse pipeline GUI series
Copyright (C) 2017- Observatoire de la Côte d'Azur

Created on Tue Jul 22 2019
@author: FMi

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.

You can use, modify and/ or redistribute the software under the terms
of the CeCILL license as circulated by CEA, CNRS and INRIA at the
following URL "http://www.cecill.info". You have a copy of the licence
in the LICENCE.md file.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
"""
import matplotlib.pyplot as plt
import argparse
import numpy as np
from astropy.io import fits
import sys
#import mat_show_oifits as msoi
import libShowOifits as msoi
import os
from matplotlib.backends.backend_pdf import PdfPages


inch=1/2.54

if __name__ == '__main__':
    print("Starting...")

    #--------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Plots the transfer function of a night of MATISSE data from its oifits files.')

    #--------------------------------------------------------------------------
    parser.add_argument('in_dir', metavar='in_dir', type=str, \
    help='The path to the directory containing your oifits data.', default=None)
    parser.add_argument('wlmin',  type=float, help='Minimum wavelength', default=3.5)
    parser.add_argument('wlmax',  type=float, help='Maximum wavelength', default=4.2)
    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_ashowTransFunc.py --help to be kind with you:\033[0m\n")
        parser.print_help()
	sys.exit(0)

    list_of_dicts          = msoi.open_oi_dir(args.in_dir)
    filtered_list_of_dicts = msoi.filter_oi_list(list_of_dicts)
    
    fig=plt.figure(figsize=(29.7*inch,19*inch))

    plts={}

    print("Setting V2 windows...")
    #####################################################
    pltv2=[]
    for i in range(6):
        if i==0:
            plts['VIS2_{0}'.format(i)]=fig.add_axes([0.075, 0.15*(i+0.5), 0.35,0.15])
        else:
            plts['VIS2_{0}'.format(i)]=fig.add_axes([0.075, 0.15*(i+0.5), 0.35,0.15],sharex=plts['VIS2_0'],sharey=plts['VIS2_0'])
            plts['VIS2_{0}'.format(i)].get_xaxis().set_visible(False)
        pltv2.append(plts['VIS2_{0}'.format(i)])

    print("Setting CP windows...")
    #####################################################
    pltcp=[]
    for i in range(4):
        if i==0:
            plts['CP_{0}'.format(i)]=fig.add_axes([0.575, 0.2*(i+0.5), 0.35,0.2],sharex=plts['VIS2_0'])
        else:
            plts['CP_{0}'.format(i)]=fig.add_axes([0.575, 0.2*(i+0.5), 0.35,0.2],sharex=plts['VIS2_0'],sharey=plts['CP_0'])
            plts['CP_{0}'.format(i)].get_xaxis().set_visible(False)
        pltcp.append(plts['CP_{0}'.format(i)])
        
    print("plotting...")        
    #####################################################
    msoi.show_oi_vs_time(filtered_list_of_dicts ,[args.wlmin,args.wlmax],key="VIS2", datatype="VIS2",subplotList=pltv2,calColor='lightgray')
    
    msoi.show_oi_vs_time(filtered_list_of_dicts ,[args.wlmin,args.wlmax],key="TF2", datatype="TF2",subplotList=pltv2,calColor='blue')
    
    msoi.show_oi_vs_time(filtered_list_of_dicts ,[args.wlmin,args.wlmax],key="T3", datatype="CLOS",subplotList=pltcp)
    
    plt.show()
