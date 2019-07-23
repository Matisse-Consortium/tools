#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 2019

@author: FMi
"""
import matplotlib.pyplot as plt
import argparse
import numpy as np
from astropy.io import fits
import sys
import mat_show_oifits as msoi
import os
from matplotlib.backends.backend_pdf import PdfPages


inch=1/2.54

if __name__ == '__main__':
    print("Starting...")

    #--------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Plots the oi data as a function of spatial frequency.')

    #--------------------------------------------------------------------------
    parser.add_argument('in_dir', metavar='in_dir', type=str, \
    help='The path to the directory containing your oifits data.', default=None)

    #--------------------------------------------------------------------------
    parser.add_argument('--log', metavar='log', \
                        help='Set logscale or not', default=False)

    #--------------------------------------------------------------------------
    parser.add_argument('--showvis', metavar='showvis', type=str, \
                        help='Plot visibilities instead of V squared.', default=False)



    
    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_showOiFreq.py --help to be kind with you:\033[0m\n")
        parser.print_help()
	sys.exit(0)

    list_of_dicts          = msoi.open_oi_dir(args.in_dir)
    filtered_list_of_dicts = msoi.filter_oi_list(list_of_dicts)

    fig=plt.figure(figsize=(29.7*inch,21*inch))
    plts={}
    
    pltv2=[]
    plts['VIS2_{0}'.format(0)]=fig.add_axes([0.075, 0.55, 0.85,0.4])
    pltv2.append(plts['VIS2_{0}'.format(0)])
    
    pltcp=[]
    plts['CP_{0}'.format(0)]=fig.add_axes([0.075, 0.08, 0.85,0.4],sharex=plts['VIS2_0'])
    pltcp.append(plts['CP_{0}'.format(0)])
    
    print("plotting data...")
    for idic in list_of_dicts:
        msoi.show_oi_vs_freq(idic, log=args.log,showvis=args.showvis, subplotV2=pltv2, subplotCP=pltcp)
    
    plt.show()
