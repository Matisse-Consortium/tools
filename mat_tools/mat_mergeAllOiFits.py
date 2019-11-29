#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 13:50:34 2019

@author: ame
"""



import sys
import argparse
from libPostTools import mat_mergeByTplStart



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merging all oifits files from the same template into a single OiFits file per detector.')
    parser.add_argument('dirIn', default="",help='Input directory')
    parser.add_argument('--dirOut', default="MERGED",help='Output directory')
    try:
        args = parser.parse_args()
    except:
        print("\n\033[93mRunning mat_MergeAllOiFits.py --help to be kind with you:\033[0m\n")
        parser.print_help()
        print("\n     Example : python mat_MergeAllOiFits.py . ")
        sys.exit(0)


    mergedData=mat_mergeByTplStart(args.dirIn,save=True,dirOut=args.dirOut)


