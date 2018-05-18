import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import numpy as np

from astropy.table import Table, vstack
from astropy.io import fits
from astropy.io import ascii as asc

import yaml

import os
import line_analysis_BSNIP



FIG_DIR = '../figures'
OUTPUT_DIR = '../data/line_info'
input_file = '../data/line_vel_analysis_input.yml'
#plt.ioff()
append = False


feature_dict= yaml.load('''
ScII5239:
    name: 'ScII5239'
    wmin: 4800
    wmax: 4985
    slope: 0 #A/day
    texpl: 2457262. #JD
    smooth_param:
        deg: 2
        width: 5
    edge_param:
        binsize: 5
        binmax: 100
        concavity_binsize: 20
    directories:
        DATA_DIR_LCO: '../data/spectra/lco'
        DATA_DIR_EFOSC: '../data/spectra/EFOSC'
    files:
        - 'asassn15oz_20150916_redblu_120911.274.fits'
        - 'asassn-15oz_20150920_redblu_135034.512.fits'
        ''')

for iline in feature_dict.keys():
        for ifile in feature_dict[iline]['files']:
            for idir_key in feature_dict[iline]['directories'].keys():
                filename = os.path.join(feature_dict[iline]['directories'][idir_key], ifile)
                if os.path.exists(filename):
                    found = True
                    break
            if not found:
                print('Unable to locate file {}'.format(ifile))
            else:
                line_return = line_analysis_BSNIP.characterize_line(feature_dict[iline], filename, visualization_level=2)