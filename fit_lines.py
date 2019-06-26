import argparse
import os

from astropy.io import fits, ascii as asc
from astropy.table import Table, Column
from astropy.wcs import WCS

from matplotlib import pyplot as plt
plt.ion()

from utilities_az import spectroscopy as spec
import fit_spectral_lines

def get_spectrum(filename, ext=None):
    '''
    Requires fits binary tables to have column names that contain 'wave' and 'flux'
    
    Assumes zero indexing for assigning wavelengths
    '''
    wave = []
    flux = []
    if '.fits' in filename:
        ofile = fits.open(filename)
        data = ofile[ext].data    
        #Table data
        if isinstance(data, (fits.BinTableHDU, fits.fitsrec.FITS_rec)):
            colnames = data.colnames
            for icol in colnames:
                if 'wave' in icol.lower():
                    wave = data[icol]
                elif 'flux' in icol.lower():
                    flux = data[icol]
        else:
            try:
                hdr_wcs = WCS(ofile[0].header)
                hdr = ofile[0].header
            except:
                try:
                    hdr_wcs = WCS(ofile[ext].header)
                    hdr = ofile[ext].header
                except:
                    sys.exit('Unable to find wcs solution in header')
            if len(data.shape) == 3:
                flux = data[0,0,:]
            
            elif len(data.shape) == 1:
                flux = data
            else:
                sys.exit('Unable to read fits file format')
            pix = np.arange(len(flux))
            wave = spec.calc_wavelength(hdr, pix)
            
    else:
        tbdata = asc.read(filename)
        if tbdata.colnames[0] == 'col0':
            wave = tbdata['col0']
            flux = tbdata['col1']
        elif tbdata.colnames[0] == 'col1':
            wave = tbdata['col1']
            flux = tbdata['col2']
        else:
            for icol in tbdata.colnames:
                if 'wave' in icol:
                    wave = tbdata[icol]
                elif 'flux' in icol:
                    flux = tbdata['icol']
    
    if (len(wave)==0) or (len(flux)==0):
        sys.exit('Unable to read fits file format')
    else:
        spectrum = spec.spectrum1d(wave, flux)
    return spectrum
        


parser = argparse.ArgumentParser(description='Fit spectral lines')
#Required
parser.add_argument('filename', 
                    help='filename (and path) to spectrum that you want to fit')
#Default values
parser.add_argument('--ext', default=0, type=int,
                    help='extension of fits file that contains the spectrum')
parser.add_argument('-p', '--profile', default='a', type=str,
                    help='profile type, either "a" for absorption or "e" for emission')
parser.add_argument('-z', '--redshift', default=0, type=float,
                    help='redshift of object to be used to correct wavelengths to rest wavelengths')
parser.add_argument('--similar_widths', action='store_false',
                    help='flag to require each component in a multi-component fit to have the same width')

parser.add_argument('--fixed_offset', type=float,
                    help='array of offsets between the centers of each component of a multi-component fit')
parser.add_argument('--search_range', type=float,
                    help='the range (from the location selected by the user) of allowed values for the center of a component')
parser.add_argument('--output', type=str, 
                    help='name of output file')
parser.add_argument('--overwrite', action='store_true',
                    help='flag to overwrite rather than appended to output table')
parser.add_argument('--rest_wave', type=float,
                    help='rest wavelength of feature')
parser.add_argument('-n', '--name', type=str, 
                    help='name of feature that is being fit')
parser.add_argument('--def_fit_range', action='store_true', 
                    help='flag to define separately the continuum edges and the fit range')
                    
args = parser.parse_args()

#Read file
spectrum = get_spectrum(args.filename, args.ext)
spectrum.wave = spec.apply_redshift(spectrum.wave, args.redshift)

#Create object for fitting
spec_feat = fit_spectral_lines.spectral_feature(spectrum, os.path.basename(args.filename))

#Define type of profile being fit
if args.profile == 'a':
    absorption=True
elif args.profile == 'e':
    absorption = False
else:
    sys.exit('--profile must be either "a" or "e"')

#define offsets if any
if args.fixed_offset is not None:
    fixed_offset = True
    offsets = args.fixed_offset
else:
    fixed_offset=False
    offsets=None

#Define output filename
if args.output is not None:
    output_filename = args.output
else:
    output_filename='fit_features.csv'

#Read in previous table or define new table for writing output
if (args.overwrite == False) and (os.path.exists(output_filename)):
    tbdata = asc.read(output_filename)
else:
    tbdata = Table(names=['name', 'component','fwhm', 'absorption',
                      'min_wave', 'min_wave_err_L', 'min_wave_err_R', 
                      'pew', 'pew_err', 'filename', 'redshift'],
                      dtype=('S15', 'i8', 'f8', 'S2',
                             'f8', 'f8', 'f8', 
                             'f8', 'f8', 'S50', 'f8'))
#Add a column for rest wavelength if provided
if args.rest_wave is not None:
    tbdata.add_column(Column(name='rest_wave', dtype='f8'))
    
#define line name
if args.name is not None:
    name = args.name
else:
    name='line1'
    linenum = 1
    while name in tbdata['name']:
        name = 'line{}'.format(linenum)
        linenum+=1
        
spec_feat.define_feature(name, 
                         absorption=absorption, 
                         overwrite=True, #overwrite the yaml files with fit info
                         similar_widths=args.similar_widths, 
                         search_range=args.search_range,
                         fixed_offset=fixed_offset, 
                         offsets=offsets,
                         define_fit_range = args.def_fit_range)
                         
for indx, (fwhm, min_tup, pew_tup) in enumerate(zip(spec_feat.fwhm_list, 
                                                    spec_feat.min_list, 
                                                    spec_feat.pew_list)):
    row = [name, indx, fwhm, args.profile, 
           min_tup[0], min_tup[1], min_tup[2],
           pew_tup[0], pew_tup[1], 
           os.path.basename(args.filename), args.redshift]
    if args.rest_wave is not None:
        row.append(args.rest_wave)
    tbdata.add_row(row)
    
tbdata.write(output_filename, overwrite=True) #already appended if args.overwrite is False
spec_feat.interactive_fig.savefig(spec_feat.output_filename.replace('yaml', 'pdf'))
plt.close(spec_feat.interactive_fig)



