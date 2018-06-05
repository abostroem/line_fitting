import os

from fit_spectral_lines import define_feature
import line_analysis_BSNIP

iline = 'CaII'

#A spectrum object can be any object with a wave, flux, and error attribute
filename = 'asassn15oz_20150904_redblu_122216.314.fits'
spectrum = line_analysis_BSNIP.read_iraf_spectrum(filename)

spectrum.__setattr__('filename',os.path.basename(filename))
min_list, pew_list, fig = define_feature(spectrum, iline, binsize, absorption=True, 
                                                     similar_widths=True, fixed_offset=True, offsets = [44, 164],
                                                     input_filename='{}_input.yaml'.format(iline), 
                                                     input_append = append, interactive=True)
