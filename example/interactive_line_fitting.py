import os

from fit_spectral_lines import define_feature
import line_analysis_BSNIP

iline = 'CaII'

min_list, pew_list, fig = define_feature(spectrum, iline, binsize, absorption=True, 
                                                     similar_widths=True, fixed_offset=True, offsets = [44, 164],
                                                     input_filename='{}_input.yaml'.format(iline), 
                                                     input_append = append, interactive=True)
