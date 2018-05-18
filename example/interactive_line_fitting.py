import os

from fit_spectral_lines import fit_feature
import line_analysis_BSNIP

for iline in line_list.keys():
            spectrum = line_analysis_BSNIP.read_iraf_spectrum(filename)
            min_list, pew_list, fig = fit_feature(spectrum, iline, binsize, absorption=True, 
            similar_widths=True, fixed_offset=True, offsets = [44, 164],
            input_filename=os.path.join(OUTPUT_DIR, '{}_input.yaml'.format(iline)))
