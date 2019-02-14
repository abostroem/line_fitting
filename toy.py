import os

from matplotlib import pyplot as plt

from astropy.io import fits, ascii as asc
from astropy.table import Table, Column

from utilities_az import spectroscopy as spec
from fit_spectral_lines import define_feature

filename = '../data/2014G/2014G_2014-03-17_19-09-29_Galileo_BC-Asi_None.dat'

tbdata = asc.read(filename, names=['wave', 'flux'])
error = tbdata['flux']*0.05
spectrum = spec.spectrum1d(tbdata['wave'], tbdata['flux'], error)
line_name = 'A'

min_list, pew_list, fig = define_feature(spectrum, line_name,
                                         absorption=True, 
                                         similar_widths=True, 
                                         interactive = True,
                                         define_fit_range=True)
                                         
min_wave = [i[0] for i in min_list]
min_wave_std = [j[1].value for j in min_list]
min_col = Column(min_wave, name='min')
min_err_col = Column(min_wave_std, name='min_err')

pew = [i[0] for i in pew_list]
pew_err = [j[1] for j in pew_list]
pew_col = Column(pew, name='pew')
pew_err_col = Column(pew_err, name='pew_err')
output = Table([min_col, min_err_col, pew_col, pew_err_col])
output.write(os.path.basename(filename).split('.')[0]+'line_fit_{}.csv'.format(line_name), format='csv')

fig.savefig(os.path.basename(filename).split('.')[0]+'line_fit_{}.pdf'.format(line_name))