'''
Fit single or multi-component emission or absorption lines.
Example call:
fit = fit_spectral_lines.fit_feature(spec, similar_widths=True)
Where spec is an object with attributes wave and flux. This can be 
created with the spectrum1d class
'''
import sys
import numpy as np

from astropy.modeling import models, fitting
from matplotlib import pyplot as plt
import matplotlib as mpl

import yaml

import line_analysis_BSNIP

PROFILES = {'g': (models.Gaussian1D, {'size_arg':'amplitude', 'mean_arg': 'mean', 'width_arg': 'stddev'}),
            'l': (models.Lorentz1D, {'size_arg':'amplitude', 'mean_arg': 'x_0', 'width_arg': 'fwhm'}),
            'm': (models.Moffat1D, {'size_arg':'amplitude', 'mean_arg': 'x_0', 'width_arg': 'fwhm'}),
            }

class spectrum1d(object):
    def __init__(self, wavelength, flux, error=None):
        self.wave = wavelength
        self.flux = flux
        self.error = error

def plot_spectrum(spectrum):
    '''
    Plot a spectrum
    Input:
        spectrum: spectrum1d class object with wave anf flux attributes
    Output:
        ax: subplot object
    '''
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.set_position([0.1, 0.1, 0.8, 0.2])
    ax1.set_position([0.1, 0.35, 0.8, 0.55])
    ax1.plot(spectrum.wave, spectrum.flux)
    ax2.set_xlabel('Wavelength $\AA$')
    ax1.set_ylabel('Flux')
    ax2.set_ylabel('Residual')
    return ax1, ax2

def continuum_normalization(ax, spectrum, absorption=True):
    '''
    Allow the user to select two points which will be fit with a straight line
    and used to normalize a feature
    Input:
        ax: subplot or axis object
        spectrum: spectrum1d class object with wave and flux attributes
        absorption (optional ): if True (default) feature is treated as an absorption feature, 
                                if False it is treated as an emission feature
    Outputs:
        ax: subplot or axis object
        x_cont: the x values of the two points selected by the user
        flux_norm: the flux normalized to the continuum
    TODO: break out interactive plotting into a separate function and have this function
    take x_cont and y_cont and fit a line. This will enable it to be able to read from an
    output file or batch mode
    '''
    #Select the edges of the feature and fit a line
    input('Zoom in on line, press enter to continue')
    print('Click on left and right continuum points')
    (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
    x_cont = np.array([x1, x2])
    sort_indx = np.argsort(x_cont)
    x_cont = x_cont[sort_indx]
    y_cont = np.array([y1, y2])
    y_cont = y_cont[sort_indx]
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    continuum_val = calc_continuum(x_cont, y_cont, spectrum.wave)
    #replot the normalized spectrum
    ax.cla()
    ax.set_ylabel('Flux')
    flux_sub = spectrum.flux - continuum_val
    ax.plot(spectrum.wave, flux_sub)
    ax.set_xlim(x_cont.min()-5, x_cont.max()+5)
    #Plot the "continuum" region
    y_cont_norm = y_cont-calc_continuum(x_cont, y_cont, x_cont)
    ax.plot(x_cont, y_cont_norm, marker='o', ls='-')
    plt.draw()
    return ax, x_cont, y_cont, flux_sub
    
def calc_continuum(x_cont, y_cont, wave):
    continuum_fit = np.polyfit(x_cont, y_cont, 1)
    continuum_val = np.polyval(continuum_fit,wave)
    return continuum_val

def find_best_fit(fit_type, center_list, line_wave, line_flux, absorption=True, similar_widths=True):
    line_flux = line_flux/np.abs(np.min(line_flux))
    plt.figure()
    plt.plot(line_wave, line_flux)
    profile, args = PROFILES[fit_type]
    model = models.Const1D(0.)
    model.amplitude.fixed=True
    avg_std = (line_wave[-1]-line_wave[0])/len(center_list)
    for x_center, y_center in center_list:
        y_center = y_center/np.abs(np.min(line_flux))
        #lines.append(plt.plot(x_center, y_center, 'r+'))
        if absorption is True:
            model = model - profile(y_center, x_center, avg_std) #amplitude, mean/x_0

    #Is there a more general way to do this that doesn't limit me to triplet features?
    nmodels = model.n_submodels()
    if (nmodels > 2) and similar_widths is True:
        for i in np.arange(2, nmodels):
            model.__getattribute__('{}_{}'.format(args['width_arg'], i)).tied = tie_widths(args['width_arg'])

    fitter = fitting.LevMarLSQFitter()
    fit = fitter(model, line_wave, line_flux)
    plt.plot(line_wave, fit(line_wave))
    for model in fit[1:]:
        model.amplitude = model.amplitude*np.abs(np.min(line_flux))
    
    import pdb; pdb.set_trace()
    return fit, fit_type

    
def fit_feature(spectrum, line_name, absorption=True, similar_widths=True, input_filename=None, input_append=True):
    '''
    Fit single or multiple components to an emission or absorption feature
    Inputs:
        spectrum: spectrum1d class object with wave and flux attributes
        absorption (optional ): if True (default) feature is treated as an absorption feature, 
                        if False it is treated as an emission feature
        similar_widths: if True (default) all lines fit are required to have the same width
                        (either stddev or FWHM depending on the model)
    Outputs:
        fit: an astropy compound model. The first component is always the continuum which
             is set to 1.
        Plot:
            A plot is made of the continuum normalized flux, the region fit, the compound fit,
            and the individual components of the fit
        Print to screen:
            Model used, center, FWHM, amplitude of each component
    
    Limitations:
    * Fits either emission or absorption, not both
    * If similar widths is used then all features have exactly the same width
    #TODO: separate out the interactive from the fitting so that an input file can be passed in and fit
    '''
    ax1, ax2 = plot_spectrum(spectrum)  
    ax1, x_cont, y_cont, flux_sub = continuum_normalization(ax1, spectrum, absorption=absorption)
    #-----------------------
    #Fit the feature
    #-----------------------
    redo = 'y'
    while redo is 'y':
        #Select the region to fit
        fit_wave = np.arange(x_cont.min(), x_cont.max()+1, 0.01)
        input('Zoom in on line, press enter to continue') #Add in option to redefine continuum
        print('select left and right edges of fit')
        (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
        #Select the normalized fit spectrum
        line_wave = spectrum.wave[(spectrum.wave<max(x1, x2)) & (spectrum.wave > min(x1, x2))]
        line_flux = flux_sub[(spectrum.wave<max(x1, x2)) & (spectrum.wave > min(x1, x2))]
        #Select the line centers
        print('Select the center(s) of the line(s), press enter to continue')
        center_list = plt.ginput(-1, timeout=0, show_clicks=True)
        plt.draw()
        #Select the fitting function
        fit_type = input('What shape line would you like to fit? (g=gaussian (default), l=lorentz, m=moffat) ')
    
        fit, fit_type= find_best_fit(fit_type, center_list, line_wave, line_flux)
        lines = []
        lines.append(ax1.plot(line_wave, line_flux))
        lines.append(ax1.plot(fit_wave, fit(fit_wave)))
        
        if hasattr(fit, 'amplitude_2'): #Check that more than 1 line was fit
            for ifit in fit[1:]:
                if absorption:
                    lines.append(ax1.plot(fit_wave, fit[0](fit_wave) - ifit(fit_wave), ls=':'))
                else:
                    lines.append(ax1.plot(fit_wave, fit[0](fit_wave) + ifit(fit_wave), ls=':'))
        #Plot the residuals
        ax2.axhline(0, color = 'k')
        lines.append(ax2.plot(line_wave, line_flux - fit(line_wave)))

        plt.draw()
        redo = input('Redo the fit? (y, n(default) ')
        if redo == 'y': #Why doesn't this work?
            for iline in lines:
                iline[0].set_visible(False)
            plt.draw()
    if input_filename is not None:
        write_input(input_filename, line_name, x_cont, y_cont, (x1,y1), (x2,y2), center_list, fit, append=input_append)
    min_list, pew_list = calc_output_values(spectrum, fit, fit_wave, x_cont, y_cont)
    return min_list, pew_list

def write_input(filename, line_name, x_cont, y_cont, l_edge, r_edge, center_list, fit, append=True):
    input_dict = {line_name: 
                    {'continuum': 
                        {'left':(x_cont[0], y_cont[0]), 
                         'right':(x_cont[1], y_cont[1])},
                     'feature':
                        {'l_edge':l_edge,
                         'r_edge': r_edge,
                         'center':center_list},
                    'fit': type(fit[1]).name}}
    if append is True:
            ofile = open(filename, 'a')
    else:
        ofile = open(filename, 'w')

    ofile.write(yaml.dump(input_dict))
    ofile.close()    

def tie_widths(width_param):
    '''
    If similar_widths keyword is set in fit_feature, this function is called
    to force the widths of each feature to be the same (this is to avoid fitting
    a really broad feature)
    '''
    def inner(model):
        width = model.__getattribute__('{}_1'.format(width_param))
        return width
    return inner

def calc_min_wavelength(fit_wave, fit):
    min_list = []
    for model in fit[1:]:  #exclude constant fit
        flux = model(fit_wave)
        min_wave = fit_wave[np.argmin(flux)]
        min_wave_stddev_l, min_wave_stddev_r = line_analysis_BSNIP.calc_velocity_error(fit_wave, flux, fit, continuum=None)
        min_list.append((min_wave, min_wave_stddev_l, min_wave_stddev_r))
    return min_list

def calc_pew(spectrum, fit, x_cont, y_cont):
    continuum_l = line_analysis_BSNIP.endpoint(x_cont[0], y_cont[0], None)
    continuum_r = line_analysis_BSNIP.endpoint(x_cont[1], y_cont[1], None)
    line_indx = (spectrum.wave<=continuum_r.wave) & (spectrum.wave>=continuum_l.wave)
    pew_list = []
    for model in fit[1:]:
        fit_flux = model(spectrum.wave)
        i_pew = line_analysis_BSNIP.calc_pseudo_ew(spectrum.wave, fit_flux, continuum_l, continuum_r)
        continuum = calc_continuum(x_cont, y_cont, spectrum.wave[line_indx])
        continuum_var = line_analysis_BSNIP.calc_continuum_variance(spectrum.wave[line_indx], continuum_l, continuum_r)
        delta_wave = np.median(spectrum.wave[1:]-spectrum.wave[:-1])
        pew_var = line_analysis_BSNIP.calc_pew_variance(flux[line_indx], continuum, delta_wave, spectrum.error**2[line_indx], continuum_var)
        pew_list.append((pew, pew_var))

def calc_output_values(spectrum, fit, fit_wave, x_cont, y_cont):
    min_list = calc_min_wavelength(fit_wave, fit)
    pew_list = calc_pew(spectrum, fit, x_cont, y_cont)
    return min_list, pew_list

def print_output(fit):
    '''
    Print output nicely to screen
    '''
    if isinstance(fit[1:], models.Gaussian1D) or isinstance(fit[1:], models.Lorentz1D) or isinstance(fit[1:], models.Moffat1D):
        print('Line 1')
        print('\tModel: {}'.format(fit[1].__class__.name))
        print('\tCenter: {}'.format(fit[1].parameters[1])) #Assumes mean/x_0 is the second parameter
        if fit[1].__class__.name == 'Gaussian1D':
            print('\tFWHM: {}'.format(fit[1].fwhm))
        else:
            print('\tFWHM: {}'.format(fit[1].fwhm.value))
        print('\tAmplitude: {}/1.'.format(fit[1].amplitude.value))
    else:
        for linenum, ifit in enumerate(fit[1:]):
            print('Line {}'.format(linenum+1))
            print('\tModel: {}'.format(ifit.__class__.name))
            print('\tCenter: {}'.format(ifit.parameters[1])) #Assumes mean/x_0 is the second parameter
            if ifit.__class__.name == 'Gaussian1D':
                print('\tFWHM: {}'.format(ifit.fwhm))
            else:
                print('\tFWHM: {}'.format(ifit.fwhm.value))
            print('\tAmplitude: {}/1.'.format(ifit.amplitude.value))