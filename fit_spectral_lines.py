'''
Fit single or multi-component emission or absorption lines.
Example call:
fit = fit_spectral_lines.fit_feature(spec, similar_widths=True)
Where spec is an object with attributes wave and flux. This can be 
created with the spectrum1d class
'''
import sys
from collections import namedtuple

import numpy as np
import yaml

from astropy.modeling import models, fitting
from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.collections as collections



import line_analysis_BSNIP

PROFILES = {'g': (models.Gaussian1D, {'size_arg':'amplitude', 'mean_arg': 'mean', 'width_arg': 'stddev'}),
            'l': (models.Lorentz1D, {'size_arg':'amplitude', 'mean_arg': 'x_0', 'width_arg': 'fwhm'}),
            'm': (models.Moffat1D, {'size_arg':'amplitude', 'mean_arg': 'x_0', 'width_arg': 'fwhm'}),
            }
#Like a light weight class with attributes wavelength and flux
spectrum1d = namedtuple("spectrum1d", ['wave', 'flux'])

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
    return fig, ax1, ax2

def calc_continuum(x_cont, y_cont, wave):
    continuum_fit = np.polyfit(x_cont, y_cont, 1)
    continuum_val = np.polyval(continuum_fit, wave)
    return continuum_val

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
    '''
    #Select the edges of the feature and fit a line
    input('Zoom in on line, press enter to continue')
    print('Click on left and right continuum points')
    (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
    x_cont = np.array([x1, x2])
    y_cont = np.array([y1, y2])
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    continuum_val = calc_continuum(x_cont, y_cont, spectrum.wave)
    #replot the normalized spectrum
    #Consider removing line rather than cla
    l = ax.get_lines()
    l[0].remove()
    flux_norm = spectrum.flux/continuum_val
    ax.plot(spectrum.wave, flux_norm)
    ylim_norm = calc_continuum(x_cont, y_cont, xlim)
    #Use percent padding rather than a fixed amount
    ax.set_xlim(x_cont.min()-5, x_cont.max()+5)
    ax.set_ylim(ylim/ylim_norm)
    if absorption is True:
        ax.set_ylim(ymax=1.01)
    else:
        ax.set_ylim(ymin=0.99)
    #Plot the "continuum" region
    y_cont_norm = y_cont/calc_continuum(x_cont, y_cont, x_cont)
    ax.plot(x_cont, y_cont_norm, marker='o', ls='-')
    plt.draw()
    return ax, x_cont, y_cont, flux_norm
    
def fit_feature(spectrum, line_name, binsize, absorption=True, similar_widths=True, input_filename=None, input_append=None):
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
    '''
    interactive_fig, ax1, ax2 = plot_spectrum(spectrum)  
    ax1, x_cont, y_cont, flux_norm = continuum_normalization(ax1, spectrum, absorption=absorption)
    continuum_l = build_continuum_object(spectrum, x_cont[0], y_cont[0])
    continuum_r = build_continuum_object(spectrum, x_cont[1], y_cont[1])
    #-----------------------
    #Fit the feature
    #-----------------------
    while True:
        #Select the region to fit
        fit_wave = np.arange(x_cont.min(), x_cont.max()+1, 0.01)
        input('Zoom in on line, press enter to continue') #Add in option to redefine continuum
        print('select left and right edges of fit')
        (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
        #Select the normalized fit spectrum
        line_indx = (spectrum.wave<max(x1, x2)) & (spectrum.wave > min(x1, x2))
        line_wave = spectrum.wave[line_indx]
        line_flux = flux_norm[line_indx]
        #Select the line centers
        print('Select the center(s) of the line(s), press enter to continue')
        center_list = plt.ginput(-1, timeout=0, show_clicks=True)
        plt.draw()
        #Select the fitting function
        fit_type = input('What shape line would you like to fit? (g=gaussian (default), l=lorentz, m=moffat) ')
    
        if fit_type in PROFILES:
            profile, args = PROFILES[fit_type]
        else:
            profile, args = PROFILES['g']
        model = models.Const1D(1.)
        model.amplitude.fixed=True
        lines = []

        for x_center, y_center in center_list:
            lines.append(ax1.plot(x_center, y_center, 'r+'))
            if absorption is True:
                model = model - profile(**{args['mean_arg']: x_center, args['size_arg']: 1.-y_center})
            else:
                model = model + profile(**{args['mean_arg']: x_center, args['size_arg']: 1.-y_center})
        if (len(center_list) > 1) and similar_widths is True:
            for i in np.arange(2, len(center_list)+1):  #For blended profiles, when similar widths is true, fix the widths to be the same for each profile
                model.__getattribute__('{}_{}'.format(args['width_arg'], i)).tied = tie_widths(args['width_arg'])
    
        fitter = fitting.LevMarLSQFitter()
        fit = fitter(model, line_wave, line_flux)
        
        lines.append(ax1.plot(line_wave, line_flux))
        lines.append(ax1.plot(fit_wave, fit(fit_wave)))
        
        if (len(center_list) > 1): #Check that more than 1 line was fit
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
        if redo.lower() == 'y': #Why doesn't this work?
            for iline in lines:
                iline[0].set_visible(False)
            plt.draw()
        else:
            break
    fit_edge_l = build_continuum_object(spectrum, x1, y1)
    fit_edge_r = build_continuum_object(spectrum, x2, y2)
    if input_filename is not None:
        write_input(input_filename, line_name, continuum_l, continuum_r, fit_edge_l, fit_edge_r, center_list, fit, append=input_append)
    min_list, pew_list = calc_output_values(spectrum, fit, fit_wave, continuum_l, continuum_r, binsize, args)
    fig = final_plot(line_wave, line_flux, spectrum.error[line_indx], 
                continuum_l, continuum_r, 
                fit_edge_l, fit_edge_r, 
                fit, min_list,
                pew_list)
    plt.close(interactive_fig)
    return min_list, pew_list, fig

def plot_model(ax, wave, continuum_l, continuum_r, vel_min, vel_err, pew, pew_err, fit):
    x_cont = np.array([continuum_l.wave, continuum_r.wave])
    y_cont = np.array([continuum_l.flux, continuum_r.flux])
    ax.errorbar(vel_min, np.min(1.0-fit(vel_min)), xerr=(vel_err,), fmt='s', label='Velocity')
    pew_collection_err = collections.BrokenBarHCollection.span_where(wave, ymin=0, ymax=1, 
                                                        where= (wave>=vel_min-pew/2-np.sqrt(pew_err))&(wave<=vel_min+pew/2+np.sqrt(pew_err)),
                                                        color='k', alpha=0.1)
    pew_collection = collections.BrokenBarHCollection.span_where(wave, ymin=0, ymax=1, 
                                                        where= (wave>=vel_min-pew/2)&(wave<=vel_min+pew/2),
                                                        color='k', alpha=0.1, label = 'pEW')
    ax.add_collection(pew_collection_err)
    ax.add_collection(pew_collection)
    
def final_plot(wave, flux, flux_err, continuum_l, continuum_r, fit_edge_l, fit_edge_r, fit, min_list, pew_list):
    #plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    x_cont = np.array([continuum_l.wave, continuum_r.wave])
    y_cont = np.array([continuum_l.flux, continuum_r.flux])
    x_fit_edge = np.array([fit_edge_l.wave, fit_edge_r.wave])
    y_fit_edge = np.array([fit_edge_l.flux, fit_edge_r.flux])
    continuum = calc_continuum(x_cont, y_cont, wave)
    ax.errorbar(wave, flux, flux_err/continuum, fmt='.', label='Spectrum', linestyle='-')
    ax.errorbar(x_cont,
                y_cont/y_cont,
                np.array([continuum_l.error, continuum_r.error])/y_cont,
                fmt='o', label='Continuum Edges')
    ax.errorbar(x_fit_edge,
            y_fit_edge,
            np.array([fit_edge_l.error, fit_edge_r.error]),
            fmt='o', label='Fit Edges')
    ax.plot(wave, fit(wave), label='Velocity fit')
    if fit.n_submodels() > 2:
        for v, p, model in zip(min_list, pew_list, fit[1:]):
            vel, vel_err_l, vel_err_r = v
            pew, pew_err = p
            plot_model(ax, wave, continuum_l, continuum_r, vel, (vel_err_l, vel_err_r), pew, pew_err, model)
    else:
        vel, vel_err_l, vel_err_r = min_list[0]
        pew, pew_err = pew_list[0]
        plot_model(ax, wave, continuum_l, continuum_r, vel, (vel_err_l, vel_err_r), pew, pew_err, fit[1])
        
    ax.legend(loc='best')
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Continuum Divided flux')
    plt.ion()
    return fig

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

def build_continuum_object(spectrum, x, y):
    err_binsize = line_analysis_BSNIP.determine_error_binsize(spectrum.wave, wave_bin=100)
    edge_indx = np.argmin(np.abs(spectrum.wave - x))
    edge_err = np.median(spectrum.error[edge_indx-int(err_binsize//2):edge_indx+int(err_binsize//2)])
    continuum = line_analysis_BSNIP.endpoint(x, y, edge_err)
    return continuum

def calc_min_wavelength(fit_wave, fit, args):
    min_list = []
    if fit.n_submodels() > 2:
        for model in fit[1:]:  #exclude constant fit
            flux = model(fit_wave)
            min_wave = fit_wave[np.argmax(flux)]
            #min_wave_stddev_l, min_wave_stddev_r = line_analysis_BSNIP.calc_velocity_error(fit_wave, flux, fit, continuum=None)
            min_wave_stddev_l = model.__getattribute__(args['width_arg'])
            min_wave_stddev_r = model.__getattribute__(args['width_arg'])
            min_list.append((min_wave, min_wave_stddev_l, min_wave_stddev_r))
    else:
        flux = fit[1](fit_wave)
        min_wave = fit_wave[np.argmax(flux)]
        #min_wave_stddev_l, min_wave_stddev_r = line_analysis_BSNIP.calc_velocity_error(fit_wave, flux, fit, continuum=None)
        min_wave_stddev_l = fit[1].__getattribute__(args['width_arg']).value
        min_wave_stddev_r = fit[1].__getattribute__(args['width_arg']).value
        min_list.append((min_wave, min_wave_stddev_l, min_wave_stddev_r))
    return min_list
    
def calc_pew(spectrum, fit, continuum_l, continuum_r, binsize):
    line_indx = (spectrum.wave<=continuum_r.wave) & (spectrum.wave>=continuum_l.wave)
    pew_list = []
    continuum = calc_continuum([continuum_l.wave, continuum_r.wave], [continuum_l.flux, continuum_r.flux], spectrum.wave[line_indx])
    continuum_var = line_analysis_BSNIP.calc_continuum_variance(spectrum.wave[line_indx], continuum_l, continuum_r)
    delta_wave = np.median(spectrum.wave[1:]-spectrum.wave[:-1])
    pew_var = line_analysis_BSNIP.calc_pew_variance(spectrum.flux[line_indx], continuum, delta_wave, (spectrum.error**2)[line_indx], continuum_var)
    if fit.n_submodels() > 2:
        for model in fit[1:]:
            fit_flux = model(spectrum.wave)
            pew = line_analysis_BSNIP.calc_pseudo_ew(spectrum.wave, fit_flux, continuum_l, continuum_r)
            pew_list.append((pew, pew_var))
    else:
        fit_flux = fit[1](spectrum.wave)
        pew = line_analysis_BSNIP.calc_pseudo_ew(spectrum.wave, fit_flux, continuum_l, continuum_r)
        pew_list.append((pew, pew_var))    
    return pew_list

def calc_output_values(spectrum, fit, fit_wave, continuum_l, continuum_r, binsize, args):
    min_list = calc_min_wavelength(fit_wave, fit, args)
    pew_list = calc_pew(spectrum, fit, continuum_l, continuum_r, binsize)
    return min_list, pew_list


def write_input(filename, line_name, continuum_l, continuum_r, fit_edge_l, fit_edge_r, center_list, fit, append=True):
    input_dict = {line_name: 
                    {'continuum': 
                        {'left':(continuum_l.wave, continuum_l.flux), 
                         'right':(continuum_r.wave, continuum_r.flux)},
                     'feature':
                        {'l_edge':(fit_edge_l.wave, fit_edge_l.flux),
                         'r_edge': (fit_edge_r.wave, fit_edge_r.flux),
                         'center':center_list},
                    'fit': type(fit[1]).name}}
    if append is True:
            ofile = open(filename, 'a')
    else:
        ofile = open(filename, 'w')

    ofile.write(yaml.dump(input_dict))
    ofile.close()  

def print_output(fit):
    '''
    Print output nicely to screen
    '''
    for linenum, ifit in enumerate(fit[1:]):
        print('Line {}'.format(linenum+1))
        print('\tModel: {}'.format(ifit.__class__.name))
        print('\tCenter: {}'.format(ifit.parameters[1])) #Assumes mean/x_0 is the second parameter
        if ifit.__class__.name == 'Gaussian1D':
            print('\tFWHM: {}'.format(ifit.fwhm))
        else:
            print('\tFWHM: {}'.format(ifit.fwhm.value))
        print('\tAmplitude: {}/1.'.format(ifit.amplitude.value))
        



        
if __name__ == "__main__":
    ofile = fits.open('asassn15oz_20151006_redblu_101906.800_multi.fits')
    tbdata = ofile[1].data
    spectrum = spectrum1d(tbdata['wavelength'], tbdata['flux'])
    fit_feature(spectrum)