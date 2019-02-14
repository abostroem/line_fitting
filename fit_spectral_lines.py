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
    
def get_continuum_from_file(spec_filename, input_filename):
    with open(input_filename, 'r') as ofile:
        input_dict = yaml.load(ofile)
    continuum_l = input_dict[spec_filename]['continuum']['left'] 
    continuum_r = input_dict[spec_filename]['continuum']['right'] 
    return continuum_l, continuum_r
    
def continuum_normalization_interactive():
    #Select the edges of the feature and fit a line
    input('Zoom in on line, press enter to continue')
    print('Click on left and right continuum points')
    (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
    return (x1, y1), (x2, y2)

def id_fit_region_interactive():
    #Select the region to fit   
    input('Zoom in on line, press enter to continue') #Add in option to redefine continuum
    print('select left and right edges of fit')
    (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
    return (x1, y1), (x2, y2)

def get_id_fit_region_from_file(spec_filename, input_filename):
    with open(input_filename, 'r') as ofile:
        input_dict = yaml.load(ofile)
    l_edge = input_dict[spec_filename]['feature']['l_edge']
    r_edge = input_dict[spec_filename]['feature']['r_edge']
    return l_edge, r_edge
    
def define_line_centers_interactive():
        print('Select the center(s) of the line(s), press enter to continue')
        center_list = plt.ginput(-1, timeout=0, show_clicks=True)
        plt.draw()
        return center_list

def get_line_centers_from_file(spec_filename, input_filename):
    with open(input_filename, 'r') as ofile:
        input_dict = yaml.load(ofile)
    center_list = input_dict[spec_filename]['feature']['center']
    return center_list

def get_fit_type_from_file(spec_filename, input_filename):
    with open(input_filename, 'r') as ofile:
        input_dict = yaml.load(ofile)
    fit_type = input_dict[spec_filename]['fit']
    if fit_type == 'Gausian1D':
        return 'g'
    elif fit_type == 'Lorentz1D':
        return 'l'
    elif fit_type == 'Moffat1D':
        return 'm'

def continuum_normalization(ax, spectrum, absorption=True, interactive=True, input_filename=None):
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
    if interactive is True:
        (x1, y1), (x2, y2) = continuum_normalization_interactive()
    else:
        (x1, y1), (x2, y2) = get_continuum_from_file(spectrum.filename, input_filename)
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
    
def fit_feature(line_wave, line_flux, fit_wave, fit_type, center_list, ax1, ax2, 
                continuum_l, continuum_r,
                offsets=None, fixed_offset=False, similar_widths=True, absorption=True,
                search_range=None):
        '''
        Fit single or multiple components to an emission or absorption feature
        Inputs:
            line_wave: array like
                wavelength in spectrum that corresponding to feature. Should be same length as line_flux
            line_flux array like
                flux in spectrum that corresponds to feature. Should be same length as line_wave
            fit_wave: array like
                wavelengths for the fit to be evaluated at (for plotting purposes)
            fit_type: str (g, l, m)
                type of function to fit: g=Gaussian, l=Lorentzian, m=Moffatt. Default is 'g'
            center_list: list
                list of (x, y) pair tuples corresponding to initial guesses for the position of 
                the features being fit. The number of objects in this list defines the number
                of lines fit.
            ax1: matplotlib subplot object
                subplot to plot the compound fit and its components
            ax2: matplotlib subplot object
                subplot to plot residuals of compound fit
            continuum_l: endpoint named tuple
                endpoint object containing the wavelength, flux, and error of the left continuum point 
            continuum_r: endpoint named tuple
                endpoint object containing the wavelength, flux, and error of the right continuum point
            absorption (optional ): bool
                if True (default) feature is treated as an absorption feature, 
                False is not implemented
            similar_widths: bool
                if True (default) all lines fit are required to have the same width
                            (either stddev or FWHM depending on the model)
            fixed_offset: bool
                If True, a fixed offset between gaussian means is used. offset keyword must also be defined.
            offsets: array like
                list of offsets of lines from left most line (e.g. for Ca II (8498, 8542, 8662), offsets = [44, 164])
            search_range: int/None
                If set to a number, the allowed values for the mean of each gaussian are bounded
                by input mean (either from input_file or by clicking) +/- search_range
        Outputs:
            fit: astropy.modeling fit object
                if return_fit is True, the astropy modeling fit object is returned.
            lines: list
                list of matplotlib line2d objects
            args: dict
                dictionary with keys 'size_arg', 'mean_arg', and 'width_arg' and values that correspond to the
                parameter names in the model defined in fit_type
            ax1: same axis object input - but with things plotted on it
            ax2: same axis object input - but with things plotted on it
              
        '''
        if fit_type in PROFILES:
            profile, args = PROFILES[fit_type]
        else:
            profile, args = PROFILES['g']
        model = models.Const1D(1.)
        model.amplitude.fixed=True
        lines = []
        initial_width = np.abs(continuum_l.wave - continuum_r.wave)/len(center_list)/10. #Following IRAF splot definition
        for x_center, y_center in center_list:
            lines.append(ax1.plot(x_center, y_center, 'r+'))
            if absorption is True:
                submodel = profile(**{args['mean_arg']: x_center, args['size_arg']: 1.-y_center, args['width_arg']:initial_width})
                submodel.amplitude.min = 0.0
                if search_range is not None:
                    mean = submodel.__getattribute__(args['mean_arg']).value
                    submodel.bounds[args['mean_arg']]= (mean-search_range, mean+search_range)
                model = model - submodel
            else:
                raise NotImplemented 
        if (len(center_list) > 1) and similar_widths is True:
            for i in np.arange(2, len(center_list)+1):  #For blended profiles, when similar widths is true, fix the widths to be the same for each profile
                model.__getattribute__('{}_{}'.format(args['width_arg'], i)).tied = tie_widths(args['width_arg'])
        if (len(center_list) > 1) and fixed_offset is True:
            for ioffset, i in zip(offsets, np.arange(2, len(center_list)+1)):
                model.__getattribute__('{}_{}'.format(args['mean_arg'], i)).tied = tie_offset(ioffset, args['mean_arg'])
        fitter = fitting.LevMarLSQFitter()
        fit = fitter(model, line_wave, line_flux, maxiter=1000)
        if fitter.fit_info:
            print('fit message',fitter.fit_info['message'])
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
        return fit, lines, args, ax1, ax2
   
def define_feature(spectrum, line_name, absorption=True, 
                similar_widths=True, fixed_offset = False, offsets=None,
                input_filename=None, input_append=False, interactive=True,
                return_fit=False, search_range=None, define_fit_range=True):
    '''
    Fit single or multiple components to an emission or absorption feature
    Inputs:
        spectrum: spectrum1d object
            spectrum1d class object with wave, flux, err, and filename attributes.
            filename will be used as a unique key when the input values are recorded
        line_name: str
            name of the line you are going to fit
        absorption (optional ): bool
            if True (default) feature is treated as an absorption feature, 
            False is not implemented
        similar_widths: bool
            if True (default) all lines fit are required to have the same width
                        (either stddev or FWHM depending on the model)
        fixed_offset: bool
            If True, a fixed offset between gaussian means is used. offset keyword must also be defined.
        offsets: array like
            list of offsets of lines from left most line (e.g. for Ca II (8498, 8542, 8662), offsets = [44, 164])
        input_filename: str
            name of file that input fit parameters will be read from when interactive=False
        input_append: bool
            if True, input parameters will be appended to current informatio in input_filename, 
            if False, the file will be over written
        interactive: bool
            If True, user is asked to define fitting parameters
            If False, this information is read from input_filename file
        return_fit: bool
            If True, fit object is returned
        search_range: int/None
            If set to a number, the allowed values for the mean of each gaussian are bounded
            by input mean (either from input_file or by clicking) +/- search_range
        define_fit_range: bool
            If True, user defines fit range and continuum separately
            If False, user defined continuum and this is used as the fit range
    Outputs:
        min_list: list
            list of wavelength locations of minimum of each fit
        pew_list: list
            list of tuples (equivalent width, left error, right error) for each feature fit.
            the error is defined as the value that includes 33.3% of the total integrated flux 
            defined between continuum_l.wave and continuum_r.wave
        fig: matplotlib figure object
            a figure of the spectrum, including errors, the continuum points with errors, 
            the fit bounding points, the minima wavelengths including errors, compound
            fit, and the equivalent width (and error). This object can be saved if desired
        fit: astropy.modeling fit object
            if return_fit is True, the astropy modeling fit object is returned.
        interactive_fig: matplotlib figure object
            if return_fit is True, the figure of the compound fit and each individual element are returned

    
    Limitations:
    * Fit to emission is not yet implemented
    * If similar widths is used then all features have exactly the same width
    '''
    interactive_fig, ax1, ax2 = plot_spectrum(spectrum)  
    ax1, x_cont, y_cont, flux_norm = continuum_normalization(ax1, spectrum, input_filename=input_filename, interactive=interactive, absorption=absorption)
    continuum_l = build_continuum_object(spectrum, x_cont[0], y_cont[0])
    continuum_r = build_continuum_object(spectrum, x_cont[1], y_cont[1])
    #-----------------------
    #Fit the feature
    #-----------------------
    while True:
        #Select the region to fit
        fit_wave = np.arange(x_cont.min(), x_cont.max()+1, 0.01)
        if define_fit_range is False:
            (x1,y1) = (continuum_l.wave, continuum_l.flux)
            (x2, y2) = (continuum_r.wave, continuum_r.flux)
        else:
            if interactive is True:
                (x1, y1), (x2, y2) = id_fit_region_interactive()
            else:
                (x1, y1), (x2, y2) = get_id_fit_region_from_file(spectrum.filename, input_filename)
        #Select the normalized fit spectrum
        line_indx = (spectrum.wave<max(x1, x2)) & (spectrum.wave > min(x1, x2))
        line_wave = spectrum.wave[line_indx]
        line_flux = flux_norm[line_indx]
        #Select the line centers
        if interactive is True:
            center_list = define_line_centers_interactive()
        else:
            center_list = get_line_centers_from_file(spectrum.filename, input_filename)
        #Select the fitting function
        if interactive is True:
            fit_type = input('What shape line would you like to fit? (g=gaussian (default), l=lorentz, m=moffat) ')
        else:
            fit_type = get_fit_type_from_file(spectrum.filename, input_filename)
        fit, lines, args, ax1, ax2 = fit_feature(line_wave, line_flux, fit_wave, fit_type, center_list, ax1, ax2, 
                                                 continuum_l, continuum_r,
                                                 offsets=offsets, fixed_offset=fixed_offset, similar_widths=similar_widths,
                                                 search_range=search_range)
        if interactive is True:
            redo = input('Redo the fit? (y, n(default) ')
            if redo.lower() == 'y': 
                for iline in lines:
                    iline[0].set_visible(False)
                plt.draw()
            else:
                break
        else:
            break
    fit_edge_l = build_continuum_object(spectrum, x1, y1)
    fit_edge_r = build_continuum_object(spectrum, x2, y2)
    if (input_filename is not None) and (interactive is True):
        write_input(spectrum.filename, input_filename, line_name, continuum_l, continuum_r, fit_edge_l, fit_edge_r, center_list, fit, append=input_append)
    min_list, pew_list = calc_output_values(spectrum, fit, fit_wave, continuum_l, continuum_r, args)
    fig = final_plot(spectrum, line_wave, line_flux, spectrum.error[line_indx], 
                continuum_l, continuum_r, 
                fit_edge_l, fit_edge_r, 
                fit, min_list,
                pew_list)
    if return_fit is True:
        return min_list, pew_list, fig, fit, interactive_fig
    else:
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
    
def final_plot(spectrum, wave, flux, flux_err, continuum_l, continuum_r, fit_edge_l, fit_edge_r, fit, min_list, pew_list):
    #plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    x_cont = np.array([continuum_l.wave, continuum_r.wave])
    y_cont = np.array([continuum_l.flux, continuum_r.flux])
    x_fit_edge = np.array([fit_edge_l.wave, fit_edge_r.wave])
    y_fit_edge = np.array([fit_edge_l.flux, fit_edge_r.flux])
    continuum = calc_continuum(x_cont, y_cont, wave)
    continuum_all = calc_continuum(x_cont, y_cont, spectrum.wave)
    full_profile_indx = (spectrum.wave>x_cont[0]) & (spectrum.wave<x_cont[1])
    full_profile_wave = spectrum.wave[full_profile_indx]
    ax.errorbar(wave, flux, flux_err/continuum, fmt='.', label='Spectrum', linestyle='-')
    ax.errorbar(x_cont,
                y_cont/y_cont,
                np.array([continuum_l.error, continuum_r.error])/y_cont,
                fmt='o', label='Continuum Edges')
    ax.errorbar(x_fit_edge,
            y_fit_edge,
            np.array([fit_edge_l.error, fit_edge_r.error]),
            fmt='o', label='Fit Edges')
    ax.plot(full_profile_wave, fit(full_profile_wave), label='Velocity fit')
    if fit.n_submodels() > 2:
        for v, p, model in zip(min_list, pew_list, fit[1:]):
            vel, vel_err_l, vel_err_r = v
            pew, pew_err = p
            plot_model(ax, full_profile_wave, continuum_l, continuum_r, vel, (vel_err_l, vel_err_r), pew, pew_err, model)
    else:
        vel, vel_err_l, vel_err_r = min_list[0]
        pew, pew_err = pew_list[0]
        plot_model(ax, wave, continuum_l, continuum_r, vel, (vel_err_l, vel_err_r), pew, pew_err, fit[1])
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.plot(spectrum.wave, spectrum.flux/continuum_all, color='k', alpha=0.25, zorder=0)  
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)      
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
    
def tie_offset(offset, mean_param):
    '''
    for multiple lines of the same element, fix the offset between the lines
    '''
    def inner(model):
        new_mean = model.__getattribute__('{}_1'.format(mean_param)).value+offset
        return new_mean
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
    
def calc_pew(spectrum, fit, continuum_l, continuum_r):
    line_indx = (spectrum.wave<=continuum_r.wave) & (spectrum.wave>=continuum_l.wave)
    pew_list = []
    continuum = calc_continuum([continuum_l.wave, continuum_r.wave], [continuum_l.flux, continuum_r.flux], spectrum.wave[line_indx])
    continuum_all = calc_continuum([continuum_l.wave, continuum_r.wave], [continuum_l.flux, continuum_r.flux], spectrum.wave)
    continuum_var = line_analysis_BSNIP.calc_continuum_variance(spectrum.wave[line_indx], continuum_l, continuum_r)
    delta_wave = np.median(spectrum.wave[1:]-spectrum.wave[:-1])
    pew_var = line_analysis_BSNIP.calc_pew_variance(spectrum.flux[line_indx], continuum, delta_wave, (spectrum.error**2)[line_indx], continuum_var)
    if fit.n_submodels() > 2:
        for model in fit[1:]:
            fit_flux = (fit[0](spectrum.wave) - model(spectrum.wave))*continuum_all
            pew = line_analysis_BSNIP.calc_pseudo_ew(spectrum.wave, fit_flux, continuum_l, continuum_r)
            pew_list.append((pew, pew_var))
    else:
        fit_flux = (fit[0](spectrum.wave) - fit[1](spectrum.wave))*continuum_all
        pew = line_analysis_BSNIP.calc_pseudo_ew(spectrum.wave, fit_flux, continuum_l, continuum_r)
        pew_list.append((pew, pew_var))    
    return pew_list

def calc_output_values(spectrum, fit, fit_wave, continuum_l, continuum_r, args):
    min_list = calc_min_wavelength(fit_wave, fit, args)
    pew_list = calc_pew(spectrum, fit, continuum_l, continuum_r)
    return min_list, pew_list

class NoAliasDumper(yaml.Dumper):
    def ignore_aliases(self, data):
        return True

def write_input(spec_filename, filename, line_name, continuum_l, continuum_r, fit_edge_l, fit_edge_r, center_list, fit, append=True):
    input_dict = {spec_filename: 
                    {'continuum': 
                        {'left':[float(continuum_l.wave), float(continuum_l.flux)], 
                         'right':[float(continuum_r.wave), float(continuum_r.flux)]},
                     'feature':
                        {'l_edge':[float(fit_edge_l.wave), float(fit_edge_l.flux)],
                         'r_edge': [float(fit_edge_r.wave), float(fit_edge_r.flux)],
                         'center':[[float(tup[0]), float(tup[1])] for tup in center_list]},
                    'fit': type(fit[1]).name}}
    if append is True:
            ofile = open(filename, 'a')
    else:
        ofile = open(filename, 'w')

    ofile.write(yaml.dump(input_dict, Dumper=NoAliasDumper))
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