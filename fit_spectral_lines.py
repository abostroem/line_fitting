'''
Fit single or multi-component emission or absorption lines.
Example call:
fit = fit_spectral_lines.fit_feature(spec, similar_widths=True)
Where spec is an object with attributes wave and flux. This can be 
created with the spectrum1d class
'''
import sys
import os
from collections import namedtuple

import numpy as np
import yaml

from astropy.modeling import models, fitting
from astropy.io import fits
from astropy.table import Table

from matplotlib import pyplot as plt
import matplotlib.collections as collections



import line_analysis_BSNIP

PROFILES = {'g': (models.Gaussian1D, {'size_arg':'amplitude', 'mean_arg': 'mean', 'width_arg': 'stddev'}),
            'l': (models.Lorentz1D, {'size_arg':'amplitude', 'mean_arg': 'x_0', 'width_arg': 'fwhm'}),
            'm': (models.Moffat1D, {'size_arg':'amplitude', 'mean_arg': 'x_0', 'width_arg': 'gamma'}),
            }
            
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
    

class NoAliasDumper(yaml.Dumper):
    def ignore_aliases(self, data):
        return True

#Like a light weight class with attributes wavelength and flux
class spectrum1d(object):
    def __init__(self, wave, flux, error=None):
        self.wave = wave
        self.flux = flux
        self.error = error


class spectral_feature(object):
    def __init__(self, spectrum, spec_filename,
                 interactive=True):
        spectrum.filename = spec_filename
        self.spectrum = spectrum
        self.interactive = interactive

    def define_feature(self, name,
                       absorption=True,
                       similar_widths=True,
                       fixed_offset=False, offsets=None,
                       input_filename=None,
                       return_fit=False,
                       search_range=None,
                       define_fit_range=True,
                       output_filename=None,
                       interactive=True,
                       overwrite=False):
        '''
        Fit single or multiple components to an emission or absorption feature
        Inputs:
            name: str
                name of the line you are going to fit
            absorption (optional ): bool (True)
                if True (default) feature is treated as an absorption feature, 
                False is not implemented
            similar_widths: bool (True)
                if True (default) all lines fit are required to have the same width
                            (either stddev or FWHM depending on the model)
            fixed_offset: bool (False)
                If True , a fixed offset between gaussian means is used. offset keyword must also be defined.
            offsets: array like
                list of offsets of lines from left most line (e.g. for Ca II (8498, 8542, 8662), offsets = [44, 164])
            input_filename: str
                name of file that input fit parameters will be read from when interactive=False
            interactive: bool (True)
                If True, user is asked to define fitting parameters
                If False, this information is read from input_filename file
            return_fit: bool (False)
                If True, fit object is returned
            search_range: int/None
                If set to a number, the allowed values for the mean of each gaussian are bounded
                by input mean (either from input_file or by clicking) +/- search_range
            define_fit_range: bool (True)
                If True, user defines fit range and continuum separately
                If False, user defined continuum and this is used as the fit range
            output_filename: str (input_filename)
                The name of the file to which output should be written when interactive is True
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

        Order of definitions
            self.name = name
            self.absorption = absorption
            self.similar_widths = similar_widths
            self.fixed_offset = fixed_offset
            self.offsets = None
            self.interactive = interactive
            self.return_fit = return_fit
            self.search_range = search_range
            self.define_fit_range = define_fit_range
            self.output_filename = output_filename
            plot_spectrum
                self.interactive_fig
                self.interactive_plot
                self.interactive_residual
            continuum_normalization
                self.x_cont
                self.y_cont
                self.flux_norm
            self.continuum_r
            self.continuum_l
            self.model_wave
            self.fit_edge_wave1
            self.fit_edge_flux1
            self.fit_edge_wave2
            self.fit_edge_flux2
            self.feature_indx
            define_line_centers
                self.center_list
            self.fit_type
            fit_feature
                self.profile
                self.args
                self.fit
            self.fit_edge_l  #consider building earlier
            self.fit_edge_r
            calc_min_wavelength
                self.min_list
            calc_pew
                self.pew_list
        '''
        self.name = name
        self.absorption = absorption
        self.similar_widths = similar_widths
        self.fixed_offset = fixed_offset
        self.offsets = offsets
        self.interactive = interactive
        self.return_fit = return_fit
        self.search_range = search_range
        self.define_fit_range = define_fit_range
        self.overwrite = overwrite
        if (output_filename is None) and (interactive is True):
            self.output_filename = '{}_{}.yaml'.format('.'.join(self.spectrum.filename.split('.')[:-1]), self.name)
        elif interactive is True:
            self.output_filename = output_filename
        if (input_filename is None) and (interactive is False):
            self.input_filename = '{}_{}.yaml'.format('.'.join(self.spectrum.filename.split('.')[:-1]), self.name)
        elif (input_filename is not None):
            self.input_filename = input_filename
        if self.spectrum.error is not None:
            if np.isnan(self.spectrum.error).all():
                self.spectrum.error = None
        self.plot_spectrum()  #defines interactive_fig, interactive_plot, interactive_residual
        self.continuum_normalization() #defines x_cont, y_cont, flux_norm
        self.continuum_l = self.build_continuum_object( self.x_cont[0], self.y_cont[0])
        self.continuum_r = self.build_continuum_object(self.x_cont[1], self.y_cont[1])
        #-----------------------
        #Fit the feature
        #-----------------------
        while True:
            #Select the region to fit
            self.model_wave = np.arange(self.x_cont.min(), self.x_cont.max()+1, 0.01)
            if self.define_fit_range is False:
                (self.fit_edge_wave1, self.fit_edge_flux1) = (self.continuum_l.wave, self.continuum_l.flux)
                (self.fit_edge_wave2, self.fit_edge_flux2) = (self.continuum_r.wave, self.continuum_r.flux)
            else:
                if self.interactive is True:
                    (self.fit_edge_wave1, self.fit_edge_flux1), (self.fit_edge_wave2, self.fit_edge_flux2) = self.id_fit_region_interactive()
                else:
                    (self.fit_edge_wave1, self.fit_edge_flux1), (self.fit_edge_wave2, self.fit_edge_flux2) = self.get_id_fit_region_from_file()
            #Select the normalized fit spectrum
            self.fit_indx = (self.spectrum.wave<max(self.fit_edge_wave1, self.fit_edge_wave2)) & (self.spectrum.wave > min(self.fit_edge_wave1, self.fit_edge_wave2))
            #Select the line centers
            if self.interactive is True:
                self.define_line_centers_interactive() #defines center_list
            else:
                self.get_line_centers_from_file() #defines center_list
            if (len(self.center_list) > 1) and self.fixed_offset is True:
                assert len(self.offsets) == len(np.arange(2, len(self.center_list)+1)), 'number of offsets must match number of features being fit'

            #Select the fitting function
            if self.interactive is True:
                self.fit_type = input('What shape line would you like to fit? (g=gaussian (default), l=lorentz, m=moffat) ')
            else:
                self.fit_type = self.get_fit_type_from_file()
            lines = self.fit_feature() #defines: fit, args, profile
            if self.interactive is True:
                redo = input('Redo the fit? (y, n(default) ')
                if redo.lower() == 'y': 
                    for iline in lines:
                        iline[0].set_visible(False)
                    plt.draw()
                else:
                    break
            else:
                break
        if self.fit_edge_wave1 < self.fit_edge_wave2:
            self.fit_edge_l = self.build_continuum_object(self.fit_edge_wave1, self.fit_edge_flux1)
            self.fit_edge_r = self.build_continuum_object(self.fit_edge_wave2, self.fit_edge_flux2)
        else:
            self.fit_edge_l = self.build_continuum_object(self.fit_edge_wave2, self.fit_edge_flux2)
            self.fit_edge_r = self.build_continuum_object(self.fit_edge_wave1, self.fit_edge_flux1)
        if self.interactive is True:    
            self.write_input()
        self.calc_min_wavelength()
        self.calc_pew()
        self.get_fwhm()
        self.final_plot()
        if return_fit is True:
            return self.min_list, self.pew_list
        else:
            plt.close(self.interactive_fig)
            return self.min_list, self.pew_list
            
    def plot_spectrum(self):
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
        ax1.plot(self.spectrum.wave, self.spectrum.flux)
        ax2.set_xlabel(r'Wavelength $\AA$')
        ax1.set_ylabel('Flux')
        ax2.set_ylabel('Residual')
        self.interactive_fig = fig
        self.interactive_plot = ax1
        self.interactive_residual = ax2
        
    def continuum_normalization(self):
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
        if self.interactive is True:
            (x1, y1), (x2, y2) = self.continuum_normalization_interactive()
        else:
            (x1, y1), (x2, y2) = self.get_continuum_from_file()
        self.x_cont = np.array([x1, x2])
        self.y_cont = np.array([y1, y2])
        xlim = self.interactive_plot.get_xlim()
        ylim = self.interactive_plot.get_ylim()
        continuum_val = self.calc_continuum(self.x_cont, self.y_cont, self.spectrum.wave)
        #replot the normalized spectrum
        #Consider removing line rather than cla
        l = self.interactive_plot.get_lines()
        l[0].remove()
        self.flux_norm = self.spectrum.flux/continuum_val
        self.interactive_plot.plot(self.spectrum.wave, self.flux_norm)
        ylim_norm = self.calc_continuum(self.x_cont, self.y_cont,xlim)
        self.interactive_plot.set_xlim(self.x_cont.min()-(0.05*np.abs(self.x_cont[0]-self.x_cont[1])), self.x_cont.max()+(0.05*np.abs(self.x_cont[0]-self.x_cont[1])))
        self.interactive_plot.set_ylim(ylim/ylim_norm)
        if self.absorption is True:
            self.interactive_plot.set_ylim(ymax=1.01)
        else:
            self.interactive_plot.set_ylim(ymin=0.99)
        #Plot the "continuum" region
        y_cont_norm = self.y_cont/self.calc_continuum(self.x_cont, self.y_cont, self.x_cont)
        self.interactive_plot.plot(self.x_cont, y_cont_norm, marker='o', ls='-')
        plt.draw()
  
    def continuum_normalization_interactive(self):
        #Select the edges of the feature and fit a line
        input('Zoom in on line, press enter to continue')
        print('Click on left and right continuum points')
        (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
        return (x1, y1), (x2, y2)
    
    def get_continuum_from_file(self):
        #TODO: disentangle filename from spectrum object
        with open(self.input_filename, 'r') as ofile:
            input_dict = yaml.load(ofile)
        continuum_l = input_dict[self.spectrum.filename]['continuum']['left'] 
        continuum_r = input_dict[self.spectrum.filename]['continuum']['right'] 
        return continuum_l, continuum_r

    def calc_continuum(self, x_endpoints, y_endpoints, wave):
        continuum_fit = np.polyfit(x_endpoints, y_endpoints, 1)
        continuum_val = np.polyval(continuum_fit, wave)
        return continuum_val

    def build_continuum_object(self, x, y):
        err_binsize = line_analysis_BSNIP.determine_error_binsize(self.spectrum.wave, wave_bin=100)
        edge_indx = np.argmin(np.abs(self.spectrum.wave - x))
        if self.spectrum.error is not None:
            edge_err = np.median(self.spectrum.error[edge_indx-int(err_binsize//2):edge_indx+int(err_binsize//2)])
        else:
            edge_err = None
        continuum = line_analysis_BSNIP.endpoint(x, y, edge_err)
        return continuum
        
    def id_fit_region_interactive(self):
        #Select the region to fit   
        input('Zoom in on line, press enter to continue') #Add in option to redefine continuum
        print('select left and right edges of fit')
        (x1,y1), (x2,y2) = plt.ginput(2, timeout=0, show_clicks=True)
        return (x1, y1), (x2, y2)

    def get_id_fit_region_from_file(self):
        with open(self.input_filename, 'r') as ofile:
            input_dict = yaml.load(ofile)
        l_edge = input_dict[self.spectrum.filename]['feature']['l_edge']
        r_edge = input_dict[self.spectrum.filename]['feature']['r_edge']
        return l_edge, r_edge
    
    def define_line_centers_interactive(self):
        print('Select the center(s) of the line(s), press enter to continue')
        self.center_list = plt.ginput(-1, timeout=0, show_clicks=True)
        plt.draw()

    def get_line_centers_from_file(self):
        with open(self.input_filename, 'r') as ofile:
            input_dict = yaml.load(ofile)
        self.center_list = input_dict[self.spectrum.filename]['feature']['center']

    def get_fit_type_from_file(self):
        with open(self.input_filename, 'r') as ofile:
            input_dict = yaml.load(ofile)
        fit_type = input_dict[self.spectrum.filename]['fit']
        if fit_type == 'Gausian1D':
            return 'g'
        elif fit_type == 'Lorentz1D':
            return 'l'
        elif fit_type == 'Moffat1D':
            return 'm'

    def fit_feature(self):
            '''
            Fit single or multiple components to an emission or absorption feature
            '''
            #Build Model
            if self.fit_type in PROFILES:
                self.profile, self.args = PROFILES[self.fit_type]
            else:
                self.profile, self.args = PROFILES['g']
            model = models.Const1D(1.)
            model.amplitude.fixed=True
            lines = []
            initial_width = np.abs(self.continuum_l.wave - self.continuum_r.wave)/len(self.center_list)/10. #Following IRAF splot definition
            for x_center, y_center in self.center_list:
                lines.append(self.interactive_plot.plot(x_center, y_center, 'r+'))
                if self.absorption is True:
                    submodel = self.profile(**{self.args['mean_arg']: x_center, self.args['size_arg']: 1.-y_center, self.args['width_arg']:initial_width})
                    submodel.amplitude.min = 0.0
                    #Bound the range of possible locations of the feature if search_range is set
                    if self.search_range is not None:
                        mean = submodel.__getattribute__(self.args['mean_arg']).value
                        submodel.bounds[self.args['mean_arg']]= (mean-self.search_range, mean+self.search_range)
                    model = model - submodel
                else:
                    submodel = self.profile(**{self.args['mean_arg']: x_center, self.args['size_arg']: 1.+y_center, self.args['width_arg']:initial_width})
                    submodel.amplitude.min = 0.0
                    #Bound the range of possible locations of the feature if search_range is set
                    if self.search_range is not None:
                        mean = submodel.__getattribute__(self.args['mean_arg']).value
                        submodel.bounds[self.args['mean_arg']]= (mean-self.search_range, mean+self.search_range)
                    model = model + submodel                    
            #For blended profiles, when similar widths is true, fix the widths to be the same for each profile
            if (len(self.center_list) > 1) and self.similar_widths is True:
                for i in np.arange(2, len(self.center_list)+1):  
                    model.__getattribute__('{}_{}'.format(self.args['width_arg'], i)).tied = tie_widths(self.args['width_arg'])
            #For features with know separation, fix the separation between features (e.g. CaII triplet)
            if (len(self.center_list) > 1) and self.fixed_offset is True:
                for ioffset, i in zip(self.offsets, np.arange(2, len(self.center_list)+1)):
                    model.__getattribute__('{}_{}'.format(self.args['mean_arg'], i)).tied = tie_offset(ioffset, self.args['mean_arg'])
                
            #Fit model to spectrum
            fitter = fitting.LevMarLSQFitter()
            self.fit = fitter(model, self.spectrum.wave[self.fit_indx], self.flux_norm[self.fit_indx], maxiter=1000)
            if fitter.fit_info:
                print('fit message',fitter.fit_info['message'])
            lines.append(self.interactive_plot.plot(self.spectrum.wave[self.fit_indx], self.flux_norm[self.fit_indx]))
            lines.append(self.interactive_plot.plot(self.model_wave, self.fit(self.model_wave)))
        
            if (len(self.center_list) > 1): #Check that more than 1 line was fit
                for ifit in self.fit[1:]:
                    if self.absorption:
                        lines.append(self.interactive_plot.plot(self.model_wave, self.fit[0](self.model_wave) - ifit(self.model_wave), ls=':'))
                    else:
                        lines.append(self.interactive_plot.plot(self.model_wave, self.fit[0](self.model_wave) + ifit(self.model_wave), ls=':'))
            #Plot the residuals
            self.interactive_residual.axhline(0, color = 'k')
            lines.append(self.interactive_residual.plot(self.spectrum.wave[self.fit_indx], self.flux_norm[self.fit_indx]- self.fit(self.spectrum.wave[self.fit_indx])))
            plt.draw() 
            return lines

    def write_input(self):
        input_dict = {self.spectrum.filename: 
                        {'continuum': 
                            {'left':[float(self.continuum_l.wave), float(self.continuum_l.flux)], 
                             'right':[float(self.continuum_r.wave), float(self.continuum_r.flux)]},
                         'feature':
                            {'l_edge':[float(self.fit_edge_l.wave), float(self.fit_edge_l.flux)],
                             'r_edge': [float(self.fit_edge_r.wave), float(self.fit_edge_r.flux)],
                             'center':[[float(tup[0]), float(tup[1])] for tup in self.center_list]},
                        'fit': type(self.fit[1]).name}}
        if os.path.exists(self.output_filename):
            if self.overwrite is True:
                ofile = open(self.output_filename, 'w')
                ofile.write(yaml.dump(input_dict, Dumper=NoAliasDumper))
                ofile.close()
            else:
                overwrite = input('{} already exists, would you like to overwrite? [n], y '.format(self.output_filename))
                if overwrite == 'y':
                    ofile = open(self.output_filename, 'w')
                    ofile.write(yaml.dump(input_dict, Dumper=NoAliasDumper))
                    ofile.close()
        else:
            ofile = open(self.output_filename, 'w')
            ofile.write(yaml.dump(input_dict, Dumper=NoAliasDumper))
            ofile.close()

    def calc_min_wavelength(self):
        min_list = []
        if self.fit.n_submodels() > 2:
            for model in self.fit[1:]:  #exclude constant fit
                flux = model(self.spectrum.wave[self.fit_indx])
                min_wave = self.spectrum.wave[self.fit_indx][np.argmax(flux)]
                min_wave_fwhm_l = model.fwhm
                min_wave_fwhm_r = model.fwhm
                min_list.append((min_wave, min_wave_fwhm_l, min_wave_fwhm_r))
        else:
            flux = self.fit[1](self.spectrum.wave[self.fit_indx])
            min_wave = self.spectrum.wave[self.fit_indx][np.argmax(flux)]
            min_wave_fwhm_l = self.fit[1].fwhm
            min_wave_fwhm_r = self.fit[1].fwhm
            min_list.append((min_wave, min_wave_fwhm_l, min_wave_fwhm_r))
        self.min_list = min_list
    
    def get_fwhm(self):
        fwhm_list = []
        if self.fit.n_submodels() > 2:
            for model in self.fit[1:]:  #exclude constant fit
                fwhm_list.append(model.fwhm)
        else:
            fwhm_list = [self.fit[1].fwhm]
        self.fwhm_list = fwhm_list

    def calc_pew(self):
        pew_list = []
        continuum = self.calc_continuum([self.continuum_l.wave, self.continuum_r.wave], 
                                        [self.continuum_l.flux, self.continuum_r.flux], 
                                        self.spectrum.wave[self.fit_indx])
        continuum_all = self.calc_continuum([self.continuum_l.wave, self.continuum_r.wave], 
                                       [self.continuum_l.flux, self.continuum_r.flux], 
                                       self.spectrum.wave)
        if self.spectrum.error is None:
            continuum_var = None
        else:
            continuum_var = line_analysis_BSNIP.calc_continuum_variance(self.spectrum.wave[self.fit_indx], self.continuum_l, self.continuum_r)
        delta_wave = np.median(self.spectrum.wave[1:]-self.spectrum.wave[:-1])
        if self.spectrum.error is None:
            pew_var = None
        else:
            pew_var = line_analysis_BSNIP.calc_pew_variance(self.spectrum.flux[self.fit_indx], continuum, delta_wave, 
                                                            (self.spectrum.error**2)[self.fit_indx], continuum_var)
        if self.fit.n_submodels() > 2:
            for model in self.fit[1:]:
                fit_flux = (self.fit[0](self.spectrum.wave) - model(self.spectrum.wave))*continuum_all
                pew = line_analysis_BSNIP.calc_pseudo_ew(self.spectrum.wave, fit_flux, self.continuum_l, self.continuum_r)
                pew_list.append((pew, pew_var))
        else:
            fit_flux = (self.fit[0](self.spectrum.wave) - self.fit[1](self.spectrum.wave))*continuum_all
            pew = line_analysis_BSNIP.calc_pseudo_ew(self.spectrum.wave, fit_flux, self.continuum_l, self.continuum_r)
            pew_list.append((pew, pew_var))    
        self.pew_list = pew_list

    def final_plot(self):
        #plt.ioff()
        self.final_fig = plt.figure()
        self.final_ax = self.final_fig.add_subplot(1,1,1)
        x_fit_edge = np.array([self.fit_edge_l.wave, self.fit_edge_r.wave])
        y_fit_edge = np.array([self.fit_edge_l.flux, self.fit_edge_r.flux])
        continuum = self.calc_continuum(self.x_cont, self.y_cont, self.spectrum.wave[self.fit_indx])
        continuum_all = self.calc_continuum(self.x_cont, self.y_cont, self.spectrum.wave)
        full_profile_indx = (self.spectrum.wave>self.x_cont[0]) & (self.spectrum.wave<self.x_cont[1])
        full_profile_wave = self.spectrum.wave[full_profile_indx]
        if self.spectrum.error is None:
            self.final_ax.plot(self.spectrum.wave[self.fit_indx], self.flux_norm[self.fit_indx], '.', label='Spectrum', linestyle='-')
            self.final_ax.plot(self.x_cont,
                        self.y_cont/self.y_cont,
                        'o', label='Continuum Edges')
            self.final_ax.plot(x_fit_edge,
                    y_fit_edge,
                    'o', label='Fit Edges')
        else:
            self.final_ax.errorbar(self.spectrum.wave[self.fit_indx], 
                        self.flux_norm[self.fit_indx], 
                        self.spectrum.error[self.fit_indx]/continuum, 
                        fmt='.', label='Spectrum', linestyle='-')
            self.final_ax.errorbar(self.x_cont,
                        self.y_cont/self.y_cont,
                        np.array([self.continuum_l.error, self.continuum_r.error])/self.y_cont,
                        fmt='o', label='Continuum Edges')
            self.final_ax.errorbar(x_fit_edge,
                    y_fit_edge,
                    np.array([self.fit_edge_l.error, self.fit_edge_r.error]),
                    fmt='o', label='Fit Edges')
        self.final_ax.plot(full_profile_wave, self.fit(full_profile_wave), label='Velocity fit')
        if self.fit.n_submodels() > 2:
            for v, p, model in zip(self.min_list, self.pew_list, self.fit[1:]):
                vel, vel_err_l, vel_err_r = v
                pew, pew_err = p
                self.plot_model(vel, (vel_err_l, vel_err_r), pew, pew_err, model)
        else:
            vel, vel_err_l, vel_err_r = self.min_list[0]
            pew, pew_err = self.pew_list[0]
            self.plot_model(vel, (vel_err_l, vel_err_r), pew, pew_err, self.fit[1])
        xlim = self.final_ax.get_xlim()
        ylim = self.final_ax.get_ylim()
        self.final_ax.plot(self.spectrum.wave, self.spectrum.flux/continuum_all, color='k', alpha=0.25, zorder=0)  
        self.final_ax.set_xlim(xlim)
        self.final_ax.set_ylim(ylim)      
        self.final_ax.legend(loc='best')
        self.final_ax.set_xlabel('Wavelength')
        self.final_ax.set_ylabel('Continuum Divided flux')
        plt.ion()
    
    def plot_model(self, vel_min, vel_err, pew, pew_err, ifit):
        xerr = (vel_err,)
        if np.shape(xerr)[0] != 2:
            xerr = np.array(xerr).T
        self.final_ax.errorbar(vel_min, np.min(1.0-ifit(vel_min)), xerr=xerr, fmt='s', label='Velocity')
        if pew_err is not None:
            pew_collection_err = collections.BrokenBarHCollection.span_where(self.spectrum.wave[self.fit_indx], ymin=0, ymax=1, 
                                                            where= (self.spectrum.wave[self.fit_indx]>=vel_min-pew/2-np.sqrt(pew_err))&(self.spectrum.wave[self.fit_indx]<=vel_min+pew/2+np.sqrt(pew_err)),
                                                            color='k', alpha=0.1)
            self.final_ax.add_collection(pew_collection_err)
    
        pew_collection = collections.BrokenBarHCollection.span_where(self.spectrum.wave[self.fit_indx], ymin=0, ymax=1, 
                                                            where= (self.spectrum.wave[self.fit_indx]>=vel_min-pew/2)&(self.spectrum.wave[self.fit_indx]<=vel_min+pew/2),
                                                            color='k', alpha=0.1, label = 'pEW')
    
        self.final_ax.add_collection(pew_collection)

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
    spectrum_filename = 'ASASSN15oz_VLT_20150921sca.asci'
    tbdata = Table.read(spectrum_filename, format='ascii', names=['wave', 'flux'])
    spectrum = spectrum1d(tbdata['wave'], tbdata['flux'])
    spec_feat = spectral_feature(spectrum, spectrum_filename)
    spec_feat.define_feature('test')
    
def test_absorption():
    spectrum_filename = 'ASASSN15oz_VLT_20150921sca.asci'
    tbdata = Table.read(spectrum_filename, format='ascii', names=['wave', 'flux'])
    spectrum = spectrum1d(tbdata['wave'], tbdata['flux'])
    spec_feat = spectral_feature(spectrum, spectrum_filename)
    spec_feat.define_feature('test_absorption')
def test_emission():
    spectrum_filename = 'ASASSN15oz_VLT_20150921sca.asci'
    tbdata = Table.read(spectrum_filename, format='ascii', names=['wave', 'flux'])
    spectrum = spectrum1d(tbdata['wave'], tbdata['flux'])
    spec_feat = spectral_feature(spectrum, spectrum_filename)
    spec_feat.define_feature('test_emission', absorption=False)

def test_no_interactive():
    spectrum_filename = 'ASASSN15oz_VLT_20150921sca.asci'
    tbdata = Table.read(spectrum_filename, format='ascii', names=['wave', 'flux'])
    spectrum = spectrum1d(tbdata['wave'], tbdata['flux'])
    spec_feat = spectral_feature(spectrum, spectrum_filename)
    spec_feat.define_feature('test_emission', interactive=False, absorption=False)