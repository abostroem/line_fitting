'''
TODO: Write a function to calculate the initial flux errors (to be used in the spline
weighting) by heavily smoothing the spectrum and calculating the stddev of the points
around the smoothed flux
'''
import os
from collections import namedtuple

from astropy.io import fits
from astropy.io import ascii as asc
from astropy.table import Table
from astropy.modeling import models,fitting
from astropy.convolution import convolve, Box1DKernel
from astropy.time import Time
import numpy as np
from scipy import signal, interpolate
from matplotlib import pyplot as plt
import matplotlib.collections as collections
from matplotlib.backends.backend_pdf import PdfPages

from utilities_az import spectroscopy as spec

endpoint = namedtuple('endpoint', ['wave', 'flux', 'error'])

FIG_DIR = '../figures'

def read_iraf_spectrum(filename, redshift=0.0069):
    ofile = fits.open(filename)
    flux = ofile[0].data[0,0,:]
    err = ofile[0].data[3,0,:]
    wave = spec.calc_wavelength(ofile[0].header, np.arange(len(flux))+1)
    rest_wave = spec.apply_redshift(wave, redshift)
    return(spec.spectrum1d(rest_wave, flux, error=err))
    
def smooth_signal(flux, width, poly_deg):
    smoothed_flux = signal.savgol_filter(flux, width, poly_deg)
    return smoothed_flux
    
def find_blue_edge(wave, flux, wcenter, binsize, wmin=None):
    '''
    Calculate the slope in each bin starting at wmin, until the bin changes sign, use center for blue_edge
    '''
    wcenter_indx = np.argmin(np.abs(wave-wcenter))
    ifit = np.polyfit(wave[wcenter_indx-binsize: wcenter_indx+1],flux[wcenter_indx-binsize: wcenter_indx+1], 1)
    slope_product = 1
    #plt.plot(wave, flux)
    #plt.xlim(wmin, wcenter)
    if wmin is None:
        search_indx = np.arange(binsize,wcenter_indx+1)
    else:
        min_indx = np.argmin(np.abs(wave - wmin))
        search_indx = np.arange(min_indx, wcenter_indx)
    for indx in search_indx[::-1]:
        last_slope = ifit[0]
        if indx-binsize < 0:
            break
        ifit = np.polyfit(wave[indx-binsize:indx+1], flux[indx-binsize:indx+1], 1)
        #plt.plot(wave[indx-binsize:indx+1], flux[indx-binsize:indx+1])
        #plt.plot(wave[indx-binsize:indx+1], np.polyval(ifit, wave[indx-binsize:indx+1]))
        slope_product = last_slope*ifit[0] #if this is negative then the slope has changed sign
        if slope_product < 0:
            break  
    edge_indx = indx - binsize//2
    return edge_indx, wave[edge_indx] 

    
def find_red_edge(wave, flux, wcenter, binsize, wmax = None):
    '''
    Calculate the slope in each bin starting at wmin, until the bin changes sign, use center for red_edge
    binsize is in pixels
    '''
    wcenter_indx = np.argmin(np.abs(wave-wcenter))
    #fig = plt.figure()
    #ax1 = fig.add_subplot(2,1,1)
    #ax2 = fig.add_subplot(2,1,2, sharex=ax1)
    #ax1.plot(wave, flux)
    #ax2.plot(wave, wave-wcenter)
    ifit = np.polyfit(wave[wcenter_indx:wcenter_indx+binsize+1],flux[wcenter_indx:wcenter_indx+binsize+1], 1)
    slope_product = 1
    #plt.plot(wave, flux)
    #plt.xlim(wmin, wcenter)
    #plt.axvline(wave[wcenter_indx], color='y')
    if wmax is None:
        search_indx = np.arange(wcenter_indx, len(flux))
    else:
        max_indx = np.argmin(np.abs(wave - wmax))
        search_indx = np.arange(wcenter_indx, max_indx+1)
    #plt.plot(wave[search_indx], flux[search_indx])
    for indx in search_indx:
        last_slope = ifit[0]
        ifit = np.polyfit(wave[indx:indx+binsize+1], flux[indx:indx+binsize+1], 1)
        slope_product = last_slope*ifit[0] #if this is negative then the slope has changed sign
        #plt.plot(wave[indx:indx+binsize+1], flux[indx:indx+binsize+1])
        #plt.plot(wave[indx:indx+binsize+1], np.polyval(ifit, wave[indx:indx+binsize+1]))
        slope_product = last_slope*ifit[0] #if this is negative then the slope has changed sign
        if slope_product < 0:
            break  
    edge_indx = indx + binsize//2
    return edge_indx, wave[edge_indx] 
    
def check_max(wave, flux, edge_indx, binsize, absorption=True):
    '''
    Fit a quadratic and verify that it is the correct direction
    '''
    wmin_indx = edge_indx - binsize//2
    wmax_indx = edge_indx + binsize//2
    fit = np.polyfit(wave[wmin_indx:wmax_indx+1], flux[wmin_indx:wmax_indx+1], 2)
    if fit[0]>0:
        concavity = 'up'
    if fit[0]<0:
        concavity = 'down'
    if ((absorption is True) and (concavity is 'down')) or ((absorption is False) and (concavity is 'up')):
        good_fit = True
    else:
        good_fit = False
    return good_fit
    
def calc_rmse(data, model):
    rmse = np.sqrt(np.sum((model-data)**2)/len(data))
    print('rmse calculated over {} points'.format(len(data)))
    return rmse

def find_boundary(wave, flux, wmin, wmax, binsize, visualize=False):
    wmin_indx = np.argmin(np.abs(wave-wmin))
    wmax_indx = np.argmin(np.abs(wave-wmax))
    slope_product = []
    ifit = np.polyfit(wave[wmin_indx-binsize//2:wmin_indx+binsize//2+1], flux[wmin_indx-binsize//2:wmin_indx+binsize//2+1], 1)
    search_indx = np.arange(wmin_indx, wmax_indx+1)
    if visualize:
        fig = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2, sharex=ax1)
        ax1.set_xlim(wmin, wmax)
    for indx in search_indx:
        last_slope = ifit[0]
        ifit = np.polyfit(wave[indx-int(binsize//2):indx+int(binsize//2)+1], flux[indx-int(binsize//2):indx+int(binsize//2)+1], 1)
        slope_product.append(last_slope*ifit[0]) #if this is negative then the slope has changed sign
        if visualize:
            ax1.plot(wave[indx-int(binsize//2):indx+int(binsize//2)+1], np.polyval(ifit, wave[indx-int(binsize//2):indx+int(binsize//2)+1]))
    slope_product = np.array(slope_product)
    slope_change_indx = search_indx[slope_product<0]
    if visualize:
        ax1.set_title('Slope Plot binsize={}, fit_wmin={:4.2f}, fit_wmax={:4.2f}'.format(binsize, wave[wmin_indx-binsize//2], wave[indx+int(binsize//2)]))
        ax1.plot(wave, flux)
        ax2.plot(wave[search_indx], slope_product)
        ax2.axhline(0, color='k', linestyle=':')
    if len(slope_change_indx) == 3:
        blue_edge_indx, wcenter_indx, red_edge_indx = slope_change_indx
        return blue_edge_indx, wcenter_indx, red_edge_indx
    else:
        return None, None, None

def determine_error_binsize(wave, wave_bin=100):
    '''
    We should be calculating noise over the same wavelength range
    rather than the same number of pixels as long as one wavelength
    bin includes enough pixels. Set binsize to be 100A. If there are
    fewer than 10 pixels in 100A (dispersion is greater than 10A/pix)
    then issue a warning and make binsize 10 pixels regardless of how
    many angstroms this represents
    
    wave: array of wavelengths
    wave_bin: size of bin in angstroms
    
    outputs: binsize in pixels
    
    Note: right now this calculates the dispersion for the full wavelength range.
    For a grating/grism with a large variation in dispersion, it might make sense to
    just calculate this over the feature wavelength range.
    '''
    dispersion = np.median(wave[1:]-wave[:-1])
    binsize = np.ceil(wave_bin/dispersion)
    if binsize < 10:
        print('WARNING: dispersion = {}, \
        leading to binsize < 10 for {}$\AA$ bins, \
        setting binsize=10, making wave_bin={}'.format(dispersion, wave_bin, 10*dispersion))
        binsize=10
    return binsize

def define_continuum(wave, flux, edge_indx, binsize, err_binsize, absorption=True, visualize=False):
    '''
    Fit a quadratic and verify that it is the correct direction
    '''
    #Silverman says: "Once these two endpoints are determined, a quadratic function is 
    # fit to the data in wavelength bins centred on each endpoint." 
    #Let's start with fitting over 2 wavelength bins?
    wmin_indx = edge_indx - int(np.floor(1*binsize))
    wmax_indx = edge_indx + int(np.ceil(1*binsize))
    quad_model = models.Polynomial1D(degree=2)
    fitter = fitting.LinearLSQFitter()
    fit = fitter(quad_model, wave[wmin_indx:wmax_indx+1], flux[wmin_indx:wmax_indx+1])
    fit_extreme_wl = -fit.c1.value/(2*fit.c2.value)

    #calc rmse over edge_indx +/- 20 pixels
    wmin_rmse = edge_indx - int(err_binsize//2)
    wmax_rmse = edge_indx + int(err_binsize//2)
    rmse = calc_rmse(flux[wmin_rmse:wmax_rmse], fit(wave[wmin_rmse:wmax_rmse]))
    
    if fit.c2.value>0:
        concavity = 'up'
    if fit.c2.value<0:
        concavity = 'down'
    
    if ((absorption is True) and (concavity is 'down')) or ((absorption is False) and (concavity is 'up')):
        good_fit = True
    else:
        good_fit = False
    if visualize:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(wave, flux)
        ax.set_xlim(wave[edge_indx-4*binsize], wave[edge_indx+4*binsize])
        ax.plot(wave[wmin_rmse:wmax_rmse], fit(wave[wmin_rmse:wmax_rmse]), label='RMSE range')
        ax.plot(wave[wmin_indx:wmax_indx+1], fit(wave[wmin_indx:wmax_indx+1]), label='fit range')
        ax.axvline(wave[edge_indx], label='input continuum', color='k')
        ax.errorbar(fit_extreme_wl, fit(fit_extreme_wl), yerr=rmse, fmt='.', label='Edge w/error', zorder=10 )
        ax.set_ylim(0.9*np.min(flux[wmin_rmse:wmax_rmse]), 1.1*np.max(flux[wmin_rmse:wmax_rmse]))
        ax.legend(loc='best')
        ax.set_title('absorption={}, concavity={}, good_fit={}'.format(absorption, concavity, good_fit))
    endpt = endpoint(fit_extreme_wl, fit(fit_extreme_wl), rmse)
    return good_fit, endpt
    
def calc_pseudo_ew(wave, flux, continuum_l, continuum_r, absorption=True, visualize=False):
    '''
    wave: array
        array of wavelength (can be whole spectrum)
    flux: array
        array of fluxes (can be whole spectrum)
    
    * Create a fit to the continuum and define the continuum for each wavelength in wave
    * Use continuum wavelengths to define index location of feature
    * Calc pseudo equivalent width using flux, continuum, and delta wave as calculated from the
    wave array
    '''
    fitter = fitting.LinearLSQFitter()
    lin_mod = models.Linear1D()
    continuum_fit = fitter(lin_mod,[continuum_l.wave, continuum_r.wave], [continuum_l.flux, continuum_r.flux])
    line_indx = np.int_(np.arange(len(wave))[(wave>=continuum_l.wave)&(wave<=continuum_r.wave)])
    continuum = continuum_fit(wave[line_indx])
    delta_lambda = wave[line_indx]-wave[line_indx-1]
    if absorption is True:
        pew = np.sum(delta_lambda*(continuum - flux[line_indx])/continuum)
    else:
        pew = np.sum(delta_lambda*(flux[line_indx]-continuum)/continuum) #Check that this is true
    if visualize is True:
        fig = plt.figure()
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
        ax2.axhline(1, color='k')
        ax1.plot(wave, flux)
        ax1.plot(wave[line_indx], flux[line_indx], label='data')
        ax1.plot(wave[line_indx], continuum, label='continuum')
        ax1.set_xlim(continuum_l.wave-10, continuum_r.wave+10)
        if absorption is True:
            ax2.plot(wave[line_indx], (continuum - flux[line_indx])/continuum, label='sum for pEW')
        else:
            ax2.plot(wave[line_indx], (flux[line_indx]-continuum)/continuum, label='sum for pEW')
    return pew
    
    
def calc_continuum(wave, continuum_l, continuum_r):
    fitter = fitting.LinearLSQFitter()
    lin_mod = models.Linear1D()
    continuum_fit = fitter(lin_mod,[continuum_l.wave, continuum_r.wave], [continuum_l.flux, continuum_r.flux])
    continuum = continuum_fit(wave)
    return continuum

def find_velocity(wave, flux, error, wcenter, continuum_l, continuum_r, binsize, visualize=False):
    line_indx = np.int_(np.arange(len(wave))[(wave>=continuum_l.wave)&(wave<=continuum_r.wave)])
    windx_min = int(line_indx[0]-binsize//2)
    windx_max = int(line_indx[-1]+binsize//2)
    fitter = fitting.LinearLSQFitter()
    lin_mod = models.Linear1D()
    continuum_fit = fitter(lin_mod,[continuum_l.wave, continuum_r.wave], [continuum_l.flux, continuum_r.flux])
    continuum = continuum_fit(wave[windx_min:windx_max])
    weight = 1./error[windx_min:windx_max]
    fit = interpolate.UnivariateSpline(wave[windx_min:windx_max], flux[windx_min:windx_max]-continuum, w=weight)
    if visualize is True:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(wave[windx_min:windx_max], flux[windx_min:windx_max]-continuum)
        s = np.sum((weight * ((flux[windx_min:windx_max]-continuum)-fit(wave[windx_min:windx_max])))**2)
        ax.errorbar(wave[windx_min:windx_max], flux[windx_min:windx_max]-continuum, error[windx_min:windx_max], fmt='.', label='spectrum', zorder=1, color='b')
        ax.plot(wave[windx_min:windx_max], flux[windx_min:windx_max]-continuum, zorder=2, color='r')
        ax.plot(wave[windx_min:windx_max], fit(wave[windx_min:windx_max]), label='fit, s={:2.2f}, len(w)={:2.2f}, med(std)={:2.2e}'.format(s, len(weight), np.median(error[windx_min:windx_max])), color='gold', zorder=3)
        min_wave = wave[line_indx][np.argmin(fit(wave[line_indx]))]
        ax.axvline(min_wave)
        knots = fit.get_knots()
        ax.vlines(knots, ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1], linestyle=':')
        ax.legend(loc='best')
    return fit
    
def calc_flux_variance(data, model, err_binsize):
    kernel = Box1DKernel(err_binsize)
    errors = convolve((data-model)**2, kernel, boundary=None)
    errors = errors
    errors = np.trim_zeros(errors)
    return errors

def calc_continuum_variance(wave, continuum_l, continuum_r):
    var = (1./(continuum_l.wave - continuum_r.wave))**2 * \
          ((wave - continuum_r.wave)**2 * continuum_l.error**2 + 
           (wave - continuum_l.wave)**2 * continuum_r.error**2)
    return var

def calc_pew_variance(flux, continuum, delta_wave, flux_var, continuum_var, visualize=False, wave=None):
    '''
    Calculate the variance of the equivalent width
    Parameters:
    -----------
        flux: array
            flux values over which equivalent width is calculated
        continuum: array
            continuum values over which equivalent width is calculated
        delta_wave: int
            the wavelength bin size (in angstroms) used in the equivalent width calculation
        flux_var: array
            variance in the flux
        continuum_var: array
            variance in the continuum
    Output:
        variance in the equivalent width
    '''
    pew_var_indiv = ((flux/(continuum**2)*delta_wave)**2 * continuum_var) + \
                    ((delta_wave/continuum)**2*flux_var)
    pew_var = np.sum(pew_var_indiv)
    if visualize is True:
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        if wave is not None:
            ax.errorbar(wave, (continuum-flux)/continuum, np.sqrt(pew_err))
        else:
            print('wavelength not defined')
            wave = np.arange(len(flux))
            ax.errorbar(wave,(continuum-flux)/continuum, np.sqrt(pew_err))
    return pew_var

def calc_velocity_error(wave, flux, vel_fit, continuum = None, visualize=False):
    '''
    Calculate 1-sigma errors
    '''
    min_indx = np.argmin(vel_fit(wave))
    indx = np.argsort(wave - wave[min_indx])
    min_indx_indx = int(indx[indx==min_indx])
    cum_sum_right = 0
    if continuum is None:
        flux_sub = flux
    else:
        flux_sub = flux - continuum
    total = np.sum((flux_sub))
    for i in indx[min_indx_indx:]:
        cum_sum_right += (flux_sub)[i]
        if cum_sum_right/total > .341:
            break
    right_err = wave[i]-wave[min_indx]
    
    cum_sum_left = 0
    j=0
    for j in indx[:min_indx_indx][::-1]:
        cum_sum_left += (flux_sub)[j]
        if cum_sum_left/total > .341:
            break
    left_err = wave[min_indx]-wave[j]
    if visualize is True:
        from visualization import make_color_wheel
        colors = make_color_wheel(wave)
        plt.figure()
        for c, ind in zip(colors, indx):
            plt.plot(wave[ind], flux[ind], marker='o', ls='none', color=c)
        for k in indx[min_indx_indx:i]:
            plt.plot(wave[k], flux[k], marker='s', ls='none', color=colors[k])
        for k in indx[j:min_indx_indx]:
            plt.plot(wave[k], flux[k], marker='s', ls='none', color=colors[k])
        plt.axvline(wave[i], label='1 $\sigma$ right error')
        plt.axvline(wave[j], label='1 $\sigma$ left error')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.legend(loc='best')
    return left_err, right_err

def find_edges(spectrum, feature_dict, smooth_flux, filename, vis=False, pp=None):
    '''
    Fit incremental slopes to find where spectrum turns over
    Fit quadratic to the turnover points to get spectrum edges
    '''
    #Estimate the edges and center of the feature
    blue_edge_indx = None
    red_edge_indx = None
    good_fit_blue = False
    good_fit_red = False
    
    adjust_binsize = feature_dict['edge_param']['binsize']
    #increase the binsize until only 3 turning points are found (edges and center)
    wmin, wmax = find_wavelength_range(feature_dict, filename)
    npts_feature = len(spectrum.wave[(spectrum.wave>=wmin) & (spectrum.wave <= wmax)])
    while ((blue_edge_indx is None) or (red_edge_indx is None) or 
            (good_fit_blue is False) or (good_fit_red is False)) and \
          (adjust_binsize < 0.4*npts_feature) and \
          (adjust_binsize < feature_dict['edge_param']['binmax']):  #TODO figure out a cutoff for this
        if plt.get_fignums() is not False:
            for ifig in plt.get_fignums():
                plt.close(ifig)
        adjust_binsize += 2
        blue_edge_indx, wcenter_indx, red_edge_indx = find_boundary(spectrum.wave, 
                                                                    smooth_flux, 
                                                                    wmin, 
                                                                    wmax, 
                                                                    adjust_binsize,
                                                                    visualize=vis)
        if (blue_edge_indx is not None) and (red_edge_indx is not None):
            err_binsize = determine_error_binsize(spectrum.wave, wave_bin=100)
            #Find the feature edges and errors
            good_fit_blue, continuum_l = define_continuum(spectrum.wave, smooth_flux, blue_edge_indx, feature_dict['edge_param']['concavity_binsize'], err_binsize, absorption=True, visualize=vis)
            good_fit_red, continuum_r = define_continuum(spectrum.wave, smooth_flux, red_edge_indx, feature_dict['edge_param']['concavity_binsize'], err_binsize, absorption=True, visualize=vis)
            if continuum_l.wave > continuum_r.wave:
                print('**** WARNING: {}, left edge {} is greater than right edge {}****'.format(os.path.basename(filename), continuum_l.wave, continuum_r.wave))
                return None

    if adjust_binsize > 0.4*npts_feature:
        blue_edge_indx, wcenter_indx, red_edge_indx = find_boundary(spectrum.wave, 
                                                                    smooth_flux, 
                                                                    wmin, 
                                                                    wmax, 
                                                                    feature_dict['edge_param']['binsize'],
                                                                    visualize=True)
        print('Unable to find edges for {}, {}'.format(feature_dict['name'],os.path.basename(filename)))
        import pdb; pdb.set_trace()
        return None
    else:
        print('filename = ',os.path.basename(filename))
        print('\tinput binsize = ', feature_dict['edge_param']['binsize'])
        print('\tadjusted binsize = ', adjust_binsize)
        print('good_fit_blue={}, good_fit_red={}, combine={}'.format(good_fit_blue, good_fit_red, ((good_fit_blue is False) or (good_fit_red is False))))
        if vis is True:
            fig1 = plt.figure(1)
            pp.savefig(fig1)
            plt.close(fig1)
            fig2 = plt.figure(2)
            pp.savefig(fig2)
            plt.close(fig2)
            fig3 = plt.figure(3)
            pp.savefig(fig3)
            plt.close(fig3)
    return err_binsize, blue_edge_indx, red_edge_indx, wcenter_indx, continuum_l, continuum_r

def find_wavelength_range(feature_dict, filename):
    date = Time(fits.getval(filename, 'date-obs', 0))
    phase = date.jd - feature_dict['texpl']
    delta_wave = feature_dict['slope']*phase
    wmin = feature_dict['wmin']+delta_wave
    wmax = feature_dict['wmax']+delta_wave
    return wmin, wmax
    
    

def final_plot(wave, flux, flux_err, continuum_l, continuum_r, vel_fit, vel_min, vel_err, pew, pew_err):
    continuum = calc_continuum(wave, continuum_l, continuum_r)
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    #ax.plot(wave, flux/continuum, label='Spectrum')
    ax.errorbar(wave, flux/continuum, flux_err/continuum, fmt='.', label='Spectrum')
    ax.errorbar(np.array([continuum_l.wave, continuum_r.wave]),
                np.array([continuum_l.flux, continuum_r.flux])/calc_continuum(np.array([continuum_l.wave, continuum_r.wave]), continuum_l, continuum_r),
                np.array([continuum_l.error, continuum_r.error])/calc_continuum(np.array([continuum_l.wave, continuum_r.wave]), continuum_l, continuum_r),
                fmt='o', label='Feature Edges')
    min_continuum = calc_continuum(np.array([vel_min]), continuum_l, continuum_r)[0]
    ax.plot(wave, (vel_fit(wave)+continuum)/continuum, label='Velocity fit')
    ax.errorbar(vel_min, np.min((vel_fit(vel_min)+min_continuum)/min_continuum), xerr=(vel_err,), fmt='s', label='Velocity')
    pew_collection_err = collections.BrokenBarHCollection.span_where(wave, ymin=0, ymax=1, 
                                                        where= (wave>=vel_min-pew/2-np.sqrt(pew_err))&(wave<=vel_min+pew/2+np.sqrt(pew_err)),
                                                        color='k', alpha=0.1)
    pew_collection = collections.BrokenBarHCollection.span_where(wave, ymin=0, ymax=1, 
                                                        where= (wave>=vel_min-pew/2)&(wave<=vel_min+pew/2),
                                                        color='k', alpha=0.1, label = 'pEW')
    ax.add_collection(pew_collection_err)
    ax.add_collection(pew_collection)
    ax.legend(loc='best')
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Continuum subtracted flux')
    #ax.set_ylim(-0.05, 1.1)
    plt.ion()
    return fig
                

def characterize_line(feature_dict, filename, visualization_level=0):
    final_vis = False
    intermediate_vis = False
    pp = None
    if (visualization_level == 1) or (visualization_level ==2):
        final_vis = True
    if visualization_level == 2:
        intermediate_vis = True
        pp = PdfPages(os.path.join(FIG_DIR, 
                                 'line_fit_intermed_{}_{}.pdf'.format(feature_dict['name'], 
                                                                      os.path.basename(filename).split('.pdf')[0])))
    #Read in spectrum
    spectrum = read_iraf_spectrum(filename)
    #Remove CR and other large deviations
    smooth_flux = smooth_signal(spectrum.flux, 
                                feature_dict['smooth_param']['width'], 
                                feature_dict['smooth_param']['deg'])
    edge_results = find_edges(spectrum, feature_dict, smooth_flux, filename, vis=intermediate_vis, pp=pp)
    if edge_results is not None:
        err_binsize, blue_edge_indx, red_edge_indx, wcenter_indx, continuum_l, continuum_r = edge_results
        #Calculate the pseudo equivalent widths
        pew = calc_pseudo_ew(spectrum.wave, smooth_flux, continuum_l, continuum_r, visualize=intermediate_vis)
        if intermediate_vis is True:
            pp.savefig()
            plt.close()
        #Calculate the most common velocity
        wcenter = spectrum.wave[wcenter_indx]
        vel_fit = find_velocity(spectrum.wave, smooth_flux, spectrum.error, wcenter, continuum_l, continuum_r, err_binsize, visualize=intermediate_vis)
        if intermediate_vis is True:
            pp.savefig()
            plt.close()
        #Find the error in the pseudo equivalent width
        line_indx = np.arange(len(spectrum.wave))[(spectrum.wave>=continuum_l.wave)&(spectrum.wave<=continuum_r.wave)]
        min_indx = int(np.floor(line_indx[0]-err_binsize/2))
        max_indx = int(np.ceil(line_indx[-1]+err_binsize/2+1))
        continuum_extended = calc_continuum(spectrum.wave[min_indx:max_indx], continuum_l, continuum_r)
        flux_var = calc_flux_variance(spectrum.flux[min_indx:max_indx]-continuum_extended,
                                    vel_fit(spectrum.wave[min_indx:max_indx]), err_binsize) #These don't include the errors from the continuum subtraction yet; ok for EW calc
        if len(flux_var) > len(spectrum.flux[line_indx]):
            flux_var = flux_var[1:-1]
        continuum = calc_continuum(spectrum.wave[line_indx], continuum_l, continuum_r)
        continuum_var = calc_continuum_variance(spectrum.wave[line_indx], continuum_l, continuum_r)
        delta_wave = np.median(spectrum.wave[1:]-spectrum.wave[:-1])
        pew_var = calc_pew_variance(spectrum.flux[line_indx], continuum, delta_wave, flux_var, continuum_var, wave=spectrum.wave[line_indx])
        #Find the velocity error
        vel_err = calc_velocity_error(spectrum.wave[line_indx], spectrum.flux[line_indx], vel_fit, continuum=continuum)
        vel_min = spectrum.wave[line_indx][np.argmin(vel_fit(spectrum.wave[line_indx]))]
        if final_vis is True:
            fig = final_plot(spectrum.wave[min_indx:max_indx], spectrum.flux[min_indx:max_indx],spectrum.error[min_indx:max_indx], continuum_l, continuum_r, vel_fit, vel_min, vel_err, pew, pew_var)
            fig.suptitle(os.path.basename(filename))
        if intermediate_vis:
            pp.close()
        return pew, pew_var, vel_min, vel_err, fig