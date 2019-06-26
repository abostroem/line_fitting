This repository holds code for fitting absorption spectra, specifically aimed
at supernova spectra.

* fit\_line.py wrapper of fit\_spectral\_lines.py allowing the user to fit a line defining only the input filename
* line\_analysis\_BSNIP.py follows the Silverman et al (2012) method. 
* fit\_spectral\_lines.py fits multiple gaussian profiles simultaneously

Quickstart:
-----------
```  
$ python fit_lines.py input_filename
```
input filename must be in one of the following formats:
* an ascii file with the first column wavelength and the second column flux
* an ascii file with a column with 'wave' in the name and a column with 'flux' in a format that astropy.io.ascii.read can recognize the column names
* a fits file with an array of fluxes and the wavelength solution in the header
* a fits file with a 4x1xn array where the first (of four) columns represents the flux and the wavelength solution in the header
* a fits binary table or fits rec array with a column name with 'wave' in it and a column name with 'flux' in it

There are a variety of command line options that can be accessed with python fit\_line.py --help

Advanced example of fitting Ha emission on a fits binary table with the spectrum in the 1st extension
```
$ python fit_lines.py 1999em_19991208_2451520.5_1.fits --ext 1 -p e --output 1999em_fit.csv --rest_wave 6563 -n Ha --def_fit_range
```

To script the input and output:
-------------------------------

Example wrappers of these scripts are in the examples directory. 

* example/example\_BSNIP\_fitting.py demonstrates how to setup the line\_analysis\_BSNIP.py. 
For your specific run of this file, you should modify feature_dict to be 
your feature of interest (line 22), supernova details, data\_directories, and filenames.
The feature\_dict is passed to line\_analysis\_BSNIP.characterize\_line in line 54. 
This is what does the actual fitting, lines 44-53 wrap multiple features and deal
with linking directories to the correct filenames.

* example/interactive\_line\_fitting.py demonstrates how to setup the fit\_spectral\_lines.py analysis.
Here you define the line you are fitting in line 6, read in the spectrum and create a spectrum
object with attributes wave, flux, and error in line 8, add the attribute 'filename' in line 12
and call the interactive line fitting in line 13.

Documentation for most commonly used functions in fit\_spectral\_lines.py

    define_feature(spectrum, line_name, absorption=True, 
                similar_widths=True, fixed_offset = False, offsets=None,
                input_filename=None, input_append=False, interactive=True,
                return_fit=False, search_range=None, define_fit_range=True):
    
    Fit single or multiple components to an emission or absorption feature  
    Inputs:   
        spectrum: spectrum1d object  
            spectrum1d class object with wave and flux attributes  
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
|  

    fit_feature(line_wave, line_flux, fit_wave, fit_type, center_list, ax1, ax2,   
                continuum_l, continuum_r,  
                offsets=None, fixed_offset=False, similar_widths=True, absorption=True,  
                search_range=None):
        
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
            

              
        

