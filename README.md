This repository holds code for fitting absorption spectra, specifically aimed
at supernova spectra.

* line\_analysis\_BSNIP.py follows the Silverman et al (2012) method. 
* fit\_spectral\_lines.py fits multiple gaussian profiles simultaneously

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

