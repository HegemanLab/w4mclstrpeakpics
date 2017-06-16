# w4mclstrpeakpics

Visualize likelihoods and intensities of features among members of a cluster of W4M samples.

This is a Galaxy-oriented R package that can be used to visualize likelihoods and intensities of features among members of a cluster of samples (in data matrix and metadata files preprocessed with "XCMS", Smith *et al.*, 2006) using Workflow4Metabolomics ("W4M", Giacomoni *et al.*, 2014) data matrix and metadata by sample-class.

## vapordocs 

**TODO** For a brief introduction to this package, please see the vignette: (https://github.com/HegemanLab/w4mclstrpeakpics/blob/master/vignettes/w4mclstrpeakpics.Rmd).

**TODO** This package has been "wrapped" as a Galaxy tool
  - the "wrapping project" is here: (https://github.com/HegemanLab/w4mclstrpeakpics_galaxy_wrapper)
  - the "w4mclstrpeakpics" tool is in the toolshed at (https://toolshed.g2.bx.psu.edu/repository?repository_id=TBD)

## License

MIT License

Copyright (c) 2017 Arthur C Eschenlauer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## NEWS

### CHANGES IN VERSION 0.98.1
#### NEW FEATURES

* First release - Wrap the w4mclstrpeakpics R package that visualizes features for a group of samples.

#### INTERNAL MODIFICATIONS

* none

#### BUILD NOTES

* This package was built with R 3.3.1 under miniconda as follows:

```
  $ ~/miniconda2/bin/conda create -n r3.3.1 r-base=3.3.1
  $ source ~/miniconda2/bin/activate r3.3.1
  (r3.3.1) $ R
  > library(devtools)
  > document()
  > test()
  > check()
  > build()
  > install()

```


## Citations

Smith, Colin A. and Want, Elizabeth J. and O’Maille, Grace and Abagyan, Ruben and Siuzdak, Gary (2006). XCMS: Processing Mass Spectrometry Data for Metabolite Profiling Using Nonlinear Peak Alignment, Matching, and Identification. In Analytical Chemistry, 78 (3), pp. 779–787. [doi:10.1021/ac051437y](http://dx.doi.org/10.1021/ac051437y)


Giacomoni, F. and Le Corguille, G. and Monsoor, M. and Landi, M. and Pericard, P. and Petera, M. and Duperier, C. and Tremblay-Franco, M. and Martin, J.-F. and Jacob, D. and et al. (2014). Workflow4Metabolomics: a collaborative research infrastructure for computational metabolomics. In Bioinformatics, 31 (9), pp. 1493–1495. [doi:10.1093/bioinformatics/btu813](http://dx.doi.org/10.1093/bioinformatics/btu813)

