# BaTMAN
## Bayesian-unfolding Toolkit for Multi-foil Activation with Neutrons

Copyright (C) 2021 Davide Chiesa (University and INFN of Milano - Bicocca)  
This program comes with ABSOLUTELY NO WARRANTY

Abstract
==================

BaTMAN was developed to unfold the neutron flux from multi-foil activation data with a Bayeasian statistical approach.
This toolkit contains C++ codes to prepare the input data, the statistical model and the launch script for the JAGS (Just Another Gibbs Sampler) program, which runs Markov Chain Monte Carlo (MCMC) simulations to sample Posterior probability density functions. 
BaTMAN also provides tools for the post-processing of the JAGS output files and for the analysis of the neutron flux unfolding results.

Before using BaTMAN it is recommended to read the JAGS user manual (https://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/) and the scientific papers cited below.

If you use this software, please cite the following articles:

1. D. Chiesa, et al.
    *Bayesian statistics applied to neutron activation data for reactor flux spectrum analysis*,  
    Ann. Nucl. Energy, vol. 70, pp. 157 – 168, 2014.  
    DOI: https://doi.org/10.1016/j.anucene.2014.02.012
    
2. D. Chiesa, et al.
    *Measurement of the neutron flux at spallation sources using multi-foil activation*,   
    Nucl. Instrum. Meth. A, vol. 902, pp. 14–24, 2018  
    DOI: https://doi.org/10.1016/j.nima.2018.06.016


Prerequisites
==================

- The BaTMAN software package is meant to be used together with JAGS.   
  To download and install JAGS visit the following webpage:  
  https://mcmc-jags.sourceforge.io/

- The BaTMAN software package makes use of the ROOT Data Analysis Framework.  
  To download and install ROOT visit the following webpage:  
  https://root.cern/


Installation instructions
==================

  The simplest way to compile this package is:

  1. `cd` to the directory containing the Makefile and type `make` to compile the package.
     
     example:
     
     `cd Dowload/BaTMAN/`
     
     `make`

  2. Add the following lines to the `.bashrc` file (or `.profile` or analogous) in your home directory:
  
     `export PATH=<YOUR PATH>/BaTMAN/bin:$PATH`
     
     `export LD_LIBRARY_PATH=<YOUR PATH>/BaTMAN/lib:$LD_LIBRARY_PATH`
     
     where `<YOUR PATH>` is the absolute path to the directory containing the BaTMAN package.
     
     
Instructions for use
==================

## Preparing the input files
 
   To unfold the neutron flux with BaTMAN, you first need to prepare the following input files:
   
   1. A file (e.g. `inputData.txt`) containing the list of the activation data in the following 6 column format:
      <pre>
      CrossSection_name    SSA(Bq/g)    SSA_uncertainty(Bq/g)    Atomic_mass(g/mol)    Isotopic_Abundance    CorrectionFile
      </pre>
      
      - `CrossSection_name` is the name of the file containing the cross section of the activation reaction
      - `SSA` is the *Specific Saturation Activity* (i.e. the activation rate per unit mass) to be provided together with its absolute uncertainty
      - `Atomic_mass` is the atomic (or molecular) mass of the element (or compound) of which the activated foil is composed. 
      - `Isotopic_Abundance` is the relative proportion of the target isotope, <ins>not</ins> to be provided in percent format (e.g. use 0.5, instead of 50%)
      - The last column is used to provide files to correct for self-shielding or Cd cover effects (up to 1 MeV). These files contain the bin-by-bin ratios of the estimated activation rate in the real geometry versus the activation rate that would obtained in the absence of flux perturbations. The binning to be used for calculating the bin-by-bin ratios must be the same used for the guess spectrum. Write `NoSelfShielding` in the last column if no correction is needed.
      
      You can use `#` to comment a line.
         
   2. A guess spectrum (e.g. `guessSpectrum.txt`) in the following 3 column format:  
      <pre>
      # Energy(MeV)   Flux          Uncertainty
        1.0000e-09	0.00000E-00   0.0000
        1.1220e-09	1.12274e-08   0.1433
        ...
      </pre> 
      - First column  (E<sub>i</sub>): energy binning
      - Second column (phi<sub>i</sub>): integral flux in the range [E<sub>i-1</sub>, E<sub>i</sub>] 
      - Third column (Err<sub>i</sub>): relative uncertainty associated to the neutron flux
      
      The guess spectrum normalization has no influence on the unfolding results (only its shape is is used to calculate the group effective cross sections)
   
   3. A directory (e.g. `XS_PATH/`) containing the activation cross section files.  
      The cross sections files must be provided in the following 3 column format:  
      <pre>
      # Energy(MeV)   XS(b)         Uncertainty(b)
        1.00000e-05	4.49772e+03   0.0000
        ...
      </pre>
      
   4. A file (e.g. `EnergyGroups.txt`) containing the list of energies (in MeV units) to be used for multi-group neutron flux unfolding.  
      Example:
      <pre>
      1e-9
      0.5e-6
      0.5
      20
      </pre>
   
## Using the BaTMAN toolkit
   
   1. Create a new directory (e.g. `MyUnfolding`) to host the files produced by the BaTMAN codes
   
      `mkdir MyUnfolding`  
      `cd MyUnfolding`
   
   2. Run `JAGS_input` to prepare input files for JAGS:
   
      `JAGS_input inputData.txt guessSpectrum.txt XS_PATH/ EnergyGroups.txt SelfShielding_PATH/`
      
      You will be asked to digit the order of magnitude of the neutron flux you are attempting to unfold.  
      If you set, for example, 1e12, the neutron flux will be unfolded in units of 10<sup>12</sup> neutrons/(cm<sup>2</sup>s) and the priors for the flux groups will be set as uniform in the range [0-10<sup>14</sup>] and zero outside this range.
      
      The following files will be produced (see the JAGS user manual for more information about their format and content):
      - `Data4Jags.dat` containing the activation rate data and the group effective cross sections 
      - `Model4Jags.bug` containing the Bayesian model for neutron flux unfolding
      - `LaunchJags.cmd` which is the script file to run `jags`
      
   3. Run `jags` with the command:
   
      `jags<LaunchJags.cmd`
      
      JAGS will produce the following output files:
      - `CODAchain*.txt`: 4 files containing *Posterior* samplings performed through 4 MCMC simulations with different seeds
      - `CODAindex.txt`: a file with an index to access the `CODAchain*.txt` files 
   
   4. Run `jags2root` to process the JAGS output files to create ROOT TNtuples:
   
      `jags2root`
   
      This program must be run in the directory where the `CODAchain*.txt` files have been created and produces two files:
      - `JagsOutput.txt`: a file with the list of the mean and standard deviation values of each marginalized Posterior PDF
      - `JagsOutput.root`: a ROOT file that contains 4 `NtuChain*` TNtuples (one for each Markov Chain) and a `NtuSum` TNtuple, obtained by merging them. You can access the content of this file with:  
      `root JagsOutput.root`  
      root [1] `new TBrowser`  
      You can take a look at the marginalized *Posteriors* of each monitored variable. You expect that each Markov Chain produces the same result (within the statistical uncertainty of the MCMC simulation).
      Remember also to examine the `TracePlots` directory in the `JagsOutput.root` file, that show the sampling sequence of each variable  during the MCMC simulation. They must appear as noisy, without any trend or jump. If you notice something unusual in the Trace Plots it means that the MCMC simulation did not converge. Possible causes are too wide *Prior* ranges for some flux groups or a multi-group binning choice that does not allow to constrain all the flux groups with the activation data provided in the input.

  5. Run `PosteriorAnalysis` to get the unfolding results:
  
     `PosteriorAnalysis inputData.txt EnergyGroups.txt`
  
     You will be asked to digit the order of magnitude of the neutron flux chosen in the `JAGS_input` run. Take care of digiting the <ins>same</ins> order of magnitude, because the flux results will be normalized accordingly.
     
     This program must be run in the directory where the `JagsOutput.root` file has been created and produces two files:
     - `PosteriorAnalysis.txt`: a file with the numerical results of the unfolding
     - `PosteriorAnalysis.root`: a ROOT file with some useful plots of the marginalized *Posteriors* and with the correlation plots between the unfolding variables.
     
     
References
==================

1. D. Chiesa, et al.
    *Bayesian statistics applied to neutron activation data for reactor flux spectrum analysis*,  
    Ann. Nucl. Energy, vol. 70, pp. 157 – 168, 2014.  
    DOI: https://doi.org/10.1016/j.anucene.2014.02.012
    
2. D. Chiesa, et al.
    *Measurement of the neutron flux at spallation sources using multi-foil activation*,   
    Nucl. Instrum. Meth. A, vol. 902, pp. 14–24, 2018  
    DOI: https://doi.org/10.1016/j.nima.2018.06.016
    
3. D. Chiesa, et al. 
   *Characterization of TRIGA RC-1 neutron irradiation facilities for radiation damage testing*,  
   Eur. Phys. J. Plus, vol. 135, no. 349, 2020.  
   DOI: https://doi.org/10.1140/epjp/s13360-020-00334-7
   
4. D. Chiesa, et al. 
   *Measurements of neutron fields in a wide energy range using multi-foil activation analysis*,  
   Submitted to the IEEE TNS for possible publication.  
   https://arxiv.org/abs/2110.05073

   
