# Border-effect-corrections-for-diagonal-line-based-recurrence-quantification-analysis-measures

This repository contains all correction schemes proposed and used in the article "Border 
effect corrections for diagonal line based recurrence quantification analysis measures" 
(Kraemer & Marwan Physics Letters A, 2019). There are basically two kinds of corrections: 

1. correction schemes for counting the diagonal lines in a recurrence plot 
   (the `dl_`-functions) and 
   
2. correction schemes for the recurrence plot itself in order to suppress the effect of 
   tangential motion (the `rp_`-functions)

See the respective docstrings for further information.

To get an idea of the functionality of the correction approaches you might just run the
MATLAB-Live scripts `examples_correction_schemes_Logistic_map.mlx` and 
`examples_correction_schemes_Roessler_system.mlx`, respectively (or the plain code found in
the corresponding `.m`-files. In these scripts we exemplary show the application of the
mentioned methods to map and flow data, similar to the data used in the paper. For a 
better understanding you should run the different sections in the script individually. 
Play around with the parameters and noise level of the data! 

In order to properly run the provided code you might need to install the CRP toolbox by
Norbert Marwan (open source; [http://tocsy.pik-potsdam.de/CRPtoolbox/]
(http://tocsy.pik-potsdam.de/CRPtoolbox/)) and the Signal Processing toolbox from MathWorks 
(just for the `rp_LM2P`-function; [https://de.mathworks.com/products/signal.html]
(https://de.mathworks.com/products/signal.html)).
