# SOBPtools
A couple of python scripts I wrote to produce spread out bragg peak profiles in 2D (i.e. depth-dose curves)


The individual Bragg peaks are produced by routines in the file BortfeldBraggPeak.py, according to the method found in Bortfeld's 1997 paper
in Med. Phys. 24 (12) 2024-. This is a simple reapplication of the pencil beam algorithm to proton dose calculation, and includes the 
convolution with a straggling kernel to include straggling in the dose distribution.

Spread Out Bragg Peaks (SOBPs) are produced using routines in the file SOBP.py. These follow the equations found in 
Jette & Chen, 2011 (Phys. Med. Biol. 56 N131). I have included the lookup table specified in their paper, but it is not 
currently used because I needed a flat-topped SOBP picture for a figure. The fudge factor of p=1.7 (rather than 1.77 which is
canonical) results in the correct behaviour only for the specific case in the bottom of the SOBP.py file. In future, I will be 
switching to use the lookup table.
