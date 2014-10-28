BF-rvplotter
============

BF_python.py is a radial velocity curve extractor using the Broadening Function technique.
I assume you have some spectral observations of an SB2 binary star and would really like RV curves.

Based loosely on Rucinski's BFall_IDL.pro and uses the PyAstronomy tools, among other fun packages.

See more here: http://www.astro.utoronto.ca/~rucinski/BFdescription.html
and here: http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/svd.html

In practice, you need to run BF_python.py twice: once to do the initial BF, and then again
to properly fit the peaks of each BF with a Gaussian. You can run it a third time to manually
tweak which BF peaks go with which star.

Finally, rvplotmaker.py makes a pretty plot of your final SB2 RV curve.
