HAPRAD
======

FORTRAN program for calculation of radiative correction to semi-inclusive
hadron leptoproduction 

:Authors: I. Akushevich,
          A. Ilyichev and
          M. Osipenko

:Version: 2.0

==========================   ====================================================
``rcdat.f``                  the main program
``rcdat_test.f``             the main program for test
``fhaprad.f``                the program for parameter initialization
``fhaprad_test.f``           the program for test parameter initialization
``ihaprad.f``                the basic program of calculation
``exclusive_model.f``        the program for exclusive cross sections calculation 
``semi_inclusive_model.f``   the program for semi-inclusive structure functions
                             calculation
``h3.f``                     the program for semi-inclusive H3 structure
                             function calculation
``h4.f``                     the program for semi-inclusive H4 structure
                             function calculation  
``haprad_utils.f``           includes some subroutines and function that are
                             necessary for calculation   
``init_pdf.f``               parton distribution initialization
``pkhff.f``                  LO and NLO fragmentation functions for charged
                             pions, kaons and the inclusive sum of charged
                             hadrons 
``res.dat``                  output file
``test.dat``                 test output file
==========================   ====================================================
 
For radiative  correction calculation run::

    $ make
    $ ./haprad20.exe

Test-run reproduce relative exclusive radiative tail contribution at different
phih to the observed cross section for the kinematical points of ref [Avakian04]_.

For test calculation run::

    $ make test
    $ ./test.exe

.. [Avakian04] H. Avakian *et al.* (CLAS Collaboration), Measurement of
               beam-spin asymmetries for pi plus electroproduction above the
               baryon resonance region, *Phys. Rev. D 69* 112004 (2004).
