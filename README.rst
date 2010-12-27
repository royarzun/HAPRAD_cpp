Radiative Correction of Hadron Leptoproduction
==============================================

C++ program for calculation of radiative correction to semi-inclusive hadron
leptoproduction, based in the original FORTRAN code HAPRAD2_.

.. _HAPRAD2: https://github.com/usm-data-analysis/HAPRAD_cpp/
             tree/master/haprad2

:Authors: H. Hakobyan,
          R. Oyarzun and
          S. Mancilla

==========================  ==================================================
``TRadCor``                 the main class
``TGlobalConfig``           class storing global configurations
``TKinematicalVariables``   class storing the input data (the five kinematical
                            variables describing the cross section)
``TLorentzInvariants``      class storing the used Lorentz invariants
``THadronKinematics``       class storing non-invariants variables describing
                            kinematics of detected hadron
``TDelta``                  class that calculates the deltas according to the
                            Eq. (19) of [Aku99]_
``TBorn``                   class that calculates the Born cross section
``TStructFunctionArray``    class for calculation of structure functions
``TThetaMatrix``            class for calculation of theta matrix, according
                            to Eq. (14) and Appendix B of [Aku09]_
``TPODINL``                 ROOT based integrable class equivalent to
                            original FORTRAN function ``podinl``
``TRV2LN``                  ROOT based integrable class equivalent to
                            original FORTRAN function ``rv2ln``
``TQQTPhi``                 ROOT based integrable class equivalent to
                            original FORTRAN function ``qqtphi``
``THapradUtils``            several functions used for calculations, including
                            exclusive cross section calculation and
                            semi-inclusive structure functions calculation
==========================  ==================================================

See the wiki_ for more information.

.. _wiki: https://github.com/usm-data-analysis/HAPRAD_cpp/wiki

.. [Aku99] I. Akushevich, N. Shumeiko, A. Soroko,
           Eur. Phys. J. C 10 (1999) 681
.. [Aku09] I. Akushevich, A. Ilyichev, M. Osipenko,
           Physics Letters B 672 (2009) 35-44
