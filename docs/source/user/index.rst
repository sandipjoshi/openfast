.. _user_guide:

User Documentation
==================

old

This section contains documentation for the OpenFAST module-coupling environment and its underlying modules.
Documentation covers usage of models, underlying theory, and in some cases module verification.

We are in the process of transitioning legacy FAST v8 documentation, which can be found at https://nwtc.nrel.gov/.
Details on the transition from FAST v8 to OpenFAST may be found in :numref:`fast_to_openfast`

new

**This should be removed or rephrased**
This document is designed to guide you through some of the changes that the FAST
wind turbine multi-physics engineering tool is undergoing, until complete 
documentation is made available. FAST v8.16.00a-bjj is the latest public release 
of FAST under the new modularization framework developed at NREL. 
The architecture of FAST v8 is entirely different from FAST v7.02.00d-bjj. 
These differences are highlighted in Figure 1.

The modules of FAST (AeroDyn, HydroDyn, etc.) correspond to different physical 
domains of the coupled aero-hydro-servo-elastic solution, most of which are 
separated by spatial boundaries. Figure 2 shows the control volumes associated 
with each module for fixed-bottom offshore wind turbines. Though not shown, 
finite-element blade structural dynamics is optionally available through the 
BeamDyn module and loading from surface ice on fixed-bottom offshore wind 
turbines is optionally available through the IceFloe or IceDyn modules. For 
land-based wind turbines, the HydroDyn hydrodynamics module and ice modules 
would not be used and the SubDyn multi-member substructure structural-dynamics 
module is optional. Figure 3 shows the control volumes associated with 
each module for floating offshore wind turbines. Though not shown, 
finite-element blade structural dynamics is optionally available through the
BeamDyn module and mooring and hydrodynamics are optionally available 
through the OrcaFlexInterface module.

.. toctree::
   :maxdepth: 1

   fast_to_openfast.rst
   openfast/index.rst
   api_change.rst
   aerodyn/index.rst
   aerodyn-olaf/index.rst
   aerodyn-aeroacoustics/index.rst
   beamdyn/index.rst
   elastodyn/index.rst
   cppapi/index.rst
