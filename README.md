# wrf-gchp-pumpkin

The "Pumpkin" Chemistry Abstraction Layer for the WRF (Weather Research and Forecasting) Model.

## What is "Pumpkin"?

Pumpkin is a cleaned-up version of WRF's Chemistry component ("WRF-Chem"), built by stubbing out chemistry routines called by the WRF Model and modified to be chemistry-agnostic, so that other chemistry mechanisms can be integrated into WRF based on the original programming endpoints developed for WRF-Chem's hard coupling into WRF. This is analoguous to a pumpkin being carved so lights can shine inside, so the project is codenamed "Pumpkin".

## Why was it developed?

The Pumpkin project was created for the WRF-GCHP project, developed by Haipeng Lin at Peking University, to couple GEOS-Chem High Performance with the WRF Model. While the Pumpkin project itself is solely related to WRF, its creation was crucial to the success of the WRF-GCHP project.

## How to use?

The main chemistry bindings should be located at `chem_driver.F`; the initialization routines should be located in `chemics_init.F`. The other files are mostly stubs or remnants from the WRF project that cannot be removed for compatibility reasons.

## Contact

Please contact Haipeng Lin at `linhaipeng at pku.edu.cn`.