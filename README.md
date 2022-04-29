# Effective-screening-for-safe-opening-of-universities-under-Omicron-and-Delta-variants-of-COVID-19
For clarity, different codes were used to generate different files.

- To generate data for Figs 1, 2, 3(a,b) and 4(a,b), please run: customized.cpp, customized.h and customized_constants.h 
- To generate data for Figs 3(c,d) and 4(c,d), please run: full_customization.cpp, full_customization.h and full_customization_constants.h
Note that to move from 50% to 95%, you can edit this command: "double OMEGA_array[] = {0.5};" and replace 0.5 by 0.95 (and vice versa). Each case should be run by   itself

- To generate data for the supplementary tables, please run: omicron_code_updated.cpp, omicron_code_updated.h and omicron_code_updated_constants.h
Note that to move from a screening compliance of 75% to 100%, you can edit this command: "ETA_array[] =  {0.75};" and replace 0.75 by 1.

- To generate data for Fig 5, please run: omicron_varying_boosted.cpp, omicron_varying_boosted.h and omicron_varying_boosted_constants.h
