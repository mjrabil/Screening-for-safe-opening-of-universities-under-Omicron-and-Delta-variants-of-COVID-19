For clarity, different codes were used to generate different files.

To generate data for Figs 1, 2, 3(a,b) and 4(a,b), please run: customized.cpp, customized.h and customized_constants.h 

To generate data for Figs 3(c,d) and 4(c,d), please run: full_customization.cpp, full_customization.h and full_customization_constants.h 

To generate data for the supplementary tables, please run: omicron_code_updated.cpp, omicron_code_updated.h and omicron_code_updated_constants.h 

To generate data for Fig 5, please run: omicron_varying_boosted.cpp, omicron_varying_boosted.h and omicron_varying_boosted_constants.h

Sensitivity Analyses:
To move from 50% to 95%, please edit the command: "double OMEGA_array[] = {0.5};" by replacing 0.5 by 0.95 (and vice versa). 
To move from a screening compliance of 75% to 100%, please edit the command: "ETA_array[] = {0.75};" and replace 0.75 by 1.
To move from a vaccine coverage of 64% boosted, 18% vac, 18% unvac to a vaccine coverage of 38% boosted, 44% vac, 18% unvac, please edit the commands: double initial_U_W_S_array[] = { 4000 }; and double initial_U_W_F_array[] = { 267 }; by replacing 4000 and 267 by 10000 and 667, respectively.
To move from the base-case transmission scenario to the best- or worst-case transmission scenarios, please edit the commands double R0_S_Delta_array[] = { 6 }; and  double R0_F_Delta_array[] = { 3.2 }; by replacing 6 by 5 or 7 and 3.2 by 2.2, 4.2, respectively.

