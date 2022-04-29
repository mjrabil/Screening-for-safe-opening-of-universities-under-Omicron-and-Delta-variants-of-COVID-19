#include <iostream>
#include <ctime>
#include <cstdlib>
#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

class region
{
public:

    double population_size;
    unsigned int ID;
    double ratio_per_100_tests;
    double ratio_of_peak_per_100_tests;


    //RATES:
    double hospitalization_rate_N_S;
    double hospitalization_rate_N_F;
    double hospitalization_rate_W_S;
    double hospitalization_rate_W_F;
    double hospitalization_rate_B_S;
    double hospitalization_rate_B_F;

    double fatality_rate_N_S;
    double fatality_rate_N_F;
    double fatality_rate_W_S;
    double fatality_rate_W_F;
    double fatality_rate_B_S;
    double fatality_rate_B_F;


    //OTHERS
    double dead_S; // D_S: Dead
    double dead_F; // D_F: Dead
    double max_sum_infections_combined;
    double day_max_sum_infections_combined;
    double max_sum_hospitalized_combined;
    double day_max_sum_hospitalized_combined;
    double cumul_infections_combined;
    double cumul_hospitalization_combined;
    double cumul_deaths_combined;
    double average_screening_tests_used_combined;
    double cumul_screening_tests_used_combined;



    //NO IMMUNITY:
    double day_max_sum_infections_no_immunity;
    double day_max_sum_hospitalized_no_immunity;
    double cumul_infections_no_immunity;
    double cumul_hospitalization_no_immunity;
    double cumul_deaths_no_immunity;
    double max_sum_infections_no_immunity;
    double max_sum_hospitalized_no_immunity;
    double average_screening_tests_used_no_immunity;
    double cumul_screening_tests_used_no_immunity;


    //STUDENTS
    double  uninfected_N_S; //U_N,S: Uninfected, susceptible individuals
    double  exposed_asymptomatic_N_S; //E_N,S: Exposed, asymptomatic, non - infectious
    double  asymptomatic_infected_N_S; //A_N,S: Infected, asymptomatic
    double recovered_unknown_N_S; //RU_N,S: recovered but does not know
    double recovered_known_N_S; // RK_N,S:Immune (either immune because of recovering or from vaccine)
    double symptomatic_infected_N_S;//S_N,S: Infected, symptomatic(true) positive test result
    double true_positive_N_S;//TP_N,S: Infected, asymptomatic, (true) positive test result
    double false_positive_N_S;//FP_N,S: Uninfected, false positive result
    double false_positive_RU_N_S; // FPRU_N,S: Recovered unknown, detected, false positive
    double hospitalized_N_S; // H_N,S: Hospitalized

     // practical variables
    double previous_asymptomatic_infected_N_S; //A_N,S(t-1)
    double current_asymptomatic_N_S;  //A_N,S(t)
    double current_exposed_asymptomatic_N_S; //E_N,S(t)
    double previous_uninfected_N_S; //U_N,S(t-1)
    double current_uninfected_N_S; //U_N,S(t)
    double previous_RU_N_S; //RU_N,S(t-1)
    double current_RU_N_S; //RU_N,S(t)
    double current_RK_N_S; //RK_N,S(t)
    double current_hospitalized_N_S; //H_N,S(t)
    double current_TP_N_S; //TP_N,S(t)
    double current_symptomatic_N_S; // S_N,S (t)
    double cumul_asymptomatic_N_S;
    double max_infections_N_S;
    double day_max_infections_N_S;
    double infections_N_S;
    double cumul_infections_N_S;
    double cumul_hospitalization_N_S;
    double cumul_deaths_N_S;
    double ratio_N_S;
    double average_screening_tests_used_N_S;
    double cumul_screening_tests_used_N_S;
    double indicator_N_S;


    //FACULTY
    double  uninfected_N_F; //U_N,F: Uninfected, susceptible individuals
    double  exposed_asymptomatic_N_F; //E_N,F: Exposed, asymptomatic, non - infectious
    double  asymptomatic_infected_N_F; //A_N,F: Infected, asymptomatic
    double recovered_unknown_N_F; //RU_N,F: recovered but does not know
    double recovered_known_N_F; // RK_N,F:Immune (either immune because of recovering or from vaccine)
    double symptomatic_infected_N_F;//S_N,F: Infected, symptomatic(true) positive test result
    double true_positive_N_F;//TP_N,F: Infected, asymptomatic, (true) positive test result
    double false_positive_N_F;//FP_N,F: Uninfected, false positive result
    double false_positive_RU_N_F; // FPRU_N,F: Recovered unknown, detected, false positive
    double hospitalized_N_F; // H_N,F: Hospitalized
    double current_TP_N_F; //TP_N,F(t)
    double current_symptomatic_N_F; // S_N,F (t)

     // practical variables
    double previous_asymptomatic_infected_N_F; //A_N,F(t-1)
    double current_asymptomatic_N_F;  //A_N,F(t)
    double current_exposed_asymptomatic_N_F; //E_N,F(t)
    double previous_uninfected_N_F; //U_N,F(t-1)
    double current_uninfected_N_F; //U_N,F(t)
    double previous_RU_N_F; //RU_N,F(t-1)
    double current_RU_N_F; //RU_N,F(t)
    double current_RK_N_F; //RK_N,F(t)
    double current_hospitalized_N_F; //H_N,F(t)
    double cumul_hospitalization_N_F;
    double cumul_deaths_N_F;
    double cumul_asymptomatic_N_F;
    double max_infections_N_F;
    double day_max_infections_N_F;
    double infections_N_F;
    double cumul_infections_N_F;
    double ratio_N_F;
    double average_screening_tests_used_N_F;
    double cumul_screening_tests_used_N_F;
    double indicator_N_F;



    //WANING IMMUNITY:
    double day_max_sum_infections_waning_immunity;
    double day_max_sum_hospitalized_waning_immunity;
    double cumul_infections_waning_immunity;
    double cumul_hospitalization_waning_immunity;
    double cumul_deaths_waning_immunity;
    double max_sum_infections_waning_immunity;
    double max_sum_hospitalized_waning_immunity;
    double average_screening_tests_used_waning_immunity;
    double cumul_screening_tests_used_waning_immunity;


    //STUDENTS
    double  uninfected_W_S; //U_W,S: Uninfected, susceptible individuals
    double  exposed_asymptomatic_W_S; //E_W,S: Exposed, asymptomatic, non - infectious
    double  asymptomatic_infected_W_S; //A_W,S: Infected, asymptomatic
    double recovered_unknown_W_S; //RU_W,S: recovered but does not know
    double recovered_known_W_S; // RK_W,S:Immune (either immune because of recovering or from vaccine)
    double symptomatic_infected_W_S;//S_W,S: Infected, symptomatic(true) positive test result
    double true_positive_W_S;//TP_W,S: Infected, asymptomatic, (true) positive test result
    double false_positive_W_S;//FP_W,S: Uninfected, false positive result
    double false_positive_RU_W_S; // FPRU_W,S: Recovered unknown, detected, false positive
    double hospitalized_W_S; // H_W,S: Hospitalized

     // practical variables
    double previous_asymptomatic_infected_W_S; //A_W,S(t-1)
    double current_asymptomatic_W_S;  //A_W,S(t)
    double current_exposed_asymptomatic_W_S; //E_W,S(t)
    double previous_uninfected_W_S; //U_W,S(t-1)
    double current_uninfected_W_S; //U_W,S(t)
    double previous_RU_W_S; //RU_W,S(t-1)
    double current_RU_W_S; //RU_W,S(t)
    double current_RK_W_S; //RK_W,S(t)
    double current_hospitalized_W_S; //H_W,S(t)
    double current_TP_W_S; //TP_W,S(t)
    double current_symptomatic_W_S; // S_W,S (t)
    double cumul_asymptomatic_W_S;
    double max_infections_W_S;
    double day_max_infections_W_S;
    double infections_W_S;
    double cumul_infections_W_S;
    double cumul_hospitalization_W_S;
    double cumul_deaths_W_S;
    double ratio_W_S;
    double average_screening_tests_used_W_S;
    double cumul_screening_tests_used_W_S;
    double indicator_W_S;


    //FACULTY
    double  uninfected_W_F; //U_W,F: Uninfected, susceptible individuals
    double  exposed_asymptomatic_W_F; //E_W,F: Exposed, asymptomatic, non - infectious
    double  asymptomatic_infected_W_F; //A_W,F: Infected, asymptomatic
    double recovered_unknown_W_F; //RU_W,F: recovered but does not know
    double recovered_known_W_F; // RK_W,F:Immune (either immune because of recovering or from vaccine)
    double symptomatic_infected_W_F;//S_W,F: Infected, symptomatic(true) positive test result
    double true_positive_W_F;//TP_W,F: Infected, asymptomatic, (true) positive test result
    double false_positive_W_F;//FP_W,F: Uninfected, false positive result
    double false_positive_RU_W_F; // FPRU_W,F: Recovered unknown, detected, false positive
    double hospitalized_W_F; // H_W,F: Hospitalized
    double current_TP_W_F; //TP_W,F(t)
    double current_symptomatic_W_F; // S_W,F (t)

     // practical variables
    double previous_asymptomatic_infected_W_F; //A_W,F(t-1)
    double current_asymptomatic_W_F;  //A_W,F(t)
    double current_exposed_asymptomatic_W_F; //E_W,F(t)
    double previous_uninfected_W_F; //U_W,F(t-1)
    double current_uninfected_W_F; //U_W,F(t)
    double previous_RU_W_F; //RU_W,F(t-1)
    double current_RU_W_F; //RU_W,F(t)
    double current_RK_W_F; //RK_W,F(t)
    double current_hospitalized_W_F; //H_W,F(t)
    double cumul_hospitalization_W_F;
    double cumul_deaths_W_F;
    double cumul_asymptomatic_W_F;
    double max_infections_W_F;
    double day_max_infections_W_F;
    double infections_W_F;
    double cumul_infections_W_F;
    double ratio_W_F;
    double average_screening_tests_used_W_F;
    double cumul_screening_tests_used_W_F;
    double indicator_W_F;



    //BOOSTED IMMUNITY:
    double day_max_sum_infections_boosted_immunity;
    double day_max_sum_hospitalized_boosted_immunity;
    double cumul_infections_boosted_immunity;
    double cumul_hospitalization_boosted_immunity;
    double cumul_deaths_boosted_immunity;
    double max_sum_infections_boosted_immunity;
    double max_sum_hospitalized_boosted_immunity;
    double average_screening_tests_used_boosted_immunity;
    double cumul_screening_tests_used_boosted_immunity;


    //STUDENTS
    double  uninfected_B_S; //U_B,S: Uninfected, susceptible individuals
    double  exposed_asymptomatic_B_S; //E_B,S: Exposed, asymptomatic, non - infectious
    double  asymptomatic_infected_B_S; //A_B,S: Infected, asymptomatic
    double recovered_unknown_B_S; //RU_B,S: recovered but does not know
    double recovered_known_B_S; // RK_B,S:Immune (either immune because of recovering or from vaccine)
    double symptomatic_infected_B_S;//S_B,S: Infected, symptomatic(true) positive test result
    double true_positive_B_S;//TP_B,S: Infected, asymptomatic, (true) positive test result
    double false_positive_B_S;//FP_B,S: Uninfected, false positive result
    double false_positive_RU_B_S; // FPRU_B,S: Recovered unknown, detected, false positive
    double hospitalized_B_S; // H_B,S: Hospitalized

     // practical variables
    double previous_asymptomatic_infected_B_S; //A_B,S(t-1)
    double current_asymptomatic_B_S;  //A_B,S(t)
    double current_exposed_asymptomatic_B_S; //E_B,S(t)
    double previous_uninfected_B_S; //U_B,S(t-1)
    double current_uninfected_B_S; //U_B,S(t)
    double previous_RU_B_S; //RU_B,S(t-1)
    double current_RU_B_S; //RU_B,S(t)
    double current_RK_B_S; //RK_B,S(t)
    double current_hospitalized_B_S; //H_B,S(t)
    double current_TP_B_S; //TP_B,S(t)
    double current_symptomatic_B_S; // S_B,S (t)
    double cumul_asymptomatic_B_S;
    double max_infections_B_S;
    double day_max_infections_B_S;
    double infections_B_S;
    double cumul_infections_B_S;
    double cumul_hospitalization_B_S;
    double cumul_deaths_B_S;
    double ratio_B_S;
    double average_screening_tests_used_B_S;
    double cumul_screening_tests_used_B_S;
    double indicator_B_S;


    //FACULTY
    double  uninfected_B_F; //U_B,F: Uninfected, susceptible individuals
    double  exposed_asymptomatic_B_F; //E_B,F: Exposed, asymptomatic, non - infectious
    double  asymptomatic_infected_B_F; //A_B,F: Infected, asymptomatic
    double recovered_unknown_B_F; //RU_B,F: recovered but does not know
    double recovered_known_B_F; // RK_B,F:Immune (either immune because of recovering or from vaccine)
    double symptomatic_infected_B_F;//S_B,F: Infected, symptomatic(true) positive test result
    double true_positive_B_F;//TP_B,F: Infected, asymptomatic, (true) positive test result
    double false_positive_B_F;//FP_B,F: Uninfected, false positive result
    double false_positive_RU_B_F; // FPRU_B,F: Recovered unknown, detected, false positive
    double hospitalized_B_F; // H_B,F: Hospitalized
    double current_TP_B_F; //TP_B,F(t)
    double current_symptomatic_B_F; // S_B,F (t)

     // practical variables
    double previous_asymptomatic_infected_B_F; //A_B,F(t-1)
    double current_asymptomatic_B_F;  //A_B,F(t)
    double current_exposed_asymptomatic_B_F; //E_B,F(t)
    double previous_uninfected_B_F; //U_B,F(t-1)
    double current_uninfected_B_F; //U_B,F(t)
    double previous_RU_B_F; //RU_B,F(t-1)
    double current_RU_B_F; //RU_B,F(t)
    double current_RK_B_F; //RK_B,F(t)
    double current_hospitalized_B_F; //H_B,F(t)
    double cumul_hospitalization_B_F;
    double cumul_deaths_B_F;
    double cumul_asymptomatic_B_F;
    double max_infections_B_F;
    double day_max_infections_B_F;
    double infections_B_F;
    double cumul_infections_B_F;
    double ratio_B_F;
    double average_screening_tests_used_B_F;
    double cumul_screening_tests_used_B_F;
    double indicator_B_F;




};