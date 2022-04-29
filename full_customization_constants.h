#define MAX_UNSIGNED_SHORT	std::numeric_limits<unsigned short>::max()
#define MAX_INT				std::numeric_limits<int>::max()
#define MIN_INT				std::numeric_limits<int>::min()
#define MAX_UNSIGNED_INT	std::numeric_limits<unsigned int>::max()
#define MAX_LONG_INT		std::numeric_limits<unsigned long>::max()
#define MIN_LONG_INT		std::numeric_limits<long int>::min()
#define MAX_FLOAT			std::numeric_limits<float>::max()
#define MIN_FLOAT			std::numeric_limits<float>::min()
#define MAX_DOUBLE			std::numeric_limits<double>::max()
#define MIN_DOUBLE			std::numeric_limits<double>::min()
//#define EPSILON				1E-12
#define EPSILON_SMALL		1E-2
#define ABS(x) 				( (x < 0) ? -(x) : (x))

//#define NUM_OPOS			58
//#define NUM_REGIONS		11

//#define DEBUG				0	// setting value to 1 will print many things to the debug file
//#define TEMPDEBUG			0	// setting value to 1 will print temporary test functions
//
//#define REPORTING_INTERVAL 	7		// screen progress display
//
//#define earthRadiusKm 6371.0	// to find distance between centers given their longitudes and latitudes (Haversine formula)#pragma once


// ASSUMING CYCLES OF 1 DAY

//SIR parameters
#define THETA       1.0/(3*3.5)     //3.5 if 50/50 (mean incubation time) or rate at which exposed individuals advance to the asymptomatic, infectious compartment

#define Se			0.8   			  //70%-80%... sensitivity of the screening test
#define Sp          0.98 			  //98%-99.7% ... specificity of the screening test
//#define X           10               // 5, 10, 25 ... number of imported infections in a given exogenous shock

#define X_N_S           15.0    
#define X_N_F           15.0/15               
#define X_W_S           10.0    
#define X_W_F           10.0/15 
#define X_B_S           5.0   
#define X_B_F           5.0/15 


//Base case: Sensitivity = 0.8,   Specificity=98% and    X=10/week => jama
//Best case: Sensitivity = 0.8 ,  Specificity=99.7% and  X=5/week
//Worst case: Sensitivity = 0.7 , Specificity=98% and    X=25/week

#define RHO_N         1.0/(3*5.0)      //5 days new guideline      // (1/rho is the time to recover) rate at which infected individuals recover from disease and are removed
#define SIGMA_N       0.0286             //(rate of symptom development) rate of symptom onset for infected individuals
#define RHO_W         1.0/(3*5.0)      //5 days new guideline      // (1/rho is the time to recover) rate at which infected individuals recover from disease and are removed
#define SIGMA_W       0.0286            	 //(rate of symptom development) rate of symptom onset for infected individuals
#define RHO_B         1.0/(3*5.0)      //5 days new guideline      // (1/rho is the time to recover) rate at which infected individuals recover from disease and are removed
#define SIGMA_B       0.0286            	 //(rate of symptom development) rate of symptom onset for infected individuals

#define MU          1.0/(3*1.0)          //(1/mu is the time to false positive return) rate at which false positives are returned to the Uninfected compartment
