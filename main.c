/**
 * @mainpage Rational function approximant to matrix exponential using Carathéodory-Fejér method.
 *
 * @author <a href="mailto:niemeyer@case.edu">Kyle E. Niemeyer</a>
 * 
 * Returns poles (\f$ \theta_j \f$) and residues (\f$ \alpha_j \f$) for the rational function
 * (partial fraction) approximation to the matrix exponential:
 * \f$ \exp(A) = \sum_{j = 1}^n \frac{\alpha_j}{A - \theta_j I} \f$
 * 
 */

/** Primary file
 * \file main.c
 *
 * \author Kyle E. Niemeyer
 * \date 07/19/2012
 *
 * Contains main function (driver).
 */

#include <stdlib.h>
#include <stdio.h>

void cf ( int n, double* poles_r, double* poles_i, double* res_r, double* res_i );

/** Main function, gets poles and residuals for best rational approximant to matrix exponential.
 * 
 * Given the type of (number of terms in) the approximant, returns the
 * poles and residues for best rational approximant to the matrix exponential.
 * 
 * \param[in]		argc	command line argument count (1 or 2)
 * \param[in]		argv	command line argument vector; either empty or contains n 
 *							(type of rational approximation)
 */
int main ( int argc, char *argv[] ) {
	
	// default type of approximant, (n, n)
	int n = 10;
	
	// look for command line argument for n
	if (argc == 2) {
		n = atoi (argv[1]);
		
		if (n == 0) {
			printf ("Error in command line argument, should be a single nonzero integer.\n");
			exit (1);
		}
	}
	
	double *poles_r = (double*) calloc (n, sizeof(double));
	double *poles_i = (double*) calloc (n, sizeof(double));
	double *res_r = (double*) calloc (n, sizeof(double));
	double *res_i = (double*) calloc (n, sizeof(double));
	
	cf ( n, poles_r, poles_i, res_r, res_i );
	
	printf("Type (%i, %i) best rational function approximant to matrix exponential.\n", n, n);
	
	printf("poles: \n");
	for (int i = 0; i < n; ++i) {
		printf("%.15e + %.15e i\n", poles_r[i], poles_i[i]);
	}
	printf("\n");
	
	printf("residuals: \n");
	for (int i = 0; i < n; ++i) {
		printf("%.15e + %.15e i\n", res_r[i], res_i[i]);
	}
	printf("\n");
	
	free (poles_r);
	free (poles_i);
	free (res_r);
	free (res_i);
	
	return 0;
}