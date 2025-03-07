/*
	clean_gait3d_pelvis213_torque.c

	Program to extract clean C source code from C file generated by Autolev program gait3d_pelvis213_torque.al
	This compiles into a Matlab MEX function.
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mex.h"

#include "gait3d_pelvis213_torque_multibody.h"
#define STRLEN 100 				// source lines are never longer than 100 characters

// function to determine whether this line of code has a nonzero element for matrix "name"
int nonzero(char *line, char *name) {
	if ( (strstr(line,name) != NULL) && (strstr(line,"] = 0;") == NULL) ) {
		return 1;
	}
	else {
		return 0;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	FILE *fid1, *fid2;
	char line[STRLEN];
	int i, copying;
	char *ptr;
    
    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    char *slash="\\";
    #else
    char *slash="/";
    #endif


    //===================================================================================================
    // gait3d_pelvis213_torque: Dynamics with derivatives
    //===================================================================================================

	// initialize the non-zero counters for the various matrices
	int nnz_dfk_dq = 0;
	int nnz_dfkdot_dq = 0;
	int nnz_df_dq = 0;
	int nnz_df_dqd = 0;
	int nnz_df_dqdd = 0;
	int nnz_df_dG = 0;

	// open the file that came from Autolev
    char filename1[100];
    strcpy(filename1, "gait3d_pelvis213_torque");
    strcat(filename1, slash);
    strcat(filename1, "gait3d_pelvis213_torque_raw.c");
	if ((fid1 = fopen(filename1,"r")) == NULL) {
		printf("Could not open %s\n", filename1);
		exit(1);
	}
	
    // write the clean C file
    char filename2[100];
    strcpy(filename2, "gait3d_pelvis213_torque");
    strcat(filename2, slash);
    strcat(filename2, "gait3d_pelvis213_torque_al.c");
	if ((fid2 = fopen(filename2,"w")) == NULL) {
		printf("Could not write %s\n", filename2);
		exit(1);
	}
	fprintf(fid2,"// This file was generated by clean_gait3d_pelvis213_torque.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include <math.h>\n");
	fprintf(fid2,"#include \"gait3d_pelvis213_torque_multibody.h\"\n");
	fprintf(fid2,"void gait3d_pelvis213_torque_al(                \n");
	fprintf(fid2,"    param_struct* par,             // input: model parameters\n");
	fprintf(fid2,"    double q[NDOF],                // input: gen.coords\n");
	fprintf(fid2,"    double qd[NDOF],               // input: gen.velocities\n");
	fprintf(fid2,"    double qdd[NDOF],              // input: gen.accelerations\n");
	fprintf(fid2,"    double G[NGRF],                // input: ground reaction forces\n");
	fprintf(fid2,"    double f[NDOF],                // output: gen.forces Q = f(q,qd,qdd,G)\n");
	fprintf(fid2,"    double df_dq[NDOF][NDOF],      // output: Jacobian df/dq\n");
	fprintf(fid2,"    double df_dqd[NDOF][NDOF],     // output: Jacobian df/dqd\n");
	fprintf(fid2,"    double df_dqdd[NDOF][NDOF],    // output: Jacobian df/dqdd (mass matrix)\n");
	fprintf(fid2,"    double df_dG[NDOF][NGRF],      // output: Jacobian df/dG\n");
	fprintf(fid2,"    double fk[NFK],                // output: forward kinematics FK(q)\n");
	fprintf(fid2,"    double dfk_dq[NFK][NDOF],      // output: Jacobian dFK/dq\n");
	fprintf(fid2,"    double fkdot[NFK],             // output: dFK/dt\n");
	fprintf(fid2,"    double dfkdot_dq[NFK][NDOF]) { // output: dFK/dt/dq\n");		
	
	// generate C code to copy q, qd, qdd into scalar variables
	for (i=0; i<NDOF; i++) {
		fprintf(fid2,"\tdouble q%1d = q[%1d];\n", i+1,i);		
		fprintf(fid2,"\tdouble q%1dp = qd[%1d];\n", i+1,i);	
		fprintf(fid2,"\tdouble q%1dpp = qdd[%1d];\n", i+1,i);	
	}
	
	// generate C code to copy ground reactions G[] into scalar variables G1...G24
	for (i=0; i<NGRF; i++) {
		fprintf(fid2,"\tdouble G%1d = G[%1d];\n", i+1,i);			
	}
	
	// generate C code to declare the airspeed components and magnitude
	fprintf(fid2,"\tdouble sx,sy,sz,s;\n");			
	
	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// allocate static (non-stack) space for the  Z[] array
		}
		if ((strcmp(line, "/* Evaluate constants */\n") == 0) || (strcmp(line, "/* Evaluate constants */\r\n") == 0)) 				// Z[] code starts here
			copying = 1;
		else if ((strcmp(line, "/* Evaluate output quantities */\n") == 0) || (strcmp(line, "/* Evaluate output quantities */\r\n") == 0))	// Z[] code ends here
			copying = 0;
		else if ((strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0) || (strcmp(line, "/* Write output to screen and to output file(s) */\r\n") == 0))		// kinematics code starts here
			copying = 1;
		else if ((strcmp(line, "  Encode[0] = 0.0;\n") == 0) || (strcmp(line, "  Encode[0] = 0.0;\r\n") == 0))	// and stop here
			copying = 0;
		else if (copying) {
			// update the non-zero counters
			nnz_dfk_dq 		+= nonzero(line, "dfk_dq[");
			nnz_dfkdot_dq 	+= nonzero(line, "dfkdot_dq[");
			nnz_df_dq 		+= nonzero(line, "df_dq[");
			nnz_df_dqd 		+= nonzero(line, "df_dqd[");
			nnz_df_dqdd 	+= nonzero(line, "df_dqdd[");
			nnz_df_dG 		+= nonzero(line, "df_dG[");
	
			// convert par__ into par->
			ptr = strstr(line, "par__");
			while (ptr != NULL) {
				strncpy(ptr, "par->", 5);
				ptr = strstr(line, "par__");
			}
			fputs(line, fid2);
		}
		
	}
	
    // close the input file
	fclose(fid1);
	
    // close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);

    //===================================================================================================
    // gait3d_pelvis213_torque_NoDer: Dynamics without derivatives
    //===================================================================================================

    // open the file that came from Autolev
    char filename1NoDer[100];
    strcpy(filename1NoDer, "gait3d_pelvis213_torque");
    strcat(filename1NoDer, slash);
    strcat(filename1NoDer, "gait3d_pelvis213_torque_NoDer_raw.c");
	if ((fid1 = fopen(filename1NoDer,"r")) == NULL) {
		printf("Could not open %s\n", filename1NoDer);
		exit(1);
	}
	
    // write the clean C file
    char filename2NoDer[100];
    strcpy(filename2NoDer, "gait3d_pelvis213_torque");
    strcat(filename2NoDer, slash);
    strcat(filename2NoDer, "gait3d_pelvis213_torque_NoDer_al.c");
	if ((fid2 = fopen(filename2NoDer,"w")) == NULL) {
		printf("Could not write %s\n", filename2NoDer);
		exit(1);
	}
	fprintf(fid2,"// This file was generated by clean_gait3d_pelvis213_torque.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include <math.h>\n");
	fprintf(fid2,"#include \"gait3d_pelvis213_torque_multibody.h\"\n");
	fprintf(fid2,"void gait3d_pelvis213_torque_NoDer_al(                \n");
	fprintf(fid2,"    param_struct* par,             // input: model parameters\n");
	fprintf(fid2,"    double q[NDOF],                // input: gen.coords\n");
	fprintf(fid2,"    double qd[NDOF],               // input: gen.velocities\n");
	fprintf(fid2,"    double qdd[NDOF],              // input: gen.accelerations\n");
	fprintf(fid2,"    double G[NGRF],                // input: ground reaction forces\n");
	fprintf(fid2,"    double f[NDOF],                // output: gen.forces Q = f(q,qd,qdd,G)\n");
	fprintf(fid2,"    double fk[NFK],                // output: forward kinematics FK(q)\n");
	fprintf(fid2,"    double dfk_dq[NFK][NDOF],      // output: Jacobian dFK/dq\n");
	fprintf(fid2,"    double fkdot[NFK],             // output: dFK/dt\n");
	fprintf(fid2,"    double dfkdot_dq[NFK][NDOF]) { // output: dFK/dt/dq\n");		
	
	// generate C code to copy q, qd, qdd into scalar variables
	for (i=0; i<NDOF; i++) {
		fprintf(fid2,"\tdouble q%1d = q[%1d];\n", i+1,i);		
		fprintf(fid2,"\tdouble q%1dp = qd[%1d];\n", i+1,i);	
		fprintf(fid2,"\tdouble q%1dpp = qdd[%1d];\n", i+1,i);	
	}
	
	// generate C code to copy ground reactions G[] into scalar variables G1...G24
	for (i=0; i<NGRF; i++) {
		fprintf(fid2,"\tdouble G%1d = G[%1d];\n", i+1,i);			
	}
	
	// generate C code to declare the airspeed components and magnitude
	fprintf(fid2,"\tdouble sx,sy,sz,s;\n");			
	
	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// allocate static (non-stack) space for the  Z[] array
		}
		if ((strcmp(line, "/* Evaluate constants */\n") == 0) || (strcmp(line, "/* Evaluate constants */\r\n") == 0)) 				// Z[] code starts here
			copying = 1;
		else if ((strcmp(line, "/* Evaluate output quantities */\n") == 0) || (strcmp(line, "/* Evaluate output quantities */\r\n") == 0))	// Z[] code ends here
			copying = 0;
		else if ((strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0) || (strcmp(line, "/* Write output to screen and to output file(s) */\r\n") == 0))		// kinematics code starts here
			copying = 1;
		else if ((strcmp(line, "  Encode[0] = 0.0;\n") == 0) || (strcmp(line, "  Encode[0] = 0.0;\r\n") == 0))	// and stop here
			copying = 0;
		else if (copying) {	
			// convert par__ into par->
			ptr = strstr(line, "par__");
			while (ptr != NULL) {
				strncpy(ptr, "par->", 5);
				ptr = strstr(line, "par__");
			}
			fputs(line, fid2);
		}
		
	}
	
    // close the input file
	fclose(fid1);
	
    // close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);

    //===================================================================================================
    // gait3d_pelvis213_torque_FK: Forward kinematics
    //===================================================================================================

    // open the file that came from Autolev
    char filename1FK[100];
    strcpy(filename1FK, "gait3d_pelvis213_torque");
    strcat(filename1FK, slash);
    strcat(filename1FK, "gait3d_pelvis213_torque_FK_raw.c");
	if ((fid1 = fopen(filename1FK,"r")) == NULL) {
		printf("Could not open %s\n", filename1FK);
		exit(1);
	}

    // write the clean C file
    char filename2FK[100];
    strcpy(filename2FK, "gait3d_pelvis213_torque");
    strcat(filename2FK, slash);
    strcat(filename2FK, "gait3d_pelvis213_torque_FK_al.c");
	if ((fid2 = fopen(filename2FK,"w")) == NULL) {
		printf("Could not write %s\n", filename2FK);
		exit(1);
	}
	fprintf(fid2,"// This file was generated by clean_gait3d_pelvis213_torque.c and contains C code generated by Autolev\n\n");
	fprintf(fid2,"#include <math.h>\n");
	fprintf(fid2,"#include \"gait3d_pelvis213_torque_multibody.h\"\n");
	fprintf(fid2,"void gait3d_pelvis213_torque_FK_al(                \n");
	fprintf(fid2,"    param_struct* par,             // input: model parameters\n");
	fprintf(fid2,"    double q[NDOF],                // input: gen.coords\n");
	fprintf(fid2,"    double qd[NDOF],               // input: gen.velocities\n");
	fprintf(fid2,"    double qdd[NDOF],              // input: gen.accelerations\n");
	fprintf(fid2,"    double fk[NFK],                // output: forward kinematics FK(q)\n");
	fprintf(fid2,"    double dfk_dq[NFK][NDOF],      // output: Jacobian dFK/dq\n");
	fprintf(fid2,"    double fkdot[NFK],             // output: dFK/dt\n");
	fprintf(fid2,"    double dfkdot_dq[NFK][NDOF]) { // output: dFK/dt/dq\n");		
	
	// generate C code to copy q, qd, qdd into scalar variables
	for (i=0; i<NDOF; i++) {
		fprintf(fid2,"\tdouble q%1d = q[%1d];\n", i+1,i);		
		fprintf(fid2,"\tdouble q%1dp = qd[%1d];\n", i+1,i);	
		fprintf(fid2,"\tdouble q%1dpp = qdd[%1d];\n", i+1,i);	
	}	
	
	// copy the necessary parts of C code from fid1 to fid2
	copying = 0;
	while (feof(fid1) == 0) {
		fgets(line, STRLEN, fid1);
		if (strncmp(line, "double   Pi,DEGtoRAD,RADtoDEG,z[", 32) == 0) {
			fprintf(fid2,"\tstatic double z[%d];\n",atoi(&line[32]));		// allocate static (non-stack) space for the  Z[] array
		}
		if ((strcmp(line, "/* Evaluate constants */\n") == 0) || (strcmp(line, "/* Evaluate constants */\r\n") == 0)) 				// Z[] code starts here
			copying = 1;
		else if ((strcmp(line, "/* Evaluate output quantities */\n") == 0) || (strcmp(line, "/* Evaluate output quantities */\r\n") == 0))	// Z[] code ends here
			copying = 0;
		else if ((strcmp(line, "/* Write output to screen and to output file(s) */\n") == 0) || (strcmp(line, "/* Write output to screen and to output file(s) */\r\n") == 0))		// kinematics code starts here
			copying = 1;
		else if ((strcmp(line, "  Encode[0] = 0.0;\n") == 0) || (strcmp(line, "  Encode[0] = 0.0;\r\n") == 0))	// and stop here
			copying = 0;
		else if (copying) {
			// convert par__ into par->
			ptr = strstr(line, "par__");
			while (ptr != NULL) {
				strncpy(ptr, "par->", 5);
				ptr = strstr(line, "par__");
			}
			fputs(line, fid2);
		}
		
	}
	
    // close the input file
	fclose(fid1);
	
    // close the output file
	fprintf(fid2,"}\n");
	fclose(fid2);
		
    // report the non-zero numbers for the sparse matrices
	printf("Nonzero elements for dfk_dq:       %d\n", nnz_dfk_dq);
	printf("Nonzero elements for dfkdot_dq:    %d\n", nnz_dfkdot_dq);
	printf("Nonzero elements for df_dq:        %d\n", nnz_df_dq);
	printf("Nonzero elements for df_dqd:       %d\n", nnz_df_dqd);
	printf("Nonzero elements for df_dqdd:      %d\n", nnz_df_dqdd);
	printf("Nonzero elements for df_dG:        %d\n", nnz_df_dG);
		
}			// END OF MAIN PROGRAM
