/*
 * -----------------------------------------------------------------
 * $Revised from cvRoberts_FSA_dns.c in cvodes example$
 * $Date: 2013/11 $
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES. The problem is from chemical
 * kinetics, and consists of equations:
 *    dy/dt = f(y)
 * on the interval from t = 0.0 to t = t_end, with initial
 * conditions y(0). The problem is stiff.
 *
 * This program solves the problem with the BDF method, Newton
 * iteration with the CVODES dense linear solver.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters.
 * The sensitivity right hand side is given analytically through the
 * user routine fS (of type SensRhs1Fn).
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on TRUE or FALSE,
 * respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % cvsRoberts_FSA_dns -nosensi
 * If sensitivities are to be computed:
 *    % cvsRoberts_FSA_dns -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

/* Accessor macros */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

/* Problem Constants */

/* INSERT GENERATED CONSTANTS */
#define RTOL  RCONST(1.0e-6)   /* scalar relative tolerance            */
#define ATOL  RCONST(1.0e-9)    /* scalar absolute tolerance components */
#define MXSTEPS 10000           /* max steps before tout */
/* END OF GENERATED CONSTANTS */

#define ZERO  RCONST(0.0)

/* Type : UserData */

typedef struct {
  realtype p[NUMPAR];           /* problem parameters */
} *UserData;

/* Initialize y */

static int InitY(realtype t, N_Vector y);

/* Initialize UserData */

static int InitUserData(void *user_data);

/* Prototypes of functions by CVODES */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

// user supplied jacobian function not currently supported
//static int Jac(long int N, realtype t,
//               N_Vector y, N_Vector fy, 
//               DlsMat J, void *user_data, 
//               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

//static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
//              int iS, N_Vector yS, N_Vector ySdot, 
//              void *user_data, N_Vector tmp1, N_Vector tmp2);

static int ewt(N_Vector y, N_Vector w, void *user_data);

/* Prototypes of private functions */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth, 
                        booleantype *err_con,
                        realtype *T0, realtype *T1, int *NOUT, char** OUTPUT_FILENAME);
static void WrongArgs(char *name);
static void PrintOutputHeader(FILE *output_file, booleantype sensi, int sensi_meth, booleantype err_con, realtype *p, int *plist);
static void PrintOutput(realtype t, N_Vector y, FILE *output_file);
static void PrintOutputS(N_Vector *yS, N_Vector y, realtype *p, int *plist, FILE *output_file);
static void PrintFinalStats(void *cvode_mem, booleantype sensi);
static int check_flag(void *flagvalue, char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
  void *cvode_mem;
  UserData data;
  realtype t, tout, tinterval;
  N_Vector y;
  int iout, flag;
  int i;

  realtype T0, T1;
  int NOUT;
  char *OUTPUT_FILENAME;
  FILE *output_file;

  realtype pbar[NUMPAR_SENSI];

/* INSERT GENERATED PLIST */
/* END OF GENERATED PLIST */

  int is; 
  N_Vector *yS;
  booleantype sensi, err_con;
  int sensi_meth;

  cvode_mem = NULL;
  data      = NULL;
  y         = NULL;
  yS        = NULL;

  /* Process arguments */
  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con, &T0, &T1, &NOUT, &OUTPUT_FILENAME);

#ifdef DEBUG
  fprintf(stdout, "\nstart_time:%f\nend_time:%f\noutput_intervals:%d\noutput_file:%s\n", T0, T1, NOUT, OUTPUT_FILENAME);
#endif

  t = T0;

  /* User data structure */
  data = (UserData) malloc(sizeof *data);
  if (check_flag((void *)data, "malloc", 2)) return(1);
  InitUserData(data);

  /* Initial conditions */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  InitY(t, y);

  /* Create CVODES object */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Allocate space for CVODES */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Use private function to compute error weights */
  flag = CVodeWFtolerances(cvode_mem, ewt);
  if (check_flag(&flag, "CVodeSetEwtFn", 1)) return(1);

  /* Attach user data */
  flag = CVodeSetUserData(cvode_mem, data);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Attach linear solver */
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, "CVDense", 1)) return(1);

// user supplied jacobian function not currently supported
//  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
//  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);
//
//  printf("\n3-species chemical kinetics problem\n");

  /* Sensitivity-related settings */
  if (sensi) {

    for(i=0;i<NUMPAR_SENSI;++i){
       pbar[i] = data->p[plist[i]];
    }

    yS = N_VCloneVectorArray_Serial(NUMPAR_SENSI, y);
    if (check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0)) return(1);
    for (is=0;is<NUMPAR_SENSI;is++) N_VConst(ZERO, yS[is]);

    flag = CVodeSensInit1(cvode_mem, NUMPAR_SENSI, sensi_meth, NULL, yS);
    if(check_flag(&flag, "CVodeSensInit", 1)) return(1);

    flag = CVodeSensEEtolerances(cvode_mem);
    if(check_flag(&flag, "CVodeSensEEtolerances", 1)) return(1);

    flag = CVodeSetSensErrCon(cvode_mem, err_con);
    if (check_flag(&flag, "CVodeSetSensErrCon", 1)) return(1);

    flag = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
    if (check_flag(&flag, "CVodeSetSensParams", 1)) return(1);
  }
 
  if((output_file = fopen(OUTPUT_FILENAME, "r")) == NULL){
      output_file = fopen(OUTPUT_FILENAME, "w");
  }else{
    fclose(output_file);
    fprintf(stderr, "\nSUNDIALS_ERROR: open output file %s failed - the file already existed\n\n", OUTPUT_FILENAME);
    return(1);
  }
  if(output_file == NULL){
    fprintf(stderr, "\nSUNDIALS_ERROR: open output file %s failed - returned NULL pointer\n\n", OUTPUT_FILENAME);
    return(1);
  }
 
  /* In loop over output points, call CVode, print results, test for error */

  iout = 0; tinterval = (T1 - T0)/NOUT;
  tout = T0 + tinterval;
  PrintOutputHeader(output_file, sensi, sensi_meth, err_con, data->p, plist);
  PrintOutput(t, y, output_file);
    fprintf(output_file,"\n-----------------------------------------");
    fprintf(output_file,"------------------------------\n");

  for (iout=0, tout=T0+tinterval; iout < NOUT; iout++, tout += tinterval) {

    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (check_flag(&flag, "CVode", 1)) break;

    PrintOutput(t, y, output_file);

    if (sensi) {
      flag = CVodeGetSens(cvode_mem, &t, yS);
      if (check_flag(&flag, "CVodeGetSens", 1)) break;
      PrintOutputS(yS, y, data->p, plist, output_file);
    } 
//    fprintf(output_file,"-----------------------------------------");
//    fprintf(output_file,"------------------------------\n");

    fprintf(output_file, "\n");
  }

  /* Print final statistics */
  PrintFinalStats(cvode_mem, sensi);

  /* Free memory */

  N_VDestroy_Serial(y);                    /* Free y vector */
  if (sensi) {
    N_VDestroyVectorArray_Serial(yS, NUMPAR_SENSI);  /* Free yS vector */
  }
  free(data);                              /* Free user data */
  CVodeFree(&cvode_mem);                   /* Free CVODES memory */

  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * InitY routine. Initialize y.
 */
static int InitY(realtype t, N_Vector y)
{
 /* INSERT GENERATED INITIAL CONDITION */

 /* END OF GENERATED INITIAL CONDITION */

  return(0);
}

/*
 * InitUserData routine. Initialize user_data.
 */
static int InitUserData(void *user_data)
{
 UserData data;

 data = (UserData) user_data;

 /* INSERT GENERATED PARAMETERS LIST */

 /* END OF GENERATED PARAMETERS LIST */

  return(0);
}

/*
 * f routine. Compute f(t,y). 
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData data;

  data = (UserData) user_data;

 /* INSERT GENERATED RIGHT HAND SIDE */

 /* END OF GENERATED RIGHT HAND SIDE */

  return(0);
}


/* 
 * Jacobian routine. Compute J(t,y). 
 */
/*
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2, y3;
  UserData data;
  realtype p1, p2, p3;
 
  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  data = (UserData) user_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
 
  IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
  IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
                      IJth(J,3,2) = 2*p3*y2;

  return(0);
}
*/ 
/* 
 * fS routine. Compute sensitivity r.h.s. 
 */
/*
static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
              int iS, N_Vector yS, N_Vector ySdot, 
              void *user_data, N_Vector tmp1, N_Vector tmp2)
{
  UserData data;
  realtype p1, p2, p3;
  realtype y1, y2, y3;
  realtype s1, s2, s3;
  realtype sd1, sd2, sd3;

  data = (UserData) user_data;
  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  y1 = Ith(y,1);  y2 = Ith(y,2);  y3 = Ith(y,3);
  s1 = Ith(yS,1); s2 = Ith(yS,2); s3 = Ith(yS,3);

  sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
  sd3 = 2*p3*y2*s2;
  sd2 = -sd1-sd3;

  switch (iS) {
  case 0:
    sd1 += -y1;
    sd2 +=  y1;
    break;
  case 1:
    sd1 +=  y2*y3;
    sd2 += -y2*y3;
    break;
  case 2:
    sd2 += -y2*y2;
    sd3 +=  y2*y2;
    break;
  }
  
  Ith(ySdot,1) = sd1;
  Ith(ySdot,2) = sd2;
  Ith(ySdot,3) = sd3;

  return(0);
}
*/
/*
 * EwtSet function. Computes the error weights at the current solution.
 */

static int ewt(N_Vector y, N_Vector w, void *user_data)
{
  int i;
  realtype yy, ww, rtol, atol;

  rtol    = RTOL;
  atol    = ATOL;

  for (i=0; i<NEQ; i++) {
    yy = Ith(y,i+1);
    ww = rtol * ABS(yy) + atol;  
    if (ww <= 0.0) return (-1);
    Ith(w,i+1) = 1.0/ww;
  }

  return(0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Process and verify arguments to cvsfwddenx.
 */

static void ProcessArgs(int argc, char *argv[], 
                        booleantype *sensi, int *sensi_meth, booleantype *err_con,
                        realtype *T0, realtype *T1, int *NOUT, char** OUTPUT_FILENAME)
{
  *sensi = FALSE;
  *sensi_meth = -1;
  *err_con = FALSE;

  if (argc < 6) WrongArgs(argv[0]);

  if (strcmp(argv[1],"-nosensi") == 0)
    *sensi = FALSE;
  else if (strcmp(argv[1],"-sensi") == 0)
    *sensi = TRUE;
  else
    WrongArgs(argv[0]);
  
  if (*sensi) {

    if (argc != 8)
      WrongArgs(argv[0]);

    if (strcmp(argv[2],"sim") == 0)
      *sensi_meth = CV_SIMULTANEOUS;
    else if (strcmp(argv[2],"stg") == 0)
      *sensi_meth = CV_STAGGERED;
    else if (strcmp(argv[2],"stg1") == 0)
      *sensi_meth = CV_STAGGERED1;
    else 
      WrongArgs(argv[0]);

    if (strcmp(argv[3],"t") == 0)
      *err_con = TRUE;
    else if (strcmp(argv[3],"f") == 0)
      *err_con = FALSE;
    else
      WrongArgs(argv[0]);

    *T0 = (realtype)atof(argv[4]);
    *T1 = (realtype)atof(argv[5]);
    *NOUT = atoi(argv[6]);
    *OUTPUT_FILENAME = argv[7];

  } else {

    if (argc != 6)
      WrongArgs(argv[0]);

    *T0 = (realtype)atof(argv[2]);
    *T1 = (realtype)atof(argv[3]);
    *NOUT = atoi(argv[4]);
    *OUTPUT_FILENAME = argv[5];

  }

}

static void WrongArgs(char *name)
{
    printf("\nUsage: %s [-nosensi] [-sensi sensi_meth err_con] (start_time) (end_time) (output_intervals) (output_filename) (keep_labl)\n",name);
    printf("         choose one from -nosensi or -sensi\n");
    printf("         sensi_meth = sim, stg, or stg1\n");
    printf("         err_con    = t or f\n");
//    printf("         keep_label = 0 (donot) or 1 (do) (keep labels in output)\n");
    
    exit(0);
}

/*
 * Print Output Header (species names)
 */
static void PrintOutputHeader(FILE *output_file, booleantype sensi, int sensi_meth, booleantype err_con, realtype *p, int *plist)
{
  int i,is;
 /* INSERT GENERATED OUTPUT SPECIES NAMES */

 /* END OF GENERATED OUTPUT SPECIES NAMES */

/* INSERT GENERATED OUTPUT PARAMETERS NAMES */

/* END OF GENERATED OUTPUT PARAMETERS NAMES */

  if(sensi){
    fprintf(output_file, "Sensitivity: YES ");
    if(sensi_meth == CV_SIMULTANEOUS)   
      fprintf(output_file,"( SIMULTANEOUS +");
    else 
      if(sensi_meth == CV_STAGGERED) fprintf(output_file,"( STAGGERED +");
      else                           fprintf(output_file,"( STAGGERED1 +");   
    if(err_con) fprintf(output_file," FULL ERROR CONTROL )");
    else        fprintf(output_file," PARTIAL ERROR CONTROL )");

  } else {

    fprintf(output_file,"Sensitivity: NO ");

  }
  fprintf(output_file, "\ntime\t\t");
  for(i=0;i<OUTPUT_SPECIES_NUMBER;++i){
    fprintf(output_file, "%s\t\t", species_names[i]);
  }
  if(sensi){
    for(is=0;is<NUMPAR_SENSI;++is){
      for(i=0;i<OUTPUT_SPECIES_NUMBER;++i){
        fprintf(output_file, "%s:%s\t\t", species_names[i], parameters_names[is]);
      }
    }
  }

  fprintf(output_file, "\n");

  for(i = 0; i < NUMPAR_SENSI; i++)
    {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(output_file, "%0.4Le\t", p[plist[i]]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(output_file, "%0.4le\t", p[plist[i]]);
#else
      fprintf(output_file, "%0.4e\t", p[plist[i]]);
#endif
    }

  fprintf(output_file, "\n");
}

/*
 * Print current t, step count, order, stepsize, and solution.
 */

static void PrintOutput(realtype t, N_Vector y, FILE *output_file)
{
  int i;
 /* INSERT GENERATED OUTPUT SPECIES INDEX */

 /* END OF GENERATED OUTPUT SPECIES INDEX */

  long int nst;
  int qu, flag;
  realtype *ydata;
  
  ydata = NV_DATA_S(y);

#if defined(SUNDIALS_EXTENDED_PRECISION)
  fprintf(output_file, "%0.4Le\t", t);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  fprintf(output_file, "%0.4le\t", t);
#else
  fprintf(output_file, "%0.4e\t", t);
#endif

//  fprintf(output_file," Solution:\t");

#if defined(SUNDIALS_EXTENDED_PRECISION)
  for(i=0;i<OUTPUT_SPECIES_NUMBER;++i){
    fprintf(output_file, "%14.6Le\t", ydata[species_index[i]]);
  }
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  for(i=0;i<OUTPUT_SPECIES_NUMBER;++i){
    fprintf(output_file, "%14.6le\t", ydata[species_index[i]]);
  }
#else
  for(i=0;i<OUTPUT_SPECIES_NUMBER;++i){
    fprintf(output_file, "%14.6e\t", ydata[species_index[i]]);
  }
#endif

//  fprintf(output_file, "\n");

}

/* 
 * Print sensitivities.
*/

static void PrintOutputS(N_Vector *yS, N_Vector y, realtype *p, int *plist, FILE *output_file)
{
  realtype *sdata;
  int i, is;
 /* INSERT GENERATED OUTPUT SPECIES INDEX */

 /* END OF GENERATED OUTPUT SPECIES INDEX */

  for(is=0;is<NUMPAR_SENSI;++is){
    sdata = NV_DATA_S(yS[is]);
//    fprintf(output_file, "\t\t Sensitivity %d:\t", is+1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    for(i=0;i<OUTPUT_SPECIES_NUMBER;++i){
      fprintf(output_file, "%14.6Le\t", sdata[species_index[i]]);// * p[plist[is]] / yp[species_index[i]]
    }
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    for(i=0;i<OUTPUT_SPECIES_NUMBER;++i){
      fprintf(output_file, "%14.6le\t", sdata[species_index[i]]);// * p[plist[is]] / yp[species_index[i]]
    }
#else
    for(i=0;i<OUTPUT_SPECIES_NUMBER;++i){
      fprintf(output_file, "%14.6e\t", sdata[species_index[i]]); // * p[plist[is]] / yp[species_index[i]]
    }
#endif

//    fprintf(output_file, "\n");
  }
  
//  fprintf(output_file, "\n");
}

/* 
 * Print some final statistics from the CVODES memory.
 */

static void PrintFinalStats(void *cvode_mem, booleantype sensi)
{
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nje, nfeLS;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  if (sensi) {
    flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
    check_flag(&flag, "CVodeGetSensNumRhsEvals", 1);
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1);
    flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    check_flag(&flag, "CVodeGetSensNumLinSolvSetups", 1);
    flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
    check_flag(&flag, "CVodeGetSensNumErrTestFails", 1);
    flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
    check_flag(&flag, "CVodeGetSensNumNonlinSolvIters", 1);
    flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
    check_flag(&flag, "CVodeGetSensNumNonlinSolvConvFails", 1);
  }

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  printf("\nFinal Statistics\n\n");
  printf("nst     = %5ld\n\n", nst);
  printf("nfe     = %5ld\n",   nfe);
  printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
  printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

  if(sensi) {
    printf("\n");
    printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
    printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
    printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
  }

  printf("\n");
  printf("nje    = %5ld    nfeLS     = %5ld\n", nje, nfeLS);

}

/* 
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
