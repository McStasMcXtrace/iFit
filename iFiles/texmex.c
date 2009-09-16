/*-----------------------------program texmex.c-----------------------------*/

/*
% [strutcs] = looktxt('filename [options]') Import text data
*/

/*
 * compile with : mex -O -output looktxt texmex.c
 *           or : mex -v -argcheck -output looktxt texmex.c
 * content: C language, MEX library functions
 * tab = 2 chars
 */

#include <mex.h>	/* include MEX library for Matlab */

#define printf  mexPrintf	/* Addapt looktxt.c code to Mex syntax */
#define malloc  mxMalloc
#define realloc mxRealloc
#define calloc  mxCalloc
#define version "1.0.8 (MeX)"
/* #define free mxFree  */
#define free NoOp

#define main looktxt		/* change stand-alone source to a Matlab usable library */
#define print_stderr mexPrintf
#define argc carg
#define argv varg
#define exit(ret) { char msg[1024]; sprintf(msg, "Looktxt/mex exited with code %i", ret); if (ret) mexErrMsgTxt(msg); }
#define TEXMEX

#include "looktxt.c" 		/* makes all the job */

int NoOp(char *pointer)
{
  return 0;
}

void mexFunction(int nlhs, mxArray *plhs[],
                        int nrhs, const mxArray *prhs[])
{
  char *varg[MAX_LENGTH];  /* build a fake argv array */
  int  carg=1;
  int i;           /* rhs argument index */
  long NumFiles; /* number of files processes */
  FILE *fout;
  char has_output_file=0;

  varg[0] = (char*)mxMalloc(MAX_LENGTH);
  strcpy(varg[0],"looktxt");

  /* check in/out parameters */

  if (nlhs > 4)
	  mexErrMsgTxt("looktxt : Too many output arguments (4 max).");

  /* allocate memory */

  for (i = 0; i < nrhs; i++)
  {
    int   buflen;        /* number of arguments for call */
    char *InputString;   /* input string which is then split into arguments */
    int   status;        /* flag for getting input string */
    char  EndFlag = 0;   /* flag to exit while loop */
    char *StartLexeme, *EndLexeme;
    char *EndString;
    char  lexeme[MAX_LENGTH];  /* current argument, stored in varg */

    if (mxIsChar(prhs[i]) != 1)
    {
      mexPrintf("looktxt/mex : argument %i\n", i);
      mexErrMsgTxt("looktxt/mex : Input should be strings");
    }

    buflen      = (mxGetM(prhs[i])*mxGetN(prhs[i]))+1;
    InputString = (char*)mxMalloc(buflen+64);
    if (InputString == NULL)
    {
      mexPrintf("looktxt/mex : argument %i. Size %i\n", i, buflen);
      mexErrMsgTxt("looktxt/mex : can not allocate memory for input string\n");
    }
    status      = mxGetString(prhs[i], InputString, buflen);
    if (status != 0)
    {
      mexPrintf("looktxt/mex : argument %i. Status %i\n", i, status);
      mexErrMsgTxt("looktxt/mex : can not get input parameter\n");
    }

    /* cut input string into separated arguments for main(argc,argv) syntax */

    EndFlag = 0;
    StartLexeme = InputString;
    EndString = InputString+strlen(InputString);

    while ((EndFlag == 0) && (carg < MAX_LENGTH-1))
    {
      if (*StartLexeme == ' ')
        while (*StartLexeme == ' ' && StartLexeme < EndString)
        {       /* look for first non ' ' : StartLexeme */
          StartLexeme++;  /* pass all spaces */
        }
      /* now StartLexeme points on a non space or is at end */

      if (*StartLexeme == '\0' || StartLexeme >= EndString) EndFlag = 1;
      else
      {
        /* look for position of first next ' ' : EndLexeme */
        EndLexeme = strchr(StartLexeme+1, ' ');

        if (EndLexeme == NULL)
          EndLexeme = EndString;

        if (EndLexeme - StartLexeme > 0)
        {
          strncpy(lexeme, StartLexeme, EndLexeme - StartLexeme+1);
          lexeme[EndLexeme - StartLexeme] = '\0';

          StartLexeme = EndLexeme+1;
        }
      }

      if (strlen(lexeme) != 0 && lexeme != NULL && EndFlag == 0)
      {
        varg[carg] = (char*)mxMalloc(strlen(lexeme)+64);
        strcpy(varg[carg], lexeme);
        if (!strcmp(lexeme,"-o") || !strcmp(lexeme,"--outfile"))
          has_output_file=1;
        carg++;
      }
      else
        EndFlag = 1;
    }

  } /* end for nrhs */
  
  mexSetTrapFlag(1);
  
  char *tempname=NULL;
  if (!has_output_file) {
    mxArray *CellElement[1];
    mexCallMATLAB(1, CellElement, 0, NULL, "tempname");
    tempname=mxArrayToString(CellElement[0]);
    varg[carg] = (char*)mxMalloc(strlen(tempname)+64);
    sprintf(varg[carg], "--outfile=%s", tempname);
    mxDestroyArray(CellElement[0]);
    carg++;
  }
  for (i=0; i<carg; i++)
  	printf("%s ", varg[i]);
  printf("\n");

  /* link integer to matlab double */
  NumFiles = looktxt(carg,varg);  /* call main routine, return -1 in case of error, or the number of processes files */
  /* the output file names are stored in the global options.files_to_convert_Target array of char* */

  if (NumFiles > 0) {
    plhs[0]=mxCreateCellMatrix(1,NumFiles);
  } else plhs[0]=mxCreateCellMatrix(1,0);

  if (NumFiles > 0)
  for (i=0; i< NumFiles; i++) {
    /* get output file name */
    /* mfile=file.Source and func=file.RootName */
    struct file_struct output_file=TexMex_Target_Array[i];
    if (output_file.TargetTxt && strlen(output_file.TargetTxt)) {
      /* first method: evaluate function */
      mxArray *CellElement[1];
      struct fileparts_struct parts = fileparts(output_file.TargetTxt);
      char addpath[MAX_LENGTH];
      sprintf(addpath, "addpath('%s')", parts.Path);
      mexEvalString(addpath);
      if (!mexCallMATLAB(1, CellElement, 0, NULL, output_file.RootName)) {
        mxSetCell(plhs[0], i, CellElement[0]);
      } else {
        long filesize=0;
        char *filestr=NULL;

        /* read content */
        struct stat stfile;
        int status = stat(output_file.TargetTxt, &stfile);
        if (status) {
          mexPrintf("looktxt/mex: Warning : unable to access file %s\n", output_file.TargetTxt);
        } else {
          filesize = stfile.st_size;
          filestr  = (char*)mxMalloc(filesize+64);
          if (!filestr) {
            mexPrintf("looktxt/mex: Warning : unable to malloc %ld bytes for file (%s)\n", filesize, output_file.TargetTxt);
          } else {
            strcpy(filestr, "");
            fout = fopen(output_file.TargetTxt, "r");
            if (!fout) {
              mexPrintf("looktxt/mex: Warning : unable to open file %s for reading\n", output_file.TargetTxt);
            } else {
              if (!fread(filestr, 1, filesize, fout)) {
                mexPrintf("looktxt/mex: Warning : unable to read file %s\n", output_file.TargetTxt);
              }
              fclose(fout);
            } /* if fopen ok */
          } /* if malloc file string ok */
        } /* if stat ok */
        if (filesize && filestr && strlen(filestr)) {
          mxArray *FileString = mxCreateString(filestr);
          mxSetCell(plhs[0], i, FileString);
          free(filestr);
        }
      } /* if output_file */
    }
  } /* for i */
  
  if(NumFiles == 1) {
    mxArray *ptr;
    ptr= mxGetCell(plhs[0], 0);
    plhs[0] = mxDuplicateArray(ptr);
  }
  
  if (tempname) remove(tempname);

  for (i=0; i < carg; free(varg[i++]));
}
