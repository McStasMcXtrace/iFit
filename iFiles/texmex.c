 /*
    Looktxt: Search and export numerics in a text/ascii file
    Copyright (C) 2009  E. Farhi <farhi at ill.eu>, Institut Laue Langevin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*-----------------------------program texmex.c-----------------------------*/

/*
% [cell] = looktxt('filename', '[options]') Import text data
*/

/*
 * compile with : mex -O -output looktxt texmex.c
 *           or : mex -v -argcheck -output looktxt texmex.c
 *           or : eval(['mex -v -output looktxt texmex.c -L"' fullfile(matlabroot,'sys','lcc','lib') '" -lcrtdll' ])
 * content: C language, MEX library functions
 * tab = 2 chars
 *
 * The LCC compiler does not handle is* functions correctly. These are then 
   redefined in looktxt.c
 */


#include <mex.h>	/* include MEX library for Matlab */

#define printf  mexPrintf	/* Addapt looktxt.c code to Mex syntax */
#define malloc  mxMalloc
#define realloc mxRealloc
#define calloc  mxCalloc
#define VERSION "Looktxt 1.1 (MeX) $Revision: 1.8 $"
/* #define free mxFree  */
#define free NoOp

#define main looktxt		/* change stand-alone source to a Matlab usable library */
#define print_stderr mexPrintf
#define argc carg
#define argv varg
#define exit(ret) { char msg[1024]; sprintf(msg, "Looktxt/mex exited with code %i\n", ret); if (ret) mexErrMsgTxt(msg); }
#define TEXMEX


int NoOp(char *pointer)
{
  return 0;
}

#include "looktxt.c" 		/* makes all the job */

void mexFunction(int nlhs, mxArray *plhs[],
                        int nrhs, const mxArray *prhs[])
{
  char *varg[MAX_LENGTH];   /* build a fake argv array */
  int  carg           = 1;
  int i               = 0;  /* rhs argument index */
  long NumFiles       = 0;  /* number of files processes */
  FILE *fout          = NULL;
  char has_output_file= 0;
  char has_debug_mode = 0;
  char *tempname      = NULL;

  /* set program name argv[0] */
  varg[0] = (char*)mxMalloc(MAX_LENGTH);
  if (!varg[0]) {
    mexPrintf("looktxt/mex : argument %i. Size %i\n", 0, MAX_LENGTH);
    mexErrMsgTxt("looktxt/mex : can not allocate program name string argv[0].\n");
  }
  strcpy(varg[0],"looktxt");

  /* check in/out parameters */
  if (nlhs > 1)
	  mexErrMsgTxt("looktxt : Too many output arguments (1 max).");

  /* allocate memory and parse input arguments (tokens from string arguments) */
  for (i = 0; i < nrhs; i++)
  {
    long  buflen     = 0;       /* number of arguments for call */
    char *InputString= NULL;    /* input string which is then split into arguments */
    int   status     = 0;       /* flag for getting input string */
    char  EndFlag    = 0;       /* flag to exit while loop */
    char *StartLexeme= NULL;
    char *EndLexeme  = NULL;
    char *EndString  = NULL;
    char  lexeme[MAX_LENGTH];   /* current argument, stored in varg */

    if (mxIsChar(prhs[i]) != 1) /* must be a char */
    {
      mexPrintf("looktxt/mex : argument %i\n", i);
      mexErrMsgTxt("looktxt/mex : Input should be strings.\n");
    }
    /* read the input argument and return it as a string */
    InputString = mxArrayToString(prhs[i]);
    if (!InputString) {
      mexPrintf("looktxt/mex : argument %i. InputString=NULL\n", i);
      mexErrMsgTxt("looktxt/mex : can not get input parameter.\n");
    }

    /* cut input string into separated arguments for main(argc,argv) syntax */
    EndFlag = 0;
    StartLexeme = InputString;
    EndString   = InputString+strlen(InputString);
    
    /* the first input argument should be used as is (filename) */
    if (i==0) {
      varg[carg] = (char*)mxMalloc(strlen(InputString)+64);
      if (varg[carg] == NULL) {
        mexPrintf("looktxt/mex : argument %i. Size %i\n", carg, strlen(InputString));
        mexErrMsgTxt("looktxt/mex : can not allocate memory for input argument string.\n");
      }
      strcpy(varg[carg], InputString);
      carg++;
    } else while ((EndFlag == 0) && (carg < MAX_LENGTH-1)) /* read tokens iteratively */
    {
      /* search for begining of a word */
      if (*StartLexeme == ' ')
        while (*StartLexeme == ' ' && StartLexeme < EndString) {
         /* look for first non ' ' : StartLexeme */
          StartLexeme++;  /* pass all spaces */
        }
      /* now StartLexeme points on a non space or is at end */
      if (*StartLexeme == '\0' || StartLexeme >= EndString) EndFlag = 1; /* end of file reached, no other word to read */
      else {
        /* look for position of end of word (first next ' ') : EndLexeme */
        EndLexeme = strchr(StartLexeme+1, ' ');

        if (EndLexeme == NULL)
          EndLexeme = EndString;

        if (EndLexeme - StartLexeme > 0) {
          /* copy this word as a 'lexeme' element */
          strncpy(lexeme, StartLexeme, EndLexeme - StartLexeme+1);
          lexeme[EndLexeme - StartLexeme] = '\0';
          StartLexeme = EndLexeme+1;  /* will continue with next word following */
        }
      } /* else */

      if (strlen(lexeme) != 0 && lexeme != NULL && EndFlag == 0)
      {
        /* transfer the word into allocated varg[] */
        varg[carg] = (char*)mxMalloc(strlen(lexeme)+64);
        if (varg[carg] == NULL) {
          mexPrintf("looktxt/mex : argument %i. Size %i\n", carg, strlen(lexeme));
          mexErrMsgTxt("looktxt/mex : can not allocate memory for input argument string.\n");
        }
        strcpy(varg[carg], lexeme);
        
        /* test if the word is something special: output file manually set and debug mode */
        if (!strncmp(lexeme,"-o",2) || !strncmp(lexeme,"--outfile", 9))
          has_output_file=1;
        if (!strncmp(lexeme,"--debug", 7))
          has_debug_mode=1;
        carg++;
      }
      else
        EndFlag = 1; /* invalid word found: we end the search for tokens */
    } /* while */
    mxFree(InputString);

  } /* end for nrhs (all input string arguments) */
  
  /* set temporary output file (will be removed at end of MeX): must make sure it is unique */
  if (!has_output_file) {
    char filename[256];
    char *ret=NULL;
    long long_r;
    /* handle path name: local or pointing to /tmp */
#ifndef P_tmpdir
#ifdef _P_tmpdir
#define P_tmpdir _P_tmpdir
#else
#define P_tmpdir "."
#endif   /* def _P_tmpdir */
#endif   /* ndef P_tmpdir */
    /* create a unique temporary file of length < 32 using random number */
    long_r = rand();
    if (!strcmp(P_tmpdir, LK_PATHSEP_S) || !strcmp(P_tmpdir, "."))
        sprintf(filename, "lk_%li_XXXXXX", long_r);
    else
        sprintf(filename, "%s%clk_%li_XXXXXX", P_tmpdir, LK_PATHSEP_C, long_r);
    mktemp(filename);
    /* allocate it as a new argument to be passed to looktxt */
    varg[carg] = (char*)mxMalloc(strlen(filename)+64);
    if (varg[carg] == NULL) {
      mexPrintf("looktxt/mex : argument %i. Size %i\n", carg, strlen(filename));
      mexErrMsgTxt("looktxt/mex : can not allocate memory for output temporary file name.\n");
    }
    sprintf(varg[carg], "--outfile=%s", filename);
    carg++;
  } /* if (!has_output_file) */
  
  /* display command line */
  for (i=0; i<carg; i++)
  	mexPrintf("%s ", varg[i]);
  mexPrintf("\n");

  /* call looktxt, using argv set from parsing arguments */
  NumFiles = looktxt(carg,varg);  /* call main routine, return -1 in case of error, or the number of processes files */
  if (has_debug_mode)
    mexPrintf("looktxt: Processed %i file%s.\n", NumFiles, NumFiles > 1 ? "s" : "");
  /* the output file names are stored in the global options.files_to_convert_Target array of char* */

  if (NumFiles > 0) {
    plhs[0]=mxCreateCellMatrix(1,NumFiles);
    if (has_debug_mode)
      mexPrintf("looktxt: Created %i matri%s for Matlab function output.\n", NumFiles, NumFiles > 1 ? "ces" : "x");
  } else plhs[0]=mxCreateCellMatrix(1,0);
  
  mexSetTrapFlag(1);  /* On error, control returns to MeX file when using mexCallMATLAB */

  if (NumFiles > 0)
  for (i=0; i< NumFiles; i++) { /* loop on output files created by looktxt */
    /* get output file name from global array TexMex_Target_Array[] */
    /* mfile=file.Source and func=file.RootName */
    struct file_struct output_file=TexMex_Target_Array[i];

    /* evaluate file content when it is defined */
    if (output_file.TargetTxt && strlen(output_file.TargetTxt)) {
      mxArray *CellElement[1];
      struct fileparts_struct parts;
      char addpath[MAX_LENGTH];
      /* first method: evaluate function and store its result as a celle element in argout[file_index] */
      if (has_debug_mode)
        mexPrintf("looktxt: init matrix[%i] to store '%s'...\n", i, output_file.TargetTxt);
      
      parts = fileparts(output_file.TargetTxt);
      /* make sure that new function created by looktxt can be accessed by Matlab */
      sprintf(addpath, "addpath('%s')", parts.Path);
      mexEvalString(addpath);
      /* evaluate function and set cell element into nargout[file_index] */
      if (!mexCallMATLAB(1, CellElement, 0, NULL, output_file.RootName)) {
        if (has_debug_mode)
          mexPrintf("looktxt: Method 1: set matrix[%i] as '%s' content...\n", i, output_file.TargetTxt);
        mxSetCell(plhs[0], i, CellElement[0]);
      } else {
        /* direct function evaluation failed: 
           we read the function script (.m) and store it as a string into nargout[file_index] */
        long filesize=0;
        char *filestr=NULL;
        struct stat stfile;
        int status;

        /* read content */
        if (has_debug_mode)
          mexPrintf("looktxt: Method 2: checking existence for matrix[%i] file '%s'...\n", i, output_file.TargetTxt);
        
        status = stat(output_file.TargetTxt, &stfile);
        if (status) {
          mexPrintf("looktxt/mex: Warning : unable to access file '%s'\n", output_file.TargetTxt);
        } else { /* file exists */
          filesize = stfile.st_size;
          filestr  = (char*)mxMalloc(filesize+64);
          if (!filestr) {
            mexPrintf("looktxt/mex: Warning : unable to malloc %ld bytes for file '%s'\n", filesize, output_file.TargetTxt);
          } else { /* allocated memory for file content */
            strcpy(filestr, "");
            fout = fopen(output_file.TargetTxt, "r");
            if (!fout) {
              mexPrintf("looktxt/mex: Warning : unable to open file '%s' for reading\n", output_file.TargetTxt);
            } else {  /*  open and read file content */
              if (has_debug_mode)
                mexPrintf("looktxt: Method 2: read matrix[%i] content from file '%s'...\n", filesize, output_file.TargetTxt);
              if (!fread(filestr, 1, filesize, fout)) {
                mexPrintf("looktxt/mex: Warning : unable to read file '%s'\n", output_file.TargetTxt);
              }
              fclose(fout);
            } /* if fopen ok */
          } /* if malloc file string ok */
        } /* if stat ok */
        if (has_debug_mode)
          mexPrintf("looktxt: Method 2:  storing matrix[%i] string as content from file '%s'...\n", 
            i, output_file.TargetTxt);
        /* could we extract file content ? */
        if (filesize && filestr && strlen(filestr)) {
          mxArray *FileString = mxCreateString(filestr);
          if (FileString == NULL) {
            mexPrintf("looktxt/mex : can not create Matlab string[%i] from  %s (length %i).\n", 
              i, output_file.TargetTxt, strlen(filestr));
            mexErrMsgTxt("looktxt/mex : can not allocate memory for output string.\n");
          }
          /* set cell element and free file content copy */
          mxSetCell(plhs[0], i, FileString);
          mexPrintf("Method 2:  transfered %i matrix %s\n", filesize, output_file.TargetTxt);
          mxFree(filestr);
        }
      } /* else (if mexCall) */
      
      /* clean up temporary files if necessary (created when no --output is present in input arguments) */
      if (!has_output_file) {
        if (has_debug_mode)
          mexPrintf("looktxt: remove matrix[%i] file '%s'.\n", i, output_file.TargetTxt);
      	remove(output_file.TargetTxt);
      	if (output_file.TargetBin && strlen(output_file.TargetBin))
      		remove(output_file.TargetBin);
      }
    } /* if output_file */
    else if (NumFiles == 1) NumFiles=0;
  } /* for i */
  
  if(NumFiles == 1) {
    mxArray *ptr;
    if (has_debug_mode)
      mexPrintf("looktxt: only a single matrix[%i].\n", NumFiles);
    
    ptr= mxGetCell(plhs[0], 0);
    plhs[0] = mxDuplicateArray(ptr);
    if (has_debug_mode)
      mexPrintf("looktxt: set out matrix[%i] (not using cell array).\n", NumFiles);
  }

  for (i=0; i < carg; mxFree(varg[i++]));
  if (has_debug_mode)
      mexPrintf("looktxt: done %i files.\n", NumFiles);
}
