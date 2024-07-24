%let pgm=utl-create-and-invert-a-two-dimensional-numeric-array-from-a-sas-dataset-r-sas-fcmp;

Create  and invert a two dimensional numeric array from a sas dataset

Useful for in memory processing of small sas datasets or prototyping sas code for conversion to matrix languages

github
https://tinyurl.com/2jzh9epb
https://github.com/rogerjdeangelis/utl-create-and-invert-a-two-dimensional-numeric-array-from-a-sas-dataset-r-sas-fcmp

The limit to the size of the 'put (all)' is 32767 bytes.
You can probably extend the length to the max size for a macro variable of 65533.

    Three

         1. build array statement
         2. sas invert datastep array usin fcmp
         3. r invert datastep array
         4. utl_numary macro
         5. related repos
/*               _     _
 _ __  _ __ ___ | |__ | | ___ _ __ ___
| `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \
| |_) | | | (_) | |_) | |  __/ | | | | |
| .__/|_|  \___/|_.__/|_|\___|_| |_| |_|
|_|
*/

/**************************************************************************************************************************/
/*                                                          |                                                             */
/*                     INPUT (HAVE.SAS7BDAT)                |  OUTPUT(CREATE THIS SAS ARRAY STATEMENT) 100X30 ARRAY       */
/*    ===================================================   |  ====================================================       */
/*                                                          |                                                             */
/*    Obs  X1  X2  X3  X4  X5 ... X26  X27  X28  X29  X30   |  array ary[100,30] xy1-xy30                                 */
/*                                                          |    (                                                        */
/*      1   1   2   3   4   5 ...  26   27   28   29   30   |     1   2   3   4   5 ...  26   27   28   29   30           */
/*      2   1   2   3   4   5 ...  26   27   28   29   30   |     1   2   3   4   5 ...  26   27   28   29   30           */
/*      3   1   2   3   4   5 ...  26   27   28   29   30   |     1   2   3   4   5 ...  26   27   28   29   30           */
/*     ..                                                   |     ...                                                     */
/*     98   1   2   3   4   5 ...  26   27   28   29   30   |     1   2   3   4   5 ...  26   27   28   29   30           */
/*     99   1   2   3   4   5 ...  26   27   28   29   30   |     1   2   3   4   5 ...  26   27   28   29   30           */
/*    100   1   2   3   4   5 ...  26   27   28   29   30   |     1   2   3   4   5 ...  26   27   28   29   30           */
/*                                                          |     )                                                       */
/*                                                          |                                                             */
/*                                                          |   THE SUM THE DIAGONAL                                      */
/*                                                          |                                                             */
/*                                                          |     DIAGONAL_SUM=46500                                      */
/*                                                          |                                                             */
/*----------------------------------------------+-----------+-----------------------------+-------------------------------*/
/*                                              |                                         |                               */
/* SAS INVERT MATRIX USING FCMP                 |                                         |                               */
/* ============================                 |                                         |                               */
/*                                              |                                         |                               */
/*          INPUT                               |          PROCESS                        |   OUTPUT                      */
/*          =====                               |         =========                       |   =======                     */
/*                                              |                                         |                               */
/*  No external input needed just               |  options cmplib=work.functions;         |   INVERSE                     */
/*  generate and add                            |  proc fcmp outlib=work.functions.math;  |                               */
/*   "array A[3,3] xy1-xy9"                     |    subroutine invert_matrix(            |   7 -3 -3                     */
/*                                              |           A[*,*]                        |   -1 1 0                      */
/*  array A[3,3] xy1-xy9                        |          ,n                             |   -1 0 1                      */
/*           (1,3,3,1,4,3,1,3,4)                |          ,Ainv[*,*]);                   |                               */
/*                                              |      outargs Ainv;                      |                               */
/*                                              |      call inv(A, Ainv);                 |                               */
/*  options validvarname=upcase;                |   endsub;                               |                               */
/*  libname sd1 "d:/sd1";                       |   run;                                  |                               */
/*  data sd1.have;                              |                                         |                               */
/*    input xy1 xy2 xy3;                        |   data _null_;                          |                               */
/*  cards4;                                     |     array A[3,3] xy1-xy9                |                               */
/*  1 3 3                                       |          (%utl_numary(sd1.have));       |                               */
/*  1 4 3                                       |     array Ainv[3,3];                    |                               */
/*  1 3 4                                       |     call invert_matrix(A, 3, Ainv);     |                               */
/*  ;;;;                                        |     put "INVERSE" /;                    |                               */
/*  run;quit;                                   |     do i = 1 to 3;                      |                               */
/*                                              |      put Ainv[i,1] Ainv[i,2] Ainv[i,3]; |                               */
/*  %put "%utl_numary(sd1.have)";               |     end;                                |                               */
/*                                              |   run;quit;                             |                               */
/*    "1,3,3,1,4,3,1,3,4"                       |                                         |                               */
/*                                              |                                         |                               */
/*----------------------------------------------|-----------------------------------------|-------------------------------*/
/*                                              |                                         |                               */
/*  R INVERT MATRIX USING FCMP                  |                                         |                               */
/*  ==========================                  |                                         |                               */
/*                                              |                                         |                               */
/*  %put "ONLY INPUT IS %utl_numary(sd1.have)"; |  %utl_rbeginx;                          |  INVERTED                     *
/*                                              |  parmcards4;                            |                               */
/*  "ONLY INPUT IS 1,3,3,1,4,3,1,3,4"           |  # Define a sample matrix               |        [,1] [,2] [,3]         */
/*                                              |  A<-matrix(                             |   [1,]    7   -3   -3         */
/*  Same as SAS FCMP but no need for ARRAY      |     c(%utl_numary(sd1.have))            |   [2,]   -1    1    0         */
/*                                              |    ,nrow = 3                            |   [3,]   -1    0    1         */
/*                                              |    ,ncol = 3)                           |                               */
/*                                              |  A_inv <- t(solve(A))                   |                               */
/*                                              |  options(digits = 1)                    |                               */
/*                                              |  cat("INVERTED \n")                     |                               */
/*                                              |  print(round(A_inv))                    |                               */
/*                                              |  ;;;;                                   |                               */
/*                                              |  %utl_rendx;                            |                               */
/*                                              |                                         |                               */
/**************************************************************************************************************************/

/*   _           _ _     _
/ | | |__  _   _(_) | __| |   __ _ _ __ _ __ __ _ _   _
| | | `_ \| | | | | |/ _` |  / _` | `__| `__/ _` | | | |
| | | |_) | |_| | | | (_| | | (_| | |  | | | (_| | |_| |
|_| |_.__/ \__,_|_|_|\__,_|  \__,_|_|  |_|  \__,_|\__, |
                                                  |___/
*/

/*
 _ __  _ __ ___   ___ ___  ___ ___
| `_ \| `__/ _ \ / __/ _ \/ __/ __|
| |_) | | | (_) | (_|  __/\__ \__ \
| .__/|_|  \___/ \___\___||___/___/
|_|
*/

%macro utl_numary(_inp);

 %symdel _array / nowarn;

 %dosubl(%nrstr(
     filename clp clipbrd lrecl=64000;
     data _null_;
       file clp;
       set &_inp(keep=_numeric_);
       put (_all_) ($) @@;
     run;quit;
     data _null_;
       infile clp;
       input;
       _infile_=compbl(_infile_);
       _infile_=translate(strip(_infile_),',',' ');
       call symputx('_array',_infile_);
     run;quit;
     ))

     &_array

%mend utl_numary;

%put "array ary[100,30] col1-col3000 (%utl_numary(have))";


%symdel _ary / nowarn;

data _null_;
  retain diagonal_sum 0;
  array ary[100,30] xy1-xy3000 (%utl_numary(have));
  set have end=dne;
  do row=1 to 30;
    diagonal_sum = diagonal_sum + ary[row,row];
  end;
  if dne then do;
     put diagonal_sum=;
     * using Gauss formula (n+1)*n/2 => repeated 10 times;
     check=100*31*30/2;
     put check;
  end;
run;quit;


/**************************************************************************************************************************/
/*                                                                                                                        */
/*    MACRO CREATES THIS STATEMENT                                                                                        */
/*                                                                                                                        */
/*      array ary[100,30] col1-col3000 (                                                                                  */
/*       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1, ...                          */
/*       ,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6 ...                          */
/*       10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10 ...                          */
/*       13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13 ...                          */
/*       16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 ...                          */
/*       19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 ...                          */
/*       22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 ...                          */
/*       25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 ...                          */
/*       28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 ...                          */
/*       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1, ...                          */
/*       ,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6 ...                          */
/*       10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10 ...                          */
/*       13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13 ...                          */
/*       16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 ...                          */
/*       19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 ...                          */
/*       22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 ...                          */
/*       25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 ...                          */
/*       28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 ...                          */
/*       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1, ...                          */
/*       ,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6 ...                          */
/*       10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10 ...                          */
/*       13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13 ...                          */
/*       16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 ...                          */
/*       19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 ...                          */
/*       22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 ...                          */
/*       25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 ...                          */
/*       28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 ...                          */
/*       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1, ...                          */
/*       ,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6 ...                          */
/*       10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10 ...                          */
/*       13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13 ...                          */
/*       16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 ...                          */
/*       19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 ...                          */
/*       22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 ...                          */
/*       25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 ...                          */
/*       28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 ...                          */
/*       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1, ...                          */
/*       ,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6 ...                          */
/*       10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10 ...                          */
/*       13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13 ...                          */
/*       16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 ...                          */
/*       19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 ...                          */
/*       22,23,24,25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 ...                          */
/*       25,26,27,28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 ...                          */
/*       28,29,30,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 ...                          */
/*       1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1, ...                          */
/*       ,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,1,2,3,4,5,6 ...                          */
/*       10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)                                                  */
/*                                                                                                                        */
/*                                                                                                                        */
/*  OUTPUT(CREATE THIS SAS ARRAY STATEMENT)                                                                               */
/*  =======================================                                                                               */
/*                                                                                                                        */
/*  array ary[100,30] xy1-3000                                                                                            */
/*    (                                                                                                                   */
/*     1   2   3   4   5 ...  26   27   28   29   30                                                                      */
/*     1   2   3   4   5 ...  26   27   28   29   30                                                                      */
/*     1   2   3   4   5 ...  26   27   28   29   30                                                                      */
/*                                                                                                                        */
/*     1   2   3   4   5 ...  26   27   28   29   30                                                                      */
/*     1   2   3   4   5 ...  26   27   28   29   30                                                                      */
/*     1   2   3   4   5 ...  26   27   28   29   30                                                                      */
/*     )                                                                                                                  */
/*                                                                                                                        */
/*   THE SUM THE DIAGONAL                                                                                                 */
/*                                                                                                                        */
/*     DIAGONAL_SUM=46500                                                                                                 */
/*                                                                                                                        */
/*
/*************************************************************************************************************************

/*___                     __                        _                     _
|___ \   ___  __ _ ___   / _| ___ _ __ ___  _ __   (_)_ ____   _____ _ __| |_
  __) | / __|/ _` / __| | |_ / __| `_ ` _ \| `_ \  | | `_ \ \ / / _ \ `__| __|
 / __/  \__ \ (_| \__ \ |  _| (__| | | | | | |_) | | | | | \ V /  __/ |  | |_
|_____| |___/\__,_|___/ |_|  \___|_| |_| |_| .__/  |_|_| |_|\_/ \___|_|   \__|
                                           |_|
*/

options validvarname=upcase;
libname sd1 "d:/sd1";
data sd1.have;
  input xy1 xy2 xy3;
cards4;
1 3 3
1 4 3
1 3 4
;;;;
run;quit;

%put "%utl_numary(sd1.have)";

options cmplib=work.functions;
proc fcmp outlib=work.functions.math;
  subroutine invert_matrix(
         A[*,*]
        ,n
        ,Ainv[*,*]);
    outargs Ainv;
    call inv(A, Ainv);
  endsub;
run;

data _null_;
  array A[3,3] xy1-xy9
       (%utl_numary(sd1.have));
  array Ainv[3,3];
  call invert_matrix(A, 3, Ainv);
  put "INVERSE" /;
  do i = 1 to 3;
   put
    Ainv[i,1]  Ainv[i,2]  Ainv[i,3];
  end;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  INVERSE                                                                                                               */
/*                                                                                                                        */
/*  7 -3 -3                                                                                                               */
/*  -1 1 0                                                                                                                */
/*  -1 0 1                                                                                                                */
/*                                                                                                                        */
/**************************************************************************************************************************/

* SAME INPUT AS FCMP ;

%utl_rbeginx;
parmcards4;
library(dplyr)
library(haven)
# Define a sample matrix
A<-matrix(
   c(%utl_numary(sd1.have))
  ,nrow = 3
  ,ncol = 3)
A_inv <- solve(A)
options(digits = 1)
cat("\nInverted matrix A^(-1):\n")
print(round(A_inv))
;;;;
%utl_rendx;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  INVERTED                                                                                                              */
/*                                                                                                                        */
/*        [,1] [,2] [,3]                                                                                                  */
/*   [1,]    7   -3   -3                                                                                                  */
/*   [2,]   -1    1    0                                                                                                  */
/*   [3,]   -1    0    1                                                                                                  */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*  _           _   _
| || |    _   _| |_| |    _ __  _   _ _ __ ___   __ _ _ __ _   _
| || |_  | | | | __| |   | `_ \| | | | `_ ` _ \ / _` | `__| | | |
|__   _| | |_| | |_| |   | | | | |_| | | | | | | (_| | |  | |_| |
   |_|    \__,_|\__|_|___|_| |_|\__,_|_| |_| |_|\__,_|_|   \__, |
                    |_____|                                |___/
*/

filename ft15f001 "c:/oto/utl_numary.sas";
parmcards4;
%macro utl_numary(_inp);

 %symdel _array / nowarn;

 %dosubl(%nrstr(
     filename clp clipbrd lrecl=64000;
     data _null_;
       file clp;
       set &_inp(keep=_numeric_);
       put (_all_) ($) @@;
     run;quit;
     data _null_;
       infile clp;
       input;
       _infile_=compbl(_infile_);
       _infile_=translate(strip(_infile_),',',' ');
       call symputx('_array',_infile_);
     run;quit;
     ))

     &_array

%mend utl_numary;
;;;;
run;quit;

%put "%utl_numary(sd1.have)";

/*___             _       _           _
| ___|   _ __ ___| | __ _| |_ ___  __| |  _ __ ___ _ __   ___  ___
|___ \  | `__/ _ \ |/ _` | __/ _ \/ _` | | `__/ _ \ `_ \ / _ \/ __|
 ___) | | | |  __/ | (_| | ||  __/ (_| | | | |  __/ |_) | (_) \__ \
|____/  |_|  \___|_|\__,_|\__\___|\__,_| |_|  \___| .__/ \___/|___/
                                                  |_|
*/

REPO
-----------------------------------------------------------------------------------------------------------------------------------------
https://github.com/rogerjdeangelis/utl-attempt-to-call-fcmp-containing-pokelong-from-sql
https://github.com/rogerjdeangelis/utl-calculate-regression-coeficients-in-base-sas-fcmp-proc-reg-r-and-python
https://github.com/rogerjdeangelis/utl-datastep-in-memory-matrix-and-submatrix-row-and-column-reductions-summations-fcmp
https://github.com/rogerjdeangelis/utl-extracting-sas-meta-data-using-sas-macro-fcmp-and-dosubl
https://github.com/rogerjdeangelis/utl-flip-a-file-upside-down-using-recent-fcmp-dynamic-array-by-Bart-Jablonski
https://github.com/rogerjdeangelis/utl-implementation-and-benchmarks-for-heapsort-mergesort-quicksort-sortn-fcmpsort
https://github.com/rogerjdeangelis/utl-many-interfaces-to-python-open-code-and-within-datastep-fcmp-wps
https://github.com/rogerjdeangelis/utl-numeric-to-numeric-formats-using-fcmp-and-proc-format
https://github.com/rogerjdeangelis/utl-proof-of-concept-using-dosubl-to-create-a-fcmp-like-function-for-a-rolling-sum-of-size-three
https://github.com/rogerjdeangelis/utl-push-and-pop-words-off-a-string-sas-fcmp
https://github.com/rogerjdeangelis/utl-round-up-to-the-nearest-mutiple-of-n-using-fcmp-code-fx-equals-ceil-x-divided-by-n-times-n
https://github.com/rogerjdeangelis/utl-sas-array-macro-fcmp-or-dosubl-take-your-choice
https://github.com/rogerjdeangelis/utl-sas-fcmp-hash-stored-programs-python-r-functions-to-find-common-words
https://github.com/rogerjdeangelis/utl-simple-datastep-replacement-for-fcmp-use-stored-program-
https://github.com/rogerjdeangelis/utl-very-fast-elegant-examples-of-recursion-using-fcmp-by-Bartosz-Jablonski
https://github.com/rogerjdeangelis/utl_sort_summarize_set_merge_using_functions_in_formats_groupformat_fcmp

/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/
