/**********************************************************************************************
***********************************************************************************************
Macro: MULTIPLE_LINREG
Created Date/Author: March 18, 2013/Dana Nickleach
Last Update Date/Person: April 17, 2013/Dana Nickleach
Current Version: V2
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu 

Purpose:  To produce a summary table in a Microsoft Word Document from a multivariable general 
linear model.  The modeling must be done before using this macro.  

Notes: The regression must be conducted using PROC GLM before calling this macro.  Three 
tables must be created using ODS output when running PROC GLM (nobs, type3, and estimate) as in
the example below.  The options SOLUTION and CLPARM must also be used in the model statement.  
Variable names must not be more than 20 characters.  This macro is not set up to handle models 
including interaction terms and will produce incorrect results.   

Parameters: 

DSN            The name of the data set to be analyzed.        
OUTPATH        File path for output table to be stored.  
FNAME          File name for output table.  
FOOTNOTE       Text of footnote to include in the table (optional).  It should be in quotes.  
               This footnote will appear below the footnote containing the number of 
               observations.  Leave this field blank if not including an additional footnote. 
TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
* For more details, please see the related documentation.
**********************************************************************************************/

%macro MULTIPLE_LINREG(DSN=,OUTPATH=,FNAME=,FOOTNOTE=,TYPE3=T,DEBUG=F);

   %let debug = %UPCASE(&debug);
   %let type3 = %UPCASE(&type3);

   %local NOBS_READ OUTLAB OUTNAME NOBS_USED;

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

   /* Calculate number of observations read and used and get dependent variable name */
   PROC SQL noprint;
      select NObsRead, NObsUsed, DependentVariables 
      into :nobs_read, :nobs_used, :outName
      from nobs (obs=1);
   QUIT;

   /* Get outcome label */
   DATA _NULL_;
      set &dsn (obs=1);
      call symput('outLab', VLABEL(&outname));
   RUN;

   /* Get variable labels */
   PROC CONTENTS DATA = &dsn out=_cont (rename=(name=parameter)) noprint;
   RUN;

   /* Prepare for merge */
   data _typ3;
       length parameter $32;
      set type3 (rename=(ProbF = p_type3 source=parameter)); 
   run;

   /* Parameter Estimates */   
   DATA _est2;
      length parameter $32;
      set estimate (rename=(parameter=temp));
       /* Preserve original variable order */
       order = _n_;
       /* Split variable name and order */
       parameter = SCAN(temp,1);
       level = LEFT(SUBSTR(temp,INDEX(temp,' ')));
       DROP temp;

       /* Recode reference levels as missing values */
       if StdErr = . then estimate = .;

       /* Drop intercept */
       if parameter = 'Intercept' then delete;
   RUN;

   PROC SORT DATA = _est2;
      by parameter;
   RUN;
   proc sort data=_typ3;
      by parameter;
   run;
   PROC SORT DATA = _cont;
      by parameter;
   RUN;

   /* Merge on type III p-values */
   data _MG; 
      length plabel $256;
      merge _typ3 (in=a keep=parameter p_type3) _est2 (in=b) _cont (keep=parameter label); 
      by parameter;
       if a or b;
       if label ~= ' ' then plabel = label;
       else plabel = parameter;
   run; 

   /* Return to original order */
   proc sort data=_mg;
      by order;
   run; 

   *---- table template -----;  
   ODS PATH WORK.TEMPLAT(UPDATE)
   SASUSR.TEMPLAT(UPDATE) SASHELP.TMPLMST(READ);

   PROC TEMPLATE;
      DEFINE STYLE STYLES.TABLES;
      NOTES "MY TABLE STYLE"; 
      PARENT=STYLES.MINIMAL;

       STYLE SYSTEMTITLE /FONT_SIZE = 12pt   FONT_FACE = "TIMES NEW ROMAN";

       STYLE HEADER /
           FONT_FACE = "TIMES NEW ROMAN"
            CELLPADDING=8
            JUST=C
            VJUST=C
            FONT_SIZE = 10pt
           FONT_WEIGHT = BOLD; 

       STYLE TABLE /
            FRAME=HSIDES            /* outside borders: void, box, above/below, vsides/hsides, lhs/rhs */
           RULES=GROUP              /* internal borders: none, all, cols, rows, groups */
           CELLPADDING=6            /* the space between table cell contents and the cell border */
            CELLSPACING=6           /* the space between table cells, allows background to show */
            JUST=C
            FONT_SIZE = 10pt
           BORDERWIDTH = 0.5pt;  /* the width of the borders and rules */

       STYLE DATAEMPHASIS /
           FONT_FACE = "TIMES NEW ROMAN"
           FONT_SIZE = 10pt
           FONT_WEIGHT = BOLD;

       STYLE DATA /
           FONT_FACE = "TIMES NEW ROMAN" 
           FONT_SIZE = 10pt;

       STYLE SYSTEMFOOTER /FONT_SIZE = 9pt FONT_FACE = "TIMES NEW ROMAN" JUST=C;
   END;

   RUN; 

   **** PRINT THE TABLE ************;
   OPTIONS ORIENTATION=PORTRAIT MISSING = "-" NODATE;
   ODS rtf STYLE=TABLES file= "&OUTPATH.&FNAME &SYSDATE..DOC"; 

   PROC REPORT DATA=_mg HEADLINE HEADSKIP CENTER SPANROWS LS=256
      STYLE(REPORT)={JUST=CENTER} SPLIT='~'   nowd; 
      COLUMN plabel Level ("&outLab" '--------------------------------------------------------'
          (estimate LowerCL UpperCL Probt %if &type3 = T %then %do; p_type3 %end;)); 
      DEFINE plabel/order order=data "Covariate" STYLE(COLUMN) = {JUST = L}; 
      DEFINE Level/DISPLAY "Level" STYLE(COLUMN) = {JUST = L CellWidth=20%}; 
      DEFINE estimate/DISPLAY "B" STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=10.2; 
       DEFINE LowerCL/DISPLAY '95%CI Low' STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=10.2;
      DEFINE UpperCL/DISPLAY '95%CI Up' STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=10.2;
       DEFINE Probt/DISPLAY "B~P-value" STYLE(COLUMN)= {JUST = C CellWidth=8%} 
          FORMAT=PVALUE8.3; 

      COMPUTE Probt; 
         IF Probt <0.05 THEN CALL DEFINE("Probt", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
      ENDCOMP; 

       %if &type3 = T %then %do; 
          DEFINE p_type3/ORDER "Type3 P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%} 
               FORMAT=PVALUE8.3; 
          COMPUTE p_type3; 
            IF p_type3 <0.05 THEN CALL DEFINE("p_type3", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         ENDCOMP;   
       %end;

       compute after plabel; 
         line '';   
      endcomp;      

       compute after _page_;
          length text1 - text2 $80.;
           text1 = "*  Number of observations in the original data set = %TRIM(&nobs_read). ";
           text2 = "   Number of observations used = %TRIM(&nobs_used).";
                   
           line @0 text1 $80.;
           line @0 text2 $80.;
           %if &footnote ~= %STR() %then %do;
              line @0 &footnote;
           %end;
      ENDCOMP; 

   RUN;

   ODS RTF CLOSE;

   /* Reload original options that were in use before running the macro */
   PROC OPTLOAD data=_options;
   RUN;

   %if &debug = F %then %do;
      /* Delete temporary datasets */
      PROC DATASETS lib=work;
         DELETE _est2 _typ3 _mg _cont _options;
      QUIT;
   %end;

%mend MULTIPLE_LINREG;


