/**********************************************************************************************
***********************************************************************************************
Macro: MULTIPLE_GENMOD
Created Date/Author: Oct. 26, 2012/Dana Nickleach 
Last Update Date/Person: Dec. 17, 2013/Dana Nickleach 
Current Version: V8
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu 

Purpose:  To produce a summary table in a Microsoft Word Document from a multivariable 
generalized linear regression model using the log or logit link function with a binomial 
distribution; a log link with a Poisson or negative binomial distribution; or an identity link
with a normal distribution.  Relative risk is reported for a binomial distribution using a log
link, odds ratio for a binomial distribution using a logit link; rate ratio for Poisson or
negative binomial models; and slope for normal models.  The modeling must be done before using
this macro.  The macro can also handle GEE models that were specified using the REPEATED 
statement.
 
Notes: The regression must be conducted using PROC GENMOD before calling this macro.  Four 
tables must be created using ODS output when running PROC GENMOD (modelinf, nobs, type3, and 
estimate) using the ODS statement "ODS OUTPUT ModelInfo=modelinf NObs = nobs Type3 = type3 
ParameterEstimates = estimate".  Note that if a REPEATED statement was used then the estimate 
table needs to be created using “GEEEmpPEst = estimate”.  The option TYPE3 must also be 
specified in the MODEL statement.  Variable names must not be more than 20 characters.  This 
macro is not set up to handle models including interaction terms and will produce incorrect 
results.

Parameters: 

OUTPATH        File path for output table to be stored.  
FNAME          File name for output table.  
FOOTNOTE       Text of footnote to include in the table.  It should be in quotes.  This 
               footnote will appear below the footnote containing the number of observations.  
               Leave this field blank if not including an additional footnote. 
TYPE3          Set to F to suppress type III p-values from being reported in the table.  The 
               default value is T.
CONTRAST       Set to T if reporting contrast estimates from the ESTIMATE statement as well.  
               The contrast must be specified using an ESTIMATE statement in PROC
               GENMOD.  The EXP option should NOT be used.  An ODS OUTPUT data set must be 
               created using Estimates=contrast.  The default value is F. 
DEBUG          Set to T to run in debug mode.  Work datasets will not be deleted in debug 
               mode.  This is useful if you are editing the code or want to further manipulate 
               the resulting data sets.  The default value is F.

***********************************************************************************************
* For more details, please see the related documentation.
**********************************************************************************************/

%macro MULTIPLE_GENMOD(OUTPATH=,FNAME=,FOOTNOTE=,TYPE3=T,CONTRAST=F,DEBUG=F);

   %local classexist dsn dist link __Macro_Err length;

   %let debug = %UPCASE(&debug);
   %let contrast = %UPCASE(&contrast);
   %let type3 = %UPCASE(&type3);

   /* Initialize error variable */
   %let __Macro_Err = 0;

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

   data _null_;
      dsid=open('estimate');
      /* Check to see if class variables were used - s/b Level1 variable */
      check=varnum(dsid,'Level1');
      /* 0 = does not exist, 1 = exists */
      exist = MIN(check, 1);
      CALL SYMPUT('classExist',PUT(exist,1.));

      /* Check to see if parameter variable is called "parm" or "parameter" */
      check=varnum(dsid,'Parameter');
      /* 0 = does not exist, 1 = exists */
      exist = MIN(check, 1);
      /* Determine which variable name nees to be used */
      if exist = 1 then pname = 'Parameter';
      else pname = 'Parm';
      CALL SYMPUT('parameter',pname);

      /* Check to see names used for CI */
      check = varnum(dsid,'LowerWaldCL');
      /* 0 = does not exist, 1 = exists */
      exist = MIN(check, 1);
      /* Determine which variable name nees to be used */
      if exist = 1 then do;
         LCL = 'LowerWaldCL';
         UCL = 'UpperWaldCL';
      end;
      else do;
         LCL = 'LowerCL';
         UCL = 'UpperCL';
      end;
      CALL SYMPUT('LCL',LCL);
      CALL SYMPUT('UCL',UCL);

      /* Check to see names used for p-value */
      check = varnum(dsid,'ProbChiSq');
      /* 0 = does not exist, 1 = exists */
      exist = MIN(check, 1);
      /* Determine which variable name nees to be used */
      if exist = 1 then pval = 'ProbChiSq';
      else pval = 'probz';
      CALL SYMPUT('pval',pval);
       
   run;

   /* Calculate number of observations read and used */
   PROC SQL noprint;
      select NObsRead, NObsUsed 
      into :nobs_read, :nobs_used
      from nobs (obs=1);
   QUIT;

   /* Get distribution used and dependent variable and dataset name */
   PROC SQL noprint;
      select cvalue1 into: dist
      from modelinf
      where label1 = 'Distribution';
      select cvalue1 into: outName
      from modelinf
      where label1 = 'Dependent Variable';
      select cvalue1 into: dsn
      from modelinf
      where label1 = 'Data Set';
      /* Get link function */
      select cvalue1 into: link
      from modelinf
      where label1 = 'Link Function';
   QUIT;

   /* Make sure that link and distribution are expected */
   %if &dist ~= Normal AND &dist ~= Poisson AND &dist ~= Negative Binomial AND &dist ~= Binomial %then %do;
      %put ERROR: Distribution must be normal, Poisson, negative binomial, or binomial.; 
      %let __Macro_Err=1;
   %end;
   %else %if (&dist = Normal and &link ~= Identity) OR 
     ((&dist = Negative Binomial OR &dist = Poisson) AND &link ~= Log) OR
     (&dist = Binomial AND (&link ~= Log AND &link ~= Logit)) %then %do; 
      %put ERROR: Link and distribution combination cannot be handled by this macro.; 
      %let __Macro_Err=1;
   %end;

   /* If there is an error in the parameters supplied then exit */
   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   /* Get outcome label */
   DATA _NULL_;
      set &dsn (obs=1);
      call symput('outLab', VLABEL(&outname));
   RUN;

   /* Get variable labels */
   PROC CONTENTS DATA = &dsn out=_cont (rename=(name=&parameter)) noprint;
   RUN;

   /* Prepare for merge */
   data _typ3;
       length &parameter $32;
      set type3 (rename=(ProbChiSq = p_type3 Source=&parameter)); 
   run;

   /* If reporting contrast estimates as well */
   %if &contrast = T %then %do;
      
      /* Create break row with header */     
      DATA _temp;
          &parameter = 'Contrasts:';
      RUN;

      /* Concatenate on contrasts */
      PROC SQL;
         create table _estimate
         as select * from estimate
         outer union CORRESPONDING
         select * from _temp
         outer union CORRESPONDING
         select * from contrast (KEEP=label LBetaEstimate LBetaLowerCL LBetaUpperCL 
          ProbChiSq rename=(label = &parameter LBetaEstimate = estimate LBetaLowerCL=&LCL 
          LBetaUpperCL=&UCL ProbChiSq=&pval));
      QUIT;

   %end;
   %else %do;
      DATA _estimate;
         set estimate;
      RUN;
   %end;

   /* Preserve original variable order */
   DATA _est2;
      set _estimate;
      order = _n_;
      pval = &pval;
      /* Calculate RR and CI */
      if estimate ~= 0 then do;
         /* For distributions that use log link - take exp of estimates */
         %if &dist = Negative Binomial OR &dist = Poisson OR &dist = Binomial %then %do;
            RR = exp(estimate);
            LCL = exp(&LCL);
            UCL = exp(&UCL);
         %end;
         /* For normal linear regression don't take exp */
         %else %if &dist = Normal %then %do;
            RR = estimate;
            LCL = &LCL;
            UCL = &UCL;
         %end;
      end;
      if &parameter in ('Intercept' 'Scale' 'Dispersion') then delete;

   RUN;

   PROC SORT DATA = _est2;
      by &parameter;
   RUN;
   proc sort data=_typ3;
      by &parameter;
   run;
   PROC SORT DATA = _cont;
      by &parameter;
   RUN;

   /* Get length of parameters in estimate file */
   DATA _NULL_;
      set _est2 (obs=1);
      CALL SYMPUT('length',PUT(VLENGTH(&parameter),best12.));
   RUN;

   /* Length should be at least 32 as it is in the other files that will be merged on */
   %let length = %sysfunc(MAX(&length,32));

   data _MG; 
      length plabel $256 &parameter $&length;
      merge _typ3 (in=a keep=&parameter p_type3) _est2 (in=b) _cont (keep=&parameter label); 
      by &parameter;
       if a or b;
       if label ~= ' ' then plabel = label;
       else plabel = &parameter;

       /* If no class variables then level is blank */
       %if &classexist = 0 %then %do;
          level1 = ' ';
       %end;
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
      COLUMN plabel Level1 ("&outLab" '--------------------------------------------------------'
          (RR LCL UCL pval %if &type3 = T %then %do; p_type3 %end;)); 
      DEFINE plabel/ order order=data "Covariate"  STYLE(COLUMN) = {JUST = L}; 
      DEFINE Level1/ DISPLAY "Level" STYLE(COLUMN) = {JUST = L CellWidth=20%}; 
      DEFINE RR/DISPLAY STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=10.2
          %if &dist = Negative Binomial OR &dist = Poisson %then %do; "Rate Ratio" %end;
          %else %if &dist = Binomial and &link = Logit %then %do; "Odds Ratio" %end;
          %else %if &dist = Binomial and &link = Log %then %do; "Relative Risk" %end;
          %else %if &dist = Normal %then %do; "B" %end;; 
      DEFINE LCL /DISPLAY '95%CI Low' STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=10.2;
      DEFINE UCL /DISPLAY '95%CI Up' STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=10.2;
      DEFINE pval/DISPLAY STYLE(COLUMN)= {JUST = C CellWidth=8%} FORMAT=PVALUE8.3
         %if &dist = Normal %then %do; "B P-Value" %end;
         %else %if &link = Log %then %do; "RR P-value" %end;
         %else %if &link = Logit %then %do; "OR P-value" %end;; 

      COMPUTE pval; 
         IF pval <0.05 THEN CALL DEFINE("pval", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
      ENDCOMP; 

      %if &type3 = T %then %do;
         DEFINE p_type3/ORDER "Type3 P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%} 
               FORMAT=PVALUE8.3 MISSING;
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
         DELETE _est2 _typ3 _mg _cont _options _temp;
      QUIT;
   %end;

%mend MULTIPLE_GENMOD;


