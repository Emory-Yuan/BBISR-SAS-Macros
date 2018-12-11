/**********************************************************************************************
***********************************************************************************************
Macro: UNI_GENMOD
Created Date/Author: May 14, 2013/Dana Nickleach 
Last Update Date/Person: May 7, 2014/Dana Nickleach 
Current Version: V7
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu

Purpose:  To conduct univariate regression analysis for each variable in the dataset using a 
GEE model via PROC GENMOD.  A non-GEE model can also be used.  The normal, binomial, Poisson, 
or negative binomial distributions can be used.  An identity link will be used with a normal
distribution, a logit link with binomial, and a log link with Poisson or negative binomial.  
The parameter estimates with 95% CI, p-value, and type III p-values are presented.  For 
categorical variables, the reference group will be shown along with the number of observation 
in each category.

Notes: 
Parameters: 

DSN            The name of the data set to be analyzed.
OUTCOME        The name of the continuous or binary outcome variable.
DESC           Set to T to reverse the order of a binary outcome (optional).  The default value 
               is F. 
CLIST          List of categorical variables, separated by empty space.
CREFLIST       If you want to change the reference groups of the categorical variables then 
               repeat CLIST here and follow each variable name by the (desc) option, where 
               needed.  Separate variable terms with an asterisk instead of a space (optional).     
NLIST          List of numerical variables, separated by empty space. 
DIST           Probability distribution to use for the model (optional).  Valid values are 
               NORMAL, BINOMIAL, POISSON, and NEGBIN.  The default value is NORMAL.  
SUBJECT        Subject-effect to be specified in REPEATED statement in order to use a GEE 
               model.  Leave blank if not using a GEE model.
WITHINSUB      Within-subject-effect to be specified in REPEATED statement.  This is optional
               if using a GEE model.  Leave blank if not using a GEE model.   
TYPE           Correlation structure keyword to be specified in REPEATED statement.  This needs 
               to be specified if SUBJECT is specified.  Leave blank if not using a GEE model.
TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.
OUTPATH        Path for output table to be stored.
FNAME          File name for output table.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
* For more details, please see the related documentation.
**********************************************************************************************/
%macro UNI_GENMOD(DSN=, OUTCOME=, DESC=F, CLIST=, CREFLIST=, NLIST=, SUBJECT=, WITHINSUB=, 
   TYPE=, DIST=NORMAL, TYPE3=T, OUTPATH=, FNAME=, DEBUG = F); 

   %local label nlist_cnt clist_cnt var cvar nobs_used outlab __Macro_Err repeated distname
      refvar;

   /* Capitalize */
   %let debug = %UPCASE(&debug);
   %let type3 = %UPCASE(&type3);
   %let desc = %UPCASE(&desc);

   /* Initialize error flag */
   %let __Macro_Err = 0;

   /* Count number of categorical variables */
   %if &clist = %STR() %then %do;
      %let clist_cnt = 0;
   %end;
   %else %do;
      %let clist_cnt = %sysfunc(countw(&CLIST));
   %end;

   /* Count number of numeric variables */
   %if &nlist = %STR() %then %do;
      %let nlist_cnt = 0;
   %end;
   %else %do;
      %let nlist_cnt = %sysfunc(countw(&NLIST));
   %end;

   /* For GEE - both subject and type should be specified */
   %if (&subject = %STR() AND &type ~= %STR()) OR 
     (&subject ~= %STR() AND &type = %STR())%then %do;
      %put ERROR: When using a GEE model both SUBJECT and TYPE need to be specified.;
      %let __Macro_Err = 1;
   %end;

   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   /* Subject and correlation type must be specified to use repeated statement */
   %if &subject = %STR() %then %let repeated = F;
   %else %let repeated = T;

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

   /* Get outcome label */
   DATA _NULL_;
      set &dsn (obs=1);
       CALL SYMPUT('outlab',VLABEL(&outcome));
   RUN;

   %do i = 1 %to %EVAL(&clist_cnt+&nlist_cnt);

      %if &i <= &clist_cnt %then %do;
         %let var = %SCAN(&CLIST, &i);
         %let cvar = T;
         %let refvar = %SCAN(&CREFLIST,&i,*);

         /* If no reference categories were provided then use the default */
         %if %sysevalf(%superq(REFVAR)=,boolean) %then %do;
            %let refvar = &var;
         %end;
      %end;
      %else %do;
         %let var = %SCAN(&NLIST, %EVAL(&i-&clist_cnt));
         %let cvar = F;
      %end;

      ODS OUTPUT ModelInfo=_modelinf NObs = _nobs Type3 = _typ3 
          %if &repeated = T %then %do; GEEEmpPEst = _est (rename=(probz=pval)) %end;
          %else %do; ParameterEstimates=_est (rename=(parameter=parm ProbChiSq=pval LowerWaldCL=lowerCL upperwaldcl=upperCL)) %end;;
      PROC GENMOD DATA = &dsn %if &desc = T %then %do; DESC %end; namelen=32;
         class &subject &withinsub %if &cvar = T %then %do; &refvar %end; /ORDER=internal;
         model &outcome = &var/type3 wald dist=&dist;
         %if &repeated = T %then %do; 
            REPEATED SUBJECT=&subject/type=&type 
               %if &withinsub ~= %STR() %then %do; withinsubject=&withinsub %end;;
         %end;
      RUN;

      /* Get distribution used (in case abbrevations are used) */
      PROC SQL noprint;
         select cvalue1 into: distname
         from _modelinf
         where label1 = 'Distribution';
      QUIT;

      /* Calculate number of observations used */
      PROC SQL noprint;
         select NObsUsed 
         into :nobs_used
         from _nobs (obs=1);
      QUIT;

      /* Get variable label */
      DATA _NULL_;
         set &dsn (obs=1);
         CALL SYMPUT('label',VLABEL(&var));
      RUN;

      /* Merge estimates and type 3 p-values */
      DATA _mg;
         length parm $32 label $256;
         merge _est (where=(parm not in ('Intercept','Scale','Dispersion'))) 
               _typ3 (keep=source ProbChiSq rename=(source=Parm ProbChiSq=p_type3));
         by Parm;
         label = "&label";

         /* Calculate RR for NB or Poisson and CI */
         if estimate ~= 0 then do;
            /* For distributions that use log link - take exp of estimates */
            %if &distname = Negative Binomial OR &distname = Poisson OR &distname = Binomial %then %do;
               est = exp(estimate);
               LCL = exp(LowerCL);
               UCL = exp(uppercl);
            %end;
            /* For normal linear regression don't take exp */
            %else %if &distname = Normal %then %do;
               est = estimate;
               LCL = lowerCL;
               UCL = upperCL;
            %end;
         end;

         /* Save order */
         order1 = &i;
         order2 = _n_;
      RUN;

      %if &cvar = T %then %do;
         ODS OUTPUT 'One-Way Frequencies' = _freq; 
         PROC FREQ DATA = &dsn;
            tables &var;
            where MISSING(&outcome) = 0 AND MISSING(&var) = 0 %if &repeated = T %then %do; AND MISSING(&subject) = 0 %end;;
         RUN;
         DATA _freq;
            set _freq;
            level1 = LEFT(F_&var);
         RUN;
         PROC SORT DATA = _freq;
            by level1;
         RUN;
         PROC SORT DATA = _mg;
            by level1;
         RUN;
         /* Merge on N's */
         DATA _mg&i;
            length level1 $256;
            merge _mg _freq (keep=level1 Frequency);
            by level1;
            rename frequency = n;
         RUN;
      %end;
      %else %do;
         DATA _mg&i;
            set _mg;
            n = &nobs_used;
            level1 = ' ';
         RUN;
      %end;

      /* Return to original order */
      PROC SORT DATA = _mg&i;
         by order1 order2;
      RUN;

   %end;

   /* Combine all results */
   DATA _summ;
      set %do i = 1 %to %EVAL(&clist_cnt+&nlist_cnt); _mg&i %end;;

      /* Estimate and 95% CI */
      if est ~= . then est_CI = TRIM(LEFT(PUT(est,8.2))) || " (" || 
         TRIM(LEFT(PUT(LCL,8.2))) || "-" || 
         TRIM(LEFT(PUT(UCL,8.2))) || ")";
      else est_ci = '-';
   RUN;

   *---- table template -----;  
   ODS PATH WORK.TEMPLAT(UPDATE) SASUSR.TEMPLAT(UPDATE) SASHELP.TMPLMST(READ);

   PROC TEMPLATE;
   DEFINE STYLE STYLES.TABLES;
   NOTES "MY TABLE STYLE"; 
   PARENT=STYLES.MINIMAL;

     STYLE SYSTEMTITLE /FONT_SIZE = 12pt     FONT_FACE = "TIMES NEW ROMAN";

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

   OPTIONS ORIENTATION=PORTRAIT MISSING = "-" NODATE;
   ODS rtf STYLE=TABLES file= "&OUTPATH.&FNAME &SYSDATE..DOC"; 

   PROC REPORT DATA=_summ HEADLINE HEADSKIP CENTER SPANROWS LS=256 STYLE(REPORT)={JUST=CENTER} 
     SPLIT='~' nowd; 
      COLUMN order1 label Level1 N ("&outLab" '----------------------------------------------'
          (est_CI pval %if &type3 = T %then %do; p_type3 %end;)); 
      DEFINE order1/order order=internal noprint;
      DEFINE label/order order=data "Covariate"  STYLE(COLUMN) = {JUST = L}; 
      DEFINE Level1/ DISPLAY "Level" STYLE(COLUMN) = {JUST = L};
      DEFINE N/DISPLAY "N"; 
      DEFINE est_CI/DISPLAY STYLE(COLUMN) = {JUST = C}
          %if &distname = Negative Binomial OR &distname = Poisson %then %do; "Rate Ratio (95% CI)" %end;
          %else %if &distname = Binomial %then %do; "Odds Ratio (95% CI)" %end;
          %else %if &distname = Normal %then %do; "B (95% CI)" %end;; 
      DEFINE pval/DISPLAY STYLE(COLUMN)= {JUST = C} FORMAT=PVALUE8.3
          %if &distname = Normal %then %do; "B~P-value" %end;
          %else %if &distname = Binomial %then %do; "OR~P-value" %end;
          %else %do; "RR~P-value" %end;;; 

      COMPUTE pval; 
         IF pval <0.05 THEN CALL DEFINE("pval", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
      ENDCOMP; 

      %if &type3 = T %then %do;
         DEFINE p_type3/ORDER "Type3 P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%} 
          FORMAT=PVALUE8.3;
         COMPUTE p_type3; 
            IF p_type3 <0.05 THEN CALL DEFINE("p_type3", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         ENDCOMP;   
      %end;

      compute after label; 
         line '';   
      endcomp;      

   RUN;

   ODS RTF CLOSE;

   /* Reload original options that were in use before running the macro */
   PROC OPTLOAD data=_options;
   RUN;

   %if &debug = F %then %do;
      /* Delete temporary datasets */
      PROC DATASETS lib=work;
         DELETE _est _typ3 _mg _options _nobs _summ _modelinf _mg: 
               %if &clist ~= %STR() %then %do; _freq %end;;
      QUIT;
   %end;


%mend uni_genmod;


