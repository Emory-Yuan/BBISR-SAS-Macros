/**********************************************************************************************
***********************************************************************************************
Macro Name: UNI_NUM
Created Date/Author: Feb. 2012/Yuan Liu
Last Update Date/Person: Oct. 3, 2013/Dana Nickleach 
Current Version: V9
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu

Purpose:  To produce a summary table of the univariate association of a numerical outcome with
all other covariates in the data individually. For categorical covariates, the sample size, 
mean and median along with ANOVA test (parametric p-value) or Kruskal-Wallis test 
(non-parametric p-value) will be carried out. For numerical covariates, both Pearson 
correlation coefficient (parametric p-value) and Spearman correlation coefficient 
(non-parametric p-value) will be calculated.

Notes: 1) The order of variables in the summary table is the same as the input order. For the 
best results, you may want to put the demographic variables together and also clinical 
characteristics variables together; 2) The biostatistician may need to help investigator decide 
which statistics (parametric or non-parametric p-value) is more appropriate to the data.

Parameters: 

DATASET        The name of the data set to be analyzed.
OUTCOME        The numerical outcome variable. Each variable name must be less than 24 
               characters long.
CLIST          List of categorical variables, separated by empty space.
NLIST          List of numerical variables, separated by empty space.
NONPAR         Set to F to suppress non-parametric statistics.  The default value is T.
OUTPATH        Path for output table to be stored.
FNAME          File name for output table.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
For more details, please see the related documentation
**********************************************************************************************/

%macro UNI_NUM(dataset=, outcome=, clist=, nlist=, nonpar = T, outpath=, fname=descript, 
     debug=F);

   %let outcome = %UPCASE(&outcome);
   %let debug = %UPCASE(&debug);
   %let nonpar = %UPCASE(&nonpar);

   /* Count number of categorical variables */
   %if &clist ~= %STR() %then %let cvar_cnt = %sysfunc(countw(&clist));
   %else %let cvar_cnt = 0;

   /* Initialize error variable */
   %let __Macro_Err = 0;

   /* Check for variable names that are too long */
   %IF %LENGTH(&outcome) > 23 %then %do;
         %put ERROR: Variable name &outcome must be less than 24 characters.; 
           %let __Macro_Err=1;
   %END;

   /* Make sure that each categorical variable has at least two non-missing values */
   %do i = 1 %to &cvar_cnt;
      PROC SQL noprint;
         select count(distinct %SCAN(&clist, &i)) into :check
         from &dataset
         where MISSING(%SCAN(&clist, &i)) = 0;
      QUIT;

       %if &check <= 1 %then %do;
          %put ERROR: The variable %SCAN(&clist, &i) has less than two non-missing levels.  Please remove from CLIST.; 
           %let __Macro_Err=1;
       %end;
   %end;

   /* If there is an error in the parameters supplied then exit */
   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   /* Get list of data sets in work library to avoid deletion later */
   ODS EXCLUDE members Directory;
   ODS OUTPUT Members(nowarn)=DataSetList;
   PROC DATASETS lib=work memtype=data;
   QUIT;

   /* If there are data sets in the work library */
   %if %sysfunc(exist(DataSetList)) %then %do;
      PROC SQL noprint;
         select Name
         into :work_sets separated by ' '
         from DataSetList;
      quit;
   %end;
   %else %do;
      %let work_sets =;
   %end;

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

   /* Format missing values consistently */
   OPTIONS MISSING = " ";

   %IF &CLIST NE  %THEN %DO; 
      %let j = 1;

      %do %while (%scan(&clist, &j, " ") ne )    ;

         %let cvar = %scan(&clist, &j, " ");

         ODS OUTPUT "Summary statistics" = summary;
         PROC MEANS DATA=&dataset N MEAN MEDIAN; 
            VAR &outcome; 
            CLASS &cvar;
         RUN;

         ods output "ANOVA" = anova "Kruskal-Wallis Test" = nonpartest;
         proc npar1way data=&dataset; 
            class &cvar; var &outcome; run;

         data anova;set anova; format ProbF;run;
         data nonpartest;set nonpartest; format  nValue1;run;

         proc sql noprint; 
            select  nValue1 into :kwp  from nonpartest where Label1 = "Pr > Chi-Square"; 
            select  ProbF into :anovap from anova where ProbF > 0; quit;

         Data summary&j;
            length    Level $20. covariate $256.;
            set summary;
            Covariate = label(&cvar);
            Anova_P = &anovap;
            KW_P = &kwp;
            Level = strip(vvalue(&cvar));
            drop &cvar;
         run;

         data summary&j;set summary&j;if substr(level, 1,1) not in ( " " "."); run;

         %let j = %eval(&j + 1) ;
      %end;

      %let j=   %eval(&j - 1) ;

      Data summmaryCategorical; set summary1 - summary&j;run;
   %END;


   %IF &NLIST NE  %THEN %DO; 

          %let i = 1;
         %do %while (%scan(&nlist, &i, " ") ne ); 
               %let nvar = %scan(&nlist, &i, " ");
                    ods output "Spearman Correlations" = spearman  "Pearson Correlations" = pearson "Simple Statistics" = simples;
                    proc corr data=&dataset spearman pearson;
                    var &outcome &nvar;run;
                    
                    data _NULL_; set &dataset; call symput('vv', put(label(&nvar),$96.));run;
                    data spearman; set spearman; format  P&outcome P&nvar;run;
                    data pearson; set pearson; format  P&outcome P&nvar;run;

                    proc sql noprint;
                    select NObs into: nobs separated by " " from simples;  quit;


                    proc sql noprint;
                    select &nvar into: pear_coef  from pearson where UPCASE(Variable) = "&outcome ";
                    select P&nvar into: pear_p from pearson where  UPCASE(Variable) = "&outcome ";
                    select &nvar into: spear_coef from spearman where UPCASE(Variable) = "&outcome ";
                    select P&nvar into: spear_p from spearman where  UPCASE(Variable) = "&outcome ";
               
                    
                    %if  %scan(&nobs, 1, " ") ne %scan(&nobs, 2, " ") %then %do  ;
                    select N&nvar into: pear_n from pearson where UPCASE(Variable) = "&outcome "; 
                    %end;
                    %else %do  ;
                    select NObs into: pear_n from simples where UPCASE(Variable) = "&outcome "; 
                    %end;

                    quit; 

                    %put &pear_coef;

                    Data coefsummary&i;
                    length covariate $256.;
                    Covariate =  "&vv";
                    N = &pear_n;
                    SpearmanCC = &spear_coef;
                    SpearmanP = &spear_p;
                    PearsonCC = &pear_coef;
                    PearsonP = &pear_P;
                    run;
          %let i = %eval(&i + 1) ;
          %end;
          %let i=   %eval(&i - 1); 

          Data summmaryNumerical; set coefsummary1 - coefsummary&i;run;
%END;


*---- table template -----;  

ODS PATH WORK.TEMPLAT(UPDATE)
SASUSR.TEMPLAT(UPDATE) SASHELP.TMPLMST(READ);

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

   *------- build the table -----;
   data _NULL_; 
      set &dataset; 
      call symput('outcomevv', put(label(&outcome),$256.));
   run;

   OPTIONS ORIENTATION=PORTRAIT MISSING = "-" NODATE;

   ODS RTF STYLE=TABLES FILE= "&OUTPATH.&FNAME &SYSDATE..DOC" startpage=never; 

   %IF &CLIST NE  %THEN %DO; 

      PROC REPORT DATA=summmarycategorical HEADLINE HEADSKIP CENTER STYLE(REPORT)={JUST=CENTER}
       SPLIT='~' nowd LS=256 SPANROWS;
         COLUMNS Covariate  Level &outcome._N (" &outcomevv " 
               %if &NONPAR = T %then %do; '______________________________________' %end; 
               %else %do; '_____________________' %end;
                    ( &outcome._Mean 
               %if &NONPAR = T %then %do; &outcome._Median %end; Anova_P 
               %if &NONPAR = T %then %do; KW_p %end;)); 
          DEFINE covariate/order  order=data    "Variable"  STYLE(COLUMN) = {JUST = L} ;
          DEFINE Level/DISPLAY "Level"     STYLE(COLUMN) = {JUST = L};
           DEFINE &outcome._N/DISPLAY "N" STYLE(COLUMN) = {JUST = C}  format=8.;
           DEFINE &outcome._Mean/DISPLAY "Mean"  STYLE(COLUMN) = {JUST = C}  format=8.2;
           DEFINE Anova_P/ORDER "ANOVA P-value"  STYLE(COLUMN) = {JUST = C} format=PVALUE5.3 ;

           %if &NONPAR = T %then %do;
              DEFINE &outcome._Median/DISPLAY "Median" STYLE(COLUMN) = {JUST = C}  format=8.2;
              DEFINE KW_P/ORDER "Kruskal-Wallis P-value" STYLE(COLUMN) = {JUST = C} format=PVALUE5.3;
               COMPUTE KW_P; 
                 IF KW_P <0.05 THEN CALL DEFINE("KW_P", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
              ENDCOMP; 
           %end;

         COMPUTE Anova_P; 
              IF Anova_P <0.05 THEN CALL DEFINE("Anova_P", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
           ENDCOMP; 

           compute after covariate; line ''; endcomp; 
      RUN; 

   %end;

   %IF &NLIST NE  %THEN %DO; 

      PROC REPORT DATA=summmaryNumerical HEADLINE HEADSKIP CENTER STYLE(REPORT)={JUST=CENTER} 
          SPLIT='~' nowd LS=256 SPANROWS;
         COLUMNS Covariate  N 
                (" &outcomevv " '______________________________________'(PearsonCC PearsonP 
               %if &NONPAR = T %then %do; SpearmanCC SpearmanP %end;)); 
          DEFINE covariate/DISPLAY "Variable"  STYLE(COLUMN) = {JUST = L} ;
           DEFINE N/DISPLAY "N" STYLE(COLUMN) = {JUST = C}  format=5.;
           DEFINE PearsonCC/DISPLAY "Pearson CC"  STYLE(COLUMN) = {JUST = C} format=6.3;
           DEFINE PearsonP/DISPLAY "Pearson P-value" STYLE(COLUMN) = {JUST = C} format=PVALUE5.3 ;
         %if &NONPAR = T %then %do;
              DEFINE SpearmanCC/DISPLAY "Spearman CC" STYLE(COLUMN) = {JUST = C} format=6.3;
              DEFINE SpearmanP/DISPLAY "Spearman P-value" STYLE(COLUMN) = {JUST = C} format=PVALUE5.3;

               COMPUTE SpearmanP; 
                 IF SpearmanP <0.05 THEN CALL DEFINE("SpearmanP", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
              ENDCOMP;
           %end;

           COMPUTE PearsonP; 
              IF PearsonP <0.05 THEN CALL DEFINE("PearsonP", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
           ENDCOMP; 

      RUN; 

   %end;

   ODS RTF CLOSE; 

   /* Reload original options that were in use before running the macro */
   PROC OPTLOAD data=_options;
   RUN;

   /* Only delete files if not in debug mode */
   %if &debug ~= T %then %do;

      /* If there are work data sets that should not be deleted */
      %if %sysevalf(%superq(work_sets)~=,boolean) %then %do;
         /* DELETE ALL TEMPORARY DATASETS that were created */
         proc datasets lib=work memtype=data noprint;  
            save &work_sets;
         quit;  
      %end;
      %else %do;
         proc datasets lib=work kill memtype=data noprint;  
         quit; 
      %end;
   %end;

%mend; 


