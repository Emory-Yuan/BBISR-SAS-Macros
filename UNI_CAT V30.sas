/**********************************************************************************************
***********************************************************************************************
Macro Name: UNI_CAT
Created Date/Author/Contact: Feb. 2012/Yuan Liu/YLIU31@emory.edu
Last Update Date/Person/Contact: Oct 6, 2015/Yuan Liu/YLIU31@emory.edu
Current Version: V30

Other Contributors: Dana Nickleach
Working Environment: SAS 9.3 English version

Purpose:  To conduct a univariate analysis for a categorical outcome with a list of covariates,
individually.  For categorical covariates, a contingency table along with the Chi-square test 
(parametric p-value) or Fisher’s exact test (non-parametric p-value) can be produced. For 
numerical covariates, the sample size, mean and median along with ANOVA test (parametric 
p-value) or Kruskal-Wallis test (non-parametric p-value) can be produced.

Notes: 1) The order of variables in the summary table is the same as the input order. For the 
best results, you may want to put the demographic variables together and also clinical 
characteristics variables together; 2) The biostatistician may need to help investigator to 
decide which statistics (parametric or non-parametric p-value) is more appropriate for the data.

Parameters: 

DATASET        The name of the data set to be analyzed.
OUTCOME        Categorical variable to be associated with CLIST and NLIST variables.  More than
               one variable can be listed separated by empty space.  However, these variables 
               appear in the table header and too many variables will cause the table to wrap 
               due to the document page width limitations, producing undesirable results.  Each 
               variable name must not be more than 30 characters long.                                                                                     
CLIST          List of categorical variables, separated by empty space.
NLIST          List of numerical variables, separated by empty space.
NONPAR         Specify a value of F, T, or A to indicate whether to conduct non-parametric 
               tests.  If the value is T then both parametric and non-parametric tests will be 
               conducted.  If the value is F then only parametric tests will be conducted.  A 
               value of A means that for categorical variables, the appropriate test statistic, 
               non-parametric or parametric, will be automatically chosen based on whether the 
               chi-square test is invalid, but for numerical covariates only the parametric test 
               will be calculated.  Option A is only available for SAS V9.3 or later.  The 
               default value is F.
SPREAD         Set to T to also report standard deviation, min, and max for numerical variables.  
               The default value is F. 
BY             A separate analysis will be conducted for each value of the variable specified 
               here.  
MHC            Set to T to report p-values from Mantel-Haenszel chi-square tests instead of 
               Pearson chi-square tests.  The default value is F.   
DOC            Set to F to suppress creation of the RTF file.  The default value is T.
OUTPATH        File path for output table to be stored.
FNAME          File name for output table.
ROWPERCENT     Set to F to report column percentages instead of row percentages from the 
               contingency table.  The default value is T.
ORIENTATION    Value of PORTRAIT or LANDSCAPE to indicate the page layout of the report.  The 
               default value is PORTRAIT.
WEIGHT         Weight variable to use in a WEIGHT statement.  Weights will not be normalized by
               the macro.  The reported N will be the sum of the weights. This option will not 
               work with NONPAR = T or A.  Leave blank if not weighting.
MATCHID        If your data is from a matched sample, indicate the id variable that links 
               matched pairs.  If not, then leave this blank.  If you specify a MATCHID then
               McNemar's test for two level categorical variables and Bowker's test of symmetry
               for more than two level categorical variables will be conducted instead of a 
               chi-square test.  A paired t-test as opposed to ANOVA for numerical variables 
               will be conducted.  This option is currently not set up to conduct 
               non-parametric tests.  This option is also only appropriate for 1-1 matching.
               Note that the data set should be in the format of one observation per subject, 
               not one observation per match.  The data will be transformed as needed in order
               to conduct the necessary tests.  The OUTCOME variable should correspond to the
               variable that identifies the repeated measurement.  For example, if a patient
               had one measurement on their right foot and one on their left then the OUTCOME
               variable should be foot.  The OUTCOME variable should have two categories.                  
DEBUG          Set to T if running in debug mode.  Work datasets will not be deleted in debug 
               mode.  This is useful if you are editing the code or want to further manipulate 
               the resulting data sets.  The default value is F.

***********************************************************************************************
For more details, please see the related documentation
**********************************************************************************************/

%MACRO UNI_CAT(DATASET=, outcome=, CLIST=, NLIST=, NONPAR=F, SPREAD=F, BY=, DOC=T, OUTPATH=,
     FNAME =, ROWPERCENT = T, ORIENTATION = PORTRAIT , WEIGHT=, MHC = F, MATCHID=, DEBUG=F); 

   /* Make sure these vars are local */
   %local OUTVAR NUM_OUT i STATL CVAR TVAR j __MACRO_ERR WORK_SETS COLABEL m TABLESUM 
     NONP_CTEST n COL OUT_CNT NVAR CTEST VV cvar_cnt check pair_var;

   /* Count number of outcomes */
   %let out_cnt = %sysfunc(countw(&outcome));
   %if &clist ~= %STR() %then %let cvar_cnt = %sysfunc(countw(&clist));
   %else %let cvar_cnt = 0;

   /* Prevent case sensitivity */
   %let nonpar = %UPCASE(&nonpar);
   %let spread = %UPCASE(&spread);
   %let rowpercent = %UPCASE(&rowpercent);
   %let debug = %UPCASE(&debug);
   %let MHC = %UPCASE(&MHC);
   %let doc = %UPCASE(&doc);
   %let outcome = %UPCASE(&outcome);
   %let clist = %UPCASE(&clist);
   %let by = %UPCASE(&by);

   /* Initialize error flag */
   %let __Macro_Err = 0;

   %if &weight ~= %STR() AND (&NONPAR = TRUE OR &NONPAR = T OR &NONPAR = A) %then %do;
      %put ERROR: Cannot use weights with nonpar=T.;
       %let __Macro_Err = 1;
   %end;

   %if &SYSVER < 9.3 AND &NONPAR=A %then %do;
      %put ERROR: Cannot use nonpar=A unless you have SAS V9.3 or greater;
       %let __Macro_Err = 1;
   %end;

   /* Check for outcome variable names that are too long */
   %DO i = 1 %to &out_cnt; 
      %IF %LENGTH(%SCAN(&outcome, &i)) > 30 %then %do;
         %put ERROR: Variable name %SCAN(&outcome, &i) must be less than 31 characters.; 
           %let __Macro_Err=1;
       %END;
   %END;

   /* Make sure that outcome variables are also not listed in CLIST */
   %do i = 1 %to &out_cnt;
     %do j = 1 %to &cvar_cnt;
         %if %SCAN(&outcome, &i) = %SCAN(&clist, &j) %then %do;
             %put ERROR: Outcome %SCAN(&outcome, &i) cannot appear in CLIST as well.; 
             %let __Macro_Err=1;
          %end;
     %end;
   %end;

   /* Make sure that BY variable is also not listed in CLIST */
   %do i = 1 %to &cvar_cnt;
      %if &by = %SCAN(&clist, &i) %then %do;
             %put ERROR: BY variable &by cannot appear in CLIST as well.; 
             %let __Macro_Err=1;
          %end;
   %end;

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

   /* Some options cannot be used with MATCHID */
   %if &matchid ~= %STR() %then %do;
      %if &nonpar = T %then %do;
         %put ERROR: MATCHID cannot be used in conjunction with NONPAR=T.; 
         %let __Macro_Err=1;
      %end;
      %else %if &nonpar = A %then %do;
         %put ERROR: MATCHID cannot be used in conjunction with NONPAR=A.; 
         %let __Macro_Err=1;
      %end;
      %if &weight ~= %STR() %then %do;
         %put ERROR: MATCHID cannot be used in conjunction with WEIGHT.; 
         %let __Macro_Err=1;
      %end;
      %if &MHC = T %then %do;
         %put ERROR: MATCHID cannot be used in conjunction with MHC=T.; 
         %let __Macro_Err=1;
      %end;
      %if &BY ~= %STR() %then %do;
         %put ERROR: MATCHID cannot be used in conjunction with BY.; 
         %let __Macro_Err=1;
      %end;

      /* Can't use with more than one outcome */
      %if &out_cnt > 1 %then %do;
         %put ERROR: MATCHID can only be used with 1 OUTCOME variable.; 
         %let __Macro_Err=1;
      %end;

      /* Count number of outcome levels */
      /* Don't use PROC SQL b/c select distinct counts unique values not formats and PROC
      TRANSPOSE will transpose based on the formatted values of the ID variable */
     
      PROC FREQ DATA = &dataset noprint;
         tables &outcome/out=_check;
      RUN;
      PROC SQL noprint;
         select count(*) into :check
         from _check;
      QUIT;
      /* The outcome should not have more than two levels */
      %if &check ~= 2 %then %do; 
         %put ERROR: When MATCHID is used the OUTCOME should only have 2 levels.  %TRIM(&outcome) has %TRIM(&check) levels.; 
         %let __Macro_Err=1;
      %end;
   %end;

   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   /* Get list of data sets in work library to avoid deletion later */
   PROC SQL noprint;
      select memname
      into :work_sets separated by ' '
      from sashelp.vmember
      where libname = 'WORK' and memtype = 'DATA';
   quit;

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

   /* Format missing values consistently */
   OPTIONS MISSING = " ";

   %IF &ROWPERCENT = T %THEN %DO; 
      %LET TABLESUM = RowPercent; 
      %LET STATL = N (Row %); 
   %end; 
   %ELSE  %DO; 
       %LET TABLESUM = ColPercent; 
       %LET STATL = N (Col %);  
   %end;

   /* Test name to use in footnotes */
   %if &MHC = T %then %do;
      %let ctest = Mantel-Haenszel chi-square;
      %let nonp_ctest = exact Mantel-Haenszel chi-square;
   %end;
   %else %do;
      %let ctest = chi-square;
      %let nonp_ctest = %STR(Fisher%'s exact);
   %end;

   /* If using a by variable sort data set */
   %if &by ~= %STR() %then %do;
      PROC SORT DATA = &dataset out = _data;
          by &by;
       RUN;
   %end;
   %else %do;
      DATA _data;
         set &dataset;
      RUN;
   %end;


   * Character Variables;
   %IF &CLIST NE  %THEN %DO; 

      %let m = 1;
      %let outvar = %SCAN(&outcome, &m);

      /* Repeat for each outcome variable */
      %DO %UNTIL (&outvar = %STR());

         %let N = 1; 
         %let cvar = %SCAN(&CLIST, &N);

         %DO %UNTIL (&cvar =  ); 

            /* Parametrics and Non-parametrics */
            %IF &NONPAR = TRUE OR  &NONPAR = T OR &NONPAR=A %THEN %DO; 
				ODS SELECT NONE;
               /*ODS EXCLUDE "Cross-Tabular Freq Table" "Chi-Square Tests" "Fisher's Exact Test"
                  %if &MHC = T %then %do; MHChiSq %end;;*/
               ODS OUTPUT  "Cross-Tabular Freq Table" = _freq "Chi-Square Tests" = _chi 
                  "Fisher's Exact Test" = _fisher %if &MHC = T %then %do; MHChiSq = _mchi %end; ;
               PROC FREQ DATA=_data;  
                  TABLE &CVAR*&outvar/exact nopercent            
                         %if &NONPAR=A %then %do; chisq(warn=output) %end;;
                      /* Request exact test for Mantel-Haenszel Chi-Square */
                  %if &MHC=T %then %do; EXACT MHCHI; %end; 
                  %if &by ~= %STR() %then %do; BY &by; %end;
               RUN; 
			   ODS SELECT ALL;

               data _freqt; 
                  length measure $96. outv $96.;
                  set _freq;
                  measure = catt(Frequency, " (", round(&TABLESUM, 0.01), ")");
                  outv = strip(vvalue(&outvar));

                  /* Concanate outcome number and category */
                  outcat = "O&m" || outv;

                  /* Get rid of total rows that won't be used */
                  where _TYPE_ not in ('01' '00' '10');
               run;

               data _freqt;
                  set _freqt; 
                  where substr(outv, 1,1) not in ( " " "."); 
               run; 

               proc sort data=_freqt; 
                  by &by &CVAR;
               run;

               proc transpose data=_freqt out=_freqtt; 
                  var measure;
                  by &by &CVAR; 
                  id outcat;
               run;

               /* Create list of transposed variable names */
               PROC CONTENTS DATA = _freqtt (drop=_NAME_ &by &CVAR) out=_cont noprint;
               RUN;

               PROC SQL noprint;
                  select name into :tvar separated by ' '
                  from _cont;
               QUIT;

               /* Fill in emtpy cells with zeros */
               DATA _freqtt;
                  set _freqtt;
                  array tran &tvar;
                  do over tran;
                     %if &ROWPERCENT = T %then %do;
                        if tran = ' ' then tran = '0 (0)';
                     %end;
                     /* Column % can't be calculated */
                     %else %do;
                        if tran = ' ' then tran = '0 (NA)';
                     %end;
                   end;
               RUN;

               %if &nonpar = A %then %do;
                  PROC SORT DATA = _chi;
                     by &by table;
                  RUN;

                  /* Replicate warning indicator for all type of test (including M-H) */
                  /* The indicator is only present for pearson chi-square */
                  DATA _chi (rename=(warn2=warning) drop=warning);
                     set _chi;
                     by &by table;
                     retain Warn2;

                    if first.table then warn2 = warning;
                 RUN;
               %end;

               /* Merge test results */
               DATA _test;
                  merge _chi (keep=&by prob statistic  %if &NONPAR=A %then %do; warning %end;
                         /* Mantel-Haenszel or Pearson Chi-Square */
                         %if &MHC = T %then %do; where=(Statistic = "Mantel-Haenszel Chi-Square")) %end;
                         %else %do; where=(Statistic = "Chi-Square")) %end;
                         /* Merge on exact tests */
                         %if &MHC=T %then %do; 
                            _mchi (keep=&by nvalue1 name1 where=(Name1='XP_MHCHI')) 
                    %end;
                         %else %do;
                            _fisher (keep=&by nvalue1 name1 where=(Name1='XP2_FISH'))
                         %end;;

                      %if &by ~= %STR() %then %do; BY &by; %end;

                      drop statistic name1;
                      FORMAT Prob nvalue1;
               RUN;

               data _freq&N;
                  length covariate $256. statistics $96. level $96.;
                  %if &by ~= %STR() %then %do;
                     merge _freqtt _test;
                     BY &by; 
                  %end;
                  %else %do;
                     set _freqtt;
                     if _n_ = 1 then set _test;
                  %end;

                  Covariate= label(&CVAR);
                  statistics = "&STATL";
                  Level =  strip(vvalue(&CVAR));
                  par_p = Prob;
                  nonpar_p = nValue1;
                  /* If auto selection use chi-square if valid otherwise Fishers */
                  %if &NONPAR=A %then %do; 
                     if warning = 1 then par_p = nonpar_p;
                     drop warning;
                  %end;
                  drop _NAME_ &CVAR;
               run;

               data _freq&N;
                  set _freq&N;  
                  if substr(level, 1,1) not in ( " " "."); 
                  /* Order variable - keep original order */
                  _order = &n;

                  rename par_p=par_p&m nonpar_p = nonpar_p&m;
               run;
            %END;


            /* Parametrics */
            %IF &NONPAR = FALSE OR &NONPAR = F %THEN %DO; 

               /* For matched data - calculate p-value */
               %if &matchid ~= %STR() %then %do;

                  PROC SORT DATA = _data;
                     by &matchid;
                  RUN;
                  /* Transpose to one observation per pair */
                  PROC TRANSPOSE DATA = _data out = _tran;
                     by &matchid;
                     id &outvar;
                     var &CVAR;
                  RUN;
                  /* Get transposed column names */
                  PROC SQL noprint;
                     select a.name into :pair_var separated by '*'
                     from sashelp.vcolumn a
                     where a.memtype = 'DATA' and a.libname = 'WORK' and a.memname = "_TRAN"
                     and UPCASE(a.name) not in ('_NAME_' '_LABEL_' "%UPCASE(&matchid)");
                  QUIT;

                  PROC FREQ DATA = _tran;
                     tables &pair_var/agree;
                     output out=_agree agree;
                  RUN;
               %end;

               ODS SELECT NONE;
               ODS OUTPUT  "Cross-Tabular Freq Table" = _freq 
                  %if &matchid = %STR() %then %do; "Chi-Square Tests" = _chi %end;;
               PROC FREQ DATA=_data;  
                  TABLE &CVAR*&outvar/nopercent CHISQ %if &SYSVER >= 9.3 %then %do;(warn=output)%end;;
                  %if &weight ~= %STR() %then %do;
                     weight &weight;
                  %end;
                  %if &by ~= %STR() %then %do; BY &by; %end;
               RUN; 
			   ODS SELECT ALL;

               /* For matched data - don't use chi-square so don't need to check for warning */
               %if &matchid = %STR() %then %do;
                  /* For V9.3 and above print warning message to log if chi-square is invalid */
                  %if &SYSVER >= 9.3 %then %do;
                     %if &by = %STR() %then %do;
                        PROC SQL noprint;
                           select count(*) into :nwarn
                           from _chi
                           where Warning = 1;
                        QUIT;

                        %if &nwarn ~= 0 %then 
                           %PUT WARNING: Chi-Square may not be a valid test for &cvar * &outvar..;
                     %end;
                     %else %do;
                        PROC SQL noprint;
                           select count(*), &by into :nwarn, :bygrp separated by ' '
                           from _chi
                           where Warning = 1;
                        QUIT;

                        %if &nwarn ~= 0 %then 
                              %PUT WARNING: Chi-Square may not be a valid test for &cvar * &outvar for &by: &bygrp..;

                     %end;

                  %end;
               %end;

               data _freqt; 
                  length measure $96. outv $96.;
                  set _freq;
                  /* Round freq's - they will only be decimals if weighting is used */
                  Frequency = ROUND(Frequency,1);
                  measure = catt(Frequency, " (", round(&TABLESUM, 0.01), ")");
                  outv = strip(vvalue(&outvar));
                  /* Concanate outcome number and category */
                  outcat = "O&m" || outv;
                  /* Get rid of total rows that won't be used */
                  where _TYPE_ not in ('01' '00' '10');
               run;

               data _freqt;
                  set _freqt; 
                  where substr(outv, 1,1) not in ( " " "."); 
               run; 

               proc sort data=_freqt; 
                  by &by &CVAR;
               run;

               proc transpose data=_freqt out=_freqtt; 
                  var measure;
                  by &by &CVAR; 
                  id outcat;
               run;

               /* Create list of transposed variable names */
               PROC CONTENTS DATA = _freqtt (drop=_NAME_ &by &CVAR) out=_cont noprint;
               RUN;

               PROC SQL noprint;
                  select name into :tvar separated by ' '
                  from _cont;
               QUIT;

               /* Fill in emtpy cells with zeros */
               DATA _freqtt;
                  set _freqtt;
                  array tran &tvar;
                  do over tran;
                     %if &ROWPERCENT = T %then %do;
                         if tran = ' ' then tran = '0 (0)';
                     %end;
                     /* Column % can't be calculated */
                     %else %do;
                        if tran = ' ' then tran = '0 (NA)';
                     %end;                     
                  end;
               RUN;

               /* For matched data */
               %if &matchid ~= %STR() %then %do;
                  /* Test results */
                  DATA _test;
                     set _agree;
                     /* Depending on how many categories there are mcnemar's test or Bowker's
                     test of symmetry will be used */
                     if P_TSYMM ~= . then prob = P_TSYMM;
                     else prob = P_MCNEM;

                     format Prob;
                     keep prob;
                  RUN;
               %end;
               %else %do;
                  /* Test results */
                  DATA _test;
                     set _chi (keep=&by prob statistic 
                     /* Mantel-Haenszel or Pearson Chi-Square */
                     %if &MHC = T %then %do; where=(Statistic = "Mantel-Haenszel Chi-Square")) %end;
                     %else %do; where=(Statistic = "Chi-Square")) %end;;

                     format Prob;
                     drop statistic;
                  RUN;
               %end;

               data _freq&N;
                  length covariate $256. statistics $96. level $96.;
                  %if &by ~= %STR() %then %do;
                     merge _freqtt _test;
                     BY &by; 
                  %end;
                  %else %do;
                     set _freqtt;
                     if _n_ = 1 then set _test;
                  %end;
                  Covariate= label(&CVAR);
                  statistics = "&STATL";
                  Level =  strip(vvalue(&CVAR));
                  par_p = prob;
                  drop _NAME_ &CVAR;
               run;

               data _freq&N;
                  set _freq&N;  
                  if substr(level, 1,1) not in ( " " "."); 
                  /* Order variable - keep original order */
                  _order = &n;
                  rename par_p=par_p&m;
               run;
            %END;

            %LET N = %EVAL(&N+1);
            %let cvar = %SCAN(&CLIST, &N);

         %end; 

         /* Combine categorical results */
         %LET N = %EVAL(&N-1);
         DATA _freq_all&m; 
            SET _freq1-_freq&N;
            /* Order variable - keep original order */
            _order2 = _N_;
         RUN; 

         /* Sort for merge */
         PROC SORT DATA = _freq_all&m; 
            by &by covariate statistics level;
         RUN;

         %let m = %EVAL(&m+1);
         %let outvar = %SCAN(&outcome, &m);

      %end;

      /* Merge categorical results from multiple outcome vars */
      %let m = %EVAL(&m-1);
      DATA _freq_all;
         merge _freq_all1 - _freq_all&m;
         by &by covariate statistics level;
         /* Order variable - to keep categorical first */
         _orderTyp = 1;
      RUN;

       /* Put back into original order */
      PROC SORT DATA = _freq_all;
         by &by _order _order2;
      RUN;

   %END;

   *NUMERIC VARIABLES ;
   %IF &NLIST NE  %THEN %DO; 

      %let m = 1;
      %let outvar = %SCAN(&outcome, &m);

      /* Repeat for each outcome variable */
      %DO %UNTIL (&outvar = %STR());

         %let N = 1; 
         %let nvar = %SCAN(&NLIST, &N);

         %DO %UNTIL (&nvar =  ); 

            PROC SORT DATA=_data;
               BY &by &outvar;
            RUN;

            PROC MEANS DATA=_data noprint; 
               var &NVAR; 
               class &by &outvar;  
               output out=_summary mean=mean2 median=median2 std=std2 min=min2 max=max2
                 %if &weight ~= %STR() %then %do; 
                  SUMWGT=n2 
               %end;
               %else %do; n=n2 %end;; 
               %if &weight ~= %STR() %then %do;
                    weight &weight;
               %end;
            RUN;

            data _summary;
               set _summary; 
               N = LEFT(PUT(ROUND(N2,1),best12.)); 
               Mean = LEFT(PUT(round(mean2,0.01),best12.));
               Median = LEFT(PUT(round(median2, 0.01),best12.));
               Std = LEFT(PUT(round(std2,.01),best12.));
               Min = LEFT(PUT(round(min2,.01),best12.));
               Max = LEFT(PUT(round(max2,.01),best12.));

               /* Concanate outcome number and category */
               outcat = "O&m" || strip(vvalue(&outvar));

               /* Labeling mean so that _LABEL_ will exist in transposed file regardless of */
               /* spread */
               LABEL mean = 'Mean'
                    std = 'Std Dev';
            run;

            proc sort data=_summary; 
               by &by &outvar;
            run;
            proc transpose data=_summary out=_summaryt; 
               var N mean median %if &spread = T %then min max std;; 
               id outcat;
                  %if &by ~= %STR() %then %do; BY &by; %end;
               where MISSING(&outvar) ~= 1;
            run;

            /* Create list of transposed variable names */
            PROC CONTENTS DATA = _summaryt (drop=_NAME_ _LABEL_ &by ) out=_cont noprint;
            RUN;

            PROC SQL noprint;
               select name into :tvar separated by ' '
               from _cont;
            QUIT;

            /* Fill in emtpy cells with zeros */
            DATA _summaryt;
               set _summaryt;
               array tran &tvar;
               do over tran;
                  if _NAME_ = 'N' then do;
                      if tran = ' ' then tran = '0';
                  end;
                  else if tran = ' ' then tran = 'NA';
               end;
            RUN;

            /* If using weights can't use NPAR1WAY */
            %if &weight ~= %STR() %then %do;
               ODS OUTPUT OverallANOVA=_anova;
               PROC GLM DATA = _data PLOTS=NONE;
                  class &outvar;
                  model &nvar = &outvar;
                  weight &weight;
                      %if &by ~= %STR() %then %do; BY &by; %end;
               RUN;
               QUIT;
            %end;
            %else %do;
               /* For matched data */
               %if &matchid ~= %STR() %then %do;

                  PROC SORT DATA = _data;
                     by &matchid;
                  RUN;
                  /* Transpose to one observation per pair */
                  PROC TRANSPOSE DATA = _data out = _tran;
                     by &matchid;
                     id &outvar;
                     var &nVAR;
                  RUN;
                  /* Get transposed column names */
                  PROC SQL noprint;
                     select a.name into :pair_var separated by '-'
                     from sashelp.vcolumn a
                     where a.memtype = 'DATA' and a.libname = 'WORK' and a.memname = "_TRAN"
                     and UPCASE(a.name) not in ('_NAME_' '_LABEL_' "%UPCASE(&matchid)");
                  QUIT;

                   /* Calculate the difference between variables */
                  DATA _tran;
                     set _tran;
                     diff = &pair_var;
                  RUN;
                  /* Calculate paired t-test */
                  ODS OUTPUT TestsForLocation=_pairedT;
                  PROC UNIVARIATE DATA = _tran;
                     var diff;
                  RUN;
               %end;
               %else %do;
                  proc npar1way data=_data PLOTS=NONE ANOVA WILCOXON noprint; 
                     class &outvar; 
                     var &NVAR; 
                     %if &by ~= %STR() %then %do; BY &by; %end;
                     output out=_anova (rename=(P_F=probf)) anova wilcoxon;
                  run;
               %end;
          
            %end;

            /* Test results */
            DATA _test;
               /* Can't calculate nonpar with weights */
               %if &weight = %STR() and &matchid = %STR() %then %do;
                  set _anova (keep=&by probf P_KW);
                   format P_KW;
               %end;
               /* For matched data */
               %else %if &matchid ~= %STR() %then %do;
                  set _pairedT (where=(test="Student's t") RENAME=(pvalue=probf)) ;
                  keep probf;
               %end;
               %else %do;
                  set _anova (keep=&by probf where=(ProbF > 0));
               %end;
               format ProbF;
            RUN;

            data _NULL_; 
               set _data; 
               call symput('vv', put(label(&NVAR),$256.));
            run;

            Data _summary&N;
               length covariate $256. statistics $96. level $96.;
               %if &by ~= %STR() %then %do;
                  merge _summaryt _test;
                  BY &by; 
               %end;
               %else %do;
                  set _summaryt;
                  if _n_ = 1 then set _test;
               %end;

               Covariate = "&vv";
               if _LABEL_ ~= ' ' then statistics = _LABEL_;
               else statistics = _NAME_ ;
               level = ' ';
               par_p = probf;
               /* Can't calculate nonpar with weights */
               /* Nonpar not calculated with matched data */
               %if &weight ~= %STR() or &matchid ~= %STR() %then %do;
                  nonpar_P = .;
               %end;
               %else %do;
                  nonpar_P = P_KW;   
               %end;

               /* Order variable - keep original order */
               _order = &n;

               drop  _NAME_;
               rename par_p=par_p&m nonpar_p = nonpar_p&m;
            run;

            %LET N = %EVAL(&N+1);
            %let nvar = %SCAN(&NLIST, &N);

         %end; 

         /* Combine numerical variable results */
         %LET N = %EVAL(&N-1);
         DATA _Summary_all&m; 
            SET _summary1 - _summary&N; 
            /* Order variable - keep original order */
            _order2 = _N_;

            %IF &NONPAR = FALSE OR &NONPAR = F OR &NONPAR = A %THEN %DO; DROP nonpar_P&m; %end;
         RUN;

         /* Sort for merge */
         PROC SORT DATA = _Summary_all&m; 
            by &by covariate statistics level;
         RUN;

         %let m = %EVAL(&m+1);
         %let outvar = %SCAN(&outcome, &m);
       
      %end;

      /* Merge numerical results from multiple outcome vars */
      %let m = %EVAL(&m-1);
      DATA _Summary_all;
         merge _Summary_all1 - _Summary_all&m;
         by &by covariate statistics level;
          /* Order variable - to keep categorical first */
          _orderTyp = 2;
      RUN;

      /* Put back into original order */
      PROC SORT DATA = _Summary_all;
         by &by _order _order2;
      RUN;
   %END;


   /* Combine categorical and numerical results */
   DATA _report;
      set %if &CLIST ~= %STR() %then %do; _freq_all %end; 
          %if &NLIST ~= %STR() %then %do; _summary_all; %end;;
   RUN;

   /* For each outcome */
   %let m = 1;
   %let outvar = %SCAN(&outcome, &m);
   %do %until (&outvar = %STR());

      /* Get outcome categories & N for header row */
      /* If BY statement is used this will still be the total N for the dataset */
      PROC FREQ DATA=_data noprint;  
         TABLE &outvar/PLOTS=NONE out=_onefreq; 
         %if &weight ~= %STR() %then %do;
             weight &weight;
         %end;
         where MISSING(&outvar) = 0;
      RUN; 

      data _onefreq;
         set _onefreq; 
          
         catelabel = catt(strip(vvalue(&outvar)), " N=", ROUND(COUNT,1));

         /* Concanate outcome name and category */
         outcat = "O&m" || strip(vvalue(&outvar));
      run;

      /* Get outcome variable categories and N into a macro var */
      proc sql noprint; 
         select catelabel
         into :catelabel&m separated by "*"
         from _onefreq; 
      quit;

      /* Get outcome variable names as they appear in the report data set */
      PROC TRANSPOSE DATA = _onefreq out=_tran (drop=_NAME_);
         id outcat;
      RUN;

      proc contents data=_tran out=_vname noprint;
      run;

      proc sort data=_vname; 
         by varNum;
      run;

      proc sql noprint; 
         select name into: categories&m separated by "  " 
         from _vname where name not in ("_LABEL_"); 
      quit;

      /* Get outcome variable label into a macro var */
      data _NULL_; 
         set _data (obs=1); 
         call symput("outcomevv&m", put(label(&outvar),$256.));
      run;

      %let m = %EVAL(&m+1);
      %let outvar = %SCAN(&outcome, &m);
   %end;

   %let num_out = %EVAL(&m-1);

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

   *------- build the table -----;
   PROC FORMAT;
      value _pval
          0-high=[PVALUE5.3]
          .='NA';
   RUN;

   OPTIONS ORIENTATION=&ORIENTATION MISSING = "-" NODATE;
   %if &doc = T %then %do;
      ODS rtf STYLE=tables FILE= "&OUTPATH.&FNAME &SYSDATE..doc"; 
   %end;

   PROC REPORT DATA=_report HEADLINE CENTER STYLE(REPORT)={JUST=CENTER} SPLIT='~' nowd 
     SPANROWS LS=256;
      COLUMNS &by _ordertyp _order Covariate  statistics level 
          %do i = 1 %to &num_out; 
                  ("&&outcomevv&i" '_____________________________'( &&categories&i)) par_p&i 
                    %if &NONPAR = T %then %do; nonpar_p&i %end;
         %end;; 

      %if &by ~= %STR() %then %do;
          DEFINE &by/order order=internal STYLE(COLUMN) = {JUST = L};
       %end;
       DEFINE _ordertyp/order order=internal noprint;
       DEFINE _order/order order=internal noprint;
       DEFINE Covariate/order order=data "Covariate" STYLE(COLUMN) = {JUST = L};
       DEFINE statistics/DISPLAY "Statistics" STYLE(COLUMN) = {JUST = C};
       DEFINE level/DISPLAY "Level" STYLE(COLUMN) = {JUST = C};

      %do j = 1 %to &num_out;
         %LET I = 1; 
           %DO %UNTIL (%SCAN(&&categories&j, &I) =   ); 
             %LET col = %SCAN(&&categories&j, &I); 
              %LET colabel = %SCAN(%BQUOTE(&&catelabel&j), &I, *); 
              DEFINE &col/DISPLAY  "&colabel" STYLE(COLUMN) = {JUST = C} ;
              %LET I = %EVAL(&I+1);
           %END; 

           %IF &NONPAR = T %THEN  %DO;
              DEFINE par_p&j/order missing "Parametric P-value*" STYLE(COLUMN) = {JUST = C} 
                    format=_PVAL.;
              DEFINE nonpar_p&j/order missing "Non-Parametric P-value**" STYLE(COLUMN) = {JUST = C}
                    format=_PVAL.;

               COMPUTE nonpar_p&j; 
                  IF . < nonpar_p&j <0.05 THEN 
                    CALL DEFINE("nonpar_p&j", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]");
               ENDCOMP; 
         %end;
           %else %if &NONPAR = A %then %do;
              DEFINE par_p&j/order missing "P-value*" STYLE(COLUMN)={JUST = C} format=_PVAL.;
         %end;
           %else %if &NONPAR = F %then %do;
              DEFINE par_p&j/order missing "Parametric P-value*" STYLE(COLUMN)={JUST=C} 
                    format=_PVAL.;
           %end;

           COMPUTE par_p&j; 
              IF . < par_p&j < 0.05 THEN 
                    CALL DEFINE("par_p&j", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]");
          ENDCOMP; 

      %end;

      compute after covariate; line ''; endcomp; 

      compute after _page_;
              
      %if &NONPAR = T %then %do;

          /* Both numerical and categorical covariates */
           %if &CLIST ~= %STR() and &NLIST ~= %STR() %then %do;
              line @0 "*  The parametric p-value is calculated by ANOVA for numerical covariates and &ctest";
              line @0 "   test for categorical covariates.";
              line @0 "** The non-parametric p-value is calculated by the Kruskal-Wallis test for";
              line @0 "   numerical covariates and &nonp_ctest test for categorical covariates.";  
           %end;
           /* Only categorical covariates */
           %if &CLIST ~= %STR() AND &NLIST = %STR() %then %do; 
              line @0 "*  The parametric p-value is calculated by &ctest test.";
              line @0 "** The non-parametric p-value is calculated by &nonp_ctest test.";  
           %end;
           /* Only numerical covariates */
           %if &CLIST = %STR() and &NLIST ~= %STR() %then %do;
              line @0 "*  The parametric p-value is calculated by ANOVA.";
              line @0 "** The non-parametric p-value is calculated by the Kruskal-Wallis test.";  
           %end;

      %end;    
      %else %if &NONPAR = A %then %do;

          /* Both numerical and categorical covariates */
           %if &CLIST ~= %STR() and &NLIST ~= %STR() %then %do;
              line @0 "*  The p-value is calculated by ANOVA for numerical covariates; and &ctest";
               line @0 "   test or &nonp_ctest for categorical covariates, where appropriate.";
           %end;
           /* Only categorical covariates */
           %if &CLIST ~= %STR() AND &NLIST = %STR() %then %do; 
              line @0 "*  The p-value is calculated by &ctest test or &nonp_ctest, where";
               line @0 "   appropriate.";
           %end;
           /* Only numerical covariates */
           %if &CLIST = %STR() and &NLIST ~= %STR() %then %do;
              line @0 "*  The p-value is calculated by ANOVA.";
               line @0 "   ";
           %end;
          
      %end;
       %else %if &NONPAR = F %then %do;
         %if &matchid = %STR() %then %do;
             /* Both numerical and categorical covariates */
              %if &CLIST ~= %STR() and &NLIST ~= %STR() %then %do;
                 line @0 "*  The parametric p-value is calculated by ANOVA for numerical covariates ";
                  line @0 "   and &ctest test for categorical covariates.";
              %end;
              /* Only categorical covariates */
              %if &CLIST ~= %STR() AND &NLIST = %STR() %then %do; 
                 line @0 "*  The parametric p-value is calculated by &ctest test.";
                  line @0 "   ";
              %end;
              /* Only numerical covariates */
              %if &CLIST = %STR() and &NLIST ~= %STR() %then %do;
                 line @0 "*  The parametric p-value is calculated by ANOVA. ";
                  line @0 "   ";
             %end;
         %end;
         /* Matched tests */
         %else %if &matchid ~= %STR() %then %do;
            /* Both numerical and categorical covariates */
            %if &CLIST ~= %STR() and &NLIST ~= %STR() %then %do;
                 line @0 "*  The parametric p-value is calculated by a paired t-test for numerical covariates, ";
                 line @0 "   McNemar's test for 2-level categorical covariates, and Bowker's test of symmetry for";
                 line @0 "   categorical covariates with more than 2 levels.";
            %end;
            /* Only categorical covariates */
            %if &CLIST ~= %STR() AND &NLIST = %STR() %then %do; 
               line @0 "*  The parametric p-value is calculated by McNemar's test for 2-level categorical covariates";
               line @0 "   and Bowker's test of symmetry for categorical covariates with more than 2 levels.";
            %end;
            /* Only numerical covariates */
            %if &CLIST = %STR() and &NLIST ~= %STR() %then %do;
                 line @0 "*  The parametric p-value is calculated by a paired t-test. ";
                  line @0 "   ";
            %end;
         %end;
       %end;

      ENDCOMP;
   RUN; 

   %if &doc = T %then %do;
      ODS RTF CLOSE; 
   %end;

   /* Reload original options that were in use before running the macro */
   /* This is mainly to reset the orientation */
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
