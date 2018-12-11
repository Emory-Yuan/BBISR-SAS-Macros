/**********************************************************************************************
* updates from last version. remove * from cvar. remove macro varaible CREFLIST; add stratum macro variable
* updated FORTH option - not using R to calculate CI.
***********************************************************************************************
Macro Name: UNI_PHREG
Created Date/Author: Feb. 2012/Yuan Liu
Last Update Date/Person: Oct, 2016/Yuan Liu
List of contributors: Dana Neckleach

Current Version: V26

Contact: Dr. Yuan Liu YLIU31@emory.edu, BBISR at Winship Cancer Institute and Emory University

Working Environment: SAS 9.4 English version

Purpose:  To conduct univariate survival analysis for each variable in the dataset. The hazard 
ratio with 95% CI is presented along with the log rank test p-value.  For categorical 
variables, the reference group will be shown along with the number of observation in each 
category.  The proportional hazard assumption can optionally be checked.Allow FORTH’s bias correction.

Notes: 1) The order of variables in the summary table is the same as the input order. For the 
best results, you may want to put the demographic variables together and also clinical 
characteristics variables together; 2) if the sample size is too big, the 
proportional hazard assumption check option may not work well.  In that case, please disable 
the option.3) 3) a list of variables with type3 p-value < PSELECT will be written to the log, 
which is stored in two global macro variables, &_select_cvar (list of categorical variables) 
and &_select_nvar (list of numerical variables) can be referred in the future multivariable modeling.
4)user should also refer to SAS help and document for technique details for some options.


Parameters: 

DATASET        The name of the data set to be analyzed.

EVENT          The variable name for time to event outcome.

START          Instead of a single failure-time variable (EVENT) a pair of failure time 
               variables can be specified.  This should be the start point of the interval for
               which the subject is at risk.  The label of this variable will be used in the
               table header.

STOP           Instead of a single failure-time variable (EVENT) a pair of failure time 
               variables can be specified.  This should be the end point of the interval for
               which the subject is at risk.  

CENSOR         The variable name for censoring indicator. If EVENTCODE not specified, it is required
			   that 1 is used for the event and 0 for censored cases. 

EVENTCODE      Specifies the number that represents the event of interest for the competing-risks 
			   analysis of Fine and Gray (1999).  

CLIST          List of categorical variables, generally separated by empty space. However, if need 
			   to change the reference level, you can follow each variable name by (DESC) or 
			   by (ref = “Ref level in formatted value”) where needed and separate terms by *.  
			   Also see code example.
 
NLIST          List of numerical variables, separated by empty space.

PHA            Value of F or T to indicate whether to run the proportional hazard assumption 
               check.

TYPE3          Value of F or T to indicate whether to print the type 3 Wald p-value in the 
               table (optional).  The default value is F.

LOGRANK        Value of F or T to indicate whether to print the log-rank p-value in the table 
               (optional).  The default value is T. This is calcualted by PROC LIFETEST separately .
 
FIRTH          Set to T to apply Firth’s modification for maximum likelihood estimation.See SAS Help and Documentation for 
			   technique details. The default value is F. 

WEIGHT         The name of the variable to include in the weight statement in PHREG (optional).
               Note that weights will not be normalized.  All N's reported will be 
               non-weighted.

ID             ID variable as ID statement in PHREG (optional). See SAS Help and Documentation for proper usage of this option.  

STRATA		   STRATA variable as STRATA statement in PHREG (optional).  

COVSAGG        Set to T to use the COVS(AGGREGATE) option in PHREG (optional). The default 
               value is F. 

PSELECT		   P-value threshold to select variable list as  type 3 p-vale < &PSELECT. The list is stored 
			   in two macro variable: &_select_var  and &_select_cvar, and  can be referred in future multivariable model.

ORIENTATION    Value of PORTRAIT or LANDSCAPE to indicate the paper layout of the report.

DOC            Set to T to create a RTF file containing the output or F to suppress creation of 
               the RTF file (optional).  The default value is T.

OUTPATH        Path for output table to be stored.

FNAME          File name for output table.

DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
For more details, please see the related documentation
**********************************************************************************************/

%MACRO UNI_PHREG(DATASET=, EVENT=, CENSOR=, START=, STOP=, EVENTCODE=, CLIST =,  NLIST=, DOC=T, 
   OUTPATH=, FNAME=, PHA=F, LOGRANK=T, TYPE3=F, FIRTH=F, WEIGHT=, ID=, STRATA = , COVSAGG=F, 
   ORIENTATION=PORTRAIT, PSELECT = , DEBUG=F); 

   %local I  WORK_SETS CHECK REFVAR EVENTVV LOGRANKP C assump cref format;
   %local cnt cnt_ pos pos_ L1 L2 clist_ _cnt_ clist_cnt i cnt_bef  cvar_ref cnt2 len2 cntw pos2 pos3 bef cnt_bef cvar_;

   %global _select_nvar  _select_cvar;

   /* Capitalize */
   %let pha = %UPCASE(&pha);
   %let debug = %UPCASE(&debug);
   %let logrank = %UPCASE(&logrank);
   %let doc = %UPCASE(&doc);
   %let covsagg = %UPCASE(&covsagg);
   %let type3 = %UPCASE(&type3);
   %let firth = %UPCASE(&firth);

   %IF &PHA = TRUE %then %let PHA = T;

   /*Set allowance of length of variable name*/
   %let length = 64;

   /* Count number of variables */

  %let nlist_cnt = 0;
  %if &nlist ~= %str( ) %then %let nlist_cnt = %sysfunc(countw(&NLIST,' '));

  %if %superq(clist) = %str( ) %then %let clist_cnt = 0;
  %if %superq(clist) ~= %str( ) %then %do;
  		
        %let cnt = %qsysfunc(countw(%superq(clist),'*'));
		%let cnt_= %qsysfunc(countw(%superq(clist),' '));

        %let pos=%qsysfunc(countc(%superq(clist), '(' ));
        %let pos_ = %qsysfunc(countc(%superq(clist), '*' ));

		%if %superq(pos) = 1 %then %do;
		%let L1=%qsysfunc(findc(%superq(clist),  '(' ));
	    %let L2=%qsysfunc(findc(%superq(clist),  ')' ));
		%let clist_=%qsysfunc(SUBSTR(%superq(clist), %superq(L1),%eval(%superq(L2) - %superq(L1) +1)));
	    %let _cnt_ = %qsysfunc(countw(%superq(clist_),' '));
        %end;

        %if %superq(cnt) = 1 and  %superq(cnt_) = 1 %then %let clist_cnt = 1;
		%else %if %superq(pos_) > 0 %then %let clist_cnt = %superq(cnt);
		%else %if %superq(pos_) = 0 and %superq(pos) = 0 %then %let clist_cnt = %superq(cnt_);
		%else %if (%superq(pos) = 1 and %superq(pos_) = 0) and %superq(cnt_) = %superq(_cnt_) %then  %let clist_cnt = 1;
		%else %do; %put ERROR: The categorical variables in CLIST should be separated by * if you specify the reference level.; 
          		   %goto exit;%end;

  %end;

/*Build up categorical variable list without reference level specified*/
  %let cvarlist =; 

  %if &clist_cnt = 1 and &pos = 1 %then %let cvarlist = %sysfunc(SUBSTR(&clist, 1,%eval(&L1-1)));
  %else %if &clist_cnt = 1 and &pos = 0 %then %let cvarlist = &clist;

  %if &clist_cnt > 1  and &pos_ = 0 %then  %let cvarlist = &clist ; 

  %if &clist_cnt > 1  and &pos_ > 0 %then %do; 
 			
      %do i = 1 %to &clist_cnt; 

				%let cvar_ref = %SCAN(&CLIST, &i, '*'); 

				%let cnt2=%qsysfunc(countc(&cvar_ref, '(' ));
				
				%let len2=%LENGTH(&cvar_ref);
				%let cntw=%qsysfunc(countw(&cvar_ref,' '));

				%put cvar_ref = &cvar_ref cnt2= &cnt2  len2=&len2 cntw =&cntw ;

				%if %superq(cnt2) = 1 %then %do;

	                %let pos2=%sysfunc(findc(&cvar_ref, '(' ));
					%let pos3=%sysfunc(findc(&cvar_ref, ')' ));

					%let bef = %sysfunc(SUBSTR(&cvar_ref, 1,%eval(&pos2-1)));
			        %let cnt_bef = %qsysfunc(countw(&bef,' '));

					%if &cnt_bef = 1 and  &pos3 = &len2 %then %let cvar_=&bef;
					%else %do; %put ERROR: The categorical variables in CLIST should be separated by * if you specify the reference level.; 
          		                %goto exit;%end;
					%end;
					
                %if &cnt2 = 0 and &cntw  = 1 %then %LET cvar_ = &cvar_ref;

				%if &cnt2 > 1 or ( &cntw > 1 and &cnt2 = 0) %then %do; %put ERROR: The categorical variables in CLIST should be separated by * if you specify the reference level.; 
          		%goto exit;  %end;	 

			%let cvarlist = &cvarlist &cvar_; 
 %end;
		 
 %end;

/* Check categorical variables that less than two non-missing levels; update varaible length */
/* %if &nlist_cnt > 0 %then %do;
 	%do i = 1 %to &nlist_cnt; 
	   %let nvar_ = %SCAN(%superq(nlist), &i, ' '); 
	    %If %LENGTH(%superq(nvar_)) > &length %then %do;
         %put ERROR: Variable name %superq(nvar_) must be less than &length characters.; 
           %goto exit; %end;
%end;
%end;*/


 %if &nlist_cnt > 0  %then %do; 
	%do i = 1 %to &nlist_cnt; 
	   %let nvar_ = %SCAN(%superq(nlist), &i, ' '); 
	    %If %LENGTH(&nvar_) > &length %then %let length = %LENGTH(&nvar_);
 %end; * end of  if nlist_cnt > 0;
 %end; * end of do i= 1 to nlist_cnt;
  
 %if &clist_cnt > 0  %then %do; 
	%do i = 1 %to &clist_cnt; 
	   %let cvar_ = %SCAN(%superq(cvarlist), &i, ' '); 

	    %If %LENGTH(&cvar_) > &length %then %let length = %LENGTH(&cvar_);

         PROC SQL noprint;
         select count(distinct &cvar_) into :check
         from &DATASET
         where MISSING(&cvar_) = 0;
         QUIT;

         %if %superq(check) <= 1 %then %do;
         %put ERROR: The variable %SCAN(&cvar_, &i) has less than two non-missing levels.  Please remove from CLIST.; 
          %goto exit; %end;
	
    %end;* end of  do i = 1 to &clist_cnt;;
  %end;* end of if &clist_cnt > 0 ;

   /* Can't use weights with the log-rank test */
   %if &weight ~= %STR() and &logrank = T %then %do;
      %put ERROR: Weights cannot be used with the log-rank test.  Please set LOGRANk=F.; 
       %goto exit;
   %end;

   /* Can't use weights with the PH test */
   %if &weight ~= %STR() and &pha = T %then %do;
      %put ERROR: The PH assumption cannot be tested when weights are used.  Please set PHA=F.; 
      %goto exit;
   %end;

   /* Should not specify both EVENT and START/STOP */
   %if &event ~= %STR() and &start ~= %STR() and &stop ~= %STR() %then %do;
      %put ERROR: Either EVENT and START/STOP, but not both should be specified.; 
      %goto exit;
   %end;

   /* Can't use log-rank test with start/stop format */
   %if &start ~= %STR() and &stop ~= %STR() and &logrank=T %then %do;
      %put ERROR: Log-rank p-values cannot be reported when a START/STOP time is specified.  Please set LOGRANK=F.; 
     %goto exit;
   %end;

   /* Can't use PL method with start/stop format */
   %if (&start ~= %STR() OR &stop ~= %STR()) and &firth=T %then %do;
      %put ERROR: Penalized likelihood cannot be used when a START/STOP time is specified.  Please set PL=F.; 
      %goto exit;
   %end;

   /* Can't use PL method with COVSAGG */
   %if &COVSAGG = T and &firth=T %then %do;
      %put ERROR: Penalized likelihood cannot be used in conjunction with COVSAGG.  Please set PL=F.; 
      %goto exit;
   %end;

   /*Can't use EVENTCODE with START and STOP, as GRARY-FINE model cannot take counting process data format in current PROC PHREG of SAS 9.4 */
	%if &EVENTCODE ~= %STR() and (&START ~=%str() or &logrank = T or &stop ~= %STR() ) %then %do;
      %put ERROR: EVENTCODE was set, which means you intent to fit the data into a competing risk model by the Gray-fine model, in current version of
PROC PHREG in SAS 9.4, it does not allow couting process data format. You may need to specify EVENT and CENSOR parameter in this macro instead.
Also you need to turn off LOGRANK and set it as F; 
     %goto exit;
   %end;


   /* Get list of data sets in work library to avoid deletion later */
   ods select none;
   ODS OUTPUT Members(nowarn)=_DataSetList;
   PROC DATASETS lib=work memtype=data ;   
   QUIT;
   ods select all;


   /* If there are data sets in the work library */
   %if %sysfunc(exist(_DataSetList)) %then %do;
      PROC SQL noprint;
         select Name
         into :work_sets separated by ' '
         from _DataSetList;
      quit;
   %end;
   %else %do;
      %let work_sets =;
   %end;

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

 
   *CHARACTER VARIABLES ;

   %IF &clist_cnt > 0 %THEN %DO; 

         /* Oct 2016, this option was turned down, if the data is cluster, ID statement will be used. The number of sample size would represent # of person time.*/
         /*For datasets where records does not represent the number of observations */
         /* Return to observation level for frequency calculations */
         /*%if &id ~= %STR() %then %do;
            PROC SORT DATA =&dataset out = _uni nodupkey;
               by &id;
            RUN;
         %end;
         %else %do;
            DATA _uni;
               set &dataset;
            RUN;
         %end;*/


   	%do C = 1 %to &clist_cnt;
   		%if &pos > 0 %then %let cvar_ref = %SCAN(%superq(CLIST), &C, '*'); 
		%else  %LET cvar_ref = %SCAN(%superq(CLIST), &C, ' ');
		%let cvar = %SCAN(%superq(cvarlist), &C, ' ');

		
        ods select none;
        ODS OUTPUT 'One-Way Frequencies' = _FREQ;   
    	PROC FREQ DATA=&dataset ORDER=internal ; 
			TABLE %superq(cvar);
            %if &event ~= %STR() %then %do;
               WHERE &EVENT ~= . and &censor ~= .;
            %end;
            %else %do;
               WHERE &START ~= . and &stop ~= . and &censor ~= .;
            %end;
         RUN;  
		 ods select all;
		 
         DATA _FREQC; 
            SET _FREQ ; 
            LENGTH Covariate Name $256. Level $256.;
            Covariate = label(%superq(cvar)); 
			Name = "&cvar";
            level = %sysfunc(TRIM(LEFT(F_%superq(cvar))));
            N =  Frequency;
            KEEP Covariate Name level Frequency N; 
         RUN;

         %if &logrank = T %then %do; 
            /* Calculate log-rank test */
		    ods select none;
            ODS OUTPUT HomTests=_score;
            PROC LIFETEST DATA = &dataset NOTABLE;
               time &EVENT*&CENSOR(0);
               strata %superq(cvar);
            RUN;
			ods select all;

            data _score;
               set _score;
               format ProbChiSq;
            run;

            /* Get log-rank p-value */
            proc sql noprint; 
               select ProbChiSq into: logrankp 
               from _score where test ="Log-Rank";
            quit;
         %end;
         %else %let logrankp=.;

         /* Run this even when penalized likelihood is used in order to get likelihood ratio 
         test because it is difficult to pull from R */
		 ODS select none;
         ODS OUTPUT %IF &PHA = T %then %do; "Supremum Test for Proportional Hazards Assumption" = _phac %end;
               "Maximum Likelihood Estimates of Model Parameters" = _mlec 
               "Number of observations" = _nobs Type3 = _typ3  GlobalTests=_global;
          proc phreg data=&dataset namelen=&length %if &covsagg=T %then %do; COVS(AGGREGATE) %end;;
           class &cvar_ref/order=internal %if &firth ~= T %then %do; param=glm %end;;
            %if &event ~= %STR() %then %do;
               model &EVENT*&CENSOR(0) =&CVAR/ %if &firth = T %then %do; firth risklimits=pl %end;%else %do; rl %end;
			   								%if &EVENTCODE ~=%str() %then %do; eventcode = &eventcode %end; ;
            %end;
            %else %do;
               model (&start,&stop)*&CENSOR(0) = &CVAR/%if &firth = T %then %do; firth risklimits=pl; %end;%else %do; rl; %end;
            %end;
            %if &PHA = T %then %do; 
               assess ph/resample seed=1;
            %end;
            %if &weight ~= %STR() %then %do;
               weight &weight;
            %end;
            %if &id ~= %STR() %then %do;
               id &id;
			%end;
			%if &strata ~= %STR() %then %do;
               strata &strata;
			%end;
         Run;
		 ODS select all;

         data _mlec;
            set _mlec;
            /* Preserve order */
            order2 = _n_;
            format ProbChiSq;
			%if &firth = T %then %do;
			rename HRLowerPLCL =HRLowerCL HRUpperPLCL=HRUpperCL;
			%end;
         run;

         %IF &PHA = T %THEN %DO;
              
            /* Merge on PHA check */
            DATA _freqc;
               if _n_ = 1 then set _phac;
               set _freqc;
               format pvalue;
            RUN;

         %end;

         PROC SORT DATA = _freqc;
            by level;
         RUN;
         PROC SORT DATA =  _mlec (rename=(ClassVal0=level));
            by level;
         RUN;

         data _freqcc; 
            merge _freqc _mlec;
            by level;
            /* Make sure reference level comes last */
            if order2 = . then order2 = 999;
         run;

		 /*proc sql;
		 create table _freqcc as
		 select a.*, b.N
		 from _mlec as a left join _freqc as b
		 on a.ClassVal0 = b.level;quit;*/


         /* Return to order */
         PROC SORT DATA = _freqcc;
            by order2; 
         RUN;

         data _freq&C; 
            /* Merge on type3 p-value */
            if _n_ = 1 then set _typ3 (RENAME=(ProbChiSq=p_type3) KEEP=ProbChiSq);
            set _freqcc;
            logrankp= &logrankp;
            /* Preserve order */
            order2 = &c;
            keep order2 Covariate Name Level Frequency logrankp HazardRatio HRLowerCL HRUpperCL 
               ProbChiSq p_type3 %IF &PHA = T %then %do; pvalue %end;;
         run;
               
      %END; 


      Data _freq_all; 
         set _FREQ1 - _FREQ%eval(&clist_cnt);
           order = 1;
      run;
   %END;

   *NUMERIC VARIABLES ;
   %IF &NLIST NE  %THEN %DO; 
          
      %LET N = 1; 

      %DO %UNTIL (%SCAN(&NLIST, &N) =   ); 
         %LET NVAR = %SCAN(&NLIST, &N); 

		    ODS select none;
            ODS OUTPUT 
             %IF &PHA = T %then %do; "Supremum Test for Proportional Hazards Assumption" = _phac %end;
             "Maximum Likelihood Estimates of Model Parameters" = _mlec 
             "Number of observations" = _nobs;
            proc phreg data=&DATASET namelen=&length  %if &covsagg=T %then %do; COVS(AGGREGATE) %end;;
               %if &event ~= %STR() %then %do;
                  model &EVENT*&CENSOR(0) = &NVAR/ %if &firth = T %then %do; firth risklimits=pl %end;%else %do; rl %end;
				  							%if &EVENTCODE ~=%str() %then %do; EVENTCODE = &EVENTCODE %end; ;
               %end;
               %else %do;
                  model (&start,&stop)*&CENSOR(0) = &NVAR/%if &firth = T %then %do; firth risklimits=pl; %end;%else %do; rl; %end;
               %end;
               %IF &PHA = T %then %do; assess ph/resample seed=1; %end;
               %if &weight ~= %STR() %then %do;
                  weight &weight;
               %end;
               %if &id ~= %STR() %then %do;
                  id &id;
               %end;
			   	%if &strata ~= %STR() %then %do;
               strata &strata;
			   %end;

            run;
			ODS select all;

            data _mlec;
               set _mlec; 
               format  ProbChiSq;
			 %if &firth = T %then %do;
			rename HRLowerPLCL =HRLowerCL HRUpperPLCL=HRUpperCL;
			%end;

            run;

            %IF &PHA = T %then %do;
               data _phac;
                  set _phac; 
                  format pvalue;
               run;
            %end;

            proc sql noprint;
               select HazardRatio into: hr from _mlec;
               select ProbChiSq into: hrp from _mlec;

   			   select HRLowerCL into: hr_lb from _mlec;
               select HRUpperCL into: hr_ub from _mlec;

			    
             %IF &PHA = T %then %do;
                  select pvalue into: assump from _phac;
               %end;

				/*sample size when ID is active.*/
               %if &id = %STR() %then %do;
                  select SumFreqsUsed into: N_obs from _nobs;
               %end;
               %else %do;
                   select SumFreqsUsed into: N_obs from _nobs;


               %if &event ~= %STR() %then %do;
                     WHERE &EVENT ~= . and &censor ~= . and &nvar ~= .;
                  %end;
                  %else %do;
                     WHERE &START ~= . and &stop ~= . and &censor ~= . and &nvar ~= .;
                  %end; 
              %end;
            quit;

         data _NULL_; 
            set &dataset; 
            call symput('VV', put(label(&NVAR),$256.));
         run;

         data _summary&N;
            length covariate name $256. Level $71.;
            covariate = "&VV";
			Name = "&NVAR";
            level = ' ';
            Frequency = &N_obs;
            /* Don't report log-rank p-value for numeric variables */
            logrankp =.;
            HazardRatio = &hr;
            HRLowerCL = &hr_lb;
            HRUpperCL = &hr_ub;
            ProbChiSq = &hrp;
            /* Type 3 p-value will be the same as the HR p-value */
            p_type3 = probchisq;
            %IF &PHA = T %THEN %DO; 
               pvalue = &assump; 
            %end;
            /* Preserve order */
            order2 = &n;

            keep order2 Covariate name level Frequency logrankp HazardRatio HRLowerCL HRUpperCL 
               ProbChiSq p_type3 %IF &PHA = T %THEN %DO; pvalue %end;;
         run;

         %LET N = %EVAL(&N+1);
      %END; 

      %LET N = %EVAL(&N-1);

      DATA _summary_all; 
         SET _summary1-_summary&N; 
         order = 2;
      RUN; 
   %END;

   /* Combine categorical and numerical covariate results */
   DATA _report; 
      set %IF &clist_cnt >0 %then %do; _freq_all %end;
          %if &nlist_cnt >0 %then %do; _summary_all %end;;
      /* HR and 95% CI */
      if HazardRatio ~= . then HR = TRIM(LEFT(PUT(HazardRatio,8.2))) || " (" || 
         TRIM(LEFT(PUT(HRLowerCL,8.2))) || "-" || 
         TRIM(LEFT(PUT(HRUpperCL,8.2))) || ")";
      else HR = '-';
   run; 

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

   **** PRINT THE TABLE ************;
   data _NULL_; 
      set &dataset; 
      %if &event ~= %STR() %then %do;
         call symput('eventvv', put(label(&EVENT),$256.));
      %end;
      %else %do;
         call symput('eventvv', put(label(&start),$256.));
      %end;
   run;

   OPTIONS ORIENTATION=&ORIENTATION MISSING = "-" NODATE;
   %if &doc = T %then %do;
      ODS RTF STYLE=TABLES file= "&OUTPATH.&FNAME &SYSDATE..DOC"; 
   %end;
      PROC REPORT DATA=_report HEADLINE HEADSKIP CENTER SPANROWS LS=256 
          STYLE(REPORT)={JUST=CENTER} SPLIT='~' nowd; 
         COLUMNS order order2 Covariate  Level frequency  
               ("&eventvv" '------------------------------------------'( HR
                ProbChiSq 
               %IF &PHA = T %then %do; pvalue %end; 
               %if &type3 = T %then %do; p_type3 %end;
               %if &logrank = T %then %do; logrankp %end;)); 
         DEFINE order/order order=internal noprint;
         DEFINE order2/order order=internal noprint;
         DEFINE Covariate/ order order=data  "Covariate"  STYLE(COLUMN) = {JUST = L}; 
         DEFINE Level/ DISPLAY   "Level"     STYLE(COLUMN) = {JUST =L}; 
         DEFINE frequency/ DISPLAY   "N"     STYLE(COLUMN) = {JUST = C};                
         DEFINE HR/DISPLAY "Hazard Ratio (95% CI)" STYLE(COLUMN) = {JUST = C CellWidth=15%};
         DEFINE ProbChiSq/DISPLAY "HR P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%}  FORMAT=PVALUE5.3; 
         %IF &PHA = T %then %do;
            DEFINE pvalue/order  "Assumption P-value"  STYLE(COLUMN) = {JUST = C} FORMAT=PVALUE5.3; 
            COMPUTE pvalue; 
               IF pvalue <0.05 THEN CALL DEFINE("pvalue", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
            ENDCOMP; 
         %end;
         %if &type3 = T %then %do;
            DEFINE p_type3/ order "Type3 P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%}  FORMAT=PVALUE5.3 MISSING; 
            COMPUTE p_type3; 
              IF p_type3 <0.05 THEN CALL DEFINE("p_type3", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
            ENDCOMP; 
         %end;
         %if &logrank = T %then %do;
            DEFINE logrankp/ order "Log-rank P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%}  FORMAT=PVALUE5.3 MISSING;
            COMPUTE logrankp; 
               IF logrankp <0.05 THEN CALL DEFINE("logrankp", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
            ENDCOMP;  
         %end;

         COMPUTE ProbChiSq; 
            IF ProbChiSq <0.05 THEN CALL DEFINE("ProbChiSq", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         ENDCOMP; 

         compute after covariate; line ''; endcomp;               

         compute after _page_;
            /* If Firth is used */
            %if &firth = T %then %do;
               line @0 "Firth’s penalized maximum likelihood estimation was used.";
            %end;
			%if &strata ~=%str( )%then %do;
               line @0 "Analysis was stratified by variable: &strata.";
            %end;
			%if &weight ~=%str( )%then %do;
               line @0 "Analysis was weighted by variable: &weight.";
            %end;
			%if &id ~=%str( )%then %do;
               line @0 "Analysis was taken the clustering effect within &id into account, and N represented number of &id-times.";
            %end;
			%if &eventcode ~=%str( )%then %do;
               line @0 "A competing risk analysis of Fine and Gray (1999) was fitted with the event of interest being &eventcode.";
            %end;

         endcomp;
      RUN;

   %if &doc = T %then %do;
      ods rtf close;
   %end;


	/** extract variable name with type 3 p-value < PSELECT in univariate analysis;*/
%if &PSELECT ~= %str() %then %do;
   proc sort data=_report out = _report_; by name; run;
   data _report_; set _report_; by name; if first.name;run;
        %let _select_cvar = ;
   		%let _select_nvar = ;		

		proc sql noprint;
		select  name into :_select_cvar separated by " "
		from _report_
		where Level ~= "" and p_type3 < &PSELECT
		order by order2;

		select  name into :_select_nvar separated by " "
		from _report_
		where Level = "" and p_type3 < &PSELECT
		order by order2;
		quit;

   %put Categorical variables selected as type-3 p-value < &PSELECT: &_select_cvar;
   %put Numerical variables selected as type-3 p-value < &PSELECT: &_select_nvar;
%end;

      /* Reload original options that were in use before running the macro */
   /* This is mainly to reset the orientation */
   PROC OPTLOAD data=_options;
   RUN;

   /* If not in debug mode then delete temporary data sets */
   %if &debug = F %then %do;
      /* If there are work data sets that should not be deleted */
      %if %sysevalf(%superq(work_sets)~=,boolean) %then %do;
         /* DELETE ALL TEMPORARY DATASETS that were created */
         proc datasets lib=work memtype=data nolist;  
            save &work_sets;
          quit;  
      %end;
      %else %do;
         proc datasets lib=work kill memtype=data nolist;  
         quit; 
      %end;
   %end;

%exit:

%MEND; 

