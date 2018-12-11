/**********************************************************************************************
***********************************************************************************************
Macro: MULTIPLE_LOGREG
Created Date/Author: June 26, 2012/Dana Nickleach
Last Update Date/Person: Mar, 2016/Yaqi Jia
Current Version: V16
Working Environment: SAS 9.4 English version

Contact: Dr. Yuan Liu yliu31@emory.edu 

Purpose:  To produce a summary table in a Microsoft Word Document from a multivariable logistic
regression model.  The modeling must be done before using this macro.  

Notes: The logistic regression must be conducted using PROC LOGISTIC before calling this macro.
A binary logit, cumulative logit, or generalized logit model can be used.  Four tables must be 
created using ODS output when running PROC LOGISTIC (nobs, type3, estimate, ci, modelinf, and
resp) as in the example below.  The options PARAM=GLM and CLPARM=WALD must also be used in the
PROC.  Variable names must not be more than 20 characters.  This macro can be set up to handle
models including interaction terms when the outcome is binary, and in the model fitting step, 
use SLICE statement to specifiy how the strafied analysis to be conducted (ODS OUTPUT NObs = nobs 
Type3 = type3 ParameterEstimates = estimate CLparmWald = CI ModelInfo=modelinf ResponseProfile=resp 
SliceDiffs = slices_diff).   

Parameters: 

EVENT          The event category for the binary response model (optional).  Specify the value 
               without quotes.  This will appear in the table header.  If using a clogit or 
               glogit model, leave this blank.  
OUTPATH        File path for output table to be stored.  
FNAME          File name for output table.  
FOOTNOTE       Text of footnote to include in the table (optional).  It should be in quotes.  
               This footnote will appear below the footnote containing the number of 
               observations.  Leave this field blank if not including an additional footnote.
EFFECT         If not empty, then the model contains a two-way interaction. It is used to specify
			   the treatment effect in the interaction.
SLICEBY		   Use in combine with EFFECT to specify the stratified variable in the interaction.
SHORTREPORT	   Use in cobmine with EFFECT and SLICEBY when there is an interaction in the model and set
			   to T to only report the stratified treatment effect.	
TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.
CLNUM		   Set to T if you want to see the number of observations for each level of covariates. 
			   The default is T.
ORIENTATION	   orientation of the output Word table. Default is portrait, can be changed to landscape.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
* For more details, please see the related documentation.
**********************************************************************************************/

%macro MULTIPLE_LOGREG(EVENT=,OUTPATH=,FNAME=,FOOTNOTE=,TYPE3=T,EFFECT=,SLICEBY=,clnum=T,shortreport=T,
ORIENTATION = portrait,DEBUG=F);

   %local classexist dsn model ref headlab __Macro_Err length length1 length2 ;

   %let debug = %UPCASE(&debug);
   %let type3 = %UPCASE(&type3);
   %let clnum= %UPCASE (&clnum);

   /* Initialize error variable */
   %let __Macro_Err = 0;

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

   /* Calculate number of observations read and used */
   PROC SQL noprint;
      select NObsRead, NObsUsed 
      into :nobs_read, :nobs_used
      from nobs;
   QUIT;

   /* Get data set name and outcome name */
   PROC SQL noprint;
      select value into :dsn
       from modelinf
       where description = 'Data Set';
       select value into :outName
       from modelinf
       where description = 'Response Variable';
       /* Model Type */
       select value into :model
       from modelinf
       where description = 'Model';
   QUIT;

   /* Make sure EVENT was not specified if glogit or clogit was used */
   %if &model ~= binary logit and %sysevalf(%superq(EVENT)~=,boolean) %then %do;
      %put ERROR: The EVENT parameter should only be specified for binary outcomes.; 
       %let __Macro_Err=1;
   %end; 

   /* If there is an error in the parameters supplied then exit */
   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   /* If using glogit get reference group */
   %if &model=generalized logit %then %do;
      PROC SQL noprint;
         select a.outcome into :ref
         from resp a
         where a.count ~= . and LEFT(a.outcome) not in (select distinct LEFT(response) from estimate);
      QUIT;
   %end;

   /* Get outcome label */
   DATA _NULL_;
      set &dsn (obs=1);
      call symput('outLab', VLABEL(&outname));
   RUN;

   /* Table header */
   /* If glogit model specify reference group in header */
   %if &model=generalized logit %then %let headlab = %TRIM(%BQUOTE(&outlab))=%TRIM(%LEFT(%BQUOTE(&ref))) as the Reference;
   /* If an event was specified then add it to the outcome label */
   %else %if %sysevalf(%superq(EVENT)~=,boolean) %then %let headlab = %TRIM(%BQUOTE(&outlab))=&event;
   %else %let headlab = &outlab;

   /* Get variable labels */
   PROC CONTENTS DATA = &dsn out=_cont (rename=(name=parameter)) noprint;
   RUN;

   /* If model is binary or glogit the classval0 variable will not exit if class variables */
   /* were not used */
   %if &model ~= cumulative logit %then %do;
      /* Check to see if class variables were used - s/b classval0 variable */
      data _null_;
         dsid=open('estimate');
         check=varnum(dsid,'classval0');
         /* 0 = class variables not used, 1 = class variables used */
         exist = MIN(check, 1);
         CALL SYMPUT('classExist',PUT(exist,1.));
      run;
   %end;
   /* For cumulative logit models the classval0 will always exist, but it will be blank except */
   /* for the intercept rows if no class variables were used */
   %else %do;
      PROC SQL;
       select count(*) into: classExist
       from estimate
       where Variable ~= 'Intercept' and classval0 ~= ' ';
    QUIT;

    /* 0 = class variables not used, 1 = class variables used */
    %let classExist = %SYSFUNC(MIN(&classExist,1));
   %end;

    %PUT number of class variable = &classexist;
   /* Preserve original variable order */
   DATA _est2;
      length parameter $32;
      set estimate (rename=(variable=parameter));
    /* Drop intercept */
    where parameter ~= 'Intercept';
      order = _n_;
   RUN;

   PROC SORT DATA = _est2;
      by parameter %if &classexist = 1 %then %do; ClassVal0 %end;;
   RUN;
   PROC SORT DATA =ci out=_ci2 ;
      by parameter %if &classexist = 1 %then %do; ClassVal0 %end;;
    /* Drop intercept */
    where parameter ~= 'Intercept';
   RUN;

   /* For glogit models need to find length of response variable before merge to prevent errors */
   %if &model=generalized logit %then %do;

      /* Get length of response variable in estimate file */
      DATA _NULL_;
         set _est2 (obs=1);
         CALL SYMPUT('length1',PUT(VLENGTH(response),best12.));
      RUN;

      /* Get length of response variable in CI file */
      DATA _NULL_;
         set _ci2 (obs=1);
         CALL SYMPUT('length2',PUT(VLENGTH(response),best12.));
      RUN;

      /* Find maximum length  */
      %let length = %sysfunc(MAX(&length1,&length2));

   %end;

   /* Merge estimates and CI */
   DATA _est2 %if &classexist = 1 %then %do; (rename=(ClassVal0=level)) %end;;
      %if &model=generalized logit %then %do; 
         /* Set length of response variable to maximum */
         length response $&length;
      %end;
	  %if &classexist = 1 %then %do; length classval0 $256.; %end;
      merge _est2 _ci2 (drop=estimate);
      by parameter %if &classexist = 1 %then %do; ClassVal0 %end; 
          %if &model=generalized logit %then %do; Response %end;;
      /* Calculate OR and CI */
      if estimate ~= 0 then do;
         OR = exp(estimate);
         LCL = exp(LowerCL);
         UCL = exp(UpperCL);
         /* OR and 95% CI */
         OR_CI = TRIM(LEFT(PUT(OR,8.2))) || " (" || 
         TRIM(LEFT(PUT(LCL,8.2))) || "-" || 
         TRIM(LEFT(PUT(UCL,8.2))) || ")";
      end;
      else OR_CI = '-';
   RUN;

   PROC SORT DATA = _cont;
      by parameter;
   RUN;

   /* If class variables were used or if generalized logit model was used */
   /* A glogit model will always have type3 p-values */
   %if &classexist = 1 OR &model=generalized logit %then %do;
      /* Prepare for merge */
      data typ3_;
         length parameter $32;
         set type3 (rename=(ProbChiSq = p_type3 effect=parameter)); 
      run;

      proc sort data=typ3_;
         by parameter;
      run;

      /* Merge on type III p-values */
      data _mg; 
         length plabel $256.;
         merge typ3_ (in=a keep=parameter p_type3) _est2 (in=b) _cont (keep=parameter label); 
         by parameter;
         if a or b;
         if label ~= ' ' then plabel = label;
         else plabel = parameter;

       /* For glogit models if class variables weren't used then created an empty level var */
       %if &classexist = 0 %then %do;
          Level = "";
       %end;

         format ProbChiSq p_type3;
      run; 

   %end;
   %else %do;

      /* Merge on labels, but not type III p-values */
      data _mg; 
         length plabel $256.;
         merge _est2 (in=b) _cont (keep=parameter label); 
         by parameter;
         if b;
         if label ~= ' ' then plabel = label;
         else plabel = parameter;

         /* Type 3 p-value */
         p_type3 = ProbChiSq;
         /* No level for continous variables */
         Level = "";

         format ProbChiSq p_type3;
      run; 

   %end;


   %if &clnum=T %then %do;
	    
       %if %sysfunc(exist(clfreq)) %then %do;
	  /*prepare to be sorted. there is missing value in class variable in clfreq, so need below data step*/
	    data clfreq %if &classexist = 1 %then %do;(rename=(value=level))%end;; 
       length parameter $256.;
       retain parameter;
       set clfreq;
       if class ~= ' ' then parameter=class;
	   else class=parameter;
       run;

       proc sort data=clfreq;by parameter %if &classexist = 1 %then %do;level %end;;run;
	   proc sort data=_mg;  by parameter %if &classexist = 1 %then %do;level %end;;run;

       data _mg;
	     length parameter $256. %if &classexist = 1 %then %do;level $256. %end; n $120.;
	     merge _mg(in=a) clfreq(in=b);
	     by parameter %if &classexist = 1 %then %do;level %end;;
		 if a;
		 n=put(total,8.);
		 if n=. then n="&nobs_used";
	   run;
      %end;
	  %else %do;
	    data _mg;
		  length n $120.;
		  set _mg;
		  n="&nobs_used";
		run;
	  %end;*end of exist else condition;
    %end;
   /* Return to original order */
   proc sort data=_mg;
      by order;
   run; 

   /* If reporting contrasts as well */
   %if &effect ~=%str() %then %do;

	/*Extract variable name,label and length*/

			/*Extract variable name,label and length*/
       %let cfivar=%sysfunc(countw(&_finalvar,' '));
        %put finalvar=&_finalvar;
        %put number of final var=&cfivar; 
        data frequency;
		  set &dsn;
		  %do i=4 %to &cfivar;
		  if missing(%scan(&_finalvar,&i,' '))=0;
		  %end;

		  where missing(&outcome)=0;
		run;
        ods select none;
		ods output list=freq_;
        proc freq data=frequency;
		   tables &sliceby*&effect/list;
	    run;
		ods select all;

	    Data freq;
	      length slice $256. effect $256.;
	      set freq_;
	      sliceby=vvalue(&sliceby);
	      effect=vvalue(&effect);
	      keep sliceby effect frequency;
	    run;

      /*extract effect and slicby values from data to prepare for the merge with freq*/ 
        data slices_diff;
	      length effect $256. _effect $256. slice $256.;
	      set slices_diff;
          effect=vvalue(&effect);
		  _effect=vvalue(_&effect);
		  sliceby=vvalue(slice);
		  rank=_n_;/*save original order*/
        run;
       
	  /* get the number of observations used for &effect */
        proc sort data=slices_diff; by effect sliceby ;run;
        proc sort data=freq; by effect sliceby;run;

        data merge_1;
           merge freq (in=a rename=(frequency=obsused1)) slices_diff(in=b);
           by effect;
           if b;
        run;

       /* get the number of observations used for _&effect*/
       proc sort data=slices_diff nodupkey; by _effect sliceby;run;

       data merge_2;
          merge freq (in=a rename=(effect=_effect frequency=obsused2)) slices_diff (in=b);
          by _effect;
          if b;
          keep sliceby _effect obsused2;
       run;

       proc sort data=merge_1;by _effect sliceby;run;
       proc sort data=merge_2;by _effect sliceby;run;
       data merge_12;
         merge merge_1(in=a) merge_2 (in=b);
         by _effect sliceby;
		 if a;
       run;

       
	  DATA merge12;
        length plabel $256. level $256. n $120.; 
        set merge_12;
		plabel = tranwrd(sliceby,scan(sliceby,1), "");
		level=cat(strip(effect)," vs. ",strip(_effect));
        n=cat(strip(obsused1)," vs. ",strip(obsused2));
		OR = ExpEstimate;
		LCL = LowerExp;
		UCL = UpperExp;
		ProbChiSq = Probz;
		 /* OR and 95% CI */
        if OR ~= . then OR_CI = TRIM(LEFT(PUT(OR,8.2))) || " (" || 
         TRIM(LEFT(PUT(LCL,8.2))) || "-" || 
         TRIM(LEFT(PUT(UCL,8.2))) || ")";
        else OR_CI = "-";
	  Run;
    
	   /* put data in original order*/
		proc sort data=merge12;by rank;run;

	  /* extract type3 p-value as interaction p-value from type3 data set*/
      
	  Data _null_; /* get the label of sliceby variable*/
	     set &dsn (obs=1);
		 call symput ('slicebylabel',vlabel(&sliceby));
		 call symput ('effectlabel',vlabel(&effect));
	  run;

	  Data _type3;
		 length plabel $256 level $256.;
		 set type3 (rename=(ProbChiSq = p_type3 effect=parameter));
         plabel = "Comparisons Stratified by &Slicebylabel :";
		 level="&effectlabel :";
		 where index(parameter,'*') ~=0;
		 keep parameter plabel level p_type3;
	  run;
	  
	  /* Combine Type3 data set and Hazardratio data set*/
      DATA _contrast;
         set _type3 merge12;
         keep plabel level n OR OR_CI probchisq p_type3;
      RUN;

	  %if &shortreport = T %then %do;

	     %let upeffect=%upcase(&effect); /* delete &effect &sliceby and &effect*&sliceby variables in above data set*/
		 %let upsliceby=%upcase(&sliceby);

	      DATA conlab;
          set _mg; 
		  if upcase(parameter) not in ("&upeffect"  "&upsliceby" "&upeffect*&upsliceby" );
          RUN;

		  proc sql noprint;
		    select distinct(plabel) into: conlab separated by ', ' from conlab;
		  quit;
	    
        Data _mg;
	      set _contrast;run; 
       %end;

	   %else %do;
	      %let upeffect=%upcase(&effect); /* delete &effect &sliceby and &effect*&sliceby variables in above data set*/
		  %let upsliceby=%upcase(&sliceby);
	      DATA _mg;
          set _mg _contrast; 
		  if upcase(parameter) not in ("&upeffect"  "&upsliceby" "&upeffect*&upsliceby" );
          RUN;
	   %end;

	%end;


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
   OPTIONS ORIENTATION=&ORIENTATION MISSING = "-" NODATE;
   ODS rtf STYLE=TABLES file= "&OUTPATH.&FNAME &SYSDATE..DOC"; 

   PROC REPORT DATA=_mg HEADLINE HEADSKIP CENTER SPANROWS LS=256
      STYLE(REPORT)={JUST=CENTER} SPLIT='~'   nowd; 
      COLUMN plabel  Level 
	      %if &clnum=T %then %do; n %end;
         /* If glogit model add column for response level */
          %if &model=generalized logit %then %do; Response %end;
          ("&headLab" '------------------------------------------'
          (OR_CI ProbChiSq  %if &type3 = T %then %do; p_type3 %end;)); 
      DEFINE plabel/order order=data "Covariate"  STYLE(COLUMN) = {JUST = L CellWidth=30%}; 
      DEFINE Level/"Level" STYLE(COLUMN) = {JUST = L CellWidth=15%}; 
        %if &clnum=T %then %do;
		DEFINE n /display 'N' STYLE(COLUMN)={JUST=C CellWidth=10%};
		%end;

       %if &type3 = T %then %do; 
         DEFINE p_type3/order MISSING "Type3 P-value" STYLE(COLUMN) = {JUST=C CellWidth=8%} 
               FORMAT=PVALUE8.3; 
       %end;

      %if &model=generalized logit %then %do;
         DEFINE response/order=data "&outlab" STYLE(COLUMN) = {JUST = L CellWidth=10%};
      %end;
      DEFINE OR_CI/DISPLAY "Odds Ratio~(95% CI)" STYLE(COLUMN) = {JUST = C CellWidth=20%}; 
      DEFINE ProbChiSq/DISPLAY "OR P-value" STYLE(COLUMN)= {JUST = C CellWidth=8%} 
               FORMAT=PVALUE8.3; 
       
      /*COMPUTE ProbChiSq; 
	   IF . < ProbChiSq <0.05 THEN CALL DEFINE("ProbChiSq", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
           if plabel = "Comparisons Stratified by &Slicebylabel (INTERACTION):" THEN CALL DEFINE("ProbChiSq", "STYLE", "STYLE=[COLOR=WHITE]");
      ENDCOMP; */
       
       %if &effect ~= %str() %then %do;
       COMPUTE plabel;
          if plabel = "Comparisons Stratified by &Slicebylabel :" THEN CALL DEFINE("plabel", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]");  
       ENDCOMP;

       COMPUTE level;
          if level = "&effectlabel :" THEN CALL DEFINE("level", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]");  
       ENDCOMP;
       
	   %end;
      %if &type3 = T %then %do; 
            COMPUTE p_type3; 
            IF .< p_type3 <0.05 THEN CALL DEFINE("p_type3", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         ENDCOMP;   
      %end;

      compute after plabel; 
         line '';   
      endcomp;      

      compute after _page_;
         length text1 - text2 $120.;
         text1 = "*  Number of observations in the original data set = %TRIM(&nobs_read). Number of observations used = %TRIM(&nobs_used).";
         line @0 text1 $120.;

         %if &footnote ~= %STR() %then %do;
            line @0 &footnote;
         %end;
		 %if &effect ~=%str() and &shortreport = T %then %do;
			line @0 "*** The estimated stratified treatement effect was controlled by: &conlab";
		%end;
      ENDCOMP; 
   RUN;

   ODS RTF CLOSE;

   /* Reload original options that were in use before running the macro */
   PROC OPTLOAD data=_options;
   RUN;

   %if &debug = F %then %do;
      /* Delete temporary datasets */
   ods select none;
      PROC DATASETS lib=work;
         DELETE _est2 _ci2 _mg _cont _options %if &classexist = 1 %then %do; typ3_ %end;
		 %if &effect = T %then %do; _contrast _temp %end;;
      QUIT;
	  ods select all;
   %end;

%mend MULTIPLE_LOGREG;


