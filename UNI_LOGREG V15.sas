/**********************************************************************************************
Macro Name: UNI_LOGREG
Created Date/Author: Aug. 23, 2012/Sungjin Kim
Last Update Date/Person: Jan, 2017/Yuan Liu
List of contributors: Dana Neckleach

Current Version: V15
Working Environment: SAS 9.4 English version

Contact: Dr. Yuan Liu YLIU31@emory.edu

PURPOSE: To conduct univariate logistic regression for each variable in the dataset.  A 
proportional odds model is fitted for a binary outcome variable; a cumulative logit function is
used for an ordinal (ordered) outcome variable with 3 or more levels.  The odds ratio with 95%
CI is presented along with the p-value.  For categorical variables, the reference group will be
shown along with the number of observation in each category.  A list of variables with 
type3 p-value < PSELECT will be written to the log. In addition two global macro variables,
_select_cvar (list of categorical variables) and _select_nvar (list of numerical variables) will be created, 
which could may be referred in the future multivariable modeling.   

NOTE:    This macro does not work with a generalized logit model with a nominal outcome. 

Parameters: 

DATASET        The name of the data set to be analyzed.
OUTCOME        The name of the outcome variable.  It must be binary or ordinal.  If ordinal
               then it should be a numeric variable.
EVENT          The event category for the binary response model (optional).  You can specify 
               the value in quotes.  This will be passed to the event= option in the model 
               statement.  Leave this blank if you have an ordinal outcome with more than 2 
               levels.
DESC           Set to T to reverse the order of an ordinal outcome (optional).  The order will be
               based on the internal order.  Only specify this if the EVENT parameter is blank.  
               The default value is F.
CLIST          List of categorical variables, separated by empty space, or by * if need to change
			   the reference level. See example.  
NLIST          List of numerical variables, separated by empty space.
STRATA		   The STRATA statement names the variables that define strata or matched sets to 
			   use in stratified logistic regression of binary response data. See SAS help Manual.
TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.  This will only have an effect if 
               categorical variables are specified in CLIST.
FIRTH          Set to T to use Firth’s penalized maximum likelihood estimation to reduce bias
               in the parameter estimates.  The default value is F.
PO             Set to T to check the proportional odds assumption by reporting the score test
               p-values.  Only set to T if using an ordinal outcome.  The default value is F.
PSELECT		   generate varaible list that with type 3 p-vale < PSELECT. Two macro variable _select_var 
				_select_cvar can be used in future multivariable model. 
DOC            Set to T to create a RTF file containing the output or F to suppress creation of 
               the RTF file (optional).  The default value is T.
OUTPATH        Path for output table to be stored.
FNAME          File name for output table.
ORIENTATION    Value of PORTRAIT or LANDSCAPE to indicate the paper layout of the report 
               (optional).  The default value is PORTRAIT.  
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
For more details, please see the related documentation
**********************************************************************************************/

%MACRO UNI_LOGREG(DATASET=, OUTCOME=, EVENT=, DESC=F, CLIST=, NLIST=, STRATA=,type3=T, 
     FIRTH=F, PO=F, PSELECT=0.2,DOC=T, OUTPATH=, FNAME=, ORIENTATION=PORTRAIT, DEBUG=F); 

   %local OUTPATH TYPE3P I CVAR OUTCOME progress DOC OR DATASET __MACRO_ERR WORK_SETS OR_LB
     EVENT REFVAR N NVNAME N_OBS CLIST OUTCOMEVV TYPE3 FNAME DEBUG CLIST_CNT CREFLIST NLIST C
     NVAR OR_UB ORP ORIENTATION out_cat_cnt MODEL MAX MIN RESP out_type;
   %local cnt cnt_ pos pos_ L1 L2 clist_ _cnt_ clist_cnt i cnt_bef  cvar_ref cnt2 len2 cntw pos2 pos3 bef cnt_bef cvar_;

   %global _select_nvar  _select_cvar;


   /* Prevent case sensitivity */
   %let DEBUG = %UPCASE(&DEBUG);
   %let type3 = %UPCASE(&type3);
   %let desc = %UPCASE(&desc);
   %let firth = %UPCASE(&firth);
   %let po = %UPCASE(&po);
   %let doc = %UPCASE(&doc);
   %let strata = %UPCASE(&strata);

   /* Initialize error flag */
   %let __Macro_Err = 0;

   /* Get list of data sets in work library to avoid deletion later */
      ods select none;
   ODS OUTPUT Members(nowarn)=_DataSetList;
   PROC DATASETS lib=work memtype=data;
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
   PROC OPTSAVE out=_options;RUN;


   /*Set allowance of length of variable name, this will be updated based on the actual data*/
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


   /* Get number of outcome levels */
   /* Have to use proc freq not proc sql because sql does not collapse formatted values */
  ods select none;
   PROC FREQ DATA = &dataset ;
      tables &outcome/out=_freq;
   RUN;
   ods select all;

   PROC SQL noprint;
      select count(*) into :out_cat_cnt
      from _freq
      /* Don't count missing values */
      where MISSING(&outcome) = 0;
   QUIT;

   /* Find out if outcome variable is character or numeric */
   DATA _NULL_;
      set &dataset;
      CALL SYMPUT('out_type',VTYPE(&outcome));
   RUN;

   /* Check for outcome variable names that are too long */
  /* %DO i = 1 %to &clist_cnt; 
      %IF %LENGTH(%SCAN(&CLIST, &i)) > 30 %then %do;
         %put ERROR: Variable name %SCAN(&CLIST, &i) must be less than 31 characters.; 
           %let __Macro_Err=1;
       %END;
   %END;*/

   /* Make sure PO was not set to T if inappropriate */
   %if &out_cat_cnt <= 2 and &po = T %then %do;
      %put ERROR: PO=T, but the outcome variable only has 2 levels.  The proportional odds test cannot be conducted for binary logit models; 
      %let __Macro_Err=1;
   %end;

   /* Make sure outcome is numeric for cumulative logit models */
   %if &out_cat_cnt > 2 and &out_type = C %then %do;
      %put ERROR: The outcome variable must be numeric for cumulative logit models.; 
      %let __Macro_Err=1;
   %end;

   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   %put &cvarlist;
   %put &clist;
   *CHARACTER Variables;
   %IF &clist_cnt > 0  %THEN %DO; 
      %do C = 1 %to &clist_cnt;

	    %if &pos > 0 %then %let REFVAR = %SCAN(%superq(CLIST), &C, '*'); 
		%else  %LET REFVAR = %SCAN(%superq(CLIST), &C, ' ');
		%let CVAR = %SCAN(%superq(cvarlist), &C, ' ');

         ODS select none;
         ODS OUTPUT  "Parameter Estimates" = _mlec "Type 3 Tests" = _type3 
          "95% Clodds=Wald" = _or  ResponseProfile=_resp ModelInfo=_info
           %if &po = T %then %do; CumulativeModelTest=_score %end;;
         proc logistic data=&dataset namelen=&length;
            /* If cumulative logit model then remove formats from outcome otherwise problems */
            /* will be caused when determining the response order */
            %if &out_cat_cnt > 2 %then %do; FORMAT &OUTCOME; %end;
            class &REFVAR / param=ref order=internal;
            /* Note that if EVENT and DESC are specified that EVENT will override DESC */
            /* EVENT will be ignored if link=clogit */
            model &OUTCOME (%if %sysevalf(%superq(event)~=,boolean) %then %do; 
               event=&event %end; order=internal %if &desc = T %then %do; desc %end;) = &CVAR/
               clodds=wald %if &firth = T %then %do; firth %end;;
			%if &strata ~=%str() %then %do; strata &strata %end;;
         run;
		 ods select all;
          
         ODS select none;
         ODS OUTPUT 'One-Way Frequencies' = FREQ; 
         PROC FREQ DATA=&dataset ORDER=data; 
            TABLE &CVAR;
            WHERE &OUTCOME is not missing;
         RUN;  
		 ods select all;

           
         DATA _FREQC; 
            length Covariate Name $256. Level $96.;
            SET FREQ ;
            Covariate = left(label(&CVAR)); 
			Name = "&CVAR";
            Level = left(F_&CVAR); 
            /*Variable = substr(Table,7,96);*/
            KEEP Covariate Name Level Frequency Variable; 
         RUN;

         data _type3_; 
            length Name $256;
            set _type3;
			Name = "&CVAR";
            rename ProbChiSq = type3p ; 
            format ProbChiSq;
            keep Name ProbChiSq;
         run;

         proc sort data=_type3_; by Name; run; 
         proc sort data=_FREQC; by Name; run; 

         data _mg; 
            merge _type3_ _FREQC; 
            by Name; 
            keep Covariate Name Level Frequency type3p; 
         run; 

         data _mlec_; 
		    length Name $256;
            set _mlec; 
			Name = "&CVAR";
            format ProbChiSq; 
            where Variable ne 'Intercept'; 
            morder = _n_; 
         run;

         data _comb; 
            merge _mlec_ (rename=(ClassVal0=Level)) _or (rename=(OddsRatioEst=OddsRatio)); 
            keep Name Level ProbChiSq OddsRatio LowerCL UpperCL morder; 
         run; 
 
         proc sort data=_comb;by Name Level;run;
         proc sort data=_mg; by Name Level; run;

         data _mergedata; 
            length Covariate $256. Level $96. ;
            merge _mg _comb; 
            by Name Level; 
            if morder = . then morder = 999;
            /* Preserve variable order */
            order2 = &c;
            keep order2 Covariate Name Level Frequency OddsRatio LowerCL UpperCL ProbChiSq morder type3p; 
         run; 
                    
         /* Merge on score test */
         %if &PO = T %then %do;
            data _mergedata; 
               set _mergedata;
               if _n_ = 1 then set _score (KEEP=ProbChiSq RENAME=(ProbChiSq=PO_Pval));
            RUN;
         %end;

         proc sort data=_mergedata out=_charr&C (DROP=morder); 
            by Covariate morder; 
         run;  
                                   
      %END; 

      data _charr_all; 
         set _charr1 - _charr%eval(&clist_cnt); 
         order = 1;
      run;
   %END;

   *NUMERIC Variables;
   %IF &nlist_cnt >0   %THEN %DO; 
          
      %do N = 1 %to &nlist_cnt;

         %LET NVAR = %SCAN(&NLIST, &N); 
         ods select none;
         ODS OUTPUT  "Parameter Estimates" = _mlec "Global Tests" = _wald "95% Clodds=Wald" = _or
          "Observations summary" = _nobs  ResponseProfile=_resp ModelInfo=_info
          %if &po = T %then %do; CumulativeModelTest=_score %end;;
         proc logistic data=&dataset namelen=&length ;
            /* If cumulative logit model then remove formats from outcome otherwise problems */
            /* will be caused when determining the response order */
            %if &out_cat_cnt > 2 %then %do; FORMAT &OUTCOME; %end;
            /* Note that if EVENT and DESC are specified that EVENT will override DESC */
            /* EVENT will be ignored if link=clogit */
            model &OUTCOME (%if %sysevalf(%superq(event)~=,boolean) %then %do; event=&event %end; 
               %if &desc = T %then %do; desc %end;
               order=internal)= &NVAR/clodds=wald %if &firth = T %then %do; firth %end;;
         run; 
		 ods select all;

         data _mlec;
		    length Name $256;
            set _mlec; 
			Name = "&CVAR";
            format ProbChiSq; 
            where Variable ne 'Intercept'; 
         run;

         data _wald;
            set _wald; 
            format ProbChiSq; 
         run;

         data _NULL_; 
            set &dataset; 
            call symput('NVNAME', put(label(&NVAR),$256.));
         run;

         proc sql noprint;
            select ProbChiSq into: type3p from _wald where test ="Wald";
            select OddsRatioEst into: or from _or;
            select LowerCL into: or_lb from _or;
            select UpperCL into: or_ub from _or;
            select ProbChiSq into: orp from _mlec;
            select SumFreqsUsed into: N_obs from _nobs;
         quit;

         data _numm&N; 
            length Covariate Name $256.;
            Covariate = "&NVNAME";
			Name = "&NVAR";
            Level = ' ';
            Frequency = &N_obs;
            type3p = &type3p;
            LowerCL = &or_lb;
            UpperCL = &or_ub;
            OddsRatio = &or;
            ProbChiSq = &orp;
            /* Preserve order */
            order2 = &n;
            keep order2 Covariate Name Frequency Level type3p OddsRatio LowerCL UpperCL ProbChiSq; 
         run;
                    
         /* Merge on score test */
         %if &PO = T %then %do;
            data _numm&N; 
               set _numm&N;
               if _n_ = 1 then set _score (KEEP=ProbChiSq RENAME=(ProbChiSq=PO_Pval));
            RUN;
         %end;

      %END; 

      data _numm_all; 
         set _numm1 - _numm%eval(&nlist_cnt);
         order = 2; 
      run; 
   %END;

   /* Combine categorical and numeric covariate results */
   DATA _report;
      set %IF &clist_cnt > 0  %then %do; _charr_all %end; 
         %IF &nlist_cnt > 0 %then %do; _numm_all %end;;

      /* OR and 95% CI */
      if OddsRatio ~= . then 
         OR_CI = TRIM(LEFT(PUT(OddsRatio,8.2))) || " (" || 
         TRIM(LEFT(PUT(lowerCL,8.2))) || "-" || 
         TRIM(LEFT(PUT(upperCL,8.2))) || ")";
      else OR_CI = '-';
   RUN;

   /* Get model type - binary logit or cumulative logit */
   PROC SQL noprint;
      select Value into :model
      from _info
      where Description = 'Model';
   QUIT;

   /* For cumulative logit models find out the direction of the outcome */
   %if &model = %QUOTE(cumulative logit) %then %do;
      PROC SQL noprint;
         select Outcome into :max
         from _resp
         where OrderedValue = 1;
         select Outcome into :min
         from _resp
         having OrderedValue = max(orderedvalue);
      QUIT;

      %if &max > &min %then %let resp = higher;
      %else %if &min > &max %then %let resp = lower;
   %end;

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

   /* Remove outer quotes for report */
   %if %sysevalf(%superq(event)~=,boolean) %then %do;
      %let event = %sysfunc(DEQUOTE(&event));
   %end;

   **** PRINT THE TABLE ************;
   data _NULL_; 
      set &dataset; 
      call symputx('outcomevv', label(&OUTCOME)); 
   run;

   OPTIONS ORIENTATION=&ORIENTATION MISSING = "-" NODATE;

   %if &doc = T %then %do;
      ODS RTF STYLE=TABLES file= "&OUTPATH.&FNAME &SYSDATE..DOC";  
   %end;

   PROC REPORT DATA=_report HEADLINE HEADSKIP CENTER STYLE(REPORT)={JUST=CENTER} SPLIT='~'
     nowd SPANROWS LS=256; 
      COLUMNS order order2 Covariate 
          %IF &clist_cnt > 0 %then %do; Level %end; 
          Frequency (
          %if %sysevalf(%superq(event)~=,boolean) %then %do; "&outcomevv=&event" %end;
          %else %do; "&outcomevv" %end;
               '------------------------------------------'(OR_CI ProbChiSq  
          %if &type3 = T %then %do; type3p %end;
          %if &PO = T %then %do; PO_pval %end;)); 
      DEFINE order/order order=internal noprint;
      DEFINE order2/order order=internal noprint;
      DEFINE Covariate/ order order=Data  "Covariate"  STYLE(COLUMN) = {JUST = L}; 
      %IF &clist_cnt > 0 %then %do; 
         DEFINE Level/ DISPLAY   "Level"     STYLE(COLUMN) = {JUST =L}; 
      %end;
      DEFINE Frequency/DISPLAY "N" STYLE(COLUMN) = {JUST = C};                
      DEFINE OR_CI/DISPLAY "Odds Ratio~(95% CI)" STYLE(COLUMN) = {JUST = C CellWidth=15%}; 
      DEFINE ProbChiSq/DISPLAY  "OR P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=PVALUE5.3; 

      COMPUTE ProbChiSq; 
         IF  . < ProbChiSq <0.05 THEN CALL DEFINE("ProbChiSq", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
      ENDCOMP; 

      %if &type3 = T %then %do;
         DEFINE type3p/order "Type3 P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=PVALUE5.3; 
         COMPUTE type3p; 
            IF  . < type3p <0.05 THEN CALL DEFINE("type3p", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         ENDCOMP; 
      %end;

      %if &PO = T %then %do;
         DEFINE PO_pval/order "Assumption P-value" STYLE(COLUMN) = {JUST = C CellWidth=9%} FORMAT=PVALUE5.3; 
         COMPUTE PO_pval; 
            IF  . < PO_pval <0.05 THEN CALL DEFINE("PO_pval", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         ENDCOMP; 
      %end;

      COMPUTE after Covariate; line ''; ENDCOMP; 

      /* For cumulative logit models add a footnote */
      %if &model = %QUOTE(cumulative logit) %then %do;
         compute after _page_;            
           line @0 "The probability of having &resp values of the outcome is being modeled.";
         ENDCOMP;
      %end;
   RUN;

   %if &doc = T %then %do;
      ODS RTF CLOSE;
   %end;

	/** extract variable name with type 3 p-value < PSELECT in univariate analysis;*/
   proc sort data=_report out = _report_; by name; run;
   data _report_; set _report_; by name; if first.name;run;

   		%let _select_cvar = ;
   		%let _select_nvar = ;		
        
		proc sql noprint;
		select  name into :_select_cvar separated by " "
		from _report_
		where Level ~= "" and type3p < &PSELECT
		order by order2;

		select  name into :_select_nvar separated by " "
		from _report_
		where Level = "" and type3p < &PSELECT
		order by order2;
		quit;

   %put Categorical variables selected at pvalue < &PSELECT: &_select_cvar;
   %put Numerical variables selected at pvalue < &PSELECT: &_select_nvar;


   /* Reload original options that were in use before running the macro */
   /* This is mainly to reset the orientation */
   PROC OPTLOAD data=_options;  RUN;

   /* Only delete files if not in debug mode */
   %if &DEBUG ~= T %then %do;

      /* If there are work data sets that should not be deleted */
      %if %sysevalf(%superq(work_sets)~=,boolean) %then %do;
         /* DELETE ALL TEMPORARY DATASETS that were created */
	     ods select none;
         proc datasets lib=work memtype=data nolist;  
            save &work_sets;
          quit;  
		  ods select all;
       %end;
       %else %do;
         proc datasets lib=work kill memtype=data nolist;  
         quit; 
      %end;
   %end;

%exit:
%MEND UNI_LOGREG; 
