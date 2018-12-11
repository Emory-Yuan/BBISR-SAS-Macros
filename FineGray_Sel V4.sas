/**********************************************************************************************
Macro: %FineGray_sel
Created Date/Author/Contact: Sep 21, 2016/Chao Zhang and Yaqi Jia;
Last Update Date/Person/Contact: Oct 2016/Yuan Liu
Current Version: V3
Working Environment: SAS 9.4 English version

Contact: Dr. Yuan Liu/email: yliu31@emory.edu

Purpose:  To conduct multivariable analysis for competing risk model. The proportional 
subdistribution hazards model as proposed by Fine and Gray (1999) was used followed by 
backward elimination. 

Notes: The model runs using PROC PHREG.  The final list of variables selected is saved in two
macro variables: &_finalcvar and &_finalvar, and is also written into the log.The macro
“MUTLIPLE_PHREG V21.sas” or later is also required.

Parameters: 

DSN            The name of the data set to be analyzed.

EVENT          Name of time to event outcome variable.   

CENSOR         Name of censoring indicator variable.  Values of 0 indicate censored. 

EVENT_CODE     The value in CENSOR that indicate event of interest, and this value will 
			   appear EVENTCODE= option.

VAR            The list of variables on interest in the initial model that would be 
			   eliminated during the backward selection procedure separated by spaces. 
			   The order of variables in this list will be preserved in the final report.

CVAR           The list of categorical variables that are in VAR and FORCEINVAR. If need to 
			   change the reference group, you can follow each variable name by (DESC) or 
			   by (ref = “Ref level in formatted value”) where needed and separate terms by *.
			   See code example. 

INC            Number of variables to include in the model (optional).  The first n variables 
               in the VAR parameter will be included in every model.  The default value is 0.   

ALPHA          the criterion for retaining a variable in the backward elimination procedure.

TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.  This only has an effect if REPORT = T.

ID             Variable to be used in the ID statement in PHREG (optional).  Refer to SAS
			   Help and Documentation for proper use of this option.


CLNUM          Set to T if you want to see the number of observations for each level of covariates. The default is T.

ORIENTATION	   orientation of the output Word table. Default is portrait, can be changed to landscape.

REPORT 		   Set it to T if a results summary table is desired. Otherwise check Log for variable 
			   selected by the backward elimination.

FILENAME       File name for output table.  This is necessary if report=T.

OUTPATH        File path for output table to be stored.  This is necessary if report=T.

DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
* For more details, please see the related documentation
**********************************************************************************************/

%macro FINEGRAY_SEL(dsn=,  event=,  censor=,  event_code=1, var=,  cVar=, inc= 0, alpha=.2,
         type3=T, id=,  clnum=T, ORIENTATION = portrait, report =T, outpath=, filename=, debug=F) ;
     
   %local cnt cnt_ pos pos_ L1 L2 clist_ _cnt_ cvar_cnt i cnt_bef  cvar_ref cnt2 len2 cntw pos2 pos3 bef cnt_bef cvar_ nn;

   /* Capitalize */
    
    %let debug      = %UPCASE(&debug);
	%let type3      = %UPCASE(&type3);
	%let report     = %UPCASE(&report);
    
   /* Macros for final variable lists */
   %global _finalvar _finalcvar;

  /* Save current options */
   PROC OPTSAVE out=_options;   RUN;

  /* Count number of variables in var */
    %let var_cnt = %sysfunc(countw(&var,' '));



  /* It could be longer if there are interaction terms */
   %let length = 64; 


   /* Time Should not be missing value*/
   
   %if &event = %STR()  %then %do;
      %put ERROR: EVENT should be specified.; 
      %goto exit;
   %end;



/*Check CVAR whether in the proper format, such as separate by * if reference level is specified*/
  /* Count number of class variables */

  %if %superq(cvar) = %str( ) %then %let cvar_cnt = 0;
  %if %superq(cvar) ~= %str( ) %then %do;
  		
        %let cnt = %qsysfunc(countw(%superq(cvar),'*'));
		%let cnt_= %qsysfunc(countw(%superq(cvar),' '));

        %let pos=%qsysfunc(countc(%superq(cvar), '(' ));
        %let pos_ = %qsysfunc(countc(%superq(cvar), '*' ));

		%if %superq(pos) = 1 %then %do;
		%let L1=%qsysfunc(findc(%superq(cvar),  '(' ));
	    %let L2=%qsysfunc(findc(%superq(cvar),  ')' ));
		%let clist_=%qsysfunc(SUBSTR(%superq(cvar), %superq(L1),%eval(%superq(L2) - %superq(L1) +1)));
	    %let _cnt_ = %qsysfunc(countw(%superq(clist_),' '));
        %end;

        %if %superq(cnt) = 1 and  %superq(cnt_) = 1 %then %let cvar_cnt = 1;
		%else %if %superq(pos_) > 0 %then %let cvar_cnt = %superq(cnt);
		%else %if %superq(pos_) = 0 and %superq(pos) = 0 %then %let cvar_cnt = %superq(cnt_);
		%else %if (%superq(pos) = 1 and %superq(pos_) = 0) and %superq(cnt_) = %superq(_cnt_) %then  %let cvar_cnt = 1;
		%else %do; %put ERROR: The categorical variables in CVAR should be separated by * if you specify the reference level.; 
          		   %goto exit; %end;

  %end;

/*Build up categorical variable list without reference level specified*/
  %let cvarlist =; 

  %if &cvar_cnt = 1 and &pos = 1 %then %let cvarlist = %sysfunc(SUBSTR(%superq(cvar), 1,%eval(&L1-1)));
  %else %if &cvar_cnt = 1 and &pos = 0 %then %let cvarlist = %superq(cvar);
  %if &cvar_cnt > 1  and %superq(pos_) = 0 %then  %let cvarlist = %superq(cvar) ; 

  %if &cvar_cnt > 1  and %superq(pos_) > 0 %then %do; 
 			
      %do i = 1 %to &cvar_cnt; 

				%let cvar_ref = %SCAN(&cvar, &i, '*'); 

				%let cnt2=%qsysfunc(countc(&cvar_ref, '(' ));
				
				%let len2=%LENGTH(&cvar_ref);
				%let cntw=%qsysfunc(countw(&cvar_ref,' '));

				/*%put cvar_ref = &cvar_ref cnt2= &cnt2  len2=&len2 cntw =&cntw ;*/

				%if %superq(cnt2) = 1 %then %do;

	                %let pos2=%sysfunc(findc(&cvar_ref, '(' ));
					%let pos3=%sysfunc(findc(&cvar_ref, ')' ));

					%let bef = %sysfunc(SUBSTR(&cvar_ref, 1,%eval(&pos2-1)));
			        %let cnt_bef = %qsysfunc(countw(&bef,' '));

					%if %superq(cnt_bef) = 1 and  %superq(pos3) = %superq(len2) %then %let cvar_=&bef;
					%else %do; %put ERROR: The categorical variables in CVAR should be separated by * if you specify the reference level.; 
          		                %goto exit;%end;
					%end;
					
                %if %superq(cnt2) = 0 and %superq(cntw) = 1 %then %LET cvar_ = &cvar_ref;

				%if %superq(cnt2) > 1 or ( %superq(cntw) > 1 and %superq(cnt2) = 0) %then %do; %put ERROR: The categorical variables in CVAR should be separated by * if you specify the reference level.; 
          		%goto exit;  %end;	 

			%let cvarlist = %superq(cvarlist) %superq(cvar_); 
 %end;
		 
 %end;

 /* Make sure that CVAR are either in VAR or FORCEINVAR*/
 %put &var   ;
 %put &cvarlist;

   %do i = 1 %to &cvar_cnt;
 		%let _check_ = 1;
     %do j = 1 %to &var_cnt;
         %if %upcase(%SCAN(&var, &j)) = %upcase(%SCAN(&cvarlist, &i)) %then %let _check_ = 0;;
     %end;
   %end;

   %if &_check_ = 1 %then %do; 
   %put ERROR: The varaible appears in CVAR should also be in VAR.; 
   %goto exit;  %end;


  /* Check for variable names that are longer than the default and increase length as needed */
   %DO i = 1 %to &cvar_cnt; 
      %IF %LENGTH(%SCAN(&var, &i,' ')) > &length %then %let length = %LENGTH(%SCAN(&var, &i,' ');
   %END;


   /*Fine&Gray model with backwards elimination */ 

   %let remlist =;

  %do %while (&var_cnt > 0);

	%if &cvar_cnt > 0 %then %let  cvar_ref = %sysfunc(tranwrd(%superq(cvar),*,%str()));
    
   ods select none;
   ODS OUTPUT  NObs = numobs ClassLevelFreq =clfreq ParameterEstimates = MLEC ModelInfo=modelinf Type3 = type3  ;
   proc phreg data=&dsn simple namelen=&length %if &id ~= %STR() %then %do; COVS(AGGREGATE) %end;;
        %if &cvar_cnt > 0 %then %do; 
        class &cvar_ref /order=internal param=glm;%end;
        model &event*&censor(0) = &var / eventcode=&event_code rl;
		%if &id ~= %STR() %then %do;
            id &id;
         %end;
    run;
	ods select all;

	/*Identifiy next candidate to be removed*/

     proc sql noprint; select count(*) into :nn from type3; quit;

	 %if &nn = &inc %then %goto exit2;
	 %else %do;

	* in Type3 data remove variables that should not be fixed in model;
    	Data typ3_n; set type3; if _N_ > &inc;run;
	* for the rest variables, see which one has the largest p value and compare it to alpha; 

	 proc sql noprint; 
	    select Effect into: var_ separated by ' ' from typ3_n ;quit; * keep the original order;

     proc sort data=typ3_n; by desending ProbChiSq;run;
     proc sql noprint;
	   select ProbChiSq into: rplist separated by ' ' from typ3_n ;
       select Effect into: rvlist separated by ' ' from typ3_n ;
     quit;


	%let premove_p = %scan(&rplist,1, ' '); 
	%let premove = %scan(&rvlist,1, ' ');

	%put &premove_p &premove;

	* if the current largest pvalue is less than alpha, then stop, else update &var and &cvar;
	%if &premove_p < &alpha %then %goto exit2;

	%else %do;
	     %let remlist = &remlist &premove;
		 /* update &var list and keep the same input order;*/
	     %let var = %sysfunc(tranwrd(%upcase(&var),%upcase(&premove),%str())); 

		 /* update &cvar;*/
		 %if &cvar_cnt >  0 %then %do;* update &cvar;

		    %let pos_=%qsysfunc(findc(%superq(cvar), '*' ));
			

		    %let newcvar = ;
		 	%do i = 1 %to &cvar_cnt;
/*				%put pos_ = &pos_ cvar= &cvar;*/
  				%if &pos_ > 0 or &cvar_cnt = 1 %then %let scanvar = %SCAN(&cvar, &i, '*');
				%else %let scanvar = %SCAN(&cvar, &i, ' ');

				%let pos = %qsysfunc(findc(%superq(scanvar), '(' ));

				%if &pos > 0 %then %let trimvar=%qsysfunc(SUBSTR(%superq(scanvar), 1,%eval(&pos-1)));
				%else %let trimvar = &scanvar;

/*				%put trimvar = &trimvar scanvar=&scanvar premove=&premove;*/

				%if %upcase(&trimvar) = %upcase(&premove) %then %let newcvar = &newcvar; 
				%else  %let newcvar = &newcvar*&scanvar; 

			%end;
			%let cvar = &newcvar;
		 %end;
	%end;

    %let var_cnt = %sysfunc(countw(&var,' ')); * update number of variable in &var; 
    %let cvar_cnt = %sysfunc(countw(&cvar,'*'));  * update number of variable in &cvar; 

%end;/*%do %if (&var_cnt > 0);*/

   /* Save in macro variables */
   %let _finalvar = &var;
   %let _finalcvar = &cvar;

%end; /*end of 	 %if &nn = &inc %then %goto exit2;  %else %do;*/

%exit2:

/* Produce report of final model */
   %if &report = T %then %do;
      %let forceInVar = .;
	  %if &inc >0 %then %do;

	    Data typ3_n; set type3; if _N_ <= &inc;run;
		proc sql noprint;
		select Effect into :inforcevar separated by ' '
		from Typ3_n;

        PROC CONTENTS DATA = &dsn (keep=&inforcevar) out=cont noprint; RUN;

		 DATA cont;
            set cont end=last;
             if label = ' ' then label = name;
             if _N_ ~= 1 and last then label = 'and ' || label;
         RUN;

         PROC SQL noprint;
            select label into :inLab separated by ', '
             from cont;
         QUIT;
		 %end;
		 %else %let inLab = None;

	     
      /* If any variables were removed from the model create a list */
      %if &remlist ~= %STR() %then %do;
         /* Get variable labels for footnote */
         PROC CONTENTS DATA = &dsn (keep=&remlist) out=cont noprint;
         RUN;

         DATA cont;
            set cont end=last;
             if label = ' ' then label = name;
             if _N_ ~= 1 and last then label = 'and ' || label;
         RUN;

         PROC SQL noprint;
            select label into :remLab separated by ', '
             from cont;
         QUIT;

           %let foottext = The following variables were forced in the model: &inLab . The following variables were removed from the model: &remLab..;

      %end;
      %else %let foottext = The following variables were forced in the model: &inLab . No variables were removed from the model.;

      %MULTIPLE_PHREG(        
          OUTPATH = &outpath, 
          FNAME = &filename,
		  ORIENTATION = &ORIENTATION,
		  clnum=&clnum,
          TYPE3 = &type3,
          debug= &debug,
          FOOTNOTE="** Backward selection with an alpha level of removal of &alpha was used.  &foottext"
        );

  %end;


   %if &debug = F %then %do;
      *--- DELETE ALL TEMPORARY DATASETS that were created; 
      proc datasets lib=work memtype=data noprint;  
           delete %if cvar ~= %STR() %then %do; freq  %end; mergedata mlec modelinf numobs type3 typ3_n keptlist Cont
                  result1 result2 result2_ result3 Ori_list Removed _options _cont ;
       quit;
    %end;
   
/* Final variables selected (get rid of double spaces */
   %put Categorical variables selected: &_finalcvar;
   %put All variables selected: &_finalvar;

%exit:
%mend finegray_sel;
