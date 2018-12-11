/**********************************************************************************************
Macro: LOGREG_SEL
Created Date/Author Oct. 19, 2012/Dana Nickleach
Last Update Date/Person Jan, 2017/Yuan Liu
Other Signficiant Contributor: Mar 2016/Yaqi Jia
Current Version: V15
Working Environment: SAS 9.4 English version

Contact: Dr. Yuan Liu yliu31@emory.edu

Purpose:  To conduct backward selection on a logistic regression model using the maximum 
possible sample size at each stage of the selection process instead of restricting to the 
sample size from the first step as SAS does when using their selection methods.  Optionally, a 
table of the resulting model can be generated.  

Notes: The model is run using PROC LOGISTIC.  A binary outcome or ordinal outcome using a
cumulative logit model can be used, but not a nominal outcome.  The final list of variables 
selected will be written to the log.  Additionally, two global macro variables, _finalvar and 
_finalcvar will be created containing the list of all variables and categorical variables 
selected, respectively.  If you are requesting a table with the model results then the macro
“MUTLIPLE_LOGREG V15” or later is also required.  Interactions can be included to obtain the
estimate of treatment effect(TRT) in each level of stratified variable (SV), and it is required 
both TRT and SV to be categorical variables. For variables selection, put TRT, SV, and TRT*SV in 
the begining of VAR; use INC = 3 to force the two main effects and their interaction in the model;
use EFFECT = TRT and SLICEBY = SV to generat the stratified treatement odds ratio.



Parameters: 

DSN            The name of the data set to be analyzed.
OUTCOME        The name of the outcome variable.  It must be binary or ordinal. 
EVENT          The event category for the binary response model.  Specify the value in quotes.  
               This is the argument that will be passed to the event= option in the model 
               statement.  Leave this blank if you have an ordinal outcome with more than 2 
               levels.
DESC           Set to T to reverse the order of an ordinal outcome (optional).  The order will 
               be based on the internal order.  Only specify this if the EVENT parameter is 
               blank.  The default value is F.
VAR            List of variables to include in the model separated by spaces.                                           
CVAR           List of categorical variables to include in the model separated by spaces.  These 
               should also appear in the var parameter.  If you want to change the reference 
               group you can follow each variable name by (desc) where needed.  However, you 
               will need to separate terms with an asterisk instead of a space. 
INC            Number of variables to include in the model (optional).  The first n variables in 
               the var parameter will be included in every model.  The default value is 0.  
SLSTAY         The significance level for removing variables from the model (optional).  The 
               default value is .05.   
WEIGHT         Variable to use in the weight statement (optional).  Weights will be normalized 
               to the original sample size using the normalize option.  Leave it blank if not 
               using weights.
REPORT         Set this to T if you want a table of the resulting model generated (optional).  
               The default value is F.  
TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.  This only has an effect if REPORT = T. Set it 
				to T when there is an interaction in the model.
EFFECT         Use to specify the treatment variable in the interaction.Use in combine with SLICEBY
			   and if not empty, VAR should contain a two-way interaction. See example.
SLICEBY		   Use in combine with EFFECT to specify the stratified variable in the interaction.

CLNUM          Set to T if you want to see the number of observations for each level of covariates. 
			   The default is T.

ORIENTATION	   orientation of the output Word table. Default is portrait, can be changed to landscape.
SHORTREPORT	   Use in cobmine with EFFECT and SLICEBY when there is an interaction in the model and set
			   to T to only report the stratified treatment effect.	
FILENAME       File name for output table.  This is necessary if report=T.
OUTPATH        File path for output table to be stored.  This is necessary if report=T.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
* For more details, please see the related documentation
**********************************************************************************************/

%macro logreg_sel(dsn=,outcome=,event=,desc=F,var=,cvar=,inc=0,slstay=.05, weight=, clnum=T,report=T,effect=,sliceby =,
shortreport=T, type3=T, outpath=,filename=,ORIENTATION = portrait, debug=F);

   /* Macros for final variable lists */
   %global _finalvar _finalcvar;

   %local FILENAME FOOTTEXT WEIGHT DSN REPORT OUTPATH CVAR_CNT CVAR I REMOVE OUTCOME CONTINUE
     __MACRO_ERR EVENT VAR DESC INC TYPE3 REMLIST DEBUG SLSTAY VAR_CNT length;

   /* Upper case T/F */
   %let report = %UPCASE(&report);
   %let debug = %UPCASE(&debug);
   %let type3 = %UPCASE(&type3);
   %let desc = %UPCASE(&desc);
   %let var = %UPCASE(&var);
   /*%let cvar = %UPCASE(&cvar);*/
   %let clnum=%UPCASE(&clnum);

   /* Count number of variables */
   %let var_cnt = %sysfunc(countw(&var,' '));

  /* Check for variable names that are longer than the default and increase length as needed */
  %let length = 64; 
      %DO i = 1 %to &var_cnt; 
      %IF %LENGTH(%SCAN(&var, &i,' ')) > &length %then %let length = %LENGTH(%SCAN(&var, &i,' ');
   %END;

 /*Check CVAR whether in the proper format, such as separate by * if reference level is specified*/
 /*Count number of class variables */

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

					%if &cnt_bef = 1 and  &pos3 = &len2 %then %let cvar_= &bef;
					%else %do; %put ERROR: The categorical variables in CVAR should be separated by * if you specify the reference level.; 
          		               %goto exit;%end;
					%end;
					
                %if &cnt2 = 0 and &cntw = 1 %then %LET cvar_ = &cvar_ref;

				%if &cnt2 > 1 or ( &cntw > 1 and &cnt2 = 0) %then %do; %put ERROR: The categorical variables in CVAR should be separated by * if you specify the reference level.; 
          		%goto exit;  %end;	 

			%let cvarlist = &cvarlist &cvar_; 
 %end;
		 
 %end;

   /* Make sure that CVAR are in VAR: to be test later oct2016*/
   %do i = 1 %to &cvar_cnt;
		%let nn= %sysfunc(findw(%upcase(&var), %upcase(%SCAN(&cvarlist, &i))));
		%put &nn;
		%if &nn = 0 %then %do; %put ERROR: The varaible appears in CVAR should also be in VAR.; %goto exit; %end;
	%end;

   /* Intialize */
   %let continue = 1;
   /* List of variables removed */
   %let remlist =;

   /* Run one step of backwards selection at a time and stop once there are no more variables*/
   /*  that should be removed. */
   %do %while(&continue=1);

   ODS SELECT NONE;
      ODS OUTPUT ClassFreq=clfreq ModelBuildingSummary=_removed NObs=nobs ParameterEstimates=estimate
          CLparmWald=CI ModelInfo=modelinf ResponseProfile=resp 
		  %if &effect ~= %STR() %then %do; SliceDiffs = slices_diff %end;
          %if %sysevalf(%superq(CVAR)~=,boolean) %then %do; Type3=type3 %end;;

      proc logistic data=&dsn namelen=&length simple;

         class %if %sysevalf(%superq(CVAR)~=,boolean) %then %do; %sysfunc(TRANSLATE(&cvar,' ','*')) %end;/order=internal param=glm;
           /* Note that if EVENT and DESC are specified that EVENT will override DESC */
         model &outcome(%if %sysevalf(%superq(event)~=,boolean) %then %do; event=&event %end; 
               order=internal %if &desc = T %then %do; desc %end;) = &var/selection=backward include=&inc 
               stop=%EVAL(&var_cnt-1) hierarchy=single
               slstay=&slstay CLPARM = wald;
           %if &weight ~= %STR() %then %do; 
              weight &weight/normalize;
           %end;
		   %if &effect ~= %STR() %then %do;    
			slice &effect * &sliceby/sliceby = &sliceby diff=control cl exp;
		   %end;
      run;
	  ODS SELECT ALL;

       /* PROC LOGISTIC creates the removed data set whether or not any variables are removed */
      /* If no variables were removed then the selection process is done */
       /* Get name of variable removed */
      PROC SQL noprint;
         select UPCASE(EffectRemoved), count(*)
         into :remove, :continue
         from _removed;
      quit;

       /* Update variable list with selected vars only */
      /* Note that the order of variables in an interaction term can be revsered so the */
      /* method used to update the categorical variable list cannot be used.  The method */
      /* below overcomes this problem. */

       /* Save order */
       DATA _est2;
          set estimate;
          order = _n_;
       RUN;

       /* Get unique list of variable names */
       PROC SORT DATA = _est2 nodupkey;
          by variable;
         where variable ~= 'Intercept';
      RUN;

       /* Return to original order */
       PROC SORT DATA = _est2;
          by order;
       RUN;

       /* Update variable list */
      PROC SQL noprint;
           select UPCASE(variable) into: var separated by ' '
           from _est2;
      QUIT;

       %if &continue = 1 %then %do;

          %let remlist = &remlist &remove;

          %if &debug = T %then %do;
            %put remlist &remlist;
           %end;

         /* Update class variable list */
		 /* update &cvar;*/
		 %if &cvar_cnt >  0 %then %do;* update &cvar;

		    %let pos_=%qsysfunc(findc(%superq(cvar), '*' ));
			
		    %let newcvar = ;

		 	%do i = 1 %to &cvar_cnt;
  				%if &pos_ > 0 or &cvar_cnt = 1 %then %let scanvar = %SCAN(&cvar, &i, '*');
				%else %let scanvar = %SCAN(&cvar, &i, ' ');

				%let pos = %qsysfunc(findc(%superq(scanvar), '(' ));

				%if &pos > 0 %then %let trimvar=%qsysfunc(SUBSTR(%superq(scanvar), 1,%eval(&pos-1)));
				%else %let trimvar = &scanvar;

				%if %upcase(&trimvar) = %upcase(&remove) %then %let newcvar = &newcvar; 
				%else  %let newcvar = &newcvar*&scanvar; 

			%end;/*end of %do i = 1 %to &cvar_cnt;*/

			%let cvar = &newcvar;
            %let cvar_cnt = %sysfunc(countw(&cvar,'*'));  * update number of variable in &cvar; 

		%end;/*%do %if (&var_cnt > 0);*/
	%end;/*end of %if &continue = 1*/
 
         /* Update variable counts */
         %let var_cnt = %EVAL(&var_cnt-1);

         /* Delete dataset */
         proc datasets lib=work noprint;
            DELETE removed;
         quit;  
  %end;


   /* Save selected vars in macro variables */
   %let _finalvar = &var;
   /* Remove asterisks */
   %if %sysevalf(%superq(CVAR)~=,boolean) %then %let _finalcvar = %sysfunc(TRANSLATE(&cvar,' ','*'));
   %else %let _finalcvar = &cvar;

   /* Produce report of final model */
   %if &report = T %then %do;

      /* If any variables were removed from the model create a list */
      %if &remlist ~= %STR() %then %do;

         /* Get variable labels for footnote */
	     ODS select none;
         PROC CONTENTS DATA = &dsn (keep=&remlist) out=_cont;
         RUN;
		 ods select all;

         DATA _cont;
            set _cont end=last;
             if label = ' ' then label = name;
             if _N_ ~= 1 and last then label = 'and ' || label;
         RUN;

         PROC SQL noprint;
            select label into :remLab separated by ', '
             from _cont;
         QUIT;

           %let foottext = The following variables were removed from the model: &remLab..;

       %end;
       %else %let foottext = No variables were removed from the model.;


      /* Remove outer quotes before macro call */
       %if %sysevalf(%superq(event)~=,boolean) %then %do;
          %let event = %sysfunc(DEQUOTE(&event));
       %end;
		
      %MULTIPLE_LOGREG(
          event=%BQUOTE(&event),
          OUTPATH = &outpath, 
          FNAME = &filename,
		  effect = &effect,
		  sliceby = &sliceby,
		  clnum=&clnum,
          FOOTNOTE="** Backward selection with an alpha level of removal of &slstay was used.  &foottext",
		  ORIENTATION = &ORIENTATION,
          TYPE3 = &type3,
		  shortreport = &shortreport,
          debug=&debug);

   %end;

   /* Only delete files if not in debug mode */
   %if &debug ~= T %then %do;
   ODS SELECT NONE;
      proc datasets lib=work memtype=data noprint;  
         delete CI estimate nobs _removed type3 resp modelinf %if &report = F %then %do;
               _est2 %end;;
      quit;  
	ODS SELECT ALL
   %end;

   /* Final variables selected (get rid of double spaces */
   %put Categorical variables selected: &_finalcvar;
   %put All variables selected: &_finalvar;


   %exit:
%mend logreg_sel;

/*********************************************************************************************/
