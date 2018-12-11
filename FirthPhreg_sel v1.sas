
/**********************************************************************************************
Macro: FirthPhreg_sel 
Created Date/Author/Contact: 04/04/2018/ Chao Zhang, and Yuan Liu 


Current Version: V1
Working Environment: SAS 9.4 English version


 Cox models probably lead to biased estimation of regression coefficients when survival dataset has too few outcome events.
 The phenomenon of monotone likelihood occurs in the fitting process of a Cox model when at least one parameter estimation 
 diverges to ±8. Firth’s penalized partial likelihood Cox regression mode to analyze survival dataset with rare events.

Purpose:  To conduct backward selection on a Firth’s penalized partial likelihood approach in Cox regression models 
, and a table of the resulting model can be generated.  


Parameters: 

DSN            The name of the data set to be analyzed.
EVENT          Time to event outcome variable.   

CENSOR         Name of censoring indicator variable.  Values of 0 indicate censored 
               observations.                             
VAR            List of variables to include in the model separated by spaces.  

CVAR           List of categorical variables to include in the model separated by spaces.  
               These should also appear in the var parameter.  If you want to change the 
               reference group you can follow each variable name by (DESC) where needed.  
               However, you will need to separate terms with an asterisk instead of a space. 

INC            Number of variables to include in the model (optional).  The first n variables 
               in the var parameter will be included in every model.  The default value is 0. 
 
ALPHA          the criterion for retaining a variable in the backward elimination procedure.

REPORT         Set this to T if you want a table of the resulting model generated (optional).  
               The default value is F.  
TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.  This only has an effect if REPORT = T.

CLNUM          Set to T if you want to see the number of observations for each level of covariates. The default is T.

ORIENTATION	   orientation of the output Word table. Default is portrait, can be changed to landscape.
	
FILENAME       File name for output table.  This is necessary if report=T.
OUTPATH        File path for output table to be stored.  This is necessary if report=T.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
* In this version, we allow the user to specify the reference level of a categorical variable.
* For more details, please see the related documentation
**********************************************************************************************/
%macro FirthPhreg_sel(dsn=, event=, censor=, var=, cvar=, inc=0, alpha=.2, type3=T, report=T,ORIENTATION = portrait,clnum=t, outpath=, filename=, debug=F);
   
  %local cnt ;
       
   /* Macros for final variable lists */
    %global _finalvar _finalcvar;
   /* Save current options */
    PROC OPTSAVE out=_options;   RUN;

   /* Upper case T/F */
   %let report = %UPCASE(&report);
   %let debug = %UPCASE(&debug);
   %let type3 = %UPCASE(&type3);
   %let var = %UPCASE(&var);
   %let cvar = %UPCASE(&cvar);

   /* Count number of variables */
   %let var_cnt = %sysfunc(countw(&var,' '));
   /* Count number of class variables */
   /* If asterisks were used as seperators or reference orders were reveresed (if only one */
   /* class then asterisks won't be used) */
   %if %INDEX(&cvar,*) > 0 OR %INDEX(&cvar,(DESC)) > 0 %then %let cvar_cnt = %sysfunc(countw(&cvar,*));
   %else %let cvar_cnt = %sysfunc(countw(&cvar,' ')); 

   /* Initialize maximum number of characters in model terms */
   /* It could be longer if there are interaction terms */
   %let length = 32; 

   /* Initialize error variable */
   %LET __Macro_Err = 0;

  /* Time Should not be missing value*/
   
   %if &event = %STR()  %then %do;
      %put ERROR: EVENT should be specified.; 
      %let __Macro_Err=1;
   %end;

  %DO i = 1 %to &var_cnt; 
      /* Check for variable names that are longer than the default and increase length as needed */
      %IF %LENGTH(%SCAN(&var, &i,' ')) > &length %then %let length = %LENGTH(%SCAN(&var, &i,' ');

       /* Check for interaction and report request */
       %if &report = T %then %do;
          %IF %INDEX(%SCAN(&var, &i, %STR( )),*) > 0 %then %do;
            %put ERROR: A report cannot be created if interaction terms are in the model.  Set REPORT=F.; 
              %let __Macro_Err=1;
          %END;
       %end;
   %END;

   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

/* Intialize */
   %let continue = 1;
   /* List of variables removed */
   %let remlist =;

   /* Run one step of backwards selection at a time and stop once there are no more variables*/
   /*  that should be removed. */
   %do %while(&continue=1);
    
   	   
    ods select none;
   ODS OUTPUT  NObs = numobs ClassLevelFreq =clfreq ParameterEstimates = MLEC ModelInfo=modelinf ModelANOVA = type3  ;
   proc phreg data=&dsn simple namelen=&length ;
         class %if %sysevalf(%superq(CVAR)~=,boolean) %then %do; 
              %sysfunc(TRANSLATE(&cvar,' ','*')) %end;/order=internal param=glm;
             
        model &event*&censor(0) = &var /firth  risklimits= pl  ;
		
    run;
	ods select all;

	/* Stop if number of variables left is equal to number of variables to include */
      %if &var_cnt = &inc %then %do;
         %let continue = 0;
      %end;
      %else %do;

	     /* See if any effects need to be removed - only consider effects that are not being in the model */
          PROC SORT DATA = type3 (firstobs=%EVAL(&inc+1)) out = _removed (where=(probchisq > &alpha));
             by DESCENDING ProbChiSq;
          RUN;

	    /* Select effect with larges p-value for removal from the model */
         /* If no variables were removed then the selection process is done */
         /* Get name of variable removed */
         PROC SQL noprint;
             select UPCASE(Effect), count(*)
             into :remove, :continue
             from _removed (OBS=1);
         quit;
      %end;

	  %if &continue = 1 %then %do;

          %let remlist = &remlist &remove;

          %if &debug = T %then %do;
             %put remlist &remlist;
          %end;

    
         /* Update variable list with selected vars only */
         /* Note that the order of variables in an interaction term can be revsered so the */
         /* method used to update the categorical variable list cannot be used.  The method */
         /* below overcomes this problem. */

         /* Save order */
          DATA _est2;
             set Mlec;
             order = _n_;
			/* Create another variable to be used for consistency */
			 parm = parameter ;
          RUN;

       /* Get unique list of variable names */
         /* Exclude removed effect */
         PROC SORT DATA = _est2 nodupkey;
             by parm;
             where  UPCASE(parm) ~= "&remove";
         RUN;

         /* Return to original order */
         PROC SORT DATA = _est2;
            by order;
         RUN;

		 /* Update variable list */
         PROC SQL noprint;
            select UPCASE(parm) into: var separated by ' '
            from _est2;
         QUIT;

      /* Update class variable list */
	    %if &cvar_cnt >0 %then %do;
         %let newcvar =;
         %do i = 1 %to &cvar_cnt;

            /* If asterisks were used as seperators */
            %if %INDEX(&cvar,*) > 0 OR %INDEX(&cvar,(DESC)) > 0 %then %do;
               /* Variable name + order */
               %let c&i = %SCAN(&cvar, &i, *);
            %end;
            %else %do;
               %let c&i = %SCAN(&cvar, &i, ' ');
            %end;                  
            /* Remove from list if it is the removed variable */
            /* Rescan to just get variable name not order as well */
            %if %BQUOTE(%SCAN(&&c&i, 1, %STR( ))) ~= %BQUOTE(&remove) %then %do;

               %if %INDEX(&cvar,*) = 0 AND %INDEX(&cvar,(DESC)) = 0 %then %do;
                  /* Recombine individual variables into a new list */
                  %let newcvar = &newcvar &&c&i; 
               %end;
               %else %do;
                  /* If first item in list don't add asterisk */
                  %if %sysevalf(%superq(NEWCVAR)=,boolean) %then %do;
                     %let newcvar = &&c&i; 
                  %end;
                  %else %do;
                     %let newcvar = &newcvar * &&c&i;
                  %end;
               %end;
            %end;

         %end;
      
         /* Reset var list */
          %let cvar = &newcvar;
         /* Trim extra blanks */
          %let cvar = %SYSFUNC(TRANWRD(&cvar,%STR(  ),%STR()));
        
         /* Update variable counts */
       
         /* Count number of class variables */
         /* If asterisks were used as seperators */
         %if %INDEX(&cvar,*) > 0 OR %INDEX(&cvar,(DESC)) > 0 %then %let cvar_cnt = %sysfunc(countw(&cvar,*));
         %else %let cvar_cnt = %sysfunc(countw(&cvar,' ')); 

       %end;
  %let var_cnt = %EVAL(&var_cnt-1);
         /* Delete dataset */
         proc datasets lib=work;
            DELETE _removed;
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
         PROC CONTENTS DATA = &dsn (keep=&remlist) out=_cont;
         RUN;

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
   %end;
  %end;

%MULTIPLE_PHREG_Firth(        
          OUTPATH = &outpath, 
          FNAME = &filename,
		  ORIENTATION = &ORIENTATION,
		  clnum=&clnum,
          TYPE3 = &type3,
          debug= &debug,
          FOOTNOTE="** Backward selection with an alpha level of removal of &alpha was used.  &foottext"
        );

 
  /* Final variables selected (get rid of double spaces */
   %put Categorical variables selected: &_finalcvar;
   %put All variables selected: &_finalvar;


%mend FirthPhreg_sel;

