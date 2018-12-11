/**********************************************************************************************
Macro: GENMOD_SEL
Created Date/Author: Dec. 11, 2013/Dana Nickleach 
Last Update Date/Person/Contact: 
Current Version: V1
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu 

Purpose:  To conduct backward selection on a multivariable generalized linear regression 
model using the maximum possible sample size at each stage of the selection process instead of 
restricting to the sample size from the first step as SAS does when using their selection 
methods.  The normal, binomial, Poisson, or negative binomial distributions can be used.  An 
identity link will be used with a normal distribution, a logit link with binomial, and a log 
link with Poisson or negative binomial. Generalized estimating equations (GEE) may also be 
used.  Optionally, a table of the resulting model can be generated.  

Notes: The model is run using PROC GENMOD.  The final list of variables selected will be 
written to the log.  Additionally, two global macro variables, _finalvar and _finalcvar will be
created containing the list of all variables and categorical variables selected, respectively.
If you are requesting a table with the model results then the macro “MUTLIPLE_GENMOD V7” or 
later is also required.  Interactions can be included if you are not requesting a table 
(REPORT=F).  However, model hierarchy will not be maintained so interactions should only be 
used with caution.  Also, the correct results might not be produced if interactions terms are
included without their main effects and the macro cannot produce a report table if the model
contains interactions.  

Parameters: 

DSN            The name of the data set to be analyzed.
OUTCOME        The name of the continuous or binary outcome variable. 
DESC           Set to T to reverse the order of a binary outcome (optional).  The default value 
               is F. 
VAR            List of variables to include in the model separated by spaces.                                            
CVAR           List of categorical variables to include in the model separated by spaces.  These 
               should also appear in the VAR parameter.  If you want to change the reference 
               group you can follow each variable name by (DESC) where needed.  However, you 
               will need to separate terms with an asterisk instead of a space. 
SUBJECT        Subject-effect to be specified in REPEATED statement.  Leave blank if not using
               a GEE model.
WITHINSUB      Within-subject-effect to be specified in REPEATED statement.  This is optional
               if using a GEE model.  Leave blank if not using a GEE model.   
DIST           Probability distribution to use for the model (optional).  Valid values are 
               NORMAL, BINOMIAL, POISSON, and NEGBIN.  The default value is NORMAL.  
TYPE           Correlation structure keyword to be specified in REPEATED statement.  Leave 
               blank if not using a GEE model.
INC            Number of variables to include in the model (optional).  The first n variables in 
               the VAR parameter will be included in every model.  The default value is 0.  
SLSTAY         The significance level for removing variables from the model (optional).  The 
               default value is .05.   
TYPE3          Set to F to suppress type III p-values from being reported in the table 
               (optional).  The default value is T.  This only has an effect if REPORT = T.
REPORT         Set this to T if you want a table of the resulting model generated (optional).  
               The default value is F.
OUTPATH        File path for output table to be stored.  This is necessary if REPORT=T. 
FNAME          File name for output table.  This is necessary if report=T.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
* For more details, please see the related documentation
**********************************************************************************************/

%macro genmod_sel(dsn=, outcome=, desc=F, var=, cvar=, subject=, withinsub=, type=, 
   dist=NORMAL, inc=0, slstay=.05, type3=T, report=T, outpath=, fname=, debug=F);

   /* Macros for final variable lists */
   %global _finalvar _finalcvar;

   %local FILENAME FOOTTEXT DSN SUBJECT REPORT OUTPATH CVAR_CNT CVAR I REMOVE OUTCOME CONTINUE
     __MACRO_ERR REMLAB VAR EVENT DESC REPEATED TYPE INC DIST TYPE3 REMLIST FNAME DEBUG NEWCVAR
     C1 SLSTAY C2 VAR_CNT;

   /* Upper case T/F */
   %let report = %UPCASE(&report);
   %let debug = %UPCASE(&debug);
   %let type3 = %UPCASE(&type3);
   %let desc = %UPCASE(&desc);
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

   /* Subject and correlation type must be specified to use repeated statement */
   %if &subject = %STR() OR &type = %STR() %then %let repeated = F;
   %else %let repeated = T;

   /* Initialize error variable */
   %LET __Macro_Err = 0;

   /* For GEE - both subject and type should be specified */
   %if (&subject = %STR() AND &type ~= %STR()) OR 
     (&subject ~= %STR() AND &type = %STR())%then %do;
      %put ERROR: When using a GEE model both SUBJECT and TYPE need to be specified.;
      %let __Macro_Err = 1;
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

   /* If there is an error in the parameters supplied then exit */
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

      ODS OUTPUT ModelInfo=modelinf NObs = nobs Type3 = type3 
          %if &repeated = T %then %do; GEEEmpPEst = estimate %end;
          %else %do; ParameterEstimates=estimate %end;;
      PROC GENMOD DATA = &dsn %if &desc = T %then %do; DESC %end; namelen=&length;
         /* Remove Asterisks from class list */
         class &subject &withinsub %if %sysevalf(%superq(CVAR)~=,boolean) %then %do; %sysfunc(TRANSLATE(&cvar,' ','*')) %end; /ORDER=internal;
         model &outcome = &var/type3 wald dist=&dist;
         %if &repeated = T %then %do; 
            REPEATED SUBJECT=&subject/type=&type
         %if &withinsub ~= %STR() %then %do; withinsubject=&withinsub %end;;
         %end;
      RUN;

      /* Stop if number of variables left is equal to number of variables to include */
      %if &var_cnt = &inc %then %do;
         %let continue = 0;
      %end;
      %else %do;
         /* See if any effects need to be removed - only consider effects that are not being */
         /* in the model */
         PROC SORT DATA = type3 (firstobs=%EVAL(&inc+1)) out = _removed (where=(probchisq > &slstay));
            by DESCENDING ProbChiSq;
         RUN;

         /* Select effect with larges p-value for removal from the model */
         /* If no variables were removed then the selection process is done */
         /* Get name of variable removed */
         PROC SQL noprint;
            select UPCASE(source), count(*)
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
            set estimate;
            order = _n_;
            %if &repeated ~= T %then %do;
               /* Create another variable to be used for consistency */
               parm = parameter;
            %end;
         RUN;

         /* Get unique list of variable names */
         /* Exclude removed effect */
         PROC SORT DATA = _est2 nodupkey;
            by parm;
            where parm not in ('Intercept' 'Scale' 'Dispersion') and UPCASE(parm) ~= "&remove";
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
         %let var_cnt = %EVAL(&var_cnt-1);
         /* Count number of class variables */
         /* If asterisks were used as seperators */
         %if %INDEX(&cvar,*) > 0 OR %INDEX(&cvar,(DESC)) > 0 %then %let cvar_cnt = %sysfunc(countw(&cvar,*));
         %else %let cvar_cnt = %sysfunc(countw(&cvar,' ')); 

         /* Delete dataset */
         proc datasets lib=work;
            DELETE _removed;
         quit;  
      %end;

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

      /* Produce report */
      %MULTIPLE_GENMOD(OUTPATH=&outpath,
          FNAME=&fname,
          FOOTNOTE="Backward selection with an alpha level of removal of &slstay was used.  &foottext",
          TYPE3=&type3,
          DEBUG=&debug);
   %end;

   /* Only delete files if not in debug mode */
   %if &debug ~= T %then %do;
      proc datasets lib=work memtype=data noprint;  
         delete estimate nobs _removed type3 modelinf %if &report = F %then %do;
               _est2 %end;;
      quit;  
   %end;

   /* Final variables selected (get rid of double spaces */
   %put Categorical variables selected: &_finalcvar;
   %put All variables selected: &_finalvar;

%mend genmod_sel;

/*********************************************************************************************/
