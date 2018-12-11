/**********************************************************************************************
***********************************************************************************************
Macro Name: MULTIPLE_PHREG
Created Date/Author: Feb. 2012/Yuan Liu
Last Update Date/Person: Oct, 2016/ Yaqi Jia
Current Version: V21
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu

Purpose:  To produce a summary table from a multivariable survival model, but the macro itself 
does not do statistical modeling automatically.  

Notes: The macro only helps to generate a nice summary table, and cannot be used for 
statistical modeling.  The modeling must be conducted using PROC PHREG before calling this 
macro.  Several tables must be created using ODS OUTPUT when running PROC PHREG (nNObs = numobs 
ParameterEstimates = MLEC ModelInfo=modelinf SliceDiffs = slices_diff(only for interaction model)
Type3 = type3). The option RL must also be used in the model statement and option PARAM= GLM in 
class statement. Variable names must not be more than 20 characters.  This macro can be set up to 
handle models including interaction terms, and in the model fitting step, use SLICE statement to specify 
how the stratified analysis to be conducted and use ods output SliceDiffs = slices_diff.   

Parameters: 

OUTPATH   	   Path for output table to be stored.
FNAME     	   File name for output table.
FOOTNOTE  	   Text of footnote to include in the table (optional).  It should be in quotes.  This 
          	   footnote will appear below the footnote containing the number of observations.  
          	   Leave this field blank if not including an additional footnote.
FOOTNOTE2 	   Text of second footnote to include in the table (optional).  It should be in quotes.  
         	   This footnote will appear below all other footnotes.  Leave this field blank if not
           	   including an additional footnote.   
TYPE3     	   Set to F to suppress type III p-values from being reported in the table (optional).
               The default value is T.
CLNUM          Set to T if you want to see the number of observations for each level of covariates. The default is T.
EFFECT         If not empty, then the model contains a two-way interaction. It is used to specify
			   the treatment effect in the interaction.
SLICEBY		   Use in combine with EFFECT to specify the stratified variable in the interaction.
SHORTREPORT	   Use in combine with EFFECT and SLICEBY when there is an interaction in the model and set
			   to T to only report the stratified treatment effect.	
ORIENTATION	   orientation of the output Word table. Default is portrait, can be changed to landscape.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted in
               debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
For more details, please see the related documentation
**********************************************************************************************/

%macro MULTIPLE_PHREG(OUTPATH=, FNAME =, FOOTNOTE=, FOOTNOTE2=, TYPE3=T, EFFECT=,SLICEBY=,CLNUM=T,STRATA=, shortreport=T,
ORIENTATION = portrait, DEBUG=F);

   %local FOOTNOTE FOOTNOTE2 NUMREAD NUMUSED OUTPATH DSN EVENT VAR CLASSEXIST N CONTRAST
     EVENTVV CLIST TYPE3 FNAME DEBUG;

   /* Upper case */
   %let type3 = %UPCASE(&type3);
   %let debug = %UPCASE(&debug);
   %let strata = %UPCASE(&strata);

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

   /* Get data set name and outcome name */
   PROC SQL noprint;
      select cvalue into :dsn
      from modelinf
      where description = 'Data Set';
      select cvalue into :event
      from modelinf
      where description = 'Dependent Variable';
   QUIT;

   /* Get outcome label */
   data _NULL_; 
      set &dsn (obs=1); 
      call symput('eventvv', VLABEL(&event));
   run;

      /* Get variable labels */
   PROC CONTENTS DATA = &dsn out=_cont (rename=(name=parameter)) noprint;
   RUN;


   /* Check to see if class variables were used - s/b classval0 variable */
   data _null_;
      dsid=open('mlec');
      check=varnum(dsid,'classval0');
      /* 0 = does not exist, 1 = exists */
      exist = MIN(check, 1);
      CALL SYMPUT('classExist',PUT(exist,1.));
   run;

   /* If class variables were used */
   /* Don't include interaction terms in categorical list */
   %if &classexist = 1 %then %do;
      /* Create list of class variables */
      PROC SQL noprint;
         select distinct Parameter into :clist separated by ' ' 
         from mlec
         where MISSING(classval0) = 0 and INDEX(parameter,'*') = 0 ;
      QUIT;
   %end;
   %else %let clist =;
   
   %put &clist;

   proc sql noprint;
      select NObsRead into :numread  from numobs; 
      select NObsUsed into :numused from numobs; 
   quit;

   %put &numread &numused;

   DATA mlec;
      set mlec;
      /* Save order */
      Num= _N_;
   RUN;

   %IF &CLIST NE %THEN %DO; 
     
       /*Merge mlec and clfreq*/

        /*prepare to be sorted. there is missing value in class variable in clfreq, so need below data step*/
       data clfreq; 
       length parameter $256.;
       retain parameter;
       set clfreq;
       if class ~= ' ' then parameter=class;
	   else class=parameter;
	   keep parameter class value freq;
       run;

       proc sort data=clfreq; by parameter value;run;
	   proc sort data=mlec;  by parameter classval0;run;

	     %if &strata ~= %str() %then %do; /*when the strata statement used in the model, the sample size is generated by strata, reduce it to total samples through strata*/
			  ods select none;
			  ods output Summary=sums_(rename=(freq_Sum=Freq));
			  proc means data=clfreq sum;
			  var Freq; by parameter Value;
			  run; ods select all; 

			  data clfreq; set sums_;run;
         %end;

       data mg;
	     length parameter $256. level $256.;
	     merge mlec(in=a rename=(classval0=level)) clfreq(in=b rename=(value=level));
	     by parameter level;
		 if a;
		 drop label; /*avoid conflict with label variable in lab data set produced by proc contents below*/
	   run;
        
	    /*Merge mg with label data set produced by proc contents*/
		ods select none;
        ods output variables=lab;
        proc contents data= &dsn ;run;
		ods select all;
        
		proc sort data=lab;by variable;run;
		proc sort data=mg;  by parameter;run;

		data mglab;
		  merge mg (in=a) lab (in=b rename=(variable=parameter) keep=variable label);
		  by parameter;
		  if a;
		run;
       
	   /*Merge type3 and mg*/
	   data type3;
           length effect $100;
           set type3; 
           order=_n_; 
           rename ProbChiSq = p_type3 effect=parameter; 
           format ProbChiSq;
           keep effect ProbChiSq order;
      run;

       proc sort data=type3;by parameter;run;
	   proc sort data=mglab;by parameter;run;

      data mergedata;
	       length parameter $256. level $256. plabel $256. n $120.;
           merge mglab type3;
           by parameter;
		   n=put(freq,8.);
           plabel=label;
           if label = "" and plabel = "" then plabel = parameter;
           format ProbChiSq;
		   if n=. then n="&numused";

           /* If coding was ref then reference parameters will not have been in the original */
         /* parameter estimates file and will have missing num - Make sure they are ordered */
           /* as the last category for a variable */
           if num = . then num = 1000;
      run;

       proc sort data=mergedata;
         by order num;
      run; 

   %END;

   %ELSE %DO;

      /* Check to see if label variable exists (it won't exist if none of the model vars have */
      /* labels */
      data _null_;
         dsid=open('mlec');
         check=varnum(dsid,'label');
           /* 0 = does not exist, 1 = exists */
           exist = MIN(check, 1);
         CALL SYMPUT('labexist',PUT(exist,1.));
      run;

      data mergedata;
         length plabel $256;
         set mlec;
         Level = "";
         p_type3 = ProbChiSq;

         /* If there are no labels use var names */
         %if &labexist ~= 0 %then %do; 
            if label ~= ' ' then plabel = label; 
               else plabel = parameter;
         %end;
         %else %do;
            plabel = Parameter;
         %end;

         format p_type3 ProbChiSq;
      run;

   %END;

 

   /* If reporting contrasts as well */
  %if &effect ~=%str() %then %do;

      	/*Extract variable name,label and length*/
       %let cfivar=%sysfunc(countw(&_finalvar,' '));
        %put finalvar=&_finalvar;
        %put number of final var=&cfivar; 
        data frequency;
		  set &dsn;
		  %do i=4 %to &cfivar;
		  if missing(%scan(&_finalvar,&i,' '))=0;
		  %end;

		  %if &event ~= %STR() %then %do;
		    where &event >= 0 and missing(&censor)=0;
		  %end;

		  %else %do;
		    where &stop >= &start and missing (&start)=0 and missing(&stop)=0 and missing(&censor)=0;
		  %end;
		run;
        
		ods output list=freq_;
        proc freq data=frequency;
		   tables &sliceby*&effect/list;
	    run;

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
		HazardRatio = ExpEstimate;
		HRLowerCL = LowerExp;
		HRUpperCL = UpperExp;
		ProbChiSq = Probz;
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
		 set type3;
         plabel = "Comparisons Stratified by &Slicebylabel :";
		 level="&effectlabel :";
		 where index(parameter,'*') ~=0;
		 keep parameter plabel level p_type3;
	  run;
	  
	  /* Combine Type3 data set and Hazardratio data set*/
      DATA _contrast;
         set _type3 merge12;
         keep plabel level n HazardRatio HRLowerCL HRUpperCL probchisq p_type3;
      RUN;

	  %if &shortreport = T %then %do;

	     %let upeffect=%upcase(&effect); /* delete &effect &sliceby and &effect*&sliceby variables in above data set*/
		 %let upsliceby=%upcase(&sliceby);

	      DATA conlab;
          set mergedata; 
		  if upcase(parameter) not in ("&upeffect"  "&upsliceby" "&upeffect*&upsliceby" );
          RUN;

          %let conlab= None;
		  proc sql noprint;
		    select distinct(plabel) into: conlab separated by ', ' from conlab;
		  quit;
	    
        Data mergedata;
	      set _contrast;run; 
       %end;

	   %else %do;
	      %let upeffect=%upcase(&effect); /* delete &effect &sliceby and &effect*&sliceby variables in above data set*/
		  %let upsliceby=%upcase(&sliceby);
	      DATA mergedata;
          set mergedata _contrast; 
		  if upcase(parameter) not in ("&upeffect"  "&upsliceby" "&upeffect*&upsliceby" );
          RUN;
	   %end;

   %end;

   DATA mergedata;
      set mergedata;
      /* HR and 95% CI */
      if HazardRatio ~= . then HR = TRIM(LEFT(PUT(HazardRatio,8.2))) || " (" || 
         TRIM(LEFT(PUT(HRLowerCL,8.2))) || "-" || 
         TRIM(LEFT(PUT(HRUpperCL,8.2))) || ")";
      else HR = '-';
      /* Number of observations read and used - save in data set */
      numused = &numused;
      numread = &numread;
   RUN;

   *---- table template -----;  
   ODS PATH WORK.TEMPLAT(UPDATE) SASUSR.TEMPLAT(UPDATE) SASHELP.TMPLMST(READ);

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

   PROC REPORT DATA=mergedata HEADLINE HEADSKIP CENTER STYLE(REPORT)={JUST=CENTER} SPLIT='~' 
     nowd SPANROWS lS=256; 
      COLUMN plabel  Level %if &clnum=T %then %do;n %end;("&eventvv" '------------------------------------------'
          (HR  ProbChiSq  %if &type3 = T %then %do; p_type3 %end;)); 
      DEFINE plabel/order order=data "Covariate"  STYLE(COLUMN) = {JUST = L CellWidth=25%}; 
      DEFINE Level/DISPLAY "Level" STYLE(COLUMN) = {JUST =L CellWidth=20%}; 
	  %if &clnum=T %then %do;
	  DEFINE n/display "N" style(column)={just=C cellwidth=10%};
	  %end;
      %if &type3 = T %then %do; 
      DEFINE p_type3/order MISSING "Type3 P-value" STYLE(COLUMN) = {JUST=C CellWidth=8%} 
               FORMAT=PVALUE8.3; 
       %end;
      DEFINE HR/DISPLAY "Hazard Ratio (95% CI)" STYLE(COLUMN) = {JUST = C CellWidth=15%}; 
      DEFINE ProbChiSq/DISPLAY "HR P-value" STYLE(COLUMN) = {JUST = C CellWidth=8%} FORMAT=PVALUE8.3; 

      COMPUTE ProbChiSq; 
         IF . < ProbChiSq <0.05 THEN CALL DEFINE("ProbChiSq", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         *if plabel = "Stratified Comparisons by &Slicebylabel :" THEN CALL DEFINE("ProbChiSq", "STYLE", "STYLE=[COLOR=WHITE]");
      ENDCOMP; 

       %if &effect ~= %str() %then %do;
       COMPUTE plabel;
          if plabel = "Comparisons Stratified by &Slicebylabel :" THEN CALL DEFINE("plabel", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]");  
       ENDCOMP;

	   COMPUTE level;
          if level = "&effectlabel :" THEN CALL DEFINE("level", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]");  
       ENDCOMP;
       %end;
       /*COMPUTE HR;
          if plabel = "Stratified Comparisons by &Slicebylabel :" THEN CALL DEFINE("HR", "STYLE", "STYLE=[COLOR=WHITE]");
       ENDCOMP;*/

       %if &type3 = T %then %do;
         COMPUTE p_type3; 
            IF . < p_type3 <0.05 THEN CALL DEFINE("p_type3", "STYLE", "STYLE=[FONT_WEIGHT=BOLD]"); 
         ENDCOMP;   
       %end;

       compute after plabel; line '';   endcomp;  

       compute after _page_;
         length text1 - text2 $512.;
         text1 = "*  Number of observations in the original data set = %TRIM(&numread). Number of observations used = %TRIM(&numused).";
         line @0 text1 $512.;

           %if &footnote ~= %STR() %then %do;
              line @0 &footnote;
           %end;
           %if &footnote2 ~= %STR() %then %do;
              line @0 &footnote2;
           %end;
		   %if &effect ~=%str() and &shortreport = T %then %do;
			line @0 "*** The estimated stratified treatement effect was controlled by: &conlab";
		   %end;
      ENDCOMP;
 
   RUN;

   ods rtf close;

   /* Reload original options that were in use before running the macro */
   PROC OPTLOAD data=_options;
   RUN;

   %if &debug = F %then %do;
      *--- DELETE ALL TEMPORARY DATASETS that were created; 
      proc datasets lib=work memtype=data noprint;  
          delete mergedata _options _cont %if &CLIST ~= %STR() %then %do; clfreq lab mglab freq mg type3 Sums_ %end;
          %if &effect ~= %str() %then %do;freq freq_ frequency slices_diff _contrast _temp merge_1 merge_2 merge_12 merge12 _type3 %end;
          %if &shortreport = T %then %do; conlab %end;;

      quit;    
   %end;

%mend;
