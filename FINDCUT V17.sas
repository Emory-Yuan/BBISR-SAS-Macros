/**********************************************************************************************
***********************************************************************************************
Macro Name: FINDCUT
Created Date/Author Feb. 2012/Yuan Liu
Last Update Date/Person Dec. 9, 2013/Dana Nickleach 
Current Version: V17
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu 

Purpose:  To search for an optimal cut point for a continuous variable that predicts time to
event outcome and produce a martingale residual plot for the variable, which will further help 
us to detect any nonlinear relationship.  

Notes:  The major part of this macro was from http://www2.sas.com/proceedings/sugi28/261-28.pdf.
Please cite properly.  The log-rank rank statistic, not log-rank chi-square statistic, is used 
to choose the optimal cut point.  The cut point is chosen based on the data being dichotomized
into groups X< cut point vs. X>= cut point.

The macro MPLOT V11 or later and the macro NO_OUT V4 or later is required if OUTLIER=T or B.  
Besides the RTF document a SAS dataset containing the optimal cut points will be created.

Parameters: 

DS             The name of the data set to be analyzed.
TIME           Time to event outcome variables separated by spaces.                        
STAT           The variable name for censoring indicators separated by spaces. The order of 
               this list should correspond to the order of the TIME parameter.  It is required 
               that 1 is used for the event and 0 for censored cases.  
CUTVAR         The numerical variables to be plotted separated by spaces.  The variable name 
               cannot be more than 26 characters.
SM             Smoothing parameter.  The default value is 0.65.
DOC            Set to T to create a RTF file containing the output or F to suppress creation of 
               the RTF file (optional).  The default value is T.
OUTPATH        Path for output table to be stored
FNAME          File name for output table.
MINOUT         Request minimal output (optional).  Set to T to minimize output, which will only 
               include the martingale residual plot, percentiles, and optimal cut point.  The 
               default value if F.  
SUBTIT         Text to be added to the title to appear after, but on the same line as the 
               first title (optional).  Text should not be enclosed in quotes.  By default no 
               additional text will be printed.
PLOT           Set to F to suppress the martingale residual plots (optional).  The default value 
               is F.
HREF           Set to T for a horizontal reference line at zero on the plot (optional).  The 
               default value is F. 
VREF           Set to T for vertical reference lines at the median and optimal cut point in blue 
               on the plot (optional).  The default value is F.
CHOSEN         Value of the chosen cut points separated by a space to display with a vertical 
               reference line in red on the plot (optional).  If multiple cut points are listed 
               they will be displayed on the same plot.  This is not recommended if processing
               more than one probe or data set since the same points will be used for all plots.  
DEC            Fixed number of decimal places to display in the percentile table (optional).  
               Note that the maximum format width is 7.  So the number of decimal places can't 
               be more than 5, but depending on the size of the numbers may be less.  The best 
               format is used by default.
OUTLIER        Set to T, F, or B (optional).  T will exclude extreme outliers from the 
               analysis.  F will not exclude any points.  B will perform the analysis both 
               ways, including all points and excluding outliers.  Note that if no outliers are
               found the second analysis will not be reported.  The default value is F.  
OUT            The name of the SAS dataset that contains the optimal cut point results.  By 
               default it will be called Cutpoint.
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted 
               in debug mode.  This is useful if you are editing the code or want to further 
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
For more details, please see the related documentation
**********************************************************************************************/

%macro FINDCUT(ds=, time=, stat=, cutvar=, sm = 0.65, DOC=T, OUTPATH=, FNAME=, MINOUT=F, 
     SUBTIT=, PLOT=T, HREF=F, VREF=F, CHOSEN=, DEC=, OUTLIER=F, OUT=Cutpoint, debug=F);

   /* Make sure that these macro variables are local and can't be changed outside this macro */
   %local i j k __MACRO_ERR CENSOR VAR EVENT NOOUT_N NOBS EVELAB OUT_CNT VV VAR_CNT tit2 median
     cutpoint color thick rlabel pattern;

   /* Initialize error variable */
   %let __Macro_Err = 0;

   /* Count number of variables */
   %let var_cnt = %sysfunc(countw(&cutvar));
   /* Count number of outcomes */
   %let out_cnt = %sysfunc(countw(&time));
   /* Count number of chosen cut points */
   %if &chosen ~= %STR() %then %do; 
      %let chosen_cnt = %sysfunc(countw(&chosen));
   %end;
   %else %let chosen_cnt = 0;

   /* Check for variable names that are too long */
   %do i = 1 %to &var_cnt;
      %let var = %SCAN(&cutvar, &i);
      %if %LENGTH(&var) > 26 %then %do;
         %put ERROR: Variable name &var must be less than 27 characters.; 
           %let __Macro_Err=1;
      %end;
   %end;

   /* If there is an error in the parameters supplied then exit */
   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   /* Capitalize */
   %let minout = %UPCASE(&minout);
   %let HREF = %UPCASE(&HREF);
   %let VREF = %UPCASE(&VREF);
   %let debug = %UPCASE(&debug);
   %let outlier = %UPCASE(&outlier);
   %let doc = %UPCASE(&doc);
   %let plot = %UPCASE(&plot);

   /* Save current options */
   PROC OPTSAVE out=_options;
   RUN;

   PROC FORMAT;
      VALUE pf 
          -99='P>0.30'
          0-<.0001='<.0001';
   RUN;

   OPTIONS NODATE;
   ODS NOPROCTITLE;
   %if &doc = T %then %do;
      ODS RTF startpage = NEVER BODYTITLE FILE= "&OUTPATH.&FNAME &SYSDATE..DOC"; 
   %end;

   /* Repeat for each cutvar */
   %do k = 1 %to &var_cnt;
      %let var = %SCAN(&cutvar, &k);

       data _dds1;
         /* Drop unnecessary variables to speed processing */
         set &ds (keep=&time &stat &var);
         where &var ne .;
      run;

      /* Get label of cutvar */
      data _NULL_; 
         set _dds1 (obs=1); 
         call symput('vv', put(label(&var),$40.));
      run;

       /* Get total number of observations */
      PROC SQL noprint;
         select count(*) into :nobs
          from _dds1;
      QUIT;
      /* Initialize */
      %let noout_n = &nobs;

       /* Remove outliers */
      %if &outlier = T OR &outlier = B %then %do;
         %no_out(dsn=_dds1,vars=&var,debug=&debug);

          /* Total number of obs excluding outliers */
          PROC SQL noprint;
             select count(*) into :noout_n
              from _dds1noout&var;
          QUIT;
       
          /* Use only one data set */
          %if &outlier = T %then %do;
             /* Dataset without outliers */
             DATA _dds1;
               set _dds1noout&var;
            RUN;
          %end;
          /* Use a second data set */
          %else %if &outlier = B %then %do;
             /* Dataset without outliers */
             DATA _dds2;
               set _dds1noout&var;
            RUN;
          %end;
      %end;

      /* Number of time to calculate optimal cut point for one probe */
      /* If want with and without outliers then repeat twice */
      /* If no outliers were removed and option was BOTH - don't run both ways */
      %if &nobs = &noout_n %then %let num = 1;
      %else %if &outlier = B %then %let num = 2;
      %else %let num = 1;

      /* For each outcome */
      %do j = 1 %to &out_cnt;
         %let event = %SCAN(&time, &j);
         %let censor = %SCAN(&stat, &j);

         /* Get outcome label */
         DATA _NULL_;
            set _dds1;
            CALL SYMPUT('eveLab',VLABEL(&event));
         RUN;

         /* May need to repeat with and without outliers */
         %do i = 1 %to &num;
            proc sort data=_dds&i; 
               by &var;

            data _difage; 
               set _dds&i; 
               by &var;
               /*title3 "step1: check no. of disticnt &var";*/
               if first.&var;

            proc sort; by &var;

            data _ttt; set _difage; by &var;
            cut=_N_;
            keep &var cut;

            proc sort; by descending &var;

            data _cut; set _ttt; by descending &var;
            if _N_=1 then do;
            nocut=cut; retain nocut;
            end;
            somecut=&var;
            drop &var;
            output;

            proc sort; by somecut;
            /* proc print n; */
            data _dim; set _cut; by somecut;
            if last.somecut; dummy=1;
            keep nocut dummy;
            proc sort; by dummy;

            *********************************************;
            * find # of &stat at each time ti *;
            * # at risk at time ti *;
            * ti are those distinct event time *;
            *********************************************;

            proc sort data=_dds&i; by descending &event &censor;

            data _step1; set _dds&i; by descending &event &censor;
            /*title3 "check step1";*/
            if _N_=1 then do;
            norisk=0; retain norisk;
            end;
            if first.&event then do;
            nodeath=0; retain nodeath;
            end;
            norisk=norisk+1;
            if &censor=1 then nodeath=nodeath+1;
            if last.&event and &censor=1 then output;
            data _step1; set _step1;
            keep &event norisk nodeath;
            proc sort; by &event;
            /* proc print n; */
            *********************************************;
            * find # of &stat at each time ti *;
            * # at risk at time ti *;
            * == with c>=&cutvar *;
            *********************************************;
            data _dummy; set _dds&i;
            dummy=1;
            proc sort; by dummy;
            data _double; merge _dummy _dim; by dummy;
            do cut=1 to nocut;
            output;
            end;
            proc sort; by cut;
            proc sort data=_cut; by cut;

            /* DCN 7/18/12 - deleted obsno from keep statement */
            data _comb; 
            merge _double _cut; 
            by cut;
            keep &event &censor &var somecut cut;
            RUN;

            proc sort; by cut descending &event &censor;
            /* proc print n; */
            proc sort data=_comb; by cut descending &event &censor;

            data _step2 (KEEP=&event cut somecut noriskc nodeathc); 
               set _comb; by cut descending &event &censor;
               /*title3 "check step2";*/
               if first.cut then do;
                  noriskc=0; retain noriskc;
               end;
               if first.&event then do;
                  nodeathc=0; retain nodeathc;
               end; 
               if &var>=somecut then noriskc=noriskc+1;
               if &censor=1 and &var>=somecut then nodeathc=nodeathc+1;
               if last.cut or (last.&event and &censor=1) then output;
            RUN;

            /* proc print n; */
            proc sort; 
               by &event;
            RUN;

            *********************************************;
            * compute Sk .. max(sk) *;
            *********************************************;
            data _step3; merge _step2 _step1; by &event;
            /*title3 "step3";*/
            sik=nodeathc-nodeath*noriskc/norisk;
            /* proc print n; */
            proc sort; by cut somecut;
            proc univariate noprint; var sik; by cut;
            output out=_step4 sum=sk;
            /* proc print n; */
            data _step4; set _step4;
            title3 "step4";
            abs_sk=abs(sk);
            dummy=1;
            proc sort; by dummy;
            /* proc print n; */
            proc univariate noprint; var abs_sk;
            output out=_step5 max=maxsk;
            *********************************************;
            * compute S**2 *;
            *********************************************;
            *********************************************;
            * figure out dim(&time where &stat=1) *;
            *********************************************;
            data _diftm; set _dds&i;
            title3 "step5";
            if &censor=1;
            proc sort; by &event;

            data _ttt; set _diftm; by &event;
            if first.&event;

            proc univariate noprint; var &event;
            output out=_deathtim N=nodeath;

               /* proc print n; */
            data _deathtim; set _deathtim;
            dummy=1;
            keep nodeath dummy;
            data _square; set _deathtim;
            do i=1 to nodeath;
            do j=1 to i;
            frac=1/(nodeath-j+1);
            keep i j frac;
            output;
            end;
            end;
            /* proc print n; */
            proc sort; by i j;
            data _sums2; set _square; by i j;
            if first.i then do;
            sumi=0; retain sumi;
            end;
            sumi=sumi+frac;
            if last.i then do;
            sumi=(1-sumi)*(1-sumi);
            output;
            end;
            /* proc print n; */
            proc univariate noprint; var sumi;
            output out=_d2 sum=s2 n=n2;
            data _step5; set _d2;
            ssquare=(1.0/(n2-1))*s2;
            dummy=1;
            proc sort; 
               by dummy;
            RUN;

            data _step5; 
               merge _step4 _step5; by dummy;
               title3 "step 5";
               q=abs_sk/(sqrt(ssquare)*sqrt(n2-1));
               if q>1 then p=2*exp(-2*q*q);
               if q<=1 then p=-99;
               format p pf.;
            RUN;

            proc sort DATA = _step5; 
               by cut;
            RUN;

            DATA _cut; 
               set _cut;
               &var=somecut;
               keep cut &var;
               /* Keep label from original data */
               LABEL &var = "%TRIM(%BQUOTE(&vv))";
            RUN;

            proc sort DATA = _cut; 
               by cut;
            RUN;

            data _temp; merge _step5 _cut; by cut;
            if _N_=1 then do;
            maxsk=0; maxcut=0; retain maxsk maxcut;
            end;
            if abs_sk>maxsk then do;
            maxsk=abs_sk;
            maxcut=cut;
            end;
            dummy=1;
            output;

            /* proc print n; */
            proc sort; by dummy;

            data _temp; set _temp; by dummy;
            cut=maxcut;
            if last.dummy then output;
            keep cut;

            proc sort; by cut;

            TITLE3 "Final Result";
            data _final&i; 
            merge _step5 _cut _temp(in=in1);
            by cut;
            pick=" ";
            if in1 then Pick="<";
            label cut="Cut point"
            sk="sk"
            abs_sk="ABS(sk)"
            Q="Q statistics"
            p="P-value";
            RUN;

            /*********************************************************************************/
            /* Output to RTF */
            /* If running with and without outliers add to subtitle */
            %if &outlier = B AND &i = 1 %then %do;
               %let tit2 = Outcome: &evelab - Using All Points;
            %end;
            %else %if &outlier = B AND &i = 2 %then %do;
               %let tit2 = Outcome: &evelab - Excluding Outliers;
            %end;
            %else %if &outlier = T OR &outlier = F %then %do;
               %let tit2 = Outcome: &evelab;
            %end;

            %if &plot = T %then %do;

               /* If drawing vertical reference lines */
               %if &VREF = T %then %do;
                  /* Calculate median value */
                  PROC UNIVARIATE DATA = _dds&i noprint;
                     var &var;
                     output out=_med MEDIAN=median;
                  RUN; 

                  PROC SQL noprint;
                     /* Get median value */
                     select median into :median
                     from _med;
                     /* Get optimal cut point value */
                     select &var into :cutpoint
                     from _final&i
                     where Pick="<";
                  QUIT;
               %end;

               /* Based on reference lines to draw set color, thickness, and label */

               %if &VREF = T %then %do;
                  %let color = blue blue;
                  %let thick = 1 1;
                  %let rlabel = 'Median' 'Optimal Cut Point';
                  %let pattern = 2 4;
                  %do m = 1 %to &chosen_cnt; 
                     %let color = &color red;
                     %let thick = &thick 2;
                     %let rlabel = &rlabel 'Chosen Cut Point';
                     %let pattern = &pattern 1;
                  %end;
               %end;
               %else %if &VREF ~= T AND &chosen ~= %STR() %then %do;
                  %let color = red;
                  %let thick = 2;
                  %let rlabel = 'Chosen Cut Point';
                  %let pattern = 1;
                  %do m = 2 %to &chosen_cnt;
                     %let color = &color red;
                     %let thick = &thick 2;
                     %let rlabel = &rlabel 'Chosen Cut Point';
                     %let pattern = &pattern 1;
                  %end;
               %end;

               /* Martingale Residual Plot */
               %MPLOT(ds = _dds&i,
                    time = &event, 
                    stat = &censor,
                    cutvar = &var, 
                    subtit2 = %BQUOTE(&tit2),
                    subtit = %BQUOTE(&subtit),
                    HREF = &HREF,
                    /* Reference lines */
                    %if &VREF = T OR &chosen ~= %STR() %then %do;
                       VREF = &median &cutpoint &chosen,
                       VREFLABEL = &rlabel,
                       VREFCOLOR = &color,
                       VREFTHICK = &thick,
                       VREFPATT = &pattern,
                    %end;
                    DEC = &dec,
                    doc=F);

            /* End PLOT=T */     
            %end;

            %if %sysevalf(%superq(subtit)~=,boolean) %then %do;
               TITLE1 "Optimal Cut Point Search - %TRIM(&subtit)";
            %end;
            %else %do;
               TITLE1 "Optimal Cut Point Search";
            %end;
            %if %sysevalf(%superq(tit2)~=,boolean) %then %do;
               TITLE2 "&tit2";
            %end; 
            proc print data=_final&i label noobs;
               /* If minimal output is requested then only print optimal cutpoint */
               %if &minout = T %then %do;
                  where Pick="<";
               %end;
               var cut &var sk abs_sk q p pick; 
            run; 
            TITLE;

         %end;

         /* Save cutpoint */
         DATA _cutpoint&k&j;
            set %do i = 1 %to &num; _final&i (in=f&i) %end;;
            length event var $32 eventL $256;
            where Pick="<";
            /* If using both with and without outliers add variable to indicate which */
            %if &outlier = B %then %do;
               length outlier $18;
               if f1 = 1 then outlier = 'All Points';
               else outlier = 'Excluding Outliers';
            %end;
            %else %if &outlier = T %then %do;
               outlier = 'Excluding Outliers';
            %end;
            %else %if &outlier = F %then %do; 
               outlier = 'All Points';
            %end;
            event = "&event";
            eventL = "&eveLab";
            var = VNAME(&var);
            VarL = VLABEL(&var);
            cut_value = &var;

            drop &var;
         RUN;

      %end;
   %end;

   %if &doc = T %then %do;
      ODS RTF CLOSE;
   %end;

   /* Merge all cut points */
   DATA &out;
      set %do i = 1 %to &var_cnt; %do k = 1 %to &out_cnt; _cutpoint&i&k %end; %end;;
   RUN;

   /* Reload original options that were in use before running the macro */
   PROC OPTLOAD data=_options;
   RUN;

   /* Only delete files if not in debug mode */
   %if &debug ~= T %then %do;
      /* DELETE ALL TEMPORARY DATASETS that were created */
      proc datasets lib=work memtype=data noprint;  
         delete _options _dds: _dummy _comb _cut _d2 _deathtim _difage _diftm _dim _double 
               _final: _square %if &plot = T and &vref = T %then %do; _med %end;
               _step1 - _step5 _sums2 _temp _ttt 
               %do i = 1 %to &var_cnt; %do k = 1 %to &out_cnt; _cutpoint&i&k %end; %end;;
      quit; 
   %end;

%mend;

