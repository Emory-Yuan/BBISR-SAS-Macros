/**********************************************************************************************
***********************************************************************************************
Macro Name: MPLOT
Created Date/Author: Feb. 2012 /Yuan Liu
Last Update Date/Person: June 7, 2014/Dana Nickleach
Current Version: V13
Working Environment: SAS 9.3 English version

Contact: Dr. Yuan Liu yliu31@emory.edu 

Purpose:  To produce martingale residual plot for a continuous variable to predict time to 
event outcome. This plot will help to detect any nonlinear relationship, as well as to choose 
the meaningful cut point to categorize a numerical variable.  

Notes:  The graph will be saved in a Word document.

Parameters: 

DS             The name of the data set to be analyzed.
TIME           Time to event outcome variable.                                                      
STAT           The variable name for censoring indicator.  It is required that 1 is used for the 
               event and 0 for censored cases.
CUTVAR         The numerical variable to be plotted.  The variable name cannot be more than 26 
               characters.  
SM             Smoothing parameter. The default value is 0.65.
TICK           X-axis major tick mark unit (optional).  Note that if the unit is very small, a 
               larger unit may be chosen that is a multiple of the supplied number to prevent
               too many tick marks.   
HREF           Set to T for a horizontal reference line at zero (optional).  The default value is 
               F. 
VREF           Values where vertical reference lines should appear separated by spaces (optional).
               This macro is not equipped to handle more than 12 reference lines.
VREFLABEL      Labels to use for vertical reference lines (optional).  Each label should be in 
               quotes and multiple labels should be separated by a space.  The order of the labels
               needs to correspond to the VREF parameter.
VREFCOLOR      Colors to use for each vertical reference lines separated by a space (optional).  
               The order of the colors needs to correspond to the VREF parameter.  By default all
               lines will be in blue.  
VREFTHICK      Line thickness in pixels to use for each vertical reference lines separated by a 
               space (optional).  The order of the thicknesses needs to correspond to the VREF 
               parameter.  
VREFFPATT      Line pattern to use for each vertical reference lines separated by a space (optional).
               The order of the patterns needs to correspond to the VREF parameter.  By default a 
               different line pattern will be used for each reference line.
SUBTIT         Text to be added to the title to appear after, but on the same line as the 
               first title (optional).  Text should not be enclosed in quotes.  By default no 
               additional text will be printed.
SUBTIT2        Text to be used for a subtitle to appear under the main title (optional).  Text 
               should not be enclosed in quotes.  By default no subtitle will be printed.  
DEC            Fixed number of decimal places to display in the percentile table (optional).  Note
               that the maximum format width is 7.  So the number of decimal places can't be more 
               than 5, but depending on the size of the numbers may be less.  The best format is 
               used by default.
DOC            Set to T to create a RTF file containing the output or F to suppress creation of 
               the RTF file (optional).  The default value is T.
OUTPATH        File path for output table to be stored.  
FNAME          File name for output table.  
DEBUG          Set to T if running in debug mode (optional).  Work datasets will not be deleted in
               debug mode.  This is useful if you are editing the code or want to further      
               manipulate the resulting data sets.  The default value is F.

***********************************************************************************************
For more details, please see the related documentation
**********************************************************************************************/

%macro MPLOT(ds=, time=, stat=, cutvar=, sm=0.65, tick=, HREF=F, VREF=, VREFLABEL=, VREFCOLOR=,
     VREFTHICK=, VREFPATT=, DEC=, OUTPATH=, FNAME=, SUBTIT=, SUBTIT2=, DOC=T, DEBUG=F);

   /* Make sure that these variables are local */
   %local min step max __Macro_Err ref_cnt i pattern color vref;

   /* Capitalize */
   %let debug = %UPCASE(&debug);
   %let href = %UPCASE(&href);
   %let doc = %UPCASE(&doc);

   /* Initialize error variable */
   %let __Macro_Err = 0;

   /* Check for variable names that are too long */
   %IF %LENGTH(&cutvar) > 26 %then %do;
         %put ERROR: Variable name &cutvar must be less than 27 characters.; 
           %let __Macro_Err=1;
   %END;

   /* If there is an error in the parameters supplied then exit */
   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   /* Save current options */
   PROC OPTSAVE out=_optionsMP;
   RUN;

   /* Count number of vertical reference lines to draw */
   %let ref_cnt = %sysfunc(countw(&VREF,' '));
   /* Line styles */
   %if &VREFPATT = %STR() %then %do;
      /* No line pattern supplied - use default */
      %let pattern = 2 4 5 8 14 15 20 26 34 35 41 42;
   %end;
   %else %let pattern = &VREFPATT;

   data dds;
      set &ds;
       /* Drop those missing outcome */
       if &time = . OR &stat = . then delete;
       /* Drop those missing cutvar */
      where &cutvar ne .;
   run;

   *******************************************************;
   * Plot Martingale Residual plot                        ;
   *******************************************************;

   PROC PHREG DATA=dds noprint;
      MODEL &time*&stat(0) = ;
      OUTPUT OUT=residuals RESMART=resmart; 
   run;

   PROC SORT DATA=residuals;
      BY &time;
   RUN;

   proc means data= dds noprint;
      var &cutvar;
       output out=mm (drop=_TYPE_ _FREQ_) min= max=/autoname; 
   run;

   data mm;
      set mm;

       /* Set min to nearest given tick unit below min */
       %if &tick ~= %STR() %then %do;
          min = FLOOR(&cutvar._Min/&tick)*&tick;
       %end;
       %else %do;
          min = &cutvar._Min;
      %end;

      del = (&cutvar._Max-min)/10; 
       /* Round increments to nearest unit given */
       %if &tick ~= %STR() %then %do;
          /* Make sure it is greater than 0 */
          del = MAX(ROUND(del,&tick),&tick);
      %end;

       /* Adjust max to nearest given tick unit above max */
       %if &tick ~= %STR() %then %do;
         max = CEIL(&cutvar._max/del)*del;
      %end;
       %else %do;
          max = &cutvar._max;
      %end;
   run;

   proc sql noprint;
      select Min into: min separated by " " from mm;
      select Max into: max separated by " " from mm;
      select del into: step separated by " " from mm;;
   quit;
   
   ODS NOPROCTITLE;
   %if &doc = T %then %do;
      OPTIONS NODATE;
      ODS RTF BODYTITLE startpage = NEVER FILE= "&OUTPATH.&FNAME &SYSDATE..DOC"; 
   %end;

   /* Calculate percentiles */
   proc univariate data=dds noprint; 
      var &cutvar;
      output out = dd n=N pctlpre=P_ pctlpts= 0 to 100 by 10;
   run;

   /* Percentiles */
   DATA dd;
      set dd;
       LABEL N='N'
          p_0 = 'Min'
          p_10 = '10th'
          p_20 = '20th'
          p_30 = '30th'
          p_40 = '40th'
          p_50 = '50th'
          p_60 = '60th'
          p_70 = '70th'
          p_80 = '80th'
          p_90 = '90th'
          p_100 = 'Max';
   RUN;

   /* Transpose so that each obersvation is a percentile */
   PROC TRANSPOSE DATA = dd out=tran;
   RUN;

   /* Annotation containing percentile labels */
   DATA anno1;
      set tran;
       length label $7;
       X1SPACE = 'GRAPHPERCENT';
      Y1SPACE = 'GRAPHPERCENT';
       label = _LABEL_;
       label = RIGHT(label);
       function = 'text';
      anchor = 'BOTTOMRIGHT';
       y1 = 5;
       x1 = 8.25*_n_;
   RUN;
   
   /* Annotation containing percentiles */
   DATA anno2;
      set tran;
       length label $7;
       X1SPACE = 'GRAPHPERCENT';
      Y1SPACE = 'GRAPHPERCENT';
       /* Show number of decimal places specified by user */
       %if &dec ~= %STR() %then %do;
          /* Don't add decimal places to the N */
          if _NAME_ = 'N' then label = PUT(col1,best7.);
          else label = PUT(col1,7.&dec);
       %end;
       /* Otherwise use best format */
       %else %do;
          label = PUT(col1,best7.);
       %end;
       label = RIGHT(label);
       function = 'text';
      anchor = 'BOTTOMRIGHT';
       y1 = 1;
       x1 = 8.25*_n_;
   RUN;

   /* Combine two annotation data sets */
   DATA anno;
      set anno1 anno2;
   RUN;

   GOPTIONS RESET=all;   

   %if %sysevalf(%superq(subtit)~=,boolean) %then %do;
      TITLE1 "Martingale Residual Plot - %TRIM(&subtit)";
   %end;
   %else %do;
      TITLE1 "Martingale Residual Plot";
   %end;
   %if %sysevalf(%superq(subtit2)~=,boolean) %then %do;
      TITLE2 "&subtit2";
   %end; 

   PROC SGPLOT DATA = residuals sganno=anno PAD=(BOTTOM=10%);
      loess y=resmart x=&cutvar/smooth=&sm MARKERATTRS=(symbol=asterisk) 
          LINEATTRS=(COLOR=blue) NOLEGFIT;
      %if &href = T %then %do; 
         REFLINE 0/LINEATTRS=(COLOR=blue);
       %end;
       /* Vertical reference lines */
       %do i = 1 %to &ref_cnt;
          /* Line color - use blue by default */
          %if %sysevalf(%superq(vrefcolor)~=,boolean) %then %let color = %SCAN(&vrefcolor,&i,%STR( ));
         %else %let color = blue;

           REFLINE %sysfunc(scan(&vref,&i,' '))/AXIS=x TRANSPARENCY=.5
               LINEATTRS=(COLOR=&color PATTERN=%sysfunc(scan(&pattern,&i,' '))
               /* Set line thickness if supplied */
               %if %sysevalf(%superq(vrefthick)~=,boolean) %then %do;
                  THICKNESS=%SCAN(&vrefthick,&i,%STR( ))
            %end;)
              %if %sysevalf(%superq(vreflabel)~=,boolean) %then %do;
                  NAME=%SCAN(&vreflabel,&i,%STR( ),Q) 
                  LEGENDLABEL=%SCAN(&vreflabel,&i,%STR( ),Q)
              %end;;
       %end;
       /* Slightly increased axis max to avoid points not displaying */
      XAXIS VALUES=(&min to %SYSEVALF(&max+0.0001) by &step);
       FOOTNOTE 'Percentiles';
   RUN;
   TITLE;
   FOOTNOTE;

   %if &doc = T %then %do;
      ods rtf close;
   %end;

   /* Did not change options if did not created RTF file */
   /* Note if calling this macro from FINDCUT then running this proc is not neccessary and */
   /* causes an insertion of a section break in the RTF document */
   %if &doc = T %then %do;
      /* Reload original options that were in use before running the macro */
      PROC OPTLOAD data=_optionsMP;
      RUN;
   %end;

   /* Only delete files if not in debug mode */
   %if &debug ~= T %then %do;

      /* Delete temporary datasets that were created */
      proc datasets lib=work memtype=data noprint;  
          delete mm residuals dds dd _optionsMP anno1 anno2 anno tran;
      quit; 

   %end;

%mend;

