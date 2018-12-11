/**********************************************************************************************
** Citation: This Macro was created by the Biostatistics and Bioinformatics Share Resource in 
** Winship Cancer Center at Emory University. Free Download at https://bbisr.winship.emory.edu/.
***********************************************************************************************

Macro Name: KM_PLOT
Created Date/Author/Contact: Oct. 5, 2012/Dana Nickleach
Last Update Date/Person/Contact: March, 2017/Yuan Liu
Current Version: V24
Working Environment: SAS 9.4 English version

Purpose:  To produce one or more Kaplan-Meier plots along with summary statistic tables.  
Multiple plots can be created using different data sets, outcomes, and strata variables with 
one macro call.  A plot will be produced for each combination of input parameters.    

Notes:  The macro "PLSurvivalTemp V5.sas" is required.     

Parameters: 

DSN            The name of the data sets to be analyzed separated by spaces.
EVENTS         List of time to event outcome variables separated by spaces.
CENSORS        List of censoring indicator variables separated by spaces.  The order must 
               correspond to the order of the EVENTS list.  Values of 0 indicate censored 
               observations.
GRPLIST        List of variables to use as strata separated by spaces.  A separate plot will be 
               produced for each strata variable.  If no variables are specified then survival 
               curves for the whole sample will be produced.
STRATA2        A second variable to stratify on in conjunction with the first.  The plots will 
               be run using two strata variables in one plot.  Note that if you specified 
               multiple variables in GRPLIST then the analysis will be run once for each 
               variable in GRPLIST stratified on both that variable and the variable in STRATA2 
               every time.  This should be empty if GRPLIST is empty.
TITLE          Title for the plot in quotes or set to DSNLABEL (no quotes) in order to title 
               each plot with the label of the data set used for the plot.  If there is no data
               set label available then the data set name will be used. For subgroup analysis, one
			   can create seperate dataset with each having a dataset lable and use DSNLABEL option
			   so that each subgroup plot will have corresponding title.
ENTRYTITLE     Specification for entrytitle (plot title).  The actual text of the title 
               should be in quotes, but other text can be included as well.  For example,
               entrytitle = HALIGN=LEFT "A" will left align the title "A".  The SAS default
               title will be used if this is not specified.  Note that this differs from the 
               TITLE parameter in that the ENTRYTITLE appears on the plot space, whereas TITLE
               appears above the plot space. Using entrytitle = "" to suppress the title.
SUBTITLE	   The second title after entrytitle. When ATRISK = T, subtitle will appear. to suppress the title
			   use subtitle = "".	
XLAB           Customize the label for x-axis. Defult value will be used when setting XLAB =,
YLAB		   Customize the label for y-axis. Defult value will be used when setting YLAB =,
WEIGHT         The name of the weight variable to be used in the WEIGHT statement. The weight should be 
			   appropriately normalized before calling this macro.  Leave this empty if not using weights.    
PAIRWISE       Set to T to request a table of pairwise comparisons.  The default value is F.  
TABLE          Set to T to request a table with N, event, censor, median survival and its 95% 
               CI.  The default is F.
TIMELIST       List of time points separated by spaces at which survival estimates and 95% CI 
               are reported in a table.  Note that this is currently only set up to work with 
               one strata variable.
UNIT           Unit of time for survival as it appears in the data, i.e. Yr, Mo, etc. to be 
               used in the table header for time point estimates.  This is only necessary if 
               TIMELIST was specified. 
JOIN           Set to T to join median survival estimates and survival point estimates into 
               one table.  This is F by default.  Note that this is not set up to work when 
               STRATA2 is specified.
PLOT           Set to F to suppress the KM plot.  The default is T.
XMAX           Maximum x-axis value to display.  By default the maximum time available in the 
               data will be used.
XTICK		   Specify x-axis tick values to display. By default the automatic even space tick 
			   values will be used. the maximum tick value should be greater than the value of 
			   XMAX. Example of useage: XTICK = (0 6 12 18 24 30 36).
YMIM 		   Specify minimum value of y-axis. The default is 0.
ATRISK		   Set up to T to generate number of subject at risk at x-axis tick values or defult tick values.
			   The defult is F.
DOC            Set to F to suppress creation of the RTF file.  The default value is T.
OUTPATH        File path for output to be stored.
FNAME          File name for output document.
STYLE          The name of the style template to use for the document.
TOC            Label for outer item in table of contents in quotes.  If an argument is 
               provided then a table of contents will be generated.  By default a table of 
               contents is not generated.  This will only have an effect when DOC=T or this 
               macro is called within an RTF statement using the CONTENTS=YES TOC_DATA options.
               Once you open the resulting document in Word, hit CTRL+A and then F9 to update
               the TOC.                 
NONCENSORED    Remove the censore indicator from the KM plot. The default value is F.
TESTNM		   Quoted text to specify the name of test used for p-value.
TESTP          The p-value from TESTNM. TESTNM and TESTP could be useful when dealing 
			   with weighted or matched data. 
FAILURE		   Generate failure curve. Under development. Default =F.
DEBUG          Set to T to run in debug mode.  Work datasets will not be deleted in debug 
               mode.  This is useful if you are editing the code or want to further manipulate
               the resulting data sets.  The default value is F.

**********************************************************************************************/

/* Create Kaplan-Meier Survival Curve */
%macro km_plot(dsn=,events=,censors=,grplist=,strata2=,title=,entrytitle="Kaplan-Meier Plot", subtitle=SECONDTITLE,xlab=, ylab=,pairwise=F,
weight=, timelist=,unit=,join=F,plot=T,table=F,xMax=MAXTIME, xTick = XTICKVALS,yMin = 0, atrisk=F, DOC=T,OUTPATH=,FNAME=,STYLE=,toc=,
DEBUG=F, noncensored=F, failure=F, testnm=, testp=  );

   %local STRATA2 subgroups DSN TIMELIST I XMAX J PLOT CENSOR TITLE EVENT GRP TIME_CNT GRPLIST STR_LAB
     CENSORS EST STRATA2TYP GRP_CNT EVENTS TIMELAB STRATA1TYP PAIRWISE TABLE EVENT_CNT DSN_CNT K
     DATA __MACRO_ERR RTF;

   /* Capitalize */
   %let pairwise = %UPCASE(&pairwise);
   %let plot = %UPCASE(&plot);
   %let table = %UPCASE(&table);
   %let doc = %UPCASE(&doc);
   %let join = %UPCASE(&join);
   %let debug = %UPCASE(&debug);
   %let noncensored=%UPCASE(&noncensored);

   /* Initialize - keep track of number of items printed to output */
   %let outNum = 0;

   /* Initialize error flag */
   %let __Macro_Err = 0;

   /* STRATA2 should be blank if GRPLIST is blank */
   %if &grplist = %STR() and &strata2 ~= %STR() %then %do;
      %put ERROR: STRATA2 should be empty if GRPLIST is empty.;
      %let __Macro_Err = 1;
   %end;

   /* If TIMELIST is specified and JOIN=T then UNIT is required */
   %if &timelist ~= %STR() and &join = T and &table = T and &unit = %STR() %then %do;
      %put ERROR: UNIT needs to be specified when using TIMELIST and JOIN=T.;
      %let __Macro_Err = 1;
   %end;

   /* Cannot specify JOIN=T with STRATA2 */
   %if &join = T and &strata2 ~= %STR() %then %do;
      %put ERROR: JOIN cannot = T if STRATA2 is specified.;
      %let __Macro_Err = 1;
   %end;

   %if &__Macro_Err. %then %do;
      data _null_;
         abort 3;
      run;
   %end;

   ODS NOPROCTITLE;

   /* Count number of data sets */
   %let dsn_cnt = %sysfunc(countw(&dsn,' '));
   /* Count number of events */
   %let event_cnt = %sysfunc(countw(&events));
   /* Count number of time estimates requested */
   %let time_cnt = %sysfunc(countw(&timelist,' '));

   /* Count number of strata variables */
   %if &grplist = %STR() %then %let grp_cnt = 1;
   %else %do;
      %let grp_cnt = %sysfunc(countw(&grplist));
   %end;   

   /* See if RTF destination is open */
   PROC SQL noprint;
      select count(*) into :RTF
      from sashelp.vdest
      where destination = 'RTF';
   QUIT;

   %if &doc = T %then %do;
      ODS RTF BODYTITLE FILE="&outpath.&fname..DOC" 
         %if &style ~= %STR() %then %do; STYLE=&style %end; 
         %if &toc ~= %STR() %then %do; CONTENTS=YES toc_data %end;; 
   %end;

   /* Repeat for each data set */
   %do k = 1 %to &dsn_cnt;
      %let data = %SCAN(&dsn,&k,%STR( ));

      /* Use data set label/name if option is specified */
      %if %QUOTE(%UPCASE(&title)) = DSNLABEL %then %do;
         /* Get data set label */
         data _null_ ;
            /* open dataset, and keep the ID */
            dsid=open("&data") ; 
            /* get the numeric attribute NOBS - number of observations */
            label=attrc(dsid,'label') ;     
            /* If the data set is not label use the data set name */
            if label = ' ' then label = "&data";
            /* write the value out to a macro variable */
            call symputx("datlab&k",label) ; 
            /* close dataset using the ID */ 
            dsid=close(dsid) ;                
         run;

         %let title_R = "&&datlab&k";
      %end;

      /* Otherwise use provided title */
      %else %let title_R = &title;

      /* Get list of variables in the data set */
      PROC CONTENTS DATA = &data out=_cont noprint;
      RUN;

      %do j = 1 %to &event_cnt;
         %let event = %SCAN(&events,&j);
         %let censor = %SCAN(&censors,&j);

         /* Check to see if outcome is available in data set */
         PROC SQL noprint;
            select count(*) into :check
            from _cont
            where UPCASE(NAME) = "%UPCASE(&event)";
         QUIT;

         %if &check = 0 %then %do;
            %put WARNING: &event was not found in the data set &data.;
         %end;
         %else %do;

            /* Get outcome label */
            data _NULL_; 
               set &data (obs=1); 
               call symput('timelab',label(&event));
            run;

            %if &grp_cnt > 0 %then %do;
               %do i = 1 %to &grp_cnt;
                  %let grp = %SCAN(&grplist,&i); 

                  %if &grp ~= %STR() %then %do;
                     /* Get label of strata variable */
                     DATA _NULL_;
                        set &data (obs=1);
                        CALL SYMPUT('str_lab',VLABEL(&grp));
                     RUN;
                  %end;

				/*When the input for entrytitle, subtitle,xMax,xTick are empty, then set to the default value*/
				%if &entrytitle = %STR() %then %do; %let entrytitle="Kaplan-Meier Plot"; %end;
				%if &subtitle = %STR() %then %do; %let subtitle=SECONDTITLE; %end;
				%if &xMax = %STR() %then %do; %let xMax=MAXTIME; %end;
				%if &xTick = %STR() %then %do; %let xTick=XTICKVALS; %end;
				%if &testnm = %STR() %then %do; %let testnm = TestName; %end;
				%if &testp = %STR() %then %do; %let testp = pValue; %end;



				  /* If only 1 strata variable then modify legend label */
                  %if &grp ~= %STR() AND &strata2 = %STR() %then %do;
                     /* Use label for legend variable */
                  %PLSurvivalTemp(legendvar="&str_lab",xViewmax=&xMax,yViewmin = &ymin, xTickv =&xTick,entrytitle=&entrytitle, subtitle = &subtitle,
									x_label = &xlab, y_label=&ylab, testnm=&testnm, testp=&testp);
                  %end;
                  %else %do;

                  /* Change entrytitle */
                  %PLSurvivalTemp(xViewmax=&xMax,xTickv =&xTick,yViewmin = &ymin, entrytitle=&entrytitle, subtitle = &subtitle, x_label = &xlab, y_label=&ylab,
					testnm=&testnm, testp=&testp);
                  %end;

                  /* Use ODS RTF commands if RTF destination is open */
                  /* This will allow for page breaks in appropriate places when this macro
                  call is placed within ODS RTF statements */
                  %if &doc = T OR &RTF = 1 %then %do;
                     ODS RTF STARTPAGE=YES;
                  %end;
                  /* Kaplan - Meier Survival Curves */
                  %if &plot = T %then %do; 
                     TITLE &title_r;
                     %if &failure ~=T %then %do; ODS SELECT SurvivalPlot %end;
					 %else %do; ODS SELECT FailurePlot %end; ;
                     /* Increment output number */
                     %let outnum = %EVAL(&outnum+1);
                  %end;
                  %else %do;
                     ODS SELECT NONE;
                  %end;
                  %if &pairwise = T %then %do;
                     ODS OUTPUT SurvDiff=_pval;
                  %end;
                  %if &pairwise = T and &strata2 ~= %STR() %then %do;
                     ODS OUTPUT legend=_legend;
                  %end;
                  %if &time_cnt > 0 %then %do; 
                     ODS OUTPUT ProductLimitEstimates=_est; 
                  %end;
                  %if &table = T %then %do;
                     ODS OUTPUT Quartiles = _quart CensoredSummary=_censor;
                  %end;
                  /* Heading 1 for Table of contents */
                  %if &toc ~= %STR() %then %do;
                     ODS PROCLABEL=&toc;
                  %end;
					/*Sort data if Subgropus analysis is needed
                  %if &subgroups ~= %STR() %then %do;
				  proc sort data=&data; by &subgroups; run;%end;*/
				  
                  proc lifetest data=&data outsurv=_estci plots=survival(test %if &atrisk = T AND &xtick ~= XTICKVALS  %then %do; atrisk = &xtick %end;
																				%else %if &atrisk = T %then %do; atrisk(outside(0.15) %end;
                                                                               %if &noncensored = T %then %do; nocensor %end;
																				%if &failure = T %then %do; failure %end;)
                         %if &time_cnt > 0 %then %do; timelist=(&timelist) %end;;
                     time &event*&censor(0);
                     %if &grp ~= %STR() %then %do;
                        strata &grp &strata2 %if &pairwise = T %then %do; /diff %end;;
                     %end;
                     %if &weight ~= %STR() %then %do;
                        WEIGHT &weight;
                     %end;
                     /*%if &subgroups ~= %STR() %then %do;
                        BY &subgroups;
						where &subgroups ~=. ;
                     %end;   */               run;
                  ODS SELECT ALL;

                  %if &plot ~= T %then %do;
                     TITLE &title_r;
                  %end;

                  /* Add Table with N, event, censor, median survival and its 95% CI. */
                  %if &table = T %then %do;
               
                     %if &grp ~= %STR() %then %do;
                        proc sort data=_quart; 
                           by &grp &strata2; 
                        run;
                        proc sort data=_censor; 
                           by &grp &strata2; 
                        run;
                     %end;

                     DATA _summ;
                        merge _quart (where=(percent=50)) _censor %if &grp ~= %STR() %then %do; (where=(MISSING(&grp)=0)) %end;;
                        %if &grp ~= %STR() %then %do; by &grp &strata2; %end;

                        per_event = failed/total;
                        per_cens = censored/total;

                        /* N & % */
                        event = PUT(Failed,8.) || ' (' || TRIM(LEFT(PUT(per_event,percent6.0))) || ')';
                        censor = PUT(censored,8.) || ' (' || TRIM(LEFT(PUT(per_cens,percent6.0))) || ')';

                        /* Median survival & 95% CI */
                        length UCL $9;
                        if UpperLimit = . then UCL = 'NA';
                        else UCL = LEFT(PUT(ROUND(UpperLimit,.1),best9.));
    
                        length LCL $9;
                        if LowerLimit = . then LCL = 'NA';
                        else LCL = LEFT(PUT(ROUND(LowerLimit,.1),best9.));

                        length est $9;
                        if estimate = . then est = 'NA';
                        else est = LEFT(PUT(ROUND(estimate,.1),best9.));

                        surv = TRIM(est) || " (" || TRIM(LCL) || ", " || TRIM(UCL) || ")";

                        /* This variable is needed for proc report */
                        dummy = 1;
                     RUN;

                     %if &join = F %then %do;
                        %if &doc = T OR &rtf = 1 %then %do;
                           ods rtf startpage=no; 
                        %end;
                        %if &plot = T %then %do; 
                           TITLE; 
                        %end;
                       
                        /* Heading 1 for Table of contents */
                        %if &toc ~= %STR() and &outnum = 0 %then %do;
                           ODS PROCLABEL=&toc;
                        %end;
                        %else %do;
                           ODS PROCLABEL=' ';
                        %end;
                        /* Increment output number */
                        %let outnum = %EVAL(&outnum+1);

                        proc report nowd data=_summ style(report)={rules=GROUPS frame=hsides} 
                                   style(header)={BACKGROUNDCOLOR=none} LS=256 CONTENTS='Survival Estimates';
                           col dummy &grp &strata2 total event Censor surv;
                           %if &grp ~= %STR() %then %do; 
                              define &grp/order style={just=l}; 
                           %end;
                           define total/style={just=c};
                           define event/style={just=c};
                           define censor/style={just=c};
                           define surv/style={just=c};
                           /* Do this to remove 3rd node level from the table of contents */
                           DEFINE dummy/order noprint;
                           break before dummy/ contents="" page;
                           LABEL total = 'No. of Subject'
                                   event = 'Event'
                                   censor = 'Censored'
                                   surv = 'Median Survival (95% CI)';
                        run;
                     %end;
               
                  %end;

                  /* If time estimates were requested */
                  %if &time_cnt > 0 %then %do;

                     %if %UPCASE(&event) = SURVIVAL %then %do;
                        %let est = survival2;
                     %end;
                     %else %do;
                        %let est = survival;
                     %end;

                     /* Merge on CI */
                     PROC SORT DATA = _est;
                        by &grp &strata2 &event &est Censor;
                     RUN;
                     PROC SORT DATA = _estci (rename=(_CENSOR_=censor));
                        by &grp &strata2 &event &est Censor;
                     RUN;
                     DATA _est;
                        merge _est (in=a) _estci (keep=&grp &strata2 &event &est Censor SDF_LCL SDF_UCL);
                        by &grp &strata2 &event &est Censor;
                        if a;

						%if &failure ~= T %then %do;
                        /* Convert estimates to character */
                        length UCL $9;
                        if SDF_UCL = . then UCL = 'NA';
                        else UCL = LEFT(PUT(SDF_UCL,percent9.1));

                        length LCL $9;
                        if SDF_LCL = . then LCL = 'NA';
                        else LCL = LEFT(PUT(SDF_LCL,percent9.1));

                        length est $9;
                        if &est = . then est = 'NA';
                        else est = LEFT(PUT(&est,percent9.1));

                        /* Combine estimate and CI */
                        rate = TRIM(est) || " (" || TRIM(LCL) || ", " || TRIM(UCL) || ")";

                        label = PUT(Timelist,8.0) || " &unit" || ' Survival';


						%end;
						%else %do;
						length UCL $9;
                        if SDF_UCL = . then UCL = 'NA';
						else UCL = LEFT(PUT(1-SDF_UCL,percent9.1));

                        length LCL $9;
                        if SDF_LCL = . then LCL = 'NA';
                        else LCL = LEFT(PUT(1-SDF_LCL,percent9.1));

                        length est $9;
                        if &est = . then est = 'NA';
                        else est = LEFT(PUT(1-&est,percent9.1));

                        /* Combine estimate and CI */
                        rate = TRIM(est) || " (" || TRIM(LCL) || ", " || TRIM(UCL) || ")";

                        label = PUT(Timelist,8.0) || " &unit" || ' Failure';

                        %end;

                        /* This variable is needed for proc report */
                        dummy = 1;

                     RUN;

                     PROC SORT DATA = _est;
                        by &grp &strata2 Timelist &event &est Censor;
                     RUN;

                     /* If joining time estimates with survival estimates */
                     %if &table = T AND &join = T %then %do;
                            
                        /* Transpose to stratum level */
                        PROC TRANSPOSE DATA = _est prefix=time out=_tran (drop=_NAME_);
                           id label;
                           idlabel label;
                           by &grp;
                           var rate;
                        RUN;

                        /* Merge onto median estimates */
                        DATA _summ;
                           merge _summ _tran;
                           by &grp;

                           /* This variable is needed for proc report */
                           dummy = 1;
                        RUN;

                        %if &doc = T OR &rtf = 1 %then %do;
                           ods rtf startpage=no; 
                        %end;
                        %if &plot = T %then %do; 
                           TITLE; 
                        %end;
                        /* Heading 1 for Table of contents */
                        %if &toc ~= %STR() and &outnum = 0 %then %do;
                           ODS PROCLABEL=&toc;
                        %end;
                        %else %do;
                           ODS PROCLABEL=' ';
                        %end;
                        /* Increment output number */
                        %let outnum = %EVAL(&outnum+1);

                        proc report nowd data=_summ style(report)={rules=GROUPS frame=hsides} 
                              style(header)={BACKGROUNDCOLOR=none} LS=256 
                              style(column)={just=c} CONTENTS='Survival Estimates';
                           col dummy &grp &strata2 total event Censor surv time:;
                           %if &grp ~= %STR() %then %do; 
                              define &grp/order style={just=l}; 
                           %end;
                           define total/style={just=c};
                           define event/style={just=c};
                           define censor/style={just=c};
                           define surv/style={just=c};
                           /* Do this to remove 3rd node level from the table of contents */
                           DEFINE dummy/order noprint;
                           break before dummy/ contents="" page;

                           LABEL total = 'No. of Subject'
                                event = 'Event'
                                censor = 'Censored'
                                surv = 'Median Survival (95% CI)';
                        run;

                     %end;
                     %else %do;

                        %if &doc = T OR &rtf = 1 %then %do;
                           ODS RTF STARTPAGE=NEVER;
                        %end;
                        %if &plot = T OR &table = T %then %do; 
                           TITLE; 
                        %end;
                        /* Heading 1 for Table of contents */
                        %if &toc ~= %STR() and &outnum = 0 %then %do;
                           ODS PROCLABEL=&toc;
                        %end;
                        %else %do;
                           ODS PROCLABEL=' ';
                        %end;
                        /* Increment output number */
                        %let outnum = %EVAL(&outnum+1);

                        proc report nowd data=_est style(report)={rules=GROUPS frame=hsides} 
                                   style(header)={BACKGROUNDCOLOR=none} LS=256 CONTENTS='Survival Estimates';
                           col dummy &grp &strata2 Timelist rate;
                           %if &grp ~= %STR() %then %do;
                              define &grp/order style={just=l};
                           %end;
                           define timelist/style={just=c};
                           define rate/style={just=c};
                           /* Do this to remove 3rd node level from the table of contents */
                           DEFINE dummy/order noprint;
                           break before dummy/ contents="" page;
                           format timelist best12.;
                           LABEL timelist = "&timelab"
                            %if &failure ~=T %then %do; rate = 'Survival Rate(95% CI)' %end;
							%else %do;rate = 'Failure Rate(95% CI)' %end;;
                        RUN;

                     %end;

                  %end;

                  /* If pairwise differences were requested */
                  %if &pairwise = T %then %do;

                     /* Pairwise output will differ depending on whether 1 or 2 strata 
                     variables were specified */
                     %if &strata2 ~= %STR() %then %do;
                        /* Get variable types */
                        DATA _NULL_;
                           set &data;
                           CALL SYMPUT('strata1typ',VTYPE(&grp));
                           CALL SYMPUT('strata2typ',VTYPE(&strata2));
                        RUN; 

                        /* Merge on strata identifiers for first level of comparision */
                        PROC SORT DATA = _pval;
                           by StratumNumber1;
                           /* Only want logrank test */
                           where Test = 'Logrank';
                        RUN;
                        DATA _pval;
                           merge _pval (in=a) _legend (rename=(stratum=StratumNumber1) drop=Symbol);
                           by StratumNumber1;
                           if a;
                           rename &grp = &grp.1 &strata2 = &strata2.1;
                        RUN;

                        /* Merge on strata identifiers for second level of comparision */
                        PROC SORT DATA = _pval;
                           by StratumNumber2;
                        RUN; 
                        DATA _pval;
                           merge _pval (in=a) _legend (rename=(stratum=StratumNumber2) drop=Symbol);
                           by StratumNumber2;
                           if a;
                           rename &grp = &grp.2 &strata2 = &strata2.2;
                        RUN;
                     
                        DATA _pval;
                           set _pval;
                           length comp1 comp2 comp $256;
                           %if &strata1typ = N %then %do;
                              comp1 = LEFT(PUTN(&grp.1,VFORMAT(&grp.1)));
                              comp2 = LEFT(PUTN(&grp.2,VFORMAT(&grp.2)));
                           %end;
                           %else %do;
                              comp1 = LEFT(TRIM(PUTC(&grp.1,VFORMAT(&grp.1))));
                              comp2 = LEFT(TRIM(PUTC(&grp.2,VFORMAT(&grp.2))));
                           %end; 
                           %if &strata2typ = N %then %do;
                              comp1 = LEFT(TRIM(comp1)) || " " || LEFT(PUTN(&strata2.1,VFORMAT(&strata2.1)));
                              comp2 = LEFT(TRIM(comp2)) || " " || LEFT(PUTN(&strata2.2,VFORMAT(&strata2.2)));
                           %end;
                           %else %do;
                              comp1 = LEFT(TRIM(comp1)) || " " || LEFT(TRIM(PUTC(&strata2.1,VFORMAT(&strata2.1))));
                              comp2 = LEFT(TRIM(comp2)) || " " || LEFT(TRIM(PUTC(&strata2.2,VFORMAT(&strata2.2))));
                           %end;

                           comp = LEFT(TRIM(comp1)) || " vs. " || LEFT(TRIM(comp2));

                           /* This variable is needed for proc report */
                           dummy = 1;
                        RUN;

                     %end;
                     /* Only one strata variable */
                     %else %do;

                        /* Get variable types */
                        DATA _NULL_;
                           set _pval;
                           /* Note in this case strata1 and strata2 type will be the same b/c
                           there is only 1 strata varaible */
                           CALL SYMPUT('strata1typ',VTYPE(Stratum1));
                        RUN; 

                        DATA _pval;
                           set _pval;
                           length comp $256;
                           %if &strata1typ = N %then %do;
                              comp1 = LEFT(TRIM(PUTN(Stratum1,VFORMAT(Stratum1))));
                              comp2 = LEFT(TRIM(PUTN(Stratum2,VFORMAT(Stratum2))));
                           %end;
                           %else %do;
                              comp1 = LEFT(TRIM(PUTC(Stratum1,VFORMAT(Stratum1))));
                              comp2 = LEFT(TRIM(PUTC(Stratum2,VFORMAT(Stratum2))));
                           %end; 

                           comp = LEFT(TRIM(comp1)) || " vs. " || LEFT(TRIM(comp2));

                           /* Only want logrank test */
                           where Test = 'Logrank';

                           /* This variable is needed for proc report */
                           dummy = 1;
                        RUN;
                     %end;

                     %if &doc = T OR &rtf = 1 %then %do;
                        ODS RTF STARTPAGE=NEVER;
                     %end;
                     /* Don't repeat title */
                     %if &plot = T OR &table = T OR &time_cnt > 0 %then %do; 
                        TITLE; 
                     %end;
                     /* Heading 1 for Table of contents */
                     %if &toc ~= %STR() and &outnum = 0 %then %do;
                        ODS PROCLABEL=&toc;
                     %end;
                     %else %do;
                        ODS PROCLABEL=' ';
                     %end;
                     /* Increment output number */
                     %let outnum = %EVAL(&outnum+1);

                     proc report nowd data=_pval style(report)={rules=GROUPS frame=hsides} 
                           style(header)={BACKGROUNDCOLOR=none} LS=256 CONTENTS='Pairwise Comparisons';
                        col dummy comp raw;
                        define comp/display style={just=l};
                        define raw/display style={just=c};
                       
                        /* Do this to remove 3rd node level from the table of contents */
                        DEFINE dummy/order noprint;
                        break before dummy/ contents="" page;

                        LABEL comp = 'Pairwise Comparison'
                                raw='Logrank p-value';
                     run;

                  %end;
                  TITLE;

                  /* Reset output number when changing outcome, data set, or strata */
                  %let outnum = 0;
               %end;

               %if &grp ~= %STR() AND &strata2 = %STR() %then %do;
                  /* Restore template */
                  proc template;
                     delete Stat.Lifetest.Graphics.ProductLimitSurvival;
                  run; 
               %end;

               /* Reset output number when changing outcome, data set, or strata */
               %let outnum = 0;

               %if &debug = F %then %do;
                  /* Delete excess data sets */
                  proc datasets lib=work memtype=data noprint; 
                     delete _estci _cont %if &table = T %then %do; _censor _quart _summ %end;
                     %if &pairwise = T %then %do; _pval %end;
                     %if &pairwise = T and &strata2 ~= %STR() %then %do; _legend %end;
                     %if &time_cnt > 0 %then %do; _est %end;
                     %if &time_cnt > 0 and &join = T %then %do; _tran %end;; 
                  quit; 
               %end;

            %end;

         %end;
      %end;
   %end;

   %if &doc = T %then %do;
      ODS RTF CLOSE;
   %end;

%mend km_plot;
