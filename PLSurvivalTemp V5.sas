/**********************************************************************************************
***********************************************************************************************
Macro Name: PLSurvivalTemp
Created Date/Author/Contact: Aug. 3, 2012/Dana Nickleach/dnickle@emory.edu
Last Update Date/Person/Contact: Mar 15, 2017/ Yuan Liu
Current Version: V5

Working Environment: SAS 9.4 English version

Purpose:  To modify the product limit survival template for KM plots.   

Notes:   

Parameters: 

LEGENDVAR      Label for strata variable to appear in legend (optional).  By defualt the 
               strata variable name is used.
XVIEWMAX       Maximum x-axis value to display (optional).  By defualt the maximum time will
               be used.
YVIEWMIN	   Minimum y-axis value to desplay (optional). Defualt is 0.
ENTRYTITLE     Specification for entrytitle (plot title).  The actual text of the entrytitle 
               should be in quotes, but other text can be included as well.  For example,
               entrytitle = HALIGN=LEFT "A" will left align the title "A".  The SAS default
               title will be used if this is not specified.
XTICKV		   Specify x-aix tick value list.
SUBTITLE	   Specify second title for entrytitle.
               

Revision History:
Date           By        Version        Description
3/26/13        DCN       V2             Added XVIEWMAX parameter.
9/2/14         DCN       V3             Added ENTRYTITLE parameter.
10/22/14       YL		 V4				Added xtickv parameter, SUBTITLE parameter and x_label/y_label parameters.
2/6/15  	   YL		 V5				Added ymax to specify the max y axis to display.
3/17/2017	   YL        V5             Added testnm and testp to specify the the test used and plvalue.

**********************************************************************************************/

%macro PLSurvivalTemp(legendvar=groupname,xViewmax=MAXTIME,xTickv = XTICKVALS, yViewmin=0,entrytitle="Kaplan-Meier Plot", 
subtitle=SECONDTITLE,testnm=TestName, testp=pValue,
x_label=, y_label=);

   /* Template for PROC LIFETEST - Survival Plot */
   PROC TEMPLATE;
      define statgraph Stat.Lifetest.Graphics.ProductLimitSurvival;
      dynamic NStrata xName plotAtRisk plotCensored plotCL plotHW plotEP
          labelCL labelHW labelEP maxTime xtickVals xtickValFitPol method
          StratumID classAtRisk plotBand plotTest GroupName yMin
          Transparency SecondTitle TestName pValue;
      BeginGraph;
         if (NSTRATA=1)
            if (EXISTS(STRATUMID))
               entrytitle &entrytitle " for " STRATUMID;
            else
               entrytitle &entrytitle;
            endif;
          
            layout overlay / xaxisopts=(%if &x_label = %STR() %then %do; %end;%else %do; label=&x_label %end; offsetmin=.05
            linearopts=(viewmax=&xViewmax tickvaluelist=&xTickv
            tickvaluefitpolicy=XTICKVALFITPOL)) yaxisopts=(%if &y_label = %STR() %then %do; label = "Survival Probability" %end;%else %do; label=&y_label %end;
			shortlabel="Survival" linearopts=(
            viewmin=&yViewmin viewmax=1 tickvaluelist=(0 .2 .4 .6 .8 1.0)));
            if (PLOTHW=1 AND PLOTEP=0)
               bandplot LimitUpper=HW_UCL LimitLower=HW_LCL x=TIME /
                  modelname="Survival" fillattrs=GRAPHCONFIDENCE name=
                  "HW" legendlabel=LABELHW;
            endif;
            if (PLOTHW=0 AND PLOTEP=1)
               bandplot LimitUpper=EP_UCL LimitLower=EP_LCL x=TIME /
                  modelname="Survival" fillattrs=GRAPHCONFIDENCE name=
                  "EP" legendlabel=LABELEP;
            endif;
            if (PLOTHW=1 AND PLOTEP=1)
               bandplot LimitUpper=HW_UCL LimitLower=HW_LCL x=TIME /
                  modelname="Survival" fillattrs=GRAPHDATA1
                  datatransparency=.55 name="HW" legendlabel=LABELHW;
               bandplot LimitUpper=EP_UCL LimitLower=EP_LCL x=TIME /
                  modelname="Survival" fillattrs=GRAPHDATA2
                  datatransparency=.55 name="EP" legendlabel=LABELEP;
            endif;
            if (PLOTCL=1)
               if (PLOTHW=1 OR PLOTEP=1)
                  bandplot LimitUpper=SDF_UCL LimitLower=SDF_LCL x=
                     TIME / modelname="Survival" display=(outline)
                     outlineattrs=GRAPHPREDICTIONLIMITS name="CL"
                     legendlabel=LABELCL;
               else
                  bandplot LimitUpper=SDF_UCL LimitLower=SDF_LCL x=
                     TIME / modelname="Survival" fillattrs=
                     GRAPHCONFIDENCE name="CL" legendlabel=LABELCL;
               endif;
            endif;
            stepplot y=SURVIVAL x=TIME / name="Survival" rolename=(
               _tip1=ATRISK _tip2=EVENT) tip=(y x Time _tip1 _tip2)
               legendlabel="Survival";
            if (PLOTCENSORED=1)
               scatterplot y=CENSORED x=TIME / markerattrs=(symbol=
                  plus) name="Censored" legendlabel="Censored";
            endif;
            if (PLOTCL=1 OR PLOTHW=1 OR PLOTEP=1)
               discretelegend "Censored" "CL" "HW" "EP" / location=
                  outside halign=center;
            else
               if (PLOTCENSORED=1)
                  discretelegend "Censored" / location=inside
                     autoalign=(topright bottomleft);
               endif;
            endif;
            if (PLOTATRISK=1)
               innermargin / align=bottom;
                  blockplot x=TATRISK block=ATRISK / repeatedvalues=
                     true display=(values) valuehalign=start
                     valuefitpolicy=truncate labelposition=left
                     labelattrs=GRAPHVALUETEXT valueattrs=
                     GRAPHDATATEXT (size=7pt) includemissingclass=
                     false;
               endinnermargin;
            endif;
         endlayout;
      else
         entrytitle &entrytitle;
         if (EXISTS(SECONDTITLE))
		  entrytitle  &subtitle / textattrs=GRAPHVALUETEXT;	
         endif;

         layout overlay / xaxisopts=(%if &x_label = %STR() %then %do; %end;%else %do; label=&x_label %end; offsetmin=.05
            linearopts=(viewmax=&xViewmax tickvaluelist= &xTickv
            tickvaluefitpolicy=XTICKVALFITPOL)) yaxisopts=(%if &y_label = %STR() %then %do; label = "Survival Probability" %end;%else %do; label=&y_label %end;
			shortlabel="Survival" linearopts=(
            viewmin=&yViewmin viewmax=1 tickvaluelist=(0 .2 .4 .6 .8 1.0)));
            if (PLOTHW)
               bandplot LimitUpper=HW_UCL LimitLower=HW_LCL x=TIME /
                  group=STRATUM index=STRATUMNUM modelname="Survival"
                  datatransparency=Transparency;
            endif;
            if (PLOTEP)
               bandplot LimitUpper=EP_UCL LimitLower=EP_LCL x=TIME /
                  group=STRATUM index=STRATUMNUM modelname="Survival"
                  datatransparency=Transparency;
            endif;
            if (PLOTCL)
               if (PLOTBAND)
                  bandplot LimitUpper=SDF_UCL LimitLower=SDF_LCL x=
                     TIME / group=STRATUM index=STRATUMNUM modelname=
                     "Survival" display=(outline);
               else
                  bandplot LimitUpper=SDF_UCL LimitLower=SDF_LCL x=
                     TIME / group=STRATUM index=STRATUMNUM modelname=
                     "Survival" datatransparency=Transparency;
               endif;
            endif;
            stepplot y=SURVIVAL x=TIME / group=STRATUM index=
               STRATUMNUM name="Survival" rolename=(_tip1=ATRISK _tip2
               =EVENT) tip=(y x Time _tip1 _tip2);
            if (PLOTCENSORED)
               scatterplot y=CENSORED x=TIME / group=STRATUM index=
                  STRATUMNUM markerattrs=(symbol=plus);
            endif;
            if (PLOTATRISK)
               innermargin / align=bottom;
                  blockplot x=TATRISK block=ATRISK / class=CLASSATRISK
                     repeatedvalues=true display=(label values)
                     valuehalign=start valuefitpolicy=truncate
                     labelposition=left labelattrs=GRAPHVALUETEXT
                     valueattrs=GRAPHDATATEXT (size=7pt)
                     includemissingclass=false;
               endinnermargin;
            endif;
            DiscreteLegend "Survival" / title=&legendvar location=
               outside;
            if (PLOTCENSORED)
               if (PLOTTEST)
                  layout gridded / rows=2 autoalign=(TOPRIGHT
                     BOTTOMLEFT TOP BOTTOM) border=true
                     BackgroundColor=GraphWalls:Color Opaque=true;
                     entry "+ Censored";
                     if (&testp < .0001)
                        entry &testnm " p " eval (
                           PUT(&testp, PVALUE6.4));
                     else
                        entry &testnm " p=" eval (
                           PUT(&testp, PVALUE6.4));
                     endif;
                  endlayout;
               else
                  layout gridded / rows=1 autoalign=(TOPRIGHT
                     BOTTOMLEFT TOP BOTTOM) border=true
                     BackgroundColor=GraphWalls:Color Opaque=true;
                     entry "+ Censored";
                  endlayout;
                  endif;
                  else if (PLOTTEST)
                         layout gridded / rows=1 autoalign=(TOPRIGHT
                         BOTTOMLEFT TOP BOTTOM) border=true
                         BackgroundColor=GraphWalls:Color Opaque=true;
                           if (&testp < .0001) 
                                   entry &testnm " p " eval (PUT(&testp, PVALUE6.4));
                           else entry &testnm " p=" eval (PUT(&testp, PVALUE6.4));
                           endif;
                        endlayout;
                     endif;
                  endif;
               endlayout;
            endif;
         EndGraph;
      end;

   RUN;

%mend PLSurvivalTemp;


