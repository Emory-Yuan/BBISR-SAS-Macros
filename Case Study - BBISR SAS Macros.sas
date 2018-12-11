/*****************************************************************************-********************
 ******************** BBISR SAS macros Case Study - Routine Data Analyses  ************************
 *************************** By: Yuan Liu, PhD*****************************************************  
 ****************************** Nov 2018***********************************************************

Citation:

Yuan Liu, Dana Nickleach, Chao Zhang,  Jeffrey Switchenko, Jeanne Kowalski. Carry Out Streamlined Routine 
Data Analyses with Reports for Observational Studies: Introduction to a Series of Generic SAS ® Macros. 
Oct 2018. Preprint. DOI:10.13140/RG.2.2.25265.84329

Purpose: 

This is a case study tutorial about conducting data analysis for a typical observational study 
by using SAS macros developed by the Biostatistics and Bioinformatics Shared Resource at Winship
Cancer Institute (BBISR). The analysis covers a streamlined routine data analysis,
including descriptive statistics, univariate associations, multivariable models, subgroup analyses,
and propensity score approaches. The macros will directly deliver the summary report in RTF format.


/*************************/
/*** 1. Initial Set up ***/
/*************************/

COMMENT Clear current SAS output and log windows; 
dm 'log;clear;output;clear;odsresults;clear';

COMMENT Load needed SAS macros. Change the path to the location where they are stored ; 
%let macrodir = H:\location you used to save all SAS macros\;
%include "&macrodir.DESCRIPTIVE V16.sas";
%include "&macrodir.UNI_CAT V30.sas";
%include "&macrodir.UNI_NUM V9.sas";
%include "&macrodir.UNI_LOGREG V15.sas";
%include "&macrodir.UNI_PHREG V26.sas";
%include "&macrodir.MULTIPLE_LOGREG V16.sas"; * This macro is required for %LOGREG_SEL to run properly;
%include "&macrodir.MULTIPLE_PHREG V21.sas";  * This macro is required for %PHREG_SEL to run properly;
%include "&macrodir.LOGREG_SEL V15.sas";
%include "&macrodir.PHREG_SEL V23.sas";
%include "&macrodir.PlSurvivalTemp V5.sas";* This macro is required for %KM_PLOT to run properly;
%include "&macrodir.KM_PLOT V23.sas";
%include "&macrodir.Table Template.sas"; 


COMMENT TIP 1: create macro variable &DIR to specify the location used to store all results. 
	    If you are conducting an updated analysis, you just need to change this path and all 
	    updated results will be saved to a new location, and you donï¿½t have to modify it in every macro.;

%let dir = H:\location you used to save all output results\;


/**********************************************************************************/
/***  Prepare for final analytic dataset    *************************************/
/**********************************************************************************/

/*** CAUTION:  1.For purpose of illustration of SAS macro usage, only a few relevant variables were selected; 
		   2.The final analytic dataset is supplied for education purpose  */

/*** load case study data and you may need to change the path to the location where you saved the data ***/
libname casedata 'H:\location you used to store data';

PROC FORMAT ;
    value fsex 	1='Male' 2='Female';
	value fapproach 
			0='No surgical procedure of primary site'
			1='Robotic assisted' 
			2='Robotic converted to open' 
			3='Laparoscopic'
			4='Laparoscopic converted to open'
			5='Open or approach unspecified'
			9='Unknown if surgery performed';
	value f30day
			0='Patient alive, or died more than 30 days after surgery performed'
			1='Patient died 30 or fewer days after surgery performed'
			9='Patient alive with fewer than 30 days of follow-up, surgery date missing, or last contact date missing';
	value fvital 0='Dead' 	1='Alive';
	value surg_app   0='Open'  1='Minimally invasive';
	value furbanrural	1-3 = "Metro" 4-7 = "Urban" 8-9 = "Rural";
	value fcombidity  2 = '2+';
  	value $ clin_T '1', '1A','1B' = '1'  '2', '2A','2B' = '2' ;
	value yesno 1= 'Yes' 0 = 'No' 9 ='Unknown';
RUN;


DATA anal; set casedata.Case_Study_Data;

	*** Define Surgical Approach: open vs. minimally invasive ;
	if 1 <= RX_HOSP_SURG_APPR_2010 <= 4 then surg_app = 1;
	else if RX_HOSP_SURG_APPR_2010 = 5 then surg_app = 0;
	format surg_app surg_app.;
	label  surg_app = "Surgical Approach";

	*** Urban/Rural;
    format UR_CD_03 furbanrural.;

	**Combidity;
	format CDCC_TOTAL fcombidity.;

	***  AJCC Clinic T stage;
	format TNM_CLIN_T clin_T.;

	*** Simplify Histology to a small number but clinical meaningful categories ;
	if  8050<=histology <=8123 then His_cat ='Squamous cell carcinomas       ' ;
	else if   histology in (8140, 8141, 8572, 8143, 8144, 8146, 8147) then His_cat = 'Adenocarcinomas';
	else if 8250 <=histology <=8323 or 8480<=histology <=8550 then His_cat = 'Adenocarcinomas';
	else His_cat = 'Other or Unknown';
	label His_cat = 'Histology';

	*** Define Outcome as survival months from date of surgery to death or last follow up;

	** OS censor: 1=dead, 0=censored *;
	if PUF_VITAL_STATUS = 0 then OS_censor = 1; 
	else OS_censor = 0;

	os_surg= ROUND(DX_LASTCONTACT_DEATH_MONTHS - DX_DEFSURG_STARTED_DAYS/(365.2425/12),0.01);
	label os_surg = 'Months (OS)';

	** 30-day Mortality;
	if PUF_30_DAY_MORT_CD = 9 then PUF_30_DAY_MORT_CD = .;
	format PUF_30_DAY_MORT_CD yesno.;

RUN;


/**************************************************/
/*** 6. Routine Data Analysis by SAS macros *******/
/**************************************************/

ods listing close;
ods html close;

/* TIP 2: create macro variables &cat_var and &num_var to store categorical and numerical variable lists 
outside macros, and reference them in all related macros. Then changes to the variable lists, only requires 
changes to those two macro variables to get all related tables updated. 

Tip 3: Also note that the order of variables in CLIST is the same as the order in which they will appear in 
the final output table. Listing similar types of variables together, e.g. demographic variables, tumor characteristic 
variables, lab values, etc., will improve the readability of the results.  */

%let cat_var = SEX UR_CD_03 CDCC_TOTAL YEAR_OF_DIAGNOSIS  TNM_CLIN_T;
%let num_var = AGE ;

TITLE 'Table 1 Descriptive Statistics for All Variables';
%DESCRIPTIVE(DATASET=anal, 
     CLIST = surg_app &cat_var,
	 NLIST = &num_var,
     OUTPATH= &dir, 
     FNAME=Table 1 Descriptive Statistics for All Variables,
     DEBUG=F); 
TITLE;

/* TIP 4: In the following example, the NONPAR option is turned off because it will 
take more time to compute Fisherï¿½s exact test on a large sample. */

data one;set anal;
if sex=1 and CDCC_TOTAL = 0 and _N_ <=1000;run;

TITLE 'Table 2 Univariate Association with Study Cohort';
%UNI_CAT (
     DATASET = one, 
     OUTCOME = surg_app, 
     CLIST = &cat_var,
     NLIST = &num_var,
     NONPAR = T, BY = his_cat,
     ROWPERCENT = F,
     SPREAD = T, OUTDATA=Tabel_2, DOC=F,
     OUTPATH = &dir,DEBUG=T,
     FNAME =Table 2 Univariate Association with Study Cohort); 
TITLE;

/* TIP 5: If you need to specify the reference level of a categorical variable, separate CLIST by ï¿½*ï¿½ and add (DESC) 
or (ref = ï¿½Reference level in formatted valueï¿½) after each desired variable name. Otherwise CLIST can be separated 
by space alone. */

%let cat_var_ref = SEX* UR_CD_03(DESC)*CDCC_TOTAL(DESC)*YEAR_OF_DIAGNOSIS*HIS_CAT(ref="Adenocarcinomas")*TNM_CLIN_T(DESC);

TITLE 'Table 3 Univariate Association with 30-day Mortality';
%UNI_LOGREG(
	DATASET = anal, 
	OUTCOME = PUF_30_DAY_MORT_CD,
	EVENT = 'Yes', 
	CLIST = surg_app* &cat_var_ref, 
	NLIST = &num_var,
	OUTPATH = &dir, 
	FNAME = Table 3 Univariate Association with 30-day Mortality);
TITLE;

TITLE 'Table 4 Univariate association with Overall Survival';
%UNI_PHREG (dataset = anal, 
	EVENT = os_surg,  
	CENSOR = os_censor,
    CLIST = surg_app * &cat_var_ref ,
	NLIST = &num_var ,
    LOGRANK=F, 
    TYPE3=T,
    PHA = F, 
	ORIENTATION = PORTRAIT,
	OUTPATH = &dir, 
	FNAME = Table 4 Univariate association with Overall Survival);
TITLE;

TITLE 'Table 5 Multivariable Logistic Regression Model for 30-day Mortality';
%LOGREG_SEL(
		DSN = anal,
   		OUTCOME = PUF_30_DAY_MORT_CD,
   		EVENT = 'Yes',
		VAR = surg_app  &cat_var &num_var,
		CVAR = surg_app*&cat_var_ref,
   		INC = 1, SLSTAY = 0.1,
		CLNUM = F,
        OUTPATH=&dir,
        FILENAME=Table 5 Multivariable Logistic Regression Model for 30-day Mortality);
TITLE;

TITLE 'Table 6 Multivariable Cox Proportional Hazard Model for Overall Survival';
%PHREG_SEL(
	DSN = anal,
	EVENT = os_surg,  
	CENSOR = os_censor,
	VAR = surg_app &cat_var &num_var,
	CVAR = surg_app*&cat_var_ref ,
	INC = 1,SLSTAY = .10,CLNUM = T,
    OUTPATH = &dir,
    FILENAME = Table 6 Multivariable Cox Proportional Hazard Model for Overall Survival);
TITLE;

TITLE 'Table 7 Multivariable Cox Model for OS Stratify by Histology';
%PHREG_SEL(
	DSN=anal,
	EVENT = os_surg,  
	CENSOR = os_censor,
	VAR = surg_app|his_cat  &cat_var &num_var,
	CVAR = surg_app  &cat_var ,
	INC = 3,
	SLSTAY = .10,
	EFFECT = surg_app,
	SLICEBY = his_cat,
    OUTPATH = &dir,
    FILENAME = Table 7 Multivariable Cox Model for OS Stratify by Histology);
TITLE;


%KM_PLOT(DSN = anal, 
	EVENTS = os_surg, 
	CENSORS = os_censor, 
	GRPLIST = surg_app, 
	TITLE = "Figure 1 KM plot of Overall Survival by Surgical Approach",	
    UNIT = Month, JOIN = T, PLOT = T, TABLE = T,
    TIMELIST = 12 60, NONCENSORED = T,
	XTICK = (0 12 24 36 48 60), XMAX = 60, ATRISK = T,
    OUTPATH = &dir,
    FNAME = Figure 1 KM plot of Overall Survival by Surgical Approaches);

