***********************************************************************************************************************
PROGRAM: sas_informative_censoring.sas
PROGRAMMER: Initial code by Paul Zivich. Edited by Chase.

Variable key -
   id:         identifier for each participant
   t:          last observed follow-up time. Includes simulated loss-to-follow-up
   ltfu:       indicator of whether the participant was lost-to-follow-up. Based on simulated data
   delta:      event indicator. 1=AIDs or death, 0=no event
   x:          measured variable related to loss to follow-up and the event. Simulated
   t_true:     observed follow-up time without simulated loss-to-follow-up. Used to create true risk function
   delta_true: event indicator had no loss-to-follow-up occurred. Used to create true risk function

Last edit: Chase Latour 10/20/2025
**********************************************************************************************************************;



%*SET FILE PATH FOR DATASET;
%let source = C:\Mac\Home\Desktop\IPCW-tutorial\;




/********************************************************************************************************************

												00 - IMPORT DATA

********************************************************************************************************************/

%*Create a global variable indicating the end of follow-up;
%let end_of_follow_up = 365;

*Read in data;
PROC IMPORT DATAFILE="&source.\actg320_sim_censor.csv" OUT = df DBMS=CSV REPLACE;
RUN;






/********************************************************************************************************************

										01 - PREPARE AND DESCRIBE COHORT

********************************************************************************************************************/

%*Create an indicator if LTFU prior to end_of_follow_up;
data df2;
set df;
	ltfu_efup = (ltfu = 1 & t < &end_of_follow_up);
run;

%*Investigate number censored due ot LTFU prior to end_of_follow_up;
proc freq data=df2;
	table ltfu_efup;
run;

%*LTFU stratified by Z;
proc freq data=df2;
	table Z*ltfu_efup;
run;







/********************************************************************************************************************

								02 - CALCULATE TRUE RISK FUNCTIONS IGNORING SIMULATED LTFU

********************************************************************************************************************/


* Data if no censoring had occurred;
PROC PHREG DATA=df2 NOPRINT;
	MODEL t_true*delta_true(0)=;
	BASELINE OUT=true_risk SURVIVAL=s/METHOD=pl;
RUN;
DATA true_risk; 
	SET true_risk; 
	r_true = 1-s; *calculating risk;
	RENAME t_true = t;
RUN;
PROC PRINT DATA=true_risk; 
	VAR r_true t; 
	TITLE "True Risk Function"; 
RUN;








/********************************************************************************************************************

								03 - CONSTRUCT IPCW USING POOLED LOGISTIC MODEL

********************************************************************************************************************/

*STEP 1: Create long dataset ---
-- From the person-level dataset, create a -long- dataset where each row corresponds to one-unit
-- (here, 1-day) of follow-up time for each person.;

DATA df_long; 
	SET df; 
	id = _n_;
	cat = ceil(t);
	delta_ind = 0; *Time-varying outcome indicator. Starts as 0 in the first interval;
	not_censored = 1; *Time-varying censoring indicator. Starts as 1 in the first interval;
	DO j=1 to cat;
		t_enter = j-1;
		t_leave = j;
		IF j=cat THEN DO; 
			IF delta = 1 THEN delta_ind = 1; 
			IF ltfu = 1 THEN not_censored = 0; 
			END;
		ELSE DO; 
			delta_ind = 0;
			not_censored = 1; 
			END;
		OUTPUT;
	END;
RUN;

*Look at the first few rows;
proc print data=df_long (obs=10);
run;



*STEP 2: Create dataset that excludes events ------
For the pooled logistic regression censoring model, we need a dataset that excludes
intervals where events occur;

data df_long2;
set df_long;
	if delta_ind = 1 then delete;
run;




*STEP 3: Estimate inverval-specific censoring probabilities ----;

*First, remove intervals where everyone is uncensored. These aren't needed for model
fitting because IPCW do not change if there is not a censoring event.;

proc sql;
	create table df_long3 as
	select *, count(*) as n_row, sum(not_censored = 1) as n_uncensored
	from df_long2
	group by t_enter
	;
	quit;

data df_long4;
set df_long3;
	where n_row ne n_uncensored;
run;


* Fitting pooled logistic model;
PROC LOGISTIC DATA=df_long4 DESC NOPRINT; * denominator model;
	CLASS t_enter;
	MODEL not_censored= z t_enter z*t_enter; 
	OUTPUT OUT=den P=d; *Get the probably of remaining uncensored in each interval;
RUN;




*STEP 4: Calculate cumulative probabiltiy of remaining uncensored----
and
STEP 5: Calculate inveral-specific weights;


*Calculate the cumulative probabilities and then calculate the IPCW. 
Merge these onto the analytic dataset;
proc sort data=df_long4;
	by id j;
run;
proc sort data=den;
	by id j;
run;
data df_long5;
	MERGE df_long4 den;
	BY id j;
	RETAIN num den;
	IF first.id THEN DO; 
		num=1; 
		den=1; 
		END;
	den = den*d; 
	ipcw = 1/den;
RUN;




/********************************************************************************************************************

										04 - CREATE KAPLAN-MEIER ESTIMATORS

********************************************************************************************************************/


******** Estimating a IPC-weighted Kaplan-Meier;
PROC PHREG DATA=df_long5 NOPRINT;
	WEIGHT ipcw;
	MODEL t_leave*delta_ind(0)= / ENTRY=t_enter;
	BASELINE OUT=ipcw_risk SURVIVAL=s / METHOD=pl;
RUN;

DATA ipcw_risk; 
	SET ipcw_risk; 
	r_ipcw = 1-s; *calculating risk;
	RENAME t_leave=t;
RUN;
PROC PRINT DATA=ipcw_risk; 
	VAR r_ipcw t; 
	TITLE "IPCW Risk Function"; 
RUN;
/* Diagnostic for IPCW
PROC MEANS DATA=df_long;
	VAR ipcw;
RUN;*/

****************************************************
* 3: naive Kaplan-Meier
***************************************************;
PROC PHREG DATA=df NOPRINT;
	MODEL t*delta(0)=;
	BASELINE OUT=obs_risk SURVIVAL=s / METHOD=pl;
RUN;
DATA obs_risk; 
	SET obs_risk; 
	r_obs = 1-s; *calculating risk;
RUN;
PROC PRINT DATA=obs_risk; 
	VAR r_obs t; 
	TITLE "Observed Risk Function"; 
RUN;

****************************************************
* Creating Plot
***************************************************;
DATA comparison; 
	MERGE true_risk bounds ipcw_risk obs_risk; 
	BY t; 
RUN;

TITLE;
ODS GRAPHICS on/  ANTIALIASMAX=410500 RESET IMAGENAME = "RiskWeighted" IMAGEFMT =PNG border=off
HEIGHT = 4in WIDTH = 6in ; ods listing style=mystyle image_dpi=150;
PROC SGPLOT DATA=comparison NOAUTOLEGEND;
	BAND X=t LOWER=r_lower UPPER=r_upper / TYPE=step FILLATTRS=(COLOR=gray TRANSPARENCY=0.75) LEGENDLABEL="Bounds" NAME="d1";
	STEP X=t Y=r_obs / LINEATTRS=(COLOR=black) LEGENDLABEL="Observed" NAME="d0";
	STEP X=t Y=r_ipcw / LINEATTRS=(COLOR=black PATTERN=3) LEGENDLABEL="IPCW" NAME="d3";
	STEP X=t Y=r_true / LINEATTRS=(COLOR=gray PATTERN=1) LEGENDLABEL="True" NAME="d4";
	KEYLEGEND "d4" "d0" "d3" "d1" / ACROSS=1 DOWN=4 LOCATION=inside NOBORDER POSITION=topleft VALUEATTRS=(SIZE=8pt);
	XAXIS LABEL="t (days)" MAX=365;
	YAXIS LABEL="Risk at t" MIN=0 MAX=0.25;
RUN;
