***********************************************************************************************************************
Some title here...

Variable key -
   id:         identifier for each participant
   t:          last observed follow-up time. Includes simulated loss-to-follow-up
   ltfu:       indicator of whether the participant was lost-to-follow-up. Based on simulated data
   delta:      event indicator. 1=AIDs or death, 0=no event
   x:          measured variable related to loss to follow-up and the event. Simulated
   t_true:     observed follow-up time without simulated loss-to-follow-up. Used to create true risk function
   delta_true: event indicator had no loss-to-follow-up occurred. Used to create true risk function

Last edit: Paul Zivich 2019/09/01
**********************************************************************************************************************;

%let source = C:\Users\zivic\Documents\Research\#PZivich\IPCW-tutorial\supplement;
%let end_of_follow_up = 365;

*Reading in data;
PROC IMPORT DATAFILE="&source./actg320_sim-censor.csv" OUT = df DBMS=CSV REPLACE;
RUN;

****************************************************
* 0: True risk functions
***************************************************;

* Data if no censoring had occurred;
PROC PHREG DATA=df NOPRINT;
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

****************************************************
* 1: Bounds on informative censoring
***************************************************;

* Upper bound: all censored have the event at their lost to follow-up time;
DATA df; 
	SET df;
	if ltfu = 1 then delta_upper = 1;
	ELSE delta_upper = delta;
RUN;
PROC PHREG DATA=df NOPRINT;
	MODEL t*delta_upper(0)=;
	BASELINE OUT=upper_risk SURVIVAL=s/METHOD=pl;
RUN;
DATA upper_risk; 
	SET upper_risk; 
	r_upper=1-s; 
	KEEP t r_upper; 
RUN;

* Lower bound: all lost to follow-up don't have the event for full study period;
DATA df; 
	SET df;
	if ltfu = 1 then t_lower = &end_of_follow_up;
	ELSE t_lower = t;
RUN;
PROC PHREG DATA=df NOPRINT;
	MODEL t_lower*delta(0)=;
	BASELINE OUT=lower_risk SURVIVAL=s/METHOD=pl;
RUN;
DATA lower_risk; 
	SET lower_risk; 
	r_lower=1-s; 
	KEEP t_lower r_lower; 
	RENAME t_lower=t;
RUN;

DATA bounds;
	MERGE lower_risk upper_risk;
	BY t;
	* Forward fills missing data;
	RETAIN _r_lower;
	IF NOT missing(r_lower) THEN _r_lower=r_lower;
	ELSE r_lower=_r_lower;
	DROP _r_lower;
RUN;
PROC PRINT DATA=bounds; 
	VAR r_upper r_lower t; 
	TITLE "Bounds for Risk Function"; 
RUN;

****************************************************
* 2: IPCW - pooled logistic model
***************************************************;

* Converting to long data set;
DATA df_long; 
	SET df; 
	id = _n_;
	cat = ceil(t);
	delta_ind = 0;
	not_censored = 1;
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

* Fitting pooled logistic model;
PROC LOGISTIC DATA=df_long DESC NOPRINT; * numerator model;
	MODEL not_censored=t_enter t_enter*t_enter t_enter*t_enter*t_enter; 
	OUTPUT OUT=num P=n;
RUN;

PROC LOGISTIC DATA=df_long DESC NOPRINT; * denominator model;
	MODEL not_censored= x t_enter t_enter*t_enter t_enter*t_enter*t_enter 
					   x*t_enter x*t_enter*t_enter x*t_enter*t_enter*t_enter; 
	OUTPUT OUT=den P=d;
RUN;

data df_long;
	MERGE df_long num den;
	BY id j;
	RETAIN num den;
	IF first.id THEN DO; 
		num=1; 
		den=1; 
		END;
	num = num*n; 
	den = den*d; 
	ipcw = num/den;
RUN;

* Estimating a IPC-weighted Kaplan-Meier;
PROC PHREG DATA=df_long NOPRINT;
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
