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

Last edit: Chase Latour April 20, 2026.
**********************************************************************************************************************;



%*SET FILE PATH FOR DATASET;
%let source = ../data/;


libname rw "&source"; * reading out data to compare to R;






/********************************************************************************************************************

												00 - IMPORT DATA

********************************************************************************************************************/


%*Read in the dataset;
PROC IMPORT DATAFILE="&source.\lau.csv" OUT = lau DBMS=CSV REPLACE;
RUN;












/********************************************************************************************************************

										01 - PREPARE AND DESCRIBE COHORT

********************************************************************************************************************/

%*Create a global variable indicating the end of follow-up (3 years);
%let end_of_fup = 3;




%*Implement administrative censoring at the end of follow-up by assigning events according to what had
occurred by the administrative censoring date.;
data lau;
set lau;
	*Event indicator;
	if t > &end_of_fup then dthev = 0;
		else dthev = dthev;

	*Censoring indicator;
	if t > &end_of_fup then artev = 0;
		else artev = artev;

	if dthev = 1 then eventtype = 2;
		else if artev = 1 then eventtype = 1;
		else eventtype = 0;

	t = int(min(t, &end_of_fup)*365.25); %*convert time to days and make into an integer.;

	%*Censor at ART initiation (create indicator);
	cens = int(eventtype = 1);

run;



*********************
	Descriptive statistics on censoring;

*Get quintiles of CD4 counts for descriptive purposes;
proc univariate data=lau noprint;
    var cd4nadir;
    output out=_quintiles 
        pctlpts=20 40 60 80 100 
        pctlpre=p_;
run;

*Assign quintile categories;
data descrip;
    if _n_ = 1 then set _quintiles;
    set lau;

    length cd4_quintiles 8;

    if cd4nadir <= p_20 then cd4_quintiles = 1;
    else if cd4nadir <= p_40 then cd4_quintiles = 2;
    else if cd4nadir <= p_60 then cd4_quintiles = 3;
    else if cd4nadir <= p_80 then cd4_quintiles = 4;
    else if cd4nadir <= p_100 then cd4_quintiles = 5;
run;

*Look at cross-tabulation ;
proc freq data=descrip;
    tables cd4_quintiles*eventtype / nocol nopercent norow;
run;

proc freq data=descrip;
    tables cd4_quintiles*eventtype / nocol norow;
run;







/********************************************************************************************************************

							02 - CREATE MACRO FOR RESTRICTED QUADRATIC SPLINE TERMS

********************************************************************************************************************/


/*
MACRO: qrspline
PURPOSE: To create columns in thedataset with the restricted quadratic spline terms
INPUT:
- datain = input dataset
- dataout = name of outputdataset
- var = varname
- probs = probabilities
- percentiles_data = dataset from which to derive the percentiles of the variable (important for time)
- percentile_var = variable to calculate percentiles from (important for time)
*/

%*initialize some necessary global variables;
%global percentiles n_knots pmax;

%macro qrspline(datain=, dataout=, var=, probs=20 40 60 80, 
			percentile_var=, percentiles_data=);
	
	%*Calculate the values of &var at those percentiles;
	proc univariate data=&percentiles_data noprint;
	    var &percentile_var;
	    output out=_percentiles 
	        pctlpts=&probs
	        pctlpre=p_;
	run;

	%*Put those variables into macro variables;
	data _null_;
	    set _percentiles;
	    call symputx('percentiles', 
	        catx(' ', p_20, p_40, p_60, p_80)
	    );
	run;
	%put CD4 at knots: &percentiles;

	%*Calculate the number of knots;
	%let n_knots = %sysfunc(countw(&probs));
	%put Number of knots: &n_knots;

	%*Merge the percentile values onto the dataset;
	data _int1;
		set &datain;
		if _n_=1 then set _percentiles;
	run;

	%*Calculate the maximum value;
	data _null_;
	    set _percentiles;

	    array pvars p_:;  /* all variables starting with p_ */
	    max_val = max(of pvars[*]);

	    call symputx('pmax', max_val);
	run;
	%put max &percentile_var at knots: &pmax;


	%*Create the restricted quadratic spline terms;
	data &dataout;
	set _int1;

		%*Calculate the upper tail restriction;
		if &var > &pmax then upper_tail_restriction = (&var - &pmax)**2;
			else upper_tail_restriction = 0;

		%*Create an array of knot values;
		array knots[&n_knots] p_:;

		%*Create an array for the restricted quadratic spline terms;
		array &var._qrs[%eval(&n_knots-1)];

		%*Calculate a restricted quadratic spline term for each person;
		do i=1 to %eval(&n_knots - 1); 
			if &var > knots[i] then &var._qrs[i] = ((&var - knots[i])**2) - upper_tail_restriction;
				else &var._qrs[i] = 0;
		end;

		drop i p_: upper_tail_restriction;

	run;

%mend;






/********************************************************************************************************************

							03 - IMPLEMENT RESTRICTED QUADRATIC SPLINES FOR CD4 COUNT AND TIME

********************************************************************************************************************/

%qrspline(datain=lau, dataout=lau2, var=cd4nadir, percentile_var=cd4nadir, percentiles_data=lau);


****Descriptive statistics: ART initiation, LTFU, and event by CD4 Count;

*Compute cutpoints;
proc univariate data=lau2 noprint;
    var cd4nadir;
    output out=cd4_quintile
        pctlpts=0 20 40 60 80 100
        pctlpre=p_;
run;

*Pull cutpoints into macro variables;
data _null_;
    set cd4_quintile;

    call symputx('p0',  p_0);
    call symputx('p20', p_20);
    call symputx('p40', p_40);
    call symputx('p60', p_60);
    call symputx('p80', p_80);
    call symputx('p100', p_100);
run;

*Create groups;
data lau_cc;
    set lau2;
    length cd4_group 8;

	*Add cd4nadir qrs term that is equal to the original value;
	cd4nadir_qrs0 = cd4nadir;

    if cd4nadir <= &p20 then cd4_group = 1;
	    else if cd4nadir <= &p40 then cd4_group = 2;
	    else if cd4nadir <= &p60 then cd4_group = 3;
	    else if cd4nadir <= &p80 then cd4_group = 4;
	    else if cd4nadir <= &p100 then cd4_group = 5;
run;

*Count outcomes by group;
proc sql;
	select cd4_group, sum(eventtype = 0) as cens_ltfu_admin,
			sum(artev = 1) as cens_art_init,
			sum(dthev = 1) as obs_evt,
			count(*) as n
	from lau_cc
	group by cd4_group
	order by cd4_group
	;
	quit;














/********************************************************************************************************************

											04 - CONVERT DATA TO LONG

Step 1: Create the long dataset. We also create restricted quadratic spline terms for time here.

********************************************************************************************************************/


*STEP 1: Create long dataset ---
-- From the person-level dataset, create a -long- dataset where each row corresponds to one-unit
-- (here, 1-day) of follow-up time for each person.;

DATA lau_long; 
	SET lau_cc; 
	id = _n_;
	cat = ceil(t);
	delta_ind = 0; *Time-varying outcome indicator. Starts as 0 in the first interval;
	not_censored = 1; *Time-varying censoring indicator. Starts as 1 in the first interval;
	DO j=1 to cat;
		t_enter = j-1;
		t_leave = j;
		IF j=cat THEN DO; 
			IF dthev = 1 THEN delta_ind = 1; 
			IF artev = 1 THEN not_censored = 0; 
			END;
		ELSE DO; 
			delta_ind = 0;
			not_censored = 1; 
			END;
		OUTPUT;
	END;
RUN;

/**Look at the first few rows;*/
/*proc print data=lau_long (obs=10);*/
/*run;*/




*Create the restricted quadratic spline terms for time at interval end (t_leave). The percentiles are 
	calculated based on the data before making the daily (long) dataset;
%qrspline(datain=lau_long, dataout=lau_long2, var=t_leave, percentile_var=t, percentiles_data=lau);

*Create a t_leave_qrs0 term for time that is the same as the original. Necessary for restricted
quadratic spline basis function;
data lau_long2;
set lau_long2;
	t_leave_qrs0 = t_leave;
run;


*Look at the distribution of time;
proc means data=lau_long2 min p25 p50 mean p75 max;
	var t_leave_qrs:;
run;














/********************************************************************************************************************

								05 - CONSTRUCT IPCW USING POOLED LOGISTIC MODEL

STEP 2: Exclude intervals where events occur.
	OPTIONAL: If time is modeled as an indicator with interactions with the rest of the model covariates, we can
	remove time points where no censoring happens for anyone. -- THIS IS NOT THE CASE IN THE MAIN TEXT, SO COMMENTED
	OUT HERE.

STEP 3: Estimate the inteval-specific probability of remaining uncensored.

STEP 4: Calculate cumulative probability of remaining uncensored

STEP 5: Calculate interval specific weights.

STEP 6: Evaluate the distribution of the weights at each time point.

********************************************************************************************************************/


*STEP 2: Exclude intervals where events occur. ------
For the pooled logistic regression censoring model, we need a dataset that excludes
intervals where events occur;
data lau_long3;
set lau_long2;
	if delta_ind = 1 then delete;
run;


/**OPTIONAL: Remove those intervals where no censoring occurs. Improves run-time.;*/
/*proc sql;*/
/*	create table n_events as*/
/*	select t_leave, sum(not_censored = 1) as n_censor, count(*) as n_rows*/
/*	from lau_long2*/
/*	group by t_leave*/
/*	;*/
/*	quit;*/
/**/
/**/
/**Subset to those rows where there are some censored persons;*/
/*data n_events2;*/
/*set n_events;*/
/*	where n_censor ne n_rows;*/
/*run;*/
/**/
/**Inner join the long dataset to only those intervals in n_events2;*/
/*proc sql;*/
/*	create table lau_long4 as*/
/*	select a.**/
/*	from lau_long3 as a*/
/*	inner join n_events2 as b*/
/*	on a.t_leave = b.t_leave*/
/*	;*/
/*	quit;*/



*STEP 3: Estimate the inteval-specific probability of remaining uncensored ----;

*We want to create all the interaction terms for the model in this data step;

*Count the number of spline terms for time;
proc sql noprint;
	select count(*) into :n_time trimmed
	from dictionary.columns
	where libname='WORK'
		and memname='LAU_LONG3'
		and upcase(name) like 'T_LEAVE_QRS%'
	;
	quit;
%put Number of spline terms for time: &n_time;

*Count the number of spline terms for CD4 count;
proc sql noprint;
	select count(*) into :n_cd4 trimmed
	from dictionary.columns
	where libname='WORK'
		and memname = 'LAU_LONG3'
		and upcase(name) like 'CD4NADIR_QRS%'
	;
	quit;

%put Number of spline terms for CD4: &n_cd4;

*Calculate the number of interaction terms -- This code is helpful if include interaction terms between the spline basis functions.;
%let nint = %eval(&n_time * &n_cd4);
%put Number of interaction terms: &nint;

data lau_long3;
set lau_long3;

	array time[&n_time] t_leave_qrs:;
	array cd4[&n_cd4] cd4nadir_qrs:;
	array interaction[&nint] int1-int&nint;

	do a=1 to &n_time;
		do b=1 to &n_cd4;
			c = ((a-1)*&n_time)+b;
			interaction[c] = time[a]*cd4[b];
		end;
	end;

	drop a b c;

run;

* Fit pooled logistic model;
PROC LOGISTIC DATA=lau_long3 DESC noprint; 
	MODEL not_censored= t_leave_qrs: cd4nadir_qrs: t_leave*cd4nadir; *int:; 
	OUTPUT OUT=den P=d; *Get the probably of remaining uncensored in each interval;
RUN;
*Look at distribution of interval-specific probabilities;
proc means data=den min p25 p50 mean p75 max missing;
	var d;
run;



*STEP 4: Calculate cumulative probabiltiy of remaining uncensored----
and
STEP 5: Calculate interval-specific weights;

*Merge the predicted probabilities onto the original dataset with all event times;
proc sql;
	create table lau_long4 as
	select a.*, b.d
	from lau_long as a
	left join den as b
	on a.id = b.id and a.j=b.j
	;
	quit;

*Replace missing predicted probabilities with 1;
data lau_long4;
set lau_long4;
	if d = . then d = 1;
		else d = d;
run;

*Calculate the cumulative probability by person ID and the lagged cumulative probability;
proc sort data=lau_long4;
	by id j;
run;

data lau_long5;
set lau_long4;
	by id j;
	
	retain num den lag_den;
	if first.id then do;
		num=1;
		den=1;
		end;
	lag_den = den; *Lag the probabilities so that first row is 1 for all persons.;
	den = den*d;
	ipcw = 1/lag_den; *IPCW are the inverse of the lagged probabilities.;
run;


*STEP 6: Evaluate the distribution of the weights at each time point;
proc means data=lau_long5 min p25 p50 mean p75 max missing;
	var ipcw;
run;

/**Confirm first interval has a weight of 1 for all persons;*/
/*proc freq data=lau_long5 (where = (j=1));*/
/*	table den lag_den ipcw / missing;*/
/*run;*/











/********************************************************************************************************************

						06 - ESTIMATE RISK WITH UNWEIGHTED AND WEIGHTED KAPLAN-MEIER ESTIMATOR

Estimate the risks and output a plot of each risk function.

********************************************************************************************************************/




****************************************************
Unweighted Kaplan-Meier
***************************************************;
PROC PHREG DATA=lau NOPRINT;
	MODEL t*dthev(0)=;
	BASELINE OUT=obs_risk SURVIVAL=s / METHOD=pl;
RUN;
DATA obs_risk; 
	SET obs_risk; 
	r = 1-s; *calculating risk;
	estimator = "Unweighted";
RUN;
PROC PRINT DATA=obs_risk; 
	VAR r t; 
	TITLE "Observed Risk Function"; 
RUN;






****************************************************
Weighted Kaplan-Meier
***************************************************;

**Estimate the IPC-weighted Kaplan-Meier;
PROC PHREG DATA=lau_long5 NOPRINT;
	WEIGHT ipcw;
	MODEL t_leave*delta_ind(0)= / ENTRY=t_enter;
	BASELINE OUT=ipcw_risk SURVIVAL=s / METHOD=pl;
RUN;
DATA ipcw_risk; 
	SET ipcw_risk; 
	r = 1-s; *calculating risk;
	estimator = "Weighted";
	RENAME t_leave=t;
RUN;
PROC PRINT DATA=ipcw_risk; 
	VAR r t; 
	TITLE "IPCW Risk Function"; 
RUN;


*Stack the weighted and unweighted risks on top of each other and look at difference;
data stacked;
set obs_risk ipcw_risk;
run;
proc means data=stacked max;
	class estimator;
	var r;
run;
	






****************************************************
* Creating Plot
***************************************************;

DATA comparison; 
	MERGE ipcw_risk (rename = (r = r_ipcw)) obs_risk (rename = (r = r_obs)); 
	BY t; 
	keep t r_ipcw r_obs; 
RUN;

TITLE;
ODS GRAPHICS on/  ANTIALIASMAX=410500 RESET IMAGENAME = "RiskWeighted" IMAGEFMT =PNG border=off
HEIGHT = 4in WIDTH = 6in ; ods listing style=mystyle image_dpi=150;
PROC SGPLOT DATA=comparison NOAUTOLEGEND;
	STEP X=t Y=r_obs / LINEATTRS=(COLOR=black) LEGENDLABEL="Unweighted" NAME="d0";
	STEP X=t Y=r_ipcw / LINEATTRS=(COLOR=black PATTERN=3) LEGENDLABEL="IPCW" NAME="d3";
	KEYLEGEND "d0" "d3" / ACROSS=1 DOWN=4 LOCATION=inside NOBORDER POSITION=topleft VALUEATTRS=(SIZE=8pt);
	XAXIS LABEL="t (days)" MAX=1096;
	YAXIS LABEL="Risk at t" MIN=0 MAX=0.31;
RUN;

