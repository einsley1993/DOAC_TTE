
LIBNAME bt "D:\bytime datasets";
LIBNAME boot "D:\Aim 2 Bootstrap";

%include 'Z:\RPLATT\YI\Aim 2\gformula4.0_Nov2023.sas';


data bytime_baseC_tvC;
  set bt.bytime_primary_g30_baseC_tvC_se;
  drop Apixaban_random Rivaroxaban_random Edoxaban_random;

  /* Define censor for censor due to LTFU: */
  /* LTFU variable was created in index_dates_final, so is 1 for all period if that person is LTFU. */
  by id period;
  censor_LTFU = (LTFU = 1 and last.id and Y ne 1 and death ne 1);
run;

proc means data = bytime_baseC_tvC sum n mean;
    var Rivaroxaban Apixaban Edoxaban censor_LTFU;
run;

data bytime_baseC2;
    set bytime_baseC_tvC;
    by ID period;
    if missing(BMI_cat_b)    then BMI_cat_b    = 'Missing';
    if missing(Smoker_cat_b) then Smoker_cat_b = 'Missing';
    index_year = year(index_date);

    if index_year >= 2019;
run;

/* For spline knots: */
proc univariate data = bytime_baseC2;
var age_b;
output out = age_b_knots
       pctlpts = 5 27.5 50 72.5 95
	   pctlpre = p_
	   pctlname = pct05 pct275 pct50 pct725 pct95;
run;
proc print data = age_b_knots noobs; run;

proc univariate data = bytime_baseC2;
var period;
output out = period_knots
       pctlpts = 5 27.5 50 72.5 95
	   pctlpre = p_
	   pctlname = pct05 pct275 pct50 pct725 pct95;
run;
proc print data = period_knots noobs; run;


proc means data = bytime_baseC2 mean max median p75 p95;
var period;
class index_drug;
run;

data bytime_baseC3;
  set bytime_baseC2;
  %RCSPLINE (period, 1 7 16 27 49);
  %makecat (index_year, 2020 2022, 3);
  %RCSPLINE (age_b, 67 76 80 84 91);
  %makecat(n_hosp,  2 5 9 15 , 5);
  %makecat(n_distinct_meds,  4 7 10 14 , 5);

  sum_drug = sum(Apixaban, Rivaroxaban, Edoxaban);
  if ((
       (index_drug = "Apixaban"    and Apixaban    = 1) or
       (index_drug = "Rivaroxaban" and Rivaroxaban = 1) or
       (index_drug = "Edoxaban"    and Edoxaban    = 1))
       and sum_drug = 1) or excused = 1
       then adhere = 1;
  else adhere = 0;
  censor_deviate = (adhere = 0);
run;

/* bytime_ipcw: one row per person-period, truncated at first deviation
   (the deviation row itself is still included; removed later with censor_deviate=0) */
data bytime_ipcw;
  set bytime_baseC3;
  by id period;
  retain deviated_before;
  if first.id then deviated_before = 0;
  if deviated_before = 0 then output;
  if censor_deviate = 1 then deviated_before = 1;
run;


*****************************************************************************
***  BOOTSTRAP MACRO                                                       ***
***  sample=0  : original data (numberhits=1 for all)  -> point estimates ***
***  sample=1+ : bootstrap resamples                   -> CI distribution  ***
***                                                                        ***
***  Bootstrap is at the PATIENT (ID) level.                              ***
***  numberhits counts how many times each patient was drawn.             ***
***  freq numberhits  in model fitting  = equivalent to row duplication   ***
***  weight numberhits in proc means   = weighted standardization         ***
***  w_final = w_overall * numberhits  in outcome model weight statement  ***
*****************************************************************************;

%macro bootstrap_iptw(nboot = );

    /* ------------------------------------------------------------------ */
    /* Step 1: Unique patient IDs                                          */
    /* ------------------------------------------------------------------ */
    data justids;
        set bytime_ipcw (keep = ID);
        by ID;
        if first.ID;
    run;

    /* ------------------------------------------------------------------ */
    /* Step 2: Expand each ID for sample 0 through nboot                  */
    /* ------------------------------------------------------------------ */
    data justids_boot;
        set justids;
        do sample = 0 to &nboot;
            output;
        end;
    run;

    proc sort data = justids_boot; by sample ID; run;

    /* ------------------------------------------------------------------ */
    /* Step 3: Resample with replacement (URS), stratified by sample      */
    /*         outall keeps un-selected patients with numberhits=0        */
    /* ------------------------------------------------------------------ */
    proc surveyselect data = justids_boot
         method  = urs
         n       = 115160
         seed    = 3104
         out     = justids_boot (drop = selected samplingweight expectedhits)
         outall
         noprint;
        strata sample;
    run;

    /* sample=0 is the original analysis: every patient appears exactly once */
    data justids_boot;
        set justids_boot;
        if sample = 0 then numberhits = 1;
    run;

    proc sort data = justids_boot; by sample ID; run;

    /* ================================================================== */
    /* Bootstrap loop: sample=0 gives point estimates,                    */
    /*                 sample=1 to nboot give the bootstrap distribution  */
    /* ================================================================== */
    %do sample = 0 %to &nboot;

	    %put *** Bootstrap sample &sample of &nboot - starting... ***;

        /* -------------------------------------------------------------- */
        /* Merge this sample's IDs with the full longitudinal data        */
        /* Patients with numberhits=0 (not selected) are dropped          */
        /* -------------------------------------------------------------- */
        data stop_boot;
            merge bytime_ipcw
                  justids_boot (where = (sample = &sample));
            by ID;
            if numberhits = 0 then delete;
        run;

        /* -------------------------------------------------------------- */
        /* IPTW: multinomial logistic at baseline (period = 0) only       */
        /* Denominator model only (nonstabilized weights)                 */
        /* -------------------------------------------------------------- */
        proc logistic data = stop_boot (where = (period = 0))
                      outmodel = model_denom desc noprint;
            ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
            class index_drug   (ref = 'Rivaroxaban')
                  BMI_cat_b    (ref = 'Underweight')
                  Smoker_cat_b (ref = 'Never')
                  AF_type      (ref = 'Others') / param = ref;
            model index_drug =
                age_b_spl1 age_b_spl2 age_b_spl3 
                male white index_year_1 index_year_2
                n_hosp_1 n_hosp_2 n_hosp_3 n_hosp_4
				n_distinct_meds_1 n_distinct_meds_2 n_distinct_meds_3 n_distinct_meds_4
                BMI_cat_b Smoker_cat_b alcoholism_b AF_type
                stroke_b TIA_b HTN_b HF_b CAD_b MI_b PAD_b DM_b CKD_b COPD_b
                Cancer_b dementia_b GI_b Bleed_b Obesity_b Peptic_b
                AntiHTN_b Antiplatelet_b Statins_b PPI_b H2B_b NSAID_b
                SSRI_b Corticosteroids_b Antidiabetic_b Insulin_b
                / link = glogit;
            freq numberhits;
        run;

        proc logistic inmodel = model_denom;
            score data = stop_boot (where = (period = 0))
                  out  = denom_pred_boot;
        run;

        data denom_pred_boot;
            set denom_pred_boot;
            if index_drug = "Apixaban"    then denom_predict = P_Apixaban;
            else if index_drug = "Rivaroxaban" then denom_predict = P_Rivaroxaban;
            else if index_drug = "Edoxaban"    then denom_predict = P_Edoxaban;
            keep ID index_drug period denom_predict;
        run;

        /* -------------------------------------------------------------- */
        /* IPCW: binary logistic for adherence at period > 0              */
        /* Denominator model only (nonstabilized weights)                 */
        /* -------------------------------------------------------------- */
        proc logistic data = stop_boot (where = (period > 0))
                      outmodel = cmodel_denom noprint;
            ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
            class index_drug   (ref = 'Rivaroxaban')
                  BMI_cat_b    (ref = 'Underweight')
                  Smoker_cat_b (ref = 'Never')
                  AF_type      (ref = 'Others') / param = ref;
            model censor_deviate =
                index_drug
                age_b_spl1 age_b_spl2 age_b_spl3 
                male white index_year_1 index_year_2
				n_hosp_1 n_hosp_2 n_hosp_3 n_hosp_4
				n_distinct_meds_1 n_distinct_meds_2 n_distinct_meds_3 n_distinct_meds_4
                BMI_cat_b Smoker_cat_b alcoholism_b AF_type
                stroke_b TIA_b HTN_b HF_b CAD_b MI_b PAD_b DM_b CKD_b COPD_b
                Cancer_b dementia_b GI_b Bleed_b Obesity_b Peptic_b
                AntiHTN_b Antiplatelet_b Statins_b PPI_b H2B_b NSAID_b
                SSRI_b Corticosteroids_b Antidiabetic_b Insulin_b
                HTN HF CAD MI DM PAD CKD COPD Cancer Dementia Obesity Peptic
                period_spl1 period_spl2 period_spl3;
            freq numberhits;
        run;

        proc logistic inmodel = cmodel_denom;
            score data = stop_boot (where = (period > 0))
                  out  = cdenom_pred_boot;
        run;

        data cdenom_pred_boot;
            set cdenom_pred_boot;
            cdenom_predict = 1 - P_1;   /* P(censor_deviate=0) = P(adhere) */
            keep ID index_drug period cdenom_predict;
        run;

        /* -------------------------------------------------------------- */
        /* IPCW: binary logistic for LTFU                                 */
        /* Fit among rows still uncensored by deviation (censor_deviate=0)*/
        /* so LTFU is modeled conditional on remaining adherent.          */
        /* Period > 0 only; LTFU weight factor is set to 1 at period 0.   */
        /* Same confounders as the deviation IPCW.                        */
        /* Denominator model only (nonstabilized weights)                 */
        /* -------------------------------------------------------------- */
        proc logistic data = stop_boot (where = (period > 0 and censor_deviate = 0))
                      outmodel = lmodel_denom noprint;
            ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
            class index_drug   (ref = 'Rivaroxaban')
                  BMI_cat_b    (ref = 'Underweight')
                  Smoker_cat_b (ref = 'Never')
                  AF_type      (ref = 'Others') / param = ref;
            model censor_LTFU =
                index_drug
                age_b_spl1 age_b_spl2 age_b_spl3 
                male white index_year_1 index_year_2
				n_hosp_1 n_hosp_2 n_hosp_3 n_hosp_4
				n_distinct_meds_1 n_distinct_meds_2 n_distinct_meds_3 n_distinct_meds_4
                BMI_cat_b Smoker_cat_b alcoholism_b AF_type
                stroke_b TIA_b HTN_b HF_b CAD_b MI_b PAD_b DM_b CKD_b COPD_b
                Cancer_b dementia_b GI_b Bleed_b Obesity_b Peptic_b
                AntiHTN_b Antiplatelet_b Statins_b PPI_b H2B_b NSAID_b
                SSRI_b Corticosteroids_b Antidiabetic_b Insulin_b
                HTN HF CAD MI DM PAD CKD COPD Cancer Dementia Obesity Peptic
                period_spl1 period_spl2 period_spl3;
            freq numberhits;
        run;

        proc logistic inmodel = lmodel_denom;
            score data = stop_boot (where = (period > 0 and censor_deviate = 0))
                  out  = ldenom_pred_boot;
        run;

        data ldenom_pred_boot;
            set ldenom_pred_boot;
            ldenom_predict = 1 - P_1;   /* P(censor_LTFU=0) = P(remain in follow-up) */
            keep ID index_drug period ldenom_predict;
        run;

        /* -------------------------------------------------------------- */
        /* Merge predictions back and compute nonstabilized weight        */
        /* -------------------------------------------------------------- */
        proc sort data = stop_boot;        by ID index_drug period; run;
        proc sort data = denom_pred_boot;  by ID index_drug period; run;
        proc sort data = cdenom_pred_boot; by ID index_drug period; run;
        proc sort data = ldenom_pred_boot; by ID index_drug period; run;

        data stop_boot2;
            merge stop_boot (in=a)
                  denom_pred_boot
                  cdenom_pred_boot
                  ldenom_pred_boot;
            by ID index_drug period;
            if a;

            /* IPTW (baseline only) */
            if period = 0 then wa_temp = 1 / denom_predict;

            /* IPCW for deviation (period > 0; 1 at period 0) */
            if period = 0                     then wc_temp = 1;
            else if period > 0 and adhere = 1 then wc_temp = 1 / cdenom_predict;
            else if period > 0 and adhere = 0 then wc_temp = 0;

            /* IPCW for LTFU (period > 0, among rows not censored by deviation) */
            if censor_deviate = 1      then wl_temp = 0;
            else if censor_LTFU = 1    then wl_temp = 0;
            else if period = 0         then wl_temp = 1;
            else                            wl_temp = 1 / ldenom_predict;
        run;

        /* -------------------------------------------------------------- */
        /* Cumulative product of weights across periods, within patient   */
        /* w_overall = IPTW x cumulative IPCW(deviation) x cumulative     */
        /*             IPCW(LTFU)                                         */
        /* -------------------------------------------------------------- */
        proc sort data = stop_boot2; by ID period; run;

        data stop_boot3;
            set stop_boot2;
            by ID;
            retain w_ac;
            if first.ID then w_ac = wa_temp * wc_temp * wl_temp;   /* period 0: IPTW only (wc_temp = wl_temp = 1) */
            else             w_ac = w_ac * wc_temp * wl_temp;      /* period > 0: accumulate both IPCWs */
            w_overall = w_ac;
        run;

        /* -------------------------------------------------------------- */
        /* Truncate w_overall at weighted 99th percentile                 */
        /* Percentile computed on uncensored rows only (the rows that     */
        /* enter the outcome model), so zero weights on censored rows do  */
        /* not pull the cutoff down                                       */
        /* -------------------------------------------------------------- */
        proc univariate data = stop_boot3 (where = (censor_deviate = 0 and censor_LTFU = 0)) noprint;
            var w_overall;
            output out = w_pct pctlpts = 99 pctlpre = P_;
            weight numberhits;
        run;

        data stop_boot4;
            if _n_ = 1 then set w_pct;
            set stop_boot3;
            if w_overall > P_99 then w_overall = P_99;
            drop P_99;
        run;

        /* -------------------------------------------------------------- */
        /* Censored dataset: remove the deviation row (censor_deviate=1)  */
        /* and the LTFU row (censor_LTFU=1)                               */
        /* w_final = w_overall * numberhits combines causal weight with   */
        /* bootstrap resampling count                                     */
        /* -------------------------------------------------------------- */
        data stop_boot_cens;
            set stop_boot4;
            where censor_deviate = 0 and censor_LTFU = 0;
            length strategy $12;
            strategy = index_drug;
            w_final = w_overall * numberhits;
        run;

        /* -------------------------------------------------------------- */
        /* Outcome model: cause-specific hazard for stroke/SE             */
        /* -------------------------------------------------------------- */
        proc logistic data = stop_boot_cens outmodel = boot_model desc noprint;
            ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
            class strategy    (ref = 'Rivaroxaban')
                  BMI_cat_b   (ref = 'Underweight')
                  Smoker_cat_b (ref = 'Never')
                  AF_type     (ref = 'Others') / param = ref;
            model Y =
                strategy
                age_b_spl1 age_b_spl2 age_b_spl3 male  
                HF_b CKD_b CHA2DS2VASc

                period_spl1 period_spl2 period_spl3
                strategy*period_spl1 strategy*period_spl2 strategy*period_spl3;
            weight w_final;
        run;

        /* -------------------------------------------------------------- */
        /* Outcome model: cause-specific hazard for the competing event   */
        /* (death). Same risk set (stop_boot_cens already truncates at    */
        /* whichever of Y=1 / death=1 / deviation occurs first), same     */
        /* covariates and weight as the stroke/SE model above.            */
        /* -------------------------------------------------------------- */
        proc logistic data = stop_boot_cens outmodel = boot_model_death desc noprint;
            ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
            class strategy    (ref = 'Rivaroxaban')
                  BMI_cat_b   (ref = 'Underweight')
                  Smoker_cat_b (ref = 'Never')
                  AF_type     (ref = 'Others') / param = ref;
            model death =
                strategy
                age_b_spl1 age_b_spl2 age_b_spl3 male  
                HF_b CKD_b CHA2DS2VASc

                period_spl1 period_spl2 period_spl3
                strategy*period_spl1 strategy*period_spl2 strategy*period_spl3;
            weight w_final;
        run;

        /* -------------------------------------------------------------- */
        /* Standardization: expand each baseline person across            */
        /* period 0-60 x 5 strategies, score both cause-specific hazards, */
        /* then compute the competing-risks cumulative incidence (CIF)    */
        /* -------------------------------------------------------------- */
        /* Baseline population taken from stop_boot4 (not stop_boot_cens) */
        /* so persons censored by LTFU at period 0 remain in the          */
        /* standardization population                                     */
        data boot_pred;
            set stop_boot4 (where = (period = 0));
            length strategy $12;
            do period = 0 to 60;
                %RCSPLINE(period, 1 7 16 27 49);
				%makecat (index_year, 2020 2022, 3);
				%RCSPLINE (age_b, 67 76 80 84 91);
				%makecat(n_hosp,  2 5 9 15 , 5);
				%makecat(n_distinct_meds,  4 7 10 14 , 5);

                strategy = "Apixaban";    output;
                strategy = "Rivaroxaban"; output;
                strategy = "Edoxaban";    output;
            end;
        run;

        proc logistic inmodel = boot_model;
            score data = boot_pred
                  out  = boot_pred_stroke (keep = ID period strategy numberhits P_1
                                           rename = (P_1 = hazard_stroke));
        run;

        proc logistic inmodel = boot_model_death;
            score data = boot_pred
                  out  = boot_pred_death (keep = ID period strategy P_1
                                          rename = (P_1 = hazard_death));
        run;

        proc sort data = boot_pred_stroke; by strategy ID period; run;
        proc sort data = boot_pred_death;  by strategy ID period; run;

        data boot_pred;
            merge boot_pred_stroke boot_pred_death;
            by strategy ID period;
        run;

        /* Discrete-time Aalen-Johansen: S(k) = S(k-1)*(1-h_stroke(k)-h_death(k)),
           CIF_stroke(k) = CIF_stroke(k-1) + S(k-1)*h_stroke(k), same for death */
        data boot_pred;
            set boot_pred;
            by strategy ID;
            retain surv risk_stroke risk_death;
            if period = 0 then do;
                surv        = 1;
                risk_stroke = 0;
                risk_death  = 0;
            end;
            risk_stroke = risk_stroke + surv * hazard_stroke;
            risk_death  = risk_death  + surv * hazard_death;
            surv        = surv * (1 - hazard_stroke - hazard_death);
        run;

        proc means data = boot_pred noprint;
            class strategy period;
            var risk_stroke risk_death;
            types strategy*period;
            output out = boot_mean (keep = strategy period risk_stroke risk_death)
                   mean = risk_stroke risk_death;
            weight numberhits;
        run;

        data boot_A boot_R boot_E;
            set boot_mean;
            if strategy = "Apixaban"    then output boot_A;
            if strategy = "Rivaroxaban" then output boot_R;
            if strategy = "Edoxaban"    then output boot_E;
            keep period risk_stroke risk_death;
        run;

        data boot_results;
            merge
                  boot_A (rename = (risk_stroke = risk_A   risk_death = riskD_A))
                  boot_R (rename = (risk_stroke = risk_R   risk_death = riskD_R))
                  boot_E (rename = (risk_stroke = risk_E   risk_death = riskD_E));
            by period;
            sample = &sample;   /* 0 = original data, 1+ = bootstrap */
        run;

        %if &sample = 0 %then %do;
            data boot_results_combined;
                set boot_results;
            run;
        %end;
        %else %do;
            proc append base = boot_results_combined data = boot_results;
            run;
        %end;

        /* Save the output dataset for plotting */
		data boot.boot_primary_post2019_LTFU;
		  set boot_results_combined;
		run;

    %end; /* end bootstrap loop */

	%put *** Completed bootstrap sample &sample of &nboot ***;

    proc sort data = boot_results_combined; by sample period; run;

%mend bootstrap_iptw;


options nomprint nonotes nospool;
%bootstrap_iptw(nboot = 500);
options mprint notes spool;

/* Point estimates */
data point_est;
    set boot.boot_primary_post2019_LTFU (where = (sample = 0));
	keep period risk_R risk_A risk_E;
run;


/* 95% CI from bootstrap  */
proc sort data = boot.boot_primary_post2019_LTFU; by period; run;
