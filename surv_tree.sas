%let main_dir = C:\Temp\covid;
%let total_procs = 10;

LIBNAME SURVLIB "&main_dir";

options AUTOSIGNON sascmd='!sascmd'; /* nosource nosource2 nonotes;*/

ods graphics on;

ods trace off;

ods pdf file="&main_dir./res.pdf";
%include "&main_dir/check_var.sas";
%include "&main_dir/build_tree.sas";
%macro add_trends(dts_in,dts_out,ivars);
%local ids n i s vars ;


data ttt;
set &dts_in;
tm=datdif(datepart(atime),datepart(RE_DATE),'ACT/ACT');
run;

%let vars=;

proc sql;
create table &dts_out (

%let n=%sysfunc(countw(&ivars,%str(|)));
%put n=&n vars=&ivars;
%do i=1 %to &n;
	%let s=%scan(&ivars,&i,%str(|));
	%let vars= &vars trnd_&s.;
	trnd_&s. float,
%end;
 patient_id int);
quit;

proc sql;
select distinct patient_id into: ids separated by '|' 
from &dts_in;
quit;

%let n=%sysfunc(countw(&ids,%str(|)));
%put n=&n vars=&ids;

%do i=1 %to &n;
	%let s=%scan(&ids,&i,%str(|));
	ods exclude all;
	ods output PearsonCorr=corrs;
	proc corr data=ttt;
		where patient_id = &s.;
		var %sysfunc(TRANWRD(&ivars,%str(|),%str( )));
		with tm;
	run;

	proc sql;
	insert into &dts_out select %sysfunc(TRANWRD(&ivars,%str(|),%str(,))), &s. from corrs;
	quit;

%end;
%mend;


%macro prep_patients(dts_in, dts_out, ivars, survtm, cns);
%local n i s s0 labels;

%let n=%sysfunc(countw(&ivars,%str(|)));
%put n=&n vars=&ivars;

%do i=1 %to &n;
	%let s=%scan(&ivars,&i,%str(|));
	data _NULL_;
	set &dts_in;
	call symput("lbl",vlabel(&s));
	run;
	%let labels = &labels | &lbl;
%end;

%let labels = %substr(&labels,3);

%put labels = &labels;
	
data &dts_out;
set &dts_in;
retain &survtm;
%do i=1 %to &n;
	%let s=%scan(&ivars,&i,%str(|));
	%let s0 = &s0 av_&s mn_&s mx_&s c_&s;
	label av_&s = "AVG of" %scan(&labels,&i,%str(|));
	label mx_&s = "MAX of" %scan(&labels,&i,%str(|));
	label mn_&s = "MIN of" %scan(&labels,&i,%str(|));
	label c_&s = "CNT of" %scan(&labels,&i,%str(|));
%end;
%put s0=&s0;
retain &s0;

by PATIENT_ID;
if first.PATIENT_ID then do;
&survtm = datdif(datepart(atime),datepart(dtime),'ACT/ACT');
put atime= dtime= &survtm=;

%do i=1 %to &n;
	%let s=%scan(&ivars,&i,%str(|));
	av_&s = &s;
	mn_&s = &s;
	mx_&s = &s;
	c_&s = (&s ne .);
%end;

end;
else; do;
%do i=1 %to &n;
	%let s=%scan(&ivars,&i,%str(|));
	if &s ne . then do;
		c_&s + 1;
		av_&s + &s;
		if (mn_&s eq .) or (mn_&s > &s) then mn_&s = &s;
		if (mx_&s eq .) or (mx_&s < &s) then mx_&s = &s;
	end;
%end;
end;

if last.PATIENT_ID then do;

%do i=1 %to &n;
	%let s=%scan(&ivars,&i,%str(|));
	if (c_&s) then av_&s = av_&s/c_&s;
%end;

output;
end;
run;

%add_trends(&dts_in,trnds,&ivars);

data &dts_out;
merge &dts_out trnds;
by PATIENT_ID;
run;

proc sort data=&dts_out out=samp_pats;
by &cns;
run;

proc surveyselect data=samp_pats method=srs out=samp_pats outall samprate=0.5;
strata &cns;
run;

data trn_&dts_out. hld_pats;
set samp_pats;
if selected then output trn_&dts_out.; else output hld_pats;
drop selected;
run;

proc surveyselect data=hld_pats method=srs out=samp_pats outall samprate=0.5;
strata &cns;
run;

data vld_&dts_out. tst_&dts_out.;
set samp_pats;
if selected then output vld_&dts_out.; else output tst_&dts_out.;
drop selected;
run;

proc sort data=trn_&dts_out.;
by PATIENT_ID;
run;

proc sort data=vld_&dts_out.;
by PATIENT_ID;
run;

proc sort data=tst_&dts_out.;
by PATIENT_ID;
run;

/*

proc sql;
cteate table tst_&dts_out. as 
select * from &dts_out where (atime>'10FEB20:00:00:00'dt);
cteate table vld_&dts_out. as 
select * from &dts_out 
where (atime<='10FEB20:00:00:00'dt) and (atime>'05FEB20:00:00:00'dt);
cteate table trn_&dts_out. as 
select * from &dts_out 
where (atime<='05FEB20:00:00:00'dt);
quit;

proc surveyselect data=hld_pats method=srs out=samp_pats outall samprate=0.75;
strata &cns;
run;

data vld_&dts_out. trn_&dts_out.;
set samp_pats;
if selected then output trn_&dts_out.; else output vld_&dts_out.;
drop selected;
run;*/

proc stdize data=&dts_out. out=imp_&dts_out. reponly;
run;

proc sql;
select PATIENT_ID into : vld_ids separated by ',' from vld_&dts_out.;
select PATIENT_ID into : tst_ids separated by ',' from tst_&dts_out.;
quit;

data imp_&dts_out.;
set imp_&dts_out.;
length chunk $ 6;
true_&survtm. = &survtm;
if PATIENT_ID in (&vld_ids) then do;
chunk='valid'; &survtm = .;
end; else if PATIENT_ID in (&tst_ids) then do;
chunk='test'; &survtm = .;
end; else chunk='train';
run;

%mend;

%macro best_phreg(dts,tgt, cns, vcns, cvars,ivars,phreg_steps);
%local phreg_effects s n i e;


%let cvars=%sysfunc(TRANWRD(&cvars,%str(|),%str( )));
%let ivars=%sysfunc(TRANWRD(&ivars,%str(|),%str( )));

%put best_phreg cvars="&cvars" ivars="&ivars";

ods output ModelBuildingSummary=phreg_Summary;
proc phreg data=&dts plots=roc rocoptions(IAUC at=2 to 30 by 5);
   class &cvars;
   model &tgt*&cns(&vcns)= &cvars &ivars / selection=forward maxstep=100;
run;

proc sql;
select EffectEntered into: phreg_effects separated by ' ' from phreg_Summary;
quit;

proc sql;
create table &phreg_steps (step int, IAUC float, vars varchar(1000));
quit;

%let n=%sysfunc(countw(&phreg_effects,%str( )));
%let s =;

data tmp_data;
set &dts;
run;

%do i=1 %to &n;
	%let s = &s %scan(&phreg_effects,&i,%str( ));

	ods exclude all;
	proc phreg data=tmp_data plots=none;
    	class &cvars;
   		model &tgt*&cns(&vcns)= &s;
   		output out=tmp_data xbeta=xb_&i.;	
	run;

	ods output IAUC=iauc;
	proc phreg data=tmp_data plots=none rocoptions(IAUC);
		where chunk eq 'valid';
   		model true_&tgt*&cns(&vcns)= &s / nofit;
		roc pred=xb_&i.;
	run;

	proc sql;
		select Estimate into : e from iauc;
		insert into &phreg_steps values (&i, &e, "&s");
	quit;
%end;
ods select all;
proc sgplot data=&phreg_steps;
	series x=step y=IAUC;
run;

%mend;
%macro main();


%local cvars ivars vcns total_vars forest_sz;

%let forest_sz=10;
%let cvars = gender | CoV_nacid ;
%let base_ivars = hctrop | hemoglobin | schloride | ptime | procalc | eosinoph | Interl2 | alkaline | albumin | basophil | Inter10 | totbilirubin | platelet | mono | antithrombin | Inter8 | indbilirubin | RBCD | neutro | totprotein | qtrep | prothromb | HBsAg | mcorpvol | hematocrit | WBCC | tumor | MCHC | fibrinogen | Inter1 | Urea | LYM | PH | RBCC | eosinophil | ccalcium | spotassium | glucose | NC | dbilirubin | MPW | ferritin | RBCSD | Thromtime | lymphocyte | HCV_antibody | Ddimer | totcholesterol | aspamin | uricacid | HCO3 | calcium | NTproBNP | Lactdehydr | PLCR | Inter6 | Fibrindeg | monocytes | PLTD | globulin | glutamyl | ISR | basophilc | MCH | apttime | CRP | HIV_antibody | ssodium | thrombocytocrit | ESR | QPT | eGFR | creatinine ;

%let vcns = 0;

data patients0 (drop= id);
set SURVLIB.patients;
retain id;
if PATIENT_ID eq . then PATIENT_ID = id; else id = PATIENT_ID;
if RE_DATE eq . then RE_DATE=atime;
run;

%prep_patients(patients0,patients,&base_ivars,surv_time,outcome);


%let total_vars=%sysevalf(%sysfunc(countw(&base_ivars,%str(|)))*5+%sysfunc(countw(&cvars,%str(|))));

%let p_skip = %sysevalf(0.75);

%local i n s s0;
%let s0=;
%let n=%sysfunc(countw(&base_ivars,%str(|)));
%do i=1 %to &n;
	%let s=%scan(&base_ivars,&i,%str(|));
	/*%let s0 = &s0 | av_&s | mn_&s | mx_&s | trnd_&s;*/
	
	%let s0 = &s0 | av_&s | trnd_&s;
%end;
%let ivars= %substr(&s0,3) | age;


%best_phreg(imp_patients,surv_time,outcome,&vcns,&cvars,&ivars,phreg_steps);

proc sql;
select vars into: best_phreg_vars from phreg_steps having iauc=max(iauc);
quit;

proc phreg data=imp_patients plots=none;
   	model surv_time*outcome(&vcns)= &best_phreg_vars;
	output out=ttt_rg_best xbeta=phreg_best;
run;


proc sql;
select vars into: last_phreg_vars from phreg_steps having step=max(step);
quit;

proc phreg data=imp_patients plots=none;
   	model surv_time*outcome(&vcns)= &last_phreg_vars;
	output out=ttt_rg_last xbeta=phreg_last;
run;

proc sql;
create table ttt_reg as 
select tb2.*, tb1.phreg_last
from ttt_rg_last as tb1 join ttt_rg_best as tb2 on tb1.PATIENT_ID=tb2.PATIENT_ID
where tb1.chunk='test' and tb2.chunk='test';
quit;

proc sql;
create table ttt_final as 
select tb1.PATIENT_ID, tb1.outcome, tb1.surv_time, tb2.phreg_last, tb2.phreg_best
from tst_patients as tb1 join ttt_reg as tb2 on tb1.PATIENT_ID=tb2.PATIENT_ID;
quit;

%local all_tests i n tst;

%let all_tests =PETO WILCOXON LOGRANK TARONE; /* MODPETO;*/

%let n=%sysfunc(countw(&all_tests,%str( )));

%let rocs=;
%local bonf;
%let bonf=1;

%do i=1 %to &n;
	%let tst=%scan(&all_tests,&i,%str( ));

	%do bonf=1 %to 1;
	
		%build_forest(SURVLIB.&tst.&bonf._forest,trn_patients,outcome,&vcns,surv_time, &ivars, &cvars, &tst, 20, 15, 15, 15, 10, &bonf, &p_skip, &forest_sz,2,&main_dir); 
	
		%score_forest(SURVLIB.&tst.&bonf._forest,t_&tst.&bonf._forest, tst_patients,ttt_&tst.&bonf._forest,&forest_sz);

		%build_tree(SURVLIB.&tst.&bonf._tree,trn_patients,outcome,&vcns,surv_time, &ivars, &cvars, &tst, 20, 15, 15, 15, 10, &bonf, 0, &total_procs, &main_dir, 0); 
		
		%prune_tree(SURVLIB.&tst.&bonf._tree,vld_patients,trn_patients,outcome,&vcns,surv_time,info_prune,&tst);

		proc sql;
			select "if " || trim(t1.longcond) || " then do; _node_=" || put(t1.id,best.) || ";t=" || put(t1.t,best.) || ";xb=" || put(t1.xb,best.) || ";end;" 
			into : code_pruned separated by ';' 
			from SURVLIB.&tst.&bonf._tree_pruned as t1 
			where not exists select * from SURVLIB.&tst.&bonf._tree_pruned as t2 where t1.id=t2.prnt;
		quit;

		proc sql;
			select "if " || trim(t1.longcond) || " then do; _node_=" || put(t1.id,best.) || ";t=" || put(t1.t,best.) || ";xb=" || put(t1.xb,best.) || ";end;" 
			into : code_rules separated by ';' 
			from SURVLIB.&tst.&bonf._tree_rules as t1 
			where not exists select * from SURVLIB.&tst.&bonf._tree_rules as t2 where t1.id=t2.prnt;
		quit;

		data ttt_tree_&tst.;
		set tst_patients;
			&code_pruned;
			t_&tst.&bonf._pruned=xb;
			&code_rules;
			t_&tst.&bonf._rules=xb;
		run;

		data ttt_final;
		merge ttt_final ttt_tree_&tst.(keep=PATIENT_ID t_&tst.&bonf._pruned t_&tst.&bonf._rules);
		by PATIENT_ID;
		run;

		
		data ttt_final;
		merge ttt_final ttt_&tst.&bonf._forest(keep=PATIENT_ID t_&tst.&bonf._forest);
		by PATIENT_ID;
		run;

		
		%let rocs = roc %str(%')Pruned &tst. Tree Bonf=&bonf.%str(%') pred=t_&tst.&bonf._pruned %str(|) 
		roc %str(%')&tst. Tree Bonf=&bonf.%str(%') pred=t_&tst.&bonf._rules %str(|)
		roc %str(%')&tst. Forest%str(%') pred=t_&tst.&bonf._forest %str(|)
		&rocs;

   	%end;
%end;

%let rocs=%sysfunc(TRANWRD(&rocs,%str(|),%str(;)));
%put rocs="&rocs";

ods select all;
ods output iauc=iauc;
proc phreg data=ttt_final plots=roc rocoptions(IAUC at=1 to 35 by 1);
   model surv_time*outcome(&vcns)=phreg_last phreg_best / nofit;
   roc 'Pruned PH Reg' pred=phreg_best;
   roc 'PH Reg' pred=phreg_last;
   &rocs.
run;
%mend;

%macro simul();
%local i;
proc sql; 
create table SURVLIB.iaucs (ii int, est float, src varchar (1000));
quit;
%do i=1 %to 10;
	%main();
proc sql; 
insert into SURVLIB.iaucs select &i, Estimate, Source from iauc;
quit;
%end;
proc sgplot data=survlib.iaucs;
vbox est / group=src;
run;
%mend;

%macro build_forest(name, dts, cns, vcns, trgt, ivars, cvars, tst, nbins, mdp, mcnt, mpv, smooth, bonf, p_skip, sz, bag_sz, main_dir);
%local i tmp_tp pname;

%let pname = ft;

%do i=1 %to &sz;

	%put build tree N&i. in the forest &name;

	proc surveyselect data=&dts. method=urs out=samp_&i. reps=&bag_sz samprate=1;
	run;

	
	data samp_&i.;
	set samp_&i.;
	run;

	signon &pname._&i.;
		%syslput r_dts=samp_&i. /remote=&pname._&i.;
		%syslput r_cns=&cns /remote=&pname._&i.; 
		%syslput r_mdp=&mdp /remote=&pname._&i.; 
		%syslput r_vcns=&vcns /remote=&pname._&i.; 
		%syslput r_trgt=&trgt /remote=&pname._&i.; 
		%syslput r_cvars=&cvars /remote=&pname._&i.; 
		%syslput r_ivars=&ivars /remote=&pname._&i.; 
		%syslput r_i=&i /remote=&pname._&i.;
		%syslput r_tst=&tst /remote=&pname._&i.; 
		%syslput r_nbins=&nbins /remote=&pname._&i.; 
		%syslput r_mcnt=&mcnt /remote=&pname._&i.; 
		%syslput r_mpv=&mpv /remote=&pname._&i.; 
		%syslput r_smooth=&smooth /remote=&pname._&i.; 
		%syslput r_main_dir=&main_dir /remote=&pname._&i.;
		%syslput r_bonf=&bonf /remote=&pname._&i.;
		%syslput r_pskip = &p_skip /remote=&pname._&i.;

		rsubmit &pname._&i. wait=no; /* log=purge;*/

		proc upload;
		run;

		data _NULL_;
			call symput('pwork',getoption('work'));
		run; 
				
		proc upload
			infile= "&r_main_dir./check_var.sas"
			outfile="&pwork./zzz_&r_i..sas";
		run;

		
		proc upload
			infile= "&r_main_dir./build_tree.sas"
			outfile="&pwork./xxx_&r_i..sas";
		run;

		%include "&pwork./zzz_&r_i..sas";
		%include "&pwork./xxx_&r_i..sas";

		data _null_;
			call sleep(1,1);
		run;

		%build_tree(new_tree, &r_dts, &r_cns, &r_vcns, &r_trgt, &r_ivars, &r_cvars, &r_tst, &r_nbins, &r_mdp, &r_mcnt, &r_mpv, &r_smooth, &r_bonf, &r_pskip, 1, &r_main_dir, &r_i);

		data new_tree_f&r_i._rules;
		set new_tree_rules;
		run;

		proc download;
		run;

	endrsubmit;

	/*
	%build_tree(&name._f&i., samp_&i., &cns, &vcns, &trgt, &ivars, &cvars, &tst, &nbins, &mdp, &mcnt, &mpv, &smooth, &bonf, &p_skip);

	proc sql;
	drop table samp_&i.;
	quit;*/
%end;


%do i=1 %to &sz;
	waitfor &pname._&i.;
	signoff &pname._&i.;

	data &name._f&i._rules;
	set new_tree_f&i._rules;
	run;

	proc sql;
		drop table samp_&i.;
		drop new_tree_f&i.;
	quit;
%end;

%put forest &name completed;

%mend;

%macro score_forest(tree_name, var_name, dts_in, dts_out, sz);

%local i;

%put score forest &tree_name on &dts_in into &dts_out;

%do i=1 %to &sz;
	proc sql;
		select "if " || trim(t1.longcond) || " then do; _node_=" || put(t1.id,best.) || ";t=" || put(t1.t,best.) || ";xb=" || put(t1.xb,best.) || ";end;" 
		into : code_f&i._rules separated by ';' 
		from &tree_name._f&i._rules as t1 
		where not exists select * from &tree_name._f&i._rules as t2 where t1.id=t2.prnt;
	quit;
%end;

data &dts_out.;
set &dts_in.;
&var_name. = 0;
%do i=1 %to &sz;
	&&code_f&i._rules;
	&var_name=&var_name+xb;
%end;
&var_name=&var_name/&sz;
run;
%mend;


%simul();
/*%main();*/




