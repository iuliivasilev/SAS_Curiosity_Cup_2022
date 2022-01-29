/*%include "&main_dir/check_var.sas";*/

%macro prune_tree(name,dts,dts_trn,cns,vcns,trgt,info_prune,tst);
%local code i pruned n s e try2prune emax smax step;

data temp_tree;
set &name._rules;
run;

proc sql;
	create table &name._pruned as 
	select * from &name._rules;
quit;

proc sql;
	create table &info_prune (step int, IAUC float, pruned varchar(1000));
quit;

%let step=0;
%let n=1;

%tree_xb_calc(temp_tree,&dts_trn,&cns,&vcns,&trgt);

proc sql;
	select "if " || trim(t1.longcond) || " then do; _node_=" || put(t1.id,best.) || ";t=" || put(t1.t,best.) || ";xb=" || put(t1.xb,best.) || ";tm=" || put(t1.tm,best.) || ";end;" into : code separated by ';' 
	from temp_tree as t1 where not exists 
	select * from temp_tree as t2 where t1.id=t2.prnt;
quit;

data ttt;
set &dts;
&code.;
/*t0=-tm;*/
run;

ods output IAUC=iauc;
proc phreg data=ttt plots=none rocoptions(IAUC);
   	model &trgt*&cns(&vcns)=xb / nofit;
   	/*roc pred=t0;*/
	roc pred=xb;
run;

proc sql;
	select Estimate into : e from iauc;
	drop table iauc;
quit;

proc sql;
	insert into &info_prune values(0,&e," ");
quit;

%do %while (&n);

	%let step=%sysevalf(&step+1);

	%let try2prune=;

	proc sql;
		select prnt, count(*) as cnt into : try2prune separated by '|' 
		from temp_tree as t1 where 
		id in (select id from temp_tree as t3 where 
		not exists select * from temp_tree as t2 
		where t3.id=t2.prnt) group by prnt
		having cnt=2;
	quit;

	%put try2prune=&try2prune;

	%let n=%sysfunc(countw(&try2prune,%str(|)));
	%let s =;
	%let emax=0;

	%do i=1 %to &n;
		%let s = %scan(&try2prune,&i,%str(|));
		proc sql;
			create table ttemp_tree as select * from temp_tree where prnt ne &s;
		quit;

		%tree_xb_calc(ttemp_tree,&dts_trn,&cns,&vcns,&trgt);

		proc sql;
			select "if " || trim(t1.longcond) || " then do; _node_=" || put(t1.id,best.) || ";t=" || put(t1.t,best.) || ";xb=" || put(t1.xb,best.) || ";tm=" || put(t1.tm,best.) || ";end;" 
			into : code_pruned separated by ';' 
			from ttemp_tree as t1 where not exists 
			select * from ttemp_tree as t2 where t1.id=t2.prnt;
		quit;

		data ttt;
			set &dts;
			&code_pruned;
			/*t0=-tm;*/
		run;

		ods exclude all;
		ods output IAUC=iauc;
		proc phreg data=ttt plots=auc rocoptions(IAUC);
   			model &trgt*&cns(&vcns)=xb / nofit;
   			/*roc pred=t0;*/
			roc pred=xb;
		run;

		%let e=0;

		proc sql;
			select Estimate into : e from iauc;
			drop table iauc;
		quit;
		%put s=&s e=&e;
		%if %sysevalf(&e>&emax) %then %do;
			%let emax=&e;
			%let smax=&s;
			%put emax=&emax smax=&smax;
		%end;
	%end;
	%if %sysevalf(&emax > 0) %then %do;
		proc sql;
			delete from temp_tree where prnt=&smax;
		quit;
		%let pruned = &pruned &smax;
		%put pruned=&pruned;
		proc sql;
			insert into &info_prune values(&step,&emax,"%sysfunc(TRANWRD(&pruned,%str( ),%str(,)))");
		quit;
	%end;
%end;

proc sql;
	select pruned into: pruned from &info_prune having iauc=max(iauc);
quit;

%put best prune=&pruned;
%if %length(&pruned)>0 %then %do;
proc sql;
	create table &name._pruned as 
	select * from &name._rules 
	where prnt not in (&pruned); 
quit;
%end;

%tree_xb_calc(&name._pruned,&dts_trn,&cns,&vcns,&trgt);

ods select all;
proc sgplot data=&info_prune;
series x=step y=iauc;
run;

proc sql;
	select "if " || trim(t1.longcond) || " then do; _node_=" || put(t1.id,best.) || "; end;" into : code separated by ';' 
	from &name._pruned as t1 where not exists 
	select * from &name._pruned as t2 where t1.id=t2.prnt;
quit;

data ttt;
	set &dts;
	&code;
run;

ods select all;
ods exclude ProductLimitEstimates Quartiles;
proc lifetest data=ttt outsurv=&name._pruned_surv
            plots=(survival (nocensor cb=hw test));
   time &trgt*&cns(&vcns);
   strata _node_ / test = &tst;
run;

%mend;

%macro tree_xb_calc(tree,dts,cns,vcns,tgt);
%put tree_xb_calc started for &tree and based on &dts;

proc sql;
	select "if " || trim(t1.longcond) || " then do; _node_=" || put(t1.id,best.) || ";t=" || put(t1.t,best.) || ";p=" || put(t1.p,best.) || ";end;" 
	into : code separated by ';' from &tree as t1 
	where not exists select * from &tree as t2 where t1.id=t2.prnt;
quit;

%put code="&code";

data ggg0;
set &dts;
&code.;
run;

ods exclude all;
ods output ParameterEstimates=params;
proc phreg data=ggg0 plots=none;
   class _node_;
   model &trgt*&cns(&vcns)=_node_;
   /*output out=ggg1  xbeta=xb;*/
run;
/*
proc sql;
	update &tree set xb=0;
	update &tree set xb=(select avg(xb) from ggg1 where 
	ggg1._node_=id group by ggg1._node_);
quit;*/

proc sql;
	update &tree set xb=0;
	update &tree set xb=(select Estimate from params
	where input(ClassVal0,best.)=id);
	update &tree set xb=0 where xb eq .;
quit;

%put tree_xb_calc completed;

%mend;

%macro build_tree(name, dts, cns, vcns, trgt, ivars, cvars, tst, nbins, mdp, mcnt, mpv,smooth, bonf, p_skip, total_procs, main_dir, tree_id);

proc sql;
create table &name._rules (id int, prnt int, cond varchar(1024), longcond varchar(4096), n int, p float, t float, tm float, xb float);
quit;

data tmp_data;
set &dts;
run;

%split_node(tmp_data,&name._rules,&cns,&vcns,&trgt, &ivars, &cvars, &tst, &nbins, &mcnt, &mpv,&smooth, &mdp,0,5000, &bonf, &p_skip, &total_procs, &main_dir, &tree_id); 

%tree_xb_calc(&name._rules,tmp_data,&cns,&vcns,&trgt);

proc sql;
drop table tmp_data;
quit;

proc sql;
select 
t1.id, trim(t1.longcond)
from &name._rules as t1 where not exists select * from &name._rules as t2 where t1.id=t2.prnt;
quit;

%mend;

%macro split_node(dts, tree, cns, vcns, trgt, ivars, cvars, tst, nbins, mcnt, mpv, smooth, mdp, prnt, smpsz, bonf, p_skip, total_procs, main_dir, tree_id);

%local i j n s maxcond maxid lid rid ln rn lp rp lt rt prntcond sn_conds pname;

%let sn_conds=sn_conds_&bonf._&prnt.;

proc sql;
create table &sn_conds (cond varchar(1000), Ch2 float, N int, ltm float, rtm float);
quit;

%let pname = cv_&tree_id.;

%put split_node start for parent id &prnt proc name &pname;

%local maxvars k p vrs ltm rtm;

%let p=0;

%do j=0 %to 1;
	%if &j %then %let vars = &ivars; %else %let vars = &cvars; 
	%let n=%sysfunc(countw(&vars,%str(|)));
	%put n=&n vars=&vars;
	%let maxvars = 0;
	%if %sysevalf(&total_procs>1) %then %do;
		%let maxvars = %sysevalf(&n/(&total_procs-1));
		%let maxvars = %sysfunc(floor(&maxvars));
	%end;
	%if not &maxvars %then %let maxvars = &n;
	%put maxvars=&maxvars;
	%let k=0;
	%let vrs=;
	%do i=1 %to &n;
		%let s=%scan(&vars,&i,%str(|));
		%let k=%sysevalf(&k+1);
		%let vrs=&vrs &s;
		%if (&k eq &maxvars) or (&i eq &n) %then %do;

			%let p=%sysevalf(&p.+1);

			%if %sysevalf(&total_procs>1) %then %do;

				%put Starting proc &p. for &vrs.;

				%put temp table=&sn_conds._&p.;
				signon &pname._&p.;
				%syslput r_dts=dts_&p. /remote=&pname._&p.;
				%syslput r_sn_conds=&sn_conds._&p. /remote=&pname._&p.; 
				%syslput r_cns=&cns /remote=&pname._&p.; 
				%syslput r_vcns=&vcns /remote=&pname._&p.; 
				%syslput r_trgt=&trgt /remote=&pname._&p.; 
				%syslput r_vrs=&vrs /remote=&pname._&p.; 
				%syslput r_j=&j /remote=&pname._&p.; 
				%syslput r_i=&i /remote=&pname._&p.;
				%syslput r_p=&p /remote=&pname._&p.; 
				%syslput r_tst=&tst /remote=&pname._&p.; 
				%syslput r_nbins=&nbins /remote=&pname._&p.; 
				%syslput r_mcnt=&mcnt /remote=&pname._&p.; 
				%syslput r_mpv=&mpv /remote=&pname._&p.; 
				%syslput r_smooth=&smooth /remote=&pname._&p.; 
				%syslput r_main_dir=&main_dir /remote=&pname._&p.;
				%syslput r_retcode=retcode_&p. /remote=&pname._&p.;
				%syslput r_pskip = &p_skip /remote=&pname._&p.;

				data dts_&p.;
				set &dts;
				run;

				rsubmit &pname._&p. wait=no log=purge;

				proc upload;
				run;

				data _NULL_;
					call symput('pwork',getoption('work'));
				run; 
			
				proc upload
					infile= "&r_main_dir./check_var.sas"
					outfile="&pwork./zzz_&r_p..sas";
				run;

				%include "&pwork./zzz_&r_p..sas";

				proc sql;
					create table conds (cond varchar(1000), Ch2 float, N int, ltm float, rtm float);
				quit;

				data _null_;
					call sleep(1,1);
				run;
			
				%check_vars(&r_dts, conds,  &r_cns, &r_vcns, &r_trgt, &r_vrs, &r_j, &r_tst, &r_nbins, &r_mcnt, &r_mpv, &r_smooth, &r_pskip); 
	
				data &r_sn_conds;
					set conds;
				run;

				proc download;
				run;

				endrsubmit;
			%end;%else %do;

				proc sql;
					create table &sn_conds._&p. (cond varchar(1000), Ch2 float, N int, ltm float, rtm float);
				quit;

				%check_vars(&dts,&sn_conds._&p.,&cns,&vcns,&trgt,&vrs,&j,&tst,&nbins,&mcnt,&mpv,&smooth,&p_skip); 

			%end;

			%let k=0;
			%let vrs=;

		%end;
	%end;
%end;

%do i=1 %to &p;
	%if %sysevalf(&total_procs>1) %then %do;
		waitfor &pname._&i.;
		signoff &pname._&i.;
		proc sql;
			drop table dts_&i.;
		quit;
	%end;
	proc sql;
		insert into &sn_conds select * from &sn_conds._&i.;
		drop table &sn_conds._&i.;
	quit;
%end;

%let rtm=.;
%let ltm=.;

data _NULL_;
set &sn_conds end=final;
length mcond $ 100;
retain mcond;
retain mworth 0;
pval = 1-CDF('CHISQUARE',ch2,1);
worth = .;
/*if n*pval then worth = -Log(n*pval); */
%if &bonf %then %do;
if pval>0 then worth = -Log(n*pval); 
%end; %else %do; 
if pval>0 then worth = -Log(pval);
%end;
if (_N_= 1) or (worth > mworth) or (pval = 0) or (worth eq .) 
then do; 
	mcond = cond;
	mworth = worth;
end;
put pval= n= worth= mworth= cond= ch2=;
if final or (worth eq .) then do;
	call symput('maxcond',mcond);
	call symput('ltm',ltm);
	call symput('rtm',rtm);
	stop;
end;
run;

proc sql;
	drop table &sn_conds;
quit;

%let mdp = %eval(&mdp-1);

%if (&mdp) and (%length(%qtrim(&maxcond))>0) %then %do;

	%put mdp=&mdp maxcond=&maxcond;

	proc sql;
		select max(id) into : maxid from &tree;
	quit;

	proc sql;
		select longcond into : longcond from &tree where id=&prnt;
	quit;

	%put maxid=&maxid longcond=&longcond;

	%if %sysfunc(notdigit(%qtrim(&maxid))) %then %do;
		%let maxid=0;
		%let longcond=(1=1);
	%end;

	%let longcond=%qtrim(&longcond); 

	%let lid = %eval(&maxid+1);
	%let rid = %eval(&maxid+2);

	%put lid=&lid rid=&rid maxid=&maxid;
	%put longcond=&longcond;

	data 
	ldts_&prnt (drop= ln rn lp rp lt rt) 
	rdts_&prnt (drop= ln rn lp rp lt rt);
	retain ln rn lp rp lt rt;
	set &dts end=final;
	if (&maxcond) then do;
		output ldts_&prnt;
		ln+1;
		lt+&trgt;
		if (&cns = &vcns) then lp+1;
	end; else do;
		output rdts_&prnt;
		rn+1;
		rt+&trgt;
		if (&cns = &vcns) then rp+1;
	end;
	if final then do;
		call symput('ln',ln);
		call symput('rn',rn);
		lt = lt / ln;
		rt = rt / rn;
		lp = lp / ln;
		rp = rp / rn;
		call symput('lt',lt);
		call symput('rt',rt);
		call symput('lp',lp);
		call symput('rp',rp);
	end;
	run;

	%put compare predictions ltm=&ltm rtm=&rtm lt=&lt rt=&rt;
	
	%if %sysevalf((&lt-&rt)*(&ltm-&rtm)<0) 
		%then %put Check it twice!!!;

	%if &ltm eq . %then %do;
		%let ltm=&lt;
	%end;

	%if &rtm eq . %then %do;
		%let rtm=&rt;
	%end;

	%put insert into &tree values(&lid, &prnt, "%qtrim(&maxcond)", "((%qtrim(&longcond)) and (%qtrim(&maxcond))", &ln, &lp, &lt,&ltm,0);
	%put insert into &tree values(&rid, &prnt, "not (%qtrim(&maxcond))", "(%qtrim(&longcond)) and (not (%qtrim(&maxcond)))", &rn, &rp, &rt,&rtm,0);  
	
	proc sql;
		insert into &tree values(&lid, &prnt, "%qtrim(&maxcond)", "(%qtrim(&longcond)) and (%qtrim(&maxcond))", &ln, &lp, &lt,&ltm,0);
		insert into &tree values(&rid, &prnt, "not (%qtrim(&maxcond))", "(%qtrim(&longcond)) and (not (%qtrim(&maxcond)))", &rn, &rp, &rt,&rtm,0); 
	quit;

	%put split_node stop for parent id &prnt cond found "&maxcond";

	%split_node(ldts_&prnt, &tree, &cns, &vcns, &trgt, &ivars, &cvars, &tst, &nbins, &mcnt, &mpv, &smooth, &mdp, &lid, &smpsz, &bonf, &p_skip, &total_procs, &main_dir, &tree_id);

	proc sql;
		drop table ldts_&prnt;
	quit;

	%split_node(rdts_&prnt, &tree, &cns, &vcns, &trgt, &ivars, &cvars, &tst, &nbins, &mcnt, &mpv, &smooth, &mdp, &rid, &smpsz, &bonf, &p_skip, &total_procs, &main_dir, &tree_id);

	proc sql;
		drop table rdts_&prnt;
	quit;
%end;
%mend;