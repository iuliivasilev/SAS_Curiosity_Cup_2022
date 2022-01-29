/*options  NOSYNTAXCHECK; ERRORABEND; */
%global lifetest_fail;

%macro check_var(dts, conds, cns, vcns, trgt, vr, ivr, tst, nbins, mcnt, mpv, smooth);
%local i j n mn rks besti betsj best vl bestch2 vr0 pval ch2 vl mncnt cond mxvr mnvr props tm besttm;
%let n = 0;
%let rks = 0;
%let besti =-1;
%let bestj =-1;
%let bestvl =-1;
%let besttm =-1;
%let bestch2 = 0;

%put check_var start for &vr;

%put check_var dataset &dts;


%let vr0=&vr;
data cv_ttt0 (keep = &cns &trgt &vr0 strt);
set &dts;
strt = .;
run;

%if &SYSERR %then %do;
	%let lifetest_fail=1;
	options nosyntaxcheck obs=max replace;
	proc sql exec; quit;
	%put set lifetest_fail and go to next try ...;
	%goto exit;
%end;
%let lifetest_fail=0;
 
%if &ivr eq 0 %then %do;

	proc sql;
		select mean(&cns) into :props from cv_ttt0;
	quit;

	%if %sysevalf(&props eq 0) %then %let props=0.00001;

	%if %sysevalf(&props eq 1) %then %let props=0.99999;

	%put props=&props;

	ods exclude all;
	proc means data=cv_ttt0;
  		class &vr / missing;
  		var &cns;
		output out=ttt_stats0(rename=(_freq_=BIN_N) where=(_type_=1)) sum=BIN_T;
	run;
	%if &SYSERR %then %do;
	%let lifetest_fail=1;
	options nosyntaxcheck obs=max replace;
	proc sql exec; quit;
	%put set lifetest_fail and go to next try ...;
	%goto exit;
	%end;
	%let lifetest_fail=0;

	data ttt_stats1 (keep=&vr woe_&vr);
	set ttt_stats0;
	woe_&vr=log((BIN_T+&props*&smooth)/(BIN_N-BIN_T+(1-&props)*&smooth));
  	run;

	proc sql;
	create table cv_tttm1 as select &cns, &trgt, ttt_stats1.&vr, woe_&vr
		from cv_ttt0 join ttt_stats1 on ttt_stats1.&vr=cv_ttt0.&vr;
	quit;
  
	%let vr0=woe_&vr;
	
	data cv_ttt0;
	set cv_tttm1;
	run;

%end;

ods exclude all;
proc rank data=cv_ttt0 out=cv_ttt1 groups=&nbins;
var &vr0;
ranks rk;
run;

%if &SYSERR %then %do;
	%let lifetest_fail=1;
	options nosyntaxcheck obs=max replace;
	proc sql exec; quit;
	%put set lifetest_fail and go to next try ...;
	%goto exit;
%end;
%let lifetest_fail=0;

proc sql;
select distinct rk into : rks separated by '|' from cv_ttt1 
order by rk;
quit;

%let n=%sysfunc(countw(&rks,%str(|)));

%local has_miss;
%let has_miss = %sysevalf(%index(&rks,.)>0); 
%put rks=&rks;
%put n=&n has_miss=&has_miss idx=%index(&rks,.);

%do i=2 %to %eval(&n);
	
	%let vl = %scan(&rks,&i,%str(|));

	%put vl= &vl;

	/*j=0 => missings in false cond*/
	%do j=0 %to &has_miss;

		%let mn=0;
		data cv_ttt2;
		retain ln 0 rn 0;
		set cv_ttt1 end=final;
		strt = (rk<&vl);
		if (rk eq .) then strt = &j;
		if strt then ln+1; else rn+1;
		if final then do;
			call symput('mn',min(ln,rn));
			put ln= rn=;
		end;
		run;
		%if &SYSERR %then %do;
			%let lifetest_fail=1;
			options nosyntaxcheck obs=max replace;
			proc sql exec; quit;
			%put set lifetest_fail and go to next try ...;
			%goto exit;
		%end;
		%let lifetest_fail=0;

		%put mn=&mn;

		%if %sysevalf(&mn>&mcnt) %then %do;
			ods exclude all;
			ods output  
			CensoredSummary=cv_ttt4_&i._&j.
			Means=cv_ttt5_&i._&j.
			Homtests=cv_ttt3_&i._&j.;
			proc lifetest data=cv_ttt2 plots=none;
				time &trgt*&cns(&vcns);
				strata strt / test = &tst;
			run;
			%put PROCESS ERROR CODE=&SYSERR;
			%if &SYSERR %then %do;
			%let lifetest_fail=1;
			options nosyntaxcheck obs=max replace;
			proc sql exec; quit;
			%put set lifetest_fail and go to next try ...;
			%goto exit;
			%end;
			%let lifetest_fail=0;

			%let mncnt=0;
			%let tm=.,.;
			proc sql;
				select min(Total) into :mncnt from cv_ttt4_&i._&j.;
				select Mean into: tm separated by ',' from cv_ttt5_&i._&j. order by strt desc;
				drop table cv_ttt4_&i._&j.;
				drop table cv_ttt5_&i._&j.;
			quit;
	

			%let pval=0;
			%let ch2=0;
			proc sql;
				select ChiSq, ProbChiSq into :ch2, :pval from cv_ttt3_&i._&j.;
				drop table cv_ttt3_&i._&j.;
			quit;


			%if "<.0001" = "&pval" %then %let pval=0.0;
			%if "." = "%qtrim(&pval)" %then %let pval=1.0;
			%put mncnt=&mncnt PVAL=&pval ch2=&ch2 tm=&tm has_miss=&has_miss;

			
			
			%put (&pval < %sysevalf(&mpv/100.0)) and (&mncnt > &mcnt) and (&ch2 > &bestch2);

			%if (%sysevalf(&mpv/100.0>&pval)) and (%sysevalf(&mncnt > &mcnt)) and (%sysevalf(&ch2 > &bestch2)) %then %do;
				%put TRUE!!!;
				%let bestch2 = &ch2;
				%let besti = &i;
				%let bestj = &j;
				%let bestvl = &vl;
				%let besttm=&tm;
				%put bestch2=&bestch2 bestj=&bestj besti=&besti bestvl=&bestvl besttm=&besttm has_miss=&has_miss;
			%end; %else %do; 
				%put FALSE!!!; 
			%end;
		%end;
	%end;
%end; 

%if %sysevalf(&besti>=0) %then
%do;
	%put final bestch2=&bestch2 bestj=&bestj besti=&besti bestvl=&bestvl besttm=&besttm has_miss=&has_miss;

	%let mxvr=;
	%let mnvr=;

	proc sql;
		select min(&vr0) into :mxvr from cv_ttt1 where (&vr0 ne .) and (rk=&bestvl);
		select max(&vr0) into :mnvr from cv_ttt1 where (&vr0 ne .) and &vr0<(&mxvr-0.00001);
	quit;


	%put mxvr=&mxvr mnvr=&mnvr;

	%let c_p=.;
	%let c_p=%sysevalf((&mxvr+&mnvr)/2.0);
	%put c_p=&c_p mxvr=&mxvr mnvr=&mnvr;

	%if %qtrim(&c_p) ne %str(.) %then %do;
		%let cond = (&vr0<&c_p);
		%if &bestj %then %let cond = ((&vr0 eq .) or &cond);
		%else %let cond = ((&vr0 ne .) and &cond);
	/*%end; %else %let cond = (&vr0 eq .);*/

		%if not &ivr %then %do;

			proc sql;
				select distinct &vr into: c_values separated by '|' from cv_ttt0 where &cond;
			quit;

			%put before c_values=&c_values;

			%let c_values =%str(%')%sysfunc(tranwrd(&c_values,%str(|),%str(',')))%str(%');

			%put after c_values=&c_values;

			%let cond = (&vr in (&c_values));
		%end;

		%put cond=&cond;

		proc sql;
			insert into &conds values ("&cond", &bestch2, &n, &besttm);
		quit;
	%end;
%end;

%put check_var stop for &vr cond inserted "&cond";
%exit:
%mend;

%macro check_vars(dts, conds, cns, vcns, trgt, vrs, ivr, tst, nbins, mcnt, mpv, smooth, p_skip);
%local n i s t;
%let n=%sysfunc(countw(&vrs,%str( )));
%do i=1 %to &n;
	%if %sysevalf(%sysfunc(ranuni(0))>&p_skip) %then %do;
		%let s = %scan(&vrs,&i,%str( ));
		%let t=3;
		%do %while (&t);
			%put left tries &t to run check_var for var=&s;
			%check_var(&dts, &conds, &cns, &vcns, &trgt, &s, &ivr, &tst, &nbins, &mcnt, &mpv, &smooth);
			%if (&lifetest_fail) %then %do;
				%put lifetest_fail is true;
				%let t=%sysevalf(&t-1);
				data _null_;
					call sleep(1,1);
				run;
			%end; %else %do; 
				%put lifetest_fail is false;
				%let t=0; 
			%end;
		%end;
	%end;
%end;
%mend;



