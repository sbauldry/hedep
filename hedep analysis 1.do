*** Purpose: produce main hedep analysis
*** Author: S Bauldry
*** Date: September 22, 2014

use "hedep data 1", replace

*** estimating propensity scores and identifying strata
qui probit w4sc w4age female r2-r4 native fs1 p1-p4 ilninc mlninc gpa ah_pvt ///
  colasp srh smoke drink bmi depr cns [pw = gswgt4_2]
eststo m1
qui probit w4ba w4age female r2-r4 native fs1 p1-p4 ilninc mlninc gpa ah_pvt ///
  colasp srh smoke drink bmi depr cns [pw = gswgt4_2]
eststo m2
esttab m1 m2 using psms.csv, replace b(%9.3f) ci(%9.3f) pr2(%9.3f) star nogaps wide

pscore w4sc w4age female r2-r4 native fs1 p1-p4 ilninc mlninc gpa ah_pvt ///
  colasp srh smoke drink bmi depr cns [pw = gswgt4_2], comsup level(0.001) ///
  pscore(psc) blockid(bsc)

pscore w4ba w4age female r2-r4 native fs1 p1-p4 ilninc mlninc gpa ah_pvt ///
  colasp srh smoke drink bmi depr cns [pw = gswgt4_2], comsup level(0.001) ///
  pscore(pba) blockid(bba)

*** estimating measurement models within each strata
capture program drop MM
program MM
	args b
	forval i = 1/14 {
		qui sem (Dep -> h4mh22 h4mh19 h4mh24 h4mh26) if `b' == `i' 
		local bic = e(chi2_ms) - e(df_ms)*ln( e(N) )
		dis "`i': " as res %5.2f e(chi2_ms) " " as res %5.2f e(df_ms) " " ///
	     as res %5.3f e(p_ms) " " as res %5.2f `bic'
	}
end

MM bsc
MM bba		

*** estimating average treatment effects
sem (Depress -> h4mh22 h4mh19 h4mh24 h4mh26) (w4sc psc -> Depress) ///
  [pw = gswgt4_2] if !mi(bsc)
  
sem (Depress -> h4mh22 h4mh19 h4mh24 h4mh26) (w4ba pba -> Depress) ///
  [pw = gswgt4_2] if !mi(bba)
  
*** program for stratification-multilevel estimator
capture program drop wste
program wste
	args iv out lb ub tit gn 
	preserve
	
	postutil clear
	tempfile wste
	postfile wste strata b v using `wste', replace
	
	forval i = 1/14 {
		qui sem (Depress -> h4mh22 h4mh19 h4mh24 h4mh26) ///
		    (w4`iv' -> Depress) [pw = gswgt4_2] if b`iv' == `i'	
		mat b = e(b)
		mat v = e(V)
		local b = b[1,1]
		local v = v[1,1]
		post wste (`i') (`b') (`v')
	}
	postclose wste
	
	use `wste', replace
	
	qui gen se = sqrt(v)
	vwls b strata, sd(se)
	
	qui gen ub = b + 1.96*se
	qui gen lb = b - 1.96*se
	qui gen pval = 2*( 1 - normal( abs(b/se) ) )
	qui gen star = "***" if pval < 0.001
	qui replace star = "**" if pval < 0.01 & mi(star)
	qui replace star = "*" if pval < 0.05 & mi(star)
	qui replace star = "+" if pval < 0.10 & mi(star)
	format b lb ub %5.2f
	list strata b star lb ub, clean
	
	format b %5.1f	
	twoway (scatter b strata if strata >= `out') ///
	  (rcap ub lb strata if strata >= `out', lc(black)) ///
	  (lfit b strata if strata >= `out', lp(dash)), scheme(s2mono) yline(0) legend(off) ///
	  xlab(1(1)14) ylab(`lb'(0.2)`ub') xtit("propensity score strata") ///
	  ytit("protective effect") title("`tit'") saving(`gn', replace)
end

*** program for matching-smoothing estimator
capture program drop mse
program mse
	args iv out lb ub tit gn
	preserve
	
	psmatch2 w4`iv', out(w4depr) pscore(p`iv') neighbor(5) caliper(3) common
	gen DepD = w4depr - _w4depr
	format DepD %5.2f
	lpoly DepD p`iv' if w4`iv' == 1, kernel(epanechnikov) bwidth(0.2) ///
	  degree(1) noscatter ci pwidth(0.3) scheme(s2mono) ylab(`lb'(0.1)`ub') ///
	  ytit("protective effect") legend(off) xtit("propensity score") note("") ///
      title("`tit'") saving(`gn', replace) yline(0)
end

*** obtaining estimates and preparing figure
wste sc 2 -0.4 0.2 "A: some college" g1
mse sc 2 -0.3 0.1 "B: some college" g3
wste ba 1 -0.4 0.2 "C: 4-year degree" g2
mse ba 1 -0.3 0.1 "D: 4-year degree" g4  
   
****** Figure 1 ******
graph combine g1.gph g2.gph g3.gph g4.gph, scheme(s2mono)
graph export Fig1.pdf, replace




