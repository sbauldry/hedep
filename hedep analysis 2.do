*** Purpose: produce auxiliary analysis requested by reviewer 3
*** Author: S Bauldry
*** Date: September 22, 2014

*** Loading data
use "hedep data 1", replace

*** constructing w4 depression scale using all items
alpha h4mh*, gen(w4depr2)

pscore w4sc w4age female r2-r4 native fs1 p1-p4 ilninc mlninc gpa ah_pvt ///
  colasp srh smoke drink bmi depr cns [pw = gswgt4_2], comsup level(0.001) ///
  pscore(psc) blockid(bsc)

pscore w4ba w4age female r2-r4 native fs1 p1-p4 ilninc mlninc gpa ah_pvt ///
  colasp srh smoke drink bmi depr cns [pw = gswgt4_2], comsup level(0.001) ///
  pscore(pba) blockid(bba)

*** estimating average treatment effects
regress w4depr w4sc psc [pw = gswgt4_2] if !mi(bsc)
regress w4depr2 w4sc psc [pw = gswgt4_2] if !mi(bsc)

regress w4depr w4ba pba [pw = gswgt4_2] if !mi(bba)
regress w4depr2 w4ba pba [pw = gswgt4_2] if !mi(bba)

*** program for stratification-multilevel estimator
capture program drop wste
program wste
	args dv iv out lb ub tit gn 
	preserve
	
	postutil clear
	tempfile wste
	postfile wste strata b v using `wste', replace
	
	forval i = 1/14 {
		qui regress `dv' w4`iv' [pw = gswgt4_2] if b`iv' == `i'
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

wste w4depr sc 2 -0.4 0.2 "A: 4-item, some college" ag1
wste w4depr2 sc 2 -0.4 0.2 "B: 10-item, some college" ag3
wste w4depr ba 1 -0.4 0.2 "C: 4-item, 4-year degree" ag2
wste w4depr2 ba 1 -0.4 0.2 "D: 10-item, 4-year degree" ag4

****** Figure 1 ******
graph combine ag1.gph ag2.gph ag3.gph ag4.gph, scheme(s2mono)
graph export aFig1.pdf, replace




