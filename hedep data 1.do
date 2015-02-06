*** Purpose: prepare data for analysis
*** Author: S Bauldry
*** Date: September 22, 2014

****** Extracting data from Wave 1
use aid bio_sex h1gi4 h1gi6a-h1gi6e h1gi11 pa12 pb8 h1rm1 h1rf1 h1ed11-h1ed14 ///
  h1ee1 h1ee2 pa55 h1gh1 ah_pvt h1to3 h1fs* h1gh59a h1gh59b h1gh60 h1to15 ///
  h1hs3 using "$w1", replace
  
*** setting missing values
recode bio_sex h1gi4 h1gi6a-h1gi6e (6 7 8 9 = .)
recode h1gi11 (7 = 1) (6 8 = .)
recode pa12 pb8 h1rm1 h1rf1 (11 12 96/99 = .)
recode h1ed11-h1ed14 (5 6 96/99 = .)
recode h1ee1 h1ee2 (6 8 9 = .)
recode pa55 (9996 = .)
recode h1gh1 (6 8 = .)
recode h1to3 (7 = 0) (6 8 = .)
recode h1fs* (6 8 9 = .)
recode h1gh59a h1gh59b (96 98 99 = .)
recode h1gh60 (996 998 999 = .)
recode h1to15 (97 = 7) (96 98 = .)
recode h1hs3 (6 8 = .)

*** preparing covariates
gen female = ( bio_sex == 2 ) if !mi(bio_sex)
lab var female "female"
drop bio_sex

gen race = 3 if h1gi4 == 1
replace race = 2 if h1gi6b == 1 & mi(race)
replace race = 4 if (h1gi6c == 1 | h1gi6d == 1 | h1gi6e == 1) & mi(race)
replace race = 1 if h1gi6a == 1 & mi(race)
lab def r 1 "white" 2 "black" 3 "hispanic" 4 "other"
lab val race r
lab var race "race"
qui tab race, gen(r)
drop h1gi4 h1gi6*

rename h1gi11 native
lab var native "born in US"

recode pa12 pb8 h1rm1 h1rf1 (10 = 1)
gen paredu = max(pa12,pb8)
gen cparedu = max(h1rm1,h1rf1)
replace paredu = cparedu if missing(paredu)
recode paredu (1 2 3 = 1) (4 5 = 2) (6 7 = 3) (8 = 4) (9 = 5)
lab var paredu "parent education"
qui tab paredu, gen(p)
drop pa12 pb8 h1rm1 h1rf1 cparedu

gen lninc = log(pa55 + 1)
lab var lninc "parent income (logged)"
drop pa55

recode h1ed11-h1ed14 (1 = 4) (2 = 3) (3 = 2) (4 = 1)
egen gpa = rowmean(h1ed11-h1ed14)
lab var gpa "gpa"
drop h1ed11-h1ed14

rename (h1ee1 h1ee2) (colasp colexp)
lab var colasp "college aspirations"
lab var colexp "college expectations"

rename (h1gh1 h1to3) (srh smoke)
lab var srh "self-rated health"
lab var smoke "ever smoked regularly"

gen drink = 7 - h1to15
lab var drink "drinking"
drop h1to15

gen bmi = (h1gh60*703)/( (12*h1gh59a + h1gh59b)^2 )
lab var bmi "BMI"
drop h1gh59a h1gh59b h1gh60

rename h1hs3 cns
lab var cns "received counseling in past 12 months"

*** saving temporary data file for merging
sort aid
tempfile d1
save `d1'


****** Extracting data from Wave 4
use aid h4ed2 imonth4 iday4 iyear4 h4od1m h4od1y h4ed6 h4mh18-h4mh27 ///
  h4hs9 using "$w4", replace
  
*** setting missing values
recode h4ed2 (96 98 = .)
recode h4ed6 (6 8 = .)
recode h4mh18-h4mh27 h4hs9 (6 8 = .)

*** preparing covariates
gen w4age = floor( ( mdy(imonth4, iday4, iyear4) - mdy(h4od1m, 15, h4od1y) )/364.25 )
lab var w4age "w4 age"
drop i* h4od*

recode h4ed2 (1/6 = 0) (7/13 = 1), gen(w4ba)
recode h4ed2 (1/5 = 0) (6/13 = 1), gen(w4sc)
lab var w4ba "w4 obtained college degree"
lab var w4sc "w4 some college"
drop h4ed2

rename h4ed6 w4incol
lab var w4incol "w4 in college"

rename h4hs9 w4cns
lab var w4cns "w4 received counseling in past 12 months"

*** saving temporary data file
sort aid
tempfile d2
save `d2', replace


****** Extracting family structure variable
use aid famst5 using "$fs", replace

*** creating indicators
lab var famst5 "family structure"
lab def fs 1 "2 bio parents" 2 "other 2 parent" 3 "single mother" 4 "single father" 5 "other"
lab val famst5 fs
qui tab famst5, gen(fs)

*** saving temporary file
sort aid
tempfile d3
save `d3', replace


****** Extracting weights and region
use aid gswgt4_2 region using "$wt", replace

*** creating region indicators
lab def rg 1 "west" 2 "midwest" 3 "south" 4 "northeast"
lab val region rg
lab var region "region"
qui tab region, gen(rg)

*** saving temporary file
sort aid
tempfile d4
save `d4', replace


****** Merging data
use `d1', replace
merge 1:1 aid using `d2'
drop _merge

merge 1:1 aid using `d3'
drop _merge

merge 1:1 aid using `d4'
drop _merge


*** selecting analysis sample
drop if mi(gswgt4_2)
dis _N

drop if mi(w4ba)
dis _N

drop if mi(h4mh19, h4mh22, h4mh24, h4mh26)
dis _N

drop if mi(ah_pvt, gpa, bmi, paredu, colasp, drink, smoke, cns, race, srh)
dis _N


*** constructing factor scores for adolescent and young adult depression
qui sem (Dep -> h1fs3 h1fs6 h1fs11 h1fs16)
predict depr, latent

qui sem (Dep -> h4mh19 h4mh22 h4mh24 h4mh26)
predict w4depr, latent

*** imputing family income
gen mlninc = ( mi(lninc) )
impute lninc native srh cns smoke-female r2-r4 paredu-bmi fs2-fs5 w4ba, gen(ilninc)
lab var mlninc "missing parent income"
lab var ilninc "imputed parent income (logged)"

*** Saving data for analysis
save "hedep data 1", replace
