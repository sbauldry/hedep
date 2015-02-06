*** Purpose: produce hedep descriptives
*** Author: S Bauldry
*** Date: September 22, 2014

use "hedep data 1", replace
svyset [pw = gswgt4_2]

*** Descriptive statistics for text
sum lninc ilninc

tab w4sc
svy: prop w4sc

tab w4ba
svy: prop w4ba

tab w4incol w4ba

*** Descriptive statistics for Table 1
foreach x of varlist h4mh19 h4mh22 h4mh24 h4mh26 w4depr w4age female r1-r4 ///
  native fs1 p1-p5 ilninc mlninc gpa ah_pvt colasp srh bmi smoke drink depr cns {
	qui svy: mean `x'
	mat b1 = e(b)
	local b1 = b1[1,1]

	qui svy: mean `x', over(w4sc)
	mat b2 = e(b)
	qui test [`x']0 - [`x']1 = 0
	local b20 = b2[1,1]
	local b21 = b2[1,2]
	local p2 = r(p)

	qui svy: mean `x', over(w4ba)
	mat b3 = e(b)
	qui test [`x']0 - [`x']1 = 0
	local b30 = b3[1,1]
	local b31 = b3[1,2]
	local p3 = r(p)

	dis "`x': " as res %5.3f `b1' " " as res %5.3f `b20' " " as res %5.3f `b21' ///
      " " as res %5.3f `p2' " " as res %5.3f `b30' " " as res %5.3f `b31' ///
      " " as res %5.3f `p3'
}














