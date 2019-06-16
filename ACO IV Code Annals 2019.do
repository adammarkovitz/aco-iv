*** Stata Program 
*** used for "Performance in the Medicare Shared Savings Program After Accounting for Nonrandom Exit" 
*** by Adam A. Markovitz, Annals of Internal Medicine, June 2019

global iv "E:/amarkov"
global rif "$iv/rif"
global npi "$iv/npi"
global bene "$iv/bene"
global medpar "$iv/medpar"
global inst "$iv/inst"
global lat "$iv/inst/lat"
global out "C:/Users/amarkov/Google Drive/thesis/iv/out_ghost"
global mat "$out/mat"
global log "$out/log"
global fig "$out/figure"
capture log close

********************************************************************************
**								RHS macros									  **
********************************************************************************
 
global x age female i.race dual hcc dib esrd poverty hs college i.season i.year	// adjusted longitudinal model
global xm age i.race dual hcc dib esrd poverty hs college i.season i.year	// adjusted longitudinal model without female sex
global z age female black hispanic race_other dual hcc dib esrd poverty hs college season2 season3 season4 year2009 year2010 year2011 year2012 year2013 year2014	// instrumental variable model
global zm age black hispanic race_other dual hcc dib esrd poverty hs college season2 season3 season4 year2009 year2010 year2011 year2012 year2013 year2014	// instrumental variable model without female sex
global h age female i.race dual hcc dib esrd poverty hs college i.season	// market-year fixed effect model
global c_x c_age c_female c_black c_hispanic c_race_other c_hcc c_dual c_dib c_esrd c_poverty c_hs c_college c_year2009 c_year2010 c_year2011 c_year2012 c_year2013 c_year2014	// ACO + market-year fixed effect model
global b age dual hcc poverty hs college i.season i.year	// beneficiary fixed effect model

********************************************************************************
**							basic descriptives								  **
********************************************************************************

use "$rif/descriptive_100pct", replace

* outcomes
sum pay hip
sum pay, d

* MSSP participation
sum inACO tr_bene diff50
sum inACO tr_bene diff50 if qtr >= 18 

** unique participants
count
foreach x in id npi tin {
	unique `x'
}

** unique MSSP participants
keep if inACO == 1 
foreach x in id npi tin ACO_NUM {
	unique `x'
}

********************************************************************************
*		Table 1: Characteristics of beneficiaries and clinicians across 	   *
* 				MSSP ACO participation and MSSP supply						   *
********************************************************************************

use "$rif/descriptive_100pct", clear

** keep pre-period data for beneficiaries who appear in post-period (i.e., with calculated MSSP supply)
keep if qtr < 18 & t50 != .
count
unique id
unique npi

** covariate balance
loc x age female white black hispanic race_other dual hcc dib esrd poverty hs college ma pay hip
covbal tr_bene `x'
mat a = r(table)
mat a = a[1..16,4], a[1..16,1], a[1..16,7]
tab tr_bene, matcell(c)
mat c = (c[1,1], c[2,1], 0)
mat a = a \ c
mat li a

covbal t50 `x'
mat b = r(table)
mat b = b[1..16,4], b[1..16,1], b[1..16,7]
tab t50, matcell(c)
mat c = (c[1,1], c[2,1], 0)
mat b = b \ c
mat li b

mat c = a , b
mat li c

putexcel set "$mat/table1_100pct", replace
putexcel a1 = mat(c)

********************************************************************************
* Figure 1: Change in spending in adjusted longitudinal model and IV model     *
********************************************************************************

**********************************************
**		Adjusted longitudinal model			**
**********************************************

mat m = J(17,7,0)
loc n = 3

foreach y in pay ip op ptb snf  {

	**********************************************
	**		Adjusted longitudinal model			**
	**********************************************
	use "$rif/fe_100pct", clear
	loc x age female race dual hcc dib esrd poverty hs college season year
	keep `y' `x' inACO hrr_npi
	areg `y' inACO $x, absorb(hrr_npi) vce(cluster hrr_npi)
	mat a = _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	margins, at(inACO=(0 1))
	mat b = r(table)
	sca b = b[1,1]
	mat m[`n',1] = b, b + _b[inACO], a
	loc n = `n' + 1
	
	**********************************************
	**					IV model				**
	**********************************************
	use "$rif/iv_100pct", clear
	keep `y' $z diff50 inACO hrr_npi
	xtset hrr_npi
	xtivreg2 `y' $z (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin
	mat m[`n',1] = b, b + _b[inACO], _b[inACO], _b[inACO] - invnormal(0.975)*_se[inACO], _b[inACO] + invnormal(0.975)*_se[inACO], 2*normal(-abs( _b[inACO]/_se[inACO])), e(N)	
	loc n = `n' + 2
	putexcel set "$mat/figure1_100pct", replace
	putexcel a1 = mat(m)
	mat li m
}

******************************************************
**						F-test						**
******************************************************

use "$rif/iv_robust_100pct", clear
keep pay hip $z diff50 inACO hrr_npi
xtset hrr_npi
loc y pay
xtivreg2 `y' $z (inACO = diff50), cluster(hrr_npi) fe robust nocollin

********************************************************************************
**		Figure 2: Changes in spending and hospitalization for hip fracture	  **
********************************************************************************

foreach y in hip pay { 

	mat m = J(8,7,0)
	loc n = 3
	
	**	Adjusted longitudinal
	use "$rif/robust_100pct", clear
	loc x age dual hcc poverty hs college season year female race dib esrd 
	keep pay hip `x' inACO hrr_npi npi id hrr_year ma
	areg `y' inACO $x, absorb(hrr_npi) vce(cluster hrr_npi)
	mat a = _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	margins, at(inACO=(0 1))
	mat b = r(table)
	sca b = b[1,1]
	mat m[`n',1] = b, b + _b[inACO], a
	loc n = `n' + 1
	
	**	Market-year FE
	areg `y' inACO $h, absorb(hrr_year) vce(cluster hrr_npi)
	mat m[`n',1] = b, b + _b[inACO], _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	loc n = `n' + 1
	mat li m
	drop hrr_year

	**	Beneficiary FE
	areg `y' inACO $b, absorb(id) vce(cluster hrr_npi)
	mat m[`n',1] = b, b + _b[inACO], _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	loc n = `n' + 1
	mat li m
	drop id
	
	**  Market-year + ACO FE
	use "$rif/aco_demean_100pct", clear
	keep c_hip c_pay c_inACO $c_x hrr_npi hrr_year
	areg c_`y' c_inACO $c_x, absorb(hrr_year) vce(cluster hrr_npi)
	mat m[`n',1] = b, b + _b[c_inACO], _b[c_inACO], _b[c_inACO] - invttail(e(df_r),0.025)*_se[c_inACO], _b[c_inACO] + invttail(e(df_r),0.025)*_se[c_inACO], 2*ttail(e(df_r),abs(_b[c_inACO]/_se[c_inACO])), e(N)
	loc n = `n' + 1
	mat li m
	
	**	Clinician FE	
	use "$rif/fe_100pct", clear
	loc x age dual hcc poverty hs college season year female race dib esrd 
	keep hip pay `x' inACO hrr_npi npi
	areg `y' inACO $x, absorb(npi) vce(cluster hrr_npi)
	mat m[`n',1] =  b, b + _b[inACO], _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	loc n = `n' + 1
	mat li m

	* Instrumental variable
	use "$rif/iv_100pct", clear
	keep hip pay $z inACO hrr_npi diff50
	xtset hrr_npi
	xtivreg2 `y' $z (inACO = diff50), cluster(hrr_npi) fe robust endog(inACO) nocollin noid
	mat m[`n',1] = b, b + _b[inACO], _b[inACO], _b[inACO] - invnormal(0.975)*_se[inACO], _b[inACO] + invnormal(0.975)*_se[inACO], 2*normal(-abs( _b[inACO]/_se[inACO])), e(N)
	loc n = `n' + 1
	mat li m
	putexcel set "$mat/fig2_`y'_100pct", replace

}

********************************************************************************
* Figure 3: Change in quality in adjusted longitudinal model and IV model      *
********************************************************************************

mat m = J(17,8,0)
loc n = 3
foreach y in proc1 proc2 proc3 proc4 {

	**********************************************
	**		Adjusted longitudinal model			**
	**********************************************
	use "$rif/quality_100pct", clear
	loc x age female race dual hcc dib esrd poverty hs college season year
	keep `y' `x' inACO hrr_npi
	keep if `y' != .
	areg `y' inACO $x, absorb(hrr_npi) vce(cluster hrr_npi)
	mat a = _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	margins, at(inACO=(0 1))
	mat b = r(table)
	sca b = b[1,1]
	mat m[`n',1] = b, b + _b[inACO], a
	loc n = `n' + 1
	
	**********************************************
	**					IV model				**
	**********************************************
	use "$rif/quality_100pct", clear
	keep `y' $z diff50 inACO hrr_npi
	xtset hrr_npi
	xtivreg2 `y' $z (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin
	mat m[`n',1] = b, b + _b[inACO], _b[inACO], _b[inACO] - invnormal(0.975)*_se[inACO], _b[inACO] + invnormal(0.975)*_se[inACO], 2*normal(-abs( _b[inACO]/_se[inACO])), e(N), e(estatp)
	loc n = `n' + 2
	putexcel set "$mat/figure3_100pct", replace
	putexcel a1 = mat(m)
	mat li m
}

* mammography
foreach y in proc5 {

	**********************************************
	**		Adjusted longitudinal model			**
	**********************************************
	use "$rif/quality_100pct", clear
	loc x age race dual hcc dib esrd poverty hs college season year
	keep `y' `x' inACO hrr_npi
	keep if `y' != .
	areg `y' inACO $xm, absorb(hrr_npi) vce(cluster hrr_npi)
	mat a = _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	margins, at(inACO=(0 1))
	mat b = r(table)
	sca b = b[1,1]
	mat m[`n',1] = b, b + _b[inACO], a
	loc n = `n' + 1
	
	**********************************************
	**					IV model				**
	**********************************************
	use "$rif/quality_100pct", clear
	keep `y' $z diff50 inACO hrr_npi
	xtset hrr_npi
	xtivreg2 `y' $zm (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin
	mat m[`n',1] = b, b + _b[inACO], _b[inACO], _b[inACO] - invnormal(0.975)*_se[inACO], _b[inACO] + invnormal(0.975)*_se[inACO], 2*normal(-abs( _b[inACO]/_se[inACO])), e(N), e(estatp)
	loc n = `n' + 2
	putexcel set "$mat/figure3_100pct", replace
	putexcel a1 = mat(m)
	mat li m

}

**************************************************************
**				Figure 4: Clinician entry/exit				**
**************************************************************

use "$rif/npi_exit4_100pct", clear
tab exiting
tab entering

* Step 1): predicted (d) and adjusted (g) matrices
mat d = (1\5\10\25\50\75\90\95\99)	// predicted matrix
mat g50 = J(50,1,0)					// adjusted matrix
forvalues i=1/50 {
	sca j = (`i'*2)-1
	mat g50[`i',1] = j
}

loc z = 1

foreach y in exiting entering {

	use "$rif/npi_exit4_100pct", clear
	
	* average beneficiary characteristics of clinician patient panel across prior 3 years
	loc x3 "l3_hcc l3_age l3_female l3_black l3_hispanic l3_race_other l3_dual l3_dib l3_esrd l3_poverty l3_college l3_hs"
	loc n i.year
	keep `y' paylead3 hrr year `x3'
	keep if `y' != .

	* Step 1: estimate regression and predicted margins
	areg `y' c.(paylead3)##c.(paylead3)##c.(paylead3) `n' `x3', absorb(hrr) vce(cluster hrr)
	di "`e(cmdline)'"
	keep if e(sample) == 1
	est sto a

	_pctile paylead3, p(1 5 10 25 50 75 90 95 99)
	forvalues i=1/9 {
		sca p`i' = r(r`i')
		loc p`i' = p`i'
	}

	margins, at(paylead3=(`p1' `p2' `p3' `p4' `p5' `p6' `p7' `p8' `p9')) post
	mat a = r(table)'
	mat b`z' = a[1..9,1], a[1..9,5..6]

	* Step 2: estimate adjusted margins across spending bins
	xtile bin = paylead3, nq(50)
	qui areg `y' i.bin `n' `x3', absorb(hrr) vce(cluster hrr)
	qui margins, at(bin=(1(1)50)) post
	mat a = r(table)'
	mat e`z' = a[1..50,1], a[1..50,5..6]
	loc z = `z' + 1

}

* Step 3: combine predicted and adjusted margins
mat a1 = d,b1,b2
mat a2 = g50,e1,e2
mat a = a2\a1

putexcel set "$mat/fig4_100pct", replace
putexcel a1 = mat(a)


********************************************************************************
** 				Appendix Table 1. Non-linear estimation of 					  **
**					discount factor used in MSSP supply						  **
********************************************************************************

* estimate pr(belonging to same ACO) as a function of distance between two clinicians
use "$inst/npi/discount_npi1", clear
sum
keep aco_same dist2 ruca

forvalues i = 1/4 {
	nl (aco_same = dist2^(-1*{b1})) if ruca == `i'
}

********************************************************************************
** 			Appendix Table 3. Formal test of statistical differences	  	  **
**			between instrumental variable and adjusted longitudinal model	  **
********************************************************************************

* Difference-in-sargan:
mat m = J(11,1,0)
loc n = 1
loc y pay hip
foreach y in `y' {
	use "$rif/iv_100pct", clear
	keep `y' $z diff50 inACO hrr_npi
	keep if `y' != .
	xtset hrr_npi
	xtivreg2 `y' $z (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin endog(inACO)
	mat m[`n',1] = e(estatp)
	loc n = `n' + 1	
	putexcel set "$mat/sargan_100pct", replace
	putexcel a1 = mat(m)
	mat li m
}

loc y proc1 proc2 proc3 proc4
foreach y in `y' {
	use "$rif/quality_100pct", clear
	keep `y' $z diff50 inACO hrr_npi
	keep if `y' != .
	xtset hrr_npi
	xtivreg2 `y' $z (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin endog(inACO)
	mat m[`n',1] = e(estatp)
	loc n = `n' + 1	
	putexcel set "$mat/sargan_100pct", replace
	putexcel a1 = mat(m)
	mat li m
}

loc y proc5
use "$rif/quality_100pct", clear
keep `y' $zm diff50 inACO hrr_npi	
xtset hrr_npi
xtivreg2 `y' $zm (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin endog(inACO)
mat m[`n',1] = e(estatp)
loc n = `n' + 1	
putexcel set "$mat/sargan_100pct", replace
putexcel a1 = mat(m)
mat li m

loc y admission acsc readmission ed
foreach y in `y' {
	use "$rif/iv_100pct", clear
	keep `y' $z diff50 inACO hrr_npi
	keep if `y' != .
	xtset hrr_npi
	xtivreg2 `y' $z (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin endog(inACO)
	mat m[`n',1] = e(estatp)
	loc n = `n' + 1	
	putexcel set "$mat/sargan_100pct", replace
	putexcel a1 = mat(m)
	mat li m
}

********************************************************************************
**			Appendix Table 4. Sensitivity analyses of the association 		  **
**			between the MSSP and change in hip fracture and spending		  **
********************************************************************************

mat m = J(18,5,0)
loc n = 1
	
foreach y in hip pay {

	******************************************************
	**		Fixed effects specification tests			**
	******************************************************

	* Adjusted longitudinal + MA penetration
	use "$rif/robust_100pct", clear
	loc w age dual hcc poverty hs college season year female race dib esrd 
	keep pay hip `w' inACO hrr_npi hrr_id npi id hrr_year ma itt_id itt_npi balanced_id
	areg `y' inACO ma $x, absorb(hrr_npi) vce(cluster hrr_npi)
	mat m[`n',1] = _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	loc n = `n' + 1

	* Beneficiary fixed effects + intention-to-treat
	loc i id
	loc x itt_`i'
	areg `y' `x' $b, absorb(`i') vce(cluster hrr_`i')
	mat m[`n',1] = _b[`x'], _b[`x'] - invttail(e(df_r),0.025)*_se[`x'], _b[`x'] + invttail(e(df_r),0.025)*_se[`x'], 2*ttail(e(df_r),abs(_b[`x']/_se[`x'])), e(N)
	loc n = `n' + 1

	* Beneficiary fixed effects + intention-to-treat + balanced panel
	areg `y' `x' $b if balanced_id == 1, absorb(`i') vce(cluster hrr_`i')
	mat m[`n',1] = _b[`x'], _b[`x'] - invttail(e(df_r),0.025)*_se[`x'], _b[`x'] + invttail(e(df_r),0.025)*_se[`x'], 2*ttail(e(df_r),abs(_b[`x']/_se[`x'])), e(N)
	loc n = `n' + 1
	
	* Clinician fixed effects + intention-to-treat
	loc i npi
	loc x itt_`i'
	areg `y' `x' $x, absorb(`i') vce(cluster hrr_`i')
	mat m[`n',1] = _b[`x'], _b[`x'] - invttail(e(df_r),0.025)*_se[`x'], _b[`x'] + invttail(e(df_r),0.025)*_se[`x'], 2*ttail(e(df_r),abs(_b[`x']/_se[`x'])), e(N)
	loc n = `n' + 1

	* Clinician fixed effects + intention-to-treat + balanced panel
	keep if balanced_id == 1
	areg `y' `x' $x, absorb(`i') vce(cluster hrr_`i')
	mat m[`n',1] = _b[`x'], _b[`x'] - invttail(e(df_r),0.025)*_se[`x'], _b[`x'] + invttail(e(df_r),0.025)*_se[`x'], 2*ttail(e(df_r),abs(_b[`x']/_se[`x'])), e(N)
	loc n = `n' + 1
	
	******************************************************
	**				IV specification tests				**
	******************************************************

	* Instrumental variable + MA penetration
	use "$rif/iv_robust_100pct", clear
	keep `y' $z inACO diff50 hrr_npi hrr_year ma
	xtset hrr_npi
	xtivreg2 `y' $z ma (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin
	mat m[`n',1] = _b[inACO], _b[inACO] - invnormal(0.975)*_se[inACO], _b[inACO] + invnormal(0.975)*_se[inACO], 2*normal(-abs( _b[inACO]/_se[inACO])), e(N)	
	loc n = `n' + 1
	drop ma

	* Instrumental variable + market-year fixed effects
	xtset hrr_year
	loc x age female black hispanic race_other dual hcc dib esrd poverty hs college season2 season3 season4
	xtivreg2 `y' `x' (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin
	mat m[`n',1] = _b[inACO], _b[inACO] - invnormal(0.975)*_se[inACO], _b[inACO] + invnormal(0.975)*_se[inACO], 2*normal(-abs( _b[inACO]/_se[inACO])), e(N)	
	loc n = `n' + 1
	
	* Instrumental variable + market-year + ACO fixed effects
	use "$rif/aco_demean_100pct", clear
	loc c_x c_age c_female c_black c_hispanic c_race_other c_hcc c_dual c_dib c_esrd c_poverty c_hs c_college c_year2009 c_year2010 c_year2011 c_year2012 c_year2013 c_year2014 
	keep c_hip c_pay c_inACO c_diff50 `c_x' hrr_npi hrr_year
	xtset hrr_year
	xtivreg2 c_`y' `c_x' (c_inACO = c_diff50), cluster(hrr_npi) fe robust noid nocollin
	mat m[`n',1] = _b[c_inACO], _b[c_inACO] - invnormal(0.975)*_se[c_inACO], _b[c_inACO] + invnormal(0.975)*_se[c_inACO], 2*normal(-abs( _b[c_inACO]/_se[c_inACO])), e(N)	
	loc n = `n' + 1
	
	putexcel set "$mat/robust_100pct", replace
	putexcel a1 = mat(m)
	mat li m
	loc n = `n' + 1
}

********************************************************************************
**	Appendix Table 5. Estimates of the association between clinician prior	  **
**		 spending and probability of joining or exiting an MSSP ACO			  **
********************************************************************************

loc x3 "l3_hcc l3_age l3_female l3_black l3_hispanic l3_race_other l3_dual l3_dib l3_esrd l3_poverty l3_college l3_hs"

loc i = 1

foreach y in entering exiting {

	use "$rif/npi_exit4_100pct", clear
	keep `y' paylead3 hrr year `x3'
	keep if `y' != .
	areg `y' c.(paylead3)##c.(paylead3)##c.(paylead3) i.year `x3', absorb(hrr) vce(cluster hrr)
	di "`e(cmdline)'"
	est sto a
	keep if e(sample) == 1
	margins, at((p1) paylead3) at((p5) paylead3) at((p10) paylead3) at((p25) paylead3) at((p50) paylead3) at((p75) paylead3) at((p90) paylead3) at((p95) paylead3) at((p99) paylead3) post
	mat r = r(table)
	mat r = r'
	mat a`i' = r[1..9,1], r[1..9,5..6], r[1..9,4]
	mat li a`i'

	* absolute differences
	mat b`i' = J(9,2,0)
	loc n = 1
	forvalues j=1/9 {
		est resto a
		margins, at((p1) paylead3) at((p5) paylead3) at((p10) paylead3) at((p25) paylead3) at((p50) paylead3) at((p75) paylead3) at((p90) paylead3) at((p95) paylead3) at((p99) paylead3) post
		lincom _b[`j'._at]-_b[5._at]
		mat b`i'[`n',1] = r(estimate), 2*normal(-abs(r(estimate)/r(se)))
		loc n = `n' + 1
	}

	* relative differences (risk ratio)
	mat d`i' = J(9,2,0)
	loc n = 1
	forvalues j=1/9 {
		est resto a
		margins, at((p1) paylead3) at((p5) paylead3) at((p10) paylead3) at((p25) paylead3) at((p50) paylead3) at((p75) paylead3) at((p90) paylead3) at((p95) paylead3) at((p99) paylead3) post
		nlcom _b[`j'._at]/_b[5._at], post
		mat d`i'[`n',1] = _b[_nl_1], 2*normal(-abs( _b[_nl_1]/_se[_nl_1]))
		loc n = `n' + 1
	}
	loc i = `i' + 1
}

mat li a1
mat li b1
mat li d1
mat a1 = a1, b1, d1
mat a2 = a2, b2, d2

mat j = J(1,8,0)

mat a = j \ a1 \ j \ a2
mat li a

putexcel set "$mat/npi_bene_exit_table_100pct", replace
putexcel a1 = mat(a)

********************************************************************************
**		Appendix Table 6. Probability of beneficiary exit from the MSSP		  **
**		depending on whether attributed clinicians remains in MSSP or exits	  **
********************************************************************************

use "$rif/npi_bene_exit_pay_100pct", clear
keep if exiting != .

* beneficiary exit as function of prior year spending, RHS, and whether NPI exited MSSP
loc x1 l1_hcc l1_age l1_female l1_black l1_hispanic l1_race_other l1_dual l1_dib l1_esrd l1_poverty l1_college l1_hs l1_ma
areg exiting i.aco_npi_2013##c.paylead1 `x1', absorb(hrr) vce(cluster hrr)
est sto a
keep if e(sample) == 1

* average probability of exit if clinician remained vs. exited
margins, at(aco_npi_2013=(1 0)) post coefl	// aco_npi_2013 = 1 if clinician remains in MSSP in 2014
mat r = r(table)'
mat c1 = r[1,1], r[1,5..6], r[1,4]
mat c2 = r[2,1], r[2,5..6], r[2,4]
lincom _b[2._at]-_b[1._at]
mat c3 = r(estimate), r(estimate) - invttail(e(df_r),0.025)*r(se), r(estimate) + invttail(e(df_r),0.025)*r(se), 2*ttail(e(df_r),abs(r(estimate)/r(se)))

est resto a
_pctile paylead1 if e(sample) == 1, p(1 5 10 25 50 75 90 95 99)
forvalues i=1/9 {
	sca p`i' = r(r`i')
	loc p`i' = p`i'
}
margins, at(aco_npi_2013=(1 0) paylead1=(`p1' `p2' `p3' `p4' `p5' `p6' `p7' `p8' `p9')) post
mat r = r(table)'
mat a1 = r[1..9,1], r[1..9,5..6], r[1..9,4]
mat a2 = r[10..18,1], r[10..18,5..6], r[10..18,4]

* (A) exit as function of spending among beneficiaries whose clinicians remain in ACO
est resto a
loc i = 1

* absolute differences
margins, at(aco_npi_2013=(1) paylead1=(`p1' `p2' `p3' `p4' `p5' `p6' `p7' `p8' `p9')) post
mat b`i' = J(9,2,0)
loc n = 1
forvalues j=1/9 {
	lincom _b[`j'._at]-_b[5._at]
	mat b`i'[`n',1] = r(estimate),  2*ttail(r(df),abs(r(estimate)/r(se)))
	loc n = `n' + 1
}

* relative differences (risk ratio)
mat d`i' = J(9,2,0)
loc n = 1
forvalues j=1/9 {
	est resto a
	margins, at(aco_npi_2013=(1) paylead1=(`p1' `p2' `p3' `p4' `p5' `p6' `p7' `p8' `p9')) post
	nlcom _b[`j'._at]/_b[5._at], post
	mat d`i'[`n',1] = _b[_nl_1], 2*normal(-abs( _b[_nl_1]/_se[_nl_1]))
	loc n = `n' + 1
}

est resto a
loc i = `i' + 1

* (B) exit as function of spending among beneficiaries whose clinicians exit ACO
* absolute differences
margins, at(aco_npi_2013=(0) paylead1=(`p1' `p2' `p3' `p4' `p5' `p6' `p7' `p8' `p9')) post
mat b`i' = J(9,2,0)
loc n = 1
forvalues j=1/9 {
	lincom _b[`j'._at]-_b[5._at]
	mat b`i'[`n',1] = r(estimate),  2*ttail(r(df),abs(r(estimate)/r(se)))
	loc n = `n' + 1
}

* relative differences (risk ratio)
mat d`i' = J(9,2,0)
loc n = 1
forvalues j=1/9 {
	est resto a
	margins, at(aco_npi_2013=(0) paylead1=(`p1' `p2' `p3' `p4' `p5' `p6' `p7' `p8' `p9')) post
	nlcom _b[`j'._at]/_b[5._at], post
	mat d`i'[`n',1] = _b[_nl_1], 2*normal(-abs( _b[_nl_1]/_se[_nl_1]))
	loc n = `n' + 1
}

mat z1 = a1, b1, d1
mat z3 = a2, b2, d2

mat j = J(1,8,0)
mat c = J(1,4,0)

mat a = j \ z1 \ (c1,c) \ j \ z3 \ (c2,c) \ (c3,c)
mat li a

putexcel set "$mat/npi_bene_exit_table_100pct", replace
putexcel a1 = mat(a)

********************************************************************************
** 		Appendix Table 7: Simulated changes in average MSSP beneficiary   	  **
**		spending with attrition of high-cost clinicians from MSSP ACOs		  **
********************************************************************************

clear all
use "$rif/npi_exit4_100pct", clear
sort npi year

* restrict to ACO clinicians in 2013 or 2014 who are in ACOs that formed in 2012 or 2013
keep if year >= 2013 & aco_start != 2014 & aco == 1

* calculate average spending if clinicians in >= 95th or >= 99th spending percentile omitted
xtile pay_20 = pay, nq(20)
xtile pay_100 = pay, nq(100)

mat a = J(3,2,0)
sum pay [aw=n]
mat a = r(mean), 0, r(N), r(sum_w)

sum pay if pay_20 < 20 [aw=n]
mat b = r(mean), 0, r(N), r(sum_w)

sum pay if pay_100 < 100 [aw=n]
mat c = r(mean), 0, r(N), r(sum_w)

sca ba = 100*((b[1,1]-a[1,1])/a[1,1])
sca li ba

sca ca = 100*((c[1,1]-a[1,1])/a[1,1])
sca li ca

mat a = (.,.,.,.) \ a \ c[1,1], ca, c[1,3..4] \  b[1,1], ba, b[1,3..4] 
mat li a

svmat a

gen spend = strofreal(round(a1, 1), "%4.0fc")
replace spend = "Average MSSP ACO spending per beneficiary-year ($)" in 1

gen diff = strofreal(round(a2, .1), "%4.1fc")
replace diff = "Difference in calculated average MSSP ACO spending (%)" in 1
replace diff = "--" in 2

gen ng = strofreal(a3)
replace ng = "Sample size (clinician-years)" in 1

gen nb = strofreal(a4)
replace nb = "Sample size (beneficiary-years)" in 1

gen sample = ""
replace sample = "Clinicians in MSSP ACOs (2013-2014)" in 1
replace sample = "All clinicians" in 2
replace sample = "Excluding clinicians with ≥ 99th of average spending" in 3
replace sample = "Excluding clinicians with ≥ 95th of average spending" in 4

keep in 1/4
loc x sample spend diff ng nb
order `x'
keep `x'
export excel using "$fig/etable7_npi_exiting_100pct", replace

********************************************************************************
**	Appendix Table 8. Simulated changes in average MSSP beneficiary spending  **
**	when including clinicians exiting MSSP and excluding clinicians entering  **
**							the MSSP (2013-2014)							  **	
********************************************************************************

clear all
use "$rif/npi_exit4_100pct", clear
sort npi year

bysort npi: egen aco_ever = max(aco)

* restrict to ACO clinicians in 2013 or 2014 who are in ACOs that formed in 2012 or 2013
drop if year < 2013 | aco_start == 2014 | aco_ever == 0

sum pay if aco == 1 [aw=n]
mat a = r(mean), 0, r(N), r(sum_w)

sum pay if aco_ever == 1 [aw=n]
mat b = r(mean), 0, r(N), r(sum_w)

sca ba = 100*((b[1,1]-a[1,1])/a[1,1])
sca li ba

mat j = J(1,4,.)
mat a = j \ a \ b[1,1], ba, b[1,3..4]
mat li a
clear
svmat a

gen spend = strofreal(round(a1, 1), "%4.0fc")
replace spend = "Average MSSP ACO spending per beneficiary-year ($)" in 1

gen diff = strofreal(round(a2, .1), "%4.1fc")
replace diff = "Difference in calculated average MSSP ACO spending (%)" in 1
replace diff = "--" in 2

gen ng = strofreal(a3)
replace ng = "Sample size (clinician-years)" in 1

gen nb = strofreal(a4)
replace nb = "Sample size (beneficiary-years)" in 1

gen sample = ""
replace sample = "Clinicians in MSSP ACOs" in 1
replace sample = "Observed sample of MSSP participants" in 2
replace sample = "Simulated sample including clinicians exiting MSSP and clinicians entering MSSP" in 3

loc x sample spend diff ng nb
order `x'
keep `x'

export excel using "$fig/etable8_npi_sim_all_100pct", replace

********************************************************************************
**								Figures									   	  **
********************************************************************************

set scheme s2color
grstyle init
grstyle graphsize x 5.5
grstyle graphsize y 4.5
grstyle color background white
grstyle anglestyle vertical_tick horizontal
grstyle yesno draw_major_hgrid yes
grstyle color major_grid gs8
grstyle linewidth major_grid thin
grstyle linepattern major_grid dot
grstyle numstyle legend_cols
grstyle linestyle legend none
grstyle set legend 10, inside
grstyle linewidth plineplot medthick
grstyle yesno grid_draw_min yes
grstyle yesno grid_draw_max yes

********************************************************************************
* Appendix Figure 2: pre-period spending trends across instrumental variable   *						   *
*									and										   *
* 		Appendix Table 2. Pre-period trends in spending across observed 
*				MSSP status and MSSP supply instrumental variable			   *
********************************************************************************

use "$rif/robust_100pct", clear

* keep beneficiaries who appear in post-period (i.e., with calculated MSSP supply)
keep if t50 != .

gen in_qtr = qtr if inACO == 1
bysort id: egen in_qtr_min = min(in_qtr)
gen byte cohort = 0
replace cohort = 1 if in_qtr_min == 18
replace cohort = 2 if in_qtr_min == 19
replace cohort = 3 if in_qtr_min == 21
replace cohort = 4 if in_qtr_min == 25
drop in_qtr in_qtr_min

keep pay t50 tr_id cohort qtr age dual hcc poverty hs college female race dib esrd hrr_npi

mkspline prespline 18 postspline = qtr

loc w age dual hcc poverty hs college female i.race dib esrd
	
foreach x in t50 tr_id {
	
	*********************************
	*		Spline regression		*
	*********************************
	
	areg pay `w' i.`x'##c.(prespline postspline), absorb(hrr_npi) vce(cluster hrr_npi)
	margins, at(`x'=(0 1) prespline = (1(1)18) postspline=0) at(`x'=(0 1) prespline = 18 postspline=(1(1)10)) post

	// control, pre-period: 1/18, post-period: 37/46
	mat c_sp = J(28,3,0)
	loc r = 1
	forvalues j = 1/18 {
		mat c_sp[`r',1] = _b[`j'._at], _b[`j'._at] - invnormal(0.975)*_se[`j'._at], _b[`j'._at] + invnormal(0.975)*_se[`j'._at]	
		loc r = `r' + 1
	}
	loc r = 19
	forvalues j = 37/46 {
		mat c_sp[`r',1] = _b[`j'._at], _b[`j'._at] - invnormal(0.975)*_se[`j'._at], _b[`j'._at] + invnormal(0.975)*_se[`j'._at]	
		loc r = `r' + 1
	}

	// treatment, pre-period: 15/28, post-period: 47/56
	mat t_sp = J(28,3,0)
	loc r = 1
	forvalues j = 19/36 {
		mat t_sp[`r',1] = _b[`j'._at], _b[`j'._at] - invnormal(0.975)*_se[`j'._at], _b[`j'._at] + invnormal(0.975)*_se[`j'._at]	
		loc r = `r' + 1
	}
	loc r = 19
	forvalues j = 47/56 {
		mat t_sp[`r',1] = _b[`j'._at], _b[`j'._at] - invnormal(0.975)*_se[`j'._at], _b[`j'._at] + invnormal(0.975)*_se[`j'._at]	
		loc r = `r' + 1
	}

	mat li c_sp
	mat li t_sp
	
	*****************************
	*		Adjusted scatter	*
	*****************************

	areg pay `w' i.`x'##i.qtr, absorb(hrr_npi) vce(cluster hrr_npi)
	margins, at(`x'=(0 1) qtr=(1(1)28)) post

	// control: 1..24
	mat c = J(28,3,0)
	loc i = 1
	forvalues j = 1/28 {
		mat c[`i',1] = _b[`j'._at], _b[`j'._at] - invnormal(0.975)*_se[`j'._at], _b[`j'._at] + invnormal(0.975)*_se[`j'._at]	
		loc i = `i' + 1
	}

	// treatment: 29..56
	mat t = J(28,3,0)
	loc i = 1
	forvalues j = 29/56 {
		mat t[`i',1] = _b[`j'._at], _b[`j'._at] - invnormal(0.975)*_se[`j'._at], _b[`j'._at] + invnormal(0.975)*_se[`j'._at]
		loc i = `i' + 1
	}

	mat li t
	mat li c

	mat x = J(28,1,0)
	forvalues i = 1/28 {
		mat x[`i',1] = `i'
	}

	mat a = x, c, t, c_sp, t_sp
	mat li a

	putexcel set "$mat/spline_`x'_100pct", replace
	putexcel a1 = mat(a)
}

********************************************************************************
**			Appendix Figure 3: pre-period trends in spending 				  **
**					across 4 ACO cohorts and controls						  **
********************************************************************************

areg pay `w' i.cohort##i.qtr, absorb(hrr_npi) vce(cluster hrr_npi)
margins, at(cohort=(0(1)4) qtr=(1(1)28)) post

loc j = 1
loc k = 28
forvalues i = 1/5 {
	mat s_`i' = J(28,3,0)
	loc r = 1
	forvalues l = `j'/`k' {
		mat s_`i'[`r',1] = _b[`l'._at], _b[`l'._at] - invnormal(0.975)*_se[`l'._at], _b[`l'._at] + invnormal(0.975)*_se[`l'._at]	
		loc r = `r' + 1
	}
	loc j = `j' + 28
	loc k = `k' + 28
}

mat x = J(28,1,0)
forvalues i = 1/28 {
	mat x[`i',1] = `i'
}

mat a = x, s_1, s_2, s_3, s_4, s_5
mat li a
putexcel set "$mat/spline_cohort_100pct", replace
putexcel a1 = mat(a)

********************************************************************************
** 				Appendix Figure 4. Changes in total Medicare spending 		  **
**					according to ACOs' year of entry into the MSSP			  **
********************************************************************************

use "$rif/subgroup_100pct", clear
xtset hrr_npi

mat m = J(6,5,0)
loc n = 1

forvalues i = 1/3 {
	
	**********************************
	**		Adjusted Longitudinal	**
	**********************************
	
	areg pay inACO_`i' $x, absorb(hrr_npi) vce(cluster hrr_npi)
	mat m[`n',1] = _b[inACO_`i'], _b[inACO_`i'] - invttail(e(df_r),0.025)*_se[inACO_`i'], _b[inACO_`i'] + invttail(e(df_r),0.025)*_se[inACO_`i'], 2*ttail(e(df_r),abs(_b[inACO_`i']/_se[inACO_`i'])), e(N)
	putexcel set "$mat/subgroup_100pct", replace
	putexcel a1 = mat(m)
	loc n = `n' + 1

	**********************************
	**				IV				**
	**********************************
	
	xtivreg2 pay $z (inACO_`i' = diff50), cluster(hrr_npi) fe robust noid nocollin
	mat m[`n',1] = _b[inACO_`i'], _b[inACO_`i'] - invnormal(0.975)*_se[inACO_`i'], _b[inACO_`i'] + invnormal(0.975)*_se[inACO_`i'], 2*normal(-abs( _b[inACO_`i']/_se[inACO_`i'])), e(N)	
	putexcel set "$mat/subgroup_100pct", replace
	putexcel a1 = mat(m)
	loc n = `n' + 1
		
	mat li m
}

********************************************************************************
** 			Appendix Figure 5. Robustness tests of MSSP supply 				  **
**						as an instrumental variable							  **
********************************************************************************

use "$rif/iv_robust_100pct", clear
keep pay hip $z inACO diff* hrr_npi
xtset hrr_npi

* IV models using MSSP supply restricted radius of 10, 25, 50, 75, or 150 miles, with or without discounting
** diff10 = discounted supply within 10mi, diff11 = non-discounted supply within 10mi, diff50 = discounted supply within 50mi, ...,
foreach y in hip pay {
	mat m = J(12,5,0)
	loc n = 1
	foreach d in 10 11 25 26 50 51 75 76 150 151 {
		xtivreg2 `y' $z (inACO = diff`d'), cluster(hrr_npi) fe robust noid nocollin
		mat m[`n',1] = _b[inACO], _b[inACO] - invnormal(0.975)*_se[inACO], _b[inACO] + invnormal(0.975)*_se[inACO], 2*normal(-abs( _b[inACO]/_se[inACO])), e(N)
		loc n = `n' + 1
		putexcel set "$mat/iv_robust_`y'_100pct", replace
		putexcel a1 = mat(m)
	}
}

********************************************************************************
**			Appendix Figure 6: County-level MA and ACO penetration 			  **
********************************************************************************

**********************************************
**		First-differences scatter plot		**
**********************************************

use "$iv/rif/aco_ma_100pct", clear

keep if year == 2011 | year == 2014
replace year = 1 if year == 2011
replace year = 2 if year == 2014
xtset county year

loc x ma aco
foreach x in `x' {
	gen d_`x' = d.`x'
}
global dx
loc x age female black hispanic race_other dual hcc dib esrd poverty hs college
foreach x in `x' {
	gen d_`x' = d.`x'
	global dx $dx d_`x'
}

keep if year == 2

reg d_aco d_ma $dx
predict hat
graph twoway (scatter d_aco d_ma, m(triangle) msize(tiny) mcolor(orange)) ///
(lfit hat d_ma, lcolor(gs8) lwidth(medium)), ///
graphregion(style(none) color(gs16)) legend(off) ///
ytitle("Change in MSSP penetration (%)", size(medium)) ///
xtitle("Change in MA penetration (%)", size(medium)) ///
yscale(titlegap(*10)) xscale(titlegap(*10)) ///
title("", span size(small) color(black))
graph export "$fig/efig6_scatter_penetration.png", replace as(png) width(1000)

********************************************************************************
* 	Appendix Figure 7: Change in hospital use in adjusted and IV model 		   *
********************************************************************************

mat m = J(14,7,0)
loc n = 3
foreach y in admission acsc readmission ed {

	**********************************************
	**		Adjusted longitudinal model			**
	**********************************************
	use "$rif/fe_100pct", clear
	loc x age female race dual hcc dib esrd poverty hs college season year
	keep `y' `x' inACO hrr_npi
	
	areg `y' inACO $x, absorb(hrr_npi) vce(cluster hrr_npi)
	mat a = _b[inACO], _b[inACO] - invttail(e(df_r),0.025)*_se[inACO], _b[inACO] + invttail(e(df_r),0.025)*_se[inACO], 2*ttail(e(df_r),abs(_b[inACO]/_se[inACO])), e(N)
	margins, at(inACO=(0 1))
	mat b = r(table)
	sca b = b[1,1]
	mat m[`n',1] = b, b + _b[inACO], a
	loc n = `n' + 1
	
	**********************************************
	**					IV model				**
	**********************************************
	use "$rif/iv_100pct", clear
	keep `y' $z diff50 inACO hrr_npi
	xtset hrr_npi
	xtivreg2 `y' $z (inACO = diff50), cluster(hrr_npi) fe robust noid nocollin
	mat m[`n',1] = b, b + _b[inACO], _b[inACO], _b[inACO] - invnormal(0.975)*_se[inACO], _b[inACO] + invnormal(0.975)*_se[inACO], 2*normal(-abs( _b[inACO]/_se[inACO])), e(N)	
	loc n = `n' + 2
	
	putexcel set "$mat/efigure7_100pct", replace
	putexcel a1 = mat(m)
	mat li m
}

********************************************************************************
** 			Appendix Figure 8. Association between provider group spending 	  **
**			performance and probability of exiting or entering an MSSP ACO    **
********************************************************************************

use "$rif/tin_exit2_100pct", clear
tab exiting
tab entering

* Step 1): predicted (d) and adjusted (g) matrices
mat d = (1\5\10\25\50\75\90\95\99)	// predicted matrix
mat g50 = J(50,1,0)					// adjusted matrix
forvalues i=1/50 {
	sca j = (`i'*2)-1
	mat g50[`i',1] = j
}

loc z = 1
foreach y in exiting entering {

	use "$rif/tin_exit2_100pct", clear
	
	* average beneficiary characteristics of clinician patient panel across prior 3 years
	loc x3 "l3_hcc l3_age l3_female l3_black l3_hispanic l3_race_other l3_dual l3_dib l3_esrd l3_poverty l3_college l3_hs"
	loc n i.year
	keep `y' paylead3 hrr year `x3'
	keep if `y' != .

	* Step 2: estimate regression and save predicted margins
	areg `y' c.(paylead3)##c.(paylead3)##c.(paylead3) `n' `x3', absorb(hrr) vce(cluster hrr)
	di "`e(cmdline)'"
	keep if e(sample) == 1
	est sto a

	_pctile paylead3, p(1 5 10 25 50 75 90 95 99)
	forvalues i=1/9 {
		sca p`i' = r(r`i')
		loc p`i' = p`i'
	}

	margins, at(paylead3=(`p1' `p2' `p3' `p4' `p5' `p6' `p7' `p8' `p9')) post
	mat a = r(table)'
	mat b`z' = a[1..9,1], a[1..9,5..6]

	* Step 3: save adjusted margins across spending bins
	xtile bin = paylead3, nq(50)
	qui areg `y' i.bin `n' `x3', absorb(hrr) vce(cluster hrr)
	qui margins, at(bin=(1(1)50)) post
	mat a = r(table)'
	mat e`z' = a[1..50,1], a[1..50,5..6]
	
	loc z = `z' + 1

}

* Step 3: combine predicted and adjusted margins
mat a1 = d,b1,b2
mat a2 = g50,e1,e2
mat a = a2\a1

putexcel set "$mat/efig8_tin_100pct", replace
putexcel a1 = mat(a)

********************************************************************************
**		Appendix Figure 9. Association between beneficiary spending and		  **
**		probability of beneficiary exiting the MSSP among beneficiaries		  **
**				whose clinicians exited or remained in the MSSP				  **
********************************************************************************

use "$rif/npi_bene_exit_pay_100pct", clear
keep if exiting != .

** Step 1: predicted beneficiary exit as function of prior year spending, RHS, and whether attributed clinician remained in or exited MSSP
loc x1 l1_hcc l1_age l1_female l1_black l1_hispanic l1_race_other l1_dual l1_dib l1_esrd l1_poverty l1_college l1_hs l1_ma
areg exiting i.aco_npi_2013##c.paylead1 `x1', absorb(hrr) vce(cluster hrr)
est sto a
keep if e(sample) == 1

* mean difference in predicted beneficiary exit if attributed clinician remains in MSSP vs. exits MSSP
margins, at(aco_npi_2013=(1 0))	// aco_npi_2013 = 1 if attributed clinician (in 2013) remains in MSSP in 2014

* predicted margins across percentiles of beneficiary spending
_pctile paylead1 if e(sample) == 1, p(1 5 10 25 50 75 90 95 99)
forvalues i=1/9 {
	sca p`i' = r(r`i')
	loc p`i' = p`i'
}

* predicted margins across percentiles of beneficiary spending and whether attributed clinician remained in or exited MSSP 
foreach c in 0 1 {
	margins, at(paylead1=(`p1' `p2' `p3' `p4' `p5' `p6' `p7' `p8' `p9') aco_npi=(`c')) post
	mat a = r(table)'
	mat a = a[1..9,1], a[1..9,5..6]
	loc i = `c' + 1
	mat d`i' = (1\5\10\25\50\75\90\95\99)
	mat d`i' = d`i', a
	est resto a
}

* step 2: adjusted beneficiary exit across spending bins

* bin beneficiary spending into 50 spending quantiles
xtile paylead1_bin = paylead1 if e(sample) == 1, nq(50)

* adjusted spending across 50 spending quantiles and whether attributed clinician remained in or exited MSSP 
areg exiting `x1' i.(paylead1_bin)##i.(aco_npi_2013), absorb(hrr) vce(cluster hrr)
est sto b
mat c = J(50,1,0)
forvalues i=1/50 {
	sca j = (`i'*2)-1
	mat c[`i',1] = j
}
foreach c in 0 1 {	
	margins, at(aco_npi_2013=(`c') paylead1_bin=(1(1)50)) post
	mat b = r(table)'
	loc i = `c' + 1
	mat b`i' = b[1..50,1], b[1..50,5..6]
	est resto b
}
mat c = (c\c),(b1\b2)

* step 3: combine predicted and adjusted margins
mat a = c\d1\d2
mat li a
putexcel set "$mat/efig9_npi_bene_exit_figure_100pct", replace
putexcel a1 = mat(a)
