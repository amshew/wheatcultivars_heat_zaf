*Replication Code for main results in: "Yield reduction under climate warming varies among wheat cultivars in South Africa";
*Authors: Shew, Tack, Nalley, and Chaminuka
*Contact: ashew@astate.edu
*2020-08-10
*Data for replication: https://doi.org/10.7910/DVN/8Y6Q7F

#delimit;
clear all;
set more off;
set maxvar 32000;
set matsize 8000;


*Directories: create and define your path to data and output folders;
global dd = "xx/data";  
global do = "xx/output";  


********************************************************************************;
*Figure 2 Marginal Effects of Temperature Bins on Wheat Yields;
********************************************************************************;
use "$dd\RegressionDataFinal.dta", clear;
reg ly i.si i.cu i.year prec c.prec#c.prec fr Bin*, cluster(ProYear) ;

gen temp = _n - 6;
gen ILT0 = 0;
replace ILT0 = 1 if temp < 0;
forvalues i = 0(5)45 {;
  local ip1 = `i' + 1;
  local ip4 = `i' + 4;
  display "`i',`ip4'";
  gen I`i'_`ip4' = 0;
  replace I`i'_`ip4' = 1 if temp >= `i' & temp <= `ip4';
};
predictnl mePref = ILT0*_b[fr] + 
                   I0_4*_b[BinExp0_4] +
                   I5_9*_b[BinExp5_9] +
                   I10_14*_b[BinExp10_14] +
                   I15_19*_b[BinExp15_19] +
                   I20_24*_b[BinExp20_24] +
                   I25_29*_b[BinExp25_29] +
                   I30_34*_b[BinExp30_inf]                    
				   , se(SEmePref);
drop if temp >= 35;
egen avg = mean(mePref);
replace mePref = mePref - avg;
gen LmePref = mePref - 1.96*SEmePref;
gen HmePref = mePref + 1.96*SEmePref;
twoway (line mePref temp, lc(black)) (line LmePref temp, lp(dash) lc(black)) 
       (line HmePref temp, lp(dash) lc(black)) , 
  legend(order(1 2) label(1 "Marginal Effect") label(2 "95% Confidence Band") pos(7) ring(0) row(3)) graphregion(color(white))
  ytitle("log-yield (bu. per acre)", size(large)) ysca(titlegap(0)) xtitle("temperature ({c 176}C)", size(large)) 
  ylabel(-.15(.05)0.05, nogrid labsize(large)) xlabel(-5(5)35, labsize(large)) name(MargEffect, replace); 
graph export "$do\Figure2.emf", replace;


********************************************************************************;
*Figure 3 Wheat Yield Impacts by +1 to +3C Warming Scenarios;
********************************************************************************;
use "$dd\RegressionDataFinal.dta", clear;
reg ly i.si i.cu i.year prec c.prec#c.prec fr Bin*, cluster(ProYear) ;
est save "$do\Regression", replace;

use "$dd\WarmingScenarios.dta", clear;
forvalues s = 1/3 {;
  est use "$do\Regression";
  nlcom exp(
   (frS`s'-frS0)*_b[fr] +
   (BinExp0_4S`s'-BinExp0_4S0)*_b[BinExp0_4] +
   (BinExp5_9S`s'-BinExp5_9S0)*_b[BinExp5_9] +
   (BinExp10_14S`s'-BinExp10_14S0)*_b[BinExp10_14] +
   (BinExp15_19S`s'-BinExp15_19S0)*_b[BinExp15_19] +
   (BinExp20_24S`s'-BinExp20_24S0)*_b[BinExp20_24] +
   (BinExp25_29S`s'-BinExp25_29S0)*_b[BinExp25_29] +
   (BinExp30_infS`s'-BinExp30_infS0)*_b[BinExp30_inf]) - 1, post;
  parmest, format(estimate min95 max95 %10.0g p %8.le) saving("$do\ImpactS`s'.dta",replace);
};
use "$do\ImpactS1.dta", clear;
gen shift = 1;
append using "$do\ImpactS2.dta";
replace shift = 2 if shift == .;
append using "$do\ImpactS3.dta";
replace shift = 3 if shift == .;
gen b = 100*estimate;
gen hib=100*max95;
gen lwb=100*min95;
gen order = _n;
graph twoway 
  (bar b order if shift==1, color(red)) 
  (bar b order if shift==2, color(red)) 
  (bar b order if shift==3, color(red)) 
  (rcap hib lwb order, color(blue)),
  xlabel(1 "+1{char 0176}C"  2 "+2{char 0176}C" 3 "+3{char 0176}C", noticks labsize(large))
  ytitle("yield impact (%)", size(large)) yscale() xtitle("warming scenario", size(large)) 
  graphregion(color(white)) ylabel(-40(10)0, nogrid labsize(large)) legend(off) ;
graph export "$do\Figure3.emf", replace;


********************************************************************************;
*Figure 4 Warming Estimates from Alternative Model that Includes Interaction of Heat and Precipitation Variables;
********************************************************************************;
use "$dd\RegressionDataFinal.dta", clear;
summ prec, detail;
gen lowP10 = (prec < 56);
reg ly i.si i.cu i.year prec c.prec#c.prec fr Bin* 1.lowP10#c.BinExp30_inf, cluster(ProYear) ;
est save "$do\Regression_LP10", replace;

use "$dd\WarmingScenarios.dta", clear;
forvalues s = 1/3 {;
  est use "$do\Regression_LP10";
  nlcom exp(
   (frS`s'-frS0)*_b[fr] +
   (BinExp0_4S`s'-BinExp0_4S0)*_b[BinExp0_4] +
   (BinExp5_9S`s'-BinExp5_9S0)*_b[BinExp5_9] +
   (BinExp10_14S`s'-BinExp10_14S0)*_b[BinExp10_14] +
   (BinExp15_19S`s'-BinExp15_19S0)*_b[BinExp15_19] +
   (BinExp20_24S`s'-BinExp20_24S0)*_b[BinExp20_24] +
   (BinExp25_29S`s'-BinExp25_29S0)*_b[BinExp25_29] +
   (BinExp30_infS`s'-BinExp30_infS0)*(_b[BinExp30_inf] + _b[1.lowP10#c.BinExp30_inf])) - 1, post;
  parmest, format(estimate min95 max95 %10.0g p %8.le) saving("$do\ImpactLP10S`s'.dta",replace);
};
use "$do\ImpactLP10S1.dta", clear;
gen shift = 1;
forvalues s = 2(1)3 {;
  append using "$do\ImpactLP10S`s'.dta";
  replace shift = `s' if shift == .;
};
gen model = "LP10"; 
forvalues s = 1(1)3 {;
  append using "$do\ImpactS`s'.dta";
  replace shift = `s' if shift == .;
};
replace model = "Preferred" if model == "";
gen b = 100*estimate;
gen hib=100*max95;
gen lwb=100*min95;
gen mod = 0;
replace mod = 1 if model == "Preferred";
replace mod = 2 if model == "LP10";
sort shift mod;
egen g = group(shift mod);
gen order = _n;
gen order2 = order;
replace order2 = order2 + 1 if order > 2;
replace order2 = order2 + 1 if order > 4;
graph twoway 
  (bar b order2 if g==1, color(red)) 
  (bar b order2 if g==2, color(orange)) 
  (bar b order2 if g==3, color(red)) 
  (bar b order2 if g==4, color(orange)) 
  (bar b order2 if g==5, color(red)) 
  (bar b order2 if g==6, color(orange)) 
  (rcap hib lwb order2),
  legend(order (1 "Preferred Model" 2 "Low Precipitation") pos(7) ring(0) rows(5) size(small) )
  xlabel(1.5 "+1{char 0176}C"  4.5 "+2{char 0176}C" 7.5 "+3{char 0176}C", noticks labsize(large))
  ytitle("yield impact (%)", size(large)) yscale() xtitle("warming scenario", size(large)) 
  graphregion(color(white)) ylabel(-40(10)0, nogrid labsize(large)) name(WarmingLowP);
graph export "$do\Figure4.emf", replace;


********************************************************************************;
*Figure 5 Mean Yields and the Effect of Heat across Wheat Cultivar Release Years;
********************************************************************************;
use "$dd\RegressionDataFinal.dta", clear;
by Cultivar year, sort: gen nvals = _n == 1;
by Cultivar: replace nvals = sum(nvals);
by Cultivar: replace nvals = nvals[_N];
keep if nvals >= 5;
save "$do\RegressionDataFinalN5.dta", replace;

*cultivar heat effects;
use "$do\RegressionDataFinalN5.dta", clear;
mixed ly i.cu i.si i.year prec c.prec#c.prec fr BinExp0_4 BinExp5_9 BinExp10_14 
         BinExp15_19 BinExp20_24 BinExp25_29 BinExp30_inf || cu: BinExp30_inf, noconstant ;
predict u0, reffects;
summ u*;
generate b_exp30_inf = _b[BinExp30_inf] + u0;
collapse (mean) b_exp30_inf, by(nvals cu Cultivar);
save "$do\CultivarHeatEffect.dta", replace;

*cultivar mean yields;
use "$do\RegressionDataFinalN5.dta", clear;
mixed ly i.cu i.si i.year prec c.prec#c.prec fr BinExp0_4 BinExp5_9 BinExp10_14 
         BinExp15_19 BinExp20_24 BinExp25_29 BinExp30_inf || cu: BinExp30_inf, noconstant ;
predict res, res;
gen eres = exp(res);
summ eres;
local Eeres = r(mean); 
summ cu; 
local min = r(min); 
local max = r(max);
margins, at(cu=(`min'(1)`max')) atmeans vsquish post coeflegend force nose;
gen lyhat = .;
replace lyhat =  _b[1bn._at] if cu == 1;
forvalues i = 2/`max' {;
  replace lyhat =  _b[`i'._at] if cu == `i';
};
gen yhat = exp(lyhat)*`Eeres';
collapse (mean) lyhat yhat, by(cu Cultivar);
save "$do\CultivarMeanYields.dta", replace;

*graphs;
use "$do\RegressionDataFinalN5.dta", clear;
merge m:1 Cultivar using "$do\CultivarHeatEffect.dta", nogen;
merge m:1 Cultivar using "$do\CultivarMeanYields.dta", nogen;
collapse (mean) b_exp30_inf yhat Release_Year, by(cu Cultivar);

twoway (scatter yhat Release_Year, m(o)) (lfit yhat Release_Year),
  legend(off) title(A, position(11) ring(0)) ytitle("mean yield") yscale() 
  xlabel(1987(5)2012,labsize(medsmall)) xtitle("release year") 
  graphregion(color(white)) ylabel(, nogrid) ;
graph save "$do\MeanRelease", replace;

twoway (scatter b_exp30_inf Release_Year, m(o)) (lfit b_exp30_inf Release_Year),
  legend(off) title(B, position(11) ring(0)) ytitle("heat effect") yscale() 
  xlabel(1987(5)2012,labsize(medsmall)) xtitle("release year") 
  graphregion(color(white)) ylabel(, nogrid) ;
graph save "$do\HeatRelease", replace;

gen n_heat = b_exp30_inf/yhat;
twoway (scatter n_heat Release_Year, m(o)) (lfit n_heat Release_Year),
  legend(off) title(C, position(11) ring(0)) ytitle("ratio: heat effect over mean yield") yscale() 
  xlabel(1987(5)2012,labsize(medsmall)) xtitle("release year") 
  graphregion(color(white)) ylabel(, nogrid) ;
graph save "$do\NormHeatRelease5", replace;

gr combine "$do\MeanRelease" "$do\HeatRelease" "$do\NormHeatRelease5", 
  graphregion(color(white)) row(1) imargin(0 2 0 0);
graph export "$do\Figure5.emf", replace;



********************************************************************************;
*Figure 6 Comparison of Warming Impacts on Wheat Yields for the Cultivars with the Largest and Smallest Heat Effects;
********************************************************************************;

use "$do\CultivarHeatEffect.dta", clear;
summ b* ;
local BinExp30_inf_large = r(min);
local BinExp30_inf_small = r(max);

*smallest heat effect;
use "$dd\WarmingScenarios.dta", clear;
forvalues s = 1/3 {;
  est use "$do\Regression";
  nlcom exp(
   (frS`s'-frS0)*_b[fr] +
   (BinExp0_4S`s'-BinExp0_4S0)*_b[BinExp0_4] +
   (BinExp5_9S`s'-BinExp5_9S0)*_b[BinExp5_9] +
   (BinExp10_14S`s'-BinExp10_14S0)*_b[BinExp10_14] +
   (BinExp15_19S`s'-BinExp15_19S0)*_b[BinExp15_19] +
   (BinExp20_24S`s'-BinExp20_24S0)*_b[BinExp20_24] +
   (BinExp25_29S`s'-BinExp25_29S0)*_b[BinExp25_29] +
   (BinExp30_infS`s'-BinExp30_infS0)*(`BinExp30_inf_small')) - 1, post;
  parmest, format(estimate min95 max95 %10.0g p %8.le) saving("$do\ImpactSmallS`s'.dta",replace);
};

*largest heat effect;
use "$dd\WarmingScenarios.dta", clear;
forvalues s = 1/3 {;
  est use "$do\Regression";
  nlcom exp(
   (frS`s'-frS0)*_b[fr] +
   (BinExp0_4S`s'-BinExp0_4S0)*_b[BinExp0_4] +
   (BinExp5_9S`s'-BinExp5_9S0)*_b[BinExp5_9] +
   (BinExp10_14S`s'-BinExp10_14S0)*_b[BinExp10_14] +
   (BinExp15_19S`s'-BinExp15_19S0)*_b[BinExp15_19] +
   (BinExp20_24S`s'-BinExp20_24S0)*_b[BinExp20_24] +
   (BinExp25_29S`s'-BinExp25_29S0)*_b[BinExp25_29] +
   (BinExp30_infS`s'-BinExp30_infS0)*(`BinExp30_inf_large')) - 1, post;
  parmest, format(estimate min95 max95 %10.0g p %8.le) saving("$do\ImpactLargeS`s'.dta",replace);
};

*graph;
use "$do\ImpactSmallS1.dta", clear;
gen shift = 1;
forvalues s = 2(1)3 {;
  append using "$do\ImpactSmallS`s'.dta";
  replace shift = `s' if shift == .;
};
gen model = "Small"; 
forvalues s = 1(1)3 {;
  append using "$do\ImpactLargeS`s'.dta";
  replace shift = `s' if shift == .;
};
replace model = "Large" if model == "";
gen b = 100*estimate;
gen hib=100*max95;
gen lwb=100*min95;
gen mod = 0;
replace mod = 1 if model == "Large";
replace mod = 2 if model == "Small";
sort shift mod;
egen g = group(shift mod);
gen order = _n;
gen order2 = order;
replace order2 = order2 + 1 if order > 2;
replace order2 = order2 + 1 if order > 4;
graph twoway 
  (bar b order2 if g==1, color(red)) 
  (bar b order2 if g==2, color(orange)) 
  (bar b order2 if g==3, color(red)) 
  (bar b order2 if g==4, color(orange)) 
  (bar b order2 if g==5, color(red)) 
  (bar b order2 if g==6, color(orange)) 
  (rcap hib lwb order2),
  legend(order (1 "Largest Heat Effect Cultivar" 2 "Smallest Heat Effect Cultivar") pos(7) ring(0) rows(2))
  xlabel(1.5 "+1{char 0176}C"  4.5 "+2{char 0176}C" 7.5 "+3{char 0176}C" , noticks labsize(large))
  ytitle("yield impact (%)", size(large)) yscale() xtitle("warming scenario", size(large)) 
  graphregion(color(white)) ylabel(, nogrid labsize(large));
graph export "$do\Figure6.emf", replace;


*END Code: for further inquiries regarding robustness checks, please contact the authors: ashew@astate.edu.