ClearAll["Global`*"];

n = 10000;
seed = 1;
SeedRandom[seed];

minPop = 1;
wp = 50;    

sampleInBins[min_?NumericQ, max_?NumericQ, n_Integer?Positive] := Module[
  {minp, maxp, delta},
  minp = SetPrecision[min, wp];
  maxp = SetPrecision[max, wp];
  If[!(NumericQ[minp] && NumericQ[maxp] && maxp >= minp), Return[$Failed]];
  If[maxp == minp, Return[ConstantArray[minp, n]]];
  delta = (maxp - minp)/n;
  RandomSample @ Table[
    RandomReal[{minp + (i - 1) delta, minp + i delta}, WorkingPrecision -> wp],
    {i, 1, n}]];

sampleInBinsLog[min_?NumericQ, max_?NumericQ, n_Integer?Positive] :=
  Power[10, sampleInBins[Log10[min], Log10[max], n]];

logParamNames = {"beta", "muV", "c", "gamma", "rL", "f", "ee", "TLt"};

rankVec[v_?VectorQ] := Module[{ord, sorted, groups, ranks, pos = 1, m, avg, idxs},
  ord = Ordering[v]; sorted = v[[ord]];
  groups = SplitBy[Transpose[{ord, sorted}], Last];
  ranks = ConstantArray[0., Length[v]];
  Do[m = Length[g]; avg = (pos + (pos + m - 1))/2.0;
    idxs = g[[All, 1]]; ranks[[idxs]] = avg; pos += m;, {g, groups}];
  ranks];

safeCorr[a_, b_] := Module[{sa = StandardDeviation[a], sb = StandardDeviation[b]},
  If[sa == 0 || sb == 0, Missing["Undefined"], Correlation[a, b]]];

residuals[design_?MatrixQ, vec_?VectorQ] := Module[{coef},
  coef = LeastSquares[design, vec]; vec - design . coef];

prccOne[XR_?MatrixQ, yR_?VectorQ, j_Integer] :=
  Module[{nobs, p, xj, keepCols, others, D, rx, ry},
    nobs = Length[yR]; p = Dimensions[XR][[2]];
    If[p < 2, Return[Missing["NeedAtLeast2Parameters"]]];
    xj = XR[[All, j]];
    If[StandardDeviation[xj] == 0, Return[Missing["ConstantParameter"]]];
    keepCols = Complement[Range[p], {j}];
    others = XR[[All, keepCols]];
    D = Join[ConstantArray[1, {nobs, 1}], others, 2];
    rx = residuals[D, xj]; ry = residuals[D, yR];
    safeCorr[rx, ry]];

cloudPlot[paramVec_?VectorQ, yVec_?VectorQ, name_, yLabel_, threshold_] :=
  Module[{pts, xmin, xmax, ep, opts},
    pts = Transpose[{paramVec, yVec}];
    xmin = Min[paramVec]; xmax = Max[paramVec];
    ep = If[threshold === None, {},
      {Red, Dashed, Line[{{xmin, threshold}, {xmax, threshold}}]}];
    opts = If[MemberQ[logParamNames, name],
      {ScalingFunctions -> {"Log10", Identity}}, {}];
    ListPlot[pts, PlotRange -> All, AxesLabel -> {name, yLabel},
      PlotLegends -> None, Epilog -> ep, ImageSize -> Large,
      Sequence @@ opts]];

LHSStats[Xraw_?MatrixQ, yraw_?VectorQ, paramNamesRaw_List, label_String,
    threshold_: 0, maxClouds_: 6] :=
  Module[{validIdx, y, X, paramNames, sdCols, varyIdx, Xv, namesV,
    spearman, spearmanTable, XR, yR, prcc, prccTable, score, order,
    idxCloud, clouds},
    validIdx = Select[Range[Length[yraw]], NumericQ[yraw[[#]]] &];
    If[validIdx === {},
      Print["[", label, "] No numeric y values. Skipping."]; Return[$Failed]];
    y = N @ yraw[[validIdx]]; X = N @ Xraw[[validIdx, All]];
    paramNames = paramNamesRaw;
    If[StandardDeviation[y] == 0,
      Print["\n===================================================="];
      Print["LHS STATS for y = ", label];
      Print["(Rows used = ", Length[y], ")"];
      Print["y is constant. Spearman/PRCC undefined. Skipping."];
      Return[<|"validIdx" -> validIdx, "paramNamesUsed" -> {},
        "Spearman" -> {}, "PRCC" -> {}|>]];
    sdCols = StandardDeviation /@ Transpose[X];
    varyIdx = Flatten @ Position[sdCols, _?(# > 0 &)];
    If[varyIdx === {},
      Print["[", label, "] All params constant. Skipping."]; Return[$Failed]];
    Xv = X[[All, varyIdx]]; namesV = paramNames[[varyIdx]];
    spearman = Table[safeCorr[rankVec[Xv[[All, j]]], rankVec[y]],
      {j, 1, Length[namesV]}];
    spearman = N @ spearman;
    spearmanTable = Transpose[{namesV, spearman}];
    Print["\n===================================================="];
    Print["LHS STATS for y = ", label];
    Print["(Rows used = ", Length[y], ", Parameters used = ", Length[namesV], ")"];
    If[threshold =!= None, Print["Cloud threshold line at y = ", threshold, "."]];
    Print["Spearman rank correlations (parameter vs y):"];
    Print[Grid[Prepend[spearmanTable, {"Parameter", "SpearmanRho"}], Frame -> All]];
    XR = Transpose[rankVec /@ Transpose[Xv]]; yR = rankVec[y];
    prcc = Table[prccOne[XR, yR, j], {j, 1, Length[namesV]}];
    prcc = N @ prcc;
    prccTable = Transpose[{namesV, prcc}];
    Print["PRCC (parameter vs y, controlling for all other parameters):"];
    Print[Grid[Prepend[prccTable, {"Parameter", "PRCC"}], Frame -> All]];
    score = Abs[prcc] /. Missing[_] -> -Infinity;
    order = Reverse @ Ordering[score];
    idxCloud = If[maxClouds === All, Range[Length[namesV]],
      Take[order, UpTo[maxClouds]]];
    clouds = Table[cloudPlot[Xv[[All, j]], y, namesV[[j]], label, threshold],
      {j, idxCloud}];
    Print["Cloud plots (showing ", Length[idxCloud], "):"];
    Print /@ clouds;
    <|"validIdx" -> validIdx, "paramNamesUsed" -> namesV,
      "Spearman" -> spearman, "PRCC" -> prcc|>];

rMin = 7.98*10^-3;    rMax = 4.00*10^-2;          
KMin = 10^9;          KMax = 10^9;
betaMin = 5.67*10^-12; betaMax = 1.00*10^-8;       
muVMin = 6.94*10^-4;  muVMax = 1.84*10^-2;         
cMin = 7.80*10^-6;    cMax = 3.33*10^-3;           

alphaMin = 7.65*10^-3; alphaMax = 3.16*10^-2;      
gammaMin = 5;         gammaMax = 500;              

beMin = 2.55*10^-2;   beMax = 1.02*10^-1;          
omegaMin = 0.05;      omegaMax = 0.27;             
riMin = 3.99*10^-3;   riMax = 6.00*10^-2;          

rLMin = 7.98*10^-4;   rLMax = 5.33*10^-2;          
fMin = 1.00*10^-3;    fMax = 9.00*10^-1;           
eeMin = 1.00*10^-9;   eeMax = 1.00*10^-4;          

TLtMin = 10;          TLtMax = 100;                
TeMin = 5;            TeMax = 20;                  

r     = sampleInBins[rMin, rMax, n];
K     = sampleInBins[KMin, KMax, n];
beta  = sampleInBinsLog[betaMin, betaMax, n];
muV   = sampleInBinsLog[muVMin, muVMax, n];
c     = sampleInBinsLog[cMin, cMax, n];
alpha = sampleInBins[alphaMin, alphaMax, n];
gamma = sampleInBinsLog[gammaMin, gammaMax, n];
be    = sampleInBins[beMin, beMax, n];
omega = sampleInBins[omegaMin, omegaMax, n];
ri    = sampleInBins[riMin, riMax, n];       
rL    = sampleInBinsLog[rLMin, rLMax, n];
f     = sampleInBinsLog[fMin, fMax, n];
ee    = sampleInBinsLog[eeMin, eeMax, n];

TLt   = sampleInBinsLog[TLtMin, TLtMax, n];
Te    = sampleInBins[TeMin, TeMax, n];

TL    = 1/alpha;    
Tat   = 1/alpha;    

alphaEffL = alpha;
SstarLAll = muV*(alphaEffL + c)/(beta*gamma*alphaEffL);
IstarLAll = (r - c - r*SstarLAll/K)/(r/K + beta*gamma*alphaEffL/muV);
VstarLAll = gamma*alphaEffL*IstarLAll/muV;
feasLAll  = MapThread[Boole[#1 > 0 && #2 > 0 && #3 > 0] &,
              {SstarLAll, IstarLAll, VstarLAll}];

chronicEquilAll = Table[
  Module[{qq, aa, Aq, Bq, Cq, roots, Svalid, Sc, xval, Ic, Vc, feas},
    qq = beta[[ii]]*be[[ii]]/muV[[ii]];
    aa = (1 - omega[[ii]])*ri[[ii]];
    Aq = qq*(K[[ii]]*qq*ri[[ii]] + aa*(r[[ii]] - ri[[ii]]));
    Bq = K[[ii]]*aa*qq*(c[[ii]] + ri[[ii]]) - 2*K[[ii]]*c[[ii]]*qq*ri[[ii]]
         + aa*c[[ii]]*(ri[[ii]] - r[[ii]]);
    Cq = K[[ii]]*c[[ii]]*(aa - c[[ii]])*(aa - ri[[ii]]);
    roots = S /. Solve[Aq*S^2 + Bq*S + Cq == 0, S];
    roots = Select[N[roots], Im[#] == 0 &];
    roots = Re /@ roots;
    Svalid = Select[roots, (# > 0 && 0 < (c[[ii]] - qq*#)/aa < 1) &];
    If[Length[Svalid] == 0,
      {0., 0., 0., 0., 0},
      Sc = Min[Svalid];
      xval = (c[[ii]] - qq*Sc)/aa;   
      Ic = K[[ii]]*(1 - xval) - Sc;
      Vc = be[[ii]]*Ic/muV[[ii]];
      feas = Boole[Sc > 0 && Ic > 0 && Vc > 0 && 0 < xval < 1];
      {Sc, Ic, Vc, xval, feas}]],
  {ii, 1, n}];

SstarCAll  = chronicEquilAll[[All, 1]];
IstarCAll  = chronicEquilAll[[All, 2]];
VstarCAll  = chronicEquilAll[[All, 3]];
LogFacCAll = chronicEquilAll[[All, 4]];
feasCAll   = chronicEquilAll[[All, 5]];

temperateEquilAll = Table[
  Module[{r0, K0, c0, beta0, alpha0, gamma0, muV0, f0, ee0, rl0,
          aa, q, cc, d, Smax, A0, B0, Khost, sCand, sGood, result,
          numI, denI, poly},
    r0 = r[[ii]]; K0 = K[[ii]]; c0 = c[[ii]];
    beta0 = beta[[ii]]; alpha0 = alpha[[ii]];  
    gamma0 = gamma[[ii]]; muV0 = muV[[ii]];
    f0 = f[[ii]]; ee0 = ee[[ii]]; rl0 = rL[[ii]];

    aa = alpha0 + c0;
    q  = (alpha0 * gamma0) / muV0;
    cc = beta0 * q;
    d  = (1 - f0) * cc;

    result = {0., 0., 0., 0., 0};

    If[d != 0,
      Smax  = aa / d;
      A0    = (1/cc) * (r0 * (c0 + ee0)/rl0 - c0);
      B0    = (r0 * f0 * ee0) / rl0;
      Khost = K0 * (1 - c0/r0);

      If[Khost > 0,
        numI = A0 * aa - S * (A0 * d + B0);
        denI = aa - d * S;
        poly = Expand[(1 + (K0 * cc) / r0) * numI + denI * (S + numI / ee0 - Khost)];

        Quiet[sCand = S /. NSolve[poly == 0, S, WorkingPrecision -> wp]];
        sGood = Select[sCand, (NumericQ[#] && Abs[Im[#]] < 10^-20 && 0 < Re[#] < Smax) &];
        sGood = Re /@ sGood;

        Do[
          If[result[[5]] == 0,
            Module[{s0, denv, Iv, Lv, Vv},
              s0 = sg;
              denv = aa - d * s0;
              If[denv > 0,
                Iv = (A0 * aa - s0 * (A0 * d + B0)) / denv;
                If[Iv > 0,
                  Lv = denv / ee0 * Iv;
                  Vv = (alpha0 * gamma0 / muV0) * Iv;
                  If[s0 > minPop && Lv > minPop && Iv > minPop && Vv > minPop,
                    result = {s0, Lv, Iv, Vv, 1}]]]]],
          {sg, sGood}]]];
    result],
  {ii, 1, n}];

SstarTAll = temperateEquilAll[[All, 1]];
LstarTAll = temperateEquilAll[[All, 2]];
IstarTAll = temperateEquilAll[[All, 3]];
VstarTAll = temperateEquilAll[[All, 4]];
feasTAll  = temperateEquilAll[[All, 5]];

LLstarLAll = 1 - (SstarLAll + IstarLAll)/K;
LLstarCAll = 1 - (SstarCAll + IstarCAll)/K;
LLstarTAll = 1 - (SstarTAll + LstarTAll + IstarTAll)/K;

Print["Fraction feasible lytic resident:     ", N[Mean[feasLAll]]];
Print["Fraction feasible chronic resident:   ", N[Mean[feasCAll]]];
Print["Fraction feasible temperate resident: ", N[Mean[feasTAll]]];

matChronicInv = {
  {lam + mc - (1 - om)*rri*LL,    -bt*SS*Exp[-mc*Tdel]*(1 - lam*Tdel)},
  {-bu,                             lam + mv}};

matLyticInv = {
  {lam + mc,  -bt*SS + bt*SS*Exp[-mc*Tdel]*(1 - lam*Tdel)},
  {0,          lam + mv - gm*bt*SS*Exp[-mc*Tdel]*(1 - lam*Tdel)}};

matTemperateInv = {
  {lam - rLs*gg + (mc + epi),
   -ff*bt*SS*Exp[-mc*Tdl]*(1 - lam*Tdl)},
  {-gm*epi*Exp[-mc*Tda]*(1 - lam*Tda),
   lam + mv - gm*(1 - ff)*bt*SS*Exp[-mc*Tda]*(1 - lam*Tda)}};

Print["\nSymbolic Chronic-invader roots:"];
Print[Simplify /@ (lam /. Solve[Det[matChronicInv] == 0, lam])];
Print["\nSymbolic Lytic-invader roots:"];
Print[Simplify /@ (lam /. Solve[Det[matLyticInv] == 0, lam])];
Print["\nSymbolic Temperate-invader roots:"];
Print[Simplify /@ (lam /. Solve[Det[matTemperateInv] == 0, lam])];

invasionLambda[mat_, subs_] :=
  Max[Re[N[lam /. Solve[Det[mat /. subs] == 0, lam]]]];
  
  invasionLambdaLyticReal[bt_, SS_, mc_, Tdel_, gm_, mv_] :=
  Max[Re[N[lam /. Solve[
    lam + mv - gm * bt * SS * Exp[-mc * Tdel] * (1 - lam * Tdel) == 0,
    lam
  ]]]];

lambdaCinvLAll = Table[
  If[feasLAll[[ii]] == 0,
    Missing["InfeasibleLyticResident"],
    invasionLambda[matChronicInv,
      {om -> omega[[ii]], rri -> ri[[ii]],
       bt -> beta[[ii]], SS -> SstarLAll[[ii]], mc -> c[[ii]],
       Tdel -> Te[[ii]], bu -> be[[ii]], mv -> muV[[ii]],
       LL -> LLstarLAll[[ii]]}]],
  {ii, 1, n}];

lambdaCinvTAll = Table[
  If[feasTAll[[ii]] == 0,
    Missing["InfeasibleTemperateResident"],
    invasionLambda[matChronicInv,
      {om -> omega[[ii]], rri -> ri[[ii]],
       bt -> beta[[ii]], SS -> SstarTAll[[ii]], mc -> c[[ii]],
       Tdel -> Te[[ii]], bu -> be[[ii]], mv -> muV[[ii]],
       LL -> LLstarTAll[[ii]]}]],
  {ii, 1, n}];

lambdaLinvCAll = Table[
  If[feasCAll[[ii]] == 0,
    Missing["InfeasibleChronicResident"],
    invasionLambdaLyticReal[
      beta[[ii]], SstarCAll[[ii]], c[[ii]],
      TL[[ii]], gamma[[ii]], muV[[ii]]
    ]
  ],
  {ii, 1, n}
];

  lambdaLinvTAll = Table[
  If[feasTAll[[ii]] == 0,
    Missing["InfeasibleTemperateResident"],
    invasionLambdaLyticReal[
      beta[[ii]], SstarTAll[[ii]], c[[ii]],
      TL[[ii]], gamma[[ii]], muV[[ii]]
    ]
  ],
  {ii, 1, n}
];

lambdaTinvLAll = Table[
  If[feasLAll[[ii]] == 0,
    Missing["InfeasibleLyticResident"],
    invasionLambda[matTemperateInv,
      {rLs -> rL[[ii]], gg -> LLstarLAll[[ii]], mc -> c[[ii]],
       epi -> ee[[ii]], ff -> f[[ii]], bt -> beta[[ii]],
       SS -> SstarLAll[[ii]], gm -> gamma[[ii]], mv -> muV[[ii]],
       Tdl -> TLt[[ii]], Tda -> Tat[[ii]]}]],
  {ii, 1, n}];

lambdaTinvCAll = Table[
  If[feasCAll[[ii]] == 0,
    Missing["InfeasibleChronicResident"],
    invasionLambda[matTemperateInv,
      {rLs -> rL[[ii]], gg -> LLstarCAll[[ii]], mc -> c[[ii]],
       epi -> ee[[ii]], ff -> f[[ii]], bt -> beta[[ii]],
       SS -> SstarCAll[[ii]], gm -> gamma[[ii]], mv -> muV[[ii]],
       Tdl -> TLt[[ii]], Tda -> Tat[[ii]]}]],
  {ii, 1, n}];

invBool[x_] := If[MissingQ[x], Missing[], If[x > 0, 1, 0]];

invCinvL = invBool /@ lambdaCinvLAll;
invCinvT = invBool /@ lambdaCinvTAll;
invLinvC = invBool /@ lambdaLinvCAll;
invLinvT = invBool /@ lambdaLinvTAll;
invTinvL = invBool /@ lambdaTinvLAll;
invTinvC = invBool /@ lambdaTinvCAll;

classify[fwd_, rev_, xWins_, yWins_, coexist_] :=
  If[MissingQ[fwd] || MissingQ[rev], "NA",
    Which[
      fwd == 1 && rev == 1, coexist,
      fwd == 1 && rev == 0, xWins,
      fwd == 0 && rev == 1, yWins,
      fwd == 0 && rev == 0, "bistable",
      True, "other"]];

outcomeCL = Table[classify[invCinvL[[ii]], invLinvC[[ii]], "C_wins", "L_wins", "LC_coexist"], {ii, 1, n}];
outcomeTL = Table[classify[invTinvL[[ii]], invLinvT[[ii]], "T_wins", "L_wins", "TL_coexist"], {ii, 1, n}];
outcomeCT = Table[classify[invCinvT[[ii]], invTinvC[[ii]], "C_wins", "T_wins", "CT_coexist"], {ii, 1, n}];

paramNamesAll = {"r", "K", "beta", "muV", "c", "alpha", "gamma",
                 "be", "omega", "ri",
                 "rL", "f", "ee",
                 "TLt", "Te"};
Xall = N @ Transpose[{r, K, beta, muV, c, alpha, gamma, be, omega, ri,
                       rL, f, ee, TLt, Te}];

LHSStats[Xall, lambdaCinvLAll, paramNamesAll, "lambda_CinvL", 0, 8];
LHSStats[Xall, lambdaLinvCAll, paramNamesAll, "lambda_LinvC", 0, 8];
LHSStats[Xall, lambdaCinvTAll, paramNamesAll, "lambda_CinvT", 0, 8];
LHSStats[Xall, lambdaLinvTAll, paramNamesAll, "lambda_LinvT", 0, 8];
LHSStats[Xall, lambdaTinvLAll, paramNamesAll, "lambda_TinvL", 0, 8];
LHSStats[Xall, lambdaTinvCAll, paramNamesAll, "lambda_TinvC", 0, 8];

safeCSV[x_] := If[MissingQ[x], "NA", x];

csvHeader = {"draw",
             "r", "K", "beta", "muV", "c", "alpha", "gamma",
             "be", "omega", "ri",
             "rL", "f", "ee",
             "TL", "Tat", "TLt", "Te",
             "SstarL", "IstarL", "VstarL", "feasL",
             "SstarC", "IstarC", "VstarC", "LogFacC", "feasC",
             "SstarT", "LstarT", "IstarT", "VstarT", "feasT",
             "lambda_CinvL", "lambda_LinvC",
             "lambda_CinvT", "lambda_LinvT",
             "lambda_TinvL", "lambda_TinvC",
             "inv_CinvL", "inv_LinvC",
             "inv_CinvT", "inv_LinvT",
             "inv_TinvL", "inv_TinvC",
             "outcome_CL", "outcome_TL", "outcome_CT"};

csvRows = Table[
  {ii, r[[ii]], K[[ii]], beta[[ii]], muV[[ii]], c[[ii]], alpha[[ii]], gamma[[ii]],
   be[[ii]], omega[[ii]], ri[[ii]],
   rL[[ii]], f[[ii]], ee[[ii]],
   TL[[ii]], Tat[[ii]], TLt[[ii]], Te[[ii]],
   SstarLAll[[ii]], IstarLAll[[ii]], VstarLAll[[ii]], feasLAll[[ii]],
   SstarCAll[[ii]], IstarCAll[[ii]], VstarCAll[[ii]], LogFacCAll[[ii]], feasCAll[[ii]],
   SstarTAll[[ii]], LstarTAll[[ii]], IstarTAll[[ii]], VstarTAll[[ii]], feasTAll[[ii]],
   safeCSV[lambdaCinvLAll[[ii]]], safeCSV[lambdaLinvCAll[[ii]]],
   safeCSV[lambdaCinvTAll[[ii]]], safeCSV[lambdaLinvTAll[[ii]]],
   safeCSV[lambdaTinvLAll[[ii]]], safeCSV[lambdaTinvCAll[[ii]]],
   safeCSV[invCinvL[[ii]]], safeCSV[invLinvC[[ii]]],
   safeCSV[invCinvT[[ii]]], safeCSV[invLinvT[[ii]]],
   safeCSV[invTinvL[[ii]]], safeCSV[invTinvC[[ii]]],
   outcomeCL[[ii]], outcomeTL[[ii]], outcomeCT[[ii]]},
  {ii, 1, n}];

csvPath = FileNameJoin[{Directory[], "Mutual_invasibility_LCT_delays_V5matched_results.csv"}];
Export[csvPath, Prepend[csvRows, csvHeader], "CSV"];
Print["\nCSV written to: ", csvPath];

Print["\n===================================================="];
Print["Mutual invasibility LCT delays V5matched - summary"];
Print["===================================================="];
Print["n                              = ", n];
Print["Fraction feasible lytic        = ", N[Mean[feasLAll]]];
Print["Fraction feasible chronic      = ", N[Mean[feasCAll]]];
Print["Fraction feasible temperate    = ", N[Mean[feasTAll]]];

Print["\nOutcome CL (chronic vs lytic):"];
Print[Grid[Prepend[Tally[outcomeCL], {"Outcome", "Count"}], Frame -> All]];
Print["\nOutcome TL (temperate vs lytic):"];
Print[Grid[Prepend[Tally[outcomeTL], {"Outcome", "Count"}], Frame -> All]];
Print["\nOutcome CT (chronic vs temperate):"];
Print[Grid[Prepend[Tally[outcomeCT], {"Outcome", "Count"}], Frame -> All]];

Print["\nAlpha vs. T ranges (matched mean latency):"];
Print["  alpha  = [", alphaMin, ", ", alphaMax, "] /min (V5 convention, linear)"];
Print["  TL=Tat = 1/alpha -> [", N[1/alphaMax], ", ", N[1/alphaMin], "] min (derived)"];
Print["\nDDE-only delay ranges (no V5 analog):"];
Print["  TLt (temperate lys commitment) = [", TLtMin, ", ", TLtMax, "] min, log-sampled"];
Print["  Te  (chronic eclipse)          = [", TeMin, ", ", TeMax, "] min, linear-sampled"];
