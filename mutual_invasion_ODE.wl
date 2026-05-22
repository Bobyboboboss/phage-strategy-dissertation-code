ClearAll["Global`*"];

n = 10000;                  
wp = 50;                    
seed = 1;                   
SeedRandom[seed];
tolStab = 10^-10;          
minPop = 1;                

ClearAll[sampleInBins];
sampleInBins[min_?NumericQ, max_?NumericQ, n_Integer?Positive] := Module[
  {minp, maxp, \[CapitalDelta]},
  minp = SetPrecision[min, wp];
  maxp = SetPrecision[max, wp];
  If[!(NumericQ[minp] && NumericQ[maxp] && maxp >= minp), Return[$Failed]];
  If[maxp == minp, Return[ConstantArray[minp, n]]];
  \[CapitalDelta] = (maxp - minp)/n;
  RandomSample@Table[
    RandomReal[{minp + (i - 1) \[CapitalDelta], minp + i \[CapitalDelta]},
      WorkingPrecision -> wp],
    {i, 1, n}
  ]
];

ClearAll[sampleInBinsLog];
sampleInBinsLog[min_?NumericQ, max_?NumericQ, n_Integer?Positive] :=
  Power[10, sampleInBins[Log10[min], Log10[max], n]];

ClearAll[rankVec, safeCorr, residuals, prccOne, cloudPlot];

rankVec[v_?VectorQ] := Module[
  {ord, sorted, groups, ranks, pos = 1, m, avg},
  ord = Ordering[v];
  sorted = v[[ord]];
  groups = SplitBy[Transpose[{ord, sorted}], Last];
  ranks = ConstantArray[0., Length[v]];
  Do[
    m = Length[g];
    avg = (pos + (pos + m - 1)) / 2.;
    Do[ranks[[g[[j, 1]]]] = avg, {j, 1, m}];
    pos += m,
    {g, groups}
  ];
  ranks
];

safeCorr[a_, b_] := Module[{sa, sb},
  sa = StandardDeviation[N[a]];
  sb = StandardDeviation[N[b]];
  If[sa == 0 || sb == 0, Missing["Undefined"], Correlation[N[a], N[b]]]
];

residuals[design_?MatrixQ, vec_?VectorQ] := Module[{coef},
  coef = LeastSquares[design, vec];
  vec - design . coef
];

prccOne[XR_?MatrixQ, yR_?VectorQ, j_Integer] := Module[
  {nobs, p, xj, keepCols, others, D, rx, ry},
  nobs = Length[yR]; p = Dimensions[XR][[2]];
  xj = XR[[All, j]];
  If[StandardDeviation[N[xj]] == 0, Return[Missing["ConstantParameter"]]];
  keepCols = Complement[Range[p], {j}];
  others = XR[[All, keepCols]];
  D = Join[ConstantArray[1, {nobs, 1}], others, 2];
  rx = residuals[D, xj];
  ry = residuals[D, yR];
  safeCorr[rx, ry]
];

cloudPlot[paramVec_?VectorQ, yVec_?VectorQ, name_String, threshold_: 0] := Module[
  {pts, xmin, xmax},
  pts = Transpose[{paramVec, yVec}];
  xmin = Min[paramVec]; xmax = Max[paramVec];
  ListPlot[pts,
    PlotRange -> All,
    AxesLabel -> {name, Subscript[\[Lambda], "max"]},
    Epilog -> {Red, Dashed, Thick, Line[{{xmin, threshold}, {xmax, threshold}}]},
    ImageSize -> Large,
    PlotStyle -> {PointSize[Small], Opacity[0.6]},
    PlotLabel -> name <> " vs \[Lambda]_max"
  ]
];

rMin = 7.98*10^-3;           rMax = 4.00*10^-2;
KMin = 10^9;                 KMax = 10^9;         

betaMin = 5.67*10^-12;       betaMax = 1.00*10^-8;
muVMin = 6.94*10^-4;         muVMax = 1.84*10^-2;
cMin = 7.80*10^-6;           cMax = 3.33*10^-3;

alphaMin = 7.65*10^-3;       alphaMax = 3.16*10^-2;
gammaMin = 5;                gammaMax = 500;

beMin = 2.55*10^-2;          beMax = 1.02*10^-1;
omegaMin = 0.05;             omegaMax = 0.27;
riMin = 3.99*10^-3;          riMax = 6.00*10^-2;

rlMin = 7.98*10^-4;          rlMax = 5.33*10^-2;
fMin = 0.001;                fMax = 0.9;
eeMin = 10^-9;               eeMax = 10^-4;

Print["Sampling ", n, " parameter sets via LHS..."];
Print["V5: Using log-scale sampling for beta, muV, c, gamma, rl, f, ee"];

rSamp     = sampleInBins[rMin, rMax, n];            
KSamp     = sampleInBins[KMin, KMax, n];            
betaSamp  = sampleInBinsLog[betaMin, betaMax, n];   
muVSamp   = sampleInBinsLog[muVMin, muVMax, n];     
cSamp     = sampleInBinsLog[cMin, cMax, n];          

alphaSamp = sampleInBins[alphaMin, alphaMax, n];     
gammaSamp = sampleInBinsLog[gammaMin, gammaMax, n];  

beSamp    = sampleInBins[beMin, beMax, n];           
omegaSamp = sampleInBins[omegaMin, omegaMax, n];    
riSamp    = sampleInBins[riMin, riMax, n];           

rlSamp    = sampleInBinsLog[rlMin, rlMax, n];        
fSamp     = sampleInBinsLog[fMin, fMax, n];          
eeSamp    = sampleInBinsLog[eeMin, eeMax, n];        

Print["Sampling complete."];

lyticEqForDraw[ii_] := Module[
  {r, K, c, beta, alpha, gamma, muV, Sstar, Istar, Vstar, R0},

  r = rSamp[[ii]]; K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; alpha = alphaSamp[[ii]];
  gamma = gammaSamp[[ii]]; muV = muVSamp[[ii]];

  R0 = beta * alpha * gamma * K * (1 - c/r) / (muV * (alpha + c));
  Sstar = muV * (alpha + c) / (beta * alpha * gamma);
  Istar = Sstar * (r * (1 - Sstar/K) - c) / (r * Sstar/K + (alpha + c));
  Vstar = gamma * alpha * Istar / muV;

  If[Sstar > 0 && Istar > 0 && Vstar > 0,
    <|"S" -> Sstar, "I" -> Istar, "V" -> Vstar, "R0" -> R0, "feasible" -> True|>,
    <|"S" -> Null, "I" -> Null, "V" -> Null, "R0" -> R0, "feasible" -> False|>
  ]
];

chronicEqForDraw[ii_] := Module[
  {r, K, c, beta, muV, be, omega, ri,
   q, a, Aq, Bq, Cq, delta, rootsS, makeEq, feasibleQ, eqs},

  r = rSamp[[ii]]; K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; muV = muVSamp[[ii]];
  be = beSamp[[ii]]; omega = omegaSamp[[ii]]; ri = riSamp[[ii]];

  q = beta * be / muV;
  a = (1 - omega) * ri;

  If[!(NumericQ[a] && a > 0 && NumericQ[q] && q != 0),
    Return[<|"S" -> Null, "I" -> Null, "V" -> Null, "feasible" -> False|>]
  ];

  Aq = q * (K * q * ri + a * (r - ri));
  Bq = K * a * q * (c + ri) - 2 * K * c * q * ri + a^2 * (c - c) + a * c * (ri - r);
  Cq = K * c * (a - c) * (a - ri);

  delta = Bq^2 - 4 * Aq * Cq;
  If[!(NumericQ[delta] && delta >= 0 && Aq != 0),
    Return[<|"S" -> Null, "I" -> Null, "V" -> Null, "feasible" -> False|>]
  ];

  rootsS = N[(-Bq + {-Sqrt[delta], +Sqrt[delta]}) / (2 * Aq), wp];

  makeEq[sStar_] := Module[{Istar, Vstar},
    Istar = K * (1 - (c - q * sStar) / a) - sStar;
    Vstar = (be / muV) * Istar;
    {sStar, Istar, Vstar}
  ];

  feasibleQ[{sStar_, Istar_, Vstar_}] := Module[{x = c - q * sStar},
    (sStar > 0) && (Istar > 0) && (Vstar > 0) && (0 < x < a)
  ];

  eqs = Select[makeEq /@ rootsS, feasibleQ];
  If[Length[eqs] == 0,
    <|"S" -> Null, "I" -> Null, "V" -> Null, "feasible" -> False|>,
    <|"S" -> eqs[[1, 1]], "I" -> eqs[[1, 2]], "V" -> eqs[[1, 3]], "feasible" -> True|>
  ]
];

temperateEqForDraw[ii_] := Module[
  {r0, K0, c0, beta0, alpha0, gamma0, muV0, f0, ee0, rl0,
   a, q, cc, d, Smax, A0, B0, Khost, S,
   numI, denI, poly, sCand, sGood, result},

  r0 = rSamp[[ii]]; K0 = KSamp[[ii]]; c0 = cSamp[[ii]];
  beta0 = betaSamp[[ii]]; alpha0 = alphaSamp[[ii]];
  gamma0 = gammaSamp[[ii]]; muV0 = muVSamp[[ii]];
  f0 = fSamp[[ii]]; ee0 = eeSamp[[ii]]; rl0 = rlSamp[[ii]];

  a = alpha0 + c0;
  q = (alpha0 * gamma0) / muV0;
  cc = beta0 * q;
  d = (1 - f0) * cc;

  If[d == 0, Return[<|"S" -> Null, "L" -> Null, "I" -> Null, "V" -> Null, "feasible" -> False|>]];

  Smax = a / d;
  A0 = (1/cc) * (r0 * (c0 + ee0)/rl0 - c0);
  B0 = (r0 * f0 * ee0) / rl0;
  Khost = K0 * (1 - c0/r0);

  If[Khost <= 0, Return[<|"S" -> Null, "L" -> Null, "I" -> Null, "V" -> Null, "feasible" -> False|>]];

  numI = A0 * a - S * (A0 * d + B0);
  denI = a - d * S;

  poly = Expand[
    (1 + (K0 * cc) / r0) * numI + denI * (S + numI / ee0 - Khost)
  ];

  Quiet[
    sCand = S /. NSolve[poly == 0, S, WorkingPrecision -> Min[wp, 30]];
  ];
  sGood = Select[sCand, (NumericQ[#] && Abs[Im[#]] < 10^-20 && 0 < Re[#] < Smax) &];
  sGood = Re /@ sGood;

  result = <|"S" -> Null, "L" -> Null, "I" -> Null, "V" -> Null, "feasible" -> False|>;

  Do[
    Module[{s0, denv, Iv, Lv, Vv},
      s0 = sg;
      denv = a - d * s0;
      If[denv <= 0, Continue[]];
      Iv = (A0 * a - s0 * (A0 * d + B0)) / denv;
      If[Iv <= 0, Continue[]];
      Lv = denv / ee0 * Iv;
      Vv = (alpha0 * gamma0 / muV0) * Iv;
      If[s0 > minPop && Lv > minPop && Iv > minPop && Vv > minPop,
        result = <|"S" -> s0, "L" -> Lv, "I" -> Iv, "V" -> Vv, "feasible" -> True|>;
        Break[];
      ];
    ],
    {sg, sGood}
  ];

  result
];

lyticStableForDraw[ii_, eq_Association] := Module[
  {r, K, c, beta, alpha, gamma, muV, S, I, V,
   sV, iV, vV, Ls, fSs, fIs, fVs, J, eigs},

  r = rSamp[[ii]]; K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; alpha = alphaSamp[[ii]];
  gamma = gammaSamp[[ii]]; muV = muVSamp[[ii]];
  {S, I, V} = eq /@ {"S", "I", "V"};

  Ls = 1 - (sV + iV)/K;
  fSs = r * sV * Ls - beta * sV * vV - c * sV;
  fIs = beta * sV * vV - (alpha + c) * iV;
  fVs = gamma * alpha * iV - muV * vV;

  J = N[D[{fSs, fIs, fVs}, {{sV, iV, vV}}] /. {sV -> S, iV -> I, vV -> V}, wp];
  eigs = Eigenvalues[J];
  Max[Re[eigs]] < -tolStab
];

chronicStableForDraw[ii_, eq_Association] := Module[
  {r, K, c, beta, muV, be, omega, ri, S, Ic, Vc,
   sV, iV, vV, Ls, fSs, fIs, fVs, J, eigs},

  r = rSamp[[ii]]; K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; muV = muVSamp[[ii]];
  be = beSamp[[ii]]; omega = omegaSamp[[ii]]; ri = riSamp[[ii]];
  {S, Ic, Vc} = eq /@ {"S", "I", "V"};

  Ls = 1 - (sV + iV)/K;
  fSs = r * sV * Ls - beta * sV * vV - c * sV + omega * ri * iV * Ls;
  fIs = beta * sV * vV - c * iV + (1 - omega) * ri * iV * Ls;
  fVs = be * iV - muV * vV;

  J = N[D[{fSs, fIs, fVs}, {{sV, iV, vV}}] /. {sV -> S, iV -> Ic, vV -> Vc}, wp];
  eigs = Eigenvalues[J];
  Max[Re[eigs]] < -tolStab
];

temperateStableForDraw[ii_, eq_Association] := Module[
  {r, K, c, beta, alpha, gamma, muV, f, ee, rl,
   S, L, It, Vt,
   sV, lV, iV, vV, Ns, gs, fSs, fLs, fIs, fVs, J, eigs},

  r = rSamp[[ii]]; K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; alpha = alphaSamp[[ii]];
  gamma = gammaSamp[[ii]]; muV = muVSamp[[ii]];
  f = fSamp[[ii]]; ee = eeSamp[[ii]]; rl = rlSamp[[ii]];
  {S, L, It, Vt} = eq /@ {"S", "L", "I", "V"};

  Ns = sV + lV + iV;
  gs = 1 - Ns/K;
  fSs = r * sV * gs - c * sV - beta * sV * vV;
  fLs = f * beta * sV * vV - (c + ee) * lV + rl * lV * gs;
  fIs = (1 - f) * beta * sV * vV + ee * lV - (alpha + c) * iV;
  fVs = alpha * gamma * iV - muV * vV;

  J = N[D[{fSs, fLs, fIs, fVs}, {{sV, lV, iV, vV}}] /.
    {sV -> S, lV -> L, iV -> It, vV -> Vt}, wp];
  eigs = Eigenvalues[J];
  Max[Re[eigs]] < -tolStab
];

temperateMaxReEigForDraw[ii_, eq_Association] := Module[
  {r, K, c, beta, alpha, gamma, muV, f, ee, rl,
   S, L, It, Vt,
   sV, lV, iV, vV, Ns, gs, fSs, fLs, fIs, fVs, J, eigs},

  r = rSamp[[ii]]; K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; alpha = alphaSamp[[ii]];
  gamma = gammaSamp[[ii]]; muV = muVSamp[[ii]];
  f = fSamp[[ii]]; ee = eeSamp[[ii]]; rl = rlSamp[[ii]];
  {S, L, It, Vt} = eq /@ {"S", "L", "I", "V"};

  Ns = sV + lV + iV;
  gs = 1 - Ns/K;
  fSs = r * sV * gs - c * sV - beta * sV * vV;
  fLs = f * beta * sV * vV - (c + ee) * lV + rl * lV * gs;
  fIs = (1 - f) * beta * sV * vV + ee * lV - (alpha + c) * iV;
  fVs = alpha * gamma * iV - muV * vV;

  J = N[D[{fSs, fLs, fIs, fVs}, {{sV, lV, iV, vV}}] /.
    {sV -> S, lV -> L, iV -> It, vV -> Vt}, wp];
  eigs = Eigenvalues[J];
  Max[Re[eigs]]
];

chronicInvadesLyticEig[ii_, lyticEq_Association] := Module[
  {r, K, c, beta, muV, be, omega, ri, S, I, Lterm, J, eigs},

  r = rSamp[[ii]]; K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; muV = muVSamp[[ii]];
  be = beSamp[[ii]]; omega = omegaSamp[[ii]]; ri = riSamp[[ii]];

  S = lyticEq["S"]; I = lyticEq["I"];
  Lterm = 1 - (S + I)/K;

  J = N[{
    {(1 - omega) * ri * Lterm - c,    beta * S},
    {be,                               -muV}
  }, wp];

  eigs = Eigenvalues[J];
  Max[Re[eigs]]
];

temperateInvadesLyticEig[ii_, lyticEq_Association] := Module[
  {r, K, c, beta, alpha, gamma, muV, f, ee, rl, S, I, Lterm, J, eigs},

  r = rSamp[[ii]]; K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; alpha = alphaSamp[[ii]];
  gamma = gammaSamp[[ii]]; muV = muVSamp[[ii]];
  f = fSamp[[ii]]; ee = eeSamp[[ii]]; rl = rlSamp[[ii]];

  S = lyticEq["S"]; I = lyticEq["I"];
  Lterm = 1 - (S + I)/K;

  J = N[{
    {rl * Lterm - c - ee,     0,             f * beta * S},
    {ee,                      -(alpha + c),   (1 - f) * beta * S},
    {0,                       alpha * gamma,  -muV}
  }, wp];

  eigs = Eigenvalues[J];
  Max[Re[eigs]]
];

lyticInvadesChronicEig[ii_, chronicEq_Association] := Module[
  {beta, alpha, gamma, muV, c, S, J, eigs},

  beta = betaSamp[[ii]]; alpha = alphaSamp[[ii]];
  gamma = gammaSamp[[ii]]; muV = muVSamp[[ii]]; c = cSamp[[ii]];

  S = chronicEq["S"];

  J = N[{
    {-(alpha + c),    beta * S},
    {gamma * alpha,   -muV}
  }, wp];

  eigs = Eigenvalues[J];
  Max[Re[eigs]]
];

lyticInvadesTemperateEig[ii_, temperateEq_Association] := Module[
  {beta, alpha, gamma, muV, c, S, J, eigs},

  beta = betaSamp[[ii]]; alpha = alphaSamp[[ii]];
  gamma = gammaSamp[[ii]]; muV = muVSamp[[ii]]; c = cSamp[[ii]];

  S = temperateEq["S"];

  J = N[{
    {-(alpha + c),    beta * S},
    {gamma * alpha,   -muV}
  }, wp];

  eigs = Eigenvalues[J];
  Max[Re[eigs]]
];

chronicInvadesTemperateEig[ii_, temperateEq_Association] := Module[
  {K, c, beta, muV, be, omega, ri, S, L, It, Lterm, J, eigs},

  K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; muV = muVSamp[[ii]];
  be = beSamp[[ii]]; omega = omegaSamp[[ii]]; ri = riSamp[[ii]];

  S = temperateEq["S"]; L = temperateEq["L"]; It = temperateEq["I"];
  Lterm = 1 - (S + L + It)/K;

  J = N[{
    {(1 - omega) * ri * Lterm - c,    beta * S},
    {be,                               -muV}
  }, wp];

  eigs = Eigenvalues[J];
  Max[Re[eigs]]
];

temperateInvadesChronicEig[ii_, chronicEq_Association] := Module[
  {K, c, beta, alpha, gamma, muV, f, ee, rl, S, Ic, Lterm, J, eigs},

  K = KSamp[[ii]]; c = cSamp[[ii]];
  beta = betaSamp[[ii]]; alpha = alphaSamp[[ii]];
  gamma = gammaSamp[[ii]]; muV = muVSamp[[ii]];
  f = fSamp[[ii]]; ee = eeSamp[[ii]]; rl = rlSamp[[ii]];

  S = chronicEq["S"]; Ic = chronicEq["I"];
  Lterm = 1 - (S + Ic)/K;

  J = N[{
    {rl * Lterm - c - ee,     0,             f * beta * S},
    {ee,                      -(alpha + c),   (1 - f) * beta * S},
    {0,                       alpha * gamma,  -muV}
  }, wp];

  eigs = Eigenvalues[J];
  Max[Re[eigs]]
];

Print["\n==================================================================="];
Print["RUNNING LHS INVASION ANALYSIS V5 (n = ", n, ")"];
Print["==================================================================="];

lambdaChrInvLyt = ConstantArray[Null, n];    
lambdaTmpInvLyt = ConstantArray[Null, n];    
lambdaLytInvChr = ConstantArray[Null, n];    
lambdaLytInvTmp = ConstantArray[Null, n];    
lambdaChrInvTmp = ConstantArray[Null, n];    
lambdaTmpInvChr = ConstantArray[Null, n];    

lyticFeasible     = ConstantArray[False, n];
chronicFeasible   = ConstantArray[False, n];
temperateFeasible = ConstantArray[False, n];

lyticStable       = ConstantArray[False, n];
chronicStable     = ConstantArray[False, n];
temperateStable   = ConstantArray[False, n];
temperateMaxReEig = ConstantArray[Null, n];

lyticEqS = ConstantArray[Null, n];    lyticEqI = ConstantArray[Null, n];    lyticEqV = ConstantArray[Null, n];
chronicEqS = ConstantArray[Null, n];  chronicEqI = ConstantArray[Null, n];  chronicEqV = ConstantArray[Null, n];
temperateEqS = ConstantArray[Null, n]; temperateEqL = ConstantArray[Null, n];
temperateEqI = ConstantArray[Null, n]; temperateEqV = ConstantArray[Null, n];

Monitor[
  Do[
    Module[{leq, ceq, teq},

      
      leq = lyticEqForDraw[ii];
      ceq = chronicEqForDraw[ii];
      teq = temperateEqForDraw[ii];

      lyticFeasible[[ii]] = leq["feasible"];
      chronicFeasible[[ii]] = ceq["feasible"];
      temperateFeasible[[ii]] = teq["feasible"];

      
      If[leq["feasible"],
        lyticEqS[[ii]] = leq["S"]; lyticEqI[[ii]] = leq["I"]; lyticEqV[[ii]] = leq["V"]];
      If[ceq["feasible"],
        chronicEqS[[ii]] = ceq["S"]; chronicEqI[[ii]] = ceq["I"]; chronicEqV[[ii]] = ceq["V"]];
      If[teq["feasible"],
        temperateEqS[[ii]] = teq["S"]; temperateEqL[[ii]] = teq["L"];
        temperateEqI[[ii]] = teq["I"]; temperateEqV[[ii]] = teq["V"]];

      
      
      
      
      
      

      
      If[leq["feasible"],
        lyticStable[[ii]] = lyticStableForDraw[ii, leq];
        lambdaChrInvLyt[[ii]] = chronicInvadesLyticEig[ii, leq];
        lambdaTmpInvLyt[[ii]] = temperateInvadesLyticEig[ii, leq];
      ];

      
      If[ceq["feasible"],
        chronicStable[[ii]] = chronicStableForDraw[ii, ceq];
        lambdaLytInvChr[[ii]] = lyticInvadesChronicEig[ii, ceq];
        lambdaTmpInvChr[[ii]] = temperateInvadesChronicEig[ii, ceq];
      ];

      
      If[teq["feasible"],
        temperateMaxReEig[[ii]] = temperateMaxReEigForDraw[ii, teq];
        temperateStable[[ii]] = temperateStableForDraw[ii, teq];
        lambdaLytInvTmp[[ii]] = lyticInvadesTemperateEig[ii, teq];
        lambdaChrInvTmp[[ii]] = chronicInvadesTemperateEig[ii, teq];
      ];
    ],
    {ii, 1, n}
  ],
  ProgressIndicator[ii, {1, n}]
];

Print["Computation complete."];

Print["\n--- FEASIBILITY & STABILITY SUMMARY ---"];
Print["Lytic equilibrium feasible: ", Count[lyticFeasible, True], " / ", n,
      " (", N[100. * Count[lyticFeasible, True] / n, 4], "%)"];
Print["Lytic equilibrium stable: ", Count[lyticStable, True], " / ", n,
      " (", N[100. * Count[lyticStable, True] / n, 4], "%)"];
Print["Chronic equilibrium feasible: ", Count[chronicFeasible, True], " / ", n,
      " (", N[100. * Count[chronicFeasible, True] / n, 4], "%)"];
Print["Chronic equilibrium stable: ", Count[chronicStable, True], " / ", n,
      " (", N[100. * Count[chronicStable, True] / n, 4], "%)"];

Module[{nFeas, nStable, nUnstable},
  nFeas = Count[temperateFeasible, True];
  nStable = Count[temperateStable, True];
  nUnstable = nFeas - nStable;
  Print["Temperate equilibrium feasible: ", nFeas, " / ", n,
        " (", N[100. * nFeas / n, 4], "%)"];
  Print["  Feasible AND stable: ", nStable,
        " (", If[nFeas > 0, N[100. * nStable / nFeas, 4], 0], "% of feasible)"];
  Print["  Feasible but UNSTABLE (oscillations): ", nUnstable,
        " (", If[nFeas > 0, N[100. * nUnstable / nFeas, 4], 0], "% of feasible)"];
  Print["  Infeasible (no positive roots with populations > ", minPop, "): ", n - nFeas];
];

Module[{lytPct, tmpUnstPct},
  lytPct = N[100. * Count[lyticStable, True] / Max[Count[lyticFeasible, True], 1], 4];
  tmpUnstPct = If[Count[temperateFeasible, True] > 0,
    N[100. * (Count[temperateFeasible, True] - Count[temperateStable, True]) /
      Count[temperateFeasible, True], 4], 0];
  Print["\n  V5 NOTE: Invasion eigenvalues are evaluated at feasible fixed points"];
  Print["  regardless of stability. For lytic (", lytPct, "% stable) and most temperate"];
  Print["  equilibria (", tmpUnstPct, "% unstable), eigenvalues at the fixed point serve as"];
  Print["  qualitative predictors of invasion success. Jensen's inequality"];
  Print["  caveat applies \[LongDash] validate numerically with RK4 simulations."];
  Print["  V5 uses log-scale LHS sampling for parameters spanning >1 order"];
  Print["  of magnitude, ensuring uniform coverage across all decades."];
];

Print["\n==================================================================="];
Print["MUTUAL INVASIBILITY CLASSIFICATION"];
Print["==================================================================="];

Print["\n  V5: All comparisons require FEASIBILITY only (not stability)."];
Print["  Log-scale sampling used for beta, muV, c, gamma, rl, f, ee."];

Print["\n--- Chronic vs Lytic ---"];
chrLytOutcome = ConstantArray["NA", n];
Module[{validIdx, nValid, coexist, chrWins, lytWins, bistable,
        fwdVals, revVals, nChrStable, nLytStable},

  
  validIdx = Select[Range[n],
    (lyticFeasible[[#]] && chronicFeasible[[#]] &&
     NumericQ[lambdaChrInvLyt[[#]]] && NumericQ[lambdaLytInvChr[[#]]]) &];
  nValid = Length[validIdx];
  nChrStable = Count[chronicStable[[validIdx]], True];
  nLytStable = Count[lyticStable[[validIdx]], True];
  Print["Valid draws (both feasible): ", nValid, " / ", n];
  Print["  Of which: chronic stable = ", nChrStable, ", lytic stable = ", nLytStable];

  If[nValid > 0,
    fwdVals = lambdaChrInvLyt[[validIdx]];
    revVals = lambdaLytInvChr[[validIdx]];

    coexist  = Count[MapThread[#1 > 0 && #2 > 0 &, {fwdVals, revVals}], True];
    chrWins  = Count[MapThread[#1 > 0 && #2 <= 0 &, {fwdVals, revVals}], True];
    lytWins  = Count[MapThread[#1 <= 0 && #2 > 0 &, {fwdVals, revVals}], True];
    bistable = Count[MapThread[#1 <= 0 && #2 <= 0 &, {fwdVals, revVals}], True];

    Do[
      Module[{idx = validIdx[[k]], fwd = fwdVals[[k]], rev = revVals[[k]]},
        chrLytOutcome[[idx]] = Which[
          fwd > 0 && rev > 0, "coexist",
          fwd > 0 && rev <= 0, "chronic_wins",
          fwd <= 0 && rev > 0, "lytic_wins",
          True, "bistable"
        ];
      ],
      {k, 1, nValid}
    ];

    Print["  Coexistence (mutual invasibility): ", coexist,
          " (", N[100. * coexist / nValid, 4], "%)"];
    Print["  Chronic wins (unidirectional):     ", chrWins,
          " (", N[100. * chrWins / nValid, 4], "%)"];
    Print["  Lytic wins (excludes chronic):     ", lytWins,
          " (", N[100. * lytWins / nValid, 4], "%)"];
    Print["  Bistable (neither invades):        ", bistable,
          " (", N[100. * bistable / nValid, 4], "%)"];
  ];
];

Print["\n--- Temperate vs Lytic ---"];
tmpLytOutcome = ConstantArray["NA", n];
Module[{validIdx, nValid, coexist, tmpWins, lytWins, bistable,
        fwdVals, revVals, nTmpStable, nLytStable},

  
  validIdx = Select[Range[n],
    (lyticFeasible[[#]] && temperateFeasible[[#]] &&
     NumericQ[lambdaTmpInvLyt[[#]]] && NumericQ[lambdaLytInvTmp[[#]]]) &];
  nValid = Length[validIdx];
  nTmpStable = Count[temperateStable[[validIdx]], True];
  nLytStable = Count[lyticStable[[validIdx]], True];
  Print["Valid draws (both feasible): ", nValid, " / ", n];
  Print["  Of which: temperate stable = ", nTmpStable, ", lytic stable = ", nLytStable];

  If[nValid > 0,
    fwdVals = lambdaTmpInvLyt[[validIdx]];
    revVals = lambdaLytInvTmp[[validIdx]];

    coexist  = Count[MapThread[#1 > 0 && #2 > 0 &, {fwdVals, revVals}], True];
    tmpWins  = Count[MapThread[#1 > 0 && #2 <= 0 &, {fwdVals, revVals}], True];
    lytWins  = Count[MapThread[#1 <= 0 && #2 > 0 &, {fwdVals, revVals}], True];
    bistable = Count[MapThread[#1 <= 0 && #2 <= 0 &, {fwdVals, revVals}], True];

    Do[
      Module[{idx = validIdx[[k]], fwd = fwdVals[[k]], rev = revVals[[k]]},
        tmpLytOutcome[[idx]] = Which[
          fwd > 0 && rev > 0, "coexist",
          fwd > 0 && rev <= 0, "temperate_wins",
          fwd <= 0 && rev > 0, "lytic_wins",
          True, "bistable"
        ];
      ],
      {k, 1, nValid}
    ];

    Print["  Coexistence (mutual invasibility): ", coexist,
          " (", N[100. * coexist / nValid, 4], "%)"];
    Print["  Temperate wins (unidirectional):   ", tmpWins,
          " (", N[100. * tmpWins / nValid, 4], "%)"];
    Print["  Lytic wins (excludes temperate):   ", lytWins,
          " (", N[100. * lytWins / nValid, 4], "%)"];
    Print["  Bistable (neither invades):        ", bistable,
          " (", N[100. * bistable / nValid, 4], "%)"];
  ];
];

Print["\n--- Chronic vs Temperate ---"];
chrTmpOutcome = ConstantArray["NA", n];
Module[{validIdx, nValid, coexist, chrWins, tmpWins, bistable,
        fwdVals, revVals, nChrStable, nTmpStable},

  
  validIdx = Select[Range[n],
    (chronicFeasible[[#]] && temperateFeasible[[#]] &&
     NumericQ[lambdaChrInvTmp[[#]]] && NumericQ[lambdaTmpInvChr[[#]]]) &];
  nValid = Length[validIdx];
  nChrStable = Count[chronicStable[[validIdx]], True];
  nTmpStable = Count[temperateStable[[validIdx]], True];
  Print["Valid draws (both feasible): ", nValid, " / ", n];
  Print["  Of which: chronic stable = ", nChrStable, ", temperate stable = ", nTmpStable];

  If[nValid > 0,
    fwdVals = lambdaChrInvTmp[[validIdx]];
    revVals = lambdaTmpInvChr[[validIdx]];

    coexist  = Count[MapThread[#1 > 0 && #2 > 0 &, {fwdVals, revVals}], True];
    chrWins  = Count[MapThread[#1 > 0 && #2 <= 0 &, {fwdVals, revVals}], True];
    tmpWins  = Count[MapThread[#1 <= 0 && #2 > 0 &, {fwdVals, revVals}], True];
    bistable = Count[MapThread[#1 <= 0 && #2 <= 0 &, {fwdVals, revVals}], True];

    Do[
      Module[{idx = validIdx[[k]], fwd = fwdVals[[k]], rev = revVals[[k]]},
        chrTmpOutcome[[idx]] = Which[
          fwd > 0 && rev > 0, "coexist",
          fwd > 0 && rev <= 0, "chronic_wins",
          fwd <= 0 && rev > 0, "temperate_wins",
          True, "bistable"
        ];
      ],
      {k, 1, nValid}
    ];

    Print["  Coexistence (mutual invasibility): ", coexist,
          " (", N[100. * coexist / nValid, 4], "%)"];
    Print["  Chronic wins (unidirectional):     ", chrWins,
          " (", N[100. * chrWins / nValid, 4], "%)"];
    Print["  Temperate wins (unidirectional):   ", tmpWins,
          " (", N[100. * tmpWins / nValid, 4], "%)"];
    Print["  Bistable (neither invades):        ", bistable,
          " (", N[100. * bistable / nValid, 4], "%)"];
  ];
];

Print["\n==================================================================="];
Print["SUMMARY TABLES (V5.2)"];
Print["==================================================================="];

safeMedian[arr_, idx_] := Module[{vals},
  If[Length[idx] == 0, Return["NA"]];
  vals = Select[arr[[idx]], NumericQ];
  If[Length[vals] == 0, "NA", N[Median[vals]]]
];

fmtNum[val_] := Module[{av, ex, mant},
  If[!NumericQ[val], Return[ToString[val]]];
  av = Abs[N[val]];
  If[av == 0, Return["0"]];
  If[0.001 <= av < 10000, Return[ToString[NumberForm[N[val], {6, 4}]]]];
  ex = Floor[Log10[av]];
  mant = N[val / 10^ex];
  StringJoin[ToString[NumberForm[mant, {4, 2}]], "e",
    If[ex >= 0, "+", ""], ToString[ex]]
];

paramArrays = {rSamp, betaSamp, muVSamp, cSamp, alphaSamp, gammaSamp,
               beSamp, omegaSamp, riSamp, rlSamp, fSamp, eeSamp, KSamp};
paramLabels = {"r", "beta", "muV", "c", "alpha", "gamma",
               "bee", "omega", "ri", "rl", "f", "ee", "K"};

lyticStableIdx     = Flatten[Position[lyticStable, True]];
lyticFeasibleIdx   = Flatten[Position[lyticFeasible, True]];
lyticUnstableIdx   = Complement[lyticFeasibleIdx, lyticStableIdx];

chronicStableIdx   = Flatten[Position[chronicStable, True]];

temperateFeasibleIdx = Flatten[Position[temperateFeasible, True]];
temperateStableIdx   = Flatten[Position[temperateStable, True]];
temperateUnstableIdx = Complement[temperateFeasibleIdx, temperateStableIdx];

Print["\n--- Table Set 1: Median steady states by stability status ---"];
Print["Strategy     | Status   | n     | S*          | I*          | V*          | L*"];
Print["-------------|----------|-------|-------------|-------------|-------------|-------------"];

printEqRow[strat_String, status_String, idx_List, sArr_, iArr_, vArr_, lArr_: None] :=
  Module[{nn, sM, iM, vM, lM},
    nn = Length[idx];
    sM = safeMedian[sArr, idx]; iM = safeMedian[iArr, idx]; vM = safeMedian[vArr, idx];
    lM = If[lArr === None, "--", safeMedian[lArr, idx]];
    Print[StringPadRight[strat, 13], "| ", StringPadRight[status, 9], "| ",
      StringPadRight[ToString[nn], 6], "| ",
      StringPadRight[If[NumericQ[sM], fmtNum[sM], sM], 12], "| ",
      StringPadRight[If[NumericQ[iM], fmtNum[iM], iM], 12], "| ",
      StringPadRight[If[NumericQ[vM], fmtNum[vM], vM], 12], "| ",
      If[NumericQ[lM], fmtNum[lM], ToString[lM]]];
  ];

printEqRow["Lytic", "Stable", lyticStableIdx, lyticEqS, lyticEqI, lyticEqV];
printEqRow["Lytic", "Unstable", lyticUnstableIdx, lyticEqS, lyticEqI, lyticEqV];
printEqRow["Chronic", "Stable", chronicStableIdx, chronicEqS, chronicEqI, chronicEqV];
printEqRow["Temperate", "Stable", temperateStableIdx, temperateEqS, temperateEqI, temperateEqV, temperateEqL];
printEqRow["Temperate", "Unstable", temperateUnstableIdx, temperateEqS, temperateEqI, temperateEqV, temperateEqL];

Print["\n--- Table Set 2: Median parameter values by stability status ---"];

printParamTable[label_String, stableIdx_List, unstableIdx_List] := Module[{},
  Print["\n  ", label, ":"];
  Print["  Parameter  | Stable (n=", Length[stableIdx], ")  | Unstable (n=", Length[unstableIdx], ")"];
  Print["  -----------|----------------------|---------------------"];
  Do[
    Module[{arr = paramArrays[[j]], name = paramLabels[[j]], sM, uM},
      sM = safeMedian[arr, stableIdx];
      uM = safeMedian[arr, unstableIdx];
      Print["  ", StringPadRight[name, 11], "| ",
        StringPadRight[If[NumericQ[sM], fmtNum[sM], sM], 21], "| ",
        If[NumericQ[uM], fmtNum[uM], uM]];
    ],
    {j, 1, Length[paramLabels]}
  ];
];

printParamTable["Lytic", lyticStableIdx, lyticUnstableIdx];
printParamTable["Temperate", temperateStableIdx, temperateUnstableIdx];
printParamTable["Chronic", chronicStableIdx, {}];  

Print["\n--- Table Set 3: Median parameter values + steady states by outcome group ---"];

printOutcomeTable[compLabel_String, outcomeArr_List, stratLabel1_String,
  eqArrays1_List, eqLabels1_List, stratLabel2_String, eqArrays2_List, eqLabels2_List] :=
  Module[{outcomes, groups, groupIdx},
    outcomes = DeleteCases[DeleteDuplicates[outcomeArr], "NA"];
    Print["\n  === ", compLabel, " ==="];
    Do[
      groupIdx = Flatten[Position[outcomeArr, outcome]];
      Print["\n  Outcome: ", outcome, " (n=", Length[groupIdx], ")"];
      Print["  --- Parameters ---"];
      Do[
        Module[{arr = paramArrays[[j]], name = paramLabels[[j]], mM},
          mM = safeMedian[arr, groupIdx];
          Print["    ", StringPadRight[name, 10], ": ",
            If[NumericQ[mM], fmtNum[mM], mM]];
        ], {j, 1, Length[paramLabels]}];
      Print["  --- ", stratLabel1, " equilibrium ---"];
      Do[
        Module[{arr = eqArrays1[[j]], label = eqLabels1[[j]], mM},
          mM = safeMedian[arr, groupIdx];
          Print["    ", StringPadRight[label, 10], ": ",
            If[NumericQ[mM], fmtNum[mM], mM]];
        ], {j, 1, Length[eqLabels1]}];
      Print["  --- ", stratLabel2, " equilibrium ---"];
      Do[
        Module[{arr = eqArrays2[[j]], label = eqLabels2[[j]], mM},
          mM = safeMedian[arr, groupIdx];
          Print["    ", StringPadRight[label, 10], ": ",
            If[NumericQ[mM], fmtNum[mM], mM]];
        ], {j, 1, Length[eqLabels2]}];
    , {outcome, outcomes}];
  ];

printOutcomeTable["Chronic vs Lytic (CL)", chrLytOutcome,
  "Lytic", {lyticEqS, lyticEqI, lyticEqV}, {"S*", "I*", "V*"},
  "Chronic", {chronicEqS, chronicEqI, chronicEqV}, {"S*", "I*", "V*"}];

printOutcomeTable["Temperate vs Lytic (TL)", tmpLytOutcome,
  "Lytic", {lyticEqS, lyticEqI, lyticEqV}, {"S*", "I*", "V*"},
  "Temperate", {temperateEqS, temperateEqL, temperateEqI, temperateEqV}, {"S*", "L*", "I*", "V*"}];

printOutcomeTable["Chronic vs Temperate (CT)", chrTmpOutcome,
  "Chronic", {chronicEqS, chronicEqI, chronicEqV}, {"S*", "I*", "V*"},
  "Temperate", {temperateEqS, temperateEqL, temperateEqI, temperateEqV}, {"S*", "L*", "I*", "V*"}];

Print["\n--- Exporting summary tables CSV ---"];

Module[{summaryRows = {}, appendRow},
  appendRow[tableType_, stratOrComp_, group_, variable_, value_] :=
    AppendTo[summaryRows, {tableType, stratOrComp, group, variable,
      If[NumericQ[value], N[value], ToString[value]]}];

  
  Do[
    Module[{strat, status, idx, sA, iA, vA, lA},
      {strat, status, idx, sA, iA, vA, lA} = entry;
      appendRow["stability_eq", strat, status, "S_star", safeMedian[sA, idx]];
      appendRow["stability_eq", strat, status, "I_star", safeMedian[iA, idx]];
      appendRow["stability_eq", strat, status, "V_star", safeMedian[vA, idx]];
      If[lA =!= None,
        appendRow["stability_eq", strat, status, "L_star", safeMedian[lA, idx]]];
    ],
    {entry, {
      {"lytic", "stable", lyticStableIdx, lyticEqS, lyticEqI, lyticEqV, None},
      {"lytic", "unstable", lyticUnstableIdx, lyticEqS, lyticEqI, lyticEqV, None},
      {"chronic", "stable", chronicStableIdx, chronicEqS, chronicEqI, chronicEqV, None},
      {"temperate", "stable", temperateStableIdx, temperateEqS, temperateEqI, temperateEqV, temperateEqL},
      {"temperate", "unstable", temperateUnstableIdx, temperateEqS, temperateEqI, temperateEqV, temperateEqL}
    }}
  ];

  
  Do[
    Module[{strat, stIdx, uIdx},
      {strat, stIdx, uIdx} = grp;
      Do[
        appendRow["stability_param", strat, "stable", paramLabels[[j]], safeMedian[paramArrays[[j]], stIdx]];
        If[Length[uIdx] > 0,
          appendRow["stability_param", strat, "unstable", paramLabels[[j]], safeMedian[paramArrays[[j]], uIdx]]
        ];
      , {j, 1, Length[paramLabels]}];
    ],
    {grp, {
      {"lytic", lyticStableIdx, lyticUnstableIdx},
      {"temperate", temperateStableIdx, temperateUnstableIdx},
      {"chronic", chronicStableIdx, {}}  
    }}
  ];

  
  Do[
    Module[{comp, outcomeArr, eqArrs, eqLbls, outcomes, groupIdx},
      {comp, outcomeArr, eqArrs, eqLbls} = spec;
      outcomes = DeleteCases[DeleteDuplicates[outcomeArr], "NA"];
      Do[
        groupIdx = Flatten[Position[outcomeArr, outcome]];
        
        Do[
          appendRow["outcome_param", comp, outcome, paramLabels[[j]], safeMedian[paramArrays[[j]], groupIdx]];
        , {j, 1, Length[paramLabels]}];
        
        Do[
          appendRow["outcome_eq", comp, outcome, eqLbls[[j]], safeMedian[eqArrs[[j]], groupIdx]];
        , {j, 1, Length[eqLbls]}];
      , {outcome, outcomes}];
    ],
    {spec, {
      {"CL", chrLytOutcome,
        {lyticEqS, lyticEqI, lyticEqV, chronicEqS, chronicEqI, chronicEqV},
        {"lytic_S", "lytic_I", "lytic_V", "chronic_S", "chronic_I", "chronic_V"}},
      {"TL", tmpLytOutcome,
        {lyticEqS, lyticEqI, lyticEqV, temperateEqS, temperateEqL, temperateEqI, temperateEqV},
        {"lytic_S", "lytic_I", "lytic_V", "temp_S", "temp_L", "temp_I", "temp_V"}},
      {"CT", chrTmpOutcome,
        {chronicEqS, chronicEqI, chronicEqV, temperateEqS, temperateEqL, temperateEqI, temperateEqV},
        {"chronic_S", "chronic_I", "chronic_V", "temp_S", "temp_L", "temp_I", "temp_V"}}
    }}
  ];

  
  Module[{summaryPath, summaryHeader},
    summaryPath = DirectoryName[$InputFileName] <> "Mutual_invasibility_LHS_V5.2_summary_tables.csv";
    summaryHeader = {"table_type", "strategy_or_comparison", "group", "variable", "median"};
    Export[summaryPath, Prepend[summaryRows, summaryHeader], "CSV"];
    Print["Summary tables exported to: ", summaryPath];
    Print["  ", Length[summaryRows], " rows"];
  ];
];

Print["\n==================================================================="];
Print["SENSITIVITY ANALYSIS"];
Print["==================================================================="];

ClearAll[runSensitivity];
runSensitivity[eigVals_List, validIdx_List, paramArrays_List,
               paramNames_List, titlePrefix_String] := Module[
  {y, X, XR, yR, spearman, prcc, sortedIdx},

  If[Length[validIdx] < 100,
    Print["\n--- ", titlePrefix, " ---"];
    Print["Too few valid draws (", Length[validIdx], ") for sensitivity analysis"];
    Return[];
  ];

  y = N[eigVals[[validIdx]]];
  X = N[Transpose[paramArrays[[All, validIdx]]]];

  Print["\n--- ", titlePrefix, " (n=", Length[validIdx], ") ---"];

  
  XR = Transpose[rankVec /@ Transpose[X]];
  yR = rankVec[y];
  spearman = Table[safeCorr[XR[[All, j]], yR], {j, 1, Length[paramNames]}];

  Print["\nSpearman rank correlations:"];
  Do[Print["  ", paramNames[[j]], ": ", N[spearman[[j]], 4]],
    {j, 1, Length[paramNames]}];

  
  prcc = Table[prccOne[XR, yR, j], {j, 1, Length[paramNames]}];

  Print["\nPRCC:"];
  Do[Print["  ", paramNames[[j]], ": ", N[prcc[[j]], 4]],
    {j, 1, Length[paramNames]}];

  
  
  
  
  Print["\nCloud plots:"];
  numericIdx = Flatten[Position[prcc, _?NumericQ]];
  sortedIdx = numericIdx[[Reverse[Ordering[Abs[N /@ prcc[[numericIdx]]]]]]];
  Do[
    Module[{j = sortedIdx[[k]]},
      Print[cloudPlot[X[[All, j]], y, titlePrefix <> " — " <> paramNames[[j]], 0]];
    ],
    {k, 1, Min[4, Length[sortedIdx]]}
  ];

  
  Print["\nHistogram of \[Lambda]_max:"];
  Print[Histogram[y, 40, "Probability",
    Frame -> True,
    FrameLabel -> {"\[Lambda]_max", "Probability"},
    Epilog -> {Red, Dashed, Thick, Line[{{0, 0}, {0, 1}}]},
    ImageSize -> Large,
    PlotLabel -> titlePrefix
  ]];

  Print["Fraction with \[Lambda]_max > 0: ",
    N[Count[Select[y, # > 0 &], _] / Length[y], 4]];
];

sharedParams = {rSamp, betaSamp, muVSamp, cSamp, alphaSamp, gammaSamp};
sharedNames = {"r", "\[Beta]", "\[Mu]V", "c", "\[Alpha]", "\[Gamma]"};

Module[{validIdx, pArrays, pNames},
  validIdx = Select[Range[n], (lyticFeasible[[#]] && NumericQ[lambdaChrInvLyt[[#]]]) &];
  pArrays = Join[sharedParams, {beSamp, omegaSamp, riSamp}];
  pNames = Join[sharedNames, {"bee", "\[Omega]", "ri"}];
  runSensitivity[lambdaChrInvLyt, validIdx, pArrays, pNames, "Chronic \[Rule] Lytic"];
];

Module[{validIdx, pArrays, pNames},
  validIdx = Select[Range[n], (lyticFeasible[[#]] && NumericQ[lambdaTmpInvLyt[[#]]]) &];
  pArrays = Join[sharedParams, {fSamp, eeSamp, rlSamp}];
  pNames = Join[sharedNames, {"f", "ee", "rL"}];
  runSensitivity[lambdaTmpInvLyt, validIdx, pArrays, pNames, "Temperate \[Rule] Lytic"];
];

Module[{validIdx, pArrays, pNames},
  validIdx = Select[Range[n], (chronicFeasible[[#]] && NumericQ[lambdaLytInvChr[[#]]]) &];
  pArrays = Join[sharedParams, {beSamp, omegaSamp, riSamp}];
  pNames = Join[sharedNames, {"bee", "\[Omega]", "ri"}];
  runSensitivity[lambdaLytInvChr, validIdx, pArrays, pNames, "Lytic \[Rule] Chronic"];
];

Module[{validIdx, pArrays, pNames},
  validIdx = Select[Range[n], (temperateFeasible[[#]] && NumericQ[lambdaLytInvTmp[[#]]]) &];
  pArrays = Join[sharedParams, {fSamp, eeSamp, rlSamp}];
  pNames = Join[sharedNames, {"f", "ee", "rL"}];
  runSensitivity[lambdaLytInvTmp, validIdx, pArrays, pNames, "Lytic \[Rule] Temperate"];
];

Module[{validIdx, pArrays, pNames},
  validIdx = Select[Range[n], (temperateFeasible[[#]] && NumericQ[lambdaChrInvTmp[[#]]]) &];
  pArrays = Join[{rSamp, betaSamp, muVSamp, cSamp, alphaSamp, gammaSamp},
                 {beSamp, omegaSamp, riSamp, fSamp, eeSamp, rlSamp}];
  pNames = Join[sharedNames, {"bee", "\[Omega]", "ri", "f", "ee", "rL"}];
  runSensitivity[lambdaChrInvTmp, validIdx, pArrays, pNames, "Chronic \[Rule] Temperate"];
];

Module[{validIdx, pArrays, pNames},
  validIdx = Select[Range[n], (chronicFeasible[[#]] && NumericQ[lambdaTmpInvChr[[#]]]) &];
  pArrays = Join[{rSamp, betaSamp, muVSamp, cSamp, alphaSamp, gammaSamp},
                 {beSamp, omegaSamp, riSamp, fSamp, eeSamp, rlSamp}];
  pNames = Join[sharedNames, {"bee", "\[Omega]", "ri", "f", "ee", "rL"}];
  runSensitivity[lambdaTmpInvChr, validIdx, pArrays, pNames, "Temperate \[Rule] Chronic"];
];

Print["\n==================================================================="];
Print["2D PARAMETER SPACE MAPS"];
Print["==================================================================="];

ClearAll[paramSpaceScatter];
paramSpaceScatter[xVals_List, yVals_List, outcomes_List,
                  xLab_String, yLab_String, title_String] := Module[
  {validMask, xV, yV, oV, pts, coexPts, win1Pts, win2Pts, biPts,
   cats, colors, legends},

  validMask = MapThread[#1 =!= "NA" && NumericQ[#2] && NumericQ[#3] &, {outcomes, xVals, yVals}];
  xV = Pick[xVals, validMask]; yV = Pick[yVals, validMask]; oV = Pick[outcomes, validMask];

  If[Length[xV] < 10, Print["Too few points for ", title]; Return[Null]];

  cats = Union[oV];
  colors = <|"coexist" -> Blue, "chronic_wins" -> Red, "lytic_wins" -> Darker[Green],
              "temperate_wins" -> Orange, "bistable" -> Gray|>;
  legends = <|"coexist" -> "Coexist", "chronic_wins" -> "Chronic wins",
              "lytic_wins" -> "Lytic wins", "temperate_wins" -> "Temperate wins",
              "bistable" -> "Bistable"|>;

  ListPlot[
    Table[
      Pick[Transpose[{xV, yV}], oV, cat],
      {cat, cats}
    ],
    PlotStyle -> Table[{PointSize[Small], Opacity[0.6],
                        Lookup[colors, cat, Black]}, {cat, cats}],
    PlotLegends -> (Lookup[legends, #, #] & /@ cats),
    AxesLabel -> {xLab, yLab},
    PlotLabel -> title,
    ImageSize -> Large,
    ScalingFunctions -> {"Log10", "Log10"}
  ]
];

Print["\n--- Chronic vs Lytic parameter space ---"];
Print[paramSpaceScatter[N[betaSamp], N[beSamp], chrLytOutcome,
  "\[Beta]", "bee", "Chronic vs Lytic: \[Beta] vs bee"]];
Print[paramSpaceScatter[N[betaSamp], N[omegaSamp], chrLytOutcome,
  "\[Beta]", "\[Omega]", "Chronic vs Lytic: \[Beta] vs \[Omega]"]];
Print[paramSpaceScatter[N[cSamp], N[beSamp], chrLytOutcome,
  "c", "bee", "Chronic vs Lytic: c vs bee"]];
Print[paramSpaceScatter[N[betaSamp], N[riSamp], chrLytOutcome,
  "\[Beta]", "ri", "Chronic vs Lytic: \[Beta] vs ri"]];

Print["\n--- Temperate vs Lytic parameter space ---"];
Print[paramSpaceScatter[N[betaSamp], N[fSamp], tmpLytOutcome,
  "\[Beta]", "f", "Temperate vs Lytic: \[Beta] vs f"]];
Print[paramSpaceScatter[N[betaSamp], N[eeSamp], tmpLytOutcome,
  "\[Beta]", "ee", "Temperate vs Lytic: \[Beta] vs ee"]];
Print[paramSpaceScatter[N[fSamp], N[eeSamp], tmpLytOutcome,
  "f", "ee", "Temperate vs Lytic: f vs ee"]];

Print["\n--- Chronic vs Temperate parameter space ---"];
Print[paramSpaceScatter[N[beSamp], N[fSamp], chrTmpOutcome,
  "bee", "f", "Chronic vs Temperate: bee vs f"]];
Print[paramSpaceScatter[N[omegaSamp], N[eeSamp], chrTmpOutcome,
  "\[Omega]", "ee", "Chronic vs Temperate: \[Omega] vs ee"]];
Print[paramSpaceScatter[N[betaSamp], N[beSamp], chrTmpOutcome,
  "\[Beta]", "bee", "Chronic vs Temperate: \[Beta] vs bee"]];

Print["\n==================================================================="];
Print["EXPORTING RESULTS"];
Print["==================================================================="];

Module[{header, data, outPath},
  header = {"draw", "r", "K", "beta", "muV", "c",
            "alpha", "gamma", "bee", "omega", "ri",
            "f", "ee", "rl",
            "lyticFeasible", "chronicFeasible", "temperateFeasible",
            "lyticStable", "chronicStable", "temperateStable",
            "lambda_chr_inv_lyt", "lambda_tmp_inv_lyt",
            "lambda_lyt_inv_chr", "lambda_lyt_inv_tmp",
            "lambda_chr_inv_tmp", "lambda_tmp_inv_chr",
            "CL_outcome", "TL_outcome", "CT_outcome"};

  data = Table[{
    ii,
    N[rSamp[[ii]]], N[KSamp[[ii]]], N[betaSamp[[ii]]], N[muVSamp[[ii]]], N[cSamp[[ii]]],
    N[alphaSamp[[ii]]], N[gammaSamp[[ii]]], N[beSamp[[ii]]], N[omegaSamp[[ii]]], N[riSamp[[ii]]],
    N[fSamp[[ii]]], N[eeSamp[[ii]]], N[rlSamp[[ii]]],
    If[lyticFeasible[[ii]], 1, 0],
    If[chronicFeasible[[ii]], 1, 0],
    If[temperateFeasible[[ii]], 1, 0],
    If[lyticStable[[ii]], 1, 0],
    If[chronicStable[[ii]], 1, 0],
    If[temperateStable[[ii]], 1, 0],
    If[NumericQ[lambdaChrInvLyt[[ii]]], N[lambdaChrInvLyt[[ii]]], "NA"],
    If[NumericQ[lambdaTmpInvLyt[[ii]]], N[lambdaTmpInvLyt[[ii]]], "NA"],
    If[NumericQ[lambdaLytInvChr[[ii]]], N[lambdaLytInvChr[[ii]]], "NA"],
    If[NumericQ[lambdaLytInvTmp[[ii]]], N[lambdaLytInvTmp[[ii]]], "NA"],
    If[NumericQ[lambdaChrInvTmp[[ii]]], N[lambdaChrInvTmp[[ii]]], "NA"],
    If[NumericQ[lambdaTmpInvChr[[ii]]], N[lambdaTmpInvChr[[ii]]], "NA"],
    chrLytOutcome[[ii]],
    tmpLytOutcome[[ii]],
    chrTmpOutcome[[ii]]
  }, {ii, 1, n}];

  outPath = DirectoryName[$InputFileName] <> "Mutual_invasibility_LHS_V5.2_results.csv";

  Export[outPath,
    Prepend[data, header],
    "CSV"
  ];
  Print["Results exported to: ", outPath];
];

Print["\n=== LHS MUTUAL INVASIBILITY ANALYSIS V5.2 COMPLETE ==="];
