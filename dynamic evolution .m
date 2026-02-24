(* ::Package:: *)

(* ::Section:: *)
(*The dynamical evolution of the EMS system in flat spacetime*)


(* ::Subsubsection:: *)
(*the initial RN solution    parameter setup  finite\:2011difference matrix   the KO dissipation matrix *)


(* ::Input:: *)
(*Off[Solve::ratnz,General::munfl]; *)
(*{M=1.,Q=9/10,b0=-8,NN=2001};(*parameters  Grid size: 2^11+1*) *)
(*amp=0.0001;*)
(*\[Zeta]RNBH[r_]=Sqrt[(2M)/r-Q^2/r^2]/Sqrt[1/r];RNBH[r_]=1-(Sqrt[1/r]\[Zeta]RNBH[r])^2;(*the RN black hole,for comparison at the initial time*) *)
(*sol=Sort[r/.Solve[RNBH[r]==0,r],#1>#2&];*)
(*{rh,rin}={M,1.3M};(*Inner boundary:slightly inside the horizon*) *)
(*{zh,zin}=r/(M+r)/.{r->{rh,rin}};(*change of coordinates*) *)
(*zout=1;(*Outer boundary:at infinity*) *)
(*dz=(zout-zin)/(NN-1); *)
(*z=Range[zin,zout,dz];(*a spatial grid with uniformly spaced points*)*)
(*zc=z;zc[[-1]]=zout(1-10^-15);(*avoid 1/0*)*)
(*{dz0,dz1}=NDSolve`FiniteDifferenceDerivative[Derivative[#],z,"DifferenceOrder"->4]["DifferentiationMatrix"]& /@ Range[0,1]; (*the fourth\:2011order accurate derivative matrix for the grid*)*)
(*kodissMatrix=NDSolve`FiniteDifferenceDerivative[Derivative[6],Range[Length[z]],"DifferenceOrder"->1]["DifferentiationMatrix"];(*KO dissipation*)*)
(*z2=Drop[z,-1];(*Separate treatment for outer boundary points*)*)
(*dt=6/10 Min[Differences[z2/(1-z2)]];*)
(*eps=1/128 dt/dz; (*the prefactor of the KO dissipation term*)*)
(**)
(*{{zin,zh,zout},{(z2/(1-z2))[[-2]],(z2/(1-z2))[[-1]]},{NN,dt}}*)
(**)
(*metricplot1=ListLinePlot[Thread[{r,RNBH[r]}/.r->z2/(1-z2)],PlotRange->All];*)
(*metricplot2=Plot[RNBH[r],{r,rin,z[[NN-1]]/(1-z[[NN-1]])},PlotStyle->Orange]; *)
(*Show[metricplot2,metricplot1]*)


(* ::Subsubsection:: *)
(*The initial conditions of the fields*)


(* ::Input:: *)
(*ic[r_,ru_,rl_]:=amp \!\(\**)
(*TagBox[GridBox[{*)
(*{"\[Piecewise]", GridBox[{*)
(*{"0", *)
(*RowBox[{"r", "<=", "rl"}]},*)
(*{*)
(*RowBox[{" ", *)
(*RowBox[{*)
(*SuperscriptBox["E", *)
(*RowBox[{*)
(*RowBox[{"-", *)
(*FractionBox["1", *)
(*RowBox[{"ru", "-", "r"}]]}], "-", *)
(*FractionBox["1", *)
(*RowBox[{*)
(*RowBox[{"-", "rl"}], "+", "r"}]]}]], " ", *)
(*SuperscriptBox[*)
(*RowBox[{"(", *)
(*RowBox[{"ru", "-", "r"}], ")"}], "2"], " ", *)
(*SuperscriptBox[*)
(*RowBox[{"(", *)
(*RowBox[{*)
(*RowBox[{"-", "rl"}], "+", "r"}], ")"}], "2"]}]}], *)
(*RowBox[{"rl", "<", "r", "<", "ru"}]},*)
(*{"0", *)
(*RowBox[{"r", ">=", "ru"}]}*)
(*},*)
(*AllowedDimensions->{2, Automatic},*)
(*Editable->True,*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},*)
(*Selectable->True]}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*"Piecewise",*)
(*DeleteWithContents->True,*)
(*Editable->False,*)
(*SelectWithContents->True,*)
(*Selectable->False,*)
(*StripWrapperBoxes->True]\)*)
(*\[Phi]0in=ic[#,9,4]&/@(z2/(1-z2));*)
(*\[Phi]0=Join[\[Phi]0in,{0}];(*Initial conditions*)*)
(*P0=(-1+z)^2 dz1 . \[Phi]0;*)
(*\[CapitalPhi]0=(-1+z)^2 dz1 . \[Phi]0;*)
(*f[\[Phi]_]:=1/(1+b0 \[Phi]^2);(*the coupling function*) *)
(*ListLinePlot[\[Phi]0,AxesOrigin->{rh,0},PlotRange->All]*)
(*ListLinePlot[Transpose[{z2/(1-z2),\[Phi]0in}],PlotRange->{{0,25},All}]*)


(* ::Subsection:: *)
(*The constraint equations and their solution*)


(* ::Subsubsection:: *)
(*Outer\:2011boundary \|01d701 constraint equation (nonlinear; Newton\[Dash]Raphson; used only initially)*)


(* ::Input:: *)
(*(*Constraint equation for the metric \[Alpha]*)*)


(* ::Input:: *)
(**)
(*\[Zeta]dzIn=(-P (1/(1-x))^(7/2) x^(3/2) \[CapitalPhi]-(Q^2/(x^2 ff[\[Phi]])+(x^2 (P^2+\[CapitalPhi]^2))/(-1+x)^4)/(2 \[Zeta][x])+Derivative[1][\[Zeta]][x]);(* \[Zeta]Constraint equations (symbolic form),evaluated only at interior grid points*)*)
(**)
(*d\[Zeta]dzIn=D[\[Zeta]dzIn/.\[Zeta]->Function[x,\[Zeta][x]+\[Epsilon] d\[Zeta][x]],\[Epsilon]]/.\[Epsilon]->0;(*\[Zeta] Variational linear function;Jacobian matrix*)*)
(**)
(*dz0in=Most[dz0];dz1in=Most[dz1];(*Derivative matrices evaluated at interior grid points*)*)
(*(*PQx4/(-1+x)^4~~(x^2 (P^2+\[CapitalPhi]^2))/(-1+x)^4*)*)
(*\[Zeta]dzOut=(-((Q^2/ff[\[Phi]]+PQx4)/(2 \[Zeta][x]))+Derivative[1][\[Zeta]][x]); (*\[Zeta]Variational linear function; Jacobian matrix; using L\[CloseCurlyQuote]H\[OHat]pital\[CloseCurlyQuote]s rule*)*)
(*d\[Zeta]dzOut=D[\[Zeta]dzOut/.\[Zeta]->Function[x,\[Zeta][x]+\[Epsilon] d\[Zeta][x]],\[Epsilon]]/.\[Epsilon]->0;*)
(*dz41=(dz1 . dz1 . dz1 . dz1)[[-1]];*)
(*dz01=dz0[[-1]];*)
(*dz11=dz1[[-1]];*)
(*dz0first=dz0[[1]];*)
(*solve\[Zeta][\[Zeta]0t_,P0t_,\[CapitalPhi]0t_,\[Phi]0t_]:=((*Newton\[Dash]Raphson iteration for solving the constraint*)*)
(*\[Zeta]0temp=\[Zeta]0t;\[Delta]\[Zeta]=\[Zeta]0temp;*)
(*{P0tin,\[CapitalPhi]0tin,\[Phi]0tin,fin}={Most[P0t],Most[\[CapitalPhi]0t],Most[\[Phi]0t],Most[f[\[Phi]0t]]};(*Extract values at interior grid points*)*)
(*PQx4v=dz41 . (z^2 (P0t^2+\[CapitalPhi]0t^2))/24;(*L\[CloseCurlyQuote]H\[OHat]pital\[CloseCurlyQuote]s rule*)*)
(*ff1=f[\[Phi]0t[[-1]]];*)
(*While[Norm[\[Delta]\[Zeta]]>10^-12,*)
(*pairsIn={d\[Zeta]dzIn,-\[Zeta]dzIn}/.{P->P0tin,\[CapitalPhi]->\[CapitalPhi]0tin,\[Phi]->\[Phi]0tin,ff[\[Phi]]->fin,x->z2,d\[Zeta][x]->dz0in,Derivative[1][d\[Zeta]][x]->dz1in,\[Zeta][x]->Most[\[Zeta]0temp], \[Zeta]'[x]->dz1in . \[Zeta]0temp};(*Jacobian matrix and equation values at interior grid points*)*)
(*pairsOut={d\[Zeta]dzOut,-\[Zeta]dzOut}/.{PQx4->PQx4v,ff[\[Phi]]->ff1,\[Zeta][x]->\[Zeta]0temp[[-1]], d\[Zeta][x]->dz01,Derivative[1][d\[Zeta]][x]->dz11,\[Zeta]'[x]->dz11 . \[Zeta]0temp};(*Jacobian matrix and equation values at the asymptotic boundary*)*)
(*pairsIn[[1,1]]=dz0first;*)
(*pairsIn[[2,1]]=0;*)
(*pairs={Join[pairsIn[[1]],{pairsOut[[1]]}],Join[pairsIn[[2]],{pairsOut[[2]]}]};(*Append the outer\:2011boundary contribution*)*)
(*\[Delta]\[Zeta]=LinearSolve@@pairs;(*Iteratively solve for the update at each step*)*)
(*\[Zeta]0temp=\[Zeta]0temp+\[Delta]\[Zeta];*)
(*];*)
(*\[Zeta]0temp (*return the solution obtained*)*)
(*) *)


(* ::Subsubsection:: *)
(*The \[Alpha]\:2011constraint equation and its solution (linear)*)


(* ::Input:: *)
(*(*Constraint equation for the metric \[Alpha]*)*)
(*inv1z=Join[1/(1-z2),{0}];(*Boundary treated separately*) *)
(*onelast=SparseArray[NN->1,NN];(*Boundary handling*)*)
(*lnadz=dz1;*)
(*lnadz[[-1]]=onelast;(*the boundary condition a[NN]=1*)*)
(*lnadzInv=Inverse[lnadz];*)
(*C0=inv1z^(7/2) z^(3/2);*)
(*solvea[\[Zeta]_,P_,\[CapitalPhi]_]:=((*Linear equation used to solve for the metric \[Alpha]*)*)
(*Exp[-lnadzInv . ((P \[CapitalPhi] )/\[Zeta] C0)]*)
(*); *)


(* ::Subsubsection:: *)
(*Solving \[Alpha] and \[CurlyPhi] in the spatial direction (required at every step)*)


(* ::Input:: *)
(*C2=(-1+z)^2;*)
(*solveSpace[\[Zeta]_,P_,\[Phi]_]:=((*\[Zeta] fixed\[RightArrow]solve \[Alpha] and \[CapitalPhi] for later times*)*)
(*\[CapitalPhi]temp0=C2 dz1 . \[Phi]; *)
(*{\[CapitalPhi]temp0,solvea[\[Zeta],P,\[CapitalPhi]temp0]}*)
(*); *)


(* ::Subsubsection:: *)
(*Metric initialization and initialization of all variables*)


(* ::Input:: *)
(*solveSpace0[\[Zeta]_,P_,\[Phi]_]:=((*Treating \[Zeta] as a trial guess,solve for \[Zeta] to obtain the precise initial metric*)*)
(*\[CapitalPhi]temp0=C2 dz1 . \[Phi];*)
(*\[Zeta]temp0=solve\[Zeta][\[Zeta],P,\[CapitalPhi]temp0,\[Phi]];*)
(*atemp0=solvea[\[Zeta]temp0,P,\[CapitalPhi]temp0];*)
(*{\[CapitalPhi]temp0,\[Zeta]temp0,atemp0}*)
(*);*)
(**)
(*\[Zeta]0=Join[\[Zeta]RNBH[z2/(1-z2)],{Sqrt[2M]}];(*the initial trial value of \[Zeta] *) *)
(*{\[CapitalPhi]t,\[Zeta]t,at}=solveSpace0[\[Zeta]0,P0,\[Phi]0];(*Metric \[Alpha] and \[Zeta] at the initial time*)*)
(*{Pt,\[Phi]t}={P0,\[Phi]0};(*Store the initial matter fields P,\[Phi],\[CapitalPhi] and the metric \[Zeta],\[Alpha]*) *)
(*(\[Zeta]t[[-1]]-\[Zeta]0[[-1]])/\[Zeta]0[[-1]]*)


(* ::Input:: *)
(*ListLinePlot[{Most[at^2 (1-(Sqrt[(1-z)/z]\[Zeta]t)^2)],RNBH[z2/(1-z2)]},PlotRange->All]*)
(*\[Zeta]dzConstraint[\[Zeta]All_,\[Zeta]_,P_,\[CapitalPhi]_,\[Phi]_]:=((dz1in . \[Zeta]All)-Q^2/(2 z2^2 f[\[Phi]]\[Zeta])-z2^(3/2)/(1-z2)^(7/2) P \[CapitalPhi]-(z2^2 (P^2+\[CapitalPhi]^2))/(2 (-1+z2)^4 \[Zeta]));(*the constraint equation for \[Zeta]*)*)
(*ListLinePlot[\[Zeta]dzConstraint[\[Zeta]t,Most[\[Zeta]t],Most[Pt],Most[\[CapitalPhi]t],Most[\[Phi]t]]//Rest//Abs//Log,PlotRange->All](*the*)
(*\|01d701-constraint imposed at the initial time*)*)


(* ::Subsection:: *)
(*The evolution equations*)


(* ::Input:: *)
(*C12=Sqrt[(1-z)/z];*)
(*C1=z^2 inv1z^2;*)
(*C3=(2(1-z))/z;*)
(*C4=((Q^2) ((-1+z)^4) )/(2 (z^4) );*)
(*evolveTime[\[Alpha]_,\[Zeta]_,P_,\[CapitalPhi]_,\[Phi]_]:=((*perform the evolution calculation for d\[Zeta]t,dPt and d\[Phi]t*)*)
(*{AA,BB}={\[Alpha](P \[Zeta] C12 +\[CapitalPhi]),\[Alpha](P+C12 \[Zeta] \[CapitalPhi])};(*intermediate vars*)*)
(*{d\[Zeta]t,dPt,d\[Phi]t}={C1 (AA BB)/(\[Zeta] \[Alpha]),C2 dz1 . AA+ C3 AA+C4 (\[Alpha] Derivative[1][f][\[Phi]])/f[\[Phi]]^2,BB}*)
(*); *)


(* ::Subsection:: *)
(*the RK4 time\:2011evolution loop*)


(* ::Input:: *)
(*TStep=50;(*Save data every 100 steps*)*)
(*Offset[General::munfl]; *)
(*t0=1;tn=t0;*)


(* ::Input:: *)
(*T0=Floor[500/dt];(*Tt[0]=0;*)*)
(*Monitor[Timing[*)
(*While[tn<T0+1,(*the termination condition of the iteration*)*)
(*diss=eps*kodissMatrix . #&/@{\[Zeta]t,Pt,\[Phi]t};(*KO dissipation*)*)
(**)
(*(*RK4 step 1: Using values of a,\[Zeta],P,\[CapitalPhi],\[Phi] at t, compute RK1 derivatives of \[Zeta],P,\[Phi] at t*)*)
(*RK1=evolveTime[at,\[Zeta]t,Pt,\[CapitalPhi]t,\[Phi]t];*)
(**)
(*(*RK4 step 2: use RK1 to get \[Zeta],P,\[Phi] at t+dt/2; update \[CapitalPhi],a; compute R2 derivatives of \[Zeta],P,\[Phi] at t+dt/2*) *)
(*{\[Zeta]temp,Ptemp,\[Phi]temp}={\[Zeta]t,Pt,\[Phi]t}+RK1*dt/2+1/2 diss;*)
(*{\[CapitalPhi]temp,atemp}=solveSpace[\[Zeta]temp,Ptemp,\[Phi]temp]; *)
(*RK2=evolveTime[atemp,\[Zeta]temp,Ptemp,\[CapitalPhi]temp,\[Phi]temp];*)
(**)
(*(*RK4 step 3: use RK2 to get \[Zeta],P,\[Phi] at t+dt/2; update \[CapitalPhi],a; compute R3 derivatives of \[Zeta],P,\[Phi] at t+dt/2*)*)
(*{\[Zeta]temp,Ptemp,\[Phi]temp}={\[Zeta]t,Pt,\[Phi]t}+RK2*dt/2+1/2 diss;*)
(*{\[CapitalPhi]temp,atemp}=solveSpace[\[Zeta]temp,Ptemp,\[Phi]temp]; *)
(*RK3=evolveTime[atemp,\[Zeta]temp,Ptemp,\[CapitalPhi]temp,\[Phi]temp];*)
(**)
(*(*RK4 step 4: use RK3 to get \[Zeta],P,\[Phi] at t+dt; update \[CapitalPhi],a; compute RK4 derivatives of \[Zeta],P,\[Phi] at t+dt*)*)
(*{\[Zeta]temp,Ptemp,\[Phi]temp}={\[Zeta]t,Pt,\[Phi]t}+RK3*dt+diss;*)
(*{\[CapitalPhi]temp,atemp}=solveSpace[\[Zeta]temp,Ptemp,\[Phi]temp]; *)
(*RK4=evolveTime[atemp,\[Zeta]temp,Ptemp,\[CapitalPhi]temp,\[Phi]temp];*)
(**)
(*(*RK4 step 5:use RK1234 to compute \[Zeta],P,\[Phi] derivatives at t+dt,then update \[CapitalPhi] and \|01d44e*)*)
(*{\[Zeta]temp,Ptemp,\[Phi]temp}={\[Zeta]t,Pt,\[Phi]t}+(RK1+2RK2+2RK3+RK4)*dt/6+diss;*)
(*{\[CapitalPhi]temp,atemp}=solveSpace[\[Zeta]temp,Ptemp,\[Phi]temp]; *)
(*{Ptemp[[-1]]=0;\[CapitalPhi]temp[[-1]]=0;\[Phi]temp[[-1]]=0};(*Boundary conditions for matter fields*)*)
(*{at,\[Zeta]t,Pt,\[CapitalPhi]t,\[Phi]t}={atemp,\[Zeta]temp,Ptemp,\[CapitalPhi]temp,\[Phi]temp};(*save all variables at t+dt for the next step*)*)
(**)
(*(*End of each iteration:save data every TStep steps*) *)
(*If[Mod[tn,TStep]==1,*)
(*{att[tn],\[Zeta]tt[tn],Ptt[tn],\[CapitalPhi]tt[tn],\[Phi]tt[tn],\[Rho]tt[tn]}={at,\[Zeta]t,Pt,\[CapitalPhi]t,\[Phi]t,*)
(*1/2*(Pt^2+\[CapitalPhi]t^2)+(Q^2*(1+b0*\[Phi]t^2))/2 *(1-z)^4/z^4};*)
(*\[Zeta]cst[tn]=\[Zeta]dzConstraint[\[Zeta]t,Most[\[Zeta]t],Most[Pt],Most[\[CapitalPhi]t],Most[\[Phi]t]](*save constraints every TStep steps*)*)
(*];*)
(*tn=tn+1;(*record iteration numbers*)*)
(*]],{tn,tn*dt,\[Phi]t[[1]]}]//AbsoluteTiming*)


(* ::Input:: *)
(*{dt tn,tn}*)


(* ::Section:: *)
(* Analysis of numerical results*)


(* ::Input:: *)
(*tns=Range[1,Floor[tn/TStep-1]TStep+1,TStep];(*The iteration steps associated with all saved data*)*)
(*times=dt tns;(*The coordinate-time values associated with all saved data*)*)


(* ::Subsubsection:: *)
(*The  \|01d701 constraint violation*)


(* ::Input:: *)
(*(*The final\:2011time \|01d701 constraint violation*) *)
(*ListLinePlot[Rest[\[Zeta]cst[tns[[-1]]]]//Abs//Log,PlotRange->All];*)
(*ListLinePlot[(Rest[\[Zeta]cst[#]]//Abs//Log)&/@tns[[1;;-1;;Floor[Length[tns]/6]]],Frame->True,PlotRange->All]*)


(* ::Subsection:: *)
(*Metric and scalar\:2011field profiles*)


(* ::Subsubsection:: *)
(*Metric profile*)


(* ::Input:: *)
(*ListLinePlot[Transpose[{z,att[#]}]&/@tns[[-1;;1;;-Floor[Length[tns]/8]]],Frame->True,LabelStyle->{16,FontFamily->Times},FrameLabel->{"z","\[Alpha]"}](*Snapshots of the evolution of the metric function*)
(*\|01d6fc*)  *)


(* ::Input:: *)
(*\[Zeta]z[t_]:=Sqrt[(1-z)/z]\[Zeta]tt[t];(*Metric \[Zeta]/Sqrt[r]*)*)
(*gttz[t_]:=att[t]^2 (1-\[Zeta]z[t]^2);(*Metric coefficient of dt^2*) *)
(*gttr[tn0_]:=Transpose[{z2/(1-z2),Most[gttz[tn0]]}];(*Metric coefficient of dt^2*) *)
(*ListLinePlot[Transpose[{z,\[Zeta]z[#]}]&/@tns[[-1;;1;;-Floor[Length[tns]/8]]],Frame->True,LabelStyle->{15,FontFamily->Times},FrameLabel->{"z","\[Zeta]"},PlotRange->{All,{0,3}},PlotLegends->{"t9","t8","t7","t6","t5","t4","t3","t2","t1"}](*Snapshots of the evolution of the metric function*)
(*\|01d701*)  *)
(**)


(* ::Input:: *)
(**)
(**)


(* ::Subsubsection:: *)
(*The location of the outer horizon*)


(* ::Input:: *)
(*metric\[Alpha]rF[t_]:=Interpolation[gttr[t]];*)
(*irrBHrn=(x/.FindRoot[metric\[Alpha]rF[#][x]==0,{x,1.4rh}])&/@tns;*)


(* ::Input:: *)
(*irrBHtrn=Transpose[{times,1/2 irrBHrn}]; *)
(*irrBHtrf=Interpolation[irrBHtrn];*)
(*ListLinePlot[irrBHtrn,PlotRange->All,PlotStyle->Blue,Frame->True,LabelStyle->{15,FontFamily->Times},FrameLabel->{"t","\!\(\*SubscriptBox[\(r\), \(h\)]\)"}]*)
(*Plot[Log[irrBHtrf'[t]],{t,times[[1]],times[[-1]]}];*)


(* ::Subsection:: *)
(*The scalar field at the outer horizon*)


(* ::Subsubsection:: *)
(*Scalar-field value at the horizon*)


(* ::Input:: *)
(*\[Phi]tr[t_]:=Transpose[{z2/(1-z2),\[Phi]tt[t]//Most}];*)
(*\[Phi]tz[t_]:=Transpose[{z2,\[Phi]tt[t]//Most}];*)
(*\[Phi]trF[t_]:=Interpolation[\[Phi]tr[t]];*)
(*\[Phi]RH=\[Phi]trF[tns[[#]]][irrBHrn[[#]]]&/@Range[Length[tns]];*)


(* ::Input:: *)
(*\[Phi]RHtrn=Transpose[{times,\[Phi]RH}]; *)
(*\[Phi]RHtrf=Interpolation[\[Phi]RHtrn];*)
(*ListLinePlot[\[Phi]RHtrn,PlotRange->All,Frame->True,LabelStyle->{15,FontFamily->Times},FrameLabel->{"t","\!\(\*SubscriptBox[\(\[Phi]\), \(h\)]\)"}](*Scalar-field value at the horizon*)*)
(**)


(* ::Subsubsection:: *)
(*Scalar curvature*)


(* ::Input:: *)
(*Rtz[t_]:=Ptt[t]^2-\[CapitalPhi]tt[t]^2;(*Metric coefficient of dt^2*) *)
(*Rtr[tn0_]:=Transpose[{z2/(1-z2),Most[Rtz[tn0]]}];(*Metric coefficient of dt^2*) *)
(*ListLinePlot[Transpose[{z,Rtz[#]}]&/@tns[[-1;;1;;-Floor[Length[tns]/6]]],Frame->True,PlotRange->All]*)


(* ::Input:: *)
(*RtrF[t_]:=Interpolation[Rtr[t]];*)
(*RRH=RtrF[tns[[#]]][irrBHrn[[#]]]&/@Range[Length[tns]];*)


(* ::Input:: *)
(*ListLinePlot[Transpose[{times,irrBHrn/(1+irrBHrn)}],PlotRange->All,PlotStyle->Blue]*)
(*ListPlot[Transpose[{times,RRH}]];*)


(* ::Subsubsection:: *)
(*Scalar-field profile*)


(* ::Input:: *)
(*ListLinePlot[\[Phi]tr[#]&/@tns[[1;;-1;;Floor[Length[tns]/10]]],Frame->True,AxesOrigin->{rh,0},PlotRange->{{0,100},All}]*)
(*(*Snapshots of the evolution of the scalar field as a function of r*)*)


(* ::Input:: *)
(*ListLinePlot[\[Phi]tz[#]&/@tns[[1;;-1;;Floor[Length[tns]/10]]],Frame->True,PlotRange->All] (*Snapshots of the evolution of the scalar field as a function of z*)*)


(* ::Subsection:: *)
(*Mass distribution and the horizon*)


(* ::Subsubsection:: *)
(*Misner-Sharp mass*)


(* ::Input:: *)
(*misnerSharpz[t_]:=\[Zeta]tt[t]^2/2;(*Misner-Sharp mass*)*)
(*misnerSharpr[t_]:=Transpose[{z2/(1-z2),misnerSharpz[t]//Most}];(*Misner-Sharp mass*)*)
(*ListLinePlot[Transpose[{z,misnerSharpz[#]}]&/@tns[[1;;-1;;Floor[Length[tns]/6]]],Frame->True,PlotRange->All](*Snapshots of the evolution of the Misner\[Dash]Sharp mass*)*)
(**)


(* ::Subsection:: *)
(*Kretschmann scalar and Energy density*)


(* ::Subsubsection:: *)
(*Kretschmann scalar*)


(* ::Input:: *)
(*rr=Join[z2/(1-z2),{10^8}];*)
(*RuuuuRdddd[t_]:=8 (Ptt[t]^2-\[CapitalPhi]tt[t]^2)^2-(8 (\[Zeta]z[t]^2) )/rr^2 (Ptt[t]^2-\[CapitalPhi]tt[t]^2)+((12 \[Zeta]z[t]^4)/(rr^4))+((20 Q^4)/(rr^8 f[\[Phi]tt[t]]^2))+((16 Q^2 (Ptt[t]^2-\[CapitalPhi]tt[t]^2))/(rr^4 f[\[Phi]tt[t]]))-((24 Q^2 \[Zeta]z[t]^2)/(rr^6 f[\[Phi]tt[t]]))(*The formula for the Kretschmann scalar*)*)


(* ::Input:: *)
(*ListLinePlot[Transpose[{z,RuuuuRdddd[#]}]&/@tns[[-1;;1;;-Floor[Length[tns]/8]]],PlotRange->All,Frame->True,FrameLabel->{"z","K(z)"},LabelStyle->{15,FontFamily->Times},PlotLegends->{"t9","t8","t7","t6","t5","t4","3","t2","t1"}](*Snapshots of the evolution of the Kretschmann scalar*)*)


(* ::Input:: *)
(*RRtn=RuuuuRdddd[#][[1]]&/@tns;*)


(* ::Input:: *)
(*ListPlot[Transpose[{times,RRtn}],Frame->True,FrameLabel->{"t","K(\!\(\*SubscriptBox[\(z\), \(h\)]\))"},PlotRange->All,LabelStyle->{15,FontFamily->Times},PlotStyle->Darker[Green]](*Time evolution of the Kretschmann scalar inside the outer horizon*)*)


(* ::Subsubsection:: *)
(*Energy density*)


(* ::Input:: *)
(*ListLinePlot[Transpose[{z,\[Rho]tt[#]}]&/@tns[[-1;;1;;-Floor[Length[tns]/10]]],Frame->True,LabelStyle->{15,FontFamily->Times},FrameLabel->{"z","\[Rho]"},PlotLegends->{"t9","t8","t7","t6","t5","t4","t3","t2","t1"},PlotRange->All](*Snapshots of the evolution of the energy density*)*)


(* ::Input:: *)
(*\[Rho]ttn=\[Rho]tt[#][[1]]&/@tns;*)
(*ListPlot[Transpose[{times,\[Rho]ttn}],Frame->True,LabelStyle->{15,FontFamily->Times},FrameLabel->{"t","\!\(\*SubscriptBox[\(\[Rho]\), \(h\)]\)"},PlotStyle->Darker[Green],PlotRange->All](*Time evolution of the energy density inside the outer horizon*)*)
