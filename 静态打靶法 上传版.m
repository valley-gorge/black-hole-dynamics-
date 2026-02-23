(* ::Package:: *)

(* ::Section:: *)
(*static solution*)


(* ::Subsection:: *)
(*the function for solving the static equations*)


(* ::Subsubsection:: *)
(*Static equations and the asymptotic expansion at the horizon*)


(* ::Input:: *)
(*b=-35.0;Q=0.448;*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*f[x_]=1/(1+b x^2);*)
(*eqd\[Zeta]dr=Derivative[1][\[Zeta]][r]-(Q^2/(2 r^3 f[\[Phi][r]] \[Zeta][r])-\[Zeta][r] /(2 r)+r/2 (1 /\[Zeta][r]-\[Zeta][r] )Derivative[1][\[Phi]][r]^2);*)
(*eqd\[Phi]dr=(\[Phi]^\[Prime]\[Prime])[r]-(-((Q^2 Derivative[1][f][\[Phi][r]])/(2 r^4 f[\[Phi][r]]^2))/(1-\[Zeta][r]^2)+(2 \[Zeta][r] Derivative[1][\[Zeta]][r] Derivative[1][\[Phi]][r])/(1-\[Zeta][r]^2)- (2/r +r Derivative[1][\[Phi]][r]^2) Derivative[1][\[Phi]][r]);(* The constraint equations used to solve for \[Zeta] and \[Phi] *)*)
(**)
(**)
(*BCH\[Phi]r[\[Phi]h_,r_]=\[Phi]h+(Q^2 Derivative[1][f][\[Phi]h] (r-rh))/(2 Q^2 rh f[\[Phi]h]-2 rh^3 f[\[Phi]h]^2)+(Q^2 Derivative[1][f][\[Phi]h] ((Q^2-2 rh^2 f[\[Phi]h]) (4 f[\[Phi]h]^2 (Q^2-rh^2 f[\[Phi]h])+Q^2 Derivative[1][f][\[Phi]h]^2)+Q^2 f[\[Phi]h] (-Q^2+rh^2 f[\[Phi]h]) (f^\[Prime]\[Prime])[\[Phi]h]) (r-rh)^2)/(16 rh^2 f[\[Phi]h]^3 (-Q^2+rh^2 f[\[Phi]h])^3);(*scalar field near\:2011horizon expansion*)*)
(**)
(**)
(*BCH\[Zeta]r[\[Phi]h_,r_]=1+((-rh^2+Q^2/f[\[Phi]h]) (r-rh))/(2 rh^3)+((6 rh^4+(Q^2 (-2 f[\[Phi]h] (Q^2+6 rh^2 f[\[Phi]h])-(3 Q^2 rh^2 Derivative[1][f][\[Phi]h]^2)/(Q^2-rh^2 f[\[Phi]h])))/f[\[Phi]h]^3) (r-rh)^2)/(16 rh^6);(*metric near\:2011horizon expansion*)*)


(* ::Subsubsection:: *)
(*The function for solving the static equations*)


(* ::Input:: *)
(*\[Epsilon]=10^-7;INF=10^7;rh\[Epsilon]=rh+\[Epsilon];nn=30;*)
(*Off[NDSolve::precw];*)
(*\[Zeta]\[Phi]r[\[Phi]h_?NumberQ]:={\[Zeta],\[Phi]}/.NDSolve[{eqd\[Zeta]dr==0, eqd\[Phi]dr==0,*)
(*\[Zeta][rh\[Epsilon]]==BCH\[Zeta]r[\[Phi]h,rh\[Epsilon]],*)
(*(*\[Zeta]'[rh\[Epsilon]]\[Equal]D[BCH\[Zeta]r[\[Phi]h,r],r]/.r->rh\[Epsilon],*)*)
(*\[Phi][rh\[Epsilon]]==BCH\[Phi]r[\[Phi]h,rh\[Epsilon]],*)
(*\[Phi]'[rh\[Epsilon]] ==D[BCH\[Phi]r[\[Phi]h,r],r]/.r->rh\[Epsilon]},*)
(* {\[Zeta],\[Phi]},{r, rh\[Epsilon],INF},WorkingPrecision->nn,PrecisionGoal->nn][[1]]*)
(**)
(*\[Phi]r[\[Phi]h_?NumberQ]:=\[Phi]/.NDSolve[{eqd\[Zeta]dr==0, eqd\[Phi]dr==0,*)
(*\[Zeta][rh\[Epsilon]]==BCH\[Zeta]r[\[Phi]h,rh\[Epsilon]],*)
(*(*\[Zeta]'[rh\[Epsilon]]\[Equal]D[BCH\[Zeta]r[\[Phi]h,r],r]/.r->rh\[Epsilon],*)*)
(*\[Phi][rh\[Epsilon]]==BCH\[Phi]r[\[Phi]h,rh\[Epsilon]],*)
(*\[Phi]'[rh\[Epsilon]] == D[BCH\[Phi]r[\[Phi]h,r],r]/.r->rh\[Epsilon]},*)
(* {\[Zeta],\[Phi]},{r, rh\[Epsilon],INF},WorkingPrecision->nn,PrecisionGoal->nn][[1]]*)


(* ::Subsubsection:: *)
(*The shooting\:2011method routine and its testing*)


(* ::Input:: *)
(**)
(*\[Phi]INF[\[Phi]h_]:=\[Phi]r[\[Phi]h][INF];*)
(*\[Phi]Hshoot[phihSeed_] :=*)
(*\[Phi]h/.FindRoot[\[Phi]INF[\[Phi]h]==10^-7, {\[Phi]h,phihSeed}, WorkingPrecision->nn,PrecisionGoal->nn]*)


(* ::Subsection:: *)
(*Computing the charge\:2011to\:2011mass ratio using the shooting method*)


(* ::Input:: *)
(*rhlist={};entropylist={};phiHlist={};ADMlist={};Qslist={};qlist={};Slist={};*)


(* ::Subsubsection:: *)
(*Solving the equations in an iterative loop*)


(* ::Input:: *)
(*\[Phi]h0=-10944.8/10000; *)
(*targetPhi=1/Sqrt[-b]; (*evaluating \[Phi] along the divergent line*)*)
(**)
(*Monitor[For[rh=15960/10000,rh>120000/100000,rh=rh-1/1000,*)
(*{phiH=Rationalize[\[Phi]Hshoot[\[Phi]h0]];*)
(*\[Zeta]\[Phi]rsol=\[Zeta]\[Phi]r[phiH];*)
(*ADM=INF/2 (\[Zeta]\[Phi]rsol[[1]])[INF]^2;*)
(*Qs=-(\[Zeta]\[Phi]rsol[[2]])'[INF] INF^2;*)
(*q=Q/ADM;*)
(**)
(*If[Abs[phiH+targetPhi]<0.01,Print["q at \!\(\*SuperscriptBox[\(b\), \(\(-1\)/2\)]\) = ",1. q]];(*check if the current \[Phi]h is sufficiently close to the target*)*)
(*AppendTo[rhlist,1. rh];*)
(*AppendTo[entropylist,1. \[Pi] rh^2];*)
(*AppendTo[phiHlist,phiH];*)
(*AppendTo[ADMlist,ADM];*)
(*AppendTo[Qslist,Qs];*)
(*AppendTo[qlist,1. q];*)
(*AppendTo[Slist,1. (\[Pi] rh^2)/(4\[Pi] ADM^2)]; *)
(*\[Phi]h0=phiH;}],1.{rh,ADM,(\[Pi] rh^2)/(4\[Pi] ADM^2),phiH,Qs,q}] *)


(* ::Section:: *)
(*Analysis of the static\:2011solution data*)


(* ::Input:: *)
(*rhlist;(*horizon radius*)*)


(* ::Input:: *)
(*phiHlist;(*the scalar\:2011field value*)*)


(* ::Input:: *)
(*qlist(*the dimensionless charge\:2011to\:2011mass ratio*)*)


(* ::Input:: *)
(*rhlist[[57]]*)


(* ::Input:: *)
(*qlist[[57]]*)


(* ::Input:: *)
(*1.54/(1+1.54)(*change of coordinates*) *)


(* ::Input:: *)
(*rh=rhlist[[57]];*)
(*rh\[Epsilon]=rh+\[Epsilon];*)
(*{\[Zeta]rsol,\[Phi]rsol}=\[Zeta]\[Phi]r[phiHlist[[57]]];*)
(*Plot[\[Phi]rsol[z/(1-z)],{z,0.6062993125984251`,0.9999745727088383},Frame->True,PlotRange->All](*scalar\:2011field profile*)*)
(*Plot[\[Zeta]rsol[z/(1-z)],{z,0.6062993125984251`2,0.9999745727088383},Frame->True,PlotRange->{All,{0,1}}](*metric\:2011function profile*)*)
(**)
(**)



