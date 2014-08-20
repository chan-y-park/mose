(* ::Package:: *)

(* ::Subsection:: *)
(*Preamble*)


BeginPackage["ConstructKWalls`"]

(*<<"/Users/pietro/Documents/Rutgers University/Work for Greg/k3 metrics/k3_notes/Software/K3Constructor_v2.m";
<<"/Users/pietro/Documents/Rutgers University/Work for Greg/k3 metrics/k3_notes/Software/CurveIntersection.m";
<<"/Users/pietro/Documents/Rutgers University/Work for Greg/k3 metrics/k3_notes/Software/KSWCF.m";*)

<<"Combinatorica`"//Quiet

BuildKWalls::usage="BuildKWalls[n_] where n is the number of iterations 
for generating children walls. 
To plot use the command PlotKWalls. 
See also the commands trajs, paths, ints, genealogy, 
singpts, PlotMSWalls, MakeVideo, ScanTheta,
AllParticles, BPSspectrum, MSpts"
gg2::usage=""
gg3::usage=""
\[Theta]::usage=""
step::usage=""
\[CapitalLambda]0::usage=""
\[CapitalLambda]::usage=""
MaxRank::usage=""
singpts::usage="List of singularities on the Coulomb branch"
trajs::usage="Organized as follows {{ {u's}, {\[Eta]1's}, \[Gamma] , \[CapitalOmega] },...}"
paths::usage="Organized as follows {{u's},...}"
PlotKWalls::usage=""
ints::usage="Organized as follows {{tr1,poin1},{tr2,pt2},u},...} 
For example {{{1,1161},{2,1139},5.742-6.706 \[ImaginaryI]},...}"
genealogy::usage="Keeps info on who are the parents of a K-wall, 
the primary walls have no parents.
For example {{},{},{},{1,2},{1,3}}"
verb::usage="True or False, toggles flags during run."
pltrng::usage="The plot range."
MSpts::usage="Array containing lists of points of MS walls"
ScanTheta::usage="Use as ScanTheta[steps, filename, generations] 
by giving the positive integer steps, the string containing the filename, 
the number of generations." 
PlotMSWalls::usage="Use only after ScanTheta"
MakeVideo::usage="Use only after ScanTheta"
NPrimKWalls::usage="Number of primary K-walls"
AllParticles::usage="Tells the particle content of the K-walls at given \[Theta]"
BPSspectrum::usage="Spectrum[u_,range_] Tells the spectrum at u, must be used 
after ScanTheta. range is a parameter that tells how fine the graining of 
CB should be"
EvolutionData::usage="Stores data about evolution with \[Theta]"
showlabels::usage="toggles labels in plot"
PlotSingularities::usage="Plots the singularities"
AllCentralCharges::usage="After scanning \[Theta], this stores info on central charges 
along PRIMARY K-walls in the format {\!\(\*SubscriptBox[\(list\), \(1\)]\),\!\(\*SubscriptBox[\(list\), \(2\)]\),...,\!\(\*SubscriptBox[\(list\), \(n\)]\)} (i-th entry 
for i-th charge), where \!\(\*SubscriptBox[\(list\), \(1\)]\) is for \!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\), ect. and \!\(\*SubscriptBox[\(list\), \(k\)]\)={{\!\(\*SubscriptBox[\(u\), \(1\)]\),Z(\!\(\*SubscriptBox[\(u\), \(1\)]\))},...}"
PlotCentralCharge::usage="Use as PlotCentralCharge[n] where n=1,...,NPrimKWalls/2. 
This plots the real and imaginary parts of \!\(\*SubscriptBox[\(Z\), \(n\)]\)."


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Kwalls-constructor*)


\[CapitalDelta]=u\[Function]Evaluate[Expand[gg2[u]^3-27 gg3[u]^2]];
PlotSingularities:=Module[{z,points},
points=z/.NSolve[\[CapitalDelta][z]==0,z];
Print["Singularities at  ",points//Column];
ListPlot[(X\[Function]{X//Re,X//Im})/@points,PlotRange->All]
];

\[Delta]=u\[Function]Evaluate[Expand[3g3[u]g2'[u]-2g2[u]g3'[u]]];

Pair=X\[Function](R \[Function]{X[[R[[2]]]],X[[R[[3]]]]})[Module[{f},Sort[Flatten[Array[(j \[Function]Array[(i \[Function]f[i,j]),Length[X]-j,j+1]),Length[X]-1]]/.{f -> ({x,y}\[Function]{Abs[X[[x]]-X[[y]]],x,y})}, #1[[1]]<=#2[[1]]&]][[1]]];

\[Eta]1A=u\[Function]Module[{e1,e2,e3,x,sl,S},
sl=x/.Solve[4x^3-g2[u]x-g3[u]==0,x];
(X\[Function](e1=X[[1]];e2=X[[2]];e3=Complement[sl,X][[1]]))[Pair[sl]];
2/(e3-e1)^(1/2) EllipticK[(e2-e1)/(e3-e1)]
];

\[Eta]1B=u\[Function]Module[{e1,e2,e3,x,sl},
sl=x/.NSolve[4x^3-g2[u]x-g3[u]==0,x];
(X\[Function](e1=X[[1]];e2=X[[2]];e3=Complement[sl,X][[1]]))[Pair[sl]];
(2I)/(e3-e1)^(1/2) EllipticK[(e3-e2)/(e3-e1)]
];

ShortDistance={u,\[CurlyTheta],L,d0,\[Alpha]}\[Function]Module[{S,N},
N=IntegerPart[L/d0];
S[0]={u,\[Alpha] \[Eta]1A[u]};
Do[ 
S[k]=(X\[Function]{X,If[Re[\[Eta]1A[X]/S[k-1][[2]]]>0,1,-1]\[Eta]1A[X]})[S[k-1][[1]]+E^(I(\[CurlyTheta]-Arg[S[k-1][[2]]])) d0];
,{k,1,N}];
{Array[S[#][[1]]&,N],Array[S[#][[2]]&,N]}
];

(*PicardFuchs={\[Eta]0,\[Chi]0,\[CurlyTheta],u0,L}\[Function]Module[{F,\[Eta],\[Chi],z,u,t,\[Eta]R,\[Eta]I,\[Chi]R,\[Chi]I,uR,uI,sl,Eqq},
F={{0,1},{-((3 g2[z] \[Delta][z]^2)/(16 \[CapitalDelta][z]^2))+(Derivative[1][\[Delta]][z] Derivative[1][\[CapitalDelta]][z])/(12 \[Delta][z] \[CapitalDelta][z])+Derivative[1][\[CapitalDelta]][z]^2/(144 \[CapitalDelta][z]^2)-Derivative[2][\[CapitalDelta]][z]/(12 \[CapitalDelta][z]),Derivative[1][\[Delta]][z]/\[Delta][z]-Derivative[1][\[CapitalDelta]][z]/\[CapitalDelta][z]}};
Eqq=Flatten[(X\[Function]{Evaluate[ComplexExpand[Re[Expand[Evaluate[X]]]]]==0,Evaluate[ComplexExpand[Im[Expand[Evaluate[X]]]]]==0})/@(Flatten[{Numerator[Factor[\[Eta][t]\!\(
\*SubscriptBox[\(\[PartialD]\), \(t\)]\({\[Eta][t], \[Chi][t]}\)\)-E^(I \[CurlyTheta]) Abs[\[Eta][t]](F/.{z->u[t]}).{\[Eta][t],\[Chi][t]}]],\[Eta][t]u'[t]-Abs[\[Eta][t]]E^(I \[CurlyTheta]),\[Eta][0]-\[Eta]0,\[Chi][0]-\[Chi]0,u[0]-u0}]/.{\[Eta]->(t\[Function]\[Eta]R[t]+I \[Eta]I[t]),\[Chi]->(t\[Function]\[Chi]R[t]+I \[Chi]I[t]),u->(t\[Function]uR[t]+I uI[t])})];
(*Print[Eqq];*)
Print["Eqq is ready"];
sl=NDSolve[Eqq,{\[Eta]R,\[Eta]I,\[Chi]R,\[Chi]I,uR,uI},{t,0,L}];
(X\[Function]Array[X[[2#-1]]&,Length[X]/2])/@{((uR/.sl)[[1,4,3]])+I((uI/.sl)[[1,4,3]]),((\[Eta]R/.sl)[[1,4,3]])+I ((\[Eta]I/.sl)[[1,4,3]])}
];*)

(*PicardFuchs={\[Eta]0,\[Chi]0,\[CurlyTheta],u0,L}\[Function]Module[{F,\[Eta],\[Chi],z,u,t,\[Eta]R,\[Eta]I,\[Chi]R,\[Chi]I,uR,uI,sl},
F={{0,1},{-((3 g2[z] \[Delta][z]^2)/(16 \[CapitalDelta][z]^2))+(Derivative[1][\[Delta]][z] Derivative[1][\[CapitalDelta]][z])/(12 \[Delta][z] \[CapitalDelta][z])+Derivative[1][\[CapitalDelta]][z]^2/(144 \[CapitalDelta][z]^2)-(\[CapitalDelta]'')[z]/(12 \[CapitalDelta][z]),Derivative[1][\[Delta]][z]/\[Delta][z]-Derivative[1][\[CapitalDelta]][z]/\[CapitalDelta][z]}};
(*Print[F/.{z->u0}];*)
sl=NDSolve[Flatten[(X\[Function]{ComplexExpand[Re[X]]==0,ComplexExpand[Im[X]]==0})/@(Flatten[{Numerator[Factor[\[Eta][t]\!\(
\*SubscriptBox[\(\[PartialD]\), \(t\)]\({\[Eta][t], \[Chi][t]}\)\)-E^(I \[CurlyTheta]) Abs[\[Eta][t]](F/.{z->u[t]}).{\[Eta][t],\[Chi][t]}]],\[Eta][t]u'[t]-Abs[\[Eta][t]]E^(I \[CurlyTheta]),\[Eta][0]-\[Eta]0,\[Chi][0]-\[Chi]0,u[0]-u0}]/.{\[Eta]->(t\[Function]\[Eta]R[t]+I \[Eta]I[t]),\[Chi]->(t\[Function]\[Chi]R[t]+I \[Chi]I[t]),u->(t\[Function]uR[t]+I uI[t])})],{\[Eta]R,\[Eta]I,\[Chi]R,\[Chi]I,uR,uI},{t,0,L}(*,AccuracyGoal->60,WorkingPrecision->80,PrecisionGoal->60,*)];
(X\[Function]Array[X[[2#-1]]&,Length[X]/2])/@{((uR/.sl)[[1,4,3]])+I((uI/.sl)[[1,4,3]]),((\[Eta]R/.sl)[[1,4,3]])+I ((\[Eta]I/.sl)[[1,4,3]])}
];*)

PicardFuchs={\[Eta]0,\[Chi]0,\[CurlyTheta],u0,L}\[Function]Module[{F,\[Eta],\[Chi],z,u,t,sl},
F={{0,1},{-((3 g2[z] \[Delta][z]^2)/(16 \[CapitalDelta][z]^2))+(Derivative[1][\[Delta]][z] Derivative[1][\[CapitalDelta]][z])/(12 \[Delta][z] \[CapitalDelta][z])+Derivative[1][\[CapitalDelta]][z]^2/(144 \[CapitalDelta][z]^2)-(\[CapitalDelta]'')[z]/(12 \[CapitalDelta][z]),Derivative[1][\[Delta]][z]/\[Delta][z]-Derivative[1][\[CapitalDelta]][z]/\[CapitalDelta][z]}};
(*Print[F/.{z->u0}];*)
sl=NDSolve[Flatten[{u'[t]==E^(I  \[CurlyTheta]) Abs[\[Eta][t]]/\[Eta][t],({X,Y}\[Function]{X[[1]]==Y[[1]],X[[2]]==Y[[2]]})[\!\(
\*SubscriptBox[\(\[PartialD]\), \(t\)]\({\[Eta][t], \[Chi][t]}\)\),E^(I  \[CurlyTheta]) Abs[\[Eta][t]]/\[Eta][t] (F/.{z->u[t]}).{\[Eta][t],\[Chi][t]}],u[0]==u0,\[Eta][0]== \[Eta]0,\[Chi][0]==\[Chi]0}],{u,\[Eta],\[Chi]},{t,0,L}];
(X\[Function]Array[X[[2#-1]]&,Length[X]/2])/@{(u/.sl)[[1,4,3]],(\[Eta]/.sl)[[1,4,3]]}
];

ConstructPrimKWalls={\[CurlyTheta],d0,l}\[Function]Module[{U,u,S},
U=DeleteDuplicates[N[u/.NSolve[\[CapitalDelta][u]==0,u]],Abs[#1-#2]<l&];
(*U=N[u/.NSolve[\[CapitalDelta][u]==0,u]];*)
If[verb,Print[U]];
Do[
S[2k-1]=ShortDistance[U[[k]],\[CurlyTheta],l,d0,1];
S[2k]=ShortDistance[U[[k]],\[CurlyTheta],l,d0,-1];
,{k,1,Length[U]}];
Array[S[#]&,2Length[U]]
];

Grow={\[CurlyTheta],L,Traj}\[Function]Module[{S},
Do[
(*Print["Trajectory #",k," will grow"];*)
S[k]=(Y\[Function]{Join[Traj[[k,1]],Y[[1]]],Join[Traj[[k,2]],Y[[2]]]})[(X\[Function]PicardFuchs[X[[2,Length[X[[2]]]]],(X[[2,Length[X[[2]]]]]-X[[2,Length[X[[2]]]-1]])/(X[[1,Length[X[[1]]]]]-X[[1,Length[X[[1]]]-1]]),\[CurlyTheta],X[[1,Length[X[[1]]]]],L])[Traj[[k]]]];
(*Print["Trajectory #",k," is grown"];*)
,{k,1,Length[Traj]}];
Array[S[#]&,Length[Traj]]
];

DecomposePrs={u,\[Eta]}\[Function]Module[{a,b,X},
X={a,b}/.Solve[(X\[Function]{ComplexExpand[Re[X]]==0,ComplexExpand[Im[X]]==0})[\[Eta]-a \[Eta]1A[u]-b \[Eta]1B[u]],{a,b}][[1]];
(*Print["before round: ",X];*)
(*Print["Periods are ", {\[Eta]1A[u],\[Eta]1B[u]}]*)
X//Round
];



(* ::Subsection::Closed:: *)
(*Curve Intersections*)


Swtch=X\[Function]Module[{N,h,t},
N=0;
h=Sign[X[[2]]-X[[1]]];
Do[If[Sign[X[[k+1]]-X[[k]]]!=h,h=Sign[X[[k+1]]-X[[k]]];N=N+1;t[N]=k; ];
,{k,2,Length[X]-1}];
{1}\[Union]Array[t,N]\[Union]{Length[X]}];

DoX=\[Gamma]\[Function]Module[{T,N,f,\[Tau]},
T=Swtch[(X\[Function]Re[X])/@\[Gamma]];
N=0;
Do[N=N+1;
f[N]=Interpolation[Table[{Re[\[Gamma][[t]]],Im[\[Gamma][[t]]]},{t,T[[k]]+1,T[[k+1]]}]];
\[Tau][N]=Interpolation[Table[{Re[\[Gamma][[t]]],t},{t,T[[k]]+1,T[[k+1]]}]]
,{k,1,Length[T]-1}];
{Array[f,N],Array[\[Tau],N]}
];

FindIntersections={f1,f2}\[Function]Module[{xx,Mn,Mx,S,x},
S={};
Mn=Max[f1[[1,1,1]],f2[[1,1,1]]];
Mx=Min[f1[[1,1,2]],f2[[1,1,2]]];
While[Mn+Pr<Mx,
xx=FindMinimum[(f1[x]-f2[x])^2,{x,Mn+Pr,Mn,Mx}];
If[xx[[1]]<Pr,S=S\[Union]{x/.xx[[2]]}];
Mn=x/.xx[[2]];
If[Mn+Pr<Mx,
xx=FindMinimum[-(f1[x]-f2[x])^2,{x,Mn+Pr,Mn,Mx}];
Mn=x/.xx[[2]]];
];
S
];

FindIntersection={\[Gamma]1,\[Gamma]2}\[Function]Module[{T1,T2,S,C},
S={};
T1=DoX[\[Gamma]1];
T2=DoX[\[Gamma]2];
f1=T1[[1]];
f2=T2[[1]];
Do[Do[
C=FindIntersections[f1[[k]],f2[[m]]];
Do[
S=S\[Union]{{IntegerPart[T1[[2,k]][C[[l]]]],IntegerPart[T2[[2,m]][C[[l]]]]}};
,{l,1,Length[C]}];
,{m,1,Length[f2]}],{k,1,Length[f1]}];
If[Length[S]>0,S=(X\[Function]{Round[X[[1]]],Round[X[[2]]],\[Gamma]1[[Round[X[[1]]]]]})/@(X\[Function]{t1,t2}/.FindRoot[{Re[ListInterpolation[\[Gamma]1][t1]-ListInterpolation[\[Gamma]2][t2]]==0,Im[ListInterpolation[\[Gamma]1][t1]-ListInterpolation[\[Gamma]2][t2]]==0},{t1,X[[1]],1,Length[\[Gamma]1]},{t2,X[[2]],1,Length[\[Gamma]2]}])/@S];
S
];



(* ::Subsection::Closed:: *)
(*KSWCF*)


Unprotect[Wedge];
\[Alpha]_\[Wedge]\[Beta]_:=\[Alpha].{{0,-m},{m,0}}.\[Beta];
Protect[Wedge];

Unprotect[Times,Power];
Subscript[X, \[Alpha]_List] Subscript[X, \[Beta]_List]:=(-1)^(\[Alpha]\[Wedge]\[Beta]) Subscript[X, \[Alpha]+\[Beta]];
\!\(
\*SubsuperscriptBox[\(X\), \(\[Alpha]_List\), \(a_Integer\)] := 
\*SubscriptBox[\(X\), \(a\ \[Alpha]\)]\);
Subscript[K, \[Alpha]_List]V_:=V/.{Subscript->({X,\[Beta]}\[Function](1-Subscript[X, \[Alpha]])^(\[Beta]\[Wedge]\[Alpha]) Subscript[X, \[Beta]])};
\!\(
\*SubsuperscriptBox[\(K\), \(\[Alpha]_List\), \(\[CapitalOmega]_Integer\)] V_\):=V/.{Subscript->({X,\[Beta]}\[Function](1-Subscript[X, \[Alpha]])^(\[CapitalOmega] \[Beta]\[Wedge]\[Alpha]) Subscript[X, \[Beta]])};
H[R_]:=R/.{Subscript->({X,\[Alpha]}\[Function]










\!\(\*SuperscriptBox[\(t\), \(
\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(Length[\[Alpha]]\)]Abs[\[Alpha]\[LeftDoubleBracket]i\[RightDoubleBracket]]\)]\) Subscript[X, \[Alpha]])};
Subscript[\[Kappa], \[Alpha]_List][h_Integer]V_:=Expand[Normal[Series[H[V/.{Subscript->({X,\[Beta]}\[Function](1-Subscript[X, \[Alpha]])^(\[Beta]\[Wedge]\[Alpha]) Subscript[X, \[Beta]])}],{t,0,h}]]/.{t->1}];
\!\(
\(\*SubsuperscriptBox[\(\[Kappa]\), \(\[Alpha]_List\), \(\[CapitalOmega]_Integer\)]\)[h_Integer]\)V_:=Expand[Normal[Series[H[V/.{Subscript->({X,\[Beta]}\[Function](1-Subscript[X, \[Alpha]])^(\[CapitalOmega] \[Beta]\[Wedge]\[Alpha]) Subscript[X, \[Beta]])}],{t,0,h}]]/.{t->1}];

Subscript[pK, \[Alpha]_]V_:=V/.{Subscript->({X,\[Beta]}\[Function](-(\[Beta]\[Wedge]\[Alpha]) Subscript[X, \[Alpha]])Subscript[X, \[Beta]])};
\!\(
\*SubsuperscriptBox[\(pK\), \(\[Alpha]_\), \(\[CapitalOmega]_\)] V_\):=V/.{Subscript->({X,\[Beta]}\[Function](-(\[CapitalOmega] \[Beta]\[Wedge]\[Alpha]) Subscript[X, \[Alpha]])Subscript[X, \[Beta]])};

\[CapitalPi]K[S_,V_,h_]:=Module[{m1,H},
H=V;
Do[ m1=Length[S]-k+1;
H=\!\(
\(\*SubsuperscriptBox[\(\[Kappa]\), \(S[\([m1, 1]\)]\), \(S[\([m1, 2]\)]\)]\)[h]\)H;
,{k,1,Length[S]}];
H
];

Protect[Times,Power];

Needs["Combinatorica`"]//Quiet;

ConstParts[n_,k_]:=Module[{f,\[CapitalTheta],P},
\[CapitalTheta]=x\[Function]If[x>0,1,0];
P=Partitions[n,k];
Flatten[(X\[Function]f/@Permutations[X])/@((X\[Function]Array[i\[Function]\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(j = 1\), \(Length[X]\)]\(\[CapitalTheta][X[\([j]\)] - i + 1]\)\),k])/@P)]/.{f->(X\[Function]X)}
];

Transl[f_]:=Module[{N,A},
N=0;
f/.{Subscript->({x,y}\[Function](N=N+1;A[N]=y;Subscript[x, y]))};
Table[A[k],{k,1,N}]
];

Spec[S1_,S2_,V_,h_]:=Module[{f,R,\[CapitalOmega],S},
R=(\[CapitalPi]K[S1,V,h]-\[CapitalPi]K[S2,V,h]);
S=DeleteDuplicates[Flatten[(X\[Function](Y\[Function]f[Y])/@X)/@Transl/@Expand/@(R((X\[Function]X^-1)/@V))]]/.{f->(X\[Function]X)};
S=Array[k\[Function]{S[[k]],\[CapitalOmega][k]},Length[S]];
R=R-(Y\[Function]\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(i = 1\), \(Length[S]\)]\(
\*SubsuperscriptBox[\(pK\), \(S[\([i, 1]\)]\), \(S[\([i, 2]\)]\)] Y\)\))/@V;
S/.Solve[(H\[Function]H==0)/@Flatten[(E\[Function]Coefficient[E,(y\[Function]Subscript[X, y[[1]]])/@DeleteDuplicates[Flatten[(X\[Function](Y\[Function]f[Y])/@X)/@Transl/@R]]/.{f->(X\[Function]X)}])/@R],(X\[Function]X[[2]])/@S][[1]]
];

Spectra[Zs_,h_,\[CapitalOmega]1_,\[CapitalOmega]2_]:=Module[{S={{{1,0},\[CapitalOmega]1},{{0,1},\[CapitalOmega]2}},S1,H,Z,r,W},
Z=\[Gamma]\[Function]\[Gamma].Zs;
S1={};
W=(y\[Function]Subscript[X, y])/@ConstParts[1,Length[Zs]];
Do[H=Spec[S,S1,W,r];
S1=Sort[Join[S1,H],Arg[Z[#1[[1]]]]>Arg[Z[#2[[1]]]]&];
,{r,2,h}];
S1//DeleteDuplicates
];

KSWCF[DSZ_,MaxRank_,\[CapitalOmega]1_,\[CapitalOmega]2_]:=Module[{Zs={1,I}},
m=DSZ;
Spectra[Zs,MaxRank+1,\[CapitalOmega]1,\[CapitalOmega]2]
];


(* ::Subsection::Closed:: *)
(*Global structure*)


(*ConstructKWalls subpackage begins here*)

uPt[{tr_Integer,pt_Integer}]:=trajs[[tr,1,pt]];
\[Eta]1[{tr_Integer,pt_Integer}]:=trajs[[tr,2,pt]];
\[CapitalOmega][tr_Integer]:=trajs[[tr,Length[trajs[[tr]]]]];
TrajCharge[tr_Integer]:=trajs[[tr,3]];
AllParticles:=Table[{trajs[[i,3]],trajs[[i,4]]},{i,1,Length[trajs]}];

LabelPlot:=Table[Graphics[Text[trajs[[i,3]],(trajs[[i,1,Floor[Length[trajs[[i,1]]]/7]]]//(X\[Function]{X//Re,X//Im}))]],{i,1,Length[trajs]}];
PlotKWalls:=Join[(Z\[Function]ListLinePlot[Z,PlotRange->pltrng,Axes->False,ImageSize->500,AspectRatio->1])/@(X\[Function](Y\[Function]{Y//Re,Y//Im})/@X)/@paths,Table[Graphics[{Red,PointSize[Large],Point[{singpts[[i]]//Re,singpts[[i]]//Im}]}],{i,1,Length[singpts]}],If[showlabels,LabelPlot,{}]]//Show;

initialize:=Module[{traj={},uu},
g2[x_]:=gg2[x];
g3[x_]:=gg3[x];
traj=If[verb,ConstructPrimKWalls[\[Theta],step,\[CapitalLambda]0],ConstructPrimKWalls[\[Theta],step,\[CapitalLambda]0]//Quiet];
traj=If[verb,Grow[\[Theta],\[CapitalLambda],traj],Grow[\[Theta],\[CapitalLambda],traj]//Quiet];
NPrimKWalls=Length[traj];
If[verb,Print["there are ",NPrimKWalls," primary walls"]];
singpts=uu/.NSolve[\[CapitalDelta][uu]==0,uu];
Table[traj[[2i-1,1,1]],{i,1,Length[traj]/2}];
trajs=Table[Join[traj[[i]]//DeleteDuplicates,{(-1)^(i+1) Subscript[\[Gamma], If[OddQ[i],(i+1)/2,i/2]],1}],{i,1,Length[traj]}];
paths=(X\[Function]X[[1]])/@trajs;

Subscript[\[Gamma], i_]:=Array[KroneckerDelta[#,i]&,NPrimKWalls/2];

Pr=10^(-10);

nT=0;
nInt=0; (* number of intersections that have been computed *)

ints={};
genealogy=Array[{}&,Length[trajs]]; 
];


siblings[i_,j_]:=TrueQ[(genealogy[[i]]==genealogy[[j]])&&(genealogy[[j]]!={})];
parentchild[i_,j_]:=(MemberQ[genealogy[[i]],j]||MemberQ[genealogy[[j]],i]);
related[i_,j_]:=(parentchild[i,j]||siblings[i,j]);


npts=8;
RefineInts[a_]:=Module[{Mn1=Max[1,a[[1,2]]-npts],Mx1=Min[Length[trajs[[a[[1,1]],1]]],a[[1,2]]+npts],Mn2=Max[1,a[[2,2]]-npts],Mx2=Min[Length[trajs[[a[[2,1]],1]]],a[[2,2]]+npts],f1x,f1y,f2x,f2y,pts1,pts2,refint,r0=Re[a[[3]]],t0,t3},
pts1=(X\[Function]{X//Re,X//Im})/@Table[trajs[[a[[1,1]],1,i]],{i,Mn1,Mx1}];
pts2=(X\[Function]{X//Re,X//Im})/@Table[trajs[[a[[2,1]],1,i]],{i,Mn2,Mx2}];
f1x=x\[Function]Evaluate[Fit[(X\[Function]X[[1]])/@pts1,{1,x,x^2,x^3},x]];
f1y=x\[Function]Evaluate[Fit[(X\[Function]X[[2]])/@pts1,{1,x,x^2,x^3},x]];
f2x=x\[Function]Evaluate[Fit[(X\[Function]X[[1]])/@pts2,{1,x,x^2,x^3},x]];
f2y=x\[Function]Evaluate[Fit[(X\[Function]X[[2]])/@pts2,{1,x,x^2,x^3},x]];
r0=a[[3]]//Re;
t0=Module[{x},x/.FindRoot[f1x[x]-r0,{x,0}]];
t3=Module[{x},x/.FindRoot[Abs[f1x[x]+I f1y[x]-(f2x[x]+I f2y[x])],{x,t0}]];
refint=(f1x[t3]+f2x[t3]+I (f1y[t3]+f2y[t3]))/2;
ReplacePart[a,3->refint]
]


iteration:=Module[{traj1,traj2,newcharges,newints,\[CapitalOmega]1,\[CapitalOmega]2,p,q,\[Gamma]1,\[Gamma]2,\[Omega]1,\[Omega]1prime,\[CapitalDelta]\[Omega]1(*,\[Omega]2*),u0,newtraj={},m=0,new\[CapitalOmega]},
newints=(RefineInts/@Select[Flatten[Table[Table[(X\[Function]{{i,X[[1]]},{j,X[[2]]},X[[3]],0})/@If[!related[i,j],FindIntersection[paths[[i]],paths[[j]]]//If[verb,(Z\[Function]Z),Quiet],{}],{j,Max[i+1,nT+1],Length[paths]}],{i,1,Length[paths]}],2],#!={}&])//(If[verb,X\[Function]X,Quiet]);
ints=Join[ints,newints];

nT=Length[trajs];

Do[
traj1=ints[[k,1,1]];
traj2=ints[[k,2,1]];
\[CapitalOmega]1=\[CapitalOmega][traj1];
\[CapitalOmega]2=\[CapitalOmega][traj2];

\[Gamma]1=DecomposePrs[uPt[ints[[k,1]]//(X\[Function]{X[[1]],X[[2]]-1})],\[Eta]1[ints[[k,1]]//(X\[Function]{X[[1]],X[[2]]-1})]];
(*\[Gamma]2=DecomposePrs[uPt[ints[[k,2]]//(X\[Function]{X[[1]],X[[2]]-1})],\[Eta]1[ints[[k,2]]//(X\[Function]{X[[1]],X[[2]]-1})]];*)
\[Gamma]2=DecomposePrs[uPt[ints[[k,1]]//(X\[Function]{X[[1]],X[[2]]-1})],\[Eta]1[ints[[k,2]]//(X\[Function]{X[[1]],X[[2]]-1})]];
m=Abs[\[Gamma]1.{{0,-1},{1,0}}.\[Gamma]2];
ints[[k,4]]=m;

If[verb,Print["New intersections: ",ints[[k]]," with m=",m]];

newcharges=Select[KSWCF[m,MaxRank,\[CapitalOmega]1,\[CapitalOmega]2],#[[1]]!={0,1}&&#[[1]]!={1,0}&];
(*Print[newcharges];*)
genealogy=Join[genealogy,Table[{newints[[k-nInt,1,1]],newints[[k-nInt,2,1]]},{i,1,Length[newcharges]}]];

Do[
p=newcharges[[ch,1,1]];
q=newcharges[[ch,1,2]];
new\[CapitalOmega]=newcharges[[ch,2]];

\[Omega]1=p \[Eta]1[ints[[k,1]]]+ q \[Eta]1[ints[[k,2]]];
\[Omega]1prime=p \[Eta]1[ints[[k,1]]//(X\[Function]{X[[1]],X[[2]]-1})]+ q \[Eta]1[ints[[k,2]]//(X\[Function]{X[[1]],X[[2]]-1})];

\[CapitalDelta]\[Omega]1=(\[Omega]1prime-\[Omega]1)/(trajs[[ints[[k,1,1]],1,ints[[k,1,2]]]]-trajs[[ints[[k,1,1]],1,(ints[[k,1,2]]-1)]]);

u0=ints[[k,3]];

newtraj=Append[newtraj,Join[PicardFuchs[\[Omega]1,\[CapitalDelta]\[Omega]1,\[Theta],u0,\[CapitalLambda]],{p TrajCharge[traj1]+ q  TrajCharge[traj2],new\[CapitalOmega]}]]//If[verb,(X\[Function]X),Quiet];

,{ch,1,Length[newcharges]}];

,{k,1+nInt,Length[ints]}];

nInt=Length[ints];
trajs=Join[trajs,newtraj];
paths=(X\[Function]X[[1]])/@trajs;
]

BuildKWalls[n_]:=Module[{},
initialize//If[!verb,Quiet];
If[verb,Print["Created primary walls"]];
Do[iteration,{i,1,n}];
];


(* ::Subsection:: *)
(*Scanning \[Theta]*)


ScanTheta[steps_,filename_,gen_]:=Module[{string},
MSpts={};
EvolutionData={};
string[j_]:=filename<>"_"<>ToString[j]<>".png";

Do[
Print["Preparing frame ",i," of ",steps];
\[Theta]=N[(i-1)\[Pi]/steps];
BuildKWalls[gen];
If[i==1,AllCentralCharges=Array[{}&,NPrimKWalls]];
Do[
AllCentralCharges[[j]]=Join[AllCentralCharges[[j]],CentralCharge[j],{trajs[[j,1,NPrimKWalls]]}]
,{j,1,NPrimKWalls}];
EvolutionData=Append[EvolutionData,trajs];
Export[string[i],PlotKWalls];
MSpts=Join[MSpts,(X\[Function]X[[3]])/@ints];
,{i,1,steps}];

AllCentralCharges=Cleanup/@AllCentralCharges;

MakeVideo:=Module[{frames={},vstring},
Print["Preparing video"];
frames=Table[Import[string[i]],{i,1,steps}];
vstring=filename<>"video.flv";
Export[vstring,frames,"FrameRate"->20];
Print["Video exported to "<>vstring]
];

];


CentralCharge[kwall_]:=Module[{temp={}},
Do[temp=Append[temp,{trajs[[kwall,1,i]],(trajs[[kwall,1,i+1]]-trajs[[kwall,1,i]])(trajs[[kwall,2,i]]+trajs[[kwall,2,i+1]])/2+If[i>1,temp[[i-1,2]],0]}];
,{i,1,Length[trajs[[kwall,1]]]-1}];
temp
];

Cleanup[list_]:=Select[list,#!={}&];

PlotReZ[kwallN_]:=Module[{options={PlotLabel->Subscript["Re Z", (2 kwallN-1)],AxesLabel->{"Re u","Im u"},ImageSize->500,PlotRange->Join[pltrng,{All}]}},
{ListPointPlot3D[(X\[Function]{X[[1]]//Re,X[[1]]//Im,X[[2]]//Re})/@AllCentralCharges[[2kwallN-1]],options],ListPointPlot3D[(X\[Function]{X[[1]]//Re,X[[1]]//Im,X[[2]]//Re})/@AllCentralCharges[[2kwallN]],options],
ListPlot3D[(X\[Function]{X[[1]]//Re,X[[1]]//Im,X[[2]]//Re})/@AllCentralCharges[[2kwallN-1]],options],ListPlot3D[(X\[Function]{X[[1]]//Re,X[[1]]//Im,X[[2]]//Re})/@AllCentralCharges[[2kwallN]],options]}//Show
];

PlotImZ[kwallN_]:=Module[{options={PlotLabel->Subscript["Im Z", (2 kwallN-1)],AxesLabel->{"Re u","Im u"},ImageSize->500,PlotRange->Join[pltrng,{All}]}},
{ListPointPlot3D[(X\[Function]{X[[1]]//Re,X[[1]]//Im,X[[2]]//Im})/@AllCentralCharges[[2kwallN-1]],options],ListPointPlot3D[(X\[Function]{X[[1]]//Re,X[[1]]//Im,X[[2]]//Im})/@AllCentralCharges[[2kwallN]],options],
ListPlot3D[(X\[Function]{X[[1]]//Re,X[[1]]//Im,X[[2]]//Im})/@AllCentralCharges[[2kwallN-1]],options],ListPlot3D[(X\[Function]{X[[1]]//Re,X[[1]]//Im,X[[2]]//Im})/@AllCentralCharges[[2kwallN]],options]}//Show
];

PlotCentralCharge[kwallN_]:=Column[{PlotReZ[kwallN],PlotImZ[kwallN]}]


PlotMSWalls[folder_]:=Module[{graph},
graph:=Show[Join[{ListPlot[(X\[Function]{X//Re,X//Im})/@MSpts,PlotRange->pltrng,PlotMarkers->{Automatic,5},AspectRatio->1]},Table[Graphics[{Red,PointSize[Large],Point[{singpts[[i]]//Re,singpts[[i]]//Im}]}],{i,1,Length[singpts]}]],Axes->False,PlotRange->(0.65 pltrng)];
Export[folder<>"MSwalls.pdf",graph];
graph
];


(* ::Subsection::Closed:: *)
(*Extracting BPS degeneracies*)


BPSspectrum[u_,folder_,range_]:=Module[{points,picks,spec={},graph},
Do[
graph=Join[(Z\[Function]ListLinePlot[Z,PlotRange->pltrng,Axes->False,ImageSize->500,AspectRatio->1])/@(X\[Function](Y\[Function]{Y//Re,Y//Im})/@X)/@((Z\[Function]Z[[1]])/@EvolutionData[[frame]]),
Table[Graphics[{Red,PointSize[Large],Point[{singpts[[i]]//Re,singpts[[i]]//Im}]}],{i,1,Length[singpts]}],
{Graphics[{Blue,Point[{u//Re,u//Im}]}],
Graphics[{Blue,Circle[{u//Re,u//Im},range]}]}]//Show;
Export[folder<>"spectrum-at-u-"<>ToString[frame]<>".png",graph];
Do[
points=EvolutionData[[frame,tr,1]];
picks=Select[points,Abs[#-u]<=range&];
If[Length[picks]!=0,spec=Append[spec,{EvolutionData[[frame,tr,3]],EvolutionData[[frame,tr,4]]}]]
,{tr,1,Length[EvolutionData[[frame]]]}];
,{frame,1,Length[EvolutionData]}];
Print["total counts ",spec//Length];
spec//DeleteDuplicates
];


(* ::Subsection:: *)
(*Closing*)


End[ ]

EndPackage[ ]


