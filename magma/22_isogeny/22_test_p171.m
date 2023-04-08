/*
In this example we compute a  (2^87,2^87)-isogeny starting at a predefined abelian surface.
This abelian surface is the Jacobian of a hyperelliptic curve defined by a Type-2 equation.
The (2^87,2^87)-group defining this isogeny is randomly selected.
 */
clear;
load "strategy.m";
load "22_symplectic_basis.m";
load "22_richelot_isogeny_chain.m";

/* typical G2SIDH parameters  with a 171-bit prime (expected 128-bit security)*/
p := 2^87*3^49*5^3-1;
Fp2<z2> := FiniteField(p^2);

type2_invariants :=  [
3446798252521358094517040152788188703411472634458348*z2 +
3934251785801114798773063613093606362278332002469273,
14442707041689001196121281746286164295977084155961*z2 +
2109612034623761207119749551680295829371348991735167,
812627276298732058211077115034347230811721036500113*z2 +
772992546050249836417619075409507366069466860817407,
4149184472552449328486684911846181278057937206823244*z2 +
76217450280973533167781361110947786357072578100390
];

A,B,C,E := Explode(type2_invariants);

R<x> := PolynomialRing(Fp2);
Y := HyperellipticCurve((x^2-1)*(x^2-A)*(E*x^2-B*x+C));
J := Jacobian(Y);

/*
A special symplectic basis was precomputed to speed-up the computation.
Alternatively it can be generated using the functions in "setup.m".
*/

BB := [ J![x^2 + (4228467102363626818406821291187317314387242721137631*z2 +
3804540158078856260408798951914626839615188809113971)*x +
4275797638838160436292168115677773824418883521347906*z2 +
1137321469212884956492608367539267130850696607667781,
(4610281969133590881952445078108252720655017745022532*z2 +
3973686199211383520298404807386383773962007438550548)*x +
648258542931357692710711566778124079589290457474914*z2 +
215955646752530011869381596473856112953124938376677],
J! [x^2 +
(4338911201117717228685979630682819315648743416178581*z2 +
164719356863546810047096688631436661008130416921123)*x +
682832746418814828949968317152011216794405052776050*z2 +
3227347167670773728077792441891694104440549509869116,
(3338502728879829955316606942134401994124769140969765*z2 +
4592275636898486852119197503541440469250199985025556)*x +
2757504835057898567366283741674487382402327863649680*z2 +
3009521334788140091333540237394896005087239208970452],
J![x^2 +
(2555738952316265353748094437339089294680814610537684*z2 +
2176640765922335894797016044059634690560411512928873)*x +
2458179105692663378496708693061286763272601074742808*z2 +
1624818658770637954475518067668700878207038193447942,
(4088437804309543062900709175013586494723123742583908*z2 +
2807979495610579340461439943471276151550973301745513)*x +
955764076954805838607110282177328845041233230865666*z2 +
564011195494660755723258516073314725149407100733060],
J![x^2 +
(3764355748013851996734964908004966244986492767609594*z2 +
4552431470848504784464885569352125742862744106050259)*x +
3689157838657674380945494042451165029342922168026393*z2 +
2127491357759698556236561733691614246364243090239406,
(916752860733121206178059110704899351395203541469843*z2 +
1183439552501690970580496543748515900250690676492295)*x +
128074129624897914327260173767029947581651278749621*z2 +
1819697913691524214918490466769661929366215392418827] ];

/*
In order to simulate the first round of a G2SIDH key exchange, where Alice needs
to compute the images of the 3^m-torsion, we generate 4 random points on the Jacobian.
*/

PP := [];
for i := 1 to 4 do
P := Random(J);
a0 := Coefficient(P[1],0);
a1 := Coefficient(P[1],1);
a2 := Coefficient(P[1],2);
b0 := Coefficient(P[2],0);
b1 := Coefficient(P[2],1);
PP := Append(PP, [a0,a1,a2,b0,b1]);
end for;

T1 := 0.0;
T2 := 0.0;
T3 := 0.0;
T4 := 0.0;
for iters:=1 to 25 do
ker := RandomSymplecticGroup(BB,87);

printf "ker";
ker;
t1 := Cputime();
Ynew1 := RichelotChain(type2_invariants, ker, 87: P_list:=PP, last_step:=true);
T1 +:= Cputime(t1);
printf "CPU time for Richelot chain (Round 1): %o\n", Cputime(t1);
t2 := Cputime();
Ynew2 := RichelotChain(type2_invariants, ker, 87: last_step:=true);
T2 +:= Cputime(t2);
printf "CPU time for Richelot chain (Round 2): %o\n", Cputime(t2);

//strategy := [87 - k : k in [1 .. 87]];
st := balanced_strategy(87);

t3 := Cputime();
Ynew3 := RichelotChain22_strategy(type2_invariants, ker, 87 : P_list:=PP, S:=st);
T3 +:= Cputime(t3);
printf "CPU time for Richelot chain using strategies (Round 1): %o\n", Cputime(t3);

t4 := Cputime();
Ynew4 := RichelotChain22_strategy(type2_invariants, ker, 87 : P_list:=[], S:=st);
T4 +:= Cputime(t4);
printf "CPU time for Richelot chain using strategies (Round 2): %o\n", Cputime(t4);


H1 := getCurve(Ynew1, x);
assert Jacobian(H1)!0 eq (p+1) * Random(Jacobian(H1));
H2 := getCurve(Ynew2, x);
assert Jacobian(H2)!0 eq (p+1) * Random(Jacobian(H2));
H3 := getCurve(Ynew3, x);
assert Jacobian(H3)!0 eq (p+1) * Random(Jacobian(H3));
H4 := getCurve(Ynew4, x);
assert Jacobian(H4)!0 eq (p+1) * Random(Jacobian(H4));


// Important fact: all the curves coincide!
assert IsIsomorphic(H1, H2);
assert IsIsomorphic(H1, H3);
assert IsIsomorphic(H1, H4);
assert IsIsomorphic(H2, H4);
assert IsIsomorphic(H3, H4);
assert IgusaInvariants(H1) eq IgusaInvariants(H2);
assert IgusaInvariants(H1) eq IgusaInvariants(H3);
assert IgusaInvariants(H1) eq IgusaInvariants(H4);
assert IgusaInvariants(H2) eq IgusaInvariants(H4);
assert IgusaInvariants(H3) eq IgusaInvariants(H4);
end for;

printf "Average (Round 1): %o/%o = %o;\n", T1, T2, T1 / T2;
printf "Average (Round 2): %o/%o = %o;\n", T3, T4, T3 / T4;

exit;