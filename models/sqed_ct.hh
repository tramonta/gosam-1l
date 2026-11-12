* vim: ts=3:sw=3
*---#[ Symbol Definitions :
Symbols [field.em], [field.ep],
	[field.mum], [field.mup],
	[field.taum], [field.taup];
Symbols [field.A],
	[field.phip], [field.phim];

Symbols [field.Cx]; 

Symbols mW, me, mmu, mtau;
Symbols wtau;
Symbols e;
Symbols gAAPP;
Symbols gAPP;
Symbols gctAAPP;
Symbols gctAPP;
Symbols gctA, gctphi, gctphi1;

#Do f={e,mu,tau}
Symbols gctAl`f', gctAr`f';
Symbols gctCL`f', gctCR`f', gctCP`f', gctCM`f';
#EndDo

*---#] Symbol Definitions :
*---#[ Procedure VertexConstants :
#Procedure VertexConstants
* All equation numbers refer to Boehm, Denner, Joos: Gauge Theories
*
*---#[   flavour quantum numbers of the fermions :
*----#[   Leptons :
#define Qe    "(-1)"
#define Qmu   "(-1)"
#define Qtau  "(-1)"
*----#]   Leptons :
*---#]   flavour quantum numbers of the fermions :

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%                                             %
*%                   EW-CT                     %
*%                                             %
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

* VV
Id C1([field.A], [field.A], [field.Cx]) = PREFACTOR(gctA);
Id C2([field.A], [field.A], [field.Cx]) = 0;

* SS
Id C1([field.phip], [field.phim], [field.Cx]) = PREFACTOR(gctphi1);
Id C2([field.phip], [field.phim], [field.Cx]) = PREFACTOR(gctphi);

* FF
#Do f={e,mu,tau}
Id CL([field.`f'p], [field.`f'm], [field.Cx]) = PREFACTOR(gctCL`f');
Id CR([field.`f'p], [field.`f'm], [field.Cx]) = PREFACTOR(gctCR`f');
Id CP([field.`f'p], [field.`f'm], [field.Cx]) = PREFACTOR(gctCP`f');
Id CM([field.`f'p], [field.`f'm], [field.Cx]) = PREFACTOR(gctCM`f');
#EndDo

*---#[ VVSS (A.2.14) :
Id C([field.A], [field.A], [field.phip], [field.phim], [field.Cx]) = PREFACTOR(e^2 * gctAAPP);
*---#] VVSS (A.2.14) :
*---#[ VSS (A.2.17) :
Id C([field.A], [field.phim], [field.phip], [field.Cx]) = PREFACTOR(e * gctAPP);
*---#] VSS (A.2.17) :

*---#[ VFF (A.2.21) :
*---#[   AFF :
#Do f={e,mu,tau}
	Id CL([field.`f'p], [field.`f'm], [field.A], [field.Cx]) = PREFACTOR( e * gctAl`f');
	Id CR([field.`f'p], [field.`f'm], [field.A], [field.Cx]) = PREFACTOR( e * gctAr`f');
#EndDo
*---#]   AFF :
*---#] VFF (A.2.21) :

*---#]   EW quantum numbers of the fermions :
*---#[ VVSS (A.2.14) :
Id C([field.A], [field.A], [field.phip], [field.phim]) = PREFACTOR(e^2 * gAAPP);
*---#] VVSS (A.2.14) :
*---#[ VSS (A.2.17) :
Id C([field.A], [field.phim], [field.phip]) = PREFACTOR(e * gAPP);
*---#] VSS (A.2.17) :
*---#[ VFF (A.2.21) :
*---#[   AFF :
#Do f={e,mu,tau}
	Id CL([field.`f'p], [field.`f'm], [field.A]) = PREFACTOR(- (`Q`f'') * e);
	Id CR([field.`f'p], [field.`f'm], [field.A]) = PREFACTOR(- (`Q`f'') * e);
#EndDo
*---#]   AFF :
*---#] VFF (A.2.21) :

#EndProcedure
*---#] Procedure VertexConstants :
