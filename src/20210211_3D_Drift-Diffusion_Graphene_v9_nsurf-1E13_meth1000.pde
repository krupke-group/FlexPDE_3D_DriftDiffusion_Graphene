TITLE 'GFET_DD-model_in_charged_layered_dielectric'
{Drift-diffusion in graphene formalism:
Biegel et al., Nasa Technical Report 1998
Selberherr, Analysis and simulation of semiconductor devices, page 19/20
Ancona, IEEE Transactions on Electronic Devices 57, 681 (2010)
Doan et al., Z. Angew. Math. Phys. 70, 55 (2019)
Champlian, J. Appl. Phys. 109, 084515 (2011)
Fang et la., Appl. Phys. Lett. 91, 092109 (2007)
Farrell et al., Challenges in Drift-Diffusion Semiconductor Simulations: A COMPARATIVE STUDY OF DIFFERENT DISCRETIZATION PHILOSOPHIES 2018
---
u:   		div(eps*grad(u))=0
phin:	-div(q*n*mun*grad(phin))=0
phip:	-div(q*p*mup*grad(phip))=0
---
n(phin,u)=NG*F1((q*(u-phin)-EDirac)/(k*T))
p(phip,u)=NG*F1((q*(phip-u)+EDirac)/(k*T))
}

COORDINATES cartesian3
VARIABLES
	u			{potential}
	phin		{quasi-Fermi level / chemical potential electrons}
	phip		{quasi-Fermi level / chemical potential holes}
!	jn			{integrated electron current}
!	jp			{integrated hole current} 

SELECT
	!autostage=off
	!ngrid = 50
	errlim=1e-2
	changelim=1e-5
	NGRID=20
	REGRID=OFF
	THREADS = 8
	SMOOTHINIT=ON
	SPECTRAL_COLORS = ON
	CONTOURS = 1e6
	PAINTED

DEFINITIONS
{ Fermi-Dirac-Integral: Aymerich-Humet et al., J. Appl. Phys. 54, 2850 (1983) }
i=1
a = ( 1 + 15/4*(i+1) + 1/40*(i+1)^2 )^0.5
b = 1.8 + 0.61*i
c = 2 + (2-sqrt(2))*2^(-i)
aa = (i+1)*2^(i+1)
bb(eta) = ((ABS(eta-b))^c + a^c)^(1/c)
F1(eta) = ( aa/(b+eta+bb(eta))^(i+1) + (EXP(-eta))/GAMMAF(i+1) )^(-1)

{fundamental constants and dimensions}
nm = 1e-9						! meters per nanometer
q = 1.602E-19					! [C]
h = 6.626E-34 					! [m² kg /s ]
k = 1.38064852e-23 		! [m² kg /s² /K]
hbar=h/(2*pi)
edge=25*nm
Lx=25*nm+edge	
eps0=8.854e-12				! [F/m]
epsr eps=epsr*eps0
!Ex=-dx(u) Ey=-dy(u) Ez=-dz(u) E=-grad(u) Em=magnitude(E)
!D=eps*E Dm=magnitude(D)
!delta(z)=(sign(z+0.1*nm)+1)/2*(sign(-z+0.1*nm)+1)/2

{graphene properties} 
T = 300							! [K ]
vF = 1E6							! [m/s]
Eg = 0								! no bandgap
mun mup
mung = 1000*1e-4			! [cm²/Vs] converted in [m²/Vs]
mupg = 1000*1e-4			! [cm²/Vs] converted in [m²/Vs]
EDirac = 0						! [eV]

{graphene carrier concentration: Champlian, J. Appl. Phys. 109, 084515 (2011)}
ni = (pi/6)*((k*T)/(hbar*vF))^2
NG=(2/pi)*((k*T)/(hbar*vF))^2
!n(phin,u)=NG*F1((q*(u-phin)-EDirac)/(k*T))
!p(phip,u)=NG*F1((q*(phip-u)+EDirac)/(k*T))
n(phin,u)=NG*F1((q*(u-phin)-EDirac)/(k*T))
p(phip,u)=NG*F1((q*(phip-u)+EDirac)/(k*T))
!n(phin,u)=NG*F1((phin-EDirac+q*u)/(k*T))
!p(phip,u)=NG*F1((-(phip-EDirac)-u)/(k*T))

{conditions}
VL=0.005		phinL=VL		phipL=VL		! potential and chemical potential left
VR=-0.005		phinR=VR	phipR=VR	! potential and chemical potential right
!VG=20.0													! back-gate potential
!VG=STAGED(-2,-1,0,1,2)
!VG=STAGED(-20,-10,-5,-2.5,-1,-0.5,0,0.5,1,2.5,5,10,20)
VGmax=20
VGmin=5
deltaVG=0.25
VG=STAGED(VGmax BY -deltaVG TO VGmin)
!n_surface = 0											! charge on alumina surface
n_surface = -q*1e13*1e4		! [/m²]

{currents}
jn=-q*mun*n(phin,u)*grad(phin)
jp=-q*mup*p(phip,u)*grad(phip)
jnx=XCOMP(jn)	jpx=XCOMP(jp)
jny=YCOMP(jn)	jpy=YCOMP(jp)
jnz=ZCOMP(jn)	jpz=ZCOMP(jp)

INITIAL VALUES
	u = 0
	phin = 0
	phip = 0

EQUATIONS
	u:			div(eps*grad(u))=0
	{2D charges introduced by boundary conditions at the global surface 'Charges" and LIMITED REGION "center" SURFACE	"graphene" }
	phin:	-div(q*mun*n(phin,u)*grad(phin))=0
	phip:	-div(q*mup*p(phip,u)*grad(phip))=0
!THEN
!	jn:		dy(jn)=q*mun*n(phin,u)*dx(phin)
!	jp:		dy(jp)=-q*mup*p(phip,u)*dx(phip)

! CONSTRAINTS    { Integral constraints }

EXTRUSION
	SURFACE	"Bottom"										Z=-30*nm
	LAYER		"SiO2"
	SURFACE	"graphene_bottom"						Z=-0.5*nm
	LAYER		"graphene"							
 	SURFACE	"graphene_top"								Z=+0.5*nm
	LAYER		"Al2O3"
	SURFACE	"Charges"										Z=+2.5*nm
	LAYER		"Methanol"
	SURFACE	"Methanol - Air"							Z=+10*nm
	LAYER		"Air"	
	SURFACE	"Top"											Z=+30*nm

BOUNDARIES
	SURFACE	"Top"							
	!SURFACE	"Methanol - Air"			
	SURFACE	"Charges"						natural(u)=n_surface
	SURFACE	"graphene_top"	
	!SURFACE	"graphene_bottom"		
	SURFACE	"Bottom"						value(u)=VG

	REGION "global"
		INACTIVE (phin,phip)
		LAYER	'Air'								epsr=1 mun=0 mup=0
		LAYER	'Methanol'						epsr=1000 mun=0 mup=0
		LAYER	'Al2O3'							epsr=8 mun=0 mup=0
		LAYER 'graphene'						epsr=4 mun=mung mup=mupg
		LAYER	'SiO2'							epsr=3.9 mun=0 mup=0
			START (-1.1*Lx,-1.1*Lx)
				LINE TO (1.1*Lx,-1.1*Lx) LINE TO (1.1*Lx,1.1*Lx) LINE TO (-1.1*Lx,1.1*Lx)  LINE TO CLOSE

	LIMITED REGION "left"
		!INACTIVE (phin,phip)
		SURFACE	"graphene_top"				value(u)=VL value(phin)=phinL value(phip)=phipL
		SURFACE	"graphene_bottom"		value(u)=VL value(phin)=phinL value(phip)=phipL
	 	LAYER	'graphene'					epsr=4 mun=mung mup=mupg
			START (-Lx, -Lx+edge)					value(u)=VL value(phin)=phinL value(phip)=phipL		!front(y<<0)
			LINE TO (-Lx+edge, -Lx+edge)	value(u)=VL	 value(phin)=phinL value(phip)=phipL		!right(>>0)
			LINE TO (-Lx+edge, Lx-edge)		value(u)=VL	 value(phin)=phinL value(phip)=phipL		!back(y>>0)
			LINE TO (-Lx, Lx-edge)					value(u)=VL	 value(phin)=phinL value(phip)=phipL		!left(x<<0)
			LINE TO CLOSE

	LIMITED REGION "center" 
		!INACTIVE (phin,phip)
		SURFACE	"graphene_top"					natural(u)=q*(p(phip,u)-n(phin,u))/2
		SURFACE	"graphene_bottom"			natural(u)=q*(p(phip,u)-n(phin,u))/2
	 	LAYER	'graphene'					epsr=4 mun=mung mup=mupg
			START (-Lx+edge, -Lx+edge)			natural(u)=0 natural(phin)=0 natural(phip)=0			!front
			LINE TO (Lx-edge, -Lx+edge) 			value(u)=VR value(phin)=phinR value(phip)=phipR	!right
			LINE TO (Lx-edge, Lx-edge)				natural(u)=0 natural(phin)=0 natural(phip)=0			!back
			LINE TO (-Lx+edge, Lx-edge)			value(u)=VL	 value(phin)=phinL value(phip)=phipL	!left
			LINE TO CLOSE

	LIMITED REGION "right"
		!INACTIVE (phin,phip)
		SURFACE	"graphene_top"					value(u)=VR value(phin)=phinR value(phip)=phipR
		SURFACE	"graphene_bottom"			value(u)=VR value(phin)=phinR value(phip)=phipR
	 	LAYER	'graphene'					epsr=4 mun=mung mup=mupg
			START (Lx-edge, -Lx+edge)			value(u)=VR value(phin)=phinR value(phip)=phipR		!front(y<<0)
			LINE TO (Lx, -Lx+edge)				value(u)=VR value(phin)=phinR value(phip)=phipR		!right(x>>0)
			LINE TO (Lx, Lx-edge)					value(u)=VR value(phin)=phinR value(phip)=phipR		!back(y>>0)
			LINE TO (Lx-edge, Lx-edge)			value(u)=VR value(phin)=phinR value(phip)=phipR		!left(x<<0)
			LINE TO CLOSE

		
! TIME 0 TO 1    { if time dependent }
MONITORS
contour(u) on z=+0.0*nm as "potential u / V"
contour(phin) on z=+0.0*nm as 'phin / eV'
contour(phip) on z=+0.0*nm as 'phip / eV'
contour((n(phin,u))*1e-4) on z=+0.0*nm as 'carrier concentration n /cm^2'
contour((p(phip,u))*1e-4) on z=+0.0*nm as 'carrier concentration p /cm^2'
contour((p(phip,u)-n(phin,u))*1e-4) on z=+0.0*nm as 'carrier concentration p-n /cm^2'


PLOTS            { save result displays }
	grid(x,y) on z=0
	grid(y,z) on x=0  zoom (-Lx,-20*nm,2*Lx,40*nm)
	grid(x,z) on y=0  zoom (-Lx,-20*nm,2*Lx,40*nm)
	contour(u) on z=0 as "potential u / V"
	contour(u) on x=0
	contour(u) on y=0
	elevation(u) from (-Lx,0,0) to (Lx,0,0)
	elevation(u) from (0,0,-300*nm) to (0,0,300*nm) report val (u,0,0,0)
	contour(phin) on z=+0.0*nm as 'phin / eV'
	contour(phip) on z=+0.0*nm as 'phip/ eV'
	elevation(phin,phip) from (-Lx,0,0) to (Lx,0,0)
	elevation(phin,phip) from (0,0,-300*nm) to (0,0,300*nm) report val (u,0,0,0)
	contour((p(phip,u))*1e-4) log on z=0 as 'carrier concentration p /cm^2'
	contour((p(phip,u))*1e-4) log on x=0 zoom (-Lx,-20*nm,2*Lx,40*nm) as 'carrier concentration p /cm^2'
	contour((p(phip,u))*1e-4) log on y=0 zoom (-Lx,-20*nm,2*Lx,40*nm) as 'carrier concentration p /cm^2'
	contour((n(phin,u))*1e-4) log on z=0 as 'carrier concentration n /cm^2'
	contour((n(phin,u))*1e-4) log on x=0 zoom (-Lx,-20*nm,2*Lx,40*nm) as 'carrier concentration n /cm^2'
	contour((n(phin,u))*1e-4) log on y=0 zoom (-Lx,-20*nm,2*Lx,40*nm) as 'carrier concentration n /cm^2'
	contour((p(phip,u)-n(phin,u))*1e-4) on z=0 as 'carrier concentration p-n /cm^2'
	contour((p(phip,u)-n(phin,u))*1e-4) on x=0 zoom (-Lx,-20*nm,2*Lx,40*nm) as 'carrier concentration p-n /cm^2'
	contour((p(phip,u)-n(phin,u))*1e-4) on y=0 zoom (-Lx,-20*nm,2*Lx,40*nm) as 'carrier concentration p-n /cm^2'
		report F1(0) report ni*1e-4 as 'ni/cm^2' report NG*1e-4 as 'NG/cm^2'
	!report delta(0) report delta(0.11*nm)
	!contour(Em) painted contours = 10000 on z=0
	!contour(Em) painted contours = 10000 on x=0 zoom (-Lx,-10,2*Lx,20)
	!contour(Em) painted contours = 10000 on y=0 zoom (-Lx,-10,2*Lx,20)
	elevation((p(phip,u)-n(phin,u))*1e-4) from (-Lx,0,0*nm) to (Lx,0,0*nm) as 'carrier concentration p-n /cm^2'
		report val (p(phip,u)*1e-4,0,0,0) report val (n(phin,u)*1e-4,0,0,0)
	elevation((p(phip,u)-n(phin,u))*1e-4) from (0,0,-5*nm) to (0,0,5*nm) as 'carrier concentration p-n /cm^2'
		report val((p(phip,u)-n(phin,u))*1e-4,0,0,0) as 'center carrier concent. p-n /cm^2' !report val((p+n)*1e-4,0,0,0) as 'center carrier concentration p+n /cm^2'
		report SURF_INTEGRAL((p(phip,u)-n(phin,u))*1e-4,"graphene_top","center","graphene")/SURF_INTEGRAL(1,"graphene_top","center","graphene") as 'mean center carrier concent. p-n /cm^2'
	elevation((p(phip,u))*1e-4,(n(phin,u))*1e-4) log from (-Lx,0,0*nm) to (Lx,0,0*nm) as 'carrier concentration /cm^2'
	elevation((p(phip,u))*1e-4,(n(phin,u))*1e-4) log from (0,0,-5*nm) to (0,0,5*nm) as 'carrier concentration /cm^2'
	elevation((p(phip,u)+n(phin,u))*1e-4) log from (-Lx,0,0*nm) to (Lx,0,0*nm) as 'carrier concentration p+n /cm^2'
	elevation(jnx,jpx,jnx+jpx) from (-Lx,0,0*nm) to (Lx,0,0*nm) as 'current A/m'
	
	contour(magnitude(jn)) on z=0 as 'jn'
	contour(magnitude(jp)) on z=0 as 'jp'
	contour(magnitude(jn+jp)) on z=0 as 'jn+jp'
	contour(jnx/magnitude(jn)) on z=0 as 'jnx/jn'
	contour(jny/magnitude(jn)) on z=0 as 'jny/jn'
	contour(jnz/magnitude(jn)) on z=0 as 'jnz/jn'
	contour(jnx+jpx) on x=0 zoom (-Lx,-20*nm,2*Lx,40*nm) as 'jnx+jpx'
	contour(jnx+jpx) on y=0 zoom (-Lx,-20*nm,2*Lx,40*nm) as 'jnx+jpx'


histories 
	history(val(p(phip,u)*1e-4,0,0,0), val(p(phip,u)*1e-4,0,0.9*(Lx-edge),0), val(n(phin,u)*1e-4,0,0,0), val(n(phin,u)*1e-4,0,0.9*(Lx-edge),0)) log 
			vs VG as "carrier concentration p,n /cm^2 in graphene vs VG / V" report (n_surface*1e-4) as "n_surface / cm^2" !EXPORT FILE="p,n.txt"
	history(val((p(phip,u)+n(phin,u))*1e-4,0,0,0)) log {val((p(phip,u)+n(phin,u))*1e-4,0,0.9*(Lx-edge),0)}
			vs VG as "carrier concentration p+n /cm^2 in graphene vs VG / V" report (n_surface*1e-4) as "n_surface / cm^2" !EXPORT FILE="p+n.txt"
	history(val(jnx+jpx,0,0,0), val(jnx+jpx,0,0.9*(Lx-edge),0))
			vs VG as "current density jnx+jpx" report (n_surface*1e-4) as "n_surface / cm^2" !EXPORT FILE="jnx+jpx.txt"

{
	history(val(p(phip,u)*1e-4,0,0,0), val(n(phin,u)*1e-4,0,0,0)) log vs VG as "carrier concentration p,n /cm^2 in graphene vs VG / V" report (n_surface*1e-4) as "n_surface / cm^2"
	history(val((p(phip,u)+n(phin,u))*1e-4,0,0,0)) log vs VG as "carrier concentration p+n /cm^2 in graphene vs VG / V" report (n_surface*1e-4) as "n_surface / cm^2"
	history(val(jnx+jpx,0,0,0)) vs VG as "current density jnx+jpx" report (n_surface*1e-4) as "n_surface / cm^2"

	history(val(p(phip,u)*1e-4,0,0.9*(Lx-edge),0), val(n(phin,u)*1e-4,0,0.9*(Lx-edge),0)) log vs VG as "carrier concentration p,n /cm^2 in graphene vs VG / V" report (n_surface*1e-4) as "n_surface / cm^2"
	history(val((p(phip,u)+n(phin,u))*1e-4,0,0.9*(Lx-edge),0)) log vs VG as "carrier concentration p+n /cm^2 in graphene vs VG / V" report (n_surface*1e-4) as "n_surface / cm^2"
	history(val(jnx+jpx,0,0.9*(Lx-edge),0)) vs VG as "current density jnx+jpx" report (n_surface*1e-4) as "n_surface / cm^2"
}
END