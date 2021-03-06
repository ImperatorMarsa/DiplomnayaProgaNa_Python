#!/usr/bin/env python3
from PIL import ImageGrab, ImageDraw, ImageFont
from vpython import *
import random, os, time, copy

scene=canvas(width=748, height=748)
MassivDlyaPereschetaPeriudGranic=[[0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1], [0,0,0]]
#####__Fundomental_Postoyannie__#################################################################################################################################################################################################################################################
Vyazkost=2.15e-3 # Paskal*sekunda

Time=1e-10 # sekund
kT=273.16*1.38e-23#1.38e-23
u0=4e-7*pi # Genri/metr

graniciVselennoy=[7.3e-8]*3
graniciVidemogo=[7.3e-8]*3

rate(1/Time)
#################################################################################################################################################################################################################################################################################

#####__Chastica__################################################################################################################################################################################################################################################################
class Chastica(object):
    """ento chastica magnitnoy jidkosti"""

    Radiuse=6.66e-9  # metrov
    Plotnost=5000 # kilogramm/metr^3
    Obyom=4/3*pi*Radiuse**3 # metrov^3
    Massa=Obyom*Plotnost # kilogramm
    MagMom=4.78e5*Obyom # Amper*metr^2 ((namagnichenost' nasisheniya=4.78*10^5 Amper/metr))
    print("Dipol'niy moment chastici ", MagMom)
    def __init__(self, coord, ugol):
        self.Otobrajenie=0

        self.color=vector(1,.5,0)

        self.radius=Chastica.Radiuse
        self.mas=Chastica.Massa

        self.pos=coord
        self.skor=vector(0, 0, 0)
        self.forse=vector(0, 0, 0)

        self.axis=ugol*(Chastica.MagMom)
        self.omega=vector(0, 0, 0)
        self.moment=vector(0, 0, 0)

        self.sfera=sphere(
            radius=self.radius,
            color=self.color,
            make_trail=False,
            trail_color=vector(0, 1, 0),
            opacity=1,
            pos=self.pos)
        self.strela=arrow(
            pos=self.sfera.pos,
            color=vector(1, 1, 0),
            opacity=1,
            axis=norm(self.axis)*(Chastica.Radiuse*3/2))

    #@!@!@!@!@!__Osnovnaya_Proga__!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
    def Kinematika(self): # https://www.desmos.com/calculator/bhjmf8p0pf
        #&$&$&$&$&$&__Algoritm_Algoritm_Varle__&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
        pom=self.StahostSmesh()
        self.skor=self.pos
        self.pos=self.pos+self.forse/(6.0*pi*self.radius*Vyazkost)*Time+pom[0]
        self.skor=(self.pos-self.skor)/Time

        self.omega=self.axis
        DeltaAlfa=self.moment/(8.0*pi*self.radius**3*Vyazkost)*Time+pom[1]
        self.axis=self.axis.rotate(angle=DeltaAlfa.mag, axis=DeltaAlfa)
        self.omega=(self.axis-self.omega)/Time

        self.PorvrkaGrani()
        self.Otobrajenie+=1
        self.sfera.color=self.color
        self.sfera.pos=self.pos
        self.strela.pos=self.sfera.pos
        self.strela.axis=norm(self.axis)*(Chastica.Radiuse*5/3)

        self.forse=vector(0, 0, 0)
        self.moment=vector(0, 0, 0)
        #&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
    #@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!

    #@!@!@!@!@!__Dopolnitelnie_Progi__!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
    def PorvrkaGrani(self):
        if self.pos.x>graniciVselennoy[0]: self.pos=vector(self.pos.x-2*graniciVselennoy[0], self.pos.y, self.pos.z)
        elif self.pos.x<-graniciVselennoy[0]: self.pos=vector(self.pos.x+2*graniciVselennoy[0], self.pos.y, self.pos.z)
        if self.pos.y>graniciVselennoy[1]: self.pos=vector(self.pos.x, self.pos.y-2*graniciVselennoy[1], self.pos.z)
        elif self.pos.y<-graniciVselennoy[1]: self.pos=vector(self.pos.x, self.pos.y+2*graniciVselennoy[1], self.pos.z)
        if self.pos.z>graniciVselennoy[2]: self.pos=vector(self.pos.x, self.pos.y, self.pos.z-2*graniciVselennoy[2])
        elif self.pos.z<-graniciVselennoy[2]: self.pos=vector(self.pos.x, self.pos.y, self.pos.z+2*graniciVselennoy[2])

    def StahostSmesh(self):
        pom=[]
        difuz=kT/(6.0*pi*self.radius*Vyazkost)
        pom.append(vector.random().norm()*((2*difuz*Time)**.5*random.gauss(0, 1)))

        difuz=kT/(8.0*pi*self.radius**3*Vyazkost)
        pom.append(vector.random().norm()*((2*difuz*Time)**.5*random.gauss(0, 1)))

        return pom
    #@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
#################################################################################################################################################################################################################################################################################

#####__Sili_Deystvuyushie_V_Sisteme__############################################################################################################################################################################################################################################
def MehMomDipolya(y, x): # https://en.wikipedia.org/wiki/Magnetic_dipole
    M1=vector(0,0,0)
    M2=vector(0,0,0)

    r=y.pos-x.pos
    if (abs(r.x)-graniciVselennoy[0]<0 and abs(r.y)-graniciVselennoy[1]<0 and abs(r.z)-graniciVselennoy[2]<0):
        B=u0/(4*pi)*(dot(x.axis, r)*3*r/(r.mag)**5-x.axis/(r.mag)**3)
        M1=cross(y.axis, B)
        r=-r
        B=u0/(4*pi)*(dot(x.axis, r)*3*r/(r.mag)**5-x.axis/(r.mag)**3)
        M2=cross(x.axis, B)

    return M1, M2

def MehSilaDipolya(y, x): # https://en.wikipedia.org/wiki/Magnetic_dipole
    F1=vector(0,0,0)
    F2=vector(0,0,0)
    r=y.pos-x.pos
    if (abs(r.x)-graniciVselennoy[0]<0 and abs(r.y)-graniciVselennoy[1]<0 and abs(r.z)-graniciVselennoy[2]<0):
        F1=3*u0/(4*pi*r.mag**5)*(dot(x.axis, r)*y.axis+dot(y.axis, r)*x.axis+dot(x.axis, y.axis)*r-5*dot(x.axis, r)*dot(y.axis, r)*r/r.mag**2)
        r=-r
        F2=3*u0/(4*pi*r.mag**5)*(dot(x.axis, r)*y.axis+dot(y.axis, r)*x.axis+dot(x.axis, y.axis)*r-5*dot(x.axis, r)*dot(y.axis, r)*r/r.mag**2)

    return F1, F2

def VneshPole(x):
    H=73e3 # Amper/metr
    B=vector(1,0,0)*H*u0
    M=cross(x.axis, B)
    return M

def SteerOttalk(x, y): # https://www.desmos.com/calculator/ddxmffkqrj
    pom=x.pos-y.pos
    dist=mag(pom)
    F_ster=vector(0, 0, 0)

    # Kakieto koefficienti #
    A=31.3 #################
    B=73.0 #################
    # ######################

    M=x.axis.mag # Magnitniy moment chastici
    q=2e-9 # Dlina volosni v metrah
    a=2.0*(x.radius+q) # Diametr chastici s volosney
    if dist<d:
        e=exp(-B*(dist/a-1))
        F_ster=A*3*u0*M**2/(4*pi*a**4)*e*pom

    return F_ster, -F_ster
#################################################################################################################################################################################################################################################################################

#####__Storonnie_Funkcii__#######################################################################################################################################################################################################################################################
def NewCoord(x, A):
    a=vector(copysign(1, -x.pos.x)*(2*graniciVselennoy[0]*A[0]-abs(x.pos.x)), copysign(1, -x.pos.y)*(2*graniciVselennoy[1]*A[1]-abs(x.pos.y)), copysign(1, -x.pos.z)*(2*graniciVselennoy[2]*A[2]-abs(x.pos.z)))
    return a
#################################################################################################################################################################################################################################################################################
#       Ya zdes' risuyu kletku      #
gray=color.gray(0.7)
d=graniciVidemogo[0]
r=0.005
boxbottom=curve(color=gray)
boxbottom.append([vector(-d,-d,-d), vector(-d,-d,d), vector(d,-d,d), vector(d,-d,-d), vector(-d,-d,-d)])
boxtop=curve(color=gray)
boxtop.append([vector(-d,d,-d), vector(-d,d,d), vector(d,d,d), vector(d,d,-d), vector(-d,d,-d)])
vert1=curve(color=gray)
vert2=curve(color=gray)
vert3=curve(color=gray)
vert4=curve(color=gray)
vert1.append([vector(-d,-d,-d), vector(-d,d,-d)])
vert2.append([vector(-d,-d,d), vector(-d,d,d)])
vert3.append([vector(d,-d,d), vector(d,d,d)])
vert4.append([vector(d,-d,-d), vector(d,d,-d)])
#       Vsyo narisoval      #

Chastici=[]
Chastici.append(Chastica(vector(0.1e-9, 0.1e-9, 0.1e-9), vector(1, 0, 0)))
#Chastici.append(Chastica(vector(-graniciVselennoy[0]*1/4, -1e-9, -1e-9), vector(1, 0, 0)))
#Chastici.append(Chastica(vector(graniciVselennoy[0]*1/4, -1e-9, -1e-9), vector(1, 0, 0)))
#Chastici.append(Chastica(vector((6.66+2.2)*2e-9, 0, 0), vector(1, 0, 0)))
#Chastici.append(Chastica(vector(-(6.66+2)*2e-9, 0, 0), vector(-sin(pi/10), cos(pi/10), 0)))
N=42 # Chislo chastic
print("Obyomnaya Koncentraciya ", (N+1)*Chastici[-1].Obyom/graniciVselennoy[0]**3)
for A in range(N):
    coord=vector.random()*(graniciVselennoy[0]-Chastici[-1].Radiuse)
    C=len(Chastici)
    while C>0:
        C=len(Chastici)
        for x in Chastici:
            r=abs(mag(coord-x.pos))
            if r<(x.Radiuse+2.5e-9)*2: coord=vector.random()*graniciVselennoy[0]
            else: C-=1
    ugol=vector.random().norm()
    Chastici.append(Chastica(coord, ugol))#'''

time.sleep(4.2)#input("-__Pognali__-")#
for AA in range(1000000): # while True:
    Hronos=time.time()
    for odna in Chastici:
        #&$&$&$&$&$__Summa_Vseh_Sil__&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
        odna.moment+=VneshPole(odna)
        for viteta in Chastici[Chastici.index(odna):]:
            if odna!=viteta:
                for A in MassivDlyaPereschetaPeriudGranic:
                    prizrac=copy.copy(viteta)
                    prizrac.pos=NewCoord(prizrac, A)
                    pom=SteerOttalk(odna, prizrac)
                    odna.forse+=pom[0]; viteta.forse+=pom[1]
                    pom=MehSilaDipolya(odna, prizrac)
                    odna.forse+=pom[0]; viteta.forse+=pom[1]
                    pom=MehMomDipolya(odna, prizrac)
                    odna.moment+=pom[0]; viteta.moment+=pom[1]

        odna.Kinematika()
        ignayt=73
        #'''
        if AA%ignayt==0:#
            Freame=ImageGrab.grab([10,10,748,748])#grab([0,0,300,300])
            puti='D:/filma/hot-%07d-.tiff'%(AA//ignayt)
            Freame.save(puti,'TIFF')#'''
        #&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

    print(AA, time.time()-Hronos)
#################################################################################################################################################################################################################################################################################
print("__ Vel dane __")
