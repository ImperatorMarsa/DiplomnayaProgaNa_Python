#!/usr/bin/env python3
from numpy import *
import random, os, time, copy

MassivDlyaPereschetaPeriudGranic=[[0.0,0.0,0.0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]]
#####__Fundomental_Postoyannie__#################################################################################################################################################################################################################################################
pi=3.141592653589793238462643

Vyazkost=2.15e-3 # Paskal*sekunda

Time=1e-10 # sekund
kT=273.16*1.38e-23#1.38e-23
u0=4e-7*pi # Genri/metr

graniciVselennoy=[1e-7]*3
graniciVidemogo=[1e-7]*3
#################################################################################################################################################################################################################################################################################

#####__Chastica__################################################################################################################################################################################################################################################################
class Chastica(object):
    """ento chastica magnitnoy jidkosti"""

    Radiuse=6.66e-9  # metrov
    Plotnost=5000 # kilogramm/metr^3
    Obyom=4/3*pi*Radiuse**3 # metrov^3
    Massa=Obyom*Plotnost # kilogramm
    MagMom=4.78e5*Obyom # Amper*metr^2 ((namagnichenost' nasisheniya = 4.78*10^5 Amper/metr))
    print("Dipol'niy moment chastici ", MagMom)
    def __init__(self, coord, ugol):
        self.Otobrajenie=0

        self.color=array([1,.5,0])

        self.radius=Chastica.Radiuse
        self.mas=Chastica.Massa

        self.pos=coord
        self.skor=array([0, 0, 0])
        self.uskorNew=array([0, 0, 0])
        self.uskorOld=array([0, 0, 0])

        self.axis=ugol*(Chastica.MagMom)
        self.omega=array([0, 0, 0])
        self.ksiNew=array([0, 0, 0])
        self.ksiOld=array([0, 0, 0])

    #@!@!@!@!@!__Osnovnaya_Proga__!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
    def Kinematika(self, forse, moment):
        #&$&$&$&$&$&__Algoritm_Algoritm_Varle__&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
        pom=self.StahostSmesh()
        self.skor=copy.copy(self.pos)
        self.pos=self.pos+forse/self.mas*Time**2+pom[0]
        self.skor=(self.pos-self.skor)/Time

        self.omega=copy.copy(self.axis)
        DeltaAlfa=moment/(2/5*self.radius**2*self.mas)*Time**2+pom[1]
        self.axis=rotate(vector=self.axis, angle=mag(DeltaAlfa), axis=DeltaAlfa)
        self.omega=(self.axis-self.omega)/Time

        self.PorvrkaGrani()
        #&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
    #@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!

    #@!@!@!@!@!__Dopolnitelnie_Progi__!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
    def PorvrkaGrani(self):
        if self.pos[0]>graniciVselennoy[0]: self.pos=array([self.pos[0]-2*graniciVselennoy[0], self.pos[1], self.pos[2]])
        elif self.pos[0]<-graniciVselennoy[0]: self.pos=array([self.pos[0]+2*graniciVselennoy[0], self.pos[1], self.pos[2]])
        if self.pos[1]>graniciVselennoy[1]: self.pos=array([self.pos[0], self.pos[1]-2*graniciVselennoy[1], self.pos[2]])
        elif self.pos[1]<-graniciVselennoy[1]: self.pos=array([self.pos[0], self.pos[1]+2*graniciVselennoy[1], self.pos[2]])
        if self.pos[2]>graniciVselennoy[2]: self.pos=array([self.pos[0], self.pos[1], self.pos[2]-2*graniciVselennoy[2]])
        elif self.pos[2]<-graniciVselennoy[2]: self.pos=array([self.pos[0], self.pos[1], self.pos[2]+2*graniciVselennoy[2]])

    def StahostSmesh(self):
        pom=[]
        difuz=kT/(6.0*pi*self.radius*Vyazkost)
        haotichVec=array([random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)])
        haotichVec/=mag(haotichVec)
        pom.append(haotichVec*((2*difuz*Time)**.5*random.gauss(0, 1)))

        difuz=kT/(8.0*pi*self.radius**3*Vyazkost)
        haotichVec=array([random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)])
        haotichVec/=mag(haotichVec)
        pom.append(haotichVec*((2*difuz*Time)**.5*random.gauss(0, 1)))
        
        return pom
    #@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
#################################################################################################################################################################################################################################################################################

#####__Sili_Deystvuyushie_V_Sisteme__############################################################################################################################################################################################################################################
def MehMomDipolya(y, x): # https://en.wikipedia.org/wiki/Magnetic_dipole
    M=array([0.0,0.0,0.0])
    r=x.pos-y.pos
    if (abs(r[0])-graniciVselennoy[0]<0 and abs(r[1])-graniciVselennoy[1]<0 and abs(r[2])-graniciVselennoy[2]<0):
        B=u0/(4*pi)*(dot(x.axis, r)*3*r/(mag(r))**5-x.axis/(mag(r))**3)
        F_vremennoe=cross(y.axis, B)
        F_vremennoe=array([copysign(F_vremennoe[0], x.pos[0]), copysign(F_vremennoe[1], x.pos[1]), copysign(F_vremennoe[2], x.pos[2])])
        M+=F_vremennoe

    return M

def MehSilaDipolya(y, x): # https://en.wikipedia.org/wiki/Magnetic_dipole
    F=array([0.0,0.0,0.0])
    r=x.pos-y.pos
    if mag(r)-graniciVselennoy[0]<0: #(abs(r[0])-graniciVselennoy[0]<0 and abs(r[1])-graniciVselennoy[1]<0 and abs(r[2])-graniciVselennoy[2]<0):
        F_vremennoe=3*u0/(4*pi*mag(r)**5)*(dot(x.axis, r)*y.axis+dot(y.axis, r)*x.axis+dot(x.axis, y.axis)*r-5*dot(x.axis, r)*dot(y.axis, r)*r/mag(r)**2)
        F_vremennoe=array([copysign(F_vremennoe[0], x.pos[0]), copysign(F_vremennoe[1], x.pos[1]), copysign(F_vremennoe[2], x.pos[2])])
        F+=F_vremennoe
            
    return F

def VneshPole(x):
    B=array([1,0,0])*4.2e-2
    M=cross(x.axis, B)
    return M

def SteerOttalk(x, y): # https://www.desmos.com/calculator/ist2tenoi8
    pom=x.pos-y.pos
    dist=mag(pom)
    N=1e18*4*pi*x.radius**2 # Koncentraciya volosni (N=10^18 1/metr^2)
    d=2.0*x.radius # Diametr chastici
    q=2e-9 # Dlina volosni
    t=2.0*q/d
    if dist<(d+2*q):#x.radius:
        F_ster=2.0*pi*N*kT*(2*d/t*log(d*(1+t)/dist))
        pom=x.pos.norm()*F_ster
        return pom

    return array([0.0,0.0,0.0])

def VyazkTrenLin(x, skor):
    f=-6.0*pi*x.radius*Vyazkost*skor*1e-1
    return f

def VyazkTrenVrash(x, omega):
    f=-8.0*pi*x.radius**3*Vyazkost*omega
    return f
#################################################################################################################################################################################################################################################################################

#####__Storonnie_Funkcii__#######################################################################################################################################################################################################################################################
def NewCoord(x, A): return array([copysign(1, -x.pos[0])*(2*graniciVselennoy[0]*A[0]-abs(x.pos[0])), copysign(1, -x.pos[1])*(2*graniciVselennoy[1]*A[1]-abs(x.pos[1])), copysign(1, -x.pos[2])*(2*graniciVselennoy[2]*A[2]-abs(x.pos[2]))])

def rotate(vector, axis, angle):
    axis = asarray(axis)
    axis = axis/sqrt(dot(axis, axis))
    a = cos(angle/2.0)
    b, c, d = -axis*sin(angle/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    vec=array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)], [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)], [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
    return dot(vec, vector)

def mag(a): return sqrt(a[0]**2+a[1]**2+a[2]**2)
#################################################################################################################################################################################################################################################################################
Chastici=[]
Chastici.append(Chastica(array([0, 0, 0]), array([0, 1, 0])))
N=42 # Chislo chastic
print("Obyomnaya Koncentraciya ", (N+1)*Chastici[-1].Obyom/graniciVselennoy[0]**3)
for A in range(N):
    haotichVec=array([random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)])
    haotichVec/=mag(haotichVec)
    coord=haotichVec*(graniciVselennoy[0]-Chastici[-1].Radiuse)
    C=len(Chastici)
    while C>0:
        C=len(Chastici)
        for x in Chastici:
            r=abs(mag(coord-x.pos))
            if r<(x.Radiuse+2.5e-9)*2:
                haotichVec=array([random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)])
                haotichVec/=mag(haotichVec)
                coord=haotichVec*(graniciVselennoy[0]-Chastici[-1].Radiuse)

            else: C-=1

    haotichVec=array([random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)])
    haotichVec/=mag(haotichVec)
    ugol=haotichVec
    Chastici.append(Chastica(coord, ugol))#'''

input("-__Pognali__-")#time.sleep(4.2)
for AA in range(1000000): # while True:
    Hronos=time.time()
    for odna in Chastici:
        #&$&$&$&$&$__Summa_Vseh_Sil__&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
        forse=array([0.0,0.0,0.0])
        moment=array([0.0,0.0,0.0])#1e-3

        forse=forse+VyazkTrenLin(odna, odna.skor)
        moment=moment+VyazkTrenVrash(odna, odna.omega)
        #moment=moment+VneshPole(odna)
       
        for viteta in Chastici:
            if odna!=viteta:
                for A in MassivDlyaPereschetaPeriudGranic:
                    prizrac=copy.copy(odna)
                    prizrac.pos=NewCoord(prizrac, A)
            
                    forse=forse+SteerOttalk(prizrac, viteta)
                    moment=moment+MehMomDipolya(prizrac, viteta)
                    forse=forse+MehSilaDipolya(prizrac, viteta)

        odna.Kinematika(forse, moment)
        #&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

    print(AA, time.time()-Hronos)
#################################################################################################################################################################################################################################################################################
print("__ Vel dane __")