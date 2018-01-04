#!/usr/bin/env python3
from PIL import ImageGrab, ImageDraw, ImageFont
from vpython import *
import random, os, time, copy

scene=canvas(width=1346, height=748)
MassivDlyaPereschetaPeriudGranic=[[0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1], [0,0,0]]
#####__Fundomental_Postoyannie__#################################################################################################################################################################################################################################################
Vyazkost=2.15e-3 # Paskal*sekunda

Time=1e-10 # sekund
kT=273.16*1.38e-23#1.38e-23
u0=4e-7*pi # Genri/metr

graniciVselennoy=[1e-7]*3
graniciVidemogo=[1e-7]*3

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
        self.uskorNew=vector(0, 0, 0)
        self.uskorOld=vector(0, 0, 0)

        self.axis=ugol*(Chastica.MagMom)
        self.omega=vector(0, 0, 0)
        self.ksiNew=vector(0, 0, 0)
        self.ksiOld=vector(0, 0, 0)
            
        self.sfera=sphere(
            radius=self.radius, 
            color=self.color, 
            make_trail=False, 
            trail_color=vector(0, 1, 0), 
            opacity=.3, 
            pos=self.pos)
        self.strela=arrow(
            pos=self.sfera.pos, 
            #shaftwidth= 
            color=vector(1, 1, 0), 
            opacity=.3, 
            axis=norm(self.axis)*(Chastica.Radiuse*3/2))          

    #@!@!@!@!@!__Osnovnaya_Proga__!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!@!
    def Kinematika(self, forse, moment): # https://www.desmos.com/calculator/bhjmf8p0pf
        #&$&$&$&$&$&__Algoritm_Algoritm_Varle__&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
        pom=self.StahostSmesh()
        self.skor=copy.copy(self.pos)
        self.pos=self.pos+forse/(6.0*pi*self.radius*Vyazkost)*Time+0*pom[0]
        self.skor=(self.pos-self.skor)/Time

        self.omega=copy.copy(self.axis)
        DeltaAlfa=moment/(8.0*pi*self.radius**3*Vyazkost)*Time+0*pom[1]
        self.axis=self.axis.rotate(angle=DeltaAlfa.mag, axis=DeltaAlfa)
        self.omega=(self.axis-self.omega)/Time

        self.PorvrkaGrani()
        self.Otobrajenie+=1
        self.sfera.color=self.color
        self.sfera.pos=self.pos
        self.strela.pos=self.sfera.pos
        self.strela.axis=norm(self.axis)*(Chastica.Radiuse*5/3)


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
    M=vector(0,0,0)
    r=y.pos-x.pos
    if (abs(r.x)-graniciVselennoy[0]<0 and abs(r.y)-graniciVselennoy[1]<0 and abs(r.z)-graniciVselennoy[2]<0):
        B=u0/(4*pi)*(dot(x.axis, r)*3*r/(r.mag)**5-x.axis/(r.mag)**3)
        M=cross(y.axis, B)

    return M

def MehSilaDipolya(y, x): # https://en.wikipedia.org/wiki/Magnetic_dipole
    F=vector(0,0,0)
    r=y.pos-x.pos
    if (abs(r.x)-graniciVselennoy[0]<0 and abs(r.y)-graniciVselennoy[1]<0 and abs(r.z)-graniciVselennoy[2]<0):
        F=3*u0/(4*pi*r.mag**5)*(dot(x.axis, r)*y.axis+dot(y.axis, r)*x.axis+dot(x.axis, y.axis)*r-5*dot(x.axis, r)*dot(y.axis, r)*r/r.mag**2)
        
        n=Chastici.index(y)
        Streli[n][0].pos=y.pos
        Streli[n][0].axis=F.norm()*y.radius

    return F

def VneshPole(x):
    B=vector(1,0,0)*4.2e-2
    M=cross(x.axis, B)
    return M

def SteerOttalk(x, y): # https://www.desmos.com/calculator/ist2tenoi8
    pom=x.pos-y.pos
    dist=mag(pom)
    N=1e18*4*pi*x.radius**2 # Koncentraciya volosni (N=10^18 1/metr^2)
    d=2.0*x.radius # Diametr chastici
    q=2e-9 # Dlina volosni v metrah
    t=2.0*q/d
    if dist<(d+2*q):#x.radius:
        F_ster=2.0*pi*N*kT*(2*d/t*log(d*(1+t)/dist))
        pom=pom.norm()*F_ster

        n=Chastici.index(x)
        Streli[n][1].pos=x.pos
        Streli[n][1].axis=pom.norm()*y.radius

        return pom

    return vector(0,0,0)
#################################################################################################################################################################################################################################################################################

#####__Storonnie_Funkcii__#######################################################################################################################################################################################################################################################
def NewCoord(x, A): 
    a=vector(copysign(1, -x.pos.x)*(2*graniciVselennoy[0]*A[0]-abs(x.pos.x)), copysign(1, -x.pos.y)*(2*graniciVselennoy[1]*A[1]-abs(x.pos.y)), copysign(1, -x.pos.z)*(2*graniciVselennoy[2]*A[2]-abs(x.pos.z)))
    #sphere(pos=a, radius=6.66e-9)
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

Streli=[[arrow(pos=vector(0,0,0), axis=vector(0,0,0), color=vector(1,0,0)), arrow(pos=vector(0,0,0), axis=vector(0,0,0), color=vector(0,1,0))], [arrow(pos=vector(0,0,0), axis=vector(0,0,0), color=vector(1,0,0)), arrow(pos=vector(0,0,0), axis=vector(0,0,0), color=vector(0,1,0))]]

Chastici=[]
Chastici.append(Chastica(vector(0.1e-9, 0.1e-9, 0.1e-9), vector(0, 1, 0)))
#Chastici.append(Chastica(vector(-graniciVselennoy[0]*1/4, -1e-9, -1e-9), vector(1, 0, 0)))
#Chastici.append(Chastica(vector(graniciVselennoy[0]*1/4, -1e-9, -1e-9), vector(1, 0, 0)))
Chastici.append(Chastica(vector((6.66+2.2)*2e-9, 0, 0), vector(0, -1, 0)))
#Chastici.append(Chastica(vector(-(6.66+2)*2e-9, 0, 0), vector(-sin(pi/10), cos(pi/10), 0)))
N=0 # Chislo chastic
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

input("-__Pognali__-")#time.sleep(4.2)
for AA in range(1000000): # while True:
    Hronos=time.time()
    for odna in Chastici:
        #&$&$&$&$&$__Summa_Vseh_Sil__&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
        forse=vector(0,0,0)
        moment=vector(0,0,0)

        #moment=moment+VneshPole(odna)       
        for viteta in Chastici:
            if odna!=viteta:
                for A in MassivDlyaPereschetaPeriudGranic:
                    prizrac=copy.copy(viteta)
                    prizrac.pos=NewCoord(prizrac, A)
                    
                    pom1, pom2=SteerOttalk(odna, prizrac), MehSilaDipolya(odna, prizrac)
                    if pom2.mag>0:
                        print("Strich ottalk = {:.4E} \t Dipol-Dpol = {:.4E}".format(pom1.mag, pom2.mag))
                    forse+=pom1+pom2

                    #forse+=SteerOttalk(odna, prizrac)+MehSilaDipolya(odna, prizrac)
                    moment+=MehMomDipolya(odna, prizrac)

        odna.Kinematika(forse, moment)
        '''if AA%100==0:#
            Freame=ImageGrab.grab()#grab([0,0,300,300])      
            puti='D:/filma/hot-%07d-.tiff'%(AA//100)
            Freame.save(puti,'TIFF')#'''
        #&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

    #print(AA, time.time()-Hronos)
#################################################################################################################################################################################################################################################################################
print("__ Vel dane __")