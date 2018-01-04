#!/usr/bin/env python3
# https://www.desmos.com/calculator/oj4vow1hj1
from vpython import *
import random, copy

#####__Fundomental_Postoyannie__#################################################################################################################################################################################################################################################
Radiuse=6.66e-9  # metrov
Plotnost=5000 # kilogramm/metr^3
Massa=4/3*pi*Radiuse**3*Plotnost # kilogramm
MagMom=.73e-9

Vyazkost=2.15e-3 # Paskal*sekunda

Time=.000000000001 # sekund
kT=273.16*1.38e-23#1.38e-23

graniciVselennoy=[7.3e-8]*3
graniciVidemogo=[7.3e-8]*3

rate(1/Time)

Sborshik={'sredSmesh':vector(0,0,0)}
#################################################################################################################################################################################################################################################################################

#####__Sili_Deystvuyushie_V_Sisteme__############################################################################################################################################################################################################################################
def StahostSila(x): # https://en.wikipedia.org/wiki/Brownian_dynamics
    pom=[]
    vek=vector.random().norm()
    rand=random.gauss(0, 1)
    pom.append(vek*((2*Vyazkost*kT)**.5*rand))

    pom.append(vector.random().norm()*((2*Vyazkost*kT)**.5*random.gauss(0, 1)))

    #print('StahostSila\t', pom[0].mag, pom[1].mag)
    return pom

def SteerOttalk(x, y):
    dist=mag(x.pos-y.pos)

    N=66600000 # Koncentraciya volosni
    d=2.0*x.radius # Diametr chastici
    q=42 # Dlina volosni
    t=2.0*q/d
    r=dist
    if dist<13.171:#x.radius:
        pom=x.pos-y.pos
        l=2*(r-d)/d
        ln=log((1+t)/(1+l/2.0))
        F_ster=2.0*pi*N*kT*(d**2.0)/t*(ln-2)
        #print('SteerOttalk\t', F_ster)
        pom=x.pos.norm()*F_ster
        return pom
        
    return vector(0,0,0)

def VyazkTrenLin(x, skor):
    f=-6.0*pi*x.radius*Vyazkost*skor
    #print('VyazkTrenLin\t{:+.4E}'.format(f.mag))
    return f

def VyazkTrenVrash(x, omega):
    f=-8.0*pi*x.radius**3*Vyazkost*omega
    #print('VyazkTrenVrash\t', f.mag)
    return f
#################################################################################################################################################################################################################################################################################

#####__Storonnie_Funkcii__#######################################################################################################################################################################################################################################################
def PorvrkaGrani(x):
    if x.pos.x>graniciVselennoy[0]: x.pos=vector(x.pos.x-2*graniciVselennoy[0], x.pos.y, x.pos.z)
    elif x.pos.x<-graniciVselennoy[0]: x.pos=vector(x.pos.x+2*graniciVselennoy[0], x.pos.y, x.pos.z)
    if x.pos.y>graniciVselennoy[1]: x.pos=vector(x.pos.x, x.pos.y-2*graniciVselennoy[1], x.pos.z)
    elif x.pos.y<-graniciVselennoy[1]: x.pos=vector(x.pos.x, x.pos.y+2*graniciVselennoy[1], x.pos.z)
    if x.pos.z>graniciVselennoy[2]: x.pos=vector(x.pos.x, x.pos.y, x.pos.z-2*graniciVselennoy[2])
    elif x.pos.z<-graniciVselennoy[2]: x.pos=vector(x.pos.x, x.pos.y, x.pos.z+2*graniciVselennoy[2])
#################################################################################################################################################################################################################################################################################

gray = color.gray(0.7)
d = graniciVidemogo[0]
r = 0.005
boxbottom = curve(color=gray)
boxbottom.append([vector(-d,-d,-d), vector(-d,-d,d), vector(d,-d,d), vector(d,-d,-d), vector(-d,-d,-d)])
boxtop = curve(color=gray)
boxtop.append([vector(-d,d,-d), vector(-d,d,d), vector(d,d,d), vector(d,d,-d), vector(-d,d,-d)])
vert1 = curve(color=gray)
vert2 = curve(color=gray)
vert3 = curve(color=gray)
vert4 = curve(color=gray)
vert1.append([vector(-d,-d,-d), vector(-d,d,-d)])
vert2.append([vector(-d,-d,d), vector(-d,d,d)])
vert3.append([vector(d,-d,d), vector(d,d,d)])
vert4.append([vector(d,-d,-d), vector(d,d,-d)])

Chastici=[]
for A in range(1):
    coord=vector(0, 0, 0)#vector.random()*graniciVselennoy[0]
    ugol=vector.random().norm()
    Chastici.append((
        sphere(
            radius=Radiuse, 
            color=vector(1, .5, 0), 
            make_trail=True, 
            opacity=0.3, 
            mas=Massa, 

            pos=coord, 
            skor=vector(0, 0, 0), 
            uskorNew=vector(0, 0, 0), 
            uskorOld=vector(0, 0, 0)),

        arrow(pos=coord, shaftwidth=MagMom, color=vector(1, 1, 0), opacity=0.3, 
            axis=ugol*(MagMom+2*Radiuse), 
            omega=vector(0, 0, 0), 
            ksiNew=vector(0, 0, 0), 
            ksiOld=vector(0, 0, 0))))

for A in range(1000000): # while True:
    for odna in Chastici:
        #&$&$&$&$&$&__Algoritm_Algoritm_Varle__&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
        buf=copy.copy(odna[0].pos)
        odna[0].pos=odna[0].pos+odna[0].skor*Time+odna[0].uskorNew*Time**2/2
        odna[0].skor=odna[0].skor+(odna[0].uskorNew+odna[0].uskorOld)*Time/2
        odna[1].axis=norm(odna[1].axis+odna[1].omega*Time+odna[1].ksiNew*Time**2/2)*(MagMom+2*Radiuse)
        odna[1].omega=odna[1].omega+(odna[1].ksiNew+odna[1].ksiOld)*Time/2

        buf=odna[0].pos-buf
        PorvrkaGrani(odna[0])        
        Sborshik['sredSmesh']+=buf
        #print("{:+.4E}\t{:+.4E}\t{:+.4E}\t\t{:+.4E}".format(Sborshik['sredSmesh'].x, Sborshik['sredSmesh'].y, Sborshik['sredSmesh'].z, Sborshik['sredSmesh'].mag))

        odna[1].pos=odna[0].pos
        #&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$

        #&$&$&$&$&$__Summa_Vseh_Sil__&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
        forse=vector(0,0,0)
        moment=vector(0,0,0)

        #forse=forse+VyazkTrenLin(odna[0], copy.copy(odna[0].skor))
        #moment=moment+VyazkTrenVrash(odna[0], copy.copy(odna[1].omega))
        pom=StahostSila(odna[0])
        forse=forse+pom[0]
        moment=moment+pom[1]

        for viteta in Chastici:
            if odna[0]!=viteta[0]:
                forse=forse+SteerOttalk(odna[0], viteta[0])
        
        odna[0].uskorOld=copy.copy(odna[0].uskorNew)
        odna[0].uskorNew=forse/odna[0].mas
        odna[1].ksiOld=copy.copy(odna[1].ksiNew)
        odna[1].ksiNew=moment/(2/5*odna[0].radius**2*odna[0].mas)
        #&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$&$
