from vpython import *

u0=4e-7*pi # Genri/metr

scene=canvas(width=2000, height=2000)

sfera=sphere(radius=1, color=vector(1, .5, 0), make_trail=False, trail_color=vector(0, 1, 0), opacity=1, pos=vector(0,0,0), axis=vector(1,0,0))

def dipolFil(pos): # https://en.wikipedia.org/wiki/Magnetic_dipole
	p=vector(1,0,0)
	E=vector(0,0,0)
	r=E-pos
	E=u0/(4*pi)*(dot(p, r)*3*r/(r.mag)**5-p/(r.mag)**3)
	#arrow(pos=pos, shaftwidth=.1, color=(pos+vector(20,20,0))/20, opacity=.5, axis=E*7300)
	#sphere(radius=.2, color=vector(1, 0, 0), pos=pos)

def Krasota():
	A=73
	r=10
	for x in range(0, A+1):
		B=x
		x/=A
		for y in range(-abs(B), abs(B)+1):
			y+=.1
			y/=(B+.1)
			dipolFil(vector(r*cos(y*pi)*cos(x*pi), r*sin(y*pi)*cos(x*pi), r*sin(x*pi)))

