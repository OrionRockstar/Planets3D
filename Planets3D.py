import datetime
import math
import random
import threading
import time
from math import sin, cos, tan, atan, pi, sqrt, floor

import numpy as np
from mayavi import mlab

class Planets:
    def __init__(self, name, mass , eccen, perihelion, inclination, argument_of_periapsis,longitude_of_the_ascending_node,mean_longitude):
        self.name=name
        self.m=mass
        self.e=eccen
        self.per=perihelion*1000
        #self.aph=aphelion*1000
        self.aph=floor(perihelion*((1+eccen)/(1-eccen)))*1000
        self.a=(self.per+self.aph)/2
        self.n=math.sqrt(mu/self.a**3)
        self.i=inclination*(pi/180)
        self.w=argument_of_periapsis*(pi/180)
        self.o=longitude_of_the_ascending_node*(pi/180)
        self.l=mean_longitude*(pi/180)
        self.T=math.ceil((2*pi*sqrt((self.a**3)/mu))/86400)
        self.x=[]
        self.y=[]
        self.z=[]
        self.current_x=0
        self.current_y=0
        self.current_z=0
        sys.planetlist.append(self)

    def eccenanom(self): #finds eccentric anomaly
        if sys.orbit:
            M=self.n*(sys.dt)

        if sys.position:
            Mo=self.l-(self.o+self.w)
            M=Mo+(self.n*((dt.d*86400)+dt.my_time))

        guess=M
        for i in range(4):
            E=guess-(M-guess+self.e*sin(guess))/(-1+self.e*cos(guess))
            guess=E
        return E

    def trueanom(self): #finds true anomaly
        E=self.eccenanom()
        v=2*atan(sqrt((1+self.e)/(1-self.e))*tan(E/2))
        return v

    def twodpos(self): #find 2-dimensional position
        v=self.trueanom()
        r=((self.a*(1-self.e**2))/(1+self.e*cos(v)))
        x=r*cos(v)
        y=r*sin(v)
        posvec=np.array([[x],[y],[0]])
        return posvec

    def threedpos(self): #rotates the orbit into the 3-dimensional space
        tran=self.transformation()
        twodposvec=self.twodpos()
        pos=np.transpose(tran@twodposvec)
        magpos=sqrt(pos[0][0]**2+pos[0][1]**2+pos[0][2]**2)
        return pos

    def velocity(self):
        tran=self.transformation()
        v=self.trueanom()
        s=sqrt(Sun.mu/(self.a*(1-self.e**2)))
        velvec=np.transpose(tran@np.array([[-s*math.sin(v)],[s*(self.e+cos(v))],[0]]))
        velmag=math.sqrt(velvec[0][0]**2+velvec[0][1]**2+velvec[0][2]**2)
        return velmag

    def transformation(self): #transformation/rotation matrix
        wtran = np.array([[cos(self.w),sin(self.w),0],[-sin(self.w),cos(self.w),0],[0,0,1]])
        itran = np.array([[1,0,0],[0,cos(self.i),sin(self.i)],[0,-sin(self.i),cos(self.i)]])
        otran = np.array([[cos(self.o),sin(self.o),0],[-sin(self.o),cos(self.o),0],[0,0,1]])

        tran=np.transpose(wtran@itran@otran)
        return tran

class Helios:
    def __init__(self,name,Mass):
        self.name=name
        self.M=Mass
        self.mu=(6.674*10**-11)*self.M

class System:
    def __init__(self):
        self.planetlist=[]
        self.dt=0
        self.orbit=False
        self.position=False

    def orbits(self): #plots the orbits
        self.orbit=True
        for planet in self.planetlist:
            for day in range(planet.T):
                pos=planet.threedpos()
                planet.x.append(pos[0][0]/1000)
                planet.y.append(pos[0][1]/1000)
                planet.z.append(pos[0][2]/1000)
                self.dt+=86400
            planet.x.append(planet.x[0])#connects to start of orbit
            planet.y.append(planet.y[0])
            planet.z.append(planet.z[0])
        self.orbit=False

    def current_pos(self): #finds instantanious postion
        self.position=True
        for planet in self.planetlist:
            curpos=planet.threedpos()
            planet.current_x=curpos[0][0]/1000
            planet.current_y=curpos[0][1]/1000
            planet.current_z=curpos[0][2]/1000
        self.position=False

sys=System()
Sun=Helios('Sun',1.989*10**30)
mu=Sun.mu

Mercury=Planets('Mercury',3.285*10**23,0.2056,46001009,3.38,29.124,48.331,252.25084)
Venus=Planets('Venus',4.867*10**24,0.0068,107476170,3.386,54.884,76.680,181.97973)
Earth=Planets('Earth',5.972*10**24,0.0167,147098291,7.155,114.208,-11.261,100.46435)
Mars=Planets('Mars',6.39*10**23,0.0934,206655215,5.65,286.502,49.558,355.45332)
Jupiter=Planets('Jupiter',1.8982*10**27,0.0489,740679835,6.09,273.867,100.464,34.40438)
Saturn=Planets('Saturn',5.6834*10**26,0.0565,1349823615,5.51,339.392,113.665,49.94432)
Uranus=Planets('Uranus',8.6810*10**25,0.046381,2734998229,6.48,69.99857,74.006,313.23218)
Neptune=Planets('Neptune',1.02413*10**26,0.008678,4459753056,6.43,276.336,131.784,304.88003)

fig=mlab.figure(bgcolor=(0,0,0),size=(2560,1440))

mlab.points3d(0,0,0,color=(1,1,0),scale_factor=(1600000))

def orbitplot():
    mlab.plot3d(Mercury.x,Mercury.y,Mercury.z,color=(.9,.9,.9), representation='wireframe', opacity=.05,line_width=.5)
    mlab.plot3d(Venus.x,Venus.y,Venus.z,representation='wireframe',color=(1,1,0),opacity=.05,line_width=.5)
    mlab.plot3d(Earth.x,Earth.y,Earth.z,color=(0,0,1),representation='wireframe',opacity=.05,line_width=.5)
    mlab.plot3d(Mars.x,Mars.y,Mars.z,color=(1,0,0),representation='wireframe',opacity=.05,line_width=.5)
    mlab.plot3d(Jupiter.x,Jupiter.y,Jupiter.z,color=(1,.5,0),representation='wireframe',opacity=.05,line_width=.5)
    mlab.plot3d(Saturn.x,Saturn.y,Saturn.z,color=(.75,.75,0),representation='wireframe',opacity=.05,line_width=.5)
    mlab.plot3d(Uranus.x,Uranus.y,Uranus.z,representation='wireframe', opacity=.05,line_width=.5)
    mlab.plot3d(Neptune.x,Neptune.y,Neptune.z,color=(0,.75,1),representation='wireframe',opacity=.05,line_width=.5)
t=time.time()
sys.orbits()
print(time.time()-t)
orbitplot()

class daytime:
    def __init__(self): #creates a new class with attributes that can be called and altered in another class method
        t=datetime.datetime.now()
        my_date=datetime.date(t.year,t.month,t.day).toordinal()+1721424.5
        self.d=my_date-2451544.5 #full days since J2000, so positions must my calculate
                  # between now and then ( add time to d)
        self.my_time=(t.hour*3600)+(t.minute*60)

dt=daytime()

class star:
    def __init__(self):
        self.size=random.randrange(1000000,500000000,500000)
        self.color=random.choice([(1,1,1),(1,1,0)])

    def pos(self):
        setattr(self,'position',(random.randint(-1000000000000,1000000000000),random.randint(-1000000000000,1000000000000),random.randint(-1000000000000,1000000000000)))
        while True:
            #setattr(self,'position',(random.randint(-100000000000000,100000000000000),random.randint(-100000000000000,100000000000000),random.randint(-100000000000000,100000000000000)))
            if -100000000<=self.position[0]<=100000000 and -100000000<=self.position[1]<=100000000 and -100000000<=self.position[2]<=100000000:
                setattr(self,'position',(random.randint(-1000000000000,1000000000000),random.randint(-1000000000000,1000000000000),random.randint(-1000000000000,1000000000000)))
            else:
                break

@mlab.animate

def updateAnimation():
    while True:
        t=datetime.datetime.now()
        my_date=datetime.date(t.year,t.month,t.day).toordinal()+1721424.5
        dt.d=my_date-2451544.5
        dt.my_time=(t.hour*3600)+(t.minute*60)
        sys.current_pos()

        mercury.mlab_source.set(x=Mercury.current_x,y=Mercury.current_y,z=Mercury.current_z)
        venus.mlab_source.set(x=Venus.current_x,y=Venus.current_y,z=Venus.current_z)
        earth.mlab_source.set(x=Earth.current_x,y=Earth.current_y,z=Earth.current_z)
        mars.mlab_source.set(x=Mars.current_x,y=Mars.current_y,z=Mars.current_z)
        jupiter.mlab_source.set(x=Jupiter.current_x,y=Jupiter.current_y,z=Jupiter.current_z)
        saturn.mlab_source.set(x=Saturn.current_x,y=Saturn.current_y,z=Saturn.current_z)
        uranus.mlab_source.set(x=Uranus.current_x,y=Uranus.current_y,z=Uranus.current_z)
        neptune.mlab_source.set(x=Neptune.current_x,y=Neptune.current_y,z=Neptune.current_z)
        yield

sys.current_pos()

mercury=mlab.points3d(Mercury.current_x,Mercury.current_y,Mercury.current_z,color=(.9,.9,.9),scale_factor=300000)
venus=mlab.points3d(Venus.current_x,Venus.current_y,Venus.current_z,color=(1,1,0),scale_factor=1000000)
earth=mlab.points3d(Earth.current_x,Earth.current_y,Earth.current_z,color=(0,0,1),scale_factor=1000000)
mars=mlab.points3d(Mars.current_x,Mars.current_y,Mars.current_z,color=(1,0,0),scale_factor=1000000)
jupiter=mlab.points3d(Jupiter.current_x,Jupiter.current_y,Jupiter.current_z,color=(1,.5,0),scale_factor=5000000)
saturn=mlab.points3d(Saturn.current_x,Saturn.current_y,Saturn.current_z,color=(.75,.75,0),scale_factor=5000000)
uranus=mlab.points3d(Uranus.current_x,Uranus.current_y,Uranus.current_z,scale_factor=5000000)
neptune=mlab.points3d(Neptune.current_x,Neptune.current_y,Neptune.current_z,color=(0,.75,1),scale_factor=5000000)

'''for i in range(100):
    st=star()
    st.pos()
    mlab.points3d(st.position[0],st.position[1],st.position[2],color=st.color,scale_factor=st.size)'''
mlab.view(focalpoint=(0,0,0))

updateAnimation()

mlab.show()
