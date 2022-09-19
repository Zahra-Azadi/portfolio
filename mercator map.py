import matplotlib.pyplot as plt
import numpy as np
fig = plt.figure()
ax = plt.axes()
f=[-1.3962,-1.047,-0.6981,-0.3490,0,0.3490,0.6981,1.047,1.3962]
#dar khat bala arzjoghrafiai bar hasb radian az -80 ta 80 daraje amade ast. dar 90 be binahayat miravad.
w=10.0
#tool naghshe 10 entekhab shode
x = np.linspace(0, 10, 1000)
y = np.linspace(-180, 180, 1000)
for i in range(len(f)):
 ax.hlines(y=np.log10(1/np.cos(f[i])+np.tan(f[i]))*w/(2*3.14),xmin=-5,xmax=5,lw=2);
t = [-180,-160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180]
#dar khat bala tooljoghrafiai bar hasb daraje az 0 ta 360 amade ast
for i in range(len(t)):
 ax.vlines(x=w*t[i]/360,ymin=-1.7,ymax=1.7,lw=2);
ti=np.array([[44.338,48.11,51.53,56.56,61.08,61.76,63.24,61.77,56.66,51.96,48.49,45.70],
             [0.685,0.6908,0.639,0.666,0.637,0.542,0.47,0.435,0.471,0.484,0.523,0.586]])
#ti mokhtasat 12 noghte marzi iran ast
for n in range(len(ti)):
  plt.plot(w*ti[0:1]/360,np.log10(1/np.cos(ti[1:2])+np.tan(ti[1:2]))*w/(2*3.14),".");
#great circles from a ponit
t00=float(input("latitude in degree"))
f00=float(input("longitude in degree"))
t0=(90-t00)*3.14/180
f0=f00*3.14/180
x0=np.sin(t0)*np.cos(f0)
y0=np.sin(t0)*np.sin(f0)
z0=np.cos(t0)
x1=np.linspace(-3,3,1000)
b=[]
xx1=[]
for i in range(len(x1)):
  p = [1-((x0/2)**2),(-2*(y0/2)*((y0/2)**2))+2*(x0/2)*(y0/2)*(x1[i]-(x0/2))-2*(z0/2)*(z0/2)*(y0/2),-(z0/2)**2+((x0/2)**2)*(x1[i]**2+(x0/2)**2-2*(x0/2)*x1[i])+((z0/2)**2)*(x1[i]**2)+((y0/2)**2)*((y0/2)**2)+2*(x0/2)*(y0/2)*((x0/2)*(y0/2)-x1[i]*(y0/2))+((z0/2)**2)*((z0/2)**2)-2*(z0/2)*(z0/2)*((x0/2)*(x1[i]-(x0/2))-(y0/2)*(y0/2))]
  a=np.roots(p)
  if np.isreal(max(a))==True:
   xx1.append(x1[i])
   b.append([max(a)])
e=round(len(xx1)/12)
xx11=[]
bb=[]
for i in range(len(xx1)):
 if i%e==0:
  xx11.append(xx1[i])
  bb.append(b[i])
y1=np.array(bb,dtype=float)
zb=[]
for i,j in zip(xx11,y1):
  zb.append(-((x0/2)*((i)-(x0/2))+(y0/2)*((j)-(y0/2)))/(z0/2))
z1=np.array(zb,dtype=float)
ew=np.column_stack((xx11,y1,z1))
a0=(x0,y0,z0)
tita=[]
for i in range(len(xx11)):
 c0=np.cross(a0,ew[i])
 tita.append(c0/((np.dot(c0,c0))**0.5))
tita=np.array(tita)
for i in range(len(tita)):
 n1=tita[i,0]
 n2=tita[i,1]
 n3=tita[i,2]
 q=y0*n2+n1*x0+n3*z0
 b=4*((n3**2)+(n2**2))*(-2*n1)*q
 c=4*((n3**2)+(n2**2))*((n1**2)*(x0**2)+(n2**2)*(y0**2)+(n3**2)*(z0**2)+2*n1*n2*x0*y0+2*z0*n3*(n1*x0+n2*y0)-(n3**2));
 d=4*(n2**2)*(n1**2)-4*((n3**2)+(n2**2))*((n3**2)+(n1**2));
 e=-2*q*n1-b
 f=(q**2)-c
 p1=[d,e,f]
 xmax=max(np.roots(p1))
 xmin=min(np.roots(p1))
 x1=np.linspace(xmin,xmax,100)
 b1=[]
 d1=[]
 for i in range(len(x1)):
  p = [1-(n1**2),(-2*y0*(n2**2))+2*n1*n2*(x1[i]-x0)-2*z0*n3*n2,-n3**2+((n1**2)*(x1[i]**2+x0**2-2*x0*x1[i]))+(n3**2)*(x1[i]**2)+(n2**2)*(y0**2)+2*n1*n2*(x0*y0-x1[i]*y0)+(z0**2)*(n3**2)-2*z0*n3*(n1*(x1[i]-x0)-n2*y0)]
  a=np.roots(p)
  b1.append([max(a)])
  d1.append([min(a)])
 y1=np.array(b1,dtype=float)
 y2=np.array(d1,dtype=float)
 zb=[]
 zd=[]
 for i,j in zip(x1,y1):
  zb.append(-(n1*(i)+n2*(j))/n3)
 for i,j in zip(x1,y2):
  zd.append(-(n1*(i)+n2*(j))/n3)
 z1=np.array(zb,dtype=float)
 z2=np.array(zd,dtype=float)
 mokh=np.column_stack((x1,b1,z1))
 mokh2=np.column_stack((x1,d1,z2))
 fi11=[]
 teta11=[]
 fi12=[]
 teta12=[]
 for i in range(0,49):
  fi11.append(np.arctan(mokh[i,1]/mokh[i,0])*180/3.14)
  teta11.append(90-(np.arccos(mokh[i,2]/(mokh[i,0]**2+mokh[i,1]**2+mokh[i,2]**2)**0.5)*180/3.14))
  teta12.append(90-(np.arccos(mokh[50+i,2]/(mokh[50+i,0]**2+mokh[50+i,1]**2+mokh[50+i,2]**2)**0.5)*180/3.14))
  fi12.append(-180+(np.arctan(mokh[50+i,1]/mokh[50+i,0])*180/3.14))
 fi21=[]
 fi22=[]
 teta21=[]
 teta22=[]
 for i in range(0,49):
  fi21.append(np.arctan(mokh2[i,1]/mokh2[i,0])*180/3.14)
  teta21.append(90-(np.arccos(mokh2[i,2]/(mokh2[i,0]**2+mokh2[i,1]**2+mokh2[i,2]**2)**0.5)*180/3.14))
  fi22.append(180+(np.arctan(mokh2[50+i,1]/mokh2[50+i,0])*180/3.14))
  teta22.append(90-(np.arccos(mokh2[50+i,2]/(mokh2[50+i,0]**2+mokh2[50+i,1]**2+mokh2[50+i,2]**2)**0.5)*180/3.14))
 xam21=[]
 yam21=[]
 for i,j in zip(fi21,teta21):
  xam21.append(w*(i)/360-5)
  yam21.append(np.log10(1/np.cos((j)*3.14/180.0)+np.tan((j)*3.14/180.0))*w/(2*3.14))
 for i in range(len(xam21)):
  if xam21[i]<-5:
   xam21[i]=xam21[i]+10
 plt.plot(xam21,yam21,".")
 xam22=[]
 yam22=[]
 for i,j in zip(fi22,teta22):
  xam22.append(w*(i)/360-5)
  yam22.append(np.log10(1/np.cos((j)*3.14/180.0)+np.tan((j)*3.14/180.0))*w/(2*3.14))
 plt.plot(xam22,yam22,".")
 xam11=[]
 yam11=[]
 for i,j in zip(fi11,teta11):
  xam11.append(w*(i)/360+5)
  yam11.append(np.log10(1/np.cos((j)*3.14/180.0)+np.tan((j)*3.14/180.0))*w/(2*3.14))
 for i in range(len(xam11)):
  if xam11[i]>5:
   xam11[i]=xam11[i]-10
 plt.plot(xam11,yam11,".")
 xam12=[]
 yam12=[]
 for i,j in zip(fi12,teta12):
  xam12.append(w*(i)/360+5)
  yam12.append(np.log10(1/np.cos((j)*3.14/180.0)+np.tan((j)*3.14/180.0))*w/(2*3.14))
 plt.plot(xam12,yam12,".")
