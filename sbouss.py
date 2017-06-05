import time
import numpy as np
import matplotlib.pyplot as plt
from scitools.std import *
import os

#Create 'indat'-file for sbouss
def indat(Nx,T,dt,dx, typee, ic, ip, redfac, a, readeta=None, terrain=None, dfile=None):
    ida = open('indat', 'w')
    if terrain == 'read' or terrain == 'readfromfile':
        ida.write('!give bathimetry option/readfromfile/\n')
        ida.write('read\n')
        ida.write('!give depthfile/depth.dat/\n')
        ida.write('%s\n' %dfile)
    else:
        ida.write('!give bathimetry option/readfromfile/\n')
        ida.write('h\n')
        ida.write('! give bathimetry option/readfromfile/\n')
        ida.write('flat\n')
        ida.write('! give total length/  100./\n')
        ida.write(' %i\n' %Nx)
    ida.write('!Give number of grid points/ 100/\n')
    ida.write(' %s\n' %(int(round(float(Nx)/dx))))
    ida.write('! give equation type/LSW/\n')
    ida.write(' %s\n' %typee)
    ida.write('!Employ discrete correction term\n')
    ida.write(' no\n')
    ida.write('!give initial condition/readfromfile/\n')
    ida.write('%s\n'%ic)
    if ic== 'read' or ic=='readfromfile':
        ida.write('!give eta-data/eta.in/\n')
        ida.write(' %s \n' %readeta)  
        ida.write('!give u1-data/u1.in/ \n')
        ida.write('stopp\n')
    else:
        ida.write('!give a/h/  0.200000003/\n')
        ida.write(' %s\n' %a)
        ida.write('!give initial position/  50./\n')
        ida.write('%s \n' %ip)
        ida.write('! Prop. towards decreasing x/ja/\n')
        ida.write('no\n')
    ida.write('!Number of cycles\n')
    ida.write('%s\n' %T)
    ida.write('!give time interval\n')
    ida.write('%s\n'%dt)
    ida.write('!gi reduction factor for dt/dx\n')
    ida.write('%s\n' %redfac)
    ida.write('!times for printing\n')
    ida.write('all\n')
    ida.close()



def animation(Nx, T, a , ex, solver, multiplot=False):
    xv = []
    yv = [] 
    for i in range(T+1):
        eta = 'eta'+str(i)
        infile = open(eta, 'r')
        x = []
        y = []
        for row in infile.readlines():
            word = infile.readline()
            A, b = row.split()
            x.append(float(A))
            y.append(float(b))
            
        if multiplot==True:
            xv.append(x)
            yv.append(y)
        else:
            plot(x,y,'-', 
                 axis=[-2, int(Nx)+2, -0.05, float(a)*2],
                 xlabel='x',
                 ylabel=r'$\eta (x)$',
                 title='Exercise: %s \n t=%s' %(ex,i),
                 legend='%s' %solver,
                 grid='on') 
            time.sleep(0.1)
        infile.close()
    if multiplot == True:
        return xv, yv
    raw_input('Press enter to proceed:')

def plotatT(Nx, T, a , ex, solver):
    infile = open('eta'+str(T), 'r')
    x = []
    y = []
    for row in infile.readlines():
        #word = infile.readline()
        a, b = row.split()
        x.append(float(a))
        y.append(float(b))
    infile.close()
    plt.plot(x,y, '-')
    plt.legend('%s' %solver)
    plt.xlabel('x')
    plt.ylabel('$\eta (x)$')
    plt.title('Exercise: %s \n t= %s' %(ex,T))
    plt.grid('on')
    
    raw_input('Press enter to proceed:')
   

def maks(dt, Nx, T):
    dx = [2, 1, 0.5]
    maksy = []
    for i in dx:
        indat(Nx,T,dt,i,'Bouss', 'sol', 20, 0.5, 0.1)
        os.system('rm eta*')
        os.system('rm u*')
        os.system('./sbouss < indat')
        infile = open('eta'+str(T), 'r')
        x = []
        y = []
        for row in infile.readlines():
            word = infile.readline()
            a, b = row.split()
            x.append(float(a))
            y.append(float(b))
        maksy.append(max(y))
        infile.close()
        plt.plot(x,y, '-')
        plt.hold('on')
        plt.legend('dx = %s' %i)
        plt.xlabel('x')
        plt.ylabel('$\eta (x)$')
        plt.title('Exercise: 3 \nt= %s' %T)
        plt.grid('on')
    
    plt.xlabel('x')
    plt.ylabel(r'$\eta (x)$')
    plt.show()
    raw_input('Press enter to proceed:')
    plt.hold('off')
    return dx, maksy
    


def oppgave2():
    Nx = 120
    T = 70
    dt = 1
    dx = 2
    a = 0.1
    
    indat(Nx,T,dt,dx,'Bouss', 'sol', 20, 0.5, a)
    os.system('rm eta*')
    os.system('rm u*')
    os.system('./sbouss < indat')
    animation(Nx, T,0.1, 2, 'Bouss')
    
    dx, maksy = maks(dt=1, Nx=120, T=70) #oppgave2
    print '---------------------------------Max eta--------------------------------'
    for i in range(len(dx)):
        print 'Max eta (dx=%s): %0.5f' %(dx[i], maksy[i])
    

def oppgave3(): #stabilitetsproblem
    T =40
    typ = 'LSW'
    dt = 1
    dx = 0.5
    Nx = 200
    indat(Nx,T,dt,dx, typ, 'sol', 20, dt/dx, 0.1)
    os.system('rm eta*')
    os.system('rm u*')
    os.system('./sbouss < indat')    
    animation(Nx, T, 0.1, 3, typ)

def oppgave4():
    T =80
    typa = 'LSW' #task a
    typb = 'LBouss' #task b
    typc = 'Bouss'
    typd = 'NLSW'
    dt = 1
    dx = 0.1
    Nx = 100
    x0 = 0
    lamb = 20
    read = 'iceta.dat'
    a = 0.1
    def eta0(x):
        diff = x- x0
        if diff > -0.5*lamb and diff < 0.5*lamb:
            return 2*a*np.cos(np.pi*diff/float(lamb))**2
        else:
            return 0
            
    x = np.linspace(-dx/2., Nx+dx/2., int(round(float(Nx)/dx))+2)
    eta = open(read, 'w')
    for i in range(len(x)):
        eta.write('%.10f  %7.10f \n' %(x[i], eta0(x[i])) )
    eta.close()
    
    #Task a    
    indat(Nx,T,dt,dx, typa, 'read', x0, 0.5, a, read)
    os.system('rm eta*')
    os.system('rm u*')
    os.system('./sbouss < indat')    
    animation(Nx, T, a, '4a', typa)
    
    #Task b
    indat(Nx,T,dt,dx, typb, 'read', x0, 0.5, a, read)
    os.system('rm eta*')
    os.system('rm u*')
    os.system('./sbouss < indat') 
    #animation(Nx, T, a, '4b', typb)
    plotatT(Nx, T, a , '4b', typb)
    
    #Task c
    plt.hold('on')
    indat(Nx,T,dt,dx, typc, 'read', x0, 0.5, a, read)
    os.system('rm eta*')
    os.system('rm u*')
    os.system('./sbouss < indat') 
    plotatT(Nx, T, a , '4c', typc)
    plt.hold('off')
    
    #Task d
    indat(Nx,T,dt,dx, typc, 'read', x0, 0.5, a, read)
    os.system('rm eta*')
    os.system('rm u*')
    os.system('./sbouss < indat') 
    animation(Nx, T, a, '4d', typd)
    
    #Task e
    plotatT(Nx, 30, a , '4e', 't=30')
    plt.hold('on')
    plotatT(Nx, 40, a , '4e', 't=40')
    plotatT(Nx, 45, a , '4e', 't=45')
    plotatT(Nx, 50, a , '4e', 't=50')
    plotatT(Nx, 60, a , '4e', 't=60')
    plt.hold('off')
    
def oppgave5(exercise):
    T =200
    typa1 = 'NLSW' #task a
    typa2 = 'Bouss'
    dt = 1
    dx = 0.1
    Nx = 1000
    x0 = 0
    lamb = 100
    read = 'iceta.dat'
    a = 0.1
    def eta0(x):
        diff = x- x0
        if diff > -0.5*lamb and diff < 0.5*lamb:
            return 2*a*np.cos(np.pi*diff/float(lamb))**2
        else:
            return 0
            
    x = np.linspace(-dx/2., Nx+dx/2., int(round(float(Nx)/dx))+2)
    eta = open(read, 'w')
    for i in range(len(x)):
        eta.write('%.10f  %7.10f \n' %(x[i], eta0(x[i])) )
    eta.close()
    
    if exercise == 'a':
        indat(Nx,T,dt,dx, typa1, 'read', x0, 0.5, a, read)
        os.system('rm eta*')
        os.system('rm u*')
        os.system('./sbouss < indat') 
        x1, y1 = animation(Nx, T, a, '5a', typa1, True) 
        indat(Nx,T,dt,dx, typa2, 'read', x0, 0.5, a, read)   
        os.system('rm eta*')
        os.system('rm u*')
        os.system('./sbouss < indat')
        x2, y2 = animation(Nx, T, a, '5a', typa2, True) 
        
        for i in range(len(x1)):
            plot(x1[i], y1[i], '-',
                 x2[i], y2[i], '-',
                 axis=[-2, int(Nx)+2, -0.02, float(a)*2],
                 xlabel='x',
                 ylabel=r'$\eta (x)$',
                 title='Exercise: 5a \n t=%s' %i,
                 legend=('%s' %typa1, '%s' %typa2),
                 grid='on') 
            time.sleep(0.1)
        raw_input('Press enter to proceed:')
        
    if exercise=='b':
        T = 800
        indat(Nx,T,dt,dx, typa2, 'read', x0, 0.5, a, read)
        os.system('rm eta*')
        os.system('rm u*')
        os.system('./sbouss < indat') 
        animation(Nx, T, a, '5b', typa2) 
    
    if exercise=='c':
        a = 0.05
        Nx = 150
        T = 150 
        
        hfile = 'depth.dat'
        h = open(hfile, 'w')
        tis = lambda x: 0.5*x**2
        for i in range(3,154):
            '''    if i < 40:
                h.write('%.10f      1.0000000\n' %i)
                if i > 50:
                h.write('%.10f      0.2000000\n' %i)
            '''
            h.write('%.10f      %5.10f\n' %(i,tis(i)))
        h.close()
        indat(Nx,T,dt,dx, 'LSW', 'sol', x0, 0.5, a, None, 'read', hfile)    
        os.system('rm eta*')
        os.system('rm u*')
        os.system('./sbouss < indat')
        animation(Nx, T, a, '5c', typa2)

def oppgave6():
    T =220
    typa = 'Bouss'
    typb = 'LSW'
    Nx = 130
    x0 = 0
    lamb = 16
    read = 'iceta.dat'
    a = 0.1
    def eta0(x):
        diff = x- x0
        if diff > -0.5*lamb and diff < 0.5*lamb:
            return 2*a*np.cos(np.pi*diff/float(lamb))**2
        else:
            return 0

    #Task i
    dt = 0.5
    dx = 0.5            
    x = np.linspace(-dx/2., Nx+dx/2., int(round(float(Nx)/dx))+2)
    eta = open(read, 'w')
    for i in range(len(x)):
        eta.write('%.10f  %7.10f \n' %(x[i], eta0(x[i])) )
    eta.close()
    

    indat(Nx,T,dt,dx, typa, 'read', x0, dt/dx, a, read)
    os.system('rm eta*')
    os.system('rm u*')
    os.system('./sbouss < indat')
    plotatT(Nx, T, a, '6 i', typa)
    plt.hold('on')
    
    #Task ii
    
    dt = 0.5
    dx = 0.5*np.sqrt(17)
    Nx = 130
    
    x = np.linspace(-dx/2., Nx+dx/2., int(round(float(Nx)/dx))+2)
    eta = open(read, 'w')
    for i in range(len(x)):
        eta.write('%.10f  %7.10f \n' %(x[i], eta0(x[i])) )
    eta.close()
    
    indat(Nx,T,dt,dx, typb, 'read', x0, dt/dx, a, read)
    os.system('rm eta*')
    os.system('rm u*')
    os.system('./sbouss < indat')    
    plotatT(Nx, T, a, '6 ii', typb)
oppgave2()
oppgave3()
oppgave4()
oppgave5('a')
oppgave5('b')
oppgave5('c')
oppgave6()

