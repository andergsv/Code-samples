'''
refinement test: solitary wave propagation throgh L-bend

globoussrefine & points: Details shared with GloBouss and FEnics
soltest: generate solitary wave shape and velocity
 
'''
from dolfin import *
import numpy as np
import copy as cp
import time as time
import math
import os
from mshr import *
import matplotlib.pyplot as plt

set_log_active(False)                                   
parameters['reorder_dofs_serial'] = False


eps = 1              #Dispersion
H = 1                #Non-dimensional Depth
alpha = Constant(1)  #Non - liearity
h = Constant(H)      #FEniCS Depth 

T = 40               #Endtime


    
def globousrefine(nx, ny):
    fenx, feny = int(point2[0]*nx), int(point2[1]*ny)
    glodx = float(1)/nx
    glody = float(1)/ny
    return fenx, feny, glodx, glody 


def point(b1):
    point1 = (0,36-b1)
    point2 = (36, 56)
    return point1, point2, b1

def soltest( fil, slag, s0, dt):
    
    os.system('/home/geirkp/bin/soliprod -%s -%s -ds 2000  -topp %s -a 0.1 > %s' \
                                                   %(fil, slag, 15+dt, slag))
    r = open(fil, 'r')
    x = []
    y = []
    for line in r:
        one, two = line.split() 
        x.append(float(one))
        y.append(float(two))
    r.close()
    return x, y


###############################################################################
##                               FEM-SOLVER                                  ##
###############################################################################

def RecMesh(Nx, Ny, kappa, mesh, aprox=False):
    
    if Nx > 8:
        mesh = refine(mesh)
    
    dtcc = mesh.hmin()
    NT = math.ceil((T-dtcc/2.)/dtcc)
    dt = (T-dtcc/2.)/NT
    print dt*NT +dt/2, dt, mesh.hmin()
    
    xu0, nu0 = soltest('t', 'hast',0, 0)
    x0, neta0 = soltest('t', 'eta',0, dt/2.)
    
    
    V = VectorFunctionSpace(mesh, 'CG', 2)
    Q = FunctionSpace(mesh, 'CG', 1) 
    
    ######################### APPROX SOLITON ###################################
    if aprox == 'True':
        c = lambda a: np.sqrt(((1+a)*np.log(1+a)-a)*6*(1+a)**2/(a**2*(3+2*a)))
        beta = lambda a: np.sqrt(3*a/(4*(1+0.68*a)))

        amp = 0.1
        xs0 = 15
        eta00 = Expression('a*4*pow(cosh(B*(x[0]-x0)),2)/pow(cosh(2*B*(x[0]-x0)) + 1,2)/(1+a*pow(tanh(B*(x[0]-x0)),2))',\
                                 a = 0.1, B = beta(amp), x0 = xs0+dt/2., element=FiniteElement("Lagrange", triangle, 1))

        from sympy.utilities.codegen import ccode
        from sympy import symbols
        import sympy as sp

        x = symbols('x[0]')

        phi = c(amp)*amp/(beta(amp)*(1+amp))*sp.tanh(beta(amp)*(x-xs0))
        dphi = phi.diff(x,1)

        u0 = interpolate(Expression((ccode(dphi),'0'), degree=5), V)
        eta0 = interpolate(eta00,Q)
    ######################### SOLITON ##########################################
    else:
        u_ = interpolate(Expression(('x[0]','x[1]'), degree=5),V)
        #eta_ = interpolate(Expression('x[0]', degree=5),Q)
        
        coords = mesh.coordinates()
        vtd = dof_to_vertex_map(Q)
        eta_ = np.zeros(len(vtd))

        #rr = np.zeros(len(u_.vector().array())/2))
        rr = np.zeros((len(u_.vector().array())/2, 3))
        for i in range(len(u_.vector().array())/2):
            tt = [u_.vector().array()[i+len(u_.vector().array())/2], \
                                         u_.vector().array()[i], i]
            #rr[i]=tt
            rr[i][0], rr[i][1], rr[i][2] = tt[0], tt[1], int(tt[2])
            
        uup = np.zeros(len(u_.vector().array())) 
        for pos in rr:
            if pos[1]<30:
                i = int(pos[2])
                x = pos[1]
                uinp = np.interp(x, xu0, nu0)
                uup[i] = uinp

        for i in range(len(coords)):
            if coords[i][0]<30:
                x = coords[i][0]
                y = np.interp(x, x0, neta0)
                eta_[vtd[i]] = y

        eta0 = Function(Q)       
        eta0.vector()[:] = eta_
        u0 = Function(V)
        u0.vector()[:] = uup
        uup = 0   

    ###########################################################################
    
    v = TestFunction(V)
    q = TestFunction(Q) 

    u = TrialFunction(V)
    eta = TrialFunction(Q)
       
    n = FacetNormal(mesh)
    
    F1 = Constant(1/dt)*dot((u-u0),v)*dx \
        + eps/(Constant(3)*dt)*h*h*inner(div(u-u0),div(v))*dx \
        + dot(grad(eta0),v)*dx 
    
    F1 += -eps/Constant(3*dt)*h*h*dot(v,n)*div(u-u0)*ds \
          -eps/Constant(3*dt)*h*h*(dot(u,n)-dot(u0,n))*div(v)*ds \
          +Constant(kappa)*eps/Constant(3*dt)*dot(v,n)*(dot(u0,n)-dot(u,n))*ds
    
    F2 = Constant(1/dt)*(eta-eta0)*q*dx \
        -dot(u0, grad(q))*dx + inner(dot(u0,n), q)*ds
    
    A1 = assemble(lhs(F1))
    
    A2 = assemble(lhs(F2))
    
    F3 = Constant(alpha)*inner(0.5*(dot(u0,grad(u))+dot(u,grad(u0))),v)*dx
    
    
    F4 =- Constant(alpha)*dot(eta*u0, grad(q))*dx 
    
    coords = mesh.coordinates()
    vtd = dof_to_vertex_map(Q)

    c = 0
    t = dt
    file = File('para/b'+str(b)+str(int(NT))+'eta.pvd')
    while t<T:
        uu = Function(V)     
        b1 = assemble(rhs(F1)+rhs(F3))
        M1 =assemble(lhs(F3))
        solver1 = LUSolver(A1+M1)
        solver1.solver_paremeters={'linear_solver': 'cg', 'preconditioner':'lu'}
        solver1.solve(uu.vector(),b1)
        u0.assign(uu)
        M2 = assemble(lhs(F4))
        
        solver2 = LUSolver(A2+M2)
        solver2.solver_paremeters={'linear_solver': 'cg', 'preconditioner':'lu'}
        etat = Function(Q)    
        b2 = assemble(rhs(F2)+rhs(F4))
        
        solver2.solve(etat.vector(),b2)
        eta0.assign(etat)    
        print  t/T*100," percent complete         \r",
        if t > T-dt: 
            print t, dt, c
            file << eta0
            c += 1
            plot(eta0, interactive=True)
        t += dt
        
    return etat, mesh.num_vertices(), mesh.hmin(), dt, mesh

Nf = [8, 16, 32, 64, 128, 256,512,1024*2]

B = [1]

for b in B:
    rec1 = Rectangle(Point(0, 0), Point(36, b))
    rec2 = Rectangle(Point(36-b, b), Point(36,56))
    domain = rec1 + rec2    
    mesh = generate_mesh(domain, 8) 
    
    out = open('datweakpicardASSEMBLETESTapro'+str(b), 'w')
    out.write('fe  hh dt vert N time\n')
    out.close()
    point1, point2, b1 = point(b) 
    print 'b=', b
    for i in Nf:
        Ny = 20    
        gi = i
        fenx, feny, gdx, gdy = globousrefine(gi, gi)
        start = time.time()
        eta, numvert, hm, DT, mesh = RecMesh(i,Ny, 800, mesh, True)
        
        print np.amax(eta.vector().array())
        
        stop = time.time()
        print stop - start
        out = open('datweakpicardASSEMBLETESTapro'+str(b), 'r+')
        out.readlines()        
        out.write('%s %s %s %s %s %s\n' %(np.amax(eta.vector().array()), hm, DT, numvert, i, stop-start))
        out.close()
    
