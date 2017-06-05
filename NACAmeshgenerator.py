'''
Mesh generator of NACA airfoil in OpenFOAM. 
'''
import matplotlib.pyplot as plt
import numpy as np
import copy as cp

def naca(plot=False, meshdict=False, AoA = None): #AoA in radians
    print 'Write the four-digit NACA type you want to plot:'
    Naca = raw_input()
    dgt = [int(i) for i in str(Naca)]
    if len(dgt) > 4 or len(dgt)< 4:
        print 'NACA number needs to be a 4-digit number!'
    else:
        md = float(dgt[0])
        pd = float(dgt[1])
        lst_digits = float(str(dgt[2])+str(dgt[3]))

        print 'Write desired chord length:'
        c = float(raw_input())
        
        if AoA == None:
            alpha = 0
            #AoA = 0
        else:
            alpha=AoA
        
        
        #symetric NACA airfoil
        def yt(lst_digits,c,x):
            t = lst_digits/100.
            return 5*t*c*(0.2969*np.sqrt(x/c)+ (-0.1260)*(x/c) + \
                (-0.3516)*(x/c)**2 + \
                0.2843*(x/c)**3 + (-0.1036)*(x/c)**4) #(-0.1015) siste

                
        def ytAoA(lst_digits, c, x):
            y = yt(lst_digits, c, x)
            uplow = np.zeros([2,2]) 
            uplow[0,0] = np.cos(alpha)
            uplow[0,1] = np.sin(alpha)
            uplow[1,0] = -np.sin(alpha)
            uplow[1,1] = np.cos(alpha)
            upper = np.zeros([2,len(x)])
            lower = np.zeros([2,len(x)])
            for i in range(len(x)):
                upper[0,i] = uplow[0,0] * x[i] + uplow[0,1] * y[i]
                upper[1,i] = uplow[1,0] * x[i] + uplow[1,1] * y[i]
                lower[0,i] = uplow[0,0] * x[i] + uplow[0,1] * -y[i]
                lower[1,i] = uplow[1,0] * x[i] + uplow[1,1] * -y[i]
            xU = upper[0,:]
            yU = upper[1,:]
            xL = lower[0,:]
            yL = lower[1,:]
                
            return xU, xL, yU, yL

        #Mean camber line    
        def yc(md,pd, c, x):
            m = md/100.
            p = pd/10.
            if x >= 0 and x <= p*c:
                return m* x/p**2 *(2*p - x/c)
            else:
                return m * (c-x)/(1-p)**2 * (1 + x/c-2*p)
            
        #cambered airfoil           
        def xyUL(x, md, pd, lst_digits, c):
            xU = np.zeros(len(x))
            xL = np.zeros(len(x))
            yU = np.zeros(len(x))
            yL = np.zeros(len(x))
            def theta(x, md, pd, c):
                m = md/100.
                p = pd/10.
        
                if x >= 0 and x<= p*c:
                    return np.arctan(2* m / p**2 * (p-x/c))
                else: 
                    return np.arctan(2*m/(1-p)**2 * (p-x/c))
                
            for i in range(len(x)):
                xU[i] = x[i] - yt(lst_digits, c, x[i]) * np.sin(theta(x[i],md,pd,c))
                xL[i] = x[i] + yt(lst_digits, c, x[i]) * np.sin(theta(x[i],md,pd,c))
                yU[i] = yc(md,pd,c,x[i]) + yt(lst_digits, c, x[i])*np.cos(theta(x[i],md,pd,c))
                yL[i] = yc(md,pd,c,x[i]) - yt(lst_digits, c, x[i])*np.cos(theta(x[i],md,pd,c))
            
            if AoA == None:
                return xU, xL, yU, yL
            
            else:
                uplow = np.zeros([2,2]) 
                uplow[0,0] = np.cos(alpha)
                uplow[0,1] = np.sin(alpha)
                uplow[1,0] = -np.sin(alpha)
                uplow[1,1] = np.cos(alpha)
                upper = np.zeros([2,len(xU)])
                lower = np.zeros([2,len(xU)])
                for i in range(len(xU)):
                    upper[0,i] = uplow[0,0] * xU[i] + uplow[0,1] * yU[i]
                    upper[1,i] = uplow[1,0] * xU[i] + uplow[1,1] * yU[i]
                    lower[0,i] = uplow[0,0] * xL[i] + uplow[0,1] * yL[i]
                    lower[1,i] = uplow[1,0] * xL[i] + uplow[1,1] * yL[i]
                xU = upper[0,:]
                yU = upper[1,:]
                xL = lower[0,:]
                yL = lower[1,:]
                
                return xU, xL, yU, yL
            
                
        x = np.linspace(0,c,101)
        if plot==True:
            #x = np.linspace(0,c,1001)
            if md > 0.2 and pd > 0.2:
                xu, xl, yu, yl = xyUL(x,md,pd,lst_digits, c)
                plt.plot(xu,yu, 'b')
                plt.plot(xl,yl, 'b')
                plt.title('NACA'+Naca+' \nAoA =' +str(AoA) )
                plt.xlabel('x')
                plt.ylabel('y')
            else:
                if AoA == None:
                    plt.plot(x,yt(lst_digits, c, x),'b') 
                    plt.plot(x,-yt(lst_digits, c, x),'b')
                    plt.title('NACA'+Naca+' \nAoA =' +str(AoA) )
                    plt.xlabel('x')
                    plt.ylabel('y')
                else:
                    xU, xL, yU, yL = ytAoA(lst_digits, c, x)
                    plt.plot(xU, yU ,'b') 
                    plt.plot(xL, yL,'b')
                    plt.title('NACA'+Naca+' \nAoA =' +str(AoA) )
                    plt.xlabel('x')
                    plt.ylabel('y')

            plt.axis('equal')
            plt.savefig(Naca+'_'+str(alpha)+'.png')
            plt.show()

        if meshdict==True:
            #globals
            global pts1
            global pts2
            global pts3
            global pts4
            global pts5
            global pts6
            global pts7
            global pts8
            # Mesh dimensions
            scale = 1            # Scaling factor
            H = 8                # *Half* height of channel
            #W = 0.5              # *Half* depth of foil (z-direction)
            W = 1
            D = 16               # Length of downstream section
            
            # Mesh resolution parameters
            Ni = 400             # Number of interpolation points along the foil
            Nx = 250             # Number of mesh cells along the foil
            ND = 150             # Number of cells in the downstream direction
            NT = 100             # Number of cells the transverse direction
            NW = 1               # Number of cells in the y-direction (along the foil axis)

            # Expansion rates
            ExpT = 500           # Expansion rate in transverse direction
            ExpD = 100           # Expansion rate in the downstream direction
            ExpArc = 50          # Expansion rate along the inlet arc
            
            #beta = np.linspace(0,np.pi,Ni)
            #x = c*(0.5*(1-np.cos(beta))
            sym_AoA = False
            if md < 0.2 and pd < 0.2:
                if AoA != None:
                    sym_AoA = True
                else:
                    sym_AoA = False
                    
            if md + pd > 0.2 or sym_AoA == True:
                if sym_AoA ==True:
                    xu, xl, yu, yl = ytAoA(lst_digits, c, x)
                else:
                    xu, xl, yu, yl = xyUL(x,md,pd,lst_digits, c)
                #x for max camber
                if sym_AoA ==True:
                    yt_ = 0
                    cmax_idx = 0
                    for i in range(len(x)):
                        yt_n = yt(lst_digits,c,x[i])
                        if yt_n > yt_:
                            yt_ = yt_n
                            cmax_idx = i
                else:
                    yc_ = 0
                    cmax_idx = 0
                    for i in range(len(x)):
                        yc_n = yc(md,pd,c,x[i])
                        if yc_n > yc_:
                            yc_ = yc_n
                            cmax_idx = i
                        
                NoseX = (-H + xu[cmax_idx])*np.cos(alpha)
                NoseY = -(-H + xu[cmax_idx])*np.sin(alpha)
                fvertices = np.zeros([12,3])
            
                #Front            
                fvertices[0,:] =  [NoseX,         NoseY,          W]
                fvertices[1,:] =  [xu[cmax_idx],  H,              W]
                fvertices[2,:] =  [xu[-1],        H,              W]
                fvertices[3,:] =  [D,             H,              W]
                fvertices[4,:] =  [0,             0,              W]
                fvertices[5,:] =  [xu[cmax_idx],  yu[cmax_idx],   W]
                fvertices[6,:] =  [xl[cmax_idx],  yl[cmax_idx],   W]
                fvertices[7,:] =  [xu[-1],        yu[-1],         W]
                fvertices[8,:] =  [D,             yu[-1],         W]
                fvertices[9,:] =  [xl[cmax_idx],  -H,             W]
                fvertices[10,:] = [xu[-1],        -H,             W]
                fvertices[11,:] = [D,             -H,             W] 
               
                #Back
                bvertices = cp.copy(fvertices)
                bvertices[:,2] = 0
            
                #Edge 4-5, 16-17
                pts1 = [xu[1:cmax_idx], yu[1:cmax_idx], W*np.ones(len(xu[1:cmax_idx]))]
                pts5 = [xu[1:cmax_idx], yu[1:cmax_idx], np.zeros(len(xu[1:cmax_idx]))]

                #Edge 5-7, 17-19
                pts2 = [xu[cmax_idx+1:-1], yu[cmax_idx+1:-1], W*np.ones(len(xu[cmax_idx+1:-1]))]
                pts6 = [xu[cmax_idx+1:-1], yu[cmax_idx+1:-1], np.zeros(len(xu[cmax_idx+1:-1]))]

                #Edge 4-6, 16-18
                pts3 = [xl[1:cmax_idx], yl[1:cmax_idx], W*np.ones(len(xu[1:cmax_idx]))]
                pts7 = [xl[1:cmax_idx], yl[1:cmax_idx], np.zeros(len(xu[1:cmax_idx]))]
                
                #Edge 6-7, 18-19
                pts4 = [xl[cmax_idx+1:-1], yl[cmax_idx+1:-1], W*np.ones(len(xu[cmax_idx+1:-1]))]
                pts8 = [xl[cmax_idx+1:-1], yl[cmax_idx+1:-1], np.zeros(len(xu[cmax_idx+1:-1]))]
                
                #Edge 0-1, 12-13
                pts9 = [-H*np.cos(np.pi/4.)+xu[cmax_idx], H*np.sin(np.pi/4.), W]
                pts11 = [-H*np.cos(np.pi/4.)+xu[cmax_idx], H*np.sin(np.pi/4.), 0]
                
                #Edge 0-9, 12-21
                pts10 = [-H*np.cos(np.pi/4.)+xu[cmax_idx], -H*np.sin(np.pi/4.), W]
                pts12 = [-H*np.cos(np.pi/4.)+xu[cmax_idx], -H*np.sin(np.pi/4.), 0]
                
            else: #elif AoA = None and md/pd == 0:
                y = [yt(lst_digits, c, x[i]) for i in range(len(x))]
                x = np.asarray(x)
                y = np.asarray(y)
                yt_ = 0
                cmax_idx = 0
                for i in range(len(x)):
                    yt_n = yt(lst_digits,c,x[i])
                    if yt_n > yt_:
                        yt_ = yt_n
                        cmax_idx = i
                        
                NoseX = (-H + x[cmax_idx])*np.cos(alpha)
                NoseY = -(-H + x[cmax_idx])*np.sin(alpha)
                fvertices = np.zeros([12,3])
                
                #Front
                fvertices[0,:] =  [NoseX,         NoseY,          W]
                fvertices[1,:] =  [x[cmax_idx],   H,              W]
                fvertices[2,:] =  [x[-1],         H,              W]
                fvertices[3,:] =  [D,             H,              W]
                fvertices[4,:] =  [0,             0,              W]
                fvertices[5,:] =  [x[cmax_idx],   y[cmax_idx],    W]
                fvertices[6,:] =  [x[cmax_idx],   -y[cmax_idx],   W]
                fvertices[7,:] =  [x[-1],         y[-1],          W]
                fvertices[8,:] =  [D,             y[-1],          W]
                fvertices[9,:] =  [x[cmax_idx],   -H,             W]
                fvertices[10,:] = [x[-1],         -H,             W]
                fvertices[11,:] = [D,             -H,             W]          
                
                #Back
                bvertices = cp.copy(fvertices)
                bvertices[:,2] = 0
                
                #Edge 4-5, 16-17
                pts1 = [x[1:cmax_idx], y[1:cmax_idx], W*np.ones(len(x[1:cmax_idx]))]
                pts5 = [x[1:cmax_idx], y[1:cmax_idx], np.zeros(len(x[1:cmax_idx]))]

                #Edge 5-7, 17-19
                pts2 = [x[cmax_idx+1:-1], y[cmax_idx+1:-1], W*np.ones(len(x[cmax_idx+1:-1]))]
                pts6 = [x[cmax_idx+1:-1], y[cmax_idx+1:-1], np.zeros(len(x[cmax_idx+1:-1]))]

                #Edge 4-6, 16-18
                pts3 = [x[1:cmax_idx], -y[1:cmax_idx], W*np.ones(len(x[1:cmax_idx]))]
                pts7 = [x[1:cmax_idx], -y[1:cmax_idx], np.zeros(len(x[1:cmax_idx]))]
                
                #Edge 6-7, 18-19
                pts4 = [x[cmax_idx+1:-1], -y[cmax_idx+1:-1], W*np.ones(len(x[cmax_idx+1:-1]))]
                pts8 = [x[cmax_idx+1:-1], -y[cmax_idx+1:-1], np.zeros(len(x[cmax_idx+1:-1]))]
                
                #Edge 0-1, 12-13
                pts9 = [-H*np.cos(np.pi/4.)+x[cmax_idx], H*np.sin(np.pi/4.), W]
                pts11 = [-H*np.cos(np.pi/4.)+x[cmax_idx], H*np.sin(np.pi/4.), 0]
                
                #Edge 0-9, 12-21
                pts10 = [-H*np.cos(np.pi/4.)+x[cmax_idx], -H*np.sin(np.pi/4.), W]
                pts12 = [-H*np.cos(np.pi/4.)+x[cmax_idx], -H*np.sin(np.pi/4.), 0]
            
            edges = [4, 5, 5, 7, 4, 6, 6, 7, 16, 17, 17, 19, 16, 18, 18, 19]
            #Write to blocMeshDict
            bd = open('constant/polyMesh/blockMeshDict','w')
            bd.write('''/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.4.0                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * / \n''')
            bd.write('\nconvertToMeters 1;\n')
            bd.write('\nvertices\n')
            bd.write('(\n    //front\n')
            for i in range(12):
                bd.write('    (%s %s %s)\n' %(fvertices[i,0], fvertices[i,1], fvertices[i,2]) )
            bd.write('\n    //back\n')
            for j in range(12):
                bd.write('    (%s %s %s)\n' %(bvertices[j,0], bvertices[j,1], bvertices[j,2]) )
            bd.write(');\n')
            bd.write('\nblocks\n')
            bd.write('(\n')
            bd.write('    hex (4 5 1 0  16 17 13 12)      (100 100 1) edgeGrading (1 0.02 0.02 1 500 500 500 500 1 1 1 1)\n') #integers 
            bd.write('    hex (5 7 2 1  17 19 14 13)      (150 100 1) simpleGrading (1 500 1)\n')
            bd.write('    hex (7 8 3 2  19 20 15 14)      (150 100 1) simpleGrading (100 500 1)\n')
            bd.write('    hex (16 18 21 12  4 6 9 0)      (100 100 1) edgeGrading (1 0.02 0.02 1 500 500 500 500 1 1 1 1)\n')
            bd.write('    hex (18 19 22 21  6 7 10 9)     (150 100 1) simpleGrading (1 500 1)\n')
            bd.write('    hex (19 20 23 22  7 8 11 10)    (150 100 1) simpleGrading (100 500 1)\n')
            bd.write(');\n')
            bd.write('\nedges\n')
            bd.write('(\n')
            #print len(pts6[0]), len(pts6[2])
            #print pts1[0][1]#], pts1[0,1], pts1[0,2]
            for i in range(8):
                pts = globals()['pts'+str(i+1)]
                bd.write('    spline %s %s\n' %(edges[2*i], edges[2*(i+1)-1]))
                bd.write('        (\n')
                for j in range(len(pts[0])):
                    bd.write('            (%s %s %s)\n' %(pts[0][j], pts[1][j], pts[2][j]))
                bd.write('        )\n')
                #bd.write('\n')
            #bd.write('        (\n')
            bd.write('    arc 0 1 (%s %s %s) \n'  %(pts9[0],pts9[1],pts9[2]))
            bd.write('    arc 0 9 (%s %s %s) \n' %(pts10[0],pts10[1],pts10[2]))
            bd.write('    arc 12 13 (%s %s %s) \n' %(pts11[0],pts11[1],pts11[2]))
            bd.write('    arc 12 21 (%s %s %s) \n' %(pts12[0],pts12[1],pts12[2]))
            bd.write(');\n')
            bd.write('\n')
            bd.write('''boundary 
( 
    inlet 
    { 
        type patch; 
        faces 
        ( 
            (1 0 12 13) 
            (0 9 21 12) 
        ); 
    } 

    outlet 
    { 
        type patch; 
        faces 
        ( 
            (11 8 20 23) 
            (8 3 15 20) 
        ); 
    } 

    topAndBottom 
    { 
        type patch; 
        faces 
        ( 
            (3 2 14 15) 
            (2 1 13 14) 
            (9 10 22 21) 
            (10 11 23 22) 
        ); 
    } 

    airfoil 
    { 
        type wall; 
        faces 
        ( 
            (5 4 16 17) 
            (7 5 17 19) 
            (4 6 18 16) 
            (6 7 19 18) 
        ); 
    } 
); 
 
mergePatchPairs 
( 
); 
 
// ************************************************************************* //''')
            bd.close()      
          
naca(plot=True, meshdict=False, AoA = 0)
                



