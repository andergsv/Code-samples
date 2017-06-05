'''
Who are you in the NBA?

A comparison of Norberg 2's indvidual averaged stat line 2015/16
 to every NBA player's stat line in the regular season 2015/16.
 
 Ten most compareable players based on the ratio of average and best 
 noted average in respective league. Smallest sum of differences deemd most
 compatible. 
'''
import numpy as np
import copy as cp

nor = open('nordberg', 'r')

nor.readline()

maxNOR = np.zeros(7) #MAX Nordberg 
norlen = 0
for line in nor:
    stat = line.split()[1:]
    for i in range(len(stat)):
        if float(stat[i]) > maxNOR[i]:
            maxNOR[i] = float(stat[i])*1.2
        else:
            pass
            
    norlen += 1
nor.close()    
nba = open('nba', 'r')

pN, rN, aN, sN, bN, FTN, PFN = 0, 0, 0, 0, 0, 0, 0 #MAX NBA 
length=1
nba.readline()
for line in nba:    
    stats = line.split()
    if len(stats) > 10:
        pp = float(stats[25]) #Poeng
        rr = float(stats[17]) #Returer
        aa = float(stats[18]) #assist
        ss = float(stats[20]) #steals
        bb = float(stats[21]) #blocks
        ft = float(stats[14]) #FT%
        pf = float(stats[22]) #PF
        if pp > pN:
            pN = pp
        if rr > rN:
            rN = rr
        if aa > aN:
            aN = aa
        if ss > sN:
            sN = ss
        if bb > bN:
            bN = bb
        if ft > FTN:
            FTN = ft
        if pf > PFN:
            PFN = pf
    length += 1

nba.close()
print pN, rN, aN, sN, bN, FTN, PFN
print maxNOR
nor = open('nordberg', 'r')

nor.readline()
def comp():
    for i in range(norlen):
        nordbergspiller = nor.readline().split()
        
        for j in range(7):
            nordbergspiller[j+1] = float(nordbergspiller[j+1])*1.2

        nba = open('nba', 'r')
        nba.readline()
        diff = 10000000
        diff1 = 10000000
        diff2 = 10000000
        diff3 = 10000000
        diff4 = 10000000
        diff5 = 10000000
        diff6 = 10000000
        diff7 = 10000000
        diff8 = 10000000
        diff9 = 10000000
        
        change1 = cp.copy(nordbergspiller)

        change1[1] = change1[1]/maxNOR[0] #Poeng
        change1[2] = change1[2]/maxNOR[1] #Returer
        change1[3] = change1[3]/maxNOR[2] #Assist
        change1[4] = change1[4]/maxNOR[3] #Steals
        change1[5] = change1[5]/maxNOR[4] #Blocks
        change1[6] = change1[6]/maxNOR[5] #FT%
        change1[7] = change1[7]/maxNOR[6] #PF
       
        
        name = 0
        name1 = 0
        name2 = 0
        name3 = 0
        name4 = 0
        
        stat = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        stat1 = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        stat2 = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        stat3 = [0,0, 0, 0 ,0 ,0 ,0 ,0] 
        stat4 = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        stat5 = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        stat6 = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        stat7 = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        stat8 = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        stat9 = [0,0, 0, 0 ,0 ,0 ,0 ,0]
        
        for k in range(int(length/2.)-1):
            nbasp = []
            
            navn = nba.readline().split()
            nbasp.append(navn[0] +' ' + navn[-1])
            
            stats = nba.readline().split()
            nbasp.append(float(stats[25])) #Poeng
            nbasp.append(float(stats[17])) #Returer
            nbasp.append(float(stats[18])) #assist
            nbasp.append(float(stats[20])) #steals
            nbasp.append(float(stats[21])) #blocks
            nbasp.append(float(stats[14])) #FT%
            nbasp.append(float(stats[22])) #PF
            
            change = cp.copy(nbasp)
            change[1] = change[1]/pN #Poeng
            change[2] = change[2]/rN #Returer
            change[3] = change[3]/aN #Assist
            change[4] = change[4]/sN #Steals
            change[5] = change[5]/bN #Blocks
            change[6] = change[6]/FTN #FT%
            change[7] = change[7]/PFN #PF
            

            
            ab = [abs(change[r+1] - change1[r+1]) for r in range(len(nbasp)-1)]
            su = sum(ab)
            
            if su < diff:
                diff1 = su
                diff = su
                
                name4 = name3
                name3 = name2
                name2 = name1
                name1 = name
                    
                stat9 = stat8
                stat8 = stat7
                stat7 = stat6
                stat6 = stat5
                stat5 = stat4    
                stat4 = stat3
                stat3 = stat2
                stat2 = stat1
                stat1 = stat
                stat = nbasp
                name = nbasp[0]
                
            elif su < diff1:
                diff2 = diff1
                diff1 = su
                stat9 = stat8
                stat8 = stat7
                stat7 = stat6
                stat6 = stat5
                stat5 = stat4    
                stat4 = stat3
                stat3 = stat2
                stat2 = stat1 
                stat1 = nbasp
            elif su < diff2:
                diff3 = diff2
                diff2 = su 
                stat9 = stat8
                stat8 = stat7
                stat7 = stat6
                stat6 = stat5
                stat5 = stat4    
                stat4 = stat3
                stat3 = stat2
                
                stat2 = nbasp
            elif su < diff3:
                diff4 = diff3
                diff3 = su 
                stat9 = stat8
                stat8 = stat7
                stat7 = stat6
                stat6 = stat5
                stat5 = stat4    
                stat4 = stat3
                 
                stat3 = nbasp
            elif su < diff4:
                diff5 = diff4
                diff4 = su
                stat9 = stat8
                stat8 = stat7
                stat7 = stat6
                stat6 = stat5
                stat5 = stat4    
                  
                stat4 = nbasp
            elif su < diff5:
                diff6 = diff5
                diff5 = su 
                stat9 = stat8
                stat8 = stat7
                stat7 = stat6
                stat6 = stat5
                
                stat5 = nbasp
            elif su < diff6:
                diff7 = diff6
                diff6 = su 
                stat9 = stat8
                stat8 = stat7
                stat7 = stat6
                
                stat6 = nbasp
            elif su < diff7:
                diff8 = diff7
                diff7 = su 
                stat9 = stat8
                stat8 = stat7
                
                stat7 = nbasp
            elif su < diff8:
                diff9 = diff8
                diff8 = su
                stat9 = stat8
                  
                stat8 = nbasp
            elif su < diff9:
                diff9 = su 
                 
                stat9 = nbasp   
                
                
                
        print '                    Navn           Poeng     Returer    Assist     Steals      Blocks     FT%          PF'
        print'--------------------------------------------------------------------------------------------------------------'
        print '{:>28} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(nordbergspiller[0],   \
                        nordbergspiller[1],nordbergspiller[2],nordbergspiller[3],nordbergspiller[4], \
                        nordbergspiller[5],nordbergspiller[6],nordbergspiller[7])
        print '1:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat[0], stat[1],   \
                        stat[2], stat[3], stat[4], stat[5], stat[6], stat[7])
        print '2:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat1[0], stat1[1], \
                        stat1[2], stat1[3], stat1[4], stat1[5], stat1[6], stat1[7])
        print '3:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat2[0], stat2[1], \
                        stat2[2], stat2[3], stat2[4], stat2[5], stat2[6], stat2[7])
        print '4:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat3[0], stat3[1], \
                        stat3[2], stat3[3], stat3[4], stat3[5], stat3[6], stat3[7])
        print '5:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat4[0], stat4[1], \
                        stat4[2], stat4[3], stat4[4], stat4[5], stat4[6], stat4[7])
        print '6:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat5[0], stat5[1], \
                        stat5[2], stat5[3], stat5[4], stat5[5], stat5[6], stat5[7])
        print '7:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat6[0], stat6[1], \
                        stat6[2], stat6[3], stat6[4], stat6[5], stat6[6], stat6[7])
        print '8:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat7[0], stat7[1], \
                        stat7[2], stat7[3], stat7[4], stat7[5], stat7[6], stat7[7])
        print '9:{:>26} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}'.format(stat8[0], stat8[1], \
                        stat8[2], stat8[3], stat8[4], stat8[5], stat8[6], stat8[7])
        print '10:{:>25} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format(stat9[0], stat9[1], \
                        stat9[2], stat9[3], stat9[4], stat9[5], stat9[6], stat9[7])
        
        nba.close()

comp()
