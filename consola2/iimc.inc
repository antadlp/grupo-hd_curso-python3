c    ** global parameters **
        integer     n
        common /global/ n

c    ** common /cofig/ Positions
        doubleprecision rx(2000), ry(2000), rz(2000)
        doubleprecision gdp(900)
        common /cofig/  rx, ry, rz, gdp

c    ** common /params/ Simulation control parameters
        doubleprecision rcut,drmax,dboxmx,box,boxz,fsxy(1000),fsz(1000)
        integer         nstep, iratio, nblock, nrad,iratb,epsw
        character       disnam*20, rundat*20,dennam*20
        common /params/ rcut, drmax, dboxmx, box,boxz,fsxy,fsz,
     :                  nstep, iratio, nblock, nrad,iratb,epsw,
     :                  disnam, rundat,dennam

c    ** common /values/ Current values of energy etc.
        doubleprecision eng, pot, vir, temp, press,presure
        doubleprecision maxrad,ro,virx,viry,virz,pressx,pressy,pressz
        common /values/ eng, pot, vir, temp, press,presure,
     :                  maxrad,ro,virx,viry,virz,pressx,pressy,pressz
 
 
c    ** common /blkavs/ cycles average values of energy etc.
        doubleprecision acp, acv, runpot, runvir,detz,np
        doubleprecision runrad(900),runpden(900),pden(900)
        common /blkavs/ acp, acv, runpot, runvir,detz,np,
     :                  runrad,runpden,pden
 
 
c    ** common /comran/ random numbers
        integer         ix1, ix2, ix3
        doubleprecision rando(97)
        common /comran/ rando, ix1, ix2, ix3 
