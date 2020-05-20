c    ** global parameters **
        integer     n
        common /global/ n

c    ** common /cofig/ Positions, velocities and forces
        doubleprecision rx(3000), ry(3000), rz(3000)
        doubleprecision vx(3000), vy(3000), vz(3000)
        doubleprecision fx(3000), fy(3000), fz(3000)
        common /cofig/ rx, ry, rz, vx, vy, vz, fx, fy, fz


c    ** common /params/ Simulation control parameters
        doubleprecision dt, rcut, req_temp
	doubleprecision siga, epsa, rcuta, sigb, epsb, rcutb
	doubleprecision boxx, boxy, boxz,epsw
        integer         nblock, nstep, nrad
	logical         fix_temp
        character       disnam*20, rundat*20
	common /params/ dt, rcut, req_temp, 
     :                  siga, epsa, rcuta, sigb, epsb, rcutb,
     :                  boxx, boxy, boxz,epsw,
     :                  nblock, nstep, fix_temp,
     :                  disnam, rundat, nrad

c    ** common /values/ Current values of energy etc.
        doubleprecision eng, kin, pot, vir, temp, press, ro, maxrad
	common /values/ eng, kin, pot, vir, temp, press, ro, maxrad
	

c    ** common /runavs/ Run average values of energy etc.
        doubleprecision runeng, runkin, runpot, runvir
	common /runavs/ runeng, runkin, runpot, runvir
	
c    ** common /blkavs/ Block average values of energy etc.
        doubleprecision blkeng,blkkin,blkpot,blkvir,runrad(200),gdp(200)
	common /blkavs/ blkeng,blkkin,blkpot,blkvir,runrad, gdp
