C    *******************************************************************
        program md
	implicit none
	include 'md.inc'
C    *******************************************************************

C    *******************************************************************

        call config
        call start
	call run
	call finish
	
	stop
	end
	
	
C    *******************************************************************
        subroutine start
	implicit none
	include 'md.inc'
C    *******************************************************************
C    ** Local variables **
        integer         i
	doubleprecision dens
C    *******************************************************************

	print *, 'md version 0.1'
	print *, 'Molecular Dynamics of Lennard-Jones atoms'
	print *, 'Periodic Boundaries in all directions'
	
        print *, 'Enter file to keep last running data'
        read  (*,'(a)') rundat
        write (*,'('' file     '',a)') rundat
        print *,'Enter file to keep g(r) data'
        read  (*,'(a)') disnam
        write (*,'('' file      '',a)') disnam
	print *, 'Enter number of blocks'
	read  *, nblock
	print *, 'Number of blocks = ', nblock
	print *, 'Enter number of steps per block'
	read  *, nstep
	print *, 'Number of steps  = ', nstep
	print *, 'Enter timestep'
	read  *, dt
	print *, 'Timestep         = ', dt
	print *, 'Enter constant temperature option'
	read  *, fix_temp
	print *, 'Constant temperature option = ', fix_temp
	print *, 'Enter required temperature'
	read  *, req_temp
	print *, 'Required temperature = ', req_temp
	print *, 'Enter potential cutoff distance'
	read  *, rcut
	
	dens = dble(n)/(boxx*boxy*boxz)
        ro   = boxx/dble(200)
        nrad = int(anint(boxx/(2.0d0*ro)))

	call force
	call kinetic
	
	eng   = kin + pot
	temp  = 2.0d0*kin/3.0d0
 	press = dens * ( temp + vir )
	
	print *, 'Initial values'
	print *, 'Energy      = ', eng
	print *, 'Kinetic     = ', kin
	print *, 'Potential   = ', pot
	print *, 'Temperature = ', temp
	print *, 'Pressure    = ', press
	print *, 'Density     = ', dens
	return
	end
	
	
C    *******************************************************************
        subroutine run
	implicit none
	include 'md.inc'
C    *******************************************************************
C    ** Local variables **
        integer         j, block, step
	doubleprecision dens, asqpot, asqprxy, flp, flpx
C    *******************************************************************

C    ** Zero run accumulators **
        runeng = 0.0d0
        runkin = 0.0d0
        runpot = 0.0d0
        runvir = 0.0d0
        asqpot = 0.0d0
        asqprxy= 0.0d0
        do j=1, nrad
         runrad(j) = 0.0d0
        enddo


        write(*,'(1x,''Block'',''    Energy'',''   Kinetic'',
     :            '' Potential'',''      Temp'',''     Press'')')
     	
        do block = 1, nblock

C       ** Zero run accumulators **
           blkeng = 0.0d0
           blkkin = 0.0d0
           blkpot = 0.0d0
           blkvir = 0.0d0

	   do step = 1, nstep
	      call movvel
	      call movpos
	      call force
              eng = kin + pot
	      runeng = runeng + eng
	      runkin = runkin + kin
	      runpot = runpot + pot
	      runvir = runvir + vir
	      blkeng = blkeng + eng
	      blkkin = blkkin + kin
	      blkpot = blkpot + pot
	      blkvir = blkvir + vir
              asqpot = asqpot + pot*pot

	      call radial 
              do j=1, nrad
               runrad(j) = runrad(j)+gdp(j)
              enddo 
	   enddo


C       ** Normalize block accumulators **
           blkeng = blkeng / dble(nstep)
           blkkin = blkkin / dble(nstep)
           blkpot = blkpot / dble(nstep)
           blkvir = blkvir / dble(nstep) 

	   dens   = dble(n)/(boxx*boxy*boxz)
	   temp   = 2.0d0*blkkin/3.0d0
	   press  = dens * ( temp + blkvir )
           asqprxy= asqprxy + press*press

	   write(*,'(1x,i5,5f10.4)') 
     :          block, blkeng, blkkin, blkpot, temp, press
	enddo

        do j=1, nrad
         runrad(j) = runrad(j)/dble(nstep*nblock)
        enddo 

	
C    ** Normalize run accumulators **
        runeng = runeng / dble(nstep*nblock)
        runkin = runkin / dble(nstep*nblock)
        runpot = runpot / dble(nstep*nblock)
        runvir = runvir / dble(nstep*nblock)
        asqpot = (asqpot/dble(nstep*nblock))  - runpot**2

	dens   = dble(n)/(boxx*boxy*boxz)
	temp   = 2.0d0*runkin/3.0d0
	press  = dens * ( temp + runvir )
	asqprxy= (asqprxy/nblock) - press**2
        if (asqpot .gt. 0.0d0)  flp  = sqrt(asqpot)
        if (asqprxy .gt. 0.0d0) flpx = sqrt(asqprxy)


	write(*,'(1x,55(''-''))')
	write(*,'(1x,5x,5f10.4)') 
     :       runeng, runkin, runpot, temp, press
	write(*,'(1x,55(''-''))')
        write(*,'(1x,'' ** fluctation values ** '')')
        write(*,'(1x,5x,8f10.4)') flp, flpx

	
	return
	end


C    *******************************************************************
        subroutine finish
	implicit none
	include 'md.inc'
C    *******************************************************************
C    ** Local variables **
	integer         i
	doubleprecision dens
C    *******************************************************************
	
        dens = dble(n)/(boxx*boxy*boxz)
        eng   = kin + pot
        temp = 2.0d0*kin/3.0d0
        press = dens * ( temp + vir )
         
	print *, 'Final values'
	print *, 'Energy      = ', eng
	print *, 'Kinetic     = ', kin
	print *, 'Potential   = ', pot
	print *, 'Temperature = ', temp
	print *, 'Pressure    = ', press

         do i=1, n
          rz(i) = rz(i) + boxx/2.0d0
         enddo
	
	open(unit=3,file=rundat,status='unknown')
	write(3,'(1x,i5)') n
	write(3,'(1x,3f10.5)') boxx, boxy, boxz
        do i = 1, n
	   write(3,'(1x,6f10.5)') rx(i),ry(i),rz(i),vx(i),vy(i),vz(i)
	enddo
	close(unit=3)

        open(unit=3,file=disnam,status='unknown')
C        write(3,'(1x,''   k'',''   r'',''   g(2)'')')
        do i=1, nrad
          write(3,'(1x,i5,2f10.5)') i,i*ro-ro/2d0,runrad(i)
        enddo
        close(unit=3) 


        return
	end


C    *******************************************************************
        subroutine movvel
	implicit none
	include 'md.inc'
C    *******************************************************************
C    ** Local variables **
	integer         i
	doubleprecision dt2, factor
C    *******************************************************************
	
	dt2 = dt / 2.0d0
	kin = 0.0d0

C    ** First half of velocity move **
C    ** Calculate kinetic energy    **
	do i = 1, n
	   vx(i) = vx(i) + dt2*fx(i)
	   vy(i) = vy(i) + dt2*fy(i)
	   vz(i) = vz(i) + dt2*fz(i)
	   kin = kin + vx(i)**2 + vy(i)**2 + vz(i)**2
	enddo
	kin = kin / 2.0d0
	kin = kin / dble(n)
	temp = 2.0d0 * kin / 3.0d0

C    ** Scale velocities if required **	
	if ( fix_temp ) then
	   factor = sqrt(req_temp/temp)
	   kin = 0.0d0
	   do i = 1, n
	      vx(i) = vx(i) * factor
	      vy(i) = vy(i) * factor
	      vz(i) = vz(i) * factor
	      kin = kin + vx(i)**2 + vy(i)**2 + vz(i)**2
	   enddo
	   kin = kin / 2.0d0
	   kin = kin / dble(n)
	   temp = 2.0d0 * kin / 3.0d0
	endif

C    ** Second half of velocity move **	
	do i = 1, n
	   vx(i) = vx(i) + dt2*fx(i)
	   vy(i) = vy(i) + dt2*fy(i)
	   vz(i) = vz(i) + dt2*fz(i)
	enddo

        return
	end

	
C    *******************************************************************
        subroutine movpos
	implicit none
	include 'md.inc'
C    *******************************************************************
C    ** Local variables **
	integer i
C    *******************************************************************

        do i = 1, n
	   rx(i) = rx(i) + dt*vx(i)
	   ry(i) = ry(i) + dt*vy(i)
	   rz(i) = rz(i) + dt*vz(i)
	   rx(i) = rx(i) - anint(rx(i)/boxx)*boxx
	   ry(i) = ry(i) - anint(ry(i)/boxy)*boxy
           rz(i) = rz(i) - anint(rz(i)/boxz)*boxz
	enddo
	
	return
	end


C    *******************************************************************
        subroutine force
	implicit none
	include 'md.inc'
C    *******************************************************************
C    ** Local variables **
	integer         i, j
	doubleprecision rxi, ryi, rzi, fxi, fyi, fzi
	doubleprecision rxij, ryij, rzij, rijsq, rcutsq, vcut
	doubleprecision r2ij, r6ij, r12ij, fij, vij
	doubleprecision fxij, fyij, fzij
C    *******************************************************************

        rcutsq = rcut*rcut
	r2ij   = 1.0d0/rcutsq
	r6ij   = r2ij*r2ij*r2ij
	r12ij  = r6ij*r6ij
	vcut   = r12ij - r6ij
	
C    ** Zero energies and forces **
        pot = 0.0d0
	vir = 0.0d0
	do i = 1, n
	   fx(i) = 0.0d0
	   fy(i) = 0.0d0
	   fz(i) = 0.0d0
	enddo

C    ** Outer loop over i **
        do i = 1, n-1
	   rxi = rx(i)
	   ryi = ry(i)
	   rzi = rz(i)
	   fxi = fx(i)
	   fyi = fy(i)
	   fzi = fz(i)
	   
C       ** Inner loop over j **
           do j = i+1, n
	      rxij = rxi - rx(j)
	      ryij = ryi - ry(j)
	      rzij = rzi - rz(j)
	      rxij = rxij - anint(rxij/boxx)*boxx
	      ryij = ryij - anint(ryij/boxy)*boxy
              rzij = rzij - anint(rzij/boxz)*boxz
              rijsq= rxij*rxij + ryij*ryij + rzij*rzij
	      if ( rijsq .lt. rcutsq ) then
	       r2ij  = 1.0d0/rijsq
	       r6ij  = r2ij*r2ij*r2ij
	       r12ij = r6ij*r6ij
	       vij   = r12ij - r6ij
	       fij   = (vij + r12ij)*r2ij
	       vij   = vij - vcut
	       pot   = pot + vij
	       fxij  = fij*rxij
	       fyij  = fij*ryij
	       fzij  = fij*rzij
	       vir   = vir + rxij*fxij + ryij*fyij + rzij*fzij
	       fxi   = fxi + fxij
	       fyi   = fyi + fyij
	       fzi   = fzi + fzij
	       fx(j) = fx(j) - fxij
	       fy(j) = fy(j) - fyij
	       fz(j) = fz(j) - fzij    
	      endif	 
	   enddo
           fx(i) = fxi
	   fy(i) = fyi
           fz(i) = fzi
	enddo

C    ** Incorporate correct numerical factors **
        pot = pot * 4.0d0
	vir = vir * 24.0d0/3.0d0
        do i = 1, n
	   fx(i) = fx(i)*24.0d0
	   fy(i) = fy(i)*24.0d0
	   fz(i) = fz(i)*24.0d0
	enddo

        pot = pot/dble(n)
	vir = vir/dble(n)


        return
	end




	
C    *******************************************************************
        subroutine kinetic
	implicit none
	include 'md.inc'
C    *******************************************************************
C    ** Local variables **
	integer         i
	doubleprecision dt2, vxi, vyi, vzi
C    *******************************************************************
	
	dt2 = dt / 2.0d0
	kin = 0.0d0

C    ** First half of velocity move **
C    ** Calculate kinetic energy    **
	do i = 1, n
	   vxi = vx(i) + dt2*fx(i)
	   vyi = vy(i) + dt2*fy(i)
	   vzi = vz(i) + dt2*fz(i)
	   kin = kin + vxi**2 + vyi**2 + vzi**2
	enddo
	kin = kin / 2.0d0
	kin = kin / dble(n)
	return
	end


C    *******************************************************************
        subroutine  radial
        implicit none
        include 'md.inc'
C    *******************************************************************
        integer         i, k, m, j
        doubleprecision xij, yij, zij, rsq, pi, rij
        doubleprecision cte, rlower, rupper,nideal
        parameter       (pi=3.14159265)
C    *******************************************************************

        maxrad = boxx/2.0d0
        ro     = maxrad/dble(nrad)

        do j=1, nrad
          gdp(j) = 0.0d0
        enddo

        do i=1, n-1
         do j = i+1, n
             xij = rx(i)-rx(j)
             yij = ry(i)-ry(j)
             zij = rz(i)-rz(j)
             xij = xij - boxx*anint(xij/boxx)
             yij = yij - boxy*anint(yij/boxy)
             zij = zij - boxz*anint(zij/boxz)
             rsq = xij*xij + yij*yij + zij*zij
             rij = sqrt(rsq)
             k   = int(rij/ro)+1
             if (k .le. nrad) then
              gdp(k) = gdp(k)+2
             endif
          enddo
         enddo
         cte = 4.0d0*pi*(dble(n)/(boxx*boxy*boxz))/3.0d0
         do m=1, nrad
          rlower = dble(m-1)*ro
          rupper = rlower+ro
          nideal = cte*(rupper**3 - rlower**3)
          gdp(m) = gdp(m)/(dble(n)*nideal)
         enddo
                   
        return
        end 


C    *******************************************************************
        subroutine  config
	implicit none
	include 'md.inc'
C    *******************************************************************
C    ** Local variables **
        integer         i
        character       name*20
C    *******************************************************************

        print *,'enter file'
        read  (*,'(a)') name
        write(*,'('' file       '',a)') name

C    *******      read initial configuration      *******
        open (unit=3, file= name, status= 'old')
         read(3,*) n
         read(3,*)  boxx, boxy, boxz
          do i=1, n
C             read(3,*) rx(i),ry(i),rz(i)
             read(3,*) rx(i),ry(i),rz(i),vx(i),vy(i),vz(i)
          enddo
        close (unit=3)

         do i=1, n
          rz(i) = rz(i) - boxz/2.0d0
         enddo
    
         return 
         end

