C        ******************************************************
         program azar
         implicit none
C        ******************************************************
         
C        ***** this program locate random particles in a box **
C        ***** given z-length  and density it determines  *****
C        ***** x and y length. ********************************
C        ****** variables *********

         integer         n, j, i, k, dr, num
         doubleprecision rx(25000), ry(25000), rz(25000)
         doubleprecision box, hlen, ranf, dummy, box2, dens
         doubleprecision rxi, ryi, rzi, rxij, ryij, rzij, rsij
C        ******************************************************

C        *** initialize random number with g05ccf (nag library)

         print *,'enter No. of particles'
         read  *, n
C         print *,'enter lenght in z'
C         read  *, hlen
         print *,'enter density'
         read  *, dens

         box  = (dble(n)/dens)**(1.0/3.0)
         box2 = box/2.0d0

C        *** put particles in a rectangular box ***
         do i=1, n
           rx(i)=box*ranf(dummy)-box2
           ry(i)=box*ranf(dummy)-box2
           rz(i)=box*ranf(dummy)
         enddo

C        *** check no-overlap among particles ****
         print *,'enter No. of iterations to check no-overlap'
         read  *, num
         do k=1, num
          dr=0
          do i=1, n-1
            rxi=rx(i)
            ryi=ry(i)
            rzi=rz(i)
            do j=i+1, n
10             rxij=rxi-rx(j)
               ryij=ryi-ry(j)
               rzij=rzi-rz(j)
               rsij=rxij*rxij+ryij*ryij+rzij*rzij
                 if (rsij .lt. 0.7d0) then
                   rx(j)=box*ranf(dummy)-box2
                   ry(j)=box*ranf(dummy)-box2
                   rz(j)=box*ranf(dummy)
                   dr=dr+1
                   goto 10
                 endif
            enddo
          enddo
          print *,'dr=',dr
         enddo

C        **** if exito then no-overlap ****
         if (dr .eq. 0.0) then
           print *,'exito, exito, exito, exito, exito, exito'
          else
           print *,'fracaso, fracaso, fracaso, fracaso'
         endif

         open (unit=10, file='azar.dat',status ='unknown')
         write (10,'(1x,i5)') n
         write (10,'(1x,3f10.6)') box, box, box
         do i=1, n
           write (10,'(1x,3f12.5)') rx(i), ry(i), rz(i)-box2
         enddo
         close (unit=10)

        stop
        end

C       ***************************************************************
         doubleprecision function ranf(dummy)
C       ***************************************************************
C       *** variables *****          
         integer l,c,m
         parameter (l=1029, c=22159, m=1048576)
         integer ran
         doubleprecision dummy
         save  ran
         data  ran /0/
C        **************************************************************
         
         ran = mod(ran*l+c,m)
         ranf= dble(ran)/m
         
          if (ranf .lt. 0.0) stop 'error in random number'
         
         return
         end       

