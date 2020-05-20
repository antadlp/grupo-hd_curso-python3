C    *******************************************************************
        program  addvel
	implicit none
C    *******************************************************************
        integer         i,n,n1,NATOM,ia
        doubleprecision rx(5000),ry(5000),rz(5000),ex
        doubleprecision box,hlen
	doubleprecision TEMP,GAUSS,DUMMY,RTEMP,XV,YV,ZV,V2,TTRAN
	doubleprecision VX(2000),VY(2000),VZ(2000)
	doubleprecision VXT(2000),VYT(2000),VZT(2000)
	doubleprecision MOPX,MOPY,MOPZ,fact
        character       name*20,rundat*20
C    *******************************************************************

        print *,'enter file'
        read  (*,'(a)') name
        write(*,'('' file       '',a)') name
        print *, 'Enter file to keep last data'
        read  (*,'(a)') rundat
        print *, 'Enter number of data'
        read *, n
        print *,'enter number of NON-frozen particles'
        READ  *, NATOM
        PRINT *,'TEMPERATURE'
        READ  *,TEMP

      RTEMP = DSQRT ( TEMP )
      DO I=1,NATOM
        VX(I) = RTEMP*GAUSS ( DUMMY )
        VY(I) = RTEMP*GAUSS ( DUMMY )
        VZ(I) = RTEMP*GAUSS ( DUMMY )
      ENDDO

      V2 = 0.0
      DO I=1,NATOM
       V2 = V2 + VX(I)**2 + VY(I)**2 + VZ(I)**2
      ENDDO
      TTRAN = 1.0D0/3.0D0 * V2/DFLOAT(NATOM)
      fact = sqrt(temp/TTRAN)
      DO I=1,NATOM
         VX(I) = fact*VX(I)
         VY(I) = fact*VY(I)
         VZ(I) = fact*VZ(I)
      ENDDO
      XV = 0.0
      YV = 0.0
      ZV = 0.0
      DO I=1,NATOM
         VX(I) = VX(I) - XV
         VY(I) = VY(I) - YV
         VZ(I) = VZ(I) - ZV
      ENDDO

      DO I=1,NATOM
         XV = XV + VX(I)
         YV = YV + VY(I)
         ZV = ZV + VZ(I)
      ENDDO


      XV = XV/DFLOAT(NATOM)
      YV = YV/DFLOAT(NATOM)
      ZV = ZV/DFLOAT(NATOM)
c
c  set total momentum to zero
c
      DO I=1,NATOM
         VX(I) = VX(I) - XV
         VY(I) = VY(I) - YV
         VZ(I) = VZ(I) - ZV
      ENDDO
      V2 = 0.0
c
c  'kinetic energies'
c
      DO I=1,NATOM
        V2 = V2 + VX(I)**2 + VY(I)**2 + VZ(I)**2
      ENDDO
c
c  resulting temperatures
c
      TTRAN = 1.0D0/3.0D0 * V2/DFLOAT(NATOM)
       MOPX = 0.0
       MOPY = 0.0
       MOPZ = 0.0
      DO I=1,NATOM
       MOPX = MOPX + VX(I)
       MOPY = MOPY + VY(I)
       MOPZ = MOPZ + VZ(I)
      ENDDO
      PRINT *,'TEMP0=',TTRAN
      PRINT *,'PX0,PY0,PZ0=',MOPX,MOPY,MOPZ


         open (unit=3, file= name, status= 'old')
          read(3,*) n
          read(3,*) box,box,hlen
          do i=1, n
C            read(3,*) rx(i),ry(i),rz(i),ex,ex,ex
            read(3,*) rx(i),ry(i),rz(i)
          enddo
        close (unit=3)
C        if (natom .gt. n) stop 'wrong number of particle'
C        if (natom .ne. n-n1+1) stop 'wrong number of particle n1'

        ex=0.0d0
        ia=1
        do i = 1, n-natom
           vxt(ia) = ex
           vyt(ia) = ex
           vzt(ia) = ex
           ia=ia+1
        enddo
        do i = 1, natom
           vxt(ia) = vx(i)
           vyt(ia) = vy(i)
           vzt(ia) = vz(i)
           ia=ia+1
        enddo
        print *,'ia=',ia-1
        if (ia-1 .ne. n) stop 'wrong number of particle n'

      V2 = 0.0
c  'kinetic energies'
c
      DO I=N-NATOM+1,N
        V2 = V2 + VXT(I)**2 + VYT(I)**2 + VZT(I)**2
      ENDDO
c
c  resulting temperatures
c
      TTRAN = 1.0D0/3.0D0 * V2/DFLOAT(NATOM)
       MOPX = 0.0
       MOPY = 0.0
       MOPZ = 0.0
      DO I=1,N
       MOPX = MOPX + VXT(I)
       MOPY = MOPY + VYT(I)
       MOPZ = MOPZ + VZT(I)
      ENDDO
      PRINT *,'TEMP=',TTRAN
      PRINT *,'PX,PY,PZ=',MOPX,MOPY,MOPZ

        open(unit=3,file=rundat,status='unknown')
         write(3,'(1x,2i14)') n
         write(3,'(1x,3f14.5)') box,box,hlen
         do i = 1, n
          if (i .le. n-natom) then
C          write(3,'(1x,3f10.4)') rx(i),ry(i),rz(i)
           write(3,'(1x,6f10.4)') rx(i),ry(i),rz(i),ex,ex,ex
          elseif (i .gt. n-natom) then
C          write(3,'(1x,3f10.4)') rx(i),ry(i),rz(i)
           write(3,'(1x,6f10.4)') rx(i),ry(i),rz(i),vxt(i),vyt(i),vzt(i)
          endif
         enddo
        close(unit=3)


         stop
         end


C    *******************************************************************
        doubleprecision FUNCTION GAUSS ( DUMMY )
	implicit none
        doubleprecision      A1, A3, A5, A7, A9
        PARAMETER ( A1 = 3.949846138, A3 = 0.252408784 )
        PARAMETER ( A5 = 0.076542912, A7 = 0.008355968 )
        PARAMETER ( A9 = 0.029899776                   )
        doubleprecision        SUM, R, R2
        doubleprecision        RANF0, DUMMY
        INTEGER     I
C    *******************************************************************

        SUM = 0.0

        DO 10 I = 1, 12

           SUM = SUM + RANF0 ( DUMMY )

10      CONTINUE

        R  = ( SUM - 6.0 ) / 4.0
        R2 = R * R

        GAUSS = (((( A9 * R2 + A7 ) * R2 + A5 ) * R2 + A3 ) * R2 +A1 )
     :          * R

        RETURN
        END

C    *******************************************************************
        doubleprecision FUNCTION RANF0 ( DUMMY )
	implicit none
        INTEGER     L, C, M
        PARAMETER ( L = 1029, C = 221591, M = 1048576 )
        INTEGER     SEED
        doubleprecision        DUMMY
        SAVE        SEED
        DATA        SEED / 0 /
C    *******************************************************************
        SEED = MOD ( SEED * L + C, M )
        RANF0 = REAL ( SEED ) / M
        RETURN
        END

