      subroutine read_ringin(nrings, rsize, ratoms,
     &                       maxrings, maxatoms)
c
c     Reads GICForge runtime file "ringin"
c
c     OUTPUT:
c       nrings   - number of rings
c       rsize(i) - size of ring i
c       ratoms(j,i) - atom indices of ring i (1-based)
c
c     INPUT:
c       maxrings - maximum number of rings allocated
c       maxatoms - maximum ring size allocated
c
c     Notes:
c     - Atom indices are 1-based (Fortran / XYZ convention)
c     - This routine performs no chemical interpretation
c

      integer nrings
      integer rsize(maxrings)
      integer ratoms(maxatoms, maxrings)
      integer maxrings, maxatoms

      integer iu, i, j, k
      character*256 line

c     --- Initialize ---
      nrings = 0
      do 10 i = 1, maxrings
         rsize(i) = 0
         do 20 j = 1, maxatoms
            ratoms(j,i) = 0
 20      continue
 10   continue

c     --- Open file ---
      iu = 98
      open(unit=iu, file='ringin', status='old', err=900)

c     --- Read file line by line ---
 100  continue
      read(iu,'(A)',end=800) line

c     Skip comments and empty lines
      if (line(1:1) .eq. '#') go to 100
      if (line .eq. ' ') go to 100

c     Read number of rings
      if (line(1:5) .eq. 'RINGS') then
         read(line(6:),*) nrings
         if (nrings .gt. maxrings) go to 910
         go to 100
      end if

c     Read ring definition
      if (line(1:4) .eq. 'RING') then
c        Expected format:
c        RING i SIZE n ATOMS a1 a2 ... an
         nr = 0
         ns = 0

         read(line,*) dummy, i, dummy, ns, dummy,
     &                 (ratoms(k,i), k=1,ns)

         rsize(i) = ns
         go to 100
      end if

      go to 100

c     --- Normal end ---
 800  continue
      close(iu)
      return

c     --- Errors ---
 900  continue
      nrings = 0
      return

 910  continue
      write(*,*) 'ERROR: too many rings in ringin'
      stop

      end
