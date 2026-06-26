      subroutine read_xyz(fname, natoms, ian, c, IOut, IERR)
c
c     Robust XYZ reader for GICForge (Fortran77)
c
c     INPUT:
c       fname   - XYZ file name
c       IOut    - output unit for messages
c
c     OUTPUT:
c       natoms  - number of atoms read
c       ian(*)  - atomic numbers
c       c(3,*)  - Cartesian coordinates (Angstrom)
c       IERR    - error code (0 = OK, <0 = error)
c
c     ERROR CODES (IERR):
c       0   : successful read
c      -1   : cannot open XYZ file
c      -2   : error reading XYZ header (natoms or comment line)
c      -3   : unexpected end-of-file while reading atom lines
c      -4   : cannot parse an XYZ coordinate line
c      -5   : invalid atomic symbol (El2IAN returned <= 0)
c
c     NOTES:
c     - Free-format read is used (no limit on decimal digits)
c     - Atomic symbols or atomic numbers are accepted
c     - Atomic symbols are converted using El2IAN(NoCase, AtSymb)
c     - Atom indexing is 1-based (Fortran / XYZ convention)
c     - The routine never calls STOP; the caller must handle IERR
c

      character*(*) fname
      integer natoms
      integer ian(*)
      double precision c(3,*)
      integer IOut
      integer IERR

      integer iu, i
      character*256 line
      character*8 atsym
      logical nocase
      double precision x, y, z
      integer iz

      external El2IAN

c     Initialize
      IERR = 0
      natoms = 0
      nocase = .true.

c     --- Open XYZ file ---
      iu = 97
      open(unit=iu, file=fname, status='old', err=900)

c     --- Read number of atoms ---
      read(iu,*,err=910,end=910) natoms

c     --- Skip comment line ---
      read(iu,'(A)',err=910,end=910) line

c     --- Read atoms ---
      do 100 i = 1, natoms

         read(iu,'(A)',err=930,end=930) line

c        Try reading atomic number first
         read(line,*,err=110) iz, x, y, z
         ian(i) = iz
         c(1,i) = x
         c(2,i) = y
         c(3,i) = z
         go to 100

c        Otherwise read atomic symbol
 110     continue
         read(line,*,err=940) atsym, x, y, z
         iz = El2IAN(nocase, atsym)
         if (iz .le. 0) go to 950
         ian(i) = iz
         c(1,i) = x
         c(2,i) = y
         c(3,i) = z

 100  continue

      close(iu)
      return

c     --- Error handling ---
 900  continue
      write(IOut,*) 'ERROR(read_xyz): cannot open XYZ file ', fname
      IERR = -1
      return

 910  continue
      write(IOut,*) 'ERROR(read_xyz): cannot read XYZ header in ', fname
      IERR = -2
      return

 930  continue
      write(IOut,*) 'ERROR(read_xyz): unexpected end of XYZ file'
      IERR = -3
      return

 940  continue
      write(IOut,*) 'ERROR(read_xyz): cannot parse XYZ line:'
      write(IOut,'(A)') line
      IERR = -4
      return

 950  continue
      write(IOut,*) 'ERROR(read_xyz): invalid atomic symbol: ', atsym
      IERR = -5
      return

      end
