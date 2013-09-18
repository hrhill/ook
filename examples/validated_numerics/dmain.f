      program dmain

      character*10 fread, fwrite
      parameter (fread='csrch.in',fwrite='csrch.out')
      integer nread, nwrite
      parameter (nread=1,nwrite=2)
      integer npmax
      parameter (npmax=20)
c     **********
c
c     This is a sample driver for the line search algorithm
c     included in MINPACK-2.
c
c     Problem parameters are read from the data file fread named
c     in the parameter statement above.
c     Information on the performance of the algorithm is printed
c     to the data file fwrite named in the parameter statement above.
c
c     This is only a sample driver, many other drivers are possible.
c
c     Subprograms called
c
c       MINPACK-2 ... dfcn, dcsrch
c
c     MINPACK-2 Project. November 1993.
c     Argonne National Laboratory.
c     Jorge J. More'.
c
c     **********
      character*60 task
      integer i, ic, nprob, ntries, ntry, nfev
      integer isave(3)
      integer ntrys(npmax), nprobs(npmax), nfevs(npmax), infos(npmax)
      double precision f, g, ftol, gtol, stp, stp0, stpmin, stpmax, xtol
      double precision factor, g0
      double precision dsave(13)
      double precision fs(npmax), gs(npmax), stps(npmax), stp0s(npmax)

      external dfcn, dcsrch

      open (nread,file=fread,status='old')
      open (nwrite,file=fwrite)

   10 continue

c     Read in problem parameters.

      read (nread,4000) nprob, ntries, stp0, ftol, gtol, xtol
      write (*,*) nprob

      if (nprob .eq. 0) go to 50

      factor = 1.0d0
      ic = 0
      do 30 ntry = 1, ntries
         ic = ic + 1

c        Initialize the search.

         stp = 0.0d0
         call dfcn(stp,f,g,nprob)
         stp = factor*stp0
         stpmin = 0.0d0
         stpmax = 4.0d0*max(1.0d0,stp)

         nfev = 0
         g0 = g

         task = 'START'
   20    continue

c        Call the search.

         call dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax,
     +               isave,dsave)

         if (task .eq. 'FG') then
            nfev = nfev + 1
            call dfcn(stp,f,g,nprob)
            go to 20
         end if

c        Record information on the algorithm.

         ntrys(ic) = ntry
         nprobs(ic) = nprob
         stp0s(ic) = factor*stp0
         nfevs(ic) = nfev
         infos(ic) = task(1:4) .eq. 'CONV'
         stps(ic) = stp
         fs(ic) = f
         gs(ic) = g
         factor = 1.0d2*factor
   30 continue

      write (nwrite,1000) ic, xtol, ftol, gtol, g0
      write (nwrite,2000)
      do 40 i = 1, ic
         write (nwrite,3000) nprobs(i), ntrys(i), nfevs(i), infos(i),
     +      stp0s(i), stps(i), fs(i), gs(i)
   40 continue

      go to 10

   50 continue


 1000 format (//' Summary of',i3,' calls to dcsrch'//'  xtol = ',d8.2,
     +       2x,' ftol = ',d8.2,2x,' gtol = ',d8.2,2x,' g0 = ',d8.2/)
 2000 format ('  nprob  ntry  nfev  info',5x,'x0',9x,'x',9x,'f',9x,'g'/)
 3000 format (4i6,2x,4d10.2)
 4000 format (2i5,4d10.2)
      end
