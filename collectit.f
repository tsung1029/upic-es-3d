c-----------------------------------------------------------------------
      subroutine PICOLLCT(if,ig,nxp,nvp,nblok)
c this subroutine performs a parallel gather of a vector, that is:
c ig(j,k,n) = if(j,k) on processor n
c if = input data
c ig = output array
c nxp = number of data values in vector
c nvp = number of processors
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      integer if, ig
      integer nxp, nvp, nblok
      dimension if(nxp,nblok), ig(nxp,nvp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer root, ierr
      root = 0
      call MPI_GATHER(if,nxp,mint,ig,nxp,mint,root,lgrp,ierr)
      return
      end
