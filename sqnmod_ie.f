!-----------------------------------------------------------------------
!
module sqn_ie
    ! subroutines for dealing with the Strange Quark Nugget

    use p0d, only: wtimer
    implicit none
    include 'mpif.h'

    private
    public :: SQNinit, ipushSQN

contains
    subroutine SQNinit(sqn, sqnr, dxSQN, edges, nppiSQN, partiSQN, npiSQNmax)
        ! this subroutine places particles for the Strange Quark Nugget
        ! WARNING: it does not know how to deal with crossing the periodic boundaries
        ! if you do insert one that crosses a boundary, expect to get some crap!
        implicit none
        integer :: npiSQNmax, npiSQNtotal, idproc, ierr, rc, funit
        real :: sqnr, dxSQN, temp, epsilon
        real, dimension(:,:) :: sqn
        real, dimension(:,:), pointer :: edges
        real, dimension(:,:,:), pointer :: partiSQN
        integer, dimension(:), pointer :: nppiSQN

        real :: x, y, z, xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0, rtemp

        ! this segment is used for mpi computers without shared memory!

        ! determine the rank of the calling process in the communicator
        call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierr)

        ! only test for placement if we're inside a cube with side 2*r centered at the SQN's position
!        xmin = sqn(1, 1) - sqnr
!        xmax = sqn(1, 1) + sqnr
!        ymin = sqn(2, 1) - sqnr
!        ymax = sqn(2, 1) + sqnr
!        zmin = sqn(3, 1) - sqnr
!        zmax = sqn(3, 1) + sqnr

        ! new layout, the grid of particles is offset from the center by dxSQN/2
        ! so, particle centers must be within sqnr - dxSQN/2 of the center
        xmin = sqn(1, 1) - sqnr + dxSQN/2
        xmax = sqn(1, 1) + sqnr - dxSQN/2
        ymin = sqn(2, 1) - sqnr + dxSQN/2
        ymax = sqn(2, 1) + sqnr - dxSQN/2
        zmin = sqn(3, 1) - sqnr + dxSQN/2
        zmax = sqn(3, 1) + sqnr - dxSQN/2

        ! check if we'll have some particles on this partition
        if (edges(2, 1) < ymin .or. edges(1, 1) > ymax) return
        if (edges(4, 1) < zmin .or. edges(3, 1) > zmax) return

        ! calculate the mins and maxes for this partition
        ymin = max(ymin, edges(1, 1))
        ymax = min(ymax, edges(2, 1))
        zmin = max(zmin, edges(3, 1))
        zmax = min(zmax, edges(4, 1))

        ! if(idproc==3) print*, idproc,'init ymin=',ymin

        ! set the mins to land an integer number of grid spacings from the SQN's center
        ! so that we can place a particle at the min and maintain even spacing across partitions
!        if ((sqn(1, 1) - xmin) /= aint((sqn(1, 1) - xmin)/dxSQN) * dxSQN) then
!            xmin = sqn(1, 1) - aint((sqn(1, 1) - xmin)/dxSQN) * dxSQN
!        endif
!
!        if (abs(sqn(2, 1) - ymin) /= aint(abs(sqn(2, 1) - ymin)/dxSQN) * dxSQN) then
!            if (ymin < sqn(2, 1)) then
!                ymin = sqn(2, 1) - aint((sqn(2, 1) - ymin)/dxSQN) * dxSQN
!            else
!                ymin = sqn(2, 1) + aint((ymin - sqn(2, 1))/dxSQN) * dxSQN
!            endif
!        endif
!
!        if (abs(sqn(3, 1) - zmin) /= aint(abs(sqn(3, 1) - zmin)/dxSQN) * dxSQN) then
!            if (zmin < sqn(3, 1)) then
!                zmin = sqn(3, 1) - aint((sqn(3, 1) - zmin)/dxSQN) * dxSQN
!            else
!                zmin = sqn(3, 1) + aint((zmin - sqn(3, 1))/dxSQN) * dxSQN
!            endif
!        endif

        ! for new layout
        ! set the mins to land an integer number of grid spacings (minus 1/2 grid spacing)
        ! from the SQN's center so that we can place a particle at the min and maintain
        ! even spacing across partitions
        epsilon = dxSQN/100
        temp = anint(abs(sqn(2, 1) - ymin + dxSQN/2)/dxSQN) * dxSQN
        if (abs(sqn(2, 1) - ymin + dxSQN/2) < temp-epsilon .or. abs(sqn(2, 1) - ymin + dxSQN/2) > temp+epsilon) then
            if (ymin < sqn(2, 1)) then
                ymin = sqn(2, 1) - anint((sqn(2, 1) - ymin)/dxSQN) * dxSQN + dxSQN/2
            else
                ymin = sqn(2, 1) + anint((ymin - sqn(2, 1))/dxSQN) * dxSQN + dxSQN/2
            endif
        endif

        temp = anint(abs(sqn(3, 1) - zmin + dxSQN/2)/dxSQN) * dxSQN
        if (abs(sqn(3, 1) - zmin + dxSQN/2) < temp-epsilon .or.  abs(sqn(3, 1) - zmin + dxSQN/2) > temp+epsilon) then
            if (zmin < sqn(3, 1)) then
                zmin = sqn(3, 1) - anint((sqn(3, 1) - zmin)/dxSQN) * dxSQN + dxSQN/2
            else
                zmin = sqn(3, 1) + anint((zmin - sqn(3, 1))/dxSQN) * dxSQN + dxSQN/2
            endif
        endif

        ! do the same for the max
        temp = anint(abs(ymax - sqn(2, 1) + dxSQN/2)/dxSQN) * dxSQN
        if (abs(ymax - sqn(2, 1) + dxSQN/2) < temp-epsilon .or. abs(ymax - sqn(2, 1) + dxSQN/2) > temp+epsilon) then
            if (ymax > sqn(2, 1)) then
                ymax = sqn(2, 1) + anint((ymax - sqn(2, 1))/dxSQN) * dxSQN - dxSQN/2
            else
                ymax = sqn(2, 1) - anint((sqn(2, 1) - ymax)/dxSQN) * dxSQN - dxSQN/2
            endif
        endif

        temp = anint(abs(zmax - sqn(3, 1) + dxSQN/2)/dxSQN) * dxSQN
        if (abs(zmax - sqn(3, 1) + dxSQN/2) < temp-epsilon .or. abs(zmax - sqn(3, 1) + dxSQN/2) > temp+epsilon) then
            if (zmax > sqn(3, 1)) then
                zmax = sqn(3, 1) + anint((zmax - sqn(3, 1))/dxSQN) * dxSQN - dxSQN/2
            else
                zmax = sqn(3, 1) - anint((sqn(3, 1) - zmax)/dxSQN) * dxSQN - dxSQN/2
            endif
        endif


        ! make sure that partitions don't replicate each others' placements on the edges
        if (ymax==edges(2,1)) ymax = ymax - dxSQN/2
        if (zmax==edges(4,1)) zmax = zmax - dxSQN/2

!        print*, idproc, 'ymin, ymax, edges = ', ymin, ymax, edges(1,1), edges(2,1)
!        print*, idproc, 'zmin, zmax, edges = ', zmin, zmax, edges(3,1), edges(4,1)

        ! debugging checking
        if (ymin < edges(1, 1)) then
            print*, 'we got an edge problem in y with the SQN on task', idproc
            call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        endif

        if (zmin < edges(3, 1)) then
            print*, 'we got an edge problem in z with the SQN on task', idproc
            call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        endif


        ! finally, place the particles
        do x = xmin, xmax+dxSQN/10, dxSQN
            do y = ymin, ymax+dxSQN/10, dxSQN
                do z = zmin, zmax+dxSQN/10, dxSQN
                    rtemp = sqrt((x - sqn(1, 1))**2 + (y - sqn(2, 1))**2 + (z - sqn(3, 1))**2)
                    if (rtemp <= sqnr) then
                    !if (rtemp < sqnr) then
                        nppiSQN = nppiSQN + 1

                        if (nppiSQN(1) > npiSQNmax) then
                            print*, 'Error: exceeded number of particles for Strange Quark Nugget on partition', idproc
                            print*, 'nppiSQN = ', nppiSQN(1), 'npiSQNmax = ', npiSQNmax
                            print*, idproc, 'zmin, zmax = ', zmin, zmax
                            print*, idproc, 'ymin, ymax = ', ymin, ymax
                            print*, idproc, 'xmin, xmax = ', xmin, xmax
                            call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
                        endif

                        partiSQN(1, nppiSQN(1), 1) = x
                        partiSQN(2, nppiSQN(1), 1) = y
                        partiSQN(3, nppiSQN(1), 1) = z
                        partiSQN(4, nppiSQN(1), 1) = sqn(4, 1)
                        partiSQN(5, nppiSQN(1), 1) = sqn(5, 1)
                        partiSQN(6, nppiSQN(1), 1) = sqn(6, 1)

                        ! write(idproc,*) x, y, z
                    endif
                enddo
            enddo
        enddo
        ! flush(idproc)
    end subroutine SQNinit

    subroutine ipushSQN(sqn, part, fxyz, npp, numSQNpart, noff, qbm, dt, ek, tpush, nx, ny, nz, ipbc)
        ! push particles with 3d electrostatic fields, 2d partition
        
        implicit none
        integer :: nx, ny, nz, ipbc, ierr, rc, idproc
        real :: qbm, dt, ek, tpush
        real, dimension(:,:) :: sqn
        real, dimension(:,:,:), pointer :: part
        real, dimension(:,:,:,:,:), pointer :: fxyz
        integer, dimension(:), pointer :: npp
        integer, dimension(:,:), pointer :: noff
        ! local data
        integer :: mnblok, nxv, nypmx, nzpmx, nxyzp, numSQNpart
        real :: tp
        ! sum1 is for kinetic energy
        double precision :: dtime

        ! common block for parallel processing
        integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
        ! lstat = length of status array
        parameter(lstat = 10)
        ! nproc = number of real or virtual processors obtained
        ! lgrp = current communicator
        ! lworld = MPI_COMM_WORLD communicator
        common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld

        mnblok = size(part, 3)
        nxv = size(fxyz, 2); nypmx = size(fxyz, 3); nzpmx = size(fxyz, 4)
        nxyzp = nxv * nypmx * nzpmx

        ! initialize timer
        call wtimer(tp, dtime, -1)

        ! the quadratic look-ahead pusher

        call pushSQN(sqn, part, fxyz,npp,numSQNpart, noff,qbm,dt,ek,nx,ny,nz,mnblok,nxv,nypmx,nxyzp,ipbc)

        ! record time
        call wtimer(tp, dtime)
        tpush = tpush + tp
    end subroutine ipushSQN

    subroutine pushSQN(sqn, part, fxyz,npp,numSQNpart,noff,qbm,dt,ek,nx,ny,nz,mnblok,nxv,nypmx,nxyzp,ipbc)
        ! the quadratic look-ahead pusher
        ! WARNING: Only works for periodic boundary conditions!!!
        ! for 3d code, this subroutine updates particle co-ordinates and
        ! velocities using leap-frog scheme in time and second-order spline
        ! interpolation in space, for distributed data
        ! with 2D spatial decomposition
        ! scalar version using guard cells, integer conversion precalculation,
        ! and 1d addressing
        ! cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
        ! 238 flops/particle, 87 loads, 6 stores
        ! input: all, output: part, ek
        ! equations used are:
        ! vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
        ! vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
        ! vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
        ! where q/m is charge/mass, and
        ! x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
        ! z(t+dt) = z(t) + vz(t+dt/2)*dt
        ! fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
        ! are approximated by interpolation from the nearest grid points:
        ! fx(x,y,z) = (.75-dz**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l)+
        ! (.5*(.5+dx)**2)*fx(n+1,m,l)+(.5*(.5-dx)**2)*fx(n-1,m,l)) +
        ! (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l)+
        ! (.5*(.5+dx)**2)*fx(n+1,m+1,l)+(.5*(.5-dx)**2)*fx(n-1,m+1,l)) +
        ! (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l)+
        ! (.5*(.5+dx)**2)*fx(n+1,m-1,l)+(.5*(.5-dx)**2)*fx(n-1,m-1,l))) +
        ! (.5*(.5+dz)**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l+1)+
        ! (.5*(.5+dx)**2)*fx(n+1,m,l+1)+(.5*(.5-dx)**2)*fx(n-1,m,l+1)) +
        ! (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l+1)+
        ! (.5*(.5+dx)**2)*fx(n+1,m+1,l+1)+(.5*(.5-dx)**2)*fx(n-1,m+1,l+1)) +
        ! (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l+1)+
        ! (.5*(.5+dx)**2)*fx(n+1,m-1,l+1)+(.5*(.5-dx)**2)*fx(n-1,m-1,l+1)))
        ! (.5*(.5-dz)**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l-1)+
        ! (.5*(.5+dx)**2)*fx(n+1,m,l-1)+(.5*(.5-dx)**2)*fx(n-1,m,l-1)) +
        ! (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l-1)+
        ! (.5*(.5+dx)**2)*fx(n+1,m+1,l-1)+(.5*(.5-dx)**2)*fx(n-1,m+1,l-1)) +
        ! (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l-1)+
        ! (.5*(.5+dx)**2)*fx(n+1,m-1,l-1)+(.5*(.5-dx)**2)*fx(n-1,m-1,l-1)))
        ! where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
        ! and similarly for fy(x,y,z) and fz(x,y,z)
        ! part(1,n,m) = position x of particle n in partition m
        ! part(2,n,m) = position y of particle n in partition m
        ! part(3,n,m) = position z of particle n in partition m
        ! part(4,n,m) = velocity vx of particle n in partition m
        ! part(5,n,m) = velocity vy of particle n in partition m
        ! part(6,n,m) = velocity vz of particle n in partition m
        ! fxyz(1,j+1,k,l,m) = x component of force/charge at grid (j,kk,ll)
        ! fxyz(2,j+1,k,l,m) = y component of force/charge at grid (j,kk,ll)
        ! fxyz(3,j+1,k,l,m) = z component of force/charge at grid (j,kk,ll)
        ! in other words, fxyz are the convolutions of the electric field
        ! over the particle shape,
        ! where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
        ! npp(m) = number of particles in partition m
        ! noff(1,m) = lowermost global gridpoint in y in particle partition m
        ! noff(2,m) = backmost global gridpoint in z in particle partition m
        ! qbm = particle charge/mass ratio
        ! dt = time interval between successive calculations
        ! kinetic energy/mass at time t is also calculated, using
        ! ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
        ! (vz(t+dt/2)+vz(t-dt/2))**2)
        ! nx/ny/nz = system length in x/y/z direction
        ! idimp = size of phase space = 6
        ! npmax = maximum number of particles in each partition
        ! mnblok = number of particle partitions.
        ! nxv = first virtual dimension of field array, must be >= nx+3
        ! nypmx = maximum size of particle partition in y, including guard cells
        ! nxyzp = second dimension of field array, must be >= nxv*nypmx*nzpmx
        ! idds = dimensionality of domain decomposition
        ! ipbc = particle boundary condition = (0,1,2,3) =
        ! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
        implicit none
        integer :: nx, ny, nz, ipbc, ierr, rc, idproc
        real :: qbm, dt, ek
        real, dimension(:,:) :: sqn
        real, dimension(:,:,:), pointer :: part
        real, dimension(3,nxyzp,mnblok) :: fxyz
        integer, dimension(:), pointer :: npp
        integer, dimension(:,:), pointer :: noff
        ! local data
        integer :: idimp, npmax, numSQNpart, mnblok, nxv, nypmx, nzpmx, nxyzp, idds
        ! for adding up the force per unit charge
        real :: ftx, fty, ftz

        ! things that were implicitly declared in the f77 routine
        real :: edgelx, edgely, edgelz, edgerx, edgery, edgerz
        real :: amx, amy, amz, dx, dy, dz, dxn, dyn, dzn
        real :: dxl, dyl, dzl, dxp, dyp, dzp, dx1, dy1, dx2, dy2
        real :: qtm
        integer :: j, m, mm, nn, lll, mmm, nnn, nop1, ml, mn, mp, mnoff, lnoff
        integer :: nxvy

        ! common block for parallel processing
        integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
        ! lstat = length of status array
        parameter(lstat = 10)
        ! nproc = number of real or virtual processors obtained
        ! lgrp = current communicator
        ! lworld = MPI_COMM_WORLD communicator
        common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld

        call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierr)

        idimp = size(part, 1); npmax = size(part, 2)
        mnblok = size(part, 3)
        idds = size(noff, 1)

        qtm = qbm*dt
        nxvy = nxv * nypmx
        ! set boundary values
        if (ipbc .eq. 1) then
            edgelx = 0.
            edgerx = float(nx)
        else if (ipbc .eq. 2) then
            edgelx = 1.
            edgely = 1.
            edgelz = 1.
            edgerx = float(nx - 1)
            edgery = float(ny - 1)
            edgerz = float(nz - 1)
        else if (ipbc .eq. 3) then
            edgelx = 1.
            edgely = 1.
            edgerx = float(nx - 1)
            edgery = float(ny - 1)
        endif


        ftx = 0.0; fty = 0.0; ftz = 0.0;
        do m = 1, mnblok
            if (npp(m) .lt. 1) cycle
            mnoff = noff(1, m)
            lnoff = noff(2, m)
            ! begin first particle
            nnn = part(1, 1, m) + .5
            mmm = part(2, 1, m) + .5
            lll = part(3, 1, m) + .5
            dxn = part(1, 1, m) - float(nnn)
            dyn = part(2, 1, m) - float(mmm)
            dzn = part(3, 1, m) - float(lll)
            mmm = mmm - mnoff
            lll = lll - lnoff
            nop1 = npp(m) - 1

            ! find interpolation weights
            do j = 1, nop1
                nn = nnn + 1
                mm = nxv * mmm + nxvy * lll
                nnn = part(1, j + 1, m) + .5
                mmm = part(2, j + 1, m) + .5
                lll = part(3, j + 1, m) + .5
                dx = dxn
                dy = dyn
                dz = dzn
                dxn = part(1, j + 1, m) - float(nnn)
                dyn = part(2, j + 1, m) - float(mmm)
                dzn = part(3, j + 1, m) - float(lll)
                amx = .75 - dx * dx
                dxl = .5 * (.5 - dx)**2
                dxp = .5 * (.5 + dx)**2
                ml = mm + nn
                dyl = .5 * (.5 - dy)**2
                amy = .75 - dy * dy
                dyp = .5 * (.5 + dy)**2
                mn = ml + nxv
                dx1 = dxl * dyl
                dx2 = amx * dyl
                dyl = dxp * dyl
                dzl = .5 * (.5 - dz)**2
                dy1 = dxl * amy
                dy2 = amx * amy
                amy = dxp * amy
                amz = .75 - dz * dz
                dxl = dxl * dyp
                amx = amx * dyp
                dyp = dxp * dyp
                dzp = .5 * (.5 + dz)**2
                mp = mn + nxv
                mmm = mmm - mnoff
                lll = lll - lnoff

                ! find acceleration
                dx = dzl * (dx1 * fxyz(1, ml, m) + dx2 * fxyz(1, ml + 1, m) + dyl * fxyz(1, ml + 2,&
                &m) + dy1 * fxyz(1, mn, m) + dy2 * fxyz(1, mn + 1, m) + amy * fxyz(1, mn + 2, m) +&
                &dxl * fxyz(1, mp, m) + amx * fxyz(1, mp + 1, m) + dyp * fxyz(1, mp + 2, m))

                dy = dzl * (dx1 * fxyz(2, ml, m) + dx2 * fxyz(2, ml + 1, m) + dyl * fxyz(2, ml + 2,&
                &m) + dy1 * fxyz(2, mn, m) + dy2 * fxyz(2, mn + 1, m) + amy * fxyz(2, mn + 2, m) +&
                &dxl * fxyz(2, mp, m) + amx * fxyz(2, mp + 1, m) + dyp * fxyz(2, mp + 2, m))

                dz = dzl * (dx1 * fxyz(3, ml, m) + dx2 * fxyz(3, ml + 1, m) + dyl * fxyz(3, ml + 2,&
                &m) + dy1 * fxyz(3, mn, m) + dy2 * fxyz(3, mn + 1, m) + amy * fxyz(3, mn + 2, m) +&
                &dxl * fxyz(3, mp, m) + amx * fxyz(3, mp + 1, m) + dyp * fxyz(3, mp + 2, m))

                ml = ml + nxvy
                mn = mn + nxvy
                mp = mp + nxvy

                dx = dx + amz * (dx1 * fxyz(1, ml, m) + dx2 * fxyz(1, ml + 1, m) + dyl * fxyz(1,&
                &ml + 2, m) + dy1 * fxyz(1, mn, m) + dy2 * fxyz(1, mn + 1, m) + amy * fxyz(1, mn + 2,&
                &m) + dxl * fxyz(1, mp, m) + amx * fxyz(1, mp + 1, m) + dyp * fxyz(1, mp + 2, m))

                dy = dy + amz * (dx1 * fxyz(2, ml, m) + dx2 * fxyz(2, ml + 1, m) + dyl * fxyz(2,&
                &ml + 2, m) + dy1 * fxyz(2, mn, m) + dy2 * fxyz(2, mn + 1, m) + amy * fxyz(2, mn + 2,&
                &m) + dxl * fxyz(2, mp, m) + amx * fxyz(2, mp + 1, m) + dyp * fxyz(2, mp + 2, m))

                dz = dz + amz * (dx1 * fxyz(3, ml, m) + dx2 * fxyz(3, ml + 1, m) + dyl * fxyz(3,&
                &ml + 2, m) + dy1 * fxyz(3, mn, m) + dy2 * fxyz(3, mn + 1, m) + amy * fxyz(3, mn + 2,&
                &m) + dxl * fxyz(3, mp, m) + amx * fxyz(3, mp + 1, m) + dyp * fxyz(3, mp + 2, m))

                ml = ml + nxvy
                mn = mn + nxvy
                mp = mp + nxvy

                dx = dx + dzp * (dx1 * fxyz(1, ml, m) + dx2 * fxyz(1, ml + 1, m) + dyl * fxyz(1,&
                &ml + 2, m) + dy1 * fxyz(1, mn, m) + dy2 * fxyz(1, mn + 1, m) + amy * fxyz(1, mn + 2,&
                &m) + dxl * fxyz(1, mp, m) + amx * fxyz(1, mp + 1, m) + dyp * fxyz(1, mp + 2, m))

                dy = dy + dzp * (dx1 * fxyz(2, ml, m) + dx2 * fxyz(2, ml + 1, m) + dyl * fxyz(2,&
                &ml + 2, m) + dy1 * fxyz(2, mn, m) + dy2 * fxyz(2, mn + 1, m) + amy * fxyz(2, mn + 2,&
                &m) + dxl * fxyz(2, mp, m) + amx * fxyz(2, mp + 1, m) + dyp * fxyz(2, mp + 2, m))

                dz = dz + dzp * (dx1 * fxyz(3, ml, m) + dx2 * fxyz(3, ml + 1, m) + dyl * fxyz(3,&
                &ml + 2, m) + dy1 * fxyz(3, mn, m) + dy2 * fxyz(3, mn + 1, m) + amy * fxyz(3, mn + 2,&
                &m) + dxl * fxyz(3, mp, m) + amx * fxyz(3, mp + 1, m) + dyp * fxyz(3, mp + 2, m))
                
                ftx = ftx + dx
                fty = fty + dy
                ftz = ftz + dz
            enddo

            ! calculate force for the last particle
            nn = nnn + 1
            mm = nxv * mmm + nxvy * lll
            amx = .75 - dxn * dxn
            dxl = .5 * (.5 - dxn)**2
            dxp = .5 * (.5 + dxn)**2
            ml = mm + nn
            dyl = .5 * (.5 - dyn)**2
            amy = .75 - dyn * dyn
            dyp = .5 * (.5 + dyn)**2
            mn = ml + nxv
            dx1 = dxl * dyl
            dx2 = amx * dyl
            dyl = dxp * dyl
            dzl = .5 * (.5 - dzn)**2
            dy1 = dxl * amy
            dy2 = amx * amy
            amy = dxp * amy
            amz = .75 - dzn * dzn
            dxl = dxl * dyp
            amx = amx * dyp
            dyp = dxp * dyp
            dzp = .5 * (.5 + dzn)**2
            mp = mn + nxv

            ! find acceleration
            dx = dzl * (dx1 * fxyz(1, ml, m) + dx2 * fxyz(1, ml + 1, m) + dyl * fxyz(1, ml + 2,&
            &m) + dy1 * fxyz(1, mn, m) + dy2 * fxyz(1, mn + 1, m) + amy * fxyz(1, mn + 2, m) +&
            &dxl * fxyz(1, mp, m) + amx * fxyz(1, mp + 1, m) + dyp * fxyz(1, mp + 2, m))

            dy = dzl * (dx1 * fxyz(2, ml, m) + dx2 * fxyz(2, ml + 1, m) + dyl * fxyz(2, ml + 2,&
            &m) + dy1 * fxyz(2, mn, m) + dy2 * fxyz(2, mn + 1, m) + amy * fxyz(2, mn + 2, m) +&
            &dxl * fxyz(2, mp, m) + amx * fxyz(2, mp + 1, m) + dyp * fxyz(2, mp + 2, m))

            dz = dzl * (dx1 * fxyz(3, ml, m) + dx2 * fxyz(3, ml + 1, m) + dyl * fxyz(3, ml + 2,&
            &m) + dy1 * fxyz(3, mn, m) + dy2 * fxyz(3, mn + 1, m) + amy * fxyz(3, mn + 2, m) +&
            &dxl * fxyz(3, mp, m) + amx * fxyz(3, mp + 1, m) + dyp * fxyz(3, mp + 2, m))

            ml = ml + nxvy
            mn = mn + nxvy
            mp = mp + nxvy

            dx = dx + amz * (dx1 * fxyz(1, ml, m) + dx2 * fxyz(1, ml + 1, m) + dyl * fxyz(1,&
            &ml + 2, m) + dy1 * fxyz(1, mn, m) + dy2 * fxyz(1, mn + 1, m) + amy * fxyz(1, mn + 2,&
            &m) + dxl * fxyz(1, mp, m) + amx * fxyz(1, mp + 1, m) + dyp * fxyz(1, mp + 2, m))

            dy = dy + amz * (dx1 * fxyz(2, ml, m) + dx2 * fxyz(2, ml + 1, m) + dyl * fxyz(2,&
            &ml + 2, m) + dy1 * fxyz(2, mn, m) + dy2 * fxyz(2, mn + 1, m) + amy * fxyz(2, mn + 2,&
            &m) + dxl * fxyz(2, mp, m) + amx * fxyz(2, mp + 1, m) + dyp * fxyz(2, mp + 2, m))

            dz = dz + amz * (dx1 * fxyz(3, ml, m) + dx2 * fxyz(3, ml + 1, m) + dyl * fxyz(3,&
            &ml + 2, m) + dy1 * fxyz(3, mn, m) + dy2 * fxyz(3, mn + 1, m) + amy * fxyz(3, mn + 2,&
            &m) + dxl * fxyz(3, mp, m) + amx * fxyz(3, mp + 1, m) + dyp * fxyz(3, mp + 2, m))

            ml = ml + nxvy
            mn = mn + nxvy
            mp = mp + nxvy

            dx = dx + dzp * (dx1 * fxyz(1, ml, m) + dx2 * fxyz(1, ml + 1, m) + dyl * fxyz(1,&
            &ml + 2, m) + dy1 * fxyz(1, mn, m) + dy2 * fxyz(1, mn + 1, m) + amy * fxyz(1, mn + 2,&
            &m) + dxl * fxyz(1, mp, m) + amx * fxyz(1, mp + 1, m) + dyp * fxyz(1, mp + 2, m))

            dy =  dy + dzp * (dx1 * fxyz(2, ml, m) + dx2 * fxyz(2, ml + 1, m) + dyl * fxyz(2,&
            &ml + 2, m) + dy1 * fxyz(2, mn, m) + dy2 * fxyz(2, mn + 1, m) + amy * fxyz(2, mn + 2,&
            &m) + dxl * fxyz(2, mp, m) + amx * fxyz(2, mp + 1, m) + dyp * fxyz(2, mp + 2, m))

            dz = dz + dzp * (dx1 * fxyz(3, ml, m) + dx2 * fxyz(3, ml + 1, m) + dyl * fxyz(3,&
            &ml + 2, m) + dy1 * fxyz(3, mn, m) + dy2 * fxyz(3, mn + 1, m) + amy * fxyz(3, mn + 2,&
            &m) + dxl * fxyz(3, mp, m) + amx * fxyz(3, mp + 1, m) + dyp * fxyz(3, mp + 2, m))

            ftx = ftx + dx
            fty = fty + dy
            ftz = ftz + dz

        enddo

        !sum the forces across partitions
        call MPI_ALLREDUCE(ftx,dx,1,mreal,MPI_SUM,lworld,ierr)
        if (ierr /= MPI_SUCCESS) then
            rc = ierr
            print*, 'Could not add forces across tasks'
            call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        endif
        call MPI_ALLREDUCE(fty,dy,1,mreal,MPI_SUM,lworld,ierr)
        if (ierr /= MPI_SUCCESS) then
            rc = ierr
            print*, 'Could not add forces across tasks'
            call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        endif
        call MPI_ALLREDUCE(ftz,dz,1,mreal,MPI_SUM,lworld,ierr)
        if (ierr /= MPI_SUCCESS) then
            rc = ierr
            print*, 'Could not add forces across tasks'
            call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        endif

!        if (idproc==0) then
!            print*, 'fx = ', dx, 'fy = ', dy, 'fz = ', dz
!        endif

        ! new velocity
        dx = sqn(4,1) + qtm * dx / numSQNpart
        dy = sqn(5,1) + qtm * dy / numSQNpart
        dz = sqn(6,1) + qtm * dz / numSQNpart

        ! average kinetic energy
        ek = 0.125*((dx + sqn(4,1))**2 + (dy + sqn(5,1))**2 + (dz + sqn(6,1))**2)
        sqn(4,1) = dx
        sqn(5,1) = dy
        sqn(6,1) = dz

        ! new position of the entire SQN
        dx = sqn(1,1) + dx * dt
        dy = sqn(2,1) + dy * dt
        dz = sqn(3,1) + dz * dt

        ! periodic boundary conditions
!        if (ipbc .eq. 1) then
            if (dx .lt. edgelx) dx = dx + edgerx
            if (dx .ge. edgerx) dx = dx - edgerx
            if (dy .lt. edgely) dy = dy + edgery
            if (dy .ge. edgery) dy = dy - edgery
            if (dz .lt. edgelz) dz = dz + edgerz
            if (dz .ge. edgerz) dz = dz - edgerz
        ! reflecting boundary conditions
!        else if (ipbc .eq. 2) then
!            if ((dx .lt. edgelx).or.(dx .ge. edgerx)) then
!            dx = part(1, j, m)
!            part(4, j, m) = -part(4, j, m)
!            endif
!            if ((dy .lt. edgely).or.(dy .ge. edgery)) then
!            dy = part(2, j, m)
!            part(5, j, m) = -part(5, j, m)
!            endif
!            if ((dz .lt. edgelz).or.(dz .ge. edgerz)) then
!            dz = part(3, j, m)
!            part(6, j, m) = -part(6, j, m)
!            endif
!        ! mixed reflecting/periodic boundary conditions
!        else if (ipbc .eq. 3) then
!            if ((dx .lt. edgelx).or.(dx .ge. edgerx)) then
!            dx = part(1, j, m)
!            part(4, j, m) = -part(4, j, m)
!            endif
!            if ((dy .lt. edgely).or.(dy .ge. edgery)) then
!            dy = part(2, j, m)
!            part(5, j, m) = -part(5, j, m)
!            endif
!        endif

        ! set new position
        sqn(1, 1) = dx
        sqn(2, 1) = dy
        sqn(3, 1) = dz

        ! move particles in the SQN and set their velocities
        do m = 1, mnblok
            if (npp(m) .lt. 1) continue
            do j = 1, npp(m)
                ! set velocity
                part(4,j,m) = sqn(4,1)
                part(5,j,m) = sqn(5,1)
                part(6,j,m) = sqn(6,1)

                ! new position
                dx = part(1, j, m) + sqn(4,1) * dt
                dy = part(2, j, m) + sqn(5,1) * dt
                dz = part(3, j, m) + sqn(6,1) * dt

                ! periodic boundary conditions
                ! if (ipbc .eq. 1) then
                    if (dx .lt. edgelx) dx = dx + edgerx
                    if (dx .ge. edgerx) dx = dx - edgerx
                    ! reflecting boundary conditions
!                else if (ipbc .eq. 2) then
!                    if ((dx .lt. edgelx).or.(dx .ge. edgerx)) then
!                    dx = part(1, j, m)
!                    part(4, j, m) = -part(4, j, m)
!                    endif
!                    if ((dy .lt. edgely).or.(dy .ge. edgery)) then
!                    dy = part(2, j, m)
!                    part(5, j, m) = -part(5, j, m)
!                    endif
!                    if ((dz .lt. edgelz).or.(dz .ge. edgerz)) then
!                    dz = part(3, j, m)
!                    part(6, j, m) = -part(6, j, m)
!                    endif
!                    ! mixed reflecting/periodic boundary conditions
!                else if (ipbc .eq. 3) then
!                    if ((dx .lt. edgelx).or.(dx .ge. edgerx)) then
!                    dx = part(1, j, m)
!                    part(4, j, m) = -part(4, j, m)
!                    endif
!                    if ((dy .lt. edgely).or.(dy .ge. edgery)) then
!                    dy = part(2, j, m)
!                    part(5, j, m) = -part(5, j, m)
!                    endif
!                endif

                ! set new position
                part(1, j, m) = dx
                part(2, j, m) = dy
                part(3, j, m) = dz
            enddo
        enddo

        end subroutine pushSQN
end module sqn_ie
