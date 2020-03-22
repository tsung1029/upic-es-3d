
!     
! File:   pinit32mod_ie.f
! Author: ellis38
!
! Created on June 7, 2010, 3:06 PM
!

module pinit32d_ie

    implicit none
    private
    public :: smoothIon,qetest, qitest
    public :: pinput3_ie, sendnml_ie
    public :: nteste, ntesti
    public :: qmetest, qmitest, rmteste, rmtesti
    public :: fixedvel, track_teste, track_testi
    public :: sqn, rsqn, qmsqn, rmsqn
    public :: check_period
    public :: scale_test

    integer, parameter :: max_test = 50   

    real, dimension(6, max_test) :: qetest = 0
    real, dimension(6, max_test) :: qitest = 0
    integer :: fixedvel = 1
    integer :: track_teste = 0
    integer :: track_testi = 0
    integer :: nteste = 0
    integer :: ntesti = 0
    integer :: scale_test = 0    ! scale the test charge inputs by vth (really helpful when shrinking cell widths)
    real :: qmetest = -1.0
    real :: rmteste = 0.0
    real :: rmtesti = 0.0
    real :: qmitest = 1.0
    integer :: smoothIon = 0 ! for smooth neutralizing ion background
    integer :: check_period = 0 ! period for checking whether we should continue for another period
                                ! or dump a restart and quit

    real, dimension(6,1) :: sqn ! strange quark nugget inital position and velocity
    real :: rsqn = 0.0  ! strange quark nugget radius
    real :: qmsqn = 0.0 ! strange quark nugget charge
    real :: rmsqn = 0.0 ! ratio of SQN mass to electron mass
    save

    namelist /pinput3_ie/ smoothIon, nteste, qmetest, qetest, ntesti, &
        & qmitest, qitest, fixedvel, rmteste, rmtesti, sqn, rsqn, qmsqn, rmsqn, &
        & check_period, track_teste, track_testi, scale_test

contains

    subroutine sendnml_ie()
        implicit none
        integer, parameter :: lenml = 21
        double precision, dimension(lenml) :: ddata
        ddata(1) = qmetest
        ddata(2) = qmitest
        ddata(3) = smoothIon
        ddata(4) = sqn(1,1)
        ddata(5) = sqn(2,1)
        ddata(6) = sqn(3,1)
        ddata(7) = sqn(4,1)
        ddata(8) = sqn(5,1)
        ddata(9) = sqn(6,1)
        ddata(10) = rsqn
        ddata(11) = qmsqn
        ddata(12) = rmsqn
        ddata(13) = check_period
        ddata(14) = fixedvel
        ddata(15) = rmteste
        ddata(16) = rmtesti
        ddata(17) = nteste
        ddata(18) = ntesti
        ddata(19) = track_teste
        ddata(20) = track_testi
        ddata(21) = scale_test
        call PBCAST(ddata, lenml)
	     qmetest = ddata(1)
	     qmitest = ddata(2)
        smoothIon = ddata(3)
        sqn(1,1) = ddata(4)
        sqn(2,1) = ddata(5)
        sqn(3,1) = ddata(6)
        sqn(4,1) = ddata(7)
        sqn(5,1) = ddata(8)
        sqn(6,1) = ddata(9)
        rsqn = ddata(10)
        qmsqn = ddata(11)
        rmsqn = ddata(12)
        check_period = ddata(13)
        fixedvel = ddata(14)
        rmteste = ddata(15)
        rmtesti = ddata(16)
        nteste = ddata(17)
        ntesti = ddata(18)
        track_teste = ddata(19)
        track_testi = ddata(20)
        scale_test = ddata(21)
    end subroutine sendnml_ie

END MODULE pinit32d_ie
