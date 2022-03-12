!*****************************************************************************************
!>
!  This program tests codes for the least-squares solution of
!  M nonlinear equations in n variables. it consists of a driver
!  And an interface subroutine fcn. the driver reads in data,
!  Calls the nonlinear least-squares solver, and finally prints
!  Out information on the performance of the solver. this is
!  Only a sample driver, many other drivers are possible. the
!  Interface subroutine fcn is necessary to take into account the
!  Forms of calling sequences used by the function and jacobian
!  Subroutines in the various nonlinear least-squares solvers.

program test

    use minpack_module
    use iso_fortran_env, only: output_unit

    implicit none

    ! originally from file22
    integer,parameter :: ncases = 28
    integer,dimension(ncases),parameter :: nprobs  = [1,1,2,2,3,3,4,5,6,7,8,9,10,11,11,11,12,13,14,15,15,15,15,16,16,16,17,18]
    integer,dimension(ncases),parameter :: ns      = [5,5,5,5,5,5,2,3,4,2,3,4,3,6,9,12,3,2,4,1,8,9,10,10,30,40,5,11]
    integer,dimension(ncases),parameter :: ms      = [10,50,10,50,10,50,2,3,4,2,15,11,16,31,31,31,10,10,20,8,8,9,10,10,30,40,33,65]
    integer,dimension(ncases),parameter :: ntriess = [1,1,1,1,1,1,3,3,3,3,3,3,2,3,3,3,1,1,3,3,1,1,1,3,1,1,1,1]

    integer :: i, ic, info, k, ldfjac, lwa, m, n, NFEv, NJEv, NPRob, ntries, icase
    real(wp) :: factor, fnorm1, fnorm2
    integer :: ma(60), na(60), nf(60), nj(60), np(60), nx(60)
    real(wp) :: fnm(60)
    integer,dimension(:),allocatable :: iwa
    real(wp),dimension(:),allocatable :: fvec, wa, x
    real(wp),dimension(:,:),allocatable :: fjac

    integer, parameter :: nwrite = output_unit ! logical output unit
    real(wp), parameter :: one = 1.0_wp
    real(wp), parameter :: ten = 10.0_wp
    real(wp), parameter :: tol = sqrt(dpmpar(1))
    real(wp), parameter :: solution_tol = 1.0e-4_wp !! reltol for matching previously generated solutions

    ic = 0
    do icase = 1, ncases+1

        if (icase == ncases+1) then
            write (nwrite, '(A,I3,A/)') '1SUMMARY OF ', ic, ' CALLS TO LMDER1'
            write (nwrite, '(A/)')      ' NPROB   N    M   NFEV  NJEV  INFO  FINAL L2 NORM'
            do i = 1, ic
                write (nwrite, '(3I5,3I6,1X,D15.7)') np(i), na(i), ma(i), nf(i), nj(i), nx(i), fnm(i)
            end do
            stop
        else

            nprob = nprobs(icase)
            n = ns(icase)
            m = ms(icase)
            lwa = 5*n+m
            ldfjac = m

            if (allocated(fjac)) deallocate(fjac); allocate(fjac(m,n))
            if (allocated(fvec)) deallocate(fvec); allocate(fvec(m))
            if (allocated(wa))   deallocate(wa);   allocate(wa(lwa))
            if (allocated(x))    deallocate(x);    allocate(x(n))
            if (allocated(iwa))  deallocate(iwa);  allocate(iwa(n))

            ntries = ntriess(icase)

            factor = one
            do k = 1, ntries
                ic = ic + 1
                call initpt(n, x, NPRob, factor)
                call ssqfcn(m, n, x, fvec, NPRob)
                fnorm1 = enorm(m, fvec)
                write (nwrite, '(////5X,A,I5,5X,A,2I5,5X//)') ' PROBLEM', NPRob, ' DIMENSIONS', n, m
                NFEv = 0
                NJEv = 0
                call lmder1(fcn, m, n, x, fvec, fjac, ldfjac, tol, info, iwa, wa, lwa)
                call ssqfcn(m, n, x, fvec, NPRob)
                fnorm2 = enorm(m, fvec)
                np(ic) = NPRob
                na(ic) = n
                ma(ic) = m
                nf(ic) = NFEv
                nj(ic) = NJEv
                nx(ic) = info
                fnm(ic) = fnorm2
                write (nwrite, '(5X,A,D15.7//5X,A,D15.7//5X,A,I10//5X,A,I10//5X,A,18X,I10//5X,A//,*(5X,5D15.7/))') &
                            ' INITIAL L2 NORM OF THE RESIDUALS', fnorm1, &
                            ' FINAL L2 NORM OF THE RESIDUALS  ', fnorm2, &
                            ' NUMBER OF FUNCTION EVALUATIONS  ', NFEv, &
                            ' NUMBER OF JACOBIAN EVALUATIONS  ', NJEv, &
                            ' EXIT PARAMETER', info, &
                            ' FINAL APPROXIMATE SOLUTION',x(1:n)
                factor = ten*factor

                ! compare with previously generated solutions:
                if (any(abs((solution(ic) - x)/(solution(ic))) > solution_tol)) then
                    write(nwrite,'(A)') 'Failed case'
                    write(nwrite, '(//5x, a//(5x, 5d15.7))') 'Expected x: ', solution(ic)
                    write(nwrite, '(/5x, a//(5x, 5d15.7))')  'Computed x: ', x
                    error stop
                end if

            end do
        end if
    end do

contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  The calling sequence of fcn should be identical to the
!  calling sequence of the function subroutine in the nonlinear
!  least-squares solver. fcn should only call the testing
!  function and jacobian subroutines ssqfcn and ssqjac with
!  the appropriate value of problem number (nprob).

    subroutine fcn(m, n, x, Fvec, Fjac, Ldfjac, Iflag)

        implicit none

        integer, intent(in) :: m !! the number of functions.
        integer, intent(in) :: n !! the number of variables.
        integer, intent(in) :: ldfjac !! leading dimension of the array fjac.
        integer, intent(inout) :: iflag !! if iflag = 1 calculate the functions at x and
                                       !! return this vector in fvec. do not alter fjac.
                                       !! if iflag = 2 calculate the jacobian at x and
                                       !! return this matrix in fjac. do not alter fvec.
                                       !!
                                       !! the value of iflag should not be changed by fcn unless
                                       !! the user wants to terminate execution of lmder.
                                       !! in this case set iflag to a negative integer.
        real(wp), intent(in) :: x(n) !! independant variable vector
        real(wp), intent(inout) :: fvec(m) !! value of function at `x`
        real(wp), intent(inout) :: fjac(ldfjac, n) !! jacobian matrix at `x`

        select case (iflag)
        case (1)
            call ssqfcn(m, n, x, Fvec, NPRob)
            NFEv = NFEv + 1
        case (2)
            call ssqjac(m, n, x, Fjac, Ldfjac, NPRob)
            NJEv = NJEv + 1
        case default
            error stop 'invalid iflag value'
        end select

    end subroutine fcn
!*****************************************************************************************

!*****************************************************************************************
!>
!  Replaced statement function in original code.

    pure elemental function dfloat(i) result(f)
        implicit none
        integer, intent(in) :: i
        real(wp) :: f
        f = real(i, wp)
    end function dfloat
!*****************************************************************************************

!*****************************************************************************************
!>
!  Get expected `x` vectors for each case.

    pure function solution(nprob) result(x)

        implicit none

        integer,intent(in) :: nprob
        real(wp),dimension(:),allocatable :: x

        select case (nprob)
        case(1); x = [-0.1000000000000000E+01_wp,-0.1000000000000000E+01_wp,-0.1000000000000000E+01_wp,&
                      -0.1000000000000000E+01_wp,-0.1000000000000000E+01_wp]
        case(2); x = [-0.9999999999999953E+00_wp,-0.1000000000000005E+01_wp,-0.9999999999999976E+00_wp,&
                      -0.9999999999999956E+00_wp,-0.9999999999999991E+00_wp]
        case(3); x = [-0.1677968180239693E+03_wp,-0.8339840901198468E+02_wp,0.2211100430795781E+03_wp,&
                      -0.4119920450599233E+02_wp,-0.3275936360479385E+02_wp]
        case(4); x = [-0.2029999900022674E+02_wp,-0.9649999500113370E+01_wp,-0.1652451975264496E+03_wp,&
                      -0.4324999750056676E+01_wp,0.1105330585100652E+03_wp]
        case(5); x = [0.1000000000000000E+01_wp,-0.2103615324224772E+03_wp,0.3212042081132130E+02_wp,&
                      0.8113456824980642E+02_wp,0.1000000000000000E+01_wp]
        case(6); x = [0.1000000000000000E+01_wp,0.3321494858957815E+03_wp,-0.4396851914289522E+03_wp,&
                      0.1636968825825863E+03_wp,0.1000000000000000E+01_wp]
        case(7); x = [0.1000000000000000E+01_wp,0.1000000000000000E+01_wp]
        case(8); x = [0.1000000000000000E+01_wp,0.1000000000000000E+01_wp]
        case(9); x = [0.1000000000000000E+01_wp,0.1000000000000000E+01_wp]
        case(10); x = [0.1000000000000000E+01_wp,-0.6243301596789443E-17_wp,0.0000000000000000E+00_wp]
        case(11); x = [0.1000000000000000E+01_wp,0.6563910805155555E-20_wp,0.0000000000000000E+00_wp]
        case(12); x = [0.1000000000000000E+01_wp,-0.1972152263052530E-29_wp,0.0000000000000000E+00_wp]
        case(13); x = [0.1652117596168394E-16_wp,-0.1652117596168393E-17_wp,0.2643388153869468E-17_wp,&
                       0.2643388153869468E-17_wp]
        case(14); x = [0.2016745112510229E-19_wp,-0.2016745112510229E-20_wp,0.3226792180016300E-20_wp,&
                       0.3226792180016300E-20_wp]
        case(15); x = [0.3226792180016378E-17_wp,-0.3226792180016378E-18_wp,0.5162867488026213E-18_wp,&
                       0.5162867488026213E-18_wp]
        case(16); x = [0.1141248446549937E+02_wp,-0.8968279137315035E+00_wp]
        case(17); x = [0.1141300466147456E+02_wp,-0.8967960386859591E+00_wp]
        case(18); x = [0.1141278178578820E+02_wp,-0.8968051074920677E+00_wp]
        case(19); x = [0.8241057657583339E-01_wp,0.1133036653471504E+01_wp,0.2343694638941154E+01_wp]
        case(20); x = [0.8406666738183293E+00_wp,-0.1588480332595655E+09_wp,-0.1643786716535352E+09_wp]
        case(21); x = [0.8406666738676455E+00_wp,-0.1589461672055184E+09_wp,-0.1644649068577712E+09_wp]
        case(22); x = [0.1928078104762493E+00_wp,0.1912626533540709E+00_wp,0.1230528010469309E+00_wp,&
                       0.1360532211505167E+00_wp]
        case(23); x = [0.7286754737686598E+06_wp,-0.1407588031293926E+02_wp,-0.3297779778419661E+08_wp,&
                       -0.2057159419780170E+08_wp]
        case(24); x = [0.1927984063846549E+00_wp,0.1914736844615448E+00_wp,0.1230924753714115E+00_wp,&
                       0.1361509629062244E+00_wp]
        case(25); x = [0.5609636471026654E-02_wp,0.6181346346286586E+04_wp,0.3452236346241439E+03_wp]
        case(26); x = [0.1291118238019656E-10_wp,0.3388670577625951E+05_wp,0.9040289597481549E+03_wp]
        case(27); x = [-0.1572496150837821E-01_wp,0.1012434882329655E+01_wp,-0.2329917223876743E+00_wp,&
                       0.1260431011028187E+01_wp,-0.1513730313944210E+01_wp,0.9929972729184218E+00_wp]
        case(28); x = [-0.1572519013866764E-01_wp,0.1012434858601050E+01_wp,-0.2329915458438263E+00_wp,&
                       0.1260429320891618E+01_wp,-0.1513727767065737E+01_wp,0.9929957342632760E+00_wp]
        case(29); x = [-0.1572470197125875E-01_wp,0.1012434909256583E+01_wp,-0.2329919227616435E+00_wp,&
                       0.1260432929295550E+01_wp,-0.1513733204527069E+01_wp,0.9929990192232205E+00_wp]
        case(30); x = [-0.1530706441663782E-04_wp,0.9997897039345963E+00_wp,0.1476396349113067E-01_wp,&
                       0.1463423301456550E+00_wp,0.1000821094549679E+01_wp,-0.2617731120708589E+01_wp,&
                       0.4104403139437982E+01_wp,-0.3143612262365286E+01_wp,0.1052626403788335E+01_wp]
        case(31); x = [-0.1530713348497491E-04_wp,0.9997897039412341E+00_wp,0.1476396297861486E-01_wp,&
                       0.1463423348188786E+00_wp,0.1000821073213774E+01_wp,-0.2617731070847226E+01_wp,&
                       0.4104403076555856E+01_wp,-0.3143612221787109E+01_wp,0.1052626393225985E+01_wp]
        case(32); x = [-0.1530703652788392E-04_wp,0.9997897039319498E+00_wp,0.1476396369347192E-01_wp,&
                       0.1463423283001695E+00_wp,0.1000821103001098E+01_wp,-0.2617731140510708E+01_wp,&
                       0.4104403164468549E+01_wp,-0.3143612278549893E+01_wp,0.1052626408008490E+01_wp]
        case(33); x = [-0.6638060467011595E-08_wp,0.1000001644117862E+01_wp,-0.5639322103223269E-03_wp,&
                       0.3478205405036302E+00_wp,-0.1567315040922755E+00_wp,0.1052815177186520E+01_wp,&
                       -0.3247271153392391E+01_wp,0.7288434897798658E+01_wp,-0.1027184824113737E+02_wp,&
                       0.9074113646924774E+01_wp,-0.4541375466623413E+01_wp,0.1012011888539188E+01_wp]
        case(34); x = [-0.6637102237178952E-08_wp,0.1000001644117874E+01_wp,-0.5639322084840861E-03_wp,&
                       0.3478205404902133E+00_wp,-0.1567315039874129E+00_wp,0.1052815176715104E+01_wp,&
                       -0.3247271152062481E+01_wp,0.7288434895390663E+01_wp,-0.1027184823833726E+02_wp,&
                       0.9074113644905919E+01_wp,-0.4541375465802633E+01_wp,0.1012011888395745E+01_wp]
        case(35); x = [-0.6638060465675741E-08_wp,0.1000001644117862E+01_wp,-0.5639322102727471E-03_wp,&
                       0.3478205405024820E+00_wp,-0.1567315040810628E+00_wp,0.1052815177127253E+01_wp,&
                       -0.3247271153204238E+01_wp,0.7288434897423372E+01_wp,-0.1027184824066345E+02_wp,&
                       0.9074113646557050E+01_wp,-0.4541375466463498E+01_wp,0.1012011888509363E+01_wp]
        case(36); x = [0.9999999999999999E+00_wp,0.1000000000000000E+02_wp,0.1000000000000000E+01_wp]
        case(37); x = [0.2578199266367936E+00_wp,0.2578299767645596E+00_wp]
        case(38); x = [-0.1159125521567487E+02_wp,0.1320248976116189E+02_wp,-0.4035744915102742E+00_wp,&
                       0.2367363118370096E+00_wp]
        case(39); x = [-0.1159592742721295E+02_wp,0.1320418669262178E+02_wp,-0.4034173628612446E+00_wp,&
                       0.2367711433960093E+00_wp]
        case(40); x = [-0.1159025975515195E+02_wp,0.1320206290836587E+02_wp,-0.4036882104139269E+00_wp,&
                       0.2366649366417166E+00_wp]
        case(41); x = [0.5000000000000000E+00_wp]
        case(42); x = [0.9817314924683995E+00_wp]
        case(43); x = [0.9817314852933997E+00_wp]
        case(44); x = [0.4315366485873882E-01_wp,0.1930916378432737E+00_wp,0.2663285938126948E+00_wp,&
                       0.4999993346289514E+00_wp,0.5000006653710486E+00_wp,0.7336714061873052E+00_wp,&
                       0.8069083621567262E+00_wp,0.9568463351412612E+00_wp]
        case(45); x = [0.4420534613578274E-01_wp,0.1994906723098809E+00_wp,0.2356191084710600E+00_wp,&
                       0.4160469078925979E+00_wp,0.5000000000000001E+00_wp,0.5839530921074019E+00_wp,&
                       0.7643808915289401E+00_wp,0.8005093276901190E+00_wp,0.9557946538642172E+00_wp]
        case(46); x = [0.5962026717536065E-01_wp,0.1667087838059440E+00_wp,0.2391710188135143E+00_wp,&
                       0.3988852903458887E+00_wp,0.3988836678710648E+00_wp,0.6011163321289351E+00_wp,&
                       0.6011147096541113E+00_wp,0.7608289811864857E+00_wp,0.8332912161940560E+00_wp,&
                       0.9403797328246394E+00_wp]
        case(47); x = [0.9794303033498624E+00_wp,0.9794303033498624E+00_wp,0.9794303033498624E+00_wp,&
                       0.9794303033498624E+00_wp,0.9794303033498624E+00_wp,0.9794303033498624E+00_wp,&
                       0.9794303033498624E+00_wp,0.9794303033498624E+00_wp,0.9794303033498624E+00_wp,&
                       0.1205696966501376E+01_wp]
        case(48); x = [0.9794303033498634E+00_wp,0.9794303033498634E+00_wp,0.9794303033498634E+00_wp,&
                       0.9794303033498634E+00_wp,0.9794303033498634E+00_wp,0.9794303033498634E+00_wp,&
                       0.9794303033498634E+00_wp,0.9794303033498634E+00_wp,0.9794303033498634E+00_wp,&
                       0.1205696966501366E+01_wp]
        case(49); x = [0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,&
                       0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,&
                       0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,&
                       0.9999999999999840E+00_wp]
        case(50); x = [0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,&
                       0.9977542164428221E+00_wp,0.9977542164428221E+00_wp,0.1067373506715322E+01_wp]
        case(51); x = [0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,0.1000000000000006E+01_wp,&
                       0.9999999999997748E+00_wp]
        case(52); x = [0.3754100492440257E+00_wp,0.1935846545431133E+01_wp,-0.1464686767487213E+01_wp,&
                       0.1286753391104403E-01_wp,0.2212270118130737E-01_wp]
        case(53); x = [0.1309976638100963E+01_wp,0.4315524807599996E+00_wp,0.6336612616028594E+00_wp,&
                       0.5994285609916951E+00_wp,0.7541797682724486E+00_wp,0.9043000823785187E+00_wp,&
                       0.1365799495210074E+01_wp,0.4823731997481072E+01_wp,0.2398684751048711E+01_wp,&
                       0.4568875547914517E+01_wp,0.5675342062730520E+01_wp]
        case default
            error stop 'invalid case'
        end select

    end function solution
!*****************************************************************************************

!*****************************************************************************************
!
!  This subroutine defines the jacobian matrices of eighteen
!  nonlinear least squares problems. The problem dimensions are
!  as described in the prologue comments of [[ssqfcn]].

    subroutine ssqjac(m, n, x, Fjac, Ldfjac, Nprob)
        implicit none

        integer, intent(in) :: m !! positive integer input variables. n must not exceed m.
        integer, intent(in) :: n !! positive integer input variables. n must not exceed m.
        integer, intent(in) :: ldfjac !! a positive integer input variable not less than m
                                   !! which specifies the leading dimension of the array fjac.
        integer, intent(in) :: nprob !! a positive integer variable which defines the
                                  !! number of the problem. nprob must not exceed 18.
        real(wp), intent(in) :: x(n) !! an input array of length n.
        real(wp), intent(out) :: fjac(ldfjac, n) !! an m by n output array which contains the jacobian
                                             !! matrix of the nprob function evaluated at x.

        real(wp), parameter :: v(11) = [4.0_wp, 2.0_wp, 1.0_wp, 5.0e-1_wp, &
                                        2.5e-1_wp, 1.67e-1_wp, 1.25e-1_wp, 1.0e-1_wp, &
                                        8.33e-2_wp, 7.14e-2_wp, 6.25e-2_wp]

        real(wp), parameter :: zero = 0.0_wp
        real(wp), parameter :: one = 1.0_wp
        real(wp), parameter :: two = 2.0_wp
        real(wp), parameter :: three = 3.0_wp
        real(wp), parameter :: four = 4.0_wp
        real(wp), parameter :: five = 5.0_wp
        real(wp), parameter :: eight = 8.0_wp
        real(wp), parameter :: ten = 10.0_wp
        real(wp), parameter :: c14 = 14.0_wp
        real(wp), parameter :: c20 = 20.0_wp
        real(wp), parameter :: c29 = 29.0_wp
        real(wp), parameter :: c45 = 45.0_wp
        real(wp), parameter :: c100 = 100.0_wp

        integer :: i, ivar, j, k, mm1, nm1
        real(wp) :: div, dx, prod, s2, temp, ti, tmp1, tmp2, tmp3, tmp4, tpi

        Fjac(1:m, 1:n) = zero

        ! JACOBIAN ROUTINE SELECTOR.

        select case (Nprob)
        case (2)

            ! LINEAR FUNCTION - RANK 1.

            do j = 1, n
                do i = 1, m
                    Fjac(i, j) = dfloat(i)*dfloat(j)
                end do
            end do

        case (3)

            ! LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.

            do j = 1, n
                do i = 1, m
                    Fjac(i, j) = zero
                end do
            end do
            nm1 = n - 1
            mm1 = m - 1
            if (nm1 >= 2) then
                do j = 2, nm1
                    do i = 2, mm1
                        Fjac(i, j) = dfloat(i - 1)*dfloat(j)
                    end do
                end do
            end if

        case (4)

            ! ROSENBROCK FUNCTION.

            Fjac(1, 1) = -c20*x(1)
            Fjac(1, 2) = ten
            Fjac(2, 1) = -one
            Fjac(2, 2) = zero

        case (5)

            ! HELICAL VALLEY FUNCTION.

            tpi = eight*atan(one)
            temp = x(1)**2 + x(2)**2
            tmp1 = tpi*temp
            tmp2 = sqrt(temp)
            Fjac(1, 1) = c100*x(2)/tmp1
            Fjac(1, 2) = -c100*x(1)/tmp1
            Fjac(1, 3) = ten
            Fjac(2, 1) = ten*x(1)/tmp2
            Fjac(2, 2) = ten*x(2)/tmp2
            Fjac(2, 3) = zero
            Fjac(3, 1) = zero
            Fjac(3, 2) = zero
            Fjac(3, 3) = one

        case (6)

            ! POWELL SINGULAR FUNCTION.

            do j = 1, 4
                do i = 1, 4
                    Fjac(i, j) = zero
                end do
            end do
            Fjac(1, 1) = one
            Fjac(1, 2) = ten
            Fjac(2, 3) = sqrt(five)
            Fjac(2, 4) = -Fjac(2, 3)
            Fjac(3, 2) = two*(x(2) - two*x(3))
            Fjac(3, 3) = -two*Fjac(3, 2)
            Fjac(4, 1) = two*sqrt(ten)*(x(1) - x(4))
            Fjac(4, 4) = -Fjac(4, 1)

        case (7)

            ! FREUDENSTEIN AND ROTH FUNCTION.

            Fjac(1, 1) = one
            Fjac(1, 2) = x(2)*(ten - three*x(2)) - two
            Fjac(2, 1) = one
            Fjac(2, 2) = x(2)*(two + three*x(2)) - c14

        case (8)

            ! BARD FUNCTION.

            do i = 1, 15
                tmp1 = dfloat(i)
                tmp2 = dfloat(16 - i)
                tmp3 = tmp1
                if (i > 8) tmp3 = tmp2
                tmp4 = (x(2)*tmp2 + x(3)*tmp3)**2
                Fjac(i, 1) = -one
                Fjac(i, 2) = tmp1*tmp2/tmp4
                Fjac(i, 3) = tmp1*tmp3/tmp4
            end do

        case (9)

            ! KOWALIK AND OSBORNE FUNCTION.

            do i = 1, 11
                tmp1 = v(i)*(v(i) + x(2))
                tmp2 = v(i)*(v(i) + x(3)) + x(4)
                Fjac(i, 1) = -tmp1/tmp2
                Fjac(i, 2) = -v(i)*x(1)/tmp2
                Fjac(i, 3) = Fjac(i, 1)*Fjac(i, 2)
                Fjac(i, 4) = Fjac(i, 3)/v(i)
            end do

        case (10)

            ! MEYER FUNCTION.

            do i = 1, 16
                temp = five*dfloat(i) + c45 + x(3)
                tmp1 = x(2)/temp
                tmp2 = exp(tmp1)
                Fjac(i, 1) = tmp2
                Fjac(i, 2) = x(1)*tmp2/temp
                Fjac(i, 3) = -tmp1*Fjac(i, 2)
            end do

        case (11)

            ! WATSON FUNCTION.

            do i = 1, 29
                div = dfloat(i)/c29
                s2 = zero
                dx = one
                do j = 1, n
                    s2 = s2 + dx*x(j)
                    dx = div*dx
                end do
                temp = two*div*s2
                dx = one/div
                do j = 1, n
                    Fjac(i, j) = dx*(dfloat(j - 1) - temp)
                    dx = div*dx
                end do
            end do
            do j = 1, n
                do i = 30, 31
                    Fjac(i, j) = zero
                end do
            end do
            Fjac(30, 1) = one
            Fjac(31, 1) = -two*x(1)
            Fjac(31, 2) = one

        case (12)

            ! BOX 3-DIMENSIONAL FUNCTION.

            do i = 1, m
                temp = dfloat(i)
                tmp1 = temp/ten
                Fjac(i, 1) = -tmp1*exp(-tmp1*x(1))
                Fjac(i, 2) = tmp1*exp(-tmp1*x(2))
                Fjac(i, 3) = exp(-temp) - exp(-tmp1)
            end do

        case (13)

            ! JENNRICH AND SAMPSON FUNCTION.

            do i = 1, m
                temp = dfloat(i)
                Fjac(i, 1) = -temp*exp(temp*x(1))
                Fjac(i, 2) = -temp*exp(temp*x(2))
            end do

        case (14)

            ! BROWN AND DENNIS FUNCTION.

            do i = 1, m
                temp = dfloat(i)/five
                ti = sin(temp)
                tmp1 = x(1) + temp*x(2) - exp(temp)
                tmp2 = x(3) + ti*x(4) - cos(temp)
                Fjac(i, 1) = two*tmp1
                Fjac(i, 2) = temp*Fjac(i, 1)
                Fjac(i, 3) = two*tmp2
                Fjac(i, 4) = ti*Fjac(i, 3)
            end do

        case (15)

            ! CHEBYQUAD FUNCTION.

            dx = one/dfloat(n)
            do j = 1, n
                tmp1 = one
                tmp2 = two*x(j) - one
                temp = two*tmp2
                tmp3 = zero
                tmp4 = two
                do i = 1, m
                    Fjac(i, j) = dx*tmp4
                    ti = four*tmp2 + temp*tmp4 - tmp3
                    tmp3 = tmp4
                    tmp4 = ti
                    ti = temp*tmp2 - tmp1
                    tmp1 = tmp2
                    tmp2 = ti
                end do
            end do

        case (16)

            ! BROWN ALMOST-LINEAR FUNCTION.

            prod = one
            do j = 1, n
                prod = x(j)*prod
                do i = 1, n
                    Fjac(i, j) = one
                end do
                Fjac(j, j) = two
            end do
            do j = 1, n
                temp = x(j)
                if (temp == zero) then
                    temp = one
                    prod = one
                    do k = 1, n
                        if (k /= j) prod = x(k)*prod
                    end do
                end if
                Fjac(n, j) = prod/temp
            end do

        case (17)

            ! OSBORNE 1 FUNCTION.

            do i = 1, 33
                temp = ten*dfloat(i - 1)
                tmp1 = exp(-x(4)*temp)
                tmp2 = exp(-x(5)*temp)
                Fjac(i, 1) = -one
                Fjac(i, 2) = -tmp1
                Fjac(i, 3) = -tmp2
                Fjac(i, 4) = temp*x(2)*tmp1
                Fjac(i, 5) = temp*x(3)*tmp2
            end do

        case (18)

            ! OSBORNE 2 FUNCTION.

            do i = 1, 65
                temp = dfloat(i - 1)/ten
                tmp1 = exp(-x(5)*temp)
                tmp2 = exp(-x(6)*(temp - x(9))**2)
                tmp3 = exp(-x(7)*(temp - x(10))**2)
                tmp4 = exp(-x(8)*(temp - x(11))**2)
                Fjac(i, 1) = -tmp1
                Fjac(i, 2) = -tmp2
                Fjac(i, 3) = -tmp3
                Fjac(i, 4) = -tmp4
                Fjac(i, 5) = temp*x(1)*tmp1
                Fjac(i, 6) = x(2)*(temp - x(9))**2*tmp2
                Fjac(i, 7) = x(3)*(temp - x(10))**2*tmp3
                Fjac(i, 8) = x(4)*(temp - x(11))**2*tmp4
                Fjac(i, 9) = -two*x(2)*x(6)*(temp - x(9))*tmp2
                Fjac(i, 10) = -two*x(3)*x(7)*(temp - x(10))*tmp3
                Fjac(i, 11) = -two*x(4)*x(8)*(temp - x(11))*tmp4
            end do

        case default

            ! LINEAR FUNCTION - FULL RANK.

            temp = two/dfloat(m)
            do j = 1, n
                do i = 1, m
                    Fjac(i, j) = -temp
                end do
                Fjac(j, j) = Fjac(j, j) + one
            end do

        end select

    end subroutine ssqjac
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine specifies the standard starting points for the
!  functions defined by subroutine [[ssqfcn]]. the subroutine returns
!  in x a multiple (factor) of the standard starting point. for
!  the 11th function the standard starting point is zero, so in
!  this case, if factor is not unity, then the subroutine returns
!  the vector  x(j) = factor, j=1,...,n.

    subroutine initpt(n, x, Nprob, Factor)
        implicit none

        integer, intent(in) :: n !! a positive integer input variable.
        integer, intent(in) :: nprob !! a positive integer input variable which defines the
                                  !! number of the problem. nprob must not exceed 18.
        real(wp), intent(in) :: factor !! an input variable which specifies the multiple of
                                    !! the standard starting point. if factor is unity, no
                                    !! multiplication is performed.
        real(wp), intent(out) :: x(n) !! an output array of length n which contains the standard
                                   !! starting point for problem nprob multiplied by factor.

        real(wp), parameter :: zero = 0.0_wp
        real(wp), parameter :: half = 0.5_wp
        real(wp), parameter :: one = 1.0_wp
        real(wp), parameter :: two = 2.0_wp
        real(wp), parameter :: three = 3.0_wp
        real(wp), parameter :: five = 5.0_wp
        real(wp), parameter :: seven = 7.0_wp
        real(wp), parameter :: ten = 10.0_wp
        real(wp), parameter :: twenty = 20.0_wp
        real(wp), parameter :: twntf = 25.0_wp

        real(wp), parameter :: c1 = 1.2_wp
        real(wp), parameter :: c2 = 2.5e-1_wp
        real(wp), parameter :: c3 = 3.9e-1_wp
        real(wp), parameter :: c4 = 4.15e-1_wp
        real(wp), parameter :: c5 = 2.0e-2_wp
        real(wp), parameter :: c6 = 4.0e3_wp
        real(wp), parameter :: c7 = 2.5e2_wp
        real(wp), parameter :: c8 = 3.0e-1_wp
        real(wp), parameter :: c9 = 4.0e-1_wp
        real(wp), parameter :: c10 = 1.5_wp
        real(wp), parameter :: c11 = 1.0e-2_wp
        real(wp), parameter :: c12 = 1.3_wp
        real(wp), parameter :: c13 = 6.5e-1_wp
        real(wp), parameter :: c14 = 7.0e-1_wp
        real(wp), parameter :: c15 = 6.0e-1_wp
        real(wp), parameter :: c16 = 4.5_wp
        real(wp), parameter :: c17 = 5.5_wp

        integer :: ivar, j
        real(wp) :: h

        x(1:n) = zero

        ! SELECTION OF INITIAL POINT.

        select case (Nprob)
        case (4)

            ! ROSENBROCK FUNCTION.

            x(1) = -c1
            x(2) = one

        case (5)

            ! HELICAL VALLEY FUNCTION.

            x(1) = -one
            x(2) = zero
            x(3) = zero

        case (6)

            ! POWELL SINGULAR FUNCTION.

            x(1) = three
            x(2) = -one
            x(3) = zero
            x(4) = one

        case (7)

            ! FREUDENSTEIN AND ROTH FUNCTION.

            x(1) = half
            x(2) = -two

        case (8)

            ! BARD FUNCTION.

            x(1) = one
            x(2) = one
            x(3) = one

        case (9)

            ! KOWALIK AND OSBORNE FUNCTION.

            x(1) = c2
            x(2) = c3
            x(3) = c4
            x(4) = c3

        case (10)

            ! MEYER FUNCTION.

            x(1) = c5
            x(2) = c6
            x(3) = c7

        case (11)

            ! WATSON FUNCTION.

            do j = 1, n
                x(j) = zero
            end do

        case (12)

            ! BOX 3-DIMENSIONAL FUNCTION.

            x(1) = zero
            x(2) = ten
            x(3) = twenty

        case (13)

            ! JENNRICH AND SAMPSON FUNCTION.

            x(1) = c8
            x(2) = c9

        case (14)

            ! BROWN AND DENNIS FUNCTION.

            x(1) = twntf
            x(2) = five
            x(3) = -five
            x(4) = -one

        case (15)

            ! CHEBYQUAD FUNCTION.

            h = one/dfloat(n + 1)
            do j = 1, n
                x(j) = dfloat(j)*h
            end do

        case (16)

            ! BROWN ALMOST-LINEAR FUNCTION.

            do j = 1, n
                x(j) = half
            end do

        case (17)

            ! OSBORNE 1 FUNCTION.

            x(1) = half
            x(2) = c10
            x(3) = -one
            x(4) = c11
            x(5) = c5

        case (18)

            ! OSBORNE 2 FUNCTION.

            x(1) = c12
            x(2) = c13
            x(3) = c13
            x(4) = c14
            x(5) = c15
            x(6) = three
            x(7) = five
            x(8) = seven
            x(9) = two
            x(10) = c16
            x(11) = c17

        case default

            ! LINEAR FUNCTION - FULL RANK OR RANK 1.

            do j = 1, n
                x(j) = one
            end do

        end select

        ! COMPUTE MULTIPLE OF INITIAL POINT.

        if (Factor /= one) then
            if (Nprob == 11) then
                do j = 1, n
                    x(j) = Factor
                end do
            else
                do j = 1, n
                    x(j) = Factor*x(j)
                end do
            end if
        end if

    end subroutine initpt
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine defines the functions of eighteen nonlinear
!  least squares problems. the allowable values of (m,n) for
!  functions 1,2 and 3 are variable but with m >= n.
!  for functions 4,5,6,7,8,9 and 10 the values of (m,n) are
!  (2,2),(3,3),(4,4),(2,2),(15,3),(11,4) and (16,3), respectively.
!  function 11 (watson) has m = 31 with n usually 6 or 9.
!  however, any n, n = 2,...,31, is permitted.
!  functions 12,13 and 14 have n = 3,2 and 4, respectively, but
!  allow any m >= n, with the usual choices being 10,10 and 20.
!  function 15 (chebyquad) allows m and n variable with m >= n.
!  function 16 (brown) allows n variable with m = n.
!  for functions 17 and 18, the values of (m,n) are
!  (33,5) and (65,11), respectively.

    subroutine ssqfcn(m, n, x, Fvec, Nprob)

        implicit none

        integer,intent(in) :: m !! positive integer input variable. n must not exceed m.
        integer,intent(in) :: n !! positive integer input variable. n must not exceed m.
        integer,intent(in) :: nprob !! a positive integer input variable which defines the
                                    !! number of the problem. nprob must not exceed 18.
        real(wp),intent(in) :: x(n) !! an input array of length n.
        real(wp),intent(out) :: fvec(m) !! an output array of length m which contains the nprob
                                        !! function evaluated at x.

        real(wp),parameter :: zero  = 0.0_wp
        real(wp),parameter :: zp25  = 0.25_wp
        real(wp),parameter :: zp5   = 0.5_wp
        real(wp),parameter :: one   = 1.0_wp
        real(wp),parameter :: two   = 2.0_wp
        real(wp),parameter :: five  = 5.0_wp
        real(wp),parameter :: eight = 8.0_wp
        real(wp),parameter :: ten   = 10.0_wp
        real(wp),parameter :: c13   = 13.0_wp
        real(wp),parameter :: c14   = 14.0_wp
        real(wp),parameter :: c29   = 29.0_wp
        real(wp),parameter :: c45   = 45.0_wp

        integer :: i, iev, ivar, j, nm1
        real(wp) :: div, dx, prod, sum, s1, s2, temp, ti, &
                    tmp1, tmp2, tmp3, tmp4, tpi

        real(wp),parameter :: v(11) = [ &
             4.0_wp, 2.0_wp, 1.0_wp, 5.0e-1_wp, 2.5e-1_wp, 1.67e-1_wp, &
             1.25e-1_wp, 1.0e-1_wp, 8.33e-2_wp, 7.14e-2_wp, 6.25e-2_wp]
        real(wp),parameter :: y1(15) = [ &
             1.4e-1_wp, 1.8e-1_wp, 2.2e-1_wp, 2.5e-1_wp, 2.9e-1_wp, 3.2e-1_wp, &
             3.5e-1_wp, 3.9e-1_wp, 3.7e-1_wp, 5.8e-1_wp, 7.3e-1_wp, 9.6e-1_wp, &
             1.34_wp, 2.1_wp, 4.39_wp]
        real(wp),parameter :: y2(11) = [1.957e-1_wp, 1.947e-1_wp, &
             1.735e-1_wp, 1.6e-1_wp, 8.44e-2_wp, 6.27e-2_wp, 4.56e-2_wp, 3.42e-2_wp,  &
             3.23e-2_wp, 2.35e-2_wp, 2.46e-2_wp ]
        real(wp),parameter :: y3(16) = [ &
             3.478e4_wp, 2.861e4_wp, 2.365e4_wp, 1.963e4_wp, &
             1.637e4_wp, 1.372e4_wp, 1.154e4_wp, 9.744e3_wp, 8.261e3_wp, 7.03e3_wp,   &
             6.005e3_wp, 5.147e3_wp, 4.427e3_wp, 3.82e3_wp, 3.307e3_wp, 2.872e3_wp ]
        real(wp),parameter :: y4(33) = [ 8.44e-1_wp,&
             9.08e-1_wp, 9.32e-1_wp, 9.36e-1_wp, 9.25e-1_wp, 9.08e-1_wp, 8.81e-1_wp,  &
             8.5e-1_wp, 8.18e-1_wp, 7.84e-1_wp, 7.51e-1_wp, 7.18e-1_wp, 6.85e-1_wp,   &
             6.58e-1_wp, 6.28e-1_wp, 6.03e-1_wp, 5.8e-1_wp, 5.58e-1_wp, 5.38e-1_wp,   &
             5.22e-1_wp, 5.06e-1_wp, 4.9e-1_wp, 4.78e-1_wp, 4.67e-1_wp, 4.57e-1_wp,   &
             4.48e-1_wp, 4.38e-1_wp, 4.31e-1_wp, 4.24e-1_wp, 4.2e-1_wp, 4.14e-1_wp,   &
             4.11e-1_wp, 4.06e-1_wp ]
        real(wp),parameter :: y5(65) = [ 1.366_wp, &
             1.191_wp, 1.112_wp, 1.013_wp, 9.91e-1_wp, 8.85e-1_wp, 8.31e-1_wp,  &
             8.47e-1_wp, 7.86e-1_wp, 7.25e-1_wp, 7.46e-1_wp, 6.79e-1_wp, 6.08e-1_wp,  &
             6.55e-1_wp, 6.16e-1_wp, 6.06e-1_wp, 6.02e-1_wp, 6.26e-1_wp, 6.51e-1_wp,  &
             7.24e-1_wp, 6.49e-1_wp, 6.49e-1_wp, 6.94e-1_wp, 6.44e-1_wp, 6.24e-1_wp,  &
             6.61e-1_wp, 6.12e-1_wp, 5.58e-1_wp, 5.33e-1_wp, 4.95e-1_wp, 5.0e-1_wp,   &
             4.23e-1_wp, 3.95e-1_wp, 3.75e-1_wp, 3.72e-1_wp, 3.91e-1_wp, 3.96e-1_wp,  &
             4.05e-1_wp, 4.28e-1_wp, 4.29e-1_wp, 5.23e-1_wp, 5.62e-1_wp, 6.07e-1_wp,  &
             6.53e-1_wp, 6.72e-1_wp, 7.08e-1_wp, 6.33e-1_wp, 6.68e-1_wp, 6.45e-1_wp,  &
             6.32e-1_wp, 5.91e-1_wp, 5.59e-1_wp, 5.97e-1_wp, 6.25e-1_wp, 7.39e-1_wp,  &
             7.1e-1_wp, 7.29e-1_wp, 7.2e-1_wp, 6.36e-1_wp, 5.81e-1_wp, 4.28e-1_wp,    &
             2.92e-1_wp, 1.62e-1_wp, 9.8e-2_wp, 5.4e-2_wp ]

        ! initialize:
        fvec(1:m) = zero

        ! FUNCTION ROUTINE SELECTOR.

        select case (Nprob)
        case (2)

            ! LINEAR FUNCTION - RANK 1.

            sum = zero
            do j = 1, n
                sum = sum + dfloat(j)*x(j)
            end do
            do i = 1, m
                Fvec(i) = dfloat(i)*sum - one
            end do
        case (3)

            ! LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.

            sum = zero
            nm1 = n - 1
            if (nm1 >= 2) then
                do j = 2, nm1
                    sum = sum + dfloat(j)*x(j)
                end do
            end if
            do i = 1, m
                Fvec(i) = dfloat(i - 1)*sum - one
            end do
            Fvec(m) = -one
        case (4)

            ! ROSENBROCK FUNCTION.

            Fvec(1) = ten*(x(2) - x(1)**2)
            Fvec(2) = one - x(1)
        case (5)

            ! HELICAL VALLEY FUNCTION.

            tpi = eight*atan(one)
            tmp1 = sign(zp25, x(2))
            if (x(1) > zero) tmp1 = atan(x(2)/x(1))/tpi
            if (x(1) < zero) tmp1 = atan(x(2)/x(1))/tpi + zp5
            tmp2 = sqrt(x(1)**2 + x(2)**2)
            Fvec(1) = ten*(x(3) - ten*tmp1)
            Fvec(2) = ten*(tmp2 - one)
            Fvec(3) = x(3)
        case (6)

            ! POWELL SINGULAR FUNCTION.

            Fvec(1) = x(1) + ten*x(2)
            Fvec(2) = sqrt(five)*(x(3) - x(4))
            Fvec(3) = (x(2) - two*x(3))**2
            Fvec(4) = sqrt(ten)*(x(1) - x(4))**2
        case (7)

            ! FREUDENSTEIN AND ROTH FUNCTION.

            Fvec(1) = -c13 + x(1) + ((five - x(2))*x(2) - two)*x(2)
            Fvec(2) = -c29 + x(1) + ((one + x(2))*x(2) - c14)*x(2)
        case (8)

            ! BARD FUNCTION.

            do i = 1, 15
                tmp1 = dfloat(i)
                tmp2 = dfloat(16 - i)
                tmp3 = tmp1
                if (i > 8) tmp3 = tmp2
                Fvec(i) = y1(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
            end do
        case (9)

            ! KOWALIK AND OSBORNE FUNCTION.

            do i = 1, 11
                tmp1 = v(i)*(v(i) + x(2))
                tmp2 = v(i)*(v(i) + x(3)) + x(4)
                Fvec(i) = y2(i) - x(1)*tmp1/tmp2
            end do
        case (10)

            ! MEYER FUNCTION.

            do i = 1, 16
                temp = five*dfloat(i) + c45 + x(3)
                tmp1 = x(2)/temp
                tmp2 = exp(tmp1)
                Fvec(i) = x(1)*tmp2 - y3(i)
            end do
        case (11)

            ! WATSON FUNCTION.

            do i = 1, 29
                div = dfloat(i)/c29
                s1 = zero
                dx = one
                do j = 2, n
                    s1 = s1 + dfloat(j - 1)*dx*x(j)
                    dx = div*dx
                end do
                s2 = zero
                dx = one
                do j = 1, n
                    s2 = s2 + dx*x(j)
                    dx = div*dx
                end do
                Fvec(i) = s1 - s2**2 - one
            end do
            Fvec(30) = x(1)
            Fvec(31) = x(2) - x(1)**2 - one
        case (12)

            ! BOX 3-DIMENSIONAL FUNCTION.

            do i = 1, m
                temp = dfloat(i)
                tmp1 = temp/ten
                Fvec(i) = exp(-tmp1*x(1)) - exp(-tmp1*x(2)) &
                        & + (exp(-temp) - exp(-tmp1))*x(3)
            end do
        case (13)

            ! JENNRICH AND SAMPSON FUNCTION.

            do i = 1, m
                temp = dfloat(i)
                Fvec(i) = two + two*temp - exp(temp*x(1)) - exp(temp*x(2))
            end do
        case (14)

            ! BROWN AND DENNIS FUNCTION.

            do i = 1, m
                temp = dfloat(i)/five
                tmp1 = x(1) + temp*x(2) - exp(temp)
                tmp2 = x(3) + sin(temp)*x(4) - cos(temp)
                Fvec(i) = tmp1**2 + tmp2**2
            end do
        case (15)

            ! CHEBYQUAD FUNCTION.

            do i = 1, m
                Fvec(i) = zero
            end do
            do j = 1, n
                tmp1 = one
                tmp2 = two*x(j) - one
                temp = two*tmp2
                do i = 1, m
                    Fvec(i) = Fvec(i) + tmp2
                    ti = temp*tmp2 - tmp1
                    tmp1 = tmp2
                    tmp2 = ti
                end do
            end do
            dx = one/dfloat(n)
            iev = -1
            do i = 1, m
                Fvec(i) = dx*Fvec(i)
                if (iev > 0) Fvec(i) = Fvec(i) + one/(dfloat(i)**2 - one)
                iev = -iev
            end do
        case (16)

            ! BROWN ALMOST-LINEAR FUNCTION.

            sum = -dfloat(n + 1)
            prod = one
            do j = 1, n
                sum = sum + x(j)
                prod = x(j)*prod
            end do
            do i = 1, n
                Fvec(i) = x(i) + sum
            end do
            Fvec(n) = prod - one
        case (17)

            ! OSBORNE 1 FUNCTION.

            do i = 1, 33
                temp = ten*dfloat(i - 1)
                tmp1 = exp(-x(4)*temp)
                tmp2 = exp(-x(5)*temp)
                Fvec(i) = y4(i) - (x(1) + x(2)*tmp1 + x(3)*tmp2)
            end do
        case (18)

            ! OSBORNE 2 FUNCTION.

            do i = 1, 65
                temp = dfloat(i - 1)/ten
                tmp1 = exp(-x(5)*temp)
                tmp2 = exp(-x(6)*(temp - x(9))**2)
                tmp3 = exp(-x(7)*(temp - x(10))**2)
                tmp4 = exp(-x(8)*(temp - x(11))**2)
                Fvec(i) = y5(i) - (x(1)*tmp1 + x(2)*tmp2 + x(3)*tmp3 + x(4)*tmp4)
            end do
        case default

            ! LINEAR FUNCTION - FULL RANK.

            sum = zero
            do j = 1, n
                sum = sum + x(j)
            end do
            temp = two*sum/dfloat(m) + one
            do i = 1, m
                Fvec(i) = -temp
                if (i <= n) Fvec(i) = Fvec(i) + x(i)
            end do
        end select

    end subroutine ssqfcn
!*****************************************************************************************

!*****************************************************************************************
    end program test
!*****************************************************************************************