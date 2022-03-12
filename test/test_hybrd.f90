
!*****************************************************************************************
!>
!  This program tests codes for the solution of n nonlinear
!  equations in n variables. it consists of a driver and an
!  interface subroutine fcn. the driver reads in data, calls the
!  nonlinear equation solver, and finally prints out information
!  on the performance of the solver. this is only a sample driver,
!  many other drivers are possible. the interface subroutine fcn
!  is necessary to take into account the forms of calling
!  sequences used by the function subroutines in the various
!  nonlinear equation solvers.

program test

    use minpack_module
    use iso_fortran_env, only: output_unit

    implicit none

    ! originally from file21
    integer,parameter :: ncases = 22
    integer,dimension(ncases),parameter :: nprobs = [1,2,3,4,5,6,6,7,7,7,7,7,8,8,8,9,10,10,11,12,13,14]
    integer,dimension(ncases),parameter :: ns = [2,4,2,4,3,6,9,5,6,7,8,9,10,30,40,10,1,10,10,10,10,10]
    integer,dimension(ncases),parameter :: ntriess = [3,3,2,3,3,2,2,3,3,3,1,1,3,1,1,3,3,3,3,3,3,3]

    integer :: i, ic, info, k, lwa, n, NFEv, NPRob, ntries, icase
    integer :: na(60), nf(60), np(60), nx(60)
    real(wp) :: fnm(60)
    real(wp) :: factor, fnorm1, fnorm2
    real(wp),allocatable :: fvec(:), wa(:), x(:)

    integer, parameter :: nwrite = output_unit !! logical output unit
    real(wp), parameter :: one = 1.0_wp
    real(wp), parameter :: ten = 10.0_wp
    real(wp), parameter :: tol = sqrt(dpmpar(1))
    real(wp), parameter :: solution_reltol = 1.0e-4_wp !! reltol for matching previously generated solutions

    ic = 0
    do icase = 1, ncases+1

        if (icase == ncases+1) then
            write (nwrite, '(A,I3,A/)') '1SUMMARY OF ', ic, ' CALLS TO HYBRD1'
            write (nwrite, '(A/)')      ' NPROB   N    NFEV  INFO  FINAL L2 NORM'
            do i = 1, ic
                write (nwrite, '(i4, i6, i7, i6, 1x, d15.7)') np(i), na(i), nf(i), nx(i), fnm(i)
            end do
            stop
        else
            nprob = nprobs(icase)
            n = ns(icase)
            lwa = (n*(3*n+13))/2
            ntries = ntriess(icase)

            if (allocated(fvec)) deallocate(fvec)
            if (allocated(wa)) deallocate(wa)
            if (allocated(x)) deallocate(x)
            allocate(fvec(n))
            allocate(wa(lwa))
            allocate(x(n))

            factor = one
            do k = 1, ntries
                ic = ic + 1
                call initpt(n, x, NPRob, factor)
                call vecfcn(n, x, fvec, NPRob)
                fnorm1 = enorm(n, fvec)
                write (nwrite, '(////5x,A,I5,5X,A,I5,5X//)') ' PROBLEM', NPRob, ' DIMENSION', n
                NFEv = 0
                call hybrd1(fcn, n, x, fvec, tol, info, wa, lwa)
                fnorm2 = enorm(n, fvec)
                np(ic) = NPRob
                na(ic) = n
                nf(ic) = NFEv
                nx(ic) = info
                fnm(ic) = fnorm2
                write (nwrite, '(5X,A,D15.7//5X,A,D15.7//5X,A,I10//5X,A,18X,I10//5X,A//,*(5X,5D15.7/))') &
                               ' INITIAL L2 NORM OF THE RESIDUALS', fnorm1, &
                               ' FINAL L2 NORM OF THE RESIDUALS  ', fnorm2, &
                               ' NUMBER OF FUNCTION EVALUATIONS  ', NFEv,   &
                               ' EXIT PARAMETER', info, &
                               ' FINAL APPROXIMATE SOLUTION', x(1:n)
                factor = ten*factor

                ! compare with previously generated solutions:
                if (any(abs(solution(ic) - x)>tol) .and. &
                    any(abs((solution(ic) - x)/(solution(ic))) > solution_reltol)) then
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
!  equation solver. fcn should only call the testing function
!  subroutine vecfcn with the appropriate value of problem
!  number (nprob).

    subroutine fcn(n, x, Fvec, Iflag)
        implicit none

        integer, intent(in) :: n !! the number of variables.
        real(wp), intent(in) :: x(n) !! independant variable vector
        real(wp), intent(out) :: fvec(n) !! value of function at `x`
        integer, intent(inout) :: iflag !! set to <0 to terminate execution

        call vecfcn(n, x, Fvec, NPRob)
        NFEv = NFEv + 1

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
        case (1); x = [0.1000000000000000E+01_wp,0.1000000000000000E+01_wp]
        case (2); x = [0.1000000000000000E+01_wp,0.1000000000000000E+01_wp]
        case (3); x = [0.1000000000000000E+01_wp,0.1000000000000054E+01_wp]
        case (4); x = [-0.1717232736993598E-17_wp,0.1717232736993598E-18_wp,0.4885791716489701E-17_wp,&
                        0.4885791716489701E-17_wp]
        case (5); x = [0.1137338840805565E-17_wp,-0.1137338840805565E-18_wp,0.1509231876962185E-17_wp,&
                        0.1509231876962185E-17_wp]
        case (6); x = [0.2071578632741476E-17_wp,-0.2071578632741476E-18_wp,0.3365608460018520E-17_wp,&
                        0.3365608460018520E-17_wp]
        case (7); x = [0.1098159327798559E-04_wp,0.9106146740037904E+01_wp]
        case (8); x = [0.1098159288127784E-04_wp,0.9106146743611500E+01_wp]
        case (9); x = [-0.9679740249513952E+00_wp,0.9471391408446626E+00_wp,-0.9695163103174519E+00_wp,&
                        0.9512476657649955E+00_wp]
        case (10); x = [-0.9679740249362032E+00_wp,0.9471391408151544E+00_wp,-0.9695163103329795E+00_wp,&
                        0.9512476657950136E+00_wp]
        case (11); x = [-0.9679740247487649E+00_wp,0.9471391404515958E+00_wp,-0.9695163105216826E+00_wp,&
                        0.9512476661606586E+00_wp]
        case (12); x = [0.1000000000000010E+01_wp,-0.1612103012704913E-13_wp,-0.8125674485239315E-34_wp]
        case (13); x = [0.1000000000004274E+01_wp,0.1388633807362216E-10_wp,0.0000000000000000E+00_wp]
        case (14); x = [0.9999999999999840E+00_wp,0.1238719384388127E-14_wp,0.0000000000000000E+00_wp]
        case (15); x = [-0.1572508640134011E-01_wp,0.1012434869369118E+01_wp,-0.2329916259567960E+00_wp,&
                        0.1260430087800365E+01_wp,-0.1513728922723441E+01_wp,0.9929964324318560E+00_wp]
        case (16); x = [-0.1572508639053690E-01_wp,0.1012434869369907E+01_wp,-0.2329916259623205E+00_wp,&
                        0.1260430087870332E+01_wp,-0.1513728922830754E+01_wp,0.9929964324984680E+00_wp]
        case (17); x = [-0.1530703652147214E-04_wp,0.9997897039319488E+00_wp,0.1476396369354227E-01_wp,&
                        0.1463423282995085E+00_wp,0.1000821103004075E+01_wp,-0.2617731140517609E+01_wp,&
                        0.4104403164477184E+01_wp,-0.3143612278555425E+01_wp,0.1052626408009917E+01_wp]
        case (18); x = [-0.1530703713913476E-04_wp,0.9997897039319459E+00_wp,0.1476396369267279E-01_wp,&
                        0.1463423283016577E+00_wp,0.1000821102994390E+01_wp,-0.2617731140495362E+01_wp,&
                        0.4104403164446781E+01_wp,-0.3143612278533631E+01_wp,0.1052626408003180E+01_wp]
        case (19); x = [0.8375125649983552E-01_wp,0.3127292952224503E+00_wp,0.5000000000008663E+00_wp,&
                        0.6872707047760241E+00_wp,0.9162487435008237E+00_wp]
        case (20); x = [0.6872707047770207E+00_wp,0.9162487435004409E+00_wp,0.4999999999996554E+00_wp,&
                        0.8375125649942240E-01_wp,0.3127292952234606E+00_wp]
        case (21); x = [0.5000000001403082E+00_wp,0.6872707047172968E+00_wp,0.8375125651152422E-01_wp,&
                        0.3127292951388482E+00_wp,0.9162487434920225E+00_wp]
        case (22); x = [0.6687659094768218E-01_wp,0.3666822990106460E+00_wp,0.2887406733471217E+00_wp,&
                        0.7112593271119373E+00_wp,0.6333177005294922E+00_wp,0.9331234090531204E+00_wp]
        case (23); x = [0.9331234090539151E+00_wp,0.2887406731191550E+00_wp,0.6687659094608109E-01_wp,&
                        0.7112593268807791E+00_wp,0.3666822992420550E+00_wp,0.6333177007580146E+00_wp]
        case (24); x = [0.3666822992460977E+00_wp,0.6333177007631426E+00_wp,0.7112593268766961E+00_wp,&
                        0.6687659094559376E-01_wp,0.9331234090527033E+00_wp,0.2887406731157666E+00_wp]
        case (25); x = [0.5806914977912385E-01_wp,0.2351716115667845E+00_wp,0.3380440957994889E+00_wp,&
                        0.4999999991484015E+00_wp,0.6619559063505691E+00_wp,0.7648283868041619E+00_wp,&
                        0.9419308505514702E+00_wp]
        case (26); x = [0.3380440947502161E+00_wp,0.4999999997939603E+00_wp,0.2351716123793344E+00_wp,&
                        0.7648283876237569E+00_wp,0.9419308503662596E+00_wp,0.5806914962097134E-01_wp,&
                        0.6619559054655015E+00_wp]
        case (27); x = [-0.3267366079625573E+02_wp,-0.2996209843926172E+02_wp,-0.8587775264169514E+02_wp,&
                        0.2222113097994968E+02_wp,0.5957249137089175E+02_wp,-0.1038025653217158E+01_wp,&
                        0.8600842862942351E+02_wp]
        case (28); x = [0.4985640222318974E-01_wp,0.1986351285003365E+00_wp,0.2698288337443381E+00_wp,&
                        0.4992723176748156E+00_wp,0.5007277255753518E+00_wp,0.7301712224978171E+00_wp,&
                        0.8013649159179719E+00_wp,0.9501436601762751E+00_wp]
        case (29); x = [0.4420534615691015E-01_wp,0.1994906721904692E+00_wp,0.2356191086780681E+00_wp,&
                        0.4160469078466623E+00_wp,0.4999999996232831E+00_wp,0.5839530926716184E+00_wp,&
                        0.7643808911925417E+00_wp,0.8005093278089829E+00_wp,0.9557946538314642E+00_wp]
        case (30); x = [0.1000000000000008E+01_wp,0.1000000000000008E+01_wp,0.1000000000000008E+01_wp,&
                        0.1000000000000008E+01_wp,0.1000000000000008E+01_wp,0.1000000000000008E+01_wp,&
                        0.1000000000000008E+01_wp,0.1000000000000008E+01_wp,0.1000000000000008E+01_wp,&
                        0.9999999999999193E+00_wp]
        case (31); x = [0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,&
                        0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,&
                        0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,0.1000000000000002E+01_wp,&
                        0.9999999999999805E+00_wp]
        case (32); x = [0.1000000000000126E+01_wp,0.1000000000000126E+01_wp,0.1000000000000126E+01_wp,&
                        0.1000000000000126E+01_wp,0.1000000000000126E+01_wp,0.1000000000000126E+01_wp,&
                        0.1000000000000126E+01_wp,0.1000000000000126E+01_wp,0.1000000000000126E+01_wp,&
                        0.9999999999987337E+00_wp]
        case (33); x = [0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,&
                        0.1000000000000116E+01_wp,0.1000000000000116E+01_wp,0.9999999999965015E+00_wp]
        case (34); x = [0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,0.9999999999999820E+00_wp,&
                        0.1000000000000726E+01_wp]
        case (35); x = [-0.4316498251876486E-01_wp,-0.8157715653538729E-01_wp,-0.1144857143805310E+00_wp,&
                        -0.1409735768625996E+00_wp,-0.1599086961819857E+00_wp,-0.1698772023127759E+00_wp,&
                        -0.1690899837812081E+00_wp,-0.1552495352218312E+00_wp,-0.1253558916789345E+00_wp,&
                        -0.7541653368589182E-01_wp]
        case (36); x = [-0.4316498251881876E-01_wp,-0.8157715653546944E-01_wp,-0.1144857143805966E+00_wp,&
                        -0.1409735768626190E+00_wp,-0.1599086961819499E+00_wp,-0.1698772023126901E+00_wp,&
                        -0.1690899837811062E+00_wp,-0.1552495352217907E+00_wp,-0.1253558916789970E+00_wp,&
                        -0.7541653368596339E-01_wp]
        case (37); x = [-0.4316498254524553E-01_wp,-0.8157715658128553E-01_wp,-0.1144857144024714E+00_wp,&
                        -0.1409735768723229E+00_wp,-0.1599086963003020E+00_wp,-0.1698772022538506E+00_wp,&
                        -0.1690899837944776E+00_wp,-0.1552495352060321E+00_wp,-0.1253558916432041E+00_wp,&
                        -0.7541653366609002E-01_wp]
        case (38); x = [-0.1528138835625800E+00_wp]
        case (39); x = [-0.1528138835625801E+00_wp]
        case (40); x = [-0.1528138835625801E+00_wp]
        case (41); x = [-0.4316498251876487E-01_wp,-0.8157715653538730E-01_wp,-0.1144857143805310E+00_wp,&
                        -0.1409735768625996E+00_wp,-0.1599086961819857E+00_wp,-0.1698772023127759E+00_wp,&
                        -0.1690899837812081E+00_wp,-0.1552495352218312E+00_wp,-0.1253558916789344E+00_wp,&
                        -0.7541653368589175E-01_wp]
        case (42); x = [-0.4316498251881876E-01_wp,-0.8157715653546946E-01_wp,-0.1144857143805966E+00_wp,&
                        -0.1409735768626190E+00_wp,-0.1599086961819498E+00_wp,-0.1698772023126901E+00_wp,&
                        -0.1690899837811062E+00_wp,-0.1552495352217907E+00_wp,-0.1253558916789970E+00_wp,&
                        -0.7541653368596334E-01_wp]
        case (43); x = [-0.4316498251876519E-01_wp,-0.8157715653538752E-01_wp,-0.1144857143805303E+00_wp,&
                        -0.1409735768625980E+00_wp,-0.1599086961819844E+00_wp,-0.1698772023127748E+00_wp,&
                        -0.1690899837812073E+00_wp,-0.1552495352218307E+00_wp,-0.1253558916789341E+00_wp,&
                        -0.7541653368589159E-01_wp]
        case (44); x = [0.5526115715943968E-01_wp,0.5695713693095779E-01_wp,0.5889020336090119E-01_wp,&
                        0.6113556183519356E-01_wp,0.6377799292665667E-01_wp,0.6700414432471043E-01_wp,&
                        0.2079417258113421E+00_wp,0.1642681193175131E+00_wp,0.8643704095817571E-01_wp,&
                        0.9133212907808361E-01_wp]
        case (45); x = [0.3439628896235289E-01_wp,0.3503231575416022E-01_wp,0.3571919583574593E-01_wp,&
                        0.3646522422001942E-01_wp,0.3728091174083566E-01_wp,0.3817986258974846E-01_wp,&
                        0.3918014109819012E-01_wp,0.4030650261419996E-01_wp,0.1797201916815169E+00_wp,&
                        0.1562408814749922E+00_wp]
        case (46); x = [0.1888395221036672E+02_wp,0.2516777354434988E+02_wp,0.1888527511739392E+02_wp,&
                        0.1888602114554983E+02_wp,0.1888683683452867E+02_wp,0.1888773578345333E+02_wp,&
                        0.1888873606083097E+02_wp,0.1888986242379029E+02_wp,0.1902927611240455E+02_wp,&
                        0.1900579680335179E+02_wp]
        case (47); x = [0.9999999999999992E+00_wp,0.9999999999999984E+00_wp,0.9999999999999977E+00_wp,&
                        0.9999999999999969E+00_wp,0.9999999999999961E+00_wp,0.9999999999999953E+00_wp,&
                        0.9999999999999946E+00_wp,0.9999999999999938E+00_wp,0.9999999999999930E+00_wp,&
                        0.9999999999999922E+00_wp]
        case (48); x = [0.1000000000000000E+01_wp,0.1000000000000000E+01_wp,0.1000000000000000E+01_wp,&
                        0.1000000000000000E+01_wp,0.1000000000000001E+01_wp,0.1000000000000001E+01_wp,&
                        0.1000000000000001E+01_wp,0.1000000000000001E+01_wp,0.1000000000000001E+01_wp,&
                        0.1000000000000001E+01_wp]
        case (49); x = [0.1000000000000000E+01_wp,0.1000000000000000E+01_wp,0.1000000000000000E+01_wp,&
                        0.1000000000000000E+01_wp,0.1000000000000000E+01_wp,0.1000000000000000E+01_wp,&
                        0.1000000000000000E+01_wp,0.1000000000000000E+01_wp,0.1000000000000000E+01_wp,&
                        0.1000000000000000E+01_wp]
        case (50); x = [-0.5707221307212121E+00_wp,-0.6818069509055232E+00_wp,-0.7022100775689857E+00_wp,&
                        -0.7055106309936168E+00_wp,-0.7049061557572888E+00_wp,-0.7014966060124587E+00_wp,&
                        -0.6918893211477919E+00_wp,-0.6657965141985400E+00_wp,-0.5960351099566767E+00_wp,&
                        -0.4164122574358191E+00_wp]
        case (51); x = [-0.5707221320932939E+00_wp,-0.6818069495100820E+00_wp,-0.7022100764111258E+00_wp,&
                        -0.7055106298493696E+00_wp,-0.7049061556844529E+00_wp,-0.7014966070294095E+00_wp,&
                        -0.6918893223739674E+00_wp,-0.6657965143753547E+00_wp,-0.5960351092038981E+00_wp,&
                        -0.4164122574142932E+00_wp]
        case (52); x = [-0.5707221320171143E+00_wp,-0.6818069499829604E+00_wp,-0.7022100760171542E+00_wp,&
                        -0.7055106298955310E+00_wp,-0.7049061557301967E+00_wp,-0.7014966070327222E+00_wp,&
                        -0.6918893223590803E+00_wp,-0.6657965144072677E+00_wp,-0.5960351090088830E+00_wp,&
                        -0.4164122575177334E+00_wp]
        case (53); x = [-0.4283028636053099E+00_wp,-0.4765964242962535E+00_wp,-0.5196524638125549E+00_wp,&
                        -0.5580993246169652E+00_wp,-0.5925061569509362E+00_wp,-0.6245036821428087E+00_wp,&
                        -0.6232394714478015E+00_wp,-0.6213938418388717E+00_wp,-0.6204535966122983E+00_wp,&
                        -0.5864692707477792E+00_wp]
        case (54); x = [-0.4283028634881692E+00_wp,-0.4765964236396713E+00_wp,-0.5196524642776766E+00_wp,&
                        -0.5580993248351936E+00_wp,-0.5925061568131795E+00_wp,-0.6245036817962692E+00_wp,&
                        -0.6232394720687789E+00_wp,-0.6213938417874499E+00_wp,-0.6204535965224117E+00_wp,&
                        -0.5864692707287930E+00_wp]
        case (55); x = [-0.4283028635608067E+00_wp,-0.4765964243232715E+00_wp,-0.5196524637037395E+00_wp,&
                        -0.5580993248328234E+00_wp,-0.5925061568292707E+00_wp,-0.6245036822076749E+00_wp,&
                        -0.6232394714256790E+00_wp,-0.6213938418143937E+00_wp,-0.6204535966527651E+00_wp,&
                        -0.5864692707189498E+00_wp]

        case default
            error stop 'invalid case'
        end select

    end function solution
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine defines fourteen test functions. the first
!  five test functions are of dimensions 2,4,2,4,3, respectively,
!  while the remaining test functions are of variable dimension
!  n for any n greater than or equal to 1 (problem 6 is an
!  exception to this, since it does not allow n = 1).

    subroutine vecfcn(n, x, Fvec, Nprob)
        implicit none

        integer, intent(in) :: n !! a positive integer input variable.
        integer, intent(in) :: nprob !! a positive integer input variable which defines the
                                     !! number of the problem. nprob must not exceed 14.
        real(wp), intent(in) :: x(n) !! an input array of length n.
        real(wp), intent(out) :: fvec(n) !! an output array of length n which contains the nprob
                                         !! function vector evaluated at x.

        real(wp), parameter :: zero = 0.0_wp
        real(wp), parameter :: one = 1.0_wp
        real(wp), parameter :: two = 2.0_wp
        real(wp), parameter :: three = 3.0_wp
        real(wp), parameter :: five = 5.0_wp
        real(wp), parameter :: eight = 8.0_wp
        real(wp), parameter :: ten = 10.0_wp
        real(wp), parameter :: c1 = 1.0e4_wp
        real(wp), parameter :: c2 = 1.0001_wp
        real(wp), parameter :: c3 = 2.0e2_wp
        real(wp), parameter :: c4 = 2.02e1_wp
        real(wp), parameter :: c5 = 1.98e1_wp
        real(wp), parameter :: c6 = 1.8e2_wp
        real(wp), parameter :: c7 = 2.5e-1_wp
        real(wp), parameter :: c8 = 5.0e-1_wp
        real(wp), parameter :: c9 = 2.9e1_wp

        integer :: i, iev, ivar, j, k, k1, k2, kp1, ml, mu
        real(wp) :: h, prod, sum, sum1, sum2, temp, temp1, &
                    temp2, ti, tj, tk, tpi

        fvec(1:n) = zero

        ! PROBLEM SELECTOR.

        select case (Nprob)
        case (2)
            ! POWELL SINGULAR FUNCTION.
            Fvec(1) = x(1) + ten*x(2)
            Fvec(2) = sqrt(five)*(x(3) - x(4))
            Fvec(3) = (x(2) - two*x(3))**2
            Fvec(4) = sqrt(ten)*(x(1) - x(4))**2
        case (3)
            ! POWELL BADLY SCALED FUNCTION.
            Fvec(1) = c1*x(1)*x(2) - one
            Fvec(2) = exp(-x(1)) + exp(-x(2)) - c2
        case (4)
            ! WOOD FUNCTION.
            temp1 = x(2) - x(1)**2
            temp2 = x(4) - x(3)**2
            Fvec(1) = -c3*x(1)*temp1 - (one - x(1))
            Fvec(2) = c3*temp1 + c4*(x(2) - one) + c5*(x(4) - one)
            Fvec(3) = -c6*x(3)*temp2 - (one - x(3))
            Fvec(4) = c6*temp2 + c4*(x(4) - one) + c5*(x(2) - one)
        case (5)
            ! HELICAL VALLEY FUNCTION.
            tpi = eight*atan(one)
            temp1 = sign(c7, x(2))
            if (x(1) > zero) temp1 = atan(x(2)/x(1))/tpi
            if (x(1) < zero) temp1 = atan(x(2)/x(1))/tpi + c8
            temp2 = sqrt(x(1)**2 + x(2)**2)
            Fvec(1) = ten*(x(3) - ten*temp1)
            Fvec(2) = ten*(temp2 - one)
            Fvec(3) = x(3)
        case (6)
            ! WATSON FUNCTION.
            do k = 1, n
                Fvec(k) = zero
            end do
            do i = 1, 29
                ti = dfloat(i)/c9
                sum1 = zero
                temp = one
                do j = 2, n
                    sum1 = sum1 + dfloat(j - 1)*temp*x(j)
                    temp = ti*temp
                end do
                sum2 = zero
                temp = one
                do j = 1, n
                    sum2 = sum2 + temp*x(j)
                    temp = ti*temp
                end do
                temp1 = sum1 - sum2**2 - one
                temp2 = two*ti*sum2
                temp = one/ti
                do k = 1, n
                    Fvec(k) = Fvec(k) + temp*(dfloat(k - 1) - temp2)*temp1
                    temp = ti*temp
                end do
            end do
            temp = x(2) - x(1)**2 - one
            Fvec(1) = Fvec(1) + x(1)*(one - two*temp)
            Fvec(2) = Fvec(2) + temp
        case (7)
            ! CHEBYQUAD FUNCTION.
            do k = 1, n
                Fvec(k) = zero
            end do
            do j = 1, n
                temp1 = one
                temp2 = two*x(j) - one
                temp = two*temp2
                do i = 1, n
                    Fvec(i) = Fvec(i) + temp2
                    ti = temp*temp2 - temp1
                    temp1 = temp2
                    temp2 = ti
                end do
            end do
            tk = one/dfloat(n)
            iev = -1
            do k = 1, n
                Fvec(k) = tk*Fvec(k)
                if (iev > 0) Fvec(k) = Fvec(k) + one/(dfloat(k)**2 - one)
                iev = -iev
            end do
        case (8)
            ! BROWN ALMOST-LINEAR FUNCTION.
            sum = -dfloat(n + 1)
            prod = one
            do j = 1, n
                sum = sum + x(j)
                prod = x(j)*prod
            end do
            do k = 1, n
                Fvec(k) = x(k) + sum
            end do
            Fvec(n) = prod - one
        case (9)
            ! DISCRETE BOUNDARY VALUE FUNCTION.
            h = one/dfloat(n + 1)
            do k = 1, n
                temp = (x(k) + dfloat(k)*h + one)**3
                temp1 = zero
                if (k /= 1) temp1 = x(k - 1)
                temp2 = zero
                if (k /= n) temp2 = x(k + 1)
                Fvec(k) = two*x(k) - temp1 - temp2 + temp*h**2/two
            end do
        case (10)
            ! DISCRETE INTEGRAL EQUATION FUNCTION.
            h = one/dfloat(n + 1)
            do k = 1, n
                tk = dfloat(k)*h
                sum1 = zero
                do j = 1, k
                    tj = dfloat(j)*h
                    temp = (x(j) + tj + one)**3
                    sum1 = sum1 + tj*temp
                end do
                sum2 = zero
                kp1 = k + 1
                if (n >= kp1) then
                    do j = kp1, n
                        tj = dfloat(j)*h
                        temp = (x(j) + tj + one)**3
                        sum2 = sum2 + (one - tj)*temp
                    end do
                end if
                Fvec(k) = x(k) + h*((one - tk)*sum1 + tk*sum2)/two
            end do
        case (11)
            ! TRIGONOMETRIC FUNCTION.
            sum = zero
            do j = 1, n
                Fvec(j) = cos(x(j))
                sum = sum + Fvec(j)
            end do
            do k = 1, n
                Fvec(k) = dfloat(n + k) - sin(x(k)) - sum - dfloat(k)*Fvec(k)
            end do
        case (12)
            ! VARIABLY DIMENSIONED FUNCTION.
            sum = zero
            do j = 1, n
                sum = sum + dfloat(j)*(x(j) - one)
            end do
            temp = sum*(one + two*sum**2)
            do k = 1, n
                Fvec(k) = x(k) - one + dfloat(k)*temp
            end do
        case (13)
            ! BROYDEN TRIDIAGONAL FUNCTION.
            do k = 1, n
                temp = (three - two*x(k))*x(k)
                temp1 = zero
                if (k /= 1) temp1 = x(k - 1)
                temp2 = zero
                if (k /= n) temp2 = x(k + 1)
                Fvec(k) = temp - temp1 - two*temp2 + one
            end do
        case (14)
            ! BROYDEN BANDED FUNCTION.
            ml = 5
            mu = 1
            do k = 1, n
                k1 = max(1, k - ml)
                k2 = min(k + mu, n)
                temp = zero
                do j = k1, k2
                    if (j /= k) temp = temp + x(j)*(one + x(j))
                end do
                Fvec(k) = x(k)*(two + five*x(k)**2) + one - temp
            end do
        case default
            ! ROSENBROCK FUNCTION.
            Fvec(1) = one - x(1)
            Fvec(2) = ten*(x(2) - x(1)**2)
        end select

    end subroutine vecfcn
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine specifies the standard starting points for
!  the functions defined by subroutine vecfcn. the subroutine
!  returns in x a multiple (factor) of the standard starting
!  point. for the sixth function the standard starting point is
!  zero, so in this case, if factor is not unity, then the
!  subroutine returns the vector  x(j) = factor, j=1,...,n.

    subroutine initpt(n, x, Nprob, Factor)

        implicit none

        integer, intent(in) :: n !! a positive integer input variable.
        real(wp), intent(out) :: x(n) !! an output array of length n which contains the standard
                                      !! starting point for problem nprob multiplied by factor.
        integer, intent(in) :: Nprob !! a positive integer input variable which defines the
                                     !! number of the problem. nprob must not exceed 14.
        real(wp), intent(in) :: Factor !! an input variable which specifies the multiple of
                                       !! the standard starting point. if factor is unity, no
                                       !! multiplication is performed.

        integer :: ivar, j
        real(wp) :: h, tj

        real(wp), parameter :: zero = 0.0_wp
        real(wp), parameter :: half = 0.5_wp
        real(wp), parameter :: one = 1.0_wp
        real(wp), parameter :: three = 3.0_wp
        real(wp), parameter :: c1 = 1.2_wp

        x(1:n) = zero

        ! selection of initial point.

        select case (nprob)
        case (2)
            ! powell singular function.
            x(1) = three
            x(2) = -one
            x(3) = zero
            x(4) = one
        case (3)
            ! powell badly scaled function.
            x(1) = zero
            x(2) = one
        case (4)
            ! wood function.
            x(1) = -three
            x(2) = -one
            x(3) = -three
            x(4) = -one
        case (5)
            ! helical valley function.
            x(1) = -one
            x(2) = zero
            x(3) = zero
        case (6)
            ! watson function.
            do j = 1, n
                x(j) = zero
            end do
        case (7)
            ! chebyquad function.
            h = one/dfloat(n + 1)
            do j = 1, n
                x(j) = dfloat(j)*h
            end do
        case (8)
            ! brown almost-linear function.
            do j = 1, n
                x(j) = half
            end do
        case (9, 10)
            ! discrete boundary value and integral equation functions.
            h = one/dfloat(n + 1)
            do j = 1, n
                tj = dfloat(j)*h
                x(j) = tj*(tj - one)
            end do
        case (11)
            ! trigonometric function.
            h = one/dfloat(n)
            do j = 1, n
                x(j) = h
            end do
        case (12)
            ! variably dimensioned function.
            h = one/dfloat(n)
            do j = 1, n
                x(j) = one - dfloat(j)*h
            end do
        case (13, 14)
            ! broyden tridiagonal and banded functions.
            do j = 1, n
                x(j) = -one
            end do
        case default
            ! rosenbrock function.
            x(1) = -c1
            x(2) = one
        end select

        ! compute multiple of initial point.

        if (factor /= one) then
            if (nprob == 6) then
                do j = 1, n
                    x(j) = factor
                end do
            else
                do j = 1, n
                    x(j) = factor*x(j)
                end do
            end if
        end if

    end subroutine initpt
!*****************************************************************************************

!*****************************************************************************************
end program test
!*****************************************************************************************