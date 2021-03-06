program main

use utilities
use uncertainty
!use asa183
use tIntegration

implicit none

integer :: passos = 100, adjustGrid = 700
real*8 :: t0 = 0.d0, tf = 406.d0, y0 = 13.18d0, R = 440.0d3, C = 5.6d-4
!real*8 :: t0 = 0.d0, tf = 10.d0, y0 = 1.0d0, R = 1.0d0, C = 3.d0
real*8 :: adjustErr = 0.2d0, Cminus, Cplus, Rminus, Rplus
real*8 :: nominalImpR = 1.d-2, measurementImpR = 3.d3
real*8 :: nominalImpC = 0.22d0
real*8 :: bleftR, tcenterR, brightR
real*8 :: bleftC, tcenterC, brightC
character(len=50) :: nameData='data.dat', nameR= 'R', nameC='C'
real*8, allocatable :: dataSet(:,:), sol(:)
integer :: i
real*8 :: ran
real*4 :: rteste
real*8, allocatable :: rteste1(:), rteste2(:)
real*8 :: xteste(2), vecWCS(2,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALIBRACAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call readDataFile(nameData, dataSet, 2, 111)
call fittingModel(dataSet, R, C, adjustErr, adjustGrid, 2, imethod=2, itol=5.d-3)
call solverRK4(sol, t0, tf, y0, R, C, passos, 0); !stop
!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WORST CASE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Rminus = R*nominalImpR
Rplus = R*nominalImpR
C = 5.6d-4
Cminus = C*nominalImpC
Cplus = C*nominalImpC
vecWCS(1,1) = R; vecWCS(1,2) = Rminus; vecWCS(1,3) = Rplus
vecWCS(2,1) = C; vecWCS(2,2) = Cminus; vecWCS(2,3) = Cplus
!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FUZZY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bleftR = R-Rminus; tcenterR = R; brightR = R+Rplus
bleftC = C-Cminus; tcenterC = C; brightC = C+Cplus
call writeTriangular(bleftR,tcenterR,brightR,inpts=201,iname=nameR); !stop
call writeTriangular(bleftC,tcenterC,brightC,inpts=201,iname=nameC); !stop
call fuzzyTwoParameters(vecWCS,method=2,parGrid=10)
!stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RANDOM GENERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call randomUniformTwoParameters(vecWCS,method=2,parGrid=4999)
vecWCS(:,3) = vecWCS(:,3)/2.d0
call randomNormalTwoParameters(vecWCS,method=2,parGrid=4999)
stop

end program
