module uncertainty

implicit none

contains

    !> @param dataSet       Matriz do conjunto de dados experimentais
    !! @param parFixed      Parametro de valor conhecido
    !! @param fitPar        Valor inicial do parametro a ser ajustado
    !! @param percErr       Percentual do fitPar que sera gerado o grid de busca do valor ajustado
    !! @param gridNumber    Numero de intervalos gerados para avaliar a funcao custo
    !! @param col           Coluna dos resultados experimentais em dataSet que sera utilizado no ajuste
    !! @param imethod       (Opcional) Metodo de integracao temporal: (1) Euler, (2) RK4
    !! @param itol          (Opcional) Tolerancia do criterio de parada do ajuste
    !! @param ikmax         (Opcional) Numero maximo de iteracoes do laco de ajuste
    subroutine fittingModel(dataSet, parFixed, fitPar, percErr, gridNumber, col, imethod, itol, ikmax)

        use utilities, only: findIndex
        use tIntegration

        implicit none

        real*8, intent(in) :: percErr, parFixed
        real*8, intent(inout) :: fitPar
        real*8, intent(in) :: dataset(:,:)
        integer, intent(in) :: gridNumber, col
        real*8, intent(in), optional :: itol
        integer, intent(in), optional :: ikmax, imethod
        real*8 :: gridSteps, varPar, t0, tf, y0
        real*8, allocatable :: gridVar(:), sol(:), error(:), errf(:)
        integer :: i, j, k, id, flag, counter=0, gridNodes, iOpt
        character(len=50) :: fileName, fileLine
        real*8 :: sol0, solf, e0, ef, tol=0.01
        integer :: nsteps = 200, iderr = 55555, kmax = 500, method = 1

        if (present(imethod)) method = imethod
        if (present(itol)) tol = itol
        if (present(ikmax)) kmax = ikmax

        !open(unit=iderr, file='errors.dat')
        gridNodes = gridNumber + 1
        allocate(gridVar(gridNodes));
        allocate(errf(gridNodes));

        write(*,'(/2x,(a)/)') "********* Begin fitting *********"

        do k=1,kmax
            write(*,'(2x,(a),1x,(i0))') "Iteration", k
            varPar = percErr*fitPar
            gridSteps = (2.d0*varPar)/gridNumber
            gridVar = 0.d0
            errf = 0.d0

            ! Creating range for looking parameter's value minimizer of error function
            do i=1,gridNodes
                gridVar(i) = fitPar-varPar + (i-1)*gridSteps
                y0 = dataSet(1,1)
                allocate(error(size(dataSet(:,col)))); error=0.d0
                do j=1,size(dataSet(:,col))-1
                    t0 = dataSet(j,col)
                    tf = dataSet(j+1,col)
                    if (method .eq. 1) call solverEuler(sol, t0, tf, y0, parFixed, gridVar(i), nsteps)
                    if (method .eq. 2) call solverRK4(sol, t0, tf, y0, parFixed, gridVar(i), nsteps)
                    sol0 = sol(1); solf = sol(nsteps+1)
                    if (j.eq.1) then
                        error(j) = dabs(dataSet(j,1)-sol0)
                        error(j+1) = dabs(dataSet(j+1,1)-solf)
                    else
                        error(j+1) = dabs(dataSet(j+1,1)-solf)
                    endif
                    y0 = solf
                    deallocate(sol)
                enddo
                !errf(i) = sum(error)
                errf(i) = sum(error)/dabs(sum(dataSet(:,1)))
                !write(iderr,*) gridVar(i), errf(i)
                deallocate(error)
            enddo

            call findIndex(errf,minval(errf),iOpt)
            fitPar = gridVar(iOpt)

            if (minval(errf).le.tol) exit

        enddo

        if (k.gt.kmax) write(*,'(2x,(a))') "Fitting iteration number exceeds."
        write(*,'(/2(2x,(a),1x,es15.5))') "fitPar =", fitPar, "relative error =", minval(errf) 
        if (method .eq. 1) write(*,'(/2x,(a))') "Numerical integration method: Implicit Euler"
        if (method .eq. 2) write(*,'(/2x,(a))') "Numerical integration method: Runge-Kutta 4"
        write(*,'(/2x,(a)/)') "********* End fitting *********"
        deallocate(gridVar); deallocate(errf)

        !close(iderr)

    endsubroutine

    !> Incompleto.
    !! @param parGrid       (Optional) [in] Número de intervalos para criação do grid de varredura
    !! @param imethod       (Optional) [in] Método de integração numérica
    subroutine worstCaseTwoParameters(parVal,parGrid,method)

        use utilities, only: linspace
        use tIntegration
        implicit none

        real*8, intent(in) :: parVal(2,3)
        integer, intent(in), optional :: parGrid, method
        real*8, allocatable :: gridWCS(:,:)
        integer :: parGridDefault = 6, gridNodes, methodDefault = 1
        integer :: i, id, passos = 100, naturalNumbering, j
        character(len=50) :: filename, fname, char_id
        real*8, allocatable :: sol(:)
        real*8 :: t0 = 0.d0, tf = 406.d0, y0 = 13.18d0

        if (present(parGrid)) parGridDefault = parGrid
        if (present(method)) methodDefault = method

        gridNodes = parGridDefault + 1

        allocate(gridWCS(2,gridNodes))
        gridWCS(1,:) = linspace(parVal(1,1)-parVal(1,2), parVal(1,1)+parVal(1,3),gridNodes)
        gridWCS(2,:) = linspace(parVal(2,1)-parVal(2,2), parVal(2,1)+parVal(2,3),gridNodes)

        do i=1,gridNodes
            do j=1,gridNodes
                naturalNumbering = (i-1)*gridNodes + j
                call solverRK4(sol, t0, tf, y0, gridWCS(1,i), gridWCS(2,j), passos,naturalNumbering,irefresh=1)
            enddo
        enddo

    endsubroutine

    !> Incompleto.
    !! @param parGrid       (Optional) [in] Número de intervalos para criação do grid de varredura
    !! @param imethod       (Optional) [in] Método de integração numérica
    subroutine fuzzyTwoParameters(parVal,parGrid,method)

        use utilities, only: linspace
        use tIntegration
        implicit none

        real*8, intent(in) :: parVal(2,3)
        integer, intent(in), optional :: parGrid, method
        real*8, allocatable :: gridWCS(:,:)
        integer :: parGridDefault = 6, gridNodes, methodDefault = 1
        integer :: i, id, passos = 100, naturalNumbering, j
        character(len=50) :: filename, fname, char_id
        real*8, allocatable :: a_cuts(:), R_cuts(:,:), C_cuts(:,:), trange(:)
        real*8, allocatable :: lsol(:), rsol(:), loutdata(:,:), routdata(:,:)
        real*8 :: t0 = 0.d0, tf = 406.d0, y0 = 13.18d0

        if (present(parGrid)) parGridDefault = parGrid
        if (present(method)) methodDefault = method

        ! Opening solution output file
        ! Making an id char variable
        id = 896762
        write(filename, '(a)') "fuzzy"
        ! Concatenate strings filename, char_id and ".dat"
        fname = trim(filename)//".dat"
        ! Opening file with unit as id and file name as attributed in fname
        open(unit=id,file=fname)

        gridNodes = parGridDefault + 1
        allocate(a_cuts(gridNodes))
        allocate(R_cuts(gridNodes,2))
        allocate(C_cuts(gridNodes,2))
        allocate(loutdata(gridNodes,passos+1)); loutdata = 0.d0
        allocate(routdata(gridNodes,passos+1)); routdata = 0.d0
        a_cuts = linspace(0.0d0,1.0d0,gridNodes)
        do i=1,gridNodes
            R_cuts(i,:) = findAlphaTri(a_cuts(i),parVal(1,1)-parVal(1,2),parVal(1,1),parVal(1,1)+parVal(1,3))
            C_cuts(i,:) = findAlphaTri(a_cuts(i),parVal(2,1)-parVal(2,2),parVal(2,1),parVal(2,1)+parVal(2,3))
            call solverRK4(lsol, t0, tf, y0, R_cuts(i,1), C_cuts(i,1), passos)
            loutdata(i,:) = lsol
            deallocate(lsol)
            call solverRK4(rsol, t0, tf, y0, R_cuts(i,2), C_cuts(i,2), passos)
            routdata(i,:) = rsol
            deallocate(rsol)
        enddo

        trange = linspace(t0,tf,passos+1)
        do i=1,gridNodes
            ! Opening solution output file
            ! Making an id char variable
            write(char_id, '(i0)') i
            write(filename, '(a)') "fuzzy"
            ! Concatenate strings filename, char_id and ".dat"
            fname = trim(filename)//trim(char_id)//".dat"
            ! Opening file with unit as id and file name as attributed in fname
            open(unit=i,file=fname)
            do j=1,(passos+1)
                write(i,*) trange(j), a_cuts(i), R_cuts(i,1), C_cuts(i,1), loutdata(i,j), R_cuts(i,2), C_cuts(i,2), routdata(i,j)
            enddo
            !write(id,*)
            close(i)
        enddo

        deallocate(a_cuts); deallocate(R_cuts); deallocate(C_cuts)
        close(id)

    endsubroutine

    !> Gera dos parâmetros de forma aleatória em distribuição uniforme resultados de forma a
    !! quantificar as incertezas. 
    !! @param parGrid       (Optional) [in] Número de intervalos para criação do grid de varredura
    !! @param imethod       (Optional) [in] Método de integração numérica
    subroutine randomUniformTwoParameters(parVal,parGrid,method)

        use utilities, only: linspace
        use random, only: rngArray
        use tIntegration
        implicit none

        real*8, intent(in) :: parVal(2,3)
        integer, intent(in), optional :: parGrid, method
        real*8, allocatable :: gridUniform1(:), gridUniform2(:)
        integer :: parGridDefault = 6, gridNodes, methodDefault = 1
        integer :: i, id, passos = 100
        character(len=50) :: filename, fname, char_id, nameR='R', namec='C'
        real*8, allocatable :: sol(:)
        real*8 :: t0 = 0.d0, tf = 406.d0, y0 = 13.18d0

        if (present(parGrid)) parGridDefault = parGrid
        if (present(method)) methodDefault = method

        gridNodes = parGridDefault + 1

        gridUniform1 = rngArray(parVal(1,1)-parVal(1,2),parVal(1,1)+parVal(1,3),idimArray=gridNodes,iname=nameR)
        gridUniform2 = rngArray(parVal(2,1)-parVal(2,2),parVal(2,1)+parVal(2,3),idimArray=gridNodes,iname=nameC)

        do i=1,gridNodes
            call solverRK4(sol, t0, tf, y0, gridUniform1(i), gridUniform2(i), passos,i,irefresh=1)
        enddo

    endsubroutine

    !> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! @param parGrid       (Optional) [in] Número de intervalos para criação do grid de varredura
    !! @param imethod       (Optional) [in] Método de integração numérica
    subroutine randomNormalTwoParameters(parVal,parGrid,method)

        use utilities, only: linspace
        use random, only: rngNormalArray
        use tIntegration
        implicit none

        real*8, intent(in) :: parVal(2,3)
        integer, intent(in), optional :: parGrid, method
        real*8, allocatable :: gridNormal1(:), gridNormal2(:)
        integer :: parGridDefault = 6, gridNodes, methodDefault = 1
        integer :: i, id, passos = 100
        character(len=50) :: filename, fname, char_id, nameR='R', nameC='C'
        real*8, allocatable :: sol(:)
        real*8 :: t0 = 0.d0, tf = 406.d0, y0 = 13.18d0

        if (present(parGrid)) parGridDefault = parGrid
        if (present(method)) methodDefault = method

        gridNodes = parGridDefault + 1

        gridNormal1 = rngNormalArray(parVal(1,1),parVal(1,3),idimArray=gridNodes,iname=nameR)
        gridNormal2 = rngNormalArray(parVal(2,1),parVal(2,3),idimArray=gridNodes,iname=nameC)

        do i=1,gridNodes
            call solverRK4(sol, t0, tf, y0, gridNormal1(i), gridNormal2(i), passos,i,irefresh=1)
        enddo

    endsubroutine

    !> Essa sub-rotina executa uma analise (e ajuste!) de worst case scenario para a incerteza associada a
    !! dois parametros: parVal (mensurado) e fitPar (estimado). parVal é conhecido, advindo de experimento,
    !! e tem uma margem de erro conhecida limitada por parErrorMinus (inf) e parErrorPlus (sup). fitPar é o
    !! parâmetro que será estimado para cada valor de parVal, dentro da sua margem de erro. Essa sub-rotina
    !! retorna os valores de fitErrorMinus (inf) e fitErrorPlus (sup) que limitam o erro para o fitPar.
    !! @param dataSet       [in] Matriz com os dados experimentais
    !! @param col           [in] Coluna que será utilizada para a leitura das medidas experimentais
    !! @param parVal        [in] O valor mensurado de parVal
    !! @param parErrorMinus [in] O erro inferior associado a parVal
    !! @param parErrorPlus  [in] O erro superior associado a parVal
    !! @param fitPar        [in/out] Estimativa inicial do valor de fitPar
    !! @param fitErrorMinus [out] Erro inferior calculado do fitPar
    !! @param fitErrorPlus  [out] Erro superior calculado do fitPar
    !! @param fitErr        (Optional) [in] Percentual do valor de fitPar para criação de grid de varredura
    !! @param parGrid       (Optional) [in] Número de intervalos para criação do grid de varredura
    !! @param imethod       (Optional) [in] Método de integração numérica
    subroutine worstCaseTwoParametersFitting(dataSet,col,parVal,parErrorMinus,parErrorPlus, &
            fitPar,fitErrorMinus,fitErrorPlus,fitErr,parGrid,method)

        implicit none

        real*8, intent(in) :: parVal, parErrorMinus, parErrorPlus, dataSet(:,:)
        real*8, intent(inout) :: fitPar
        real*8, intent(out) :: fitErrorMinus, fitErrorPlus
        integer, intent(in) :: col
        real*8, intent(in), optional :: fitErr
        integer, intent(in), optional :: parGrid, method
        real*8 :: fitErrDefault = 2.d-1, gridSteps
        real*8, allocatable :: parValues(:,:)
        integer :: parGridDefault = 500, gridNodes, methodDefault = 1
        integer :: i, id
        character(len=50) :: filename, fname, char_id

        if (present(fitErr)) fitErrDefault = fitErr
        if (present(parGrid)) parGridDefault = parGrid
        if (present(method)) methodDefault = method

        ! Opening solution output file
        ! Making an id char variable
        id = 11111
        write(filename, '(a)') "RC"
        ! Concatenate strings filename, char_id and ".dat"
        fname = trim(filename)//".dat"
        ! Opening file with unit as id and file name as attributed in fname
        open(unit=id,file=fname)

        gridSteps = (parErrorPlus+parErrorMinus)/parGridDefault
        gridNodes = parGridDefault + 1

        allocate(parValues(gridNodes,2)); parValues = 0.d0
        parValues(:,2) = fitPar

        do i=1,gridNodes
            parValues(i,1) = parVal - parErrorMinus + (i-1)*gridSteps
            call fittingModel(dataSet, parValues(i,1), parValues(i,2), fitErrDefault, parGridDefault, col, imethod=methodDefault)
            write(id,*) parValues(i,1), parValues(i,2)
        enddo

        fitErrorPlus = maxval(parValues(:,2)); fitErrorMinus = minval(parValues(:,2))

        close(id)

    endsubroutine

    !> Função que retorna o valor da função de pertinencia triangular de um dado x.
    !! @param x         [in] Valor de entrada
    !! @param bleft     [in] Valor da coordenada que define o limite esquerdo onde a função é 0
    !! @param tleft     [in] Valor da coordenada que define o limite esquerdo onde a função é 1
    !! @param tright    [in] Valor da coordenada que define o limite direito onde a função é 1
    !! @param bright    [in] Valor da coordenada que define o limite direito onde a função é 0
    function triangularMembership(x,bleft,tcenter,bright)

        implicit none

        real*8, intent(in) :: x, bleft, tcenter, bright
        real*8 :: triangularMembership
        real*8 :: a, b

        if ((x.ge.bleft).and.(x.lt.tcenter)) then
            a = 1.d0/(tcenter-bleft); b = -a*bleft
            triangularMembership = a*x + b
        else if ((x.gt.tcenter).and.(x.le.bright)) then
            a = 1.d0/(tcenter-bright); b = -a*bright
            triangularMembership = a*x + b
        else
            triangularMembership = 1.d0
        endif

        return

    endfunction

    !> Função que retorna o valor da função de pertinencia trapezoidal de um dado x.
    !! @param x         [in] Valor de entrada
    !! @param bleft     [in] Valor da coordenada que define o limite esquerdo onde a função é 0
    !! @param tleft     [in] Valor da coordenada que define o limite esquerdo onde a função é 1
    !! @param tright    [in] Valor da coordenada que define o limite direito onde a função é 1
    !! @param bright    [in] Valor da coordenada que define o limite direito onde a função é 0
    function trapezoidalMembership(x,bleft,tleft,tright,bright)

        implicit none

        real*8, intent(in) :: x, bleft, tleft, tright, bright
        real*8 :: trapezoidalMembership
        real*8 :: a, b

        if ((x.ge.bleft).and.(x.lt.tleft)) then
            a = 1.d0/(tleft-bleft); b = -a*bleft
            trapezoidalMembership = a*x + b
        else if ((x.gt.tright).and.(x.le.bright)) then
            a = 1.d0/(tright-bright); b = -a*bright
            trapezoidalMembership = a*x + b
        else
            trapezoidalMembership = 1.d0
        endif

        return

    endfunction

    subroutine writeTriangular(bleft,tcenter,bright,inpts,iname)

        use utilities

        implicit none

        real*8, intent(in) :: bleft, tcenter, bright
        integer, intent(in), optional :: inpts
        character(len=*), intent(in), optional :: iname
        integer :: npts = 100, i, id
        real*8, allocatable :: xrange(:)
        character(len=50) :: filename, fname, char_id, idname

        if (present(inpts)) npts = inpts
        if (present(iname)) idname = trim(iname)

        ! Opening solution output file
        ! Making an id char variable
        id = 8463294
        write(filename, '(a)') "alphaTri"
        ! Concatenate strings filename, char_id and ".dat"
        if(present(iname)) then 
            fname = trim(filename)//trim(idname)//".dat"
        else
            fname = trim(filename)//".dat"
        endif
        ! Opening file with unit as id and file name as attributed in fname
        open(unit=id,file=fname)

        xrange = linspace(bright,bleft,npts)

        do i=1,size(xrange)
            write(id,'(2x,2(es15.5))') xrange(i), triangularMembership(xrange(i),bleft,tcenter,bright)
        enddo

        close(id)

    endsubroutine

    subroutine writeTrapezoidal(bleft,tleft,tright,bright,inpts,iname)

        use utilities

        implicit none

        real*8, intent(in) :: bleft, tleft, tright, bright
        integer, intent(in), optional :: inpts
        character(len=*), intent(in), optional :: iname
        integer :: npts = 100, i, id, ilen
        real*8, allocatable :: xrange(:)
        character(len=50) :: filename, fname, char_id, idname

        if (present(inpts)) npts = inpts
        if (present(iname)) idname = trim(iname)

        ! Opening solution output file
        ! Making an id char variable
        id = 736381
        write(filename, '(a)') "alphaTrap"
        ! Concatenate strings filename, char_id and ".dat"
        if(present(iname)) then 
            fname = trim(filename)//trim(idname)//".dat"
        else
            fname = trim(filename)//".dat"
        endif
        ! Opening file with unit as id and file name as attributed in fname
        open(unit=id,file=fname)

        xrange = linspace(bright,bleft,npts)

        do i=1,size(xrange)
            write(id,'(2x,2(es15.5))') xrange(i), trapezoidalMembership(xrange(i),bleft,tleft,tright,bright)
        enddo

        close(id)

    endsubroutine

    function findAlphaTri(aval,bleft,tcenter,bright)

        implicit none

        real*8, intent(in) :: aval, bleft, tcenter, bright
        real*8 :: findAlphaTri(2)
        real*8 :: a, b

        if (aval .ne. 1.d0) then        
            a = 1.d0/(tcenter-bleft); b = -a*bleft
            findAlphaTri(1) = (aval-b)/a 

            a = 1.d0/(tcenter-bright); b = -a*bright
            findAlphaTri(2) = (aval-b)/a 
        else
            findAlphaTri(1) = tcenter 
            findAlphaTri(2) = tcenter
        endif

        return

    endfunction

end module
