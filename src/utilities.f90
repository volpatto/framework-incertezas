module utilities

    implicit none

    contains

    !> @param fileName      Nome do arquivo de leitura dos dados
    !! @param dataSet       Matriz onde sera armazenado os dados lidos
    !! @param cols          Quantas colunas de dados estao no arquivo
    !! @param id            Id de abertura de arquivo para o read do Fortran
    subroutine readDataFile(fileName, dataSet, cols, id)

        implicit none

        character(len=50), intent(in) :: fileName
        integer, intent(in) :: cols, id
        real*8, allocatable :: dataSet(:,:)
        character(len=200) :: lineContent
        integer :: flag, counter=0, i, j

        open(unit=id, file=fileName)

        do
            read(id, "(A)", iostat=flag) lineContent
            if (flag .ne. 0) exit
            counter =  counter + 1
        enddo

        allocate(dataSet(counter,cols)); dataSet=0.d0

        rewind(id)
        do i=1,counter
            read(id,*) (dataSet(i,j), j=1,cols)
            !write(*,*) (dataSet(i,j), j=1,cols)
        enddo

        close(id)

    endsubroutine

    !> @param vet       Vetor de entrada
    !! @param val       Valor que sera buscado no vetor
    !! @param idx       Indice do valor buscado
    subroutine findIndex(vet, val, idx)

        implicit none

        real*8, allocatable :: vet(:)
        real*8, intent(in) :: val
        integer, intent(inout) :: idx
        integer :: imax, i

        imax = size(vet)

        do i=1,imax
            if (vet(i).eq.val) then
                idx = i; exit
            endif
            if (i.eq.imax) stop "Index not found"
        enddo

    endsubroutine

    !> Usa a rotina geradora de número aleatório uniformes r8_uni da asa183
    !! para criar um vetor com valores aleatorios entre minRange e maxRange.
    !! @param minRange      Valor minimo que podera ser gerado
    !! @param maxRange      Valor maximo que podera ser gerado
    !! @param rArray        Vetor com entradas pseudoaleatorias
    !! @param idimArray     (Opcional) Tamanho do vetor pseudoaleatorio
    subroutine rngArray(minRange, maxRange, rArray, idimArray)

        use asa183, only: r8_uni

        implicit none

        real*8, intent(in) :: minRange, maxRange
        real*8, allocatable, intent(inout) :: rArray(:)
        integer, intent(in), optional :: idimArray
        integer :: dimArray = 10, i
        integer :: rseed1, rseed2, n
        !integer, parameter :: seed = 82547
        integer, allocatable :: seed(:)

        if (present(idimArray)) dimArray = idimArray

        allocate(rArray(dimArray)); rArray = 0.d0

        call random_seed(size = n)
        allocate(seed(n))
        call random_seed(put=seed)
        rseed1 = irand(); rseed2 = irand()

        do i=1,dimArray
            rArray(i) = r8_uni(rseed1, rseed2)*(maxRange-minRange) + minRange
            print*, rArray(i)
        enddo

    endsubroutine

    function linspace(a,b,inpts)

        implicit none

        real*8, intent(in) :: a, b
        integer, intent(in), optional :: inpts
        real*8, allocatable :: linspace(:)
        integer :: npts = 100, i
        real*8 :: h

        if (present(inpts)) npts = inpts

        h = (b-a)/(npts-1)

        allocate(linspace(npts)); linspace = 0.d0

        do i=1,npts
            linspace(i) = (i-1)*h + a
            !print*, linspace(i)
        enddo

        return

    endfunction

    subroutine writeMatrixData(data_matrix, iname, id)

        implicit none

        real*8, intent(in) :: data_matrix(:,:)
        character(len=*), intent(in) :: iname
        integer, intent(in) :: id
        integer :: i, j
        character(len=50) :: filename, fname, char_id

        ! Opening solution output file
        ! Making an id char variable
        write(char_id, '(i0)') id
        write(filename, '(a)') trim(iname)
        ! Concatenate strings filename, char_id and ".dat"
        fname = trim(filename)//trim(char_id)//".dat"
        ! Opening file with unit as id and file name as attributed in fname
        open(unit=id,file=fname)

        do i = 1,size(data_matrix(:,1))
            write(id,*) (data_matrix(i,j), j=1,size(data_matrix(1,:)))
        enddo

        close(id)
    endsubroutine

endmodule
