module random

    implicit none

    contains 

    !> Usa a rotina geradora de número aleatório uniformes r8_uni da asa183
    !! para criar um vetor com valores aleatorios entre minRange e maxRange.
    !! @param minRange      Valor minimo que podera ser gerado
    !! @param maxRange      Valor maximo que podera ser gerado
    !! @param rArray        Vetor com entradas pseudoaleatorias
    !! @param idimArray     (Opcional) Tamanho do vetor pseudoaleatorio
    !! @param iname         (Opcional) Nome do arquivo para escrita
    function rngArray(minRange, maxRange, idimArray, iname)

        implicit none

        real*8, intent(in) :: minRange, maxRange
        real*8, allocatable :: rngArray(:)
        integer, intent(in), optional :: idimArray
        character(len=*), intent(in), optional :: iname
        integer :: dimArray = 10, i, id
        integer :: rseed1, rseed2, n
        integer, allocatable :: seed(:)
        character(len=50) :: filename, fname, char_id, idname

        if (present(idimArray)) dimArray = idimArray
        if (present(iname)) then 
            idname = trim(iname)
            ! Opening solution output file
            ! Making an id char variable
            id = 7825374
            write(filename, '(a)') "uniDist"
            ! Concatenate strings filename, char_id and ".dat"
            fname = trim(filename)//trim(idname)//".dat"
            ! Opening file with unit as id and file name as attributed in fname
            open(unit=id,file=fname)
        endif

        !allocate(rArray(dimArray)); rArray = 0.d0
        allocate(rngArray(dimArray)); rngArray = 0.d0

        call random_seed(size = n)
        allocate(seed(n))
        call random_seed(put=seed)
        rseed1 = irand(); rseed2 = irand()

        do i=1,dimArray
            rngArray(i) = r8_uni(rseed1, rseed2)*(maxRange-minRange) + minRange
            if (present(iname)) write(id,*) rngArray(i)
        enddo

        if (present(iname)) close(id)

    endfunction

    !>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! @param mean          Média
    !! @param sd            Desvio-padrão
    !! @param idimArray     (Opcional) Tamanho do vetor pseudoaleatorio
    !! @param iname         (Opcional) Nome do arquivo para escrita
    function rngNormalArray(mean, sd, idimArray, iname)

        implicit none

        real*8, intent(in) :: mean, sd
        character(len=*), intent(in), optional :: iname
        real*8, allocatable :: rngNormalArray(:)
        integer, intent(in), optional :: idimArray
        integer :: dimArray = 10, i, id
        integer :: rseed1, n
        integer, allocatable :: seed(:)
        character(len=50) :: filename, fname, char_id, idname

        if (present(idimArray)) dimArray = idimArray
        if (present(iname)) then 
            idname = trim(iname)
            ! Opening solution output file
            ! Making an id char variable
            id = 7492642
            write(filename, '(a)') "normDist"
            ! Concatenate strings filename, char_id and ".dat"
            fname = trim(filename)//trim(idname)//".dat"
            ! Opening file with unit as id and file name as attributed in fname
            open(unit=id,file=fname)
        endif

        allocate(rngNormalArray(dimArray)); rngNormalArray = 0.d0

        call random_seed(size = n)
        allocate(seed(n))
        call random_seed(put=seed)
        rseed1 = irand(); !rseed2 = irand()

        do i=1,dimArray
            rngNormalArray(i) = r8_normal_ab(mean, sd, rseed1)
            if (present(iname)) write(id,*) rngNormalArray(i)
        enddo

        if (present(iname)) close(id)

    endfunction

    function r8_normal_ab ( a, b, seed )

        !*****************************************************************************80
        !
        !! R8_NORMAL_AB returns a scaled pseudonormal R8.
        !
        !  Discussion:
        !
        !    The normal probability distribution function (PDF) is sampled,
        !    with mean A and standard deviation B.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    06 August 2013
        !    29 September 2017 by Diego T. Volpatto (volpatto@lncc.br)
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) A, the mean of the PDF.
        !
        !    Input, real ( kind = 8 ) B, the standard deviation of the PDF.
        !
        !    Input/output, integer ( kind = 4 ) SEED, a seed for the random
        !    number generator.
        !
        !    Output, real ( kind = 8 ) R8_NORMAL_AB, a sample of the normal PDF.
        !
        implicit none

        real ( kind = 8 ) a
        real ( kind = 8 ) b
        real ( kind = 8 ) r1
        real ( kind = 8 ) r2
        real ( kind = 8 ) r8_normal_ab
        real ( kind = 8 ), parameter :: r8_pi = 4.d0*datan(1.d0)
        !real ( kind = 8 ) r8_uniform_01
        integer ( kind = 4 ) seed
        real ( kind = 8 ) x

        r1 = r8_uniform_01 ( seed )
        r2 = r8_uniform_01 ( seed )
        x = dsqrt ( - 2.0D+00 * dlog ( r1 ) ) * dcos ( 2.0D+00 * r8_pi * r2 )

        r8_normal_ab = a + b * x

        return
    end
    function r8_uniform_01 ( seed )

        !*****************************************************************************80
        !
        !! R8_UNIFORM_01 returns a unit pseudorandom R8.
        !
        !  Discussion:
        !
        !    This routine implements the recursion
        !
        !      seed = 16807 * seed mod ( 2^31 - 1 )
        !      r8_uniform_01 = seed / ( 2^31 - 1 )
        !
        !    The integer arithmetic never requires more than 32 bits,
        !    including a sign bit.
        !
        !    If the initial seed is 12345, then the first three computations are
        !
        !      Input     Output      R8_UNIFORM_01
        !      SEED      SEED
        !
        !         12345   207482415  0.096616
        !     207482415  1790989824  0.833995
        !    1790989824  2035175616  0.947702
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    31 May 2007
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Reference:
        !
        !    Paul Bratley, Bennett Fox, Linus Schrage,
        !    A Guide to Simulation,
        !    Second Edition,
        !    Springer, 1987,
        !    ISBN: 0387964673,
        !    LC: QA76.9.C65.B73.
        !
        !    Bennett Fox,
        !    Algorithm 647:
        !    Implementation and Relative Efficiency of Quasirandom
        !    Sequence Generators,
        !    ACM Transactions on Mathematical Software,
        !    Volume 12, Number 4, December 1986, pages 362-376.
        !
        !    Pierre L'Ecuyer,
        !    Random Number Generation,
        !    in Handbook of Simulation,
        !    edited by Jerry Banks,
        !    Wiley, 1998,
        !    ISBN: 0471134031,
        !    LC: T57.62.H37.
        !
        !    Peter Lewis, Allen Goodman, James Miller,
        !    A Pseudo-Random Number Generator for the System/360,
        !    IBM Systems Journal,
        !    Volume 8, 1969, pages 136-143.
        !
        !  Parameters:
        !
        !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
        !    should NOT be 0.
        !    On output, SEED has been updated.
        !
        !    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
        !    strictly between 0 and 1.
        !
        implicit none

        integer ( kind = 4 ) k
        real ( kind = 8 ) r8_uniform_01
        integer ( kind = 4 ) seed

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
            seed = seed + 2147483647
        end if
        !
        !  Although SEED can be represented exactly as a 32 bit integer,
        !  it generally cannot be represented exactly as a 32 bit real number!
        !
        r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

        return

        end
    function r8_uni ( s1, s2 )

        !*****************************************************************************80
        !
        !! R8_UNI returns a pseudorandom number between 0 and 1.
        !
        !  Discussion:
        !
        !    This function generates uniformly distributed pseudorandom numbers
        !    between 0 and 1, using the 32-bit generator from figure 3 of
        !    the article by L'Ecuyer.
        !
        !    The cycle length is claimed to be 2.30584E+18.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    08 July 2008
        !
        !  Author:
        !
        !    Original Pascal original version by Pierre L'Ecuyer
        !    FORTRAN90 version by John Burkardt
        !
        !  Reference:
        !
        !    Pierre LEcuyer,
        !    Efficient and Portable Combined Random Number Generators,
        !    Communications of the ACM,
        !    Volume 31, Number 6, June 1988, pages 742-751.
        !
        !  Parameters:
        !
        !    Input/output, integer ( kind = 4 ) S1, S2, two values used as the
        !    seed for the sequence.  On first call, the user should initialize
        !    S1 to a value between 1 and 2147483562;  S2 should be initialized
        !    to a value between 1 and 2147483398.
        !
        !    Output, real ( kind = 8 ) R8_UNI, the next value in the sequence.
        !
        implicit none

        integer ( kind = 4 ) k
        real ( kind = 8 ) r8_uni
        integer ( kind = 4 ) s1
        integer ( kind = 4 ) s2
        integer ( kind = 4 ) z

        k = s1 / 53668
        s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
        if ( s1 < 0 ) then
            s1 = s1 + 2147483563
        end if

        k = s2 / 52774
        s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
        if ( s2 < 0 ) then
            s2 = s2 + 2147483399
        end if

        z = s1 - s2
        if ( z < 1 ) then
            z = z + 2147483562
        end if

        r8_uni = real ( z, kind = 8 ) / 2147483563.0D+00

        return
    end
    endmodule
