module tIntegration

    implicit none

    contains

    !> @param ysol  Vetor solucao 
    !! @param t0    Tempo inicial
    !! @param tf    Tempo final
    !! @param y0    Valor inicial
    !! @param cteR  Resistencia
    !! @param cteC  Capacitancia
    !! @param steps Quantidade de passos de tempo
    !! @param id    (Opcional) Id do arquivo. Deve ser individual para cada caso.
    subroutine solverEuler (ysol, t0, tf, y0, cteR, cteC, steps, id)

        implicit none

        real*8, intent(in) :: t0, y0, cteR, cteC, tf
        integer, optional :: id
        real*8, allocatable :: ysol(:)
        real*8 :: yprev, tprev, t, dt
        integer, intent(in) :: steps
        integer :: i, nnodes
        character(len=50) :: filename, fname, char_id

        nnodes = steps+1
        allocate(ysol(nnodes)); ysol=0.d0

        yprev = y0; tprev = t0

        if (present(id)) then
            ! Opening solution output file
            ! Making an id char variable
            write(char_id, '(i0)') id
            write(filename, '(a)') "solutionE"
            ! Concatenate strings filename, char_id and ".dat"
            fname = trim(filename)//trim(char_id)//".dat"
            ! Opening file with unit as id and file name as attributed in fname
            open(unit=id,file=fname)

            ! Writing initial condition
            write(id,2000) t0,y0; 
        endif

        ! Computing the time step
        dt = (tf-t0)/steps

        ! Solver Euler
        ysol(1) = y0
        do i=2,nnodes
            t=tprev+dt; 
            ysol(i)=yprev/(1.d0+dt/(cteR*cteC))
            yprev = ysol(i) 
            tprev = t
            if (present(id)) write(id,2000) t, ysol(i)
        end do

        if (present(id)) close(id); !deallocate(ysol)

1000    format((f10.4),4(4x,(f20.8)))
2000    format((f10.4),1(4x,(f20.8)))

    end subroutine

    !> @param ysol      Vetor solucao 
    !! @param t0        Tempo inicial
    !! @param tf        Tempo final
    !! @param y0        Valor inicial
    !! @param cteR      Resistencia
    !! @param cteC      Capacitancia
    !! @param steps     Quantidade de passos de tempo
    !! @param id        (Opcional) Id do arquivo. Deve ser individual para cada caso.
    !! @param irefresh  (Opcional) Flag de desalocação do vetor solução.
    subroutine solverRK4 (ysol, t0, tf, y0, cteR, cteC, steps, id, irefresh)

        implicit none

        real*8, intent(in) :: t0, y0, cteR, cteC, tf
        integer, optional :: id, irefresh
        real*8, allocatable :: ysol(:)
        real*8 :: yprev, tprev, t, dt, k1, k2, k3, k4
        integer, intent(in) :: steps
        integer :: i, nnodes, refresh=0
        character(len=50) :: filename, fname, char_id

        if (present(irefresh)) refresh = irefresh

        nnodes = steps+1
        allocate(ysol(nnodes)); ysol=0.d0

        yprev = y0; tprev = t0

        if (present(id)) then
            ! Opening solution output file
            ! Making an id char variable
            write(char_id, '(i0)') id
            write(filename, '(a)') "solutionRK"
            ! Concatenate strings filename, char_id and ".dat"
            fname = trim(filename)//trim(char_id)//".dat"
            ! Opening file with unit as id and file name as attributed in fname
            open(unit=id,file=fname)

            ! Writing initial condition
            write(id,2000) t0,y0; 
        endif

        ! Computing the time step
        dt = (tf-t0)/steps

        ! Solver Euler
        ysol(1) = y0
        do i=2,nnodes
            t=tprev+dt; 
            k1 = -yprev/(cteR*cteC)
            k2 = -(yprev+dt*k1/2.d0)/(cteR*cteC)
            k3 = -(yprev+dt*k2/2.d0)/(cteR*cteC)
            k4 = -(yprev+dt*k3)/(cteR*cteC)
            ysol(i)=yprev + dt*(k1 + 2.d0*k2 + 2.d0*k3 + k4)/6.d0
            yprev = ysol(i) 
            tprev = t
            if (present(id)) write(id,2000) t, ysol(i)
        end do

        if (present(id)) close(id); 
        if (refresh .ne. 0) deallocate(ysol)

1000    format((f10.4),4(4x,(f20.8)))
2000    format((f10.4),1(4x,(f20.8)))

    end subroutine

endmodule
