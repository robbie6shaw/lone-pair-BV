! FILE: calc.f90

subroutine distance(point1, point2, rcut, result)
!
! Calculates the distance between two points. If the the distance on one axis exceeds the cutoff distance, only that axes distance is returned. 
!       
    implicit none
    real, dimension(3), intent(in) :: point1
    real, dimension(3), intent(in) :: point2
    real, intent(in) :: rcut
    real, intent(out) :: result

    integer :: i
    real :: delta
    real :: sum
    integer :: cut

    cut = 0
    sum = 0
    do i = 1,3
        delta = abs(point2(i) - point1(i))
        if (delta > rcut) then
            result = delta
            cut = 1
            exit 
        else
            sum = sum + delta ** 2
        endif
    enddo

    if (cut == 1) then
        result = sqrt(sum)
    endif

end subroutine distance
