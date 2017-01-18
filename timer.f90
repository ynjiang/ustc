module timer
    implicit none

    integer, parameter :: dbl = selected_real_kind(p=14)
    type, public :: timer_class
        private
        real(kind=dbl) :: saved_time ! saved time in ms
    contains
        procedure, public :: start_timer => start_timer_sub
        procedure, public :: elapsed_time => elapsed_time_fn
    endtype timer_class

    ! restrict control to the real subroutine names
    private :: start_timer_sub, elapsed_time_fn

contains
    ! ==============================================
    ! subroutine to get and save the initial time
    ! ==============================================
    subroutine start_timer_sub(this)
    implicit none
    class(timer_class) :: this
    integer,dimension(8) :: value
    
    ! Get time
    call date_and_time (values=value)
    this%saved_time = 86400.D0 * value(3) + 3600.D0 * value(5) &
                + 60.D0 * value(6) + value(7) + 0.001D0 * value(8)
    end subroutine start_timer_sub

    ! ====================================
    ! function to calculate elapsed time
    ! ====================================
    real function elapsed_time_fn(this)
    implicit none
    class(timer_class) :: this
    integer,dimension(8) :: value
    real(kind=dbl) :: current_time ! Current time (ms)
    
    ! Get time
    call date_and_time (values=value)
    current_time = 86400.D0 * value(3) + 3600.D0 * value(5) &
                + 60.D0 * value(6) + value(7) + 0.001D0 * value(8)
    ! Get elapsed ime in seconds
    elapsed_time_fn = current_time - this%saved_time
    end function elapsed_time_fn
end module timer
