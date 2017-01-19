module wavefield
    implicit none
    type, public :: wavefield_class
        integer :: nx, ny, nz
        real :: dx, dy, dz
        complex,allocatable,dimension(:,:,:) :: dp
        logical :: dp_allocated = .false.
    contains
        procedure,public :: init => init_sub
    endtype
    private :: init_sub
contains
    ! ======================================================
    ! subroutine to initialize the wavefield
    ! ======================================================
    subroutine init_sub(this) ! {
        implicit none
        class (wavefield_class) :: this
        integer :: istat

        if (this%dp_allocated) then
            deallocate (this%dp)
            this%dp_allocated = .false.
        endif

        allocate (this%dp(this%nx,this%ny,this%nz),stat=istat)
        if (istat/=0) then
            print*, 'wavefield: Error in allocate dp.'
            return
        endif
        this%dp_allocated = .true.
        this%dp = cmplx(0.0,0.0)
    end subroutine
    ! }
end module
