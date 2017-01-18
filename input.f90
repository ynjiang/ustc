module input
    implicit none
    
    ! define of type input_class {
    type,public :: input_class
        character(len=200) :: parfile
        integer :: verbose
        integer :: nx, ny, nz
        real :: dx, dy, dz
        integer :: nt, ntfft
        real :: dt
        real :: fmin, fmax
        real :: fp, amp, t0
        integer :: padding
        character(len=200) :: velfile
        character(len=200) :: datpref
        character(len=200) :: infofile
        character(len=200) :: snappref
        character(len=200) :: imgpref
    contains
        procedure, public :: read_input => read_input_sub
    endtype
    ! }
    private :: read_input_sub

contains
    ! ==============================================
    ! subroutine to read parameters from partable
    ! ==============================================
    subroutine read_input_sub (this)
        implicit none
        class (input_class) :: this
        integer :: error
        character(len=200) :: chartmp

        open (11,file=trim(this%parfile),status="old",iostat=error)
        if (error/=0) then
            print*, 'error in opening ', trim(this%parfile)
            close (11)
            stop
        endif
        read (11,*) chartmp
        read (11,*) this%verbose
        read (11,*) chartmp
        read (11,*) this%nx, this%ny, this%nz
        read (11,*) chartmp
        read (11,*) this%dx, this%dy, this%dz
        read (11,*) chartmp
        read (11,*) this%nt, this%ntfft, this%dt
        read (11,*) chartmp
        read (11,*) this%fmin, this%fmax
        read (11,*) chartmp
        read (11,*) this%fp, this%amp, this%t0
        read (11,*) chartmp
        read (11,*) this%padding
        read (11,*) chartmp
        read (11,'(A100)') this%velfile
        read (11,*) chartmp
        read (11,'(A100)') this%datpref
        read (11,*) chartmp
        read (11,'(A100)') this%infofile
        read (11,*) chartmp
        read (11,'(A100)') this%snappref
        read (11,*) chartmp
        read (11,'(A100)') this%imgpref
        close (11)

        print*, '      verbose'
        write (*,'(I11)') this%verbose
        print*, '        nx        ny         nz'
        write (*,'(3I10)') this%nx, this%ny, this%nz
        print*, '        dx        dy         dz'
        write (*,'(3F10.1)') this%dx, this%dy, this%dz
        print*, '        nt       ntfft       dt'
        write (*,'(2I10,F6.3)') this%nt, this%ntfft, this%dt
        print*, '       fmin       fmax'
        write (*,'(2F10.1)') this%fmin, this%fmax
        print*, '        fp         amp      t0'
        write (*,'(3F10.1)') this%fp, this%amp, this%t0
        print*, '      padding'
        write (*,'(I10)') this%padding
        print*, "velfile:"
        write (*,*) trim(this%velfile)
        print*, "datpref:"
        write (*,*) trim(this%datpref)
        print*, "infofile:"
        write (*,*) trim(this%infofile)
        print*, "snappref:"
        write (*,*) trim(this%snappref)
        print*, "imgpref:"
        write (*,*) trim(this%imgpref)
    end subroutine
end module
