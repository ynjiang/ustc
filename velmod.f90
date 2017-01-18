module velmod
    implicit none
    type, public :: velmod_class
        character(len=200) :: velfile ! only used in init_sub
        integer :: nx, ny, nz
        real :: dx, dy, dz ! (m)
        real,allocatable,dimension(:,:,:) :: vel
        logical :: vel_allocated = .false.
        ! shot location (single shot)
        integer :: sx, sy
        ! receiver locations (multiple receivers)
        integer :: num_of_rec=0
        integer,allocatable,dimension(:) :: rx, ry
        logical :: rx_allocated = .false.
        logical :: ry_allocated = .false.
        real,allocatable,dimension(:) :: vb
        logical :: vb_allocated
    contains
        procedure, public :: init_from_file => init_from_file_sub
        procedure, public :: init_rec_location => init_rec_location_sub
        procedure, public :: init_expand => init_expand_sub
        procedure, public :: get_vb_min => get_vb_min_sub
        procedure, public :: get_vb_max => get_vb_max_sub
        procedure, public :: get_vb_avg => get_vb_avg_sub
    endtype
    ! }
    private :: init_from_file_sub
    private :: init_rec_location_sub
    private :: init_expand_sub
    private :: get_vb_min_sub
    private :: get_vb_max_sub
    private :: get_vb_avg_sub

contains
    ! =======================================================
    ! subroutine to initialize the velocity model from file
    ! =======================================================
    subroutine init_from_file_sub(this) ! {
        implicit none
        class (velmod_class) :: this
        integer :: istat, irec, ix, iy, iz

        if (this%vel_allocated) then
            deallocate (this%vel)
            this%vel_allocated = .false.
        endif

        allocate (this%vel(this%nx,this%ny,this%ny),stat=istat)
        if (istat/=0) then
            print*, 'Error in allocating vel.'
            stop
        endif
        this%vel_allocated = .true.

        open (13,file=trim(this%velfile),status='old',recl=this%nz,form='unformatted',access='direct',iostat=istat)
        if (istat/=0) then
            print*, trim(this%velfile)
            print*, 'Error in opening velfile.'
            deallocate (this%vel)
            stop
        endif
        irec = 0
         do ix = 1,this%nx
            do iy = 1,this%ny
                irec = irec+1
                read (13,rec=irec) (this%vel(ix,iy,iz),iz=1,this%nz)
            enddo
        enddo
        close (13)
    end subroutine
    ! }

    ! =====================================================
    ! subroutine to initialize the receivers' locations
    ! =====================================================
    subroutine init_rec_location_sub(this, rx, ry) ! {
        implicit none
        class (velmod_class) :: this
        integer,dimension(:),intent(in) :: rx, ry
        
        if (size(rx,1)/=size(ry,1)) then
            print*, 'Input rx and ry is not in the same size.'
            return
        endif
        this%num_of_rec = size(rx,1)
        
        if (this%rx_allocated) then
            deallocate (this%rx)
            this%rx_allocated = .false.
        endif
        
        if (this%ry_allocated) then
            deallocate (this%ry)
            this%ry_allocated = .false.
        endif

        allocate (this%rx(this%num_of_rec))
        this%rx_allocated = .true.
        
        allocate (this%ry(this%num_of_rec))
        this%ry_allocated = .true.

        this%rx = rx
        this%ry = ry

    end subroutine
    ! }

    ! ===================================================
    ! subroutine to initialize the velocity by expanding
    ! ===================================================
    subroutine init_expand_sub (this, vm, padding) ! {
        ! the padding is only used in x-y plane, not in z-direction
        implicit none
        class (velmod_class) :: this
        class (velmod_class), intent(in) :: vm
        integer, intent(in) :: padding
        integer :: pad_nx, pad_ny
        integer :: pad_px, pad_py
        integer :: ix, iy, iz, i
        integer :: istat

        pad_nx = padding
        pad_ny = padding
        pad_px = padding
        pad_py = padding

        this%nx = vm%nx + pad_nx + pad_px
        this%ny = vm%ny + pad_ny + pad_py
        this%nz = vm%nz

        this%dx = vm%dx
        this%dy = vm%dy
        this%dz = vm%dz

        this%sx = vm%sx + pad_nx
        this%sy = vm%sy + pad_ny

        if (this%rx_allocated) then
            deallocate (this%rx)
            this%rx_allocated = .false.
        endif
        if (this%ry_allocated) then
            deallocate (this%ry)
            this%ry_allocated = .false.
        endif

        this%num_of_rec = vm%num_of_rec
        if (vm%rx_allocated) then
            allocate (this%rx(this%num_of_rec), stat=istat)
            if (istat/=0) then
                print*, 'velmod: Error in allocate rx in expand.'
                return
            endif
            this%rx_allocated = .true.
        endif
        if (vm%ry_allocated) then
            allocate (this%ry(this%num_of_rec), stat=istat)
            if (istat/=0) then
                print*, 'velmod: Error in allocate ry in expand.'
                return
            endif
            this%ry_allocated = .true.
        endif

        print*, 'velmod: number of receivers: ', this%num_of_rec
        !do i=1,this%num_of_rec
        !    this%rx(i) = vm%rx(i) + pad_nx
        !    this%ry(i) = vm%ry(i) + pad_nx
        !enddo

        if (this%vel_allocated) then
            deallocate (this%vel)
            this%vel_allocated = .false.
        endif
        allocate (this%vel(this%nx,this%ny,this%nz),stat=istat)
        if (istat/=0) then
            print*, 'velmod: Error in allocate vel.'
            return
        endif
        this%vel_allocated = .true.
        do iz=1,this%nz
            do iy=pad_ny+1,pad_ny+vm%ny
                do ix=1,pad_nx
                    this%vel(ix,iy,iz)=vm%vel(1,iy-pad_ny,iz)
                enddo
                do ix=pad_nx+1,pad_nx+vm%nx
                    this%vel(ix,iy,iz)=vm%vel(ix-pad_nx,iy-pad_ny,iz)
                enddo
                do ix=pad_nx+vm%nx+1,this%nx
                    this%vel(ix,iy,iz)=vm%vel(vm%nx,iy-pad_ny,iz)
                enddo
            enddo
            do ix=1,this%nx
                do iy=1,pad_ny
                    this%vel(ix,iy,iz)=this%vel(ix,pad_ny+1,iz)
                enddo
                do iy=pad_ny+vm%ny+1,this%ny
                    this%vel(ix,iy,iz)=this%vel(ix,pad_ny+vm%ny,iz)
                enddo
            enddo
        enddo
    end subroutine
    ! }

    ! ==================================================
    ! subroutine to get background velocity - min
    ! ==================================================
    subroutine get_vb_min_sub (this) ! {
        implicit none
        class (velmod_class) :: this
        integer :: istat, ix, iy, iz
        real :: tmp

        if (this%vb_allocated) then
            deallocate (this%vb)
            this%vb_allocated = .false.
        endif

        allocate (this%vb(this%nz), stat=istat)
        if (istat/=0) then
            print*, 'Error in allocating vb_min.'
            return
        endif

        do iz=1,this%nz
            tmp = 100000.0
            do ix=1,this%nx
                do iy=1,this%ny
                    if (this%vel(ix,iy,iz) < tmp) tmp=this%vel(ix,iy,iz)
                enddo
            enddo
            this%vb(iz) = tmp
        enddo
    end subroutine
    ! }
    
    ! ==================================================
    ! subroutine to get background velocity - max
    ! ==================================================
    subroutine get_vb_max_sub (this) ! {
        implicit none
        class (velmod_class) :: this
        integer :: istat, ix, iy, iz
        real :: tmp

        if (this%vb_allocated) then
            deallocate (this%vb)
            this%vb_allocated = .false.
        endif

        allocate (this%vb(this%nz), stat=istat)
        if (istat/=0) then
            print*, 'Error in allocating vb_max.'
            return
        endif

        do iz=1,this%nz
            tmp = -100000.0
            do ix=1,this%nx
                do iy=1,this%ny
                    if (this%vel(ix,iy,iz) > tmp) tmp=this%vel(ix,iy,iz)
                enddo
            enddo
            this%vb(iz) = tmp
        enddo

    end subroutine
    ! }

    ! ==================================================
    ! subroutine to get background velocity - avg
    ! ==================================================
    subroutine get_vb_avg_sub (this) ! {
        implicit none
        class (velmod_class) :: this
        integer :: istat, ix, iy, iz
        real :: tmp

        if (this%vb_allocated) then
            deallocate (this%vb)
            this%vb_allocated = .false.
        endif

        allocate (this%vb(this%nz), stat=istat)
        if (istat/=0) then
            print*, 'Error in allocating vb_avg.'
            return
        endif

        do iz=1,this%nz
            tmp = 0.0
            do ix=1,this%nx
                do iy=1,this%ny
                    tmp = tmp + this%vel(ix,iy,iz)
                enddo
            enddo
            this%vb(iz) = tmp/this%nx/this%ny
        enddo

    end subroutine
    ! }

end module
