! =============================================
! 
! Prestack depth migration
! 
! 3-d superwide-angle one-way propagator
!
!=============================================

program main
    use timer
    use input
    use velmod
    use dataset
    use source
    use wavefield
    implicit none 

    type (timer_class) :: time
    type (input_class) :: reader
    type (velmod_class) :: vm_raw, vm_glo, vm_loc(4), vmr_loc
    type (dataset_class) :: ur
    type (source_class) :: us
    type (wavefield_class) :: wf

    integer :: ishot, iw, iz, idirec
    integer*8 :: plan, planinv, ploc(4), plocinv(4), prloc(4), prlocinv(4)

    print*, "======================================="
    print*, "main: program begins."
    call time % start_timer

    ! input parameters
    reader % parfile = 'partable'
    call reader % read_input
    print*, 'main: read parameters from partable successfully.'
    print*, '---------------------------------------'

    ! initialize the source
    us % fmin = reader % fmin
    us % fmax = reader % fmax
    us % nt   = reader % nt
    us % dt   = reader % dt
    us % fp   = reader % fp
    us % amp  = reader % amp
    us % t0   = reader % t0
    call us % ricker
    print*, 'main: initialize the source successfully.'
    print*, '---------------------------------------'

    ! read velocity model from velfile
    vm_raw % nx = reader % nx
    vm_raw % ny = reader % ny
    vm_raw % nz = reader % nz
    vm_raw % dx = reader % dx
    vm_raw % dy = reader % dy
    vm_raw % dz = reader % dz
    vm_raw % velfile = reader % velfile
    call vm_raw % init_from_file
    print*, 'main: read raw velmod from velfile successfully.'
    print*, '---------------------------------------'

    ! generate expanded veloctity from raw velocity
    ! for each z direction, padding is only in the x-y plane.
    call vm_glo % init_expand (vm_raw, reader % padding)
    print*, 'main: expand raw velmod to global velmod successfully.'
    print*, '---------------------------------------'

    if (vm_raw % vel_allocated) then
        deallocate (vm_raw % vel)
        vm_raw % vel_allocated = .false.
    endif

    call fftw2d_f77_create_plan (plan, vm_glo % nx, vm_glo % ny, -1, 1)
    call fftw2d_f77_create_plan (planinv, vm_glo % nx, vm_glo % ny, 1, 1)

    wf % nx = vm_glo % nx
    wf % ny = vm_glo % ny
    wf % nz = vm_glo % nz
    wf % dx = vm_glo % dx
    wf % dy = vm_glo % dy
    wf % dz = vm_glo % dz

    call wf % init
    frequency_loop: do iw = 1, us % nwnew
        if (mod(iw,20) == 0) print*, 'main: iw = ', iw
        wf % dp = cmplx(0.0,0.0)
    enddo frequency_loop
    
    call fftwnd_f77_destroy_plan (planinv)
    call fftwnd_f77_destroy_plan (plan)

    if (vm_glo % vel_allocated) then
        deallocate (vm_glo % vel)
        vm_glo % vel_allocated = .false.
    endif
    if (wf % dp_allocated) then
        deallocate (wf % dp)
        wf % dp_allocated = .false.
    endif

    print*, 'main: Done!'
    print*, 'main: time total: ', time % elapsed_time(), ' s'
    print*, "======================================="
end program
