module source
    implicit none
    ! define of type source_class {
    type, public :: source_class
        real :: dt,fmin,fmax
        integer :: nt,nwmin,nwmax,nwnew
        real :: fp, amp, t0
        real,allocatable,dimension(:) :: w_time
        logical :: w_time_allocated = .false.
        complex,allocatable,dimension(:) :: w_freq
        logical :: w_freq_allocated = .false.
    contains
        procedure, public :: ricker => ricker_sub
    endtype
    ! }
    private :: ricker_sub

contains
    ! ========================================
    ! subroutine to generate ricker serial
    ! ========================================
    subroutine ricker_sub(this)
    implicit none
    class (source_class) :: this                                                                             
    integer :: it0,i,istat,it
    real    :: pi,a,wmax,coef,t0s
    complex cnum,im
    integer*8 plansource
    real :: df,dw,fwmin,fwmax
    complex,allocatable,dimension(:) :: work,workin
    character(len=100) :: srcfile_t = '~/data/1/srcfile_time.dat'
    character(len=100) :: srcfile_f = '~/data/1/srcfile_freq.dat'

    print*, ''
    print*, 'source: members that need to be initialized:'
    print*, 'source: fmin ', this%fmin, ', fmax ', this%fmax
    print*, 'source: nt ', this%nt, ', dt ', this%dt
    print*, 'source: fp ', this%fp, ', amp ', this%amp, ', t0 ', this%t0

    ! ----generate time-damain sequence----
    if (this%w_time_allocated) then
        deallocate (this%w_time)
        this%w_time_allocated = .false.
    endif
    allocate (this%w_time(this%nt),stat=istat)
    if (istat/=0) then
        print*, 'Error when allocate w_time'
        stop
    endif
    this%w_time_allocated = .true.
                       
    pi=3.1415926
    t0s=1.0/this%fp
    it0=int(t0s/this%dt+1.e-6)+1
    t0s=(it0-1.0)*this%dt

    if (this%nt.lt.2*it0) then
        write(*,*) 'Error: nt should be >=',2*it0,'!!!'
        stop
    endif
                                                                                
    wmax=0.0
    do i=1,this%nt
        a=pi*this%fp*((i-1)*this%dt-t0s)
        a=a*a
        this%w_time(i)=(1.-2.*a)*exp(-a)
        wmax=max(abs(this%w_time(i)),wmax)
    enddo
    if (abs(wmax)<1.e-6) then
        print*, "Error: wmax should not be 0"
        stop
    endif
    coef=this%amp/wmax
    do i=1,this%nt
        this%w_time(i)=coef*this%w_time(i)
    enddo

    open (11,file=srcfile_t,status='replace')
    do it=1,this%nt
        write(11,*) this%w_time(it)
    enddo
    close (11)
    
    ! ----generate frequency-domain sequence----
    df = 1.0/(this%nt*this%dt)
    dw = 2.0*pi*df
    this%nwmin = max(2,ifix(this%fmin/df)+1)
    this%nwmax = min(this%nt/2,ifix(this%fmax/df)+1)
    fwmin = (this%nwmin-1)*df
    fwmax = (this%nwmax-1)*df
    this%nwnew = this%nwmax-this%nwmin+1
    if (this%w_freq_allocated) then
        deallocate (this%w_freq)
        this%w_freq_allocated = .false.
    endif
    allocate (this%w_freq(this%nwnew),stat=istat)
    if (istat/=0) then
        print*, 'Error when allocate w_freq'
        stop
    endif
    this%w_freq_allocated = .true.
    
    allocate (work(this%nt))
    allocate (workin(this%nt))

    do it=1,this%nt
        workin(it) = cmplx(this%w_time(it),0.0)
    enddo
    call fftw_f77_create_plan (plansource,this%nt,1,0)
    call fftw_f77_one (plansource,workin,work)
    call fftw_f77_destroy_plan (plansource)
    
    im = cmplx(0.0,1.0)
    if (t0s/=0.0) then
        print*, 'source: t0s ', t0s
        do it=1,this%nt/2+1
            cnum=cexp(-im*(it-1)*dw*t0s)
            work(it)=work(it)*cnum
        enddo
    endif
    if (this%t0/=0.0) then
        print*, 'source: t0 ',this%t0
        do it=1,this%nt/2+1
            cnum=cexp(im*(it-1)*dw*this%t0)
            work(it)=work(it)*cnum
        enddo
    endif
    
    open (11,file=srcfile_f,status='replace')
    do it=1,this%nt
        write(11,*) work(it), df*(it-1), abs(work(it))
    enddo
    close (11)

    this%w_freq(1:this%nwnew)=work(this%nwmin:this%nwmax)
    deallocate (work, workin)

    print*, ''
    print*, 'source: menbers calculated:'
    print*, 'source: nwmin ', this%nwmin, ', nwmax ', this%nwmax, ', nwnew', this%nwnew
    print*, 'source: time-domain source is stored in ', trim(srcfile_t)
    print*, 'source: frequency-domain source is stored in ', trim(srcfile_f)
    end subroutine 
end module
