module dataset
    implicit none
    
    ! define of type trace_class {
    type,public :: trace_class
    endtype
    ! }

    ! define of type dataset_class {
    type,public :: dataset_class
        character(len=200) :: datafile
        integer :: ns ! number of samples
        real :: ds ! interval of samples
    contains
        procedure, public :: read_data => read_data_sub
    endtype
    ! }

    !define of type dataset_class with su header {
    type,public :: dataset_class_su
    endtype
    ! }
    
    private :: read_data_sub
contains
    ! =======================================
    ! subroutine to read data from datfile
    ! =======================================
    subroutine read_data_sub (this)
        implicit none
        class (dataset_class) :: this
        
    end subroutine
end module
