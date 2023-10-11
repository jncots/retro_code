module traj_bulk_m
!================================================
! Loading calculated data
! Example of usage in "test_load" subroutine
!================================================
 use ode_solver_m, only : ode_solution, ode_solver
 use search_data_m, only : search_data
 implicit none
 private
 save

 public :: traj_bulk
! public :: test_calc

 type traj_bulk
  logical :: no_fdata=.true.
  logical :: not_load=.true.
  logical :: change_file=.false.
  character(500) :: dname, prefix
  character(500) :: fname=''
  integer :: ifile=0, itraj=0
  integer :: nfile=0, ntraj=0
  type(search_data) :: sdata
  type(ode_solution), allocatable :: sol(:)
 contains
  procedure :: set=>set_tb, del=>del_tb
  procedure :: read_next, load_traj, find_files=>find_tb
 end type traj_bulk



 contains


subroutine del_tb(this)
 class(traj_bulk) :: this

 this%no_fdata=.true.
 this%not_load=.true.
 call this%sdata%del
 if (allocated(this%sol)) deallocate(this%sol)
 this%ifile=0 
 this%itraj=0
 this%nfile=0 
 this%ntraj=0
 
end subroutine del_tb


subroutine find_tb(this)
 class(traj_bulk) :: this

 call this%del

 this%sdata%search_dir=trim(this%dname)
 this%sdata%prefix=trim(this%prefix)
 call this%sdata%search
 this%nfile=this%sdata%nfiles
 this%no_fdata=.false.

end subroutine find_tb


subroutine set_tb(this)
 class(traj_bulk) :: this

 if (this%no_fdata) call this%find_files
 this%not_load=.false.

 this%ifile=1 
 this%itraj=1

 this%fname=this%sdata%dfile(this%ifile)
 call read_traj(this%fname,this%sol)
 this%change_file=.true.
 this%ntraj=size(this%sol)

end subroutine set_tb




function load_traj(this,sol)
 class(traj_bulk) :: this
 type(ode_solution) :: sol
 logical :: load_traj, no_next

 if (this%not_load) then
  call this%set
 else
  call make_step(this,no_next)
  if (no_next) then
   load_traj=.false.
   call sol%del
   return
  end if
 end if

 if (this%sol(this%itraj)%np>0) then
  call sol%copy_from(this%sol(this%itraj))
 else
  do
   call make_step(this,no_next)
   
   if (no_next) then
    load_traj=.false.
    call sol%del
    return
   end if
   
   if (this%sol(this%itraj)%np>0) then
    call sol%copy_from(this%sol(this%itraj))
    exit  
   end if
  
  end do
 end if
 
 load_traj=.true. 

end function load_traj


subroutine read_next(this,sol,no_next)
!===========================================
! Old subroutine, use load_traj !!!!
!===========================================
 class(traj_bulk) :: this
 type(ode_solution) :: sol
 logical, intent(out) :: no_next


 if (this%not_load) call this%set

 sol%np=0
 if (this%sol(this%itraj)%np>0) then
  call sol%copy_from(this%sol(this%itraj))
  call make_step(this,no_next)
 else
  do
   if (this%sol(this%itraj)%np>0) then
    call sol%copy_from(this%sol(this%itraj))
    call make_step(this,no_next)
    exit
   end if
   call make_step(this,no_next)
   if (no_next) exit
  end do 
 end if
 

end subroutine read_next


subroutine make_step(this,no_next)
 class(traj_bulk) :: this
 logical, intent(out) :: no_next

 this%change_file=.false.
 this%itraj=this%itraj+1
 
 if (this%itraj>this%ntraj) then
  this%ifile=this%ifile+1
  
  if (this%ifile>this%nfile) then
   no_next=.true.
   return
  else
   this%fname=this%sdata%dfile(this%ifile)
   call read_traj(this%fname,this%sol)
   this%change_file=.true.
   this%ntraj=size(this%sol)
   this%itraj=1
  end if
 end if

 no_next=.false.

end subroutine make_step



subroutine read_traj(fname,sol)
 character(500), intent(in) :: fname
 type(ode_solution), allocatable :: sol(:)
 integer :: ndev, nsol, i, iost

 open(newunit=ndev,file=fname,form='unformatted',access='stream',status='old')
 read(ndev,iostat=iost) nsol

 if (iost<0) return

 if (allocated(sol)) deallocate(sol)
 allocate(sol(nsol))
 do i=1,nsol
  call sol(i)%read(ndev)
 end do

end subroutine read_traj


!subroutine test_load
! type(traj_bulk) :: tbulk
! type(ode_solution) :: sol
! logical :: no_next
! integer :: j

! tbulk%dname='/home/anton/work/project/halo/chp_prop/in_bs_field/plane_z20r20/bfs/bfs_17_cs510_50/'
! tbulk%prefix='traj_part'

! j=0
! do
!  call tbulk%read_next(sol,no_next)
!  j=j+1
!!  write(*,*) j
!!  if (sol%np==0) then
!!   write(*,*) trim(tbulk%fname)
!!   write(*,*) tbulk%itraj
!!   read(*,*)

!!  end if

! 
!   if (tbulk%change_file) then  
!    write(*,*) trim(tbulk%fname)
!   end if  
!  if (no_next) exit
! end do

! write(*,*) j

!end subroutine test_load


!subroutine test_calc
! 
!! call load_traj

! call test_load
! 
!end subroutine test_calc


end module traj_bulk_m



!program main
! use traj_bulk_m, only : test_calc

! call test_calc

!end program main


