!!-------------------------------------------------------
!!---- Crystallographic Fortran Modules Library (CrysFML)
!!-------------------------------------------------------
!!---- The CrysFML project is distributed under LGPL. In agreement with the
!!---- Intergovernmental Convention of the ILL, this software cannot be used
!!---- in military applications.
!!----
!!---- Copyright (C) 1999-2012  Institut Laue-Langevin (ILL), Grenoble, FRANCE
!!----                          Universidad de La Laguna (ULL), Tenerife, SPAIN
!!----                          Laboratoire Leon Brillouin(LLB), Saclay, FRANCE
!!----
!!---- Authors: Juan Rodriguez-Carvajal (ILL)
!!----          Javier Gonzalez-Platas  (ULL)
!!----
!!---- Contributors: Laurent Chapon     (ILL)
!!----               Marc Janoschek     (Los Alamos National Laboratory, USA)
!!----               Oksana Zaharko     (Paul Scherrer Institute, Switzerland)
!!----               Tierry Roisnel     (CDIFX,Rennes France)
!!----               Eric Pellegrini    (ILL)
!!----
!!---- This library is free software; you can redistribute it and/or
!!---- modify it under the terms of the GNU Lesser General Public
!!---- License as published by the Free Software Foundation; either
!!---- version 3.0 of the License, or (at your option) any later version.
!!----
!!---- This library is distributed in the hope that it will be useful,
!!---- but WITHOUT ANY WARRANTY; without even the implied warranty of
!!---- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!!---- Lesser General Public License for more details.
!!----
!!---- You should have received a copy of the GNU Lesser General Public
!!---- License along with this library; if not, see <http://www.gnu.org/licenses/>.
!!----
!!----
!!---- MODULE: CFML_IO_Formats
!!----   INFO: Creation/Conversion for several formats
!!----
!!---- HISTORY
!!----    Update: 07/03/2011
!!----
!!----
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!----    ERR_FORM
!!----    ERR_FORM_MESS
!!----    INTERVAL_TYPE
!!----    JOB_INFO_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       FILE_TO_FILELIST
!!----       GET_JOB_INFO
!!----       GET_PHASES_FILE
!!--++       GET_NPHASES_CIFFILE
!!--++       GET_NPHASES_PCRFILE
!!----       INIT_ERR_FORM
!!----       READ_ATOM
!!----       READ_CELL
!!----       READ_CIF_ATOM
!!----       READ_CIF_CELL
!!----       READ_CIF_CHEMICALNAME
!!----       READ_CIF_CONT
!!----       READ_CIF_HALL
!!----       READ_CIF_HM
!!----       READ_CIF_LAMBDA
!!----       READ_CIF_SYMM
!!----       READ_CIF_TITLE
!!----       READ_CIF_Z
!!----       READ_FILE_ATOM
!!--++       READ_FILE_ATOMLIST             [Overloaded]
!!--++       READ_FILE_POINTLIST            [Overloaded]
!!----       READ_FILE_CELL
!!--++       READ_FILE_CELLc                [Overloaded]
!!--++       READ_FILE_CELLt                [Overloaded]
!!----       READ_FILE_LAMBDA
!!----       READ_FILE_RNGSINTL
!!----       READ_FILE_SPG
!!----       READ_FILE_TRANSF
!!----       READ_SHX_ATOM
!!----       READ_SHX_CELL
!!----       READ_SHX_CONT
!!----       READ_SHX_FVAR
!!----       READ_SHX_LATT
!!----       READ_SHX_SYMM
!!----       READ_SHX_TITL
!!----       READ_UVALS
!!--++       READN_SET_XTAL_CFL             [Private]
!!--++       READN_SET_XTAL_CFL_MOLEC       [Private]
!!--++       READN_SET_XTAL_CIF             [Private]
!!--++       READN_SET_XTAL_PCR             [Private]
!!--++       READN_SET_XTAL_SHX             [Private]
!!----       READN_SET_XTAL_STRUCTURE
!!--++       READN_SET_XTAL_STRUCTURE_MOLCR [Overloaded]
!!--++       READN_SET_XTAL_STRUCTURE_SPLIT [Overloaded]
!!----       WRITE_CIF_POWDER_PROFILE
!!----       WRITE_CIF_TEMPLATE
!!----       WRITE_SHX_TEMPLATE
!!----
!!
 Module CFML_IO_Formats

    !---- Use modules ----!
    Use CFML_GlobalDeps,                only: cp,sp,pi,eps
    Use CFML_Math_General,              only: sind
    Use CFML_String_Utilities
    Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Set_Crystal_Cell, Convert_U_Betas, &
                                              Convert_B_Betas, U_Equiv, Convert_Betas_U
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, Set_SpaceGroup, Get_Multip_Pos
    Use CFML_Atom_TypeDef,              only: Atom_Type, Init_Atom_Type,atom_list_type,         &
                                              Allocate_atom_list, Deallocate_atom_list
    Use CFML_Molecular_Crystals,        only: Err_Molec, Err_Molec_Mess,Molecular_Crystal_Type, &
                                              Read_Molecule, Set_Euler_Matrix, Write_Molecule
    Use CFML_Geometry_Calc,             only: Point_List_Type, Get_Euler_from_Fract

    !---- Variables ----!
    implicit none

    private


    !---- List of public functions ----!

    !---- List of public subroutines ----!
    public :: Init_Err_Form, Read_Atom, Read_Cell, Read_Cif_Atom, Read_Cif_Cell,                 &
              Read_Cif_Cont, Read_Cif_Hall, Read_Cif_Hm, Read_Cif_Lambda, Read_Cif_Symm,         &
              Read_Cif_Title, Read_Cif_Z, Read_File_Atom, Read_File_Spg, Read_Cif_ChemicalName,  &
              Read_File_Transf, Read_Shx_Atom, Read_Shx_Cell, Read_Shx_Cont, Read_Shx_Fvar,      &
              Read_Shx_Latt, Read_Shx_Symm, Read_Shx_Titl, Read_Uvals, Write_Cif_Powder_Profile, &
              Write_Cif_Template, Write_Shx_Template, Read_File_rngSINTL, Read_File_Lambda,      &
              Get_job_info, File_To_FileList, Get_Phases_File

    !---- List of public overloaded procedures: subroutines ----!
    public :: Read_File_Cell, Readn_Set_Xtal_Structure, Write_Atoms_CFL, Write_CFL

    !---- List of private functions ----!

    !---- List of private subroutines ----!
    private:: Read_File_Cellc, Read_File_Cellt, Read_File_Atomlist,Read_File_Pointlist,               &
              Readn_Set_Xtal_CFL, Readn_Set_Xtal_CIF, Readn_Set_Xtal_PCR,Readn_Set_Xtal_SHX,          &
              Readn_Set_Xtal_CFL_Molec, Readn_Set_Xtal_Structure_Split,                               &
              Readn_Set_Xtal_Structure_Molcr, Get_NPhases_CIFFile,Get_NPHases_PCRFile,                &
              Write_CFL_Molcrys, Write_CFL_Atom_List_Type, Write_Atoms_CFL_ATM, Write_Atoms_CFL_MOLX, &
              Write_Atoms_CFL_MOLX_orig

    !---- Definitions ----!


    !!----
    !!---- ERR_FORM
    !!----    logical, public :: err_form
    !!----
    !!----    Logical Variable indicating an error in CFML_IO_Formats
    !!----
    !!---- Update: February - 2005
    !!
    logical, public :: err_form

    !!----
    !!---- ERR_FORM_MESS
    !!----    character(len=150), public :: ERR_Form_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150),       public  :: ERR_Form_Mess

    !!----
    !!---- EPSS
    !!----    real(kind=cp), parameter, private :: epss=1.0e-5_cp
    !!----
    !!----    Private small real number for floating point comparisons
    !!----
    !!---- Update: February - 2011
    !!
    real(kind=cp), parameter, private :: epss=1.0e-5_cp

    !!----
    !!---- TYPE :: INTERVAL_TYPE
    !!--..
    !!---- Type, public :: interval_type
    !!----    real(kind=cp) :: mina  !low limit
    !!----    real(kind=cp) :: maxb  !high limit
    !!---- End Type interval_type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: interval_type
       real(kind=cp) :: mina  !low limit
       real(kind=cp) :: maxb  !high limit
    End Type interval_type

    !!----
    !!---- TYPE :: JOB_INFO_TYPE
    !!--..
    !!---- Type, public :: Job_Info_type
    !!----    character(len=120)                            :: Title          ! Title
    !!----    integer                                       :: Num_Phases     ! Number of phases
    !!----    integer                                       :: Num_Patterns   ! Number of patterns
    !!----    integer                                       :: Num_cmd        ! Number of command lines
    !!----    character(len=16),  dimension(:), allocatable :: Patt_typ       ! Type of Pattern
    !!----    character(len=128), dimension(:), allocatable :: Phas_nam       ! Name of phases
    !!----    character(len=128), dimension(:), allocatable :: cmd            ! Command lines: text for actions
    !!----    type(interval_type),dimension(:), allocatable :: range_stl      ! Range in sinTheta/Lambda
    !!----    type(interval_type),dimension(:), allocatable :: range_q        ! Range in 4pi*sinTheta/Lambda
    !!----    type(interval_type),dimension(:), allocatable :: range_d        ! Range in d-spacing
    !!----    type(interval_type),dimension(:), allocatable :: range_2theta   ! Range in 2theta-spacing
    !!----    type(interval_type),dimension(:), allocatable :: range_Energy   ! Range in Energy
    !!----    type(interval_type),dimension(:), allocatable :: range_tof      ! Range in Time of Flight
    !!----    type(interval_type),dimension(:), allocatable :: Lambda         ! Lambda
    !!----    real(kind=cp)      ,dimension(:), allocatable :: ratio          ! ratio lambda2/lambda1
    !!----    real(kind=cp)      ,dimension(:), allocatable :: dtt1,dtt2      ! d-to-TOF coefficients
    !!---- End Type Job_Info_type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Job_Info_type
       character(len=120)                            :: Title
       integer                                       :: Num_Phases
       integer                                       :: Num_Patterns
       integer                                       :: Num_cmd
       character(len=16),  dimension(:), allocatable :: Patt_typ
       character(len=128), dimension(:), allocatable :: Phas_nam
       character(len=128), dimension(:), allocatable :: cmd
       type(interval_type),dimension(:), allocatable :: range_stl
       type(interval_type),dimension(:), allocatable :: range_q
       type(interval_type),dimension(:), allocatable :: range_d
       type(interval_type),dimension(:), allocatable :: range_2theta
       type(interval_type),dimension(:), allocatable :: range_Energy
       type(interval_type),dimension(:), allocatable :: range_tof
       type(interval_type),dimension(:), allocatable :: Lambda
       real(kind=cp)      ,dimension(:), allocatable :: ratio
       real(kind=cp)      ,dimension(:), allocatable :: dtt1,dtt2
    End Type Job_Info_type

    !!----
    !!---- TYPE :: FILE_LIST_TYPE
    !!--..
    !!---- Type,public :: File_List_Type
    !!----    integer                                       :: nlines ! Number of lines in the file
    !!----    character(len=256), allocatable, dimension(:) :: line   ! Content of the lines
    !!---- End Type file_list_type
    !!----
    !!---- Updated: February - 2005, November 2012
    !!
    Type,public :: File_List_Type
       integer                                       :: nlines
       character(len=256), allocatable, dimension(:) :: line
    End Type File_List_Type


    !---- Interfaces - Overloaded procedures--!
    Interface  Read_File_Cell
       Module Procedure Read_File_Cellc  !Last Output Argument Vector Of Six Component With The Cell Parameters
       Module Procedure Read_File_Cellt  !Last output argument object of type Crystal_cell_type
    End interface

    Interface Read_File_Atom
       Module Procedure Read_File_Atomlist   !Last Output Argument of type Atom_list_type
       Module Procedure Read_File_Pointlist  !Last output argument of type Point_list_type
    End Interface

    Interface Readn_Set_Xtal_Structure
       Module Procedure Readn_Set_Xtal_Structure_Molcr ! For Molecular Crystal Type
       Module Procedure Readn_Set_Xtal_Structure_Split ! For Cell, Spg, A types
    End Interface

    Interface Write_CFL
       Module Procedure Write_CFL_Molcrys        ! For Molecular Crystal Type
       Module Procedure Write_CFL_Atom_List_Type ! For Cell, Spg, A Types
    End Interface

    Interface Write_Atoms_CFL
       Module Procedure Write_Atoms_CFL_MOLX ! For Molecular Crystal Type
       Module Procedure Write_Atoms_CFL_ATM  ! For Cell, Spg, A Types
    End Interface

 Contains

    !---- Functions ----!

    !---- Subroutines ----!

    !!----
    !!---- Subroutine File_To_FileList(File_dat,File_list)
    !!----   character(len=*),     intent( in) :: file_dat  !Input data file
    !!----   type(file_list_type), intent(out) :: file_list !File list structure
    !!----
    !!----    Charge an external file to an object of File_List_Type.
    !!----
    !!---- Update: August - 2008
    !!
    Subroutine File_To_FileList(File_dat,File_list)
       !---- Arguments ----!
       character(len=*),      intent( in) :: file_dat
       type(file_list_type),  intent(out) :: file_list

       !---- Local Variables ----!
       integer                           :: nlines

       !---- Number of Lines in the input file ----!
       call Number_Lines(trim(File_dat), nlines)

       if (nlines==0) then
          err_form=.true.
          ERR_Form_Mess="The file "//trim(File_dat)//" contains nothing"
          return
       else
          file_list%nlines=nlines
          if (allocated(file_list%line)) deallocate(file_list%line)
          allocate(file_list%line(nlines))
          call reading_Lines(trim(File_dat),nlines,file_list%line)
       end if

       return
    End Subroutine File_To_FileList

    !!----
    !!---- Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
    !!----   character(len=*), dimension(:), intent( in) :: file_dat     !Lines of text (content of a file)
    !!----   integer,                        intent( in) :: i_ini,i_end  !Lines to explore
    !!----   type(job_info_type),            intent(out) :: Job_info     !Object to be constructed here
    !!----
    !!----
    !!----    Constructor of the object Job_info. The arrary of strings file_dat
    !!----    have to be provided as input. It contains lines corresponding to the
    !!----    input control file. The analysis of the command lines is not given here.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Job_Info(file_dat,i_ini,i_end,Job_info)
       !---- Arguments ----!
       character(len=*), dimension(:), intent( in) :: file_dat
       integer,                        intent( in) :: i_ini,i_end
       type(job_info_type),            intent(out) :: Job_info

       !---- Local Variables ----!
       integer                           :: i,nphas, ncmd,n_pat,ier, j
       integer, dimension(i_end-i_ini+1) :: ip,ic,ipt
       real(kind=sp)                     :: a1,a2,a3,a4,a5
       character(len=120)                :: line, fmtfields, fmtformat

       !--- Initialize FindFMT
       call Init_FindFMT(i_ini)
       nphas=0
       ncmd=0
       n_pat=0
       ip=i_end
       ic=0
       ipt=0
       Job_info%title=" General Job: CrysFML"
       Job_info%Num_Patterns=1

       do i=i_ini,i_end
          line=u_case(adjustl(file_dat(i)))
          if (line(1:5) == "TITLE") Job_info%title=line(7:)
          if (line(1:5) == "NPATT") then
             read(unit=line(7:), fmt=*,iostat=ier) Job_info%Num_Patterns
             if (ier /= 0) Job_info%Num_Patterns=1
          end if
          if (line(1:6) == "PHASE_") then
             nphas=nphas+1
             ip(nphas)=i
          end if
          if (line(1:4) == "CMDL") then
             ncmd=ncmd+1
             ic(ncmd)=i
          end if
          if (line(1:5) == "PATT_") then
             n_pat=n_pat+1
             ipt(n_pat)=i
          end if
       end do

       if (nphas == 0) then
          nphas=1
          ip(nphas)=0
       end if
       if (n_pat == 0) then
          n_pat=1
          ipt(n_pat) = 0
       end if

       if (Job_info%Num_Patterns /= n_pat) Job_info%Num_Patterns = n_pat
       Job_info%Num_Phases=nphas
       Job_info%Num_Cmd=ncmd

       if (allocated(Job_Info%Patt_typ)) deallocate(Job_Info%Patt_typ)
       allocate(Job_Info%Patt_typ(n_pat))

       if (allocated(Job_Info%Phas_nam)) deallocate(Job_Info%Phas_nam)
       allocate(Job_Info%Phas_nam(nphas))

       if (allocated(Job_Info%range_stl)) deallocate(Job_Info%range_stl)
       allocate(Job_Info%range_stl(n_pat))

       if (allocated(Job_Info%range_q)) deallocate(Job_Info%range_q)
       allocate(Job_Info%range_q(n_pat))

       if (allocated(Job_Info%range_d)) deallocate(Job_Info%range_d)
       allocate(Job_Info%range_d(n_pat))

       if (allocated(Job_Info%range_2theta)) deallocate(Job_Info%range_2theta)
       allocate(Job_Info%range_2theta(n_pat))

       if (allocated(Job_Info%range_energy)) deallocate(Job_Info%range_energy)
       allocate(Job_Info%range_energy(n_pat))

       if (allocated(Job_Info%range_tof)) deallocate(Job_Info%range_tof)
       allocate(Job_Info%range_tof(n_pat))

       if (allocated(Job_Info%lambda)) deallocate(Job_Info%lambda)
       allocate(Job_Info%lambda(n_pat))

       if (allocated(Job_Info%ratio)) deallocate(Job_Info%ratio)
       allocate(Job_Info%ratio(n_pat))

       if (allocated(Job_Info%dtt1)) deallocate(Job_Info%dtt1)
       allocate(Job_Info%dtt1(n_pat))

       if (allocated(Job_Info%dtt2)) deallocate(Job_Info%dtt2)
       allocate(Job_Info%dtt2(n_pat))

       !---- Initialize all variables
       Job_Info%Patt_typ    =" "
       Job_Info%Phas_nam    =" "
       Job_Info%range_stl%mina=0.0
       Job_Info%range_stl%maxb=0.0
       Job_Info%range_q%mina=0.0
       Job_Info%range_q%maxb=0.0
       Job_Info%range_d%mina=0.0
       Job_Info%range_d%maxb=0.0
       Job_Info%range_2theta%mina=0.0
       Job_Info%range_2theta%maxb=0.0
       Job_Info%range_Energy%mina=0.0
       Job_Info%range_Energy%maxb=0.0
       Job_Info%range_tof%mina=0.0
       Job_Info%range_tof%maxb=0.0
       Job_Info%Lambda%mina=0.0
       Job_Info%Lambda%maxb=0.0
       Job_Info%ratio = 0.0
       Job_Info%dtt1 = 0.0
       Job_Info%dtt2 = 0.0
       if (ncmd > 0) then
          if (allocated(Job_Info%cmd)) deallocate(Job_Info%cmd)
          allocate(Job_Info%cmd(ncmd))
          Job_Info%cmd=" "
       end if

       !---- Fill the different fields of Job_Info
       !---- Start with patterns
       fmtfields = "9fffff"

       !---- First asks if there is a PATT_ card, if not a standard is taken
       if (ipt(1) /= 0) then
          do n_pat=1, Job_info%Num_Patterns
             i=ipt(n_pat)
             line=u_case(adjustl(file_dat(i)))
             line=line(8:)
             call findfmt(0,line,fmtfields,fmtformat)
             read(unit=line,fmt=fmtformat) Job_Info%Patt_typ(n_pat), a1,a2,a3,a4,a5
             if (ierr_fmt /= 0) return
             line=u_case(Job_Info%Patt_typ(n_pat))

             select case(line(1:9))
                case("XRAY_2THE","NEUT_2THE","XRAY_SXTA","NEUT_SXTA")
                   if ( a1 <= 0.000001) a1=1.5405
                   if ( a2 <= 0.000001) then
                      a2=a1
                      a3=0.0
                   end if
                   if (a5 <= a4) a5=120.0
                   Job_Info%Lambda(n_pat)%mina=a1
                   Job_Info%Lambda(n_pat)%maxb=a2
                   Job_Info%ratio(n_pat)=a3
                   Job_Info%range_2theta(n_pat)%mina=a4
                   Job_Info%range_2theta(n_pat)%maxb=a5
                   a4=sind(0.5*a4)/a1
                   a5=sind(0.5*a5)/a2
                   Job_Info%range_stl(n_pat)%mina=a4
                   Job_Info%range_stl(n_pat)%maxb=a5
                   Job_Info%range_q(n_pat)%mina=a4*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
                   Job_Info%range_d(n_pat)%mina=0.5/a5
                   Job_Info%range_d(n_pat)%maxb=0.5/a4

                case("NEUT_TOF ")
                   if (a1 <= 0.000001) a1=1000.0
                   if (a4 <= a3) a4=2.0*abs(a3)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=a2
                   Job_Info%range_tof(n_pat)%mina=a3
                   Job_Info%range_tof(n_pat)%maxb=a4
                   Job_Info%range_d(n_pat)%mina=0.5*(-1.0+sqrt(1.0+4.0*a2*a3/a1/a1))
                   Job_Info%range_d(n_pat)%maxb=0.5*(-1.0+sqrt(1.0+4.0*a2*a4/a1/a1))
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

                case("XRAY_ENER")
                   if (a1 <= 0.000001) a1=12.4 !(=hc(keV.Angstr.)
                   Job_Info%dtt1(n_pat)=a1
                   Job_Info%dtt2(n_pat)=0.0
                   Job_Info%range_energy(n_pat)%mina=a3
                   Job_Info%range_energy(n_pat)%maxb=a4
                   if (a3 <= 0.00001) a3=0.01
                   if (a4 <= 0.00001) a4=2.00
                   Job_Info%range_d(n_pat)%mina=a1/a4
                   Job_Info%range_d(n_pat)%maxb=a1/a3
                   Job_Info%range_stl(n_pat)%mina=0.5/Job_Info%range_d(n_pat)%maxb
                   Job_Info%range_stl(n_pat)%maxb=0.5/Job_Info%range_d(n_pat)%mina
                   Job_Info%range_q(n_pat)%mina=Job_Info%range_stl(n_pat)%mina*4.0*pi
                   Job_Info%range_q(n_pat)%maxb=Job_Info%range_stl(n_pat)%maxb*4.0*pi

             end select
          end do

       else
          n_pat=1
          a1=1.5405
          a2=a1
          a3=0.0
          a4=0.0
          a5=120.0
          Job_Info%Patt_typ(n_pat)="XRAY_2THE"
          Job_Info%Lambda(n_pat)%mina=a1
          Job_Info%Lambda(n_pat)%maxb=a2
          Job_Info%ratio(n_pat)=a3
          Job_Info%range_2theta(n_pat)%mina=a4
          Job_Info%range_2theta(n_pat)%maxb=a5
          a4=sind(0.5*a4)/a1
          a5=sind(0.5*a5)/a2
          Job_Info%range_stl(n_pat)%mina=a4
          Job_Info%range_stl(n_pat)%maxb=a5
          Job_Info%range_q(n_pat)%mina=a4*4.0*pi
          Job_Info%range_q(n_pat)%maxb=a5*4.0*pi
          Job_Info%range_d(n_pat)%mina=0.5/a5
          Job_Info%range_d(n_pat)%maxb=0.5/a4
       end if

       !---- Phase names
       if (ip(1) /= 0) then
          do i=1,nphas
             j=ip(i)
             line=adjustl(file_dat(j))
             Job_Info%Phas_nam(i)=line(8:)
          end do
       else
          Job_Info%Phas_nam(1)= Job_info%title
       end if

       !---- Command Lines, stored but not analysed here
       do i=1,ncmd
          j=ic(i)
          line=adjustl(file_dat(j))
          Job_Info%cmd(i)=line(8:)
       end do

       return
    End Subroutine Get_Job_Info

    !!----
    !!---- Subroutine Init_Err_Form()
    !!----
    !!----    Initialize Errors Variable for this module
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Form()

       err_form=.false.
       ERR_Form_Mess=" "

       return
    End Subroutine Init_Err_Form

    !!----
    !!---- Subroutine Read_Atom(Line,Atomo)
    !!----    character(len=*), intent(in out ) :: line    !  In -> Input String with ATOM directive
    !!----    Type (Atom_Type), intent(out)     :: Atomo   ! Out -> Parameters on variable Atomo
    !!----
    !!----    Subroutine to read the atom parameters from a given "line"
    !!----    it construct the object Atomo of type Atom.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Atom(line,Atomo)
       !---- Arguments ----!
       character(len=*), intent(in out ) :: line
       Type (Atom_Type), intent(out)     :: Atomo

       !---- Local variables -----!
       integer                           :: iv, nlong1,n,ier,q
       real(kind=cp), dimension (10)     :: vet1
       real(kind=cp), dimension (10)     :: vet2
       character(len=4)                  :: dire
       character(len=5)                  :: label
       character(len=132), dimension(1)  :: filevar
       character(len=*), parameter       :: digpm="0123456789+-"

       !---- Init ----!
       call init_err_form()
       call init_atom_type(Atomo)
       q=0
       iv=index(line,"#")
       if(iv /= 0) atomo%AtmInfo=line(iv+1:)

       call cutst(line,nlong1,dire)
       if (u_case(dire) /= "ATOM") then
          err_form=.true.
          ERR_Form_Mess=" Error reading the ATOM keyword"
          return
       end if

       !---- Atom Label ----!
       call cutst(line,nlong1,label)
       atomo%lab=label(1:5)

       !---- Atom Type (Chemical symbol & Scattering Factor) ----!
       call cutst(line,nlong1,label)
       n=index(digpm,label(2:2))
       if (n /=0) then
         atomo%chemsymb=U_case(label(1:1))
       else
         atomo%chemsymb=U_case(label(1:1))//L_case(label(2:2))
       end if
       atomo%SfacSymb=label(1:4)

       !---- Parameters ----!
       filevar(1)="atm "//trim(line)

       n=1
       call Read_Key_ValueSTD(filevar,n,n,"atm",vet1,vet2,iv)
      ! call getnum(line,vet,ivet,iv)
       if (iv <= 0) then
          err_form=.true.
          ERR_Form_Mess= "Error reading parameters of atom:"//atomo%lab
          return
       end if

       !---- Coordinates  ----!
       if (iv < 3) then
          err_form=.true.
          ERR_Form_Mess= "Error reading Coordinates of atom:"//atomo%lab
          return
       end if

       atomo%x(:)=vet1(1:3)
       atomo%x_std(:)=vet2(1:3)

       !---- Biso ----!
       if (iv > 3) then
         atomo%biso=vet1(4)
         atomo%biso_std=vet2(4)
       end if

       !---- Occ ----!
       if (iv > 4) then
          atomo%occ=vet1(5)
          atomo%occ_std=vet2(5)
       end if

       !---- Moment ----!
       if (iv > 5) atomo%moment=vet1(6)

       !---- Charge ----!
       if (iv > 6) atomo%charge=vet1(7)

       !Attempt to get the oxidation state from "Label"
       if(abs(atomo%charge) < eps) then
         iv=index(label,"+")
         Select Case(iv)
           Case(0) !No + sign
             n=index(label,"-")
             Select Case(n)
               Case(2) !Element with a single character symbol F-1
                  read(unit=label(3:),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
               Case(3) !Element in the form: F1- or Br-1
                  read(unit=label(2:2),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) then
                        read(unit=label(4:4),fmt="(i1)",iostat=ier)  q
                        if (ier /= 0) q=0
                  end if
               Case(4) !Element in the form: Br1-
                  read(unit=label(3:3),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
             End Select
             q=-q   !anions
           Case(2) !Element with a single character symbol C+4
                  read(unit=label(3:),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
           Case(3) !Element in the form: C4+ or Fe+3
                  read(unit=label(2:2),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) then
                        read(unit=label(4:4),fmt="(i1)",iostat=ier)  q
                        if (ier /= 0) q=0
                  end if
           Case(4) !Element in the form: Fe3+
                  read(unit=label(3:3),fmt="(i1)",iostat=ier)  q
                  if (ier /= 0) q=0
         End Select
         atomo%charge=real(q)
       end if
       return
    End Subroutine Read_Atom

    !!----
    !!---- Subroutine Read_Cell(Line,Celda)
    !!----    character(len=*),          intent(in out ) :: line   !  In -> Input String with CELL Directive
    !!----    real(kind=cp),dimension(6),intent(out)     :: Celda  !  In -> Parameters on Celda Variable
    !!----
    !!----    Subroutine to read the cell parameters from a given "line"
    !!----    it construct the object Celda of type Crystal_Cell.
    !!----    Assumes the string "line" has been read from a file and
    !!----    starts with the word "cell", that is removed before reading
    !!----    the values of the parameters.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cell(line,Celda)
       !---- Arguments ----!
       character(len=*),          intent(in out ) :: line
       real(kind=cp),dimension(6),intent(out)     :: Celda

       !---- Local variables -----!
       integer, dimension (6)               :: ivet
       real(kind=cp), dimension (6)         :: vet
       integer                              :: nlong1,iv
       character(len=4)                     :: dire

       call init_err_form()

       call cutst(line,nlong1,dire)
       if (u_case(dire) /= "CELL") then
          err_form=.true.
          ERR_Form_Mess=" Error reading the CELL keyword"
          return
       end if

       call getnum(line,vet,ivet,iv)
       if (iv /= 6 ) then
          err_form=.true.
          ERR_Form_Mess=" Error reading the Cell Parameters"
          return
       else
          celda=vet
       end if

       return
    End Subroutine Read_Cell

    !!----
    !!---- Subroutine Read_Cif_Atom(Filevar,Nline_Ini,Nline_End,N_Atom,Atm_List)
    !!----    character(len=*),dimension(:), intent(in)     :: filevar    !  In -> Input strings information
    !!----    integer,                       intent(in out) :: nline_ini  !  In -> Line to beginning search
    !!----                                                                   Out -> Current line on Filevar
    !!----    integer,                       intent(in)     :: nline_end  !  In -> Line to the End search
    !!----    integer,                       intent(out)    :: n_atom     ! Out -> Actual number of atom
    !!----    type (atom_list_type),        intent(out)    :: Atm_List   ! Out -> Atom list
    !!----
    !!----    Obtaining Atoms parameters from Cif file. A control error is present.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Atom(filevar,nline_ini,nline_end,n_atom,Atm_List)
       !---- Arguments ----!
       character(len=*), dimension(:),   intent(in)      :: filevar
       integer,                          intent(in out)  :: nline_ini
       integer,                          intent(in)      :: nline_end
       integer,                          intent(out)     :: n_atom
       type (atom_list_type),            intent(out)     :: Atm_List

       !---- Local Variables ----!
       character(len=len(filevar(1)))               :: string,cp_str
       character(len=20),dimension(15)              :: label

       integer                         :: i, j, nc, nct, nline, iv
       !integer, dimension(1)           :: ivet
       integer, dimension( 7)          :: lugar   !   1 -> label
                                                  !   2 -> Symbol
                                                  ! 3-5 -> coordinates
                                                  !   6 -> occupancy
                                                  !   7 -> Uequi
       real(kind=cp), dimension(1)     :: vet1,vet2
       type(atom_list_type)            :: Atm

       !---- Estimacion Inicial ----!
       lugar=0
       call allocate_atom_list(nline_end-nline_ini+1,Atm)

       n_atom=0
       call Read_Key_StrVal(filevar,nline_ini,nline_end,"_atom_site_",string)

       j=0
       do i=nline_ini,nline_end
          string=adjustl(filevar(i))
          if ("_atom_site_label" == string(1:16)) then
             j=j+1
             lugar(1)=j
             cycle
          end if
          if ("_atom_site_type_symbol" == string(1:22)) then
             j=j+1
             lugar(2)=j
             cycle
          end if
          if ("_atom_site_fract_x" == string(1:18)) then
             j=j+1
             lugar(3)=j
             cycle
          end if
          if ("_atom_site_fract_y" == string(1:18)) then
             j=j+1
             lugar(4)=j
             cycle
          end if
          if ("_atom_site_fract_z" == string(1:18)) then
             j=j+1
             lugar(5)=j
             cycle
          end if
          if ("_atom_site_occupancy" == string(1:20)) then
             j=j+1
             lugar(6)=j
             cycle
          end if
          if ("_atom_site_U_iso_or_equiv" == string(1:25)) then
             j=j+1
             lugar(7)=j
             cycle
          end if
          if ("_atom_site_" == string(1:11)) then
             j=j+1
             cycle
          end if

          if ("_oxford_atom_site_" == string(1:18)) then
             j=j+1
             cycle
          end if

          nline=i
          exit
       end do

       if (any(lugar(3:5) == 0)) then
          err_form=.true.
          ERR_Form_Mess=" Error reading atoms"
          return
       end if
       nct=count(lugar > 0)

       nline_ini=nline
       string=" "
       do i=nline_ini,nline_end
          string=adjustl(trim(string)//" "//filevar(i))
          if (string(1:1) == "#" .or. string(1:1) == "?") cycle
          if (len_trim(string) == 0) exit
          if (string(1:1) == "_" .or. string(1:5) == "loop_") exit
          cp_str=string
          call getword(cp_str,label,nc)
          if (nc < nct) cycle

          n_atom=n_atom+1

          ! _atom_site_label
          atm%atom(n_atom)%lab=label(lugar(1))

          ! _atom_site_type_symbol
          if (lugar(2) /= 0) then
             atm%atom(n_atom)%SfacSymb=label(lugar(2))(1:4)
             if(index("1234567890+-",label(lugar(2))(2:2)) /= 0 ) then
                atm%atom(n_atom)%chemSymb=U_case(label(lugar(2))(1:1))
             else
                atm%atom(n_atom)%chemSymb=U_case(label(lugar(2))(1:1))//L_case(label(lugar(2))(2:2))
             end if
           !  call getnum(label(lugar(2))(2:2),vet1,ivet,iv)
           !  if (iv <= 0) then
           !     atm%atom(n_atom)%chemSymb=label(lugar(2))(1:2)
           !  else
           !     atm%atom(n_atom)%chemSymb=label(lugar(2))(1:1)
           !  end if
          else
             if(index("1234567890+-",label(lugar(1))(2:2)) /= 0 ) then
                atm%atom(n_atom)%chemSymb=U_case(label(lugar(1))(1:1))
             else
                atm%atom(n_atom)%chemSymb=U_case(label(lugar(1))(1:1))//L_case(label(lugar(1))(2:2))
             end if
             atm%atom(n_atom)%SfacSymb=atm%atom(n_atom)%chemSymb

           !   call getnum(label(lugar(1))(2:2),vet1,ivet,iv)
           !  if (iv <= 0) then
           !     atm%atom(n_atom)%chemSymb=label(lugar(1))(1:2)
           !  else
           !     atm%atom(n_atom)%chemSymb=label(lugar(1))(1:1)
           !  end if
          end if

          call getnum_std(label(lugar(3)),vet1,vet2,iv)    ! _atom_site_fract_x
          atm%atom(n_atom)%x(1)=vet1(1)
          atm%atom(n_atom)%x_std(1)=vet2(1)
          call getnum_std(label(lugar(4)),vet1,vet2,iv)    ! _atom_site_fract_y
          atm%atom(n_atom)%x(2)=vet1(1)
          atm%atom(n_atom)%x_std(2)=vet2(1)
          call getnum_std(label(lugar(5)),vet1,vet2,iv)    ! _atom_site_fract_z
          atm%atom(n_atom)%x(3)=vet1(1)
          atm%atom(n_atom)%x_std(3)=vet2(1)

          ! _atom_site_occupancy
          if (lugar(6) /= 0) then
             call getnum_std(label(lugar(6)),vet1,vet2,iv)
          else
             vet1=1.0
          end if
          atm%atom(n_atom)%occ=vet1(1)
          atm%atom(n_atom)%occ_std=vet2(1)

          if (lugar(7) /= 0) then
             call getnum_std(label(lugar(7)),vet1,vet2,iv)    ! _atom_site_Uiso_or_equiv
          else
             vet1=0.0
          end if
          atm%atom(n_atom)%ueq=vet1(1)
          atm%atom(n_atom)%Biso=vet1(1)*78.95683521     !If anisotropic they
          atm%atom(n_atom)%Biso_std=vet2(1)*78.95683521 !will be put to zero

          atm%atom(n_atom)%utype="u_ij"
          string=" "

       end do
       nline=i

       !---- Anisotropic parameters ----!
       nline_ini=nline
       lugar=0
       call Read_Key_StrVal(filevar,nline_ini,nline_end,"_atom_site_aniso_",string)

       j=0
       do i=nline_ini,nline_end
          string=adjustl(filevar(i))
          if ("_atom_site_aniso_label" == string(1:22)) then
             j=j+1
             lugar(1)=j
             cycle
          end if
          if ("_atom_site_aniso_U_11" == string(1:21)) then
             j=j+1
             lugar(2)=j
             cycle
          end if
          if ("_atom_site_aniso_U_22" == string(1:21)) then
             j=j+1
             lugar(3)=j
             cycle
          end if
          if ("_atom_site_aniso_U_33" == string(1:21)) then
             j=j+1
             lugar(4)=j
             cycle
          end if
          if ("_atom_site_aniso_U_12" == string(1:21)) then
             j=j+1
             lugar(5)=j
             cycle
          end if
          if ("_atom_site_aniso_U_13" == string(1:21)) then
             j=j+1
             lugar(6)=j
             cycle
          end if
          if ("_atom_site_aniso_U_23" == string(1:21)) then
             j=j+1
             lugar(7)=j
             cycle
          end if

          if ("_atom_site_aniso" == string(1:16) ) then
             j=j+1
             cycle
          end if

          nline=i
          exit
       end do

       if (all(lugar > 0)) then
          nct=count(lugar > 0)
          nline_ini=nline
          string=" "
          do i=nline_ini,nline_end
             string=adjustl(trim(string)//" "//filevar(i))
             if (string(1:1) == "#" .or. string(1:1) =="?") cycle
             if (len_trim(string) == 0) exit

             cp_str=string
             call getword(cp_str,label,nc)
             if (nc < nct) cycle

             do j=1,n_atom
                if (atm%atom(j)%lab(1:4) /= label(lugar(1))(1:4)) cycle

                call getnum_std(label(lugar(2)),vet1,vet2,iv)    ! _atom_site_aniso_U_11
                atm%atom(j)%u(1)=vet1(1)
                atm%atom(j)%u_std(1)=vet2(1)
                call getnum_std(label(lugar(3)),vet1,vet2,iv)    ! _atom_site_aniso_U_22
                atm%atom(j)%u(2)=vet1(1)
                atm%atom(j)%u_std(2)=vet2(1)
                call getnum_std(label(lugar(4)),vet1,vet2,iv)    ! _atom_site_aniso_U_33
                atm%atom(j)%u(3)=vet1(1)
                atm%atom(j)%u_std(3)=vet2(1)
                call getnum_std(label(lugar(5)),vet1,vet2,iv)    ! _atom_site_aniso_U_12
                atm%atom(j)%u(4)=vet1(1)
                atm%atom(j)%u_std(4)=vet2(1)
                call getnum_std(label(lugar(6)),vet1,vet2,iv)    ! _atom_site_aniso_U_13
                atm%atom(j)%u(5)=vet1(1)
                atm%atom(j)%u_std(5)=vet2(1)
                call getnum_std(label(lugar(7)),vet1,vet2,iv)    ! _atom_site_aniso_U_23
                atm%atom(j)%u(6)=vet1(1)
                atm%atom(j)%u_std(6)=vet2(1)

                atm%atom(j)%thtype="aniso"
                atm%atom(j)%Biso=0.0
                atm%atom(j)%Biso_std=0.0
                exit
             end do
             nline=i
             string=" "
          end do

       end if
       nline_ini=nline

       !---- Adjusting ... ----!
       if (n_atom > 0) then
          call allocate_atom_list(n_atom,Atm_list)
          atm_list%natoms=n_atom
          do i=1,n_atom
             atm_list%atom(i)=atm%atom(i)
          end do
       end if
       call Deallocate_atom_list(atm)

       return
    End Subroutine Read_Cif_Atom

    !!----
    !!---- Subroutine Read_Cif_Cell(Filevar,Nline_Ini,Nline_End,Celda,Stdcelda)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar      !  In -> String vector input
    !!----    integer,                        intent(in out) :: nline_ini    !  In -> Line to start the search
    !!----                                                                      Out -> Current line on Filevar
    !!----    integer,                        intent(in)     :: nline_end    !  In -> Line to finish the search
    !!----    real(kind=cp),dimension(6),     intent (out)   :: Celda        ! Out -> Cell variable
    !!----    real(kind=cp),dimension(6),     intent (out)   :: StdCelda     ! Out -> Cell variable
    !!----
    !!----    Read Cell Parameters from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Cell(Filevar,Nline_Ini,Nline_End,Celda,StdCelda)
       !---- Arguments ----!
       character(len=*),  dimension(:),     intent(in)     :: filevar
       integer,                             intent(in out) :: nline_ini
       integer,                             intent(in)     :: nline_end
       real(kind=cp),dimension(6),          intent(out)    :: Celda
       real(kind=cp),dimension(6),optional, intent(out)    :: StdCelda

       !---- Local Variables ----!
       integer                     :: iv,initl
       real(kind=cp), dimension(1) :: vet1,vet2
       real(kind=cp), dimension(6) :: a

       !---- Valores iniciales ----!
       celda=(/1.0,1.0,1.0,90.0,90.0,90.0/)
       a=0.0
       if (present(stdcelda)) stdcelda=0.0

       !---- Celda ----!
       initl=nline_ini  !Preserve initial line => some CIF files have random order for cell parameters
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_length_a",vet1,vet2,iv)
       if (iv == 1) then
          Celda(1)   =vet1(1)
          a(1)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_length_b",vet1,vet2,iv)
       if (iv == 1) then
          Celda(2)   =vet1(1)
          a(2)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_length_c",vet1,vet2,iv)
       if (iv == 1) then
          Celda(3)   =vet1(1)
         a(3)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_angle_alpha",vet1,vet2,iv)
       if (iv == 1) then
          Celda(4)   =vet1(1)
          a(4)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_angle_beta",vet1,vet2,iv)
       if (iv == 1) then
          Celda(5)   =vet1(1)
          a(5)=vet2(1)
       end if

       nline_ini=initl
       call read_key_valueSTD(filevar,nline_ini,nline_end,"_cell_angle_gamma",vet1,vet2,iv)
       if (iv == 1) then
          Celda(6)   =vet1(1)
          a(6)=vet2(1)
       end if
       if (present(stdcelda)) stdcelda=a

       return
    End Subroutine Read_Cif_Cell

    !!----
    !!---- Subroutine Read_Cif_ChemicalName(Filevar,Nline_Ini,Nline_End,ChemName)
    !!----    character(len=*),  dimension(:), intent(in) :: filevar      !  In -> String vector
    !!----    integer,           intent(in out)           :: nline_ini    !  In -> Line to start the search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,           intent(in)               :: nline_end    !  In -> Line to finish the search
    !!----    character(len=*),  intent(out)              :: ChemName     ! Out -> Title string
    !!----
    !!----    Obtaining Chemical Name from Cif file
    !!----
    !!---- Update: March - 2009
    !!
    Subroutine Read_Cif_ChemicalName(Filevar,Nline_Ini,Nline_End,ChemName)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       character(len=*),  intent(out)              :: ChemName

       !---- Local variables ----!
       integer :: np1, np2

       ChemName=" "
       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_chemical_name_common",ChemName)

       if (len_trim(chemname) == 0) then
          call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_chemical_name_systematic",ChemName)
       end if

       if (len_trim(chemname) > 0) then
          if (trim(chemname) =="; ?" .or. trim(chemname)=="#") chemname=" "
          np1=index(chemname,"'")
          np2=index(chemname,"'",back=.true.)
          if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
             chemname=chemname(np1+1:np2-1)
          else
             np1=index(chemname,'"')
             np2=index(chemname,'"',back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                chemname=chemname(np1+1:np2-1)
             end if
          end if
       end if

       return
    End Subroutine Read_Cif_ChemicalName

    !!----
    !!---- Subroutine Read_Cif_Cont(Filevar,Nline_Ini,Nline_End,N_Elem_Type,Elem_Type,N_Elem)
    !!----    character(len=*), dimension(:),      intent(in)      :: filevar       !  In -> String vector input
    !!----    integer,                             intent(in out)  :: nline_ini     !  In -> Line to start the search
    !!----                                                                             Out -> Actual line on Filevar
    !!----    integer,                             intent(in)      :: nline_end     !  In -> Line to finish the search
    !!----    integer,                             intent(out)     :: n_elem_type   ! Out -> N. of different elements
    !!----    character(len=*), dimension(:),      intent(out)     :: elem_type     ! Out -> String for Element type
    !!----    real(kind=cp), dimension(:),optional,intent(out)     :: n_elem        ! Out -> Number of elements
    !!----
    !!----    Obtaining the chemical contents from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Cont(Filevar,Nline_Ini,Nline_End,N_Elem_Type,Elem_Type,N_Elem)
       !---- Arguments ----!
       character(len=*), dimension(:),      intent(in)      :: filevar
       integer,                             intent(in out)  :: nline_ini
       integer,                             intent(in)      :: nline_end
       integer,                             intent(out)     :: n_elem_type
       character(len=*), dimension(:),      intent(out)     :: elem_type
       real(kind=cp), dimension(:),optional,intent(out)     :: n_elem

       !---- Local  variables ----!
       character(len=len(filevar(1)))      :: string
       character(len=10),dimension(15)     :: label

       integer                    :: iv
       integer                    :: i,np1,np2,nlabel,nlong
       integer, dimension(1)      :: ivet

       real(kind=cp),dimension(1) :: vet

       n_elem_type = 0
       elem_type   = " "
       if (present(n_elem)) n_elem = 0.0

       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_chemical_formula_sum",string)
       if (len_trim(string) ==0) string=filevar(nline_ini+1)
       string=adjustl(string)
       if (string(1:1) == "?") return
       np1=index(string,"'")
       np2=index(string,"'",back=.true.)
       nlabel=0
       if (np1 /= 0 .and. np2 /= 0 .and. np2 > np1) then
          call getword(string(np1+1:np2-1),label,nlabel)
       end if
       if (nlabel /=0) then
          n_elem_type = nlabel
          do i=1,nlabel
             nlong=len_trim(label(i))
             select case (nlong)
                 case (1)
                    elem_type(i)=label(i)(1:1)
                    if (present(n_elem)) n_elem(i)   = 1.0

                 case (2)
                    call getnum(label(i)(2:),vet,ivet,iv)
                    if (iv == 1) then
                       elem_type(i)=label(i)(1:1)
                       if (present(n_elem)) n_elem(i)   =vet(1)
                    else
                       elem_type(i)=label(i)(1:2)
                       if (present(n_elem)) n_elem(i)   = 1.0
                    end if

                 case (3:)
                    call getnum(label(i)(2:),vet,ivet,iv)
                    if (iv == 1) then
                       elem_type(i)=label(i)(1:1)
                       if (present(n_elem)) n_elem(i)   =vet(1)
                    else
                       call getnum(label(i)(3:),vet,ivet,iv)
                       if (iv == 1) then
                          elem_type(i)=label(i)(1:2)
                          if (present(n_elem)) n_elem(i)   =vet(1)
                       else
                          elem_type(i)=label(i)(1:2)
                          if (present(n_elem)) n_elem(i)   = 1.0
                       end if

                    end if

             end select
          end do
       end if

       return
    End Subroutine Read_Cif_Cont

    !!----
    !!---- Subroutine Read_Cif_Hall(Filevar,Nline_Ini,Nline_End,Spgr_Ha)
    !!----    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String vector input
    !!----    integer,          intent(in out)           :: nline_ini    !  In -> Line to start the search
    !!----                                                                 Out -> Actual line on Filevar
    !!----    integer,          intent(in)               :: nline_end    !  In -> Line to finish the search
    !!----    character(len=*), intent(out)              :: spgr_ha      ! Out -> Hall symbol
    !!----
    !!----    Obtaining the Hall symbol of the Space Group
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Hall(Filevar,Nline_Ini,Nline_End,Spgr_Ha)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: filevar
       integer,          intent(in out)           :: nline_ini
       integer,          intent(in)               :: nline_end
       character(len=*), intent(out)              :: spgr_ha

       !---- Local variables ----!
       integer :: np1, np2

       spgr_ha=" "
       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_symmetry_space_group_name_Hall",spgr_ha)
       if (len_trim(spgr_ha)==0) spgr_ha=adjustl(filevar(nline_ini+1))
       ! TR  feb. 2015 .(re-reading the same item with another name)
       if(len_trim(spgr_ha) == 0) then
        spgr_ha=" "
        call Read_Key_StrVal(filevar,nline_ini,nline_end, "_space_group_name_Hall",spgr_ha)
        if (len_trim(spgr_ha)==0) spgr_ha=adjustl(filevar(nline_ini+1))
       end if

       if (spgr_ha =="?" .or. spgr_ha=="#") then
          spgr_ha=" "
       else
          np1=index(spgr_ha,"'")
          np2=index(spgr_ha,"'",back=.true.)
          if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
             spgr_ha=spgr_ha(np1+1:np2-1)
          else
             np1=index(spgr_ha,'"')
             np2=index(spgr_ha,'"',back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                spgr_ha=spgr_ha(np1+1:np2-1)
             else
                spgr_ha=" "
             end if
          end if
       end if

       return
    End Subroutine Read_Cif_Hall

    !!----
    !!---- Subroutine Read_Cif_Hm(Filevar,Nline_Ini,Nline_End,Spgr_Hm)
    !!----    character(len=*),  dimension(:), intent(in) :: filevar     !  In -> String vector
    !!----    integer,           intent(in out)           :: nline_ini   !  In -> Line to start the search
    !!----                                                                 Out -> Actual Line on Filevar
    !!----    integer,           intent(in)               :: nline_end   !  In -> Line to finish the search
    !!----    character(len=*),  intent(out)              :: spgr_hm     ! Out -> Hermann-Mauguin symbol
    !!----
    !!----    Obtaining the Herman-Mauguin symbol of Space Group
    !!----
    !!---- Update: March - 2010
    !!
    Subroutine Read_Cif_Hm(Filevar,Nline_Ini,Nline_End,Spgr_Hm)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       character(len=*),  intent(out)              :: spgr_hm

       !---- Local variables ----!
       character(len=1) :: csym, csym2
       integer          :: np1, np2

       spgr_hm=" "
       np1=nline_ini
       call Read_Key_Str(filevar,nline_ini,nline_end, &
                            "_symmetry_space_group_name_H-M",spgr_hm)
       !if (len_trim(spgr_hm) ==0 ) spgr_hm=adjustl(filevar(nline_ini+1))
       nline_ini=np1
       ! TR  feb. 2015 .(re-reading the same item with another name)
       if(len_trim(spgr_hm) == 0) then
        spgr_hm = " "
        call Read_Key_Str(filevar,nline_ini,nline_end, "_space_group_name_H-M_alt",spgr_hm)
        if (len_trim(spgr_hm) ==0 ) spgr_hm=adjustl(filevar(nline_ini+1))
       end if

       if (spgr_hm =="?" .or. spgr_hm=="#") then
          spgr_hm=" "
       else
          np1=index(spgr_hm,"'")
          np2=index(spgr_hm,"'",back=.true.)
          if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
             spgr_hm=spgr_hm(np1+1:np2-1)
          else
             np1=index(spgr_hm,'"')
             np2=index(spgr_hm,'"',back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                spgr_hm=spgr_hm(np1+1:np2-1)
             else
                spgr_hm=" "
             end if
          end if
       end if

       !---- Adapting Nomenclature from ICSD to our model ----!
       np1=len_trim(spgr_hm)
       if (np1 > 0) then
          csym=u_case(spgr_hm(np1:np1))
          select case (csym)
             case("1")
                csym2=u_case(spgr_hm(np1-1:np1-1))
                if (csym2 == "Z" .or. csym2 =="S") then
                   spgr_hm=spgr_hm(:np1-2)//":1"
                end if

             case("S","Z")
                csym2=u_case(spgr_hm(np1-1:np1-1))
                select case (csym2)
                   case ("H")
                      spgr_hm=spgr_hm(:np1-2)
                   case ("R")
                      spgr_hm=spgr_hm(:np1-2)//":R"
                   case default
                      spgr_hm=spgr_hm(:np1-1)
                end select

             case("R")
                csym2=u_case(spgr_hm(np1-1:np1-1))
                if (csym2 == "H" ) then
                   spgr_hm=spgr_hm(:np1-2)
                else
                   spgr_hm=spgr_hm(:np1-1)//":R"
                end if
          end select
       end if

       return
    End Subroutine Read_Cif_Hm

    !!----
    !!---- Subroutine Read_Cif_Lambda(Filevar,Nline_Ini,Nline_End,Lambda)
    !!----    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String vector
    !!----    integer,           intent(in out)          :: nline_ini    !  In -> Line to start of search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,           intent(in)              :: nline_end    !  In -> Line to finish the search
    !!----    real(kind=cp),     intent(out)             :: lambda       !  Out -> lamda value
    !!----
    !!----    Radiation length
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Lambda(Filevar,Nline_Ini,Nline_End,Lambda)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       real(kind=cp),     intent(out)              :: lambda

       !---- Local Variables ----!
       integer                    :: iv
       integer,dimension(1)       :: ivet
       real(kind=cp), dimension(1):: vet

       lambda=0.71073    ! Mo

       call read_key_value(filevar,nline_ini,nline_end, &
                           "_diffrn_radiation_wavelength",vet,ivet,iv)
       if (iv == 1) then
          lambda=vet(1)
       end if

       return
    End Subroutine Read_Cif_Lambda

    !!----
    !!---- Subroutine Read_Cif_Symm(Filevar,Nline_Ini,Nline_End,N_Oper,Oper_Symm)
    !!----    character(len=*), dimension(:), intent(in) :: filevar       !  In -> String vector
    !!----    integer,          intent(in out)           :: nline_ini     !  In -> Line to start the search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,          intent(in)               :: nline_end     !  In -> Line to finish the search
    !!----    integer,          intent(out)              :: n_oper        ! Out -> Number of Operators
    !!----    character(len=*), dimension(:),intent(out) :: oper_symm     ! Out -> Vector with Symmetry Operators
    !!----
    !!----    Obtaining Symmetry Operators from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Symm(Filevar,Nline_Ini,Nline_End,N_Oper,Oper_Symm)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: filevar
       integer,          intent(in out)           :: nline_ini
       integer,          intent(in)               :: nline_end
       integer,          intent(out)              :: n_oper
       character(len=*), dimension(:),intent(out) :: oper_symm

       !---- Local variables ----!
       character(len=len(filevar(1))) :: string
       integer                        :: i,np1,np2

       n_oper=0
       oper_symm=" "
       np1=nline_ini
       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_symmetry_equiv_pos_as_xyz",string)
       nline_ini=np1
       ! TR  feb. 2015 .(re-reading the same item with another name)
       if(len_trim(string) == 0) then
        call Read_Key_StrVal(filevar,nline_ini,nline_end, "_space_group_symop_operation_xyz",string)
       end if

       if (len_trim(string) /=0) then
          string=adjustl(string)

          if (string(1:1) /="#" .and. string(1:1) /= "?") then      ! Comentario
             np1=index(string,"'")
             np2=index(string,"'",back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                n_oper=n_oper+1
                oper_symm(n_oper)=string(np1+1:np2-1)
             else
                np1=index(string,'"')
                np2=index(string,'"',back=.true.)
                if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                   n_oper=n_oper+1
                   oper_symm(n_oper)=string(np1+1:np2-1)
                end if
             end if
          end if
       end if

       do i=nline_ini+1,nline_end
          string=adjustl(filevar(i))
          if (len_trim(string) /=0) then
             if (string(1:1) /="#" .and. string(1:1) /= "?") then      ! Comentario o Vacio
                np1=index(string,"'")
                np2=index(string,"'",back=.true.)
                if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                   n_oper=n_oper+1
                   oper_symm(n_oper)=string(np1+1:np2-1)
                else
                   np1=index(string,'"')
                   np2=index(string,'"',back=.true.)
                   if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                      n_oper=n_oper+1
                      oper_symm(n_oper)=string(np1+1:np2-1)
                   end if
                end if
             end if
          else
             nline_ini=i+1
             exit
          end if
       end do

       return
    End Subroutine Read_Cif_Symm

    !!----
    !!---- Subroutine Read_Cif_Title(Filevar,Nline_Ini,Nline_End,Title)
    !!----    character(len=*),  dimension(:), intent(in) :: filevar      !  In -> String vector
    !!----    integer,           intent(in out)           :: nline_ini    !  In -> Line to start the search
    !!----                                                                  Out -> Actual line on Filevar
    !!----    integer,           intent(in)               :: nline_end    !  In -> Line to finish the search
    !!----    character(len=*),  intent(out)              :: title        ! Out -> Title string
    !!----
    !!----    Obtaining Title from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Title(Filevar,Nline_Ini,Nline_End,title)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       character(len=*),  intent(out)              :: title

       !---- Local variables ----!
       integer :: np, np1, np2

       title=" "
       call Read_Key_StrVal(filevar,nline_ini,nline_end, &
                            "_publ_section_title",title)

       if (len_trim(title) ==0 ) title=adjustl(filevar(nline_ini+1))
       if (title =="; ?" .or. title=="#") then
          title=" "
       else
          np=len_trim(title)
          if (np <= 3) title=adjustl(filevar(nline_ini+2))
          np1=index(title,"'")
          np2=index(title,"'",back=.true.)
          if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
             title=title(np1+1:np2-1)
          else
             np1=index(title,'"')
             np2=index(title,'"',back=.true.)
             if (np1 > 0 .and. np2 > 0 .and. np2 > np1) then
                title=title(np1+1:np2-1)
             end if
          end if
       end if

       return
    End Subroutine Read_Cif_Title

    !!----
    !!---- Subroutine Read_Cif_Z(Filevar,Nline_Ini,Nline_End,Z)
    !!----    character(len=*), dimension(:), intent(in) :: filevar     !  In -> String vector
    !!----    integer,           intent(in out)          :: nline_ini   !  In -> Line to start the search
    !!----                                                                Out -> Actual line on Filevar
    !!----    integer,           intent(in)              :: nline_end   !  In -> Line to finish the search
    !!----    integer,           intent(out)             :: Z           ! Out -> Z value
    !!----
    !!----    Unit formula from Cif file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Cif_Z(filevar,nline_ini,nline_end,z)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar
       integer,           intent(in out)           :: nline_ini
       integer,           intent(in)               :: nline_end
       integer,           intent(out)              :: z

       !---- Local Variables ----!
       integer                     :: iv
       integer,dimension(1)        :: ivet
       real(kind=cp), dimension(1) :: vet

       z=0
       call read_key_value(filevar,nline_ini,nline_end, &
                           "_cell_formula_units_Z",vet,ivet,iv)
       if (iv == 1) then
          z=ivet(1)
       end if

       return
    End Subroutine Read_Cif_Z

    !!----
    !!---- Subroutine Read_File_Atom(Filevar,Nline_Ini,Nline_End,Atomos)
    !!----    character(len=*),dimension(:), intent(in)       :: filevar     !  In -> String vector
    !!----    integer,                       intent(in)       :: nline_ini   !  In -> Line to start the search
    !!----                                                                     Out -> Actual line on Filevar
    !!----    integer,                       intent(in)       :: nline_end   !  In -> Line to finish the search
    !!----    type (atom_list_type),        intent(out)      :: Atomos      ! Out -> Atom list
    !!----           or
    !!----    type (Point_list_Type),        intent(out)      :: Atomos      ! Out -> point list
    !!----
    !!----     Subroutine to read an atom (or point) list from a file. Atomos should be previously allocated.
    !!----     Control of error is present.
    !!----
    !!---- Update: June - 2005
    !!

    !!--++
    !!--++ Subroutine Read_File_Atomlist(Filevar,Nline_Ini,Nline_End,Atomos)
    !!--++    character(len=*),dimension(:), intent(in)       :: filevar     !  In -> String vector
    !!--++    integer,                       intent(in)       :: nline_ini   !  In -> Line to start the search
    !!--++                                                                     Out -> Actual line on Filevar
    !!--++    integer,                       intent(in)       :: nline_end   !  In -> Line to finish the search
    !!--++    type (atom_list_type),        intent(out)       :: Atomos      ! Out -> Atom list
    !!--++
    !!--++     Subroutine to read an atom list from a file. Atomos should be previously allocated.
    !!--++     Control of error is present
    !!--++
    !!--++ Update: June - 2005
    !!
    Subroutine Read_File_Atomlist(filevar,nline_ini,nline_end,Atomos)
       !---- Arguments ----!
       character(len=*), dimension(:),   intent(in)      :: filevar
       integer,                          intent(in out)  :: nline_ini
       integer,                          intent(in)      :: nline_end
       type (atom_list_type),            intent(in out)  :: Atomos

       !---- Local variables -----!
       character(len=len(filevar(1))) :: line
       character(len=4)               :: dire
       integer                        :: i,na
       type (Atom_Type)               :: Atomo

       !---- Initial Values ----!
       na=0
       do i=nline_ini,nline_end
          dire=adjustl(u_case(filevar(i)(1:4)))
          if (dire /= "ATOM") cycle
          line=adjustl(filevar(i))
          call read_atom(line,atomo)
          if (err_form) cycle

          !---- Trial to read anisotropic thermal parameters ----!
          if( i < size(filevar) ) then
           line=adjustl(filevar(i+1))
           select case (u_case(line(1:4)))
             case ("U_IJ")
                call read_uvals(line,atomo, "u_ij")
             case ("B_IJ")
                call read_uvals(line,atomo, "b_ij")
             case ("BETA")
                call read_uvals(line,atomo, "beta")
           end select
           if (err_form) cycle
          end if
          na=na+1
          Atomos%atom(na)=atomo
       end do

       Atomos%natoms=na

       return
    End Subroutine Read_File_Atomlist

    !!----
    !!---- Subroutine Read_File_PointList(Filevar,Nline_Ini,Nline_End,Atomos)
    !!----    character(len=*),dimension(:), intent(in)       :: filevar     !  In -> String vector
    !!----    integer,                       intent(in)       :: nline_ini   !  In -> Line to start the search
    !!----                                                                     Out -> Actual line on Filevar
    !!----    integer,                       intent(in)       :: nline_end   !  In -> Line to finish the search
    !!----    type (Point_List_Type),        intent(out)      :: Atomos      ! Out -> point list
    !!----
    !!----     Subroutine to read an point list from a file. Atomos should be previously allocated.
    !!----     Control of error is present
    !!----
    !!---- Update: June - 2005
    !!
    Subroutine Read_File_PointList(filevar,nline_ini,nline_end,Atomos)
       !---- Arguments ----!
       character(len=*), dimension(:),   intent(in)      :: filevar
       integer,                          intent(in out)  :: nline_ini
       integer,                          intent(in)      :: nline_end
       type (Point_List_Type),           intent(in out)  :: Atomos

       !---- Local variables -----!
       character(len=len(filevar(1))) :: line
       character(len=4)               :: dire
       integer                        :: i,na
       type (Atom_Type)               :: Atomo

       !---- Initial Values ----!
       na=0

       do i=nline_ini,nline_end
          dire=adjustl(u_case(filevar(i)(1:4)))
          if (dire /= "ATOM") cycle
          line=adjustl(filevar(i))
          call read_atom(line,atomo)
          if (err_form) cycle
          na=na+1
          Atomos%x(:,na) =atomo%x(:)
          Atomos%p(na)   = 0
          Atomos%nam(na) = atomo%lab
       end do

       Atomos%np=na

       return
    End Subroutine Read_File_PointList

    !!----
    !!---- Subroutine Read_File_Cell(Filevar,Nline_Ini,Nline_End,Celda)
    !!----    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String Vector
    !!----    integer,           intent(in out)          :: nline_ini    !  In -> Line to start the search
    !!----                                                                 Out -> Atual line on Filevar
    !!----    integer,           intent(in)              :: nline_end    !  In -> line to finish the search
    !!----
    !!----    real(kind=cp),dimension(6), intent (out)   :: Celda        ! Out -> Cell variable
    !!----                            or
    !!----    type (Crystal_Cell_Type), intent (out)     :: Celda        ! Out -> Cell variable
    !!----
    !!----    Read Cell Parameters from file. Control error is present
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Read_File_Cellc(Filevar,Nline_Ini,Nline_End,Celda)
    !!--++    character(len=*), dimension(:), intent(in) :: filevar      !  In -> String Vector
    !!--++    integer,           intent(in out)          :: nline_ini    !  In -> Line to start the search
    !!--++                                                                 Out -> Atual line on Filevar
    !!--++    integer,           intent(in)              :: nline_end    !  In -> line to finish the search
    !!--++    real(kind=cp),dimension(6), intent (out)   :: Celda        ! Out -> Cell variable
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Read Cell Parameters from file. Control error is present
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_File_Cellc(filevar,nline_ini,nline_end,Celda)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in)     :: filevar
       integer,                         intent(in)     :: nline_ini
       integer,                         intent(in)     :: nline_end
       real(kind=cp),dimension(6),      intent(out)    :: Celda

       !---- Local Variables ----!
       integer                     :: iv, i,j
       integer, dimension(6)       :: ivet
       real(kind=cp), dimension(6) :: vet

       !---- Valores iniciales ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       !---- Celda ----!
       call read_key_value(filevar,i,j,"cell",vet,ivet,iv)
       if (iv /=6) then
          err_form=.true.
          ERR_Form_Mess=" Bad Cell Parameters..."
          return
       else
          celda=vet(:)
       end if

       return
    End Subroutine Read_File_Cellc

    !!--++
    !!--++ Subroutine Read_File_Cellt(Filevar,Nline_Ini,Nline_End,Celda,CFrame)
    !!--++    character(len=*),  dimension(:), intent(in)     :: filevar     !  In -> String Vector
    !!--++    integer,                         intent(in)     :: nline_ini   !  In -> Line to start the search
    !!--++    integer,                         intent(in)     :: nline_end   !  In -> line to finish the search
    !!--++    type (Crystal_Cell_Type),        intent(out)    :: Celda       ! Out -> Cell structure
    !!--++    character(len=*),  optional,     intent(in)     :: CFrame      !  Cartesian Frame "A" or "C" (if absent -> "A")
    !!--++          ! Out -> Cell variable
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Read Cell Parameters from file. Control error is present
    !!--++    The object Celda is constructed just after reading the cell parameters.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_File_Cellt(filevar,nline_ini,nline_end,Celda,CFrame)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in)     :: filevar
       integer,                         intent(in)     :: nline_ini
       integer,                         intent(in)     :: nline_end
       type (Crystal_Cell_Type),        intent(out)    :: Celda
       character(len=*),  optional,     intent(in)     :: CFrame

       !---- Local Variables ----!
       integer                     :: iv, i,j
       real(kind=cp), dimension(6) :: vet1,vet2

       !---- Valores iniciales ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       !---- Celda ----!

       call read_key_valueSTD(filevar,i,j,"cell",vet1,vet2,iv)
       if (iv /=6) then
          err_form=.true.
          ERR_Form_Mess=" Bad Cell Parameters..."
          return
       end if
       if(present(CFrame)) then
         call Set_Crystal_Cell(vet1(1:3),vet1(4:6),Celda,CFrame)
       else
         call Set_Crystal_Cell(vet1(1:3),vet1(4:6),Celda,"A")
       end if
       celda%cell_std=vet2(1:3)
       celda%ang_std=vet2(4:6)

       return
    End Subroutine Read_File_Cellt

    !!----
    !!---- Subroutine Read_File_lambda(Filevar,Nline_Ini,Nline_End,v1,v2,v3)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar   !  In -> String Vector
    !!----    integer,                        intent(in out) :: nline_ini !  In -> Line to start the search
    !!----                                                                  Out -> Atual line on Filevar
    !!----    integer,                        intent(in)     :: nline_end !  In -> line to finish the search
    !!----    real(kind=cp),                  intent(   out) :: v1,v2,v3  ! Out -> Lambda1,lambda2,ratio
    !!----
    !!----    Read wavelengths and ratio.
    !!----    If no value is read, Lambda1=Lambda2=1.54056 Angstroms, ratio=0.0
    !!----    If only one value is read Lambda1=Lambda2=v1, ratio=0
    !!----    If only two values iare read Lambda1=v1, Lambda2=v2, ratio=0.5
    !!----    In other cases Lambda1=v1, Lambda2=v2, ratio=v3
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_File_Lambda(Filevar,Nline_Ini,Nline_End,v1,v2,v3)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)     :: filevar
       integer,                        intent(in out) :: nline_ini
       integer,                        intent(in)     :: nline_end
       real(kind=cp),                  intent(   out) :: v1,v2,v3

       !---- Local Variables ----!
       integer                    :: iv, i,j
       integer, dimension(3)      :: ivet
       real(kind=cp), dimension(3):: vet

       !---- Valores iniciales ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       v3=0.0
       v1=1.54056
       !---- Read Lambda ----!
       call read_key_value(filevar,i,j,"wave",vet,ivet,iv)
       if      (iv == 0) then
         v2=1.54056
       else if (iv == 1) then
         v1=vet(1)
         v2=vet(1)
       else if (iv == 2) then
         v1=vet(1)
         v2=vet(2)
         v3=0.5
       else if (iv == 3) then
         v1=vet(1)
         v2=vet(2)
         v3=vet(3)
       end if

       return
    End Subroutine Read_File_Lambda

    !!----
    !!---- Subroutine Read_File_RngSintL(Filevar,Nline_Ini,Nline_End,v1,v2)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar   !  In -> String Vector
    !!----    integer,                        intent(in out) :: nline_ini !  In -> Line to start the search
    !!----                                                                  Out -> Atual line on Filevar
    !!----    integer,                        intent(in)     :: nline_end !  In -> line to finish the search
    !!----    real(kind=cp),                  intent(   out) :: v1,v2     ! Out -> Interval [v1,v2] in sinT/Lambda
    !!----
    !!----    Read range for sintheta/lambda.
    !!----    If only one value is read v1=0 and v2= read value
    !!----    If the keyword RNGSL is not given in the file, the default
    !!----    values are v1=0.0, v2=1.0
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_File_RngSintL(Filevar,Nline_Ini,Nline_End,v1,v2)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)     :: filevar
       integer,                        intent(in out) :: nline_ini
       integer,                        intent(in)     :: nline_end
       real(kind=cp),                  intent(   out) :: v1,v2

       !---- Local Variables ----!
       integer                     :: iv, i,j
       integer,       dimension(2) :: ivet
       real(kind=cp), dimension(2) :: vet

       !---- Valores iniciales ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       !---- Range in sinTheta/Lambda ----!
       call read_key_value(filevar,i,j,"rngsl",vet,ivet,iv)
       if      (iv == 0) then
         v1=0.0
         v2=1.0
       else if (iv == 1) then
         v1=0.0
         v2=vet(1)
       else if (iv == 2) then
         v1=vet(1)
         v2=vet(2)
       end if

       return
    End Subroutine Read_File_RngSintL

    !!----
    !!---- Subroutine Read_File_Spg (Filevar,Nline_Ini,Nline_End,Spg,Sub)
    !!----    character(len=*),  dimension(:), intent(in) :: filevar       !  In -> String vector
    !!----    integer,           intent(in)               :: nline_ini     !  In -> Line to start the search
    !!----                                                                   Out -> Actual line on Filevar
    !!----    integer,           intent(in)               :: nline_end     !  In -> Line to Finish the search
    !!----    character(len=*),  intent(out)              :: spg           ! Out -> Space Group symbol
    !!----    character(len=*),  intent(in ),optional     :: sub           ! in  -> The space sroup symbol is a subgroup
    !!----                                                                 !        of an already given space group
    !!----    Reads the cards "SPGR", "SPACEG" or "SUBG" in filvar. Control of error is present
    !!----
    !!---- Update: February - 2011
    !!
    Subroutine Read_File_Spg(filevar,nline_ini,nline_end,Spg,sub)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in) :: filevar   ! Variable
       integer,           intent(in)               :: nline_ini
       integer,           intent(in)               :: nline_end
       character(len=*),  intent(out)              :: spg
       character(len=*),  intent(in),  optional    :: sub

       !--Local variables--!
       integer  :: i

       call init_err_form()
       i=nline_ini
       if(present(sub)) then
         call Read_Key_StrVal(filevar,i,nline_end, "subg",spg)
       else
         call Read_Key_StrVal(filevar,i,nline_end, "spgr",spg)
       end if
       if (len_trim(spg) == 0 ) then
         call Read_Key_StrVal(filevar,i,nline_end, "spaceg",spg)
         if (len_trim(spg) == 0 ) then
           err_form=.true.
           ERR_Form_Mess=" Problems reading the Space Group symbol/number"
           return
         end if
       end if

       return
    End Subroutine Read_File_Spg

    !!----
    !!---- Read_File_Transf(Filevar,Nline_Ini,Nline_End,Transf,Orig)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar      !  In -> String Vector
    !!----    integer,                        intent(in out) :: nline_ini    !  In -> Line to start the search
    !!----                                                                     Out -> Atual line on Filevar
    !!----    integer,                        intent(in)     :: nline_end    !  In -> line to finish the search
    !!----    real(kind=cp),dimension(3,3),   intent(out)    :: transf       ! Out -> Cell variable
    !!----    real(kind=cp),dimension(3  ),   intent(out)    :: orig
    !!----
    !!----    Read transformation matrix for changing the space group or cell setting.
    !!----    First the matrix M is read row by row and then the origin in the old setting
    !!----    is finally read. A single line with 12 real numbers should be given.
    !!--<<
    !!----    e.g.: TRANS  m11 m12 m13  m21 m22 m33  m31 m32 m33   o1 o2 o3
    !!----
    !!----    That means       a'=m11 a + m12 b + m13 c
    !!----                     b'=m21 a + m22 b + m23 c
    !!----                     c'=m31 a + m32 b + m33 c
    !!----
    !!----                     X' = inv(Mt) (X-O)
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_File_transf(filevar,nline_ini,nline_end,trans,orig)
       !---- Arguments ----!
       character(len=*),  dimension(:), intent(in)     :: filevar
       integer,                         intent(in)     :: nline_ini
       integer,                         intent(in)     :: nline_end
       real(kind=cp),dimension(3,3),    intent(out)    :: trans
       real(kind=cp),dimension(3  ),    intent(out)    :: orig

       !---- Local Variables ----!
       integer                      :: iv, i,j
       integer,       dimension(12) :: ivet
       real(kind=cp), dimension(12) :: vet
       character(len=80)            :: transf_key

       !---- Initial values ----!
       call init_err_form()

       i=nline_ini
       j=nline_end

       !---- transformation matrix ----!
       call read_key_value(filevar,i,j,"trans",vet,ivet,iv,"#",transf_key)
       if (iv /= 12 .or. err_string) then
          !Try to read the transformation from transf_key
          if(len_trim(transf_key) /= 0) then
            call Get_Transf(transf_key,trans,orig)
            if(err_string) then
               err_form=.true.
               ERR_Form_Mess=" Bad matrix/origin setting in string: "//trim(transf_key)//" -> "//trim(Err_String_Mess)
               return
            end if
          else
               err_form=.true.
               ERR_Form_Mess=" Bad matrix/origin setting..."
               return
          end if

       else
          trans(1,1:3)=vet(1:3)
          trans(2,1:3)=vet(4:6)
          trans(3,1:3)=vet(7:9)
          orig(1:3) = vet(10:12)
       end if

       return
    End Subroutine Read_File_transf

    !!----
    !!---- Subroutine Read_Shx_Atom(Filevar,Nline_Ini,Nline_End,N_Fvar,Fvar,Elem_Type,Celda,Atm_List)
    !!----    character(len=*), dimension(:), intent(in)      :: filevar        !  In -> String vector
    !!----    integer,                        intent(in out)  :: nline_ini      !  In -> Line to start the search
    !!----                                                                         Out -> Actual line on Filevar
    !!----    integer,                        intent(in)      :: nline_end      !  In -> Line to finish the search
    !!----    integer,                        intent(in)      :: n_fvar         !  In -> Number of parameters on FVAR
    !!----    real(kind=cp), dimension(:),    intent(in)      :: fvar           !  In -> Values for FVAR
    !!----    character(len=*), dimension(:), intent(in)      :: elem_type      !  In -> type of elements
    !!----    type (Crystal_Cell_Type),       intent(in)      :: Celda          !  In -> Cell type variable
    !!----    type (Atom_list_type),          intent(out)     :: Atm_List       ! Out -> number of atoms
    !!----         ! Out -> Atom List
    !!----
    !!----    Obtaining Atoms parameters from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Atom(filevar,nline_ini,nline_end,n_fvar,fvar,elem_type,celda,Atm_List)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)      :: filevar
       integer,                        intent(in out)  :: nline_ini
       integer,                        intent(in)      :: nline_end
       integer,                        intent(in)      :: n_fvar
       real(kind=cp), dimension(:),    intent(in)      :: fvar
       character(len=*), dimension(:), intent(in)      :: elem_type
       type (Crystal_Cell_Type),       intent(in)      :: Celda
       type (Atom_list_type),          intent(out)     :: Atm_List

       !---- Local Variables ----!
       character(len=80)               :: string
       character(len=30),dimension(15) :: label
       character(len=2)                :: el
       integer                         :: i, nc, iv
       integer                         :: j, n_atom
       integer, dimension(15)          :: ivet
       real(kind=cp)                   :: x, p, u
       real(kind=cp), dimension(15)    :: vet
       type(atom_list_type)            :: Atm

       call allocate_atom_list(nline_end-nline_ini+1,Atm)
       n_atom=0

       do i=nline_ini,nline_end
          string=filevar(i)
          if (len_trim(string) == 0) cycle
          call getword(string,label,nc)
          select case (nc)
             case (5) ! Atomname Sfac X Y Z
                call getnum(label(2),vet,ivet,iv)   ! Is Sfac integer?
                if (iv /= 1) cycle
                call getnum(label(3),vet,ivet,iv)   ! Is X real?
                if (iv /= 1) cycle
                call getnum(label(4),vet,ivet,iv)   ! Is Y real?
                if (iv /= 1) cycle
                call getnum(label(5),vet,ivet,iv)   ! Is Z real?
                if (iv /= 1) cycle

                n_atom=n_atom+1
                atm%atom(n_atom)%lab=label(1)(1:4)
                call getnum(label(2),vet,ivet,iv)
                el=elem_type(ivet(1))
                atm%atom(n_atom)%chemSymb=U_case(el(1:1))//L_case(el(2:2))
                call getnum(label(3),vet,ivet,iv)
                atm%atom(n_atom)%x(1)=vet(1)
                call getnum(label(4),vet,ivet,iv)
                atm%atom(n_atom)%x(2)=vet(1)
                call getnum(label(5),vet,ivet,iv)
                atm%atom(n_atom)%x(3)=vet(1)
                atm%atom(n_atom)%utype="u_ij"

             case (6) ! Atomname Sfac X Y Z Occ
                call getnum(label(2),vet,ivet,iv)   ! Is Sfac integer?
                if (iv /= 1) cycle
                call getnum(label(3),vet,ivet,iv)   ! Is X real?
                if (iv /= 1) cycle
                call getnum(label(4),vet,ivet,iv)   ! Is Y real?
                if (iv /= 1) cycle
                call getnum(label(5),vet,ivet,iv)   ! Is Z real?
                if (iv /= 1) cycle
                call getnum(label(6),vet,ivet,iv)   ! Is Occ real?
                if (iv /= 1) cycle

                n_atom=n_atom+1
                atm%atom(n_atom)%lab=label(1)(1:4)
                call getnum(label(2),vet,ivet,iv)
                el=elem_type(ivet(1))
                atm%atom(n_atom)%chemSymb=U_case(el(1:1))//L_case(el(2:2))
                call getnum(label(3),vet,ivet,iv)
                atm%atom(n_atom)%x(1)=vet(1)
                call getnum(label(4),vet,ivet,iv)
                atm%atom(n_atom)%x(2)=vet(1)
                call getnum(label(5),vet,ivet,iv)
                atm%atom(n_atom)%x(3)=vet(1)
                call getnum(label(6),vet,ivet,iv)
                atm%atom(n_atom)%occ=vet(1)
                atm%atom(n_atom)%utype="u_ij"

             case (7,8) ! Atomname Sfac X Y Z Occ Uiso   (TR: item 8 can be electronic density created by SHELXS)
                call getnum(label(2),vet,ivet,iv)   ! Is Sfac integer?
                if (iv /= 1) cycle
                call getnum(label(3),vet,ivet,iv)   ! Is X real?
                if (iv /= 1) cycle
                call getnum(label(4),vet,ivet,iv)   ! Is Y real?
                if (iv /= 1) cycle
                call getnum(label(5),vet,ivet,iv)   ! Is Z real?
                if (iv /= 1) cycle
                call getnum(label(6),vet,ivet,iv)   ! Is Occ real?
                if (iv /= 1) cycle
                call getnum(label(7),vet,ivet,iv)   ! Is Uiso real?
                if (iv /= 1) cycle

                n_atom=n_atom+1
                atm%atom(n_atom)%lab=label(1)(1:4)
                call getnum(label(2),vet,ivet,iv)
                el=elem_type(ivet(1))
                atm%atom(n_atom)%chemSymb=U_case(el(1:1))//L_case(el(2:2))
                call getnum(label(3),vet,ivet,iv)
                atm%atom(n_atom)%x(1)=vet(1)
                call getnum(label(4),vet,ivet,iv)
                atm%atom(n_atom)%x(2)=vet(1)
                call getnum(label(5),vet,ivet,iv)
                atm%atom(n_atom)%x(3)=vet(1)
                call getnum(label(6),vet,ivet,iv)
                atm%atom(n_atom)%occ=vet(1)
                call getnum(label(7),vet,ivet,iv)
                atm%atom(n_atom)%ueq=vet(1)
                atm%atom(n_atom)%utype="u_ij"
                atm%atom(n_atom)%thtype="isotr"

          case (9) ! Atomname Sfac X Y Z Occ U11 U22 = U33 U23 U13 U12
                call getnum(label(2),vet,ivet,iv)   ! Is Sfac integer?
                if (iv /= 1) cycle
                call getnum(label(3),vet,ivet,iv)   ! Is X real?
                if (iv /= 1) cycle
                call getnum(label(4),vet,ivet,iv)   ! Is Y real?
                if (iv /= 1) cycle
                call getnum(label(5),vet,ivet,iv)   ! Is Z real?
                if (iv /= 1) cycle
                call getnum(label(6),vet,ivet,iv)   ! Is Occ real?
                if (iv /= 1) cycle
                call getnum(label(7),vet,ivet,iv)   ! Is U11 real?
                if (iv /= 1) cycle
                call getnum(label(8),vet,ivet,iv)   ! Is U22 real?
                if (iv /= 1) cycle
                call getnum(filevar(i+1),vet,ivet,iv) ! Are U33 U23 U13 U12?
                if (iv /= 4) cycle

                n_atom=n_atom+1
                atm%atom(n_atom)%lab=label(1)(1:4)
                call getnum(label(2),vet,ivet,iv)
                el=elem_type(ivet(1))
                atm%atom(n_atom)%chemSymb=U_case(el(1:1))//L_case(el(2:2))
                call getnum(label(3),vet,ivet,iv)
                atm%atom(n_atom)%x(1)=vet(1)
                call getnum(label(4),vet,ivet,iv)
                atm%atom(n_atom)%x(2)=vet(1)
                call getnum(label(5),vet,ivet,iv)
                atm%atom(n_atom)%x(3)=vet(1)
                call getnum(label(6),vet,ivet,iv)
                atm%atom(n_atom)%occ=vet(1)
                !---- U11 U22 U33 U12 U13 U23 Order ----!
                call getnum(label(7),vet,ivet,iv)
                atm%atom(n_atom)%u(1)=vet(1)
                call getnum(label(8),vet,ivet,iv)
                atm%atom(n_atom)%u(2)=vet(1)
                call getnum(filevar(i+1),vet,ivet,iv)
                atm%atom(n_atom)%u(3)=vet(1)
                atm%atom(n_atom)%u(4)=vet(4)
                atm%atom(n_atom)%u(5)=vet(3)
                atm%atom(n_atom)%u(6)=vet(2)
                atm%atom(n_atom)%utype="u_ij"
                atm%atom(n_atom)%thtype="aniso"
             case default
                cycle
          end select
       end do

       !---- Adjusting ... ----!
       call allocate_atom_list(n_atom,Atm_list)
       do i=1,n_atom
          atm_list%atom(i)=atm%atom(i)
       end do
       call Deallocate_atom_list(atm)

       !---- Tratamiento de Datos del Shelx ----!
       do i=1,n_atom
          !---- coordinates ----!
          if (atm_list%atom(i)%x(1) >= 10.0) atm_list%atom(i)%x(1)=atm_list%atom(i)%x(1)-10.0
          if (atm_list%atom(i)%x(2) >= 10.0) atm_list%atom(i)%x(2)=atm_list%atom(i)%x(2)-10.0
          if (atm_list%atom(i)%x(3) >= 10.0) atm_list%atom(i)%x(3)=atm_list%atom(i)%x(3)-10.0

          !---- ocupancy ----!
          if (abs(atm_list%atom(i)%occ)  > 10.0) then
             x=atm_list%atom(i)%occ
             if (x > 10.0) then
                atm_list%atom(i)%occ=x-10.0
             else
                x=abs(atm_list%atom(i)%occ)
                do j=2,n_fvar
                   if (x > 10.0*real(j) .and. x < 10.0*real(j+1)) then
                      p=x-10.0*real(j)
                      if (atm_list%atom(i)%occ > 0.0) then
                         atm_list%atom(i)%occ=p*fvar(j)
                      else
                         atm_list%atom(i)%occ=p*(fvar(j)-1.0)
                      end if
                   end if
                end do
             end if
          end if

          !---- Thermal factors ----!
          if (atm_list%atom(i)%thtype == "aniso") then
             atm_list%atom(i)%ueq=U_Equiv(celda,atm_list%atom(i)%u(1:6))  ! Uequi
             atm_list%atom(i)%biso= atm_list%atom(i)%ueq*78.95683521
          else
             if (atm_list%atom(i)%ueq < 0.0) then
                u=-atm_list%atom(i)%ueq
                if (u <= 5.0 .and. u >= 0.5) then
                   do j=i-1,1,-1
                      if (atm_list%atom(j)%ChemSymb == "H " .or. atm_list%atom(j)%ChemSymb == "h " ) cycle
                      atm_list%atom(i)%ueq=u*U_Equiv(celda,atm_list%atom(j)%u(1:6))  ! Uequi
                      atm_list%atom(i)%biso= atm_list%atom(i)%ueq*78.95683521
                   end do
                end if
             end if
          end if

       end do

       return
    End Subroutine Read_Shx_Atom

    !!----
    !!---- Subroutine Read_Shx_Cell(Filevar,Nline_Ini,Nline_End,Celda,Stdcelda,Lambda,Z)
    !!----    character(len=*), dimension(:), intent(in)     :: filevar       !  In -> String vector
    !!----    integer,                        intent(in out) :: nline_ini     !  In -> Line to start the search
    !!----                                                                      Out -> Actual line on Filevar
    !!----    integer,                        intent(in)     :: nline_end     !  In -> Line to finish the search
    !!----    real(kind=cp),dimension(6),     intent(out)    :: celda         ! Out -> Cell Parameters
    !!----    real(kind=cp),dimension(6),     intent(out)    :: Stdcelda      ! Out -> Std Cell Parameters
    !!----    real(kind=cp),                  intent(out)    :: lambda        ! Out -> Lambda
    !!----    integer,                        intent(out)    :: Z             ! Out -> Z
    !!----
    !!----    Obtaining Cell Parameter from Shelx file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Cell(filevar,nline_ini,nline_end,Celda,StdCelda,lambda,z)
       !---- Arguments ----!
       character(len=*), dimension(:),     intent(in)     :: filevar
       integer,                            intent(in out) :: nline_ini
       integer,                            intent(in)     :: nline_end
       real(kind=cp),dimension(6),         intent(out)    :: Celda
       real(kind=cp),dimension(6),optional,intent(out)    :: StdCelda
       real(kind=cp),             optional,intent(out)    :: lambda
       integer,          optional,         intent(out)    :: z

       !---- Local Variables ----!
       integer                      :: iv,z_shx
       integer, dimension(10)       :: ivet
       real(kind=cp), dimension(10) :: vet
       real(kind=cp)                :: lambda_shx
       real(kind=cp),dimension(6)   :: std

       !---- Valores iniciales ----!
       celda=0.0
       if (present(stdcelda)) stdcelda=0.0
       if (present(Lambda))   lambda=0.0
       if (present(z))        z=0

       !---- CELL ----!
       call read_key_value(filevar,nline_ini,nline_end,"CELL",vet,ivet,iv)
       if (iv == 7) then
          lambda_shx = vet(1)
          celda      = vet(2:7)
       end if

       !---- Z, STD ----!
       call read_key_value(filevar,nline_ini,nline_end,"ZERR",vet,ivet,iv)
       if (iv == 7) then
          z_shx= ivet(1)
          std  = vet(2:7)
       end if

       if (present(stdcelda)) stdcelda=std
       if (present(lambda)) lambda=lambda_shx
       if (present(z)) z=z_shx

       return
    End Subroutine Read_Shx_Cell

    !!----
    !!---- Subroutine Read_Shx_Cont(Filevar,Nline_Ini,Nline_End,N_Elem_Type,Elem_Type,N_Elem)
    !!----    character(len=*),  dimension(:),    intent(in)    :: filevar       !  In -> String Vector
    !!----    integer,                            intent(in out):: nline_ini     !  In -> Line to start the search
    !!----                                                                         Out -> Actual Line on Filevar
    !!----    integer,                            intent(in)    :: nline_end     !  In -> Line to finish the search
    !!----    integer,                            intent(out)   :: n_elem_type   ! Out -> N. of different species
    !!----    character(len=*), dimension(:),     intent(out)   :: elem_type     ! Out -> Character to identify the specie
    !!----    real(kind=cp),dimension(:),optional,intent(out)   :: n_elem        ! Out -> Number of elements into the same species
    !!----
    !!----    Obtaining Chemical contents from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Cont(filevar,nline_ini,nline_end,n_elem_type,elem_type,n_elem)
       !---- Arguments ----!
       character(len=*), dimension(:),     intent(in)      :: filevar
       integer,                            intent(in out)  :: nline_ini
       integer,                            intent(in)      :: nline_end
       integer,                            intent(out)     :: n_elem_type
       character(len=*), dimension(:),     intent(out)     :: elem_type
       real(kind=cp),dimension(:),optional,intent(out)     :: n_elem

       !---- Local  variables ----!
       character(len=len(filevar(1)))      :: string
       integer                     :: iv
       integer,      dimension(15) :: ivet
       real(kind=cp),dimension(15) :: vet

       n_elem_type = 0
       elem_type   = " "
       if (present(n_elem)) n_elem = 0.0

       call Read_Key_StrVal(filevar,nline_ini,nline_end,"SFAC",string)
       if (len_trim(string) /=0) then
          call getword(string,elem_type,n_elem_type)
       end if

       if (present(n_elem)) then
          call read_key_value(filevar,nline_ini,nline_end,"UNIT",vet,ivet,iv)
          if (iv /= 0) n_elem=vet
       end if

       return
    End Subroutine Read_Shx_Cont

    !!----
    !!---- Subroutine Read_Shx_Fvar(Filevar,Nline_Ini,Nline_End,N_Fvar,Fvar)
    !!----    character(len=*), dimension(:), intent(in)    :: filevar       !  In -> String vector
    !!----    integer,                        intent(in out):: nline_ini     !  In -> Line to start the search
    !!----                                                                   ! Out -> Actual line on Filevar
    !!----    integer,                        intent(in)    :: nline_end     !  In -> Line to finish the search
    !!----    integer,                        intent(out)   :: n_fvar        ! Out -> N. of parameters on FVAR
    !!----    real(kind=cp), dimension(:),    intent(out)   :: fvar          ! Out -> values of FVAR
    !!----
    !!----    Obtaining Fvar parameters from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Read_Shx_Fvar(filevar,nline_ini,nline_end,n_fvar,fvar)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)    :: filevar
       integer,                        intent(in out):: nline_ini
       integer,                        intent(in)    :: nline_end
       integer,                        intent(out)   :: n_fvar
       real(kind=cp), dimension(:),    intent(out)   :: fvar

       !---- Local  variables ----!
       integer                      :: iv
       integer,       dimension(15) :: ivet
       real(kind=cp), dimension(15) :: vet

       n_fvar = 1
       fvar   = 1.0

       call read_key_value(filevar,nline_ini,nline_end,"FVAR",vet,ivet,iv)
       if (iv /= 0) then
          n_fvar=iv
          fvar=vet
       end if

       return
    End Subroutine Read_Shx_Fvar

    !!----
    !!---- Subroutine Read_Shx_Latt(Filevar,Nline_Ini,Nline_End,Latt)
    !!----    character(len=*), dimension(:), intent(in) :: filevar     !  In -> String Vector
    !!----    integer,           intent(in out)          :: nline_ini   !  In -> Line to start the search
    !!----                                                                Out -> Actual line on Filevar
    !!----    integer,           intent(in)              :: nline_end   !  In -> Line to finish the search
    !!----    integer,           intent(out)             :: latt        ! Out -> Lattice number
    !!----
    !!----    Obtaining lattice from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Latt(filevar,nline_ini,nline_end,latt)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: filevar
       integer,           intent(in out)          :: nline_ini
       integer,           intent(in)              :: nline_end
       integer,           intent(out)             :: latt

       !---- Local Variables ----!
       integer                     :: iv
       integer,       dimension(2) :: ivet
       real(kind=cp), dimension(2) :: vet

       latt=1
       call read_key_value(filevar,nline_ini,nline_end,"LATT",vet,ivet,iv)
       if (iv == 1) latt = ivet(1)

       return
    End Subroutine Read_Shx_Latt

    !!----
    !!---- Subroutine Read_Shx_Symm(Filevar,Nline_Ini,Nline_End,N_Oper,Oper_Symm)
    !!----    character(len=*), dimension(:), intent(in) :: filevar       !  In -> String Vector
    !!----    integer,           intent(in out)          :: nline_ini     !  In -> Line to start the search
    !!----                                                                  Out -> Actual Line on Filevar
    !!----    integer,           intent(in)              :: nline_end     !  In -> Line to finish the search
    !!----    integer,           intent(out)             :: n_oper        ! Out -> Number of Operators
    !!----    character(len=*), dimension(:),intent(out) :: oper_symm     ! Out -> String for Symmetry Operators
    !!----
    !!----    Obtaining Symmetry Operators from Shelx file (.ins or .res)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Symm(filevar,nline_ini,nline_end,n_oper,oper_symm)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in) :: filevar
       integer,          intent(in out)           :: nline_ini
       integer,          intent(in)               :: nline_end
       integer,          intent(out)              :: n_oper
       character(len=*), dimension(:),intent(out) :: oper_symm

       !---- Local variables ----!
       character(len=80) :: string
       integer           :: nline

       n_oper=0
       oper_symm=" "

       do
          call Read_Key_StrVal(filevar,nline_ini,nline_end,"SYMM",string)
          if (len_trim(string) /=0) then
             n_oper=n_oper+1
             oper_symm(n_oper)=string
             nline_ini=nline_ini+1
             nline=nline_ini
          else
             exit
          end if
       end do
       nline_ini=nline

       return
    End Subroutine Read_Shx_Symm

    !!----
    !!---- Subroutine Read_Shx_Titl(Filevar,Nline_Ini,Nline_End,Title)
    !!----    character(len=*),dimension(:), intent(in)     :: filevar      !  In -> String Vector
    !!----    integer,                       intent(in out) :: nline_ini    !  In -> Line to start the search
    !!----                                                                    Out -> Actual Line on Filevar
    !!----    integer,                       intent(in)     :: nline_end    !  In -> Line to finish the search
    !!----    character(len=*),              intent(out)    :: title        ! Out -> Title
    !!----
    !!----    Obtaining Title from Shelx file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Shx_Titl(filevar,nline_ini,nline_end,Title)
       !---- Arguments ----!
       character(len=*),dimension(:), intent(in)     :: filevar
       integer,                       intent(in out) :: nline_ini
       integer,                       intent(in)     :: nline_end
       character(len=*),              intent(out)    :: title

       call Read_Key_StrVal(filevar,nline_ini,nline_end,"TITL",title)

       return
    End Subroutine Read_Shx_Titl

    !!----
    !!---- Subroutine Read_Uvals(Line,Atomo,Ulabel)
    !!----    character(len=*),  intent(in out)  :: line      !  In -> String
    !!----    Type (Atom_Type),  intent(in out)  :: Atomo     !  In -> Atomo variable
    !!----                                                      Out ->
    !!----    character(len=4),  intent(in)      :: ulabel    !  In -> u_ij, b_ij, beta
    !!----
    !!----    Subroutine to read the anisotropic thermal parameters from a given Line
    !!----    it complets the object Atomo of type Atom.
    !!----    Assumes the string Line has been read from a file and
    !!----    starts with one of the words (u_ij, b_ij or beta), that is removed before reading
    !!----    the values of the parameters.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Uvals(Line,Atomo,Ulabel)
       !---- Arguments ----!
       character(len=*),  intent(in )     :: line
       Type (Atom_Type),  intent(in out)  :: Atomo
       character(len=4),  intent(in)      :: ulabel

       !---- Local variables -----!
       character(len=len(line)),dimension(1):: line2
       real(kind=cp), dimension (6)         :: vet1,vet2
       integer                              :: iv,n

       call init_err_form()

       atomo%utype    = ulabel
       line2(1)=line
       n=1
       call cutst(line2(1))
       line2(1)="Uval "//line2(1)
       call Read_Key_ValueSTD(line2,n,n,"Uval",vet1,vet2,iv)

        if (iv /= 6) then
          err_form=.true.
          ERR_Form_Mess="  Error reading the anisotropic thermal parameters of atom:"//atomo%lab
          return
       end if
       atomo%U(1:6)=vet1(1:6)
       atomo%U_std(1:6)=vet2(1:6)
       atomo%thtype="aniso"

       return
    End Subroutine Read_Uvals

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++    Type(Job_Info_type), optional,intent(out)  :: Job_Info
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CFL File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CFL(file_dat,nlines,Cell,SpG,A,CFrame,NPhase,Job_Info)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       character(len=*),    optional,intent(in)   :: CFrame
       Integer,             optional,intent( in)  :: Nphase
       Type(Job_Info_type), optional,intent(out)  :: Job_Info

       !---- Local variables ----!
       character(len=132)               :: line
       character(len= 20)               :: Spp
       character(len= 40),dimension(192):: gen
       integer                          :: i, nauas, ndata, iph, n_ini,n_end,ngen,k,nsym
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip

       real(kind=cp),dimension(3):: vet

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!
       iph=1
       if (present(nphase)) iph=nphase
       if (present(Job_Info)) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Get_Job_Info(file_dat,n_ini,n_end,Job_info)
       end if

       !---- Reading Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if(present(CFrame)) then
         call read_File_Cell(file_dat,n_ini,n_end,Cell,CFrame) !Read and construct Cell
       else
         call read_File_Cell(file_dat,n_ini,n_end,Cell) !Read and construct Cell
       end if
       if (err_form) return

       !---- Reading Space Group Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call read_File_Spg (file_dat,n_ini,n_end,Spp)
       if (err_form) then !Try to read symmetry operators or generators
         ngen=0
         nsym=0
         do i=n_ini, n_end
           line=l_case(adjustl(file_dat(i)))
           if(line(1:4) == "symm") nsym=nsym+1
           if(line(1:3) == "gen")  ngen=ngen+1
         end do
         if(ngen > 0) then
           k=0
           do i=n_ini, n_end
             line=l_case(adjustl(file_dat(i)))
             if(line(1:3) == "gen")  then
              k=k+1
              gen(k)=adjustl(line(5:))
             end if
           end do
           call Set_SpaceGroup(" ",SpG,gen,ngen,"gen")   !Construct the space group from generators
         else if (nsym > 0) then
           k=0
           do i=n_ini, n_end
             line=l_case(adjustl(file_dat(i)))
             if(line(1:4) == "symm")  then
              k=k+1
              gen(k)=adjustl(line(6:))
             end if
           end do
           call Set_SpaceGroup(" ",SpG,gen,nsym,"fix")  !Construct the space group from fixed symmetry elements
         else
           return
         end if
       else
          call Set_SpaceGroup(Spp,SpG) !Construct the space group
       end if
       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)

       !---- Calculating number of Atoms in the Phase ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (l_case(line(1:4)) == "atom")  nauas=nauas+1
       end do

       if (nauas > 0) then
          call Allocate_atom_list(nauas,A)  !allocation space for Atom list
          call read_File_Atom(file_dat,n_ini,n_end,A)
          if (err_form) return

          do i=1,A%natoms
             vet=A%atom(i)%x
             A%atom(i)%Mult=Get_Multip_Pos(vet,SpG)
             if(A%atom(i)%occ < epss) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)
             if (A%atom(i)%thtype == "aniso") then
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) =  Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"
             end if
          end do
       end if

       return
    End Subroutine Readn_Set_XTal_CFL

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CFL_Molec(file_dat, nlines, Molcrys)
    !!--++    character(len=*),dimension(:),  intent(in)     :: file_dat
    !!--++    integer,                        intent(in)     :: nlines
    !!--++    Type (Molecular_Crystal_Type),  intent(in out) :: Molcrys
    !!--++
    !!--++ (Private)
    !!--++ Read Molecule Information in a CFL
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CFL_Molec(file_dat, nlines, Molcrys, Nphase)
       !---- Arguments ----!
       character(len=*),dimension(:),  intent(in)     :: file_dat
       integer,                        intent(in)     :: nlines
       type (Molecular_Crystal_Type),  intent(in out) :: Molcrys
       Integer, optional,              intent(in)     :: Nphase

       !---- Local variables ----!
       character(len=132)            :: line
       integer                       :: i,n,nmol,npos,n_ini,n_end,ierr,nauas, iph, ndata
       integer, parameter               :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)        :: ip
       real(kind=cp)                 :: theta,phi,chi
       real(kind=cp), dimension(3)   :: x1f,x2f,x3f
       real(kind=cp), dimension(3,3) :: EuM

       !---- Standard CrysFML file *.CFL ----!
       nauas=0
       ndata=0
       ip=nlines
       ip(1)=1

       !---- Calculating number of Phases ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:6)) == "phase_")  then
             ndata=ndata+1
             ip(ndata)=i
          end if
       end do

       !---- Reading Phase Information ----!

       if (present(nphase)) then
           iph=nphase
       else
           iph=1
       end if

       n_ini=ip(iph)
       n_end=ip(iph+1)

       !---- Detecting the Molecules defined in the file ----!
       nmol=0
       do i=n_ini,n_end
          line=u_case(adjustl(file_dat(i)))
          if (line(1:1) == " ") cycle
          if (line(1:1) == "!") cycle
          npos=index(line,"MOLE")
          if (npos /= 0) nmol=nmol+1
       end do
       if (nmol==0) return

       !---- Allocating Memory for all molecules ----!
       if (allocated(molcrys%mol)) deallocate(molcrys%mol)
       molcrys%n_mol=nmol
       allocate(molcrys%mol(nmol))

       !---- Reading Molecules ----!

       do n=1,nmol
          !---- Read ----!
          do i=n_ini,n_end
             line=u_case(adjustl(file_dat(i)))
             if (line(1:1) == " ") cycle
             if (line(1:1) == "!") cycle
             npos=index(line,"MOLE")
             if (npos == 0) cycle
             call read_molecule(file_dat,n_ini,n_end,molcrys%mol(n))
             err_form=err_molec
             ERR_Form_Mess=err_molec_mess
             if (err_form) then
                molcrys%n_mol=n-1
                return
             end if
             exit
          end do

          !---- Search for three points (fractional coordinates) ----!
          !---- defining a Cartesian frame                       ----!
          do
             if (n_ini > n_end) exit
             line=adjustl(file_dat(n_ini))
             if (u_case(line(1:9)) == "XYZ_FRAME") then
                read(unit=line(10:),fmt=*,iostat=ierr) x1f,x2f,x3f
                if (ierr == 0) then
                   call get_euler_from_fract(x1f,x2f,x3f,molcrys%Cell%Cr_Orth_cel,phi,theta,chi,EuM, Code="D")
                   molcrys%mol(n)%orient(1)= phi
                   molcrys%mol(n)%orient(2)= theta
                   molcrys%mol(n)%orient(3)= chi
                   molcrys%mol(n)%xcentre= x3f
                   call Set_euler_matrix(molcrys%mol(n)%rot_type, phi,theta,chi,EuM)
                   molcrys%mol(n)%Euler=EuM
                   molcrys%mol(n)%is_EulerMat=.true.
                   molcrys%mol(n)%in_Xtal=.true.
                end if
                n_ini=n_ini+1
                exit
             else
                if (u_case(line(1:4)) =="MOLE") exit
                n_ini=n_ini+1
             end if
          end do

       end do

       return
    End Subroutine Readn_Set_XTal_CFL_Molec

    !!--++
    !!--++ Subroutine Readn_Set_XTal_CIF(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    Character(len=*),    optional,intent( in)  :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a CIF File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_CIF(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       Character(len=*),    optional,intent( in)  :: CFrame
       Integer,             optional,intent( in)  :: Nphase

       !---- Local Variables ----!
       character(len=132)                :: line
       character(len= 20)                :: Spp
       character(len=60), dimension(192) :: symm_car

       integer                   :: i, nauas, ndata, iph, n_ini,n_end,noper
       integer, parameter        :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph) :: ip

       real(kind=cp),dimension(6):: vet,vet2

       ip=nlines
       ip(1)=1

       !---- First determine if there is more than one structure ----!
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:5)) == "data_" .and. l_case(line(1:11)) /= "data_global" )  then
             n_ini=i
             ip(1)=i
             exit
          end if
       end do

       ndata=0
       do i=n_ini,nlines
          line=adjustl(file_dat(i))
          if (l_case(line(1:5)) == "data_")  then
             ndata=ndata+1
             if (ndata > maxph-1) then
                err_form=.true.
                ERR_Form_Mess=" => Too many phases in this file "
                return
             end if
             ip(ndata)=i   !Pointer to the number of the line starting a single phase
          end if
       end do

       iph=1
       if (present(nphase)) iph=nphase

       !---- Read Cell Parameters ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Cell(file_dat,n_ini,n_end,vet,vet2)
       if (err_form) return
       if(present(CFrame)) then
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,CFrame,vet2(1:3),vet2(4:6))
       else
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,"A",vet2(1:3),vet2(4:6))
       end if
       !---- Read Atoms Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Atom(file_dat,n_ini,n_end,nauas,A)
       if (err_form) return

       !---- SpaceGroup Information ----!
       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       call Read_Cif_Hm(file_dat,n_ini,n_end,Spp)

       n_ini=ip(iph)           !Updated values to handle non-conventional order
       n_end=ip(iph+1)
       if (len_trim(Spp) == 0) call Read_Cif_Hall(file_dat,n_ini,n_end,Spp)

       if (len_trim(Spp) == 0) then
          n_ini=ip(iph)           !Updated values to handle non-conventional order
          n_end=ip(iph+1)
          call Read_Cif_Symm(file_dat,n_ini,n_end,noper,symm_car)

          if (noper ==0) then
             err_form=.true.
             ERR_Form_Mess=" => No Space Group/No Symmetry information in this file "
             return
          else
             call Set_SpaceGroup("  ",SpG,symm_car,noper,"GEN")
          end if
       else
          call Set_SpaceGroup(Spp,SpG) !Construct the space group
       end if

       !---- Modify occupation factors and set multiplicity of atoms
       !---- in order to be in agreement with the definitions of Sfac in CrysFML
       !---- Convert Us to Betas and Uiso to Biso
       do i=1,A%natoms
          vet(1:3)=A%atom(i)%x
          A%atom(i)%Mult=Get_Multip_Pos(vet(1:3),SpG)
          A%atom(i)%Occ=A%atom(i)%Occ*real(A%atom(i)%Mult)/real(SpG%Multip)
          if(A%atom(i)%occ < epss) A%atom(i)%occ=real(A%atom(i)%Mult)/real(SpG%Multip)

          select case (A%atom(i)%thtype)
             case ("isotr")
                A%atom(i)%biso= A%atom(i)%ueq*78.95683521

             case ("aniso")
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) = Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"

             case default
                A%atom(i)%biso = A%atom(i)%ueq*78.95683521
                A%atom(i)%thtype = "isotr"
          end select
       end do

       return
    End Subroutine Readn_Set_XTal_CIF

    !!--++
    !!--++ Subroutine Readn_Set_XTal_PCR(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (atom_list_type),        intent(out)  :: A
    !!--++    character(len=*),    optional,intent(in)   :: CFrame
    !!--++    Integer,             optional,intent( in)  :: Nphase
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a PCR File
    !!--++
    !!--++ Update: 17/05/2010
    !!
    Subroutine Readn_Set_XTal_PCR(file_dat, nlines, Cell, Spg, A, CFrame, NPhase)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (atom_list_type),        intent(out)  :: A
       character(len=*),    optional,intent(in)   :: CFrame
       Integer,             optional,intent(in)   :: Nphase

       !---- Local Variables ----!
       logical                           :: multi,ask_phase,is_codewords
       character(len=132)                :: line
       character(len= 20)                :: Spp, label
       integer                           :: i,j, k,iv, nauas, ndata, iph, n_ini,n_end, nlong1
       integer, parameter                :: maxph=21  !Maximum number of phases "maxph-1"
       integer, dimension(maxph)         :: ip
       integer, dimension(30)            :: ivet

       real(kind=cp),dimension(30)       :: vet

       ip=nlines
       ip(1)=1

       !> Simple / Multi format
       multi=.false.
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (line(1:1) =='!' .or. line(1:1)==' ') cycle
          if (index(line,'NPATT ') <=0) cycle
          multi=.true.
       end do

       !> Number of Phases
       if (.not. multi) then
          do i=2,nlines
             line=adjustl(file_dat(i))
             if (line(1:1) =='!' .or. line(1:1)==' ') cycle
             call getnum(line,vet,ivet,iv)
             if (iv > 3) then
                iph=ivet(3)
                exit
             end if
          end do

       else
          do i=1,nlines
             line=adjustl(file_dat(i))
             if (line(1:4) /='!Nph') cycle

             line=adjustl(file_dat(i+1))
             call getnum(line,vet,ivet,iv)
             if (iv > 1) then
                iph=ivet(1)
                exit
             end if
          end do
       end if
       if (iph == 0) then
          err_form=.true.
          ERR_Form_Mess=" No Phase information was found in this PCR file. Please, check it! "
          return
       end if

       !> Locate where begin each Phase
       k=0
       ask_phase=.true.

       do i=1,nlines
          line=adjustl(file_dat(i))
          if (ask_phase) then
             if (index(line,'Data for PHASE') <= 0) cycle
          else
             if (line(1:1) /='!') then
                k=k+1
                ip(k)=i
                if (k == iph) exit

                ask_phase=.true.
             end if
             cycle
          end if
          ask_phase=.false.
       end do
       if (iph /= k) then
          err_form=.true.
          ERR_Form_Mess=" Locating Phases failed in this PCR. Please, check it!"
          return
       end if

       !> Select the Phase
       iph=1
       if (present(nphase)) iph=nphase
       n_ini=ip(iph)
       n_end=ip(iph+1)

       !---- Read Cell Parameters ----!
       do i=n_ini,n_end
          if (index(file_dat(i),'alpha') /=0 .and. index(file_dat(i),'gamma') /=0) then
             do j=i+1,n_end
                line=adjustl(file_dat(j))
                if (line(1:1) == '!' .or. line(1:1) == ' ') cycle
                iv=index(line,'#')
                if (iv > 1) line=line(1:iv-1)

                call getnum(line, vet, ivet,iv)
                if (iv /= 6) then
                   err_form=.true.
                   ERR_Form_Mess=" => Problems reading Cell Parameters on PCR file "
                   return
                end if
                if(present(CFrame)) then
                  call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,CFrame)
                else
                  call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell)
                end if
                exit
             end do
             exit
          end if
       end do

       !---- SpaceGroup Information ----!
       Spp=' '
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (line(1:1) == '!' .or. line(1:1)==' ') cycle
          if (index(file_dat(i),'<--Space') /=0) then
             j=index(file_dat(i),'<--Space')
             Spp=adjustl(file_dat(i)(1:j-1))
             if (len_trim(Spp) <= 0) then
                err_form=.true.
                ERR_Form_Mess=" => Problems reading Space group on PCR file "
                return
             end if
             call Set_SpaceGroup(Spp,SpG) !Construct the space group
             exit
          end if
       end do

       !---- Read Atoms Information ----!
       do i=n_ini,n_end
          line=adjustl(file_dat(i))
          if (line(1:4) /= '!Nat') cycle
          do j=i+1,n_end
             line=adjustl(file_dat(j))
             if (line(1:1) == '!' .or. line(1:1)==' ') cycle
             call getnum(line(1:5),vet,ivet,iv)
             ndata=ivet(1)
             exit
          end do
          exit
       end do

       if (ndata > 0) then
          call allocate_atom_list(ndata,A)

          is_codewords=.false.
          nauas=0

          do i=n_ini,n_end
             line=adjustl(file_dat(i))
             if (index(line,'!Atom') == 0 .or. index(line,'Typ') == 0) cycle

             do j=i+1,n_end
                line=adjustl(file_dat(j))
                if (line(1:1) == '!' .or. line(1:1)==' ') cycle
                if (is_codewords) then
                   is_codewords=.false.
                   cycle
                end if

                iv=index(line,'#')
                if (iv > 1) line=line(1:iv-1)

                nauas=nauas+1
                ! Atom Label
                call cutst(line,nlong1,label)
                A%atom(nauas)%lab=trim(label)

                ! Atom Type
                call cutst(line,nlong1,label)
                A%Atom(nauas)%chemsymb=U_case(label(1:1))//L_case(label(2:2))

                ! Atom Coordinates,Biso and Occ
                call getnum(line,vet,ivet,iv)
                if (iv < 5) then    !Line reading for the second time anisotropic temperature factors
                   nauas = nauas -1 !see below
                   is_codewords=.true.
                   cycle
                end if

                A%atom(nauas)%x=vet(1:3)
                A%atom(nauas)%Mult=Get_Multip_Pos(vet(1:3),SpG)
                A%atom(nauas)%biso=vet(4)
                A%atom(nauas)%occ=vet(5)
                A%atom(nauas)%thtype='isotr'
                A%atom(nauas)%Utype="beta"
                if (ivet(8) == 2) then    ! Anisotropic reading
                   A%atom(nauas)%thtype='aniso'
                   call getnum(file_dat(j+2),vet,ivet,iv)
                   A%atom(nauas)%u(1:6)=vet(1:6)
                end if
                is_codewords=.true.
                if (nauas == ndata) exit
             end do
             exit
          end do
       end if

       return
    End Subroutine Readn_Set_XTal_PCR

    !!--++
    !!--++ Subroutine Readn_Set_XTal_SHX(file_dat,nlines,Cell,SpG,A,CFrame)
    !!--++    character(len=*),dimension(:),intent(in)   :: file_dat
    !!--++    integer,                      intent(in)   :: nlines
    !!--++    Type (Crystal_Cell_Type),     intent(out)  :: Cell
    !!--++    Type (Space_Group_Type),      intent(out)  :: SpG
    !!--++    Type (Atom_list_type),        intent(out)  :: A
    !!--++    Character(len=*), optional,   intent(in)   :: CFrame
    !!--++
    !!--++ (Private)
    !!--++ Read and Set Crystal Information in a Shelx File
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_XTal_SHX(file_dat,nlines,Cell,SpG,A,CFrame)
       !---- Arguments ----!
       character(len=*),dimension(:),intent(in)   :: file_dat
       integer,                      intent(in)   :: nlines
       Type (Crystal_Cell_Type),     intent(out)  :: Cell
       Type (Space_Group_Type),      intent(out)  :: SpG
       Type (Atom_list_type),        intent(out)  :: A
       Character(len=*), optional,   intent(in)   :: CFrame

       !---- Local Variables ----!
       character(len=60), dimension(192) :: symm_car
       character(len=2),  dimension(15)  :: elem_atm
       integer                           :: i,n_ini, n_end, nl, noper
       integer                           :: n_elem_atm, n_fvar
       real(kind=cp), dimension(6)       :: vet,vet2
       real(kind=cp), dimension(10)      :: fvar

       n_ini=1
       n_end=nlines

       !---- CELL / ZERR ----!
       call Read_Shx_Cell(file_dat,n_ini,n_end,vet,vet2)
       if(present(CFrame)) then
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,CFrame,vet2(1:3),vet2(4:6))
       else
         call Set_Crystal_Cell(vet(1:3),vet(4:6),Cell,"A",vet2(1:3),vet2(4:6))
       end if

       !---- OBTAIN SPACE GROUP (LATT / SYMM) ----!
       call Read_Shx_Latt(file_dat,n_ini,n_end,nl)
       call Read_Shx_Symm(file_dat,n_ini,n_end,noper,symm_car)
       if (nl > 0) then
          noper=noper+1
          symm_car(noper)="-X,-Y,-Z"
       end if
       select case (abs(nl))
          case (2) ! I
             noper=noper+1
             symm_car(noper)="X+1/2,Y+1/2,Z+1/2"
          case (3) ! Rom, Hex
             noper=noper+1
             symm_car(noper)="X+2/3,Y+1/3,Z+1/3"
             noper=noper+1
             symm_car(noper)="X+1/3,Y+2/3,Z+2/3"
          case (4) ! F
             noper=noper+1
             symm_car(noper)="X,Y+1/2,Z+1/2"
          case (5) ! A
             noper=noper+1
             symm_car(noper)="X,Y+1/2,Z+1/2"
             noper=noper+1
             symm_car(noper)="X+1/2,Y,Z+1/2"
             noper=noper+1
             symm_car(noper)="X+1/2,Y+1/2,Z"
          case (6) ! B
             noper=noper+1
             symm_car(noper)="X+1/2,Y,Z+1/2"
          case (7) ! C
             noper=noper+1
             symm_car(noper)="X+1/2,Y+1/2,Z"
       end select ! nl
       call set_spacegroup(" ",SPG,symm_car,noper,"gen")

       !---- ATOMS ----!
       call Read_Shx_Cont(file_dat,n_ini,n_end,n_elem_atm,elem_atm)
       call Read_Shx_Fvar(file_dat,n_ini,n_end,n_fvar,fvar)
       call Read_Shx_Atom(file_dat,n_ini,n_end,n_fvar,fvar,elem_atm,cell,A)
       if (err_form) return

       !---- Convert Us to Betas and Uiso to Biso
       do i=1,A%natoms
          vet(1:3)=A%atom(i)%x
          A%atom(i)%Mult=Get_Multip_Pos(vet(1:3),SpG)

          select case (A%atom(i)%thtype)
             case ("isotr")
                A%atom(i)%biso= A%atom(i)%ueq*78.95683521

             case ("aniso")
                A%atom(i)%ueq=U_Equiv(cell,a%atom(i)%u(1:6))  ! Uequi
                A%atom(i)%biso= A%atom(i)%ueq*78.95683521
                select case (A%atom(i)%Utype)
                   case ("u_ij")
                      A%atom(i)%u(1:6) =  Convert_U_Betas(A%atom(i)%u(1:6),Cell)
                   case ("b_ij")
                      A%atom(i)%u(1:6) = Convert_B_Betas(A%atom(i)%u(1:6),Cell)
                end select
                A%atom(i)%Utype="beta"

             case default
                A%atom(i)%ueq=0.05
                A%atom(i)%biso = A%atom(i)%ueq*78.95683521
                A%atom(i)%thtype = "isotr"
          end select
       end do

       return
    End Subroutine Readn_Set_XTal_SHX

    !!--++
    !!--++ Subroutine Readn_Set_Xtal_Structure_Molcr(filenam,Molcrys,Mode,Iphase, Job_Info, file_list,CFrame)
    !!--++    character(len=*),              intent( in)     :: filenam  ! In -> Name of the file
    !!--++    Type (Molecular_Crystal_Type), intent(out)     :: Molcrys  ! Molecular crytal
    !!--++    Character(len=*),    optional, intent( in)     :: Mode     ! In -> if Mode="CIF" filenam
    !!--++                                                                       is of CIF type format
    !!--++    Integer,             optional, intent( in)     :: Iphase   ! Number of the phase.
    !!--++    Type(Job_Info_type), optional, intent(out)     :: Job_Info ! Diffaction conditions
    !!--++    Type(file_list_type),optional, intent(in out)  :: file_list! Complete file to be used by
    !!--++                                                              the calling program or other procedures
    !!--++    Character(len=*),    optional, intent(in)      :: CFrame
    !!--++    Overloaded
    !!--++    Subroutine to read and input file and construct the crystal structure
    !!--++    in terms of the ofjects Cell, SpG and A. The optional argument Iphase is an integer
    !!--++    telling to the program to read the phase number Iphase in the case of the presence
    !!--++    of more than one phase. If absent only the first phase is read.
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_Xtal_Structure_Molcr(filenam,Molcrys,Mode,Iphase,Job_Info,file_list,CFrame)
       !---- Arguments ----!
       character(len=*),              intent( in)     :: filenam
       Type (Molecular_Crystal_Type), intent(out)     :: Molcrys
       Character(len=*),     optional,intent( in)     :: Mode
       Integer,              optional,intent( in)     :: Iphase
       Type(Job_Info_type),  optional,intent(out)     :: Job_Info
       Type(file_list_type), optional,intent(in out)  :: file_list
       Character(len=*),     optional,intent(in)      :: CFrame
       !---- Local variables -----!
       Type (Atom_list_type)                         :: A
       character(len=132), allocatable, dimension(:) :: file_dat
       character(len=3)                              :: modec
       integer                                       :: i,nlines


       call init_err_form()

       nlines=0
       if (present(file_list)) nlines=file_list%nlines

       !---- Number of Lines in the input file ----!
       if(nlines == 0) then
           call Number_Lines(trim(filenam), nlines)
           if (nlines==0) then
              err_form=.true.
              ERR_Form_Mess="The file "//trim(filenam)//" contains nothing"
              return
           else
              if (allocated(file_dat)) deallocate( file_dat)
              allocate( file_dat(nlines))
              call reading_Lines(trim(filenam),nlines,file_dat)
           end if
           if (present(file_list)) then
              file_list%nlines=nlines
              if (allocated(file_list%line)) deallocate(file_list%line)
              allocate(file_list%line(nlines))
              file_list%line=file_dat
           end if
       else
           if (allocated(file_dat)) deallocate( file_dat)
           allocate( file_dat(nlines))
           file_dat=file_list%line
       end if


       !---- Define the type of file: CIF, CFL, RES,... ----!
       modec=" "
       if (present(mode)) modec=l_case(mode(1:3))

       select case(modec)
           case("cif")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                 end if
              end if

           case("pcr")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                 end if
              end if

           case("shx")
              if(present(CFrame)) then
                call readn_set_xtal_shx(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
              else
                call readn_set_xtal_shx(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
              end if
           case default
              !---- CFL Format ----!
              if (present(Job_Info)) then
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,NPhase=IPhase,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,NPhase=IPhase,Job_Info=Job_Info)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,Job_Info=Job_Info)
                    end if
                 end if
              else
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame,NPhase=IPhase)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,NPhase=IPhase)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A,CFrame)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,molcrys%Cell,molcrys%Spg,A)
                    end if
                 end if
              end if
              !---- Reading molecules ----!
              if (present(iphase)) then
                call readn_set_xtal_cfl_molec(file_dat,nlines,molcrys,NPhase=IPhase)
              else
                call readn_set_xtal_cfl_molec(file_dat,nlines,molcrys)
              end if

       end select
       if (err_form) return

       !---- Passing from Atom_List_Type -> Molcrys ----!
       molcrys%n_free=A%natoms
       if (A%natoms > 0) then
          if (allocated(molcrys%Atm)) deallocate(molcrys%Atm)
          allocate(molcrys%Atm(A%natoms))
          molcrys%Atm=A%Atom
       end if

       call deallocate_atom_list(A)

       !---- Testing if Xtal was defined ----!
       if (all(molcrys%cell%cell > 0.0)) then
          do i=1,molcrys%n_mol
             if (.not. molcrys%mol(i)%in_xtal) then
                 molcrys%mol(i)%in_xtal=.true.
             end if
          end do
       end if

       return
    End Subroutine Readn_Set_Xtal_Structure_Molcr

    !!--++
    !!--++ Subroutine Readn_Set_Xtal_Structure_Split(filenam,Cell,SpG,A,Mode,Iphase,Job_Type,File_List,CFrame)
    !!--++    character(len=*),              intent( in)     :: filenam  ! In -> Name of the file
    !!--++    Type (Crystal_Cell_Type),      intent(out)     :: Cell     ! Out -> Cell object
    !!--++    Type (Space_Group_Type),       intent(out)     :: SpG      ! Out -> Space Group object
    !!--++    Type (atom_list_type),         intent(out)     :: A        ! Out -> Atom_List object
    !!--++    Character(len=*),    optional, intent( in)     :: Mode     ! In -> if Mode="CIF" filenam
    !!--++                                                                       is of CIF type format
    !!--++    Integer,             optional, intent( in)     :: Iphase   ! Number of the phase.
    !!--++    Type(Job_Info_type), optional, intent(out)     :: Job_Info ! Diffaction conditions
    !!--++    Type(file_list_type),optional, intent(in out)  :: file_list! Complete file to be used by
    !!--++                                                                 the calling program or other procedures
    !!--++    Character(len=*),    optional, intent( in)     :: CFrame   !Cartesian Frame
    !!--++
    !!--++    Overloaded
    !!--++    Subroutine to read and input file and construct the crystal structure
    !!--++    in terms of the ofjects Cell, SpG and A. The optional argument Iphase is an integer
    !!--++    telling to the program to read the phase number Iphase in the case of the presence
    !!--++    of more than one phase. If absent only the first phase is read.
    !!--++
    !!--++ Update: April - 2005
    !!
    Subroutine Readn_Set_Xtal_Structure_Split(filenam,Cell,SpG,A,Mode,Iphase,Job_Info,file_list,CFrame)
       !---- Arguments ----!
       character(len=*),             intent( in)     :: filenam
       Type (Crystal_Cell_Type),     intent(out)     :: Cell
       Type (Space_Group_Type),      intent(out)     :: SpG
       Type (atom_list_type),        intent(out)     :: A
       Character(len=*),    optional,intent( in)     :: Mode
       Integer,             optional,intent( in)     :: Iphase
       Type(Job_Info_type), optional,intent(out)     :: Job_Info
       Type(file_list_type),optional,intent(in out)  :: file_list
       Character(len=*),    optional,intent( in)     :: CFrame

       !---- Local variables -----!
       character(len=132), allocatable, dimension(:) :: file_dat
       character(len=3)                              :: modec
       integer                                       :: nlines

       call init_err_form()

       nlines=0
       if (present(file_list)) nlines=file_list%nlines

       !---- Number of Lines in the input file ----!
       if(nlines == 0) then
           call Number_Lines(trim(filenam), nlines)
           if (nlines==0) then
              err_form=.true.
              ERR_Form_Mess="The file "//trim(filenam)//" contains nothing"
              return
           else
              if (allocated(file_dat)) deallocate( file_dat)
              allocate( file_dat(nlines))
              call reading_Lines(trim(filenam),nlines,file_dat)
           end if
           if (present(file_list)) then
              file_list%nlines=nlines
              if (allocated(file_list%line)) deallocate(file_list%line)
              allocate(file_list%line(nlines))
              file_list%line=file_dat
           end if
       else
           if (allocated(file_dat)) deallocate( file_dat)
           allocate( file_dat(nlines))
           file_dat=file_list%line
       end if

       !---- Define the type of file: CIF, CFL, RES,... ----!
       modec=" "
       if (present(mode)) modec=l_case(mode(1:3))

       select case(modec)
           case("cif")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,Cell,Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,Cell,Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_cif(file_dat,nlines,Cell,Spg,A,CFrame)
                 else
                   call readn_set_xtal_cif(file_dat,nlines,Cell,Spg,A)
                 end if
              end if

           case("pcr")
              if (present(iphase)) then
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg, A,CFrame,NPhase=IPhase)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg, A,NPhase=IPhase)
                 end if
              else
                 if(present(CFrame)) then
                   call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg,A,CFrame)
                 else
                   call readn_set_xtal_pcr(file_dat,nlines,Cell,Spg,A)
                 end if
              end if

           case("shx")
              if(present(CFrame)) then
                call readn_set_xtal_shx(file_dat,nlines,Cell,Spg,A,CFrame)
              else
                call readn_set_xtal_shx(file_dat,nlines,Cell,Spg,A)
              end if
           case default
              !---- CFL Format ----!
              if (present(Job_Info)) then
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,CFrame,NPhase=IPhase,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,NPhase=IPhase,Job_Info=Job_Info)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,CFrame,Job_Info=Job_Info)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,Job_Info=Job_Info)
                    end if
                 end if
              else
                 if (present(iphase)) then
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,CFrame,NPhase=IPhase)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,NPhase=IPhase)
                    end if
                 else
                    if(present(CFrame)) then
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A,CFrame)
                    else
                      call readn_set_xtal_cfl(file_dat,nlines,Cell,Spg,A)
                    end if
                 end if
              end if

       end select

       return
    End Subroutine Readn_Set_Xtal_Structure_Split

    !!----
    !!---- Subroutine Write_Cif_Powder_Profile(Filename,Code)
    !!----    character(len=*), intent(in) :: filename     !  In -> Name of File
    !!----    integer,     intent(in)      :: code         !  In -> 0 Shelxs-Patterson
    !!----                                                          1 Shelxs-Direct Methods
    !!----                                                          2 Shelxl-Refinement
    !!----
    !!----    Write a Cif Powder Profile file
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Cif_Powder_Profile(filename,code)
       !---- Arguments ----!
       character(len=*), intent(in) :: filename
       integer,          intent(in) :: code

       !---- Local Variables ----!
       logical :: info

       integer :: iunit !,nlong

       !---- Inicializacion de variables ----!
       info=.false.
       iunit=0

       !---- Esta abierto este Fichero? ----!
       inquire(file=filename,opened=info)
       if (info) then
          inquire(file=filename,number=iunit)
          close(unit=iunit)
       end if

       !---- Escritura ----!
       if (iunit==0) iunit=61
       open(unit=iunit,file=filename,status="unknown",action="write")
       rewind(unit=iunit)

       !---- Head ----!
       write(unit=iunit,fmt="(a)") "data_profile"

       write(unit=iunit,fmt="(a)") " "
       if (code == 0) then
          write(unit=iunit,fmt="(a)")     "_pd_block_id      ?"
       else
          write(unit=iunit,fmt="(a,i3)")  "_pd_block_id       ",code
       end if

       !---- Profile ----!
       write(unit=iunit,fmt="(a)") " "

       write(unit=iunit,fmt="(a)") "loop_"
       write(unit=iunit,fmt="(a)") "_pd_proc_point_id"
       write(unit=iunit,fmt="(a)") "_pd_proc_2theta_corrected             # one of "
       write(unit=iunit,fmt="(a)") "_pd_proc_energy_incident              # these "
       write(unit=iunit,fmt="(a)") "_pd_proc_d_spacing                    # three"
       write(unit=iunit,fmt="(a)") "_pd_proc_intensity_net"
       write(unit=iunit,fmt="(a)") "_pd_calc_intensity_net "
       write(unit=iunit,fmt="(a)") "_pd_proc_ls_weight      "
       write(unit=iunit,fmt="(a)") "?     ?     ?     ?     ?     ?     ?"

       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "# The following lines are used to test the character set of files sent by     "
       write(unit=iunit,fmt="(a)") "# network email or other means. They are not part of the CIF data set.        "
       write(unit=iunit,fmt="(a)") "# abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789              "
       write(unit=iunit,fmt="(a)") "# !@#$%^&*()_+{}:"//""""//"~<>?|\-=[];'`,./ "

       close (unit=iunit)

       return
    End Subroutine Write_Cif_Powder_Profile

    !!----
    !!---- Subroutine Write_Cif_Template(filename,type_data,code,cell,SpG,A)
    !!----    character(len=*),        intent(in) :: filename   !  In -> Filename
    !!----    integer,                 intent(in) :: type_data  !  In -> 0: Single Crystal, 1: Powder Data, 2:Only structural data
    !!----    character(len=*),        intent(in) :: code       !  In -> Code name for the data set
    !!----    Type (Crystal_Cell_Type),intent(in) :: Cell       ! Cell type to be output
    !!----    Type (Space_Group_Type), intent(in) :: SpG        ! Space group type to be output
    !!----    Type (Atom_List_Type),   intent(in) :: A          ! Atom list type to be output
    !!----
    !!----    Write a Cif File
    !!----
    !!---- Updated: February - 2005, January 2015
    !!
    Subroutine Write_Cif_Template(filename,type_data,code,cell,SpG,A)
       !---- Arguments ----!
       character(len=*),        intent(in) :: filename
       integer,                 intent(in) :: type_data
       character(len=*),        intent(in) :: code
       Type (Crystal_Cell_Type),intent(in) :: Cell
       Type (Space_Group_Type), intent(in) :: SpG
       Type (Atom_List_Type),   intent(in) :: A

       !---- Local Variables ----!
       logical                           :: info,aniso
       character(len=132)                :: line
       character(len=1), parameter       :: qmark='?'
       character(len=30)                 :: comm,adptyp
       character(len=30),dimension(6)    :: text
       real(kind=cp)                     :: u,su, ocf
       real(kind=cp), dimension(6)       :: Ua,sua,aux
       real(kind=cp), dimension(A%natoms):: occup,soccup
       integer                           :: iunit,i, j

       !---- Initialization of variables ----!
       info=.false.
       iunit=0

       !---- Is this file opened? ----!
       inquire(file=filename,opened=info)
       if (info) then
          inquire(file=filename,number=iunit)
          close(unit=iunit)
       end if

       !---- Writing ----!
       if (iunit==0) iunit=61
       open(unit=iunit,file=filename,status="unknown",action="write")
       rewind(unit=iunit)

       !---- Head Information ----!
       if(type_data == 0) then
           write(unit=iunit,fmt="(a)") "##############################################################################"
           write(unit=iunit,fmt="(a)") "###    CIF submission form for molecular structure report (Acta Cryst. C)  ###"
           write(unit=iunit,fmt="(a)") "##############################################################################"
           write(unit=iunit,fmt="(a)") " "
           write(unit=iunit,fmt="(a)") "#============================================================================="
           write(unit=iunit,fmt="(a)") "data_global"
           write(unit=iunit,fmt="(a)") "#============================================================================="
           write(unit=iunit,fmt="(a)") " "
          else if(type_data > 1) then
           write(unit=iunit,fmt="(a)") "##################################################################"
           write(unit=iunit,fmt="(a)") "###    CIF file from CrysFML, contains only structural data    ###"
           write(unit=iunit,fmt="(a)") "##################################################################"

       end if

       !---- Processing Summary ----!
       if(type_data < 2) then
         write(unit=iunit,fmt="(a)") "# PROCESSING SUMMARY (IUCr Office Use Only)"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_journal_data_validation_number      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_journal_date_recd_electronic        ?"
         write(unit=iunit,fmt="(a)") "_journal_date_to_coeditor            ?"
         write(unit=iunit,fmt="(a)") "_journal_date_from_coeditor          ?"
         write(unit=iunit,fmt="(a)") "_journal_date_accepted               ?"
         write(unit=iunit,fmt="(a)") "_journal_date_printers_first         ?"
         write(unit=iunit,fmt="(a)") "_journal_date_printers_final         ?"
         write(unit=iunit,fmt="(a)") "_journal_date_proofs_out             ?"
         write(unit=iunit,fmt="(a)") "_journal_date_proofs_in              ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_name               ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_code               ?"
         write(unit=iunit,fmt="(a)") "_journal_coeditor_notes"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_journal_techeditor_code             ?"
         write(unit=iunit,fmt="(a)") "_journal_techeditor_notes"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_journal_coden_ASTM                  ?"
         write(unit=iunit,fmt="(a)") "_journal_name_full                   ?"
         write(unit=iunit,fmt="(a)") "_journal_year                        ?"
         write(unit=iunit,fmt="(a)") "_journal_volume                      ?"
         write(unit=iunit,fmt="(a)") "_journal_issue                       ?"
         write(unit=iunit,fmt="(a)") "_journal_page_first                  ?"
         write(unit=iunit,fmt="(a)") "_journal_page_last                   ?"
         write(unit=iunit,fmt="(a)") "_journal_paper_category              ?"
         write(unit=iunit,fmt="(a)") "_journal_suppl_publ_number           ?"
         write(unit=iunit,fmt="(a)") "_journal_suppl_publ_pages            ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Submission details ----!
         write(unit=iunit,fmt="(a)") "# 1. SUBMISSION DETAILS"
         write(unit=iunit,fmt="(a)") " "

         write(unit=iunit,fmt="(a)") "_publ_contact_author_name            ?   # Name of author for correspondence"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_address             # Address of author for correspondence"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_email           ?"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_fax             ?"
         write(unit=iunit,fmt="(a)") "_publ_contact_author_phone           ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_contact_letter"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_requested_journal              ?"
         write(unit=iunit,fmt="(a)") "_publ_requested_coeditor_name        ?"
         write(unit=iunit,fmt="(a)") "_publ_requested_category             ?   # Acta C: one of CI/CM/CO/FI/FM/FO"

         write(unit=iunit,fmt="(a)") "#=============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Title  and Author List ----!
         write(unit=iunit,fmt="(a)") "# 3. TITLE AND AUTHOR LIST"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_section_title"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_title_footnote"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The loop structure below should contain the names and addresses of all "
         write(unit=iunit,fmt="(a)") "# authors, in the required order of publication. Repeat as necessary."

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _publ_author_name"
         write(unit=iunit,fmt="(a)") "    _publ_author_footnote"
         write(unit=iunit,fmt="(a)") "    _publ_author_address"
         write(unit=iunit,fmt="(a)") "?                                   #<--'Last name, first name' "
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Text ----!
         write(unit=iunit,fmt="(a)") "# 4. TEXT"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_publ_section_synopsis"
         write(unit=iunit,fmt="(a)") ";  ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_abstract"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";          "
         write(unit=iunit,fmt="(a)") "_publ_section_comment"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_exptl_prep      # Details of the preparation of the sample(s)"
         write(unit=iunit,fmt="(a)") "                              # should be given here. "
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_exptl_refinement"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_references"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_figure_captions"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_publ_section_acknowledgements"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Identifier ----!
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") "# If more than one structure is reported, the remaining sections should be "
         write(unit=iunit,fmt="(a)") "# completed per structure. For each data set, replace the '?' in the"
         write(unit=iunit,fmt="(a)") "# data_? line below by a unique identifier."
       end if !type_data < 2

       write(unit=iunit,fmt="(a)") " "
       if (len_trim(code) == 0) then
          write(unit=iunit,fmt="(a)") "data_?"
       else
          write(unit=iunit,fmt="(a)") "data_"//code(1:len_trim(code))
       end if
       write(unit=iunit,fmt="(a)") " "
       if(type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Chemical Data ----!
         write(unit=iunit,fmt="(a)") "# 5. CHEMICAL DATA"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_chemical_name_systematic"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
         write(unit=iunit,fmt="(a)") "_chemical_name_common             ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_moiety          ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_structural      ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_analytical      ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_iupac           ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_sum             ?"
         write(unit=iunit,fmt="(a)") "_chemical_formula_weight          ?"
         write(unit=iunit,fmt="(a)") "_chemical_melting_point           ?"
         write(unit=iunit,fmt="(a)") "_chemical_compound_source         ?       # for minerals and "
         write(unit=iunit,fmt="(a)") "                                          # natural products"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "loop_"
         write(unit=iunit,fmt="(a)") "    _atom_type_symbol               "
         write(unit=iunit,fmt="(a)") "    _atom_type_description          "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_dispersion_real "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_dispersion_imag "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_source          "
         write(unit=iunit,fmt="(a)") "    _atom_type_scat_length_neutron       # include if applicable"
         write(unit=iunit,fmt="(a)") "    ?    ?    ?    ?    ?      ?    "

       end if !type_data < 2
       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "#============================================================================="
       write(unit=iunit,fmt="(a)") " "
       !---- Crystal Data ----!
       select case (type_data)
          case (0,2) ! Single Crystal or structural data only
             write(unit=iunit,fmt="(a)") "# 6. CRYSTAL DATA"
          case (1) ! Powder Data + Crystal Data
             write(unit=iunit,fmt="(a)") "# 6. POWDER SPECIMEN AND CRYSTAL DATA"
       end select

       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "_symmetry_cell_setting               ?"
       line=SpG%SPG_Symb
       write(unit=iunit,fmt="(a)") "_symmetry_space_group_name_H-M       '"//trim(line)//"'"
       line=SpG%Hall
       write(unit=iunit,fmt="(a)") "_symmetry_space_group_name_Hall      '"//trim(line)//"'"

       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "loop_"
       write(unit=iunit,fmt="(a)") "    _symmetry_equiv_pos_as_xyz   #<--must include 'x,y,z'"
       do i=1,SpG%multip
          line="'"//trim(SpG%SymopSymb(i))//"'"
          write(iunit,'(a)') trim(line)
       end do

       write(unit=iunit,fmt="(a)") " "
       do i=1,3
          call setnum_std(Cell%cell(i),Cell%cell_std(i),text(i))
       end do
       do i=1,3
          call setnum_std(Cell%ang(i),Cell%ang_std(i),text(i+3))
       end do
       write(iunit,'(a)') "_cell_length_a                       "//trim(adjustl(text(1)))
       write(iunit,'(a)') "_cell_length_b                       "//trim(adjustl(text(2)))
       write(iunit,'(a)') "_cell_length_c                       "//trim(adjustl(text(3)))
       write(iunit,'(a)') "_cell_angle_alpha                    "//trim(adjustl(text(4)))
       write(iunit,'(a)') "_cell_angle_beta                     "//trim(adjustl(text(5)))
       write(iunit,'(a)') "_cell_angle_gamma                    "//trim(adjustl(text(6)))

       write(unit=iunit,fmt="(a,f14.4)") "_cell_volume                   ",Cell%CellVol
       if(type_data < 2) then
         write(unit=iunit,fmt="(a)") "_cell_formula_units_Z                ?"
         write(unit=iunit,fmt="(a)") "_cell_measurement_temperature        ?"
         write(unit=iunit,fmt="(a)") "_cell_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"
       end if

       select case (type_data)
          case (0) ! Single Crystal
             write(unit=iunit,fmt="(a)") "_cell_measurement_reflns_used        ?"
             write(unit=iunit,fmt="(a)") "_cell_measurement_theta_min          ?"
             write(unit=iunit,fmt="(a)") "_cell_measurement_theta_max          ?"

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "_exptl_crystal_description           ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_colour                ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_size_max              ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_size_mid              ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_size_min              ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_size_rad              ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_density_diffrn        ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_density_meas          ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_density_method        ?"
             write(unit=iunit,fmt="(a)") "_exptl_crystal_F_000                 ?"

          case (1) ! Powder Data
             write(unit=iunit,fmt="(a)") "# The next three fields give the specimen dimensions in mm.  The equatorial"
             write(unit=iunit,fmt="(a)") "# plane contains the incident and diffracted beam."

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "_pd_spec_size_axial               ?       # perpendicular to "
             write(unit=iunit,fmt="(a)") "                                          # equatorial plane"

             write(unit=iunit,fmt="(a)") "_pd_spec_size_equat               ?       # parallel to "
             write(unit=iunit,fmt="(a)") "                                          # scattering vector"
             write(unit=iunit,fmt="(a)") "                                          # in transmission"
             write(unit=iunit,fmt="(a)") "_pd_spec_size_thick               ?       # parallel to "
             write(unit=iunit,fmt="(a)") "                                          # scattering vector"
             write(unit=iunit,fmt="(a)") "                                          # in reflection"

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "# The next five fields are character fields that describe the specimen."

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "_pd_spec_mounting                         # This field should be"
             write(unit=iunit,fmt="(a)") "                                          # used to give details of the "
             write(unit=iunit,fmt="(a)") "                                          # container."
             write(unit=iunit,fmt="(a)") "; ?"
             write(unit=iunit,fmt="(a)") ";"
             write(unit=iunit,fmt="(a)") "_pd_spec_mount_mode               ?       # options are 'reflection'"
             write(unit=iunit,fmt="(a)") "                                          # or 'transmission'"
             write(unit=iunit,fmt="(a)") "_pd_spec_shape                    ?       # options are 'cylinder' "
             write(unit=iunit,fmt="(a)") "                                          # 'flat_sheet' or 'irregular'"
             write(unit=iunit,fmt="(a)") "_pd_char_particle_morphology      ?"
             write(unit=iunit,fmt="(a)") "_pd_char_colour                   ?       # use ICDD colour descriptions"

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "# The following three fields describe the preparation of the specimen."
             write(unit=iunit,fmt="(a)") "# The cooling rate is in K/min.  The pressure at which the sample was "
             write(unit=iunit,fmt="(a)") "# prepared is in kPa.  The temperature of preparation is in K.        "

             write(unit=iunit,fmt="(a)") " "
             write(unit=iunit,fmt="(a)") "_pd_prep_cool_rate                ?"
             write(unit=iunit,fmt="(a)") "_pd_prep_pressure                 ?"
             write(unit=iunit,fmt="(a)") "_pd_prep_temperature              ?"
       end select
       if(type_data < 2) then
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The next four fields are normally only needed for transmission experiments."
         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_coefficient_mu        ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_type       ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_process_details       ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_T_min      ?"
         write(unit=iunit,fmt="(a)") "_exptl_absorpt_correction_T_max      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Experimental Data ----!
         write(unit=iunit,fmt="(a)") "# 7. EXPERIMENTAL DATA"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_exptl_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         if (type_data == 1) then
            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# The following item is used to identify the equipment used to record "
            write(unit=iunit,fmt="(a)") "# the powder pattern when the diffractogram was measured at a laboratory "
            write(unit=iunit,fmt="(a)") "# other than the authors' home institution, e.g. when neutron or synchrotron"
            write(unit=iunit,fmt="(a)") "# radiation is used."

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_instr_location"
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"
            write(unit=iunit,fmt="(a)") "_pd_calibration_special_details           # description of the method used"
            write(unit=iunit,fmt="(a)") "                                          # to calibrate the instrument"
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"
         end if

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_diffrn_ambient_temperature          ?"
         write(unit=iunit,fmt="(a)") "_diffrn_radiation_type               ?"
         write(unit=iunit,fmt="(a)") "_diffrn_radiation_wavelength         ?"
         write(unit=iunit,fmt="(a)") "_diffrn_radiation_source             ?"
         write(unit=iunit,fmt="(a)") "_diffrn_source                       ?"
         write(unit=iunit,fmt="(a)") "_diffrn_source_target                ?"
         write(unit=iunit,fmt="(a)") "_diffrn_source_type                  ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_diffrn_radiation_monochromator      ?"
         write(unit=iunit,fmt="(a)") "_diffrn_measurement_device_type      ?"
         write(unit=iunit,fmt="(a)") "_diffrn_measurement_method           ?"
         write(unit=iunit,fmt="(a)") "_diffrn_detector_area_resol_mean     ?   # Not in version 2.0.1"
         write(unit=iunit,fmt="(a)") "_diffrn_detector                     ?"
         write(unit=iunit,fmt="(a)") "_diffrn_detector_type                ?   # make or model of detector"
         if (type_data == 1) then
            write(unit=iunit,fmt="(a)") "_pd_meas_scan_method                 ?   # options are 'step', 'cont',"
            write(unit=iunit,fmt="(a)") "                                         # 'tof', 'fixed' or"
            write(unit=iunit,fmt="(a)") "                                         # 'disp' (= dispersive)"
            write(unit=iunit,fmt="(a)") "_pd_meas_special_details"
            write(unit=iunit,fmt="(a)") ";  ?"
            write(unit=iunit,fmt="(a)") ";"
         end if

         select case (type_data)
            case (0)
               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_number                ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_av_R_equivalents      ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_av_sigmaI/netI        ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_min             ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_max             ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_theta_full            ?   # Not in version 2.0.1"
               write(unit=iunit,fmt="(a)") "_diffrn_measured_fraction_theta_max  ?   # Not in version 2.0.1"
               write(unit=iunit,fmt="(a)") "_diffrn_measured_fraction_theta_full ?   # Not in version 2.0.1"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_h_min           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_h_max           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_k_min           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_k_max           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_l_min           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_limit_l_max           ?"
               write(unit=iunit,fmt="(a)") "_diffrn_reflns_reduction_process     ?"

               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_diffrn_standards_number             ?"
               write(unit=iunit,fmt="(a)") "_diffrn_standards_interval_count     ?"
               write(unit=iunit,fmt="(a)") "_diffrn_standards_interval_time      ?"
               write(unit=iunit,fmt="(a)") "_diffrn_standards_decay_%            ?"
               write(unit=iunit,fmt="(a)") "loop_"
               write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_h"
               write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_k"
               write(unit=iunit,fmt="(a)") "    _diffrn_standard_refln_index_l"
               write(unit=iunit,fmt="(a)") "?   ?   ?"

            case (1)
               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "#  The following four items give details of the measured (not processed)"
               write(unit=iunit,fmt="(a)") "#  powder pattern.  Angles are in degrees."

               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_pd_meas_number_of_points         ?"
               write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_min         ?"
               write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_max         ?"
               write(unit=iunit,fmt="(a)") "_pd_meas_2theta_range_inc         ?"

               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "# The following three items are used for time-of-flight measurements only."

               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_pd_instr_dist_src/spec           ?"
               write(unit=iunit,fmt="(a)") "_pd_instr_dist_spec/detc          ?"
               write(unit=iunit,fmt="(a)") "_pd_meas_2theta_fixed             ?"

         end select

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "#============================================================================="
         write(unit=iunit,fmt="(a)") " "

         !---- Refinement Data ----!
         write(unit=iunit,fmt="(a)") "# 8. REFINEMENT DATA"

         write(unit=iunit,fmt="(a)") " "

         write(unit=iunit,fmt="(a)") "_refine_special_details"
         write(unit=iunit,fmt="(a)") "; ?"
         write(unit=iunit,fmt="(a)") ";"

         if (type_data == 1) then
            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# Use the next field to give any special details about the fitting of the"
            write(unit=iunit,fmt="(a)") "# powder pattern."

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "_pd_proc_ls_special_details"
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"

            write(unit=iunit,fmt="(a)") " "
            write(unit=iunit,fmt="(a)") "# The next three items are given as text."
            write(unit=iunit,fmt="(a)") " "

            write(unit=iunit,fmt="(a)") "_pd_proc_ls_profile_function      ?"
            write(unit=iunit,fmt="(a)") "_pd_proc_ls_background_function   ?"
            write(unit=iunit,fmt="(a)") "_pd_proc_ls_pref_orient_corr"
            write(unit=iunit,fmt="(a)") "; ?"
            write(unit=iunit,fmt="(a)") ";"
         end if

         select case (type_data)
            case (0)
               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_reflns_number_total                 ?"
               write(unit=iunit,fmt="(a)") "_reflns_number_gt                    ?  # Not in version 2.0.1"
               write(unit=iunit,fmt="(a)") "_reflns_threshold_expression         ?  # Not in version 2.0.1"

            case (1)
               write(unit=iunit,fmt="(a)") " "
               write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_R_factor         ?"
               write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_wR_factor        ?"
               write(unit=iunit,fmt="(a)") "_pd_proc_ls_prof_wR_expected      ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# The following four items apply to angular dispersive measurements."
              write(unit=iunit,fmt="(a)") "# 2theta minimum, maximum and increment (in degrees) are for the "
              write(unit=iunit,fmt="(a)") "# intensities used in the refinement."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_min         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_max         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_2theta_range_inc         ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_wavelength               ?"

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_block_diffractogram_id        ?  # The id used for the block containing"
              write(unit=iunit,fmt="(a)") "                                     # the powder pattern profile (section 11)."

              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "# Give appropriate details in the next two text fields."
              write(unit=iunit,fmt="(a)") " "
              write(unit=iunit,fmt="(a)") "_pd_proc_info_excluded_regions    ?"
              write(unit=iunit,fmt="(a)") "_pd_proc_info_data_reduction      ?"

         end select

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "_refine_ls_structure_factor_coef     ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_matrix_type               ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_R_I_factor                ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_R_Fsqd_factor             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_R_factor_all              ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_R_factor_gt               ?   # Not in version 2.0.1"
         write(unit=iunit,fmt="(a)") "_refine_ls_wR_factor_all             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_wR_factor_ref             ?   # Not in version 2.0.1"
         write(unit=iunit,fmt="(a)") "_refine_ls_goodness_of_fit_all       ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_goodness_of_fit_ref       ?   # Not in version 2.0.1"
         write(unit=iunit,fmt="(a)") "_refine_ls_restrained_S_all          ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_restrained_S_obs          ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_number_reflns             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_number_parameters         ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_number_restraints         ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_number_constraints        ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_hydrogen_treatment        ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_weighting_scheme          ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_weighting_details         ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_shift/su_max              ?   # Not in version 2.0.1"
         write(unit=iunit,fmt="(a)") "_refine_ls_shift/su_mean             ?   # Not in version 2.0.1"
         write(unit=iunit,fmt="(a)") "_refine_diff_density_max             ?"
         write(unit=iunit,fmt="(a)") "_refine_diff_density_min             ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_extinction_method         ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_extinction_coef           ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_details     ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_Flack       ?"
         write(unit=iunit,fmt="(a)") "_refine_ls_abs_structure_Rogers      ?"

         write(unit=iunit,fmt="(a)") " "
         write(unit=iunit,fmt="(a)") "# The following items are used to identify the programs used."
         write(unit=iunit,fmt="(a)") " "

         write(unit=iunit,fmt="(a)") "_computing_data_collection           ?"
         write(unit=iunit,fmt="(a)") "_computing_cell_refinement           ?"
         write(unit=iunit,fmt="(a)") "_computing_data_reduction            ?"
         write(unit=iunit,fmt="(a)") "_computing_structure_solution        ?"
         write(unit=iunit,fmt="(a)") "_computing_structure_refinement      ?"
         write(unit=iunit,fmt="(a)") "_computing_molecular_graphics        ?"
         write(unit=iunit,fmt="(a)") "_computing_publication_material      ?"

       End if  !(type_data < 2) then
       write(unit=iunit,fmt="(a)") " "
       write(unit=iunit,fmt="(a)") "#============================================================================="
       write(unit=iunit,fmt="(a)") " "
       !---- Atomic Coordinates and Displacement Parameters ----!
       write(unit=iunit,fmt="(a)") "# 9. ATOMIC COORDINATES AND DISPLACEMENT PARAMETERS"

       write(unit=iunit,fmt="(a)") " "

       write(unit=iunit,fmt="(a)") "loop_"
       write(unit=iunit,fmt='(a)') "    _atom_site_label"
       write(unit=iunit,fmt='(a)') "    _atom_site_fract_x"
       write(unit=iunit,fmt='(a)') "    _atom_site_fract_y"
       write(unit=iunit,fmt='(a)') "    _atom_site_fract_z"
       write(unit=iunit,fmt='(a)') "    _atom_site_U_iso_or_equiv"
       write(unit=iunit,fmt='(a)') "    _atom_site_occupancy"
       write(unit=iunit,fmt='(a)') "    _atom_site_adp_type              # Not in version 2.0.1"
       write(unit=iunit,fmt='(a)') "    _atom_site_type_symbol"

       !Calculation of the factor corresponding to the occupation factor provided in A
       do i=1,A%natoms
         occup(i)=A%Atom(i)%occ/(real(A%Atom(i)%mult)/real(SpG%multip))
         soccup(i)=A%Atom(i)%occ_std/(real(A%Atom(i)%mult)/real(SpG%multip))
       end do
       ocf=sum(abs(A%atom(1)%x-A%atom(2)%x))
       if( ocf < 0.001) then
         ocf=occup(1)+occup(2)
       else
         ocf=occup(1)
       end if
       occup=occup/ocf; soccup=soccup/ocf
       aniso=.false.
       do i=1,A%natoms
          line(1:132)=" "
          line(2:)= A%Atom(i)%Lab
           ! _atom_site_fract_x, _atom_site_fract_y, _atom_site_fract_z
          do j=1,3
            comm=" "
            call setnum_std(A%Atom(i)%x(j),A%Atom(i)%x_std(j),comm)
            line=trim(line)//" "//trim(comm)
          end do
          ! _atom_site_U_iso_or_equiv
            comm=" "
          if(A%Atom(i)%thtype == "isotr") then
            adptyp='Uiso'
            u=A%Atom(i)%Biso/(8.0*pi*pi)
            su=A%Atom(i)%Biso_std/(8.0*pi*pi)
            call setnum_std(u,su,comm)

          else if(A%Atom(i)%thtype == "aniso") then
            aniso=.true.
            adptyp='Uani'
            if(A%atom(i)%Utype == "beta") then
               aux=A%atom(i)%u
               ua=convert_betas_u(aux,cell)
               aux=A%atom(i)%u_std
               sua=convert_betas_u(aux,cell)
            else if(A%atom(i)%Utype == "u_ij") then
               ua=A%atom(i)%u
               sua=A%atom(i)%u_std
            end if
            u=(ua(1)+ua(2)+ua(3))/3.0
            su=(ua(1)+ua(2)+ua(3))/3.0
            call setnum_std(u,su,comm)
          else
            adptyp='.'
          end if
            line=trim(line)//" "//trim(comm)

           !_atom_site_occupancy
            comm=" "
            call setnum_std(occup(i),soccup(i),comm)
            line=trim(line)//" "//trim(comm)

          WRITE(iunit,"(a)") trim(line)//" "//trim(adptyp)//" "//A%atom(i)%SfacSymb
       end do
       if(aniso) then
          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_label "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_11  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_22  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_33  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_12  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_13  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_U_23  "
          write(unit=iunit,fmt="(a)") "    _atom_site_aniso_type_symbol"
          do i=1,A%natoms
             if(A%Atom(i)%thtype /= "aniso") cycle
             line(1:132)=" "
             line(2:)= A%Atom(i)%Lab
             if(A%atom(i)%Utype == "beta") then
                aux=A%atom(i)%u
                ua=convert_betas_u(aux,cell)
                aux=A%atom(i)%u_std
                sua=convert_betas_u(aux,cell)
             else if(A%atom(i)%Utype == "u_ij") then
                ua=A%atom(i)%u
                sua=A%atom(i)%u_std
             end if
             do j=1,6
               comm=" "
               call setnum_std(ua(j),sua(j),comm)
               line=trim(line)//" "//trim(comm)
             end do
              WRITE(iunit,"(a)") trim(line)//"  "//A%atom(i)%SfacSymb
          end do
       end if

       if(type_data < 2) then
          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "# Note: if the displacement parameters were refined anisotropically"
          write(unit=iunit,fmt="(a)") "# the U matrices should be given as for single-crystal studies."

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "#============================================================================="
          write(unit=iunit,fmt="(a)") " "

          !---- Molecular Geometry ----!
          write(unit=iunit,fmt="(a)") "# 10. MOLECULAR GEOMETRY"

          write(unit=iunit,fmt="(a)") " "


          write(unit=iunit,fmt="(a)") "_geom_special_details                ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "    _geom_bond_atom_site_label_1  "
          write(unit=iunit,fmt="(a)") "    _geom_bond_atom_site_label_2  "
          write(unit=iunit,fmt="(a)") "    _geom_bond_site_symmetry_1    "
          write(unit=iunit,fmt="(a)") "    _geom_bond_site_symmetry_2    "
          write(unit=iunit,fmt="(a)") "    _geom_bond_distance           "
          write(unit=iunit,fmt="(a)") "    _geom_bond_publ_flag          "
          write(unit=iunit,fmt="(a)") "    ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "    _geom_contact_atom_site_label_1 "
          write(unit=iunit,fmt="(a)") "    _geom_contact_atom_site_label_2 "
          write(unit=iunit,fmt="(a)") "    _geom_contact_distance          "
          write(unit=iunit,fmt="(a)") "    _geom_contact_site_symmetry_1   "
          write(unit=iunit,fmt="(a)") "    _geom_contact_site_symmetry_2   "
          write(unit=iunit,fmt="(a)") "    _geom_contact_publ_flag         "
          write(unit=iunit,fmt="(a)") "    ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_1 "
          write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_2 "
          write(unit=iunit,fmt="(a)") "_geom_angle_atom_site_label_3 "
          write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_1   "
          write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_2   "
          write(unit=iunit,fmt="(a)") "_geom_angle_site_symmetry_3   "
          write(unit=iunit,fmt="(a)") "_geom_angle                   "
          write(unit=iunit,fmt="(a)") "_geom_angle_publ_flag         "
          write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_1 "
          write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_2 "
          write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_3 "
          write(unit=iunit,fmt="(a)") "_geom_torsion_atom_site_label_4 "
          write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_1   "
          write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_2   "
          write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_3   "
          write(unit=iunit,fmt="(a)") "_geom_torsion_site_symmetry_4   "
          write(unit=iunit,fmt="(a)") "_geom_torsion                   "
          write(unit=iunit,fmt="(a)") "_geom_torsion_publ_flag         "
          write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "loop_"
          write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_D "
          write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_H "
          write(unit=iunit,fmt="(a)") "_geom_hbond_atom_site_label_A "
          write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_D   "
          write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_H   "
          write(unit=iunit,fmt="(a)") "_geom_hbond_site_symmetry_A   "
          write(unit=iunit,fmt="(a)") "_geom_hbond_distance_DH       "
          write(unit=iunit,fmt="(a)") "_geom_hbond_distance_HA       "
          write(unit=iunit,fmt="(a)") "_geom_hbond_distance_DA       "
          write(unit=iunit,fmt="(a)") "_geom_hbond_angle_DHA         "
          write(unit=iunit,fmt="(a)") "_geom_hbond_publ_flag         "
          write(unit=iunit,fmt="(a)") "?   ?   ?   ?   ?   ?   ?   ?   ?   ?   ?"

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "#============================================================================="
          write(unit=iunit,fmt="(a)") " "


          !---- Final Informations ----!
          write(unit=iunit,fmt="(a)") "#============================================================================="
          write(unit=iunit,fmt="(a)") "# Additional structures (last six sections and associated data_? identifiers) "
          write(unit=iunit,fmt="(a)") "# may be added at this point.                                                 "
          write(unit=iunit,fmt="(a)") "#============================================================================="

          write(unit=iunit,fmt="(a)") " "
          write(unit=iunit,fmt="(a)") "# The following lines are used to test the character set of files sent by     "
          write(unit=iunit,fmt="(a)") "# network email or other means. They are not part of the CIF data set.        "
          write(unit=iunit,fmt="(a)") "# abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789              "
          write(unit=iunit,fmt="(a)") "# !@#$%^&*()_+{}:"//""""//"~<>?|\-=[];'`,./ "
       end if

       close(unit=iunit)

       return
    End Subroutine Write_Cif_Template

    !!----
    !!---- Subroutine Write_Shx_Template(Filename,Code,Title,Lambda,Z,Celda,Space,Atomos)
    !!----    character(len=*),        intent(in) :: filename  !  In -> Filename
    !!----    integer,                 intent(in) :: code      !  In -> 0 Shelxs-Patterson
    !!----                                                              1 Shelxs-Direct Methods
    !!----                                                              2 Shelxl-Refinement
    !!----    character(len=*),        intent(in) :: title     !  In -> Title
    !!----    real(kind=cp),           intent(in) :: lambda    !  In -> Lambda
    !!----    integer,                 intent(in) :: z         !  In -> Z
    !!----    type(Crystal_cell_Type), intent(in) :: celda     !  In -> Cell variable
    !!----    type(Space_Group_Type),  intent(in) :: Space     !  In -> SpaceGroup variable
    !!----    type(atom_list_type),    intent(in) :: atomos    !  In -> Atom List
    !!----
    !!----    Write a Shelx File
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Shx_Template(filename,code,title,lambda,z,celda,space,atomos)
       !---- Arguments ----!
       character(len=*),        intent(in) :: filename
       integer,                 intent(in) :: code
       character(len=*),        intent(in) :: title
       real(kind=cp),           intent(in) :: lambda
       integer,                 intent(in) :: z
       type(Crystal_cell_Type), intent(in) :: celda
       type(Space_Group_Type),  intent(in) :: Space
       type(atom_list_type),    intent(in) :: atomos

       !---- Local Variables ----!
       logical                :: info

       integer                :: i,j,k,nc,iunit !,nlong
       integer                :: nlat
       integer, dimension(15) :: z_cont

       !---- Inicializacion de variables ----!
       info=.false.
       iunit=0
       z_cont=0
       nc=0  !this depends on scattering factor?

       !---- Esta abierto este Fichero? ----!
       inquire(file=filename,opened=info)
       if (info) then
          inquire(file=filename,number=iunit)
          close(unit=iunit)
       end if

       !---- Escritura ----!
       if (iunit == 0) iunit=61
       open(unit=iunit,file=filename,status="unknown",action="write")
       rewind(unit=iunit)

       !---- Title ----!
       write(unit=iunit,fmt="(a)") "TITL "//title(1:len_trim(title))

       !---- Lambda, Cell ----!
       write(unit=iunit,fmt="(a,f8.5,3f9.4,3f7.3)") "CELL ",lambda,celda%cell,celda%ang

       !---- Z, Std ----!
       write(unit=iunit,fmt="(a,i3,a,3f8.4,3f7.3)") "ZERR ",z,"     ",celda%cell_std,celda%ang_std

       !---- Latt ----!
       nlat=1
       select case (space%centred)
          case (0) ! Centric

          case (1) ! Acentric
             nlat=-1

          case (2) ! Not used in Shelx
             write(unit=iunit,fmt="(a)") " ERROR: Origin not at -1 "
             close(unit=iunit)
             return

       end select
       select case (space%spg_lat)
          case ("P")

          case ("I")
             nlat=2*nlat

          case ("R")
             nlat=3*nlat

          case ("F")
             nlat=4*nlat

          case ("A")
             nlat=5*nlat

          case ("B")
             nlat=6*nlat

          case ("C")
             nlat=7*nlat

       end select
       write(unit=iunit,fmt="(a,i2)") "LATT ",nlat

       !---- Symm ----!
       do i=2,space%numops
          write(unit=iunit,fmt="(a)") "SYMM "//u_case(space%symopsymb(i))
       end do

       !---- Sfac ----!
       j=0
       do i=1,atomos%natoms
          if (j == 0) then
             j=1
             z_cont(j)=atomos%atom(i)%z
          else
             do k=1,j
                if (z_cont(k) == atomos%atom(i)%z) exit
             end do
             if (z_cont(k) /= atomos%atom(i)%z) then
                j=j+1
                z_cont(j)=atomos%atom(i)%z
             end if
          end if
       end do


       write(unit=iunit,fmt="(a)") "SFAC "

       !---- Unit ----!
       write(unit=iunit,fmt="(a)") "UNIT "

       select case (code)
          case (0) ! Shelxs - Patterson
             write(unit=iunit,fmt="(a)") "PATT "

          case (1) ! Shelxs - Direct Methods
             write(unit=iunit,fmt="(a)") "TREF "

          case (2) ! Shelxl - Refinement
             !---- L.S. ----!
             write(unit=iunit,fmt="(a)") "L.S. 10"

             !---- Fvar ----!
             write(unit=iunit,fmt="(a)") "FVAR 1.0"

             !---- Weight ----!
             write(unit=iunit,fmt="(a)") "WGHT 0.2"

             !---- Fmap ----!
             write(unit=iunit,fmt="(a)") "FMAP 2"

             !---- Atoms ----!
             do i=1,atomos%natoms
                write(unit=iunit,fmt="(a4,i3,4f11.5)") &
                     atomos%atom(i)%lab, nc, atomos%atom(i)%x, atomos%atom(i)%occ+10.0
             end do
       end select

       !---- Format ----!
       write(unit=iunit,fmt="(a)") "HKLF 4"

       !---- End ----!
       write(unit=iunit,fmt="(a)") "END "

       return
    End Subroutine Write_Shx_Template

    !!----
    !!---- Subroutine Get_Phases_File(filecode, Nphas, PhasesName,ILines)
    !!----    character(len=*),                intent(in)   :: filecode
    !!----    Integer,                         intent(out)  :: Nphas
    !!----    Character(len=80), dimension(:), intent(out)  :: PhasesName
    !!----    Integer,dimension(2,:),          intent(out)  :: ILines
    !!----
    !!---- Determine how many phases there are in a CIF or PCR file and
    !!---- give the lines to locate
    !!----
    !!---- Update: 01/05/2013
    !!
    Subroutine Get_Phases_File(filecode, NPhas, PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),             intent(in)   :: filecode
       integer,                      intent(out)  :: Nphas
       character(len=*),dimension(:),intent(out)  :: PhasesName
       integer,dimension(:,:),       intent(out)  :: ILines

       !---- Local Variables ----!
       character(len=3) :: ext
       integer          :: npos

       !> Error
       call init_err_form()

       !> Init
       Nphas=0
       PhasesName=' '
       Ilines=0

       !> PCR or CIF file
       npos=index(filecode,'.',back=.true.)
       if (npos <=0) then
          err_form=.true.
          err_form_mess='No extension was found in the name of the file!'
          return
       end if

       ext=filecode(npos+1:)
       ext=u_case(ext)
       select case (ext)
          case ('CIF')
             call get_nphases_ciffile(filecode, NPhas, PhasesName,ILines)
          case ('PCR')
             call get_nphases_pcrfile(filecode, NPhas, PhasesName,ILines)
          case default
             err_form=.true.
             err_form_mess='Extension for this file not valid!'
       end select

       return
    End Subroutine Get_Phases_File

    !!--++
    !!--++ Subroutine Get_NPhases_CIFFile(Filecode,NPhas,PhasesName,ILines)
    !!--++    character(len=*),                 intent(in)  :: Filecode    ! Filename
    !!--++    integer,                          intent(out) :: NPhas       ! Number of Phases in the file
    !!--++    character(len=*), dimension(:),   intent(out) :: PhasesName     ! Name of Phases in the file
    !!--++    integer,          dimension(:,:), intent(out) :: ILines        ! Index for lines for each Phase
    !!--++
    !!--++ Determine the number of phases are included into the file
    !!--++
    !!--++ Date: 01/05/2013
    !!
    Subroutine Get_NPhases_CIFFile(Filecode,NPhas,PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),                 intent(in)  :: Filecode    ! Filename
       integer,                          intent(out) :: NPhas       ! Number of Phases in the file
       character(len=*), dimension(:),   intent(out) :: PhasesName     ! Name of Phases in the file
       integer,          dimension(:,:), intent(out) :: ILines        ! Index for lines for each Phase

       !---- Local Variables ----!
       character(len=150), dimension(:), allocatable :: filen
       character(len=150)                            :: line
       integer                                       :: i,j,nl

       !> Error
       call init_err_form()

       !> Initialize
       NPhas=0
       PhasesName=' '
       ILines=0

       !> Reading file
       nl=0
       call number_lines(trim(filecode),nl)
       if (nl <=0) then
          err_form=.true.
          err_form_mess='No lines were readen for '//trim(filecode)//' !!'
          return
       end if
       allocate(filen(nl))
       call reading_lines(trim(filecode),nl,filen)

       !> Number of Phases
       do i=1,nl
          line=adjustl(filen(i))

          !> empty line
          if (len_trim(line) <= 0) cycle

          !> comment line
          if (line(1:1) =='#') cycle

          !> No data_global
          j=index(line,'data_global')
          if (j > 0) cycle

          !> Just only lines beginning with data...
          j=index(line,'data_')
          if (j /= 1) cycle

          nphas=nphas+1
          ILines(1,Nphas)=i
          PhasesName(nphas)=trim(line(j+5:))
          if (nphas > 1) ILines(2,nphas-1)=i-1
       end do
       if (nphas > 0 .and. ILines(2,nphas)==0) ILines(2,nphas)=nl

       if (allocated(filen)) deallocate(filen)

       return
    End Subroutine  Get_NPhases_CIFFile

    !!--++
    !!--++ Subroutine Get_NPhases_PCRFile(filecode, Nphas,PhasesName,ILines)
    !!--++    character(len=*),                intent(in)   :: filecode
    !!--++    Integer,                         intent(out)  :: Nphas
    !!--++    Character(len=80), dimension(:), intent(out)  :: PhasesName
    !!--++    Integer,dimension(2,:),          intent(out)  :: ILines
    !!--++
    !!--++ Determine how many phases and where there in a PCR file
    !!--++
    !!--++ Update: 01/05/2013
    !!
    Subroutine Get_NPhases_PCRFile(filecode, NPhas, PhasesName,ILines)
       !---- Arguments ----!
       character(len=*),             intent(in)   :: filecode
       integer,                      intent(out)  :: Nphas
       character(len=*),dimension(:),intent(out)  :: PhasesName
       integer,dimension(:,:),       intent(out)  :: ILines

       !---- Local Variables ----!
       logical                                      :: multi, ask_phase
       character(len=80), dimension(:), allocatable :: file_dat
       character(len=80)                            :: line
       integer                                      :: i,k,iv,nlines
       integer, dimension(30)                       :: ivet
       real(kind=cp), dimension(30)                 :: vet

       !> Err
       call init_err_form()

       !> Init
       NPhas=0
       PhasesName=' '
       ILines=0

       !> Reading file
       nlines=0
       call number_lines(trim(filecode),nlines)
       if (nlines <=0) then
          err_form=.true.
          err_form_mess='No lines were readen for '//trim(filecode)//' !!'
          return
       end if
       allocate(file_dat(nlines))
       call reading_lines(trim(filecode),nlines,file_dat)

       ILines(1,:)=1
       ILines(2,:)=nlines

       !> Simple / Multi format
       multi=.false.
       do i=1,nlines
          line=adjustl(file_dat(i))
          if (line(1:1) =='!' .or. line(1:1)==' ') cycle
          if (index(line,'NPATT ') <=0) cycle
          multi=.true.
       end do

       !> Number of Phases
       if (.not. multi) then
          do i=2,nlines
             line=adjustl(file_dat(i))
             if (line(1:1) =='!' .or. line(1:1)==' ') cycle
             call getnum(line,vet,ivet,iv)
             if (iv > 3) then
                NPhas=ivet(3)
                exit
             end if
          end do

       else
          do i=1,nlines
             line=adjustl(file_dat(i))
             if (line(1:4) /='!Nph') cycle

             line=adjustl(file_dat(i+1))
             call getnum(line,vet,ivet,iv)
             if (iv > 1) then
                NPhas=ivet(1)
                exit
             end if
          end do
       end if

       if (NPhas == 0) then
          err_form=.true.
          err_form_mess=" No Phase information was found in this PCR file. Please, check it! "
          return
       end if

       !> Locate where begin each Phase
       k=0
       ask_phase=.true.

       do i=1,nlines
          line=adjustl(file_dat(i))
          if (ask_phase) then
             if (index(line,'Data for PHASE') <= 0) cycle
          else
             if (line(1:1) /='!') then
                k=k+1
                ILines(1,k)=i
                PhasesName(k)=trim(adjustl(line))
                if (k == NPhas) exit

                ask_phase=.true.
             end if
             cycle
          end if
          ask_phase=.false.
       end do

       if (NPhas /= k) then
          err_form=.true.
          err_form_mess=" Locating Phases failed in this PCR. Please, check it!"
          return
       end if

       do i=1,Nphas
          if (nphas > 1) then
             ilines(2,i)=ilines(1,i+1)-1
          end if
       end do

       return
    End Subroutine Get_NPhases_PCRFile
    !!----
    !!---- Subroutine Write_CFL(lun,Cel,SpG,Atm,comment)
    !!----    integer,                  intent(in)    :: lun
    !!----    type (Space_Group_Type),  intent(in)    :: SpG
    !!----    type (Crystal_Cell_Type), intent(in)    :: Cel
    !!----    type (atom_list_type),    intent(in)    :: Atm
    !!----    character(len=*),optional,intent(in)    :: comment
    !!----
    !!----    (OVERLOADED)
    !!----
    !!----    Write a CFL-file with atom_list_type
    !!----
    !!---- Update: July - 2014
    !!
    Subroutine Write_CFL_Atom_List_Type(lun,Cel,SpG,Atm,comment)
       !---- Arguments ----!
       integer,                  intent(in)    :: lun
       type (Space_Group_Type),  intent(in)    :: SpG
       type (Crystal_Cell_Type), intent(in)    :: Cel
       type (atom_list_type),    intent(in)    :: Atm
       character(len=*),optional,intent(in)    :: comment

       !----- Local variables -----!
       integer                         :: j !,loc
       real(kind=cp), dimension(6)     :: a,sa
       character(len=30), dimension(6) :: text

       if(present(comment)) write(unit=lun,fmt="(a)") "TITLE "//trim(comment)
       write(unit=lun,fmt="(a)") "!  Automatically generated CFL file (Write_CFL)"

       a(1:3)=Cel%Cell
       a(4:6)=Cel%ang
       sa(1:3)=Cel%Cell_std
       sa(4:6)=Cel%ang_std
       do j=1,6
          call SetNum_Std(a(j), sa(j), text(j))
       end do
       write(unit=lun,fmt="(a)") "!         a               b               c            alpha           beta            gamma"
       write(unit=lun,fmt="(a,6a16)") "Cell ",text
       write(unit=lun,fmt="(a,i3)")"!     Space Group # ",SpG%NumSpg
       write(unit=lun,fmt="(a,a)") "Spgr  ",SpG%SPG_Symb
       call Write_Atoms_CFL(Atm,Lun,cel)

       return
    End Subroutine Write_CFL_Atom_List_Type
    !!----
    !!---- Subroutine Write_CFL(lun,Molx,comment)
    !!----    integer,                       intent(in) :: lun
    !!----    type (Molecular_Crystal_Type), intent(in) :: Molx
    !!----    character(len=*),optional,     intent(in) :: comment
    !!----
    !!----    (OVERLOADED)
    !!----
    !!----    Write a CFL-file with molecular_crystal_type
    !!----
    !!---- Update: July - 2014
    !!
    Subroutine Write_CFL_Molcrys(lun,Molx,comment)
       !---- Arguments ----!
       integer,                       intent(in) :: lun
       type (Molecular_Crystal_Type), intent(in) :: Molx
       character(len=*),optional,     intent(in) :: comment

       !----- Local variables -----!
       integer                         :: j !,loc
       real(kind=cp), dimension(6)     :: a,sa
       character(len=30), dimension(6) :: text

       if(present(comment)) write(unit=lun,fmt="(a)") "TITLE "//trim(comment)
       write(unit=lun,fmt="(a)") "!  Automatically generated CFL file (Write_CFL)"

       a(1:3)=molx%cell%Cell
       a(4:6)=molx%cell%ang
       sa(1:3)=molx%cell%Cell_std
       sa(4:6)=molx%cell%ang_std
       do j=1,6
          call SetNum_Std(a(j), sa(j), text(j))
       end do
       write(unit=lun,fmt="(a)") "!         a               b               c            alpha           beta            gamma"
       write(unit=lun,fmt="(a,6a16)") "Cell ",text
       write(unit=lun,fmt="(a,i3)")"!     Space Group # ",molx%spg%NumSpg
       write(unit=lun,fmt="(a,a)") "Spgr  ",molx%spg%SPG_Symb
       call Write_Atoms_CFL(Molx,Lun)

       return
    End Subroutine Write_CFL_Molcrys
    !!----
    !!---- Subroutine Write_Atoms_CFL(Ats,Lun,Cell)
    !!----    Type (atom_list_type),dimension(:),  intent(in) :: Ats     !  In -> Atom List
    !!----    integer, optional,                   intent(in) :: lun     !  In -> Unit to write
    !!----    Type(Crystal_Cell_Type), optional,   intent(in) :: Cell    !  In -> Transform to thermal parameters
    !!----
    !!----    Write the atoms in the asymmetric unit for a CFL file
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Write_Atoms_CFL_ATM(Ats,Lun,cell)
       !---- Arguments ----!
       type (atom_list_type),            intent(in) :: Ats
       integer, optional,                intent(in) :: Lun
       Type(Crystal_Cell_Type), optional,intent(in) :: Cell

       !---- Local Variables ----!
       character(len=30),dimension(6) :: text
       character(len=36)              :: forma,fom
       integer                        :: i, j, iunit, leng, maxl,ish
       real(kind=cp), dimension(6)    :: u,bet,sb

       iunit=6
       if (present(lun)) iunit=lun

       if(ats%natoms == 0) then
         write (unit=iunit,fmt="(a)") "!  No atoms ..."
         return
       end if
       !Determine the maximum length of the atom labels
       maxl=0
       do i=1,ats%natoms
         leng=len_trim(ats%atom(i)%lab)
         if(leng > maxl) maxl=leng
       end do
       maxl=max(maxl,4)+1
       ish=maxl-4
       fom   ="(a,tr  ,a)"
       Select Case(ish)
          Case(:9)
            write(unit=fom(6:6),fmt="(i1)") ish
          Case(10:)
            write(unit=fom(6:7),fmt="(i2)") ish
       End Select
       forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"
       Select Case(maxl)
         Case(:9)
             write(unit=forma(5:5),fmt="(i1)") maxl
         Case(10:)
             write(unit=forma(5:6),fmt="(i2)") maxl
       End Select
       write (unit=iunit,fmt=fom) "!     ", &
             "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"
       do i=1,ats%natoms

          do j=1,3
             call SetNum_Std(ats%atom(i)%x(j), ats%atom(i)%x_std(j), text(j))
          end do
          call SetNum_Std(ats%atom(i)%Biso, ats%atom(i)%Biso_std, text(4))
          call SetNum_Std(ats%atom(i)%Occ, ats%atom(i)%Occ_std, text(5))

          write (unit=iunit,fmt=forma) &
                "Atom   ",trim(ats%atom(i)%lab),ats%atom(i)%chemsymb, (text(j),j=1,5), &
                 ats%atom(i)%moment,ats%atom(i)%charge,"# "//ats%atom(i)%AtmInfo

          if (ats%atom(i)%thtype == "aniso") then

             if (ats%atom(i)%utype == "beta") then
                bet=ats%atom(i)%u(1:6)
                sb=ats%atom(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(bet(j), sb(j), text(j))
                end do
                write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                if (present(Cell)) then
                   u=convert_betas_u(bet,cell)
                   sb=convert_betas_u(ats%atom(i)%u_std,cell)
                   do j=1,6
                      call SetNum_Std(u(j), sb(j), text(j))
                   end do
                   write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
                end if

             else if(ats%atom(i)%thtype == "u_ij") then
                u=ats%atom(i)%u(1:6)
                sb=ats%atom(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(u(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "U_ij  ", text
                if (present(Cell)) then
                   bet=convert_u_betas(u,cell)
                   sb=convert_u_betas(ats%atom(i)%u_std,cell)
                   do j=1,6
                      call SetNum_Std(bet(j), sb(j), text(j))
                   end do
                   write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
                end if
             end if

          end if
       end do

       return
    End Subroutine Write_Atoms_CFL_ATM

    !!----
    !!---- Subroutine Write_Atoms_CFL(Ats,Lun,Cell)
    !!----    Type (atom_list_type),dimension(:),  intent(in) :: Ats     !  In -> Atom List
    !!----    integer, optional,                   intent(in) :: lun     !  In -> Unit to write
    !!----    Type(Crystal_Cell_Type), optional,   intent(in) :: Cell    !  In -> Transform to thermal parameters
    !!----
    !!----    Write the atoms in the asymmetric unit for a CFL file
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Write_Atoms_CFL_MOLX(Molx,Lun)
        !---- Arguments ----!
        type (Molecular_Crystal_Type), intent(in) :: Molx
        integer, optional,             intent(in) :: Lun

        !---- Local Variables ----!
        character(len=30),dimension(6) :: text
        character(len=36)              :: forma,fom
        integer                        :: i, j, iunit, leng, maxl,ish
        real(kind=cp), dimension(6)    :: u,bet,sb

        iunit=6
        if (present(lun)) iunit=lun

        if(molx%n_free > 0) then
            !Determine the maximum length of the atom labels
            maxl=0
            do i=1,molx%n_free
                leng=len_trim(molx%atm(i)%lab)
                if(leng > maxl) maxl=leng
            end do
            maxl=max(maxl,4)+1
            ish=maxl-4
            fom   ="(a,tr  ,a)"
            Select Case(ish)
                Case(:9)
                    write(unit=fom(6:6),fmt="(i1)") ish
                Case(10:)
                    write(unit=fom(6:7),fmt="(i2)") ish
            End Select
            forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"
            Select Case(maxl)
                Case(:9)
                    write(unit=forma(5:5),fmt="(i1)") maxl
                Case(10:)
                    write(unit=forma(5:6),fmt="(i2)") maxl
            End Select
            write (unit=iunit,fmt=fom) "!     ", &
                  "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"
            do i=1,molx%n_free

                do j=1,3
                   call SetNum_Std(molx%atm(i)%x(j), molx%atm(i)%x_std(j), text(j))
                end do
                call SetNum_Std(molx%atm(i)%Biso, molx%atm(i)%Biso_std, text(4))
                call SetNum_Std(molx%atm(i)%Occ, molx%atm(i)%Occ_std, text(5))

                write (unit=iunit,fmt=forma) &
                      "Atom   ",trim(molx%atm(i)%lab),molx%atm(i)%chemsymb, (text(j),j=1,5), &
                       molx%atm(i)%moment,molx%atm(i)%charge,"# "//molx%atm(i)%AtmInfo

                if (molx%atm(i)%thtype == "aniso") then

                    if (molx%atm(i)%utype == "beta") then
                        bet=molx%atm(i)%u(1:6)
                        sb=molx%atm(i)%u_std(1:6)
                        do j=1,6
                            call SetNum_Std(bet(j), sb(j), text(j))
                        end do
                        write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                        u=convert_betas_u(bet,molx%cell)
                        sb=convert_betas_u(molx%atm(i)%u_std,molx%cell)
                        do j=1,6
                            call SetNum_Std(u(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
                    else if(molx%atm(i)%thtype == "u_ij") then
                        u=molx%atm(i)%u(1:6)
                        sb=molx%atm(i)%u_std(1:6)
                        do j=1,6
                            call SetNum_Std(u(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "U_ij  ", text
                        bet=convert_u_betas(u,molx%cell)
                        sb=convert_u_betas(molx%atm(i)%u_std,molx%cell)
                        do j=1,6
                            call SetNum_Std(bet(j), sb(j), text(j))
                        end do
                        write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
                    end if
                end if
            end do ! i=1,molx%n_free
        end if ! molx%n_free > 0

        if (molx%n_mol > 0) then
            do i=1,molx%n_mol
                write(unit=iunit,fmt="(/,a,tr2,i3,tr2,a,tr2,a)") &
                     "MOLEX",molx%mol(i)%natoms,trim(molx%mol(i)%Name_mol),molx%mol(i)%coor_type
                write(unit=iunit,fmt="(a)") &
                     "!    Xc         Yc          Zc        Phi        Theta      Chi     TypeAngles TypeThermal"
                write(unit=iunit,fmt="(6f11.5,tr6,a,tr10,a)") &
                     molx%mol(i)%xcentre,molx%mol(i)%orient,molx%mol(i)%rot_type,molx%mol(i)%therm_type
                write(unit=iunit,fmt="(t1,6i10,tr2,a)") &
                     molx%mol(i)%lxcentre,molx%mol(i)%lorient," ! Refinemencodes"

                select case (molx%mol(i)%coor_type)
                    case ("C","c")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type        XC          YC          ZC    N1  N2  N3      Biso        Occ "
                    case ("F","f")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type        X           Y           Z     N1  N2  N3      Biso        Occ "
                    case ("S","s")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type    distance      Theta       Phi     N1  N2  N3      Biso        Occ "
                    case ("Z","z")
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type    distance  Bond-Angle Torsion-Ang  N1  N2  N3      Biso        Occ "
                    case default
                        write(unit=iunit,fmt="(a)") &
                        "!Atom   Type      Coor1       Coor2       Coor3   N1  N2  N3      Biso        Occ "
                end select ! molx%mol(i)%coor_type

                do j=1,molx%mol(i)%natoms
                    write(unit=iunit,fmt="(a,tr2,a,3f12.5,3i4,2f12.5)")  &
                          molx%mol(i)%AtName(j), molx%mol(i)%AtSymb(j),molx%mol(i)%I_Coor(:,j),  &
                          molx%mol(i)%Conn(:,j), molx%mol(i)%Biso(j),  molx%mol(i)%Occ(j)
                end do ! j = molx%mol(i)%natoms
            end do ! i = 1,molx%n_mol
        end if ! molx%n_mol > 0
        return
    End Subroutine Write_Atoms_CFL_MOLX
    !!----
    !!---- Subroutine Write_Atoms_CFL(Ats,Lun,Cell)
    !!----    Type (atom_list_type),dimension(:),  intent(in) :: Ats     !  In -> Atom List
    !!----    integer, optional,                   intent(in) :: lun     !  In -> Unit to write
    !!----    Type(Crystal_Cell_Type), optional,   intent(in) :: Cell    !  In -> Transform to thermal parameters
    !!----
    !!----    Write the atoms in the asymmetric unit for a CFL file
    !!----
    !!---- Update: February - 2003
    !!
    Subroutine Write_Atoms_CFL_MOLX_orig(Molx,Lun)
       !---- Arguments ----!
       type (Molecular_Crystal_Type), intent(in) :: Molx
       integer, optional,             intent(in) :: Lun

       !---- Local Variables ----!
       character(len=30),dimension(6) :: text
       character(len=36)              :: forma,fom
       integer                        :: i, j, iunit, leng, maxl,ish
       real(kind=cp), dimension(6)    :: u,bet,sb

       iunit=6
       if (present(lun)) iunit=lun

       if(molx%n_free == 0) then
         write (unit=iunit,fmt="(a)") "!  No atoms ..."
         return
       end if
       !Determine the maximum length of the atom labels
       maxl=0
       do i=1,molx%n_free
         leng=len_trim(molx%atm(i)%lab)
         if(leng > maxl) maxl=leng
       end do
       maxl=max(maxl,4)+1
       ish=maxl-4
       fom   ="(a,tr  ,a)"
       Select Case(ish)
          Case(:9)
            write(unit=fom(6:6),fmt="(i1)") ish
          Case(10:)
            write(unit=fom(6:7),fmt="(i2)") ish
       End Select
       forma="(a,a  ,tr2,a,tr3,5a14,2f8.2,tr3,a)"
       Select Case(maxl)
         Case(:9)
             write(unit=forma(5:5),fmt="(i1)") maxl
         Case(10:)
             write(unit=forma(5:6),fmt="(i2)") maxl
       End Select
       write (unit=iunit,fmt=fom) "!     ", &
             "Atom  Type     x/a           y/b           z/c           Biso          Occ           Spin    Charge    Info"
       do i=1,molx%n_free

          do j=1,3
             call SetNum_Std(molx%atm(i)%x(j), molx%atm(i)%x_std(j), text(j))
          end do
          call SetNum_Std(molx%atm(i)%Biso, molx%atm(i)%Biso_std, text(4))
          call SetNum_Std(molx%atm(i)%Occ, molx%atm(i)%Occ_std, text(5))

          write (unit=iunit,fmt=forma) &
                "Atom   ",trim(molx%atm(i)%lab),molx%atm(i)%chemsymb, (text(j),j=1,5), &
                 molx%atm(i)%moment,molx%atm(i)%charge,"# "//molx%atm(i)%AtmInfo

          if (molx%atm(i)%thtype == "aniso") then

             if (molx%atm(i)%utype == "beta") then
                bet=molx%atm(i)%u(1:6)
                sb=molx%atm(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(bet(j), sb(j), text(j))
                end do
                write (unit=iunit,fmt="(a,tr1,6a14)") "Beta  ", text
                u=convert_betas_u(bet,molx%cell)
                sb=convert_betas_u(molx%atm(i)%u_std,molx%cell)
                do j=1,6
                    call SetNum_Std(u(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "!U_ij  ", text
             else if(molx%atm(i)%thtype == "u_ij") then
                u=molx%atm(i)%u(1:6)
                sb=molx%atm(i)%u_std(1:6)
                do j=1,6
                   call SetNum_Std(u(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "U_ij  ", text
                bet=convert_u_betas(u,molx%cell)
                sb=convert_u_betas(molx%atm(i)%u_std,molx%cell)
                do j=1,6
                    call SetNum_Std(bet(j), sb(j), text(j))
                end do
                write(unit=iunit,fmt="(a,6a14)") "!Beta  ", text
             end if
          end if
       end do

       return
    End Subroutine Write_Atoms_CFL_MOLX_orig

 End Module CFML_IO_Formats

