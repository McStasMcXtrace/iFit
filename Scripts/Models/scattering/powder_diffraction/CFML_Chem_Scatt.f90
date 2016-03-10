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
!!---- MODULE: CFML_Scattering_Chemical_Tables
!!----
!!----   INFO: Tabulated information about atomic chemical and scattering data. A set of fortran
!!----         TYPEs and variables are defined. Tables are declared as allocatable arrays of
!!----         types and they are charged only if the setting (initialising) procedures are called.
!!----         It is convenient in a particular program using this moduled to call the "removing"
!!----         procedures (making a deallocation) to liberate memory after the required information
!!----         is found and stored in user-defined variables.
!!----
!!---- HISTORY
!!----    Updated: 04/03/2011
!!----
!!---- DEPENDENCIES
!!--++    Use CFML_GlobalDeps,       only: Cp
!!--++    Use CFML_String_Utilities, only: L_Case, U_Case
!!----
!!---- VARIABLES
!!----    ANOMALOUS_SC_TYPE
!!----    ANOMALOUS_SCFAC
!!----    CHEM_INFO_TYPE
!!----    CHEM_INFO
!!----    MAGNETIC_FORM_TYPE
!!----    MAGNETIC_FORM
!!----    MAGNETIC_J2
!!----    MAGNETIC_J4
!!----    MAGNETIC_J6
!!----    NUM_CHEM_INFO
!!----    NUM_DELTA_FP
!!----    NUM_MAG_FORM
!!----    NUM_MAG_J2
!!----    NUM_MAG_J4
!!----    NUM_MAG_J6
!!----    NUM_XRAY_FORM
!!----    XRAY_FORM_TYPE
!!----    XRAY_FORM
!!----    XRAY_WAVELENGTH_TYPE
!!----    XRAY_WAVELENGTHS
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       GET_ATOMIC_MASS
!!----       GET_ATOMIC_VOL
!!----       GET_CHEMSYMB
!!----       GET_COVALENT_RADIUS
!!----       GET_FERMI_LENGTH
!!----       GET_ABS_XS
!!----       GET_INC_XS
!!----       GET_IONIC_RADIUS
!!----       REMOVE_CHEM_INFO
!!----       REMOVE_DELTA_FP_FPP
!!----       REMOVE_MAGNETIC_FORM
!!----       REMOVE_XRAY_FORM
!!----       SET_CHEM_INFO
!!----       SET_DELTA_FP_FPP
!!----       SET_MAGNETIC_FORM
!!----       SET_XRAY_FORM
!!----
!!
 Module CFML_Scattering_Chemical_Tables
    !---- Use Modules ----!
    Use CFML_GlobalDeps,       only: Cp
    Use CFML_String_Utilities, only: U_Case, L_Case

    implicit none

    private

    !---- List of public subroutines ----!
    public :: Get_Atomic_Mass, Get_Atomic_Vol, Get_ChemSymb, Get_Covalent_radius, Get_Ionic_radius
    public :: Get_Fermi_Length, Get_Abs_Xs, Get_Inc_Xs
    public :: Remove_Chem_Info, Remove_Delta_Fp_Fpp, Remove_Magnetic_Form, Remove_Xray_Form
    public :: Set_Chem_Info, Set_Delta_Fp_Fpp, Set_Magnetic_Form, Set_Xray_Form

    !---- Definitions ----!

    !!----
    !!---- TYPE, PUBLIC :: ANOMALOUS_SC_TYPE
    !!--..
    !!---- Type, public :: Anomalous_Sc_Type
    !!----    character (len= 2)           :: Symb  ! Symbol of the Chemical species
    !!----    real(kind=cp), dimension(5)  :: Fp    ! Delta Fp
    !!----    real(kind=cp), dimension(5)  :: Fpp   ! Delta Fpp
    !!---- End Type Anomalous_Sc_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Anomalous_Sc_Type
       character(len= 2)            :: Symb
       real(kind=cp), dimension(5)  :: Fp
       real(kind=cp), dimension(5)  :: Fpp
    End Type Anomalous_Sc_Type

    !!----
    !!---- ANOMALOUS_SCFAC
    !!----    Type(Anomalous_Sc_Type), allocatable, dimension(:), public :: Anomalous_ScFac
    !!----
    !!----    Table of Delta-Fp and Delta-Fpp for 5 common radiations.
    !!----    The order is the following:
    !!--<<
    !!----                          1         2         3          4          5
    !!----        Wavelenghts:     Cr        Fe        Cu         Mo         Ag
    !!----             Lambda   2.28962   1.93597   1.54051    0.70926    0.556363
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Type(Anomalous_Sc_Type), allocatable, dimension(:), public :: Anomalous_ScFac

    !!----
    !!---- TYPE, PUBLIC :: CHEM_INFO_TYPE
    !!--..
    !!---- Type, public :: Chem_Info_Type
    !!----    character (len= 2)         :: Symb     ! Symbol of the Element
    !!----    character (len=12)         :: Name     ! Name of the Element
    !!----    integer                    :: Z        ! Atomic Number
    !!----    real(kind=cp)              :: AtWe     ! Atomic weight
    !!----    real(kind=cp)              :: RCov     ! Covalent Radio
    !!----    real(kind=cp)              :: RWaals   ! van der Waals Radio
    !!----    real(kind=cp)              :: VAtm     ! Atomic volumen
    !!----    integer, dimension(5)      :: Oxid     ! Oxidation State
    !!----    real(kind=cp), dimension(5):: Rion     ! Ionic Radio (depending of the oxidation)
    !!----    real(kind=cp)              :: SctF     ! Scattering length Fermi
    !!----    real(kind=cp)              :: SedInc   ! Incoherent Scattering Neutron cross-section (barns -> [10**(-24) cm**2] )
    !!----    real(kind=cp)              :: Sea      ! Neutron Absorption cross-section ( barns, for v= 2200m/s, l(A)=3.95/v (km/s) )
    !!---- End Type Chem_Info_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Chem_Info_Type
       character (len= 2)         :: Symb          ! Symbol of the Element
       character (len=12)         :: Name          ! Name of the Element
       integer                    :: Z             ! Atomic Number
       real(kind=cp)              :: AtWe          ! Atomic weight
       real(kind=cp)              :: RCov          ! Covalent Radius
       real(kind=cp)              :: RWaals        ! van der Waals Radius
       real(kind=cp)              :: VAtm          ! Atomic volumen
       integer, dimension(5)      :: Oxid          ! Oxidation State
       real(kind=cp), dimension(5):: Rion          ! Ionic Radius (depending of the oxidation)
       real(kind=cp)              :: SctF          ! Fermi length [10**(-12) cm]
       real(kind=cp)              :: SedInc        ! Incoherent Scattering Neutron cross-section (barns -> [10**(-24) cm**2] )
       real(kind=cp)              :: Sea           ! Neutron Absorption cross-section ( barns, for v= 2200m/s, l(A)=3.95/v (km/s) )
    End Type Chem_Info_Type

    !!----
    !!---- CHEM_INFO
    !!----    Type (Chem_Info_Type), allocatable, dimension(:), public :: Chem_Info
    !!----
    !!----    Tabulated chemical data according to the items specified in the definition of Chem_Info_Type.
    !!----
    !!---- Update: February - 2005
    !!
    Type(Chem_Info_Type), allocatable, dimension(:), public :: Chem_Info

    !!----
    !!---- TYPE :: MAGNETIC_FORM_TYPE
    !!--..
    !!---- Type, public :: Magnetic_Form_Type
    !!----    character (len= 4)          :: Symb   ! Symbol of the Chemical species
    !!----    real(kind=cp), dimension(7) :: SctM   ! Scattering Factors
    !!---- End Type Magnetic_Form_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Magnetic_Form_Type
       character (len= 4)         :: Symb         ! Symbol of the Chemical species
       real(kind=cp), dimension(7):: SctM
    End Type Magnetic_Form_Type

    !!----
    !!---- MAGNETIC_FORM
    !!----    Type (Magnetic_Form_Type), allocatable, dimension(:), public :: Magnetic_Form
    !!----
    !!----    Tabulated magnetic form factor data
    !!----
    !!---- Update: February - 2005
    !!
    Type(Magnetic_Form_Type), allocatable, dimension(:), public :: Magnetic_Form

    !!----
    !!---- MAGNETIC_J2
    !!----    Type (Magnetic_Form_Type), allocatable, dimension(:), public :: Magnetic_j2
    !!----
    !!----    Tabulated magnetic form factor J2
    !!----
    !!---- Update: February - 2005
    !!
    Type(Magnetic_Form_Type), allocatable, dimension(:), public :: Magnetic_j2

    !!----
    !!---- MAGNETIC_J4
    !!----    Type (Magnetic_Form_Type), allocatable, dimension(:), public :: Magnetic_J4
    !!----
    !!----    Tabulated magnetic form factor J4
    !!----
    !!---- Update: February - 2005
    !!
    Type(Magnetic_Form_Type), allocatable, dimension(:), public :: Magnetic_j4

    !!----
    !!---- MAGNETIC_J6
    !!----    Type (Magnetic_Form_Type), allocatable, dimension(:), public :: Magnetic_J6
    !!----
    !!----    Tabulated magnetic form factor J6
    !!----
    !!---- Update: February - 2005
    !!
    Type(Magnetic_Form_Type), allocatable, dimension(:), public :: Magnetic_j6

    !!----
    !!---- NUM_CHEM_INFO
    !!----    integer, parameter, public :: Num_Chem_Info = 108
    !!----
    !!----    Number of total Chem_info Data
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Num_Chem_Info = 108

    !!----
    !!---- NUM_DELTA_FP
    !!----    integer, parameter, public :: Num_Delta_Fp  = 98
    !!----
    !!----    Number of total Delta (Fp,Fpp) Data
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Num_Delta_Fp  = 98

    !!----
    !!---- NUM_MAG_FORM
    !!----    integer, parameter, public :: Num_Mag_Form  = 117
    !!----
    !!----    Number of total Magnetic_Form Data
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Num_Mag_Form  = 119

    !!----
    !!---- NUM_MAG_J2
    !!----    integer, parameter, public :: Num_Mag_J2 = 97
    !!----
    !!----    Number of <j2> Magnetic_Form Data
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Num_Mag_j2  = 97

    !!----
    !!---- NUM_MAG_J4
    !!----    integer, parameter, public :: Num_Mag_J4 = 97
    !!----
    !!----    Number of <j4> Magnetic_Form Data
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Num_Mag_j4  = 97

    !!----
    !!---- NUM_MAG_J6
    !!----    integer, parameter, public :: Num_Mag_J6 = 39
    !!----
    !!----    Number of <j5> Magnetic_Form Data
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Num_Mag_j6  = 39

    !!----
    !!---- NUM_XRAY_FORM
    !!----    integer, parameter, public :: Num_Xray_Form = 214
    !!----
    !!----    Number of total Xray_Form Data
    !!----
    !!---- Update: February - 2005
    !!
    integer, parameter, public :: Num_Xray_Form = 214

    !!----
    !!---- TYPE :: XRAY_FORM_TYPE
    !!--..
    !!---- Type, public :: Xray_Form_Type
    !!----    character (len= 4)         :: Symb  ! Symbol of the Chemical species
    !!----    integer                    :: Z     ! Atomic Number
    !!----    real(kind=cp), dimension(4):: a     ! Coefficients for calculating the X-ray scattering factors
    !!----    real(kind=cp), dimension(4):: b     ! f(s) = Sum_{i=1,4} { a(i) exp(-b(i)*s^2) } + c
    !!----    real(kind=cp)              :: c     ! s=sinTheta/Lambda
    !!---- End Type Xray_Form_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Xray_Form_Type
       character (len= 4)         :: Symb
       integer                    :: Z
       real(kind=cp), dimension(4):: a
       real(kind=cp), dimension(4):: b
       real(kind=cp)              :: c
    End Type Xray_Form_Type

    !!----
    !!---- XRAY_FORM
    !!----    Type (Xray_Form_Type), allocatable, dimension(:), public :: Xray_Form
    !!----
    !!----    Tabulated Xray scattering factor coefficients
    !!----
    !!---- Update: February - 2005
    !!
    Type(Xray_Form_Type), allocatable, dimension(:), public :: Xray_Form

    !!----
    !!---- TYPE :: XRAY_WAVELENGTH_TYPE
    !!--..
    !!---- Type, public :: Xray_Wavelength_Type
    !!----    character (len= 2)                :: Symb  ! Symbol of the Chemical species
    !!----    real(kind=cp), dimension(2)       :: Kalfa ! K-Serie for X-ray
    !!----    real(kind=cp)                     :: Kbeta ! K-Serie for X-ray
    !!---- End Type Xray_Wavelength_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Xray_Wavelength_Type
       character (len= 2)         :: Symb
       real(kind=cp), dimension(2):: Kalfa
       real(kind=cp)              :: Kbeta
    End Type Xray_Wavelength_Type

    !!----
    !!---- XRAY_WAVELENGTHS
    !!----    Type (Xray_Wavelength_Type), dimension(7), public :: Xray_Wavelengths
    !!----
    !!----    Tabulated K-Series for Xray
    !!----
    !!---- Update: February - 2005
    !!
    Type(Xray_Wavelength_Type), dimension(7), public :: Xray_Wavelengths =(/                            &
                                                Xray_Wavelength_type("CR",(/2.28988,2.29428/),2.08480), &
                                                Xray_Wavelength_type("FE",(/1.93631,1.94043/),1.75650), &
                                                Xray_Wavelength_type("CU",(/1.54059,1.54431/),1.39220), &
                                                Xray_Wavelength_type("MO",(/0.70932,0.71360/),0.63225), &
                                                Xray_Wavelength_type("AG",(/0.55942,0.56380/),0.49708), &
                                                Xray_Wavelength_type("CO",(/1.78919,1.79321/),1.62083), &
                                                Xray_Wavelength_type("NI",(/1.65805,1.66199/),1.50017)  /)

 Contains

    !---------------------!
    !---- Subroutines ----!
    !---------------------!

    !!----
    !!---- Subroutine Get_Atomic_Mass(Atm,Mass)
    !!----    character(len=2), intent(in)  :: Atm
    !!----    real(kind=cp),    intent(out) :: Mass
    !!----
    !!----    Provides the atomic mass given the chemical symbol of the element
    !!----    In case of problems the returned mass is ZERO.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Atomic_Mass(atm,mass)
       !---- Arguments ----!
       character(len=2), intent (in) :: atm
       real(kind=cp),    intent(out) :: Mass

       !---- Local variables ----!
       character(len=2) :: atm_car
       integer :: i

       mass=0.0
       atm_car=u_case(atm)
       if (.not. allocated(chem_info) ) call set_chem_info()

       do i=1,Num_Chem_Info
          if (index(atm_car,chem_info(i)%Symb) /=0) then
             mass=chem_info(i)%AtWe
             exit
          end if
       end do

       return
    End Subroutine Get_Atomic_Mass

    !!----
    !!---- Subroutine Get_Atomic_Vol(Atm,Vol)
    !!----    character(len=2), intent(in)  :: Atm
    !!----    real(kind=cp),    intent(out) :: Vol
    !!----
    !!----    Provides the atomic volume given the chemical symbol of the element
    !!----    In case of problems the returned Volume is ZERO.
    !!----
    !!---- Update: March- 2013
    !!
    Subroutine Get_Atomic_Vol(atm,vol)
       !---- Arguments ----!
       character(len=2), intent (in) :: atm
       real(kind=cp),    intent(out) :: Vol

       !---- Local variables ----!
       character(len=2) :: atm_car
       integer :: i

       vol=0.0
       atm_car=u_case(atm)
       if (.not. allocated(chem_info) ) call set_chem_info()

       do i=1,Num_Chem_Info
          if (index(atm_car,chem_info(i)%Symb) /=0) then
             vol=chem_info(i)%VAtm
             exit
          end if
       end do

       return
    End Subroutine Get_Atomic_Vol

    !!----
    !!---- Subroutine Get_ChemSymb(Label, ChemSymb, Z)
    !!----   character(len=*),  intent(in) :: Label    ! Label
    !!----   character(len=*),  intent(out):: ChemSymb ! Chemical Symbol
    !!----   integer, optional, intent(out):: Z        ! Atomic number
    !!----
    !!----  Subroutine to get the chemical symbol from label and optionally
    !!----  the atomic number
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_ChemSymb(Label, ChemSymb, Z)
       !---- Argument ----!
       character(len=*),  intent(in) :: Label    ! Label
       character(len=*),  intent(out):: ChemSymb ! Chemical Symbol
       integer, optional, intent(out):: Z        ! Atomic number

       !---- Local variables ----!
       character(len=*),  parameter :: parcar="1234567890+-."
       character(len=2)             :: car
       integer                      :: npos

       ChemSymb="**"
       car=adjustl(label)
       npos=index(parcar,car(2:2))
       if (npos /=0) car(2:2)=" "
       car=u_case(car)
       car(2:2)=l_case(car(2:2))
       ChemSymb=car

       if (present(z)) then
          if (.not. allocated(chem_info) ) call set_chem_info()
          car=u_case(chemsymb)
          do npos=1,num_chem_info
             if (car == Chem_Info(npos)%Symb) then
                Z=Chem_Info(npos)%Z
                exit
             end if
          end do
       end if

       return
    End Subroutine Get_ChemSymb

    !!----
    !!---- Subroutine Get_Covalent_Radius(nam,rad)
    !!----    character(len=*), intent (in) :: nam
    !!----    real(kind=cp),    intent(out) :: rad
    !!----
    !!----    Provides the covalent radius given the chemical symbol of the element
    !!----    In case of problems the returned radius is 1.4 angstroms.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Covalent_Radius(nam,rad)
       !---- Arguments ----!
       character(len=*), intent (in) :: nam
       real(kind=cp),    intent(out) :: rad

       !---- Local variables ----!
       character(len=2) :: atm_car
       integer          :: i

       rad=1.4
       atm_car=u_case(nam(1:2))
       if (atm_car(2:2) > "Z" .or. atm_car(2:2) < "A") atm_car(2:2)=" "
       if (.not. allocated(chem_info) ) call set_chem_info()
       do i=1,Num_Chem_Info
          if (index(atm_car,chem_info(i)%Symb) /=0) then
             rad=chem_info(i)%RCov
             exit
          end if
       end do

       return
    End Subroutine Get_Covalent_Radius

    !!----
    !!---- Subroutine Get_Fermi_Length(nam,b)
    !!----    character(len=*), intent (in) :: nam
    !!----    real(kind=cp),    intent(out) :: b
    !!----
    !!----    Provides the Fermi length (in 10-12 cm) given the chemical
    !!----    symbol of the element. In case of problems the returned Fermi
    !!----    length is 0.0 10-12 cm.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Fermi_Length(nam,b)
       !---- Arguments ----!
       character(len=*), intent (in) :: nam
       real(kind=cp),    intent(out) :: b

       !---- Local variables ----!
       character(len=2) :: atm_car
       integer          :: i

       b=0.0
       atm_car=u_case(nam(1:2))
       if (atm_car(2:2) > "Z" .or. atm_car(2:2) < "A") atm_car(2:2)=" "
       if (.not. allocated(chem_info) ) call set_chem_info()
       do i=1,Num_Chem_Info
          if (index(atm_car,chem_info(i)%Symb) /=0) then
             b=chem_info(i)%SctF
             exit
          end if
       end do

       return
    End Subroutine Get_Fermi_Length

    !!----
    !!---- Subroutine Get_Inc_Xs(nam,u)
    !!----    character(len=*), intent (in) :: nam
    !!----    real(kind=cp),    intent(out) :: u
    !!----
    !!----    Provides incoherent scattering neutron cross-section (barns -> [10**(-24) cm**2] )
    !!----    for given chemical symbol of the element. In case of problems the returned value is 0.0.
    !!----
    !!----
    !!---- Update: Mai - 2013
    !!

    Subroutine Get_Inc_Xs(nam,u)
       !---- Arguments ----!
       character(len=*), intent (in) :: nam
       real(kind=cp),    intent(out) :: u

       !---- Local variables ----!
       character(len=2) :: atm_car
       integer          :: i

       u=0.0
       atm_car=u_case(nam(1:2))
       if (atm_car(2:2) > "Z" .or. atm_car(2:2) < "A") atm_car(2:2)=" "
       if (.not. allocated(chem_info) ) call set_chem_info()
       do i=1,Num_Chem_Info
          if (index(atm_car,chem_info(i)%Symb) /=0) then
             u=chem_info(i)%SedInc
             exit
          end if
       end do

       return
    End Subroutine Get_Inc_Xs

    !!----
    !!---- Subroutine Get_Abs_Xs(nam,u)
    !!----    character(len=*), intent (in) :: nam
    !!----    real(kind=cp),    intent(out) :: u
    !!----
    !!----    Provides the absorption cross-section ( barns, for v= 2200m/s, l(A)=3.95/v (km/s) )
    !!----    for given chemical symbol of the element. In case of problems the returned value is 0.0.
    !!----
    !!----
    !!---- Update: April - 2013
    !!

    Subroutine Get_Abs_Xs(nam,u)
       !---- Arguments ----!
       character(len=*), intent (in) :: nam
       real(kind=cp),    intent(out) :: u

       !---- Local variables ----!
       character(len=2) :: atm_car
       integer          :: i

       u=0.0
       atm_car=u_case(nam(1:2))
       if (atm_car(2:2) > "Z" .or. atm_car(2:2) < "A") atm_car(2:2)=" "
       if (.not. allocated(chem_info) ) call set_chem_info()
       do i=1,Num_Chem_Info
          if (index(atm_car,chem_info(i)%Symb) /=0) then
             u=chem_info(i)%Sea
             exit
          end if
       end do

       return
    End Subroutine Get_Abs_Xs

    !!----
    !!---- Subroutine Get_Ionic_Radius(nam,valence,rad)
    !!----    character(len=*), intent (in) :: nam
    !!----    integer,          intent (in) :: valence
    !!----    real(kind=cp),    intent(out) :: rad
    !!----
    !!----    Provides the ionic radius given the chemical symbol of the element
    !!----    and the valence as an integer. In case of problems the returned radius is 0.0 angstroms.
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Get_Ionic_Radius(nam,valence,rad)
       !---- Arguments ----!
       character(len=*), intent (in) :: nam
       integer,          intent (in) :: valence
       real(kind=cp),    intent(out) :: rad

       !---- Local variables ----!
       character(len=2) :: atm_car
       integer          :: i,j

       rad=0.0
       atm_car=u_case(nam(1:2))
       if (atm_car(2:2) > "Z" .or. atm_car(2:2) < "A") atm_car(2:2)=" "
       if (.not. allocated(chem_info) ) call set_chem_info()
       do i=1,Num_Chem_Info
          if (index(atm_car,chem_info(i)%Symb) /=0) then
             do j=1,5
                if (valence == chem_info(i)%oxid(j)) then
                   rad=chem_info(i)%Rion(j)
                   exit
                end if
             end do
          end if
       end do

       return
    End Subroutine Get_Ionic_Radius

    !!----
    !!---- Subroutine Remove_Chem_Info()
    !!----
    !!----    Deallocate Chem_Info Table
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Remove_Chem_Info()

       if (allocated(chem_info)) deallocate(chem_info)

       return
    End Subroutine Remove_Chem_Info

    !!----
    !!---- Subroutine Remove_Delta_Fp_Fpp()
    !!----
    !!----    Deallocate Anomalous_ScFac Table
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Remove_Delta_Fp_Fpp()

       if (allocated(Anomalous_ScFac)) deallocate(Anomalous_ScFac)

       return
    End Subroutine Remove_Delta_Fp_Fpp

    !!----
    !!---- Subroutine Remove_Magnetic_Form()
    !!----
    !!----    Deallocate Magnetic_Form Table
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Remove_Magnetic_Form()

       if (allocated(Magnetic_Form)) deallocate(Magnetic_Form)
       if (allocated(Magnetic_j2))   deallocate(Magnetic_j2)
       if (allocated(Magnetic_j4))   deallocate(Magnetic_j4)
       if (allocated(Magnetic_j6))   deallocate(Magnetic_j6)

       return
    End Subroutine Remove_Magnetic_form

    !!----
    !!---- Subroutine Remove_Xray_Form()
    !!----
    !!----    Deallocate Xray_Form Table
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Remove_Xray_Form()

       if (allocated(Xray_Form)) deallocate(Xray_Form)

       return
    End Subroutine Remove_Xray_form

    !!----
    !!---- Subroutine Set_Chem_Info()
    !!----    Allocates and loads the table  chem_info(num_chem_info):
    !!--<<
    !!----        1: Symbol of the Element
    !!----        2: Name of the Element
    !!----        3: Atomic Number
    !!----        4: Atomic weight
    !!----        5: Covalent Radius
    !!----        6: Van der Waals Radius
    !!----        7: Atomic volumen
    !!----        8: Oxidation State (5 states)
    !!----        9: Ionic Radius (depending of the oxidation)
    !!----       10: Fermi lenght [10**(-12) cm]
    !!----       11: Incoherent Scattering Neutron cross-section (barns -> [10**(-24) cm**2] )
    !!----       12: Neutron Absorption cross-section ( barns, for v= 2200m/s, l(A)=3.95/v (km/s) )
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Chem_Info()

       if (.not. allocated(chem_info)) allocate(chem_info(num_chem_info))

       !  Symb , Name, Z , AtWe  , RCov , RWaals, VAtm, Oxid(5), Rion(5), b=SctF, SedInc, Sea
       chem_info( 1:10) = (/  &
                          chem_info_type("H ","Hydrogen    ",  1,  1.00797, 0.320, 1.33, 14.1, (/-1, 1, 0, 0, 0/)  ,  &
                                                           (/ 2.08, 0.00, 0.00, 0.00, 0.00/),-0.3739,80.2600,  0.33260    ) ,  &
                          chem_info_type("HE","Helium      ",  2,  4.00260, 0.930, 1.50, 31.8, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.3260, 0.0000,  0.00747    ) ,  &
                          chem_info_type("LI","Lithium     ",  3,  6.94100, 1.230, 1.78, 13.1, (/ 1, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.60, 0.00, 0.00, 0.00, 0.00/),-0.1900, 0.9200, 70.50000    ) ,  &
                          chem_info_type("BE","Beryllium   ",  4,  9.01218, 0.900, 1.45,  5.0, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.31, 0.00, 0.00, 0.00, 0.00/), 0.7790, 0.0018,  0.00760    ) ,  &
                          chem_info_type("B ","Boron       ",  5, 10.81000, 0.820, 1.93,  4.6, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.20, 0.00, 0.00, 0.00, 0.00/), 0.5300, 1.7000,767.00000    ) ,  &
                          chem_info_type("C ","Carbon      ",  6, 12.01100, 0.770, 1.70,  5.3, (/ 2,-4, 4, 0, 0/)  ,  &
                                                           (/ 2.60, 0.15, 0.00, 0.00, 0.00/), 0.6646, 0.0010,  0.0035     ) ,  &
                          chem_info_type("N ","Nitrogen    ",  7, 14.00670, 0.750, 1.70, 17.3, (/ 2,-3, 3, 4, 5/)  ,  &
                                                           (/ 0.00, 1.71, 0.00, 0.00, 0.11/), 0.9360, 0.5000,  1.9000     ) ,  &
                          chem_info_type("O ","Oxygen      ",  8, 15.99940, 0.730, 1.50, 14.0, (/-2, 6, 0, 0, 0/)  ,  &
                                                           (/ 1.40, 0.09, 0.00, 0.00, 0.00/), 0.5803, 0.0000,  0.00019    ) ,  &
                          chem_info_type("F ","Fluorine    ",  9, 18.99840, 0.720, 1.47, 17.1, (/-1, 7, 0, 0, 0/)  ,  &
                                                           (/ 1.36, 0.07, 0.00, 0.00, 0.00/), 0.5654, 0.0008,  0.0096     ) ,  &
                          chem_info_type("NE","Neon        ", 10, 20.17900, 0.710, 1.50, 16.8, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.4566, 0.0080,  0.0039     ) /)

       chem_info(11:20) = (/  &
                          chem_info_type("NA","Sodium      ", 11, 22.98977, 1.540, 2.07, 23.7, (/ 1, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.95, 0.00, 0.00, 0.00, 0.00/), 0.3630, 1.6200,  0.5300     ) ,  &
                          chem_info_type("MG","Magnesium   ", 12, 24.30500, 1.360, 2.20, 14.0, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.65, 0.00, 0.00, 0.00, 0.00/), 0.5375, 0.0800,  0.063      ) ,  &
                          chem_info_type("AL","Aluminum    ", 13, 26.98154, 1.180, 2.45, 10.0, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.50, 0.00, 0.00, 0.00, 0.00/), 0.3449, 0.0082,  0.231      ) ,  &
                          chem_info_type("SI","Silicon     ", 14, 28.08600, 1.310, 2.30, 12.1, (/-1, 4, 0, 0, 0/)  ,  &
                                                           (/ 2.71, 0.41, 0.00, 0.00, 0.00/), 0.4149, 0.0040,  0.171       ) , &
                          chem_info_type("P ","Phosphorus  ", 15, 30.97376, 1.060, 2.15, 17.0, (/-3, 3, 4, 5, 0/)  ,  &
                                                           (/ 2.12, 0.00, 0.00, 0.34, 0.00/), 0.5130, 0.0050,  0.172      ) ,  &
                          chem_info_type("S ","Sulfur      ", 16, 32.06000, 1.020, 1.74, 15.5, (/-2, 2, 4, 6, 0/)  ,  &
                                                           (/ 1.84, 0.29, 0.00, 0.00, 0.00/), 0.2847, 0.0070,  0.530      ) ,  &
                          chem_info_type("CL","Chlorine    ", 17, 35.45300, 0.990, 1.76, 18.7, (/-1, 1, 3, 5, 7/)  ,  &
                                                           (/ 1.81, 0.00, 0.00, 0.00, 0.26/), 0.9577, 5.3000, 33.500      ) ,  &
                          chem_info_type("AR","Argon       ", 18, 39.94800, 0.980, 2.00, 24.2, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.1909, 0.2250,  0.675      ) ,  &
                          chem_info_type("K ","Potassium   ", 19, 39.09800, 2.030, 2.43, 45.3, (/ 1, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.33, 0.00, 0.00, 0.00, 0.00/), 0.3670, 0.2700,  2.100      ) ,  &
                          chem_info_type("CA","Calcium     ", 20, 40.08000, 1.740, 2.09, 29.9, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.99, 0.00, 0.00, 0.00, 0.00/), 0.4700, 0.0500,  0.430      ) /)

       chem_info(21:30) = (/  &
                          chem_info_type("SC","Scandium    ", 21, 44.95590, 1.440, 2.54, 15.0, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.81, 0.00, 0.00, 0.00, 0.00/), 1.2290, 4.5000, 27.500      ) ,  &
                          chem_info_type("TI","Titanium    ", 22, 47.8670, 1.320, 2.57, 10.6, (/ 2, 4, 0, 0, 0/)  ,  &
                                                           (/ 0.90, 0.68, 0.00, 0.00, 0.00/),-0.3438, 2.8700,  6.090      ) ,  &
                          chem_info_type("V ","Vanadium    ", 23, 50.94140, 1.320, 2.43,  8.3, (/ 2, 3, 4, 5, 0/)  ,  &
                                                           (/ 0.00, 0.74, 0.00, 0.59, 0.00/),-0.0382, 5.0800,  5.080      ) ,  &
                          chem_info_type("CR","Chromium    ", 24, 51.99600, 1.180, 2.45,  7.2, (/ 2, 3, 6, 0, 0/)  ,  &
                                                           (/ 0.00, 0.69, 0.52, 0.00, 0.00/), 0.3635, 1.8370,  3.050      ) ,  &
                          chem_info_type("MN","Manganese   ", 25, 54.93800, 1.170, 2.45,  7.4, (/ 2, 3, 4, 6, 7/)  ,  &
                                                           (/ 0.80, 0.72, 0.53, 0.46, 0.46/),-0.3730, 0.4000, 13.300      ) ,  &
                          chem_info_type("FE","Iron        ", 26, 55.84700, 1.170, 2.44,  7.1, (/ 2, 3, 0, 0, 0/)  ,  &
                                                           (/ 0.76, 0.64, 0.00, 0.00, 0.00/), 0.9450, 0.4000,  2.560      ) ,  &
                          chem_info_type("CO","Cobalt      ", 27, 58.93320, 1.160, 2.43,  6.7, (/ 2, 3, 0, 0, 0/)  ,  &
                                                           (/ 0.74, 0.63, 0.00, 0.00, 0.00/), 0.2490, 4.8000, 37.180      ) ,  &
                          chem_info_type("NI","Nickel      ", 28, 58.70000, 1.160, 2.60,  6.6, (/ 2, 3, 0, 0, 0/)  ,  &
                                                           (/ 0.74, 0.63, 0.00, 0.00, 0.00/), 1.0300, 5.2000,  4.490      ) ,  &
                          chem_info_type("CU","Copper      ", 29, 63.54600, 1.170, 2.62,  7.1, (/ 1, 2, 0, 0, 0/)  ,  &
                                                           (/ 0.96, 0.69, 0.00, 0.00, 0.00/), 0.7718, 0.5500,  3.780      ) ,  &
                          chem_info_type("ZN","Zinc        ", 30, 65.38000, 1.250, 2.55,  9.2, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.74, 0.00, 0.00, 0.00, 0.00/), 0.5680, 0.0770,  1.110      ) /)

       chem_info(31:40) = (/  &
                          chem_info_type("GA","Gallium     ", 31, 69.72000, 1.260, 2.32, 11.8, (/ 1, 3, 0, 0, 0/)  ,  &
                                                           (/ 1.13, 0.62, 0.00, 0.00, 0.00/), 0.7288, 0.1600,  2.750      ) ,  &
                          chem_info_type("GE","Germanium   ", 32, 72.59000, 1.220, 2.27, 13.6, (/ 2, 4, 0, 0, 0/)  ,  &
                                                           (/ 0.93, 0.53, 0.00, 0.00, 0.00/), 0.8185, 0.1700,  2.200      ) ,  &
                          chem_info_type("AS","Arsenic     ", 33, 74.92160, 1.200, 2.11, 13.1, (/-3, 3, 5, 0, 0/)  ,  &
                                                           (/ 2.22, 0.00, 0.47, 0.00, 0.00/), 0.6580, 0.0600,  4.500      ) ,  &
                          chem_info_type("SE","Selenium    ", 34, 78.96000, 1.160, 2.32, 16.5, (/-2, 2, 4, 6, 0/)  ,  &
                                                           (/ 1.98, 0.00, 0.00, 0.42, 0.00/), 0.7970, 0.3200, 11.700      ) ,  &
                          chem_info_type("BR","Bromine     ", 35, 79.90400, 1.140, 1.85, 23.5, (/-1, 1, 3, 5, 7/)  ,  &
                                                           (/ 1.95, 0.00, 0.00, 0.00, 0.39/), 0.6795, 0.1000,  6.900      ) ,  &
                          chem_info_type("KR","Krypton     ", 36, 83.80000, 1.120, 2.50, 32.2, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.7810, 0.0100, 25.000      ) ,  &
                          chem_info_type("RB","Rubidium    ", 37, 85.46780, 2.160, 2.57, 55.9, (/ 1, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.48, 0.00, 0.00, 0.00, 0.00/), 0.7090, 0.5000,  0.380      ) ,  &
                          chem_info_type("SR","Strontium   ", 38, 87.62000, 1.910, 2.22, 33.7, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.13, 0.00, 0.00, 0.00, 0.00/), 0.7020, 0.0600,  1.280      ) ,  &
                          chem_info_type("Y ","Yttrium     ", 39, 88.90590, 1.620, 2.88, 19.8, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.93, 0.00, 0.00, 0.00, 0.00/), 0.7750, 0.1500,  1.280      ) ,  &
                          chem_info_type("ZR","Zirconium   ", 40, 91.22000, 1.450, 2.66, 14.1, (/ 4, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.80, 0.00, 0.00, 0.00, 0.00/), 0.7160, 0.0200,  0.185      ) /)

       chem_info(41:50) = (/  &
                          chem_info_type("NB","Niobium     ", 41, 92.90640, 1.340, 2.58, 10.8, (/ 3, 5, 0, 0, 0/)  ,  &
                                                           (/ 0.70, 0.00, 0.00, 0.00, 0.00/), 0.7054, 0.0024,  1.150      ) ,  &
                          chem_info_type("MO","Molybdenum  ", 42, 95.94000, 1.300, 2.57,  9.4, (/ 2, 3, 4, 5, 6/)  ,  &
                                                           (/ 0.00, 0.00, 0.68, 0.00, 0.62/), 0.6715, 0.0400,  2.480      ) ,  &
                          chem_info_type("TC","Technetium  ", 43, 97.00000, 1.270, 2.45,  0.0, (/ 7, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.6800, 0.5000, 20.000      ) ,  &
                          chem_info_type("RU","Ruthenium   ", 44,101.07000, 1.250, 2.50,  8.3, (/ 2, 3, 4, 6, 8/)  ,  &
                                                           (/ 0.00, 0.69, 0.67, 0.00, 0.00/), 0.7030, 0.4000,  2.560      ) ,  &
                          chem_info_type("RH","Rhodium     ", 45,102.90550, 1.250, 2.55,  8.3, (/ 2, 3, 4, 0, 0/)  ,  &
                                                           (/ 0.86, 0.00, 0.00, 0.00, 0.00/), 0.5880, 0.3000,144.800      ) ,  &
                          chem_info_type("PD","Palladium   ", 46,106.40000, 1.280, 2.60,  8.9, (/ 2, 4, 0, 0, 0/)  ,  &
                                                           (/ 0.86, 0.00, 0.00, 0.00, 0.00/), 0.5910, 0.0930,  6.900      ) ,  &
                          chem_info_type("AG","Silver      ", 47,107.86800, 1.340, 2.69, 10.3, (/ 1, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.26, 0.00, 0.00, 0.00, 0.00/), 0.5922, 0.5800, 63.300      ) ,  &
                          chem_info_type("CD","Cadmium     ", 48,112.40000, 1.480, 2.79, 13.1, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.97, 0.00, 0.00, 0.00, 0.00/), 0.4870, 3.4600,2520.00      ) ,  &
                          chem_info_type("IN","Indium      ", 49,114.82000, 1.440, 2.73, 15.7, (/ 1, 3, 0, 0, 0/)  ,  &
                                                           (/ 1.32, 0.81, 0.00, 0.00, 0.00/), 0.4065, 0.5400,193.800      ) ,  &
                          chem_info_type("SN","Tin         ", 50,118.69000, 1.410, 2.56, 16.3, (/ 2, 4, 0, 0, 0/)  ,  &
                                                           (/ 1.12, 0.71, 0.00, 0.00, 0.00/), 0.6225, 0.0220,  0.626      ) /)

       chem_info(51:60) = (/  &
                          chem_info_type("SB","Antimony    ", 51,121.75000, 1.400, 2.56, 18.4, (/-3, 3, 5, 0, 0/)  ,  &
                                                           (/ 2.45, 0.00, 0.62, 0.00, 0.00/), 0.5570, 0.0000,  4.910      ) ,  &
                          chem_info_type("TE","Tellurium   ", 52,127.60000, 1.360, 2.57, 20.5, (/-2, 2, 4, 6, 0/)  ,  &
                                                           (/ 2.21, 0.00, 0.00, 0.56, 0.00/), 0.5800, 0.0900,  4.700      ) ,  &
                          chem_info_type("I ","Iodine      ", 53,126.90450, 1.330, 1.98, 25.7, (/-1, 1, 3, 5, 7/)  ,  &
                                                           (/ 2.16, 0.00, 0.00, 0.00, 0.50/), 0.5280, 0.3100,  6.150      ) ,  &
                          chem_info_type("XE","Xenon       ", 54,131.30000, 1.310, 2.50, 42.9, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.4920, 0.0000, 23.900      ) ,  &
                          chem_info_type("CS","Cesium      ", 55,132.90540, 2.350, 2.77, 70.0, (/ 1, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.69, 0.00, 0.00, 0.00, 0.00/), 0.5420, 0.2100, 29.000      ) ,  &
                          chem_info_type("BA","Barium      ", 56,137.34000, 1.980, 2.44, 39.0, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.35, 0.00, 0.00, 0.00, 0.00/), 0.5070, 0.1500,  1.100      ) ,  &
                          chem_info_type("LA","Lanthanum   ", 57,138.90550, 1.690, 2.97, 22.5, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.15, 0.00, 0.00, 0.00, 0.00/), 0.8240, 1.1300,  8.970      ) ,  &
                          chem_info_type("CE","Cerium      ", 58,140.12000, 1.650, 2.93, 21.0, (/ 3, 4, 0, 0, 0/)  ,  &
                                                           (/ 1.11, 1.01, 0.00, 0.00, 0.00/), 0.4840, 0.0000,  0.630      ) ,  &
                          chem_info_type("PR","Praseodymium", 59,140.90770, 1.650, 2.92, 20.8, (/ 3, 4, 0, 0, 0/)  ,  &
                                                           (/ 1.09, 0.92, 0.00, 0.00, 0.00/), 0.4580, 0.0150, 11.500      ) ,  &
                          chem_info_type("ND","Neodymium   ", 60,144.24000, 1.640, 2.91, 20.6, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.08, 0.00, 0.00, 0.00, 0.00/), 0.7690, 9.2000, 50.500      ) /)

       chem_info(61:70) = (/  &
                          chem_info_type("PM","Promethium  ", 61,145.00000, 1.630, 2.90,  0.0, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.06, 0.00, 0.00, 0.00, 0.00/), 1.2600, 1.3000,168.400      ) ,  &
                          chem_info_type("SM","Samarium    ", 62,150.40000, 1.620, 2.90, 19.9, (/ 2, 3, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 1.04, 0.00, 0.00, 0.00/), 0.8000,39.0000,5922.00      ) ,  &
                          chem_info_type("EU","Europium    ", 63,151.96000, 1.850, 2.90, 28.9, (/ 2, 3, 0, 0, 0/)  ,  &
                                                           (/ 1.12, 0.00, 0.00, 0.00, 0.00/), 0.7220, 2.5000,4530.00      ) ,  &
                          chem_info_type("GD","Gadolinium  ", 64,157.25000, 1.610, 2.89, 19.9, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.02, 0.00, 0.00, 0.00, 0.00/), 0.6500, 0.0000,49700.0      ) ,  &
                          chem_info_type("TB","Terbium     ", 65,158.92540, 1.590, 2.86, 19.2, (/ 3, 4, 0, 0, 0/)  ,  &
                                                           (/ 1.00, 0.00, 0.00, 0.00, 0.00/), 0.7380, 0.0040, 23.400      ) ,  &
                          chem_info_type("DY","Dysprosium  ", 66,162.50000, 1.590, 2.85, 19.0, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.99, 0.00, 0.00, 0.00, 0.00/), 1.6900,54.4000,994.000      ) ,  &
                          chem_info_type("HO","Holmium     ", 67,164.93040, 1.580, 2.84, 18.7, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.97, 0.00, 0.00, 0.00, 0.00/), 0.8010, 0.3600, 64.700      ) ,  &
                          chem_info_type("ER","Erbium      ", 68,167.26000, 1.570, 2.83, 18.4, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.96, 0.00, 0.00, 0.00, 0.00/), 0.7790, 1.1000,159.000      ) ,  &
                          chem_info_type("TM","Thulium     ", 69,168.93420, 1.560, 2.82, 18.1, (/ 2, 3, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.95, 0.00, 0.00, 0.00/), 0.7070, 0.1000,100.000      ) ,  &
                          chem_info_type("YB","Ytterbium   ", 70,173.04000, 0.000, 3.04, 24.8, (/ 2, 3, 0, 0, 0/)  ,  &
                                                           (/ 1.13, 0.94, 0.00, 0.00, 0.00/), 1.2430, 4.0000, 34.800      ) /)

       chem_info(71:80) = (/  &
                          chem_info_type("LU","Lutetium    ", 71,174.97000, 1.560, 2.82, 17.8, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.93, 0.00, 0.00, 0.00, 0.00/), 0.7210, 0.7000, 74.000      ) ,  &
                          chem_info_type("HF","Hafnium     ", 72,178.49000, 1.440, 2.67, 13.6, (/ 4, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.81, 0.00, 0.00, 0.00, 0.00/), 0.7770, 2.6000, 74.000      ) ,  &
                          chem_info_type("TA","Tantalum    ", 73,180.94790, 1.440, 2.53, 10.9, (/ 5, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.73, 0.00, 0.00, 0.00, 0.00/), 0.6910, 0.0100, 20.600      ) ,  &
                          chem_info_type("W ","Tungsten    ", 74,183.85000, 1.300, 2.47,  9.5, (/ 2, 3, 4, 5, 6/)  ,  &
                                                           (/ 0.00, 0.00, 0.68, 0.00, 0.64/), 0.4860, 1.6300, 18.300      ) ,  &
                          chem_info_type("RE","Rhenium     ", 75,186.20700, 1.280, 2.45,  8.8, (/ 1, 2, 4, 6, 7/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.9200, 0.9000, 89.700      ) ,  &
                          chem_info_type("OS","Osmium      ", 76,190.20000, 1.260, 2.47,  8.4, (/ 2, 3, 4, 6, 8/)  ,  &
                                                           (/ 0.00, 0.00, 0.69, 0.00, 0.00/), 1.0700, 0.3000, 16.000      ) ,  &
                          chem_info_type("IR","Iridium     ", 77,192.22000, 1.270, 2.42,  8.5, (/ 2, 3, 4, 6, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.66, 0.00, 0.00/), 1.0600, 0.0000,425.000      ) ,  &
                          chem_info_type("PT","Platinum    ", 78,195.09000, 1.300, 2.60,  9.1, (/ 2, 4, 0, 0, 0/)  ,  &
                                                           (/ 0.96, 0.00, 0.00, 0.00, 0.00/), 0.9600, 0.1300, 10.300      ) ,  &
                          chem_info_type("AU","Gold        ", 79,196.96650, 1.340, 2.60, 10.2, (/ 1, 3, 0, 0, 0/)  ,  &
                                                           (/ 1.37, 0.00, 0.00, 0.00, 0.00/), 0.7630, 0.4300, 98.650      ) ,  &
                          chem_info_type("HG","Mercury     ", 80,200.59000, 1.490, 2.80, 14.8, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.10, 0.00, 0.00, 0.00, 0.00/), 1.2692, 6.6000,372.300      ) /)

       chem_info(81:90) = (/  &
                          chem_info_type("TL","Thallium    ", 81,204.37000, 1.480, 2.65, 17.2, (/ 1, 3, 0, 0, 0/)  ,  &
                                                           (/ 1.40, 0.95, 0.00, 0.00, 0.00/), 0.8776, 0.2100,  3.430      ) ,  &
                          chem_info_type("PB","Lead        ", 82,207.20000, 1.470, 2.64, 18.3, (/ 2, 4, 0, 0, 0/)  ,  &
                                                           (/ 1.20, 0.84, 0.00, 0.00, 0.00/), 0.9405, 0.0030,  0.171      ) ,  &
                          chem_info_type("BI","Bismuth     ", 83,208.98040, 1.460, 2.64, 21.3, (/ 3, 5, 0, 0, 0/)  ,  &
                                                           (/ 1.20, 0.74, 0.00, 0.00, 0.00/), 0.8532, 0.0084,  0.034      ) ,  &
                          chem_info_type("PO","Polonium    ", 84,209.00000, 1.460, 2.60, 22.7, (/ 2, 4, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("AT","Astatine    ", 85,210.00000, 0.000, 2.60,  0.0, (/-1, 1, 3, 5, 7/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("RN","Radon       ", 86,222.00000, 0.000, 2.60,  0.0, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("FR","Francium    ", 87,223.00000, 0.000, 3.00,  0.0, (/ 1, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.76, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("RA","Radium      ", 88,226.02540, 0.000, 3.00, 45.0, (/ 2, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.40, 0.00, 0.00, 0.00, 0.00/), 1.0000, 0.0000, 12.800      ) ,  &
                          chem_info_type("AC","Actinium    ", 89,227.00000, 0.000, 2.98,  0.0, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 1.18, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("TH","Thorium     ", 90,232.03810, 1.650, 2.89, 19.9, (/ 3, 4, 0, 0, 0/)  ,  &
                                                           (/ 1.14, 0.95, 0.00, 0.00, 0.00/), 1.0310, 0.0000,  7.370      ) /)

       chem_info(91:100)= (/  &
                          chem_info_type("PA","Protactinium", 91,231.03597, 0.000, 2.71, 15.0, (/ 3, 4, 5, 0, 0/)  ,  &
                                                           (/ 1.12, 0.98, 0.00, 0.00, 0.00/), 0.9100, 0.1000,200.600      ) ,  &
                          chem_info_type("U ","Uranium     ", 92,238.02900, 1.420, 2.68, 12.5, (/ 3, 4, 5, 6, 0/)  ,  &
                                                           (/ 1.11, 0.97, 0.00, 0.00, 0.00/), 0.8417, 0.0050,  7.570      ) ,  &
                          chem_info_type("NP","Neptunium   ", 93,237.04820, 0.000, 2.65, 21.1, (/ 3, 4, 5, 6, 0/)  ,  &
                                                           (/ 1.09, 0.95, 0.00, 0.00, 0.00/), 1.0550, 0.5000,175.900      ) ,  &
                          chem_info_type("PU","Plutonium   ", 94,244.00000, 0.000, 2.43,  0.0, (/ 3, 4, 5, 6, 0/)  ,  &
                                                           (/ 1.07, 0.93, 0.00, 0.00, 0.00/), 1.4100, 0.0000,558.000      ) ,  &
                          chem_info_type("AM","Americium   ", 95,243.00000, 0.000, 2.61, 20.8, (/ 3, 4, 5, 6, 0/)  ,  &
                                                           (/ 1.06, 0.92, 0.00, 0.00, 0.00/), 0.8300, 0.3000, 75.300      ) ,  &
                          chem_info_type("CM","Curium      ", 96,247.00000, 0.000, 2.60,  0.0, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.9500, 0.0000, 16.200      ) ,  &
                          chem_info_type("BK","Berkelium   ", 97,247.00000, 0.000, 2.60,  0.0, (/ 3, 4, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("CF","Californium ", 98,251.00000, 0.000, 2.60,  0.0, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("ES","Einsteinium ", 99,254.00000, 0.000, 0.00,  0.0, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("FM","Fermium     ",100,257.00000, 0.000, 0.00,  0.0, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) /)

       chem_info(101:108)=(/  &
                          chem_info_type("MD","Mendelevium ",101,258.00000, 0.000, 0.00,  0.0, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("NO","Nobelium    ",102,255.00000, 0.000, 0.00,  0.0, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("LR","Lawrencium  ",103,260.00000, 0.000, 0.00,  0.0, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("KU","            ",104,261.00000, 0.000, 0.00,  0.0, (/ 4, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("HA","            ",105,262.00000, 0.000, 0.00,  0.0, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.000      ) ,  &
                          chem_info_type("BS","Boron-11    ",  5, 10.81000, 0.820, 1.93,  4.6, (/ 3, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.20, 0.00, 0.00, 0.00, 0.00/), 0.6650, 1.7000,767.00000    ) ,  &
                          chem_info_type("ZE","Zero-scatter",  1,  0.00000, 0.820, 1.93,  4.6, (/ 0, 0, 0, 0, 0/)  ,  &
                                                           (/ 0.00, 0.00, 0.00, 0.00, 0.00/), 0.0000, 0.0000,  0.00000    ) ,  &
                          chem_info_type("D ","Deuterium   ",  1,  2.00797, 0.320, 1.33, 14.1, (/-1, 1, 0, 0, 0/)  ,  &
                                                           (/ 2.08, 0.00, 0.00, 0.00, 0.00/),0.6671, 0.0000,  0.00000    ) /)
       !  Symb , Name, Z , AtWe  , RCov , RWaals, VAtm, Oxid(5), Rion(5), b=SctF, SedInc, Sea
       return
    End Subroutine Set_Chem_Info

    !!----
    !!---- Subroutine Set_Delta_Fp_Fpp()
    !!--<<
    !!----    Wavelenghts:     Cr        Fe        Cu         Mo         Ag
    !!----         Lambda   2.28962   1.93597   1.54051    0.70926    0.556363
    !!-->>
    !!----    Set values for Delta-fp & Delta-fpp for the above wavelengths
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Delta_Fp_Fpp()

       if (.not. allocated(anomalous_ScFac)) allocate(anomalous_ScFac(Num_Delta_Fp))

       Anomalous_ScFac( 1)=Anomalous_Sc_Type("h ", (/   0.000,   0.000,   0.000,   0.000,   0.000/), &
                                                   (/   0.000,   0.000,   0.000,   0.000,   0.000/)  )
       Anomalous_ScFac( 2)=Anomalous_Sc_Type("he", (/   0.000,   0.000,   0.000,   0.000,   0.000/), &
                                                   (/   0.000,   0.000,   0.000,   0.000,   0.000/)  )
       Anomalous_ScFac( 3)=Anomalous_Sc_Type("li", (/   0.002,   0.002,   0.001,   0.000,   0.000/), &
                                                   (/   0.001,   0.001,   0.000,   0.000,   0.000/)  )
       Anomalous_ScFac( 4)=Anomalous_Sc_Type("be", (/   0.008,   0.005,   0.003,   0.000,  -0.001/), &
                                                   (/   0.003,   0.002,   0.001,   0.000,   0.000/)  )
       Anomalous_ScFac( 5)=Anomalous_Sc_Type("b ", (/   0.018,   0.013,   0.008,   0.000,   0.000/), &
                                                   (/   0.009,   0.007,   0.004,   0.001,   0.000/)  )
       Anomalous_ScFac( 6)=Anomalous_Sc_Type("c ", (/   0.035,   0.026,   0.017,   0.002,   0.000/), &
                                                   (/   0.021,   0.015,   0.009,   0.002,   0.001/)  )
       Anomalous_ScFac( 7)=Anomalous_Sc_Type("n ", (/   0.059,   0.044,   0.029,   0.004,   0.001/), &
                                                   (/   0.042,   0.029,   0.018,   0.003,   0.002/)  )
       Anomalous_ScFac( 8)=Anomalous_Sc_Type("o ", (/   0.090,   0.069,   0.047,   0.008,   0.003/), &
                                                   (/   0.073,   0.052,   0.032,   0.006,   0.004/)  )
       Anomalous_ScFac( 9)=Anomalous_Sc_Type("f ", (/   0.129,   0.100,   0.069,   0.014,   0.006/), &
                                                   (/   0.119,   0.085,   0.053,   0.010,   0.006/)  )
       Anomalous_ScFac(10)=Anomalous_Sc_Type("ne", (/   0.174,   0.138,   0.097,   0.021,   0.011/), &
                                                   (/   0.184,   0.132,   0.083,   0.016,   0.010/)  )
       Anomalous_ScFac(11)=Anomalous_Sc_Type("na", (/   0.223,   0.180,   0.129,   0.030,   0.016/), &
                                                   (/   0.270,   0.195,   0.124,   0.025,   0.015/)  )
       Anomalous_ScFac(12)=Anomalous_Sc_Type("mg", (/   0.272,   0.224,   0.165,   0.042,   0.023/), &
                                                   (/   0.381,   0.277,   0.177,   0.036,   0.022/)  )
       Anomalous_ScFac(13)=Anomalous_Sc_Type("al", (/   0.318,   0.269,   0.204,   0.056,   0.032/), &
                                                   (/   0.522,   0.381,   0.246,   0.052,   0.031/)  )
       Anomalous_ScFac(14)=Anomalous_Sc_Type("si", (/   0.355,   0.311,   0.244,   0.072,   0.042/), &
                                                   (/   0.693,   0.509,   0.330,   0.071,   0.043/)  )
       Anomalous_ScFac(15)=Anomalous_Sc_Type("p ", (/   0.377,   0.347,   0.283,   0.090,   0.055/), &
                                                   (/   0.900,   0.664,   0.434,   0.095,   0.058/)  )
       Anomalous_ScFac(16)=Anomalous_Sc_Type("s ", (/   0.374,   0.370,   0.319,   0.110,   0.068/), &
                                                   (/   1.142,   0.847,   0.557,   0.124,   0.076/)  )
       Anomalous_ScFac(17)=Anomalous_Sc_Type("cl", (/   0.335,   0.375,   0.348,   0.132,   0.084/), &
                                                   (/   1.423,   1.061,   0.702,   0.159,   0.099/)  )
       Anomalous_ScFac(18)=Anomalous_Sc_Type("ar", (/   0.243,   0.352,   0.366,   0.155,   0.101/), &
                                                   (/   1.747,   1.309,   0.872,   0.201,   0.125/)  )
       Anomalous_ScFac(19)=Anomalous_Sc_Type("k ", (/   0.070,   0.286,   0.365,   0.179,   0.118/), &
                                                   (/   2.110,   1.589,   1.066,   0.250,   0.156/)  )
       Anomalous_ScFac(20)=Anomalous_Sc_Type("ca", (/  -0.221,   0.163,   0.341,   0.203,   0.137/), &
                                                   (/   2.514,   1.904,   1.286,   0.306,   0.193/)  )
       Anomalous_ScFac(21)=Anomalous_Sc_Type("sc", (/  -0.717,  -0.038,   0.285,   0.226,   0.156/), &
                                                   (/   2.968,   2.256,   1.533,   0.372,   0.235/)  )
       Anomalous_ScFac(22)=Anomalous_Sc_Type("ti", (/  -1.683,  -0.357,   0.189,   0.248,   0.175/), &
                                                   (/   3.470,   2.643,   1.807,   0.446,   0.283/)  )
       Anomalous_ScFac(23)=Anomalous_Sc_Type("v ", (/  -3.841,  -0.896,   0.035,   0.267,   0.194/), &
                                                   (/   0.459,   3.070,   2.110,   0.530,   0.338/)  )
       Anomalous_ScFac(24)=Anomalous_Sc_Type("cr", (/  -2.161,  -1.973,  -0.198,   0.284,   0.213/), &
                                                   (/   0.548,   3.533,   2.443,   0.624,   0.399/)  )
       Anomalous_ScFac(25)=Anomalous_Sc_Type("mn", (/  -1.639,  -3.367,  -0.568,   0.295,   0.229/), &
                                                   (/   0.650,   0.481,   2.808,   0.729,   0.468/)  )
       Anomalous_ScFac(26)=Anomalous_Sc_Type("fe", (/  -1.339,  -2.095,  -1.179,   0.301,   0.244/), &
                                                   (/   0.764,   0.566,   3.204,   0.845,   0.545/)  )
       Anomalous_ScFac(27)=Anomalous_Sc_Type("co", (/  -1.124,  -1.623,  -2.464,   0.299,   0.256/), &
                                                   (/   0.893,   0.662,   3.608,   0.973,   0.630/)  )
       Anomalous_ScFac(28)=Anomalous_Sc_Type("ni", (/  -0.956,  -1.343,  -2.956,   0.285,   0.261/), &
                                                   (/   1.036,   0.769,   0.509,   1.113,   0.724/)  )
       Anomalous_ScFac(29)=Anomalous_Sc_Type("cu", (/  -0.795,  -1.129,  -2.019,   0.263,   0.265/), &
                                                   (/   1.196,   0.888,   0.589,   1.266,   0.826/)  )
       Anomalous_ScFac(30)=Anomalous_Sc_Type("zn", (/  -0.684,  -0.978,  -1.612,   0.222,   0.260/), &
                                                   (/   1.373,   1.021,   0.678,   1.431,   0.938/)  )
       Anomalous_ScFac(31)=Anomalous_Sc_Type("ga", (/  -0.570,  -0.841,  -1.354,   0.163,   0.249/), &
                                                   (/   1.569,   1.168,   0.777,   1.609,   1.059/)  )
       Anomalous_ScFac(32)=Anomalous_Sc_Type("ge", (/  -0.462,  -0.717,  -1.163,   0.081,   0.228/), &
                                                   (/   1.786,   1.331,   0.886,   1.801,   1.190/)  )
       Anomalous_ScFac(33)=Anomalous_Sc_Type("as", (/  -0.365,  -0.607,  -1.011,  -0.030,   0.196/), &
                                                   (/   2.022,   1.508,   1.006,   2.007,   1.332/)  )
       Anomalous_ScFac(34)=Anomalous_Sc_Type("se", (/  -0.273,  -0.503,  -0.879,  -0.178,   0.152/), &
                                                   (/   2.283,   1.704,   1.139,   2.223,   1.481/)  )
       Anomalous_ScFac(35)=Anomalous_Sc_Type("br", (/  -0.198,  -0.413,  -0.767,  -0.374,   0.090/), &
                                                   (/   2.563,   1.916,   1.283,   2.456,   1.643/)  )
       Anomalous_ScFac(36)=Anomalous_Sc_Type("kr", (/  -0.130,  -0.328,  -0.665,  -0.652,   0.008/), &
                                                   (/   2.872,   2.149,   1.439,   2.713,   1.820/)  )
       Anomalous_ScFac(37)=Anomalous_Sc_Type("rb", (/  -0.082,  -0.256,  -0.574,  -1.044,  -0.099/), &
                                                   (/   3.201,   2.398,   1.608,   2.973,   2.003/)  )
       Anomalous_ScFac(38)=Anomalous_Sc_Type("sr", (/  -0.012,  -0.161,  -0.465,  -1.657,  -0.230/), &
                                                   (/   3.608,   2.709,   1.820,   3.264,   2.203/)  )
       Anomalous_ScFac(39)=Anomalous_Sc_Type("y ", (/   0.006,  -0.106,  -0.386,  -2.951,  -0.406/), &
                                                   (/   4.002,   3.009,   2.025,   3.542,   2.411/)  )
       Anomalous_ScFac(40)=Anomalous_Sc_Type("zr", (/   0.007,  -0.061,  -0.314,  -2.965,  -0.639/), &
                                                   (/   4.422,   3.329,   2.245,   0.560,   2.630/)  )
       Anomalous_ScFac(41)=Anomalous_Sc_Type("nb", (/  -0.013,  -0.028,  -0.248,  -2.197,  -0.957/), &
                                                   (/   4.876,   3.676,   2.482,   0.621,   2.860/)  )
       Anomalous_ScFac(42)=Anomalous_Sc_Type("mo", (/  -0.063,  -0.012,  -0.191,  -1.825,  -1.416/), &
                                                   (/   5.353,   4.043,   2.735,   0.688,   3.103/)  )
       Anomalous_ScFac(43)=Anomalous_Sc_Type("tc", (/  -0.153,  -0.017,  -0.145,  -1.590,  -2.205/), &
                                                   (/   5.862,   4.434,   3.005,   0.759,   3.353/)  )
       Anomalous_ScFac(44)=Anomalous_Sc_Type("ru", (/  -0.270,  -0.039,  -0.105,  -1.420,  -5.524/), &
                                                   (/   6.406,   4.854,   3.296,   0.836,   3.651/)  )
       Anomalous_ScFac(45)=Anomalous_Sc_Type("rh", (/  -0.424,  -0.083,  -0.077,  -1.287,  -2.649/), &
                                                   (/   6.984,   5.300,   3.605,   0.919,   0.596/)  )
       Anomalous_ScFac(46)=Anomalous_Sc_Type("pd", (/  -0.639,  -0.157,  -0.059,  -1.177,  -2.128/), &
                                                   (/   7.594,   5.773,   3.934,   1.007,   0.654/)  )
       Anomalous_ScFac(47)=Anomalous_Sc_Type("ag", (/  -0.924,  -0.259,  -0.060,  -1.085,  -1.834/), &
                                                   (/   8.235,   6.271,   4.282,   1.101,   0.717/)  )
       Anomalous_ScFac(48)=Anomalous_Sc_Type("cd", (/  -1.303,  -0.416,  -0.079,  -1.005,  -1.637/), &
                                                   (/   8.912,   6.800,   4.653,   1.202,   0.783/)  )
       Anomalous_ScFac(49)=Anomalous_Sc_Type("in", (/  -1.788,  -0.626,  -0.126,  -0.936,  -1.493/), &
                                                   (/   9.627,   7.356,   5.045,   1.310,   0.854/)  )
       Anomalous_ScFac(50)=Anomalous_Sc_Type("sn", (/  -2.401,  -0.888,  -0.194,  -0.873,  -1.378/), &
                                                   (/  10.380,   7.943,   5.459,   1.424,   0.930/)  )
       Anomalous_ScFac(51)=Anomalous_Sc_Type("sb", (/  -3.194,  -1.214,  -0.287,  -0.816,  -1.284/), &
                                                   (/  11.166,   8.557,   5.894,   1.546,   1.010/)  )
       Anomalous_ScFac(52)=Anomalous_Sc_Type("te", (/  -4.267,  -1.630,  -0.418,  -0.772,  -1.212/), &
                                                   (/  11.995,   9.203,   6.352,   1.675,   1.096/)  )
       Anomalous_ScFac(53)=Anomalous_Sc_Type("i ", (/  -5.852,  -2.147,  -0.579,  -0.726,  -1.144/), &
                                                   (/  12.850,   9.885,   6.835,   1.812,   1.187/)  )
       Anomalous_ScFac(54)=Anomalous_Sc_Type("xe", (/  -8.133,  -2.812,  -0.783,  -0.684,  -1.084/), &
                                                   (/  11.933,  10.608,   7.348,   1.958,   1.284/)  )
       Anomalous_ScFac(55)=Anomalous_Sc_Type("cs", (/ -10.742,  -3.652,  -1.022,  -0.644,  -1.029/), &
                                                   (/  12.919,  11.382,   7.904,   2.119,   1.391/)  )
       Anomalous_ScFac(56)=Anomalous_Sc_Type("ba", (/ -11.460,  -4.832,  -1.334,  -0.613,  -0.983/), &
                                                   (/   9.981,  12.164,   8.460,   2.282,   1.500/)  )
       Anomalous_ScFac(57)=Anomalous_Sc_Type("la", (/ -12.135,  -6.683,  -1.716,  -0.588,  -0.942/), &
                                                   (/   3.565,  12.937,   9.036,   2.452,   1.615/)  )
       Anomalous_ScFac(58)=Anomalous_Sc_Type("ce", (/  -9.574,  -8.388,  -2.170,  -0.564,  -0.904/), &
                                                   (/   3.843,  11.953,   9.648,   2.632,   1.735/)  )
       Anomalous_ScFac(59)=Anomalous_Sc_Type("pr", (/  -7.817, -12.457,  -2.939,  -0.530,  -0.859/), &
                                                   (/   4.130,   6.285,  10.535,   2.845,   1.873/)  )
       Anomalous_ScFac(60)=Anomalous_Sc_Type("nd", (/  -7.486, -11.016,  -3.431,  -0.535,  -0.842/), &
                                                   (/   4.427,   9.874,  10.933,   3.018,   1.995/)  )
       Anomalous_ScFac(61)=Anomalous_Sc_Type("pm", (/  -6.891, -12.122,  -4.357,  -0.530,  -0.818/), &
                                                   (/   4.741,   3.627,  11.614,   3.225,   2.135/)  )
       Anomalous_ScFac(62)=Anomalous_Sc_Type("sm", (/  -6.429,  -9.616,  -5.696,  -0.533,  -0.798/), &
                                                   (/   5.073,   3.883,  12.320,   3.442,   2.281/)  )
       Anomalous_ScFac(63)=Anomalous_Sc_Type("eu", (/  -6.050,  -8.352,  -7.718,  -0.542,  -0.782/), &
                                                   (/   5.416,   4.149,  11.276,   3.669,   2.435/)  )
       Anomalous_ScFac(64)=Anomalous_Sc_Type("gd", (/  -5.779,  -7.565,  -9.242,  -0.564,  -0.774/), &
                                                   (/   5.773,   4.427,  11.946,   3.904,   2.595/)  )
       Anomalous_ScFac(65)=Anomalous_Sc_Type("tb", (/  -5.525,  -6.980,  -9.498,  -0.591,  -0.767/), &
                                                   (/   6.153,   4.721,   9.242,   4.151,   2.764/)  )
       Anomalous_ScFac(66)=Anomalous_Sc_Type("dy", (/  -5.250,  -6.492, -10.423,  -0.619,  -0.761/), &
                                                   (/   6.549,   5.026,   9.748,   4.410,   2.940/)  )
       Anomalous_ScFac(67)=Anomalous_Sc_Type("ho", (/  -5.040,  -6.112, -12.255,  -0.666,  -0.765/), &
                                                   (/   6.958,   5.343,   3.704,   4.678,   3.124/)  )
       Anomalous_ScFac(68)=Anomalous_Sc_Type("er", (/  -4.878,  -5.810,  -9.733,  -0.723,  -0.773/), &
                                                   (/   7.387,   5.675,   3.937,   4.958,   3.316/)  )
       Anomalous_ScFac(69)=Anomalous_Sc_Type("tm", (/  -4.753,  -5.565,  -8.488,  -0.795,  -0.790/), &
                                                   (/   7.833,   6.022,   4.181,   5.248,   3.515/)  )
       Anomalous_ScFac(70)=Anomalous_Sc_Type("yb", (/  -4.652,  -5.361,  -7.701,  -0.884,  -0.815/), &
                                                   (/   8.291,   6.378,   4.432,   5.548,   3.723/)  )
       Anomalous_ScFac(71)=Anomalous_Sc_Type("lu", (/  -4.580,  -5.190,  -7.133,  -0.988,  -0.847/), &
                                                   (/   8.759,   6.745,   4.693,   5.858,   3.937/)  )
       Anomalous_ScFac(72)=Anomalous_Sc_Type("hf", (/  -4.592,  -5.088,  -6.715,  -1.118,  -0.890/), &
                                                   (/   9.277,   7.148,   4.977,   6.185,   4.164/)  )
       Anomalous_ScFac(73)=Anomalous_Sc_Type("ta", (/  -4.540,  -4.948,  -6.351,  -1.258,  -0.937/), &
                                                   (/   9.811,   7.565,   5.271,   6.523,   4.399/)  )
       Anomalous_ScFac(74)=Anomalous_Sc_Type("w ", (/  -4.499,  -4.823,  -6.048,  -1.421,  -0.993/), &
                                                   (/  10.364,   7.996,   5.577,   6.872,   4.643/)  )
       Anomalous_ScFac(75)=Anomalous_Sc_Type("re", (/  -4.483,  -4.719,  -5.790,  -1.598,  -1.048/), &
                                                   (/  10.929,   8.439,   5.891,   7.232,   4.894/)  )
       Anomalous_ScFac(76)=Anomalous_Sc_Type("os", (/  -4.503,  -4.647,  -5.581,  -1.816,  -1.127/), &
                                                   (/  11.520,   8.903,   6.221,   7.605,   5.156/)  )
       Anomalous_ScFac(77)=Anomalous_Sc_Type("ir", (/  -4.527,  -4.578,  -5.391,  -2.066,  -1.216/), &
                                                   (/  12.140,   9.389,   6.566,   7.990,   5.427/)  )
       Anomalous_ScFac(78)=Anomalous_Sc_Type("pt", (/  -4.584,  -4.535,  -5.233,  -2.352,  -1.319/), &
                                                   (/  12.787,   9.895,   6.925,   8.388,   5.708/)  )
       Anomalous_ScFac(79)=Anomalous_Sc_Type("au", (/  -4.668,  -4.510,  -5.096,  -2.688,  -1.438/), &
                                                   (/  13.451,  10.418,   7.297,   8.798,   5.998/)  )
       Anomalous_ScFac(80)=Anomalous_Sc_Type("hg", (/  -4.803,  -4.523,  -4.990,  -3.084,  -1.576/), &
                                                   (/  14.143,  10.963,   7.686,   9.223,   6.299/)  )
       Anomalous_ScFac(81)=Anomalous_Sc_Type("tl", (/  -4.945,  -4.532,  -4.883,  -3.556,  -1.730/), &
                                                   (/  14.860,  11.528,   8.089,   9.659,   6.610/)  )
       Anomalous_ScFac(82)=Anomalous_Sc_Type("pb", (/  -5.161,  -4.596,  -4.818,  -4.133,  -1.910/), &
                                                   (/  15.595,  12.108,   8.505,  10.102,   6.930/)  )
       Anomalous_ScFac(83)=Anomalous_Sc_Type("bi", (/  -5.420,  -4.688,  -4.776,  -4.861,  -2.116/), &
                                                   (/  16.341,  12.700,   8.930,  10.559,   7.258/)  )
       Anomalous_ScFac(84)=Anomalous_Sc_Type("po", (/  -5.742,  -4.817,  -4.756,  -5.924,  -2.353/), &
                                                   (/  17.139,  13.331,   9.383,  11.042,   7.600/)  )
       Anomalous_ScFac(85)=Anomalous_Sc_Type("at", (/  -6.132,  -4.992,  -4.772,  -7.444,  -2.630/), &
                                                   (/  17.942,  13.969,   9.843,   9.961,   7.949/)  )
       Anomalous_ScFac(86)=Anomalous_Sc_Type("rn", (/  -6.545,  -5.173,  -4.787,  -8.862,  -2.932/), &
                                                   (/  18.775,  14.629,  10.317,  10.403,   8.307/)  )
       Anomalous_ScFac(87)=Anomalous_Sc_Type("fr", (/  -7.052,  -5.402,  -4.833,  -7.912,  -3.285/), &
                                                   (/  19.615,  15.299,  10.803,   7.754,   8.674/)  )
       Anomalous_ScFac(88)=Anomalous_Sc_Type("ra", (/  -7.614,  -5.659,  -4.898,  -7.620,  -3.702/), &
                                                   (/  20.461,  15.977,  11.296,   8.105,   9.047/)  )
       Anomalous_ScFac(89)=Anomalous_Sc_Type("ac", (/  -8.318,  -5.976,  -4.994,  -7.725,  -4.192/), &
                                                   (/  21.327,  16.668,  11.799,   8.472,   9.428/)  )
       Anomalous_ScFac(90)=Anomalous_Sc_Type("th", (/  -9.150,  -6.313,  -5.091,  -8.127,  -4.784/), &
                                                   (/  22.240,  17.397,  12.330,   8.870,   9.819/)  )
       Anomalous_ScFac(91)=Anomalous_Sc_Type("pa", (/ -10.382,  -6.695,  -5.216,  -8.960,  -5.555/), &
                                                   (/  23.161,  18.140,  12.868,   9.284,  10.227/)  )
       Anomalous_ScFac(92)=Anomalous_Sc_Type("u ", (/ -10.930,  -7.126,  -5.359, -10.673,  -6.735/), &
                                                   (/  23.121,  18.879,  13.409,   9.654,  10.637/)  )
       Anomalous_ScFac(93)=Anomalous_Sc_Type("np", (/ -12.152,  -7.624,  -5.529, -11.158,  -7.842/), &
                                                   (/  24.097,  19.642,  13.967,   4.148,   9.570/)  )
       Anomalous_ScFac(94)=Anomalous_Sc_Type("pu", (/ -12.280,  -8.187,  -5.712,  -9.725,  -8.473/), &
                                                   (/  23.658,  20.425,  14.536,   4.330,   6.999/)  )
       Anomalous_ScFac(95)=Anomalous_Sc_Type("am", (/ -12.771,  -8.872,  -5.930,  -8.926,  -7.701/), &
                                                   (/  24.607,  21.173,  15.087,   4.511,   7.296/)  )
       Anomalous_ScFac(96)=Anomalous_Sc_Type("cm", (/ -13.513,  -9.743,  -6.176,  -8.416,  -7.388/), &
                                                   (/  25.540,  21.896,  15.634,   4.697,   7.589/)  )
       Anomalous_ScFac(97)=Anomalous_Sc_Type("bk", (/ -14.827, -10.539,  -6.498,  -7.990,  -7.485/), &
                                                   (/  26.801,  21.942,  16.317,   4.908,   7.931/)  )
       Anomalous_ScFac(98)=Anomalous_Sc_Type("ze", (/   0.000,   0.000,   0.000,   0.000,   0.000/), &
                                                   (/   0.000,   0.000,   0.000,   0.000,   0.000/)  )
       return
    End Subroutine Set_Delta_Fp_Fpp

    !!----
    !!---- Subroutine Set_Magnetic_Form()
    !!----
    !!----    Magnetic form factors <j0> P.J. Brown, ILL prep. SP.88BR5016
    !!----    (March 1988)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Magnetic_Form()

       if (.not. allocated(magnetic_form)) allocate(magnetic_form(num_mag_form))
       if (.not. allocated(magnetic_j2))   allocate(magnetic_j2(num_mag_j2))
       if (.not. allocated(magnetic_j4))   allocate(magnetic_j4(num_mag_j4))
       if (.not. allocated(magnetic_j6))   allocate(magnetic_j6(num_mag_j6))

       Magnetic_Form(  1) = Magnetic_Form_Type("MSC0", &
                                              (/  0.251200, 90.029602,  0.329000, 39.402100,  0.423500, 14.322200, -0.004300/) )
       Magnetic_Form(  2) = Magnetic_Form_Type("MSC1", &
                                              (/  0.488900, 51.160301,  0.520300, 14.076400, -0.028600,  0.179200,  0.018500/) )
       Magnetic_Form(  3) = Magnetic_Form_Type("MSC2", &
                                              (/  0.504800, 31.403500,  0.518600, 10.989700, -0.024100,  1.183100,  0.000000/) )
       Magnetic_Form(  4) = Magnetic_Form_Type("MTI0", &
                                              (/  0.465700, 33.589802,  0.549000,  9.879100, -0.029100,  0.323200,  0.012300/) )
       Magnetic_Form(  5) = Magnetic_Form_Type("MTI1", &
                                              (/  0.509300, 36.703300,  0.503200, 10.371300, -0.026300,  0.310600,  0.011600/) )
       Magnetic_Form(  6) = Magnetic_Form_Type("MTI2", &
                                              (/  0.509100, 24.976299,  0.516200,  8.756900, -0.028100,  0.916000,  0.001500/) )
       Magnetic_Form(  7) = Magnetic_Form_Type("MTI3", &
                                              (/  0.357100, 22.841299,  0.668800,  8.930600, -0.035400,  0.483300,  0.009900/) )
       Magnetic_Form(  8) = Magnetic_Form_Type("MV0 ", &
                                              (/  0.408600, 28.810900,  0.607700,  8.543700, -0.029500,  0.276800,  0.012300/) )
       Magnetic_Form(  9) = Magnetic_Form_Type("MV1 ", &
                                              (/  0.444400, 32.647900,  0.568300,  9.097100, -0.228500,  0.021800,  0.215000/) )
       Magnetic_Form( 10) = Magnetic_Form_Type("MV2 ", &
                                              (/  0.408500, 23.852600,  0.609100,  8.245600, -0.167600,  0.041500,  0.149600/) )
       Magnetic_Form( 11) = Magnetic_Form_Type("MV3 ", &
                                              (/  0.359800, 19.336399,  0.663200,  7.617200, -0.306400,  0.029600,  0.283500/) )
       Magnetic_Form( 12) = Magnetic_Form_Type("MV4 ", &
                                              (/  0.310600, 16.816000,  0.719800,  7.048700, -0.052100,  0.302000,  0.022100/) )
       Magnetic_Form( 13) = Magnetic_Form_Type("MCR0", &
                                              (/  0.113500, 45.199001,  0.348100, 19.493099,  0.547700,  7.354200, -0.009200/) )
       Magnetic_Form( 14) = Magnetic_Form_Type("MCR1", &
                                              (/ -0.097700,  0.047000,  0.454400, 26.005400,  0.557900,  7.489200,  0.083100/) )
       Magnetic_Form( 15) = Magnetic_Form_Type("MCR2", &
                                              (/  1.202400, -0.005500,  0.415800, 20.547501,  0.603200,  6.956000, -1.221800/) )
       Magnetic_Form( 16) = Magnetic_Form_Type("MCR3", &
                                              (/ -0.309400,  0.027400,  0.368000, 17.035500,  0.655900,  6.523600,  0.285600/) )
       Magnetic_Form( 17) = Magnetic_Form_Type("MCR4", &
                                              (/ -0.232000,  0.043300,  0.310100, 14.951800,  0.718200,  6.172600,  0.204200/) )
       Magnetic_Form( 18) = Magnetic_Form_Type("MMN0", &
                                              (/  0.243800, 24.962900,  0.147200, 15.672800,  0.618900,  6.540300, -0.010500/) )
       Magnetic_Form( 19) = Magnetic_Form_Type("MMN1", &
                                              (/ -0.013800,  0.421300,  0.423100, 24.667999,  0.590500,  6.654500, -0.001000/) )
       Magnetic_Form( 20) = Magnetic_Form_Type("MMN2", &
                                              (/  0.422000, 17.684000,  0.594800,  6.005000,  0.004300, -0.609000, -0.021900/) )
       Magnetic_Form( 21) = Magnetic_Form_Type("MMN3", &
                                              (/  0.419800, 14.282900,  0.605400,  5.468900,  0.924100, -0.008800, -0.949800/) )
       Magnetic_Form( 22) = Magnetic_Form_Type("MMN4", &
                                              (/  0.376000, 12.566100,  0.660200,  5.132900, -0.037200,  0.563000,  0.001100/) )
       Magnetic_Form( 23) = Magnetic_Form_Type("MMN5", &
                                              (/  0.74050,  5.07409,    0.29237,  11.66547,  -1.78834,   0.00593,   1.75568 /) )
       Magnetic_Form( 24) = Magnetic_Form_Type("MFE0", &
                                              (/  0.070600, 35.008499,  0.358900, 15.358300,  0.581900,  5.560600, -0.011400/) )
       Magnetic_Form( 25) = Magnetic_Form_Type("MFE1", &
                                              (/  0.125100, 34.963299,  0.362900, 15.514400,  0.522300,  5.591400, -0.010500/) )
       Magnetic_Form( 26) = Magnetic_Form_Type("MFE2", &
                                              (/  0.026300, 34.959702,  0.366800, 15.943500,  0.618800,  5.593500, -0.011900/) )
       Magnetic_Form( 27) = Magnetic_Form_Type("MFE3", &
                                              (/  0.397200, 13.244200,  0.629500,  4.903400, -0.031400,  0.349600,  0.004400/) )
       Magnetic_Form( 28) = Magnetic_Form_Type("MFE4", &
                                              (/  0.378200, 11.380000,  0.655600,  4.592000, -0.034600,  0.483300,  0.000500/) )
       Magnetic_Form( 29) = Magnetic_Form_Type("MCO0", &
                                              (/  0.413900, 16.161600,  0.601300,  4.780500, -0.151800,  0.021000,  0.134500/) )
       Magnetic_Form( 30) = Magnetic_Form_Type("MCO1", &
                                              (/  0.099000, 33.125198,  0.364500, 15.176800,  0.547000,  5.008100, -0.010900/) )
       Magnetic_Form( 31) = Magnetic_Form_Type("MCO2", &
                                              (/  0.433200, 14.355300,  0.585700,  4.607700, -0.038200,  0.133800,  0.017900/) )
       Magnetic_Form( 32) = Magnetic_Form_Type("MCO3", &
                                              (/  0.390200, 12.507800,  0.632400,  4.457400, -0.150000,  0.034300,  0.127200/) )
       Magnetic_Form( 33) = Magnetic_Form_Type("MCO4", &
                                              (/  0.351500, 10.778500,  0.677800,  4.234300, -0.038900,  0.240900,  0.009800/) )
       Magnetic_Form( 34) = Magnetic_Form_Type("MNI0", &
                                              (/ -0.017200, 35.739201,  0.317400, 14.268900,  0.713600,  4.566100, -0.014300/) )
       Magnetic_Form( 35) = Magnetic_Form_Type("MNI1", &
                                              (/  0.070500, 35.856098,  0.398400, 13.804200,  0.542700,  4.396500, -0.011800/) )
       Magnetic_Form( 36) = Magnetic_Form_Type("MNI2", &
                                              (/  0.016300, 35.882599,  0.391600, 13.223300,  0.605200,  4.338800, -0.013300/) )
       Magnetic_Form( 37) = Magnetic_Form_Type("MNI3", &
                                              (/ -0.013400, 35.867699,  0.267800, 12.332600,  0.761400,  4.236900, -0.016200/) )
       Magnetic_Form( 38) = Magnetic_Form_Type("MNI4", &
                                              (/ -0.009000, 35.861401,  0.277600, 11.790400,  0.747400,  4.201100, -0.016300/) )
       Magnetic_Form( 39) = Magnetic_Form_Type("MCU0", &
                                              (/  0.090900, 34.983799,  0.408800, 11.443200,  0.512800,  3.824800, -0.012400/) )
       Magnetic_Form( 40) = Magnetic_Form_Type("MCU1", &
                                              (/  0.074900, 34.965599,  0.414700, 11.764200,  0.523800,  3.849700, -0.012700/) )
       Magnetic_Form( 41) = Magnetic_Form_Type("MCU2", &
                                              (/  0.023200, 34.968601,  0.402300, 11.564000,  0.588200,  3.842800, -0.013700/) )
       Magnetic_Form( 42) = Magnetic_Form_Type("MCU3", &
                                              (/  0.003100, 34.907398,  0.358200, 10.913800,  0.653100,  3.827900, -0.014700/) )
       Magnetic_Form( 43) = Magnetic_Form_Type("MCU4", &
                                              (/ -0.013200, 30.681700,  0.280100, 11.162600,  0.749000,  3.817200, -0.016500/) )
       Magnetic_Form( 44) = Magnetic_Form_Type("MY0 ", &
                                              (/  0.591500, 67.608101,  1.512300, 17.900400, -1.113000, 14.135900,  0.008000/) )
       Magnetic_Form( 45) = Magnetic_Form_Type("MZR0", &
                                              (/  0.410600, 59.996101,  1.054300, 18.647600, -0.475100, 10.540000,  0.010600/) )
       Magnetic_Form( 46) = Magnetic_Form_Type("MZR1", &
                                              (/  0.453200, 59.594799,  0.783400, 21.435699, -0.245100,  9.036000,  0.009800/) )
       Magnetic_Form( 47) = Magnetic_Form_Type("MNB0", &
                                              (/  0.394600, 49.229698,  1.319700, 14.821600, -0.726900,  9.615600,  0.012900/) )
       Magnetic_Form( 48) = Magnetic_Form_Type("MNB1", &
                                              (/  0.457200, 49.918201,  1.027400, 15.725600, -0.496200,  9.157300,  0.011800/) )
       Magnetic_Form( 49) = Magnetic_Form_Type("MMO0", &
                                              (/  0.180600, 49.056801,  1.230600, 14.785900, -0.426800,  6.986600,  0.017100/) )
       Magnetic_Form( 50) = Magnetic_Form_Type("MMO1", &
                                              (/  0.350000, 48.035400,  1.030500, 15.060400, -0.392900,  7.479000,  0.013900/) )
       Magnetic_Form( 51) = Magnetic_Form_Type("MTC0", &
                                              (/  0.129800, 49.661098,  1.165600, 14.130700, -0.313400,  5.512900,  0.019500/) )
       Magnetic_Form( 52) = Magnetic_Form_Type("MTC1", &
                                              (/  0.267400, 48.956600,  0.956900, 15.141300, -0.238700,  5.457800,  0.016000/) )
       Magnetic_Form( 53) = Magnetic_Form_Type("MRU0", &
                                              (/  0.106900, 49.423801,  1.191200, 12.741700, -0.317600,  4.912500,  0.021300/) )
       Magnetic_Form( 54) = Magnetic_Form_Type("MRU1", &
                                              (/  0.441000, 33.308601,  1.477500,  9.553100, -0.936100,  6.722000,  0.017600/) )
       Magnetic_Form( 55) = Magnetic_Form_Type("MRH0", &
                                              (/  0.097600, 49.882500,  1.160100, 11.830700, -0.278900,  4.126600,  0.023400/) )
       Magnetic_Form( 56) = Magnetic_Form_Type("MRH1", &
                                              (/  0.334200, 29.756399,  1.220900,  9.438400, -0.575500,  5.332000,  0.021000/) )
       Magnetic_Form( 57) = Magnetic_Form_Type("MPD0", &
                                              (/  0.200300, 29.363300,  1.144600,  9.599300, -0.368900,  4.042300,  0.025100/) )
       Magnetic_Form( 58) = Magnetic_Form_Type("MPD1", &
                                              (/  0.503300, 24.503700,  1.998200,  6.908200, -1.524000,  5.513300,  0.021300/) )
       Magnetic_Form( 59) = Magnetic_Form_Type("MCE2", &
                                              (/  0.295300, 17.684601,  0.292300,  6.732900,  0.431300,  5.382700, -0.019400/) )
       Magnetic_Form( 60) = Magnetic_Form_Type("MND2", &
                                              (/  0.164500, 25.045300,  0.252200, 11.978200,  0.601200,  4.946100, -0.018000/) )
       Magnetic_Form( 61) = Magnetic_Form_Type("MND3", &
                                              (/  0.054000, 25.029301,  0.310100, 12.102000,  0.657500,  4.722300, -0.021600/) )
       Magnetic_Form( 62) = Magnetic_Form_Type("MSM2", &
                                              (/  0.090900, 25.203199,  0.303700, 11.856200,  0.625000,  4.236600, -0.020000/) )
       Magnetic_Form( 63) = Magnetic_Form_Type("MSM3", &
                                              (/  0.028800, 25.206800,  0.297300, 11.831100,  0.695400,  4.211700, -0.021300/) )
       Magnetic_Form( 64) = Magnetic_Form_Type("MEU2", &
                                              (/  0.075500, 25.296000,  0.300100, 11.599300,  0.643800,  4.025200, -0.019600/) )
       Magnetic_Form( 65) = Magnetic_Form_Type("MEU3", &
                                              (/  0.020400, 25.307800,  0.301000, 11.474400,  0.700500,  3.942000, -0.022000/) )
       Magnetic_Form( 66) = Magnetic_Form_Type("MGD2", &
                                              (/  0.063600, 25.382299,  0.303300, 11.212500,  0.652800,  3.787700, -0.019900/) )
       Magnetic_Form( 67) = Magnetic_Form_Type("MGD3", &
                                              (/  0.018600, 25.386700,  0.289500, 11.142100,  0.713500,  3.752000, -0.021700/) )
       Magnetic_Form( 68) = Magnetic_Form_Type("MTB2", &
                                              (/  0.054700, 25.508600,  0.317100, 10.591100,  0.649000,  3.517100, -0.021200/) )
       Magnetic_Form( 69) = Magnetic_Form_Type("MTB3", &
                                              (/  0.017700, 25.509501,  0.292100, 10.576900,  0.713300,  3.512200, -0.023100/) )
       Magnetic_Form( 70) = Magnetic_Form_Type("MDY2", &
                                              (/  0.130800, 18.315500,  0.311800,  7.664500,  0.579500,  3.146900, -0.022600/) )
       Magnetic_Form( 71) = Magnetic_Form_Type("MDY3", &
                                              (/  0.115700, 15.073200,  0.327000,  6.799100,  0.582100,  3.020200, -0.024900/) )
       Magnetic_Form( 72) = Magnetic_Form_Type("MHO2", &
                                              (/  0.099500, 18.176100,  0.330500,  7.855600,  0.592100,  2.979900, -0.023000/) )
       Magnetic_Form( 73) = Magnetic_Form_Type("MHO3", &
                                              (/  0.056600, 18.317600,  0.336500,  7.688000,  0.631700,  2.942700, -0.024800/) )
       Magnetic_Form( 74) = Magnetic_Form_Type("MER2", &
                                              (/  0.112200, 18.122299,  0.346200,  6.910600,  0.564900,  2.761400, -0.023500/) )
       Magnetic_Form( 75) = Magnetic_Form_Type("MER3", &
                                              (/  0.058600, 17.980200,  0.354000,  7.096400,  0.612600,  2.748200, -0.025100/) )
       Magnetic_Form( 76) = Magnetic_Form_Type("MTM2", &
                                              (/  0.098300, 18.323601,  0.338000,  6.917800,  0.587500,  2.662200, -0.024100/) )
       Magnetic_Form( 77) = Magnetic_Form_Type("MTM3", &
                                              (/  0.058100, 15.092200,  0.278700,  7.801500,  0.685400,  2.793100, -0.022400/) )
       Magnetic_Form( 78) = Magnetic_Form_Type("MYB2", &
                                              (/  0.085500, 18.512300,  0.294300,  7.373400,  0.641200,  2.677700, -0.021300/) )
       Magnetic_Form( 79) = Magnetic_Form_Type("MYB3", &
                                              (/  0.041600, 16.094900,  0.284900,  7.834100,  0.696100,  2.672500, -0.022900/) )
       Magnetic_Form( 80) = Magnetic_Form_Type("MU3 ", &
                                              (/  0.505800, 23.288200,  1.346400,  7.002800, -0.872400,  4.868300,  0.019200/) )
       Magnetic_Form( 81) = Magnetic_Form_Type("MU4 ", &
                                              (/  0.329100, 23.547501,  1.083600,  8.454000, -0.434000,  4.119600,  0.021400/) )
       Magnetic_Form( 82) = Magnetic_Form_Type("MU5 ", &
                                              (/  0.365000, 19.803801,  3.219900,  6.281800, -2.607700,  5.301000,  0.023300/) )
       Magnetic_Form( 83) = Magnetic_Form_Type("MNP3", &
                                              (/  0.515700, 20.865400,  2.278400,  5.893000, -1.816300,  4.845700,  0.021100/) )
       Magnetic_Form( 84) = Magnetic_Form_Type("MNP4", &
                                              (/  0.420600, 19.804600,  2.800400,  5.978300, -2.243600,  4.984800,  0.022800/) )
       Magnetic_Form( 85) = Magnetic_Form_Type("MNP5", &
                                              (/  0.369200, 18.190001,  3.151000,  5.850000, -2.544600,  4.916400,  0.024800/) )
       Magnetic_Form( 86) = Magnetic_Form_Type("MNP6", &
                                              (/  0.292900, 17.561100,  3.486600,  5.784700, -2.806600,  4.870700,  0.026700/) )
       Magnetic_Form( 87) = Magnetic_Form_Type("MPU3", &
                                              (/  0.384000, 16.679300,  3.104900,  5.421000, -2.514800,  4.551200,  0.026300/) )
       Magnetic_Form( 88) = Magnetic_Form_Type("MPU4", &
                                              (/  0.493400, 16.835501,  1.639400,  5.638400, -1.158100,  4.139900,  0.024800/) )
       Magnetic_Form( 89) = Magnetic_Form_Type("MPU5", &
                                              (/  0.388800, 16.559200,  2.036200,  5.656700, -1.451500,  4.255200,  0.026700/) )
       Magnetic_Form( 90) = Magnetic_Form_Type("MPU6", &
                                              (/  0.317200, 16.050699,  3.465400,  5.350700, -2.810200,  4.513300,  0.028100/) )
       Magnetic_Form( 91) = Magnetic_Form_Type("MAM2", &
                                              (/  0.474300, 21.776100,  1.580000,  5.690200, -1.077900,  4.145100,  0.021800/) )
       Magnetic_Form( 92) = Magnetic_Form_Type("MAM3", &
                                              (/  0.423900, 19.573900,  1.457300,  5.872200, -0.905200,  3.968200,  0.023800/) )
       Magnetic_Form( 93) = Magnetic_Form_Type("MAM4", &
                                              (/  0.373700, 17.862499,  1.352100,  6.042600, -0.751400,  3.719900,  0.025800/) )
       Magnetic_Form( 94) = Magnetic_Form_Type("MAM5", &
                                              (/  0.295600, 17.372499,  1.452500,  6.073400, -0.775500,  3.661900,  0.027700/) )
       Magnetic_Form( 95) = Magnetic_Form_Type("MAM6", &
                                              (/  0.230200, 16.953300,  1.486400,  6.115900, -0.745700,  3.542600,  0.029400/) )
       Magnetic_Form( 96) = Magnetic_Form_Type("MAM7", &
                                              (/  0.360100, 12.729900,  1.964000,  5.120300, -1.356000,  3.714200,  0.031600/) )
       Magnetic_Form( 97) = Magnetic_Form_Type("MPR3", &
                                              (/  0.050400, 24.998900,  0.257200, 12.037700,  0.714200,  5.003900, -0.021900/) )
       Magnetic_Form( 98) = Magnetic_Form_Type("MO1", &
                                              (/  0.115285, 85.197300,  0.556229, 25.252200,  0.332476,  6.362070, -0.00460676/) )
       Magnetic_Form( 99) = Magnetic_Form_Type("MXX1", &
                                              (/  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/) ) !for future use
       Magnetic_Form(100) = Magnetic_Form_Type("MXX2", &
                                              (/  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/) ) !for future use
       Magnetic_Form(101) = Magnetic_Form_Type("JCE2", &
                                              (/  0.031972,  8.926222,  0.265792,  7.678510,  0.682151,  2.329783,  0.020578/) )
       Magnetic_Form(102) = Magnetic_Form_Type("JCE3", &
                                              (/  0.051183,  6.115375,  0.277738,  7.952485,  0.654079,  2.287000,  0.016355/) )
       Magnetic_Form(103) = Magnetic_Form_Type("JPR3", &
                                              (/  0.023288,  0.582954,  0.349391,  5.601756,  0.615363,  1.932779,  0.011454/) )
       Magnetic_Form(104) = Magnetic_Form_Type("JND2", &
                                              (/  0.089354,  2.282004,  0.206157,  1.708607,  0.669916,  2.297662,  0.048390/) )
       Magnetic_Form(105) = Magnetic_Form_Type("JND3", &
                                              (/  0.073287,  4.412361,  0.371485,  4.019648,  0.539459,  1.557985,  0.017335/) )
       Magnetic_Form(106) = Magnetic_Form_Type("JGD3", &
                                              (/  0.060537, 10.775218,  0.271475, 13.097898,  0.665241,  3.162837,  0.001566/) )
       Magnetic_Form(107) = Magnetic_Form_Type("JTB2", &
                                              (/  0.049801, 18.734161,  0.277437, 10.084129,  0.661194,  2.745624,  0.010774/) )
       Magnetic_Form(108) = Magnetic_Form_Type("JTB3", &
                                              (/  0.049792, 15.112189,  0.270644,  9.158312,  0.679388,  2.880260, -0.000131/) )
       Magnetic_Form(109) = Magnetic_Form_Type("JDY2", &
                                              (/  0.175586,  5.938148,  0.228867, 11.464046,  0.583298,  2.167554,  0.011186/) )
       Magnetic_Form(110) = Magnetic_Form_Type("JDY3", &
                                              (/  0.146536, 12.639305,  0.375822,  5.511785,  0.515731,  2.090789,  0.093576/) )
       Magnetic_Form(111) = Magnetic_Form_Type("JHO2", &
                                              (/  0.023234,  0.703240,  0.270745,  9.993475,  0.677581,  2.521403,  0.027101/) )
       Magnetic_Form(112) = Magnetic_Form_Type("JHO2", &
                                              (/  0.023234,  0.703240,  0.270745,  9.993475,  0.677581,  2.521403,  0.027101/) )
       Magnetic_Form(113) = Magnetic_Form_Type("JHO3", &
                                              (/  0.043204,  0.910121,  0.279392,  8.683387,  0.668537,  2.417518,  0.008207/) )
       Magnetic_Form(114) = Magnetic_Form_Type("JER2", &
                                              (/  0.037734,  6.081446,  0.256447,  9.598846,  0.679204,  2.139296,  0.025543/) )
       Magnetic_Form(115) = Magnetic_Form_Type("JER3", &
                                              (/  0.038871,  5.311772,  0.259781,  8.173226,  0.678414,  2.082836,  0.022169/) )
       Magnetic_Form(116) = Magnetic_Form_Type("JTM2", &
                                              (/  0.037670,  4.455198,  0.254184,  9.151058,  0.677308,  2.021746,  0.029718/) )
       Magnetic_Form(117) = Magnetic_Form_Type("JTM3", &
                                              (/  0.028279,  2.291633,  0.265583,  7.776700,  0.675720,  2.018924,  0.029883/) )
       Magnetic_Form(118) = Magnetic_Form_Type("JYB3", &
                                           (/  0.092380,  2.046342,  0.258408,  7.471918,  0.609716,  1.913869,  0.038824/) )
       Magnetic_Form(119) = Magnetic_Form_Type("JO1 ", &
                                              (/  0.115285, 85.197300,  0.556229, 25.252200,  0.332476,  6.362070,-0.00460676/) )

       !---- <j2> Coefficients ----!
       Magnetic_j2(  1) = Magnetic_Form_Type("SC0 ",(/10.8172,54.327, 4.7353,14.847, 0.6071, 4.218, 0.0011/))
       Magnetic_j2(  2) = Magnetic_Form_Type("SC1 ",(/ 8.5021,34.285, 3.2116,10.994, 0.4244, 3.605, 0.0009/))
       Magnetic_j2(  3) = Magnetic_Form_Type("SC2 ",(/ 4.3683,28.654, 3.7231,10.823, 0.6074, 3.668, 0.0014/))
       Magnetic_j2(  4) = Magnetic_Form_Type("TI0 ",(/ 4.3583,36.056, 3.8230,11.133, 0.6855, 3.469, 0.0020/))
       Magnetic_j2(  5) = Magnetic_Form_Type("TI1 ",(/ 6.1567,27.275, 2.6833, 8.983, 0.4070, 3.052, 0.0011/))
       Magnetic_j2(  6) = Magnetic_Form_Type("TI2 ",(/ 4.3107,18.348, 2.0960, 6.797, 0.2984, 2.548, 0.0007/))
       Magnetic_j2(  7) = Magnetic_Form_Type("TI3 ",(/ 3.3717,14.444, 1.8258, 5.713, 0.2470, 2.265, 0.0005/))
       Magnetic_j2(  8) = Magnetic_Form_Type("V0  ",(/ 3.8099,21.347, 2.3295, 7.409, 0.4333, 2.632, 0.0015/))
       Magnetic_j2(  9) = Magnetic_Form_Type("V1  ",(/ 4.7474,23.323, 2.3609, 7.808, 0.4105, 2.706, 0.0014/))
       Magnetic_j2( 10) = Magnetic_Form_Type("V2  ",(/ 3.4386,16.530, 1.9638, 6.141, 0.2997, 2.267, 0.0009/))
       Magnetic_j2( 11) = Magnetic_Form_Type("V3  ",(/ 2.3005,14.682, 2.0364, 6.130, 0.4099, 2.382, 0.0014/))
       Magnetic_j2( 12) = Magnetic_Form_Type("V4  ",(/ 1.8377,12.267, 1.8247, 5.458, 0.3979, 2.248, 0.0012/))
       Magnetic_j2( 13) = Magnetic_Form_Type("CR0 ",(/ 3.4085,20.127, 2.1006, 6.802, 0.4266, 2.394, 0.0019/))
       Magnetic_j2( 14) = Magnetic_Form_Type("CR1 ",(/ 3.7768,20.346, 2.1028, 6.893, 0.4010, 2.411, 0.0017/))
       Magnetic_j2( 15) = Magnetic_Form_Type("CR2 ",(/ 2.6422,16.060, 1.9198, 6.253, 0.4446, 2.372, 0.0020/))
       Magnetic_j2( 16) = Magnetic_Form_Type("CR3 ",(/ 1.6262,15.066, 2.0618, 6.284, 0.5281, 2.368, 0.0023/))
       Magnetic_j2( 17) = Magnetic_Form_Type("CR4 ",(/ 1.0293,13.950, 1.9933, 6.059, 0.5974, 2.346, 0.0027/))
       Magnetic_j2( 18) = Magnetic_Form_Type("MN0 ",(/ 2.6681,16.060, 1.7561, 5.640, 0.3675, 2.049, 0.0017/))
       Magnetic_j2( 19) = Magnetic_Form_Type("MN1 ",(/ 3.2953,18.695, 1.8792, 6.240, 0.3927, 2.201, 0.0022/))
       Magnetic_j2( 20) = Magnetic_Form_Type("MN2 ",(/ 2.0515,15.556, 1.8841, 6.063, 0.4787, 2.232, 0.0027/))
       Magnetic_j2( 21) = Magnetic_Form_Type("MN3 ",(/ 1.2427,14.997, 1.9567, 6.118, 0.5732, 2.258, 0.0031/))
       Magnetic_j2( 22) = Magnetic_Form_Type("MN4 ",(/ 0.7879,13.886, 1.8717, 5.743, 0.5981, 2.182, 0.0034/))
       Magnetic_j2( 23) = Magnetic_Form_Type("MN5 ",(/-0.11904,6.59893,-0.23941,10.73086, 0.35048,1.49116,0.00776/))
       Magnetic_j2( 24) = Magnetic_Form_Type("FE0 ",(/ 1.9405,18.473, 1.9566, 6.323, 0.5166, 2.161, 0.0036/))
       Magnetic_j2( 25) = Magnetic_Form_Type("FE1 ",(/ 2.6290,18.660, 1.8704, 6.331, 0.4690, 2.163, 0.0031/))
       Magnetic_j2( 26) = Magnetic_Form_Type("FE2 ",(/ 1.6490,16.559, 1.9064, 6.133, 0.5206, 2.137, 0.0035/))
       Magnetic_j2( 27) = Magnetic_Form_Type("FE3 ",(/ 1.3602,11.998, 1.5188, 5.003, 0.4705, 1.991, 0.0038/))
       Magnetic_j2( 28) = Magnetic_Form_Type("FE4 ",(/ 1.5582, 8.275, 1.1863, 3.279, 0.1366, 1.107,-0.0022/))
       Magnetic_j2( 29) = Magnetic_Form_Type("CO0 ",(/ 1.9678,14.170, 1.4911, 4.948, 0.3844, 1.797, 0.0027/))
       Magnetic_j2( 30) = Magnetic_Form_Type("CO1 ",(/ 2.4097,16.161, 1.5780, 5.460, 0.4095, 1.914, 0.0031/))
       Magnetic_j2( 31) = Magnetic_Form_Type("CO2 ",(/ 1.9049,11.644, 1.3159, 4.357, 0.3146, 1.645, 0.0017/))
       Magnetic_j2( 32) = Magnetic_Form_Type("CO3 ",(/ 1.7058, 8.859, 1.1409, 3.309, 0.1474, 1.090,-0.0025/))
       Magnetic_j2( 33) = Magnetic_Form_Type("CO4 ",(/ 1.3110, 8.025, 1.1551, 3.179, 0.1608, 1.130,-0.0011/))
       Magnetic_j2( 34) = Magnetic_Form_Type("NI0 ",(/ 1.0302,12.252, 1.4669, 4.745, 0.4521, 1.744, 0.0036/))
       Magnetic_j2( 35) = Magnetic_Form_Type("NI1 ",(/ 2.1040,14.866, 1.4302, 5.071, 0.4031, 1.778, 0.0034/))
       Magnetic_j2( 36) = Magnetic_Form_Type("NI2 ",(/ 1.7080,11.016, 1.2147, 4.103, 0.3150, 1.533, 0.0018/))
       Magnetic_j2( 37) = Magnetic_Form_Type("NI3 ",(/ 1.1612, 7.700, 1.0027, 3.263, 0.2719, 1.378, 0.0025/))
       Magnetic_j2( 38) = Magnetic_Form_Type("NI4 ",(/ 1.1612, 7.700, 1.0027, 3.263, 0.2719, 1.378, 0.0025/))
       Magnetic_j2( 39) = Magnetic_Form_Type("CU0 ",(/ 1.9182,14.490, 1.3329, 4.730, 0.3842, 1.639, 0.0035/))
       Magnetic_j2( 40) = Magnetic_Form_Type("CU1 ",(/ 1.8814,13.433, 1.2809, 4.545, 0.3646, 1.602, 0.0033/))
       Magnetic_j2( 41) = Magnetic_Form_Type("CU2 ",(/ 1.5189,10.478, 1.1512, 3.813, 0.2918, 1.398, 0.0017/))
       Magnetic_j2( 42) = Magnetic_Form_Type("CU3 ",(/ 1.2797, 8.450, 1.0315, 3.280, 0.2401, 1.250, 0.0015/))
       Magnetic_j2( 43) = Magnetic_Form_Type("CU4 ",(/ 0.9568, 7.448, 0.9099, 3.396, 0.3729, 1.494, 0.0049/))
       Magnetic_j2( 44) = Magnetic_Form_Type("Y0  ",(/14.4084,44.658, 5.1045,14.904,-0.0535, 3.319, 0.0028/))
       Magnetic_j2( 45) = Magnetic_Form_Type("ZR0 ",(/10.1378,35.337, 4.7734,12.545,-0.0489, 2.672, 0.0036/))
       Magnetic_j2( 46) = Magnetic_Form_Type("ZR1 ",(/11.8722,34.920, 4.0502,12.127,-0.0632, 2.828, 0.0034/))
       Magnetic_j2( 47) = Magnetic_Form_Type("NB0 ",(/ 7.4796,33.179, 5.0884,11.571,-0.0281, 1.564, 0.0047/))
       Magnetic_j2( 48) = Magnetic_Form_Type("NB1 ",(/ 8.7735,33.285, 4.6556,11.605,-0.0268, 1.539, 0.0044/))
       Magnetic_j2( 49) = Magnetic_Form_Type("MO0 ",(/ 5.1180,23.422, 4.1809, 9.208,-0.0505, 1.743, 0.0053/))
       Magnetic_j2( 50) = Magnetic_Form_Type("MO1 ",(/ 7.2367,28.128, 4.0705, 9.923,-0.0317, 1.455, 0.0049/))
       Magnetic_j2( 51) = Magnetic_Form_Type("TC0 ",(/ 4.2441,21.397, 3.9439, 8.375,-0.0371, 1.187, 0.0066/))
       Magnetic_j2( 52) = Magnetic_Form_Type("TC1 ",(/ 6.4056,24.824, 3.5400, 8.611,-0.0366, 1.485, 0.0044/))
       Magnetic_j2( 53) = Magnetic_Form_Type("RU0 ",(/ 3.7445,18.613, 3.4749, 7.420,-0.0363, 1.007, 0.0073/))
       Magnetic_j2( 54) = Magnetic_Form_Type("RU1 ",(/ 5.2826,23.683, 3.5813, 8.152,-0.0257, 0.426, 0.0131/))
       Magnetic_j2( 55) = Magnetic_Form_Type("RH0 ",(/ 3.3651,17.344, 3.2121, 6.804,-0.0350, 0.503, 0.0146/))
       Magnetic_j2( 56) = Magnetic_Form_Type("RH1 ",(/ 4.0260,18.950, 3.1663, 7.000,-0.0296, 0.486, 0.0127/))
       Magnetic_j2( 57) = Magnetic_Form_Type("PD0 ",(/ 3.3105,14.726, 2.6332, 5.862,-0.0437, 1.130, 0.0053/))
       Magnetic_j2( 58) = Magnetic_Form_Type("PD1 ",(/ 4.2749,17.900, 2.7021, 6.354,-0.0258, 0.700, 0.0071/))
       Magnetic_j2( 59) = Magnetic_Form_Type("CE2 ",(/ 0.9809,18.063, 1.8413, 7.769, 0.9905, 2.845, 0.0120/))
       Magnetic_j2( 60) = Magnetic_Form_Type("ND2 ",(/ 1.4530,18.340, 1.6196, 7.285, 0.8752, 2.622, 0.0126/))
       Magnetic_j2( 61) = Magnetic_Form_Type("ND3 ",(/ 0.6751,18.342, 1.6272, 7.260, 0.9644, 2.602, 0.0150/))
       Magnetic_j2( 62) = Magnetic_Form_Type("SM2 ",(/ 1.0360,18.425, 1.4769, 7.032, 0.8810, 2.437, 0.0152/))
       Magnetic_j2( 63) = Magnetic_Form_Type("SM3 ",(/ 0.4707,18.430, 1.4261, 7.034, 0.9574, 2.439, 0.0182/))
       Magnetic_j2( 64) = Magnetic_Form_Type("EU2 ",(/ 0.8970,18.443, 1.3769, 7.005, 0.9060, 2.421, 0.0190/))
       Magnetic_j2( 65) = Magnetic_Form_Type("EU3 ",(/ 0.3985,18.451, 1.3307, 6.956, 0.9603, 2.378, 0.0197/))
       Magnetic_j2( 66) = Magnetic_Form_Type("GD2 ",(/ 0.7756,18.469, 1.3124, 6.899, 0.8956, 2.338, 0.0199/))
       Magnetic_j2( 67) = Magnetic_Form_Type("GD3 ",(/ 0.3347,18.476, 1.2465, 6.877, 0.9537, 2.318, 0.0217/))
       Magnetic_j2( 68) = Magnetic_Form_Type("TB2 ",(/ 0.6688,18.491, 1.2487, 6.822, 0.8888, 2.275, 0.0215/))
       Magnetic_j2( 69) = Magnetic_Form_Type("TB3 ",(/ 0.2892,18.497, 1.1678, 6.797, 0.9437, 2.257, 0.0232/))
       Magnetic_j2( 70) = Magnetic_Form_Type("DY2 ",(/ 0.5917,18.511, 1.1828, 6.747, 0.8801, 2.214, 0.0229/))
       Magnetic_j2( 71) = Magnetic_Form_Type("DY3 ",(/ 0.2523,18.517, 1.0914, 6.736, 0.9345, 2.208, 0.0250/))
       Magnetic_j2( 72) = Magnetic_Form_Type("HO2 ",(/ 0.5094,18.515, 1.1234, 6.706, 0.8727, 2.159, 0.0242/))
       Magnetic_j2( 73) = Magnetic_Form_Type("HO3 ",(/ 0.2188,18.516, 1.0240, 6.707, 0.9251, 2.161, 0.0268/))
       Magnetic_j2( 74) = Magnetic_Form_Type("ER2 ",(/ 0.4693,18.528, 1.0545, 6.649, 0.8679, 2.120, 0.0261/))
       Magnetic_j2( 75) = Magnetic_Form_Type("ER3 ",(/ 0.1710,18.534, 0.9879, 6.625, 0.9044, 2.100, 0.0278/))
       Magnetic_j2( 76) = Magnetic_Form_Type("TM2 ",(/ 0.4198,18.542, 0.9959, 6.600, 0.8593, 2.082, 0.0284/))
       Magnetic_j2( 77) = Magnetic_Form_Type("TM3 ",(/ 0.1760,18.542, 0.9105, 6.579, 0.8970, 2.062, 0.0294/))
       Magnetic_j2( 78) = Magnetic_Form_Type("YB2 ",(/ 0.3852,18.550, 0.9415, 6.551, 0.8492, 2.043, 0.0301/))
       Magnetic_j2( 79) = Magnetic_Form_Type("YB3 ",(/ 0.1570,18.555, 0.8484, 6.540, 0.8880, 2.037, 0.0318/))
       Magnetic_j2( 80) = Magnetic_Form_Type("U3  ",(/ 4.1582,16.534, 2.4675, 5.952,-0.0252, 0.765, 0.0057/))
       Magnetic_j2( 81) = Magnetic_Form_Type("U4  ",(/ 3.7449,13.894, 2.6453, 4.863,-0.5218, 3.192, 0.0009/))
       Magnetic_j2( 82) = Magnetic_Form_Type("U5  ",(/ 3.0724,12.546, 2.3076, 5.231,-0.0644, 1.474, 0.0035/))
       Magnetic_j2( 83) = Magnetic_Form_Type("NP3 ",(/ 3.7170,15.133, 2.3216, 5.503,-0.0275, 0.800, 0.0052/))
       Magnetic_j2( 84) = Magnetic_Form_Type("NP4 ",(/ 2.9203,14.646, 2.5979, 5.559,-0.0301, 0.367, 0.0141/))
       Magnetic_j2( 85) = Magnetic_Form_Type("NP5 ",(/ 2.3308,13.654, 2.7219, 5.494,-0.1357, 0.049, 0.1224/))
       Magnetic_j2( 86) = Magnetic_Form_Type("NP6 ",(/ 1.8245,13.180, 2.8508, 5.407,-0.1579, 0.044, 0.1438/))
       Magnetic_j2( 87) = Magnetic_Form_Type("PU3 ",(/ 2.0885,12.871, 2.5961, 5.190,-0.1465, 0.039, 0.1343/))
       Magnetic_j2( 88) = Magnetic_Form_Type("PU4 ",(/ 2.7244,12.926, 2.3387, 5.163,-0.1300, 0.046, 0.1177/))
       Magnetic_j2( 89) = Magnetic_Form_Type("PU5 ",(/ 2.1409,12.832, 2.5664, 5.152,-0.1338, 0.046, 0.1210/))
       Magnetic_j2( 90) = Magnetic_Form_Type("PU6 ",(/ 1.7262,12.324, 2.6652, 5.066,-0.1695, 0.041, 0.1550/))
       Magnetic_j2( 91) = Magnetic_Form_Type("AM2 ",(/ 3.5237,15.955, 2.2855, 5.195,-0.0142, 0.585, 0.0033/))
       Magnetic_j2( 92) = Magnetic_Form_Type("AM3 ",(/ 2.8622,14.733, 2.4099, 5.144,-0.1326, 0.031, 0.1233/))
       Magnetic_j2( 93) = Magnetic_Form_Type("AM4 ",(/ 2.4141,12.948, 2.3687, 4.945,-0.2490, 0.022, 0.2371/))
       Magnetic_j2( 94) = Magnetic_Form_Type("AM5 ",(/ 2.0109,12.053, 2.4155, 4.836,-0.2264, 0.027, 0.2128/))
       Magnetic_j2( 95) = Magnetic_Form_Type("AM6 ",(/ 1.6778,11.337, 2.4531, 4.725,-0.2043, 0.034, 0.1892/))
       Magnetic_j2( 96) = Magnetic_Form_Type("AM7 ",(/ 1.8845, 9.161, 2.0746, 4.042,-0.1318, 1.723, 0.0020/))
       Magnetic_j2( 97) = Magnetic_Form_Type("PR3 ",(/ 0.8734,18.9876,1.5594,6.0872, 0.8142,2.4150, 0.0111/))

       !---- <j4> Coefficients ----!
       Magnetic_j4(  1) = Magnetic_Form_Type("SC0 ",(/ 1.3420,10.200, 0.3837, 3.079, 0.0468, 0.118,-0.0328/))
       Magnetic_j4(  2) = Magnetic_Form_Type("SC1 ",(/ 7.1167,15.487,-6.6671,18.269, 0.4900, 2.992, 0.0047/))
       Magnetic_j4(  3) = Magnetic_Form_Type("SC2 ",(/-1.6684,15.648, 1.7742, 9.062, 0.4075, 2.412, 0.0042/))
       Magnetic_j4(  4) = Magnetic_Form_Type("TI0 ",(/-2.1515,11.271, 2.5149, 8.859, 0.3555, 2.149, 0.0045/))
       Magnetic_j4(  5) = Magnetic_Form_Type("TI1 ",(/-1.0383,16.190, 1.4699, 8.924, 0.3631, 2.283, 0.0044/))
       Magnetic_j4(  6) = Magnetic_Form_Type("TI2 ",(/-1.3242,15.310, 1.2042, 7.899, 0.3976, 2.156, 0.0051/))
       Magnetic_j4(  7) = Magnetic_Form_Type("TI3 ",(/-1.1117,14.635, 0.7689, 6.927, 0.4385, 2.089, 0.0060/))
       Magnetic_j4(  8) = Magnetic_Form_Type("V0  ",(/-0.9633,15.273, 0.9274, 7.732, 0.3891, 2.053, 0.0063/))
       Magnetic_j4(  9) = Magnetic_Form_Type("V1  ",(/-0.9606,15.545, 1.1278, 8.118, 0.3653, 2.097, 0.0056/))
       Magnetic_j4( 10) = Magnetic_Form_Type("V2  ",(/-1.1729,14.973, 0.9092, 7.613, 0.4105, 2.039, 0.0067/))
       Magnetic_j4( 11) = Magnetic_Form_Type("V3  ",(/-0.9417,14.205, 0.5284, 6.607, 0.4411, 1.967, 0.0076/))
       Magnetic_j4( 12) = Magnetic_Form_Type("V4  ",(/-0.7654,13.097, 0.3071, 5.674, 0.4476, 1.871, 0.0081/))
       Magnetic_j4( 13) = Magnetic_Form_Type("CR0 ",(/-0.6670,19.613, 0.5342, 6.478, 0.3641, 1.905, 0.0073/))
       Magnetic_j4( 14) = Magnetic_Form_Type("CR1 ",(/-0.8309,18.043, 0.7252, 7.531, 0.3828, 2.003, 0.0073/))
       Magnetic_j4( 15) = Magnetic_Form_Type("CR2 ",(/-0.8930,15.664, 0.5590, 7.033, 0.4093, 1.924, 0.0081/))
       Magnetic_j4( 16) = Magnetic_Form_Type("CR3 ",(/-0.7327,14.073, 0.3268, 5.674, 0.4114, 1.810, 0.0085/))
       Magnetic_j4( 17) = Magnetic_Form_Type("CR4 ",(/-0.6748,12.946, 0.1805, 6.753, 0.4526, 1.800, 0.0098/))
       Magnetic_j4( 18) = Magnetic_Form_Type("MN0 ",(/-0.5452,15.471, 0.4406, 4.902, 0.2884, 1.543, 0.0059/))
       Magnetic_j4( 19) = Magnetic_Form_Type("MN1 ",(/-0.7947,17.867, 0.6078, 7.704, 0.3798, 1.905, 0.0087/))
       Magnetic_j4( 20) = Magnetic_Form_Type("MN2 ",(/-0.7416,15.255, 0.3831, 6.469, 0.3935, 1.800, 0.0093/))
       Magnetic_j4( 21) = Magnetic_Form_Type("MN3 ",(/-0.6603,13.607, 0.2322, 6.218, 0.4104, 1.740, 0.0101/))
       Magnetic_j4( 22) = Magnetic_Form_Type("MN4 ",(/-0.5127,13.461, 0.0313, 7.763, 0.4282, 1.701, 0.0113/))
       Magnetic_j4( 23) = Magnetic_Form_Type("MN5 ",(/0.19236,0.32487,1.67062,6.65663,-1.82036,6.19424,-0.04334/))
       Magnetic_j4( 24) = Magnetic_Form_Type("FE0 ",(/-0.5029,19.677, 0.2999, 3.776, 0.2576, 1.424, 0.0071/))
       Magnetic_j4( 25) = Magnetic_Form_Type("FE1 ",(/-0.5109,19.250, 0.3896, 4.891, 0.2810, 1.526, 0.0069/))
       Magnetic_j4( 26) = Magnetic_Form_Type("FE2 ",(/-0.5401,17.227, 0.2865, 3.742, 0.2658, 1.424, 0.0076/))
       Magnetic_j4( 27) = Magnetic_Form_Type("FE3 ",(/-0.5507,11.493, 0.2153, 4.906, 0.3468, 1.523, 0.0095/))
       Magnetic_j4( 28) = Magnetic_Form_Type("FE4 ",(/-0.5352, 9.507, 0.1783, 5.175, 0.3584, 1.469, 0.0097/))
       Magnetic_j4( 29) = Magnetic_Form_Type("CO0 ",(/-0.4221,14.195, 0.2900, 3.979, 0.2469, 1.286, 0.0063/))
       Magnetic_j4( 30) = Magnetic_Form_Type("CO1 ",(/-0.4115,14.561, 0.3580, 4.717, 0.2644, 1.418, 0.0074/))
       Magnetic_j4( 31) = Magnetic_Form_Type("CO2 ",(/-0.4759,14.046, 0.2747, 3.731, 0.2458, 1.250, 0.0057/))
       Magnetic_j4( 32) = Magnetic_Form_Type("CO3 ",(/-0.4466,13.391, 0.1419, 3.011, 0.2773, 1.335, 0.0093/))
       Magnetic_j4( 33) = Magnetic_Form_Type("CO4 ",(/-0.4091,13.194,-0.0194, 3.417, 0.3534, 1.421, 0.0112/))
       Magnetic_j4( 34) = Magnetic_Form_Type("NI0 ",(/-0.4428,14.485, 0.0870, 3.234, 0.2932, 1.331, 0.0096/))
       Magnetic_j4( 35) = Magnetic_Form_Type("NI1 ",(/-0.3836,13.425, 0.3116, 4.462, 0.2471, 1.309, 0.0079/))
       Magnetic_j4( 36) = Magnetic_Form_Type("NI2 ",(/-0.3803,10.403, 0.2838, 3.378, 0.2108, 1.104, 0.0050/))
       Magnetic_j4( 37) = Magnetic_Form_Type("NI3 ",(/-0.3715, 8.952, 0.1211, 2.940, 0.2526, 1.105, 0.0061/))
       Magnetic_j4( 38) = Magnetic_Form_Type("NI4 ",(/-0.3509, 8.157, 0.2220, 2.106, 0.1567, 0.925, 0.0065/))
       Magnetic_j4( 39) = Magnetic_Form_Type("CU0 ",(/-0.3204,15.132, 0.2335, 4.021, 0.2312, 1.196, 0.0068/))
       Magnetic_j4( 40) = Magnetic_Form_Type("CU1 ",(/-0.3572,15.125, 0.2336, 3.966, 0.2315, 1.197, 0.0070/))
       Magnetic_j4( 41) = Magnetic_Form_Type("CU2 ",(/-0.3914,14.740, 0.1275, 3.384, 0.2548, 1.255, 0.0103/))
       Magnetic_j4( 42) = Magnetic_Form_Type("CU3 ",(/-0.3671,14.082,-0.0078, 3.315, 0.3154, 1.377, 0.0132/))
       Magnetic_j4( 43) = Magnetic_Form_Type("CU4 ",(/-0.2915,14.124,-0.1065, 4.201, 0.3247, 1.352, 0.0148/))
       Magnetic_j4( 44) = Magnetic_Form_Type("Y0  ",(/-8.0767,32.201, 7.9197,25.156, 1.4067, 6.827,-0.0001/))
       Magnetic_j4( 45) = Magnetic_Form_Type("ZR0 ",(/-5.2697,32.868, 4.1930,24.183, 1.5202, 6.048,-0.0002/))
       Magnetic_j4( 46) = Magnetic_Form_Type("ZR1 ",(/-5.6384,33.607, 4.6729,22.338, 1.3258, 5.924,-0.0003/))
       Magnetic_j4( 47) = Magnetic_Form_Type("NB0 ",(/-3.1377,25.595, 2.3411,16.569, 1.2304, 4.990,-0.0005/))
       Magnetic_j4( 48) = Magnetic_Form_Type("NB1 ",(/-3.3598,25.820, 2.8297,16.427, 1.1203, 4.982,-0.0005/))
       Magnetic_j4( 49) = Magnetic_Form_Type("MO0 ",(/-2.8860,20.572, 1.8130,14.628, 1.1899, 4.264,-0.0008/))
       Magnetic_j4( 50) = Magnetic_Form_Type("MO1 ",(/-3.2618,25.486, 2.3596,16.462, 1.1164, 4.491,-0.0007/))
       Magnetic_j4( 51) = Magnetic_Form_Type("TC0 ",(/-2.7975,20.159, 1.6520,16.261, 1.1726, 3.943,-0.0008/))
       Magnetic_j4( 52) = Magnetic_Form_Type("TC1 ",(/-2.0470,19.683, 1.6306,11.592, 0.8698, 3.769,-0.0010/))
       Magnetic_j4( 53) = Magnetic_Form_Type("RU0 ",(/-1.5042,17.949, 0.6027, 9.961, 0.9700, 3.393,-0.0010/))
       Magnetic_j4( 54) = Magnetic_Form_Type("RU1 ",(/-1.6278,18.506, 1.1828,10.189, 0.8138, 3.418,-0.0009/))
       Magnetic_j4( 55) = Magnetic_Form_Type("RH0 ",(/-1.3492,17.577, 0.4527,10.507, 0.9285, 3.155,-0.0009/))
       Magnetic_j4( 56) = Magnetic_Form_Type("RH1 ",(/-1.4673,17.957, 0.7381, 9.944, 0.8485, 3.126,-0.0012/))
       Magnetic_j4( 57) = Magnetic_Form_Type("PD0 ",(/-1.1955,17.628, 0.3183,11.309, 0.8696, 2.909,-0.0006/))
       Magnetic_j4( 58) = Magnetic_Form_Type("PD1 ",(/-1.4098,17.765, 0.7927, 9.999, 0.7710, 2.930,-0.0006/))
       Magnetic_j4( 59) = Magnetic_Form_Type("CE2 ",(/-0.6468,10.533, 0.4052, 5.624, 0.3412, 1.535, 0.0080/))
       Magnetic_j4( 60) = Magnetic_Form_Type("ND2 ",(/-0.5416,12.204, 0.3571, 6.169, 0.3154, 1.485, 0.0098/))
       Magnetic_j4( 61) = Magnetic_Form_Type("ND3 ",(/-0.4053,14.014, 0.0329, 7.005, 0.3759, 1.707, 0.0209/))
       Magnetic_j4( 62) = Magnetic_Form_Type("SM2 ",(/-0.4150,14.057, 0.1368, 7.032, 0.3272, 1.582, 0.0192/))
       Magnetic_j4( 63) = Magnetic_Form_Type("SM3 ",(/-0.4288,10.052, 0.1782, 5.019, 0.2833, 1.236, 0.0088/))
       Magnetic_j4( 64) = Magnetic_Form_Type("EU2 ",(/-0.4145,10.193, 0.2447, 5.164, 0.2661, 1.205, 0.0065/))
       Magnetic_j4( 65) = Magnetic_Form_Type("EU3 ",(/-0.4095,10.211, 0.1485, 5.175, 0.2720, 1.237, 0.0131/))
       Magnetic_j4( 66) = Magnetic_Form_Type("GD2 ",(/-0.3824,10.344, 0.1955, 5.306, 0.2622, 1.203, 0.0097/))
       Magnetic_j4( 67) = Magnetic_Form_Type("GD3 ",(/-0.3621,10.353, 0.1016, 5.310, 0.2649, 1.219, 0.0147/))
       Magnetic_j4( 68) = Magnetic_Form_Type("TB2 ",(/-0.3443,10.469, 0.1481, 5.416, 0.2575, 1.182, 0.0104/))
       Magnetic_j4( 69) = Magnetic_Form_Type("TB3 ",(/-0.3228,10.476, 0.0638, 5.419, 0.2566, 1.196, 0.0159/))
       Magnetic_j4( 70) = Magnetic_Form_Type("DY2 ",(/-0.3206,12.071, 0.0904, 8.026, 0.2616, 1.230, 0.0143/))
       Magnetic_j4( 71) = Magnetic_Form_Type("DY3 ",(/-0.2829, 9.525, 0.0565, 4.429, 0.2437, 1.066, 0.0092/))
       Magnetic_j4( 72) = Magnetic_Form_Type("HO2 ",(/-0.2976, 9.719, 0.1224, 4.635, 0.2279, 1.005, 0.0063/))
       Magnetic_j4( 73) = Magnetic_Form_Type("HO3 ",(/-0.2717, 9.731, 0.0474, 4.638, 0.2292, 1.047, 0.0124/))
       Magnetic_j4( 74) = Magnetic_Form_Type("ER2 ",(/-0.2975, 9.829, 0.1189, 4.741, 0.2116, 1.004, 0.0117/))
       Magnetic_j4( 75) = Magnetic_Form_Type("ER3 ",(/-0.2568, 9.834, 0.0356, 4.741, 0.2172, 1.028, 0.0148/))
       Magnetic_j4( 76) = Magnetic_Form_Type("TM2 ",(/-0.2677, 9.888, 0.0925, 4.784, 0.2056, 0.990, 0.0124/))
       Magnetic_j4( 77) = Magnetic_Form_Type("TM3 ",(/-0.2292, 9.895, 0.0124, 4.785, 0.2108, 1.007, 0.0151/))
       Magnetic_j4( 78) = Magnetic_Form_Type("YB2 ",(/-0.2393, 9.947, 0.0663, 4.823, 0.2009, 0.965, 0.0122/))
       Magnetic_j4( 79) = Magnetic_Form_Type("YB3 ",(/-0.2121, 8.197, 0.0325, 3.153, 0.1975, 0.884, 0.0093/))
       Magnetic_j4( 80) = Magnetic_Form_Type("U3  ",(/-0.9859,16.601, 0.6116, 6.515, 0.6020, 2.597,-0.0010/))
       Magnetic_j4( 81) = Magnetic_Form_Type("U4  ",(/-1.0540,16.605, 0.4339, 6.512, 0.6746, 2.599,-0.0011/))
       Magnetic_j4( 82) = Magnetic_Form_Type("U5  ",(/-0.9588,16.485, 0.1576, 6.440, 0.7785, 2.640,-0.0010/))
       Magnetic_j4( 83) = Magnetic_Form_Type("NP3 ",(/-0.9029,16.586, 0.4006, 6.470, 0.6545, 2.563,-0.0004/))
       Magnetic_j4( 84) = Magnetic_Form_Type("NP4 ",(/-0.9887,12.441, 0.5918, 5.294, 0.5306, 2.263,-0.0021/))
       Magnetic_j4( 85) = Magnetic_Form_Type("NP5 ",(/-0.8146,16.581,-0.0055, 6.475, 0.7956, 2.562,-0.0004/))
       Magnetic_j4( 86) = Magnetic_Form_Type("NP6 ",(/-0.6738,16.553,-0.2297, 6.505, 0.8513, 2.553,-0.0003/))
       Magnetic_j4( 87) = Magnetic_Form_Type("PU3 ",(/-0.7014,16.369,-0.1162, 6.697, 0.7778, 2.450, 0.0000/))
       Magnetic_j4( 88) = Magnetic_Form_Type("PU4 ",(/-0.9160,12.203, 0.4891, 5.127, 0.5290, 2.149,-0.0022/))
       Magnetic_j4( 89) = Magnetic_Form_Type("PU5 ",(/-0.7035,16.360,-0.0979, 6.706, 0.7726, 2.447, 0.0000/))
       Magnetic_j4( 90) = Magnetic_Form_Type("PU6 ",(/-0.5560,16.322,-0.3046, 6.768, 0.8146, 2.426, 0.0001/))
       Magnetic_j4( 91) = Magnetic_Form_Type("AM2 ",(/-0.7433,16.416, 0.3481, 6.788, 0.6014, 2.346, 0.0000/))
       Magnetic_j4( 92) = Magnetic_Form_Type("AM3 ",(/-0.8092,12.854, 0.4161, 5.459, 0.5476, 2.172,-0.0011/))
       Magnetic_j4( 93) = Magnetic_Form_Type("AM4 ",(/-0.8548,12.226, 0.3037, 5.909, 0.6173, 2.188,-0.0016/))
       Magnetic_j4( 94) = Magnetic_Form_Type("AM5 ",(/-0.6538,15.462,-0.0948, 5.997, 0.7295, 2.297, 0.0000/))
       Magnetic_j4( 95) = Magnetic_Form_Type("AM6 ",(/-0.5390,15.449,-0.2689, 6.017, 0.7711, 2.297, 0.0002/))
       Magnetic_j4( 96) = Magnetic_Form_Type("AM7 ",(/-0.4688,12.019,-0.2692, 7.042, 0.7297, 2.164,-0.0011/))
       Magnetic_j4( 97) = Magnetic_Form_Type("PR3 ",(/-0.3970,10.9919, 0.0818, 5.9897, 0.3656, 1.5021, 0.0110/))

       !---- <j6> Coefficients ----!
       Magnetic_j6(  1) = Magnetic_Form_Type("CE2 ",(/-0.1212, 7.994,-0.0639, 4.024, 0.1519, 1.096, 0.0078/))
       Magnetic_j6(  2) = Magnetic_Form_Type("ND2 ",(/-0.1600, 8.009, 0.0272, 4.028, 0.1104, 1.068, 0.0139/))
       Magnetic_j6(  3) = Magnetic_Form_Type("ND3 ",(/-0.0416, 8.014,-0.1261, 4.040, 0.1400, 1.087, 0.0102/))
       Magnetic_j6(  4) = Magnetic_Form_Type("SM2 ",(/-0.1428, 6.041, 0.0723, 2.033, 0.0550, 0.513, 0.0081/))
       Magnetic_j6(  5) = Magnetic_Form_Type("SM3 ",(/-0.0944, 6.030,-0.0498, 2.074, 0.1372, 0.645,-0.0132/))
       Magnetic_j6(  6) = Magnetic_Form_Type("EU2 ",(/-0.1252, 6.049, 0.0507, 2.085, 0.0572, 0.646, 0.0132/))
       Magnetic_j6(  7) = Magnetic_Form_Type("EU3 ",(/-0.0817, 6.039,-0.0596, 2.120, 0.1243, 0.764,-0.0001/))
       Magnetic_j6(  8) = Magnetic_Form_Type("GD2 ",(/-0.1351, 5.030, 0.0828, 2.025, 0.0315, 0.503, 0.0187/))
       Magnetic_j6(  9) = Magnetic_Form_Type("GD3 ",(/-0.0662, 6.031,-0.0850, 2.154, 0.1323, 0.891, 0.0048/))
       Magnetic_j6( 10) = Magnetic_Form_Type("TB2 ",(/-0.0758, 6.032,-0.0540, 2.158, 0.1199, 0.890, 0.0051/))
       Magnetic_j6( 11) = Magnetic_Form_Type("TB3 ",(/-0.0559, 6.031,-0.1020, 2.237, 0.1264, 1.107, 0.0167/))
       Magnetic_j6( 12) = Magnetic_Form_Type("DY2 ",(/-0.0568, 6.032,-0.1003, 2.240, 0.1401, 1.106, 0.0109/))
       Magnetic_j6( 13) = Magnetic_Form_Type("DY3 ",(/-0.0423, 6.038,-0.1248, 2.244, 0.1359, 1.200, 0.0188/))
       Magnetic_j6( 14) = Magnetic_Form_Type("HO2 ",(/-0.0725, 6.045,-0.0318, 2.243, 0.0738, 1.202, 0.0252/))
       Magnetic_j6( 15) = Magnetic_Form_Type("HO3 ",(/-0.0289, 6.050,-0.1545, 2.230, 0.1550, 1.260, 0.0177/))
       Magnetic_j6( 16) = Magnetic_Form_Type("ER2 ",(/-0.0648, 6.056,-0.0515, 2.230, 0.0825, 1.264, 0.0250/))
       Magnetic_j6( 17) = Magnetic_Form_Type("ER3 ",(/-0.0110, 6.061,-0.1954, 2.224, 0.1818, 1.296, 0.0149/))
       Magnetic_j6( 18) = Magnetic_Form_Type("TM2 ",(/-0.0842, 4.070, 0.0807, 0.849,-0.2087, 0.039, 0.2095/))
       Magnetic_j6( 19) = Magnetic_Form_Type("TM3 ",(/-0.0727, 4.073, 0.0243, 0.689, 3.9459, 0.002,-3.9076/))
       Magnetic_j6( 20) = Magnetic_Form_Type("YB2 ",(/-0.0739, 5.031, 0.0140, 2.030, 0.0351, 0.508, 0.0174/))
       Magnetic_j6( 21) = Magnetic_Form_Type("YB3 ",(/-0.0345, 5.007,-0.0677, 2.020, 0.0985, 0.549,-0.0076/))
       Magnetic_j6( 22) = Magnetic_Form_Type("U3  ",(/-0.3797, 9.953, 0.0459, 5.038, 0.2748, 1.607, 0.0016/))
       Magnetic_j6( 23) = Magnetic_Form_Type("U4  ",(/-0.1793,11.896,-0.2269, 5.428, 0.3291, 1.701, 0.0030/))
       Magnetic_j6( 24) = Magnetic_Form_Type("U5  ",(/-0.0399,11.891,-0.3458, 5.580, 0.3340, 1.645, 0.0029/))
       Magnetic_j6( 25) = Magnetic_Form_Type("NP3 ",(/-0.2427,11.844,-0.1129, 5.377, 0.2848, 1.568, 0.0022/))
       Magnetic_j6( 26) = Magnetic_Form_Type("NP4 ",(/-0.2436, 9.599,-0.1317, 4.101, 0.3029, 1.545, 0.0019/))
       Magnetic_j6( 27) = Magnetic_Form_Type("NP5 ",(/-0.1157, 9.565,-0.2654, 4.260, 0.3298, 1.549, 0.0025/))
       Magnetic_j6( 28) = Magnetic_Form_Type("NP6 ",(/-0.0128, 9.569,-0.3611, 4.304, 0.3419, 1.541, 0.0032/))
       Magnetic_j6( 29) = Magnetic_Form_Type("PU3 ",(/-0.0364, 9.572,-0.3181, 4.342, 0.3210, 1.523, 0.0041/))
       Magnetic_j6( 30) = Magnetic_Form_Type("PU4 ",(/-0.2394, 7.837,-0.0785, 4.024, 0.2643, 1.378, 0.0012/))
       Magnetic_j6( 31) = Magnetic_Form_Type("PU5 ",(/-0.1090, 7.819,-0.2243, 4.100, 0.2947, 1.404, 0.0015/))
       Magnetic_j6( 32) = Magnetic_Form_Type("PU6 ",(/-0.0001, 7.820,-0.3354, 4.144, 0.3097, 1.403, 0.0020/))
       Magnetic_j6( 33) = Magnetic_Form_Type("AM2 ",(/-0.3176, 7.864, 0.0771, 4.161, 0.2194, 1.339, 0.0018/))
       Magnetic_j6( 34) = Magnetic_Form_Type("AM3 ",(/-0.3159, 6.982, 0.0682, 3.995, 0.2141, 1.188,-0.0015/))
       Magnetic_j6( 35) = Magnetic_Form_Type("AM4 ",(/-0.1787, 7.880,-0.1274, 4.090, 0.2565, 1.315, 0.0017/))
       Magnetic_j6( 36) = Magnetic_Form_Type("AM5 ",(/-0.0927, 6.073,-0.2227, 3.784, 0.2916, 1.372, 0.0026/))
       Magnetic_j6( 37) = Magnetic_Form_Type("AM6 ",(/ 0.0152, 6.079,-0.3549, 3.861, 0.3125, 1.403, 0.0036/))
       Magnetic_j6( 38) = Magnetic_Form_Type("AM7 ",(/ 0.1292, 6.082,-0.4689, 3.879, 0.3234, 1.393, 0.0042/))
       Magnetic_j6( 39) = Magnetic_Form_Type("PR3 ",(/-0.0224, 7.9931,-0.1202, 3.9406, 0.1299, 0.8938, 0.0051/))

       return
    End Subroutine Set_Magnetic_Form

    !!----
    !!---- Subroutine Set_Xray_Form()
    !!----    Set Xray_Form Array:
    !!--<<
    !!----        1: Symbol of the Element
    !!----        2: Name of the Element
    !!----        3: a(4)
    !!----        4: b(4)
    !!----        5: c
    !!----    Coefficients for calculating the X-ray scattering factors
    !!----        f(s) = Sum_{i=1,4} { a(i) exp(-b(i)*s^2) } + c
    !!----
    !!----    where s=sinTheta/Lambda
    !!-->>
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Set_Xray_Form()

       if (.not. allocated(xray_form)) allocate(xray_form(num_xray_form))

       Xray_form( 1:10) = (/ &
                          xray_form_type("h   ",  1, (/  0.493002,   0.322912,   0.140191,   0.040810/), &
                                                     (/ 10.510900,  26.125700,   3.142360,  57.799698/),  0.003038) , &
                          xray_form_type("h-1 ",  1, (/  0.897661,   0.565616,   0.415815,   0.116973/), &
                                                     (/ 53.136799,  15.187000, 186.575989,   3.567090/),  0.002389) , &
                          xray_form_type("he  ",  2, (/  0.873400,   0.630900,   0.311200,   0.178000/), &
                                                     (/  9.103700,   3.356800,  22.927601,   0.982100/),  0.006400) , &
                          xray_form_type("li  ",  3, (/  1.128200,   0.750800,   0.617500,   0.465300/), &
                                                     (/  3.954600,   1.052400,  85.390503, 168.261002/),  0.037700) , &
                          xray_form_type("li+1",  3, (/  0.696800,   0.788800,   0.341400,   0.156300/), &
                                                     (/  4.623700,   1.955700,   0.631600,  10.095300/),  0.016700) , &
                          xray_form_type("be  ",  4, (/  1.591900,   1.127800,   0.539100,   0.702900/), &
                                                     (/ 43.642700,   1.862300, 103.483002,   0.542000/),  0.038500) , &
                          xray_form_type("be+2",  4, (/  6.260300,   0.884900,   0.799300,   0.164700/), &
                                                     (/  0.002700,   0.831300,   2.275800,   5.114600/), -6.109200) , &
                          xray_form_type("b   ",  5, (/  2.054500,   1.332600,   1.097900,   0.706800/), &
                                                     (/ 23.218500,   1.021000,  60.349800,   0.140300/), -0.193200) , &
                          xray_form_type("c   ",  6, (/  2.310000,   1.020000,   1.588600,   0.865000/), &
                                                     (/ 20.843899,  10.207500,   0.568700,  51.651199/),  0.215600) , &
                          xray_form_type("cv  ",  6, (/  2.260690,   1.561650,   1.050750,   0.839259/), &
                                                     (/ 22.690701,   0.656665,   9.756180,  55.594898/),  0.286977) /)

       Xray_form(11:20) = (/ &
                          xray_form_type("n   ",  7, (/ 12.212600,   3.132200,   2.012500,   1.166300/), &
                                                     (/  0.005700,   9.893300,  28.997499,   0.582600/),-11.528999) , &
                          xray_form_type("o   ",  8, (/  3.048500,   2.286800,   1.546300,   0.867000/), &
                                                     (/ 13.277100,   5.701100,   0.323900,  32.908897/),  0.250800) , &
                          xray_form_type("o-1 ",  8, (/  4.191600,   1.639690,   1.526730, -20.306999/), &
                                                     (/ 12.857300,   4.172360,  47.017899,  -0.014040/), 21.941200) , &
                          xray_form_type("f   ",  9, (/  3.539200,   2.641200,   1.517000,   1.024300/), &
                                                     (/ 10.282499,   4.294400,   0.261500,  26.147600/),  0.277600) , &
                          xray_form_type("f-1 ",  9, (/  3.632200,   3.510570,   1.260640,   0.940706/), &
                                                     (/  5.277560,  14.735300,   0.442258,  47.343700/),  0.653396) , &
                          xray_form_type("ne  ", 10, (/  3.955300,   3.112500,   1.454600,   1.125100/), &
                                                     (/  8.404200,   3.426200,   0.230600,  21.718399/),  0.351500) , &
                          xray_form_type("na  ", 11, (/  4.762600,   3.173600,   1.267400,   1.112800/), &
                                                     (/  3.285000,   8.842199,   0.313600, 129.423996/),  0.676000) , &
                          xray_form_type("na+1", 11, (/  3.256500,   3.936200,   1.399800,   1.003200/), &
                                                     (/  2.667100,   6.115300,   0.200100,  14.039000/),  0.404000) , &
                          xray_form_type("mg  ", 12, (/  5.420400,   2.173500,   1.226900,   2.307300/), &
                                                     (/  2.827500,  79.261101,   0.380800,   7.193700/),  0.858400) , &
                          xray_form_type("mg+2", 12, (/  3.498800,   3.837800,   1.328400,   0.849700/), &
                                                     (/  2.167600,   4.754200,   0.185000,  10.141100/),  0.485300) /)

       Xray_form(21:30) = (/ &
                          xray_form_type("al  ", 13, (/  6.420200,   1.900200,   1.593600,   1.964600/), &
                                                     (/  3.038700,   0.742600,  31.547199,  85.088600/),  1.115100) , &
                          xray_form_type("al+3", 13, (/  4.174480,   3.387600,   1.202960,   0.528137/), &
                                                     (/  1.938160,   4.145530,   0.228753,   8.285240/),  0.706786) , &
                          xray_form_type("si  ", 14, (/  6.291500,   3.035300,   1.989100,   1.541000/), &
                                                     (/  2.438600,  32.333698,   0.678500,  81.693695/),  1.140700) , &
                          xray_form_type("siv ", 14, (/  5.662690,   3.071640,   2.624460,   1.393200/), &
                                                     (/  2.665200,  38.663399,   0.916946,  93.545799/),  1.247070) , &
                          xray_form_type("si+4", 14, (/  4.439180,   3.203450,   1.194530,   0.416530/), &
                                                     (/  1.641670,   3.437570,   0.214900,   6.653650/),  0.746297) , &
                          xray_form_type("p   ", 15, (/  6.434500,   4.179100,   1.780000,   1.490800/), &
                                                     (/  1.906700,  27.157000,   0.526000,  68.164497/),  1.114900) , &
                          xray_form_type("s   ", 16, (/  6.905300,   5.203400,   1.437900,   1.586300/), &
                                                     (/  1.467900,  22.215099,   0.253600,  56.172001/),  0.866900) , &
                          xray_form_type("cl  ", 17, (/ 11.460400,   7.196400,   6.255600,   1.645500/), &
                                                     (/  0.010400,   1.166200,  18.519400,  47.778400/), -9.557400) , &
                          xray_form_type("cl-1", 17, (/ 18.291500,   7.208400,   6.533700,   2.338600/), &
                                                     (/  0.006600,   1.171700,  19.542400,  60.448601/),-16.378000) , &
                          xray_form_type("ar  ", 18, (/  7.484500,   6.772300,   0.653900,   1.644200/), &
                                                     (/  0.907200,  14.840700,  43.898300,  33.392899/),  1.444500) /)

       Xray_form(31:40) = (/ &
                          xray_form_type("k   ", 19, (/  8.218599,   7.439800,   1.051900,   0.865900/), &
                                                     (/ 12.794900,   0.774800, 213.186996,  41.684097/),  1.422800) , &
                          xray_form_type("k+1 ", 19, (/  7.957800,   7.491700,   6.359000,   1.191500/), &
                                                     (/ 12.633100,   0.767400,  -0.002000,  31.912800/), -4.997800) , &
                          xray_form_type("ca  ", 20, (/  8.626600,   7.387300,   1.589900,   1.021100/), &
                                                     (/ 10.442100,   0.659900,  85.748398, 178.436996/),  1.375100) , &
                          xray_form_type("ca+2", 20, (/ 15.634800,   7.951800,   8.437200,   0.853700/), &
                                                     (/ -0.007400,   0.608900,  10.311600,  25.990499/),-14.875000) , &
                          xray_form_type("sc  ", 21, (/  9.189000,   7.367900,   1.640900,   1.468000/), &
                                                     (/  9.021299,   0.572900, 136.108002,  51.353100/),  1.332900) , &
                          xray_form_type("sc+3", 21, (/ 14.400800,   8.027300,   1.659430,   1.579360/), &
                                                     (/  0.298540,   7.962900,  -0.286040,  16.066200/), -6.666700) , &
                          xray_form_type("ti  ", 22, (/  9.759500,   7.355800,   1.699100,   1.902100/), &
                                                     (/  7.850800,   0.500000,  35.633801, 116.104996/),  1.280700) , &
                          xray_form_type("ti+2", 22, (/  9.114230,   7.621740,   2.279300,   0.087899/), &
                                                     (/  7.524300,   0.457585,  19.536100,  61.655800/),  0.897155) , &
                          xray_form_type("ti+3", 22, (/ 17.734400,   8.738160,   5.256910,   1.921340/), &
                                                     (/  0.220610,   7.047160,  -0.157620,  15.976800/),-14.652000) , &
                          xray_form_type("ti+4", 22, (/ 19.511400,   8.234730,   2.013410,   1.520800/), &
                                                     (/  0.178847,   6.670180,  -0.292630,  12.946400/),-13.280000) /)

       Xray_form(41:50) = (/ &
                          xray_form_type("v   ", 23, (/ 10.297100,   7.351100,   2.070300,   2.057100/), &
                                                     (/  6.865700,   0.438500,  26.893799, 102.477997/),  1.219900) , &
                          xray_form_type("v+2 ", 23, (/ 10.106000,   7.354100,   2.288400,   0.022300/), &
                                                     (/  6.881800,   0.440900,  20.300400, 115.122002/),  1.229800) , &
                          xray_form_type("v+3 ", 23, (/  9.431410,   7.741900,   2.153430,   0.016865/), &
                                                     (/  6.395350,   0.383349,  15.190800,  63.969002/),  0.656565) , &
                          xray_form_type("v+5 ", 23, (/ 15.688700,   8.142080,   2.030810,  -9.576000/), &
                                                     (/  0.679003,   5.401350,   9.972780,   0.940464/),  1.714300) , &
                          xray_form_type("cr  ", 24, (/ 10.640600,   7.353700,   3.324000,   1.492200/), &
                                                     (/  6.103800,   0.392000,  20.262600,  98.739899/),  1.183200) , &
                          xray_form_type("cr+2", 24, (/  9.540340,   7.750900,   3.582740,   0.509107/), &
                                                     (/  5.660780,   0.344261,  13.307500,  32.422401/),  0.616898) , &
                          xray_form_type("cr+3", 24, (/  9.680900,   7.811360,   2.876030,   0.113575/), &
                                                     (/  5.594630,   0.334393,  12.828800,  32.876099/),  0.518275) , &
                          xray_form_type("mn  ", 25, (/ 11.281900,   7.357300,   3.019300,   2.244100/), &
                                                     (/  5.340900,   0.343200,  17.867399,  83.754303/),  1.089600) , &
                          xray_form_type("mn+2", 25, (/ 10.806100,   7.362000,   3.526800,   0.218400/), &
                                                     (/  5.279600,   0.343500,  14.343000,  41.323502/),  1.087400) , &
                          xray_form_type("mn+3", 25, (/  9.845210,   7.871940,   3.565310,   0.323613/), &
                                                     (/  4.917970,   0.294393,  10.817100,  24.128099/),  0.393974) /)

       Xray_form(51:60) = (/ &
                          xray_form_type("mn+4", 25, (/  9.962530,   7.970570,   2.760670,   0.054447/), &
                                                     (/  4.848500,   0.283303,  10.485200,  27.573000/),  0.251877) , &
                          xray_form_type("fe  ", 26, (/ 11.769500,   7.357300,   3.522200,   2.304500/), &
                                                     (/  4.761100,   0.307200,  15.353500,  76.880501/),  1.036900) , &
                          xray_form_type("fe+2", 26, (/ 11.042400,   7.374000,   4.134600,   0.439900/), &
                                                     (/  4.653800,   0.305300,  12.054600,  31.280899/),  1.009700) , &
                          xray_form_type("fe+3", 26, (/ 11.176400,   7.386300,   3.394800,   0.072400/), &
                                                     (/  4.614700,   0.300500,  11.672900,  38.556599/),  0.970700) , &
                          xray_form_type("co  ", 27, (/ 12.284100,   7.340900,   4.003400,   2.348800/), &
                                                     (/  4.279100,   0.278400,  13.535900,  71.169197/),  1.011800) , &
                          xray_form_type("co+2", 27, (/ 11.229600,   7.388300,   4.739300,   0.710800/), &
                                                     (/  4.123100,   0.272600,  10.244300,  25.646599/),  0.932400) , &
                          xray_form_type("co+3", 27, (/ 10.337999,   7.881730,   4.767950,   0.725591/), &
                                                     (/  3.909690,   0.238668,   8.355830,  18.349100/),  0.286667) , &
                          xray_form_type("ni  ", 28, (/ 12.837600,   7.292000,   4.443800,   2.380000/), &
                                                     (/  3.878500,   0.256500,  12.176300,  66.342102/),  1.034100) , &
                          xray_form_type("ni+2", 28, (/ 11.416600,   7.400500,   5.344200,   0.977300/), &
                                                     (/  3.676600,   0.244900,   8.873000,  22.162600/),  0.861400) , &
                          xray_form_type("ni+3", 28, (/ 10.780600,   7.758680,   5.227460,   0.847114/), &
                                                     (/  3.547700,   0.223140,   7.644680,  16.967300/),  0.386044) /)

       Xray_form(61:70) = (/ &
                          xray_form_type("cu  ", 29, (/ 13.337999,   7.167600,   5.615800,   1.673500/), &
                                                     (/  3.582800,   0.247000,  11.396600,  64.812599/),  1.191000) , &
                          xray_form_type("cu+1", 29, (/ 11.947500,   7.357300,   6.245500,   1.557800/), &
                                                     (/  3.366900,   0.227400,   8.662500,  25.848700/),  0.890000) , &
                          xray_form_type("cu+2", 29, (/ 11.816800,   7.111810,   5.781350,   1.145230/), &
                                                     (/  3.374840,   0.244078,   7.987600,  19.896999/),  1.144310) , &
                          xray_form_type("zn  ", 30, (/ 14.074300,   7.031800,   5.162500,   2.410000/), &
                                                     (/  3.265500,   0.233300,  10.316299,  58.709702/),  1.304100) , &
                          xray_form_type("zn+2", 30, (/ 11.971900,   7.386200,   6.466800,   1.394000/), &
                                                     (/  2.994600,   0.203100,   7.082600,  18.099499/),  0.780700) , &
                          xray_form_type("ga  ", 31, (/ 15.235400,   6.700600,   4.359100,   2.962300/), &
                                                     (/  3.066900,   0.241200,  10.780500,  61.413498/),  1.718900) , &
                          xray_form_type("ga+3", 31, (/ 12.691999,   6.698830,   6.066920,   1.006600/), &
                                                     (/  2.812620,   0.227890,   6.364410,  14.412200/),  1.535450) , &
                          xray_form_type("ge  ", 32, (/ 16.081600,   6.374700,   3.706800,   3.683000/), &
                                                     (/  2.850900,   0.251600,  11.446800,  54.762501/),  2.131300) , &
                          xray_form_type("ge+4", 32, (/ 12.917200,   6.700030,   6.067910,   0.859041/), &
                                                     (/  2.537180,   0.205855,   5.479130,  11.603000/),  1.455720) , &
                          xray_form_type("as  ", 33, (/ 16.672300,   6.070100,   3.431300,   4.277900/), &
                                                     (/  2.634500,   0.264700,  12.947900,  47.797199/),  2.531000) /)

       Xray_form(71:80) = (/ &
                          xray_form_type("se  ", 34, (/ 17.000599,   5.819600,   3.973100,   4.354300/), &
                                                     (/  2.409800,   0.272600,  15.237200,  43.816299/),  2.840900) , &
                          xray_form_type("br  ", 35, (/ 17.178900,   5.235800,   5.637700,   3.985100/), &
                                                     (/  2.172300,  16.579599,   0.260900,  41.432800/),  2.955700) , &
                          xray_form_type("br-1", 35, (/ 17.171799,   6.333800,   5.575400,   3.727200/), &
                                                     (/  2.205900,  19.334499,   0.287100,  58.153500/),  3.177600) , &
                          xray_form_type("kr  ", 36, (/ 17.355499,   6.728600,   5.549300,   3.537500/), &
                                                     (/  1.938400,  16.562300,   0.226100,  39.397202/),  2.825000) , &
                          xray_form_type("rb  ", 37, (/ 17.178400,   9.643499,   5.139900,   1.529200/), &
                                                     (/  1.788800,  17.315100,   0.274800, 164.933990/),  3.487300) , &
                          xray_form_type("rb+1", 37, (/ 17.581600,   7.659800,   5.898100,   2.781700/), &
                                                     (/  1.713900,  14.795700,   0.160300,  31.208700/),  2.078200) , &
                          xray_form_type("sr  ", 38, (/ 17.566299,   9.818399,   5.422000,   2.669400/), &
                                                     (/  1.556400,  14.098800,   0.166400, 132.376007/),  2.506400) , &
                          xray_form_type("sr+2", 38, (/ 18.087400,   8.137300,   2.565400, -34.193001/), &
                                                     (/  1.490700,  12.696300,  24.565100,  -0.013800/), 41.402500) , &
                          xray_form_type("y   ", 39, (/ 17.775999,  10.294600,   5.726290,   3.265880/), &
                                                     (/  1.402900,  12.800600,   0.125599, 104.353996/),  1.912130) , &
                          xray_form_type("y+3 ", 39, (/ 17.926800,   9.153100,   1.767950, -33.108002/), &
                                                     (/  1.354170,  11.214500,  22.659901,  -0.013190/), 40.260201) /)

       Xray_form(81:90) = (/ &
                          xray_form_type("zr  ", 40, (/ 17.876499,  10.948000,   5.417320,   3.657210/), &
                                                     (/  1.276180,  11.916000,   0.117622,  87.662697/),  2.069290) , &
                          xray_form_type("zr+4", 40, (/ 18.166800,  10.056200,   1.011180,  -2.647900/), &
                                                     (/  1.214800,  10.148300,  21.605400,  -0.102760/),  9.414539) , &
                          xray_form_type("nb  ", 41, (/ 17.614201,  12.014400,   4.041830,   3.533460/), &
                                                     (/  1.188650,  11.766000,   0.204785,  69.795700/),  3.755910) , &
                          xray_form_type("nb+3", 41, (/ 19.881199,  18.065300,  11.017700,   1.947150/), &
                                                     (/  0.019175,   1.133050,  10.162100,  28.338900/),-12.912000) , &
                          xray_form_type("nb+5", 41, (/ 17.916300,  13.341700,  10.799000,   0.337905/), &
                                                     (/  1.124460,   0.028781,   9.282060,  25.722799/), -6.393400) , &
                          xray_form_type("mo  ", 42, (/  3.702500,  17.235600,  12.887600,   3.742900/), &
                                                     (/  0.277200,   1.095800,  11.004000,  61.658401/),  4.387500) , &
                          xray_form_type("mo+3", 42, (/ 21.166401,  18.201700,  11.742300,   2.309510/), &
                                                     (/  0.014734,   1.030310,   9.536590,  26.630699/),-14.421000) , &
                          xray_form_type("mo+5", 42, (/ 21.014900,  18.099199,  11.463200,   0.740625/), &
                                                     (/  0.014345,   1.022380,   8.788090,  23.345200/),-14.316000) , &
                          xray_form_type("mo+6", 42, (/ 17.887100,  11.175000,   6.578910,   0.000000/), &
                                                     (/  1.036490,   8.480610,   0.058881,   0.000000/),  0.344941) , &
                          xray_form_type("tc  ", 43, (/ 19.130100,  11.094800,   4.649010,   2.712630/), &
                                                     (/  0.864132,   8.144870,  21.570700,  86.847198/),  5.404280) /)

       Xray_form(91:100)= (/ &
                          xray_form_type("ru  ", 44, (/ 19.267399,  12.918200,   4.863370,   1.567560/), &
                                                     (/  0.808520,   8.434669,  24.799700,  94.292801/),  5.378740) , &
                          xray_form_type("ru+3", 44, (/ 18.563801,  13.288500,   9.326019,   3.009640/), &
                                                     (/  0.847329,   8.371640,   0.017662,  22.886999/), -3.189200) , &
                          xray_form_type("ru+4", 44, (/ 18.500299,  13.178699,   4.713040,   2.185350/), &
                                                     (/  0.844582,   8.125340,   0.364950,  20.850399/),  1.423570) , &
                          xray_form_type("rh  ", 45, (/ 19.295700,  14.350100,   4.734250,   1.289180/), &
                                                     (/  0.751536,   8.217580,  25.874901,  98.606201/),  5.328000) , &
                          xray_form_type("rh+3", 45, (/ 18.878500,  14.125900,   3.325150,  -6.198900/), &
                                                     (/  0.764252,   7.844380,  21.248699,  -0.010360/), 11.867800) , &
                          xray_form_type("rh+4", 45, (/ 18.854500,  13.980600,   2.534640,  -5.652600/), &
                                                     (/  0.760825,   7.624360,  19.331699,  -0.010200/), 11.283500) , &
                          xray_form_type("pd  ", 46, (/ 19.331900,  15.501699,   5.295370,   0.605844/), &
                                                     (/  0.698655,   7.989290,  25.205200,  76.898598/),  5.265930) , &
                          xray_form_type("pd+2", 46, (/ 19.170099,  15.209600,   4.322340,   0.000000/), &
                                                     (/  0.696219,   7.555730,  22.505699,   0.000000/),  5.291600) , &
                          xray_form_type("pd+4", 46, (/ 19.249300,  14.790000,   2.892890,  -7.949200/), &
                                                     (/  0.683839,   7.148330,  17.914400,   0.005127/), 13.017400) , &
                          xray_form_type("ag  ", 47, (/ 19.280800,  16.688499,   4.804500,   1.046300/), &
                                                     (/  0.644600,   7.472600,  24.660500,  99.815598/),  5.179000) /)

       Xray_form(101:110)=(/ &
                          xray_form_type("ag+1", 47, (/ 19.181200,  15.971900,   5.274750,   0.357534/), &
                                                     (/  0.646179,   7.191230,  21.732599,  66.114700/),  5.215720) , &
                          xray_form_type("ag+2", 47, (/ 19.164299,  16.245600,   4.370900,   0.000000/), &
                                                     (/  0.645643,   7.185440,  21.407200,   0.000000/),  5.214040) , &
                          xray_form_type("cd  ", 48, (/ 19.221399,  17.644400,   4.461000,   1.602900/), &
                                                     (/  0.594600,   6.908900,  24.700800,  87.482498/),  5.069400) , &
                          xray_form_type("cd+2", 48, (/ 19.151400,  17.253500,   4.471280,   0.000000/), &
                                                     (/  0.597922,   6.806390,  20.252100,   0.000000/),  5.119370) , &
                          xray_form_type("in  ", 49, (/ 19.162399,  18.559601,   4.294800,   2.039600/), &
                                                     (/  0.547600,   6.377600,  25.849899,  92.802902/),  4.939100) , &
                          xray_form_type("in+3", 49, (/ 19.104500,  18.110800,   3.788970,   0.000000/), &
                                                     (/  0.551522,   6.324700,  17.359501,   0.000000/),  4.996350) , &
                          xray_form_type("sn  ", 50, (/ 19.188900,  19.100500,   4.458500,   2.466300/), &
                                                     (/  5.830300,   0.503100,  26.890900,  83.957100/),  4.782100) , &
                          xray_form_type("sn+2", 50, (/ 19.109400,  19.054800,   4.564800,   0.487000/), &
                                                     (/  0.503600,   5.837800,  23.375200,  62.206100/),  4.786100) , &
                          xray_form_type("sn+4", 50, (/ 18.933300,  19.713100,   3.418200,   0.019300/), &
                                                     (/  5.764000,   0.465500,  14.004900,  -0.758300/),  3.918200) , &
                          xray_form_type("sb  ", 51, (/ 19.641800,  19.045500,   5.037100,   2.682700/), &
                                                     (/  5.303400,   0.460700,  27.907400,  75.282501/),  4.590900) /)

       Xray_form(111:120)=(/ &
                          xray_form_type("sb+3", 51, (/ 18.975500,  18.932999,   5.107890,   0.288753/), &
                                                     (/  0.467196,   5.221260,  19.590200,  55.511299/),  4.696260) , &
                          xray_form_type("sb+5", 51, (/ 19.868500,  19.030199,   2.412530,   0.000000/), &
                                                     (/  5.448530,   0.467973,  14.125900,   0.000000/),  4.692630) , &
                          xray_form_type("te  ", 52, (/ 19.964399,  19.013800,   6.144870,   2.523900/), &
                                                     (/  4.817420,   0.420885,  28.528400,  70.840302/),  4.352000) , &
                          xray_form_type("i   ", 53, (/ 20.147200,  18.994900,   7.513800,   2.273500/), &
                                                     (/  4.347000,   0.381400,  27.765999,  66.877602/),  4.071200) , &
                          xray_form_type("i-1 ", 53, (/ 20.233200,  18.997000,   7.806900,   2.886800/), &
                                                     (/  4.357900,   0.381500,  29.525900,  84.930397/),  4.071400) , &
                          xray_form_type("xe  ", 54, (/ 20.293301,  19.029800,   8.976700,   1.990000/), &
                                                     (/  3.928200,   0.344000,  26.465900,  64.265800/),  3.711800) , &
                          xray_form_type("cs  ", 55, (/ 20.389200,  19.106199,  10.662000,   1.495300/), &
                                                     (/  3.569000,   0.310700,  24.387899, 213.903992/),  3.335200) , &
                          xray_form_type("cs+1", 55, (/ 20.352400,  19.127800,  10.282100,   0.961500/), &
                                                     (/  3.552000,   0.308600,  23.712799,  59.456497/),  3.279100) , &
                          xray_form_type("ba  ", 56, (/ 20.336100,  19.297001,  10.888000,   2.695900/), &
                                                     (/  3.216000,   0.275600,  20.207300, 167.201996/),  2.773100) , &
                          xray_form_type("ba+2", 56, (/ 20.180700,  19.113600,  10.905399,   0.776340/), &
                                                     (/  3.213670,   0.283310,  20.055799,  51.745998/),  3.029020) /)

       Xray_form(121:130)=(/ &
                          xray_form_type("la  ", 57, (/ 20.577999,  19.598999,  11.372700,   3.287190/), &
                                                     (/  2.948170,   0.244475,  18.772600, 133.123993/),  2.146780) , &
                          xray_form_type("la+3", 57, (/ 20.248899,  19.376301,  11.632299,   0.336048/), &
                                                     (/  2.920700,   0.250698,  17.821100,  54.945297/),  2.408600) , &
                          xray_form_type("ce  ", 58, (/ 21.167099,  19.769501,  11.851299,   3.330490/), &
                                                     (/  2.812190,   0.226836,  17.608299, 127.112999/),  1.862640) , &
                          xray_form_type("ce+3", 58, (/ 20.803600,  19.559000,  11.936900,   0.612376/), &
                                                     (/  2.776910,   0.231540,  16.540800,  43.169201/),  2.090130) , &
                          xray_form_type("ce+4", 58, (/ 20.323500,  19.818600,  12.123300,   0.144583/), &
                                                     (/  2.659410,   0.218850,  15.799200,  62.235500/),  1.591800) , &
                          xray_form_type("pr  ", 59, (/ 22.043999,  19.669701,  12.385600,   2.824280/), &
                                                     (/  2.773930,   0.222087,  16.766899, 143.643997/),  2.058300) , &
                          xray_form_type("pr+3", 59, (/ 21.372700,  19.749100,  12.132900,   0.975180/), &
                                                     (/  2.645200,   0.214299,  15.323000,  36.406502/),  1.771320) , &
                          xray_form_type("pr+4", 59, (/ 20.941299,  20.053900,  12.466800,   0.296689/), &
                                                     (/  2.544670,   0.202481,  14.813700,  45.464298/),  1.242850) , &
                          xray_form_type("nd  ", 60, (/ 22.684500,  19.684700,  12.774000,   2.851370/), &
                                                     (/  2.662480,   0.210628,  15.885000, 137.903000/),  1.984860) , &
                          xray_form_type("nd+3", 60, (/ 21.961000,  19.933899,  12.120000,   1.510310/), &
                                                     (/  2.527220,   0.199237,  14.178300,  30.871700/),  1.475880) /)

       Xray_form(131:140)=(/ &
                          xray_form_type("pm  ", 61, (/ 23.340500,  19.609501,  13.123500,   2.875160/), &
                                                     (/  2.562700,   0.202088,  15.100900, 132.720993/),  2.028760) , &
                          xray_form_type("pm+3", 61, (/ 22.552700,  20.110800,  12.067100,   2.074920/), &
                                                     (/  2.417400,   0.185769,  13.127500,  27.449100/),  1.194990) , &
                          xray_form_type("sm  ", 62, (/ 24.004200,  19.425800,  13.439600,   2.896040/), &
                                                     (/  2.472740,   0.196451,  14.399600, 128.007004/),  2.209630) , &
                          xray_form_type("sm+3", 62, (/ 23.150400,  20.259899,  11.920200,   2.714880/), &
                                                     (/  2.316410,   0.174081,  12.157100,  24.824200/),  0.954586) , &
                          xray_form_type("eu  ", 63, (/ 24.627399,  19.088600,  13.760300,   2.922700/), &
                                                     (/  2.387900,   0.194200,  13.754600, 123.173996/),  2.574500) , &
                          xray_form_type("eu+2", 63, (/ 24.006300,  19.950399,  11.803400,   3.872430/), &
                                                     (/  2.277830,   0.173530,  11.609600,  26.515600/),  1.363890) , &
                          xray_form_type("eu+3", 63, (/ 23.749699,  20.374500,  11.850900,   3.265030/), &
                                                     (/  2.222580,   0.163940,  11.311000,  22.996599/),  0.759344) , &
                          xray_form_type("gd  ", 64, (/ 25.070900,  19.079800,  13.851800,   3.545450/), &
                                                     (/  2.253410,   0.181951,  12.933100, 101.397995/),  2.419600) , &
                          xray_form_type("gd+3", 64, (/ 24.346600,  20.420799,  11.870800,   3.714900/), &
                                                     (/  2.135530,   0.155525,  10.578199,  21.702900/),  0.645089) , &
                          xray_form_type("tb  ", 65, (/ 25.897600,  18.218500,  14.316700,   2.953540/), &
                                                     (/  2.242560,   0.196143,  12.664800, 115.362000/),  3.582240) /)

       Xray_form(141:150)=(/ &
                          xray_form_type("tb+3", 65, (/ 24.955900,  20.327099,  12.247100,   3.773000/), &
                                                     (/  2.056010,   0.149525,  10.049900,  21.277300/),  0.691967) , &
                          xray_form_type("dy  ", 66, (/ 26.507000,  17.638300,  14.559600,   2.965770/), &
                                                     (/  2.180200,   0.202172,  12.189899, 111.874001/),  4.297280) , &
                          xray_form_type("dy+3", 66, (/ 25.539499,  20.286100,  11.981200,   4.500730/), &
                                                     (/  1.980400,   0.143384,   9.349720,  19.580999/),  0.689690) , &
                          xray_form_type("ho  ", 67, (/ 26.904900,  17.293999,  14.558300,   3.638370/), &
                                                     (/  2.070510,   0.197940,  11.440700,  92.656601/),  4.567960) , &
                          xray_form_type("ho+3", 67, (/ 26.129601,  20.099400,  11.978800,   4.936760/), &
                                                     (/  1.910720,   0.139358,   8.800180,  18.590799/),  0.852795) , &
                          xray_form_type("er  ", 68, (/ 27.656300,  16.428499,  14.977900,   2.982330/), &
                                                     (/  2.073560,   0.223545,  11.360400, 105.703003/),  5.920460) , &
                          xray_form_type("er+3", 68, (/ 26.722000,  19.774799,  12.150600,   5.173790/), &
                                                     (/  1.846590,   0.137290,   8.362249,  17.897400/),  1.176130) , &
                          xray_form_type("tm  ", 69, (/ 28.181900,  15.885099,  15.154200,   2.987060/), &
                                                     (/  2.028590,   0.238849,  10.997499, 102.960999/),  6.756210) , &
                          xray_form_type("tm+3", 69, (/ 27.308300,  19.332001,  12.333900,   5.383480/), &
                                                     (/  1.787110,   0.136974,   7.967780,  17.292200/),  1.639290) , &
                          xray_form_type("yb  ", 70, (/ 28.664101,  15.434500,  15.308700,   2.989630/), &
                                                     (/  1.988900,   0.257119,  10.664700, 100.417000/),  7.566720) /)

       Xray_form(151:160)=(/ &
                          xray_form_type("yb+2", 70, (/ 28.120899,  17.681700,  13.333500,   5.146570/), &
                                                     (/  1.785030,   0.159970,   8.183040,  20.389999/),  3.709830) , &
                          xray_form_type("yb+3", 70, (/ 27.891701,  18.761400,  12.607200,   5.476470/), &
                                                     (/  1.732720,   0.138790,   7.644120,  16.815300/),  2.260010) , &
                          xray_form_type("lu  ", 71, (/ 28.947599,  15.220800,  15.100000,   3.716010/), &
                                                     (/  1.901820,   9.985189,   0.261033,  84.329803/),  7.976280) , &
                          xray_form_type("lu+3", 71, (/ 28.462799,  18.121000,  12.842899,   5.594150/), &
                                                     (/  1.682160,   0.142292,   7.337270,  16.353500/),  2.975730) , &
                          xray_form_type("hf  ", 72, (/ 29.143999,  15.172600,  14.758600,   4.300130/), &
                                                     (/  1.832620,   9.599899,   0.275116,  72.028999/),  8.581540) , &
                          xray_form_type("hf+4", 72, (/ 28.813099,  18.460100,  12.728500,   5.599270/), &
                                                     (/  1.591360,   0.128903,   6.762320,  14.036600/),  2.396990) , &
                          xray_form_type("ta  ", 73, (/ 29.202400,  15.229300,  14.513500,   4.764920/), &
                                                     (/  1.773330,   9.370460,   0.295977,  63.364399/),  9.243540) , &
                          xray_form_type("ta+5", 73, (/ 29.158699,  18.840700,  12.826799,   5.386950/), &
                                                     (/  1.507110,   0.116741,   6.315240,  12.424400/),  1.785550) , &
                          xray_form_type("w   ", 74, (/ 29.081800,  15.430000,  14.432700,   5.119820/), &
                                                     (/  1.720290,   9.225900,   0.321703,  57.056000/),  9.887500) , &
                          xray_form_type("w+6 ", 74, (/ 29.493599,  19.376301,  13.054399,   5.064120/), &
                                                     (/  1.427550,   0.104621,   5.936670,  11.197200/),  1.010740) /)

       Xray_form(161:170)=(/ &
                          xray_form_type("re  ", 75, (/ 28.762100,  15.718900,  14.556400,   5.441740/), &
                                                     (/  1.671910,   9.092270,   0.350500,  52.086098/), 10.472000) , &
                          xray_form_type("os  ", 76, (/ 28.189400,  16.154999,  14.930500,   5.675890/), &
                                                     (/  1.629030,   8.979480,   0.382661,  48.164700/), 11.000500) , &
                          xray_form_type("os+4", 76, (/ 30.418999,  15.263700,  14.745800,   5.067950/), &
                                                     (/  1.371130,   6.847060,   0.165191,  18.003000/),  6.498040) , &
                          xray_form_type("ir  ", 77, (/ 27.304899,  16.729599,  15.611500,   5.833770/), &
                                                     (/  1.592790,   8.865530,   0.417916,  45.001099/), 11.472200) , &
                          xray_form_type("ir+3", 77, (/ 30.415600,  15.862000,  13.614500,   5.820080/), &
                                                     (/  1.343230,   7.109090,   0.204633,  20.325399/),  8.279030) , &
                          xray_form_type("ir+4", 77, (/ 30.705799,  15.551200,  14.232600,   5.536720/), &
                                                     (/  1.309230,   6.719830,   0.167252,  17.491100/),  6.968240) , &
                          xray_form_type("pt  ", 78, (/ 27.005899,  17.763901,  15.713100,   5.783700/), &
                                                     (/  1.512930,   8.811740,   0.424593,  38.610298/), 11.688300) , &
                          xray_form_type("pt+2", 78, (/ 29.842899,  16.722401,  13.215300,   6.352340/), &
                                                     (/  1.329270,   7.389790,   0.263297,  22.942600/),  9.853290) , &
                          xray_form_type("pt+4", 78, (/ 30.961201,  15.982900,  13.734800,   5.920340/), &
                                                     (/  1.248130,   6.608340,   0.168640,  16.939199/),  7.395340) , &
                          xray_form_type("au  ", 79, (/ 16.881901,  18.591299,  25.558201,   5.860000/), &
                                                     (/  0.461100,   8.621600,   1.482600,  36.395599/), 12.065800) /)

       Xray_form(171:180)=(/ &
                          xray_form_type("au+1", 79, (/ 28.010899,  17.820400,  14.335899,   6.580770/), &
                                                     (/  1.353210,   7.739500,   0.356752,  26.404301/), 11.229900) , &
                          xray_form_type("au+3", 79, (/ 30.688599,  16.902901,  12.780100,   6.523540/), &
                                                     (/  1.219900,   6.828720,   0.212867,  18.659000/),  9.096800) , &
                          xray_form_type("hg  ", 80, (/ 20.680901,  19.041700,  21.657499,   5.967600/), &
                                                     (/  0.545000,   8.448400,   1.572900,  38.324600/), 12.608900) , &
                          xray_form_type("hg+1", 80, (/ 25.085300,  18.497299,  16.888300,   6.482160/), &
                                                     (/  1.395070,   7.651050,   0.443378,  28.226200/), 12.020500) , &
                          xray_form_type("hg+2", 80, (/ 29.564100,  18.059999,  12.837400,   6.899120/), &
                                                     (/  1.211520,   7.056390,   0.284738,  20.748199/), 10.626800) , &
                          xray_form_type("tl  ", 81, (/ 27.544600,  19.158400,  15.538000,   5.525930/), &
                                                     (/  0.655150,   8.707510,   1.963470,  45.814899/), 13.174600) , &
                          xray_form_type("tl+1", 81, (/ 21.398500,  20.472300,  18.747799,   6.828470/), &
                                                     (/  1.471100,   0.517394,   7.434630,  28.848200/), 12.525800) , &
                          xray_form_type("tl+3", 81, (/ 30.869499,  18.384100,  11.932800,   7.005740/), &
                                                     (/  1.100800,   6.538520,   0.219074,  17.211399/),  9.802700) , &
                          xray_form_type("pb  ", 82, (/ 31.061699,  13.063700,  18.441999,   5.969600/), &
                                                     (/  0.690200,   2.357600,   8.618000,  47.257900/), 13.411800) , &
                          xray_form_type("pb+2", 82, (/ 21.788601,  19.568199,  19.140600,   7.011070/), &
                                                     (/  1.336600,   0.488383,   6.772700,  23.813200/), 12.473400) /)

       Xray_form(181:190)=(/ &
                          xray_form_type("pb+4", 82, (/ 32.124397,  18.800301,  12.017500,   6.968860/), &
                                                     (/  1.005660,   6.109260,   0.147041,  14.714000/),  8.084280) , &
                          xray_form_type("bi  ", 83, (/ 33.368900,  12.951000,  16.587700,   6.469200/), &
                                                     (/  0.704000,   2.923800,   8.793700,  48.009300/), 13.578199) , &
                          xray_form_type("bi+3", 83, (/ 21.805300,  19.502600,  19.105301,   7.102950/), &
                                                     (/  1.235600,   6.241490,   0.469999,  20.318501/), 12.471100) , &
                          xray_form_type("bi+5", 83, (/ 33.536400,  25.094601,  19.249699,   6.915550/), &
                                                     (/  0.916540,   0.390420,   5.714140,  12.828500/), -6.799400) , &
                          xray_form_type("po  ", 84, (/ 34.672600,  15.473300,  13.113800,   7.025880/), &
                                                     (/  0.700999,   3.550780,   9.556419,  47.004501/), 13.677000) , &
                          xray_form_type("at  ", 85, (/ 35.316299,  19.021099,   9.498870,   7.425180/), &
                                                     (/  0.685870,   3.974580,  11.382400,  45.471500/), 13.710800) , &
                          xray_form_type("rn  ", 86, (/ 35.563099,  21.281601,   8.003700,   7.443300/), &
                                                     (/  0.663100,   4.069100,  14.042200,  44.247299/), 13.690500) , &
                          xray_form_type("fr  ", 87, (/ 35.929901,  23.054699,  12.143900,   2.112530/), &
                                                     (/  0.646453,   4.176190,  23.105200, 150.644989/), 13.724700) , &
                          xray_form_type("ra  ", 88, (/ 35.763000,  22.906399,  12.473900,   3.210970/), &
                                                     (/  0.616341,   3.871350,  19.988701, 142.324997/), 13.621099) , &
                          xray_form_type("ra+2", 88, (/ 35.215000,  21.670000,   7.913420,   7.650780/), &
                                                     (/  0.604909,   3.576700,  12.601000,  29.843599/), 13.543100) /)

       Xray_form(191:200)=(/ &
                          xray_form_type("ac  ", 89, (/ 35.659698,  23.103199,  12.597700,   4.086550/), &
                                                     (/  0.589092,   3.651550,  18.598999, 117.019997/), 13.526600) , &
                          xray_form_type("ac+3", 89, (/ 35.173599,  22.111200,   8.192160,   7.055450/), &
                                                     (/  0.579689,   3.414370,  12.918700,  25.944300/), 13.463699) , &
                          xray_form_type("th  ", 90, (/ 35.564499,  23.421900,  12.747300,   4.807030/), &
                                                     (/  0.563359,   3.462040,  17.830900,  99.172195/), 13.431400) , &
                          xray_form_type("th+4", 90, (/ 35.100700,  22.441799,   9.785540,   5.294440/), &
                                                     (/  0.555054,   3.244980,  13.466100,  23.953300/), 13.375999) , &
                          xray_form_type("pa  ", 91, (/ 35.884701,  23.294800,  14.189100,   4.172870/), &
                                                     (/  0.547751,   3.415190,  16.923500, 105.250999/), 13.428699) , &
                          xray_form_type("u   ", 92, (/ 36.022800,  23.412800,  14.949100,   4.188000/), &
                                                     (/  0.529300,   3.325300,  16.092699, 100.612999/), 13.396600) , &
                          xray_form_type("u+3 ", 92, (/ 35.574699,  22.525900,  12.216499,   5.370730/), &
                                                     (/  0.520480,   3.122930,  12.714800,  26.339399/), 13.309200) , &
                          xray_form_type("u+4 ", 92, (/ 35.371498,  22.532600,  12.029100,   4.798400/), &
                                                    (/  0.516598,   3.050530,  12.572300,  23.458200/), 13.267099) , &
                          xray_form_type("u+6 ", 92, (/ 34.850899,  22.758400,  14.009900,   1.214570/), &
                                                     (/  0.507079,   2.890300,  13.176700,  25.201700/), 13.166500) , &
                          xray_form_type("np  ", 93, (/ 36.187401,  23.596399,  15.640200,   4.185500/), &
                                                     (/  0.511929,   3.253960,  15.362200,  97.490799/), 13.357300) /)

       Xray_form(201:210)=(/ &
                          xray_form_type("np+3", 93, (/ 35.707397,  22.612999,  12.989799,   5.432270/), &
                                                     (/  0.502322,   3.038070,  12.144899,  25.492800/), 13.254400) , &
                          xray_form_type("np+4", 93, (/ 35.510300,  22.578699,  12.776600,   4.921590/), &
                                                     (/  0.498626,   2.966270,  11.948400,  22.750200/), 13.211599) , &
                          xray_form_type("np+6", 93, (/ 35.013599,  22.728600,  14.388400,   1.756690/), &
                                                     (/  0.489810,   2.810990,  12.330000,  22.658100/), 13.113000) , &
                          xray_form_type("pu  ", 94, (/ 36.525398,  23.808300,  16.770700,   3.479470/), &
                                                     (/  0.499384,   3.263710,  14.945499, 105.979996/), 13.381200) , &
                          xray_form_type("pu+3", 94, (/ 35.840000,  22.716900,  13.580700,   5.660160/), &
                                                    (/  0.484936,   2.961180,  11.533100,  24.399200/), 13.199100) , &
                          xray_form_type("pu+4", 94, (/ 35.649300,  22.646000,  13.359500,   5.188310/), &
                                                     (/  0.481422,   2.890200,  11.316000,  21.830099/), 13.155500) , &
                          xray_form_type("pu+6", 94, (/ 35.173599,  22.718100,  14.763500,   2.286780/), &
                                                     (/  0.473204,   2.738480,  11.552999,  20.930300/), 13.058200) , &
                          xray_form_type("am  ", 95, (/ 36.670601,  24.099199,  17.341499,   3.493310/), &
                                                     (/  0.483629,   3.206470,  14.313600, 102.272995/), 13.359200) , &
                          xray_form_type("cm  ", 96, (/ 36.648800,  24.409599,  17.399000,   4.216650/), &
                                                     (/  0.465154,   3.089970,  13.434600,  88.483398/), 13.288700) , &
                          xray_form_type("bk  ", 97, (/ 36.788101,  24.773600,  17.891899,   4.232840/), &
                                                     (/  0.451018,   3.046190,  12.894600,  86.002998/), 13.275400) /)

       Xray_form(211:214)=(/ &
                          xray_form_type("cf  ", 98, (/ 36.918499,  25.199499,  18.331699,   4.243910/), &
                                                     (/  0.437533,   3.007750,  12.404400,  83.788101/), 13.267400) , &
                          xray_form_type("o-2 ",  8, (/  4.758000,   3.637000,   0.000000,   0.000000/), &
                                                     (/  7.831000,  30.049999,   0.000000,   0.000000/), 1.5940000) , &
                          xray_form_type("ze  ",  1, (/  0.000000,   0.000000,   0.000000,   0.000000/), &
                                                     (/  0.000000,   0.000000,   0.000000,   0.000000/), 0.0000000) , &
                          xray_form_type("d   ",  1, (/  0.493002,   0.322912,   0.140191,   0.040810/), &
                                                     (/ 10.510900,  26.125700,   3.142360,  57.799698/),  0.003038) /)
       return
    End  Subroutine Set_Xray_Form

 End Module CFML_Scattering_Chemical_Tables

