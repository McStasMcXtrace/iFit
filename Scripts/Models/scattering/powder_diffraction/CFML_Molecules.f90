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
!!---- MODULE: CFML_Molecular_Crystals
!!----   INFO: Module to define molecules on Crystals
!!----
!!---- HISTORY
!!----    Update: 07/03/2011
!!----
!!---- COMMENTARY
!!--..    Explanations about Eulerian angles, active and passive rotations
!!--..    Ref. Texture Analysis in Material Science, H.J. Bunge.
!!--..    Ed Butterworths London 1970?
!!--..
!!--..
!!--..    First variant:
!!--..    -------------
!!--..    Eulerian angles g={phi1,PHI,phi2}, positive rotations: anti-clockwise
!!--..
!!--..     1: Rotation around the common Z,Zm axis of an angle phi1
!!--..     2: Rotation around the new Xm axis of an angle PHI
!!--..     3: Rotation around the new Zm-axis of an angle phi2
!!--..
!!--..          g = gZm(phi2) . gXm(PHI) . gZm(Zm)
!!--..
!!--..
!!--..                   (  cosphi2  sinphi2    0  )                     (   1      0       0      )
!!--..
!!--..       gZm(phi2) =(  -sinphi2  cosphi2    0   )        gXm(PHI)  =(    0    cosPHI  sinPHI    )
!!--..
!!--..                   (    0         0       1  )                     (   0   -sinPHI  cosPHI   )
!!--..
!!--..
!!--..
!!--..                   (  cosphi1  sinphi1    0  )
!!--..
!!--..       gZm(phi1) =(  -sinphi1  cosphi1    0   )
!!--..
!!--..                   (    0         0       1  )
!!--..
!!--..
!!--..
!!--..    Second variant:
!!--..    ---------------
!!--..     Eulerian angles g={PSI,THETA,PHI}, positive rotations: anti-clockwise
!!--..
!!--..     1: Rotation around the common Z,Zm axis of an angle PSI   (PHI)
!!--..     2: Rotation around the new Ym axis of an angle THETA      (THETA) <--- FullProf
!!--..     3: Rotation around the new Zm-axis of an angle PHI        (CHI)
!!--..
!!--..       phi1=PSI+pi/2   PHI = THETA   phi2=PHI - pi/2
!!--..
!!--..
!!--..
!!--..    Rotation Axis and Rotation Angle
!!--..    --------------------------------
!!--..
!!--..    The rotation axis is given by the unit vector u represented by its polar
!!--..    coordinates (theta,phi) and the rotation angle (omega) around u, so that
!!--..    one can write the rotation g={u,omega}={theta,phi,omega}
!!--..
!!--..    Passive rotations: one looks for the coordinates of a unique point respect
!!--..                       to two rotated frames
!!--..
!!--..                             ( cosphi  sinphi)
!!--..         -----\--------  M = (               )     is the matrix relating the basis (e)=M(i)
!!--..         |\    \             (-sinphi  cosphi)
!!--..         | \    \
!!--..         |  \    r        The point r  has coordinates (x,y) and coordinates (x',y')
!!--..         |   \ /          w.r.t. the rotated axes the relation is:
!!--..         |    \
!!--..         | phi \          (x')   ( cosphi  sinphi)  (x)     x'= x cosphi + y sinphi
!!--..                          (  ) = (               )  ( )  => y'=-x sinphi + y cosphi
!!--..                          (y')   (-sinphi  cosphi)  (y)
!!--..
!!--..     Active  rotations: one looks for the new coordinates of a point respect
!!--..                       to the same frame when a rotation is applied
!!--..
!!--..
!!--..                             ( cosphi -sinphi)
!!--..         --------------  R = (               )
!!--..         |\                  ( sinphi  cosphi)
!!--..         | \
!!--..         |  \   r'           x'= x cosphi - y sinphi
!!--..         | r \               y'= x sinphi + y cosphi
!!--..         |    \
!!--..         | phi \
!!--..
!!--..
!!--..      The representative matrices are one the inverse of the other. R=Minv=Mt
!!--..
!!--..      In molecular crystals one looks for determining the position of each atom of
!!--..      the molecule in the crystallographic frame when one knows the internal coordinates
!!--..      of the atoms, the position of the origin of the internal frame in the crystallographic
!!--..      frame and the orientation (Euler or Euler-like angles) of the internal frame with
!!--..      respect to the crystallographic frame.
!!--..
!!--..      The problem is to define a simple set of orientational angles
!!--..
!!--..      We shall adopt the conventional definition of Euler angles but we will call then
!!--..      a=phi1, b=PHI, c=phi2. The above matrices correspond to passive rotations, so that
!!--..      when applied to a fixed point their product will give the coordinates of this point
!!--..      with respect to the rotated system. In our case will give the position of an external
!!--..      point (cartesian crystal  frame, CCF) w.r.t the cartesian molecular frame (CMF).
!!--..      Taking the transpose of the final matrix one obtains an active rotation matrix that
!!--..      applied to a point moves it to a new point referred to the CCF.
!!--..
!!--..
!!--..
!!--..               ( cosa cosc - sina sinc cosb     sina cosc + cosa sinc cosb    sinc sinb )
!!--..
!!--..    g(a,b,c) =( -cosa sinc - sina cosc cosb    -sina sinc + cosa cosc cosb    cosc sinb  )
!!--..
!!--..               (         sina  sinb                   -cosa sinb                 cosb   )
!!--..
!!--..
!!--..               ( cosa cosc - sina sinc cosb    -cosa sinc - sina cosc cosb    sina sinb )
!!--..
!!--..   gt(a,b,c) =(  sina cosc + cosa sinc cosb    -sina sinc + cosa cosc cosb   -cosa sinb  )
!!--..
!!--..               (      sinc sinb                      cosc sinb                   cosb   )
!!--..
!!--..
!!--..
!!--..     The matrix g applied to a point with coordinates given w.r.t. CCF, provides the coordinates
!!--..     w.r.t. CMF. If we take a point in the CMF and we apply the matrix gt we obtain the coordinates
!!--..     of this point w.r.t. CCF.
!!--..
!!--..     Orientational angles used in FullProf
!!--..     -------------------------------------
!!--..
!!--..     The molecular frame (CMF) is supposed to coincide at the begining with the Cartesian
!!--..     crystallographic frame (CCF). To position a molecule in an arbitrary position the
!!--..     total movement is decomposed in the following way:
!!--..
!!--..     1) Perform a rotation of angle CHI around the Z,Zm-axis : the rotation matrix relating
!!--..        the two unitary bases (Em and E in form of columns) is the following:
!!--..
!!--..                       (cosCHI    sinCHI   0 )            (e1)           (i)
!!--..                                                          (  )           ( )
!!--..             Rz(CHI) =(-sinCHI    cosCHI   0  )        Em=(e2) = Rz(CHI) (j) = Rz(CHI) E
!!--..                                                          (  )           ( )
!!--..                       (  0         0      1 )            (e3)           (k)
!!--..
!!--..        An active rotation is obtained transposing Rz(CHI)t = Az(CHI). This matrix is
!!--..        applied to a point in CCF and provides the new coordinates in the CCF after
!!--..        the rotation of angle CHI around Z.
!!--..
!!--..
!!--..     2) Perform a rotation of angle THE around the Y-axis : the rotation matrix relating
!!--..        the two unitary bases (Em and E in form of columns) is now the following:
!!--..
!!--..                       (cosTHE   0   -sinTHE )            (e1)           (i)
!!--..                                                          (  )           ( )
!!--..             Ry(THE) =(   0      1      0     )        Em=(e2) = Ry(THE) (j) = Ry(THE) E
!!--..                                                          (  )           ( )
!!--..                       (sinTHE   0    cosTHE )            (e3)           (k)
!!--..
!!--..        An active rotation is obtained transposing Ry(THE)t = Ay(THE). This matrix is
!!--..        applied to a point in CCF and provides the new coordinates in the CCF after
!!--..        the rotation of angle THE around Y.
!!--..
!!--..     3) Perform a rotation of angle PHI around the Z-axis : the rotation matrix relating
!!--..        the two unitary bases (Em and E in form of columns) is the following:
!!--..
!!--..                       (cosPHI    sinPHI   0 )            (e1)           (i)
!!--..                                                          (  )           ( )
!!--..             Rz(PHI) =(-sinPHI    cosPHI   0  )        Em=(e2) = Rz(PHI) (j) = Rz(PHI) E
!!--..                                                          (  )           ( )
!!--..                       (  0         0      1 )            (e3)           (k)
!!--..
!!--..        An active rotation is obtained transposing Rz(PHI)t = Az(PHI). This matrix is
!!--..        applied to a point in CCF and provides the new coordinates in the CCF after
!!--..        the rotation of angle PHI around Z.
!!--..
!!--..    With this rotational angles the interpretation of the angles (THE,PHI) correspond to
!!--..    the spherical angles of the CMF Zm-axis with respect to the CCF. The total active
!!--..    matrix to be applied to atoms of the molecule in the initial position (when the two
!!--..    frames coincide) to get the final coordinates is the following:
!!--..
!!--..
!!--..                 M = Az(PHI) . Ay(THE) . Az(CHI) =  XA(PHI,THE)  . XAp(CHI)
!!--..
!!--..
!!--..     In the initial state the Cartesian coordinates of atoms (x), in columns, are
!!--..     the same in both frames, the positions after the total rotation are given by:
!!--..
!!--..                         (x)-final =   M (x)
!!--..
!!--..     To obtain the internal coordinates of a point in the CCF one must apply the
!!--..     following formula:
!!--..
!!--..                   X-internal  =  Mt  X = XAp(CHI)t . XA(PHI,THE)t  X
!!--..
!!--..     the final expressions of the different matrices are the following:
!!--..
!!--..
!!--..                    (cosPHI cosTHE      -sinPHI      cosPHI sinTHE )
!!--..
!!--..     XA(PHI,THE) = ( sinPHI cosTHE       cosPHI      sinPHI sinTHE  )
!!--..
!!--..                    (  -sinTHE             0             cosTHE    )
!!--..
!!--..
!!--..                   (cosCHI   -sinCHI   0 )
!!--..
!!--..        XAp(CHI) =( sinCHI    cosCHI   0  )
!!--..
!!--..                   (  0         0      1 )
!!--..
!!--..
!!--..               ( cosa cosc - sina sinc cosb     sina cosc + cosa sinc cosb    sinc sinb )
!!--..
!!--..    g(a,b,c) =( -cosa sinc - sina cosc cosb    -sina sinc + cosa cosc cosb    cosc sinb  )
!!--..
!!--..               (         sina  sinb                   -cosa sinb                 cosb   )
!!--..
!!--..
!!--..               ( cosa cosc - sina sinc cosb    -cosa sinc - sina cosc cosb    sina sinb )
!!--..
!!--..   gt(a,b,c) =(  sina cosc + cosa sinc cosb    -sina sinc + cosa cosc cosb   -cosa sinb  )
!!--..
!!--..               (      sinc sinb                      cosc sinb                   cosb   )
!!--..
!!--..
!!--..
!!--..
!!--..
!!--..
!!--..  M(PHI,THE,CHI) =
!!--..
!!--..     (cosPHI cosTHE cosCHI - sinPHI sinCHI   -cosPHI cosTHE sinCHI - sinPHI cosCHI    cosPHI sinTHE)
!!--..
!!--..   =( sinPHI cosTHE cosCHI + cosPHI sinCHI   -sinPHI cosTHE sinCHI + cosPHI cosCHI    sinPHI sinTHE )
!!--..
!!--..     (       -sinTHE cosCHI                           sinTHE sinCHI                      cosTHE    )
!!--..
!!--..
!!--..
!!--..   Comparing the matrix M(THE,PHI,CHI) with the matrix gt(a,b,c)=gt(alpha,beta,gamma)=gt(phi1,PHI,phi2)
!!--..
!!--..   One can see that both matrices are identical if we take:
!!--..
!!--..        alpha=phi1=PHI+pi/2     beta=PHI=THETA   gamma=phi2=CHI-pi/2
!!--..
!!--..         (phi1=PSI+pi/2   PHI = THETA   phi2=PHI - pi/2)
!!--..
!!--..
!!--..      The angles used in FullProf correspond to the second variant of Euler angles making the
!!--..      sustitution:
!!--..
!!--..        (PSI,THETA,PHI)  --->   (PHI, THETA, CHI)
!!--..
!!--..            2nd variant   -->      FullProf
!!--..
!!--..      This is clear from the following. If we take passive rotations as for deriving the matrix
!!--..  corresponding to the Euler angles the matrix Mt should be the result
!!--..
!!--..      Mt = (  Az(PHI) . Ay(THE) . Az(CHI) )t = Rz(CHI) . Ry(THE) . Rz(PHI)
!!--..
!!--..  Then the interpretation of the rotations are strictly the same as given in the description
!!--..  of the second variant of Euler angles.
!!--..
!!--..
!!----
!!---- DEPENDENCIES
!!----
!!---- VARIABLES
!!----    ERR_MOLEC
!!----    ERR_MOLEC_MESS
!!----    MOLECULE_TYPE
!!----    MOLECULAR_CRYSTAL_TYPE
!!----
!!---- PROCEDURES
!!----    Functions:
!!----
!!----    Subroutines:
!!----       CARTESIAN_TO_FRACTIONAL
!!----       CARTESIAN_TO_SPHERICAL
!!----       CARTESIAN_TO_ZMATRIX
!!--++       CREATE_CONNECTIVITY_CARTESIAN   [Private]
!!----       EMPIRIC_FORMULA
!!--++       EMPIRIC_FORMULA_FATOM           [Overloaded]
!!--++       EMPIRIC_FORMULA_MOLCRYS         [Overloaded]
!!--++       EMPIRIC_FORMULA_MOLEC           [Overloaded]
!!----       FIX_REFERENCE
!!----       FIX_ORIENT_CARTESIAN
!!----       FRACTIONAL_TO_CARTESIAN
!!----       FRACTIONAL_TO_SPHERICAL
!!----       FRACTIONAL_TO_ZMATRIX
!!--++       GET_CARTESIAN_FROM_Z            [Private]
!!--++       GET_Z_FROM_CARTESIAN            [Private]
!!----       INIT_ERR_MOLEC
!!----       INIT_MOLECULE
!!----       MOLCRYS_TO_ATOMLIST
!!----       MOLEC_TO_ATOMLIST
!!----       READ_FREE_ATOMS
!!----       READ_MOLECULE
!!--++       READ_MOLECULE_IN_FILE           [Overloaded]
!!--++       READ_MOLECULE_IN_VAR            [Overloaded]
!!----       SET_EULER_MATRIX
!!----       SPHERICAL_TO_CARTESIAN
!!----       SPHERICAL_TO_FRACTIONAL
!!----       SPHERICAL_TO_ZMATRIX
!!----       WRITE_FREE_ATOMS
!!----       WRITE_MOLECULAR_CRYSTAL
!!----       WRITE_MOLECULE
!!----       ZMATRIX_TO_CARTESIAN
!!----       ZMATRIX_TO_FRACTIONAL
!!----       ZMATRIX_TO_SPHERICAL
!!----
!!
 Module CFML_Molecular_Crystals

    !---- Use Modules ----!
    use CFML_GlobalDeps,                only: cp, eps, to_rad
    use CFML_Math_General,              only: acosd, asind, cosd, sind
    use CFML_Math_3D,                   only: cross_product, Get_Spheric_Coord
    use CFML_Crystallographic_Symmetry, only: Space_Group_type, Write_SpaceGroup
    use CFML_Atom_TypeDef,              only: Atom_Type, Atom_List_Type, Allocate_Atom_List, Deallocate_Atom_List,&
                                              Init_Atom_Type
    use CFML_Crystal_Metrics,           only: Crystal_Cell_Type, Set_Crystal_Cell,Err_crys, Err_Crys_Mess, &
                                              Write_Crystal_Cell
    use CFML_String_Utilities,          only: u_case, l_case, getword, getnum, cutst
    use CFML_Geometry_Calc,             only: angle_dihedral,distance,Get_PhiTheChi
    use CFML_Scattering_Chemical_Tables,only: Num_Chem_Info,Chem_Info,Set_Chem_Info,Remove_Chem_Info,Get_ChemSymb

    implicit none

    private

    !---- List of public functions ----!

    !---- List of public overloaded procedures: functions ----!

    !---- List of public subroutines ----!
    public :: Init_Err_Molec, Init_Molecule, Read_Free_Atoms, Read_Molecule,             &
              Write_Molecule, Write_Molecular_Crystal, Write_Free_Atoms

    public :: Cartesian_to_Fractional, Cartesian_to_Spherical, Cartesian_to_Zmatrix,     &
              Fractional_to_Cartesian, Fractional_to_Spherical, Fractional_to_Zmatrix,   &
              Zmatrix_to_Cartesian, Zmatrix_to_Fractional, Zmatrix_to_Spherical,         &
              Spherical_to_Cartesian, Spherical_to_Zmatrix,Spherical_to_Fractional,      &
              Fix_Reference,Fix_Orient_Cartesian, Set_Euler_Matrix, Molcrys_to_AtomList, &
              Molec_to_AtomList, Empiric_Formula,Init_Mol_Crys

    !---- List of private functions ----!

    !---- List of private Subroutines ----!
    private :: Create_Connectivity_Cartesian, Get_Cartesian_From_Z, Get_Z_From_Cartesian, &
               Empiric_Formula_FAtom, Empiric_Formula_Molcrys, Empiric_Formula_Molec


    !---- Definitions ----!

    !!----
    !!---- ERR_MOLEC
    !!----    logical, public :: err_molec
    !!----
    !!----    Logical Variable indicating an error in MOLECULAR_CRYSTAL module
    !!----
    !!---- Update: February - 2005
    !!
    logical, public          :: Err_Molec

    !!----
    !!---- ERR_MOLEC_MESS
    !!----    character(len=150), public :: ERR_Molec_Mess
    !!----
    !!----    String containing information about the last error
    !!----
    !!---- Update: February - 2005
    !!
    character(len=150), public :: ERR_Molec_Mess

    !!----
    !!----  TYPE :: MOLECULE_TYPE
    !!--..
    !!----  Type, public :: Molecule_Type
    !!----     character(len=80)                               :: Name_mol     !Global name for the molecule
    !!----     integer                                         :: natoms       !Number of atoms
    !!----     logical                                         :: in_xtal      !True if global coordinates xcentre, orient are defined
    !!----     logical                                         :: is_EulerMat  !True if the Euler Matrix has been set
    !!----     logical                                         :: is_connect   !True if the connectivity is correct
    !!----     character(len=1)                                :: rot_type     !Type of rotational angles
    !!----                                                                     !"E": Conventional Euler angles (alpha,beta,gamma)
    !!----                                                                     !"P": Second variant of Euler angles (default)
    !!----                                                                     !     Polar:(theta,phi,chi)
    !!----     character(len=1)                                :: coor_type    !Type of internal coordinates
    !!----                                                                     !"Z": Z-matrix
    !!----                                                                     !"C": Cartesian
    !!----                                                                     !"S": Spherical
    !!----                                                                     !"F": Fractional coordinates (only if in_xtal = .true.)
    !!----     character(len=3)                                :: therm_type   !Type of thermal factor
    !!----                                                                     !"ISO": No collective motion
    !!----                                                                     !"T  ": Translational
    !!----                                                                     !"TL ": Translational + Librational
    !!----                                                                     !"TLS": Translational + Librational + Correlation
    !!----     real(kind=cp), dimension(3)                     :: xcentre      !Fractional coordinates of the centre
    !!----     real(kind=cp), dimension(3)                     :: mxcentre     !Refinement codes (or multipliers) of Fractional coordinates of the centre
    !!----     integer,       dimension(3)                     :: lxcentre     !Numbers of LSQ parameters for Fractional coordinates of the centre
    !!----     real(kind=cp), dimension(3)                     :: Orient       !Orientation angles (Euler angles or variant ...)
    !!----     real(kind=cp), dimension(3)                     :: mOrient      !Refinement codes (or multipliers) of Orientation angles (Euler angles or variant ...)
    !!----     integer,       dimension(3)                     :: lOrient      !Numbers of LSQ parameters for Orientation angles (Euler angles or variant ...)
    !!----     real(kind=cp), dimension(6)                     :: T_TLS        !Translational Thermal factor tensor
    !!----     real(kind=cp), dimension(6)                     :: mT_TLS       !Refinement codes (or multipliers) of Translational Thermal factor tensor
    !!----     integer,       dimension(6)                     :: lT_TLS       !Numbers of LSQ parameters for Translational Thermal factor tensor
    !!----     real(kind=cp), dimension(6)                     :: L_TLS        !Librational Thermal factor tensor
    !!----     real(kind=cp), dimension(6)                     :: mL_TLS       !Refinement codes (or multipliers) of Librational Thermal factor tensor
    !!----     integer,       dimension(6)                     :: lL_TLS       !Numbers of LSQ parameters for Librational Thermal factor tensor
    !!----     real(kind=cp), dimension(3,3)                   :: S_TLS        !TL-correlation Thermal factor
    !!----     real(kind=cp), dimension(3,3)                   :: mS_TLS       !Refinement codes (or multipliers) of TL-correlation Thermal factor
    !!----     integer,       dimension(3,3)                   :: lS_TLS       !Numbers of LSQ parameters for TL-correlation Thermal factor
    !!----     real(kind=cp), dimension(3,3)                   :: Euler        !Euler matrix
    !!----     character(len=6),  allocatable, dimension(  :)  :: AtName       !Atom Name
    !!----     character(len=4),  allocatable, dimension(  :)  :: AtSymb       !Atom species
    !!----     integer,           allocatable, dimension(  :)  :: AtZ          !Atomic Number
    !!----     integer,           allocatable, dimension(:,:)  :: Ptr          !Pointer to scat.factors (first index -> pattern)
    !!----     real(kind=cp),     allocatable, dimension(:,:)  :: I_coor       !Internal coordinates (d,ang,dang)
    !!----     real(kind=cp),     allocatable, dimension(:,:)  :: mI_coor      !Refinement codes (or multipliers) of internal coordinates
    !!----     integer,           allocatable, dimension(:,:)  :: lI_coor      !Numbers of LSQ parameters for internal coordinates
    !!----     real(kind=cp),     allocatable, dimension(  :)  :: biso         !Isotropic temperature factor
    !!----     real(kind=cp),     allocatable, dimension(  :)  :: mbiso        !Refinement codes (or multipliers) of Isotropic temperature factor
    !!----     integer,           allocatable, dimension(  :)  :: lbiso        !Numbers of LSQ parameters for Isotropic temperature factor
    !!----     real(kind=cp),     allocatable, dimension(  :)  :: occ          !Occupation factor
    !!----     real(kind=cp),     allocatable, dimension(  :)  :: mocc         !Refinement codes (or multipliers) of Occupation factor
    !!----     integer,           allocatable, dimension(  :)  :: locc         !Numbers of LSQ parameters for Occupation factor
    !!----     integer,           allocatable, dimension(  :)  :: Nb           !Number of neighbours
    !!----     integer,           allocatable, dimension(:,:)  :: inb          !Index of neighbours
    !!----     integer,           allocatable, dimension(:,:)  :: Tb           !Type of bonds
    !!----     integer,           allocatable, dimension(:,:)  :: conn         !Conectivity (N1,N2,N3)
    !!----  End Type Molecule_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Molecule_type
       character(len=80)                               :: Name_mol     !Global name for the molecule
       integer                                         :: natoms       !Number of atoms
       logical                                         :: in_xtal      !True if global coordinates xcentre, orient are defined
       logical                                         :: is_EulerMat  !True if the Euler Matrix has been set
       logical                                         :: is_connect   !True if the connectivity is given and correct
       character(len=1)                                :: rot_type     !Type of rotational angles
                                                                       !"E": Conventional Euler angles (alpha,beta,gamma)
                                                                       !"P": Second variant of Euler angles (default)
                                                                       !     Polar:(theta,phi,chi)
       character(len=1)                                :: coor_type    !Type of internal coordinates
                                                                       !"Z": Z-matrix
                                                                       !"C": Cartesian
                                                                       !"S": Spherical
                                                                       !"F": Fractional coordinates (only if in_xtal = .true.)
       character(len=3)                                :: therm_type   !Type of thermal factor
                                                                       !"ISO": No collective motion
                                                                       !"T  ": Translational
                                                                       !"TL ": Translational + Librational
                                                                       !"TLS": Translational + Librational + Correlation
       real(kind=cp), dimension(3)                     :: xcentre      !Fractional coordinates of the centre
       real(kind=cp), dimension(3)                     :: mxcentre     !Refinement codes (or multipliers) of Fractional coordinates of the centre
       integer,       dimension(3)                     :: lxcentre     !Numbers of LSQ parameters for Fractional coordinates of the centre
       real(kind=cp), dimension(3)                     :: Orient       !Orientation angles (Euler angles or variant ...)
       real(kind=cp), dimension(3)                     :: mOrient      !Refinement codes (or multipliers) of Orientation angles (Euler angles or variant ...)
       integer,       dimension(3)                     :: lOrient      !Numbers of LSQ parameters for Orientation angles (Euler angles or variant ...)
       real(kind=cp), dimension(6)                     :: T_TLS        !Translational Thermal factor tensor
       real(kind=cp), dimension(6)                     :: mT_TLS       !Refinement codes (or multipliers) of Translational Thermal factor tensor
       integer,       dimension(6)                     :: lT_TLS       !Numbers of LSQ parameters for Translational Thermal factor tensor
       real(kind=cp), dimension(6)                     :: L_TLS        !Librational Thermal factor tensor
       real(kind=cp), dimension(6)                     :: mL_TLS       !Refinement codes (or multipliers) of Librational Thermal factor tensor
       integer,       dimension(6)                     :: lL_TLS       !Numbers of LSQ parameters for Librational Thermal factor tensor
       real(kind=cp), dimension(3,3)                   :: S_TLS        !TL-correlation Thermal factor
       real(kind=cp), dimension(3,3)                   :: mS_TLS       !Refinement codes (or multipliers) of TL-correlation Thermal factor
       integer,       dimension(3,3)                   :: lS_TLS       !Numbers of LSQ parameters for TL-correlation Thermal factor
       real(kind=cp), dimension(3,3)                   :: Euler        !Euler matrix
       character(len=20), allocatable, dimension(  :)  :: AtName       !Atom Name
       character(len=4),  allocatable, dimension(  :)  :: AtSymb       !Atom species
       integer,           allocatable, dimension(  :)  :: AtZ          !Atomic Number
       integer,           allocatable, dimension(:,:)  :: Ptr          !Pointer to scat.factors (first index -> pattern)
       real(kind=cp),     allocatable, dimension(:,:)  :: I_coor       !Internal coordinates (d,ang,dang)
       real(kind=cp),     allocatable, dimension(:,:)  :: mI_Coor      !Refinement codes (or multipliers) of internal coordinates
       integer,           allocatable, dimension(:,:)  :: lI_coor      !Numbers of LSQ parameters for internal coordinates
       real(kind=cp),     allocatable, dimension(  :)  :: biso         !Isotropic temperature factor
       real(kind=cp),     allocatable, dimension(  :)  :: mbiso        !Refinement codes (or multipliers) of Isotropic temperature factor
       integer,           allocatable, dimension(  :)  :: lbiso        !Numbers of LSQ parameters for Isotropic temperature factor
       real(kind=cp),     allocatable, dimension(  :)  :: occ          !Occupation factor
       real(kind=cp),     allocatable, dimension(  :)  :: mocc         !Refinement codes (or multipliers) of Occupation factor
       integer,           allocatable, dimension(  :)  :: locc         !Numbers of LSQ parameters for Occupation factor
       integer,           allocatable, dimension(  :)  :: Nb           !Number of neighbours
       integer,           allocatable, dimension(:,:)  :: INb          !Index of neighbours
       integer,           allocatable, dimension(:,:)  :: Tb           !Type of Bonds
       integer,           allocatable, dimension(:,:)  :: Conn         !Conectivity (N1,N2,N3)
    End Type Molecule_type

    !!----
    !!----  TYPE :: MOLECULAR_CRYSTAL_TYPE
    !!--..
    !!----  Type, public :: Molecular_Crystal_Type
    !!----     integer                                              :: N_free      !Number of free atoms
    !!----     integer                                              :: N_mol       !Number of Molecules
    !!----     integer                                              :: N_species   !Number of species
    !!----     integer                                              :: Npat        !
    !!----     type(Crystal_Cell_type)                              :: Cell        !Cell Information
    !!----     type(Space_Group_type)                               :: SpG         !Space Group Information
    !!----     type(Atom_type),         allocatable, dimension(  :) :: Atm         !Free Atoms
    !!----     type(Molecule_type ),    allocatable, dimension(  :) :: Mol         !Molecules
    !!----  End type Molecular_Crystal_Type
    !!----
    !!---- Update: February - 2005
    !!
    Type, public :: Molecular_Crystal_Type
       integer                                              :: N_Free
       integer                                              :: N_Mol
       integer                                              :: N_Species
       integer                                              :: Npat
       type(Crystal_Cell_type)                              :: Cell
       type(Space_Group_type)                               :: SpG
       type(Atom_type),         allocatable, dimension(  :) :: Atm
       type(Molecule_type ),    allocatable, dimension(  :) :: Mol
    End type Molecular_Crystal_Type

    !---- Overloading Section ----!
    Interface Empiric_Formula
       Module Procedure Empiric_Formula_FAtom
       Module Procedure Empiric_Formula_Molec
       Module Procedure Empiric_Formula_Molcrys
    End Interface

    Interface Read_Molecule
       Module Procedure Read_Molecule_in_File
       Module Procedure Read_Molecule_in_Var
    End Interface

 Contains
    !---- Subroutines ----!

    !!----
    !!---- Subroutine Cartesian_to_Fractional(Molecule,Cell,NewMolecule)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Crystal_Cell_Type), intent(in)           :: Cell
    !!----    type (Molecule_type), intent(   out), optional :: Newmolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a
    !!----    molecule from cartesian coordinates to  fractional coordinates.
    !!----    If a third argument is present the subroutine creates a new
    !!----    molecule (copy of the old one) with fractional coordinates,
    !!----    preserving the input molecule in Cartesian Coordinates. Otherwise
    !!----    the input molecule is changed on output.
    !!----    Control of error is present
    !!--..       Xc= Euler.Xic  (Cartesian in the crystal frame)
    !!--..       xf= Orth_Cr_cel Xc (fractional before translating to the centre)
    !!--..       Xf = Orth_Cr_cel (Euler.Xic) + Xo (final fractional coordinates)
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Cartesian_to_Fractional(Molecule,Cell,NewMolecule)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Crystal_Cell_Type), intent(in)           :: Cell
       type (Molecule_type), intent(   out), optional :: NewMolecule

       !---- Local Variables ----!
       integer                       :: i,na
       real(kind=cp)                 :: phi,theta,chi
       real(kind=cp), dimension(3)   :: ci,xi
       real(kind=cp), dimension(3,3) :: Eu
       type (Molecule_type)          :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "C") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Cartesian_to_Fractional: the input molecule is not in Cartesian coordinates"
          return
       end if

       na=molecule%natoms
       if (na <=0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Cartesian_to_Fractional: No atoms are defined on molecule variable"
          return
       end if

       if (.not. molecule%in_xtal) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Cartesian_to_Fractional: the input molecule haven't crystal information"
          return
       end if

       !---- Step 1----!
       call init_molecule(newmol,na)
       newmol=molecule

       !---- Frame after a rotation defined by the matrix M(theta,phi,Chi)
       phi   = newmol%orient(1)
       theta = newmol%orient(2)
       chi   = newmol%orient(3)
       if (newmol%is_EulerMat) then
          Eu=newmol%Euler
       else
          call Set_Euler_matrix(newmol%rot_type,phi,theta,chi,Eu)
          newmol%Euler=Eu
          newmol%is_EulerMat=.true.
       end if

       do i=1,na
          ci=matmul(Eu,newmol%I_coor(:,i))         !Cartesian components in the Crystal Frame
          xi=matmul(cell%Orth_Cr_cel,ci)           !Fractional coordinates before translation
          newmol%I_coor(:,i) = newmol%xcentre + xi !Final fractional coordinates
       end do
       newmol%coor_type = "F"

       !---- Step 3 ----!
       if (present(newmolecule)) then
          call Init_molecule(NewMolecule,na)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Cartesian_to_Fractional: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Cartesian_to_Fractional

    !!----
    !!---- Subroutine Cartesian_to_Spherical(Molecule,NewMolecule)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Molecule_type), intent(   out), optional :: Newmolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a
    !!----    molecule from cartesian coordinates to  spherical coordinaters.
    !!----    If a second argument is present the subroutine creates a new
    !!----    molecule (copy of the old one) with spherical coordinates,
    !!----    preserving the input molecule in Cartesian Coordinates. Otherwise
    !!----    the input molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Cartesian_to_Spherical(Molecule,NewMolecule)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Molecule_type), intent(   out), optional :: NewMolecule

       !---- Local variables -----!
       integer                     :: i,na
       real(kind=cp)               :: r, theta, phi
       real(kind=cp), dimension(3) :: ri
       type (Molecule_type)        :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "C") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Cartesian_to_Spherical: the input molecule is not in Cartesian coordinates"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Cartesian_to_Spherical: No atoms are defined"
          return
       end if

       !---- Start calculations for each atom of the molecule ----!
       call init_molecule(newmol,na)
       NewMol=Molecule

       do i=1,na
          ri=Molecule%I_Coor(:,i)
          call  Get_Spheric_Coord(ri,r,theta,phi,"D")
          NewMol%I_Coor(1,i) = r
          NewMol%I_Coor(2,i) = theta
          NewMol%I_Coor(3,i) = phi
       end do
       NewMol%coor_type="S"

       if (present(newmolecule)) then
          call Init_molecule(NewMolecule,na)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Cartesian_to_Spherical: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Cartesian_to_Spherical

    !!----
    !!---- Subroutine Cartesian_to_Zmatrix(Molecule,NewMolecule,Cell, D_min,D_max)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Molecule_type), intent(   out), optional :: NewMolecule
    !!----    Type(Crystal_Cell_Type), intent(in),  optional :: Cell
    !!----    real(kind=cp),        intent(in    ), optional :: D_min
    !!----    real(kind=cp),        intent(in    ), optional :: D_max
    !!----
    !!----    Subroutine to transform the internal coordinates of a molecule
    !!----    from cartesian coordinates to  Z-matrix.
    !!----    If a second argument is present the subroutine creates a new
    !!----    molecule (copy of the old one) with Z-matrix, preserving
    !!----    the input molecule in Cartesian Coordinates. Otherwise the input
    !!----    molecule is changed on output.
    !!----    The input cartesian coordinates may be defined with respect to another
    !!----    internal frame. The final internal frame is that defined for Z-matrices:
    !!----    the x-axis is from the first to the second atom and the x-y plane is formed
    !!----    by the three first atoms. The Euler matrix and the molecular centre in the
    !!----    crystallographic system is changed in consequence.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Cartesian_to_Zmatrix(Molecule,NewMolecule,Cell,D_min,D_max)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Molecule_type), intent(   out), optional :: NewMolecule
       Type(Crystal_Cell_Type), intent(in),  optional :: Cell
       real(kind=cp),        intent(in    ), optional :: D_min
       real(kind=cp),        intent(in    ), optional :: D_max

       !---- Local variables -----!
       integer                       :: i,na,j,k,n,mode
       real(kind=cp)                 :: dist, ang, phi, theta, chi
       real(kind=cp), dimension(3)   :: ci,ri,rj,rk,rn,u1,u2,u3
       real(kind=cp), dimension(3,3) :: Mat, Eu
       type (Molecule_type)          :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "C") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Cartesian_to_Zmatrix: the input molecule is not in Cartesian coordinates"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Cartesian_to_Zmatrix: Not atoms are defined"
          return
       end if

       if (na < 3) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Cartesian_to_Zmatrix: You need at least three atoms"
          return
       end if

       !---- Call Connectivity if necessary ----!
       if (.not. molecule%is_connect) then
          mode=0
          if (present(d_min)) mode=1
          if (present(d_max)) mode=mode + 2
          select case (mode)
             case (0)
                call create_connectivity_cartesian(molecule)
             case (1)
                call create_connectivity_cartesian(molecule,dmin=d_min)
             case (2)
                call create_connectivity_cartesian(molecule,dmax=d_max)
             case (3)
                call create_connectivity_cartesian(molecule,dmin=d_min,dmax=d_max)
          end select
          if (err_molec) then
             ERR_Molec_Mess="Error in Cartesian_to_Zmatrix: the connectivity is wrong"
             return
          end if
          molecule%is_connect=.true.
       end if

       !---- Start calculations for each atom of the molecule ----!
       call init_molecule(newmol,na)
       newmol=molecule

       !---- First atom is always at origin (Z-matrix) ----!
       NewMol%I_Coor(:,1) = 0.0_cp
       NewMol%conn(:,1)   = 0

       !---- Second atom is always along "x" ----!
       ri=molecule%I_coor(:,2)-molecule%I_coor(:,1)
       dist=sqrt(dot_product(ri,ri))
       NewMol%I_Coor(1,2)   = dist
       NewMol%I_Coor(2:3,2) = 0.0_cp
       NewMol%conn(2:3,2)   = 0
       NewMol%conn(1,2)     = 1

       !---- Third atom is always in the "xy" plane ----!
       !---- A(i) d_ij  ang_ijk   dang_ijkl  j k l
       if (NewMol%conn(1,3) == 1) then
          NewMol%conn(2,3) = 2
          NewMol%conn(3,3) = 0
          ri=molecule%I_coor(:,3)-molecule%I_coor(:,1)
          rj=molecule%I_coor(:,2)-molecule%I_coor(:,1)
          dist= sqrt(dot_product(ri,ri))
          ang = acosd(dot_product(ri,rj)/dist/sqrt(dot_product(rj,rj)))
          NewMol%I_coor(1,3) = dist
          NewMol%I_coor(2,3) = ang
          NewMol%I_coor(3,3) = 0.0_cp
       else
          NewMol%conn(1,3) = 2
          NewMol%conn(2,3) = 1
          NewMol%conn(3,3) = 0
          ri=molecule%I_coor(:,3)-molecule%I_coor(:,2)
          rj=molecule%I_coor(:,1)-molecule%I_coor(:,2)
          dist= sqrt(dot_product(ri,ri))
          ang = acosd(dot_product(ri,rj)/dist/sqrt(dot_product(rj,rj)))
          NewMol%I_coor(1,3) = dist
          NewMol%I_coor(2,3) = ang
          NewMol%I_coor(3,3) = 0.0_cp
       end if

       if (Molecule%in_xtal) then    !Modify the Euler matrix, orientation angles and centre
          if (Molecule%is_EulerMat) then
             Eu=Molecule%Euler
          else
             phi=Molecule%orient(1)
             theta=Molecule%orient(2)
             chi=Molecule%orient(3)
             Call Set_Euler_matrix(Molecule%rot_type,phi,theta,chi,Eu)
          end if
          newmol%Euler=Eu
          newmol%is_EulerMat=.true.

          ri=molecule%I_coor(:,1)
          rj=molecule%I_coor(:,2)
          rk=molecule%I_coor(:,3)
          u1=rj-ri
          u1=u1/sqrt(dot_product(u1,u1))
          u3=cross_product(u1,rk-ri)
          u3=u3/sqrt(dot_product(u3,u3))
          u2=cross_product(u3,u1)
          Mat(:,1)=u1
          Mat(:,2)=u2  !Active matrix needed to get the new Euler matrix
          Mat(:,3)=u3

          newmol%Euler=matmul(Eu,Mat)  !New Euler Matrix
          call Get_PhiTheChi(newmol%Euler,Phi,Theta,Chi,"D")
          newmol%orient(1)=  phi
          newmol%orient(2)=theta
          newmol%orient(3)=  chi

          !---- New centre (?) Needs the Cell argument
          if (present(Cell)) then
             rj=Matmul(Mat,ri)
             newmol%xcentre=matmul(Cell%Orth_Cr_cel,rj)+molecule%xcentre
          else
             if (dot_product(ri,ri) > eps) then
                err_molec=.true.
                ERR_Molec_Mess="Error in Cartesian_to_Zmatrix: First atom not at the origin => a cell has to be provided "
                return
             end if
          end if
       end if

       do i=4,na                      !The result of this calculation is independent of the type of
          ri = molecule%I_coor(:,i)   !cartesian coordinates => it is not needed to transforn the input Cartesian!
          j  = molecule%conn(1,i)     !The connectivity is needed for the Z-matrix description
          k  = molecule%conn(2,i)     !If the connectivity is given it is possible to transform to
          n  = molecule%conn(3,i)     !Z-matrix if cartesian/spherical coordinates are given.
          if ( j == 0 .or. k == 0 .or. n == 0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Cartesian_to_Zmatrix: the connectivity is wrong for atom: " &
                            //molecule%Atname(i)
             return
          end if
          rj = molecule%I_coor(:,j)
          rk = molecule%I_coor(:,k)
          rn = molecule%I_coor(:,n)
          call get_Z_from_cartesian(ci,ri,rj,rk,rn)
          NewMol%I_coor(:,i) = ci
       end do
       NewMol%coor_type="Z"

       if (present(NewMolecule)) then
          call Init_molecule(NewMolecule,na)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Cartesian_to_Zmatrix: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Cartesian_to_Zmatrix

    !!--++
    !!--++ Subroutine Create_Connectivity_Cartesian(Molecule, Dmin, Dmax)
    !!--++    type (Molecule_type),          intent(in out):: Molecule
    !!--++    real(kind=cp), optional,       intent(in)    :: Dmin
    !!--++    real(kind=cp), optional,       intent(in)    :: Dmax
    !!--++
    !!--++    (PRIVATE)
    !!--++    Subroutine that create the connectivity for the molecule.
    !!--++    The coordinates must be in Cartesian system. Control of
    !!--++    error is implemented.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Create_Connectivity_Cartesian(Molecule,Dmin,Dmax)
       !---- Arguments ----!
       type (Molecule_type),          intent(in out):: Molecule
       real(kind=cp), optional,       intent(in)    :: Dmin
       real(kind=cp), optional,       intent(in)    :: Dmax

       !---- Local variables ----!
       logical                                                      :: re_order
       integer                                                      :: i,j,k,l,m,nc1,nc2,nc3
       integer, dimension(molecule%natoms,molecule%natoms)          :: T_Conn
       integer, dimension(3,molecule%natoms)                        :: T_N
       integer, dimension(molecule%natoms)                          :: T_Ind
       real(kind=cp),    dimension(molecule%natoms,molecule%natoms) :: T_Dist
       real(kind=cp)                                                :: d_min, d_max
       real(kind=cp)                                                :: dist !,ang,tors
       type (Molecule_type)                                         :: Newmol


       !---- Initialize ----!
       d_min=0.6
       d_max=3.0
       T_Conn=0
       T_N   =0
       T_Ind =0
       T_Dist=0.0
       if (present(dmin)) d_min=dmin
       if (present(dmax)) d_max=dmax

       !---- Controls ----!
       if (molecule%coor_type /= "C") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Connectivity: the input molecule is not in Cartesian coordinates"
          return
       end if

       !---- Creating Tables ----!
       do i=1,molecule%natoms
          do j=i+1,molecule%natoms
             dist=distance(molecule%I_coor(:,i),molecule%I_coor(:,j))
             if (dist < d_min .or. dist > d_max) cycle
             if (adjustl(molecule%Atsymb(i)) == "H   " .and. &
                 adjustl(molecule%Atsymb(j)) == "H   ") cycle
             T_Conn(i,j)=i
             T_Conn(j,i)=i
             T_Dist(i,j)=dist
             T_Dist(j,i)=dist
          end do
       end do

       !---- Test for reorder atoms ----!
       re_order=.false.

       do i=2,molecule%natoms
          j=count(T_conn(i,1:i-1) > 0)
          if (j==0) re_order=.true.
       end do

       if (re_order) then
          m=1
          T_ind(m)=1
          do i=1,molecule%natoms
             do j=1,molecule%natoms
                if (T_Conn(i,j) <= 0) cycle
                l=0
                do k=1,m
                     if (j == T_ind(k)) then
                        l=1
                        exit
                     end if
                end do
                if (l > 0) cycle
                m=m+1
                T_ind(m)=j
             end do
          end do

          call init_molecule(newmol,molecule%natoms)
          newmol=molecule
          do i=2,newmol%natoms
             j=T_ind(i)
             newmol%AtName(i)=   molecule%AtName(j)
             newmol%AtSymb(i)=   molecule%AtSymb(j)
             newmol%AtZ(i)=      molecule%AtZ(j)
             newmol%Ptr(:,i)=    molecule%Ptr(:,j)
             newmol%I_Coor(:,i)= molecule%I_Coor(:,j)
             newmol%mI_Coor(:,i)=molecule%mI_Coor(:,j)
             newmol%lI_Coor(:,i)=molecule%lI_Coor(:,j)
             newmol%biso(i)=     molecule%biso(j)
             newmol%mbiso(i)=    molecule%mbiso(j)
             newmol%lbiso(i)=    molecule%lbiso(j)
             newmol%occ(i)=      molecule%occ(j)
             newmol%mocc(i)=     molecule%mocc(j)
             newmol%locc(i)=     molecule%locc(j)
             newmol%nb(i)=       molecule%nb(j)
             newmol%Inb(:,i)=    molecule%Inb(:,j)
             newmol%Tb(:,i)=     molecule%Tb(:,j)
             newmol%Conn(:,i)=   molecule%Conn(:,j)
          end do
          molecule=newmol
          call init_molecule(newmol,0)

          T_Conn=0
          T_Dist=0.0
          do i=1,molecule%natoms
             do j=i+1,molecule%natoms
                dist=distance(molecule%I_coor(:,i),molecule%I_coor(:,j))
                if (dist < d_min .or. dist > d_max) cycle
                if (adjustl(molecule%Atsymb(i)) == "H   " .and. &
                    adjustl(molecule%Atsymb(j)) == "H   ") cycle
                T_Conn(i,j)=i
                T_Conn(j,i)=i
                T_Dist(i,j)=dist
                T_Dist(j,i)=dist
             end do
          end do
       end if

       !---- Connectivity Info ----!
       do i=2, molecule%natoms

          !---- Distances: Fill N1 ----!
          j=minloc(T_Dist(i,1:i-1),dim=1,mask=(T_Dist(i,1:i-1) > 0.0))
          T_N(1,i)=j

          if (j == 0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Connectivity: Some Index are zeros"
             return
          end if

          !---- Angles: Fill N2 ----!
          if (i > 2) then
             nc1=count((T_Conn(j,1:i-1) > 0 .and. T_Conn(j,1:i-1) /=j),dim=1)
             nc2=count((T_Conn(i,1:i-1) > 0 .and. T_Conn(i,1:i-1) /=j),dim=1)
             k=0
             if (nc1 > 0) then
                do
                   k=minloc(T_Dist(j,1:i-1),dim=1, mask=(T_Dist(j,1:i-1) > 0.0))
                   if (k == j) then
                      T_Dist(j,k)=-T_Dist(j,k)
                      cycle
                   else
                      exit
                   end if
                end do
             elseif (nc2 > 0) then
                do
                   k=minloc(T_Dist(i,1:i-1),dim=1, mask=(T_Dist(i,1:i-1) > 0.0))
                   if (k == j) then
                      T_Dist(i,k)=-T_Dist(i,k)
                      cycle
                   else
                      exit
                   end if
                end do
             end if
             if (k == 0) then
                !---- Elegir uno cualquiera ----!
                do l=1,i-1
                   if (l == j) cycle
                   k=l
                   exit
                end do
             end if
             T_N(2,i)=k
          end if
          T_Dist=abs(T_Dist)

          !---- Torsion ----!
          if (i > 3) then
             nc1=count((T_Conn(k,1:i-1) > 0 .and. T_Conn(k,1:i-1) /=j .and. T_Conn(k,1:i-1) /=k),dim=1)
             nc2=count((T_Conn(j,1:i-1) > 0 .and. T_Conn(j,1:i-1) /=j .and. T_Conn(j,1:i-1) /=k),dim=1)
             nc3=count((T_Conn(i,1:i-1) > 0 .and. T_Conn(i,1:i-1) /=j .and. T_Conn(i,1:i-1) /=k),dim=1)

             l=0
             if (nc1 > 0) then
                do
                   l=minloc(T_Dist(k,1:i-1),dim=1, mask=(T_Dist(k,1:i-1) > 0.0))
                   if (l == j .or. l == k) then
                      T_Dist(k,l)=-T_Dist(k,l)
                      cycle
                   else
                      exit
                   end if
                end do
             elseif (nc2 > 0) then
                do
                   l=minloc(T_Dist(j,1:i-1),dim=1, mask=(T_Dist(j,1:i-1) > 0.0))
                   if (l == j .or. l == k) then
                      T_Dist(j,l)=-T_Dist(j,l)
                      cycle
                   else
                      exit
                   end if
                end do
             elseif (nc3 > 0) then
                do
                   l=minloc(T_Dist(i,1:i-1),dim=1, mask=(T_Dist(i,1:i-1) > 0.0))
                   if (l == j .or. l == k) then
                      T_Dist(i,l)=-T_Dist(i,l)
                      cycle
                   else
                      exit
                   end if
                end do
             end if
             if (l==0) then
                !---- Elegir uno cualquiera ----!
                do m=1,i-1
                   if (m == j .or. m == k) cycle
                   l=m
                   exit
                end do
             end if
             T_N(3,i)=l
          end if
          T_Dist=abs(T_Dist)

       end do

       !---- Final Part ----!
       do i=1, molecule%natoms
          molecule%Conn(:,i)=T_N(:,i)
          select case (i)
             case (2)
                if (T_N(1,i) == 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error in Connectivity: Some Index are zeros"
                end if
             case (3)
                if (any(T_N(1:2,i) == 0)) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error in Connectivity: Some Index are zeros"
                end if
             case (4:)
                if (any(T_N(:,i) == 0)) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error in Connectivity: Some Index are zeros"
                end if
          end select
       end do

       return
    End Subroutine Create_Connectivity_Cartesian

    !!----
    !!---- Subroutine Fix_Reference(Molecule, NewMolecule, NAtom_O, NAtom_X, NAtom_XY)
    !!----    type (Molecule_type),     intent(in out)           :: Molecule
    !!----    type (Molecule_type),     intent(   out), optional :: Newmolecule
    !!----    integer,                  intent(in),     optional :: NAtom_O
    !!----    integer,                  intent(in),     optional :: NAtom_X
    !!----    integer,                  intent(in),     optional :: NAtom_XY
    !!----
    !!----    Subroutine to order the molecule choosing which atom is the origin,
    !!----    which define the X axis and which defines the XY Plane
    !!----    If the second argument is present the subroutine creates a new molecule
    !!----    preserving the input molecule in Cartesian. Otherwise the input molecule is
    !!----    changed on output.
    !!----    If Natom_0 is absent, then the first atom on the molecule will be the origin.
    !!----    If Natom_X is absent, then the second atom on the molecule will define the X axis.
    !!----    If Natom_XY is absent, then the third atom on the molecule will define the XY Plane.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Fix_Reference(Molecule, NewMolecule, NAtom_O, NAtom_X, NAtom_XY)
       !---- Arguments ----!
       type (Molecule_type),     intent(in out)           :: Molecule
       type (Molecule_type),     intent(   out), optional :: Newmolecule
       integer,                  intent(in),     optional :: NAtom_O
       integer,                  intent(in),     optional :: NAtom_X
       integer,                  intent(in),     optional :: NAtom_XY

       !---- Local variables ----!
       integer                   :: n_or, n_x, n_xy
       integer                   :: i
       type (Molecule_type)      :: Newmol,SetMol

       !---- Initialize ----!
       n_or=1
       n_x =2
       n_xy=3
       if (present(natom_O))  n_or=natom_o
       if (present(natom_x))  n_x =natom_x
       if (present(natom_xy)) n_xy=natom_xy

       call Init_Err_Molec()

       if (n_x == n_or) then
          err_molec=.true.
          ERR_Molec_Mess="The atom defining origin and X axis is the same"
          return
       end if

       if (n_xy == n_or .or. n_xy ==n_x) then
          err_molec=.true.
          ERR_Molec_Mess="The atom defining the Plane XY is equal to the origin or that define the X axis"
          return
       end if

       if (molecule%natoms > 0) call Init_Molecule(Newmol,molecule%natoms)
       Newmol=molecule

       !---- Sorting the Atom List ----!
       call init_molecule(SetMol,1)

       !---- Fix Origin ----!
       if (n_or /= 1) then
          SetMol%AtName(1)    =NewMol%AtName(n_or)
          SetMol%AtSymb(1)    =NewMol%AtSymb(n_or)
          SetMol%AtZ(1)       =NewMol%AtZ(n_or)
          SetMol%Ptr(:,1)     =NewMol%Ptr(:,n_or)
          SetMol%I_Coor(:,1)  =NewMol%I_Coor(:,n_or)
          SetMol%mI_Coor(:,1) =NewMol%mI_Coor(:,n_or)
          SetMol%lI_Coor(:,1) =NewMol%lI_Coor(:,n_or)
          SetMol%Biso(1)      =NewMol%Biso(n_or)
          SetMol%mbiso(1)     =NewMol%mbiso(n_or)
          SetMol%lBiso(1)     =NewMol%lBiso(n_or)
          SetMol%Occ(1)       =NewMol%Occ(n_or)
          SetMol%mocc(1)      =NewMol%mocc(n_or)
          SetMol%lOcc(1)      =NewMol%lOcc(n_or)
          SetMol%Nb(1)        =NewMol%Nb(n_or)
          SetMol%INb(:,1)     =NewMol%INb(:,n_or)
          SetMol%Tb(:,1)      =NewMol%Tb(:,n_or)
          SetMol%Conn(:,1)    =NewMol%Conn(:,n_or)

          NewMol%AtName(2:n_or)    =NewMol%AtName(1:n_or-1)
          NewMol%AtSymb(2:n_or)    =NewMol%AtSymb(1:n_or-1)
          NewMol%AtZ(2:n_or)       =NewMol%AtZ(1:n_or-1)
          NewMol%Ptr(:,2:n_or)     =NewMol%Ptr(:,1:n_or-1)
          NewMol%I_Coor(:,2:n_or)  =NewMol%I_Coor(:,1:n_or-1)
          NewMol%mI_Coor(:,2:n_or) =NewMol%mI_Coor(:,1:n_or-1)
          NewMol%lI_Coor(:,2:n_or) =NewMol%lI_Coor(:,1:n_or-1)
          NewMol%Biso(2:n_or)      =NewMol%Biso(1:n_or-1)
          NewMol%mbiso(2:n_or)     =NewMol%mbiso(1:n_or-1)
          NewMol%lBiso(2:n_or)     =NewMol%lBiso(1:n_or-1)
          NewMol%Occ(2:n_or)       =NewMol%Occ(1:n_or-1)
          NewMol%mocc(2:n_or)      =NewMol%mocc(1:n_or-1)
          NewMol%lOcc(2:n_or)      =NewMol%lOcc(1:n_or-1)
          NewMol%Nb(2:n_or)        =NewMol%Nb(1:n_or-1)
          NewMol%INb(:,2:n_or)     =NewMol%INb(:,1:n_or-1)
          NewMol%Tb(:,2:n_or)      =NewMol%Tb(:,1:n_or-1)
          NewMol%Conn(:,2:n_or)    =NewMol%Conn(:,1:n_or-1)

          NewMol%AtName(1)    =SetMol%AtName(1)
          NewMol%AtSymb(1)    =SetMol%AtSymb(1)
          NewMol%AtZ(1)       =SetMol%AtZ(1)
          NewMol%Ptr(:,1)     =SetMol%Ptr(:,1)
          NewMol%I_Coor(:,1)  =SetMol%I_Coor(:,1)
          NewMol%mI_Coor(:,1) =SetMol%mI_Coor(:,1)
          NewMol%lI_Coor(:,1) =SetMol%lI_Coor(:,1)
          NewMol%Biso(1)      =SetMol%Biso(1)
          NewMol%mbiso(1)     =SetMol%mbiso(1)
          NewMol%lBiso(1)     =SetMol%lBiso(1)
          NewMol%Occ(1)       =SetMol%Occ(1)
          NewMol%mocc(1)      =SetMol%mocc(1)
          NewMol%lOcc(1)      =SetMol%lOcc(1)
          NewMol%Nb(1)        =SetMol%Nb(1)
          NewMol%INb(:,1)     =SetMol%INb(:,1)
          NewMol%Tb(:,1)      =SetMol%Tb(:,1)
          NewMol%Conn(:,1)    =SetMol%Conn(:,1)

          if (Newmol%is_connect) then
             do i=1,n_or
                if (newmol%conn(1,i) == n_or) then
                   newmol%conn(1,i)=1
                else if (newmol%conn(1,i) < n_or) then
                   newmol%conn(1,i)=newmol%conn(1,i)+1
                end if

                if (newmol%conn(2,i) == n_or) then
                   newmol%conn(2,i)=1
                else if (newmol%conn(2,i) < n_or) then
                   newmol%conn(2,i)=newmol%conn(2,i)+1
                end if

                if (newmol%conn(3,i) == n_or) then
                   newmol%conn(3,i)=1
                else if (newmol%conn(3,i) < n_or) then
                   newmol%conn(3,i)=newmol%conn(3,i)+1
                end if
             end do
          end if

          if (n_x < n_or) then
             n_x=n_x+1
          end if
          if (n_xy < n_or) then
             n_xy=n_xy+1
          end if
       end if

       !---- Fix X Axis ----!
       if (n_x /= 2) then

          SetMol%AtName(1)    =NewMol%AtName(n_x)
          SetMol%AtSymb(1)    =NewMol%AtSymb(n_x)
          SetMol%AtZ(1)       =NewMol%AtZ(n_x)
          SetMol%Ptr(:,1)     =NewMol%Ptr(:,n_x)
          SetMol%I_Coor(:,1)  =NewMol%I_Coor(:,n_x)
          SetMol%mI_Coor(:,1) =NewMol%mI_Coor(:,n_x)
          SetMol%lI_Coor(:,1) =NewMol%lI_Coor(:,n_x)
          SetMol%Biso(1)      =NewMol%Biso(n_x)
          SetMol%mbiso(1)     =NewMol%mbiso(n_x)
          SetMol%lBiso(1)     =NewMol%lBiso(n_x)
          SetMol%Occ(1)       =NewMol%Occ(n_x)
          SetMol%mocc(1)      =NewMol%mocc(n_x)
          SetMol%lOcc(1)      =NewMol%lOcc(n_x)
          SetMol%Nb(1)        =NewMol%Nb(n_x)
          SetMol%INb(:,1)     =NewMol%INb(:,n_x)
          SetMol%Tb(:,1)      =NewMol%Tb(:,n_x)
          SetMol%Conn(:,1)    =NewMol%Conn(:,n_x)

          NewMol%AtName(3:n_x)    =NewMol%AtName(2:n_x-1)
          NewMol%AtSymb(3:n_x)    =NewMol%AtSymb(2:n_x-1)
          NewMol%AtZ(3:n_x)       =NewMol%AtZ(2:n_x-1)
          NewMol%Ptr(:,3:n_x)     =NewMol%Ptr(:,2:n_x-1)
          NewMol%I_Coor(:,3:n_x)  =NewMol%I_Coor(:,2:n_x-1)
          NewMol%mI_Coor(:,3:n_x) =NewMol%mI_Coor(:,2:n_x-1)
          NewMol%lI_Coor(:,3:n_x) =NewMol%lI_Coor(:,2:n_x-1)
          NewMol%Biso(3:n_x)      =NewMol%Biso(2:n_x-1)
          NewMol%mbiso(3:n_x)     =NewMol%mbiso(2:n_x-1)
          NewMol%lBiso(3:n_x)     =NewMol%lBiso(2:n_x-1)
          NewMol%Occ(3:n_x)       =NewMol%Occ(2:n_x-1)
          NewMol%mocc(3:n_x)      =NewMol%mocc(2:n_x-1)
          NewMol%lOcc(3:n_x)      =NewMol%lOcc(2:n_x-1)
          NewMol%Nb(3:n_x)        =NewMol%Nb(2:n_x-1)
          NewMol%INb(:,3:n_x)     =NewMol%INb(:,2:n_x-1)
          NewMol%Tb(:,3:n_x)      =NewMol%Tb(:,2:n_x-1)
          NewMol%Conn(:,3:n_x)    =NewMol%Conn(:,2:n_x-1)

          NewMol%AtName(2)    =SetMol%AtName(1)
          NewMol%AtSymb(2)    =SetMol%AtSymb(1)
          NewMol%AtZ(2)       =SetMol%AtZ(1)
          NewMol%Ptr(:,2)     =SetMol%Ptr(:,1)
          NewMol%I_Coor(:,2)  =SetMol%I_Coor(:,1)
          NewMol%mI_Coor(:,2) =SetMol%mI_Coor(:,1)
          NewMol%lI_Coor(:,2) =SetMol%lI_Coor(:,1)
          NewMol%Biso(2)      =SetMol%Biso(1)
          NewMol%mbiso(2)     =SetMol%mbiso(1)
          NewMol%lBiso(2)     =SetMol%lBiso(1)
          NewMol%Occ(2)       =SetMol%Occ(1)
          NewMol%mocc(2)      =SetMol%mocc(1)
          NewMol%lOcc(2)      =SetMol%lOcc(1)
          NewMol%Nb(2)        =SetMol%Nb(1)
          NewMol%INb(:,2)     =SetMol%INb(:,1)
          NewMol%Tb(:,2)      =SetMol%Tb(:,1)
          NewMol%Conn(:,2)    =SetMol%Conn(:,1)

          if (Newmol%is_connect) then
             do i=1,n_x
                if (newmol%conn(1,i) == n_x) then
                   newmol%conn(1,i)=2
                else if (newmol%conn(1,i) < n_x .and. newmol%conn(1,i) > 1) then
                   newmol%conn(1,i)=newmol%conn(1,i)+1
                end if

                if (newmol%conn(2,i) == n_x) then
                   newmol%conn(2,i)=2
                else if (newmol%conn(2,i) < n_x .and. newmol%conn(2,i) > 1) then
                   newmol%conn(2,i)=newmol%conn(2,i)+1
                end if

                if (newmol%conn(3,i) == n_x) then
                   newmol%conn(3,i)=2
                else if (newmol%conn(3,i) < n_x .and. newmol%conn(3,i) > 1) then
                   newmol%conn(3,i)=newmol%conn(3,i)+1
                end if
             end do
          end if
          if (n_xy < n_x) then
             n_xy=n_xy+1
          end if
       end if

       !---- Fix XY Plane ----!
       if (n_xy /= 3) then

          SetMol%AtName(1)    =NewMol%AtName(n_xy)
          SetMol%AtSymb(1)    =NewMol%AtSymb(n_xy)
          SetMol%AtZ(1)       =NewMol%AtZ(n_xy)
          SetMol%Ptr(:,1)     =NewMol%Ptr(:,n_xy)
          SetMol%I_Coor(:,1)  =NewMol%I_Coor(:,n_xy)
          SetMol%mI_Coor(:,1) =NewMol%mI_Coor(:,n_xy)
          SetMol%lI_Coor(:,1) =NewMol%lI_Coor(:,n_xy)
          SetMol%Biso(1)      =NewMol%Biso(n_xy)
          SetMol%mbiso(1)     =NewMol%mbiso(n_xy)
          SetMol%lBiso(1)     =NewMol%lBiso(n_xy)
          SetMol%Occ(1)       =NewMol%Occ(n_xy)
          SetMol%mocc(1)      =NewMol%mocc(n_xy)
          SetMol%lOcc(1)      =NewMol%lOcc(n_xy)
          SetMol%Nb(1)        =NewMol%Nb(n_xy)
          SetMol%INb(:,1)     =NewMol%INb(:,n_xy)
          SetMol%Tb(:,1)      =NewMol%Tb(:,n_xy)
          SetMol%Conn(:,1)    =NewMol%Conn(:,n_xy)

          NewMol%AtName(4:n_xy)    =NewMol%AtName(3:n_xy-1)
          NewMol%AtSymb(4:n_xy)    =NewMol%AtSymb(3:n_xy-1)
          NewMol%AtZ(4:n_xy)       =NewMol%AtZ(3:n_xy-1)
          NewMol%Ptr(:,4:n_xy)     =NewMol%Ptr(:,3:n_xy-1)
          NewMol%I_Coor(:,4:n_xy)  =NewMol%I_Coor(:,3:n_xy-1)
          NewMol%mI_Coor(:,4:n_xy) =NewMol%mI_Coor(:,3:n_xy-1)
          NewMol%lI_Coor(:,4:n_xy) =NewMol%lI_Coor(:,3:n_xy-1)
          NewMol%Biso(4:n_xy)      =NewMol%Biso(3:n_xy-1)
          NewMol%mbiso(4:n_xy)     =NewMol%mbiso(3:n_xy-1)
          NewMol%lBiso(4:n_xy)     =NewMol%lBiso(3:n_xy-1)
          NewMol%Occ(4:n_xy)       =NewMol%Occ(3:n_xy-1)
          NewMol%mocc(4:n_xy)      =NewMol%mocc(3:n_xy-1)
          NewMol%lOcc(4:n_xy)      =NewMol%lOcc(3:n_xy-1)
          NewMol%Nb(4:n_xy)        =NewMol%Nb(3:n_xy-1)
          NewMol%INb(:,4:n_xy)     =NewMol%INb(:,3:n_xy-1)
          NewMol%Tb(:,4:n_xy)      =NewMol%Tb(:,3:n_xy-1)
          NewMol%Conn(:,4:n_xy)    =NewMol%Conn(:,3:n_xy-1)

          NewMol%AtName(3)    =SetMol%AtName(1)
          NewMol%AtSymb(3)    =SetMol%AtSymb(1)
          NewMol%AtZ(3)       =SetMol%AtZ(1)
          NewMol%Ptr(:,3)     =SetMol%Ptr(:,1)
          NewMol%I_Coor(:,3)  =SetMol%I_Coor(:,1)
          NewMol%mI_Coor(:,3) =SetMol%mI_Coor(:,1)
          NewMol%lI_Coor(:,3) =SetMol%lI_Coor(:,1)
          NewMol%Biso(3)      =SetMol%Biso(1)
          NewMol%mbiso(3)     =SetMol%mbiso(1)
          NewMol%lBiso(3)     =SetMol%lBiso(1)
          NewMol%Occ(3)       =SetMol%Occ(1)
          NewMol%mocc(3)      =SetMol%mocc(1)
          NewMol%lOcc(3)      =SetMol%lOcc(1)
          NewMol%Nb(3)        =SetMol%Nb(1)
          NewMol%INb(:,3)     =SetMol%INb(:,1)
          NewMol%Tb(:,3)      =SetMol%Tb(:,1)
          NewMol%Conn(:,3)    =SetMol%Conn(:,1)

          if (Newmol%is_connect) then
             do i=1,n_xy
                if (newmol%conn(1,i) == n_xy) then
                   newmol%conn(1,i)=3
                else if (newmol%conn(1,i) < n_xy .and. newmol%conn(1,i) > 2) then
                   newmol%conn(1,i)=newmol%conn(1,i)+1
                end if

                if (newmol%conn(2,i) == n_xy) then
                   newmol%conn(2,i)=3
                else if (newmol%conn(2,i) < n_xy .and. newmol%conn(2,i) > 2) then
                   newmol%conn(2,i)=newmol%conn(2,i)+1
                end if

                if (newmol%conn(3,i) == n_xy) then
                   newmol%conn(3,i)=3
                else if (newmol%conn(3,i) < n_xy .and. newmol%conn(3,i) > 2) then
                   newmol%conn(3,i)=newmol%conn(3,i)+1
                end if
             end do
          end if
       end if

       if (present(Newmolecule)) then
          call Init_molecule(NewMolecule,Newmol%natoms)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Fix_Reference: The optional variable was not dimensioned!"
             return
          end if
          Newmolecule=Newmol
       else
          Molecule=Newmol
       end if

       return
    End Subroutine Fix_Reference

    !!----
    !!---- Subroutine Fix_Orient_Cartesian(Molecule, NewMolecule, NAtom_O, NAtom_X, NAtom_XY,Mat)
    !!----    type (Molecule_type),     intent(in out)           :: Molecule
    !!----    type (Molecule_type),     intent(   out), optional :: Newmolecule
    !!----    integer,                  intent(in),     optional :: NAtom_O
    !!----    integer,                  intent(in),     optional :: NAtom_X
    !!----    integer,                  intent(in),     optional :: NAtom_XY
    !!----    real(kind=cp),dimension(3,3),intent(out), optional :: Mat
    !!----
    !!----    Subroutine to transform the Cartesian coordinates of the molecule choosing
    !!----    which atom is the origin, which define the X axis and which defines the XY Plane
    !!----    If the second argument is present the subroutine creates a new molecule
    !!----    preserving the input molecule in Cartesian. Otherwise the input molecule is
    !!----    changed on output.
    !!----    If Natom_0 is absent, then the first atom on the molecule will be the origin.
    !!----    If Natom_X is absent, then the second atom on the molecule will define the X axis.
    !!----    If Natom_XY is absent, then the third atom on the molecule will define the XY Plane.
    !!----    The optional output matrix Mat is the active rotation matrix passing from the old
    !!----    Cartesian frame to the new one. The transpose matrix has served to transform the
    !!----    original Cartesian coordinates.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Fix_Orient_Cartesian(Molecule, NewMolecule, NAtom_O, NAtom_X, NAtom_XY,Mat)
       !---- Arguments ----!
       type (Molecule_type),     intent(in out)           :: Molecule
       type (Molecule_type),     intent(   out), optional :: Newmolecule
       integer,                  intent(in),     optional :: NAtom_O
       integer,                  intent(in),     optional :: NAtom_X
       integer,                  intent(in),     optional :: NAtom_XY
       real(kind=cp),dimension(3,3),intent(out), optional :: Mat

       !---- Local variables ----!
       integer                       :: n_or, n_x, n_xy
       integer                       :: i
       real(kind=cp),dimension(3)    :: u1,u2,u3
       real(kind=cp),dimension(3,3)  :: R
       type (Molecule_type)          :: Newmol

       n_or=1
       n_x =2
       n_xy=3
       if (present(natom_O))  n_or=natom_o
       if (present(natom_x))  n_x =natom_x
       if (present(natom_xy)) n_xy=natom_xy

       if (molecule%natoms > 0) call Init_Molecule(Newmol,molecule%natoms)
       call Fix_Reference(Molecule,Newmol,n_or,n_x,n_xy)
       if (err_molec) return

       !---- Traslation the Origin ----!
       do i=2,Newmol%natoms
          newmol%I_coor(:,i)=newmol%I_coor(:,i)-newmol%I_coor(:,1)
       end do
       newmol%I_coor(:,1)=0.0

       u1=Newmol%I_coor(:,2)
       u1=u1/sqrt(dot_product(u1,u1))
       u2=Newmol%I_coor(:,3)
       u3=cross_product(u1,u2)
       u3=u3/sqrt(dot_product(u3,u3))
       u2=cross_product(u3,u1)
       R(1,:)=u1
       R(2,:)=u2  !Passive matrix needed to get the new coordinates
       R(3,:)=u3  !The active matrix can be output in the optional argument
       if (present(Mat)) Mat=transpose(R)

       do i=2,Newmol%natoms
          newmol%I_coor(:,i)=matmul(R,newmol%I_coor(:,i))
       end do

       if (present(Newmolecule)) then
          if (NewMol%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Fix_Orient_Cartesian: The optional variable was not dimensioned!"
             return
          end if
          Newmolecule=Newmol
       else
          Molecule=Newmol
       end if

       return
    End Subroutine Fix_Orient_Cartesian

    !!----
    !!---- Subroutine Empiric_Formula(Atm/Molcrys/Molecule,Formula,Form_Weight)
    !!----    type(Atom_List_Type),          intent(in)  :: Atm
    !!----    or
    !!----    type(molecular_crystal_type),  intent(in)  :: Molcrys
    !!----    or
    !!----    type(molecule_type),           intent(in)  :: Molecule
    !!----    character(len=*),              intent(out) :: Formula
    !!----    real(kind=cp), optional,       intent(out) :: Form_Weight
    !!----
    !!----    Obtain the Empiric Formula from Atm/Molcrys/Molecule variable
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Empiric_Formula_FAtom(Atm,Formula,Form_Weight)
    !!--++    type(Atom_List_Type),    intent(in)  :: Atm
    !!--++    character(len=*),        intent(out) :: Formula
    !!--++    real(kind=cp), optional, intent(out) :: Form_Weight
    !!--++
    !!--++    (OVERLOADED)
    !!--++    Obtain the Empiric Formula from Atm variable
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Empiric_Formula_FAtom(Atm,Formula,Form_Weight)
       !---- Arguments ----!
       type(Atom_List_Type),     intent(in)  :: Atm
       character(len=*),         intent(out) :: Formula
       real(kind=cp), optional,  intent(out) :: Form_Weight

       !---- Local variables ----!
       character(len=2)                  :: car
       character(len=5)                  :: numcar
       integer                           :: i,j
       integer, dimension(Num_Chem_Info) :: N_PT
       real(kind=cp)                     :: weight

       !---- Init ----!
       N_PT=0
       weight=0.0

       Formula=" "
       if (Atm%natoms <= 0) then
          if (present(Form_weight)) Form_weight=0.0
          return
       end if

       !---- Set Information Table ----!
       call Set_Chem_Info()

       do i=1,atm%natoms
          car=atm%atom(i)%chemsymb
          car=u_case(car)
          do j=1,Num_Chem_Info
             if (car == Chem_Info(j)%Symb) then
                n_pt(j)=n_pt(j)+1
                exit
             end if
          end do
       end do

       if (all (n_pt ==0)) then
          if (present(Form_weight)) Form_weight=0.0
          call Remove_Chem_Info()
          return
       end if

       do i=1,Num_Chem_Info
          if (n_pt(i) == 0) cycle
          car=Chem_Info(i)%Symb
          car(2:2)=l_case(car(2:2))
          write(unit=numcar,fmt="(i5)") n_pt(i)
          Formula=trim(Formula)//trim(car)//adjustl(numcar)
          weight=weight+n_pt(i)*Chem_Info(i)%atwe
       end do

       call Remove_Chem_Info()

       if (present(Form_weight)) Form_weight=weight

       return
    End Subroutine Empiric_Formula_FAtom

    !!--++
    !!--++ Subroutine Empiric_Formula_Molcrys(Molcrys,Formula,Form_Weight)
    !!--++    type(molecular_crystal_type), intent(in)  :: Molcrys
    !!--++    character(len=*),             intent(out) :: Formula
    !!--++    real(kind=cp), optional,      intent(out) :: Form_Weight
    !!--++
    !!--++    (Overloaded)
    !!--++    Obtain the Empiric Formula from Molecule variable and
    !!--++    the Weight is the variable is present.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Empiric_Formula_Molcrys(Molcrys,Formula,Form_Weight)
       !---- Arguments ----!
       type(molecular_crystal_type), intent(in)  :: Molcrys
       character(len=*),             intent(out) :: Formula
       real(kind=cp), optional,      intent(out) :: Form_Weight

       !---- Local variables ----!
       character(len=2)                  :: car
       character(len=5)                  :: numcar
       integer                           :: i,j,k
       integer, dimension(Num_Chem_Info) :: N_PT
       real(kind=cp)                     :: weight


       !---- Init ----!
       N_PT=0
       weight=0.0

       Formula=" "

       if (molcrys%n_free <= 0 .and. molcrys%n_mol <=0) then
          if (present(Form_weight)) Form_weight=0.0
          return
       end if

       !---- Set Information Table ----!
       call Set_Chem_Info()

       do i=1,molcrys%n_free
          car=molcrys%atm(i)%chemsymb
          car=u_case(car)
          do j=1,Num_Chem_Info
             if (car == Chem_Info(j)%Symb) then
                n_pt(j)=n_pt(j)+1
                exit
             end if
          end do
       end do

       do k=1,molcrys%n_mol
          do i=1,molcrys%mol(k)%natoms
                 car=molcrys%mol(k)%atsymb(i)
             car=u_case(car)
             do j=1,Num_Chem_Info
                if (car == Chem_Info(j)%Symb) then
                   n_pt(j)=n_pt(j)+1
                   exit
                end if
             end do
          end do
       end do

       if (all (n_pt ==0)) then
          if (present(Form_weight)) Form_weight=0.0
          call Remove_Chem_Info()
          return
       end if

       do i=1,Num_Chem_Info
          if (n_pt(i) == 0) cycle
          car=Chem_Info(i)%Symb
          car(2:2)=l_case(car(2:2))
          write(unit=numcar,fmt="(i5)") n_pt(i)
          Formula=trim(Formula)//trim(car)//adjustl(numcar)
          weight=weight+n_pt(i)*Chem_Info(i)%atwe
       end do

       call Remove_Chem_Info()

       if (present(Form_weight)) Form_weight=weight

       return
    End Subroutine Empiric_Formula_Molcrys

    !!--++
    !!--++ Subroutine Empiric_Formula_Molec(Molecule,Formula,Form_Weight)
    !!--++    type(molecule_type),     intent(in)  :: Molecule
    !!--++    character(len=*),        intent(out) :: Formula
    !!--++    real(kind=cp), optional, intent(out) :: Form_Weight
    !!--++
    !!--++    (Overloaded)
    !!--++    Obtain the Empiric Formula from Molecule variable and
    !!--++    the Weight is the variable is present.
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Empiric_Formula_Molec(Molecule,Formula,Form_Weight)
       !---- Arguments ----!
       type(molecule_type),      intent(in)  :: Molecule
       character(len=*),         intent(out) :: Formula
       real(kind=cp), optional,  intent(out) :: Form_Weight

       !---- Local variables ----!
       character(len=2)                  :: car
       character(len=5)                  :: numcar
       integer                           :: i,j
       integer, dimension(Num_Chem_Info) :: N_PT
       real(kind=cp)                     :: weight

       !---- Init ----!
       N_PT=0
       weight=0.0

       Formula=" "
       if (molecule%natoms <= 0) then
          if (present(Form_weight)) Form_weight=0.0
          return
       end if

       !---- Set Information Table ----!
       call Set_Chem_Info()

       do i=1,molecule%natoms
          call Get_ChemSymb(molecule%atsymb(i),car)
          car=u_case(car)
          do j=1,Num_Chem_Info
             if (car == Chem_Info(j)%Symb) then
                n_pt(j)=n_pt(j)+1
                exit
             end if
          end do
       end do

       if (all (n_pt ==0)) then
          if (present(Form_weight)) Form_weight=0.0
          call Remove_Chem_Info()
          return
       end if

       do i=1,Num_Chem_Info
          if (n_pt(i) == 0) cycle
          car=Chem_Info(i)%Symb
          car(2:2)=l_case(car(2:2))
          write(unit=numcar,fmt="(i5)") n_pt(i)
          Formula=trim(Formula)//trim(car)//adjustl(numcar)
          weight=weight+n_pt(i)*Chem_Info(i)%atwe
       end do

       call Remove_Chem_Info()

       if (present(Form_weight)) Form_weight=weight

       return
    End Subroutine Empiric_Formula_Molec

    !!----
    !!---- Subroutine Fractional_to_Cartesian(Molecule,Cell,NewMolecule)
    !!----    type (Molecule_type),     intent(in out)           :: Molecule
    !!----    type (Crystal_Cell_Type), intent(in    )           :: Cell
    !!----    type (Molecule_type),     intent(   out), optional :: Newmolecule
    !!----
    !!----    Subroutine to transform the fractional coordinates to cartesian internal
    !!----    coordinates of a molecule.
    !!----    If Newmolecule is present the subroutine creates a new molecule
    !!----    (copy of the old one) with cartesian coordinates, preserving
    !!----    the input molecule in fractional. Otherwise the input molecule is
    !!----    changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Fractional_to_Cartesian(Molecule,Cell,NewMolecule)
       !---- Arguments ----!
       type (Molecule_type),     intent(in out)           :: Molecule
       type (Crystal_Cell_Type), intent(in    )           :: Cell
       type (Molecule_type),     intent(   out), optional :: NewMolecule

       !---- Local variables -----!
       integer                       :: i, na
       real(kind=cp)                 :: phi,theta,chi
       real(kind=cp), dimension(3)   :: ci,xi
       real(kind=cp), dimension(3,3) :: Eu

       type (Molecule_type)          :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "F") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Fractional_to_Cartesian: the input molecule is not in fractional coordinates"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Fractional_to_Cartesian: No atoms are defined"
          return
       end if

       call Init_molecule(NewMol,na)
       NewMol=Molecule

       if (molecule%in_xtal) then
          if (newmol%is_EulerMat) then
             Eu=newmol%Euler
          else
             phi=newmol%orient(1)
             theta=newmol%orient(2)
             chi=newmol%orient(3)
             Call Set_Euler_matrix(newmol%rot_type,phi,theta,chi,Eu)
             newmol%Euler=Eu
             newmol%is_EulerMat=.true.
          end if

          !---- Newmol contains fractional coordinates ----!
          do i=1,newmol%natoms
             xi=newmol%I_coor(:,i) - newmol%xcentre !Fractional coordinates after removing translation
             ci=matmul(cell%Cr_Orth_cel,xi)       !Cartesian components in the Crystal Frame
             newmol%I_coor(1:3,i) = matmul(ci,Eu)   !Final Cartesian internal coordinates (use passive matrix!)
          end do
       else
          do i=1,newmol%natoms
             newmol%I_coor(:,i)=matmul(cell%cr_orth_cel,newmol%I_coor(:,i))
          end do
          call Fix_Orient_Cartesian(newmol)  ! Select the internal frame as needed for Z-matrices
       end if
       newmol%coor_type = "C"

       if (present(NewMolecule)) then
          call Init_molecule(NewMolecule,Newmol%natoms)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Fractional to Cartesian: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Fractional_to_Cartesian

    !!----
    !!---- Subroutine Fractional_to_Spherical(Molecule,Cell,NewMolecule)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Crystal_Cell_Type), intent(in)           :: Cell
    !!----    type (Molecule_type), intent(   out), optional :: Newmolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a
    !!----    molecule from Fractional coordinates to  Spherical coordinaters.
    !!----    If a third argument is present the subroutine creates a new
    !!----    molecule (copy of the old one) with Spherical coordinates,
    !!----    preserving the input molecule in Fractional Coordinates. Otherwise
    !!----    the input molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Fractional_to_Spherical(Molecule, Cell, NewMolecule)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Crystal_Cell_Type), intent(in)           :: Cell
       type (Molecule_type), intent(   out), optional :: NewMolecule

       !---- Local Variables ----!
       integer                     :: na
       type (Molecule_type)        :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "F") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Fractional_to_Spherical: the input molecule is not in Fractional coordinates"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Fractional_to_Spherical: No atoms are defined"
          return
       end if

       !---- Step 1----!
       call init_molecule(NewMol,na)
       NewMol= Molecule
       call Fractional_to_Cartesian(NewMol,Cell)
       if (err_molec) then
          ERR_Molec_Mess="Error in Fractional_to_Spherical: Intermediate procedure fail (I)!"
          return
       end if

       !---- Step 2 ----!
       call Cartesian_to_Spherical(NewMol)
       if (err_molec) then
          ERR_Molec_Mess="Error in Fractional_to_Spherical: Intermediate procedure fail (II)!"
          return
       end if

       !---- Step 3 ----!
       if (present(newmolecule)) then
          call Init_molecule(NewMolecule,Newmol%natoms)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Fractional to Spherical: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Fractional_to_Spherical

    !!----
    !!---- Subroutine Fractional_to_Zmatrix(Molecule,Cell,NewMolecule)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Crystal_Cell_Type), intent(in)           :: Cell
    !!----    type (Molecule_type), intent(   out), optional :: Newmolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a
    !!----    molecule from Fractional coordinates to  Zmatrix coordinaters.
    !!----    If a second argument is present the subroutine creates a new
    !!----    molecule (copy of the old one) with Zmatrix coordinates,
    !!----    preserving the input molecule in Fractional Coordinates. Otherwise
    !!----    the input molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Fractional_to_Zmatrix(Molecule,Cell,NewMolecule)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Crystal_Cell_Type), intent(in)           :: Cell
       type (Molecule_type), intent(   out), optional :: NewMolecule

       !---- Local Variables ----!
       integer                     :: na
       type (Molecule_type)        :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "F") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Fractional_to_Zmatrix: the input molecule is not in Fractional coordinates"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Fractional_to_Spherical: No atoms are defined"
          return
       end if

       !---- Step 1----!
       call Init_Molecule(NewMol,na)
       NewMol=Molecule
       call Fractional_to_Cartesian(NewMol,Cell)
       if (err_molec) then
          ERR_Molec_Mess="Error in Fractional_to_Zmatrix: Intermediate procedure fail (I)!"
          return
       end if

       !---- Step 2 ----!
       call Cartesian_to_Zmatrix(NewMol, Cell=Cell)  !The cell is needed to eventually take into account
       if (err_molec) then                           !a different Cartesian frame on the input molecule
          ERR_Molec_Mess="Error in Fractional_to_Zmatrix: Intermediate procedure fail (II)!"
          return
       end if

       !---- Step 3 ----!
       if (present(newmolecule)) then
          call Init_molecule(NewMolecule,na)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Fractional to ZMatrix: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Fractional_to_Zmatrix

    !!--++
    !!--++ Subroutine Get_Cartesian_from_Z(ci,ri,rj,rk,rn)
    !!--++    real, dimension(3), intent ( in) :: ci,rj,rj,rn
    !!--++    real, dimension(3), intent (out) :: ri
    !!--++
    !!--++    Subroutine to calculate the cartesian coordinates of an atom (i)
    !!--++    when its distance (dij=ci(1)) to another atom (j), the angle (aijk=ci(2))
    !!--++    spanned with another atom (k) centred at (j), the torsion angle
    !!--++    (bijkn=ci(3)) with a fourth atom (n) and the coordinates of
    !!--++    the three atoms (jkn), rj,rk,rn are all given.
    !!--++
    !!--<<    The algorithm used to determine the Cartesian coordinates of atom (i) is the
    !!--++    following:
    !!--++       - Select a local Cartesian frame with (j) at origin, x-axis along (jk),
    !!--++         z-axis perpendicular to the plane (jkn), y-axis right-handled frame
    !!--++            e1 = rjk/djk, e2 = e3 x e1,  e3= rjk x rkn / djk/dkn
    !!--++       - The above system determine a matrix M = (e1,e2,e3), with components ei in columns
    !!--++         that serves to transform interatomic vector components back to the original system.
    !!--++       - In the above system the coordinates of atom (i) is given by
    !!--++            ri = rj + M ui
    !!--++
    !!--++         where
    !!--++            ui = d ( cos(aijk), cos(bijkn) sin(aijk), sqrt(1 - cos(aijk)^2 -(cos(bijkn) sin(aijk))^2))
    !!-->>
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Cartesian_from_Z(ci,ri,rj,rk,rn)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent ( in) :: ci,rj,rk,rn
       real(kind=cp), dimension(3), intent (out) :: ri

       !--- Local variables ---!
       real(kind=cp)                 :: ca,cb,sa
       real(kind=cp), dimension(3)   :: r,e1,e2,e3
       real(kind=cp), dimension(3,3) :: M

       ca = cosd(ci(2))                  ! cos(aijk)
       sa = sqrt(abs(1.0_cp - ca*ca))    ! sin(aijk)
       cb = cosd(ci(3))                  ! cos(bijkn)
       r(1) = ci(1) * ca                 ! Coordinates in the local system
       r(2) = ci(1)*cb*sa
       r(3) = ci(1)*sqrt(abs(1.0_cp - ca*ca - sa*sa*cb*cb )) *sign(1.0_cp,ci(3))

       e1  = rk - rj
       e1  = e1/sqrt(dot_product(e1,e1))
       e3  = cross_product( rk - rj, rn - rk)
       e3  = e3/sqrt(dot_product(e3,e3))
       e2  = cross_product( e3, e1)
       M(:,1) = e1
       M(:,2) = e2
       M(:,3) = e3

       ri = rj + matmul(M,r)

       return
    End Subroutine Get_Cartesian_from_Z


    !!--++
    !!--++ Subroutine Get_Z_from_Cartesian(ci,ri,rj,rk,rn)
    !!--++    real, dimension(3), intent ( in) :: ri,rj,rj,rn
    !!--++    real, dimension(3), intent (out) :: ci
    !!--++
    !!--++     Subroutine to calculate the distance of an atom (i)
    !!--++     (dij=ci(1)) to another atom (j), the angle (aijk=ci(2))
    !!--++     spanned with another atom (k) centred at (j) and  the torsion angle
    !!--++     (bijkn=ci(3)) with a fourth atom (n) when the cartesian coordinates are given
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Get_Z_from_Cartesian(ci,ri,rj,rk,rn)
       !---- Arguments ----!
       real(kind=cp), dimension(3), intent ( in) :: ri,rj,rk,rn
       real(kind=cp), dimension(3), intent (out) :: ci

       !--- Local variables ---!
       real(kind=cp)                 :: dji,djk
       real(kind=cp), dimension(3)   :: rji,rjk

       rji = ri-rj
       ci(1) = sqrt(dot_product(rji,rji))
       rjk = rk-rj
       dji = ci(1)
       djk = sqrt(dot_product(rjk,rjk))
       ci(2) = acosd( dot_product(rji,rjk)/dji/djk)
       ci(3) = angle_dihedral(ri,rj,rk,rn)
       if (abs(ci(3)+180.00) <= 0.001) ci(3)=180.0

       return
    End Subroutine Get_Z_from_Cartesian

    !!----
    !!---- Subroutine Init_Err_Molec()
    !!----
    !!----    Initialize Flags of Errors in this module
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Err_Molec()

       err_molec=.false.
       ERR_Molec_Mess=" "

       return
    End Subroutine Init_Err_Molec

    !!----
    !!---- Subroutine Init_Molecule(Molecule,Natm)
    !!----    type(Molecule_Type), intent(out) :: Molecule
    !!----    integer, optional,   intent(in)  :: Natm
    !!----
    !!----    Initialize the Variable Molecule
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Init_Molecule(Molecule,Natm)
       !---- Argument ----!
       type(Molecule_Type), intent(out) :: Molecule
       integer, optional,   intent(in)  :: Natm

       molecule%name_mol   =" "
       molecule%natoms     =0

       molecule%in_xtal    = .false.
       molecule%is_eulerMat= .false.
       molecule%is_connect = .false.
       molecule%rot_type   =" "
       molecule%coor_type  =" "
       molecule%therm_type =" "

       molecule%xcentre    =0.0
       molecule%mxcentre   =0.0
       molecule%lxcentre   =0

       molecule%orient     =0.0
       molecule%mOrient    =0.0
       molecule%lorient    =0

       molecule%t_tls      =0.0
       molecule%mT_TLS     =0.0
       molecule%lt_tls     =0

       molecule%l_tls      =0.0
       molecule%mL_TLS     =0.0
       molecule%ll_tls     =0

       molecule%s_tls      =0.0
       molecule%mS_TLS     =0.0
       molecule%ls_tls     =0

       molecule%euler      =0.0

       if (allocated(molecule%AtName))  deallocate(molecule%AtName)
       if (allocated(molecule%AtSymb))  deallocate(molecule%AtSymb)
       if (allocated(molecule%AtZ))     deallocate(molecule%AtZ)
       if (allocated(molecule%Ptr))     deallocate(molecule%Ptr)
       if (allocated(molecule%I_Coor))  deallocate(molecule%I_Coor)
       if (allocated(molecule%mI_Coor)) deallocate(molecule%mI_Coor)
       if (allocated(molecule%lI_Coor)) deallocate(molecule%lI_Coor)
       if (allocated(molecule%Biso))    deallocate(molecule%Biso)
       if (allocated(molecule%mbiso))   deallocate(molecule%mbiso)
       if (allocated(molecule%lBiso))   deallocate(molecule%lBiso)
       if (allocated(molecule%Occ))     deallocate(molecule%Occ)
       if (allocated(molecule%mocc))    deallocate(molecule%mocc)
       if (allocated(molecule%lOcc))    deallocate(molecule%lOcc)
       if (allocated(molecule%Nb))      deallocate(molecule%Nb)
       if (allocated(molecule%INb))     deallocate(molecule%INb)
       if (allocated(molecule%Tb))      deallocate(molecule%Tb)
       if (allocated(molecule%Conn))    deallocate(molecule%Conn)

       if (present(natm)) then
          if (natm > 0) then
             molecule%natoms=natm

             allocate(molecule%AtName(natm))
             allocate(molecule%AtSymb(natm))
             allocate(molecule%AtZ(natm))
             allocate(molecule%Ptr(2,natm))
             allocate(molecule%I_Coor(3,natm))
             allocate(molecule%mI_Coor(3,natm))
             allocate(molecule%lI_Coor(3,natm))
             allocate(molecule%Biso(natm))
             allocate(molecule%mbiso(natm))
             allocate(molecule%lBiso(natm))
             allocate(molecule%Occ(natm))
             allocate(molecule%mocc(natm))
             allocate(molecule%lOcc(natm))
             allocate(molecule%Nb(natm))
             allocate(molecule%INb(10,natm))
             allocate(molecule%Tb(10,natm))
             allocate(molecule%Conn(3,natm))

             molecule%AtName  =" "
             molecule%AtSymb  =" "
             molecule%AtZ     =0
             molecule%Ptr     =0
             molecule%I_Coor  =0.0
             molecule%mI_Coor =0.0
             molecule%lI_Coor =0
             molecule%Biso    =0.0
             molecule%mbiso   =0.0
             molecule%lBiso   =0
             molecule%Occ     =0.0
             molecule%mocc    =0.0
             molecule%lOcc    =0
             molecule%Nb      =0
             molecule%INb     =0
             molecule%Tb      =0
             molecule%Conn    =0


          end if
       end if

       return
    End Subroutine Init_Molecule
    !!----
    !!---- Subroutine Init_Molecule(Molx,Natm,Nmol)
    !!----    type(Molecular_Crystal_Type), intent(out) :: Molx
    !!----    integer, optional,            intent(in)  :: Natm
    !!----    integer, optional,            intent(in)  :: Nmol
    !!----
    !!----    Initialization for Molecular Crystal
    !!----
    !!---- Update: October - 2014
    !!
    Subroutine Init_Mol_Crys(Molx,Natm,Nmol)
        !---- Argument ----!
        type(Molecular_Crystal_Type), intent(out) :: Molx
        integer, optional,            intent(in)  :: Natm
        integer, optional,            intent(in)  :: Nmol

        integer :: i

        molx%N_Free    = 0
        molx%N_Mol     = 0
        molx%N_Species = 0
        molx%Npat      = 0

        if (allocated(molx%atm))  deallocate(molx%atm)
        if (allocated(molx%mol))  deallocate(molx%mol)

        if (present(nmol) .and. nmol > 0) then
            molx%N_Mol = nmol
            allocate(molx%mol(nmol))
            do i=1,nmol
                call init_molecule(molx%mol(i))
            end do
        end if

        if (present(natm) .and. natm > 0) then
            molx%N_Free = natm
            allocate (molx%atm(natm))
            do i=1,natm
                call init_atom_type(molx%atm(i))
            end do
        end if
        return
    End Subroutine Init_Mol_Crys
    !!----
    !!---- Subroutine Molcrys_to_AtomList(Molcrys,Atm)
    !!----    type (Molecular_Crystal_Type), intent(in)  :: Molec
    !!----    type (Atom_List_Type),         intent(out) :: Atm
    !!----
    !!---- Subroutine to pass all information from Molecular_Crystal_Type
    !!---- to Atom_List_Type
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Molcrys_to_AtomList(Molcrys,Atm)
       !---- Arguments ----!
       type (Molecular_Crystal_Type), intent(in)  :: Molcrys
       type (Atom_List_Type),         intent(out) :: Atm

       !---- Local variables ----!
       integer               :: i, n
       integer               :: Nat, NaF, NMol
       type (Atom_List_Type) :: A

       !---- Number of Atoms ----!
       NaF=molcrys%n_free
       NMol=molcrys%n_mol
       if (NMol > 0) then
          Nat=NaF+sum(molcrys%mol(1:NMol)%natoms)
       else
          Nat=NaF
       end if
       if (Nat <= 0) return

       !---- Allocating Atom_List_Type ----!
       call allocate_atom_list(Nat,Atm)

       !---- Fill information from Free atoms Part ----!
       if (naF > 0) Atm%Atom(1:NaF)=molcrys%atm(1:NaF)

       !---- Fill information from Molecules Part ----!
       n=naF
       do i=1,NMol
          if (molcrys%mol(i)%natoms <= 0) cycle
          if (.not. molcrys%mol(i)%in_xtal) cycle
          call molec_to_AtomList(molcrys%mol(i),A,"F",molcrys%cell)
          if (err_molec) return
          if (A%natoms <= 0) cycle
          Atm%Atom(n+1:n+A%natoms)=A%Atom(1:A%natoms)
          n=n+A%natoms
          call deallocate_atom_list(A)
       end do

       return
    End Subroutine Molcrys_to_AtomList

    !!----
    !!---- Subroutine Molec_to_AtomList(Molec,Atm, Coor_Type, Cell)
    !!----    type (Molecule_Type),               intent(in)  :: Molec
    !!----    type (Atom_List_Type),              intent(out) :: Atm
    !!----    character(len=*),         optional, intent(in)  :: Coor_type
    !!----    type (Crystal_Cell_type), optional, intent(in)  :: Cell
    !!----
    !!---- Subroutine to pass all information from Molecule_Type
    !!---- to Atom_List_Type. Coor_type determine the type of
    !!---- cordinates parameter in output. In general Cell if
    !!---- necessary to obtain on Output fractional coordinates or
    !!---- special case for ZMatrix.
    !!----
    !!---- Update: April - 2005
    !!
    Subroutine Molec_to_AtomList(Molec, Atm, Coor_Type, Cell)
       !---- Arguments ----!
       type (Molecule_Type),               intent(in)     :: Molec
       type (Atom_List_Type),              intent(out)    :: Atm
       character(len=*),         optional, intent(in)     :: Coor_type
       type (Crystal_Cell_type), optional, intent(in)     :: Cell

       !---- Local Variables ----!
       character(len=1)      :: car
       integer               :: i,nat
       type (Molecule_Type)  :: newmol

       !---- Number of Atoms ----!
       Nat=molec%natoms
       Atm%natoms=0
       if (Nat <= 0) return

       car="F"
       if (present(coor_type)) car=adjustl(coor_type)
       call init_molecule(newmol,nat)

       newmol=molec

       select case (car)
          case ("C")
             select case (molec%coor_type)
                case ("C")

                case ("F")
                   if (present(cell)) then
                      call Fractional_to_Cartesian(newmol,cell)
                      if (err_molec) then
                         call init_molecule(newmol)
                         return
                      end if
                   else
                      err_molec=.true.
                      ERR_Molec_Mess="You need the Cell_Type on this routine"
                      call init_molecule(newmol)
                      return
                   end if

                case ("S")
                   call Spherical_to_Cartesian(newmol)
                   if (err_molec) then
                      call init_molecule(newmol)
                      return
                   end if

                case ("Z")
                   call ZMatrix_to_Cartesian(newmol)
                   if (err_molec) then
                      call init_molecule(newmol)
                      return
                   end if

             end select

          case ("F")
             select case (molec%coor_type)
                case ("C")
                   if ( present(cell)) then
                      call Cartesian_to_Fractional(newmol,cell)
                      if (err_molec) then
                         call init_molecule(newmol)
                         return
                      end if
                   else
                      err_molec=.true.
                      ERR_Molec_Mess="You need the Cell_Type on this routine"
                      call init_molecule(newmol)
                      return
                   end if

                case ("F")

                case ("S")
                   if (present(cell)) then
                      call Spherical_to_Fractional(newmol,cell)
                      if (err_molec) then
                         call init_molecule(newmol)
                         return
                      end if
                   else
                      err_molec=.true.
                      ERR_Molec_Mess="You need the Cell_Type on this routine"
                      call init_molecule(newmol)
                      return
                   end if

                case ("Z")
                   if (present(cell)) then
                      call ZMatrix_to_Fractional(newmol,cell)
                      if (err_molec) then
                         call init_molecule(newmol)
                         return
                      end if
                   else
                      err_molec=.true.
                      ERR_Molec_Mess="You need the Cell_Type on this routine"
                      call init_molecule(newmol)
                      return
                   end if
             end select

          case ("S")
             select case (molec%coor_type)
                case ("C")
                   call Cartesian_to_Spherical(newmol)
                   if (err_molec) then
                      call init_molecule(newmol)
                      return
                   end if

                case ("F")
                   if (present(cell)) then
                      call Fractional_to_Spherical(newmol,cell)
                      if (err_molec) then
                         call init_molecule(newmol)
                         return
                      end if
                   else
                      err_molec=.true.
                      ERR_Molec_Mess="You need the Cell_Type on this routine"
                      call init_molecule(newmol)
                      return
                   end if

                case ("S")

                case ("Z")
                   call ZMatrix_to_Spherical(newmol)
                   if (err_molec) then
                      call init_molecule(newmol)
                      return
                   end if

             end select

          case ("Z")
             select case (molec%coor_type)
                case ("C")
                   if (present(cell)) then
                      call Cartesian_to_ZMatrix(newmol,cell=cell)
                   else
                      call Cartesian_to_ZMatrix(newmol)
                   end if
                   if (err_molec) then
                      call init_molecule(newmol)
                      return
                   end if

                case ("F")
                   if (present(cell)) then
                      call Fractional_to_ZMatrix(newmol,cell)
                      if (err_molec) then
                         call init_molecule(newmol)
                         return
                      end if
                   else
                      err_molec=.true.
                      ERR_Molec_Mess="You need the Cell_Type on this routine"
                      call init_molecule(newmol)
                      return
                   end if

                case ("S")
                   if (present(cell)) then
                      call Spherical_to_ZMatrix(newmol,cell=cell)
                   else
                      call Spherical_to_ZMatrix(newmol)
                   end if
                   if (err_molec) then
                      call init_molecule(newmol)
                      return
                   end if

                case ("Z")
             end select

       end select

       !---- Allocating Atom_List_Type ----!
       call allocate_atom_list(Nat,Atm)

       !---- Passing Information ----!
       Atm%Atom(1:Nat)%Lab      =Newmol%AtName(1:Nat)
       Atm%Atom(1:Nat)%SfacSymb =Newmol%AtSymb(1:Nat)
       Atm%Atom(1:Nat)%Active   =.true.
       Atm%Atom(1:Nat)%Z        =Newmol%AtZ(1:Nat)
       Atm%Atom(1:Nat)%Mult     =1
       !Atm%Atom(1:Nat)%X        =Newmol%I_Coor(:,1:Nat)
       !Atm%Atom(1:Nat)%X_Std    =0.0
       !Atm%Atom(1:Nat)%MX       =Newmol%mI_Coor(:,1:Nat)
       !Atm%Atom(1:Nat)%LX       =Newmol%lI_Coor(:,1:Nat)
       Atm%Atom(1:Nat)%Occ      =Newmol%Occ(1:Nat)
       Atm%Atom(1:Nat)%Occ_Std  =0.0
       Atm%Atom(1:Nat)%MOcc     =Newmol%mOcc(1:Nat)
       Atm%Atom(1:Nat)%LOcc     =Newmol%lOcc(1:Nat)
       Atm%Atom(1:Nat)%Biso     =Newmol%biso(1:Nat)
       Atm%Atom(1:Nat)%Biso_std =0.0
       Atm%Atom(1:Nat)%MBiso    =Newmol%mbiso(1:Nat)
       Atm%Atom(1:Nat)%LBiso    =Newmol%lbiso(1:Nat)
       Atm%Atom(1:Nat)%Utype    ="none"
       Atm%Atom(1:Nat)%ThType   ="isotr"
       !Atm%Atom(1:Nat)%U        =0.0
       !Atm%Atom(1:Nat)%U_std    =0.0
       Atm%Atom(1:Nat)%Ueq      =0.0
       !Atm%Atom(1:Nat)%MU       =0.0
       !Atm%Atom(1:Nat)%LU       =0
       Atm%Atom(1:Nat)%Charge   =0.0
       Atm%Atom(1:Nat)%Moment   =0.0
       !Atm%Atom(1:Nat)%Ind      =0
       Atm%Atom(1:Nat)%NVar     =0
       !Atm%Atom(1:Nat)%VarF     =0.0

       do i=1,Nat
          call Get_ChemSymb(Atm%Atom(i)%SfacSymb, Atm%Atom(i)%ChemSymb)
          Atm%Atom(i)%X=Newmol%I_Coor(:,i)
          Atm%Atom(i)%X_Std=0.0
          Atm%Atom(i)%mX=Newmol%mI_Coor(:,i)
          Atm%Atom(i)%lX=Newmol%lI_Coor(:,i)
          Atm%Atom(i)%U    =0.0
          Atm%Atom(i)%U_Std=0.0
          Atm%Atom(i)%mU   =0.0
          Atm%Atom(i)%lU   =0
          Atm%Atom(i)%Ind  =0
          Atm%Atom(i)%VarF =0.0
       end do

       call init_molecule(newmol)

       return
    End Subroutine Molec_to_AtomList

    !!----
    !!---- Subroutine Read_Free_Atoms(Lun,AtmF,N)
    !!----    integer,                       intent(in)   :: Lun        ! Logical unit to be rad
    !!----    type(Atom_Type), dimension(:), intent(out)  :: AtmF       ! Free atoms
    !!----    integer,                       intent(out)  :: N          ! Free atoms read
    !!----
    !!--<<    Subroutine to read a set of Free Atoms from a file.
    !!----    The format is:
    !!----        ATOMS N_Atoms
    !!----
    !!----    Internal Coordinates for Atoms (N_Atoms Lines)
    !!----        Atom_Name(6)  Atom_Specie(4)  Coordinates(3)  Biso  Occ [VARY]
    !!----
    !!----    if VARY is present as last option on the Internal Coordinates line,
    !!----    then an extra line is read
    !!----        Codes_Coordinates(3)   Code_BIso  Code_Occ
    !!-->>
    !!----
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Read_Free_Atoms(Lun, AtmF, N)
       !---- Arguments ----!
       integer,                       intent(in)   :: Lun    ! Logical unit to be rad
       type(Atom_Type), dimension(:), intent(out)  :: AtmF   ! Free atoms
       integer,                       intent(out)  :: N      ! Free atoms read

       !---- Local Variables ----!
       character(len=80)           :: line
       character(len=6)            :: label
       character(len=4)            :: var,symb
       integer                     :: i,ier,nlong,iv
       integer,       dimension(5) :: ivet
       real(kind=cp), dimension(5) :: vet

       call Init_Err_Molec()
       N=0
       do
          read(unit=lun,fmt="(a)",iostat=ier) line
          if (ier == 0) then
             line=adjustl(line)
             if (u_case(line(1:4)) /= "ATOM") cycle
          else
             err_molec=.true.
             ERR_Molec_Mess="Atoms Information not found in file! "
             return
          end if

          call cutst(line,nlong)
          call getnum(line,vet,ivet,iv)
          if (iv /= 1) then
             err_molec=.true.
             ERR_Molec_Mess="Number of Free atoms not found in file! "
             return
          end if
          N=ivet(1)
          exit
       end do

       do i=1,N
          read(unit=lun,fmt="(a)",iostat=ier) line
          if (ier /=0) then
             err_molec=.true.
             ERR_Molec_Mess="Free atoms Information was incomplete "
             return
          end if
          call cutst(line,nlong,label)
          call cutst(line,nlong,symb)

          line=u_case(line)
          var=" "
          iv=index(line,"VARY")
          if (iv /= 0) then
             line=line(1:iv-1)
             var="VARY"
          end if

          call getnum(line,vet,ivet,iv)
          select case (iv)
             case (:2)
                vet(1:3)=0.0
                vet(4)=0.0
                vet(5)=1.0
             case (3)
                vet(4)=0.0
                vet(5)=1.0
             case (4)
                vet(5)=1.0
          end select
          AtmF(i)%Lab =label
          AtmF(i)%ChemSymb=symb
          AtmF(i)%x=vet(1:3)
          AtmF(i)%biso=vet(4)
          AtmF(i)%occ =vet(5)

          if (var == "VARY") then
             do
                read(unit=lun,fmt="(a)", iostat=ier) line
                if (ier /= 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading the refinement codes of free atoms "
                   return
                end if
                line=adjustl(line)
                if (line(1:1) =="!") cycle
                exit
             end do

             call getnum(line,vet,ivet,iv)
             select case (iv)
                case (3)
                   AtmF(i)%mx =vet(1:3)

                case (5)
                   AtmF(i)%mx    =vet(1:3)
                   AtmF(i)%mbiso =vet(4)
                   AtmF(i)%mocc  =vet(5)

                case default
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading the refinement codes of free atoms  "
                   return
             end select
          end if

       end do

       return
    End Subroutine Read_Free_Atoms

    !!----
    !!---- Subroutine Read_Molecule(Lun,Molecule) or (File_Dat, N_Ini, N_End, Molecule)
    !!----    integer,              intent( in)           :: Lun         !  In -> Logical unit to be read
    !!----    or
    !!----    character(len=*), dimension(:), intent(in)  :: File_Dat
    !!----    integer,                        intent(in)  :: N_Ini
    !!----    integer,                        intent(in)  :: N_End
    !!----    and
    !!----    type (Molecule_type),           intent(out) :: Molecule    ! Out -> Molecule
    !!----
    !!--<<    Subroutine to read a molecule from a file
    !!----    The format is:
    !!----        MOLE[X] N_Atoms Molecule_Name Coordinates_Type
    !!----
    !!----    where:
    !!----        N_atoms             Number of atoms in the molecule definition
    !!----        Molecule_Name       Name for the molecule
    !!----        Coordinates_Type    C: Cartesian coordinates
    !!----                            F: Fractional coordinates
    !!----                            S: Spherical coordinates
    !!----                            Z: Z-Matrix coordinates
    !!----
    !!----    If keyword MOLEX is present, then the next line will be read (6 reals, 2 characters)
    !!----        Molecule_Centre(3), Molecule_Orient(3), Rotational_Angle Type(1), Thermal_Factor Type(1)
    !!----
    !!----    where:
    !!----        Molecule_Centre     Coordinate of Center of Molecule
    !!----        Molecule_Orient     Angles orientation
    !!----        Rotational Angle    E: Conventional Euler angles (alpha, beta, gamma)
    !!----                            P: Polar Euler angles (Phi, theta, Chi) (default)
    !!----        Thermal Factor    ISO: No collective motion
    !!----                          TLS: Traslational + Librational + Correlation
    !!----                           TL: Traslational + Librational
    !!----                            T: Traslational
    !!----
    !!----        According to Thermal Factors, next lines will be read
    !!----                          [T]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!----
    !!----                         [TL]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!----                               6 Thermal Factors (Line3) + 6 Codes Thermal Factors (Line4)
    !!----
    !!----                        [TLS]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!----                               6 Thermal Factors (Line3) + 6 Codes Thermal Factors (Line4)
    !!----                               9 Thermal Factors (Line5) + 9 Codes Thermal Factors (Line6)
    !!----
    !!----    Internal Coordinates for Atoms (N_Atoms Lines)
    !!----        Atom_Name(6)  Atom_Specie(4)  Coordinates(3)  N1  N2  N3  Biso  Occ [VARY]
    !!----
    !!----    if VARY is present as last option on the Internal Coordinates line,
    !!----    then an extra line is read
    !!----        Codes_Coordinates(3)   Code_BIso  Code_Occ
    !!-->>
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!

    !!--++
    !!--++ Subroutine Read_Molecule_in_File(Lun,Molecule)
    !!--++    integer,              intent( in)                  :: Lun         !  In -> Logical unit to be read
    !!--++    type (Molecule_type), intent(out)                  :: Molecule    ! Out -> Molecule
    !!--++
    !!--++    (Overloaded)
    !!--++    Subroutine to read a molecule from a file.
    !!--++    The format is:
    !!--++        MOLE[X] N_Atoms Molecule_Name Coordinates_Type
    !!--++
    !!--++    where:
    !!--++        N_atoms             Number of atoms in the molecule definition
    !!--++        Molecule_Name       Name for the molecule
    !!--++        Coordinates_Type    C: Cartesian coordinates
    !!--++                            F: Fractional coordinates
    !!--++                            S: Spherical coordinates
    !!--++                            Z: Z-Matrix coordinates
    !!--++
    !!--++    If keyword MOLEX is present, then the next line will be read (6 reals, 2 characters)
    !!--++        Molecule_Centre(3), Molecule_Orient(3), Rotational_Angle Type(1), Thermal_Factor Type(1)
    !!--++
    !!--++    where:
    !!--++        Molecule_Centre     Coordinate of Center of Molecule
    !!--++        Molecule_Orient     Angles orientation
    !!--++        Rotational Angle    E: Conventional Euler angles (alpha, beta, gamma)
    !!--++                            P: Polar Euler angles (Phi, theta, Chi) (default)
    !!--++        Thermal Factor    ISO: No collective motion
    !!--++                          TLS: Traslational + Librational + Correlation
    !!--++                           TL: Traslational + Librational
    !!--++                            T: Traslational
    !!--++
    !!--++        According to Thermal Factors, next lines will be read
    !!--++                          [T]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!--++
    !!--++                         [TL]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!--++                               6 Thermal Factors (Line3) + 6 Codes Thermal Factors (Line4)
    !!--++
    !!--++                        [TLS]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!--++                               6 Thermal Factors (Line3) + 6 Codes Thermal Factors (Line4)
    !!--++                               9 Thermal Factors (Line5) + 9 Codes Thermal Factors (Line6)
    !!--++
    !!--++    Internal Coordinates for Atoms (N_Atoms Lines)
    !!--++        Atom_Name(6)  Atom_Specie(4)  Coordinates(3)  N1  N2  N3  Biso  Occ [VARY]
    !!--++
    !!--++    if VARY is present as last option on the Internal Coordinates line,
    !!--++    then an extra line is read
    !!--++        Codes_Coordinates(3)   Code_BIso  Code_Occ
    !!--++
    !!--++    Control of error is present
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Molecule_in_File(Lun,Molecule)
       !---- Arguments ----!
       integer,              intent(in    ) :: lun
       type (Molecule_type), intent(   out) :: Molecule

       !---- Local variables -----!
       character(len=150)              :: line
       character(len=20),dimension(10) :: dire
       character(len=4)                :: var
       integer                         :: i,j,ic,npos,na,ier
       integer,dimension(10)           :: ivet
       real(kind=cp), dimension(10)    :: vet
       real(kind=cp),dimension(3,3)    :: Eu
       logical                         :: in_xtal,mol_found

       in_xtal=.false.
       mol_found =.false.
       call Init_Err_Molec()

       do
          read(unit=lun,fmt="(a)", iostat=ier) line
          if (ier == 0) then
             line=adjustl(line)
             if (u_case(line(1:4)) /= "MOLE") cycle
          else
             if(.not. mol_found) then
              err_molec=.true.
              ERR_Molec_Mess="Molecule not found in file! "
             end if
             return
          end if

          mol_found =.true.
          if (u_case(line(1:5)) == "MOLEX") in_xtal=.true.
          i=index(line,"!")
          if( i /= 0 ) line=line(1:i-1)

          !---- Coordinates format ----!
          call getword(line,dire,ic)
          if (ic /= 4) then
             err_molec=.true.
             ERR_Molec_Mess="Instruction: MOLE[X] N_Atoms Molecule_Name Coordinates_Type, not found in file! "
             return
          end if

          !---- Name and Number of Atoms in the molecule ----!
          read(unit=dire(2),fmt=*,iostat=ier) na
          if (na > 0) then
             call init_molecule(molecule,na)
             Molecule%Name_mol =dire(3)
          else
             err_molec=.true.
             ERR_Molec_Mess="Error reading the number of atoms in a molecule: "//trim(line)
             return
          end if

          select case (dire(4)(1:1)) ! Coordinates_Type [C,S,F,Z]
             case ("C","c")
                molecule%coor_type="C"
             case ("F","f")
                molecule%coor_type="F"
             case ("S","s")
                molecule%coor_type="S"
             case ("Z","z")
                molecule%coor_type="Z"
             case default
                err_molec=.true.
                ERR_Molec_Mess="Coordinates Type for Molecule Unknown! "
                return
          end select ! dire
          exit !The molecule has been found
       end do

       !---- Initialize the crystal part of the molecule
       Molecule%xcentre    = 0.0_cp
       Molecule%orient     = 0.0_cp
       Molecule%therm_type = "   "
       Molecule%T_TLS      = 0.0_cp
       Molecule%L_TLS      = 0.0_cp
       Molecule%S_TLS      = 0.0_cp
       Molecule%in_xtal    = .false.
       Molecule%is_EulerMat=.false.

       if (in_xtal) then
          !---- Read the global coordinates of the centre of molecule and orientational angles
          do
             read(unit=lun,fmt="(a)", iostat=ier) line
             if (ier /= 0) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading Molecule information! "
                return
             end if
             line=adjustl(line)
             if (line(1:1) =="!") cycle
             exit
          end do

          i=index(line,"!")
          if( i /= 0 ) line=line(1:i-1)

          call getword(line,dire,ic)
          if (ic /= 8) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the position and angles of the molecule: "//trim(Molecule%Name_mol)
             return
          end if

          line=trim(dire(1))//"   "//trim(dire(2))//"   "//trim(dire(3))
          call getnum(line,vet,ivet,ic)
          if (ic /= 3) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the position and angles of the molecule: "//trim(Molecule%Name_mol)
             return
          end if
          Molecule%xcentre=vet(1:3)

          line=trim(dire(4))//"   "//trim(dire(5))//"   "//trim(dire(6))
          call getnum(line,vet,ivet,ic)
          if (ic /= 3) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the position and angles of the molecule: "//trim(Molecule%Name_mol)
             return
          end if
          Molecule%orient=vet(1:3)

          Molecule%rot_type=adjustl(u_case(dire(7)))
          Molecule%therm_type=adjustl(u_case(dire(8)))

          do
             read(unit=lun,fmt="(a)", iostat=ier) line
             if (ier /= 0) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading Molecule information! "
                return
             end if
             line=adjustl(line)
             if (line(1:1) =="!") cycle
             exit
          end do
          i=index(line,"!")
          if( i /= 0 ) line=line(1:i-1)

          call getnum(line,vet,ivet,ic)
          if (ic /= 6) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the position and angles of the molecule: "//trim(Molecule%Name_mol)
             return
          end if
          Molecule%mxcentre=vet(1:3)
          Molecule%mOrient =vet(4:6)

          Molecule%in_xtal = .true.

          !---- Set the Euler Matrix
          if (Molecule%rot_type /= "E") Molecule%rot_type="P"

          call Set_euler_matrix(Molecule%rot_type,  &
                                Molecule%orient(1),Molecule%orient(2),Molecule%orient(3),Eu)
                                !    Phi/alpha          Theta/beta          Chi/gamma
          Molecule%Euler=Eu
          Molecule%is_EulerMat=.true.

          !---- Read the THERMAL PARAMETERS
          if (Molecule%therm_type(1:1) == "T") then
             do
                read(unit=lun,fmt="(a)", iostat=ier) line
                if (ier /= 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(line)
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 6) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the tensor T of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%T_TLS=vet(1:6)

             do
                read(unit=lun,fmt="(a)", iostat=ier) line
                if (ier /= 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(line)
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 6) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the codes of tensor T of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%mT_TLS=vet(1:6)
          end if

          if (Molecule%therm_type(2:2) == "L") then
             do
                read(unit=lun,fmt="(a)", iostat=ier) line
                if (ier /= 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(line)
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 6) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the tensor L of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%L_TLS=vet(1:6)

             do
                read(unit=lun,fmt="(a)", iostat=ier) line
                if (ier /= 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(line)
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 6) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the codes of the tensor L of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%mL_TLS=vet(1:6)
          end if

          if (Molecule%therm_type(3:3) == "S") then
             do
                read(unit=lun,fmt="(a)", iostat=ier) line
                if (ier /= 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(line)
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 9) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the tensor S of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%S_TLS(1,:)=vet(1:3)
             Molecule%S_TLS(2,:)=vet(4:6)
             Molecule%S_TLS(3,:)=vet(7:9)

             do
                read(unit=lun,fmt="(a)", iostat=ier) line
                if (ier /= 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(line)
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 9) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the code of tensor S of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%mS_TLS(1,:)=vet(1:3)
             Molecule%mS_TLS(2,:)=vet(4:6)
             Molecule%mS_TLS(3,:)=vet(7:9)
          end if

       end if  !(in_xtal)

       !---- Read the internal coordinates of the atoms in the molecule
       !---- Read the Z-matrix/Cartesian/spherical/Fractional coordinates of the molecule
       molecule%is_connect=.true.
       do i=1,na
          do
             read(unit=lun,fmt="(a)", iostat=ier) line
             if (ier /= 0) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading Molecule information! "
                return
             end if
             line=adjustl(line)
             if (line(1:1) =="!") cycle
             exit
          end do
          j=index(line,"!")
          if( j /= 0 ) line=line(1:j-1)

          !---- Atom Name ---!
          call Cutst(line,ic,Molecule%Atname(i))

          !---- Atom specie ----!
          call Cutst(line,ic,Molecule%Atsymb(i))

          !---- Passing Codes? ----!
          call getword(line,dire,ic)
          var=adjustl(dire(ic))
          var=u_case(var)
          if (var == "VARY") then
             ic=len_trim(line)
             npos=index(line(1:ic)," ",back=.true.)
             if (npos <=0) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading Molecule information (II)! "
                return
             end if
             line=line(1:npos)
          end if

          !---- Rest of Information ----!
          vet =0.0
          ivet=0
          call getnum(line,vet,ivet,ic)
          select case (ic)
             case (0)
                Molecule%I_Coor(:,i)=0.0
                Molecule%Conn(:,i)  =0
                Molecule%Biso(i)    =0.5
                Molecule%Occ(i)     =1.0
             case (1)
                Molecule%I_Coor(1,i)  =vet(1)
                Molecule%I_Coor(2:3,i)=0.0
                Molecule%conn(:,i)    =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (2)
                Molecule%I_Coor(1:2,i)=vet(1:2)
                Molecule%I_Coor(3,i)  =0.0
                Molecule%conn(:,i)    =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (3)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(:,i)    =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (4)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(1,i)    =ivet(4)
                Molecule%conn(2:3,i)  =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (5)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(1:2,i)  =ivet(4:5)
                Molecule%conn(3,i)    =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (6)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(:,i)    =ivet(4:6)
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (7)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(:,i)    =ivet(4:6)
                Molecule%biso(i)      =vet(7)
                Molecule%Occ(i)       =1.0

             case (8)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(:,i)    =ivet(4:6)
                Molecule%biso(i)      =vet(7)
                Molecule%Occ(i)       =vet(8)

             case default
                err_molec=.true.
                ERR_Molec_Mess="Error reading the atoms in the molecule: "//trim(Molecule%Name_mol)
                return

          end select ! ic

          if (Molecule%coor_type == "Z") then

             if (i == 2 .and. (ivet(4) ==0 .and. ivet(5) ==0 .and. ivet(6) ==0)) then
                Molecule%conn(1,i)=1
             end if
             if(Molecule%I_Coor(3,i) > 180.0) Molecule%I_Coor(3,i) = Molecule%I_Coor(3,i) -360.0
             if(Molecule%I_Coor(3,i) <-180.0) Molecule%I_Coor(3,i) = Molecule%I_Coor(3,i) +360.0

             if (ivet(4) >= i .or. ivet(5) >= i .or. ivet(6) >= i )                err_molec=.true.
             if (i == 3 .and. (ivet(4) == 0 .or. ivet(5) == 0))                    err_molec=.true.
             if (i > 3 .and. (ivet(4) == 0 .or. ivet(5) == 0 .or. ivet(6) == 0))   err_molec=.true.
             if (err_molec) then
                ERR_Molec_Mess = "The Z-matrix connectivity is wrong: "//trim(line)
                return
             end if
          else
             if (ivet(4) >= i .or. ivet(5) >= i .or. ivet(6) >= i )               molecule%is_connect=.false.
             if (i == 3 .and. (ivet(4) == 0 .or. ivet(5) == 0))                   molecule%is_connect=.false.
             if (i > 3 .and. (ivet(4) == 0 .or. ivet(5) == 0 .or. ivet(6) == 0))  molecule%is_connect=.false.
          end if

          Molecule%mI_Coor(:,i)=0.0
          Molecule%mbiso(i)  =0.0
          Molecule%mocc(i)   =0.0

          if (var == "VARY") then
             do
                read(unit=lun,fmt="(a)", iostat=ier) line
                if (ier /= 0) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading the refinement codes of atoms in the molecule: "//trim(Molecule%Name_mol)
                   return
                end if
                line=adjustl(line)
                if (line(1:1) =="!") cycle
                exit
             end do
             j=index(line,"!")
             if( j /= 0 ) line=line(1:j-1)

             call getnum(line,vet,ivet,ic)
             select case (ic)
                case (3)
                   Molecule%mI_Coor(:,i)=vet(1:3)

                case (5)
                   Molecule%mI_Coor(:,i)=vet(1:3)
                   Molecule%mbiso(i)  =vet(4)
                   Molecule%mocc(i)   =vet(5)

                case default
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading the refinement codes of atoms in the molecule: "//trim(Molecule%Name_mol)
                   return
             end select
          end if

       end do

       return
    End Subroutine Read_Molecule_in_File

    !!--++
    !!--++ Subroutine Read_Molecule_in_Var(File_Dat, N_Ini, N_End, Molecule)
    !!--++    character(len=*), dimension(:), intent(in)  :: File_Dat
    !!--++    integer,                        intent(in)  :: N_Ini
    !!--++    integer,                        intent(in)  :: N_End
    !!--++    type (Molecule_type),           intent(out) :: Molecule    ! Out -> Molecule
    !!--++
    !!--++    (Overloaded)
    !!--++    Subroutine to read a molecule from a file.
    !!--++    The format is:
    !!--++        MOLE[X] N_Atoms Molecule_Name Coordinates_Type
    !!--++
    !!--++    where:
    !!--++        N_atoms             Number of atoms in the molecule definition
    !!--++        Molecule_Name       Name for the molecule
    !!--++        Coordinates_Type    C: Cartesian coordinates
    !!--++                            F: Fractional coordinates
    !!--++                            S: Spherical coordinates
    !!--++                            Z: Z-Matrix coordinates
    !!--++
    !!--++    If keyword MOLEX is present, then the next line will be read (6 reals, 2 characters)
    !!--++        Molecule_Centre(3), Molecule_Orient(3), Rotational_Angle Type(1), Thermal_Factor Type(1)
    !!--++
    !!--++    where:
    !!--++        Molecule_Centre     Coordinate of Center of Molecule
    !!--++        Molecule_Orient     Angles orientation
    !!--++        Rotational Angle    E: Conventional Euler angles (alpha, beta, gamma)
    !!--++                            P: Polar Euler angles (Phi, theta, Chi) (default)
    !!--++        Thermal Factor    ISO: No collective motion
    !!--++                          TLS: Traslational + Librational + Correlation
    !!--++                           TL: Traslational + Librational
    !!--++                            T: Traslational
    !!--++
    !!--++        According to Thermal Factors, next lines will be read
    !!--++                          [T]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!--++
    !!--++                         [TL]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!--++                               6 Thermal Factors (Line3) + 6 Codes Thermal Factors (Line4)
    !!--++
    !!--++                        [TLS]: 6 Thermal Factors (Line1) + 6 Codes Thermal Factors (Line2)
    !!--++                               6 Thermal Factors (Line3) + 6 Codes Thermal Factors (Line4)
    !!--++                               9 Thermal Factors (Line5) + 9 Codes Thermal Factors (Line6)
    !!--++
    !!--++    Internal Coordinates for Atoms (N_Atoms Lines)
    !!--++        Atom_Name(6)  Atom_Specie(4)  Coordinates(3)  N1  N2  N3  Biso  Occ [VARY]
    !!--++
    !!--++    if VARY is present as last option on the Internal Coordinates line,
    !!--++    then an extra line is read
    !!--++        Codes_Coordinates(3)   Code_BIso  Code_Occ
    !!--++
    !!--++    Control of error is present
    !!--++
    !!--++ Update: February - 2005
    !!
    Subroutine Read_Molecule_in_Var(File_dat,n_ini,n_end,Molecule)
       !---- Arguments ----!
       character(len=*), dimension(:), intent(in)      :: File_Dat
       integer,                        intent(in out)  :: N_Ini
       integer,                        intent(in)      :: N_End
       type (Molecule_type),           intent(out)     :: Molecule

       !---- Local variables -----!
       character(len=150)              :: line
       character(len=20),dimension(10) :: dire
       character(len=4)                :: var
       integer                         :: i,j,ic,npos,na,ier
       integer,dimension(10)           :: ivet
       real(kind=cp), dimension(10)    :: vet
       real(kind=cp),dimension(3,3)    :: Eu
       logical                         :: in_xtal

       in_xtal=.false.
       call Init_Err_Molec()

       n_ini=n_ini-1

       do
          n_ini=n_ini+1
          if (n_ini > n_end) then
             err_molec=.true.
             ERR_Molec_Mess="Not found Molecule"
             return
          end if
          line=adjustl(file_dat(n_ini))
          if (u_case(line(1:4)) /= "MOLE") cycle

          if (u_case(line(1:5)) == "MOLEX") in_xtal=.true.
          i=index(line,"!")
          if( i /= 0 ) line=line(1:i-1)

          !---- Coordinates format ----!
          call getword(line,dire,ic)
          if (ic /= 4) then
             err_molec=.true.
             ERR_Molec_Mess="Instruction: MOLE[X] N_Atoms Molecule_Name Coordinates_Type, not found! "
             return
          end if

          !---- Name and Number of Atoms in the molecule ----!
          read(unit=dire(2),fmt=*,iostat=ier) na
          if(ier /= 0) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the number of atoms in a molecule: "//trim(line)
             return
          else
             if (na > 0) then
                call init_molecule(molecule,na)
                Molecule%Name_mol =dire(3)
             else
                err_molec=.true.
                ERR_Molec_Mess="Error reading the number of atoms in a molecule: "//trim(line)
                return
             end if
          end if

          select case (dire(4)(1:1)) ! Coordinates_Type [C,S,F,Z]
             case ("C","c")
                molecule%coor_type="C"
             case ("F","f")
                molecule%coor_type="F"
             case ("S","s")
                molecule%coor_type="S"
             case ("Z","z")
                molecule%coor_type="Z"
             case default
                err_molec=.true.
                ERR_Molec_Mess="Coordinates Type for Molecule Unknown! "
                return
          end select ! dire

          exit !The molecule has been found

       end do

       !---- Initialize the crystal part of the molecule
       Molecule%xcentre    = 0.0_cp
       Molecule%orient     = 0.0_cp
       Molecule%therm_type = "   "
       Molecule%T_TLS      = 0.0_cp
       Molecule%L_TLS      = 0.0_cp
       Molecule%S_TLS      = 0.0_cp
       Molecule%in_xtal    = .false.
       Molecule%is_EulerMat=.false.

       if (in_xtal) then
          !---- Read the global coordinates of the centre of molecule and orientational angles
          do
             n_ini=n_ini+1
             if (n_ini > n_end) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading Molecule information! "
                return
             end if
             line=adjustl(file_dat(n_ini))
             if (line(1:1) =="!") cycle
             exit
          end do

          i=index(line,"!")
          if( i /= 0 ) line=line(1:i-1)

          call getword(line,dire,ic)
          if (ic /= 8) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the position and angles of the molecule: "//trim(Molecule%Name_mol)
             return
          end if

          line=trim(dire(1))//"   "//trim(dire(2))//"   "//trim(dire(3))
          call getnum(line,vet,ivet,ic)
          if (ic /= 3) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the position and angles of the molecule: "//trim(Molecule%Name_mol)
             return
          end if
          Molecule%xcentre=vet(1:3)

          line=trim(dire(4))//"   "//trim(dire(5))//"   "//trim(dire(6))
          call getnum(line,vet,ivet,ic)
          if (ic /= 3) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the position and angles of the molecule: "//trim(Molecule%Name_mol)
             return
          end if
          Molecule%orient=vet(1:3)

          Molecule%rot_type=adjustl(u_case(dire(7)))
          Molecule%therm_type=adjustl(u_case(dire(8)))

          do
             n_ini=n_ini+1
             if (n_ini > n_end) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading Molecule information! "
                return
             end if
             line=adjustl(file_dat(n_ini))
             if (line(1:1) =="!") cycle
             exit
          end do
          i=index(line,"!")
          if( i /= 0 ) line=line(1:i-1)

          call getnum(line,vet,ivet,ic)
          if (ic /= 6) then
             err_molec=.true.
             ERR_Molec_Mess="Error reading the position and angles of the molecule: "//trim(Molecule%Name_mol)
             return
          end if
          Molecule%mxcentre=vet(1:3)
          Molecule%mOrient =vet(4:6)

          Molecule%in_xtal = .true.

          !---- Set the Euler Matrix
          if (Molecule%rot_type /= "E") Molecule%rot_type="P"

          call Set_euler_matrix(Molecule%rot_type,  &
                                Molecule%orient(1),Molecule%orient(2),Molecule%orient(3),Eu)
                                !    Phi/alpha          Theta/beta          Chi/gamma
          Molecule%Euler=Eu
          Molecule%is_EulerMat=.true.

          !---- Read the THERMAL PARAMETERS
          if (Molecule%therm_type(1:1) == "T") then
             do
                n_ini=n_ini+1
                if (n_ini > n_end) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(file_dat(n_ini))
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 6) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the tensor T of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%T_TLS=vet(1:6)

             do
                n_ini=n_ini+1
                if (n_ini > n_end) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(file_dat(n_ini))
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 6) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the codes of tensor T of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%mT_TLS=vet(1:6)
          end if

          if (Molecule%therm_type(2:2) == "L") then
             do
                n_ini=n_ini+1
                if (n_ini > n_end) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(file_dat(n_ini))
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 6) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the tensor L of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%L_TLS=vet(1:6)

             do
                n_ini=n_ini+1
                if (n_ini > n_end) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(file_dat(n_ini))
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 6) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the codes of the tensor L of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%mL_TLS=vet(1:6)
          end if

          if (Molecule%therm_type(3:3) == "S") then
             do
                n_ini=n_ini+1
                if (n_ini > n_end) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(file_dat(n_ini))
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 9) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the tensor S of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%S_TLS(1,:)=vet(1:3)
             Molecule%S_TLS(2,:)=vet(4:6)
             Molecule%S_TLS(3,:)=vet(7:9)

             do
                n_ini=n_ini+1
                if (n_ini > n_end) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading Molecule information! "
                   return
                end if
                line=adjustl(file_dat(n_ini))
                if (line(1:1) =="!") cycle
                exit
             end do
             i=index(line,"!")
             if( i /= 0 ) line=line(1:i-1)

             call getnum(line,vet,ivet,ic)
             if (ic /= 9) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading the code of tensor S of the molecule: "//trim(Molecule%Name_mol)
                return
             end if
             Molecule%mS_TLS(1,:)=vet(1:3)
             Molecule%mS_TLS(2,:)=vet(4:6)
             Molecule%mS_TLS(3,:)=vet(7:9)
          end if

       end if  !(in_xtal)

       !---- Read the internal coordinates of the atoms in the molecule
       !---- Read the Z-matrix/Cartesian/spherical/Fractional coordinates of the molecule
       molecule%is_connect=.true.
       do i=1,na
          do
             n_ini=n_ini+1
             if (n_ini > n_end) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading Molecule information! "
                return
             end if
             line=adjustl(file_dat(n_ini))
             if (line(1:1) =="!") cycle
             exit
          end do
          j=index(line,"!")
          if( j /= 0 ) line=line(1:j-1)

          !---- Atom Name ---!
          call Cutst(line,ic,Molecule%Atname(i))

          !---- Atom specie ----!
          call Cutst(line,ic,Molecule%Atsymb(i))

          !---- Passing Codes? ----!
          call getword(line,dire,ic)
          var=adjustl(dire(ic))
          var=u_case(var)
          if (var == "VARY") then
             ic=len_trim(line)
             npos=index(line(1:ic)," ",back=.true.)
             if (npos <=0) then
                err_molec=.true.
                ERR_Molec_Mess="Error reading Molecule information (II)! "
                return
             end if
             line=line(1:npos)
          end if

          !---- Rest of Information ----!
          vet =0.0
          ivet=0
          call getnum(line,vet,ivet,ic)
          select case (ic)
             case (0)
                Molecule%I_Coor(:,i)=0.0
                Molecule%Conn(:,i)  =0
                Molecule%Biso(i)    =0.5
                Molecule%Occ(i)     =1.0
             case (1)
                Molecule%I_Coor(1,i)  =vet(1)
                Molecule%I_Coor(2:3,i)=0.0
                Molecule%conn(:,i)    =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (2)
                Molecule%I_Coor(1:2,i)=vet(1:2)
                Molecule%I_Coor(3,i)  =0.0
                Molecule%conn(:,i)    =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (3)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(:,i)    =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (4)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(1,i)    =ivet(4)
                Molecule%conn(2:3,i)  =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (5)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(1:2,i)  =ivet(4:5)
                Molecule%conn(3,i)    =0
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (6)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(:,i)    =ivet(4:6)
                Molecule%biso(i)      =0.5
                Molecule%Occ(i)       =1.0

             case (7)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(:,i)    =ivet(4:6)
                Molecule%biso(i)      =vet(7)
                Molecule%Occ(i)       =1.0

             case (8)
                Molecule%I_Coor(:,i)  =vet(1:3)
                Molecule%conn(:,i)    =ivet(4:6)
                Molecule%biso(i)      =vet(7)
                Molecule%Occ(i)       =vet(8)

             case default
                err_molec=.true.
                ERR_Molec_Mess="Error reading the atoms in the molecule: "//trim(Molecule%Name_mol)
                return

          end select ! ic

          if (Molecule%coor_type == "Z") then

             if (i == 2 .and. (ivet(4) ==0 .and. ivet(5) ==0 .and. ivet(6) ==0)) then
                Molecule%conn(1,i)=1
             end if

             if (ivet(4) >= i .or. ivet(5) >= i .or. ivet(6) >= i )                err_molec=.true.
             if (i == 3 .and. (ivet(4) == 0 .or. ivet(5) == 0))                    err_molec=.true.
             if (i > 3 .and. (ivet(4) == 0 .or. ivet(5) == 0 .or. ivet(6) == 0))   err_molec=.true.
             if (err_molec) then
                ERR_Molec_Mess = "The Z-matrix connectivity is wrong: "//trim(line)
                return
             end if
          else
             if (ivet(4) >= i .or. ivet(5) >= i .or. ivet(6) >= i )               molecule%is_connect=.false.
             if (i == 3 .and. (ivet(4) == 0 .or. ivet(5) == 0))                   molecule%is_connect=.false.
             if (i > 3 .and. (ivet(4) == 0 .or. ivet(5) == 0 .or. ivet(6) == 0))  molecule%is_connect=.false.
          end if

          Molecule%mI_Coor(:,i)=0.0
          Molecule%mbiso(i)  =0.0
          Molecule%mocc(i)   =0.0

          if (var == "VARY") then
             do
                n_ini=n_ini+1
                if (n_ini > n_end) then
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading the refinement codes of atoms in the molecule: "//trim(Molecule%Name_mol)
                   return
                end if
                line=adjustl(file_dat(n_ini))
                if (line(1:1) =="!") cycle
                exit
             end do
             j=index(line,"!")
             if( j /= 0 ) line=line(1:j-1)

             call getnum(line,vet,ivet,ic)
             select case (ic)
                case (3)
                   Molecule%mI_Coor(:,i)=vet(1:3)

                case (5)
                   Molecule%mI_Coor(:,i)=vet(1:3)
                   Molecule%mbiso(i)  =vet(4)
                   Molecule%mocc(i)   =vet(5)

                case default
                   err_molec=.true.
                   ERR_Molec_Mess="Error reading the refinement codes of atoms in the molecule: "//trim(Molecule%Name_mol)
                   return
             end select
          end if

       end do

       return
    End Subroutine Read_Molecule_in_Var

    !!----
    !!---- Subroutine Set_Euler_Matrix(Rt,Phi,Theta,Chi,Eu)
    !!----    character(len=*),              intent ( in) :: Rt
    !!----    real(kind=cp),                 intent ( in) :: Phi,Theta,Chi
    !!----    real(kind=cp), dimension(3,3), intent (out) :: Eu
    !!----
    !!----    Subroutine to obtain the Euler active matrix to transform a point
    !!----    to another point. For instance the internal coordinates of a molecule
    !!----    can be transformed to absolute positions using columns vectors.
    !!----    If the Cartesian coordinates of an atom in the molecular frame is the
    !!----    column vector  Xm, the cartesian coordinates in the crystal frame X
    !!----    are obtained from:  X = Eu Xm
    !!----    The internal coordinates of a point are obtained from Xm = EuT X.
    !!----    The character variable "rt" indicates the type of Euler angles provided.
    !!----    If rt="E", the angles PHI,THETA,CHI correspond to the conventional
    !!----    Euler angles ALPHA, BETA, GAMMA. Otherwise, they correspond to the
    !!----    2nd setting, allowing to interpret PHI and THETA as the polar angles of
    !!----    the molecular frame Zm-axis, and CHI a rotation around Zm.
    !!----
    !!----   Update: February - 2005
    !!
    Subroutine Set_Euler_Matrix(Rt,Phi,Theta,Chi,Eu)
       !---- Arguments ----!
       character(len=*),              intent ( in) :: Rt
       real(kind=cp),                 intent ( in) :: Phi,Theta,Chi
       real(kind=cp), dimension(3,3), intent (out) :: Eu

       !---- Local Variables ----!
       real(kind=cp) :: PH,TH,CH

       TH=THETA
       if (rt(1:1) == "E") then
          PH=PHI+90.0_cp
          CH=CHI-90.0_cp
       else
          PH=PHI
          CH=CHI
       end if
       Eu(1,1) =  cosd(PH)* cosd(TH)* cosd(CH) - sind(PH)* sind(CH)
       Eu(1,2) = -cosd(PH)* cosd(TH)* sind(CH) - sind(PH)* cosd(CH)
       Eu(1,3) =  cosd(PH)* sind(TH)
       Eu(2,1) =  sind(PH)* cosd(TH)* cosd(CH) + cosd(PH)* sind(CH)
       Eu(2,2) = -sind(PH)* cosd(TH)* sind(CH) + cosd(PH)* cosd(CH)
       Eu(2,3) =  sind(PH)* sind(TH)
       Eu(3,1) = -cosd(CH)* sind(TH)
       Eu(3,2) =  sind(CH)* sind(TH)
       Eu(3,3) =            cosd(TH)

       return
    End Subroutine Set_Euler_Matrix

    !!----
    !!---- Subroutine Spherical_to_Cartesian(Molecule,NewMolecule)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Molecule_type), intent(   out), optional :: Newmolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a molecule from Spherical
    !!----    coordinates to  cartesian coordinaters.
    !!----    If a second argument is present the subroutine creates a new molecule
    !!----    (copy of the old one) with spherical coordinates, preserving
    !!----    the input molecule in Cartesian Coordinates. Otherwise the input
    !!----    molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Spherical_to_Cartesian(Molecule,NewMolecule)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Molecule_type), intent(   out), optional :: Newmolecule

       !---- Local variables -----!
       integer                     :: i,na
       real(kind=cp)               :: r, theta, phi

       type (Molecule_type)        :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "S") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Spherical_to_Cartesian: the input molecule is not in Spherical coordinates"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Spherical_to_Cartesian: No atoms are defined"
          return
       end if

       call init_molecule(newmol,na)
       NewMol=Molecule

       !---- Start calculations for each atom of the molecule ----!
       do i=1,na
          r     = Molecule%I_coor(1,i)
          theta = Molecule%I_coor(2,i)
          phi   = Molecule%I_coor(3,i)
          NewMol%I_coor(1,i) = r*sind(theta)*cosd(phi)
          NewMol%I_coor(2,i) = r*sind(theta)*sind(phi)
          NewMol%I_coor(3,i) = r*cosd(theta)
       end do
       NewMol%coor_type="C"

       if (present(NewMolecule)) then
          call Init_molecule(NewMolecule,Newmol%natoms)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Spherical to Cartesian: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Spherical_to_Cartesian

    !!----
    !!---- Subroutine Spherical_to_Fractional(Molecule,Cell,NewMolecule)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Crystal_Cell_Type), intent(in)           :: Cell
    !!----    type (Molecule_type), intent(   out), optional :: Newmolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a
    !!----    molecule from Spherical coordinates to  Fractional coordinaters.
    !!----    If a second argument is present the subroutine creates a new
    !!----    molecule (copy of the old one) with Fractional coordinates,
    !!----    preserving the input molecule in Spherical Coordinates. Otherwise
    !!----    the input molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Spherical_to_Fractional(Molecule, Cell, NewMolecule)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Crystal_Cell_Type), intent(in)           :: Cell
       type (Molecule_type), intent(   out), optional :: NewMolecule

       !---- Local Variables ----!
       integer                     :: na
       type (Molecule_type)        :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "S") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Spherical_to_Fractional: the input molecule is not in Spherical coordinates"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Spherical_to_Fractional: No atoms are defined"
          return
       end if

       !---- Step 1----!
       call init_molecule(newmol,na)
       newmol=Molecule
       call Spherical_to_Cartesian(NewMol)
       if (err_molec) then
          ERR_Molec_Mess="Error in Spherical_to_Fractional: Intermediate procedure fail (I)!"
          return
       end if

       !---- Step 2 ----!
       call Cartesian_to_Fractional(NewMol,Cell)
       if (err_molec) then
          ERR_Molec_Mess="Error in Spherical_to_Fractional: Intermediate procedure fail (II)!"
          return
       end if

       !---- Step 3 ----!
       if (present(newmolecule)) then
          call Init_molecule(NewMolecule,na)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Spherical to Fractional: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Spherical_to_Fractional

    !!----
    !!---- Subroutine Spherical_to_Zmatrix(Molecule,NewMolecule,Cell)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Molecule_type), intent(   out), optional :: Newmolecule
    !!----    Type(Crystal_Cell_Type), intent(in),  optional :: Cell
    !!----
    !!----    Subroutine to transform the internal coordinates of a
    !!----    molecule from Spherical coordinates to  Zmatrix coordinaters.
    !!----    If a second argument is present the subroutine creates a new
    !!----    molecule (copy of the old one) with Zmatrix coordinates,
    !!----    preserving the input molecule in Spherical Coordinates. Otherwise
    !!----    the input molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Spherical_to_Zmatrix(Molecule, NewMolecule,Cell)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Molecule_type), intent(   out), optional :: NewMolecule
       Type(Crystal_Cell_Type), intent(in),  optional :: Cell

       !---- Local Variables ----!
       integer                     :: na
       type (Molecule_type)        :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "S") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Spherical_to_ZMatrix: the input molecule is not in Spherical coordinates"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Spherical_to_ZMatrix: No atoms are defined"
          return
       end if

       !---- Step 1----!
       call init_molecule(newmol,na)
       newmol= Molecule
       call Spherical_to_Cartesian(NewMol)
       if (err_molec) then
          ERR_Molec_Mess="Error in Spherical_to_Zmatrix: Intermediate procedure fail (I)!"
          return
       end if

       !---- Step 2 ----!
       if(present(Cell)) then
          call Cartesian_to_Zmatrix(NewMol,Cell=Cell)
       else
          call Cartesian_to_Zmatrix(NewMol)
       end if
       if (err_molec) then
          ERR_Molec_Mess="Error in Spherical_to_Zmatrix: Intermediate procedure fail (II)!"
          return
      end if

       !---- Step 3 ----!
       if (present(newmolecule)) then
          call Init_molecule(NewMolecule,na)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in Spherical to ZMatrix: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Spherical_to_Zmatrix

    !!----
    !!---- Subroutine Write_Free_Atoms(AtmF,N,Lun)
    !!----    type (Atom_type), dimension(:), intent(in) :: AtmF
    !!----    integer,                        intent(in) :: N
    !!----    integer, optional,              intent(in) :: Lun
    !!----
    !!----    Write information about Free Atoms for Molecular Crystal
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Free_Atoms(AtmF,N,Lun)
       !---- Arguments ----!
       type (Atom_type), dimension(:), intent(in) :: AtmF
       integer,                        intent(in) :: N
       integer, optional,              intent(in) :: Lun

       !---- Local Variables ----!
       integer :: i,uni

       uni=6
       if (present(lun)) uni=lun

       write(unit=uni,fmt="(a)")     " "
       write(unit=uni,fmt="(a,i4)")  " => Number of Free Atoms: ",N
       write(unit=uni,fmt="(a)")     " "
       write (unit=uni,fmt="(T5,a)") " Atom      Chem        x/a        y/b        z/c       Occ     Biso"
       write (unit=uni,fmt="(T5,a)") "===================================================================="
       do i=1,N
          write(unit=uni,fmt="(T5,a,T16,a,T21,5f11.4)") atmF(i)%Lab,atmF(i)%chemsymb,atmF(i)%x,atmF(i)%occ,atmF(i)%biso
       end do

       return
    End Subroutine Write_Free_Atoms

    !!----
    !!---- Subroutine Write_Molecular_Crystal(MolCrys,Lun)
    !!----    type (Molecular_Crystal_type), intent(in) :: MolCrys
    !!----    integer, optional,             intent(in) :: Lun
    !!----
    !!----    Write information about Molecular Crystal
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Molecular_Crystal(MolCrys,Lun)
       !---- Arguments ----!
       type(Molecular_Crystal_Type), intent(in) :: MolCrys
       integer, optional,            Intent(in) :: Lun

       !---- Local Variables ----!
       integer :: i,uni

       uni=6
       if (present(lun)) uni=lun

       write(unit=uni,fmt="(/,/,a)") "      Molecular Crystal Information  "
       write(unit=uni,fmt="(a)")     "      ----------------------------- "

       write(unit=uni,fmt="(a)")     " "
       call Write_Crystal_Cell(MolCrys%Cell,uni)
       write(unit=uni,fmt="(a)")     " "

       write(unit=uni,fmt="(a)")     " "
       call Write_SpaceGroup(MolCrys%SPG,uni)
       write(unit=uni,fmt="(a)")     " "

       if (MolCrys%N_Free > 0) then
          write(unit=uni,fmt="(a)")     " "
          call Write_Free_Atoms(MolCrys%Atm,MolCrys%N_Free,uni)
          write(unit=uni,fmt="(a)")     " "
       end if

       if (MolCrys%N_Mol > 0) then
          do i=1,MolCrys%N_Mol
             write(unit=uni,fmt="(a)")     " "
             call Write_Molecule(MolCrys%Mol(i),uni)
             write(unit=uni,fmt="(a)")     " "
          end do
       end if

       return
    End Subroutine Write_Molecular_Crystal

    !!----
    !!---- Subroutine Write_Molecule(Molecule,Lun)
    !!----    type (Molecule_type), intent(in)           :: Molecule
    !!----    integer,              intent(in), optional :: Lun
    !!----
    !!----    Write information about molecule
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Write_Molecule(Molecule,Lun)
       !---- Arguments ----!
       type (Molecule_type), intent(in):: Molecule
       integer,optional,     intent(in):: Lun

       !---- Local variables -----!
       integer            :: i,uni,j
       character(len=4)   :: var
       real(kind=cp), dimension(3  ) :: geom_cent

       uni=6
       if (present(lun)) uni=lun

       write(unit=uni,fmt="(/,/,a,a)")    " =>  MOLECULE of name :  ",trim(Molecule%Name_mol)
       select case (molecule%coor_type)
          case ("C","c")
             write(unit=uni,fmt="(a)")     "            Type of Molecular description: CARTESIAN COORDINATES"
          case ("F","f")
             write(unit=uni,fmt="(a)")     "            Type of Molecular description: FRACTIONAL COORDINATES"
          case ("S","s")
             write(unit=uni,fmt="(a)")     "            Type of Molecular description: SPHERICAL COORDINATES"
          case ("Z","z")
             write(unit=uni,fmt="(a)")     "            Type of Molecular description: Z-MATRIX"
          case default
             write(unit=uni,fmt="(a)")     "            Type of Molecular description: UNKNOWN "
       end select ! molecule%coor_type

       write(unit=uni,fmt="(a,i3)")     "                          Number of atoms: ",  Molecule%natoms
       if (Molecule%in_xtal) then
          write(unit=uni,fmt="(a,3f11.5)")   "         Fractional coordinates of centre: ",  Molecule%xcentre
          write(unit=uni,fmt="(a,3f11.5)")   "                         Refinement codes: ",  Molecule%mxcentre

          if (Molecule%rot_type == "E") then
             write(unit=uni,fmt="(a,3f11.5,a,3f9.5,a)") &
                                 "   Orientation EULER angles (alp,bet,gam): ",  Molecule%orient,&
                                 " (radians:", Molecule%orient*to_rad,")"
          else
             write(unit=uni,fmt="(a,3f11.5,a,3f9.5,a)") &
                                 "   Orientation POLAR angles (PHI,THE,CHI): ",  Molecule%orient,&
                                 " (radians:", Molecule%orient*to_rad,")"
          end if
          write(unit=uni,fmt="(a,3f11.5)") "                         Refinement codes: ",  Molecule%mOrient

          if (Molecule%therm_type(1:1) == "T") then
             write(unit=uni,fmt="(a,6f11.5)")"       T-tensor (T11,T22,T33,T12,T13,T23): ", Molecule%T_TLS
             write(unit=uni,fmt="(a,6f11.5)")"                         Refinement codes: ", Molecule%mT_TLS
          end if

          if (Molecule%therm_type(2:2) == "L") then
             write(unit=uni,fmt="(a,6f11.5)")"       L-tensor (L11,L22,L33,L12,L13,L23): ", Molecule%L_TLS
             write(unit=uni,fmt="(a,6f11.5)")"                         Refinement codes: ", Molecule%mL_TLS
          end if

          if (Molecule%therm_type(3:3) == "S") then
             write(unit=uni,fmt="(a,3f11.5,tr5,3f11.5)")"       S-tensor             (S11,S12,S13): ", &
                                              Molecule%S_TLS(1,:), Molecule%mS_TLS(1,:)
             write(unit=uni,fmt="(a,3f11.5,tr5,3f11.5)")"     + Refinement codes     (S21,S22,S23): ", &
                                              Molecule%S_TLS(2,:), Molecule%mS_TLS(2,:)
             write(unit=uni,fmt="(a,3f11.5,tr5,3f11.5)")"                            (S31,S32,S33): ", &
                                              Molecule%S_TLS(3,:), Molecule%mS_TLS(3,:)
          end if

          select case (Molecule%coor_type)
             case ("C","c")
                write(unit=uni,fmt="(t29,a)")"Atom    Type        XC          YC          ZC    N1  N2  N3      Biso        Occ "
             case ("F","f")
                write(unit=uni,fmt="(t29,a)")"Atom    Type        X           Y           Z     N1  N2  N3      Biso        Occ "
             case ("S","s")
                write(unit=uni,fmt="(t29,a)")"Atom    Type    distance      Theta       Phi     N1  N2  N3      Biso        Occ "
             case ("Z","z")
                write(unit=uni,fmt="(t29,a)")"Atom    Type    distance  Bond-Angle Torsion-Ang  N1  N2  N3      Biso        Occ "
             case default
                write(unit=uni,fmt="(t29,a)")"Atom    Type      Coor1       Coor2       Coor3   N1  N2  N3      Biso        Occ "
          end select ! Molecule%coor_type

       else  !(Molecule%in_xtal)

          select case (Molecule%coor_type)
             case ("C","c")
                write(unit=uni,fmt="(t29,a)")"Atom    Type        XC          YC          ZC    N1  N2  N3 "
             case ("F","f")
                write(unit=uni,fmt="(t29,a)")"Atom    Type        X           Y           Z     N1  N2  N3 "
             case ("S","s")
                write(unit=uni,fmt="(t29,a)")"Atom    Type    distance      Theta       Phi     N1  N2  N3 "
             case ("Z","z")
                write(unit=uni,fmt="(t29,a)")"Atom    Type    distance  Bond-Angle Torsion-Ang  N1  N2  N3 "
             case default
                write(unit=uni,fmt="(t29,a)")"Atom    Type      Coor1       Coor2       Coor3   N1  N2  N3 "
          end select ! Molecule%coor_type

       end if  !(Molecule%in_xtal)

          geom_cent=0.0_cp

          if (Molecule%in_xtal ) then
             do i=1,Molecule%natoms
                  if(Molecule%AtSymb(i) /= "ZE") geom_cent=geom_cent + Molecule%I_Coor(:,i)
                  write(unit=uni,fmt="(t29,a,tr2,a,3f12.5,3i4,2f12.5)")  &
                       Molecule%AtName(i), Molecule%AtSymb(i),Molecule%I_Coor(:,i),  &
                       Molecule%Conn(:,i), Molecule%Biso(i),  Molecule%Occ(i)
                  var="    "
                  do j=1,3
                     if (abs(Molecule%mI_Coor(j,i)) > eps) var="VARY"
                  end do
                  if (abs(Molecule%mbiso(i)) > eps)      var="VARY"
                  if (abs(Molecule%mocc(i))  > eps)      var="VARY"
                  if (var == "VARY") then
                     write(unit=uni,fmt="(t41,3f12.5,tr12,2f12.5)")  Molecule%mI_Coor(:,i), &
                          Molecule%mbiso(i),Molecule%mocc(i)
                  end if
             end do
          else
             do i=1,Molecule%natoms
                  if(Molecule%AtSymb(i) /= "DU") geom_cent=geom_cent + Molecule%I_Coor(:,i)
                  write(unit=uni,fmt="(t29,a,tr2,a,3f12.5,3i4       )")  &
                  Molecule%Atname(i), Molecule%Atsymb(i), Molecule%I_coor(:,i),  &
                  Molecule%conn(:,i)
             end do
          end if

       if(      molecule%coor_type == "C" .or. molecule%coor_type == "c" &
           .or. molecule%coor_type == "F" .or. molecule%coor_type == "f") then
           geom_cent=geom_cent/real(Molecule%natoms)
           write(unit=uni,fmt="(//,a,3F10.5)")  "  => Geometrical centre of molecule ( "//trim(Molecule%Name_mol)//" ):", geom_cent
       end if

       write(unit=uni,fmt="(/,a)")              "  => Euler Matrix of molecule ( "//trim(Molecule%Name_mol)//" ):"
       do i=1,3
          write(unit=uni,fmt="(t29,3f10.5)")  Molecule%Euler(i,:)
       end do

       return
    End Subroutine Write_Molecule

    !!----
    !!---- Subroutine Zmatrix_to_Cartesian(Molecule,NewMolecule)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Molecule_type), intent(   out), optional :: NewMolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a molecule from
    !!----    Z-matrix to cartesian coordinates.
    !!----    If a second argument is present the subroutine creates a new molecule
    !!----    (copy of the old one) with cartesian coordinates, preserving
    !!----    the input molecule. Otherwise the input molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Zmatrix_to_Cartesian(Molecule,NewMolecule)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Molecule_type), intent(   out), optional :: NewMolecule

       !---- Local variables -----!
       integer                     :: i,na,j,k,n
       real(kind=cp)               :: dist, ang
       real(kind=cp), dimension(3) :: ci,ri,rj,rk,rn

       type (Molecule_type)        :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "Z") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Zmatrix_to_Cartesian: the input molecule is not a Z-matrix"
          return
       end if

       na= Molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Zmatrix_to_Cartesian: Not atoms are defined"
          return
       end if

       call init_molecule(newmol,na)
       NewMol=Molecule

       !---- Start calculations for each atom of the molecule ----!

       !---- First atom is always at origin (Z-matrix)
       NewMol%I_coor(:,1) = 0.0_cp
       NewMol%conn(:,1) = 0

       !---- Second atom is always along "x"
       NewMol%I_coor(2:3,2) = 0.0
       NewMol%conn(2:3,2) = 0
       NewMol%conn(1,2)   = 1

       !--- Third atom is always in the "xy" plane       !A(i) d_ij  ang_ijk   dang_ijkl  j k l
       if (NewMol%conn(1,3) == 1) then
          NewMol%conn(2,3) = 2
          NewMol%conn(3,3) = 0
          dist= NewMol%I_coor(1,3)
          ang = NewMol%I_coor(2,3)
          NewMol%I_coor(1,3) = dist * cosd(ang)
          NewMol%I_coor(2,3) = dist * sind(ang)
          NewMol%I_coor(3,3) = 0.0_cp
       else
          NewMol%conn(1,3) = 2
          NewMol%conn(2,3) = 1
          NewMol%conn(3,3) = 0
          dist= NewMol%I_coor(1,3)
          ang = NewMol%I_coor(2,3)
          NewMol%I_coor(1,3) = dist * cosd(180.0_cp-ang) +  NewMol%I_coor(1,2)
          NewMol%I_coor(2,3) = dist * sind(180.0_cp-ang)
          NewMol%I_coor(3,3) = 0.0_cp
       end if

       do i=4,na
          ci(:) = NewMol%I_coor(:,i)
          j     = NewMol%conn(1,i)         !The connectivity is needed for the Z-matrix description
          k     = NewMol%conn(2,i)         !If the connectivity is given it is possible to transform to
          n     = NewMol%conn(3,i)         !Z-matrix if cartesian/spherical coordinates are given.
          if (j == 0 .or. k == 0 .or. n == 0) cycle
          rj(:) = NewMol%I_coor(:,j)
          rk(:) = NewMol%I_coor(:,k)
          rn(:) = NewMol%I_coor(:,n)
          call get_cartesian_from_Z(ci,ri,rj,rk,rn)
          NewMol%I_coor(:,i) = ri
       end do
       NewMol%coor_type="C"

       if (present(NewMolecule)) then
          call Init_molecule(NewMolecule,na)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in ZMatrix to Cartesian: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Zmatrix_to_Cartesian

    !!----
    !!---- Subroutine Zmatrix_to_Fractional(Molecule,Cell,NewMolecule)
    !!----    type (Molecule_type),     intent(in out)           :: Molecule
    !!----    type (Crystal_Cell_Type), intent(in    )           :: Cell
    !!----    type (Molecule_type),     intent(   out), optional :: NewMolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a molecule from
    !!----    Z-matrix to fractional coordinates.
    !!----    If a third argument is present the subroutine creates a new molecule
    !!----    (copy of the old one) with fractional coordinates, preserving
    !!----    the input molecule in Z-matrix. Otherwise the input molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Zmatrix_to_Fractional(Molecule,Cell,NewMolecule)
       !---- Arguments ----!
       type (Molecule_type),     intent(in out)           :: Molecule
       type (Crystal_Cell_Type), intent(in    )           :: Cell
       type (Molecule_type),     intent(   out), optional :: NewMolecule

       !---- Local variables -----!
       integer                       :: na
       type (Molecule_type)          :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "Z") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Zmatrix_to_Fractional: the input molecule is not in Zmatrix coordinates"
          return
       end if

       na=molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Zmatrix_to_Fractional: No atoms found"
          return
       end if

       call init_molecule(newmol,na)
       newmol= molecule
       call Zmatrix_to_Cartesian(newmol)
       call Cartesian_to_Fractional(newmol,cell)

       if (present(NewMolecule)) then
          call Init_molecule(NewMolecule,Newmol%natoms)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in ZMatrix_to_Fractional: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Zmatrix_to_Fractional

    !!----
    !!---- Subroutine Zmatrix_to_Spherical(Molecule,NewMolecule)
    !!----    type (Molecule_type), intent(in out)           :: Molecule
    !!----    type (Molecule_type), intent(   out), optional :: Newmolecule
    !!----
    !!----    Subroutine to transform the internal coordinates of a
    !!----    molecule from Zmatrix coordinates to  Spherical coordinaters.
    !!----    If a second argument is present the subroutine creates a new
    !!----    molecule (copy of the old one) with Spherical coordinates,
    !!----    preserving the input molecule in Zmatrix Coordinates. Otherwise
    !!----    the input molecule is changed on output.
    !!----    Control of error is present
    !!----
    !!---- Update: February - 2005
    !!
    Subroutine Zmatrix_to_Spherical(Molecule, NewMolecule)
       !---- Arguments ----!
       type (Molecule_type), intent(in out)           :: Molecule
       type (Molecule_type), intent(   out), optional :: NewMolecule

       !---- Local Variables ----!
       integer                     :: na
       type (Molecule_type)        :: Newmol

       !---- Controls ----!
       if (molecule%coor_type /= "Z") then
          err_molec=.true.
          ERR_Molec_Mess="Error in Zmatrix_to_Spherical: the input molecule is not in Zmatrix coordinates"
          return
       end if

       na=molecule%natoms
       if (na <= 0) then
          err_molec=.true.
          ERR_Molec_Mess="Error in Zmatrix_to_Fractional: No atoms found"
          return
       end if

       !---- Step 1----!
       call init_Molecule(newmol,na)
       newmol=Molecule
       call Zmatrix_to_Cartesian(NewMol)
       if (err_molec) then
          ERR_Molec_Mess="Error in Zmatrix_to_Spherical: Intermediate procedure fail (I)!"
          return
       end if

       !---- Step 2 ----!
       call Cartesian_to_Spherical(NewMol)
       if (err_molec) then
          ERR_Molec_Mess="Error in Zmatrix_to_Spherical: Intermediate procedure fail (II)!"
          return
       end if

       !---- Step 3 ----!
       if (present(newmolecule)) then
          call Init_molecule(NewMolecule,na)
          if (NewMolecule%natoms <=0) then
             err_molec=.true.
             ERR_Molec_Mess="Error in ZMatrix to Spherical: The optional variable was not dimensioned!"
             return
          end if
          NewMolecule=newmol
       else
          Molecule=newmol
       end if

       return
    End Subroutine Zmatrix_to_Spherical

 End Module CFML_Molecular_Crystals
