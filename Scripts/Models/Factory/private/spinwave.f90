! SPINWAVE 2.2 - (c) S. Petit, LLB CEA/CNRS - Sept 2018
! http://www-llb.cea.fr/logicielsllb/SpinWave/SW.html

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	lecture de mots clefs
!
	module m_borne
	implicit none

	integer   :: mxcar, mxcle, mxlmc,mxlvc
	parameter (mxcar=132,mxcle=40,mxlmc=12,mxlvc=800)
	
	end module m_borne
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
	module m_cartes
	use m_borne
	implicit none
	contains
!
!--------------------------------------------------------------------------
!	analyseur cartes de donnees
!
	subroutine anacart(ligne,nbcle,mc,av,iv,rv)
	implicit none
	character(len=mxcar), intent(in)   :: ligne
	integer             , intent(out)  :: nbcle
	character(len=mxcar), intent(out)  :: mc(:)
	character(len=mxcar), intent(out)  :: av(:)
	integer             , intent(out)  :: iv(:)
	real                , intent(out)  :: rv(:)

	character(len=mxcar)               :: carte
	character                          :: car
	integer                            :: i,j,k,nucle,lncar,ios
	integer                            :: lmc,lvc
	integer, pointer                   :: pos(:)
!
!	-------------------------------------------------------------------
!	on initialise la lecture des mots clefs sur la ligne
!
	carte = ligne
	nbcle = 1
	do nucle = 1, mxcle
		do i = 1, mxlvc
			av(nucle)(i:i) = " "
		enddo
		do i = 1, mxlmc
			mc(nucle)(i:i) = " "
		enddo
	enddo
!
!	-------------------------------------------------------------------
!	on cherche la taille effective notee i de la ligne
!	partant de la taille maximale mxvar, on parcourt la ligne a 
!	l'envers jusqu'a ce qu'on trouve un blanc " "
!	si la ligne est blanche, on provoque une erreur
!
	do i = mxcar,1,-1
		car = carte(i:i)
		if(car.ne." ") exit
	enddo
	lncar = i
	if(lncar.eq.0) call tilt("anacart","pas de ligne blanche")
!
!	-------------------------------------------------------------------
!	on supprime les blancs internes et on recalcule la taille 
!	effective de la ligne
!
	do i = lncar,1,-1
		car = carte(i:i)
		if(car.eq." ") then
			carte(i:lncar-1)   = carte(i+1:lncar)
			carte(lncar:lncar) = car
		endif
	enddo
	do i = mxcar,1,-1
		car = carte(i:i)
		if (car.ne." ") exit
	enddo
	lncar = i
!
!	-------------------------------------------------------------------
!	recherche du nombre de mots clefs sur la ligne
!	separateur = ","
!	sauf en premiere position ce qui cree une erreur
!
	do i = 1,lncar
		car = carte(i:i)
		if(car.eq.",") then
			if((i.eq.1).or.(i.eq.lncar)) then
				call tilt("anacart",&
				"pas de virugule en debut ou fin de ligne")
			else
				if(nbcle.lt.mxcle) then
					nbcle = nbcle + 1
				else
					call tilt("anacart","trop de mot clef")
				endif
			endif
		endif
	enddo
!
!	-------------------------------------------------------------------
!	recherche des positions des mots clef
!	on renvoie un nbcle+1 mot clef fictif
!	pour le traitement particulier du dernier mot clef reel
!
	allocate(pos(nbcle+1))
	pos(1) = 1
	j      = 1
	do i = 1, lncar
		car = carte(i:i)
		if(car.eq.",") then
			j = j + 1
			pos(j)= i+1
		endif
	enddo
	pos(nbcle+1) = lncar+2
!
!	-------------------------------------------------------------------
!	lecture des mots clefs avec traitement eventuel des affectations
!
	do i = 1, nbcle
		k = 0
		do j = pos(i), pos(i+1)-2
			car = carte(j:j)
			if(car.eq."=") then
				k = j
			endif
		enddo
		
		if(k.ne.0) then
			lmc = k-1-pos(i)+1
			lvc = pos(i+1)-2-(k+1)+1
			if(lmc.gt.mxlmc) then
				call tilt("anacart",&
				"mot clef a + de mxlmc caracteres")
			endif
			if(lvc.gt.mxlvc) then
				call tilt("anacart",&
				"valeur clef a + de mxlvc caracteres")
			endif

			mc(i) = carte(pos(i):k-1)
			av(i) = carte(k+1:pos(i+1)-2)
		else
			lmc = pos(i+1)-2-pos(i)+1
			if(lmc.gt.mxlmc) then
				call tilt("anacart",&
				"mot clef a + de mxlmc caracteres")
			endif
			mc(i) = carte(pos(i):pos(i+1)-2)			
		endif
	enddo
	deallocate(pos)
!
!	--------------------------------------------------------------------
!	traduire valeur du mot clef en entier ou reel
!	
	if(nbcle.gt.0) then
		do nucle = 1, nbcle
			read(av(nucle), "(i12)"  , iostat=ios) iv(nucle)
			read(av(nucle), "(e12.0)", iostat=ios) rv(nucle)
		enddo
	endif
!
	end subroutine anacart
!
!
!--------------------------------------------------------------------------
!	mettre en majuscules
!	
	subroutine capit(carte)
	implicit none
	character(len=mxcar), intent(out)  :: carte
	integer                            :: interval,i,long

	interval = ichar("A") - ichar("a")
	do i = 1, mxcar
		if((carte(i:i).le."z").and.(carte(i:i).ge."a")) then
		carte(i:i) = char(ichar(carte(i:i)) + interval)
		endif
	enddo
	end subroutine capit
!
!
!--------------------------------------------------------------------------
!	traitement des erreurs
!
	subroutine tilt(str1,str2)
	implicit none
	character(len=*), intent(in)   :: str1,str2
	integer                        :: i
	write(*,*) "***********************************"
	write(*,*) "***********************************"
	do i = 1, 10
		write(*,*) str1,str2
	enddo
	write(*,*) "***********************************"
	write(*,*) "***********************************"
	write(*,*) 
	stop 1
	end subroutine tilt
!
!
!--------------------------------------------------------------------------
!	Affichage
!
	subroutine banner(str1)
	implicit none
	character(len=*), intent(in)   :: str1
!
	write(*,*) 
	write(*,*) "**********************************************************************"
	write(*,*) "**********************************************************************"
	write(*,*) str1
	write(*,*) "**********************************************************************"
	write(*,*) "**********************************************************************"
	write(*,*) 
	end subroutine banner
!
!
!--------------------------------------------------------------------------
!	Affichage
!
	subroutine banner2(str1)
	implicit none
	character(len=*), intent(in)   :: str1
!
	write(*,*) 
	write(*,*) str1
	write(*,*) "******************************"
	write(*,*) 
	end subroutine banner2
!
!--------------------------------------------------------------------------
!
	end module m_cartes
	
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	routines de calcul pour la description du reseau reciproque
!	et des angles associes.
!
	module m_util2p1
	implicit none
	contains
!
!	------------------------------------------------------------
!	produit vectoriel w des vecteurs u et v de dimension 3
	subroutine vct(u,v,w)
	implicit none
	real, intent(in)  :: u(3),v(3)
	real, intent(out) :: w(3)

	w(1) = u(2)*v(3) - u(3)*v(2)
	w(2) = u(3)*v(1) - u(1)*v(3) 
	w(3) = u(1)*v(2) - u(2)*v(1)
	end subroutine vct
!
!	------------------------------------------------------------
!	produit vectoriel w des vecteurs u et v de dimension 3
!	cas de coordonnees complexes
	subroutine vctx(ur,ui,vr,vi,wr,wi)
	implicit none
	real, intent(in)  :: ur(3),ui(3),vr(3),vi(3)
	real, intent(out) :: wr(3),wi(3)

	wr(1) = (ur(2)*vr(3)-ui(2)*vi(3)) - (ur(3)*vr(2)-ui(3)*vi(2))
	wi(1) = (ur(2)*vi(3)+ui(2)*vr(3)) - (ur(3)*vi(2)+ui(3)*vr(2))
	wr(2) = (ur(3)*vr(1)-ui(3)*vi(1)) - (ur(1)*vr(3)-ui(1)*vi(3))
	wi(2) = (ur(3)*vi(1)+ui(3)*vr(1)) - (ur(1)*vi(3)+ui(1)*vr(3))
	wr(3) = (ur(1)*vr(2)-ui(1)*vi(2)) - (ur(2)*vr(1)-ui(2)*vi(1))
	wi(3) = (ur(1)*vi(2)+ui(1)*vr(2)) - (ur(2)*vi(1)+ui(2)*vr(1))
	end subroutine vctx
!
!	------------------------------------------------------------
!	determinant de 3 vecteurs a,b,c de dimension 3
!
	subroutine dtr(a,b,c,s)
	implicit none
	real, intent(in)  :: a(3),b(3),c(3)
	real, intent(out) :: s

	s =     a(1)*b(2)*c(3) + b(1)*c(2)*a(3) + c(1)*a(2)*b(3)
	s = s - a(3)*b(2)*c(1) - b(3)*c(2)*a(1) - c(3)*a(2)*b(1)

	end subroutine dtr
!
!	------------------------------------------------------------
!	calcul des donnees relatives au reseau reciproque 
!	definition du plan de diffusion par i,j et k
!
	subroutine reseau(ax,ay,az,alfa,beta,gama,vol,e,ee)
	implicit none
	real, intent(in)  :: ax,ay,az,alfa,beta,gama
	real, intent(out) :: vol,e(3,3),ee(3,3)
	integer           :: i,d
	real              :: pi
!
	d      = 3
	pi     = 4.0*atan(1.0)
!
!	definition des axes du reseau reel dans un repere othornorme
!	
	e      = 0.0
	e(1,1) = 1.0
	e(2,1) = cos(gama)
	e(2,2) = sin(gama)
	e(3,1) = cos(beta)
	e(3,2) = (cos(alfa)-e(3,1)*e(2,1))/e(2,2)
	e(3,3) = sqrt(1.0-e(3,1)**2-e(3,2)**2)
	do i=1,d
		e(i,:) = e(i,:)/sqrt(sum(e(i,:)*e(i,:)))
	enddo
	e(1,:) = e(1,:)*ax
	e(2,:) = e(2,:)*ay
	e(3,:) = e(3,:)*az
!	
!	definition des axes du reseau reciproque dans le meme repere
!	
	call dtr(e(1,:),e(2,:),e(3,:),vol)
	call vct(e(1,:),e(2,:),ee(3,:))
	call vct(e(2,:),e(3,:),ee(1,:))
	call vct(e(3,:),e(1,:),ee(2,:))
	ee = 2.0*pi/vol*ee

	end subroutine reseau
!
!	------------------------------------------------------------
!	rotation formule de Rodrigues 
!	rotation dans le plan perpendiculaire au vecteur n 
!	d'un angle a. le vecteur resultat est dans u
!
	subroutine rodrigue(a,n,u,v)
	implicit none
	real, intent(in)  :: a,n(3),u(3)
	real, intent(out) :: v(3)
	real              :: m(3,3),mu(3,3)
	integer           :: i,j,d
!
	d = 3
	m = 0.0
	do i=1,d
		m(i,i) = cos(a)
	enddo
	do i=1,d
	do j=1,d
		m(i,j) = m(i,j) + (1.0-cos(a))*n(i)*n(j)
	enddo
	enddo
	mu = 0.0
	mu(1,2) = -n(3)
	mu(2,1) = +n(3)
	mu(1,3) = +n(2)
	mu(3,1) = -n(2)
	mu(2,3) = -n(1)
	mu(3,2) = +n(1)
	
	v = matmul(m+sin(a)*mu,u)
	end subroutine rodrigue
!
!	------------------------------------------------------------
!	matrice formule de Rodrigues
!	rotation dans le plan perpendiculaire au vecteur n 
!
	subroutine rodrigue_mat(n,nn,rr,ri)
	implicit none
	real, intent(in)  :: n(3)
	real, intent(out) :: nn(3,3),rr(3,3),ri(3,3)
	integer           :: i,j,d
!
	d  = 3
	rr = 0.0
	do i=1,d
		rr(i,i) = 1.0
		do j=1,d
			nn(i,j) = n(i)*n(j)
		enddo
	enddo
	rr = rr - nn
	ri = 0.0
	ri(1,2) = -n(3)
	ri(2,1) = +n(3)
	ri(1,3) = +n(2)
	ri(3,1) = -n(2)
	ri(2,3) = -n(1)
	ri(3,2) = +n(1)
	end subroutine rodrigue_mat
!
!	------------------------------------------------------------
!	caracteristiques des RE
!	
	subroutine paramRE(nom,j,gj,alf,bet,gam)
	implicit none
	character(len=8), intent(in)  :: nom
	real,             intent(out) :: j,gj,alf,bet,gam

	select case(nom)

	case("ND")
		j   = 9.0/2.0
		gj  = 8.0/11.0
		alf = -7.0/(9*121)
		bet = -8*17.0/(3*3*3*11*11*11*13)
		gam = -5*17*19.0/(3*3*3*7*11*11*11*13*13)
	case("ER")
		j   = 15./2.
		gj  = 6./5.
		alf = 4./(9*25*7)
		bet = 2./(9*5*7*11*13)
		gam = 8./(27*7*121*169)
	case("YB") 
		j   = 7./2.
		gj  = 8./7.
		alf = 2./63
		bet =-2./(9*5*7*11)
		gam = 4./(27*7*11*13)
	case("HO") 
		j   = 8.
		gj  = 5./4.
		alf =-1./(2*9*25)
		bet =-1./(2*3*5*7*11*13)
		gam =-5./(27*7*121*169)
	case("TB") 
		j   = 6.
		gj  = 3./2.
		alf =-1./99
		bet = 2./(27*5*121)
		gam =-1./(81*7*121*13)
	case("PR") 
		j   = 4.
		gj  = 4./5.
		alf =-4.*13./(9*25*11)
		bet =-4./(9*5*121)
		gam = 16.*17./(81*5*7*121*13)
	case("DY") 
		j   = 15./2.
		gj  = 4./3.
		alf =-2./(9*5*7)
		bet =-8./(27*5*7*11*13)
		gam = 4./(27*7*121*169)
	case("CE") 
		j   = 5./2.
		gj  = 6./7.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S1") 
		j   = 1.
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S2") 
		j   = 2.
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("SD2") 
		j   = 0.5
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S3D2") 
		j   = 1.5
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S5D2") 
		j   = 2.5
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S7D2") 
		j   = 3.5
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("MN3") 
		j   = 2.
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("MN4") 
		j   = 3./2.
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("FE3") 
		j   = 5./2.
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S3") 
		j   = 3
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S4") 
		j   = 4
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S5") 
		j   = 5
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S6") 
		j   = 6
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S7") 
		j   = 7
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S8") 
		j   = 8
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S9") 
		j   = 9
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case("S10") 
		j   = 10
		gj  = 2.
		alf = 1.
		bet = 1.
		gam = 1.
	case default
	end select
!
	end subroutine paramRE
!
!	------------------------------------------------------------
!	facteur de forme des RE3+
!	
	subroutine formRE(nom,fga,fa,fgb,fb,fgc,fc,fgd)
	implicit none
	character(len=8), intent(in)  :: nom
	real,             intent(out) :: fga,fa,fgb,fb,fgc,fc,fgd

	select case(nom)
	case("ER")
		fga = 0.0586
		fa  = 17.9802
		fgb = 0.3540
		fb  = 7.0964
		fgc = 0.6126
		fc  = 2.7482
		fgd = -0.0251
	case("YB") 
		fga = 0.0416
		fa  = 16.0949
		fgb = 0.2849
		fb  = 7.8341
		fgc = 0.6961
		fc  = 2.6725
		fgd = -0.0229
	case("HO") 
		fga = 0.0566
		fa  = 18.3176
		fgb = 0.3365
		fb  = 7.6880
		fgc = 0.6317
		fc  = 2.9427
		fgd = -0.0248
	case("TB") 
		fga = 0.0177
		fa  = 25.5095
		fgb = 0.2921
		fb  = 10.5769
		fgc = 0.7133
		fc  = 3.5122
		fgd = -0.0231
	case("PR") 
		fga = 0.0504
		fa  = 16.0949
		fgb = 0.2849
		fb  = 7.8341
		fgc = 0.6961
		fc  = 2.6725
		fgd = -0.0218
	case("DY") 
		fga = 0.1157
		fa  = 15.0732
		fgb = 0.3270
		fb  = 6.7991
		fgc = 0.5821
		fc  = 3.0202
		fgd = -0.0249
	case("ZZ") 
		fga = 1.0
		fa  = 0.0
		fgb = 0.0
		fb  = 0.0
		fgc = 0.0
		fc  = 0.0
		fgd = 0.0
	case default
		fga = 1.0
		fa  = 0.0
		fgb = 0.0
		fb  = 0.0
		fgc = 0.0
		fc  = 0.0
		fgd = 0.0
	end select
!
	end subroutine formRE
!
!	------------------------------------------------------------
	end module m_util2p1
	
	!-------------------------------------------------------------------------------
!	methodes d'algebre lineaire pompees dans Numerical Recipes
!	pour la resolution de systemes lineaires
!
	module m_solve
	implicit none
	contains
!
!	-------------------------------------------------------------
!	Transposition d'une matrice a(n,m)
!
	subroutine trsp(n,m,a,ta)
	implicit none
	integer, intent(in)    :: n,m
	real,    intent(in)    :: a(:,:)
	real,    intent(out)   :: ta(:,:)
	integer                :: i,j
!	
	do i=1,n
	do j=1,m
		ta(j,i) = a(i,j)
	enddo
	enddo
!
	end subroutine trsp
!
!	-------------------------------------------------------------
!	Decomposition LU
!	index est un vecteur nnant les permutations effectuees 
!	sur les lignes pendant la decomposition
!	d vaut + ou - 1, suivant que le nombre de permutations de 
!	lignes a ete pair ou impair.
!
	subroutine ludcmp(n,a,index,d)
	implicit none
	integer, intent(in)    :: n
	real,    intent(inout) :: a(:,:)
	integer, intent(out)   :: index(:)
	real,    intent(out)   :: d
!
	integer                :: i,imax,j,k
	real                   :: aamax,dum,sum,tiny
	real, pointer          :: vv(:)
!
	allocate(vv(n))
!
	tiny = 1.0e-20
	d    = 1.0
    do i = 1,n 
		aamax = 0.0
		do j = 1,n 
			if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
	 	enddo
	  	if(aamax.eq.0) then
			write(*,*) "matrice singuliere dans ludcmp"
	       		stop 1
		endif
		vv(i) = 1.0/aamax
	enddo
!
	do j = 1,n 
		do i = 1,j-1 
			sum = a(i,j)
			do k = 1,i-1 
				sum = sum-a(i,k)*a(k,j)
			enddo
			a(i,j) = sum
		enddo

		aamax = 0.0
		do i = j,n 
			sum = a(i,j)
			do k = 1,j-1 
				sum = sum-a(i,k)*a(k,j)
			enddo
		    a(i,j) = sum
			dum    = vv(i)*abs(sum)
			if(dum.ge.aamax) then
				imax  = i
				aamax = dum
			endif
		enddo
!
	  	if(j.ne.imax) then
			do k = 1,n 
				dum       = a(imax,k)
				a(imax,k) = a(j,k)
				a(j,k)    = dum
			enddo
			d = -d
			vv(imax) = vv(j)
	  	endif
		index(j) = imax
		if(a(j,j).eq.0.) a(j,j) = tiny
		if(j.ne.n) then
			dum = 1.0/a(j,j)
			do i = j+1,n 
				a(i,j) =a(i,j)*dum
			enddo
		endif
	enddo
	deallocate(vv)
	end subroutine ludcmp
!
!	-------------------------------------------------------------
!	Phase de remontee associee a la decomposition LU
!	Resout le systeme de n equations lineaires :
!	ax=b ou a est sous decomposee LU par ludcmp.
!	b contient le second membre en entree et
!	la solution x en sortie.
!
	subroutine lubksb(n,a,index,b)
  	implicit none
	integer, intent(in)    :: n
	real,    intent(in)    :: a(:,:)
	integer, intent(in)    :: index(:)
	real,    intent(inout) :: b(:)
!
	integer                :: i,ii,j,ll
	real                   :: sum

	ii = 0
	do i = 1,n 
		ll    = index(i)
		sum   = b(ll)
		b(ll) = b(i)
		if(ii.ne.0) then
			do j = ii,i-1 
				sum = sum-a(i,j)*b(j)
			enddo
		elseif(sum.ne.0) then 
			ii = i
		endif
		b(i) = sum
	enddo
	do i = n,1,-1
		sum = b(i)
		do j = i+1,n 
			sum = sum-a(i,j)*b(j)
		enddo
		b(i) = sum/a(i,i)
	enddo
	end subroutine lubksb
!
!	-------------------------------------------------------------
!	resolution systeme lineaire  ax = b
!
	subroutine solve(n,a,b,x)
	implicit none
	integer,    intent(in)  :: n
	real,       intent(in)  :: a(:,:),b(:)
	real,       intent(out) :: x(:)
!
	integer, pointer        :: index(:)
	real                    :: d     
	real, pointer           :: a2(:,:),b2(:)   

	allocate(index(n),a2(n,n),b2(n))
	a2 = a
	b2 = b
 	call ludcmp(n,a2,index,d)
	call lubksb(n,a2,index,b2)
	x = b2
	deallocate(index,a2,b2)
	end subroutine solve
!
!	-------------------------------------------------------------
!	resolution systeme lineaire  ax = b pour a et b complexes
!
	subroutine solvecplx(n,ar,ai,br,bi,xr,xi)
	implicit none
	integer,    intent(in)  :: n
	real,       intent(in)  :: ar(:,:),ai(:,:),br(:),bi(:)
	real,       intent(out) :: xr(:),xi(:)
!
	integer                 :: n2
	integer, pointer        :: index(:)
	real                    :: d     
	real, pointer           :: a2(:,:),b2(:),x2(:)

	n2 = 2*n
	allocate(index(n2),a2(n2,n2),b2(n2))
	a2(  1:n ,  1:n ) = +ar
	a2(  1:n ,n+1:n2) = -ai
	a2(n+1:n2,  1:n ) = +ai
	a2(n+1:n2,n+1:n2) = +ar
	b2(  1:n )        = br
	b2(n+1:n2)        = bi
!
 	call ludcmp(n2,a2,index,d)
	call lubksb(n2,a2,index,b2)
	xr = b2(  1:n )
	xi = b2(n+1:n2)
	deallocate(index,a2,b2)
	end subroutine solvecplx
!
!	-------------------------------------------------------------
!	calcul de l'inverse de la matrice a
!
	subroutine inverse(n,a,y)
	implicit none
	integer,    intent(in)  :: n
	real,       intent(in)  :: a(:,:)
	real,       intent(out) :: y(:,:)
!
	integer, pointer        :: index(:)
	integer                 :: i
	real                    :: d         
	real, pointer           :: a2(:,:)
!           
	y = 0.0
	do i = 1,n
		y(i,i) = 1.0
	enddo
!
 	allocate(index(n),a2(n,n))
	a2 = a
	call ludcmp(n,a2,index,d)
	do i = 1,n
		call lubksb(n,a2,index,y(:,i))
	enddo
	deallocate(index,a2)
	end subroutine inverse
!
!	-------------------------------------------------------------
!	calcul du determinant de la matrice a
!
	subroutine det(n,a,d)
	implicit none
	integer,    intent(in)  :: n
	real,       intent(in)  :: a(:,:)
	real,       intent(out) :: d
!
	integer                 :: i
	integer, pointer        :: index(:)
	real,    pointer        :: a2(:,:)
!           
 	allocate(index(n),a2(n,n))
	a2 = a
	call ludcmp(n,a2,index,d)
	do i = 1,n
		d = d * a2(i,i)
	enddo
	deallocate(index,a2)
	end subroutine det

!	-------------------------------------------------------------
!
	end module m_solve
	
	!-------------------------------------------------------------------------------
!	methodes d'algebre lineaire pompees dans Numerical Recipes
!	pour la resolution de systemes lineaires
!	modifs pour traitement des nombres complexes
!
	module m_solvex
	implicit none
	contains
!
!	-------------------------------------------------------------
!	Transpose-conjugue d'une matrice a(n,m)
!
	subroutine trspx(n,m,ar,ai,tar,tai)
	implicit none
	integer, intent(in)    :: n,m
	real,    intent(in)    :: ar(:,:),ai(:,:)
	real,    intent(out)   :: tar(:,:),tai(:,:)
	integer                :: i,j
!	
	do i=1,n
	do j=1,m
		tar(j,i) =  ar(i,j)
		tai(j,i) = -ai(i,j)
	enddo
	enddo
!
	end subroutine trspx
!
!	-------------------------------------------------------------
!	Decomposition LU
!	index est un vecteur donnant les permutations effectuees 
!	sur les lignes pendant la decomposition
!	d vaut + ou - 1, suivant que le nombre de permutations de 
!	lignes a ete pair ou impair.
!
	subroutine ludcmpx(n,ar,ai,index,d)
	implicit none
	integer, intent(in)    :: n
	real,    intent(inout) :: ar(:,:),ai(:,:)
	integer, intent(out)   :: index(:)
	real,    intent(out)   :: d
!
	integer                :: i,imax,j,k
	real                   :: aamax,tiny,v,dum
	real                   :: dumr,dumi,sumr,sumi,sr,si
	real, pointer          :: vv(:)
!
	allocate(vv(n))
!
	tiny = 1.0e-20
	d    = 1.0
    	do i = 1,n 
		aamax = 0.0
		do j = 1,n 
			v = sqrt(ar(i,j)*ar(i,j) + ai(i,j)*ai(i,j))
			if (v.gt.aamax) aamax = v
	 	enddo
	  	if(aamax.eq.0) then
			write(*,*) "matrice singuliere dans ludcmp"
	       		stop 1
		endif
		vv(i) = 1.0/aamax
	enddo
!
	do j = 1,n 
		do i = 1,j-1 
			sumr = ar(i,j)
			sumi = ai(i,j)
			do k = 1,i-1 
				sumr = sumr-(ar(i,k)*ar(k,j)-ai(i,k)*ai(k,j))
				sumi = sumi-(ar(i,k)*ai(k,j)+ai(i,k)*ar(k,j))
			enddo
			ar(i,j) = sumr
			ai(i,j) = sumi
		enddo

		aamax = 0.0
		do i = j,n 
			sumr = ar(i,j)
			sumi = ai(i,j)
			do k = 1,j-1 
				sumr = sumr-(ar(i,k)*ar(k,j)-ai(i,k)*ai(k,j))
				sumi = sumi-(ar(i,k)*ai(k,j)+ai(i,k)*ar(k,j))
			enddo
		    ar(i,j) = sumr
			ai(i,j) = sumi
			dum     = vv(i)*sqrt(sumr*sumr+sumi*sumi)
			if(dum.ge.aamax) then
				imax  = i
				aamax = dum 
			endif
		enddo
!
	  	if(j.ne.imax) then
			do k = 1,n 
				dumr       = ar(imax,k)
				dumi       = ai(imax,k)
				ar(imax,k) = ar(j,k)
				ai(imax,k) = ai(j,k)
				ar(j,k)    = dumr
				ai(j,k)    = dumi
			enddo
			d = -d
			vv(imax) = vv(j)
	  	endif
		index(j) = imax
		v = ar(j,j)*ar(j,j)+ai(j,j)*ai(j,j)
		if(v.eq.0.) then
			ar(j,j) = tiny
			ai(j,j) = tiny
		endif
		if(j.ne.n) then
			v    = ar(j,j)*ar(j,j)+ai(j,j)*ai(j,j)
			dumr = ar(j,j)/v
			dumi =-ai(j,j)/v
			do i = j+1,n 
				sr = ar(i,j)
				si = ai(i,j)
				ar(i,j) = sr*dumr-si*dumi
				ai(i,j) = sr*dumi+si*dumr
			enddo
		endif
	enddo
	deallocate(vv)
	end subroutine ludcmpx
!
!	-------------------------------------------------------------
!	Phase de remontee associee a la decomposition LU
!	Resout le systeme de n equations lineaires :
!	ax=b ou a est sous decomposee LU par ludcmp.
!	b contient le second membre en entree et
!	la solution x en sortie.
!
	subroutine lubksbx(n,ar,ai,index,br,bi)
  	implicit none
	integer, intent(in)    :: n
	real,    intent(in)    :: ar(:,:),ai(:,:)
	integer, intent(in)    :: index(:)
	real,    intent(inout) :: br(:),bi(:)
!
	integer                :: i,ii,j,ll
	real                   :: sumr,sumi,v,vr,vi

	ii = 0
	do i = 1,n 
		ll     = index(i)
		sumr   = br(ll)
		sumi   = bi(ll)
		br(ll) = br(i)
		bi(ll) = bi(i)
		if(ii.ne.0) then
			do j = ii,i-1 
				sumr = sumr-(ar(i,j)*br(j)-ai(i,j)*bi(j))
				sumi = sumi-(ar(i,j)*bi(j)+ai(i,j)*br(j))
			enddo
		else
			v = sumr*sumr+sumi*sumi
			if(v.ne.0) ii = i
		endif
		br(i) = sumr
		bi(i) = sumi
	enddo
	do i = n,1,-1
		sumr = br(i)
		sumi = bi(i)
		do j = i+1,n 
			sumr = sumr-(ar(i,j)*br(j)-ai(i,j)*bi(j))
			sumi = sumi-(ar(i,j)*bi(j)+ai(i,j)*br(j))
		enddo
		v     =  ar(i,i)*ar(i,i)+ai(i,i)*ai(i,i)
		vr    =  ar(i,i)/v
		vi    = -ai(i,i)/v
		br(i) = sumr*vr-sumi*vi
		bi(i) = sumr*vi+sumi*vr
	enddo
	end subroutine lubksbx
!
!	-------------------------------------------------------------
!	resolution systeme lineaire  ax = b
!
	subroutine solvex(n,ar,ai,br,bi,xr,xi)
	implicit none
	integer,    intent(in)  :: n
	real,       intent(in)  :: ar(:,:),ai(:,:),br(:),bi(:)
	real,       intent(out) :: xr(:),xi(:)
!
	integer, pointer        :: index(:)
	real                    :: d     
	real, pointer           :: a2r(:,:),a2i(:,:),b2r(:),b2i(:)   

	allocate(index(n),a2r(n,n),a2i(n,n),b2r(n),b2i(n))
	a2r = ar
	a2i = ai
	b2r = br
	b2i = bi
 	call ludcmpx(n,a2r,a2i,index,d)
	call lubksbx(n,a2r,a2i,index,b2r,b2i)
	xr = b2r
	xi = b2i
	deallocate(index,a2r,a2i,b2r,b2i)
	end subroutine solvex
!
!	-------------------------------------------------------------
!	multiplication matrice  complexe 
!
	subroutine matmultcplx(ar,ai,br,bi,cr,ci)
	implicit none
	real,       intent(in)  :: ar(:,:),ai(:,:),br(:,:),bi(:,:)
	real,       intent(out) :: cr(:,:),ci(:,:)
!
	cr = matmul(ar,br) - matmul(ai,bi)
	ci = matmul(ar,bi) + matmul(ai,br)
	
	end subroutine matmultcplx
!
!	-------------------------------------------------------------
!	multiplication matrice par vecteur complexe 
!
	subroutine multcplx(ar,ai,xr,xi,br,bi)
	implicit none
	real,       intent(in)  :: ar(:,:),ai(:,:),xr(:),xi(:)
	real,       intent(out) :: br(:),bi(:)
!
	br = matmul(ar,xr) - matmul(ai,xi)
	bi = matmul(ar,xi) + matmul(ai,xr)
	
	end subroutine multcplx
!
!	-------------------------------------------------------------
!	produit scalaire u A v avec matrice et vecteur complexe 
!
	subroutine pdtcplx(n,ur,ui,ar,ai,vr,vi,xr,xi)
	implicit none
	integer,    intent(in)  :: n
	real,       intent(in)  :: ur(:),ui(:),ar(:,:),ai(:,:),vr(:),vi(:)
	real,       intent(out) :: xr,xi
!
	real, pointer           :: w1(:),w2(:)
!
	allocate(w1(n),w2(n))
	w1 = matmul(ar,vr) - matmul(ai,vi)
	w2 = matmul(ar,vi) + matmul(ai,vr)
	xr = sum(ur*w1) - sum(ui*w2)
	xi = sum(ur*w2) + sum(ui*w1)
	deallocate(w1,w2)
	end subroutine pdtcplx
!
!	-------------------------------------------------------------
!	calcul de l'inverse de la matrice a
!
	subroutine inversex(n,ar,ai,yr,yi)
	implicit none
	integer,    intent(in)  :: n
	real,       intent(in)  :: ar(:,:),ai(:,:)
	real,       intent(out) :: yr(:,:),yi(:,:)
!
	integer, pointer        :: index(:)
	integer                 :: i,j
	real                    :: d         
	real, pointer           :: a2r(:,:),a2i(:,:)
!           
	yr = 0.0
	yi = 0.0
	do i = 1,n
		yr(i,i) = 1.0
	enddo
!
 	allocate(index(n),a2r(n,n),a2i(n,n))
	a2r = ar
	a2i = ai
	call ludcmpx(n,a2r,a2i,index,d)
	do i=1,n
		write(*,"(10f8.3)") (a2r(i,j),a2i(i,j),j=1,n)
	enddo
	do i = 1,n
		call lubksbx(n,a2r,a2i,index,yr(:,i),yi(:,i))
	enddo
	deallocate(index,a2r,a2i)
	end subroutine inversex
!
!	-------------------------------------------------------------
!	calcul du determinant de la matrice a
!
	subroutine detx(n,ar,ai,dr,di)
	implicit none
	integer,    intent(in)  :: n
	real,       intent(in)  :: ar(:,:),ai(:,:)
	real,       intent(out) :: dr,di
!
	integer                 :: i,j
	integer, pointer        :: index(:)
	real,    pointer        :: a2r(:,:),a2i(:,:)
	real                    :: sr,si
!           
 	allocate(index(n),a2r(n,n),a2i(n,n))
	a2r = ar
	a2i = ai
	dr  = 1.0
	di  = 0.0
	call ludcmpx(n,a2r,a2i,index,dr)
	do i = 1,n
		sr = dr*a2r(i,i) - di*a2i(i,i)
		si = dr*a2i(i,i) + di*a2r(i,i)
		dr = sr
		di = si
	enddo
	deallocate(index,a2r,a2i)
	end subroutine detx
!
!	-------------------------------------------------------------
!
	end module m_solvex
	
	!-------------------------------------------------------------------------------
!	methodes d'algebre lineaire pompees dans Numerical Recipes
!	pour la diagonalisation de matrices quelconques 
!
	module m_diag
	implicit none
	contains
!	
!-------------------------------------------------------------------------
!
!	'balancing' d'une matrice a , afin eviter les pbs d'arrondis
!	les valeurs propres ne sont pas modifiees
!	une matrice symetrique n'est pas modifiee
!
	subroutine balanc(n,a)
	implicit none
	integer, intent(in)    :: n
	real,    intent(inout) :: a(:,:)
	integer                :: i,last,ndum
	real                   :: c,f,g,r,s
	real                   :: radx,sqradx
!	
	radx   = radix(a)
	sqradx = radx**2
!
	do
		last = 1
		do i = 1,n
			c = sum(abs(a(:,i)))-a(i,i)
			r = sum(abs(a(i,:)))-a(i,i)
			if ((c .ne. 0.0).and.(r .ne. 0.0)) then
				g = r/radx
				f = 1.0
				s = c+r
				do
					if (c.ge.g) exit
					f = f*radx
					c = c*sqradx
				enddo
				g=r*radx
				do
					if (c.le.g) exit
					f = f/radx
					c = c/sqradx
				enddo
				if ((c+r)/f .lt. 0.95*s) then
					last   = 0
					g      = 1.0/f
					a(i,:) = a(i,:)*g
					a(:,i) = a(:,i)*f
				endif
			end if
		enddo
		if (last.ne.0) exit
	enddo 
	end subroutine balanc
!
!-------------------------------------------------------------------------
!
!	reduction de Hessemberg par la methode dite de l'elimination
!	en sortie, les elements a(i,j) tels que i > j+1 sont nuls
!	les valeurs propres ne sont pas modifiees
!
	subroutine elmhes(n,a)
	implicit none
	integer, intent(in)   :: n
	real,	intent(inout) :: a(:,:)
!
	integer               :: i,j,m
	real                  :: x
	real, pointer         :: y(:),u(:,:)
!
	allocate(y(n),u(n,n))
	do m = 2,n-1 
!
!		recherche du pivot x et de la ligne correspondante i
!
		x = 0.0
		i = m
		do j = m,n
			if(abs(a(j,m-1)).gt.abs(x)) then
				x = a(j,m-1)
				i = j
			endif
		enddo
!
!		echange des lignes puis des colonnes
!
		if (i.ne.m) then
			y(m-1:n)   = a(i,m-1:n)
			a(i,m-1:n) = a(m,m-1:n)
			a(m,m-1:n) = y(m-1:n)
!
			y(:)       = a(:,i)
			a(:,i)     = a(:,m)
			a(:,m)     = y(:)
		endif
!
!		elimination
!
		if (x.ne.0.0) then
		 	y(m+1:n)     = a(m+1:n,m-1)/x 
			a(m+1:n,m-1) = y(m+1:n)
			do i = m+1,n
			do j = m,n
				u(i,j) = y(i)*a(m,j)
			enddo
			enddo
			a(m+1:n,m:n) = a(m+1:n,m:n)-u(m+1:n,m:n)
			a(:,m)       = a(:,m)+matmul(a(:,m+1:n),y(m+1:n))
		endif
	enddo
	deallocate(y,u)
	end subroutine elmhes
!	
!-------------------------------------------------------------------------
!
!	algorithme QR pour matrice de Hessemberg
!
	subroutine hqr(n,a,wr,wi,marche)
	implicit none
	integer, intent(in)    :: n
	real,    intent(inout) :: a(:,:)
	real,    intent(inout) :: wr(:),wi(:)
	logical, intent(out)   :: marche
!
	integer                :: i,its,j,k,l,m,nn,mnnk
	integer                :: maxits
	real                   :: anorm,p,q,r,s,t,u,v,w,x,y,z
	real, pointer          :: pp(:)
!
	allocate(pp(n))
!
	marche = .true.
	maxits = 500
	anorm  = 0.0
	do i = 1,n
	do j = max(i-1,1),n
		anorm = anorm + abs(a(i,j))
	enddo
	enddo
!
	t  = 0.0
	nn = n
!
	do
	if (nn.lt.1) exit
	its = 0
	iterate: do
		small: do l = nn,2,-1
				s = abs(a(l-1,l-1))+abs(a(l,l))
				if (s.eq.0.0) s=anorm
				if (abs(a(l,l-1))+s.eq.s) then
					a(l,l-1) = 0.0
					exit small
				endif
			end do small
		x = a(nn,nn)
		if (l.eq.nn) then
			wr(nn) = x+t
			wi(nn) = 0.0
			nn     = nn-1
			exit iterate
		endif
!
		y = a(nn-1,nn-1)
		w = a(nn,nn-1)*a(nn-1,nn)
		if (l.eq.nn-1) then
			p = 0.5*(y-x)
			q = p**2+w
			z = sqrt(abs(q))
			x = x+t
			if (q.ge.0.0) then
				z        = p+sign(z,p)
				wr(nn)   = x+z
				wr(nn-1) = wr(nn)
				if (z.ne.0.0) wr(nn) = x-w/z
				wi(nn)   = 0.0
				wi(nn-1) = 0.0
			else
				wr(nn)   = x+p
				wr(nn-1) = wr(nn)
				wi(nn)   = z
				wi(nn-1) = -z
			endif
			nn = nn-2
			exit iterate
		endif
!
		if (its == maxits) then
			write(*,*) "Trop d iterations dans hqr"
			marche = .false.
			return
		endif
		if (its == 10 .or. its == 20) then
			t = t+x
			do i=1,nn
				a(i,i) = a(i,i) -x
			enddo
			s = abs(a(nn,nn-1)) + abs(a(nn-1,nn-2))
			x = 0.75*s
			y = x
			w = -0.4375*s**2
		endif
		its = its+1
		do m = nn-2,l,-1
			z = a(m,m)
			r = x-z
			s = y-z
			p = (r*s-w)/a(m+1,m)+a(m,m+1)
			q = a(m+1,m+1)-z-r-s
			r = a(m+2,m+1)
			s = abs(p)+abs(q)+abs(r)
			p = p/s
			q = q/s
			r = r/s
			if (m.eq.l) exit
			u = abs(a(m,m-1))*(abs(q)+abs(r))
			v = abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
			if (u+v.eq.v) exit
		enddo
		do i = m+2,nn
			a(i,i-2) = 0.0
			if (i.ne.m+2) a(i,i-3) = 0.0
		enddo
		do k = m,nn-1
			if (k.ne.m) then
				p = a(k,k-1)
				q = a(k+1,k-1)
				r = 0.0
				if (k.ne.nn-1) r = a(k+2,k-1)
				x = abs(p)+abs(q)+abs(r)
				if (x.ne.0.0) then 
					p = p/x
					q = q/x
					r = r/x
				endif
			endif
			s = sign(sqrt(p**2+q**2+r**2),p)
			if (s.ne.0.0) then 
				if (k.eq.m) then
					if (l.ne.m) a(k,k-1) = -a(k,k-1)
				else
					a(k,k-1) = -s*x
				endif
				p = p+s
				x = p/s
				y = q/s
				z = r/s
				q = q/p
				r = r/p
				pp(k:nn) = a(k,k:nn)+q*a(k+1,k:nn)
				if (k.ne.nn-1) then
				pp(k:nn)    = pp(k:nn)+r*a(k+2,k:nn)
				a(k+2,k:nn) = a(k+2,k:nn)-pp(k:nn)*z
				endif
				a(k+1,k:nn) = a(k+1,k:nn)-pp(k:nn)*y
				a(k,k:nn)   = a(k,k:nn)-pp(k:nn)*x
				mnnk        = min(nn,k+3)
				pp(l:mnnk)  = x*a(l:mnnk,k)+y*a(l:mnnk,k+1)
				if (k.ne.nn-1) then
				pp(l:mnnk)    = pp(l:mnnk)+z*a(l:mnnk,k+2)
				a(l:mnnk,k+2) = a(l:mnnk,k+2)-pp(l:mnnk)*r
				endif
				a(l:mnnk,k+1) = a(l:mnnk,k+1)-pp(l:mnnk)*q
				a(l:mnnk,k)   = a(l:mnnk,k)-pp(l:mnnk)
			endif
		enddo
	end do iterate
	enddo
	deallocate(pp)
	end subroutine hqr
!
!-------------------------------------------------------------------------------
!	tri des valeurs d'un tableau dans l'ordre croissant
!
	subroutine valpsrt(n,dr,di)
	implicit none
	integer, intent(in)    :: n
	real,    intent(inout) :: dr(:),di(:)
!
	integer                :: i,j,k
	real                   :: pr,pi
!
	do i = 1,n-1
		k = i
		pr = dr(i)
		pi = di(i)
		do j = i+1,n
			if(dr(j).le.pr) then
				k = j
				pr = dr(j)
				pi = di(j)
			endif
		enddo
		if(k.ne.i) then
			dr(k) = dr(i)
			dr(i) = pr
			di(k) = di(i)
			di(i) = pi
		endif
	enddo
!
	end subroutine valpsrt
!	
!-------------------------------------------------------------------------
!
!	diagonalisation d'une matrice a coeffcients reels
!	les valeurs propres sont rangees dans l'ordre croissant
!	de leur partie reelle
!
	subroutine diago(n,a,wr,wi,success)
	implicit none
	integer, intent(in)  :: n
	real,    intent(in)  :: a(:,:)
	real,    intent(out) :: wr(:),wi(:)
	logical, intent(out) :: success
!
	real, pointer        :: a2(:,:)
	logical              :: marche
!
	allocate(a2(n,n))
	a2 = a
	call balanc(n,a2)
	call elmhes(n,a2)
	call hqr(n,a2,wr,wi,marche)
	if(marche) then
		call valpsrt(n,wr,wi)
		success = .true.
	else
		success = .false.
	endif
	deallocate(a2)
!
	end subroutine diago
!	
!-------------------------------------------------------------------------
!	diagonalisation d'une matrice a coeffcients complexes
!	les valeurs propres sont rangees dans l'ordre croissant
!	de leur partie reelle
!
	subroutine diagocplx(n,ar,ai,wr,wi,success)
	implicit none
	integer, intent(in)  :: n
	real,    intent(in)  :: ar(:,:),ai(:,:)
	real,    intent(out) :: wr(:),wi(:)
	logical, intent(out) :: success
!
	integer              :: n2,i,j
	real                 :: eps
	real, pointer        :: a2(:,:),wr2(:),wi2(:)
	logical              :: marche
!
	eps = 0.05
	n2 = 2*n
	allocate(a2(n2,n2),wr2(n2),wi2(n2))
	a2(  1:n ,  1:n ) = +ar
	a2(  1:n ,n+1:n2) = -ai
	a2(n+1:n2,  1:n ) = +ai
	a2(n+1:n2,n+1:n2) = +ar
	call balanc(n2,a2)
	call elmhes(n2,a2)
	call hqr(n2,a2,wr2,wi2,marche)
	if(marche) then
		call valpsrt(n2,wr2,wi2)
		do i = 1,n2
		if(abs(wi2(i)).ge.eps) then
		write(*,*) "--------------------------------------"
		write(*,*) "En sortie de la routine diagocplx ...."
		write(*,*) "Attention au choix des valeurs propres"
		write(*,"(i4,2f12.5)") i,wr2(i),wi2(i)
		write(*,*) "--------------------------------------"
		endif
		enddo
!
!		on suppose que les valeurs propres sont complexes-conjuguees
!		les unes des autres et on n'en prend qu'une sur deux. 
!		valable que lorsque ces valeurs propres sont reelles.
!
		do i = 1,n
			j     = 2*(i-1)+1
			wr(i) = wr2(j)
			wi(i) = wi2(j)
		enddo
		success = .true.
	else
		success = .false.
	endif
	deallocate(a2,wr2,wi2)
!
	end subroutine diagocplx
!	
!-------------------------------------------------------------------------
!
	end module m_diag
	
	!-------------------------------------------------------------------------------
!	methodes d'algebre lineaire pompees dans Numerical Recipes
!	pour la diagonalisation de matrices symetriques reelles 
!	ou hermitiennes
!
	module m_diag_sym
	implicit none
	contains
!
!-------------------------------------------------------------------------------
!	diagonalisation selon la methode de Jacobi
!	
	subroutine jacobi(n,a_in,d,v,nrot)
	implicit none
	integer, intent(in)  :: n
	real,    intent(in)  :: a_in(:,:)
	real,    intent(out) :: d(:),v(:,:)
	integer, intent(out) :: nrot
!
	integer              :: i,ip,iq,j
	real                 :: c,g,h,s,sm,t,tau,theta,tresh
	real, pointer        :: a(:,:),b(:),z(:)
!
	allocate(a(n,n),b(n),z(n))
!
!	------------------------------------------------------------------------
!	initialisations
!	v      ->  matrice identite
!	b et d -> diagonale de a 
!
	a = a_in
	v = 0.0
	do ip = 1,n
		v(ip,ip) = 1.0
		b(ip)    = a(ip,ip)
	enddo
	d = b
	z = 0.0
	nrot = 0
!
	i = 0
	do
!		mise a jour du nombre d'iterations
		i  = i + 1
!
!		somme des elements non diagonaux
		sm = 0.0
		do ip = 1   ,n-1
		do iq = ip+1,n
			sm = sm  + abs(a(ip,iq))
		enddo
		enddo
!		exit normal de la boucle precision machine atteinte
		if(sm.eq.0) exit
!
!		choix de la valeur de tresh
		tresh  = 0.0
		if(i.lt.4) tresh = 0.2*sm/n**2
!
!		rotations
		do ip = 1   ,n-1
		do iq = ip+1,n
			g = 100.0*abs(a(ip,iq))
			if((i.gt.4) & 
			   .and.(abs(d(ip))+g.eq.abs(d(ip))) &
			   .and.(abs(d(iq))+g.eq.abs(d(iq)))) then
				a(ip,iq) = 0.0
			elseif(abs(a(ip,iq)).gt.tresh) then
				h = d(iq) - d(ip)
				if(abs(h)+g.eq.abs(h)) then
					t = a(ip,iq)/h
				else
					theta = 0.5*h/a(ip,iq)
					t = 1.0/&
					(abs(theta)+sqrt(1.0+theta**2))
					if(theta.lt.0.0) t=-t
				endif
					
				c     = 1.0/sqrt(1.0+t**2)
				s     = t*c
				tau   = s/(1.0+c)
				h     = t*a(ip,iq)
				z(ip) = z(ip) - h
				z(iq) = z(iq) + h
				d(ip) = d(ip) - h
				d(iq) = d(iq) + h
				a(ip,iq) = 0.0
!
				do j = 1,ip-1
					g       = a(j,ip)
					h       = a(j,iq)
					a(j,ip) = g - s*(h+g*tau)
					a(j,iq) = h + s*(g-h*tau)
				enddo
				do j = ip+1,iq-1
					g       = a(ip,j)
					h       = a(j,iq)
					a(ip,j) = g - s*(h+g*tau)
					a(j,iq) = h + s*(g-h*tau)
				enddo
				do j = iq+1,n
					g       = a(ip,j)
					h       = a(iq,j)
					a(ip,j) = g - s*(h+g*tau)
					a(iq,j) = h + s*(g-h*tau)
				enddo		
				do j = 1,n
					g       = v(j,ip)
					h       = v(j,iq)
					v(j,ip) = g - s*(h+g*tau)
					v(j,iq) = h + s*(g-h*tau)
				enddo
				nrot = nrot + 1
			endif
		enddo
		enddo
!
		b = b + z
		d = b
		z = 0.0
!	
!		exit si trop d'iterations
!
		if(i.gt.50) then
			write(*,*) "Nb iterations > 50 dans jacobi"
			stop 1
		endif
	enddo
	deallocate(a,b,z)
!
	end subroutine jacobi
!
!-------------------------------------------------------------------------------
!	tri des valeurs et vecteur propres
!	dans l'ordre croissant
!
	subroutine eigsrt(n,d,v)
	implicit none
	integer, intent(in)    :: n
	real,    intent(inout) :: d(:),v(:,:)
!
	integer                :: i,j,k
	real                   :: p
!
	do i = 1,n-1
		k = i
		p = d(i)
		do j = i+1,n
			if(d(j).le.p) then
				k = j
				p = d(j)
			endif
		enddo
		if(k.ne.i) then
			d(k) = d(i)
			d(i) = p
			do j = 1,n
				p      = v(j,i)
				v(j,i) = v(j,k)
				v(j,k) = p
			enddo
		endif
	enddo
	end subroutine eigsrt
!
!-------------------------------------------------------------------------------
!	reduction de Householder d'une matrice symetrique a
!	en une matrice bande
!	q matrice orthogonale effectuant la transformation
!	d elements diagonaux 
!	e elements non diagonaux avec e(1) = 0
!
	subroutine tred2(n,a_in,d,e,q)
	implicit none
	integer, intent(in)  :: n
	real,    intent(in)  :: a_in(:,:)
	real,    intent(out) :: d(:),e(:),q(:,:)
!
	integer              :: i,j,k,l
	real                 :: f,g,h,hh,scale
	real, pointer        :: a(:,:)
!
	allocate(a(n,n))
	a = a_in
!
	do i = n,2,-1
		l = i-1
		h = 0.0
		if(l.gt.1) then
			scale = sum(abs(a(i,1:l)))
			if(scale.eq.0.0) then
				e(i) = a(i,l)
			else
				a(i,1:l) = a(i,1:l)/scale
				h        = sum(a(i,1:l)**2)
				f        = a(i,l)
				g        = -sign(sqrt(h),f)
				e(i)     = scale*g
				h        = h - f*g
				a(i,l)   = f - g
				a(1:l,i) = a(i,1:l)/h
				do j = 1,l
					e(j) = (sum(a(j    ,1:j)*a(i,  1:j))+& 
					        sum(a(j+1:l,  j)*a(i,j+1:l)) &
					       )/h
				enddo
				f      = sum(e(1:l)*a(i,1:l))
				hh     = f/(h+h)
				e(1:l) = e(1:l) - hh*a(i,1:l) 
				do j = 1,l
					a(j,1:j) = a(j,1:j) - a(i,j)*e(  1:j) -&
					                      e(  j)*a(i,1:j)
				enddo
			endif
		else
			e(i) = a(i,l)
		endif
		d(i) = h
	enddo
!
	d(1) = 0.0
	e(1) = 0.0
	do i = 1,n
		l = i-1
		if (d(i).ne.0.0) then
			do j = 1,l
				g = sum(a(i,1:l)*a(1:l,j))
				do k = 1,l
					a(k,j) = a(k,j) - g*a(k,i)
				enddo
			enddo
		endif
		d(i)   = a(i,i)
		a(i,i) = 1.0
		do j = 1,l
			a(i,j) = 0.0
			a(j,i) = 0.0
		enddo
	enddo
!
	q = a
	deallocate(a)
!
	end subroutine tred2
!
!-------------------------------------------------------------------------------
!
	subroutine tqli(n,d_in,e_in,vp,z)
	integer, intent(in)    :: n
	real,    intent(in)    :: d_in(:),e_in(:)
	real,    intent(out)   :: vp(:)
	real,    intent(inout) :: z(:,:)
	integer                :: i,iter,maxiter,k,l,m
	real                   :: b,c,dd,f,g,p,r,s
	real, pointer          :: d(:),e(:),ff(:)
!
	maxiter = 200
!
	allocate(d(n),e(n),ff(n))
	d = d_in
	e = e_in
!
!	renumerotation des elements de e
	do i = 2,n
		e(i-1) = e(i)
	enddo
	e(n) = 0.0
!
	do l = 1,n
		iter = 0
		iterate: do
			do m = l,n-1
				dd = abs(d(m)) + abs(d(m+1))
				if (abs(e(m))+dd.eq.dd) exit
			enddo
			if(m.eq.l) exit iterate
			if(iter.eq.maxiter) then
				write(*,*) "Nb iterations > maxiter tqli"
				stop 1
			endif
			iter = iter + 1
			g    = (d(l+1)-d(l))/(2.0*e(l))
			r    = sqrt(g*g+1.0)
			g    = d(m) - d(l) + e(l)/(g+sign(r,g))
			s    = 1.0
			c    = 1.0
			p    = 0.0
			do i= m-1,l,-1
				f      = s*e(i)
				b      = c*e(i)
				r      = sqrt(f**2+g**2)
				e(i+1) = r
				if(r.eq.0.0) then
					d(i+1) = d(i+1) - p
					e(m)   = 0.0
					cycle iterate
				endif
				s        = f/r
				c        = g/r
				g        = d(i+1) - p
				r        = (d(i)-g)*s+2.0*c*b
				p        = s*r
				d(i+1)   = g+p
				g        = c*r-b
				ff       = z(:,i+1)
				z(:,i+1) = s*z(:,i) + c*ff
				z(:,i  ) = c*z(:,i) - s*ff
			enddo
!
			d(l) = d(l) - p
			e(l) = g
			e(m) = 0.0
		enddo iterate
	enddo
	vp = d
	deallocate(d,e,ff)
	end subroutine tqli
!
!-------------------------------------------------------------------------------
!	diagonalisation d'une matrice reelle symetrique a selon 
!	la procedure tqli
!
	subroutine diago_sym(n,a,val,vec)
	implicit none
	integer, intent(in)  :: n
	real,    intent(in)  :: a(:,:)
	real,    intent(out) :: val(:),vec(:,:)
	real, pointer        :: q(:,:),d(:),e(:)
	integer              :: i,j
!	
!	verification du caractere symetrique de a
	do i = 1,n
	do j = i+1,n
	if(a(i,j).ne.a(j,i)) then
		write(*,"(2i4,a12)") i,j,"Non symetrique" 
		write(*,"(3e12.5)") a(i,j),a(j,i),a(i,j)-a(j,i)
	endif
	enddo
	enddo
!
!	reduction sous forme tridiagonale
!	recherche des valeurs propres
	allocate(q(n,n),d(n),e(n))
	call tred2 (n,a,d,e,q)
	call tqli  (n,d,e,val,q)
	vec = q
	deallocate(q,d,e)
	call eigsrt(n,val,vec)
!
	end subroutine diago_sym
!
!-------------------------------------------------------------------------------
!	diagonalisation d'une matrice hermitienne a selon la procedure tqli
!
	subroutine diagocplx_sym(n,are,aim,val,vecre,vecim)
	implicit none
	integer, intent(in)  :: n
	real,    intent(in)  :: are(:,:),aim(:,:)
	real,    intent(out) :: val(:),vecre(:,:),vecim(:,:)
	real, pointer        :: a(:,:),q(:,:),d(:),e(:),val2(:)
	integer              :: n2,i,j
!
!	reduction sous forme tridiagonale
!	recherche des valeurs propres
	n2 = 2*n
	allocate(a(n2,n2),d(n2),e(n2),q(n2,n2),val2(n2))
!
	a(  1:n,   1: n) =  are
	a(n+1:n2,n+1:n2) =  are
 	a(  1:n, n+1:n2) = -aim
 	a(n+1:n2,  1: n) =  aim
!
	call tred2 (n2,a,d,e,q)
	call tqli  (n2,d,e,val2,q)
	call eigsrt(n2,val2,q)
	do i = 1,n
		j            = 2*(i-1)+1
		val(i)       = val2(j)
		vecre(1:n,i) = q(  1: n,j)
		vecim(1:n,i) = q(n+1:n2,j)
	enddo
	deallocate(a,d,e,q,val2)
!
	end subroutine diagocplx_sym
!
!-------------------------------------------------------------------------------
!	diagonalisation d'une matrice reelle a selon la methode de jacobi
!
	subroutine diago2(n,a,val,vec)
	implicit none
	integer, intent(in)  :: n
	real,    intent(in)  :: a(:,:)
	real,    intent(out) :: val(:),vec(:,:)
	integer              :: nrot
!
	call jacobi(n,a,val,vec,nrot)
	write(*,*) "NROT = ",nrot
	call eigsrt(n,val,vec)
!
	end subroutine diago2
!
!-------------------------------------------------------------------------------
!	diagonalisation d'une matrice reelle a selon la methode de Lanczos
!
	subroutine lanczos(n,a,u1,val,vec)
	implicit none
	integer, intent(in)  :: n
	real,    intent(in)  :: a(:,:),u1(:)
	real,    intent(out) :: val(:),vec(:,:)
!
	integer              :: i,j
	real,    pointer     :: alfa(:),beta(:),u(:,:),v(:)
!
	allocate(alfa(n),beta(n),v(n),u(n,n))
!
!	premier vecteur de la base
	i        = 1
	u(:,i)   = u1/sqrt(sum(u1*u1))
	beta(i)  = 0
!
!	calcul des vecteurs de la base par iterations successives
!	v           = A u_i
!	a_i         = u_i A ui
!	b_i+1^2     = u_i A A u_i - a_i^2 -b_i^2
!	b_i+1 u_i+1 = A u_i - a_i u_i - b_i u_i-1 
! 
	do
		v         = matmul(a,u(:,i))
		alfa(i)   = sum(u(:,i)*v)
		if(i.lt.n) then
			beta(i+1) = sum(v*v) - alfa(i)**2 - beta(i)**2
			if(beta(i+1).lt.0.0) then
				write(*,*) "i=",i,"beta_i+1 < 0"
				stop 1
			endif
			if(beta(i+1).eq.0.0) then
				exit
			else
				beta(i+1) = sqrt(beta(i+1))
				if(i.eq.1) then
				u(:,i+1)  = v-alfa(i)*u(:,i)
				else
				u(:,i+1)  = v-alfa(i)*u(:,i)-beta(i)*u(:,i-1)
				endif
				u(:,i+1)  = u(:,i+1)/beta(i+1)
				i         = i + 1
			endif
		else
			exit
		endif
	enddo
	if(i.eq.n) then
		call tqli  (n,alfa,beta,val,u)
		vec = u
	endif
	deallocate(alfa,beta,v,u)
!	
	end subroutine lanczos
!-------------------------------------------------------------------------------
!
	end module m_diag_sym
	
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	recherche des vecteurs propres connaissant les valeurs propres
!	par resolution d'un systeme lineaire
!	la methode n'est pas au point pour les racines multiples ...
!	avec normalisation pour des bosons
!
	module m_vecboson
	use m_diag
	use m_diag_sym
	use m_solvex
	implicit none
	contains
!
!	-----------------------------------------------------------------------
!	diagonalisation et calcul des vecteurs propres : Methode originale
!	
	subroutine diagboson0(n2,ar,ai,dr,ur,ui,vr,vi)
	implicit none
	integer, intent(in)  :: n2
	real,    intent(in)  :: ar(:,:),ai(:,:)
	real,    intent(out) :: dr(:),ur(:,:),ui(:,:),vr(:,:),vi(:,:)
!
	integer              :: i,j,k,n
	real,     pointer    :: di(:)
	real,     pointer    :: a2r(:,:),a2i(:,:),b2r(:),b2i(:)
	real,     pointer    :: xpr(:),xpi(:)
	real                 :: xr
	logical              :: success
!	
!	dimension
	n = n2/2
!
!	diagonalisation
	allocate(di(n2))
	call diagocplx(n2,ar,ai,dr,di,success)
	if(success.eqv..false.) then
		write(*,*) "Buuuuuuurp"
		stop
	endif
	do k = 1,n
		dr(k  ) =  dr(k+n)
		dr(k+n) = -dr(k  )
	enddo
	deallocate(di)
!
!	vecteurs propres
	allocate(a2r(n2-1,n2-1),a2i(n2-1,n2-1),b2r(n2-1),b2i(n2-1))
	allocate(xpr(n2),xpi(n2))
	do k = 1,n
		if(dr(k).eq.dr(k+1)) write(*,*) "Matrice de rang 2"
		b2r = -ar(2:n2,1)
		b2i = -ai(2:n2,1)
		a2r =  ar(2:n2,2:n2)
		a2i =  ai(2:n2,2:n2)
		xpr    = 0.0
		xpi    = 0.0
		xpr(1) = 1.0
		xpi(1) = 0.0

		do i=1,n2-1
			a2r(i,i) = a2r(i,i) - dr(k)
		enddo
		call solvex(n2-1,a2r,a2i,b2r,b2i,xpr(2:n2),xpi(2:n2))
		ur(:,k) = xpr(  1:n )
		ui(:,k) = xpi(  1:n )
		vr(:,k) = xpr(n+1:n2)
		vi(:,k) = xpi(n+1:n2)
	enddo
	deallocate(a2r,a2i,b2r,b2i)
!
!	on impose la relation de commutation canonique.
!
	do k = 1,n
		xpr(  1:n ) = ur(:,k)
		xpr(n+1:n2) = vr(:,k)
		xpi(  1:n ) = ui(:,k)
		xpi(n+1:n2) = vi(:,k)
!
		xr =    sum(xpr(  1:n )*xpr(  1:n )+xpi(  1:n )*xpi(  1:n ))
		xr = xr-sum(xpr(n+1:n2)*xpr(n+1:n2)+xpi(n+1:n2)*xpi(n+1:n2))
		if(xr.eq.0.0) then
			write(*,*) "xr = 0.0"
			cycle
		endif
		if(xr.lt.0.0) then
			xr      =  sqrt(abs(xr))
			ur(:,k) = -xpi(  1:n )/xr
			vr(:,k) = -xpi(n+1:n2)/xr
			ui(:,k) =  xpr(  1:n )/xr
			vi(:,k) =  xpr(n+1:n2)/xr
		elseif(xr.gt.0.0) then
			xr      =  sqrt(abs(xr))
			ur(:,k) =  xpr(  1:n )/xr
			vr(:,k) =  xpr(n+1:n2)/xr
			ui(:,k) =  xpi(  1:n )/xr
			vi(:,k) =  xpi(n+1:n2)/xr
		endif
	enddo
	deallocate(xpr,xpi)
	end subroutine diagboson0
!
!	-----------------------------------------------------------------------
!	calcul des vecteurs et valeurs propres
!	
	subroutine diagboson(n2,ar,ai,dr,pur,pui,pvr,pvi)
	implicit none
	integer, intent(in)  :: n2
	real,    intent(in)  :: ar(:,:),ai(:,:)
	real,    intent(out) :: dr(:),pur(:,:),pui(:,:),pvr(:,:),pvi(:,:)
!
	integer              :: i,j,k,n
	real,     pointer    :: t(:),vr(:,:),vi(:,:)
	real,     pointer    :: rr(:,:),ri(:,:),wr(:,:),wi(:,:)
	real,     pointer    :: g(:),u(:)
	real                 :: x,y,z,p
!
!	initialisation et allocation dyamique de la memoire
	allocate(t(n2),vr(n2,n2),vi(n2,n2))
	allocate(rr(n2,n2),ri(n2,n2),wr(n2,n2),wi(n2,n2))
	allocate(g(n2),u(n2))
	n         = n2/2
	g(  1:n ) =  1.0
	g(n+1:n2) = -1.0
!
!	on diagonalise directement la matrice h
!	elle doit etre definie positive (valeurs propres positives)
!	sinon on genere une erreur
!	write(*,*) "Matrice AR"
!	do i=1,n2
!		write(*,"(10f12.5)") (ar(i,j),j=1,n2)
!	enddo
!	write(*,*) "Matrice AI"
!	do i=1,n2
!		write(*,"(10f12.5)") (ai(i,j),j=1,n2)
!	enddo
!
	call diagocplx_sym(n2,ar,ai,t,vr,vi)
	do i=1,n2
		if(t(i).lt.0.0) then
!			si valeur propore negative, la diagonlaisation echoue
!			on sort de la subroutine en forcant les vecteurs propres à 0.
			write(*,"(a8,i3)") "VP < 0 ",i
			pur = 0.0
			pui = 0.0
			pvr = 0.0
			pvi = 0.0
			return
!			stop
		endif
!		write(*,"(i4,f12.5)") i,t(i)
	enddo
!	write(*,*) "Matrice VR"
!	do i=1,n2
!		write(*,"(10f12.5)") (vr(i,j),j=1,n2)
!	enddo
!	write(*,*) "Matrice VI"
!	do i=1,n2
!		write(*,"(10f12.5)") (vi(i,j),j=1,n2)
!	enddo
!	
!	calcul de la matrice R = T^1/2 V+ g V T^1/2
!	et diagonalisation
	do i=1,n2
	do j=1,n2
		z = sqrt(t(i)*t(j))
		rr(i,j) = z*(sum(vr(:,i)*g*vr(:,j))+sum(vi(:,i)*g*vi(:,j)))
		ri(i,j) = z*(sum(vr(:,i)*g*vi(:,j))-sum(vi(:,i)*g*vr(:,j)))
	enddo
	enddo
	call diagocplx_sym(n2,rr,ri,dr,wr,wi)
!	write(*,*) "Matrice RR"
!	do i=1,n2
!		write(*,"(10f12.5)") (rr(i,j),j=1,n2)
!	enddo
!	write(*,*) "Matrice RI"
!	do i=1,n2
!		write(*,"(10f12.5)") (ri(i,j),j=1,n2)
!	enddo
!	write(*,*) "Valeurs propres de R"
!	do i=1,n2
!		write(*,"(i4,f12.5)") i,dr(i)
!	enddo
!	write(*,*) "Matrice WR"
!	do i=1,n2
!		write(*,"(10f12.5)") (wr(i,j),j=1,n2)
!	enddo
!	write(*,*) "Matrice WI"
!	do i=1,n2
!		write(*,"(10f12.5)") (wi(i,j),j=1,n2)
!	enddo
!
!	tri des valeurs et vecteur propres
!	dans l'ordre decroissant
	do i = 1,n2-1
		k = i
		p = dr(i)
		do j = i+1,n2
			if(dr(j).ge.p) then
				k = j
				p = dr(j)
			endif
		enddo
		if(k.ne.i) then
			dr(k) = dr(i)
			dr(i) = p
			do j = 1,n2
				p       = wr(j,i)
				wr(j,i) = wr(j,k)
				wr(j,k) = p
				p       = wi(j,i)
				wi(j,i) = wi(j,k)
				wi(j,k) = p
			enddo
		endif
	enddo
!
!	calcul de la matrice de changement de base 
!	des ondes de spin
!	P = V T^(-1/2) W A avec A = sqrt(abs(dr))*signe(dr)
!	la moitie superieure contient pur et pui (1      < i < n2/2)
!	la seconde moitie    contient pvr et pvi (n2/2+1 < i < n)
!	on ne calcule que pour j < n
	u = 1.0/sqrt(t)
	do j=1,n
		z = sqrt(abs(dr(j)))
		do i=1,n
			x = sum(vr(i,:)*u*wr(:,j)) - sum(vi(i,:)*u*wi(:,j))
			y = sum(vr(i,:)*u*wi(:,j)) + sum(vi(i,:)*u*wr(:,j))
			pur(i,j) = x*z
			pui(i,j) = y*z
			x = sum(vr(i+n,:)*u*wr(:,j)) - sum(vi(i+n,:)*u*wi(:,j))
			y = sum(vr(i+n,:)*u*wi(:,j)) + sum(vi(i+n,:)*u*wr(:,j))
			pvr(i,j) = x*z
			pvi(i,j) = y*z
		enddo
	enddo
!
!	deallocation de la memoire
	deallocate(t,vr,vi,rr,ri,wr,wi,g,u)
	end subroutine diagboson
!
	end module m_vecboson
	
	!--------------------------------------------------------------------
!
	module m_fibo
	implicit none
	contains
!	-------------------------------------
!	construction de la suite de Fibonacci
!
	subroutine suite_fibo(n,jf,jf1)
	implicit none
	integer, intent(in)  :: n
	integer, intent(out) :: jf,jf1
	integer              :: j
	integer, pointer     :: a(:)
!
	allocate(a(n))	
	a(1) = 1
	a(2) = 1
	do j = 3,n
		a(j) = a(j-1) + a(j-2)
	enddo
	jf  = a(n)
	jf1 = a(n-1)
	deallocate(a)
	end subroutine suite_fibo
!
!	--------------------------------------
!	coordonnees des points sur un cylindre
!	
	subroutine point_fibo(jf,jf1,phi,z)
	implicit none
	integer, intent(in)  :: jf,jf1
	real,    intent(out) :: phi(:),z(:)
	integer              :: j
	real                 :: pi	
	pi = 4.0*atan(1.0)

	do j = 1,jf+1
		z  (j) = -1.0 + 2.0/jf*(j-1.0)
		phi(j) = 2.0*pi*jf1/jf*(j-1.0)
	enddo
	end subroutine point_fibo
	end module m_fibo


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	routines de calcul pour la description du couplage magnon-phonon
!
	module m_couplage
	use m_util2p1
	implicit none
	contains
!
!	----------------------------------------------------
!	retourne le couplage d'échange entre les ions i et j
!
	subroutine cvsavary(g,jpmpm,jpm,jzpm,jzz,ja,jb,jc,j4)
	implicit none
	real,             intent(in)  :: g(3,3),jpmpm,jpm,jzpm,jzz
	real,             intent(out) :: ja,jb,jc,j4
!
	real :: r2,gp,gz
!
	gz = g(3,3)
	gp = (g(1,1)+g(2,2))/2.0
	r2 = sqrt(2.0)
	ja = 4.*( jpmpm-jpm)/(3.*gp*gp) + 4.*r2*jzpm/(3.*gp*gz)+    jzz/(3.*gz*gz)
	jb = 2.*( jpmpm+jpm)/(gp*gp)
	jc = 2.*(-jpmpm+jpm)/(3.*gp*gp) + 4.*r2*jzpm/(3.*gp*gz)- 2.*jzz/(3.*gz*gz)
	j4 = 2.*( jpmpm-jpm)/(3.*gp*gp) -    r2*jzpm/(3.*gp*gz)-    jzz/(3.*gz*gz)
!
	write(*,"(a8,f12.5)") "JPMPM = ",jpmpm
	write(*,"(a8,f12.5)") "JPM   = ",jpm
	write(*,"(a8,f12.5)") "JZPM  = ",jzpm
	write(*,"(a8,f12.5)") "JZZ   = ",jzz
	write(*,"(a8,f12.5)") "JA    = ",ja
	write(*,"(a8,f12.5)") "JB    = ",jb
	write(*,"(a8,f12.5)") "JC    = ",jc
	write(*,"(a8,f12.5)") "J4    = ",j4
	write(*,"(a8,f12.5)") "gp    = ",gp
	write(*,"(a8,f12.5)") "gz    = ",gz
!
!
	end subroutine cvsavary
!
!	----------------------------------------------------
!	retourne le couplage d'échange entre les ions i et j
!
	subroutine mag(i,j,vec,pos,u1,u2,e,nz,modele,dani,dist,jech,qc,affich,jm)
	implicit none
	integer,          intent(in)  :: i,j,nz
	real,             intent(in)  :: vec(:),pos(:,:),u1(3),u2(3),e(3,3)
	real,             intent(in)  :: dani(:,:,:),dist(:),jech(:,:,:),qc(:)
	character(len=8), intent(in)  :: modele
	logical,          intent(in)  :: affich
	real,             intent(out) :: jm(:,:)
!
	integer              :: d,k,l,l1,l2
	real                 :: pi,p1,p2,p3,p4,z,dd
	real                 :: r0(3),rc(3),ss(3),matj(3,3),u(3),v(3),w(3)
!
	d  = 3
	pi = 4.0*atan(1.0)
!
	rc = pos(j,:) + vec - pos(i,:)
	rc = matmul(e,rc)
	dd = sqrt(sum(rc*rc))
	jm = 0.0
	if((i.eq.j).and.(dd.eq.0.0)) then
		do l1=1,d
		do l2=1,d
!			jm(l1,l2) = dval*dech(l1)*dech(l2)
			jm(l1,l2) = dani(i,l1,l2)
		enddo
		enddo
	endif
	if(dd.gt.0.0) then
			do k=nz,1,-1
			if(dd.le.dist(k)) jm = jech(k,:,:)
			enddo
			if((modele.eq."DIP").and.(dd.le.dist(1))) then
				z = jm(1,1)
				r0 = pos(j,:) + vec - pos(i,:)
				r0 = matmul(e,r0)/dd
				do l1=1,d
				do l2=1,d
					jm(l1,l2) = -r0(l1)*r0(l2)*z
				enddo
				enddo
			endif
			if((modele.eq."BV").and.(dd.le.dist(1))) then
				do l1=1,d
					ss(l1) = jm(l1,l1)
				enddo
				r0 = pos(j,:) + vec - pos(i,:)
				r0 = matmul(e,r0)/dd
				p1 = acos(r0(3))
				if((p1.eq.0.0).or.(p1.eq.pi)) then
					p1 = 0.0
				else
					if(1.-abs(r0(1)/sin(p1)).lt.(1e-4)) then
						if(abs(r0(1)/sin(p1)-1.).lt.(1e-4)) p2=0.0
						if(abs(r0(1)/sin(p1)+1.).lt.(1e-4)) p2=pi
					else
						if(r0(2).lt.0.0) p2 = (2*pi-acos(r0(1)/sin(p1)))
						if(r0(2).ge.0.0) p2 = (     acos(r0(1)/sin(p1)))
					endif
				endif
				matj(3,:) = r0
				matj(2,1) = +sin(p2)
				matj(2,2) = -cos(p2)
				matj(2,3) = 0.0
				call vct(matj(2,:),matj(3,:),matj(1,:))
				do l1=1,d
				do l2=1,d
					jm(l1,l2) = sum(matj(:,l1)*ss(:)*matj(:,l2))
				enddo
				enddo
			endif
			if((modele.eq."ANISOBV").and.(dd.le.dist(1))) then
				p1 = jech(1,1,1)
				p2 = jech(2,1,1)
				p3 = jech(3,1,1)
				p4 = jech(4,1,1)
!
!				convention de Pierre Bonville
!
!				j1-> lambda_a,j2-> lambda_b,j3-> lambda_c
!
!				axe c 
				u = u1-u2
				u = u/sqrt(sum(u*u))
!				axe a 
				v = u1+u2
				v = v/sqrt(sum(v*v))
!				axe b 
				call vct(u,v,w)
				w = w/sqrt(sum(w*w))
!
!				u est le long de la liaison, v perp à la liaison
!				mais dans le plan des axes de cef et w ferme le triedre
!
!				terme symetrique
				jm = 0.0
				do l1=1,d
				do l2=1,d
					jm(l1,l2) = p1*v(l1)*v(l2)+p2*w(l1)*w(l2)+p3*u(l1)*u(l2)
				enddo
				enddo
!
!				terme antisymetrique
				matj      = 0.0
				matj(1,2) = +w(3)
				matj(2,1) = -w(3)
				matj(1,3) = -w(2)
				matj(3,1) = +w(2)
				matj(2,3) = +w(1)
				matj(3,2) = -w(1)
!				jm = jm + p4*matj*sqrt(2.0)
				jm = jm - p4*matj*sqrt(2.0)
			endif
		endif
		if(affich) then
		if(sum(jm*jm).ne.0.0) then
			l2=1
			write(*,"(2i4,3f7.3,4f12.5,x,a8)") i,j,(vec(k),k=1,d),(jm(l2,l1),l1=1,d),dd,modele
			do l2=2,d
			write(*,"(8x,21x,3f12.5)") (jm(l2,l1),l1=1,d)
			enddo
		endif
		endif
	end subroutine mag
!
!	----------------------------------------------------
!	retourne le couplage elastique entre les ions i et j
!
	subroutine for(i,j,vec,pos,e,nz,modele,dist,ctef,qc,fm)
	implicit none
	integer,          intent(in)  :: i,j,nz
	real,             intent(in)  :: vec(:),pos(:,:),e(3,3)
	real,             intent(in)  :: dist(:),ctef(:,:,:),qc(:)
	character(len=3), intent(in)  :: modele
	real,             intent(out) :: fm(:,:)
!
	integer              :: d,ntt,k,l,l1,l2
	real                 :: pi,p1,p2,z,dd
	real                 :: r0(3),rc(3),jm(3,3)
!
	d  = 3
	pi = 4.0*atan(1.0)
!
	rc = pos(j,:) + vec - pos(i,:)
	rc = matmul(e,rc)
	dd = sqrt(sum(rc*rc))
	jm = 0.0
	if(dd.gt.0.0) then
		do k=nz,1,-1
			if(dd.le.dist(k)) jm = ctef(k,:,:)
		enddo
	endif
	fm = jm
	end subroutine for
!
!	----------------------------------------------------
!	retourne le couplage magnon-phonon
!	de la forme
!	(uur,uui) = sum_R (u_j-u_i+R)/|u_j-u_i+R| exp(iq R) et |u_j-u_i+R| < dist
!	uu = uur(q=0)
!	uo = \sum_R (u_j-u_i+R)/|u_j-u_i+R|.(u_j-u_i+R)
!
	subroutine cpl(i,j,vec,pos,e,nz,modele,dist,ctef,qc,uo,uu)
	implicit none
	integer,          intent(in)  :: i,j,nz
	real,             intent(in)  :: vec(:),pos(:,:),e(3,3)
	real,             intent(in)  :: dist(:),ctef(:,:,:),qc(:)
	character(len=3), intent(in)  :: modele
	real,             intent(out) :: uo,uu(3)
!
	integer              :: d,k,l,l1,l2,l3
	real                 :: pi,p1,p2,z,dd
	real                 :: r0(3),rc(3),ss(3)
!
	d  = 3
	pi = 4.0*atan(1.0)
!
	uo  = 0.0
	uu  = 0.0
	rc = pos(j,:) + vec - pos(i,:)
	rc = matmul(e,rc)
	dd = sqrt(sum(rc*rc))
	ss = 0.0
	if(dd.gt.0.0) then
		do k=nz,1,-1
			if(dd.le.dist(k)) then
				ss = (rc/dd)
			endif
		enddo
	endif
	uo  = sum(ss*rc)
	uu  = ss
!	if(i.ne.j) uu(1) = 1
	end subroutine cpl
!
	end module m_couplage
	
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	Calcul des operateurs de stevens a partir des 
!	composantes +, - et z d'un moment cinetique j
!
	module m_stevens
	implicit none
	contains
!
!	-----------------------------------------------------------------------
!
	subroutine stevens1(j,vj,o20,o22,o40,o42,o43,o44,o60,o62,o63,o64,o66)
	implicit none
	real,    intent(in)  :: j
	real,    intent(in)  :: vj(:,:,:)
	real,    intent(out) :: o20(:,:),o22(:,:)
	real,    intent(out) :: o40(:,:),o42(:,:),o43(:,:),o44(:,:)
	real,    intent(out) :: o60(:,:),o62(:,:),o63(:,:),o64(:,:),o66(:,:)
!	
	real,    pointer     :: id(:,:),a(:,:)
	real,    pointer     :: vp (:,:),vp2(:,:),vp3(:,:),vp4(:,:),vp6(:,:)
	real,    pointer     :: vm (:,:),vm2(:,:),vm3(:,:),vm4(:,:),vm6(:,:)
	real,    pointer     :: vz (:,:),vz2(:,:),vz3(:,:),vz4(:,:),vz6(:,:)
	real                 :: j1
	integer              :: nu,i
!
	nu = int(2*j+1)
	j1 = j+1.0
!
	allocate(id(nu,nu),a(nu,nu),vp(nu,nu),vm(nu,nu),vz(nu,nu))
	allocate(vp2(nu,nu),vm2(nu,nu),vz2(nu,nu))
	allocate(vp3(nu,nu),vm3(nu,nu),vz3(nu,nu))
	allocate(vp4(nu,nu),vm4(nu,nu),vz4(nu,nu))
	allocate(vp6(nu,nu),vm6(nu,nu),vz6(nu,nu))
!
	id = 0.0
	do i = 1,nu
		id(i,i) = 1.0
	enddo
!
	vp  = vj(1,:,:)
	vp2 = matmul(vp ,vp )
	vp3 = matmul(vp ,vp2)
	vp4 = matmul(vp ,vp3)
	vp6 = matmul(vp3,vp3)
!
	vm  = vj(2,:,:)
	vm2 = matmul(vm ,vm )
	vm3 = matmul(vm ,vm2)
	vm4 = matmul(vm ,vm3)
	vm6 = matmul(vm3,vm3)
!
	vz  = vj(3,:,:)
	vz2 = matmul(vz ,vz )
	vz3 = matmul(vz ,vz2)
	vz4 = matmul(vz ,vz3)
	vz6 = matmul(vz3,vz3)
!	---------------------
	o20 = 3*vz2 -j*j1*id
	o22 = (vp2+vm2)/2.0
!	---------------------
	o40 = 35.0*vz4 + (-30.0*j*j1+25.0)*vz2
	o40 = o40 + (-6.0*j*j1 + 3.0*j*j*j1*j1)*id

	a   = 7.0*vz2 -(j*j1+5.0)*id
	o42 = (matmul(a,vp2+vm2) + matmul(vp2+vm2,a))/4.0
	o43 = (matmul(vz,vp3+vm3) + matmul(vp3+vm3,vz))/4.0
	o44 = (vp4+vm4)/2.0
!	---------------------
	o60 = 231.0*vz6 + (-315.0*j*j1 + 735)*vz4
	o60 = o60 + (105.0*j*j*j1*j1 - 525.0*j*j1 + 294)*vz2
	o60 = o60 + (-5.0*j*j*j*j1*j1*j1+ 40.0*j*j*j1*j1 - 60.0*j*j1)*id
	a   = 33.0*vz4 - (18.0*j*j1+123.0)*vz2
	a   = a + (j*j*j1*j1 + 10.0*j*j1 + 102.0)*id
	o62 = (matmul(a,vp2+vm2) + matmul(vp2+vm2,a))/4.0
	a   = 11.0*vz3 - (3.0*j*j1 + 59.0)*vz
	o63 = (matmul(a,vp3+vm3) + matmul(vp3+vm3,a))/4.0
	a   = 11.0*vz2-(j*j1+38.0)*id
	o64 = (matmul(a,vp4+vm4)+matmul(vp4+vm4,a))/4.0
	o66 = (vp6+vm6)/2.0
!
	deallocate(id,a,vp,vm,vz,vp2,vm2,vz3,vp3,vm3,vz2,vp4,vm4,vz4,vz6)
	end subroutine stevens1
!
!
!	-----------------------------------------------------------------------
!
	subroutine stevens2(j,vj,pxy,pyz,pzx)
	implicit none
	real,    intent(in)  :: j
	real,    intent(in)  :: vj(:,:,:)
	real,    intent(out) :: pxy(:,:),pyz(:,:),pzx(:,:)
!	
	real,    pointer     :: vp (:,:),vm (:,:),vz (:,:)
	integer              :: nu
!
	nu = int(2*j+1)
	allocate(vp(nu,nu),vm(nu,nu),vz(nu,nu))
!
	vp  = vj(1,:,:)
	vm  = vj(2,:,:)
	vz  = vj(3,:,:)
!
!	attention, pxy et pyz sont imaginaires purs
	pxy = (matmul(vp,vp)   -matmul(vm,vm)   )/4.0
	pyz = (matmul(vp-vm,vz)+matmul(vz,vp-vm))/4.0
	pzx = (matmul(vp+vm,vz)+matmul(vz,vp+vm))/4.0
!
	deallocate(vp,vm,vz)
	end subroutine stevens2
!
!	-----------------------------------------------------------------------
!
	end module m_stevens
	
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!	Hamiltonien de champ cristallin
!
	module m_cef
	use m_stevens
	implicit none
	contains	
!	
	subroutine cef(wj,b20,b22,b2m2,b21,b2m1,b40,b42,b43,b44,b60,b62,b63,b64,b66,dd,a,h,hr,hi)
	implicit none
	real,    intent(in)  :: wj,b20,b22,b2m2,b21,b2m1,b40,b42,b43,b44,b60,b62,b63,b64,b66,dd
	real,    intent(in)  :: h(3),a(6)
	real,    intent(out) :: hr(:,:),hi(:,:)
!
	integer              :: d,i,j,nu
	real                 :: pi,mj,x,y
	real,        pointer :: vj(:,:,:),u(:,:,:)
	real,        pointer :: jx(:,:),jy(:,:),jz(:,:)
	real,        pointer :: jx2(:,:),jy2(:,:),jz2(:,:)
	real,        pointer :: jxy(:,:),jxz(:,:),jyz(:,:)
	real,        pointer :: o20(:,:),o22(:,:)
	real,        pointer :: o40(:,:),o42(:,:),o43(:,:),o44(:,:)
	real,        pointer :: o60(:,:),o62(:,:),o63(:,:),o64(:,:),o66(:,:)
!
!	-----------------------------------------------------------------------
!	constantes
!	----------
!
	d  = 3
	pi = 4.0*atan(1.0)
	nu = int(2.*wj+1) 
!
!	-----------------------------------------------------------------------!
!	degenerescence 2j+1
!	dimension base espace vectoriel local (2j+1)
!
!	le hamiltonien a diagonaliser s'ecrit dans 
!	la base des etats decrits par mj
!	avec -j < mj < j 
!	et reperes par un entier i tel que mj = i-1-j
!
!	le hamiltonien est de la forme 
!	h= Hloc + g hmag.J
!	On introduit les composantes +(1),-(2),z(3) 
!	des operateurs moments cinetiques de sorte que 
!	les pdts scalaires sont alors de la forme 
!	(a^+ b^- + a^- b^+)/2 + a^z b^z 
!
	allocate(vj(d,nu,nu))
	allocate(jx(nu,nu),jy(nu,nu),jz(nu,nu),jx2(nu,nu),jy2(nu,nu),jz2(nu,nu))
	allocate(jxy(nu,nu),jxz(nu,nu),jyz(nu,nu))
	allocate(o20(nu,nu),o22(nu,nu))
	allocate(o40(nu,nu),o42(nu,nu),o43(nu,nu),o44(nu,nu))
	allocate(o60(nu,nu),o62(nu,nu),o63(nu,nu),o64(nu,nu),o66(nu,nu))
!
!
!	J de l'atome
!	calcul des composantes +,-,z
!	de l'operateur J 
!
	vj = 0.d0
	do i = 1,nu
		mj = (wj+1)-i
		x = sqrt(wj*(wj+1.d0)-(i-wj-1.d0)*(i-wj))
		if(i.lt.nu) then
			vj(1,i,i+1) = x
			vj(2,i+1,i) = x
		endif
		vj(3,i,i) = mj
	enddo
!
!	calcul des operateurs de stevens
	call stevens1(wj,vj,o20,o22,o40,o42,o43,o44,o60,o62,o63,o64,o66)
!
!	composantes de J
!	attention jy = (J+ - J-)/2i = -i (J+ - J-)/2
	jx =  (vj(1,:,:)+vj(2,:,:))/2.0
	jy = -(vj(1,:,:)-vj(2,:,:))/2.0
	jz =   vj(3,:,:)
!
!	calcul de jxz = (Jx Jz + Jz Jx)
	jxz = matmul(jx,jz)+matmul(jz,jx)
	jyz = matmul(jy,jz)+matmul(jz,jy)
	jxy = matmul(jx,jy)+matmul(jy,jx)
!
!	calcul de jx2 = Jx Jx
	jx2 = matmul(jx,jx)
	jy2 = matmul(jy,jy)
	jz2 = matmul(jz,jz)
!
!	termes de champ cristallin
	hr =      b20*o20 + b22*o22 + b21*jxz
	hr = hr + b40*o40 + b42*o42 + b43*o43 + b44*o44
	hr = hr + b60*o60 + b62*o62 + b63*o63 + b64*o64 + b66*o66
	hi = b2m2*jxy + b2m1*jyz
!
!	attention couplage du type +(hx hy hz).J
!	soit = h^x s^x + h^y s^y + s^z h^z = hx/2 (s^+ + s^-) - i hy/2 (s^+ - s^-) + h^z s^z
	hr = hr + h(1)*jx + h(3)*jz
	hi = hi + h(2)*jy
!
!	perturbation (distorsion)
	hr = hr + dd*(a(1)*jx2+a(2)*jy2+a(3)*jz2         +a(5)*jxz         )
	hi = hi + dd*(                          +a(4)*jxy         +a(6)*jyz) 
!
	deallocate(vj,jx,jy,jz,jx2,jy2,jz2,jxy,jxz,jyz)
	deallocate(o20,o22,o40,o42,o43,o44,o60,o62,o63,o64,o66)
!
	end subroutine cef
!
!	-----------------------------------------------------------------------
!	calcul des composantes x,y,z de l'operateur J
!
	subroutine jz(nu,wj,mjzr,mjzi)
	implicit none
	integer, intent(in) :: nu
	real,    intent(in) :: wj
	real,   intent(out) :: mjzr(:,:,:),mjzi(:,:,:)
!
	integer             :: i,j,l1,l2,d
	real                :: x,y,mj
	real, pointer       :: w(:,:,:)
!
!	calcul des composantes +,-,z de l'operateur J
	d = 3
	allocate(w(d,nu,nu))
	w = 0.0
	do i = 1,nu
		mj = (wj+1)-i
		x = sqrt(wj*(wj+1.d0)-(i-wj-1.d0)*(i-wj))
		if(i.lt.nu) then
			w(1,i,i+1) = x
			w(2,i+1,i) = x
		endif
		w(3,i,i) = mj
	enddo
!
!	calcul des composantes x,y,z de l'operateur J
	mjzr = 0.0
	mjzi = 0.0
	mjzr(3,:,:) =   w(3,:,:)
	mjzr(1,:,:) =  (w(1,:,:)+w(2,:,:))/2.0
	mjzi(2,:,:) = -(w(1,:,:)-w(2,:,:))/2.0
	deallocate(w)
	end subroutine jz
!
!	-----------------------------------------------------------------------
!
	end module m_cef
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	module m_fmag
	implicit none
	contains	
	
!	------------------------------------------
!	Facteurs de structure magnetique
!	FM(q)=sum_i S(i) fc(i) exp(-w(i) q^2/4pi) exp(i 2 pi q r_i)
!	q = vecteur de diffusin en unites reduites
!	r = position en unites reduites
!
	subroutine fm(mtot,q,qc,r,s,ff,fmr,fmi)
	implicit none
	integer, intent(in)  :: mtot
	real,    intent(in)  :: q(:),qc(:)
	real,    intent(in)  :: r(:,:),s(:,:),ff
!	real,    intent(in)  :: fc(:),w(:)
	real,    intent(out) :: fmr(:),fmi(:)
!
	integer              :: i,j
	real                 :: pi
	real                 :: s1,x,y
	real                 :: sre(3),sim(3)
!
	pi  = 4.0*atan(1.0)

!	facteurs de structure magnetiques
!
	sre = 0.0
	sim = 0.0
	do i = 1,mtot
		x     = 2.0*pi*sum(q*r(i,:))
!		y     = w(i)*sum(qc*qc)/(4.0*pi)**2
		y     = 0.0
!
		s1 = sum(s(i,:)*s(i,:))
		if(s1.ne.0.0) then
!		sre   = sre + exp(-y)*cos(x)*s(i,:)*fc(i)
!		sim   = sim + exp(-y)*sin(x)*s(i,:)*fc(i)
		sre   = sre + cos(x)*s(i,:)
		sim   = sim + sin(x)*s(i,:)
		endif
	enddo
	fmr = ff*sre
	fmi = ff*sim
!
	end subroutine fm
!
	end module m_fmag
	
		program sw2
	use m_borne
	use m_cartes
	use m_util2p1
	use m_solve
	use m_solvex
	use m_vecboson
	use m_fibo
	use m_couplage
	use m_stevens
	use m_cef
	use m_fmag
	implicit none
!
	type atome
		character(len=mxcar) nam
		real                 pos(3)
	end type atome
!
	integer              :: n,m,n0,n2,d,nw,nt,ntt,nz,niter
	integer              :: np,np1,np2,np3,nq,nq2
	integer              :: ntw
	integer              :: i,j,k,l,ll
	integer              :: k1,k2,k3,k4,l1,l2,l3,l4
	integer              :: p,pr,calc
	integer              :: ii,jj
	integer              :: nfibo,jf,jf1
	integer              :: ncalc,ncalc2
	integer              :: nato
	integer              :: nbx,nby,nbz
	integer              :: nnx,nny,nnz,ntt2
	integer              :: nbk,nk,pq,pk
	integer,pointer      :: vq(:),vq2(:),vk(:),vk2(:),tab(:,:),ntab(:)
!
	real                 :: jm(3,3),jmr(3,3),jmi(3,3),jm0(3,3),jmp(3,3),vcr(3),vci(3),uu(3),uo
	real                 :: vdmr(3),vdmi(3),udm(3),uur(3),uui(3)
	real                 :: ax,ay,az,alfa,beta,gama,e(3,3),ee(3,3),mat(3,3)
	real                 :: q0(3),dq(3),dq1(3),dq2(3),dq3(3),q(3),en0,qmin,qmax,wmax,wmin
	real                 :: jpmpm,jpm,jzpm,jzz
	real                 :: fga,fa,fgb,fb,fgc,fc,fgd,s,fmag
	real                 :: pi,d2r,mub,muo,kb,ec,j2unit,k2unit
	real                 :: dd,sig,vsig,vol,unit1,ki,emax,reg1,reg2,reg3
	real                 :: ex(3),ey(3),ez(3),rc(3),qc(3),qc2(3),qc3(3),ss(3),r0(3),matj(3,3)
	real                 :: kprop(3),etao(3),zdro(3),zdio(3),v(3)
	real                 :: x,y,z,xr,xi,yr,yi,wr,wi,p1,p2,p3,p4,p5,g,gdm,gksea,gdiag
	real                 :: nn1(3,3),rr1(3,3),ri1(3,3),nn2(3,3),rr2(3,3),ri2(3,3)
	real                 :: pdtuvr,pdtuvi
	real,pointer         :: pos(:,:),pos1(:,:),sp(:),phi(:),theta(:),ma(:),listvec(:,:)
	real,pointer         :: listvec2(:,:)
	real,pointer         :: dech(:,:),dval(:),dani(:,:,:),jech(:,:,:,:,:),ctef(:,:,:,:,:)
	real,pointer         :: dm(:,:,:,:),pt(:,:,:,:)
	real,pointer         :: dist(:,:,:),disf(:,:,:)
	real,pointer         :: matl(:,:,:),matr(:,:,:),zdr(:,:),zdi(:,:),eta(:,:)
	real,pointer         :: omega(:),h(:)
	real,pointer         :: ar(:,:),ai(:,:),a2r(:,:),a2i(:,:)
	real,pointer         :: dr(:,:),ur(:,:),ui(:,:),vr(:,:),vi(:,:)
	real,pointer         :: u1(:),u2(:)
	real,pointer         :: tq(:,:),tqc(:,:),tqr(:,:),su(:),listeq(:,:)
	real,pointer         :: tq2(:,:),su2(:),tqr2(:,:)
	real,pointer         :: w(:),matc(:,:,:),fq(:,:,:),fpartiel(:,:,:)
	real,pointer         :: pou(:,:),poud(:),pous(:),poun(:),pouy(:),pouz(:)
	real,pointer         :: vphi(:),vz(:)
	real,pointer         :: vkr(:,:),vki(:,:),fn(:,:),fy(:,:),fz(:,:)
	real,pointer         :: xpr(:),xpi(:),xmr(:),xmi(:)
	real,pointer         :: uprop(:,:),eta0(:,:,:),zdr0(:,:,:),zdi0(:,:,:)
	real,pointer         :: redm(:,:)
	real,pointer         :: rtw(:,:),atw(:)
	real,pointer         :: biq(:,:,:,:),dbiq12(:,:),dbiq13(:,:)
	real,pointer         :: urkn(:,:,:),uikn(:,:,:),vrkn(:,:,:),vikn(:,:,:)
	logical              :: fform,success,affich,trouve,angle
	logical              :: mf,pho,couple,incom
	logical              :: detail,grille
	logical              :: maxime,savary,bbq
	logical              :: simpleq,tof,integ,coupe,coup1d,reduc,twin,cnt,partiel
!
	integer,pointer      :: nu(:)
	integer              :: nui
	real                 :: h0(3),hmol(3),alf,bet,gam,dis(6),temp,gl(3,3),gl2(3,3)
	real,pointer         :: wj(:),gj(:),b20(:),b22(:),b2m2(:),b21(:),b2m1(:),b40(:),b42(:),b43(:),b44(:)
	real,pointer         :: b60(:),b62(:),b63(:),b64(:),b66(:),dis0(:)
	real,pointer         :: xcef(:,:),matlx(:,:,:),matrx(:,:,:)
 	real,pointer         :: hr(:,:),hi(:,:),df(:),vfr(:,:),vfi(:,:),vm(:),nrj(:)
	real,pointer         :: jr(:,:,:),ji(:,:,:),jmy(:,:)
	real,pointer         :: u2r(:,:),u2i(:,:)
!
	character(len=8),pointer :: modele(:,:)
	character(len=mxcar)     :: fres,ffst
	character(len=mxcar)     :: nom1
	character(len=8),pointer :: nom(:)
	character(len=8)         :: unit
!
	type(atome),pointer  :: lato(:)
!
	integer              :: lu,nbcle,nucle,ios
	character(len=mxcar) :: carte
	character(len=mxcar) :: moc(mxcle)
	character(len=mxcar) :: av(mxcle)
	integer              :: iv(mxcle)
	real                 :: rv(mxcle)
!
!	-----------------------------------------------------------------------
!
	call banner("SPINWAVE 2.2 - (c) S. Petit, LLB CEA/CNRS - Sept 2018")
!
!	-----------------------------------------------------------------------
!	initilisation des constantes et des parametres par defaut
!
	pi     = 4.0*atan(1.0)
	d2r    = pi/180.0
	mub    = 9.2740154e-24
	muo    = 4.0*pi*1.0e-7
	kb     = 1.380658e-23
	ec     = 1.602177e-19
	j2unit = 1.0
!
	ax     = 6.28
	ay     = 6.28
	az     = 6.28
	alfa   = pi/2.0
	beta   = pi/2.0
	gama   = pi/2.0
!
	d      = 3
	n      = 1
	nt     = 2
	nz     = 4
	nato   = 0
!
	lu     = 5
	nw     = 100
	wmax   = 20.0
	wmin   = 0.0
	sig    = 0.5
!
	q0     = 0.0
	q      = 0.0
	qmin   = 0.0
	qmax   = 0.0
	en0    = 0.0
	ki     = 2.662
	unit1  = 2.0
!
	couple = .false.
	mf     = .false.
	pho    = .false.
	h0     = 0.0
	niter  = 30
	g      = 0.0
	gdm    = 0.0
	gksea  = 0.0
	gdiag  = 0.0
	temp   = 0.1
	unit   = "MEV"
	gl     = 0.0
	do i=1,3
		gl(i,i) = 1.0
	enddo
	gl2    = gl
!
	affich = .false.
	fform  = .false.
	tof    = .false.
	integ  = .false.
	coupe  = .false.
	coup1d = .false.
	angle  = .false.
	incom  = .false.
	kprop  = 0.0
	ex     = 0.0
	ey     = 0.0
	ez     = 0.0
	ez(3)  = 1.0
	ncalc  = 8
!
	nfibo  = 15
	calc   = 2
	reg1   = 0.01
	reg2   = 0.01
	reg3   = 0.01
!
!	facteur de forme
!	par defaut, celui du Mn3+
!
	fga    =  0.4198
	fa     = 14.2830
	fgb    =  0.6054
	fb     =  5.4689
	fgc    =  0.9241
	fc     = -0.0088
	fgd    = -0.9498
!
	nnx    = nt
	nny    = nt
	nnz    = nt
!
	ffst   = ""
!
	maxime = .false.
	savary = .false.
!
	simpleq = .true.
	reduc   = .false.
	partiel = .false.
!
	twin    = .false.
	ntw     = 0
!
	bbq     = .false.
	cnt     = .false.
	nbk     = 0
!
	np      = 0
	np1     = 0
	np2     = 0
	np3     = 0
!
!	-----------------------------------------------------------------------
!	premiere lecture des donnees
!
	call banner("LECTURE DU FICHIER DE DONNEES")
	rewind(lu)
	do
		read(lu,"(a)",iostat=ios) carte
		write(*,"(a)") carte
		if(ios.ne.0) exit
		if(carte(1:1).eq."#") cycle
		call anacart(carte,nbcle,moc,av,iv,rv)
		do nucle = 1, nbcle
			if(moc(nucle).eq."I") then
			     if(iv(nucle).gt.n) n = iv(nucle)
			endif
			if(moc(nucle).eq."J")      nato     = nato+1
			if(moc(nucle).eq."AX")     ax       = rv(nucle)
			if(moc(nucle).eq."AY")     ay       = rv(nucle)
			if(moc(nucle).eq."AZ")     az       = rv(nucle)
			if(moc(nucle).eq."ALFA")   alfa     = rv(nucle)*d2r
			if(moc(nucle).eq."BETA")   beta     = rv(nucle)*d2r
			if(moc(nucle).eq."GAMA")   gama     = rv(nucle)*d2r
!
			if(moc(nucle).eq."Q0X")    q0(1)    = rv(nucle)
			if(moc(nucle).eq."Q0Y")    q0(2)    = rv(nucle)
			if(moc(nucle).eq."Q0Z")    q0(3)    = rv(nucle)		
			if(moc(nucle).eq."DQX")    dq(1)    = rv(nucle)
			if(moc(nucle).eq."DQY")    dq(2)    = rv(nucle)
			if(moc(nucle).eq."DQZ")    dq(3)    = rv(nucle)
!
			if(moc(nucle).eq."COUPE")  coupe    = .true.
			if(moc(nucle).eq."COUP1D") coup1d   = .true.
			if(moc(nucle).eq."DQ1X")   dq1(1)   = rv(nucle)
			if(moc(nucle).eq."DQ1Y")   dq1(2)   = rv(nucle)
			if(moc(nucle).eq."DQ1Z")   dq1(3)   = rv(nucle)
			if(moc(nucle).eq."DQ2X")   dq2(1)   = rv(nucle)
			if(moc(nucle).eq."DQ2Y")   dq2(2)   = rv(nucle)
			if(moc(nucle).eq."DQ2Z")   dq2(3)   = rv(nucle)
			if(moc(nucle).eq."DQ3X")   dq3(1)   = rv(nucle)
			if(moc(nucle).eq."DQ3Y")   dq3(2)   = rv(nucle)
			if(moc(nucle).eq."DQ3Z")   dq3(3)   = rv(nucle)
			if(moc(nucle).eq."EN0")    en0      = rv(nucle)
			if(moc(nucle).eq."NP")     np       = iv(nucle)
			if(moc(nucle).eq."NP1")    np1      = iv(nucle)
			if(moc(nucle).eq."NP2")    np2      = iv(nucle)
			if(moc(nucle).eq."NP3")    np3      = iv(nucle)
!
			if(moc(nucle).eq."WMAX")   wmax     = rv(nucle)
			if(moc(nucle).eq."WMIN")   wmin     = rv(nucle)
			if(moc(nucle).eq."SIG")    sig      = rv(nucle)
			if(moc(nucle).eq."NW")     nw       = iv(nucle)
			if(moc(nucle).eq."FICH")   fres     = av(nucle)
			if(moc(nucle).eq."FFST")   ffst     = av(nucle)
			if(moc(nucle).eq."AFFICH") affich   = .true.
			if(moc(nucle).eq."REG1")   reg1     = rv(nucle)
			if(moc(nucle).eq."REG2")   reg2     = rv(nucle)
			if(moc(nucle).eq."REG3")   reg3     = rv(nucle)
			if(moc(nucle).eq."CALC")   calc     = iv(nucle)
			if(moc(nucle).eq."NT")     nt       = iv(nucle)
			if(moc(nucle).eq."FFORM")  fform    = .true.
			if(moc(nucle).eq."TOF")    tof      = .true.
			if(moc(nucle).eq."NFIBO")  nfibo    = iv(nucle)
			if(moc(nucle).eq."QMIN")   qmin     = rv(nucle)
			if(moc(nucle).eq."QMAX")   qmax     = rv(nucle)
			if(moc(nucle).eq."INTEG")  integ    = .true.
			if(moc(nucle).eq."UNIT")   unit1    = rv(nucle)
			if(moc(nucle).eq."KI")     ki       = rv(nucle)
!
			if(moc(nucle).eq."EZX")    ez(1)    = rv(nucle)
			if(moc(nucle).eq."EZY")    ez(2)    = rv(nucle)
			if(moc(nucle).eq."EZZ")    ez(3)    = rv(nucle)					
!
			if(moc(nucle).eq."fA")     fga      = rv(nucle)
			if(moc(nucle).eq."fa")     fa       = rv(nucle)
			if(moc(nucle).eq."fB")     fgb      = rv(nucle)
			if(moc(nucle).eq."fb")     fb       = rv(nucle)
			if(moc(nucle).eq."fC")     fgc      = rv(nucle)
			if(moc(nucle).eq."fc")     fc       = rv(nucle)
			if(moc(nucle).eq."fD")     fgd      = rv(nucle)
!
			if(moc(nucle).eq."PHO")    pho      = .true.
			if(moc(nucle).eq."COUPLE") couple   = .true.
			if(moc(nucle).eq."MF")     mf       = .true.
			if(moc(nucle).eq."G")      g        = rv(nucle)
			if(moc(nucle).eq."TEMP")   temp     = rv(nucle)
			if(moc(nucle).eq."NITER")  niter    = iv(nucle)
			if(moc(nucle).eq."HX")     h0(1)    = rv(nucle)
			if(moc(nucle).eq."HY")     h0(2)    = rv(nucle)
			if(moc(nucle).eq."HZ")     h0(3)    = rv(nucle)
			if(moc(nucle).eq."UNIT")   unit     = av(nucle)
			if(moc(nucle).eq."GX")	   gl(1,1)  = rv(nucle)
			if(moc(nucle).eq."GY")     gl(2,2)  = rv(nucle)
			if(moc(nucle).eq."GZ")     gl(3,3)  = rv(nucle)
			if(moc(nucle).eq."GX2")	   gl2(1,1) = rv(nucle)
			if(moc(nucle).eq."GY2")    gl2(2,2) = rv(nucle)
			if(moc(nucle).eq."GZ2")    gl2(3,3) = rv(nucle)
!
			if(moc(nucle).eq."NX ")    nnx      = iv(nucle)
			if(moc(nucle).eq."NY ")    nny      = iv(nucle)
			if(moc(nucle).eq."NZ ")    nnz      = iv(nucle)
!
			if(moc(nucle).eq."INCOM")  incom    = .true.
			if(moc(nucle).eq."KPROPX") kprop(1) = rv(nucle)
			if(moc(nucle).eq."KPROPY") kprop(2) = rv(nucle)
			if(moc(nucle).eq."KPROPZ") kprop(3) = rv(nucle)
!
			if(moc(nucle).eq."MAXIME") maxime   = .true.
			if(moc(nucle).eq."SAVARY") savary   = .true.
!
			if(moc(nucle).eq."REDM")    reduc   = .true.
			if(moc(nucle).eq."PARTIEL") partiel = .true.
			if(moc(nucle).eq."BBQ")     bbq     = .true.
			if(moc(nucle).eq."CONTINUUM") cnt   = .true.
			if(moc(nucle).eq."NBK")     nbk     = iv(nucle)
!
			if(moc(nucle).eq."TWIN") then
				twin = .true.
				ntw  = ntw+1
			endif
!
		enddo
	enddo
!
	if(coupe) then
		if(np.ne.0) then
			np1=np
			np2=np
		endif
	endif
	if(coupe)  simpleq  = .false.
	if(tof)    simpleq  = .false.
!
!	-----------------------------------------------------------------------
!	lecture de la position des spins dans la maille magnetique
!	lecture de l'axe d'orientation e3
!
!	calcul de la base locale de chacun des spins
!	e1,e2,e3 pour le calcul de mat(i,alfa,beta)
!	i indice de spin
!	alfa = x,y,z base orthogonale espace reel
!	beta = 1,2,3 ou 3 = spin statique
!
	allocate(pos(n,d),pos1(n,d),sp(n),phi(n),theta(n),ma(n))
	allocate(matl(n,d,d),matr(n,d,d),zdr(n,d),zdi(n,d),eta(n,d))
	allocate(nom(n),wj(n),gj(n),b20(n),b22(n),b2m2(n),b21(n),b2m1(n),b40(n),b42(n),b43(n),b44(n))
	allocate(b60(n),b62(n),b63(n),b64(n),b66(n),dis0(n))
	allocate(matlx(n,d,d),matrx(n,d,d),xcef(n,d),jmy(n,d))
	allocate(nu(n))
!
	pos   = 0.0
	phi   = 0.0
	theta = 0.0
	ma    = 0.0
	matl  = 0.0
	matr  = 0.0
	mat   = 0.0
!
	matlx = 0.0
	matrx = 0.0
	b20   = 0.0
	b22   = 0.0
	b2m2  = 0.0
	b21   = 0.0
	b2m1  = 0.0
	b40   = 0.0
	b42   = 0.0
	b43   = 0.0
	b44   = 0.0
	b60   = 0.0
	b62   = 0.0
	b63   = 0.0
	b64   = 0.0
	b66   = 0.0
	dis0  = 0.0
!
	rewind(lu)
	do
		read(lu,"(a)",iostat=ios) carte
		if(ios.ne.0) exit
		if(carte(1:1).eq."#") cycle
		call anacart(carte,nbcle,moc,av,iv,rv)
		trouve = .false.
		angle  = .false.
		do nucle = 1, nbcle
			if(moc(nucle).eq."I") then
				i      = iv(nucle)
				trouve = .true.
			endif
		enddo
		if(trouve) then
			do nucle = 1, nbcle
				if((moc(nucle).eq."PHI").or.(moc(nucle).eq."THETA")) then
				angle = .true.
				endif
			enddo
		endif
		ss = 0.0
		if(trouve) then
			do nucle = 1, nbcle
				if(moc(nucle).eq."X")     pos(i,1)  = rv(nucle)
				if(moc(nucle).eq."Y")     pos(i,2)  = rv(nucle)
				if(moc(nucle).eq."Z")     pos(i,3)  = rv(nucle)
				if(moc(nucle).eq."PHI")   phi(i)    = rv(nucle)
				if(moc(nucle).eq."THETA") theta(i)  = rv(nucle)
				if(moc(nucle).eq."SX")    jmy(i,1)  = rv(nucle)
				if(moc(nucle).eq."SY")    jmy(i,2)  = rv(nucle)
				if(moc(nucle).eq."SZ")    jmy(i,3)  = rv(nucle)
				if(moc(nucle).eq."MA")    ma(i)     = rv(nucle)
				if(moc(nucle).eq."NOM")   nom(i)    = av(nucle)
				if(moc(nucle).eq."WJ")    wj(i)     = rv(nucle)
				if(moc(nucle).eq."B20")   b20(i)    = rv(nucle)
				if(moc(nucle).eq."B22")   b22(i)    = rv(nucle)
				if(moc(nucle).eq."B2M2")  b2m2(i)   = rv(nucle)
				if(moc(nucle).eq."B21")   b21(i)    = rv(nucle)
				if(moc(nucle).eq."B2M1")  b2m1(i)   = rv(nucle)
				if(moc(nucle).eq."B40")   b40(i)    = rv(nucle)
				if(moc(nucle).eq."B42")   b42(i)    = rv(nucle)
				if(moc(nucle).eq."B43")   b43(i)    = rv(nucle)
				if(moc(nucle).eq."B44")   b44(i)    = rv(nucle)
				if(moc(nucle).eq."B60")   b60(i)    = rv(nucle)
				if(moc(nucle).eq."B62")   b62(i)    = rv(nucle)
				if(moc(nucle).eq."B63")   b63(i)    = rv(nucle)
				if(moc(nucle).eq."B64")   b64(i)    = rv(nucle)
				if(moc(nucle).eq."B66")   b66(i)    = rv(nucle)
				if(moc(nucle).eq."DD")    dis0(i)   = rv(nucle)
				if(moc(nucle).eq."CX")    xcef(i,1) = rv(nucle)
				if(moc(nucle).eq."CY")    xcef(i,2) = rv(nucle)
				if(moc(nucle).eq."CZ")    xcef(i,3) = rv(nucle)
			enddo
			if(angle.eqv..false.) then
				ss = jmy(i,:)
				s  = sqrt(sum(ss*ss))
				if(s.eq.0.0) then
					write(*,*) "DIRECTION NON DEFINIE POUR LE SPIN",i,ss
				else
					ss       = ss/sqrt(sum(ss*ss))
					theta(i) = acos(ss(3))
					if((theta(i).eq.0.0).or.(theta(i).eq.pi)) then
						phi(i) = 0.0
					else
						
						if(1.-abs(ss(1)/sin(theta(i))).lt.(1e-4)) then
							if(abs(ss(1)/sin(theta(i))-1.).lt.(1e-4)) phi(i)=0.0
							if(abs(ss(1)/sin(theta(i))+1.).lt.(1e-4)) phi(i)=pi/d2r
						else
							if(ss(2).lt.0.0) phi(i) = (2*pi-acos(ss(1)/sin(theta(i))))/d2r
							if(ss(2).ge.0.0) phi(i) = (     acos(ss(1)/sin(theta(i))))/d2r
						endif
					endif
					theta(i) = theta(i)/d2r
				endif
			else
				jmy(i,1) = cos(phi(i)*d2r)*sin(theta(i)*d2r)
				jmy(i,2) = sin(phi(i)*d2r)*sin(theta(i)*d2r)
				jmy(i,3) = cos(theta(i)*d2r)
			endif			
		endif
	enddo
!
!	Axes locaux de champ cristallin
!	matlx: coordonnees  cartesiennes -> base locale
!	matrx: coordonnees  base locale  -> cartesiennes
!
	do i=1,n
!		vecteur e3 de la base
		matlx(i,3,:) = xcef(i,:)/sqrt(sum(xcef(i,:)*xcef(i,:)))
!		vecteur e1 de la base
!		on choisit a priori e1 dans le plan qui contient ez et e3
!		ss est le vecteur ez
		ss    = 0.0
		ss(3) = 1.0
		y     = 1.0
		if(matlx(i,3,3).lt.0.0) y=-1.0
!		x = ez.xcef
		x     = sum(ss*matlx(i,3,:))
		if(x.ne.1.0) then
			ss    = ss - x*matlx(i,3,:)
			matlx(i,1,:) = y*ss/sqrt(sum(ss*ss))
!			e2 complete le triedre direct
			call vct(matlx(i,3,:),matlx(i,1,:),matlx(i,2,:))
		else
!		sinon, xcef est l'axe z; on prend alors le repere std
		matlx(i,1,1) = 1.0
		matlx(i,2,2) = 1.0
		endif
!		matrice inverse
		matrx(i,:,:) = transpose(matlx(i,:,:))
	enddo
!
!	parametres de champ cristallin
!
	do i=1,n
		call paramRE(nom(i),wj(i),gj(i),alf,bet,gam)
		b20(i)  = b20(i) *alf/2.
		b22(i)  = b22(i) *alf/2.
		b2m2(i) = b2m2(i)*alf/2.
		b21(i)  = b21(i) *alf/2.
		b2m1(i) = b2m1(i)*alf/2.
		b40(i)  = b40(i)*bet/8
		b43(i)  = b43(i)*bet*(-sqrt(35.)/2.)
		b60(i)  = b60(i)*gam/16.
		b63(i)  = b63(i)*gam*(-sqrt(105.)/8.)
		b66(i)  = b66(i)*gam*( sqrt(231.)/16.)
		nu(i)   = int(2*wj(i)+1)
	enddo
!
!	-----------------------------------------------------------------------
!	unites pour la temperature et le champ
!
	call banner("UNITS")
	select case(unit)
	case("K")
		j2unit = 1.0/kb
		k2unit = 1.0
	case("MEV")
		j2unit = 1000./ec
		k2unit = kb*j2unit
	case("THZ")
		j2unit = 1000./ec
		k2unit = kb*j2unit
	case default
		mub = 1.0
		muo = 1.0
	end select
	temp = temp*k2unit
	write(*,"(a8)") unit
	write(*,"(a8,e12.5)") "j2unit ",j2unit
	write(*,"(a8,e12.5)") "k2unit ",k2unit
	write(*,"(a8,e12.5)") "mub    ",mub
!
!	-----------------------------------------------------------------------
!	reseau reel et reseau reciproque
!	
	call reseau(ax,ay,az,alfa,beta,gama,vol,e,ee)
	e  = transpose(e)
	ee = transpose(ee)
!
!	-----------------------------------------------------------------------
!	vecteur z du laboratoire, perp au plan de diffusion
!
	ez = matmul(e,ez)
	ez = ez/sqrt(sum(ez*ez))
!
!	-----------------------------------------------------------------------
!	vecteurs espace reel pour sommation de Fourier
!
	ntt = (2*nt+1)**3
	allocate(listvec(ntt,d))
	do l1 = -nt,nt,1
	do l2 = -nt,nt,1
	do l3 = -nt,nt,1
		l = 1 + l3+nt + (2*nt+1)*(l2+nt + (2*nt+1)*(l1+nt))
		listvec(l,1) = l1
		listvec(l,2) = l2
		listvec(l,3) = l3
	enddo
	enddo
	enddo
!
!	-----------------------------------------------------------------------
!	lecture des atomes suppl\E9mentaires
!
	if(nato.ne.0) then
		allocate(lato(nato))
		do i=1,nato
			lato(i)%nam = ""
			lato(i)%pos = 0.0
		enddo
		j = 0
		rewind(lu)
		do
			read(lu,"(a)",iostat=ios) carte
			if(ios.ne.0) exit
			if(carte(1:1).eq."#") cycle
			call anacart(carte,nbcle,moc,av,iv,rv)
			trouve = .false.
			do nucle = 1, nbcle
				if(moc(nucle).eq."J") then
					j           = j+1
					trouve      = .true.
				endif
			enddo
			if(trouve) then
				do nucle = 1, nbcle
					if(moc(nucle).eq."J") lato(j)%nam    = av(nucle)
					if(moc(nucle).eq."X") lato(j)%pos(1) = rv(nucle)
					if(moc(nucle).eq."Y") lato(j)%pos(2) = rv(nucle)
					if(moc(nucle).eq."Z") lato(j)%pos(3) = rv(nucle)
				enddo
			endif
		enddo
	endif
!
!	-----------------------------------------------------------------------
!	lecture des interactions magnetiques
!
	allocate(dech(n,d),dval(n),dani(n,d,d),jech(n,n,nz,d,d),ctef(n,n,nz,d,d))
	allocate(dist(nz,n,n),disf(nz,n,n))
	allocate(modele(n,n))
	allocate(h(n))
	allocate(uprop(n,d))
	dech   = 0.0
	dval   = 0.0
	dani   = 0.0
	jech   = 0.0
	ctef   = 0.0
	dist   = 0.0
	disf   = 0.0
	h      = 1.0e-4
	modele = ""
	uprop  = 0.0
!
	rewind(lu)
	do
		read(lu,"(a)",iostat=ios) carte
		if(ios.ne.0) exit
		if(carte(1:1).eq."#") cycle
		call anacart(carte,nbcle,moc,av,iv,rv)
		trouve = .false.
		do nucle = 1, nbcle
			if(moc(nucle).eq."I1") trouve = .true.
		enddo
		if(trouve) then
			do nucle = 1, nbcle
				if(moc(nucle).eq."I1")  then
					if(av(nucle).eq."ALL") then
						l1  = 0.0
					else
						l1  = iv(nucle)
					endif
				endif
				if(moc(nucle).eq."I2")  then
					if(av(nucle).eq."ALL") then
						l2  = 0.0
					else
						l2  = iv(nucle)
					endif
				endif
			enddo
			do nucle = 1, nbcle
				if(moc(nucle).eq."MODELE") then
					if(l1.ne.0) then
						if(l2.ne.0) then
							modele(l1,l2) = av(nucle)
							modele(l2,l1) = av(nucle)
						else
							modele(l1,:)  = av(nucle)
							modele(:,l1)  = av(nucle)
						endif
					else
						if(l2.ne.0) then
							modele(:,l2)  = av(nucle)
							modele(l2,:)  = av(nucle)
						else
							modele(:,:)   = av(nucle)
							modele(:,:)   = av(nucle)
						endif
					endif
				endif
			enddo
			do nucle = 1, nbcle
				if(moc(nucle)(1:1).eq."J") then
					k=0
					l=0
					if(moc(nucle)(2:2).eq."1")  l=1
					if(moc(nucle)(2:2).eq."2")  l=2
					if(moc(nucle)(2:2).eq."3")  l=3
					if(moc(nucle)(2:2).eq."4")  l=4
					select case (moc(nucle)(3:4))
						case ("XX")
							k =1
							k1=1
						case ("YY")
							k =2
							k1=2
						case ("ZZ")
							k =3
							k1=3
						case ("XY")
							k =1
							k1=2
						case ("XZ")
							k =1
							k1=3
						case ("YZ")
							k =2
							k1=3
						case ("YX")
							k =2
							k1=1
						case ("ZX")
							k =3
							k1=1
						case ("ZY")
							k =3
							k1=2
						case default
							k=0
					end select
					if(l.ne.0) then
					if(k.eq.0) then
						do i=1,d
							if(l2.ne.0) then
								if(l1.ne.0) then
								jech(l1,l2,l,i,i) = rv(nucle)
								jech(l2,l1,l,i,i) = rv(nucle)
								else
								jech( :,l2,l,i,i) = rv(nucle)
								jech(l2, :,l,i,i) = rv(nucle)
								endif
							else
								if(l1.ne.0) then
								jech(l1,:,l,i,i)  = rv(nucle)
								else
								jech( :,:,l,i,i)  = rv(nucle)
								endif
							endif
						enddo
					else
						if(l2.ne.0) then
							if(l1.ne.0) then
							jech(l1,l2,l,k,k1) = rv(nucle)
							jech(l2,l1,l,k,k1) = rv(nucle)
							else
							jech( :,l2,l,k,k1) = rv(nucle)
							jech(l2, :,l,k,k1) = rv(nucle)
							endif
						else
							if(l1.ne.0) then
							jech(l1,:,l,k,k1)  = rv(nucle)
							else
							jech( :,:,l,k,k1)  = rv(nucle)
							endif
						endif
					endif
					endif
				endif
				if(moc(nucle)(1:1).eq."D") then
					l=0
					if(moc(nucle)(2:2).eq."1")  l=1
					if(moc(nucle)(2:2).eq."2")  l=2
					if(moc(nucle)(2:2).eq."3")  l=3
					if(moc(nucle)(2:2).eq."4")  l=4
					if(l.ne.0) then
						if(l2.ne.0) then
							if(l1.ne.0) then
							dist(l,l1,l2  ) = rv(nucle)
							dist(l,l2,l1  ) = rv(nucle)
							else
							dist(l, :,l2  ) = rv(nucle)
							dist(l,l2, :  ) = rv(nucle)
							endif
						else
							if(l1.ne.0) then
							dist(l,l1,:   ) = rv(nucle)
							else
							dist(l, :,:   ) = rv(nucle)
							endif
						endif
					endif
				endif
				if(moc(nucle).eq."H") then
					if(l1.ne.0) then
						h(l1) = rv(nucle)
					else
						h( :) = rv(nucle)
					endif
				endif
				if(moc(nucle).eq."DV") then
					if(l1.ne.0) then
						dval(l1) = rv(nucle)
					else
						dval( :) = rv(nucle)
					endif
				endif
				if(moc(nucle).eq."UPX") then
					if(l1.ne.0) then
						uprop(l1,1) = rv(nucle)
					else
						uprop(:,1)= rv(nucle)
					endif
				endif
				if(moc(nucle).eq."UPY") then
					if(l1.ne.0) then
						uprop(l1,2) = rv(nucle)
					else
						uprop(:,2)= rv(nucle)
					endif
				endif
				if(moc(nucle).eq."UPZ") then
					if(l1.ne.0) then
						uprop(l1,3) = rv(nucle)
					else
						uprop(:,3)= rv(nucle)
					endif
				endif
			enddo
!
			do nucle = 1, nbcle
				if(moc(nucle)(1:1).eq."F") then
					k=0
					l=0
					if(moc(nucle)(2:2).eq."1")  l=1
					if(moc(nucle)(2:2).eq."2")  l=2
					if(moc(nucle)(2:2).eq."3")  l=3
					if(moc(nucle)(2:2).eq."4")  l=4
					if(moc(nucle)(3:4).eq."XX") k=1
					if(moc(nucle)(3:4).eq."YY") k=2
					if(moc(nucle)(3:4).eq."ZZ") k=3
					if(l.ne.0) then
					if(k.eq.0) then
						do i=1,d
							if(l2.ne.0) then
								if(l1.ne.0) then
									ctef(l1,l2,l,i,i) = rv(nucle)
									ctef(l2,l1,l,i,i) = rv(nucle)
								else
									ctef( :,l2,l,i,i) = rv(nucle)
									ctef(l2, :,l,i,i) = rv(nucle)
								endif
							else
								if(l1.ne.0) then
									ctef(l1,:,l,i,i)  = rv(nucle)
								else
									ctef( :,:,l,i,i)  = rv(nucle)
								endif
							endif
						enddo
					else
						if(l2.ne.0) then
							if(l1.ne.0) then
							ctef(l1,l2,l,k,k) = rv(nucle)
							ctef(l2,l1,l,k,k) = rv(nucle)
							else
							ctef( :,l2,l,k,k) = rv(nucle)
							ctef(l2, :,l,k,k) = rv(nucle)
							endif
						else
							if(l1.ne.0) then
							ctef(l1,:,l,k,k)  = rv(nucle)
							else
							ctef( :,:,l,k,k)  = rv(nucle)
							endif
						endif
					endif
					endif
				endif
				if(moc(nucle)(1:1).eq."L") then
					l=0
					if(moc(nucle)(2:2).eq."1")  l=1
					if(moc(nucle)(2:2).eq."2")  l=2
					if(moc(nucle)(2:2).eq."3")  l=3
					if(moc(nucle)(2:2).eq."4")  l=4
					if(l.ne.0) then
						if(l2.ne.0) then
							if(l1.ne.0) then
							disf(l,l1,l2  ) = rv(nucle)
							disf(l,l2,l1  ) = rv(nucle)
							else
							disf(l, :,l2  ) = rv(nucle)
							disf(l,l2, :  ) = rv(nucle)
							endif
						else
							if(l1.ne.0) then
							disf(l,l1,:   ) = rv(nucle)
							else
							disf(l, :,:   ) = rv(nucle)
							endif
						endif
					endif
				endif
			enddo
		endif
	enddo
!
!	anisotropie
!	le vecteur qui definit l'anisotropie est confondu avec l'axe du CEF
!	O20  = 3Jz^2  + Cte
!	O22  = Jx^2-Jy^2
!	O2m2 = JxJy
!	O2m2 = 
	do i=1,n
		if(dval(i).ne.0.0) b20(i)  = dval(i)/3.0

		dani(i,1,1) =  b22(i)
		dani(i,2,2) = -b22(i)
		dani(i,3,3) = 3.0*b20(i)
		dani(i,1,2) =  b2m2(i)
		dani(i,2,1) =  b2m2(i)
		dani(i,1,3) =  b21(i)
		dani(i,3,1) =  b21(i)
		dani(i,2,3) =  b2m1(i)
		dani(i,3,2) =  b2m1(i)
!
!		on r\E9\E9crit le tenseur dans le rep\E8re global			
		dani(i,:,:) = matmul(transpose(matlx(i,:,:)),matmul(dani(i,:,:),matlx(i,:,:)))
	enddo
!
!	-----------------------------------------------------------------------
!	lecture des interactions magnetiques
!	type Dzyaloshinskii-Moryia
!
	allocate(dm(n,n,ntt,d),pt(n,n,ntt,d))
	dm = 0.0
	pt = 0.0
	rewind(lu)
	do
		read(lu,"(a)",iostat=ios) carte
		if(ios.ne.0) exit
		if(carte(1:1).eq."#") cycle
		call anacart(carte,nbcle,moc,av,iv,rv)
		do nucle = 1, nbcle
			if(moc(nucle).eq."GDM")   gdm    = rv(nucle)
			if(moc(nucle).eq."GKSEA") gksea  = rv(nucle)
			if(moc(nucle).eq."GDIAG") gdiag  = rv(nucle)
		enddo
		trouve = .false.
		do nucle = 1, nbcle
			if(moc(nucle).eq."UX")  trouve = .true.
			if(moc(nucle).eq."UY")  trouve = .true.
			if(moc(nucle).eq."UZ")  trouve = .true.
		enddo
		if(trouve) then
			udm = 0
			do nucle = 1, nbcle
				if(moc(nucle).eq."I1")  l1  = iv(nucle)
				if(moc(nucle).eq."I2")  l2  = iv(nucle)
				if(moc(nucle).eq."UX")	udm(1) = rv(nucle)
				if(moc(nucle).eq."UY")	udm(2) = rv(nucle)
				if(moc(nucle).eq."UZ")	udm(3) = rv(nucle)
			enddo
			i = 0
			j = 0
			do l=1,ntt
				x = sum((listvec(l,:)-udm)**2)
				if(x.eq.0.0) i=l
				x = sum((listvec(l,:)+udm)**2)
				if(x.eq.0.0) j=l
			enddo
			if(i.eq.0)  write(*,*) "Erreur i =0"
			if(l1.eq.0) write(*,*) "Erreur l1=0"
			if(l2.eq.0) write(*,*) "Erreur l2=0"
			do nucle = 1, nbcle
				if(moc(nucle).eq."DMX")	dm(l1,l2,i,1) = rv(nucle)
				if(moc(nucle).eq."DMY")	dm(l1,l2,i,2) = rv(nucle)
				if(moc(nucle).eq."DMZ")	dm(l1,l2,i,3) = rv(nucle)
				if(moc(nucle).eq."PT")  nom1          = av(nucle)
				k = 0
				do l=1,nato
					if(nom1.eq.(lato(l)%nam)) then
						pt(l1,l2,i,:) = lato(l)%pos
						k             = l
					endif
				enddo
			enddo
			if(k.eq.0) then
				write(*,*) "PONT NON TROUVE POUR",l1,l2,i,udm,nom1
			endif
			dm(l2,l1,j,:) = -dm(l1,l2,i,:)
			pt(l2,l1,j,:) =  pt(l1,l2,i,:)-udm
		endif
	enddo
!
!	-----------------------------------------------------------------------
!	lecture des interactions biquadratiques
!
	if(bbq) then
		allocate(biq(n,n,n,n),dbiq12(n,n),dbiq13(n,n))
		biq    = 0.0
		dbiq12 = 0.0
		dbiq13 = 0.0
		rewind(lu)
		do
			read(lu,"(a)",iostat=ios) carte
			if(ios.ne.0) exit
			if(carte(1:1).eq."#") cycle
			call anacart(carte,nbcle,moc,av,iv,rv)
			trouve = .false.
			do nucle = 1, nbcle
				if(moc(nucle).eq."BIQ") trouve = .true.
			enddo
			if(trouve) then
				p1 = 0.
				p2 = 0.
				p4 = 0.
				do nucle = 1, nbcle
					if(moc(nucle).eq."I1")  l1  = iv(nucle)
					if(moc(nucle).eq."I2")  l2  = iv(nucle)
					if(moc(nucle).eq."I3")  l3  = iv(nucle)
					if(moc(nucle).eq."I4")  l4  = iv(nucle)
					if(moc(nucle).eq."D12") p1  = rv(nucle)
					if(moc(nucle).eq."D13") p2  = rv(nucle)
					if(moc(nucle).eq."BIQ") p4  = rv(nucle)
				enddo
				biq(l1,l2,l3,l4) = p4
				biq(l1,l2,l4,l3) = p4
				biq(l2,l1,l3,l4) = p4
				biq(l2,l1,l4,l3) = p4
				biq(l3,l4,l1,l2) = p4
				biq(l4,l3,l1,l2) = p4
				biq(l3,l4,l2,l1) = p4
				biq(l4,l3,l2,l1) = p4
				dbiq12(l1,l2) = p1
				dbiq12(l2,l1) = p1
				dbiq13(l1,l3) = p2
				dbiq13(l3,l1) = p2
			endif
		enddo
	endif
!
!	-----------------------------------------------------------------------
!	normalisation dans le cas incommensurable
!	uprop doit etre donn\E9 en cartesiennes
!
	if(incom) then
		do i=1,n
			x = sqrt(sum(uprop(i,:)*uprop(i,:)))
			uprop(i,:) = uprop(i,:)/x
		enddo
	endif
!
!	-----------------------------------------------------------------------
!	lecture des twins
!
	if(twin) then
		allocate(rtw(ntw,d),atw(ntw))
		ntw = 0
		rewind(lu)
		do
			read(lu,"(a)",iostat=ios) carte
			if(ios.ne.0) exit
			if(carte(1:1).eq."#") cycle
			call anacart(carte,nbcle,moc,av,iv,rv)
			do nucle = 1, nbcle
				if(moc(nucle).eq."TWIN") ntw = ntw+1
				if(moc(nucle).eq."ATW") atw(ntw) = rv(nucle)
				if(moc(nucle).eq."TWX") rtw(ntw,1) = rv(nucle)
				if(moc(nucle).eq."TWY") rtw(ntw,2) = rv(nucle)
				if(moc(nucle).eq."TWZ") rtw(ntw,3) = rv(nucle)
			enddo
		enddo
	endif
!
!	---------------------------
!	rustine pour Maxime et MnGe
!
	if(maxime) then
		do i=1,n
		do j=1,n
			x= sum(jech(i,j,1,:,:))
			if(x.ne.0.0) then
				do l=1,ntt
					r0 = pos(j,:)+listvec(l,:)-pos(i,:)
					r0 = matmul(e,r0)
					z  = sqrt(sum(r0*r0))
					if((z.le.dist(1,i,j)).and.(z.gt.0.0)) then
						dm(j,i,l,:) = r0/sqrt(sum(r0*r0))
					endif
				enddo
			endif
		enddo
		enddo
	endif
!
!	--------------------------------------
!	rustine pour integrales echange savary
!
	if(savary) then
		do i=1,n
		do j=1,n
			modele(i,j) = "ANISOBV"
			jpmpm = jech(i,j,1,1,1)
			jpm   = jech(i,j,2,1,1)
			jzpm  = jech(i,j,3,1,1)
			jzz   = jech(i,j,4,1,1)
			call cvsavary(gl2,jpmpm,jpm,jzpm,jzz,jech(i,j,1,1,1),jech(i,j,2,1,1),jech(i,j,3,1,1),jech(i,j,4,1,1))
		enddo
		enddo
	endif
!
!	--------------------------------------------
!	recherche de la configuration du fondamental
!	H    = Hel + Hmag + V
!	Hel  = 1/2 G_ij (u_i-u_j)^2
!	Hmag = J_ij S_i.S_j + Hcef
!	V    = g K_ij.(u_i-u_j) S_i.S_j
!	Hypothese de couplage K_ij = selon vecteur norme (u_i-u_j)
!
	if(mf) then
		call banner("OPTION MF : Recherche de la configuration du fondamental")
!
		allocate(nrj(niter))
		nrj = 0.0
!
!		premiere iteration
		p    = 0
		write(*,"(a10,i4)") "Iteration ",p
		do i=1,n
			write(*,"(i4,13f8.3)") i,(jmy(i,k),k=1,d),(pos(i,k),k=1,d)
		enddo

!		boucle sur iterations
		do p=1,niter
!
!			traitement des positions atomiques
!			u_i sum_j G_ij - G_ij uj + u_i (sum_j g K_ij S_i S_j) = 0
			if(couple) then
				allocate(ar(n*d+d,n*d),u1(n*d+d))
				ar = 0.0
				u1 = 0.0
				ar = 0.0
				qc = 0.0
!				Terme elastique et terme de couplage
				do i=1,n
				do j=1,n
				do l=1,ntt
					ii = (i-1)*d
					jj = (j-1)*d
					call for(i,j,listvec(l,:),pos,e,nz,modele(i,j),disf(:,i,j),ctef(i,j,:,:,:),qc,jm)
					ar(ii+1:ii+d,jj+1:jj+d) = ar(ii+1:ii+d,jj+1:jj+d) - jm
					ar(ii+1:ii+d,ii+1:ii+d) = ar(ii+1:ii+d,ii+1:ii+d) + jm
					call cpl(i,j,listvec(l,:),pos,e,nz,modele(i,j),disf(:,i,j),ctef(i,j,:,:,:),qc,uo,uu)
!					write(*,*) i,j,uu
					u1(ii+1:ii+d) = u1(ii+1:ii+d) - g*uu*sum(jmy(i,:)*jmy(j,:))
				enddo
				enddo
				enddo
!				On impose la conservation du centre de masse pour regulariser le systeme
!				on ajoute donc une n+1 eme ligne : sum_j m_j u_j = 0
				i  = n+1
				ii = (i-1)*d
				do j=1,n
					jj  = (j-1)*d
					jm0 = 0.0
					do k=1,d
						jm0(k,k) = ma(j)
					enddo
					ar(ii+1:ii+d,jj+1:jj+d) = jm0
				enddo
				u1(ii+1:ii+d) = 0.0
				do i=1,(n+1)*d
					write(*,"(15f8.3)") (ar(i,j),j=1,n*d)
				enddo
				do i=1,(n+1)*d
					write(*,"(f8.3)") u1(i)
				enddo
!
!				resolution du systeme
!				ar = matrice de dimension (n*d+d) lignes x (n*d) colonnes
!				ai = T(ar) ar => matrice de dimension (n*d) lignes x (n*d) colonnes
!				u1 = vecteur de (n*d+d) lignes
!				u2 = T(ar) u1 = vecteur de (n*d) lignes
				allocate(ai(n*d,n*d),u2(n*d))
				do i=1,n*d
				do j=1,n*d
					ai(i,j) = sum(ar(:,i)*ar(:,j))
				enddo
				enddo
				do i=1,n*d
					u2(i) = sum(ar(:,i)*u1)
				enddo
				deallocate(u1)
				allocate(u1(n*d))
				call solve(n*d,ai,u2,u1)
				do i=1,n
					ii = (i-1)*d
					pos(i,:) = pos(i,:)+u1(ii+1:ii+d)
				enddo
				deallocate(ar,ai,u1,u2)
			endif
!
!			champ moyen sur les variables de spin
			do i=1,n
!				calcul du champ moleculaire
				ss = 0.0
				do j=1,n
				do l=1,ntt
					x  = 0.0
					uu = 0.0
					uo = 0.0
!					terme en J_ij S_i.S_j
					detail = .false.
!					x et uu remplacent dval et dech
!					seuls sont pris en compte les couplages pour des distances non nulles
!
					allocate(u1(3),u2(3))
					u1 = matrx(i,:,3)
					u2 = matrx(j,:,3)
					v  = jmy(j,:)
					if(incom) then
						x  = sum(2.0*pi*kprop*listvec(l,:))
						call rodrigue(x,uprop(j,:),matrx(j,:,3),u2)
						call rodrigue(x,uprop(j,:),jmy(j,:),v)
					endif
!					write(*,*) v
!
					x=0.0
					call mag(i,j,listvec(l,:),pos,u1,u2,e,nz,modele(i,j),&
					x*dani,dist(:,i,j),jech(i,j,:,:,:),qc,detail,jm)
					deallocate(u1,u2)
!
!					terme en V = K_ij (u_i-u_j) S_i.S_j
					call cpl(i,j,listvec(l,:),pos,e,nz,modele(i,j),&
					disf(:,i,j),ctef(i,j,:,:,:),qc,uo,uu)
!
!					on somme les termes echange et r\E9seau
					ss = ss + matmul(jm,v)+g*uo*v
!
!					terme DM
					call vct(v,dm(i,j,l,:),vdmr)
					ss = ss + gdm*vdmr

!					terme KSEA
					jm = 0.0
					do l1=1,d
					do l2=1,d
						jm(l1,l2) = gksea*dm(i,j,l,l1)*dm(i,j,l,l2)
					enddo
					enddo
					do l1=1,d
						jm(l1,l1) = 0.0
					enddo
					ss = ss + matmul(jm,v)
!
!					terme GDIAG
					jm = 0.0
					do l1=1,d
						jm(l1,l1) = 1.0
					enddo
					ss = ss + gdiag*sum(dm(i,j,l,:)**2)*matmul(jm,v)					
!
!					terme biquadratique
!
					if(bbq) then
					jm = 0.0

					allocate(u1(3),u2(3))
					rc = pos(i,:)+listvec(l,:)-pos(j,:)
					rc = matmul(e,rc)
					if(sqrt(sum(rc*rc)).le.dbiq12(i,j)) then
					do k3=1,n
					do l3=1,ntt
						rc = pos(k3,:)-pos(i,:)+listvec(l3,:)
						rc = matmul(e,rc)
						if(sqrt(sum(rc*rc)).le.dbiq13(i,k3)) then
						do k4=1,n
						do l4=1,ntt
							rc = pos(k3,:)-pos(k4,:)+listvec(l3,:)-listvec(l4,:)
							rc = matmul(e,rc)
							if(sqrt(sum(rc*rc)).le.dbiq12(k3,k4)) then
					
							u1 = jmy(k3,:)
							u2 = jmy(k4,:)
							if(incom) then
								y  = sum(2.0*pi*kprop*listvec(l3,:))
								call rodrigue(y,uprop(k3,:),jmy(k3,:),u1)
								y  = sum(2.0*pi*kprop*listvec(l4,:))
								call rodrigue(y,uprop(k4,:),jmy(k4,:),u2)
							endif
							x = 2.0*biq(i,j,k3,k4)*sum(u1*u2)
							do l1=1,d
							jm(l1,l1) = jm(l1,l1)+x
							enddo
							endif
						enddo
						enddo
						endif
					enddo
					enddo
					endif
					if(sqrt(sum(rc*rc)).le.dbiq13(i,j)) then
					do k3=1,n
					do l3=1,ntt
						rc = pos(k3,:)-pos(i,:)+listvec(l3,:)
						rc = matmul(e,rc)
						if(sqrt(sum(rc*rc)).le.dbiq12(i,k3)) then
						do k4=1,n
						do l4=1,ntt
							rc = pos(k4,:)-pos(j,:)+listvec(l4,:)-listvec(l,:)
							rc = matmul(e,rc)
							if(sqrt(sum(rc*rc)).le.dbiq12(j,k4)) then
							u1 = jmy(k3,:)
							u2 = jmy(k4,:)
							if(incom) then
								y  = sum(2.0*pi*kprop*listvec(l3,:))
								call rodrigue(y,uprop(k3,:),jmy(k3,:),u1)
								y  = sum(2.0*pi*kprop*listvec(l4,:))
								call rodrigue(y,uprop(k4,:),jmy(k4,:),u2)
							endif
							do l1=1,d
							do l2=1,d
							jm(l1,l2) = jm(l1,l2)+4.0*biq(i,k3,j,k4)*u1(l1)*u2(l2)
							enddo
							enddo
							endif
						enddo
						enddo
						endif
					enddo
					enddo
					endif
					deallocate(u1,u2)
					ss = ss + matmul(jm,v)
					endif		
!
				enddo
				enddo

!				write(*,"(6f12.5)") gj(i)*mub*j2unit*matmul(matlx(i,:,:),h0),ss
				hmol = gj(i)*mub*j2unit*matmul(gl,matmul(matlx(i,:,:),h0))+matmul(matlx(i,:,:),ss)
!
!				allocation memoire pour champ moyen
				nui = nu(i)
				allocate(hr(nui,nui),hi(nui,nui),df(nui),vfr(nui,nui),vfi(nui,nui))
				allocate(jr(d,nui,nui),ji(d,nui,nui),vm(nui))
!				champ cristallin
				dis = 0.0
				call cef(wj(i),b20(i),b22(i),b2m2(i),b21(i),b2m1(i),b40(i),b42(i),b43(i),b44(i),&
				b60(i),b62(i),b63(i),b64(i),b66(i),dis0(i),dis,hmol,hr,hi)		
!				diagonalisation
				call diagocplx_sym(nui,hr,hi,df,vfr,vfi)
!				calcul des valeurs moyennes
				z = sum(exp(-(df(:)-df(1))/temp))
				do l1=1,nui
					vm(l1) = exp(-(df(l1)-df(1))/temp)
				enddo
				vm     = vm/z
				nrj(p) =  nrj(p) + sum(df*vm)
!				calcul de la valeur moyenne de J
!				sum P*_ac P_bc exp(-Ec/T)/z
				call jz(nui,wj(i),jr,ji)
				allocate(u2r(nui,nui),u2i(nui,nui))
				do l1=1,nui
				do l2=1,nui
					u2r(l1,l2) = sum(( vfr(l1,:)*vfr(l2,:) + vfi(l1,:)*vfi(l2,:))*vm)
					u2i(l1,l2) = sum((-vfi(l1,:)*vfr(l2,:) + vfr(l1,:)*vfi(l2,:))*vm)
				enddo
				enddo
				do k=1,d
					jmy(i,k) = sum(jr(k,:,:)*u2r(:,:)-ji(k,:,:)*u2i(:,:))
				enddo
				deallocate(u2r,u2i)
!				transformation base locale -> cartesiennes
				jmy(i,:) = matmul(matrx(i,:,:),jmy(i,:))
!				deallocation de la memoire
				deallocate(hr,hi,df,vfr,vfi,jr,ji,vm)
!
			enddo
			nrj(p) = nrj(p)/n
			write(*,"(a10,i4)") "Iteration ",p
			do i=1,n
				write(*,"(i4,12f8.3)") i,(jmy(i,k),k=1,d),(pos(i,k),k=1,d)
			enddo
		enddo
!		fin des iterations
		write(*,*) "Evolution de l energie"
		do p=1,niter
			write(*,"(a10,i4,e15.5)") "Iteration ",p,nrj(p)
		enddo
		deallocate(nrj)
	endif
!
!	Mise a jour de l'\E9tat fondamental
!	---------------------------------
	do i=1,n
!
		if(mf) then
			sp(i)  = sqrt(sum(jmy(i,:)*jmy(i,:)))
			ss     = jmy(i,:)/sqrt(sum(jmy(i,:)*jmy(i,:)))
			theta(i) = acos(ss(3))
			y = 1.0e-4
			x = abs(theta(i))
			if((abs(x).le.y).or.(abs(x-pi).lt.y)) then
				phi(i) = 0.0
			else
				x = ss(1)/sin(theta(i))
				if(1.-abs(x).lt.y) then
					if(abs(x-1.).lt.y) phi(i) = 0.0
					if(abs(x+1.).lt.y) phi(i) = pi/d2r
				else
					if(ss(2).lt.0.0) phi(i) = (2*pi-acos(ss(1)/sin(theta(i))))/d2r
					if(ss(2).ge.0.0) phi(i) = (     acos(ss(1)/sin(theta(i))))/d2r
				endif
			endif
			theta(i) = theta(i)/d2r
		else
			sp(i)  = wj(i)
		endif
!
!		axes li\E9s \E0 la direction des spins d\E9finis par phi et theta
!		vecteur e3 de la base
		matl(i,3,1) = cos(phi(i)*d2r)*sin(theta(i)*d2r)
		matl(i,3,2) = sin(phi(i)*d2r)*sin(theta(i)*d2r)
		matl(i,3,3) = cos(theta(i)*d2r)
!		vecteur e1 de la base
!		on choisit e1 dans le plan perp a e3
!		tq OM.e3 = 0
!		on impose que e1 soit norme et de cote z=0
		matl(i,2,1) = +sin(phi(i)*d2r)
		matl(i,2,2) = -cos(phi(i)*d2r)
		matl(i,2,3) = 0.0
		call vct(matl(i,2,:),matl(i,3,:),matl(i,1,:))
		call inverse(d,matl(i,:,:),matr(i,:,:))
		zdr(i,:) = matr(i,:,1)
		zdi(i,:) = matr(i,:,2)
		eta(i,:) = matr(i,:,3)
!
	enddo
!
!	-----------------------------------------------------------------------
!
	call banner("DIRECTIONS DES SPINS")
	write(*,"(2a4,3a8,3a12)") "i","l","h","k","l","ETAX","ETAY","ETAZ"
	allocate(eta0(n,ntt,d),zdr0(n,ntt,d),zdi0(n,ntt,d))
	do i=1,n
		do l=1,ntt
			x  = sum(2.0*pi*kprop*listvec(l,:))
			call rodrigue(x,uprop(i,:),zdr(i,:),zdr0(i,l,:))
			call rodrigue(x,uprop(i,:),zdi(i,:),zdi0(i,l,:))
			call rodrigue(x,uprop(i,:),eta(i,:),eta0(i,l,:))
			write(*,"(2i4,3f8.3,3f12.5,f12.5)") i,l,(listvec(l,k),k=1,d),(eta0(i,l,k),k=1,d),x/pi*180.
		enddo
	enddo
!
!	trace des spin dans la maille format fullprof
!
	if(ffst.ne."") then
	call banner("SAUVEGARDE FORMAT FST")
	lu=1
	open(lu,file=ffst)
	write(lu,"(a15)") "BKG 0.9 0.9 0.9 1.0"
	write(lu,"(a5,6f12.5,a29)") "CELL ",ax,ay,az,alfa*180/pi,beta*180/pi,gama*180/pi,&
	" COLOR 0.5 0.5 0.5 1 MULTIPLE"
	write(lu,"(a4,6f12.5)") "BOX ",-0.05,1.05,-0.05,1.05,-0.05,1.05
	write(lu,"(a1)")  "{"
	write(lu,"(a9)")  "LATTICE P"
	write(lu,"(a13)") "K 0.0 0.0 0.0"
	write(lu,"(a10)") "SYMM x,y,z"
	write(lu,"(a14)") "MSYM u,v,w,0.0"
	do i=1,n
		write(lu,"(a9,3f12.5,a10)") "MATOM A A",(pos(i,j),j=1,d)," SCALE 0.5"
		write(lu,"(a7,7f12.5,a6,3f12.5,a3)") "SKP 1 1",(jmy(i,j),j=1,d),0.,0.,0.,0.," COLOR",0.,i*1./n,0.5," 1."
	enddo
	write(lu,"(a1)") "}"
	close(lu)
	endif
!
!	-----------------------------------------------------------------------
!	defintion des valeurs du vecteur d'onde
!	Liste des valeurs de Q ou on veut afficher les sections efficaces
!
	if(tof) then
!		liste 	
		allocate(listeq(np,d))
!		liste globale (suite de fibonacci)
		call suite_fibo(nfibo,jf,jf1)
		nq = np*(jf+1)
		allocate(tq(nq,d),tqc(nq,d),tqr(nq,d),vq(nq),su(nq))
		listeq = 0.0
		tq     = 0.0
!
		allocate(vphi(jf+1),vz(jf+1))
		call point_fibo(jf,jf1,vphi,vz)
		do p = 1,np
			x = qmin+(qmax-qmin)*p/np
			listeq(p,3) = x
			do j = 1,jf+1
				l = (p-1)*(jf+1)+j
				tqc(l,1) = x*sqrt(1.0-vz(j)*vz(j))*cos(vphi(j))
				tqc(l,2) = x*sqrt(1.0-vz(j)*vz(j))*sin(vphi(j))
				tqc(l,3) = x*vz(j)					
				su(l)    = 4.0*pi/jf
				vq(l)    = p
			enddo
		enddo
		call inverse(d,ee,mat)
		do p = 1,np
			listeq(p,:) = matmul(mat,listeq(p,:))
		enddo
		do p = 1,nq
			tq(p,:) = matmul(mat,tqc(p,:))
		enddo
		tqr = tqc
	endif
!
	if(coupe) then
		np  = np1*np2
		nq  = np
		allocate(listeq(np,d))
		allocate(tq(nq,d),tqc(nq,d),tqr(nq,d),vq(nq),su(nq))
		listeq = 0.0
		tq     = 0.0
		
		do i=1,np1
		do j=1,np2
			p = (i-1)*np2+j
			listeq(p,:) = q0 + (i-1)*dq1+(j-1)*dq2
			tq(p,:)     = listeq(p,:)
			tqc(p,:)    = matmul(ee,tq(p,:))
			su(p)       = 1.0
			vq(p)       = p
		enddo
		enddo
		tqr = tqc
	endif
!	
	if(simpleq) then
	
		nq  = np
		allocate(listeq(np,d))
		allocate(tq(nq,d),tqc(nq,d),tqr(nq,d),vq(nq),su(nq))
		listeq = 0.0
		tq     = 0.0
!	
		do p=1,np
			listeq(p,:) = q0 + (p-1)*dq
			tq(p,:)     = listeq(p,:)
			tqc(p,:)    = matmul(ee,tq(p,:))
			su(p)       = 1.0
			vq(p)       = p
		enddo
		tqr = tqc
	endif
!
	if(reduc) then
		np = nbk**3
		allocate(listeq(np,d))
		allocate(tq(nq,d),tqc(nq,d),tqr(nq,d),vq(nq),su(nq))
		listeq = 0.0
		tq     = 0.0
		
		l = 0
		do i = 1,nbk
		do j = 1,nbk
		do k = 1,nbk
			q(1) = i*1.0/nbk
			q(2) = j*1.0/nbk
			q(3) = k*1.0/nbk
			l = l+1
			listeq(l,:) = q
			tq(l,:)     = listeq(p,:)
			tqc(l,:)    = matmul(ee,tq(l,:))
			su(l  )     = 1.0	
			vq(l  )     = l	
		enddo
		enddo
		enddo
		tqr = tqc
	endif
!	
	if(coup1d) then
!	
		nk = np1*np2*np3
		allocate(tq2(nq + nq*nk,d),vq2(nq + nq*nk),su2(nq + nq*nk))
		tq2(1:nq,:) = tq(1:nq,:)
		vq2(1:nq)   = vq(1:nq)
		su2(1:nq)   = su(1:nq)
!	
		do p=1,nq
			l = nq + (p-1)*nk
			do i = 1,np1
			do j = 1,np2
			do k = 1,np3
				l        = l +1
				tq2(l,:) = tq(p,:)+i*dq1+j*dq2+k*dq3
				vq2(l  ) = vq(p)
				su2(l  ) = su(p)
			enddo
			enddo
			enddo
		enddo
!		
		deallocate(tq,tqc,tqr,vq,su)
		nq = nq + nq*nk
		allocate(tq(nq,d),tqc(nq,d),tqr(nq,d),vq(nq),su(nq))
		tq = tq2
		vq = vq2
		su = su2
		do p=1,nq
			tqc(p,:) = matmul(ee,tq(p,:))
		enddo
		tqr = tqc
!
		deallocate(tq2,vq2,su2)		
	endif
!
!	dans le cas ou on veut calculer le continuum
	if(cnt) then
!		
		nk = nbk**3
		allocate(tq2(nq+nk+nq*nk,d),vq2(nq+nk+nq*nk),su2(nq+nk+nq*nk),vk2(nq+nk+nq*nk))
		allocate(tqr2(nq+nk+nq*nk,d))
		tq2(1:nq,:)  = tq(1:nq,:)
		tqr2(1:nq,:) = tqr(1:nq,:)
		vq2(1:nq)    = vq(1:nq)
		su2(1:nq)    = su(1:nq)
		vk2          = 0
!
!		valeurs pour k et Q+k
		l = 0
!		call init_random_seed()
		do i = 1,nbk
		do j = 1,nbk
		do k = 1,nbk
!			call random_number(q)
			q(1) = i*1.0/nbk
			q(2) = j*1.0/nbk
			q(3) = k*1.0/nbk
			l = l+1
			tq2(nq+l,:) = q
			vq2(nq+l  ) = 0
			su2(nq+l  ) = 0.0
			do p=1,nq
				tq2(nq+p*nk+l,:)  = q+tq(p,:)
				tqr2(nq+p*nk+l,:) = tqr(p,:)
				vq2(nq+p*nk+l  )  = vq(p)
				su2(nq+p*nk+l  )  = 0.0
				vk2(nq+p*nk+l  )  = nq+l
			enddo
		enddo
		enddo
		enddo
!		
		deallocate(tq,tqc,tqr,vq,su)
		nq = nq+nk+nq*nk
		allocate(tq(nq,d),tqc(nq,d),tqr(nq,d),vq(nq),su(nq),vk(nq))
		tq  = tq2
		tqr = tqr2
		vq  = vq2
		su  = su2
		vk  = vk2
		do p=1,nq
			tqc(p,:) = matmul(ee,tq(p,:))
		enddo
		
		deallocate(tq2,tqr2,vq2,su2,vk2)
	endif
!
!	dans le cas des twins, on ajoute les Q correspondants
!	listeq n'est pas modifi\E9e
	if(twin) then
		allocate(tq2(nq + nq*ntw,d),tqr2(nq + nq*ntw,d),vq2(nq + nq*ntw),su2(nq + nq*ntw))
		tq2(1:nq,:)  = tq(1:nq,:)
		tqr2(1:nq,:) = tqr(1:nq,:)
		vq2(1:nq)    = vq(1:nq)
		su2(1:nq)    = su(1:nq)
!		
		l = nq
		do i = 1,ntw
		do p = 1,nq
			l = nq+(i-1)*nq
			qc  = matmul(ee,tq(p,:))
			ss  = matmul(e,rtw(i,:))
			ss  = ss/sqrt(sum(ss*ss))
			call rodrigue(atw(i)*pi/180,ss,qc,qc2)
			call inverse(d,ee,mat)
			tq2(l+p,:)  = matmul(mat,qc2)
			tqr2(l+p,:) = tqr(p,:)
			vq2(l+p  )  = vq(p)
			su2(l+p  )  = su(p)
		enddo		
		enddo
!	
		deallocate(tq,tqc,tqr,vq,su)
		nq = nq + nq*ntw
		allocate(tq(nq,d),tqc(nq,d),tqr(nq,d),vq(nq),su(nq))
		tq  = tq2
		tqr = tqr2
		vq  = vq2
		su  = su2
		do p=1,nq
			tqc(p,:) = matmul(ee,tq(p,:))
		enddo
		deallocate(tq2,tqr2,vq2,su2)
	endif
!
!	dans le cas incommensurable, on triple le nombre de calculs
!	pour k+kprop, k-kprop
!	ce bloc doit rester en fin de d\E9claration des Q
	if(incom) then
		allocate(tq2(nq+2*nq,d),tqr2(nq+2*nq,d),vq2(nq+2*nq),su2(nq+2*nq))
		tq2(1:nq,:)  = tq(1:nq,:)
		tqr2(1:nq,:) = tqr(1:nq,:)
		vq2(1:nq)    = vq(1:nq)
		su2(1:nq)    = su(1:nq)
		do p=1,nq
			tq2 (  nq+p,:) = tq(p,:)+kprop
			tqr2(  nq+p,:) = tqr(p,:)
			vq2 (  nq+p  ) = vq(p)
			su2 (  nq+p  ) = su(p)
			tq2 (2*nq+p,:) = tq(p,:)-kprop
			tqr2(2*nq+p,:) = tqr(p,:)
			vq2 (2*nq+p  ) = vq(p)
			su2 (2*nq+p  ) = su(p)
		enddo
!		
		deallocate(tq,tqc,tqr,vq,su)
		nq2= nq
		nq = 3*nq
		allocate(tq(nq,d),tqc(nq,d),tqr(nq,d),vq(nq),su(nq))
		tq  = tq2
		tqr = tqr2
		vq  = vq2
		su  = su2
		do p=1,nq
			tqc(p,:) = matmul(ee,tq(p,:))
		enddo
		deallocate(tq2,tqr2,vq2,su2)
	endif
!
!	tableau des associations
!
	allocate(tab(np,nq),ntab(np))
	do p=1,np
		k=0
		do ii=1,nq
			if(vq(ii).eq.p) then
				k=k+1
				tab(p,k) = ii
			endif
		enddo
		ntab(p) = k
	enddo
!
!	-----------------------------------------------------------------------
!	grille en energie pour affichage final
!
	if(coupe) nw = 1
	allocate(w(nw))
	do i = 1,nw
		w(i) = wmin+(i-1)*wmax/(nw-1)
	enddo
	if(coupe) w(1)=en0
!
!	-----------------------------------------------------------------------
!	Affichage des donnees du calcul
!
	lu = 1
	open(lu,file=fres)
	call banner("DONNEES DU CALCUL")
!
	if(tof) then
	write(*,"(a12,i12)")    "NFIBO" ,nfibo
	write(*,"(a12,i12)")    "FN"    ,jf
	write(*,"(a12,i12)")    "NP"    ,nq
	write(*,"(a12,f12.5)")  "QMIN"  ,qmin
	write(*,"(a12,f12.5)")  "QMAX"  ,qmax
	endif
	if(coupe) then
	write(*,"(a12,f12.5)")  "EN0"   ,en0
	write(*,"(a12,3f12.5)") "Q0"    ,(q0(i) ,i=1,d)
	write(*,"(a12,3f12.5)") "DQ1"   ,(dq1(i),i=1,d)
	write(*,"(a12,3f12.5)") "DQ2"   ,(dq2(i),i=1,d)
	write(*,"(a12,i12)")    "NP1"   ,np1
	write(*,"(a12,i12)")    "NP2"   ,np2
	write(*,"(a12,i12)")    "NP1"   ,np1
	write(*,"(a12,i12)")    "NP"    ,np
	endif
	if(coup1d) then
	write(*,"(a12,i12)")    "NP1"   ,np1
	write(*,"(a12,i12)")    "NP2"   ,np2
	write(*,"(a12,i12)")    "NP3"   ,np3
	write(*,"(a12,3f12.5)") "Q0"    ,(q0(i) ,i=1,d)
	write(*,"(a12,3f12.5)") "DQ1"   ,(dq1(i),i=1,d)
	write(*,"(a12,3f12.5)") "DQ2"   ,(dq2(i),i=1,d)
	write(*,"(a12,3f12.5)") "DQ3"   ,(dq3(i),i=1,d)
	endif
	if(simpleq) then
	write(*,"(a12,3f12.5)") "Q0"    ,(q0(i),i=1,d)
	write(*,"(a12,3f12.5)") "DQ"    ,(dq(i),i=1,d)
	endif
!
	write(*,"(a12,f12.5)")  "AX"    ,ax
	write(*,"(a12,f12.5)")  "AY"    ,ay
	write(*,"(a12,f12.5)")  "AZ"    ,az
	write(*,"(a12,f12.5)")  "ALFA"  ,alfa/d2r
	write(*,"(a12,f12.5)")  "BETA"  ,beta/d2r
	write(*,"(a12,f12.5)")  "GAMA"  ,gama/d2r
	write(*,"(a12,f12.5)")  "SIG"   ,sig
	write(*,"(a12,f12.5)")  "fA"    ,fga
	write(*,"(a12,f12.5)")  "fa"    ,fa
	write(*,"(a12,f12.5)")  "fB"    ,fgb
	write(*,"(a12,f12.5)")  "fb"    ,fb
	write(*,"(a12,f12.5)")  "fC"    ,fgc
	write(*,"(a12,f12.5)")  "fc"    ,fc
	write(*,"(a12,f12.5)")  "fD"    ,fgd
!
	write(*,"(a12,f12.5)")  "EZX"   ,ez(1)
	write(*,"(a12,f12.5)")  "EZY"   ,ez(2)
	write(*,"(a12,f12.5)")  "EZZ"   ,ez(3)
!
	write(*,"(a12,i12  )")  "NITER" ,niter
	write(*,"(a12,f12.5)")  "TEMP"  ,temp
	write(*,"(a12,i12  )")  "CALC"  ,calc
!
	write(*,"(a12,f12.5)")  "HX"   ,h0(1)
	write(*,"(a12,f12.5)")  "HY"   ,h0(2)
	write(*,"(a12,f12.5)")  "HZ"   ,h0(3)
	write(*,"(a12,f12.5)")  "H"    ,sqrt(sum(h0*h0))
	write(*,"(a12,f12.5)")  "GX"   ,gl(1,1)
	write(*,"(a12,f12.5)")  "GY"   ,gl(2,2)
	write(*,"(a12,f12.5)")  "GZ"   ,gl(3,3)
	write(*,"(a12,f12.5)")  "GX2"  ,gl2(1,1)
	write(*,"(a12,f12.5)")  "GY2"  ,gl2(2,2)
	write(*,"(a12,f12.5)")  "GZ2"  ,gl2(3,3)
!
	if(partiel) write(*,"(a12)")  "PARTIEL"
	if(reduc)   write(*,"(a12)")  "REDUC"
	if(twin) then
		write(*,"(a12)")  "TWIN"
		do i=1,ntw
			write(*,"(i4,4f12.5)") i,rtw(i,1),rtw(i,2),rtw(i,3),atw(i)
		enddo
	endif
!
!	-----------------------------------------------------------------------
!
	call banner("LISTE DES ATOMES MAGNETIQUES")
!
	write(*,"(a4,16a8)") "ION","NOM","rlu","rlu","rlu",&
	"X","Y","Z","MASSE","PHI","THETA","ETAX","ETAY","ETAZ","CX","CY","CZ"
	do i=1,n
		rc = matmul(e,pos(i,:))
		write(*,"(i4,a8,15f8.3)") i,nom(i),(pos(i,j),j=1,d),(rc(j),j=1,3),&
		ma(i),phi(i),theta(i),(eta(i,j),j=1,d),(xcef(i,j),j=1,d)
	enddo
	write(*,"(a4,15a8)") "ION","NOM","WJ","GJ","B20","B22",&
	"B40","B42","B43","B44","B60","B62","B63","B64","B66","DD"
	do i=1,n
		write(*,"(i4,a8,14f8.3)") i,nom(i),wj(i),gj(i),b20(i),b22(i),&
		b40(i),b42(i),b43(i),b44(i),b60(i),b62(i),b63(i),b64(i),b66(i),dis0(i)
	enddo
	write(*,"(a4,15a8)") "ION","NOM","WJ","GJ","B20","B22","B2M2","B21","B2M1"
	do i=1,n
		write(*,"(i4,a8,14f8.3)") i,nom(i),wj(i),gj(i),b20(i),b22(i),&
		b22(i),b2m2(i),b21(i),b2m1(i)
	enddo
	write(*,"(a4,6a8)") "ION","NOM","JMX","JMY","JMZ","JM","SP"
	do i=1,n
		write(*,"(i4,a8,5f8.3)") i,nom(i),(jmy(i,j),j=1,d),sqrt(sum(jmy(i,:)*jmy(i,:))),sp(i)
	enddo
!
!
	write(*,*)
	call banner("MATRICES DE PASSAGE")
	do i=1,n
		write(*,"(a5,i4)") "Site ",i
		write(*,"(a36,5x,a36)") "locales -> cartesiennes","cartesiennes -> locales"
		do l1=1,d
			write(*,"(3f12.5,5x,3f12.5)") (matr(i,l1,l2),l2=1,d),(matl(i,l1,l2),l2=1,d)
		enddo
	enddo

	write(*,*)
	call banner("MATRICES DE PASSAGE CEF")
	do i=1,n
		write(*,"(a5,i4)") "Site ",i
		write(*,"(a36,5x,a36)") "locales -> cartesiennes","cartesiennes -> locales"
		do l1=1,d
			write(*,"(3f12.5,5x,3f12.5)") (matrx(i,l1,l2),l2=1,d),(matlx(i,l1,l2),l2=1,d)
		enddo
	enddo
	if(nato.ne.0) then
		write(*,*)
		call banner("LISTE DES ATOMES PONTANTS")
		write(*,"(a4,3a12)") "ION","rlu","rlu","rlu"
		do i=1,nato
			write(*,"(a4,3f12.5)")  lato(i)%nam,(lato(i)%pos(j),j=1,d)
		enddo
	endif
	if(incom) then
		write(*,*)
		call banner("STRUCTURE INCOMMENSURABLE")
		write(*,"(a8,3f8.3)") "K prop ",(kprop(j),j=1,d)
		do i=1,n
			write(*,"(i3,3f8.3)") i,(uprop(i,j),j=1,d)		
		enddo
	else
		write(*,*)
		call banner("STRUCTURE COMMENSURABLE")
		write(*,"(a8,3f8.3)") "K prop ",(kprop(j),j=1,d)	
	endif
!
	write(*,*)
	call banner("DISTANCES INTER-ATOMIQUES")
	do i=1,n
		allocate(u1(n))
		do j=1,n
			rc    = pos(j,:) - pos(i,:)
			rc    = matmul(e,rc)
			u1(j) = sqrt(sum(rc*rc))
		enddo
		write(*,"(i4,50f8.3)") i,(u1(j),j=1,n)
		deallocate(u1)
	enddo
!			
	write(*,*) "Calcul de z, z* et eta"
	do i=1,n
!		write(*,"(i4,3a8)") i,"zr","zi","eta"
		do j=1,d
!			write(*,"(4x,3f8.3)") zdr(i,j),zdi(i,j),eta(i,j)
		enddo
	enddo
!
	write(*,*) "Calcul de z z"
!	write(*,"(4x,16i12)") (j,j=1,n)
	do i=1,n
	allocate(u1(n),u2(n))
	do j=1,n
		u1(j) = sum( zdr(i,:)*zdr(j,:)-zdi(i,:)*zdi(j,:))
		u2(j) = sum( zdr(i,:)*zdi(j,:)+zdi(i,:)*zdr(j,:))
	enddo
!	write(*,"(i4,16(f6.2,f6.2))") i,(u1(j),u2(j),j=1,n)
	deallocate(u1,u2)
	enddo
!
	write(*,*) "Calcul de z z*"
!	write(*,"(4x,16i12)") (j,j=1,n)
	do i=1,n
	allocate(u1(n),u2(n))
	do j=1,n
		u1(j) = sum( zdr(i,:)*zdr(j,:)+zdi(i,:)*zdi(j,:))
		u2(j) = sum(-zdr(i,:)*zdi(j,:)+zdi(i,:)*zdr(j,:))
	enddo
!	write(*,"(i4,16(f6.2,f6.2))") i,(u1(j),u2(j),j=1,n)
	deallocate(u1,u2)
	enddo
!
	write(*,*) "Calcul de eta eta"
!	write(*,"(4x,16i6)") (j,j=1,n)
	do i=1,n
	allocate(u1(n))
	do j=1,n
		u1(j) = sum(eta(i,:)*eta(j,:))
	enddo
!	write(*,"(i4,16f6.2)") i,(u1(j),j=1,n)
	deallocate(u1)
	enddo
!
!	write(*,*) "Calcul de eta x eta"
!	write(*,"(4x,16i6)") (j,j=1,n)
!	do i=1,n
!		do j=1,n
!			call vct(eta(i,:),eta(j,:),vdmr)
!			write(*,"(2i4,16f6.2)") i,j,(vdmr(l),l=1,d)
!		enddo
!	enddo
!
!	write(*,*) "Calcul de z x z"
!	write(*,"(4x,16i6)") (j,j=1,n)
!	do i=1,n
!		do j=1,n
!			call vctx(zdr(i,:),zdi(i,:),zdr(j,:),zdi(j,:),vdmr,vdmi)
!			write(*,"(2i4,16f6.2)") i,j,&
!			(vdmr(l),l=1,d),(vdmi(l),l=1,d)
!		enddo
!	enddo
!	write(*,*) "Calcul de z x zbar"
!	write(*,"(4x,16i6)") (j,j=1,n)
!	do i=1,n
!		do j=1,n
!			call vctx(zdr(i,:),zdi(i,:),zdr(j,:),-zdi(j,:),vdmr,vdmi)
!			write(*,"(2i4,16f6.2)") i,j,&
!			(vdmr(l),l=1,d),(vdmi(l),l=1,d)
!		enddo
!	enddo
!
!	write(*,*) "Calcul de z x eta"
!	write(*,"(4x,16i6)") (j,j=1,n)
!	do i=1,n
!		do j=1,n
!			call vct(zdr(i,:),eta(j,:),vdmr)
!			call vct(zdi(i,:),eta(j,:),vdmi)
!			write(*,"(2i4,16f6.2)") i,j,(vdmr(l),vdmi(l),l=1,d)
!		enddo
!	enddo
!
	call banner("INTERACTIONS MAGNETIQUES")
	do l=1,nz
	write(*,"(a17,i2,a8)") "Termes d echange ",l," voisins"
!	write(*,"(2a4,4a8)") "I1","I2","JXX","JYY","JZZ","D"
!	do i=1,n
!	do j=i,n
!		x = 0.0
!		do k=1,d
!			x=x+jech(i,j,l,k,k)
!		enddo
!		if(x.ne.0.0) write(*,"(2i4,4f8.3)") i,j,(jech(i,j,l,k,k),k=1,d),dist(l,i,j)
!	enddo
!	enddo
	enddo
	call banner("ANISOTROPIE H")
	write(*,"(a4,a8)") "I1","H"
	do i=1,n
		write(*,"(i4,f8.3)") i,h(i)
	enddo
	call banner("ANISOTROPIE D")
	write(*,"(a4,4a8)") "I1","DV","DXX","DYY","DZZ"
	do i=1,n
		write(*,"(i4,4f8.3)") i,dval(i),(dech(i,k),k=1,d)
	enddo

	call banner("ANISOTROPIE Bnm")
	do i=1,n
		write(*,"(a4,i4)") "I1",i
		do k=1,d
			write(*,"(3f8.3)") (dani(i,k,l),l=1,d)
		enddo
	enddo
	if((gdm.ne.0.0).or.(gksea.ne.0.0).or.(gdiag.ne.0.0)) then
	call banner("DZYALOSHINSKII")
	write(*,"(a8,f8.3)") "GDM   =",gdm
	write(*,"(a8,f8.3)") "GKSEA =",gksea
	write(*,"(a8,f8.3)") "GDIAG =",gdiag
	write(*,"(2a4,9a8)") "i","j","DMx","DMy","DMz","Ux","Uy","Uz","Ptx","Pty","Ptz"
	do i=1,n
	do j=1,n
		do l=1,ntt
			x = sum(dm(i,j,l,:)*dm(i,j,l,:))
			if(x.ne.0.0) then
			write(*,"(2i4,9f8.3)") i,j,(dm(i,j,l,k),k=1,d),(listvec(l,k),k=1,d),(pt(i,j,l,k),k=1,d)
			endif
		enddo
	enddo
	enddo
	endif
!
	if(bbq) then
	call banner("BIQUADRATIQUE")
	write(*,"(4a4,3a8)") "i","j","k","l","Value","dij","dkl"
	do l1=1,n
	do l2=1,n
	do l3=1,n
	do l4=1,n
		if(biq(l1,l2,l3,l4).ne.0.0) then
		write(*,"(4i4,3f8.3)") l1,l2,l3,l4,biq(l1,l2,l3,l4),&
		dbiq12(l1,l2),dbiq13(l3,l4)
		endif
	enddo
	enddo
	enddo
	enddo
	endif
!
	call banner("CONSTANTES DE FORCE")
	do l=1,nz
	write(*,"(a17,i2,a8)") "Ctes de force",l," voisins"
!	write(*,"(2a4,4a8)") "I1","I2","GXX","GYY","GZZ","D"
!	do i=1,n
!	do j=i,n
!		x = 0.0
!		do k=1,d
!			x=x+ctef(i,j,l,k,k)
!		enddo
!		if(x.ne.0.0) write(*,"(2i4,4f8.3)") i,j,(ctef(i,j,l,k,k),k=1,d),disf(l,i,j)
!	enddo
!	enddo
	enddo
!	
!	-----------------------------------------------------------------------
!	Affichages valeurs de Q
!
	call banner("LISTE DES VALEURS DE Q")
	call banner2("Liste reduite")
	do p=1,np
		write(*,"(i8,3f8.3,i8)") p,(listeq(p,i),i=1,d),ntab(p)
		do k=1,ntab(p)
		write(*,"(3i8,12f8.3)") p,k,tab(p,k),(tq(tab(p,k),i),i=1,d)
		enddo
	enddo
	call banner2("Liste complete")
	do p=1,nq
		write(*,"(i8,3f8.3,i8)") p,(tq(p,i),i=1,d),vq(p)
	enddo
!
!	-----------------------------------------------------------------------
!	calcul du facteur de structure magnetique pour la structure trouvee
!
	call banner("FACTEUR DE STRUCTURE MAGNETIQUE")
	if(incom.eqv..false.) then
	write(*,"(a4,3a6,5a13)") "No","h","k","l","MyMy*","MzMz*","Re(MyMz*)","Im(MyMz*)","MyMy*+MzMz*"
	
	ntt2 = (2*nnx+1)*(2*nny+1)*(2*nnz+1)
	allocate(listvec2(ntt2,d))
	l = 0
	do l1 = -nnx,nnx,1
	do l2 = -nny,nny,1
	do l3 = -nnz,nnz,1
		l = l+1
		listvec2(l,1) = l1
		listvec2(l,2) = l2
		listvec2(l,3) = l3
	enddo
	enddo
	enddo
!	
	do l=1,ntt2
		qc = matmul(ee,listvec2(l,:))
!
!		facteur de forme (attention pas au carre)
		fmag = 1.0
		if(fform) then
			s    = sqrt(sum(qc*qc))/(4.0*pi)
			fmag = fga*exp(-fa*s*s)+fgb*exp(-fb*s*s)+fgc*exp(-fc*s*s)+fgd
		endif
!
!		facteur de structure
		allocate(u1(d),u2(d))
		call fm(n,listvec2(l,:),qc,pos,jmy,fmag,u1,u2)
!
!		direction de Y et de Z
		ex = qc/sqrt(sum(qc*qc))
		call vct(ez,ex,ey)
!
		u1=u1-ex*(sum(ex*u1))
		u2=u2-ex*(sum(ex*u2))
		xr = sum(u1*ey)
		xi = sum(u2*ey)
		y  = xr*xr+xi*xi
		yr = sum(u1*ez)
		yi = sum(u2*ez)
		z  = yr*yr+yi*yi		
!
		wr = (+xr*yr+xi*yi)*2.0
		wi = (-xr*yi+xi*yr)*2.0
		write(*,"(i4,3f6.1,5f13.5)") l,(listvec2(l,k),k=1,d),y,z,wr,wi,y+z
		deallocate(u1,u2)
	enddo
	deallocate(listvec2)
	endif
!
!	-----------------------------------------------------------------------
!	on libere partiellement la memoire
!
	deallocate(phi,theta)
!
!	-----------------------------------------------------------------------
!	Calcul de la matrice de omega
!
	n0 = n
	if(couple) n = n0+d*n0
	n2 = 2*n
	write(*,"(a20,i4)") "Magnetic Matrix Dimension ",n0
	write(*,"(a20,i4)") "Total    Matrix Dimension ",n
!
	call banner("MATRICE OMEGA")
!
!	contribution Heisenberg
!
	call banner2("Contribution Heisenberg")
	write(*,"(2a4,3a7,4a12,a8)") "i","j","dx","dy","dz","Jx","Jy","Jz","Dist","Modele"
!
	allocate(omega(n0))
	detail =.true.
	qc     = 0.0
	do i=1,n0
		x = 0.0
		do j=1,n0
		do l=1,ntt
!
			call mag(i,j,listvec(l,:),pos,matrx(i,:,3),matrx(j,:,3),e,nz,modele(i,j),&
			dani,dist(:,i,j),jech(i,j,:,:,:),qc,detail,jm)
			call cpl(i,j,listvec(l,:),pos,e,nz,modele(i,j),disf(:,i,j),ctef(i,j,:,:,:),&
			qc,uo,uu)
!
			etao = eta0(j,l,:)
!
			z =   sp(j)   *sum(eta(i,:)*matmul(jm,etao))
!			if(z.ne.0.0) write(*,*) j,l,etao,z
			y = g*sp(j)*uo*sum(eta(i,:)*          etao )
			x = x + z + y
		enddo
		enddo
		ss = gj(i)*mub*j2unit*matmul(gl,matmul(matlx(i,:,:),h0))

		omega(i) = x-h(i)*sp(i)	

		write(*,"(i4,f12.5)") i,omega(i)
		write(*,*)	
	enddo
!
!	contribution DM
!
	if(gdm.ne.0.0) then
	call banner2("Contribution DM")
	write(*,"(2a4,6a7)") "i","j","dx","dy","dz","DMx","DMy","DMz"
	do i=1,n0
		x = 0.0
		do j=1,n0
		do l=1,ntt
!
			etao = eta0(j,l,:)
!
			call vct(eta(i,:),etao,vdmr)
			x = x + sp(j)*sum(dm(i,j,l,:)*vdmr)
!
			if(sum(dm(i,j,l,:)**2).ne.0.0) then
			write(*,"(2i4,6f7.3)") i,j,(listvec(l,k),k=1,d),dm(i,j,l,:)
			endif
		enddo
		enddo
		omega(i) = omega(i) + x*gdm
		write(*,"(i4,f12.5)") i,omega(i)
	enddo
	endif
!
!	contribution KSEA
!
	if(gksea.ne.0.0) then
	call banner2("Contribution KSEA")
	write(*,"(2a4,6a7)") "i","j","dx","dy","dz","DMx","DMy","DMz"
	do i=1,n0
		x = 0.0
		do j=1,n0
		do l=1,ntt
!
			etao = eta0(j,l,:)
!
			jm = 0.0
			do l1=1,d
			do l2=1,d
				jm(l1,l2) =  gksea*dm(i,j,l,l1)*dm(i,j,l,l2)
			enddo
			enddo
			do l1=1,d
				jm(l1,l1) = 0.0
			enddo
			x = x + sp(j)*sum(eta(i,:)*matmul(jm,etao))
!
			if(sum(dm(i,j,l,:)**2).ne.0.0) then
			write(*,"(2i4,6f7.3)") i,j,(listvec(l,k),k=1,d),(jm(1,l2),l2=1,d)
			write(*,"(29x,3f7.3)") (jm(2,l2),l2=1,d)
			write(*,"(29x,3f7.3)") (jm(3,l2),l2=1,d)

			endif
		enddo
		enddo
		omega(i) = omega(i) + x
		write(*,"(i4,f12.5)") i,omega(i)
	enddo
	endif
!
!	contribution GDIAG
!
	if(gdiag.ne.0.0) then
	call banner2("Contribution GDIAG")
	write(*,"(2a4,6a7)") "i","j","dx","dy","dz","DMx","DMy","DMz"
	do i=1,n0
		x = 0.0
		do j=1,n0
		do l=1,ntt
!
			etao = eta0(j,l,:)
!
			jm = 0.0
			do l1=1,d
				jm(l1,l1) = 1.0
			enddo
			jm = gdiag*sum(dm(i,j,l,:)**2)*jm
			x = x + sp(j)*sum(eta(i,:)*matmul(jm,etao))
!
			if(sum(dm(i,j,l,:)**2).ne.0.0) then
			write(*,"(2i4,6f7.3)") i,j,(listvec(l,k),k=1,d),(jm(1,l2),l2=1,d)
			write(*,"(29x,3f7.3)") (jm(2,l2),l2=1,d)
			write(*,"(29x,3f7.3)") (jm(3,l2),l2=1,d)

			endif
		enddo
		enddo
		omega(i) = omega(i) + x
		write(*,"(i4,f12.5)") i,omega(i)
	enddo
	endif
!
!	contribution biquadratique
!
	if(bbq) then
	call banner2("Contribution biquadratique")
	write(*,"(2a4,3a7,a21)") "i","j","dx","dy","dz","J_effectif"
	do i=1,n0
		x = 0.0
		do j=1,n0
		do l=1,ntt
!
			etao = eta0(j,l,:)
!			
			jm = 0.0
			allocate(u1(3),u2(3))
			rc = pos(i,:)+listvec(l,:)-pos(j,:)
			rc = matmul(e,rc)
			if(sqrt(sum(rc*rc)).le.dbiq12(i,j)) then
			do k3=1,n
			do l3=1,ntt
				rc = pos(k3,:)-pos(i,:)+listvec(l3,:)
				rc = matmul(e,rc)
				if(sqrt(sum(rc*rc)).le.dbiq13(i,k3)) then
				do k4=1,n
				do l4=1,ntt
					rc = pos(k3,:)-pos(k4,:)+listvec(l3,:)-listvec(l4,:)
					rc = matmul(e,rc)
					if(sqrt(sum(rc*rc)).le.dbiq12(k3,k4)) then					
						u1 = sp(k3)*eta0(k3,l3,:)
						u2 = sp(k4)*eta0(k4,l4,:)
						x = 2.0*biq(i,j,k3,k4)*sum(u1*u2)
						do l1=1,d
						jm(l1,l1) = jm(l1,l1)+x
						enddo
					endif
				enddo
				enddo
				endif
			enddo
			enddo
			endif
			if(sqrt(sum(rc*rc)).le.dbiq13(i,j)) then
			do k3=1,n
			do l3=1,ntt
				rc = pos(k3,:)-pos(i,:)+listvec(l3,:)
				rc = matmul(e,rc)
				if(sqrt(sum(rc*rc)).le.dbiq12(i,k3)) then
				do k4=1,n
				do l4=1,ntt
					rc = pos(k4,:)-pos(j,:)+listvec(l4,:)-listvec(l,:)
					rc = matmul(e,rc)
					if(sqrt(sum(rc*rc)).le.dbiq12(j,k4)) then
						u1 = sp(k3)*eta0(k3,l3,:)
						u2 = sp(k4)*eta0(k4,l4,:)
						do l1=1,d
						do l2=1,d
						jm(l1,l2) = jm(l1,l2)+4.0*biq(i,k3,j,k4)*u1(l1)*u2(l2)
						enddo
						enddo
					endif
				enddo
				enddo
				endif
			enddo
			enddo
			endif
			if(sum(jm*jm).ne.0.0) then
			write(*,"(2i4,6f7.3)") i,j,(listvec(l,k),k=1,d),(jm(1,l2),l2=1,d)
			do l1=2,3
			write(*,"(29x,3f7.3)") (jm(l1,l2),l2=1,d)
			enddo
			endif
			x = x + sp(j)*sum(eta(i,:)*matmul(jm,etao))
!
		enddo
		enddo
		omega(i) = omega(i) + x
		write(*,"(i4,f12.5)") i,omega(i)
	enddo
	endif
!
!	constantes de force
!
	if(couple) then
	qc = 0.0
	call banner2("Ctes de force")
	do i=1,n0
		do j=1,n0
		do l=1,ntt
			call for(i,j,listvec(l,:),pos,e,nz,modele(i,j),disf(:,i,j),ctef(i,j,:,:,:),qc,jm)
			if(sum(jm*jm).ne.0.0) then
				write(*,"(2i4,3f7.3)") i,j,(jm(k,k),k=1,d)
			endif
		enddo
		enddo
	enddo
	endif
!
	call banner2("Omega")
	do i=1,n
		write(*,"(i4,f12.5)") i,omega(i)
	enddo
!
!	-----------------------------------------------------------------------
!	On diagonalise puis calcule l'intensite
!
	allocate(matc(ncalc,d,d))
	allocate(ar(n2,n2),ai(n2,n2),dr(nq,n2))
	allocate(ur(n,n),ui(n,n),vr(n,n),vi(n,n))
	if(cnt) then
		allocate(urkn(nq,n,n),uikn(nq,n,n),vrkn(nq,n,n),vikn(nq,n,n))
	endif
!
	ncalc2 = ncalc
	if(partiel) ncalc2 = ncalc2+n0
	allocate(fq(nq,n,ncalc2))
	allocate(fpartiel(nq,n,n0))
	allocate(redm(nq,n0))
	fq       = 0.0
	redm     = 0.0
	fpartiel = 0.0
	if(couple) then
		allocate(fn(nq,n),fy(nq,n),fz(nq,n))	
		allocate(vkr(n0,n0*d),vki(n0,n0*d))
		fn  = 0.0
		fy  = 0.0
		fz  = 0.0
		vkr = 0.0
		vki = 0.0
	endif
!
!	pour regularisation
	call random_seed
!
!	boucle sur les valeurs de k
!	---------------------------
	call banner2("Boucle sur les valeurs de k")
	write(*,*) "NP  = ", np
	write(*,*) "NQ  = ", nq
	do p = 1,nq
!
!		Vecteur Q		
		qc  = tqc(p,:)
		qc2 = tqr(p,:)
		write(*,"(i8,a5,3e13.5)") p,"k=  ",(tq(p,j),j=1,d)
		write(*,"(i8,a5,3e13.5)") p,"kc= ",(tqc(p,j),j=1,d)
		write(*,"(i8,a5,3e13.5)") p,"kr= ",(tqr(p,j),j=1,d)
		write(*,"(i8,a5,3e13.5)") p,"Q=  ",(listeq(vq(p),j),j=1,d)
!
!		facteur de forme magnetique
!
		fmag = 1.0
		if(fform) then
			s    = sqrt(sum(qc*qc))/(4.0*pi)
			fmag = fga*exp(-fa*s*s)+fgb*exp(-fb*s*s)+fgc*exp(-fc*s*s)+fgd
			fmag = fmag*fmag
		endif
!
!		Matrices pour le calcul des fct de correlation
!		par convention
!		1     : facteur d'orientation
!		2,3,4 : xx yy et zz dans le repere (ax ay az) du cristal
!		5,6,7 : yy,zz et chiral dans le repere defini par (q,qperp,qz)
!
		matc = 0.0
		do i=1,d
			matc(1,i,i) = 1.0
		enddo
		do i=1,d
		do j=1,d
			matc(1,i,j) = matc(1,i,j)-qc(i)*qc(j)/sum(qc*qc)
		enddo
		enddo
		matc(2,1,1) = 1.0
		matc(3,2,2) = 1.0
		matc(4,3,3) = 1.0
!
!		repere neutrons polarises
		ex = qc/sqrt(sum(qc*qc))
		call vct(ez,ex,ey)
		do i=1,d
		do j=1,d
			matc(5,i,j) = ey(i)*ey(j)
			matc(6,i,j) = ez(i)*ez(j)
			matc(7,i,j) = ey(i)*ez(j)
			matc(8,i,j) = matc(7,i,j)
		enddo
		enddo
!
!		remplissage de la matrice a diagonaliser
!		----------------------------------------
!
		ar = 0.00
		ai = 0.00
!
!		contribution magnetique Heisenberg et DM
		do i=1,n0
		do j=1,n0
		do l=1,ntt
!
!			phase
			r0   = matmul(e,listvec(l,:))
			x    = sum(qc*r0)
			xr   = cos(x)
			xi   = sin(x)
!
			zdro = zdr0(j,l,:)
			zdio = zdi0(j,l,:)
			x    = sqrt(sp(i)*sp(j))
!
!			contributions Heisenberg et DM
!
			detail = .false.
			call mag(i,j,listvec(l,:),pos,matrx(i,:,3),matrx(j,:,3),e,nz,modele(i,j),&
			dani,dist(:,i,j),jech(i,j,:,:,:),qc,detail,jm)
			jmr = jm*xr*x/2.0
			jmi = jm*xi*x/2.0
!
!			contribution DM
!
			uur = dm(i,j,l,:)*xr*x/2.0
			uui = dm(i,j,l,:)*xi*x/2.0
!
!			contribution KSEA
!
			jm = 0.0
			do l1=1,d
			do l2=1,d
				jm(l1,l2) = gksea*dm(i,j,l,l1)*dm(i,j,l,l2)
			enddo
			enddo
			do l1=1,d
				jm(l1,l1) = 0.0
			enddo
			jmr = jmr+jm*xr*x/2.0
			jmi = jmi+jm*xi*x/2.0
!
!			contribution GDIAG
!
			jm = 0.0
			do l1=1,d
				jm(l1,l1) = 1.0
			enddo
			jm = gdiag*sum(dm(i,j,l,:)**2)*jm
			jmr = jmr+jm*xr*x/2.0
			jmi = jmi+jm*xi*x/2.0
!
!			contribution biquadratique
!
			if(bbq) then

			jm = 0.0
			allocate(u1(3),u2(3))
			rc = pos(i,:)+listvec(l,:)-pos(j,:)
			rc = matmul(e,rc)
			if(sqrt(sum(rc*rc)).le.dbiq12(i,j)) then
			do k3=1,n
			do l3=1,ntt
				rc = pos(k3,:)-pos(i,:)+listvec(l3,:)
				rc = matmul(e,rc)
				if(sqrt(sum(rc*rc)).le.dbiq13(i,k3)) then
				do k4=1,n
				do l4=1,ntt
					rc = pos(k3,:)-pos(k4,:)+listvec(l3,:)-listvec(l4,:)
					rc = matmul(e,rc)
					if(sqrt(sum(rc*rc)).le.dbiq12(k3,k4)) then
						u1 = sp(k3)*eta0(k3,l3,:)
						u2 = sp(k4)*eta0(k4,l4,:)
						x = 2.0*biq(i,j,k3,k4)*sum(u1*u2)
						do l1=1,d
						jm(l1,l1) = jm(l1,l1)+x
						enddo
					endif
				enddo
				enddo
				endif
			enddo
			enddo
			endif
			if(sqrt(sum(rc*rc)).le.dbiq13(i,j)) then
			do k3=1,n
			do l3=1,ntt
				rc = pos(k3,:)-pos(i,:)+listvec(l3,:)
				rc = matmul(e,rc)
				if(sqrt(sum(rc*rc)).le.dbiq12(i,k3)) then
				do k4=1,n
				do l4=1,ntt
					rc = pos(k4,:)-pos(j,:)+listvec(l4,:)-listvec(l,:)
					rc = matmul(e,rc)
					if(sqrt(sum(rc*rc)).le.dbiq12(j,k4)) then
						u1 = sp(k3)*eta0(k3,l3,:)
						u2 = sp(k4)*eta0(k4,l4,:)
						do l1=1,d
						do l2=1,d
						jm(l1,l2) = jm(l1,l2)+4.0*biq(i,k3,j,k4)*u1(l1)*u2(l2)
						enddo
						enddo
					endif
				enddo
				enddo
				endif
			enddo
			enddo
			endif
			jmr = jmr+jm*xr*x/2.0
			jmi = jmi+jm*xi*x/2.0
			endif
!
!			Fin des contributions
!
!			bloc en z_i zbar_j
!	
			call pdtcplx(d,zdr(i,:), zdi(i,:),jmr,jmi,zdro,-zdio,p1,p2)
			call vctx(zdr(i,:),+zdi(i,:),zdro,-zdio,vdmr,vdmi)
			ar(i  ,j  ) = ar(i  ,j  )+p1+gdm*sum(uur*vdmr-uui*vdmi)
			ai(i  ,j  ) = ai(i  ,j  )+p2+gdm*sum(uur*vdmi+uur*vdmi)
!	
!			bloc en zbar_i zbar_j
!
			call pdtcplx(d,zdr(i,:),-zdi(i,:),jmr,jmi,zdro,-zdio,p1,p2)
			call vctx(zdr(i,:),-zdi(i,:),zdro,-zdio,vdmr,vdmi)
			ar(i+n,j  ) = ar(i+n,j  )+p1+gdm*sum(uur*vdmr-uui*vdmi)
			ai(i+n,j  ) = ai(i+n,j  )+p2+gdm*sum(uur*vdmi+uur*vdmi)
!	
!			bloc en z_i z_j
!
			call pdtcplx(d,zdr(i,:), zdi(i,:),jmr,jmi,zdro, zdio,p1,p2)
			call vctx(zdr(i,:),+zdi(i,:),zdro,+zdio,vdmr,vdmi)
			ar(i  ,j+n) = ar(i  ,j+n)+p1+gdm*sum(uur*vdmr-uui*vdmi)
			ai(i  ,j+n) = ai(i  ,j+n)+p2+gdm*sum(uur*vdmi+uur*vdmi)
!	
!			bloc en zbar_i z_j
!
			call pdtcplx(d,zdr(i,:),-zdi(i,:),jmr,jmi,zdro, zdio,p1,p2)
			call vctx(zdr(i,:),-zdi(i,:),zdro,+zdio,vdmr,vdmi)
			ar(i+n,j+n) = ar(i+n,j+n)+p1+gdm*sum(uur*vdmr-uui*vdmi)
			ai(i+n,j+n) = ai(i+n,j+n)+p2+gdm*sum(uur*vdmi+uur*vdmi)
!		
		enddo
		enddo
		enddo
!
!		on incorpore le vecteur omega
		do i=1,n0
			ar(i  ,i  ) = ar(i  ,i  ) - omega(i)
			ar(i+n,i+n) = ar(i+n,i+n) - omega(i)
		enddo
!
!		contribution phonon
!
		if(couple) then
!			terme d'inertie
!			p^2/2m = (a-a^+)/sqrt(2)i (a-a^+)/sqrt(2)i = -1/2 ( a a - a a^+ - a^+ a + a^+ a^+)
!                  = (a^+ a)M(a a^+)
!			avec :
!			M = {  1/2 -1/2 }
!               { -1/2  1/2 }
!
			do i=1,n0
!				indice corrige
				ii = n0+(i-1)*d
				x = 1.0/(2*ma(i))
				do j=1,d
					ar(  ii+j,  ii+j) = + x/2.
					ar(n+ii+j,  ii+j) = - x/2.
					ar(  ii+j,n+ii+j) = - x/2.
					ar(n+ii+j,n+ii+j) = + x/2.
				enddo
			enddo
!			terme de force
!			u^2 = (a+a^+)/sqrt(2) (a+a^+)/sqrt(2) = 1/2 ( a a + a a^+ + a^+ a + a^+ a^+)
!                  = (a^+ a)M(a a^+)
!			avec :
!			M = {  1/2  1/2 }
!               {  1/2  1/2 }
			do i=1,n0
				ii = n0+(i-1)*d
				do j=1,n0
				do l=1,ntt
!
!					phase
					r0   = matmul(e,listvec(l,:))
					x    = sum(qc*r0)
					xr   = cos(x)
					xi   = sin(x)
!
					call for(i,j,listvec(l,:),pos,e,nz,modele(i,j),disf(:,i,j),ctef(i,j,:,:,:),qc,jm)
					jj = n0+(j-1)*d
					ar(  ii+1:  ii+d,  jj+1:  jj+d) = ar(  ii+1:  ii+d,  jj+1:  jj+d)-jm*xr
					ai(  ii+1:  ii+d,  jj+1:  jj+d) = ai(  ii+1:  ii+d,  jj+1:  jj+d)-jm*xi
					ar(n+ii+1:n+ii+d,  jj+1:  jj+d) = ar(n+ii+1:n+ii+d,  jj+1:  jj+d)-jm*xr
					ai(n+ii+1:n+ii+d,  jj+1:  jj+d) = ai(n+ii+1:n+ii+d,  jj+1:  jj+d)-jm*xi
					ar(  ii+1:  ii+d,n+jj+1:n+jj+d) = ar(  ii+1:  ii+d,n+jj+1:n+jj+d)-jm*xr
					ai(  ii+1:  ii+d,n+jj+1:n+jj+d) = ai(  ii+1:  ii+d,n+jj+1:n+jj+d)-jm*xi
					ar(n+ii+1:n+ii+d,n+jj+1:n+jj+d) = ar(n+ii+1:n+ii+d,n+jj+1:n+jj+d)-jm*xr
					ai(n+ii+1:n+ii+d,n+jj+1:n+jj+d) = ai(n+ii+1:n+ii+d,n+jj+1:n+jj+d)-jm*xi
				enddo
				enddo
				jm0 = 0.0
				ss  = 0.0
				do j=1,n0
				do l=1,ntt
					call for(i,j,listvec(l,:),pos,e,nz,modele(i,j),disf(:,i,j),ctef(i,j,:,:,:),ss,jm)
					jm0 = jm0+jm
				enddo
				enddo
				ar(  ii+1:  ii+d,  ii+1:  ii+d) = ar(  ii+1:  ii+d,  ii+1:  ii+d) + jm0
				ar(n+ii+1:n+ii+d,  ii+1:  ii+d) = ar(n+ii+1:n+ii+d,  ii+1:  ii+d) + jm0
				ar(  ii+1:  ii+d,n+ii+1:n+ii+d) = ar(  ii+1:  ii+d,n+ii+1:n+ii+d) + jm0
				ar(n+ii+1:n+ii+d,n+ii+1:n+ii+d) = ar(n+ii+1:n+ii+d,n+ii+1:n+ii+d) + jm0
			enddo
!
!			terme de couplage magnon phonon
!			B^+k (V_k  a_k + V_k a^+_-k ) + B_-k (U_k  a_k + U_k  a^+_-k ) +
!			a^+k (U_-k B_k + V_-k B^+_-k ) + a_-k (U_-k B_k + V_-k B^+_-k )
!			sous forme matricielle
!
!			|	 V_k      V_k  |			|	 V_k      V_k  |
!			|  U_-k     V_-k   |			|  V_k*     V_k    |
!			|    U_k      U_k  |	soit	|    V_k*     V_k* |
!			|  U_-k     V_-k   |			|  V_k*     V_k    |
!	
!			ici, on suppose V_k = V_-k
!	
!			on a bien un couplage hermitique adjoint
!			V_k ~  g_k z_i
!			U_k ~  g_-k zbar_i = V_k^*
!
			if(g.ne.0.0) then
			do i=1,n0
			do j=1,n0
			do l=1,ntt
!
!				phase
				r0   = matmul(e,listvec(l,:))
				x    = sum(qc*r0)
				xr   = cos(x)
				xi   = sin(x)
!	
				etao = eta0(j,l,:)
				zdro = zdr0(j,l,:)
				zdio = zdi0(j,l,:)
!
				call cpl(i,j,listvec(l,:),pos,e,nz,modele(i,j),disf(:,i,j),ctef(i,j,:,:,:),qc,uo,uu)
				ii = (i-1)*d
				jj = (j-1)*d
				p1 = +sum(eta(i,:)*zdro    )*sp(i)*sqrt(sp(j))
				p2 = -sum(eta(i,:)*zdio    )*sp(i)*sqrt(sp(j))
				p3 = +sum(etao    *zdr(i,:))*sp(j)*sqrt(sp(i))
				p4 = -sum(etao    *zdi(i,:))*sp(j)*sqrt(sp(i))
!				write(*,"(2i4,13f8.3)") i,j,p1,p2,p3,p4,uu,uur,uui
				vkr(i,jj+1:jj+d) = g*(p1*uu*xr-p2*uu*xi)
				vki(i,jj+1:jj+d) = g*(p1*uu*xi+p2*uu*xr)
				vkr(i,ii+1:ii+d) = vkr(i,ii+1:ii+d) + g*p3*uu
				vki(i,ii+1:ii+d) = vki(i,ii+1:ii+d) + g*p4*uu
			enddo
			enddo
			enddo
			do i=1,n0
				ar(i,n0+1:n)      = ar(i,n0+1:n)      + vkr(i,:)
				ai(i,n0+1:n)      = ai(i,n0+1:n)      + vki(i,:)
				ar(i,n+n0+1:n2)   = ar(i,n+n0+1:n2)   + vkr(i,:)
				ai(i,n+n0+1:n2)   = ai(i,n+n0+1:n2)   + vki(i,:)
	
				ar(n0+1:n,  i)    = ar(n0+1:n,  i)    + vkr(i,:)
				ai(n0+1:n,  i)    = ai(n0+1:n,  i)    - vki(i,:)
				ar(n0+1:n,n+i)    = ar(n0+1:n,n+i)    + vkr(i,:)
				ai(n0+1:n,n+i)    = ai(n0+1:n,n+i)    + vki(i,:)

				ar(n+i,n0+1:n)    = ar(n+i,n0+1:n)    + vkr(i,:)
				ai(n+i,n0+1:n)    = ai(n+i,n0+1:n)    - vki(i,:)
				ar(n+i,n+n0+1:n2) = ar(n+i,n+n0+1:n2) + vkr(i,:)
				ai(n+i,n+n0+1:n2) = ai(n+i,n+n0+1:n2) - vki(i,:)

				ar(n+n0+1:n2,i)   = ar(n+n0+1:n2,i)   + vkr(i,:)
				ai(n+n0+1:n2,i)   = ai(n+n0+1:n2,i)   - vki(i,:)
				ar(n+n0+1:n2,n+i) = ar(n+n0+1:n2,n+i) + vkr(i,:)
				ai(n+n0+1:n2,n+i) = ai(n+n0+1:n2,n+i) + vki(i,:)
			enddo
			endif	
		endif
!	
!		diagonalisation et calcul des vecteurs propres
!		----------------------------------------------
		if(affich) then
			write(*,*) "terme magnetique Re"
			do i=1,n0
			write(*,"(24f8.3)") (ar(i,j),j=1,n0)
			enddo
			write(*,*) "terme magnetique Im"
			do i=1,n0
			write(*,"(24f8.3)") (ai(i,j),j=1,n0)
			enddo
			write(*,*) "terme phononique Re"
			do i=n0+1,n
			write(*,"(24f8.3)") (ar(i,j),j=n0+1,n)
			enddo
			write(*,*) "terme phononique Im"
			do i=n0+1,n
			write(*,"(24f8.3)") (ai(i,j),j=n0+1,n)
			enddo
			write(*,*) "terme de couplage Re"
			do i=1,n0
			write(*,"(24f8.3)") (ar(n+i,j),j=n0+1,n)
			enddo
			write(*,*) "terme de couplage Im"
			do i=1,n0
			write(*,"(24f8.3)") (ai(n+i,j),j=n0+1,n)
			enddo
			write(*,*) "terme global Re"
			do i=1,n2
			write(*,"(24f8.3)") (ar(i,j),j=1,n2)
			enddo
			write(*,*) "terme global Im"
			do i=1,n2
			write(*,"(24f8.3)") (ai(i,j),j=1,n2)
			enddo
		endif
!
!		terme de regularisation
		do i=1,n2
			call random_number(x)
			ar(i,i) = ar(i,i) + reg2*x
			do j=1,i-1
				call random_number(x)
				ai(i,j) =  ai(i,j) + reg1*(x-0.5)
				ai(j,i) = -ai(i,j)
			enddo
		enddo
		if(couple) then
			do i=n0+1,n
				call random_number(x)
				ar(i,i) = ar(i,i) + reg3*x
			enddo
			do i=n+n0+1,n2
				call random_number(x)
				ar(i,i) = ar(i,i) + reg3*x
			enddo
		endif
!		pour ancienne methode de calcul		
		if(calc.eq.1) then
			do i=1,n
				ar(i+n,:) = -ar(i+n,:)
				ai(i+n,:) = -ai(i+n,:)
			enddo
			call diagboson0(n2,ar,ai,dr(p,:),ur,ui,vr,vi)
		endif
		if(calc.eq.2) call diagboson (n2,ar,ai,dr(p,:),ur,ui,vr,vi)
!

		if(cnt) then
			urkn(p,:,:) = ur
			uikn(p,:,:) = ui
			vrkn(p,:,:) = vr
			vikn(p,:,:) = vi
		endif
!
!		Pour affichage des vecteurs et valeurs propres
!		do i=1,n
!			write(*,"(i4,f12.5)") i,dr(p,i)
!			do j=1,n
!				write(*,"(2i4,4f12.5)") i,j,ur(j,i),ui(j,i),vr(j,i),vi(j,i)
!			enddo
!		enddo
!	
!		--------------------------------------------------------------------
!		Calcul de la section efficace selon (U* V*) h (U V) avec h par blocs
!	
		do ll=1,ncalc
!
!			calcul\E9 avec le qc
			mat = matc(ll,:,:)
!
!			construction de h =
!			z_i    M zbar_j  z_i    M z_j
!			zbar_i M zbar_j	 zbar_i M z_j
			allocate(a2r(n2,n2),a2i(n2,n2))
			a2r = 0.0
			a2i = 0.0
			do i=1,n0
			do j=1,n0
!
				zdro = zdr(j,:)
				zdio = zdi(j,:)
!				option incommensurable
				if(incom) then
!
					call rodrigue_mat(uprop(i,:),nn1,rr1,ri1)
					call rodrigue_mat(uprop(j,:),nn2,rr2,ri2)
!
!					W = Rr - I Ri
					if(p.le.nq2) then
						jmr = matmul(nn1,matmul(mat,nn2))
						jmi = 0.0
					endif
					if((p.le.2*nq2).and.(p.gt.nq2)) then
						jmr =  (matmul(rr1,matmul(mat,rr2))-matmul(ri1,matmul(mat,ri2)))/4.0
						jmi =  (matmul(ri1,matmul(mat,rr2))+matmul(rr1,matmul(mat,ri2)))/4.0
					endif
					if((p.le.3*nq2).and.(p.gt.2*nq2)) then
						jmr =  (matmul(rr1,matmul(mat,rr2))-matmul(ri1,matmul(mat,ri2)))/4.0
						jmi = -(matmul(ri1,matmul(mat,rr2))+matmul(rr1,matmul(mat,ri2)))/4.0
					endif
!
				else
					jmr = mat
					jmi = 0.0
				endif
!
				call pdtcplx(d,zdr(i,:), zdi(i,:),jmr,jmi,zdro,-zdio,p1,p2)
				a2r(i  ,j  ) = a2r(i  ,j  )+p1
				a2i(i  ,j  ) = a2i(i  ,j  )+p2
				call pdtcplx(d,zdr(i,:),-zdi(i,:),jmr,jmi,zdro,-zdio,p1,p2)
				a2r(i+n,j  ) = a2r(i+n,j  )+p1
				a2i(i+n,j  ) = a2i(i+n,j  )+p2
				call pdtcplx(d,zdr(i,:), zdi(i,:),jmr,jmi,zdro, zdio,p1,p2)
				a2r(i  ,j+n) = a2r(i  ,j+n)+p1
				a2i(i  ,j+n) = a2i(i  ,j+n)+p2
				call pdtcplx(d,zdr(i,:),-zdi(i,:),jmr,jmi,zdro, zdio,p1,p2)
				a2r(i+n,j+n) = a2r(i+n,j+n)+p1
				a2i(i+n,j+n) = a2i(i+n,j+n)+p2
!
			enddo
			enddo
!
			do k = 1,n
				allocate(xpr(n2),xpi(n2))
				xpr = 0.0
				xpi = 0.0		
				do i = 1,n0
!					vecteur (U,V) stocke dans (xpr + i xpi)e*{-iq u_i}
!					phase liee a la position des spins
!					ce signe - dans l'expression de p1 est correct
					rc       = matmul(e,pos(i,:))
					p1       = -sum(qc *rc)
					xpr(i)   = (ur(i,k)*cos(p1)-ui(i,k)*sin(p1))*sqrt(sp(i))
					xpr(i+n) = (vr(i,k)*cos(p1)-vi(i,k)*sin(p1))*sqrt(sp(i))
					xpi(i)   = (ur(i,k)*sin(p1)+ui(i,k)*cos(p1))*sqrt(sp(i))
					xpi(i+n) = (vr(i,k)*sin(p1)+vi(i,k)*cos(p1))*sqrt(sp(i))
!
!					calcul des partiels
					fpartiel(p,k,i) = ur(i,k)**2+ui(i,k)**2+vr(i,k)**2+vi(i,k)**2
				enddo
!
!				calcul de (U* V*) h (U V)
				call pdtcplx(n2,xpr,-xpi,a2r,a2i,xpr,xpi,xr,xi)
				deallocate(xpr,xpi)
					
				fq(p,k,ll) = xr*fmag/2.0
				if(ll.eq.ncalc) fq(p,k,ll) = xi*fmag
			enddo
!
			deallocate(a2r,a2i)
		enddo
		if(partiel) then
			do k = 1,n
			do i = 1,n0
				fq(p,k,ncalc+i) = fpartiel(p,k,i)
			enddo
			enddo
		endif
!	
!		reduction du moment par fluctuations quantiques
!		-----------------------------------------------
!
		if(reduc) then
			do i=1,n0
			redm(p,i) = 0.0
			do k=1,n
				redm(p,i) = redm(p,i) + (vr(i,k)**2+vi(i,k)**2)
			enddo
			redm(p,i) = redm(p,i)/n
			enddo
		endif
!
!		-----------------------------------------------
!		dans le cas du couplage magnon-phonon

		if(couple) then
		allocate(a2r(n2,n2),a2i(n2,n2))
		a2r = 0.0
		a2i = 0.0
		do l1=1,d
		do l2=1,d
			jm(l1,l2) = qc(l1)*qc(l2)
		enddo
		enddo
		do i = 1,n0
		do j = 1,n0
			ii = n0+(i-1)*d
			jj = n0+(j-1)*d
			a2r(  ii+1:  ii+d,  jj+1:  jj+d) = jm
			a2r(n+ii+1:n+ii+d,  jj+1:  jj+d) = jm
			a2r(  ii+1:  ii+d,n+jj+1:n+jj+d) = jm
			a2r(n+ii+1:n+ii+d,n+jj+1:n+jj+d) = jm
		enddo		
		enddo
		do k = 1,n
			allocate(xpr(n2),xpi(n2))
			xpr = 0.0
			xpi = 0.0		
			do i = 1,n0
!				vecteur (U,V) stocke dans (xpr + i xpi)e*{-iq u_i}
!				phase liee a la position des spins
!				ce signe - dans l'expression de p1 est correct
				rc = matmul(e,pos(i,:))
				p1 = -sum(qc*rc)
				ii = n0+(i-1)*d
				xpr(  ii+1:  ii+d) = ur(ii+1:ii+d,k)*cos(p1)-ui(ii+1:ii+d,k)*sin(p1)
				xpr(n+ii+1:n+ii+d) = vr(ii+1:ii+d,k)*cos(p1)-vi(ii+1:ii+d,k)*sin(p1)
				xpi(  ii+1:  ii+d) = ur(ii+1:ii+d,k)*sin(p1)+ui(ii+1:ii+d,k)*cos(p1)
				xpi(n+ii+1:n+ii+d) = vr(ii+1:ii+d,k)*sin(p1)+vi(ii+1:ii+d,k)*cos(p1)
			enddo
!
!			calcul de (U* V*) h (U V)
			call pdtcplx(n2,xpr,-xpi,a2r,a2i,xpr,xpi,xr,xi)
			deallocate(xpr,xpi)		
			fn(p,k) = xr
		enddo
		deallocate(a2r,a2i)
		endif
!
!		Affichage des resultats
!		------------------------
!
		if(couple) then
			do k=1,n
				write(*,"(13e13.5)") dr(p,k),(fq(p,k,ll),ll=1,ncalc),fn(p,k),fy(p,k),fz(p,k)
				if(partiel) write(*,"(13e13.5)") (fpartiel(p,k,i),i=1,n0)
			enddo
			if(reduc) then
				write(*,*)
				write(*,"(13e13.5)") (redm(p,i),i=1,n0)
			endif
		else
			do k=1,n
				write(*,"(10e13.5)") dr(p,k),(fq(p,k,ll),ll=1,ncalc)
				if(partiel) write(*,"(13e13.5)") (fpartiel(p,k,i),i=1,n0)
			enddo
			if(reduc) then
				write(*,*)
				write(*,"(13e13.5)") (redm(p,i),i=1,n0)
			endif
		endif
!
!		fin de la boucle sur le vct d'onde k
	enddo
!
	if(reduc) then
		write(*,*)
		write(*,*) "Reduction du moment"
		write(*,"(a4,a13)") "site","s/<S>"
		do i=1,n0
			write(*,"(i4,e13.5)") i,(1.0-sum(redm(:,i))/nq)
		enddo
	endif
!	
!	------------------------------------------------
!	Affichage (Q,omega) final
!	
	write(*,*)
	call banner("AFFICHAGE FINAL")
!	
	allocate(pou(nw,ncalc2),poud(nw),pous(nw))
	if(couple) then
		allocate(poun(nw),pouy(nw),pouz(nw))
	endif
!
!	conversion en fwhm
!	et facteur de normalisation des gaussiennes
	sig  = sig/sqrt(4.0*log(2.0))
	vsig = 1.0/(sqrt(pi)*sig)
!
	do p=1,np
		pou  = 0.0
		pous = 0.0
		poud = 0.0
!
!		continuum \E0 2 magnons (Manila)
!		Calcul des termes du continuum
!		Terme 1 seul ...
!		------------------------------
		if(cnt) then
!			
!			matrice d'orientation correspondant \E0 Q	
			qc  = matmul(ee,listeq(p,:))
			mat = 0.0
			if(sum(qc*qc).ne.0.0) then
			do i=1,d
				mat(i,i) = 1
				do j=1,d
					mat(i,j) = mat(i,j)-qc(i)*qc(j)/sum(qc*qc)
				enddo
			enddo
			endif
!
!			regroupement des valeurs de k pour Q
!			on somme sur k (pk) et on utilile k-Q (pq)
			do ii=2,ntab(p)
				pq = tab(p,ii)
				pk = vk(pq)
				write(*,"(3(a5,i5),2(a5,3f8.3))") "Q=",p,"k+Q=",pq,"pk=",pk,"k+Q=",tq(pq,:),"k=",tq(pk,:)
				do i=1,n0
				do j=1,n0
				rc = matmul(e,pos(i,:)-pos(j,:))
				p1 = -sum(qc*rc)
				p2 = sum(eta(i,:)*matmul(mat(:,:),eta(j,:)))	
				do l1=1,n
				do l2=1,n
!					Vkil UQ+kil'			
					xr = vrkn(pk,i,l1)*urkn(pq,i,l2)-vikn(pk,i,l1)*uikn(pq,i,l2)
					xi = vrkn(pk,i,l1)*uikn(pq,i,l2)+vikn(pk,i,l1)*urkn(pq,i,l2)
!
!					U*kjl V*Q+kjl'					
					yr = urkn(pk,j,l1)*vrkn(pq,j,l2)-uikn(pk,j,l1)*vikn(pq,j,l2)
					yi =-urkn(pk,j,l1)*vikn(pq,j,l2)-uikn(pk,j,l1)*vrkn(pq,j,l2)
!
!					U*Q+kjl' V*kjl									
					yr = yr + urkn(pq,j,l2)*vrkn(pk,j,l1)-uikn(pq,j,l2)*vikn(pk,j,l1)
					yi = yi - urkn(pq,j,l2)*vikn(pk,j,l1)-uikn(pq,j,l2)*vrkn(pk,j,l1)
!
					pdtuvr = xr*yr-xi*yi
					pdtuvi = xr*yi+xi*yr
!
					x = p2*(pdtuvr*cos(p1)-pdtuvi*sin(p1))/(ntab(p)-2+1)
					y = p2*(pdtuvr*sin(p1)+pdtuvi*cos(p1))/(ntab(p)-2+1)
!
					do l=1,nw
						z = dr(pk,l1)+dr(pq,l2)
						pous(l) = pous(l) + x*exp(-((w(l)-z)/sig)**2)
					enddo
				enddo
				enddo					
				enddo
				enddo
			enddo
		endif		
!	
		do i=1,nw
!			regroupement des valeur de k pour Q	
			do ii=1,ntab(p)
				pq = tab(p,ii)
				do ll=1,ncalc2	
					pou(i,ll) = pou(i,ll) + su(pq)*sum(fq(pq,:,ll)*exp(-((w(i)-dr(pq,:))/sig)**2))
				enddo
			enddo
!			
!			dispersion principale			
			pq = tab(p,1)
			poud(i) = sum(exp(-((w(i)-dr(pq,:))/sig)**2))
!			
!			couplage magnon-phonon
			if(couple) then
			poun = 0.0
			pouy = 0.0
			pouz = 0.0
			do ii=1,ntab(p)
				pq = tab(p,ii)
				poun(i) = poun(i) + sum(fn(pq,:)*exp(-((w(i)-dr(pq,:))/sig)**2))
				pouy(i) = pouy(i) + sum(fy(pq,:)*exp(-((w(i)-dr(pq,:))/sig)**2))
				pouz(i) = pouz(i) + sum(fz(pq,:)*exp(-((w(i)-dr(pq,:))/sig)**2))
			enddo
			endif	
!
!			affichage
			if(tof) then		
				q = matmul(ee,listeq(p,:))
				x = sqrt(sum(q*q))
				write(lu,"(32(e13.5))") x,w(i),pou(i,1)
			else
				if(couple) then
				write(lu,"(33(e13.5))") (listeq(p,ll),ll=1,d),w(i),&
				(pou(i,ll),ll=1,ncalc2),pous(i),poun(i),pouy(i),pouz(i),poud(i)
				else
				write(lu,"(33(e13.5))") (listeq(p,ll),ll=1,d),w(i),&
				(pou(i,ll),ll=1,ncalc2),pous(i),poud(i)
				endif
			endif
		enddo
		if(simpleq) write(lu,*)
		if(coupe) then
			if(mod(p,np2).eq.0) write(lu,*)
		endif
	enddo	
!
!	
	if(tof.and.integ) then
!
!		Int\E9gration en \E9nergie de 0 a emax (=unit*ki**2)
!		unit = 2.0 (mev <-> A-1) ki = 2.662 par defaut
		emax = unit1*(ki**2)
		j=1
		do i=1,nw
			if(w(i).lt.emax) j=i
		enddo
		do p = 1,np
!			valeur de Q puis de cos(2theta)
			q = matmul(ee,listeq(p,:))
			x = sqrt(sum(q*q))
			if(x.le.(2*ki)) then
				p1 = 1.0-0.5*(x/ki)**2
!				integration en energie jusqu'a emax
				p3 = 0.0
				do i=1,j
!					pour cos(2theta) et w(i), on calcul p2=module de Q
					p2 = ki*sqrt(2.0-w(i)/emax-2.0*p1*sqrt(1.0-w(i)/emax))
!					on cherche le module de Q le plus proche
!					et on integre
					xr = qmin
					do k=1,np
					q = matmul(ee,listeq(k,:))
					xi = sqrt(sum(q*q))
						if((p2.ge.xr).and.(p2.lt.xi)) then	
							do ii=1,ntab(k)
								pk = tab(k,ii)
								p3 = p3 + su(pk)*sum(fq(pk,:,1)*exp(-((w(i)-dr(pk,:))/sig)**2))
							enddo
						endif
						xr = xi
					enddo
				enddo
				write(lu,"(3(e13.5))") x,p3,acos(p1)*180./pi
			endif
		enddo
	endif
!
	close(lu)
!
!	------------------------------------------------
!	Fin du programme
!
	deallocate(listvec,dist,jech,dech,dval,h,dm,pt,matl,matr,ma,ctef,disf)
	if(nato.ne.0) deallocate(lato)
	deallocate(omega,ar,ai,dr,ur,ui,vr,vi)
	deallocate(pos,pos1,zdr,zdi,eta,w,fq,tq,matc,fpartiel)
	if(couple) deallocate(fn,fy,fz,vkr,vki)
	if(reduc) deallocate(redm)
	if(cnt) deallocate(urkn,uikn,vrkn,vikn)
	end program sw2
	
	
