program ising
implicit none
!size
integer, parameter :: size = 100
!numero di misure e tempo di decorrelazione
integer, parameter :: nmis = 10000
integer, parameter :: ndecorr = size*size
!il lattice andrà da 1 a size
integer lattice(size, size)
!ottimizzazione per le condizioni al bordo
integer modul(0:size+1)

real*8 magnetizzazione(nmis)
real*8 energia(nmis)
real*8 chi(nmis-1)

real*8 somma1, somma2, somma3,energia_media,magnetizzazione_media,energia_sigma,magnetizzazione_sigma,b,T,tau_integrato
integer iseed, i,j,somma,magnet,ener,idecorr,imis,nT

iseed = -11

call random_seed()

!inizializza modul (che fa il modulo)
do i=1, size
    modul(i) = i
enddo
modul(0) = size
modul(size+1) = 1

T=0

if (T .le. 2.0) then
	!setup a freddo
	lattice = 1
	magnet = 1*size*size
	!L'energia sarebbe -4J*size^2 /2
	!perchè ogni sito interagisce con i 4 intorno, ma conto ogni coppia 2 volte. J=1
	ener = -2*size*size
	!fine partenza fredda
else
	!partenza calda
	magnet = 0
	ener = 0
	!setup e calcola magnetizzazione
	do i=1, size
		do j=1,size
			lattice(i,j) = 1
			magnet = magnet+1
			call random_number(somma1)
			if (somma1 .le. 0.5) then
				lattice(i,j) = -1
				magnet = magnet-2
			endif
		enddo
	enddo
	!calcola energia
	do i=1, size
		do j=1,size
			somma = lattice(i,modul(j+1)) + lattice(i,modul(j-1)) + lattice(modul(i+1),j) + lattice(modul(i-1),j)
			ener = ener - lattice(i,j) * somma
		enddo
	enddo
	ener = ener/2
	!fine partenza calda
endif

do nT=0,40
T = real(nT)*0.1
!beta*J. critico è 0.440687
b=1./T

write(*,*) "Temperatura: ", T

!termalizzazione
do imis=1, 500
    do idecorr=1, ndecorr
        call step(lattice, magnet, ener)
    enddo
enddo

write(*,*) "Magnet", magnet

do imis=1, nmis
    do idecorr=1, ndecorr
        call step(lattice, magnet, ener)
    enddo
    magnetizzazione(imis) = real(magnet)/(size*size)
    energia(imis) = real(ener)/(size*size)
    write(10*nT+1,*) magnetizzazione(imis)!, energia(imis)
enddo

!correlazione della magnetizzazione
do i=0, nmis-1
    somma1 = 0
    somma2 = 0
    somma3 = 0
    do j=1, (nmis-i)
        somma1 = somma1 + magnetizzazione(j)*magnetizzazione(i+j)
        somma2 = somma2 + magnetizzazione(j)
        somma3 = somma3 + magnetizzazione(i+j)
    enddo
    chi(i+1) = somma1 - 1./(nmis-i) * somma2 * somma3
    write(10*nT+2,*) chi(i+1)/chi(1)
enddo

tau_integrato = 0
do i=1, nmis
    !somma solo finchè si hanno contributi significativi
    if (chi(i) < chi(1)*0.01) then
        exit
    endif
    tau_integrato = tau_integrato + chi(i)
enddo
tau_integrato = tau_integrato/chi(1)
write(3,*) "Tau integrato magnet:", tau_integrato

!ora calcoliamo media e varianza della magnetizzazione
magnetizzazione_media = sum(magnetizzazione)/nmis
magnetizzazione_sigma = sqrt((sum(magnetizzazione*magnetizzazione)/nmis - magnetizzazione_media*magnetizzazione_media)&
    *((2*tau_integrato+1)/nmis))

!rifare per energia
!do i=0, nmis-1
!    somma1 = 0
!    somma2 = 0
!    somma3 = 0
!    do j=1, (nmis-i)
!        somma1 = somma1 + energia(j)*energia(i+j)
!        somma2 = somma2 + energia(j)
!        somma3 = somma3 + energia(i+j)
!    enddo
!    chi(i+1) = somma1 - 1./(nmis-i) * somma2 * somma3
!enddo

!tau_integrato = 0
!do i=1, nmis
!    !somma solo finchè si hanno contributi significativi
!    if (chi(i) < chi(1)*0.1) then
!        exit
!    endif
!    tau_integrato = tau_integrato + chi(i)
!    !write(40,*) chi(i)/chi(1)
!enddo
!tau_integrato = tau_integrato/chi(1)
!write(*,*) "Tau integrato energia:", tau_integrato

!energia_media = sum(energia)/nmis
!energia_sigma = sqrt((sum(energia*energia)/nmis - energia_media*energia_media)*((2*tau_integrato+1)/nmis))
!end calcoli energia

!write(*,*) "Energia:", energia_media,"+/-",energia_sigma
write(4,*) "Magnetizzazione:", magnetizzazione_media,"+/-",magnetizzazione_sigma

!write(50,*) energia_media, magnetizzazione_media, energia_sigma, magnetizzazione_sigma

enddo

print *, char(7)

contains
subroutine step(lattice, magnet, ener)
    real*8 u,v,rand
    integer i,j,somma,esp,magnet,ener
    integer lattice(size, size)
    !il sito del lattice da considerare
    !notare che u,v<1 strettamente
    !u = ran(iseed)
    !v = ran(iseed)
    call random_number(u)
    call random_number(v)
    i = floor(size*u)+1
    j = floor(size*v)+1
    !calcolo energia e accettanza
    somma = lattice(i,modul(j+1)) + lattice(i,modul(j-1)) + lattice(modul(i+1),j) + lattice(modul(i-1),j)
    esp = 2*lattice(i,j)*somma !differenza di energia
    !se è negativo, accetta. altrimenti
    if (esp>0) then
        !accettanza = exp(-b*esp)
        !estrazione di numero random
        !rand = ran(iseed)
        call random_number(rand)
        if (exp(-b*esp)<rand) then!reject
            return
        endif
    endif
    !ora accetta il cambiamento del sito
    !e aggiorna la magnetizzazione
    lattice(i,j) = lattice(i,j)*(-1)
    magnet = magnet + 2*lattice(i,j)
    ener = ener + esp
end subroutine step

!“Minimal” random number generator of Park and Miller combined with a Marsaglia shiftsequence. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpointvalues). This fully portable, scalar generator has the “traditional” (not Fortran 90) callingsequence with a random deviate as the returned function value: call withidumanegativeinteger to initialize; thereafter, do not alteridumexcept to reinitialize. The period of thisgenerator is about 3.1×10E18
FUNCTION ran(idum)
    IMPLICIT NONE
    INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B), INTENT(INOUT) :: idum
    REAL :: ran
    INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    REAL, SAVE :: am
    INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
    if (idum <= 0 .or. iy < 0) then
        am=nearest(1.0,-1.0)/IM
        iy=ior(ieor(888889999,abs(idum)),1)
        ix=ieor(777755555,abs(idum))
        idum=abs(idum)+1
    end if
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/IQ
    iy=IA*(iy-k*IQ)-IR*k
    if (iy < 0) iy=iy+IM
    ran=am*ior(iand(IM,ieor(ix,iy)),1)
END FUNCTION ran
end program ising
