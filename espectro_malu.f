c       extensao dos espectros para o UV e NIR
c       laerte 060114 181014
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	parameter(nspec=2639,nl=5501,nbase=150,nlambdae=19001)
	PARAMETER (imw=12000,maxt=2000)
	REAL lambda(imw),flux(imw),x(nbase),flux0(imw),lambdac(imw)
	character*32 nomeBSc

	real spec(imw,maxt),alam(imw)
	integer ns(25)
	character*80 file0
	character*80 nomebase(6)
	real fluxbase(imw,maxt),fluxsint(imw),linear,fbase(imw,maxt)
        real q_CCM(nlambdae),lambdae(nlambdae),qext(imw)
	real lnorm,lini,lfim

	real xx(imw),yy(imw)


      call cpu_time(tempo1)

c       arquivo de saida
c       
      open(unit=10,file='uv2ir.res',status='unknown')

c       lambda de normalizacao
	lnorm=4020.

c     lei de extincao de Cardelli et al
      open(unit=1,file='ExtLaws.out',status='old')
	lini=1.e30
	lfim=-1.e30
      read(1,*)
      do i=1,nlambdae
         read(1,*)lambdae(i),q_CCM(i)
	 lini=min(lini,lambdae(i))
	 lfim=max(lfim,lambdae(i))
      enddo
      close(1) 

c       base B&C
c       identificacao da base
c       
c       nsl: no. de idades na base
c       nmetal: no. de arquivos com SEDs de B&C
c       o arquivo 'base' contem:
c       - no, idade (yr), metalicidade, indice das 221 idades
c       - nomes dos arquivos com SEDs de B&C
	nsl=25
	nmetal=6

	open(unit=1,file='base',status='old')
	do i=1,nsl
	   read(1,*)ii,ageb,zb,ns(i)
	enddo
	do i=1,nmetal
	   read(1,*)nomebase(i)
	enddo
	close(unit=1)

	nb=0
	do i=1,nmetal
c       le os espectros de B&C para uma data metalicidade da base
	   file0=nomebase(i)
	   call lbc(file0,alam,spec,ll,ks,lnorm)

c       seleciona os espectros pelas idades da base
	   do j=1,nsl
	      nb=nb+1
	      ii=ns(j)
	      kk=0
	      do k=1,ll
		 if(alam(k).ge.lini.and.alam(k).le.lfim)then
		    kk=kk+1
		    lambdac(kk)=alam(k)
		    fbase(kk,nb)=spec(k,ii)
		 endif
	      enddo
	   enddo   

	enddo

	write(10,*)kk,(lambdac(j),j=1,kk)

	do j=1,kk
   	   qext(j)=LINEAR(lambdac(j),lambdae,q_CCM,nlambdae)
	enddo

c       leitura do output do Starlight
	open(unit=1,file='list_malu.dat',status='old')
	do i=1,nspec
	   read(1,'(a32)')nomeBSc
           print*,i,trim(nomeBSc)
	   open(unit=2,file=trim(nomeBSc)//'.7xt.sc5.C11.im.BS',
     i status='old')
	   do j=1,57
	      read(2,*)
	   enddo  
	   read(2,*)v0 
	   read(2,*)vd 
	   read(2,*)Av
	   do j=1,3
	      read(2,*)
	   enddo  
	   s=0.
	   do j=1,nbase
	      read(2,*)jj,x(j),a1,a2,age,zmetal
	      s=s+x(j)
	   enddo 

	   do j=214,531
	      read(2,*)
	   enddo 

	   nll=0
	   do j=1,nl
	      read(2,*)a1,a2,a3,a4,a5
	      if(a4.ge.0..and.a1.ge.lini.and.a1.le.lfim)then
		 nll=nll+1
		 lambda(nll)=a1
		 flux(nll)=a2
		 if(lambda(nll).eq.lnorm)inorm=nll
	      endif
	   enddo  
	   do j=1,nll 
	      flux(j)=flux(j)/flux(inorm)
	   enddo

	   close(unit=2)

c       sintese  
	do k=1,kk
	   fluxsint(k)=0.
	   do j=1,nb
		 fluxsint(k)=fluxsint(k)+x(j)*fbase(k,j)
	   enddo  
	   fluxsint(k)=fluxsint(k)/s*10.**(-0.4*Av*(qext(k)-qext(inorm)))
	enddo 

	write(10,*)nomeBSc
	write(10,*)(fluxsint(k),k=1,kk)

C	   call plota2l(nll,lambda,flux,kk,lambdac,fluxsint)

	enddo
	close(unit=1)

      call cpu_time(tempo2)
      write(10,*)'cpu time: ',tempo2-tempo1
      write(*,*)'cpu time: ',tempo2-tempo1

	close(unit=10)

	END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine lbc(file0,w,spec,iw,ks,lnorm)
c	Reads galaxy s.e.d.'s = output file from program GALAXEV.

	PARAMETER (imw=12000)
	parameter (maxt=2000)
	real ta(maxt),w(imw),f(imw,maxt),spec(imw,maxt),alam(imw)
	character*80 file0
	real aux(52)
	real lnorm

c	Open file   arquivo B&C de uma dada metalicidade
	file0='./BC03/'//file0
	open (2,file=file0,status='old')

c	Read time scale
	read (2,*) ks,(ta(i),i=1,ks)
	do i=1,5
	   read(2,*)
	enddo  

c	Read wavelength scale
	read (2,*) iw,(w(i),i=1,iw)

c       identifica indice do lambda de normalizacao
	do i=1,iw
	   if(w(i).eq.lnorm)inorm=i
	enddo

c	Read spectra (para todas as idades)
	do k=1,ks 
	read (2,*) iw,(f(i,k),i=1,iw),kk,(aux(i),i=1,kk)
	enddo

c       normalizacao dos espectros da base
	do k=1,ks
	   do l=1,iw
	      spec(l,k)=f(l,k)/f(inorm,k)
	   enddo   
	enddo
   
	return
	end
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL FUNCTION LINEAR (X0,X,Y,N)
c     INTERPOLATES LINEARLY THE FUNCTION Y(X) AT X=X0.
c     minha adaptacao laerte 100210
      REAL X(N),Y(N)

      do i=1,n-1

c        IF (X0.EQ.X(I).OR.X(I).EQ.X(I+1))then corrigido 060910
        IF (X0.EQ.X(I).OR.X(I).EQ.X(I+1))then
            LINEAR=Y(I)
            RETURN
        ENDIF

        IF (X0.gt.X(I).and.X0.le.X(I+1))then
            LINEAR=Y(I) + (Y(I+1)-Y(I))*(X0-X(I))/(X(I+1)-X(I))
            RETURN
        ENDIF

      enddo

      write(*,*)'LINEAR: nao devia imprimir isso!!'
      write(*,*)n,x0
      linear=-999.
c      stop

      END
