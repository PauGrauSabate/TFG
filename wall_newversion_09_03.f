!23456789012345678901234567890123456789012345678901234567890123456789012
!
      program wall_new
      implicit none
! posiciones
      real(8) :: x(1000),y(1000)
! posiciones justo despues de un update (para calcular la velocidad instantanea)
      real(8) :: xn(1000),yn(1000)
      real(8) :: vx(1000),vy(1000)
      real(8) :: dt
! time step
      integer :: it,tt,i,j,k,sqn,npar,nit,step
! nit= numero de pasos de simu
      real(8) :: rho,ls,omega,hueco,h
      real(8) :: f0,sigma,cutoff,visc
      real(8) :: dx,dy,dist,displx, disply,fmodulo
      real(8) :: fx(1000,1000),fy(1000,1000),fxj(1000),fyj(1000)
      real(8) ::  arg1,arg2,arg3
! fx(i,j) fuerza de j sobre i en la direccion x --> f_{ij}^x
! fxj(i) fuerza total sobre la particula i --->F_i^x=\Sigma_{j\neq i}f_{ij}^x
      real(8) :: ux(1000,1000),uy(1000,1000)
c       real*8 ux(1000,1000),uy(1000,1000),inicial(1000,5)
c         real*8 vcolx(1000),vcoly(1000)
      character(5) :: numf
c declara las funciones
      real(8) :: force,fuerza, vrx,vry,vix,viy, fuerzaREP
      real(8) :: rseed
c vrx Rotlet en la direccion x
c vix Stokeslet en la direccion x

c lee parámetros al bulto
      open(14,file="input-parameters.dat")
      read(14,*) dt,nit,step,rho,sqn,omega,h
c 1.0d-3 10000000 10000 0.1 10 1.0d-1 1
      close(14)
c      write(*,*) dt,nit,ls,sqn,h,omega
c lee parámetros nuevos del potencial
      open(15,file="input-potential.dat")
c f0=epsilon L-J
      read(15,*) f0,sigma,cutoff,visc
c 1.0d-1 2.0d0 6.0d0 1.0d0
      close(15)

      npar=sqn*sqn
c numero de particulas
      ls=sqrt(1.d0*npar/rho)
      write(*,*) ls
c tamaño de la caja cuadrada
      hueco=ls/dfloat(sqn)
      write(*,*) 'hueco:', hueco
      write(*,*) 'cutoff:', cutoff
c hueco entre particulas en la configuracion reticular

c genera red con posiciones iniciales
c ya distorsionada

      i=1
      do j=1,sqn
        do k=1,sqn
	  call random_number(rseed)
          x(i)=(dfloat(j)-0.5)*ls/dfloat(sqn)+5.0d-1*hueco*(rseed-0.5)
          y(i)=(dfloat(k)-0.5)*ls/dfloat(sqn)+5.0d-1*hueco*(rseed-0.5)
          i=i+1
        enddo
      enddo

c      npar = 3
c      ls=sqrt(1.d0*npar/rho)
c      write(*,*) ls
c      x(1) = 1.0d0
c      x(2) = 1.0d0 + hueco
c      x(3) = 1.0d0 + hueco/2.0d0
c      y(1) = 1.0d0
c      y(2) = 1.0d0 
c      y(3) = 1.0d0 + 2.0d0*hueco/3.0d0


      tt=10000
      write(numf,'(i5)') tt
      open(12,file='pos'//numf//'.dat')
      do i=1,npar
        write(12,*) x(i),y(i)
      enddo
      close(12)
c hasta aquí va bien

      tt=tt+1
      do i=1,npar
        xn(i)=x(i)
        yn(i)=y(i)
      enddo

c bucle de tiempo principal
      do it=1,nit

c inicializa matriz y vector de fuerzas y posiciones nuevas
        do i=1,npar
          do j=1,npar
            fx(i,j)=0.0d0
            fy(i,j)=0.0d0
          enddo
          fxj(i)=0.0d0
          fyj(i)=0.0d0
          xn(i)=x(i)
          yn(i)=y(i)
        enddo

c primer bucle particular para rellenar matriz de fuerzas
        do i=1,npar
c	  if ((i.gt.190).and.(i.le.210)) cycle
          do j=1,npar
c	    if ((j.gt.190).and.(j.le.210)) cycle
	    if (i==j) cycle
	    dx=x(i)-x(j)
c displacement in the x-axis (can be >0 <0)
            dx=dx-ls*nint(dx/ls)
            dy=y(i)-y(j)
            dy=dy-ls*nint(dy/ls)
            dist=dsqrt(dx*dx+dy*dy)
c xy plane distance
	    if ((i.gt.195).and.(i.le.205)) then
	      if (dist.lt.cutoff) then
c calculo de la fuerza repulsiva
	        fmodulo=fuerzaREP(dist,f0,sigma)
		fx(i,j)=fmodulo*dx/dist
		fxj(i)=fxj(i)+fx(i,j)
		fy(i,j)=fmodulo*dy/dist
                fyj(i)=fyj(i)+fy(i,j)
		if (fmodulo.gt.1000.0d0) then
		  write(*,*) 'Fuerza muy alta en la particula :', i,j,tt
		endif
	      endif
	    else if ((j.gt.195).and.(j.le.205)) then
	      if (dist.lt.cutoff) then
c calculo de la fuerza repulsiva
                fmodulo=fuerzaREP(dist,f0,sigma)
                fx(i,j)=fmodulo*dx/dist
                fxj(i)=fxj(i)+fx(i,j)
                fy(i,j)=fmodulo*dy/dist
                fyj(i)=fyj(i)+fy(i,j)
		if (fmodulo.gt.1000.0d0) then
                  write(*,*) 'Fuerza muy alta en la particula :', i,j,it
                endif
              endif
	    else
              if (dist.lt.cutoff) then
c calculo de la fuerza
                fmodulo=fuerza(dist,f0,sigma)
                fx(i,j)=fmodulo*dx/dist
                fxj(i)=fxj(i)+fx(i,j)
                fy(i,j)=fmodulo*dy/dist
                fyj(i)=fyj(i)+fy(i,j)
              endif
	     endif
c            write(*,*) 'cosas'
c            write(*,*) i,j,fx(i,j),fy(i,j),fxj(i),fyj(i)
          enddo
        enddo

c segundo bucle para actualizar posición de la partícula
        do i=1,npar
          do j=1,npar
            if (i==j) cycle
c la posición se actualiza por 3 contribuciones:
c fuerza que le hacen las vecinas
c flujo debido a rotación de vecinas
c flujo debido a fuerza que siente de las vecinas
            dx=x(i)-x(j)
            dx=dx-ls*nint(dx/ls)
            dy=y(i)-y(j)
            dy=dy-ls*nint(dy/ls)
            dist=dsqrt(dx*dx+dy*dy)
	    if ((j.gt.195).and.(j.le.205)) then
              displx=(fx(i,j)/(3.0d0*visc)+0*vrx(omega,dy,dist,h)+
     +              vix(fxj(j),fyj(j),dx,dy,dist,visc,h))*dt
              disply=(fy(i,j)/(3.0d0*visc)+0*vry(omega,dx,dist,h)+
     +              viy(fxj(j),fyj(j),dx,dy,dist,visc,h))*dt
	    else
	      displx=(fx(i,j)/(3.0d0*visc)+vrx(omega,dy,dist,h)+
     +              vix(fxj(j),fyj(j),dx,dy,dist,visc,h))*dt
              disply=(fy(i,j)/(3.0d0*visc)+vry(omega,dx,dist,h)+
     +              viy(fxj(j),fyj(j),dx,dy,dist,visc,h))*dt
	    endif      
            xn(i)=xn(i)+displx
c      write(*,*) i,j,x(i),xn(i),fx(i,j)
            yn(i)=yn(i)+disply
c      write(*,*) i,j,y(i),yn(i),fy(i,j)
c de mi experiencia con el anterior código vi que era mejor
c aplicar tras cada paso de tiempo las ccppc
      if (displx.gt.0.5) then
        write(*,*) 'error :'
	write(*,*) i,j 
      write(*,*) displx
      write(*,*) (fx(i,j)/(3.0d0*visc))*dt 
      write(*,*) (vrx(omega,dy,dist,h))*dt 
      write(*,*) (vix(fxj(j),fyj(j),dx,dy,dist,visc,h))*dt 
c      stop
    
      endif
c      write(*,*) displx, disply
c   PBC
            if (xn(i).gt.ls) then
              xn(i)=xn(i)-ls
c              write(*,*) 'X-'
            elseif (xn(i).lt.0.0d0) then
              xn(i)=xn(i)+ls
c              write(*,*) 'X+'
            endif
           
            if (yn(i).gt.ls) then
              yn(i)=yn(i)-ls
c              write(*,*) 'Y-'
            elseif (yn(i).lt.0.0d0) then
              yn(i)=yn(i)+ls
c              write(*,*) 'Y+'
            endif
            if (xn(i).gt.ls) then
	      write(*,*) 'error'
	      stop
	    else if (yn(i).gt.ls) then
              write(*,*) 'error'
              stop
	    endif

          enddo
        enddo
c calculo velocidad instantanea
        do i=1,npar
c          write(*,*) xn(i), x(i)
          vx(i)=(xn(i)-x(i))/dt
          vy(i)=(yn(i)-y(i))/dt
c          write(*,*) vx(i), vy(i)
        enddo

c calculo flujo de velocidad
        if (mod(it,step).eq.0) then
          do i=1,200
            do j=1,200
              ux(i,j)=0.0d0
              uy(i,j)=0.0d0
              do k=1,npar
                dx=dfloat(i)*ls*1.0d-2-x(k)
                dx=dx-ls*nint(dx/ls)
                dy=dfloat(j)*ls*1.0d-2-y(k)
                dy=dy-ls*nint(dy/ls)
                dist=dsqrt(dx*dx+dy*dy)
        if(dist.gt.sigma/5) then
                ux(i,j)=ux(i,j)+vrx(omega,dy,dist,h)+
     +                  vix(fxj(k),fyj(k),dx,dy,dist,visc,h)
                uy(i,j)=uy(i,j)+vry(omega,dx,dist,h)+
     +                  viy(fxj(k),fyj(k),dx,dy,dist,visc,h)
        endif
        
              enddo
            enddo
          enddo
        endif      

c actualizo las posiciones fuera del bucle principal
c si hubiera guardado las distancias en una matriz

        do i=1,npar
          x(i)=xn(i)
          y(i)=yn(i)
        enddo

c escritura de resultados
        if (mod(it,step).eq.0) then
          write(numf,'(i5)')tt
          open(12,file='conf'//numf//'.dat')
          do i=1,npar
            write(12,*) x(i),y(i),vx(i),vy(i),fx(i,49),fxj(i) 
          enddo
          close(12)
c          open(13,file='ufield'//numf//'.dat')
c          do i=1,200
c        do j=1,200
c                    write(13,*) i,j,ux(i,j), uy(i,j)
c        enddo
c      enddo
c          close(13)
      tt=tt+1
        endif

      enddo

      stop
      end

c Funciones y subrutinas

c fuerza
c arg1: dist
c arg2: epsilon
c arg3: sigma

c calculo de la fuerza que deriva del potencial de interaccion
      real(8) function fuerza(arg1, arg2, arg3) result(res) ! arg1, arg2, arg3 = dist, f0, sigma
      implicit none
      real(8), intent(in) :: arg1, arg2, arg3
      real(8) :: r6, ir6, ir12, s6, s12

      if (arg1 == 0.0d0) then
        res = 0.0d0  ! Evitem divisió per zero
      return
      endif

      r6 = arg1**6
      ir6 = 1.0d0 / r6
      ir12 = ir6 * ir6
      s6 = arg3**6
      s12 = s6 * s6

c      res = 48.0d0 * arg2 * (s12 * ir12 - 0.5d0 * s6 * ir6) / arg1

      res = 48.d0 * arg2 * (s12 / arg1**13 - 0.5d0 * s6 / arg1**7)

      end function fuerza

c calculo de la fuerza que deriva del potencial de interaccion
      real(8) function fuerzaREP(arg1, arg2, arg3) result(res) ! arg1, arg2, arg3 = dist, f0, sigma
      implicit none
      real(8), intent(in) :: arg1, arg2, arg3
      real(8) :: r6, ir6, ir12, s6, s12, cut

      if (arg1 == 0.0d0) then
        res = 0.0d0  ! Evitem divisió per zero
      return
      endif

      r6 = arg1**6
      ir6 = 1.0d0 / r6
      ir12 = ir6 * ir6
      s6 = arg3**6
      s12 = s6 * s6
      
      cut = arg3*2.0d0**(1.0d0/6.0d0)
      if (arg1.le.cut) then
	      res = 48.d0 * arg2 * (s12 / arg1**13 - 0.5d0 * s6 / arg1**7)
      else
	      res = 0
      endif
      end function fuerzaREP

c vrx
c arg1: omega
c arg2: dy
c arg3: dist
c     ROTLET
       real(8) function vrx(arg1,arg2,arg3,h8) result(res)
c call: vrx(omega,dy,dist,h)
      implicit none
       real(8) ::  arg1,arg2,arg3,h8,r2

      r2=arg3*arg3+4.0d0*h8*h8
c  r2=dist^2+4*h^2
c 4???
      res=arg1*arg2*(1.0d0/dsqrt(r2*r2*r2)-1.0d0/(arg3*arg3*arg3))
c   vrx=omega*dy*(1/R^3-1/r^3)
      
      end function vrx

c vry
c arg1: omega
c arg2: dx
c arg3: dist
  
      real(8) function vry(arg1,arg2,arg3,h8) result(res)
      implicit none
       real(8) ::  arg1,arg2,arg3,h8,r2

      r2=arg3*arg3+4.0d0*h8*h8
      res=arg1*arg2*(1.0d0/(arg3*arg3*arg3)-1.0d0/dsqrt(r2*r2*r2))

      end function vry

c vix
c arg1: fx(j)
c arg2: fy(j)
c arg3: xi-xj
! arg4: yi-yj
! arg5: dist
! arg6: visc

!: STOKESLET

      real(8) function vix(arg1,arg2,arg3,arg4,arg5,arg6,h8) result(res)
 !   call: vix(fxj(k),fyj(k),dx,dy,dist,visc,h)
      implicit none
       real(8) ::  arg1,arg2,arg3,arg4,arg5,arg6,h8,r2

      r2=arg5*arg5+4.0d0*h8*h8
!   r2=dist^2+4*h^2
! 4???
      res=(arg1*(arg5*arg5+arg3*arg3)+arg2*arg3*arg4)/
     $    (4.0d0*arg6*arg5*arg5*arg5)+
! (Fx*(dist^2+dx^2)+Fy*dx*dy)/(4*visc*dist^3)
! equivalent to Fx*(1/r + r_x^2/r^3)/4*visc + Fy*r_x*r_y/(4*visc*r^3)
! pared
     $    (arg1*(2.0d0*h8*h8*(r2-3.0d0*arg3*arg3)-r2*(r2+arg3*arg3))+
!   Fx/(4*visc)*(2h*(h(R^2-R_x^2))/R^5)-Fx/(4*visc)*(1/R+R_x^2/R^3)
     $   arg2*6.0d0*h8*h8*arg3*arg4)/(4.0d0*arg6*dsqrt(r2*r2*r2*r2*r2))
!   Fy/(4*visc)*(6h^2*R_x*R_y)/R^5
! Falta R_x*R_y/R^3???

      end function vix

! viy
! arg1: fx(j)
! arg2: fy(j)
! arg3: xi-xj
! arg4: yi-yj
! arg5: dist
! arg6: visc

      real(8) function viy(arg1,arg2,arg3,arg4,arg5,arg6,h8) result(res)
!      real(8) ::function viy(arg1,arg2,arg3,arg4,arg5,arg6,h8)
      implicit none
      real(8) ::  arg1,arg2,arg3,arg4,arg5,arg6,h8,r2

      r2=arg5*arg5+4.0d0*h8*h8
      res=(arg2*(arg5*arg5+arg4*arg4)+arg1*arg3*arg4)/
     $    (4.0d0*arg6*arg5*arg5*arg5)+
! pared
     $    (arg2*(2.0d0*h8*h8*(r2-3.0d0*arg4*arg4)-r2*(r2+arg4*arg4))+
     $   arg1*6.0d0*h8*h8*arg3*arg4)/(4.0d0*arg6*dsqrt(r2*r2*r2*r2*r2))

      end function viy

