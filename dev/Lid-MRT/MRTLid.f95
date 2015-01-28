Program MRTCavity
!Convalidação do Metodo
parameter (n=100,m=100)
real f(0:8,0:n,0:m),feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8),cx(0:8),cy(0:8),u(0:n,0:m),v(0:n,0:m)
real tminv(0:8,0:8),Sm(0:8),tm(0:8,0:8),Stm(0:8,0:8),ev(0:8,0:8)
real mom(0:8,0:n,0:m),meq(0:8,0:n,0:m)
real sumb
integer l,i,j,k
open(2,file='resultCavity1')
uo=0.3                                           !10
rhoo=5.0              
dx=1.0
dy=1.0
dt=1.0
visco=0.01
Re=uo*m/visco
print*,"Re=",Re
omega=1.0/(3.0*visco+0.5)
mstep=4000                                       
w(0)=4./9.										!20
do i=1,4
  w(i)=1./9.
end do
do i=5,8
  w(i)=1./36.
end do
cx(:)=(/0.,1.,-1.,0.,0.,1.,-1.,-1.,1./)
cy(:)=(/0.,0.,0.,1.,-1.,1.,-1.,1.,-1./)
tm(0,:)=(/1.,1.,1.,1.,1.,1.,1.,1.,1./)
tm(1,:)=(/-4.,-1.,-1.,-1.,-1.,2.,2.,2.,2./)			!30
tm(2,:)=(/4.,-2.,-2.,-2.,-2.,1.,1.,1.,1./)
tm(3,:)=(/0.,1.,0.,-1.,0.,1.,-1.,-1.,1./)
tm(4,:)=(/0.,-2.,0.,2.,0.,1.,-1.,-1.,1./)
tm(5,:)=(/0.,0.,1.,0.,-1.,1.,1.,-1.,-1./)
tm(6,:)=(/0.,0.,-2.,0.,2.,1.,1.,-1.,-1./)
tm(7,:)=(/0.,1.,-1.,1.,-1.,0.,0.,0.,0./)
tm(8,:)=(/0.,0.,0.,0.,0.,1.,-1.,1.,-1./)
al=1./1.
tminv(0,:)=(/4.*al,-4.*al,4.*al,0.,0.,0.,0.,0.,0./)
tminv(1,:)=(/4.*al,-al,-2.*al,6.*al,-6.*al,0.,0.,9.*al,0./)				!40
tminv(2,:)=(/4.*al,-al,-2.*al,0.,0.,6.*al,-6.*al,-9.*al,0./)
tminv(3,:)=(/4.*al,-al,-2.*al,-6*al,6.*al,0.,0.,9.*al,0./)
tminv(4,:)=(/4.*al,-al,-2.*al,0.,0.,-6.*al,6.*al,-9.*al,0./)
tminv(5,:)=(/4.*al,2.*al,al,6.*al,3.*al,6.*al,3.*al,0.,9.*al/)
tminv(6,:)=(/4.*al,2.*al,al,-6.*al,-3.*al,6.*al,3.*al,0.,-9.*al/)
tminv(7,:)=(/4.*al,2.*al,al,-6.*al,-3.*al,-6.*al,-3.*al,0.,9.*al/)
tminv(8,:)=(/4.*al,2.*al,al,6.*al,3.*al,-6.*al,-3.*al,0.,-9.*al/)
do i=0,8
  do j=0,8
    sumcc=0.0																!50
    do l=0,8
      sumcc=sumcc+tminv(i,l)*tm(l,j)
    end do
    ev(i,j)=sumcc*1./36.
  end do
end do
tau=1./omega

Sm(:)=(/1.,1.4,1.4,1.,1.2,1.,1.2,tau,tau/)
do i=0,8																	!60
  do j=0,8
    Stm(i,j)=tminv(i,j)*Sm(j)*ev(i,j)
  end do
end do

do j=0,m
  do i=0,n
    rho(i,j)=rhoo
    u(i,j)=0.0                                      
    v(i,j)=0.0																!70
  end do
end do
do j=0,m
  do i=0,n
    do k=0,8
      f(k,i,j)=0.0
    end do
  end do
end do                                               
do i=1,n-1																	!80
  u(i,m)=uo
  
end do

!Loop Principal
do kk=1,mstep 
  do i=0,n
    do j=0,m         !Calculo do momento de equilibrio
      meq(0,i,j)=rho(i,j)                                      
      meq(1,i,j)=rho(i,j)*(-2.+3.*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))			!90
      meq(2,i,j)=rho(i,j)*(1.+3.*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
      meq(3,i,j)=rho(i,j)*u(i,j)
      meq(4,i,j)=-rho(i,j)*u(i,j)
      meq(5,i,j)=rho(i,j)*v(i,j)
      meq(6,i,j)=-rho(i,j)*v(i,j)
      meq(7,i,j)=rho(i,j)**(2)*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
      meq(6,i,j)=rho(i,j)*u(i,j)*v(i,j)
    end do
  end do
  do i=0,n  !Calculo do momento																!100
    do j=0,m
      do k=0,8
        suma=0.0
        do l=0,8
          suma=suma+tm(k,l)*f(l,i,j)
        end do
        mom(k,i,j)=suma
      end do
    end do
  end do																				!110
  do i=0,n    !Colisão no espaço de momento
    do j=0,m
      do k=0,8        
        sumb=0.0
        do l=0,8
          sumb=sumb+Stm(k,l)*(mom(l,i,j)-meq(l,i,j))
        end do
        f(k,i,j)=f(k,i,j)-sumb
      end do
    end do																					!120
  end do                                            
  do j=0,m                     !Processo de Transmissão
    do i=1,n                                                       
      f(1,n-i+1,j)=f(1,n-i,j)        !Left to Right
      f(2,i-1,j)=f(2,i,j)            !Right to Left
    end do                          
  end do
  do i=0,n
    do j=1,m
      f(3,i,m-j+1)=f(3,i,m-j)        !Bottom to Top
      f(4,i,j-1)=f(4,i,j)            !Top to Bottom         
    end do
  end do                                                          !90
  do j=1,m
    do i=1,n
      f(5,n-i+1,m-j+1)=f(5,n-i,m-j)   !Diagonal ++     
      f(6,i-1,j-1)=f(6,i,j)           !Diagonal --
    end do
  end do
  do j=1,m
    do i=1,n                                               
      f(7,i-1,m-j+1)=f(7,i,m-j)       !Diagonal -+
      f(8,n-i+1,j-1)=f(8,n-i,j)       !Diagonal +-                   !100
    end do
  end do
  do  j=0,m                  !Condições de Contorno
    f(1,0,j)=f(2,0,j)        !Metodo Bounce Back na Camada Oeste
    f(8,0,j)=f(7,0,j)
    f(5,0,j)=f(6,0,j)
    f(2,n,j)=f(1,n,j)        !Metodo Bounce Back na Camada Leste     
    f(7,n,j)=f(8,n,j)
    f(6,n,j)=f(5,n,j)                                                                                      
  end do                                                                    !110
  do i=0,n
    f(3,i,0)=f(4,i,0)        !Metodo Bounce Back na Camada Sul
    f(5,i,0)=f(6,i,0)
    f(7,i,0)=f(8,i,0)
  end do
  do i=1,n-1                     !Metodo Bounce Back na Camada Norte
    rhon=f(0,i,m)+f(1,i,m)+f(2,i,m)+2.*(f(3,i,m)+f(7,i,m)+f(5,i,m))          
    f(4,i,m)=f(3,i,m)
    f(8,i,m)=f(7,i,m)+(rhon*uo/6.0)                                   
    f(6,i,m)=f(5,i,m)-(rhon*uo/6.0)                                                !120
  end do
  do j=0,m                  !Montando resultados
    do i=0,n
      ssum=0.0
      do k=0,8
        ssum=ssum+f(k,i,j)
      end do                                                
      rho(i,j)=ssum
    end do                                                         
  end do
  do i=0,n
    do j=0,m-1
      usum=0.0
      vsum=0.0
      do k=0,8
        usum=usum+f(k,i,j)*cx(k)
        vsum=vsum+f(k,i,j)*cy(k)                                            
      end do
      u(i,j)=usum/rho(i,j)
      v(i,j)=vsum/rho(i,j)
    end do
  end do
end do
!Plotando resultados
do j=0,m
  write(2,*)(u(i,j),i=0,n)
end do
do j=0,m
  write(2,*)(v(i,j),i=0,n)
end do 

                                                                
stop
end program MRTCavity