Program lidDrivOMP
!Convalidação do Metodo

parameter (n=100,m=100)
real f(0:8,0:n,0:m),feq(0:8,0:n,0:m),rho(0:n,0:m)
real w(0:8),cx(0:8),cy(0:8),u(0:n,0:m),v(0:n,0:m)
integer i,j,k
open(2,file='lidDrivOMP.out')

uo=0.1                                           !10
!sumvelo=0.0      !não seve pra nada
rhoo=5.0              !pq?
dx=1.0
dy=1.0
dt=1.0
visco=0.01
Re=uo*m/visco
print*,"Re=",Re
omega=1.0/(3.0*visco+0.5)
mstep=40000                                       !20
w(0)=4./9.
do i=1,4
  w(i)=1./9.
end do
do i=5,8
  w(i)=1./36.
end do
cx(0)=0
cx(1)=1
cx(2)=-1                                             !30
cx(3)=0
cx(4)=0
cx(5)=1
cx(6)=-1
cx(7)=-1
cx(8)=1
cy(0)=0
cy(1)=0
cy(2)=0
cy(3)=1                                           !40
cy(4)=-1
cy(5)=1
cy(6)=-1
cy(7)=1
cy(8)=-1

do j=0,m
  do i=0,n
    rho(i,j)=rhoo
    u(i,j)=0.0                                      !50
    v(i,j)=0.0
  end do
end do
do j=0,m
  do i=0,n
    do k=0,8
      f(k,i,j)=0.0
    end do
  end do
end do                                               !60
do i=1,n-1
  u(i,m)=uo
  v(i,m)=0.0
end do

!Loop Principal
do kk=1,mstep                                       
  do i=0,n                     !Processo de Colisão
    do j=0,m
      t1=u(i,j)*u(i,j)+v(i,j)*v(i,j)                           !70
      do k=0,8
        t2=u(i,j)*cx(k)+v(i,j)*cy(k)
        feq(k,i,j)=rho(i,j)*w(k)*(1.0+3.0*t2+4.5*t2*t2-1.5*t1)
        !if(k.eq.0) feq(k,i,j)=w(k)*rho(i,j)
        f(k,i,j)=omega*feq(k,i,j)+(1.0-omega)*f(k,i,j)
      end do
    end do
  end do                                            
  do j=0,m                     !Processo de Transmissão
    do i=1,n                                                       !80
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
end program lidDrivOMP

